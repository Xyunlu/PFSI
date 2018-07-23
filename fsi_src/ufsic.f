      subroutine ufsic(iorder,iorder1,iorder0,
     & coor,bf,nodvar,uvar,
     & eu1,ur,ulap,ulapn,
     & vlap,vlapn,coor0,coorlap,
     & ua,u1a,dua,ums,
     & vms,wms,vs,acs,
     & ac,eu,apool,knode,kdgof,
     & kvar,knodeg,kdgofg,kvarg,knode1,kdgofb,neq,kcoor,
     & time,dt,it,itn_tol,cc,nsource,iblk,msend,
     & msend_tol,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      dimension nodvar(knode,kdgof),bf(knode,kdgof),coor(knode,kcoor),
     & eu(knode,kdgof),eu1(knode,kdgof),
     * uvar(neq)
      dimension apool(kvarg)
      dimension iorder(knode),iorder1(knode1),iorder0(knodeg),
     +          ur(knodeg,kdgof),
     +          ulap(knodeg,kdgof),ulapn(knodeg,kdgof),
     +          vlap(knodeg,kdgof),vlapn(knodeg,kdgof),
     +          coorlap(knodeg,kcoor),coor0(knodeg,kcoor),
     +          ua(knodeg,kdgofg),u1a(knodeg,kdgofg),dua(knodeg,kdgofg),
     +          ums(knode1),vms(knode1),wms(knode1),
     +          vs(knode1,kdgofb),acs(knode1,kdgofb),
     +          ac(knode1,kdgofb)

      dimension mnode(100),nnode(100),nmate(100),mmate(100)

      open(88,file='control.dat',status='old')
      read(88,*) emax,numstep,emax,numstep
      close(88)
6     format (1x,10i5)
7     format (1x,6e12.5)

c      call recvai(iblk,nsource,iorder,knode)
      nextern=0
      nupdate=0
      do inod=1,knode
        if(iorder(inod).eq.0) nextern=nextern+1
        if(iorder(inod).eq.1) nupdate=nupdate+1
      enddo
c      call recvai(iblk,nsource,iorder,knode)
      nextern1=0
      nupdate1=0
      do inod=1,knode1
        if(iorder1(inod).eq.0) nextern1=nextern1+1
        if(iorder1(inod).eq.1) nupdate1=nupdate1+1
      enddo
      nextern0=0
      nupdate0=0
      do inod=1,knodeg
        if(iorder0(inod).eq.0) nextern0=nextern0+1
        if(iorder0(inod).eq.1) nupdate0=nupdate0+1
      enddo

C.... OPEN time file
C.... Open coor file
c.....read bfd file
c      read (22) ((eu(i,j),i=1,kdgof),j=1,knode)
      do inod=1,knode
       do idof=1,kdgof
       eu(inod,idof)=bf(inod,idof)
       enddo
      enddo
cc      write (*,*) 'u ='
cc      write(*,7) ((eu(i,j),i=1,kdgof),j=1,knode)
c.....read nv file
c      read (23) nskp
c      read (23) ((nodvar(i,j),i=1,kdgof),j=1,knode)
cc      write (*,*) 'nodvar ='
cc      write (*,6) ((nodvar(i,j),i=1,kdgof),j=1,knode)
c
c       Lsp have modified here, for convenient, old code left here
c
c.....read u file
c      read (24) (uvar(i),i=1,neq)
c      write (*,*) 'uvar ='
c      write (*,7) (uvar(i),i=1,neq)
c      do 200 inod=1,knode
c       do 100 idfg=1,kdgof
c        n=nodvar(idfg,inod)
c        if (n.le.0) goto 100
c        eu(idfg,inod)=uvar(n)
c100     continue
c200     continue
c
c.....read u file
c      read (24) (uvar(i),i=1,neq)
c      write (*,*) 'uvar ='
c      write (*,7) (uvar(i),i=1,neq)
        k = 0
      do inod=1,knode
       do idof=1,kdgof
        n=nodvar(inod,idof)
        if (n.gt.0) then
          k = k + 1
         eu(inod,idof)=uvar(n)
        endif
       enddo
      enddo
c
c      do n=1,knode
c...... write displace file
c        write (25,1600) n,(eu(i,n),i=1,kdgof)
c      enddo
c      write(25,*)'**********************************************'

CC
CC..  update the solution ...
      ntemp = 0
      do i=1,knode
        if(iorder(i).eq.1) then
          do j=1,kdgof
            ntemp=ntemp+1
            apool(ntemp)=eu(i,j)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
      if(ntemp.ne.nupdate*kdgof) then
        print *,'Error! ntemp ne nupdate , U program !!',ntemp,nupdate
        call endjob(ierr)
      endif
CC..
      call recvar(iblk,nsource,apool,nextern*kdgof)
      ntemp=0
      do i=1,knode
        if(iorder(i).eq.0) then
          do j=1,kdgof
            ntemp=ntemp+1
            eu(i,j)=apool(ntemp)
          enddo
        endif
      enddo
      if(ntemp.ne.nextern*kdgof) then
        print *,'Error! ntemp ne nextern , U program !',ntemp, nextern
        call endjob(ierr)
      endif

      MSend_tol = 0
      ntemp = 0
      do i=1,knode
        if(iorder(i).eq.1) then
          do j=1,kdgof
            ntemp=ntemp+1
            apool(ntemp)=eu1(i,j)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
      if(ntemp.ne.nupdate*kdgof) then
        print *,'Error! Send eu1 to master In Ulaplace !!'
        print *,'Error! ntemp ne nupdate :',ntemp,nupdate
        call endjob(ierr)
      endif
C...
C==   Read id0, subdomain.idx, coor00 in master processor
C..   read id0
C..   read subdomain.idx
C..   read coor00
C=======================================================
      ntemp = 0
      do i=1,knodeg
        if(iorder0(i).eq.1) then
          do j=1,kdgof
            ntemp=ntemp+1
            apool(ntemp)=ur(i,j)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
      if(ntemp.ne.nupdate0*kdgof) then
        print *,'Error! Send ur to master In Ulaplace !!'
        print *,'Error! ntemp ne nupdate0 :',ntemp,nupdate0
        call endjob(ierr)
      endif
C...  Send ulap to master processor
      ntemp = 0
      do i=1,knodeg
        if(iorder0(i).eq.1) then
          do j=1,kdgof
            ntemp=ntemp+1
            apool(ntemp)=ulap(i,j)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
      if(ntemp.ne.nupdate0*kdgof) then
        print *,'Error! Send ulap to master In Ulaplace !!'
        print *,'Error! ntemp ne nupdate0 :',ntemp,nupdate0
        call endjob(ierr)
      endif
C...  Send ulapn to master processor
      ntemp = 0
      do i=1,knodeg
        if(iorder0(i).eq.1) then
          do j=1,kdgof
            ntemp=ntemp+1
            apool(ntemp)=ulapn(i,j)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
      if(ntemp.ne.nupdate0*kdgof) then
        print *,'Error! Send ulapn to master In Ulaplace !!'
        print *,'Error! ntemp ne nupdate0 :',ntemp,nupdate0
        call endjob(ierr)
      endif
C...
C==  Read vangle file in master processor
C...
C==  Compute coorlap and update coor and ulap in master processor
C=== Update unodd file
C..  Recv ulap from master processor
      call recvar(iblk,nsource,apool,knodeg*kdgof)
      ntemp = 0
      do i=1,knodeg
        do j=1,kdgof
          ntemp=ntemp+1
          ulap(i,j)=apool(ntemp)
        enddo
      enddo
C...  Recv vlap from master processor
      call recvar(iblk,nsource,apool,knodeg*kdgof)
      ntemp = 0
      do i=1,knodeg
        do j=1,kdgof
          ntemp=ntemp+1
          vlap(i,j)=apool(ntemp)
        enddo
      enddo

C..  Recv coor from master processor to update coor0 file
      call recvar(iblk,nsource,apool,knodeg*kcoor)
      ntemp = 0
      do i=1,knodeg
        do j=1,kdgof
          ntemp=ntemp+1
          coor0(i,j)=apool(ntemp)
        enddo
      enddo

C..  Recv coorlap from master processor to update cooreule file
      call recvar(iblk,nsource,apool,knodeg*kcoor)
      ntemp = 0
      do i=1,knodeg
        do j=1,kdgof
          ntemp=ntemp+1
          coorlap(i,j)=apool(ntemp)
        enddo
      enddo

      err = 0.0
      sum = 0.0
      do 500 inod=1,knode
      do 501 idof=1,kdgof
      if(iorder(inod) .eq. 1) then
        err = err+(eu(inod,idof)-eu1(inod,idof))**2
        sum = sum+eu1(inod,idof)**2
      endif
501   continue
500   continue
      ione = itn_tol
      call sendint(nsource,iblk,ione)
      call sendr(nsource,iblk,err)
      call sendr(nsource,iblk,sum)

      call recvint(iblk,nsource,MSend_tol)

      if (MSend_tol .eq. 1) then
C...    Update unodd file
        do i=1,knodeg
          do j=1,kdgof
            ulapn(i,j) = ulap(i,j)
            vlapn(i,j) = vlap(i,j)
          enddo
        enddo
C....   Update unoda file
        do i=1,knodeg
          do j=1,kdgofg
            u1a(i,j) = ua(i,j)
            dua(i,j) = 0.D0
          enddo
        enddo
C....  Read accel from unods to Update unodac file
        do i=1,knode1
          do j=1,kdgof
            ac(i,j) = acs(i,j)
          enddo
        enddo
C...  Send ac to master processor
c        ntemp = 0
c        do i=1,knode1
c          if(iorder1(i).eq.1) then
c            do j=1,kdgofb
c              ntemp=ntemp+1
c              apool(ntemp)=ac(i,j)
c            enddo
c          endif
c        enddo
c        call sendar(nsource,iblk,apool,ntemp)
c        if(ntemp.ne.nupdate1*kdgof) then
c          print *,'Error! Send ac to master In Ulaplace !!'
c          print *,'Error! ntemp ne nupdate1 :',ntemp,nupdate1
c          call endjob(ierr)
c        endif
      else
        ione=ione+1
        itn_tol = itn_tol+1
      endif
c      call sswitch(nsource,iblk,MSend,0)
      call sendar (nsource,iblk,eu,knode*kdgof)
      do 502 i=1,knode
      do 502 j=1,kdgof
      eu1(i,j)=eu(i,j)
502   continue
c@send

1600  format (1x,i4,6e12.5)
      return
      end
