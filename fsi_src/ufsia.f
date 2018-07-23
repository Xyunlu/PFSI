      subroutine ufsia(iorder,iorder1,iorder2,
     & coor,bf,nodvar,uvar,
     & eu,eu1,edu,coor0,
     & eum,evm,ewm,eumn,
     & evmn,ewmn,evlx,evly,
     & evlz,evlxn,evlyn,evlzn,
     & ubc,accel,ums,vms,
     & wms,vs,acs,evf,
     & ev,eue,apool,knode,knode1,
     & knode2,kdgof,kdgofb,kdgofc,kvar,neq,kcoor,time,
     & dt,it,itn,cc,nsource,iblk,msend,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      dimension nodvar(knode,kdgof),bf(knode,kdgof),coor(knode,kcoor),
     & coor0(knode,kcoor),
     & eu1(knode,kdgof),eu(knode,kdgof),ev(knode,kdgof),
     & eue(knode,kdgof),edu(knode,kdgof),evf(knode,kdgof),
     & eum(knode),evm(knode),ewm(knode),
     & eumn(knode),evmn(knode),ewmn(knode),
     & evlx(knode),evly(knode),evlz(knode),
     & evlxn(knode),evlyn(knode),evlzn(knode),
     * uvar(neq)
      dimension iorder(knode),apool(kvar)
      dimension iorder1(knode1),iorder2(knode2)
      dimension accel(knode1,kdgofb),ubc(knode2,kdgofc)
      dimension ums(knode1),vms(knode1),wms(knode1),
     & vs(knode1,kdgofb),acs(knode1,kdgofb)

      dimension mnode(100),nnode(100),nmate(100),mmate(100)

      logical filflg1, filflg
      common /rdata/ vangle,angle,penalty,pen4sd
      common /idata/ id4updown
c......
      open (41,file='control.dat',status='old')
      read (41,*) emax,numstep,skip,nskip,relax
      read (41,*) skip
      read (41,*) skip
      read (41,*) pen4p,pen4sd
      close(41)
      close(41)
c
      open (41,file='init.dat',status='old')
      read (41,*) xmax,ymax,zmax,xo,yo,zo,radius,thick,radaxis,angle,
     +            radart_mult,ratiol
      read (41,*) nx,nz,ns,npt,nct,nsc,naxis,nyi,nyo
      read (41,*) ai,bi,ao,bo
      read (41,*) u_in,visc_f,pe,pv,denf,dens,gravity,omega
      read (41,*) idpoint,timespring
      read (41,*) skip
      read (41,*) skip
      read (41,*) id4updown
      close(41)
c
c      open (41,file='vangle',form='formatted',status='old')
c      read (41,*) vangle,vangle1,angle,angle1
c      close(41)
6     format (1x,10i5)
7     format (1x,6e12.5)

c      call recvai(iblk,nsource,iorder,knode)
      nextern=0
      nupdate=0
      do inod=1,knode
        if(iorder(inod).eq.0) nextern=nextern+1
        if(iorder(inod).eq.1) nupdate=nupdate+1
      enddo

      nextern1=0
      nupdate1=0
      do inod=1,knode1
        if(iorder1(inod).eq.0) nextern1=nextern1+1
        if(iorder1(inod).eq.1) nupdate1=nupdate1+1
      enddo

      nextern2=0
      nupdate2=0
      do inod=1,knode2
        if(iorder2(inod).eq.0) nextern2=nextern2+1
        if(iorder2(inod).eq.1) nupdate2=nupdate2+1
      enddo
C.... OPEN time file
C.... Open coor file
c.....read bfd file
c      read (22) ((ev(i,j),i=1,kdgof),j=1,knode)
      do inod=1,knode
       do idof=1,kdgof
       ev(inod,idof)=bf(inod,idof)
       enddo
      enddo
cc      write (*,*) 'u ='
cc      write(*,7) ((ev(i,j),i=1,kdgof),j=1,knode)
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
c        ev(idfg,inod)=uvar(n)
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
         ev(inod,idof)=uvar(n)
        endif
       enddo
      enddo
c
c      do n=1,knode
c...... write displace file
c        write (25,1600) n,(ev(i,n),i=1,kdgof)
c      enddo
c      write(25,*)'**********************************************'

CC
CC..  update the solution ...
      ntemp = 0
      do i=1,knode
        if(iorder(i).eq.1) then
          do j=1,kdgof
            ntemp=ntemp+1
            apool(ntemp)=ev(i,j)
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
            ev(i,j)=apool(ntemp)
          enddo
        endif
      enddo
      if(ntemp.ne.nextern*kdgof) then
        print *,'Error! ntemp ne nextern , U program !',ntemp, nextern
        call endjob(ierr)
      endif

      MSend = 0
      aa = 0.0
      ab = 0.0
      bb = 0.0
      do 500 inod=1,knode
      do 501 idof=1,kdgof
        eue(inod,idof) = ev(inod,idof)-eu(inod,idof)
        if(iorder(inod) .eq. 1) then
          aa = aa+eue(inod,idof)*eue(inod,idof)
          ab = ab+eue(inod,idof)*edu(inod,idof)
          bb = bb+edu(inod,idof)*edu(inod,idof)
        endif
501   continue
500   continue
      ione = itn
      call sendint(nsource,iblk,ione)
cc@tran_aa_ab_bb
      call sendr(nsource,iblk,aa)
      call sendr(nsource,iblk,ab)
      call sendr(nsource,iblk,bb)
c      cc = 1.0
c      if( ione .gt. 1) then
c        rab = sqrt(aa)*sqrt(bb)
c        if (aa.gt.bb) cc = sqrt(bb/aa)
c        if (ab.gt.0.8*rab) cc = cc*2.0
c        if (ab.lt.0.3*rab) cc = cc*0.5
c        if (ab.lt.0.0) cc = cc*0.5
c        if (ab.lt.-0.40*rab) cc = cc*0.5
c        if (ab.lt.-0.80*rab) cc = cc*0.5
c      endif
c      if (cc.gt.1.0) cc = 1.0
      call recvr(iblk,nsource,cc)
      err = 0.0
      sum = 0.0
      do 502 inod=1,knode
      do 503 idof=1,kdgof
        if( iorder(inod) .eq. 1) then
          err = err+eue(inod,idof)**2
          sum = sum+eu(inod,idof)**2
        endif
        eue(inod,idof) = eue(inod,idof)*cc
        ev(inod,idof) = eu(inod,idof)+eue(inod,idof)
        edu(inod,idof) = 0.0
503   continue
502   continue
c      if (sum .lt. 1.d-20) then
c        err1=dsqrt(err/knode)
c      else
c        err1=dsqrt(err/sum)
c      endif
c      if (err1.lt.emax .or. ione.ge.numstep) then
c      MSend = 1
c      endif
c      call sswitch(nsource,iblk,MSend,1)
      call sendr(nsource,iblk,err)
      call sendr(nsource,iblk,sum)
      call recvint(iblk,nsource,MSend)
      if (MSend.eq.1) then
C.... Open coor0 file
        angle=datan(1.d0)*4.d0/180.d0*angle
        theta=omega*time-angle
        do 504 inod=1,knode
          if (nodvar(inod,4).eq.-1) then
            utmp=(ev(inod,1)+eu1(inod,1))*dt/2.d0+eumn(inod)
            vtmp=(ev(inod,2)+eu1(inod,2))*dt/2.d0+evmn(inod)
            wtmp=(ev(inod,3)+eu1(inod,3))*dt/2.d0+ewmn(inod)
            if (nodvar(inod,1).eq.-1.and.nodvar(inod,2).eq.-1.and.
     &          nodvar(inod,3).eq.-1.and.
     &          idpoint.ne.1.and.dabs(omega).gt.1.d-20) then
              if (id4updown.eq.0) then
                utmp=(dcos(theta)-1.d0)*(coor0(inod,1)-xo)+
     &               dsin(theta)*(coor0(inod,3)-zo)
                wtmp=-dsin(theta)*(coor0(inod,1)-xo)+
     &               (dcos(theta)-1.d0)*(coor0(inod,3)-zo)
              elseif (id4updown.eq.1) then
                utmp=(dcos(theta)-1.d0)*(coor0(inod,1)-xo)-
     &               dsin(theta)*(coor0(inod,2)-zo)
                vtmp=dsin(theta)*(coor0(inod,1)-xo)+
     &               (dcos(theta)-1.d0)*(coor0(inod,2)-zo)
              endif
            endif
            eum(inod)=eum(inod)+(utmp-eum(inod))*relax
            evm(inod)=evm(inod)+(vtmp-evm(inod))*relax
            ewm(inod)=ewm(inod)+(wtmp-ewm(inod))*relax
          endif
504     continue
        call sendar (nsource,iblk,eum,knode)
        call sendar (nsource,iblk,evm,knode)
        call sendar (nsource,iblk,ewm,knode)
        call sendar (nsource,iblk,eumn,knode)
        call sendar (nsource,iblk,evmn,knode)
        call sendar (nsource,iblk,ewmn,knode)
        call sendar (nsource,iblk,evlx,knode)
        call sendar (nsource,iblk,evly,knode)
        call sendar (nsource,iblk,evlz,knode)
        call sendar (nsource,iblk,evlxn,knode)
        call sendar (nsource,iblk,evlyn,knode)
        call sendar (nsource,iblk,evlzn,knode)
C....   Update ubc and bfdc file
CC..    Send ubc to master processor ...
        ntemp = 0
        do i=1,knode2
          if(iorder2(i).eq.1) then
            do j=1,kdgofc
              ntemp=ntemp+1
              apool(ntemp)=ubc(i,j)
            enddo
          endif
        enddo
        call sendar(nsource,iblk,apool,ntemp)
        if(ntemp.ne.nupdate2*kdgofc) then
          print *,'Error! ntemp ne nupdate2 , U program !!',ntemp,nupdate2
          call endjob(ierr)
        endif
C...
        call recvar(iblk,nsource,apool,knode2*kdgofc)
        ntemp=0
        do i=1,knode2
          do j=1,kdgofc
            ntemp=ntemp+1
            ubc(i,j)=apool(ntemp)
          enddo
        enddo
        if(ntemp.ne.knode2*kdgofc) then
          print *,'Error! ntemp ne knode2 , U program !',ntemp, knode2*kdgofc
          call endjob(ierr)
        endif
C
CC Compute the acceleration for the structure and update unods
C
C....   Read unodac
CC..    Send accel to master processor ...
        ntemp = 0
        do i=1,knode1
          if(iorder1(i).eq.1) then
            do j=1,kdgofb
              ntemp=ntemp+1
              apool(ntemp)=accel(i,j)
            enddo
          endif
        enddo
        call sendar(nsource,iblk,apool,ntemp)
        if(ntemp.ne.nupdate1*kdgofb) then
          print *,'Error! ntemp ne nupdate1 , U program !!',ntemp,nupdate1
          call endjob(ierr)
        endif
C..
C...    Send ev to master processor
C..
        call sendar (nsource,iblk,ev,knode*kdgof)
C...
        call recvar(iblk,nsource,ums,knode1)
        call recvar(iblk,nsource,vms,knode1)
        call recvar(iblk,nsource,wms,knode1)
        call recvar(iblk,nsource,apool,knode1*kdgofb)
        ntemp = 0
        do i=1,knode1
          do j=1,kdgofb
            ntemp = ntemp+1
            vs(i,j) = apool(ntemp)
          enddo
        enddo
        call recvar(iblk,nsource,apool,knode1*kdgofb)
        ntemp = 0
        do i=1,knode1
          do j=1,kdgofb
            ntemp = ntemp+1
            acs(i,j) = apool(ntemp)
          enddo
        enddo
      else
        ione=ione+1
        itn = itn + 1
      endif
      do 505 i=1,knode
      do 505 j=1,kdgof
      eu(i,j)=ev(i,j)
505   continue
      do 506 i=1,knode
      do 506 j=1,kdgof
      edu(i,j)=eue(i,j)
506   continue
      call sendar(nsource,iblk,eu,knode*kdgof)
      call sendar(nsource,iblk,eu1,knode*kdgof)
      call sendar(nsource,iblk,edu,knode*kdgof)
      do 507 inod=1,knode
        if (nodvar(inod,4).eq.-1) then
          do 508 idof=1,kdgof
            ev(inod,idof)=0.d0
508       continue
        endif
507   continue
      do 509 i=1,knode
      do 509 j=1,kdgof
      evf(i,j)=ev(i,j)
509   continue
      call sendar (nsource,iblk,evf,knode*kdgof)
c@send

1600  format (1x,i4,6e12.5)
      return
      end
