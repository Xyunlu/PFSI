      subroutine Mmufsia(filename,kkcom)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      common /jcom/ jcom(100)
      common /stop/ msstop,msend
      common /ia/ ia(150000000)
      common /ib/ ib(75000000)
      common /jc/ jc(75000000)

      nsource=0
c.....open partition.dat file
      open(21,file='partition.dat',form='formatted')
      read(21,*) numblk
      close(21)
c      write(*,*) 'numblk ===',numblk
c.....open disp0 file
      open(23,file=filename(1),form='unformatted',status='unknown')
      read(23) knode,kdgof
      close(23)
c.....open dispb0 file
      open(23,file=filename(2),form='unformatted',status='unknown')
      read(23) knode1,kdgof1
      close(23)
c.....open dispc0 file
      open(23,file=filename(3),form='unformatted',status='unknown')
      read(23) knode2,kdgof2
      close(23)
c.....open nbefile
      open(21,file='nbefile',form='formatted',status='old')
      read(21,*) numtyp_body
      close(21)

      max_lknode=knode*2/numblk+1000
      max_lknode1=knode1*2/numblk+1000
      max_lknode2=knode2*2/numblk+1000
      maxap=max_lknode*kdgof
      kend=0
      kendit=0


      kna1=knode*kdgof*2
      kna2=maxap*2
      kna3=max_lknode*1
      if (kna3/2*2 .lt. kna3) kna3=kna3+1
      knb1=numblk*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb4=knode*kdgof*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb2=numblk*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knc1=numblk*1
      if (knc1/2*2 .lt. knc1) knc1=knc1+1
      knc2=numblk*max_lknode*1
      if (knc2/2*2 .lt. knc2) knc2=knc2+1
      knc3=numblk*max_lknode*1
      if (knc3/2*2 .lt. knc3) knc3=knc3+1
      knc4=numblk*1
      if (knc4/2*2 .lt. knc4) knc4=knc4+1
      knc5=numtyp_body*1
      if (knc5/2*2 .lt. knc5) knc5=knc5+1
      knc6=numtyp_body*1
      if (knc6/2*2 .lt. knc6) knc6=knc6+1
      knc7=knode*1
      if (knc7/2*2 .lt. knc7) knc7=knc7+1
      knc8=numblk*1
      if (knc8/2*2 .lt. knc8) knc8=knc8+1
      knc9=numblk*max_lknode1*1
      if (knc9/2*2 .lt. knc9) knc9=knc9+1
      knc10=numblk*max_lknode1*1
      if (knc10/2*2 .lt. knc10) knc10=knc10+1
      knc11=numblk*1
      if (knc11/2*2 .lt. knc11) knc11=knc11+1
      knc12=numtyp_body*1
      if (knc12/2*2 .lt. knc12) knc12=knc12+1
      knc13=numtyp_body*1
      if (knc13/2*2 .lt. knc13) knc13=knc13+1
      knc14=knode1*1
      if (knc14/2*2 .lt. knc14) knc14=knc14+1
      knc15=numblk*1
      if (knc15/2*2 .lt. knc15) knc15=knc15+1
      knc16=numblk*max_lknode2*1
      if (knc16/2*2 .lt. knc16) knc16=knc16+1
      knc17=numblk*max_lknode2*1
      if (knc17/2*2 .lt. knc17) knc17=knc17+1
      knc18=numblk*1
      if (knc18/2*2 .lt. knc18) knc18=knc18+1
      knc19=numtyp_body*1
      if (knc19/2*2 .lt. knc19) knc19=knc19+1
      knc20=numtyp_body*1
      if (knc20/2*2 .lt. knc20) knc20=knc20+1
      knc21=knode2*1
      if (knc21/2*2 .lt. knc21) knc21=knc21+1
      kna4=knode*kdgof*2
      kna5=knode*kdgof*2
      kna6=knode*kdgof*2
      kna7=knode*kdgof*2
      kna8=knode*kdgof*2
      kna9=knode*kdgof*2
      kna10=knode*2
      kna11=knode*2
      kna12=knode*2
      kna13=knode*2
      kna14=knode*2
      kna15=knode*2
      kna16=knode*2
      kna17=knode*2
      kna18=knode*2
      kna19=knode*2
      kna20=knode*2
      kna21=knode*2
      jcom(1) = 1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
      kna7=kna7+kna6
      kna8=kna8+kna7
      kna9=kna9+kna8
      kna10=kna10+kna9
      kna11=kna11+kna10
      kna12=kna12+kna11
      kna13=kna13+kna12
      kna14=kna14+kna13
      kna15=kna15+kna14
      kna16=kna16+kna15
      kna17=kna17+kna16
      kna18=kna18+kna17
      kna19=kna19+kna18
      kna20=kna20+kna19
      kna21=kna21+kna20
      if (kna21-1.gt.150000000) then
      write(*,*) 'exceed memory of array ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'memory needed = ',kna21,' in prgram mufsia'
      stop 55555
      endif
      knb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      knb4=knb4+knb3
      if (knb4-1.gt.75000000) then
      write(*,*) 'exceed memory of array ib'
      write(*,*) 'memory of ib = 75000000'
      write(*,*) 'memory needed = ',knb4,' in prgram mufsia'
      stop 55555
      endif
      kcom=kkcom+1
      knc0=jcom(kcom)
      knc1=knc1+knc0
      knc2=knc2+knc1
      knc3=knc3+knc2
      knc4=knc4+knc3
      knc5=knc5+knc4
      knc6=knc6+knc5
      knc7=knc7+knc6
      knc8=knc8+knc7
      knc9=knc9+knc8
      knc10=knc10+knc9
      knc11=knc11+knc10
      knc12=knc12+knc11
      knc13=knc13+knc12
      knc14=knc14+knc13
      knc15=knc15+knc14
      knc16=knc16+knc15
      knc17=knc17+knc16
      knc18=knc18+knc17
      knc19=knc19+knc18
      knc20=knc20+knc19
      knc21=knc21+knc20
      jcom(kcom+1) = knc21
      if (knc21-1.gt.75000000) then
      write(*,*) 'exceed memory of array jc'
      write(*,*) 'memory of jc = 75000000'
      write(*,*) 'memory needed = ',knc21,' in prgram mufsia'
      stop 55555
      endif
      call mufsia(knode,knode1,knode2,kdgof,
     *kdgof1,kdgof2,numblk,max_lknode,max_lknode1,max_lknode2,
     *maxap,nsource,msend,numtyp_body,ia(kna0),ia(kna1),
     *ia(kna2),ia(kna3),ia(kna4),ia(kna5),ia(kna6),
     *ia(kna7),ia(kna8),ia(kna9),ia(kna10),ia(kna11),
     *ia(kna12),ia(kna13),ia(kna14),ia(kna15),ia(kna16),
     *ia(kna17),ia(kna18),ia(kna19),ia(kna20),ib(knb0),
     *ib(knb1),ib(knb2),ib(knb3),jc(knc0),jc(knc1),
     *jc(knc2),jc(knc3),jc(knc4),jc(knc5),jc(knc6),
     *jc(knc7),jc(knc8),jc(knc9),jc(knc10),jc(knc11),
     *jc(knc12),jc(knc13),jc(knc14),jc(knc15),jc(knc16),
     *jc(knc17),jc(knc18),jc(knc19),jc(knc20),
     *filename)
      return
      end
      subroutine mufsia(knode,knode1,knode2,kdgof,
     *kdgof1,kdgof2,numblk,max_lknode,max_lknode1,max_lknode2,
     *maxap,nsource,msend,numtyp_body,u,apool,
     *ipool,eu1,eu,ev,eue,edu,
     *evf,eum,evm,ewm,eumn,evmn,
     *ewmn,evlx,evly,evlz,evlxn,evlyn,
     *evlzn,nextern,nextern1,nertern2,nodvar,knode_iblk,
     *mlgnod,iorder,nupdate_iblk,mnodebo,nnodebo,nod_ord,
     *knode_iblk1,mlgnod1,iorder1,nupdate_iblk1,mnodebo1,nnodebo1,
     *nod_ord1,knode_iblk2,mlgnod2,iorder2,nupdate_iblk2,mnodebo2,
     *nnodebo2,nod_ord2,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
      dimension u(knode,kdgof),apool(maxap),ipool(max_lknode)
      dimension nextern(numblk),nodvar(knode,kdgof)
      dimension nextern1(numblk),nextern2(numblk)
      dimension knode_iblk(numblk),mlgnod(numblk,max_lknode),
     +          iorder(numblk,max_lknode)
      dimension nupdate_iblk(numblk)
      dimension mnodebo(numtyp_body),nnodebo(numtyp_body)
      dimension nod_ord(knode)
      dimension knode_iblk1(numblk),mlgnod1(numblk,max_lknode1),
     +          iorder1(numblk,max_lknode1)
      dimension nupdate_iblk1(numblk)
      dimension mnodebo1(numtyp_body),nnodebo1(numtyp_body)
      dimension nod_ord1(knode1)
      dimension knode_iblk2(numblk),mlgnod2(numblk,max_lknode2),
     +          iorder2(numblk,max_lknode2)
      dimension nupdate_iblk2(numblk)
      dimension mnodebo2(numtyp_body),nnodebo2(numtyp_body)
      dimension nod_ord2(knode2)
      dimension Eu1(knode,kdgof),Eu(knode,kdgof),
     & Ev(knode,kdgof),Eue(knode,kdgof),Edu(knode,kdgof),
     & Evf(knode,kdgof)
      dimension Eum(knode),Evm(knode),
     & Ewm(knode),Eumn(knode),Evmn(knode),
     & Ewmn(knode),Evlx(knode),Evly(knode),
     & Evlz(knode),Evlxn(knode),Evlyn(knode),
     & Evlzn(knode)
      logical filflg1, filflg
      real*8,allocatable::accel(:,:)
      integer,allocatable::ifluidgl(:)
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
      open (41,file='vangle',form='formatted',status='old')
      read (41,*) vangle,vangle1,angle,angle1
      close(41)

C .................................................................
C KDGOF  ..... NUMBER OF D.O.F
C KNODE  ..... NUMBER OF NODES
C .................................................................
      open(21,file='time0',form='formatted',status='old')
      read(21,*) t0,tmax,dt
      close(21)

      do iblk=1,numblk
        nextern(iblk)=0
        nupdate_i=0
        do inod=1,knode_iblk(iblk)
          if(iorder(iblk,inod).eq.0) nextern(iblk)=nextern(iblk)+1
          if(iorder(iblk,inod).eq.1) nupdate_i=nupdate_i+1
        enddo
        if(nupdate_i.ne.nupdate_iblk(iblk)) then
          print *,'nupdate nodes Error !!!!! nupdate_i,nupdate_iblk ==',
     &            nupdte_i,nupdate_iblk(iblk)
          call endjob(ierr)
        endif
      enddo

c
cc     initialize u array
cc
c      do i=1,knode
c        do j=1,kdgof
c             u(i,j)=0.0d0
c        enddo
c      enddo
c
c     receive displace at every points from iblk subdomain and set vector u
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk(iblk)*kdgof)
        ntemp=0
        do i=1,knode_iblk(iblk)
          if(iorder(iblk,i).eq.1) then
             ng_nod=mlgnod(iblk,i)
             do j=1,kdgof
               ntemp=ntemp+1
               neqmax=neqmax+1
               u(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode*kdgof) then
         print *,"error! different number of neqmax at mu program..!!"
      endif
CC
CC==============================================================CC
c
c     send displace at every points to iblk subdomain
c
      do iblk=1,numblk
        ntemp=0
        do i=1,knode_iblk(iblk)
          if(iorder(iblk,i).eq.0) then
            ng_nod=mlgnod(iblk,i)
            do j=1,kdgof
              ntemp=ntemp+1
              apool(ntemp)=u(ng_nod,j)
            enddo
          endif
        enddo
        if(ntemp.ne.nextern(iblk)*kdgof) then
          print *,'Error! ntemp ne nextern ...',ntemp,nextern(iblk)
          call endjob(ierr)
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo
      do iblk=1,numblk
        call recvint(nsource,iblk,ione)
      enddo
CC
cc@tran_aa_ab_bb
      aa_all = 0.D0
      ab_all = 0.D0
      bb_all = 0.D0
      do iblk=1,numblk
        call recvr(nsource,iblk,aa)
        call recvr(nsource,iblk,ab)
        call recvr(nsource,iblk,bb)
        aa_all = aa_all + aa
        ab_all = ab_all + ab
        bb_all = bb_all + bb
      enddo
CC
      aa = aa_all
      ab = ab_all
      bb = bb_all
      cc = 1.0D0
      if( ione .gt. 1) then
        rab = sqrt(aa)*sqrt(bb)
        if (aa.gt.bb) cc = sqrt(bb/aa)
        if (ab.gt.0.8*rab) cc = cc*2.0
        if (ab.lt.0.3*rab) cc = cc*0.5
        if (ab.lt.0.0) cc = cc*0.5
        if (ab.lt.-0.40*rab) cc = cc*0.5
        if (ab.lt.-0.80*rab) cc = cc*0.5
      endif
      if (cc.gt.1.0) cc = 1.0
      do iblk=1,numblk
        call sendr(iblk,nsource,cc)
      enddo
      write(*,'(a,3E12.5)') 'aa, ab, bb=', aa, ab, bb
      write(*,'(a,E12.5)') 'cc =', cc
C
      err_all =0.D0
      sum_all= 0.D0
      do iblk=1,numblk
        call recvr(nsource,iblk,err)
        call recvr(nsource,iblk,sum)
        err_all = err_all + err
        sum_all = sum_all + sum
      enddo
      err = err_all
      sum = sum_all
      if (sum .lt. 1.d-20) then
        err1=dsqrt(err/knode)
      else
        err1=dsqrt(err/sum)
      endif
      write(*,'(a,3E12.5,i3)') 'inner err1,err,sum,step =',err1,err,sum,ione
      if (err1.lt.emax .or. ione.ge.numstep) then
        MSend = 1
      endif
      do iblk=1,numblk
        call sendint(iblk,nsource,MSend)
      enddo
c      call mswitch(nsource,MSend)
CC
CC...
CC
      nunit = 21
      if( MSend .eq. 1) then
        open(21,file='id0',form='unformatted',status='old')
        read(21) mm,nn,((nodvar(i,j),j=1,kdgof),i=1,knode)
        close(21)
        open(nunit,file='unodd',form='unformatted',status='unknown')
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Eum,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evm,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Ewm,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Eumn,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evmn,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Ewmn,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evlx,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evly,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evlz,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evlxn,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evlyn,knode_iblk,mlgnod,iorder)
        call mrecvu1d(nsource,numblk,knode,apool,
     &                Evlzn,knode_iblk,mlgnod,iorder)
        write(nunit)
     &       (Eum(i),i=1,knode),
     &       (Evm(i),i=1,knode),
     &       (Ewm(i),i=1,knode),
     &       (Eumn(i),i=1,knode),
     &       (Evmn(i),i=1,knode),
     &       (Ewmn(i),i=1,knode),
     &       (Evlx(i),i=1,knode),
     &       (Evly(i),i=1,knode),
     &       (Evlz(i),i=1,knode),
     &       (Evlxn(i),i=1,knode),
     &       (Evlyn(i),i=1,knode),
     &       (Evlzn(i),i=1,knode)
        close(nunit)

C
C.... Read subdomain file
C
        allocate(ifluidgl(knode))
        open(41,file='subdomain.idx',form='unformatted',status='unknown')
        read(41) nn,knodef,knodes,knodefs
        read(41) (ifluidgl(i),i=1,knode)
        close(41)
C....
C.... Update bfdc file
C....
        do iblk=1,numblk
          nextern2(iblk)=0
          nupdate_i=0
          do inod=1,knode_iblk2(iblk)
            if(iorder2(iblk,inod).eq.0) nextern2(iblk)=nextern2(iblk)+1
            if(iorder2(iblk,inod).eq.1) nupdate_i=nupdate_i+1
          enddo
          if(nupdate_i.ne.nupdate_iblk2(iblk)) then
            print *,'nupdate nodes Error !!! nupdate_i,nupdate_iblk ==',
     &              nupdte_i,nupdate_iblk2(iblk)
            call endjob(ierr)
          endif
        enddo
C...
C.... Recv ubc from slave processor
C...
        neqmax2=0
        do iblk=1,numblk
          call recvar(nsource,iblk,apool,nupdate_iblk2(iblk)*kdgof2)
          ntemp2=0
          do i=1,knode_iblk2(iblk)
            if(iorder2(iblk,i).eq.1) then
               ng_nod=mlgnod2(iblk,i)
               do j=1,kdgof2
                 ntemp2=ntemp2+1
                 neqmax2=neqmax2+1
                 u(ng_nod,j)=apool(ntemp2)
              enddo
            endif
          enddo
        enddo
        if(neqmax2.ne.knode2*kdgof2) then
           print *,"error! different number of neqmax2 at mu program.!!",neqmax2, knode2*kdgof2
        endif
C
CC    Compute bfdc
CC
        do inod=1,knode
          if(nodvar(inod,4).eq.-1 .and. nodvar(inod,1).lt.-1) then
            u(ifluidgl(-nodvar(inod,1)-1),1)=eum(inod)
            u(ifluidgl(-nodvar(inod,2)-1),2)=evm(inod)
            u(ifluidgl(-nodvar(inod,3)-1),3)=ewm(inod)
          endif
        enddo
        open(21,file='bfdc',form='unformatted',status='unknown')
        write(21) ((u(i,j),j=1,kdgof2),i=1,knode2)
        close(21)
C...
C...  Send ubc to slave processor
        do iblk=1,numblk
          ntemp=0
          do i=1,knode_iblk2(iblk)
            ng_nod=mlgnod2(iblk,i)
            do j=1,kdgof2
              ntemp=ntemp+1
              apool(ntemp)=u(ng_nod,j)
            enddo
          enddo
          if(ntemp.ne.knode_iblk2(iblk)*kdgof2) then
            print *,'Error! ntemp ne nextern ...',ntemp,nextern(iblk)
            call endjob(ierr)
          endif
          call sendar(iblk,nsource,apool,ntemp)
        enddo
C
C
CC Compute the acceleration for the structure
C
C...
C.... Recv accel from slave processor
C...
        allocate(accel(knode1,kdgof1))
        neqmax1=0
        do iblk=1,numblk
          call recvar(nsource,iblk,apool,nupdate_iblk1(iblk)*kdgof1)
          ntemp1=0
          do i=1,knode_iblk1(iblk)
            if(iorder1(iblk,i).eq.1) then
              ng_nod=mlgnod1(iblk,i)
              do j=1,kdgof1
                ntemp1=ntemp1+1
                neqmax1=neqmax1+1
                accel(ng_nod,j)=apool(ntemp2)
              enddo
            endif
          enddo
        enddo
        if(neqmax1.ne.knode1*kdgof1) then
           print *,"error! different number of neqmax1 at mu program.!!",neqmax1, knode1*kdgof1
        endif
C
CCC
C
        do k=knode2+1,knode2+knode1
          m=ifluidgl(k)
          accel(m,1)=(eum(k)-eumn(k))/(dt*dt/4.d0)-dt/4.d0*evlxn(k)
     &                -accel(m,1)
          accel(m,2)=(evm(k)-evmn(k))/(dt*dt/4.d0)-dt/4.d0*evlyn(k)
     &                -accel(m,2)
          accel(m,3)=(ewm(k)-ewmn(k))/(dt*dt/4.d0)-dt/4.d0*evlzn(k)
     &                -accel(m,3)
        enddo
C....
C....  recv ev from slave processor
C...
        call mrecvu2d(nsource,numblk,knode,kdgof,apool,
     &                Eu,knode_iblk,mlgnod,iorder)
C
C.... Update unods file
C
C.... Send ums to slave
        do iblk=1,numblk
          do i=1,knode_iblk1(iblk)
            ng_nod = mlgnod1(iblk,i)
            apool(i) = eum(ng_nod+knode2)
          enddo
          call sendar(iblk,nsource,apool,knode_iblk1(iblk))
        enddo
C.... Send vms to slave
        do iblk=1,numblk
          do i=1,knode_iblk1(iblk)
            ng_nod = mlgnod1(iblk,i)
            apool(i) = evm(ng_nod+knode2)
          enddo
          call sendar(iblk,nsource,apool,knode_iblk1(iblk))
        enddo
C.... Send wms to slave
        do iblk=1,numblk
          do i=1,knode_iblk1(iblk)
            ng_nod = mlgnod1(iblk,i)
            apool(i) = ewm(ng_nod+knode2)
          enddo
          call sendar(iblk,nsource,apool,knode_iblk1(iblk))
        enddo
C.... Send vs to slave
        do iblk=1,numblk
          ntemp = 0
          do i=1,knode_iblk1(iblk)
            ng_nod = mlgnod1(iblk,i)
            do j=1,kdgof1
              ntemp = ntemp+1
              apool(ntemp) = eu(ng_nod+knode2,j)
            enddo
          enddo
          call sendar(iblk,nsource,apool,knode_iblk1(iblk)*kdgof1)
        enddo
C.... Send acs to slave
        do iblk=1,numblk
          ntemp = 0
          do i=1,knode_iblk1(iblk)
            ng_nod = mlgnod1(iblk,i)
            do j=1,kdgof1
              ntemp = ntemp+1
              apool(ntemp) = accel(i,j)
            enddo
          enddo
          call sendar(iblk,nsource,apool,knode_iblk1(iblk)*kdgof1)
        enddo
C
        open(33,file='unods',form='unformatted',status='unknown')
        write(33) (eum(i),i=knode2+1,knode2+knode1),
     &    (evm(i),i=knode2+1,knode2+knode1),
     &    (ewm(i),i=knode2+1,knode2+knode1),
     &    ((ev(i,k),i=knode2+1,knode2+knode1),k=1,kdgof1),
     &    ((accel(i,j),i=1,knode1),j=1,kdgof1)
        close(33)
C
        deallocate(ifluidgl, accel)

      endif

      open(nunit,file='unod',form='unformatted',status='unknown')
      call mrecvu2d(nsource,numblk,knode,kdgof,apool,
     & Eu,knode_iblk,mlgnod,iorder)
      call mrecvu2d(nsource,numblk,knode,kdgof,apool,
     & Eu1,knode_iblk,mlgnod,iorder)
      call mrecvu2d(nsource,numblk,knode,kdgof,apool,
     & Edu,knode_iblk,mlgnod,iorder)
      write(nunit)
     & ((Eu(i,j),i=1,knode),j=1,kdgof),
     & ((Eu1(i,j),i=1,knode),j=1,kdgof),
     & ((Edu(i,j),i=1,knode),j=1,kdgof)
      close(nunit)
      open(nunit,file='unodf',form='unformatted',status='unknown')
      call mrecvu2d(nsource,numblk,knode,kdgof,apool,
     & Evf,knode_iblk,mlgnod,iorder)
      write(nunit)
     & ((Evf(i,j),i=1,knode),j=1,kdgof)
      close(nunit)
c@mend
      END
