      subroutine Mmufsic(filename,kkcom)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      common /jcom/ jcom(100)
      common /stop/ msstop,msend,msend_tol
      common /ia/ ia(150000000)
      common /ib/ ib(75000000)
      common /jc/ jc(75000000)

      nsource=0
c.....open partition.dat file
      open(21,file='partition.dat',form='formatted')
      read(21,*) numblk
      close(21)
c      write(*,*) 'numblk ===',numblk
c.....open dispc0 file
      open(23,file=filename(1),form='unformatted',status='unknown')
      read(23) knode2,kdgof2
      close(23)
c.....open disp0 file
      open(23,file=filename(2),form='unformatted',status='unknown')
      read(23) knode,kdgof
      close(23)
c.....open dispb0 file
      open(23,file=filename(3),form='unformatted',status='unknown')
      read(23) knode1,kdgof1
      close(23)
c.....open nbefile
      open(21,file='nbefile',form='formatted',status='old')
      read(21,*) numtyp_body
      close(21)
C.....open coor0 file
      open(21,file='coor0',form='unformatted',status='old')
      read(21) knode,kcoor
      close(21)

      max_lknode=knode*2/numblk+1000
      max_lknode1=knode1*2/numblk+1000
      max_lknode2=knode2*2/numblk+1000
      maxap=max_lknode*kdgof
      kend=0
      kendit=0


      kna2=maxap*2
      kna3=max_lknode*1
      if (kna3/2*2 .lt. kna3) kna3=kna3+1
      knb1=numblk*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb4=knode*kdgof*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb2=numblk*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=numblk*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
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
      kna1=knode2*kdgof2*2
      kna14=knode2*kdgof2*2
      kna15=knode2*kdgof2*2
      kna5=knode*kcoor*2
      knb5=knode*kdgof*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      kna13=knode*1
      if (kna13/2*2 .lt. kna13) kna13=kna13+1
      kna7=knode*kdgof2*2
      kna8=knode*kdgof2*2
      kna9=knode*kdgof2*2
      kna10=knode*kdgof2*2
      kna12=knode*kcoor*2
      kna11=knode*kdgof*2
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
      if (kna15-1.gt.150000000) then
      write(*,*) 'exceed memory of array ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'memory needed = ',kna15,' in prgram mufsic'
      stop 55555
      endif
      knb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      knb4=knb4+knb3
      knb5=knb5+knb4
      if (knb5-1.gt.75000000) then
      write(*,*) 'exceed memory of array ib'
      write(*,*) 'memory of ib = 75000000'
      write(*,*) 'memory needed = ',knb5,' in prgram mufsic'
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
      write(*,*) 'memory needed = ',knc21,' in prgram mufsic'
      stop 55555
      endif
      call mufsic(kcoor,kdgof,kdgof1,kdgof2,
     *knode,knode1,knode2,numblk,max_lknode,max_lknode1,
     *max_lknode2,maxap,nsource,msend,msend_tol,numtyp_body,ia(kna0),
     *ia(kna1),ia(kna2),ia(kna3),ia(kna4),ia(kna5),
     *ia(kna6),ia(kna7),ia(kna8),ia(kna9),ia(kna10),
     *ia(kna11),ia(kna12),ia(kna13),ia(kna14),ib(knb0),
     *ib(knb1),ib(knb2),ib(knb3),ib(knb4),jc(knc0),
     *jc(knc1),jc(knc2),jc(knc3),jc(knc4),jc(knc5),
     *jc(knc6),jc(knc7),jc(knc8),jc(knc9),jc(knc10),
     *jc(knc11),jc(knc12),jc(knc13),jc(knc14),jc(knc15),
     *jc(knc16),jc(knc17),jc(knc18),jc(knc19),jc(knc20),
     *filename)
      return
      end
      subroutine mufsic(kcoor,kdgof,kdgof1,kdgof2,
     *knode,knode1,knode2,numblk,max_lknode,max_lknode1,
     *max_lknode2,maxap,nsource,msend,msend_tol,numtyp_body,
     *u,apool,ipool,urg,coor0,ud,
     *ulap,ulapn,vlap,vlapn,ua,coorlap,
     *ifluidgl,eu,eu1,nextern,nextern1,nextern2,
     *nodvar,id0,knode_iblk,mlgnod,iorder,nupdate_iblk,
     *mnodebo,nnodebo,nod_ord,knode_iblk1,mlgnod1,iorder1,
     *nupdate_iblk1,mnodebo1,nnodebo1,nod_ord1,knode_iblk2,mlgnod2,
     *iorder2,nupdate_iblk2,mnodebo2,nnodebo2,nod_ord2,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
      dimension apool(maxap),ipool(max_lknode)
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
      dimension u(knode2,kdgof2)
      dimension Eu(knode2,kdgof2),Eu1(knode2,kdgof2)
      dimension coor0(knode,kcoor),id0(knode,kdgof),ifluidgl(knode),
     +          ur(knode,kdgof2),ulap(knode,kdgof2),ulapn(knode,kdgof2),
     +          vlap(knode,kdgof2),vlapn(knode,kdgof2),
     +          coorlap(knode,kcoor),ua(knode,kdgof)

      open(88,file='control.dat',status='old')
      read(88,*) emax,numstep,emax,numstep
      close(88)

C.......OPEN TIME File
      OPEN(21,FILE='time',FORM='UNFORMATTED',STATUS='OLD')
      READ(21) TMAX,DT,TIME,IT
c      WRITE(*,*) ' TMAX,DT,TIME,IT =',TMAX,DT,TIME,IT
      CLOSE(21)
C .................................................................
C KDGOF  ..... NUMBER OF D.O.F
C KNODE  ..... NUMBER OF NODES
C .................................................................
      do iblk=1,numblk
        nextern2(iblk)=0
        nupdate2_i=0
        do inod=1,knode_iblk2(iblk)
          if(iorder2(iblk,inod).eq.0) nextern2(iblk)=nextern2(iblk)+1
          if(iorder2(iblk,inod).eq.1) nupdate2_i=nupdate2_i+1
        enddo
        if(nupdate2_i.ne.nupdate_iblk2(iblk)) then
          print *,'nupdate nodes Error !!! nupdate2_i,nupdate_iblk2 ==',
     &            nupdte2_i,nupdate_iblk2(iblk)
          call endjob(ierr)
        endif
      enddo

c
cc     initialize u array
cc
c      do i=1,knode2
c        do j=1,kdgof2
c             u(i,j)=0.0d0
c        enddo
c      enddo
c
c     receive displace at every points from iblk subdomain and set vector u
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk2(iblk)*kdgof2)
        ntemp=0
        do i=1,knode_iblk2(iblk)
          if(iorder2(iblk,i).eq.1) then
             ng_nod=mlgnod2(iblk,i)
             do j=1,kdgof2
               ntemp=ntemp+1
               neqmax=neqmax+1
               u(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode2*kdgof2) then
         print *,"error! different number of neqmax at mu program..!!"
      endif
CC
CC==============================================================CC
c
c     send displace at every points to iblk subdomain
c
      do iblk=1,numblk
        ntemp=0
        do i=1,knode_iblk2(iblk)
          if(iorder2(iblk,i).eq.0) then
            ng_nod=mlgnod2(iblk,i)
            do j=1,kdgof2
              ntemp=ntemp+1
              apool(ntemp)=u(ng_nod,j)
            enddo
          endif
        enddo
        if(ntemp.ne.nextern2(iblk)*kdgof2) then
          print *,'Error! ntemp ne nextern ...',ntemp,nextern2(iblk)
          call endjob(ierr)
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo
C=====================================================C
c
c     receive Eu1
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk2(iblk)*kdgof2)
        ntemp=0
        do i=1,knode_iblk2(iblk)
          if(iorder2(iblk,i).eq.1) then
             ng_nod=mlgnod2(iblk,i)
             do j=1,kdgof2
               ntemp=ntemp+1
               neqmax=neqmax+1
               Eu1(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode2*kdgof2) then
         print *,"error! different number of neqmax at mu program..!!"
         write(*,*) 'Recv Eu1 in Mulaplace:',neqmax, knode2*kdgof2
      endif
C
C...  Read id0,subdomain.idx,coor00 file
C
      open (31,file='id0',form='unformatted',status='old')
      read (31) mmm,nnn,((nodvar(j,i),i=1,kdgof),j=1,knode)
      close(31)
      open (41,file='subdomain.idx',form='unformatted',status='unknown')
      read (41) nnn,knodef,knodes
      read (41) (ifluidgl(i),i=1,knode)
      close(41)
      open (31,file='coor0.0',form='unformatted',status='unknown')
      read (31) mmm,nnn,((coor0(i,j),j=1,kcoor),i=1,knode)
      close(31)
      open(41,file='unodd',form='unformatted',status='old')
      read(41) ((ulap(j,i),j=1,knode),i=1,kdgof2),
     &         ((ulapn(j,i),j=1,knode),i=1,kdgof2),
     &         ((vlap(j,i),j=1,knode),i=1,kdgof2),
     &         ((vlapn(j,i),j=1,knode),i=1,kdgof2)
      close(41)
c
c     receive ur
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk(iblk)*kdgof2)
        ntemp=0
        do i=1,knode_iblk(iblk)
          if(iorder(iblk,i).eq.1) then
             ng_nod=mlgnod(iblk,i)
             do j=1,kdgof2
               ntemp=ntemp+1
               neqmax=neqmax+1
               ur(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode*kdgof2) then
         print *,"error! different number of neqmax at mu program..!!"
         write(*,*) 'Recv ur in Mulaplace:',neqmax, knode*kdgof2
      endif
c
c     receive ulap
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk(iblk)*kdgof2)
        ntemp=0
        do i=1,knode_iblk(iblk)
          if(iorder(iblk,i).eq.1) then
             ng_nod=mlgnod(iblk,i)
             do j=1,kdgof2
               ntemp=ntemp+1
               neqmax=neqmax+1
               ulap(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode*kdgof2) then
         print *,"error! different number of neqmax at mu program..!!"
         write(*,*) 'Recv ulap in Mulaplace:',neqmax, knode*kdgof2
      endif
c
c     receive ulapn
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk(iblk)*kdgof2)
        ntemp=0
        do i=1,knode_iblk(iblk)
          if(iorder(iblk,i).eq.1) then
             ng_nod=mlgnod(iblk,i)
             do j=1,kdgof2
               ntemp=ntemp+1
               neqmax=neqmax+1
               ulapn(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode*kdgof2) then
         print *,"error! different number of neqmax at mu program..!!"
         write(*,*) 'Recv ulapn in Mulaplace:',neqmax, knode*kdgof2
      endif
C..
C....  Compute coorlap, update coor and ulap
C...
      do i=1,knode
        if (i .gt. knode1+knode2) then
          print *,'111111, i, knode1, knode2 =', i,knode1,knode2
          do j=1,kcoor
            coorlap(i,j)=coor0(i,j)
          enddo
        elseif (nodvar(i,4) .eq. -1) then
          do j=1,kcoor
            coorlap(i,j)=coor0(i,j)+ulap(i,j)
          enddo
        else
          jnod=ifluidgl(i)
          do j=1,kcoor
            ulap(i,j)=u(jnod,j)+ur(i,j)
            coor0(i,j)=coor0(i,j)+ulap(i,j)
            coorlap(i,j)=coor0(i,j)
          enddo
        endif
        do j=1,kdgof2
          vlap(i,j)=(ulap(i,j)-ulapn(i,j))/dt
        enddo
      enddo

c
c     Send ulap to slave processor
c
      neqmax=0
      do iblk=1,numblk
        ntemp=0
        do i=1,knode_iblk(iblk)
          ng_nod=mlgnod(iblk,i)
          do j=1,kdgof2
            ntemp=ntemp+1
            neqmax=neqmax+1
            apool(ntemp)=ulap(ng_nod,j)
         enddo
        enddo
        if(ntemp .ne. knode_iblk(iblk)*kdgof2) then
          print *,"error! different number of ntemp at mu program..!!"
          write(*,*) 'Send ulap in Mulaplace:',iblk, ntemp, knode_iblk(iblk)*kdgof2
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo

c
c     Send vlap to slave processor
c
      neqmax=0
      do iblk=1,numblk
        ntemp=0
        do i=1,knode_iblk(iblk)
          ng_nod=mlgnod(iblk,i)
          do j=1,kdgof2
            ntemp=ntemp+1
            neqmax=neqmax+1
            apool(ntemp)=vlap(ng_nod,j)
         enddo
        enddo
        if(ntemp .ne. knode_iblk(iblk)*kdgof2) then
          print *,"error! different number of ntemp at mu program..!!"
          write(*,*) 'Send vlap in Mulaplace:',iblk, ntemp, knode_iblk(iblk)*kdgof2
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo

c
c     Send coor0 to slave processor to update coor0 file
c
      neqmax=0
      do iblk=1,numblk
        ntemp=0
        do i=1,knode_iblk(iblk)
          ng_nod=mlgnod(iblk,i)
          do j=1,kcoor
            ntemp=ntemp+1
            neqmax=neqmax+1
            apool(ntemp)=coor0(ng_nod,j)
         enddo
        enddo
        if(ntemp .ne. knode_iblk(iblk)*kcoor) then
          print *,"error! different number of ntemp at mu program..!!"
          write(*,*) 'Send coor0 in Mulaplace:',iblk, ntemp, knode_iblk(iblk)*kcoor
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo
      open(31,file='coor0',form='unformatted',status='unknown')
      write(31) knode,kcoor,((coor0(i,j),j=1,kcoor),i=1,knode)
      close(31)

c
c     Send coorlap to slave processor to update cooreule file
c
      neqmax=0
      do iblk=1,numblk
        ntemp=0
        do i=1,knode_iblk(iblk)
          ng_nod=mlgnod(iblk,i)
          do j=1,kcoor
            ntemp=ntemp+1
            neqmax=neqmax+1
            apool(ntemp)=coorlap(ng_nod,j)
         enddo
        enddo
        if(ntemp .ne. knode_iblk(iblk)*kcoor) then
          print *,"error! different number of ntemp at mu program..!!"
          write(*,*) 'Send coorlap in Mulaplace:',iblk, ntemp, knode_iblk(iblk)*kcoor
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo
      open(31,file='cooreule',form='unformatted',status='unknown')
      write(31) knode,kcoor,((coorlap(i,j),j=1,kcoor),i=1,knode)
      close(31)
CC
c@tran_aa_ab_bb
C
      err_all =0.D0
      sum_all= 0.D0
      do iblk=1,numblk
        call recvint(nsource,iblk,ione)
        call recvr(nsource,iblk,err)
        call recvr(nsource,iblk,sum)
        err_all = err_all + err
        sum_all = sum_all + sum
      enddo
      err = err_all
      sum = sum_all
      if (sum .lt. 1.d-20) then
        err1=dsqrt(err/knode2)
      else
        err1=dsqrt(err/sum)
      endif
      WRITE(*,'(a,2E12.5,i3,i6,E12.5)') 'Mesh ERR,sum,step,it,time =',ERR1,sum,ione,it,time
      if (err1 .lt. emax .or. ione .ge. numstep) then
        MSend_tol = 1
      endif
      do iblk=1,numblk
        call sendint(iblk,nsource,MSend_tol)
      enddo
c      call mswitch(nsource,MSend)
CC
CC...
CC
      if( MSend_tol .eq. 1) then
        open(41,file='unodd',form='unformatted',status='old')
        write(41) ((ulap(j,i),j=1,knode),i=1,kdgof2),
     &            ((ulap(j,i),j=1,knode),i=1,kdgof2),
     &            ((vlap(j,i),j=1,knode),i=1,kdgof2),
     &            ((vlap(j,i),j=1,knode),i=1,kdgof2)
        close(41)
        open(51,file='unod',form='unformatted',status='old')
        read(51) ((ua(j,i),j=1,knode),i=1,kdgof)
        rewind(51)
        write(51) ((ua(j,i),j=1,knode),i=1,kdgof),
     &            ((ua(j,i),j=1,knode),i=1,kdgof),
     &            ((0.d0,j=1,knode),i=1,kdgof)
        close(51)
        open(42,file='vangle',form='formatted',status='old')
        read(42,*) vangle,vangle1,angle,angle1
        rewind(42)
        write(42,*) vangle,vangle,angle,angle
        close(42)
c        open(51,file='unods',form='unformatted',status='old')
c        read(51) ((ulap(j,i),j=1,knode1),i=1,kdgof1),
c     &           ((ulap(j,i),j=1,knode1),i=1,kdgof1),
c     &           ((ulap(j,i),j=1,knode1),i=1,kdgof1)
c        close(51)
        open(51,file='unodac',form='unformatted',status='unknown')
        write(51) ((ulap(j,i),j=1,knode1),i=1,kdgof1)
        close(51)
      else
        open(41,file='unodd',form='unformatted',status='old')
        write(41) ((ulap(j,i),j=1,knode),i=1,kdgof2),
     &            ((ulapn(j,i),j=1,knode),i=1,kdgof2),
     &            ((vlap(j,i),j=1,knode),i=1,kdgof2),
     &            ((vlapn(j,i),j=1,knode),i=1,kdgof2)
        close(41)
      endif
      nunit = 21
      open(nunit,file='unodc',form='unformatted',status='unknown')
      call mrecvu2d(nsource,numblk,knode2,kdgof2,apool,
     & Eu,knode_iblk2,mlgnod2,iorder2)
      write(nunit)
     & ((Eu(i,j),i=1,knode2),j=1,kdgof2)
      close(nunit)
c@mend
      END
