      subroutine Mmufsib(filename,kkcom)
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
c.....open dispb0 file
      open(23,file=filename(1),form='unformatted',status='unknown')
      read(23) knode1,kdgof1
      close(23)
c.....open disp0 file
      open(23,file=filename(2),form='unformatted',status='unknown')
      read(23) knode,kdgof
      close(23)
c.....open dispc0 file
      open(23,file=filename(3),form='unformatted',status='unknown')
      read(23) knode2,kdgof2
      close(23)
c.....open nbefile
      open(21,file='nbefile',form='formatted',status='old')
      read(21,*) numtyp_body
      close(21)
      kcoor = 3
 
      max_lknode=knode*2/numblk+1000
      max_lknode1=knode1*2/numblk+1000
      max_lknode2=knode2*2/numblk+1000
      maxap=max_lknode*kdgof
      kend=0
      kendit=0
 
 
      kna1=maxap*2
      kna2=max_lknode*1
      if (kna2/2*2 .lt. kna2) kna2=kna2+1
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
      kna3=knode*kdgof1*2
      kna4=knode*kcoor*2
      kna5=knode2*kdgof2*2
      jcom(1) = 1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      if (kna5-1.gt.150000000) then
      write(*,*) 'exceed memory of array ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'memory needed = ',kna5,' in prgram mufsib'
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
      write(*,*) 'memory needed = ',knb4,' in prgram mufsib'
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
      write(*,*) 'memory needed = ',knc21,' in prgram mufsib'
      stop 55555
      endif
      call mufsib(knode,knode1,knode2,kdgof,
     *kdgof1,kdgof2,kcoor,numblk,max_lknode,max_lknode1,
     *max_lknode2,maxap,nsource,msend,numtyp_body,ia(kna0),ia(kna1),
     *ia(kna2),ia(kna3),ia(kna4),ib(knb0),ib(knb1),
     *ib(knb2),ib(knb3),jc(knc0),jc(knc1),jc(knc2),
     *jc(knc3),jc(knc4),jc(knc5),jc(knc6),jc(knc7),
     *jc(knc8),jc(knc9),jc(knc10),jc(knc11),jc(knc12),
     *jc(knc13),jc(knc14),jc(knc15),jc(knc16),jc(knc17),
     *jc(knc18),jc(knc19),jc(knc20),
     *filename)
      return
      end
      subroutine mufsib(knode,knode1,knode2,kdgof,
     *kdgof1,kdgof2,kcoor,numblk,max_lknode,max_lknode1,
     *max_lknode2,maxap,nsource,msend,numtyp_body,apool,
     *ipool,urg,coor,ubc,nextern,nextern1,
     *nertern2,nodvar,knode_iblk,mlgnod,iorder,nupdate_iblk,
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
      dimension urg(knode,kdgof1),coor(knode,kcoor),ubc(knode2,kdgof2)
      logical filflg1, filflg
      integer,allocatable::ifluidgl(:), id_elm(:)
      double precision,dimension(:),allocatable :: force,rotiner,rpool
      common /rdata/ vangle,angle,penalty,pen4sd
      common /idata/ id4updown
 
C..... Recv and compute the force and rotiner
 
c......
      open (40,file='init.dat',form='formatted',status='old')
      read (40,*) xmax,ymax,zmax,xo,yo,zo,radius,thick,radaxis,angle0
      read (40,*) nx,nz,ns,npt,nct,nsc,naxis,nyi,nyo
      read (40,*) skip
      read (40,*) skip
      read (40,*) skip
      read (40,*) skip
      read (40,*) skip
      read (40,*) id4updown
      read (40,*) index4pro
      close(40)
      angle0=datan(1.d0)*4.d0/180.d0*angle0
 
      open (88,file='control.dat',status='old')
      read (88,*) emax,numstep,emax,numstep,relax
      close(88)
 
C .................................................................
C KDGOF  ..... NUMBER OF D.O.F
C KNODE  ..... NUMBER OF NODES
C .................................................................
      open(21,file='time0',form='formatted',status='old')
      read(21,*) t0,tmax,dt
      close(21)
 
      open(21,file=filename(4),form='unformatted',status='old')
      open(22,file='mlgelm_no_1',form='formatted',status='old')
      read(22,*) nnn, numtyp
      do ityp=1,numtyp
        read(21) numg,nnode
        read(21) nnn,mmm
        allocate(id_elm(numg),force(numg),rotiner(numg))
        do ne=1, numg
          id_elm(ne) = 0
          force(ne) = 0.D0
          rotiner(ne) = 0.D0
        enddo
        do iblk=1,numblk
          read(22,*) numi
          read(22,*) (ipool(i),i=1,numi)
          allocate(rpool(numi))
          call recvar(nsource,iblk,apool,numi)
          call recvar(nsource,iblk,rpool,numi)
          do ie=1,numi
            ieg = ipool(ie)
            if(ieg .gt. numg) then
              write(*,*) 'Error! ieg, numg = ', ieg, numg
              call endjob(ierr)
            endif
            if( id_elm(ieg) .eq. 0) then
              id_elm(ieg) = 1
              rotiner(ieg) = apool(ie)
              force(ieg) = rpool(ie)
            else
              if(dabs(apool(ie)-rotiner(ieg)) .gt. 1.D-06) then
                write(*,*) 'Error!! apool, rotiner:', apool(ie),rotiner(ieg)
              endif
              if(dabs(rpool(ie)-force(ieg)) .gt. 1.D-06) then
                write(*,*) 'Error!! rpool, force:',rpool(ie),force(ieg)
              endif
            endif
          enddo
          deallocate(rpool)
        enddo
 
        forceg = 0.D0
        rotinerg = 0.D0
        do ie=1,numg
          forceg = forceg + force(ie)
          rotinerg = rotinerg + rotiner(ie)
        enddo
        write(*,*) 'forceg, rotinerg=',forceg, rotinerg
        deallocate(id_elm,force,rotiner)
      enddo
      close(22)
      close(21)
 
      open (40,file='vangle',form='formatted',status='old')
      read (40,*) vangle,vangle1,angle,angle1
      rewind(40)
 
      vangle=forceg/rotinerg
      angle=angle1+(vangle+vangle1)*dt/2.d0
      write(40,*) vangle,vangle1,angle,angle1
      close(40)

      write(*,*) 'vangle,vangle1 =', vangle, vangle1
      write(*,*) 'angle, angle1 = ',angle, angle1
 
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
c     receive coor0 from slave processor
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk(iblk)*kdgof1)
        ntemp=0
        do i=1,knode_iblk(iblk)
          if(iorder(iblk,i).eq.1) then
             ng_nod=mlgnod(iblk,i)
             do j=1,kcoor
               ntemp=ntemp+1
               neqmax=neqmax+1
               coor(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode*kcoor) then
         print *,"error! different number of neqmax at mu program..!!"
      endif
 
c
c     receive ur from slave processor
c
      neqmax=0
      do iblk=1,numblk
        call recvar(nsource,iblk,apool,nupdate_iblk(iblk)*kdgof1)
        ntemp=0
        do i=1,knode_iblk(iblk)
          if(iorder(iblk,i).eq.1) then
             ng_nod=mlgnod(iblk,i)
             do j=1,kdgof1
               ntemp=ntemp+1
               neqmax=neqmax+1
               urg(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode*kdgof1) then
         print *,"error! different number of neqmax at mu program..!!"
      endif
CC
CC==============================================================CC
      allocate(ifluidgl(knode))
      open (41,file='subdomain.idx',form='unformatted',status='unknown')
      read (41) nnn,knodef,knodes,knodefs
      read (41) (ifluidgl(i),i=1,knode)
      close(41)
c      write(*,*) 'In mufsib::: '
c      write(*,*) 'knodef,knodes,knodefs =', knodef,knodes,knodefs
c      write(*,*) 'knode,knode1,knode2=', knode, knode1,knode2
C..
C....  Compute urg
C..
      do inod=1,knode
         if (inod.gt.knodef+knodes) then
            urg(inod,1)=0.d0
            urg(inod,2)=0.d0
            urg(inod,3)=0.d0
         else
            if (id4updown.eq.0) then
               utmp=(dcos(angle-angle0)-1.d0)*(coor(inod,1)-xo)+
     &              dsin(angle-angle0)*(coor(inod,3)-zo)
               urg(inod,1)=urg(inod,1)+(utmp-urg(inod,1))*relax
               urg(inod,2)=0.d0
               utmp=-dsin(angle-angle0)*(coor(inod,1)-xo)+
     &              (dcos(angle-angle0)-1.d0)*(coor(inod,3)-zo)
               urg(inod,3)=urg(inod,3)+(utmp-urg(inod,3))*relax
            elseif (id4updown.eq.1) then
               utmp=(dcos(angle-angle0)-1.d0)*(coor(inod,1)-xo)-
     &              dsin(angle-angle0)*(coor(inod,2)-zo)
               urg(inod,1)=urg(inod,1)+(utmp-urg(inod,1))*relax
               urg(inod,3)=0.d0
               utmp=dsin(angle-angle0)*(coor(inod,1)-xo)+
     &              (dcos(angle-angle0)-1.d0)*(coor(inod,2)-zo)
               urg(inod,2)=urg(inod,2)+(utmp-urg(inod,2))*relax
            endif
         endif
      enddo
C...
C...  Update coor0
C...
      do i=1,knode
        do j=1,kcoor
          coor(i,j)=coor(i,j)+urg(i,j)
        enddo
      enddo
c
c     receive ubc from slave processor
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
               ubc(ng_nod,j)=apool(ntemp)
            enddo
          endif
        enddo
      enddo
      if(neqmax.ne.knode2*kdgof2) then
         print *,"error! different number of neqmax at mu program..!!"
      endif
 
      open (51,file='id0',form='unformatted',status='old')
      read (51) mmm,nnn,((nodvar(j,i),i=1,kdgof),j=1,knode)
      close(51)
c
ccc ale's a @interface = (u_s - u_r)
c
      do j=1,knode
         if (nodvar(j,4).eq.-1 .and. nodvar(j,1).lt.-1) then
               k=ifluidgl(-nodvar(j,1)-1)
               do i=1,kdgof2
                  ubc(k,i)=ubc(k,i)-urg(j,i)
               enddo
         endif
      enddo
 
C...
C...  Save bfdc and send to slave processor
C...
      open(21,file='bfdc',form='unformatted',status='unknown')
      write(21) ((ubc(i,j),j=1,kdgof2),i=1,knode2)
      close(21)
C...
C...  Send ubc to slave processor
      do iblk=1,numblk
        ntemp=0
        do i=1,knode_iblk2(iblk)
          ng_nod=mlgnod2(iblk,i)
          do j=1,kdgof2
            ntemp=ntemp+1
            apool(ntemp)=ubc(ng_nod,j)
          enddo
        enddo
        if(ntemp.ne.knode_iblk2(iblk)*kdgof2) then
          print *,'Error! ntemp ne nextern ...',ntemp,knode_iblk2(iblk)*kdgof2
          call endjob(ierr)
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo
C
C.... Update and Save unodr file
C
      open(21,file='unodr',form='unformatted',status='unknown')
      write(21) ((urg(i,j),j=1,kdgof1),i=1,knode)
      close(21)
C.... Send urg to slave
      do iblk=1,numblk
        ntemp = 0
        do i=1,knode_iblk(iblk)
          ng_nod = mlgnod(iblk,i)
          do j=1,kdgof1
            ntemp = ntemp+1
            apool(ntemp) = urg(ng_nod,j)
          enddo
        enddo
        if(ntemp.ne.knode_iblk(iblk)*kdgof1) then
          print *,'Error! ntemp ne nextern ...',ntemp,knode_iblk(iblk)*kdgof1
          call endjob(ierr)
        endif
        call sendar(iblk,nsource,apool,ntemp)
      enddo
C
c@mend
      END
