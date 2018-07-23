      subroutine Mbftm(filename,kkcom)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      common /jcom/ jcom(100)
      logical filflg
      character*12 fname1,fname2
      character*12 nvfile(100),bfdfile(100),disfile(100)
      dimension ndof(100)
      common /ia/ ia(150000000)
      common /jc/ jc(75000000)
      common /stop/ msstop
      KCOOR=3
      numnf=1
      maxdof=0
C ....... field a .................................
      open(11,file='id0',form='unformatted',status='old')
      read(11) nnode,kdgof
      close(11)
      nvfile(1)='id0'
      bfdfile(1)='bfda'
      disfile(1)='disp0'
      if (kdgof.gt.maxdof) maxdof=kdgof
      ndof(1)=kdgof
      KNODE = NNODE
      KVAR=KNODE*KDGOF
C................................................................
      open(21,file='partition.dat',form='formatted',status='old')
      read(21,*) numblk
      close(21)
c.................
      max_lknode=knode*2/numblk+1000
c.....................
      kna4=knode*maxdof*1
      if (kna4/2*2 .lt. kna4) kna4=kna4+1
      kna1=knode*kcoor*2
      kna2=knode*maxdof*2
      kna7=3*2
      kna5=knode*maxdof*1
      if (kna5/2*2 .lt. kna5) kna5=kna5+1
      kna3=maxdof*2
      kna6=knode*1
      if (kna6/2*2 .lt. kna6) kna6=kna6+1
      knb1=numblk*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb4=numblk*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb2=numblk*max_lknode*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=numblk*max_lknode*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      jcom(1) = 1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
      kna7=kna7+kna6
      if (kna7-1.gt.150000000) then
      write(*,*) 'exceed memory of array ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'memory needed = ',kna7,' in prgram bft'
      stop 55555
      endif
      kcom=kkcom+1
      knb0=jcom(kcom)
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      knb4=knb4+knb3
      jcom(kcom+1) = knb4
      if (knb4-1.gt.75000000) then
      write(*,*) 'exceed memory of array jc'
      write(*,*) 'memory of jc = 75000000'
      write(*,*) 'memory needed = ',knb4,' in prgram bft'
      stop 55555
      endif
      call bftm_1(knode,kdgof,kcoor,numblk,
     *tmax,max_lknode,time,it,maxdof,numnf,
     *ndof,nvfile,bfdfile,disfile,msstop,ia(kna0),ia(kna1),
     *ia(kna2),ia(kna3),ia(kna4),ia(kna5),ia(kna6),
     *jc(knb0),jc(knb1),jc(knb2),jc(knb3),
     *filename)
 
C ...... OPEN THE FILE TO OBTAIN GRAPH FILE NAMES
      inquire(file='plotname',exist=filflg)
      if (.NOT. filflg) then
      fname1 = 'unod'
      open(26,file='plotname',form='formatted',status='unknown')
      write(26,'(8a)') fname1
      close(26)
      endif
      open(26,file='plotname',form='formatted',status='old')
C ...... OPEN THE BATCH FILE FOR STORING THE RESULT TO GRAPHIC
      open(27,file='post.bat',form='formatted',status='unknown')
 
C ...... STORE THE RESULT FOR EACH nstep TIME STEP
      nstep = 1
      ik = it/nstep
      kk = it-ik*nstep
cc      write(*,*) 'nstep,it,ik,kk =',nstep,it,ik,kk
      if (kk.gt.0) goto 9999
9998  CONTINUE
 
C ...... GET THE GRAPHIC FILE NAME
      read(26,'(8a)',END=9999) fname1
C       WRITE(*,*) 'fname1 =',fname1
      fname2 = fname1
      call bftm_2(fname2,ik)
C ...... WRITE COPY COMMAND TO post.bat FILE FOR STORING THE RESULT
      write(27,*) 'cp -f ',fname1,' ',fname2
      GOTO 9998
9999  CONTINUE
      close(26)
      write(27,*)
      close(27)
c       write(*,*) 'it =',it
c       write(*,*) fname1
c       write(*,*) fname2
 
 
      return
      end
      subroutine bftm_1(knode,kdgof,kcoor,numblk,
     *tmax,max_lknode,time,it,maxdof,numnf,
     *ndof,nvfile,bfdfile,disfile,msstop,coor,
     *bf,bfi,nodvar,id,inode,r,
     *knodei,mlgnod,iorder,nupdatei,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
      DIMENSION NODVAR(KNODE,MAXDOF),COOR(KNODE,kcoor),
     *  BF(KNODE,MAXDOF),R(3),id(KNODE,MAXDOF),
     *  bfi(MAXDOF*knode),inode(knode)
      dimension knodei(numblk),nupdatei(numblk)
      dimension mlgnod(numblk,max_lknode)
      dimension iorder(numblk,max_lknode)
      character*12 nvfile(100),bfdfile(100),disfile(100)
 
      integer, dimension(:), allocatable :: nvartmp
 
      msstop = 0
      nsource = 0
C.......OPEN TIME File AND UPDATE THE TIME
      OPEN(21,FILE=filename(1),FORM='UNFORMATTED')
      READ(21) TMAX,DT,TIME,IT
      T = TIME+DT
      TIME = TIME+DT
      IT = IT+1
cc        WRITE(*,*) ' TMAX,DT,TIME,IT =',TMAX,DT,TIME,IT
      REWIND(21)
      WRITE(21) TMAX,DT,TIME,IT
      CLOSE(21)
 
C.......OPEN COOR file
      OPEN (21,FILE='coor0',FORM='UNFORMATTED',STATUS='OLD')
      READ (21) KNODE,NCOOR,((COOR(I,J),J=1,NCOOR),I=1,KNODE)
      CLOSE(21)
 
C.......OPEN STOP file IF THE LAST TIME IS ARRIVED
      IF (TIME-TMAX.GT.-1.0d-20) THEN
        msstop = 1
      ENDIF
 
      open(21,file='inlet.id',form='unformatted',status='old')
      read(21) mmm
      allocate(nvartmp(mmm))
      rewind(21)
      read(21) mmm,(nvartmp(i),i=1,mmm)
      close(21)
 
      open(21,file='init.dat',status='old')
      read(21,*) xmax,ymax,zmax
      read(21,*) skip
      read(21,*) skip
      read(21,*) u_in
      read(21,*) skip
      read(21,*) skip
      read(21,*) skip
      read(21,*) id4updown
      read(21,*) index4pro
      close(21)
      pi=datan(1.d0)*4.d0
 
      do 1100 nf=1,numnf
C.......  OPEN id0 file
        OPEN (22,FILE=trim(nvfile(nf)),FORM='UNFORMATTED',STATUS='OLD')
        READ (22) KNODE,KDGOF,((id(I,J),J=1,KDGOF),I=1,KNODE)
        CLOSE (22)
        open(21,file=trim(disfile(NF)),form='unformatted',status='old')
        read(21) knode,kdgof,((bf(j,i),i=1,kdgof),j=1,knode)
        close(21)
c        do j=1,kdgof
c          do i=1,knode
c            nodvar(i,j)=id(i,j)
c          enddo
c        enddo
c        neq=0
c        do i=1,knode
c          inod=nod_ord(i)
c          do j=1,kdgof
c            IF (id(inod,j).gt.0) then
c              NEQ = NEQ + 1
c              nodvar(inod,j) = NEQ
c            endif
c          enddo
c        enddo
c        neq1=neq+1
c        do i=1,knode
c          do j=1,kdgof
c            IF (nodvar(i,j).LT.-1) then
c              N = -nodvar(i,j)-1
c27            CONTINUE
c              IF (nodvar(i,j).LT.-1) THEN
c                N=-nodvar(i,j)-1
c                GOTO 27
c              ENDIF
c              nodvar(i,j) = NODVAR(n,j)
c            endif
c          enddo
c        enddo
C        WRITE(*,*) 'KNODE =',KNODE,' KDGOF =',KDGOF
C        WRITE (*,*) 'NODVAR ='
C        WRITE (*,6) ((NODVAR(I,J),J=1,KDGOF),I=1,KNODE)
       nchk = 0
 
C......COMPUTE BOUNDARY CONDITION
        DO 333 N=1,KNODE
          DO 100 J=1,NCOOR
100         R(J) = COOR(N,J)
          DO 200 J=1,KDGOF
            IDm = ID(N,J)
c
ccc Inlet ID
c
            if(nvartmp(n) .eq. 0 .and. j .lt. kdgof) idm=-2

            if( N .eq. nchk) then
              print *,'nchk,j, Idm =',nchk, j,idm
              print *,'bf =', bf(nchk,j)
            endif

            IF(IDm.eq.-2) then
              if (id4updown.eq.0) then
                if (index4pro.eq.0) then
                  BF(N,j) = bound(r,u_in,xmax,zmax,time,pi,j)
                else
                  BF(N,j) = boundpro(r,u_in,xmax,zmax,time,pi,j)
                endif
              elseif (id4updown.eq.1) then
                r(3)=r(2)
                if (index4pro.eq.0) then
                  BF(N,j) = -bound(r,u_in,xmax,zmax,time,pi,j)
                else
                  BF(N,j) = -boundpro(r,u_in,xmax,zmax,time,pi,j)
                endif
              endif
            endif
200       CONTINUE
c
          if (id4updown .eq. 1) then
            tmp=bf(n,2)
            bf(n,2)=bf(n,3)
            bf(n,3)=tmp
          endif
333     CONTINUE
CC        WRITE(*,*) ' BF = '
CC        WRITE(*,'(6F13.3)') ((BF(N,J),J=1,KDGOF),N=1,KNODE)
C.....    partition the boundary value
        do 1200 iblk=1,numblk
          call sendint(iblk,nsource,msstop)
          call sendr(iblk,nsource,tmax)
          call sendr(iblk,nsource,dt)
          call sendr(iblk,nsource,time)
          call sendint(iblk,nsource,it)
          nod_iblk=knodei(iblk)
          do i=1,nod_iblk
            do j=1,kdgof
              bfi((j-1)*nod_iblk+i)=bf(mlgnod(iblk,i),j)
            enddo
          enddo
          call sendar(iblk,nsource,bfi,nod_iblk*kdgof)
c..............................................................
1200    continue
C.......OPEN BFD file and WRITE BOUNDARY CONDITION
        OPEN (22,FILE=bfdfile(nf),FORM='UNFORMATTED',STATUS='unknown')
        WRITE(22) ((BF(I,J),I=1,KNODE),J=1,KDGOF)
        CLOSE (22)
 
1100  continue
      deallocate(nvartmp)
 
      END
 
      subroutine bftm_2(str,num)
      implicit real*8(a-h,o-z)
      character*25 str, ext
 
      write(ext,'(i5)') num
      str = trim(str) // '.' // trim(adjustl(ext))
 
      return
      end
