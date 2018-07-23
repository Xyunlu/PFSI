      subroutine Mpartition0(filename,kkcom)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      common /jcom/ jcom(100)
      common /ia/ ia(150000000)
      common /ib/ ib(75000000)
      common /jc/ jc(75000000)
      character*22 fname1,fname2
      character*1 material
 
c
c     get the lmddm control file
c
      isignpart=0
      open(21,file='partition.dat',form='formatted',status='unknown')
      read(21,*) numblk,isignpart
      read (21,*) mdivx,mdivy,mdivz
      close(21)
      nparts = numblk
c
      open(22,file='npart',form='formatted',status='unknown')
      write(22,*) numblk
      close(22)
C.......OPEN No. of DIV file
c      open (21,file='div.dat',form='formatted',status='old')
c      read (21,*) mdivx,mdivy,mdivz
c      close(21)
C      write (*,*) MDIVX,MDIVY,MDIVZ
c
c     open data file from GID preprocessing
c
C.....OPEN IoFile
      open(21,file=filename(1),form='formatted',status='old')
      read(21,*) material
c      read(21,*) numtyp
      close(21)
C.... Open Numtyp_bodyFile
      open(21,file='nbefile',form='formatted',status='old')
      read(21,*) numtyp
      close(21)
C.....OPEN ELEM0 FILE
      kelem=0
      numall=0
      OPEN (30,FILE=filename(2),FORM='UNFORMATTED',STATUS='UNKNOWN')
      do ityp=1,numtyp
        READ(30) mnode0,NNODE0
        numall=numall+mnode0
        kelem=kelem+mnode0*nnode0
        if(material.eq.'Y'.or.material.eq.'y') then
          read(30)
        endif
      enddo
C.....Open coor0 file
      open(31,file=filename(3),form='unformatted',status='old')
      read(31) knode,ncoor
      knode1=knode+1
      max_lknode=knode*2/numblk+1000
      max_lkelem=kelem*2/numblk+1000
      max_lnum=numall*2/numblk+1000
      max_lknode3=3*max_lknode
      maxa=75000000-1000-knode-1-100
      maxt=150000000/2-kelem-1000

C
      kna5=kelem*1
      if (kna5/2*2 .lt. kna5) kna5=kna5+1
      knb3=maxa*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knb1=knode1*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb2=1000*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      kna2=maxt*1
      if (kna2/2*2 .lt. kna2) kna2=kna2+1
      kna1=maxt*1
      if (kna1/2*2 .lt. kna1) kna1=kna1+1
      kna3=numtyp*1
      if (kna3/2*2 .lt. kna3) kna3=kna3+1
      kna4=numtyp*1
      if (kna4/2*2 .lt. kna4) kna4=kna4+1
      knc8=knode*1
      if (knc8/2*2 .lt. knc8) knc8=knc8+1
      knd5=numtyp*1
      if (knd5/2*2 .lt. knd5) knd5=knd5+1
      knd6=numtyp*1
      if (knd6/2*2 .lt. knd6) knd6=knd6+1
      knc1=knode*1
      if (knc1/2*2 .lt. knc1) knc1=knc1+1
      knd1=numblk*1
      if (knd1/2*2 .lt. knd1) knd1=knd1+1
      knc5=numblk*1
      if (knc5/2*2 .lt. knc5) knc5=knc5+1
      knc9=numblk*1
      if (knc9/2*2 .lt. knc9) knc9=knc9+1
      knd4=numblk*1
      if (knd4/2*2 .lt. knd4) knd4=knd4+1
      knd2=numblk*max_lknode*1
      if (knd2/2*2 .lt. knd2) knd2=knd2+1
      knc3=numblk*max_lkelem*1
      if (knc3/2*2 .lt. knc3) knc3=knc3+1
      knc4=numblk*max_lnum*1
      if (knc4/2*2 .lt. knc4) knc4=knc4+1
      knc6=numblk*numtyp*1
      if (knc6/2*2 .lt. knc6) knc6=knc6+1
      knc7=numblk*numtyp*1
      if (knc7/2*2 .lt. knc7) knc7=knc7+1
      kne1=max_lknode3*2
      knc11=knode*ncoor*2
      knc2=1000*1
      if (knc2/2*2 .lt. knc2) knc2=knc2+1
      kne2=max_lknode*1
      if (kne2/2*2 .lt. kne2) kne2=kne2+1
      kne3=max_lknode*1
      if (kne3/2*2 .lt. kne3) kne3=kne3+1
      knd3=numblk*max_lknode*1
      if (knd3/2*2 .lt. knd3) knd3=knd3+1
      knc10=knode*1
      if (knc10/2*2 .lt. knc10) knc10=knc10+1
      knd7=knode*1
      if (knd7/2*2 .lt. knd7) knd7=knd7+1
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
      write(*,*) 'memory needed = ',kna5,' in prgram partition'
      stop 55555
      endif
      knb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      if (knb3-1.gt.75000000) then
      write(*,*) 'exceed memory of array ib'
      write(*,*) 'memory of ib = 75000000'
      write(*,*) 'memory needed = ',knb3,' in prgram partition'
      stop 55555
      endif
      knc0=1
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
      if (knc11-1.gt.150000000) then
      write(*,*) 'exceed memory of array ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'memory needed = ',knc11,' in prgram partition'
      stop 55555
      endif
      kcom=kkcom+1
      knd0=jcom(kcom)
      knd1=knd1+knd0
      knd2=knd2+knd1
      knd3=knd3+knd2
      knd4=knd4+knd3
      knd5=knd5+knd4
      knd6=knd6+knd5
      knd7=knd7+knd6
      jcom(kcom+1) = knd7
      if (knd7-1.gt.75000000) then
      write(*,*) 'exceed memory of array jc'
      write(*,*) 'memory of jc = 75000000'
      write(*,*) 'memory needed = ',knd7,' in prgram partition'
      stop 55555
      endif
      kne0=1
      kne1=kne1+kne0
      kne2=kne2+kne1
      kne3=kne3+kne2
      if (kne3-1.gt.75000000) then
      write(*,*) 'exceed memory of array ib'
      write(*,*) 'memory of ib = 75000000'
      write(*,*) 'memory needed = ',kne3,' in prgram partition'
      stop 55555
      endif
      call partition0(knode,knode1,ncoor,nparts,
     *numblk,kelem,max_lknode,max_lknode3,max_lkelem,max_lnum,
     *maxa,maxt,numtyp,material,mdivx,mdivy,
     *mdivz,isignpart,ia(kna0),ia(kna1),ia(kna2),
     *ia(kna3),ia(kna4),ib(knb0),ib(knb1),ib(knb2),
     *ia(knc0),ia(knc1),ia(knc2),ia(knc3),ia(knc4),
     *ia(knc5),ia(knc6),ia(knc7),ia(knc8),ia(knc9),
     *ia(knc10),jc(knd0),jc(knd1),jc(knd2),jc(knd3),
     *jc(knd4),jc(knd5),jc(knd6),ib(kne0),ib(kne1),
     *ib(kne2),
     *filename)
 
      return
      end
      subroutine partition0(knode,knode1,ncoor,nparts,
     *numblk,kelem,max_lknode,max_lknode3,max_lkelem,max_lnum,
     *maxa,maxt,numtyp,material,mdivx,mdivy,
     *mdivz,isignpart,naj,nap,mnode1,nnode1,
     *node,numcol,lm,na,mglnod,nb_elm,
     *node_iblk,mlgelm_no,num_iblk,mnodei,nnodei,idnode,
     *kelem_iblk,nod_ord,coor,knode_iblk,mlgnod,iorder,
     *n_update,mnode,nnode,iord_nod,coor_iblk,itemp,
     *iord,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
c====================================================================
c     idnode:the index of the node belone which iblk
c     knode_iblk:the number of the node in every subdomain
c     N_update:the number of the internal nodes in every subdomain
c     num_iblk:the number of the elements in every subdomain
c     mglnod:the map of node number from global to local
c     mlgelm_no:the map of elem number from local to global
c     node_iblk:the local element array
c     mlgnod:the map of node number from local to global(include the
c            external nodes)
c=====================================================================
      character*22 fname1,fname2
      DIMENSION NODE(kelem),na(maxa),numcol(knode1),lm(1000)
      dimension nap(maxt),naj(maxt)
      dimension mnode1(numtyp),nnode1(numtyp)
      dimension IDNODE(knode)
      dimension mnode(numtyp),nnode(numtyp)
      dimension nmnode(numtyp),nnnode(numtyp)
      dimension mglnod(knode)
      dimension knode_iblk(numblk),num_iblk(numblk),kelem_iblk(numblk)
      dimension N_update(numblk),mlgnod(numblk,max_lknode)
      dimension node_iblk(numblk,max_lkelem),mlgelm_no(numblk,max_lnum)
      dimension mnodei(numblk,numtyp),nnodei(numblk,numtyp)
      dimension coor_iblk(max_lknode3),coor(knode,ncoor)
      dimension nb_elm(1000)
      dimension nknode_iblk(numblk)
      dimension nN_update(numblk),nmlgnod(numblk,max_lknode)
      dimension itemp(max_lknode),iord(max_lknode)
      dimension iorder(numblk,max_lknode)
      dimension norder(numblk,max_lknode)
      dimension nod_ord(knode),iord_nod(knode),niord_nod(knode)
      CHARACTER*1 MATERIAL
      character*22 fname,GraphFile
      character*1000 licensefile
      include 'partdata.h'
 
      dimension indx(mdivx+1),indy(mdivx*mdivy+1)
      dimension ind(mdivx*mdivy*mdivz+1),indz(mdivx*mdivy*mdivz+1)

      integer, dimension(:,:), allocatable :: id
 
6     FORMAT (1X,10I5)
7     FORMAT (1X,6E12.5)
8     FORMAT (1X,1000I12)
 
      nsource = 0
c      print *,'111  GraphFile = ',GraphFile
      if(nparts.ne.numblk) then
        write(*,*)'Error when get partition information of the dimain'
        stop 1111
      end if
 
      do ityp=1,numtyp
         mnode(ityp)=0
         nnode(ityp)=0
      enddo
 
C.......OPEN ELEM0 file
      rewind(30)
      kelem1=0
      do ityp=1,numtyp
         read(30) mnode1(ityp),nnode1(ityp),
     *            ((node(kelem1+(i-1)*nnode1(ityp)+j),
     *            j=1,nnode1(ityp)),i=1,mnode1(ityp))
c         print *,'ityp == ',ityp
c         do ne=1,mnode1(ityp)
c            print *,(node(kelem1+(ne-1)*nnode(ityp)+j),j=1,nnode(ityp))
c         enddo
         kelem1=kelem1+mnode1(ityp)*nnode1(ityp)
         if(material.eq.'y'.or.material.eq.'Y') then
           read(30)
         endif
      enddo
      close(30)
c
c     isignpart is the sign of the method of partition
c
      if(isignpart.eq.1) then
C========================  Get the Graph of the Mesh ============
        DO 350 I=1,knode1
          numcol(i)=0
350     CONTINUE
        do i=1,maxt
          nap(i)=0
          naj(i)=0
        enddo
        jna=knode
        NUMEL=0
        kelem=0
        kemate=0
        DO 2000 ITYP=1,NUMTYP
C.......INPUT ENODE
          num=mnode1(ityp)
          nod=nnode1(ityp)
          NNe=nod
          IF(MATERIAL.EQ.'Y' .OR. MATERIAL.EQ.'y') THEN
            NNE = NNE-1
          ENDIF
          DO 1000 NE=1,NUM
           L=0
           DO 700 INOD=1,NNE
            NODI=NODE(kelem+(NE-1)*NOD+INOD)
            L=L+1
            LM(L)=NODI
700        CONTINUE
           NUMEL=NUMEL+1
           CALL partition0_1(knode,NUMCOL,NAP,NAJ,L,LM,JNA)
           IF (JNA.GT.MAXT) THEN
             WRITE(*,*) 'EXCEET ARRAY LENGTH MAXT ....',MAXT,' < ',JNA
             STOP 1111
           ENDIF
1000      continue
          kelem=kelem+num*nod
2000    CONTINUE
        CALL partition0_2(knode,NAP,NAJ,NUMCOL,NA,LM)
c      write(*,*) 'nvtex,nedges == ',knode,nedges
c      write(*,*) 'na =='
c      write(*,*) (na(i),i=1,nedges)
        nedges=NUMCOL(knode+1)/2
        call pmetismain(nparts,idnode,knode,nedges,numcol,na)
      endif
      do ityp=1,numtyp
        mnode(ityp)=mnode1(ityp)
        nnode(ityp)=nnode1(ityp)
      enddo
c
c     read coor0 file .......
c
      rewind(31)
      read(31) knode,ncoor,((coor(i,j),j=1,ncoor),i=1,knode)
      close(31)
      if(isignpart.eq.2) then
        call partcoor(knode,mdivx,mdivy,mdivz,coor(1,1),coor(1,2),
     *    coor(1,3),mglnod,nod_ord,knode_iblk,idnode,indx,indy,indz,ind)
      endif
 
      open(21,file='id0',form='unformatted',status='unknown')
      read(21) knode,kdgof
      allocate(id(knode,kdgof))
      rewind(21)
      read(21) knode,kdgof,((id(i,j),j=1,kdgof),i=1,knode)
      close(21)

c      print *,'28056, id =',(id(28056,j),j=1,kdgof)
c      print *,'85819, id =',(id(85819,j),j=1,kdgof)

      do i=1,knode
        do j=1,kdgof
          N=id(i,j)
          if( N .lt. -1) then
            inod = -N-1
            idflg = idnode(inod)
            exit
          endif
        enddo
      enddo

c      print *,'inod, idflg = ',inod,idnode(inod)

      do i=1,knode
        do j=1,kdgof
          N=id(i,j)
          if( N .lt. -1) then
            inod = -N-1
c            idnode(inod) = idflg
            idnode(i) = idnode(inod)
            exit
          endif
        enddo
      enddo

c      print *,'idnode(51) = ', idnode(51)

      deallocate(id)

      open(21,file='unod_id',form='unformatted',status='unknown')
      write(21) (dble(idnode(i)+1),i=1,knode)
      close(21)
 
c
c     call partpost function now
c
      call partpost(knode,ncoor,numblk,
     *kelem,max_lknode,max_lknode3,max_lkelem,max_lnum,maxa,
     *numtyp,material,licensefile,mglnod,
     *nb_elm,node_iblk,mlgelm_no,num_iblk,mnodei,nnodei,
     *idnode,kelem_iblk,nod_ord,knode_iblk,n_update,mlgnod,
     *iorder,mnode,nnode,iord_nod,itemp,iord,node)

      open(21,file='time0',form='formatted',status='old')
      read(21,*) t0,tmax,dt
      close(21)
c.....open time file
      open(44,file='time',form='unformatted',status='unknown')
      write(44) tmax,dt,t0,0
      close(44)

c
c     get the local coor for subdomain
c
      do iblk=1,numblk
        do nc=1,ncoor
           do inod=1,knode_iblk(iblk)
              nod_g=mlgnod(iblk,inod)
c              if(nod_g .eq. 75943 .or. nod_g .eq. 85819) then
c                write(*,'(a,12i7)') 'nod_g,iblk,inod,iorder(nod_g)=',nod_g,iblk,inod,iorder(iblk,inod)
c              endif
              coor_iblk((nc-1)*knode_iblk(iblk)+inod)=coor(nod_g,nc)
           enddo
        enddo
        call sendint(iblk,nsource,numblk)
        call sendr(iblk,nsource,t0)
        call sendr(iblk,nsource,tmax)
        call sendr(iblk,nsource,dt)
        call sendint(iblk,nsource,knode_iblk(iblk))
        call sendint(iblk,nsource,ncoor)
        call sendar(iblk,nsource,coor_iblk,knode_iblk(iblk)*ncoor)
        do i=1,knode_iblk(iblk)
           itemp(i)=iorder(iblk,i)
        enddo
        call sendai(iblk,nsource,itemp,knode_iblk(iblk))
      enddo
c
c----------------------------------------------------
      return
      end
 
      SUBROUTINE partition0_1(NEQ,NUMCOL,NAP,NAJ,ND,LM,JNA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NAP(*),NAJ(*),NUMCOL(*),LM(*)
6     FORMAT(1X,5I15)
C      WRITE (*,*) 'ND= ',ND, (LM(I),I=1,ND)
      DO 400 I=1,ND
      NI = LM(I)
      DO 300 J=1,ND
      NJ = LM(J)
      IF (NJ.EQ.NI) GOTO 300
      NUMJ = NUMCOL(NI)
      IF (NUMJ.EQ.0) THEN
       JNA = JNA+1
       NAP(NI) = JNA
       NAJ(NI) = NJ
       NUMCOL(NI) = NUMCOL(NI)+1
      ELSE
       JP = NAP(NI)
       JV = NAJ(NI)
       IF (NJ.EQ.JV) GOTO 300
       DO K=1,NUMJ-1
        JV = NAJ(JP)
        JP = NAP(JP)
        IF (NJ.EQ.JV) GOTO 300
       ENDDO
       JNA = JNA+1
       NAP(JP) = JNA
       NAJ(JP) = NJ
       NUMCOL(NI) = NUMCOL(NI)+1
      ENDIF
300   CONTINUE
400   CONTINUE
      RETURN
      END
 
      SUBROUTINE partition0_2(NEQ,NAP,NAJ,NUMCOL,NA,LMI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NUMCOL(*),NAP(*),NAJ(*),NA(*),LMI(*)
C        IF(NUMEL.EQ.0) GO TO 1000
      NN = 0
      DO 600 N=1,NEQ
      JP = NAP(N)
      JV = NAJ(N)
      LI = NUMCOL(N)
      DO 500 I=1,LI
      LMI(I) = JV
      JV = NAJ(JP)
      JP = NAP(JP)
500   CONTINUE
      CALL partition0_3(LI,LMI)
      DO 550 I=1,LI
      NN = NN+1
550   NA(NN) = LMI(I)
600   CONTINUE
      DO 800 N=1,NEQ-1
800   NUMCOL(N+1) = NUMCOL(N+1)+NUMCOL(N)
      DO 850 N=1,NEQ
850   NUMCOL(NEQ-N+2) = NUMCOL(NEQ-N+1)
      NUMCOL(1) = 0
C       WRITE(*,*) 'NUMCOL ='
C       WRITE(*,6) (NUMCOL(N),N=1,NEQ+1)
C       WRITE(*,*) 'NA ='
C       WRITE(*,6) (NA(I),I=1,NN)
1000    RETURN
6       FORMAT(1X,5I15)
      RETURN
      END
 
      SUBROUTINE partition0_3(ND,LM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LM(1)
C       WRITE(*,*) '**** ORDER ****'
C       WRITE(*,*) (LM(I),I=1,ND)
      DO 200 I=1,ND
      LS=LM(I)+1
      DO 100 J=I,ND
      IF (LM(J).GT.LS) GOTO 100
      LS=LM(J)
      J0=J
100     CONTINUE
      LM(J0)=LM(I)
      LM(I)=LS
200     CONTINUE
C       WRITE(*,*) (LM(I),I=1,ND)
C       WRITE(*,*) '-----------------'
      RETURN
      END
 
      subroutine partition0_4(kdgof,knode,id)
      dimension id(kdgof,knode),ms(1000),is(1000)
      do 1000 k=1,kdgof
      m = 0
      do 800 n=1,knode
      if (id(k,n).le.1) goto 800
       j=id(k,n)
       j0=0
       if (m.gt.0) then
        do i=1,m
         if (j.eq.ms(i)) j0=is(i)
        enddo
       endif
       if (j0.eq.0) then
        m=m+1
        ms(m)=j
        is(m)=n
        id(k,n)=1
       else
        id(k,n)=-j0-1
       endif
800   continue
1000  continue
      return
      end
 
