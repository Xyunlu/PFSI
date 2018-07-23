      subroutine Mgetpart(filename,kkcom)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      common /jcom/ jcom(100)
      character*12 fname1,fname2
      character*1 material
      common /ia/ ia(150000000)
      common /ib/ ib(75000000)
      common /jc/ jc(75000000)
 
c     Get the element data file and domain partition number'
c
c     Check if the file name and the partition number are correct !
c
c      write(*,*)'The input partition filname and number of subdomains:'
c      write(*,*) filename(1), nparts
c
c     get the lmddm control file
c
c      lgio=0
      open(1,file='partition.dat',form='formatted',status='unknown')
      read(1,*) numblk
      close(1)
c
      nparts = numblk
C.....OPEN IoFile to get material and element type number ..
      open(21,file=filename(1),form='formatted',status='old')
      read(21,*) material
      read(21,*) numtyp
      close(21)
      print *,'numtyp == ',numtyp
C.... Open Numtyp_bodyFile
c
c     need to get numtyp_body here
c
      open(21,file='nbefile',form='formatted',status='old')
      read(21,*) numtyp_body
      close(21)
      print *,'numtyp_body == ',numtyp_body
c
c
      open(2,file='part',form='formatted',status='unknown')
      write(2,*) numblk
      close(2)
c
      nsource=0
C
C.....OPEN ELEM0 FILE
c     get the some constant to malloc memory
c
      OPEN (30,FILE=filename(2),FORM='UNFORMATTED',STATUS='UNKNOWN')
      kelem_body=0
      kemate_body=0
      numall=0
      do itype=1,numtyp_body
        READ(30) NUM0,NNODE0
        numall=numall+num0
        kelem_body=kelem_body+num0*nnode0
        if(material.eq.'y'.or.material.eq.'Y') then
          READ(30) MMATE0,NMATE0
          kemate_body=kemate_body+mmate0*nmate0
        endif
      enddo
      if(kemate_body.eq.0) kemate_body=1
c
      kelem_bound=0
      kemate_bound=0
c      numallbd=0
      do itype=numtyp_body+1,numtyp
        READ(30) NUM0,NNODE0
c        numallbd=numallbd+num0
        kelem_bound=kelem_bound+num0*nnode0
        if(material.eq.'y'.or.material.eq.'Y') then
          READ(30) MMATE0,NMATE0
          kemate_bound=kemate_bound+mmate0*nmate0
        endif
      enddo
      if(kemate_body.eq.0) kemate_body=1
C.....OPNE DISP0 FILE (Boundary condition file)
      OPEN (33,FILE=filename(3),FORM='UNFORMATTED',STATUS='UNKNOWN')
      READ(33) KNODE,KDGOF
c
      max_lknode=knode*2/numblk+1000
      max_lnum=numall*2/numblk+1000
      max_lkelem=kelem_body*2/numblk+1000
c
      maxip=max_lknode*kdgof
      if(maxip.lt.max_lkelem) maxip=max_lkelem
      maxap=max_lknode*kdgof
c      maxap2=maxap*2
c
      kna10=knode*kdgof*1
      if (kna10/2*2 .lt. kna10) kna10=kna10+1
      knb4=knode*kdgof*2
      kna12=knode*kdgof*1
      if (kna12/2*2 .lt. kna12) kna12=kna12+1
      kna13=max_lknode*1
      if (kna13/2*2 .lt. kna13) kna13=kna13+1
      kna9=kemate_body*2
      kna16=kelem_bound*1
      if (kna16/2*2 .lt. kna16) kna16=kna16+1
      kna11=kemate_bound*2
      knc5=numtyp_body*1
      if (knc5/2*2 .lt. knc5) knc5=knc5+1
      knc6=numtyp_body*1
      if (knc6/2*2 .lt. knc6) knc6=knc6+1
      kna14=numtyp*1
      if (kna14/2*2 .lt. kna14) kna14=kna14+1
      kna15=numtyp*1
      if (kna15/2*2 .lt. kna15) kna15=kna15+1
      knb1=maxip*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb2=maxap*2
      kna3=numblk*max_lkelem*1
      if (kna3/2*2 .lt. kna3) kna3=kna3+1
      kna4=numblk*max_lnum*1
      if (kna4/2*2 .lt. kna4) kna4=kna4+1
      knc1=numblk*1
      if (knc1/2*2 .lt. knc1) knc1=knc1+1
      knc4=numblk*1
      if (knc4/2*2 .lt. knc4) knc4=knc4+1
      knc2=numblk*max_lknode*1
      if (knc2/2*2 .lt. knc2) knc2=knc2+1
      kna6=numblk*numtyp_body*1
      if (kna6/2*2 .lt. kna6) kna6=kna6+1
      kna17=kelem_bound*1
      if (kna17/2*2 .lt. kna17) kna17=kna17+1
      kna7=numblk*numtyp_body*1
      if (kna7/2*2 .lt. kna7) kna7=kna7+1
      kna2=1000*1
      if (kna2/2*2 .lt. kna2) kna2=kna2+1
      kna5=numblk*1
      if (kna5/2*2 .lt. kna5) kna5=kna5+1
      kna8=knode*1
      if (kna8/2*2 .lt. kna8) kna8=kna8+1
      kna1=knode*1
      if (kna1/2*2 .lt. kna1) kna1=kna1+1
      knc3=numblk*max_lknode*1
      if (knc3/2*2 .lt. knc3) knc3=knc3+1
      kna18=numtyp*1
      if (kna18/2*2 .lt. kna18) kna18=kna18+1
      kna19=numtyp*1
      if (kna19/2*2 .lt. kna19) kna19=kna19+1
      kna20=numblk*numtyp*1
      if (kna20/2*2 .lt. kna20) kna20=kna20+1
      kna21=numblk*numtyp*1
      if (kna21/2*2 .lt. kna21) kna21=kna21+1
      knb3=maxap*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knc7=knode*1
      if (knc7/2*2 .lt. knc7) knc7=knc7+1
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
      write(*,*) 'memory needed = ',kna21,' in prgram getpart'
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
      write(*,*) 'memory needed = ',knb4,' in prgram getpart'
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
      jcom(kcom+1) = knc7
      if (knc7-1.gt.75000000) then
      write(*,*) 'exceed memory of array jc'
      write(*,*) 'memory of jc = 75000000'
      write(*,*) 'memory needed = ',knc7,' in prgram getpart'
      stop 55555
      endif
      call getpart(kkcom,knode,kdgof,numtyp,
     *numtyp_body,nparts,nsource,kemate_body,kelem_bound,t0,
     *tmax,dt,kemate_bound,numblk,max_lknode,max_lkelem,
     *max_lnum,material,maxip,maxap,ia(kna0),ia(kna1),
     *ia(kna2),ia(kna3),ia(kna4),ia(kna5),ia(kna6),
     *ia(kna7),ia(kna8),ia(kna9),ia(kna10),ia(kna11),
     *ia(kna12),ia(kna13),ia(kna14),ia(kna15),ia(kna16),
     *ia(kna17),ia(kna18),ia(kna19),ia(kna20),ib(knb0),
     *ib(knb1),ib(knb2),ib(knb3),jc(knc0),jc(knc1),
     *jc(knc2),jc(knc3),jc(knc4),jc(knc5),jc(knc6),
     *filename)
      return
      end
      subroutine getpart(kkcom,knode,kdgof,numtyp,
     *numtyp_body,nparts,nsource,kemate_body,kelem_bound,t0,
     *tmax,dt,kemate_bound,numblk,max_lknode,max_lkelem,
     *max_lnum,material,maxip,maxap,mglnod,nb_elm,
     *node_iblk,mlgelm_no,num_iblk,mnodeboi,nnodeboi,idnode,
     *emate_body,id,emate_bound,nodvar,inode,mmate,
     *nmate,node_bound,nodei_bound,mnode,nnode,mnodei,
     *nnodei,ipool,apool,iapool,u012,knode_iblk,
     *mlgnod,iorder,nupdatei,mnodebo,nnodebo,nod_ord,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Get id0,disp0,disp1,disp2,disp3,elem0 and lgnv information
c     for one field in every subdomain.
c     Parameter list:
c       ID:          array use to read id0 file for the whole domain.
c       U0:          array use to read disp1 file for the whole domain.
c       U1:          array use to read disp2 file for the whole domain.
c       U2:          array use to read disp3 file for the whole domain.
c       BFU:         array use to read disp0 file for the whole domain.
c       NODVAR:      array use to save global equation numbering
c                    of the whole domain.
c       inode:       temporary array use to get global node numbering
c                    of every subdomain.
c       emate_body:  array use to save material parameter of body element
c                    of the whole domain.
c       node_bound:  array use to save the global node number of boundary
c                    element for the whole domain.
c       emate_bound: array use to save material parameter of boundary element
c                    of the whole domain.
c       mnode:       array use to save the element total number of every type
c                    of the whole domain.
c       nnode:       array use to save the node number of every type element
c                    of the whole domain if not exist material parameter,
c                    the node number +1 of every type element
c                    of the whole domain if exist material parameter.
c       mmate:       array use to save the total number of material
c                    of every type element of the whole domain.
c       nmate:       array use to save the number of material parameters
c                    of every type element of the whole domain.
c       ipool:       temporary integer array to send messages.
c       apool:       temporary real*8 array to send messages.
c       knodei:      array use to save the total node number in iblk subdomain,
c                    which come from partition() program.
c       nupdatei:    array use to save the total number of update nodes
c                    in iblk subdomain,which come from partition() program.
c       mlgnod:      array use to save the global node numbers of all nodes
c                    at every subdomian, which come from partition() program.
c       node_iblk:       array use to save the local node numbers of every element
c                    at every subdomian, which come from partition() program.
c       mlgelm_no:   array use to save the global element number of all elements
c                    at every subdomian, which come from partition() program.
c       mnodei:      array use to save the element total number of every type
c                    at every subdomain.
c       nnodei:      array use to save the node number of every type element
c                    at every subdomain if not exist material parameter,
c                    the node number +1 of every type element
c                    at every subdomain if exist material parameter.
c       num_iblk:        temporary array
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      DIMENSION  ID(KNODE,KDGOF),
     &  U012(KNODE,KDGOF),NODVAR(KNODE,KDGOF),
     &  inode(max_lknode),
     &  emate_body(kemate_body),
     &  node_bound(kelem_bound),emate_bound(kemate_bound),
     &  mnodebo(numtyp_body),nnodebo(numtyp_body),
     &  mmate(numtyp),nmate(numtyp),
     &  ipool(maxip),apool(maxap)
      dimension node_iblk(numblk,max_lkelem),mlgelm_no(numblk,max_lnum)
      dimension knode_iblk(numblk),nupdatei(numblk)
      dimension mlgnod(numblk,max_lknode)
      dimension mnodeboi(numblk,numtyp_body)
      dimension nodei_bound(kelem_bound),nnodeboi(numblk,numtyp_body)
      dimension nb_elm(1000),num_iblk(numblk)
      dimension idnode(knode),mglnod(knode)
      dimension iorder(numblk,max_lknode)
c         temparary array
      dimension mnode(numtyp),nnode(numtyp)
      dimension mnodei(numblk,numtyp),nnodei(numblk,numtyp)
      dimension iapool(maxap)
      dimension nod_ord(knode)
      character*12 fname1,fname2
      character*3 uname(20)
      character*1 material
      logical filflgdisp(20),filflg
c
c     set initial value of mnode and nnode
c
      do ityp=1,numtyp_body
        mnode(ityp)=mnodebo(ityp)
      enddo
      do ityp=1,numtyp_body
        nnode(ityp)=nnodebo(ityp)
      enddo
      do ityp=1,numtyp_body
        do iblk=1,numblk
          mnodei(iblk,ityp)=mnodeboi(iblk,ityp)
        enddo
      enddo
      do ityp=1,numtyp_body
        do iblk=1,numblk
          nnodei(iblk,ityp)=nnodeboi(iblk,ityp)
        enddo
      enddo
      do ityp=numtyp_body+1,numtyp
        do iblk=1,numblk
          mnodei(iblk,ityp)=0
        enddo
      enddo
      do ityp=numtyp_body+1,numtyp
        do iblk=1,numblk
          nnodei(iblk,ityp)=0
        enddo
      enddo
c
      do ityp=1,numtyp
        mmate(ityp)=0
      enddo
      do ityp=1,numtyp
        nmate(ityp)=0
      enddo
c
c...  initial finished .........c
c
      if(nparts.ne.numblk) then
        write(*,*)'Error when get partition information of the dimain'
        stop 1111
      end if
 
      rewind(30)
c
c     read body element from elem0 file firstly
c
      kemate=0
      do ityp=1,numtyp_body
        READ (30)
        if(material.eq.'y'.or.material.eq.'Y') then
          READ (30) MMATE(ityp),NMATE(ityp),
     *        ((EMATE_body((I-1)*NMATE(ityp)+J+kemate),
     *                  J=1,NMATE(ityp)),I=1,MMATE(ityp))
c          WRITE(*,*) 'MMATE =',MMATE(ityp),' NMATE =',NMATE(ityp)
c          WRITE (*,*) 'EMATE ='
c          WRITE (*,*) ((EMATE_body((I-1)*NMATE(ityp)+J+kemate),
c     *          J=1,NMATE(ityp)),I=1,MMATE(ityp))
          kemate=kemate+mmate(ityp)*nmate(ityp)
        endif
      enddo
c
c     read boundary element from elem0 file
c
      kelem=0
      kemate=0
      do ityp=numtyp_body+1,numtyp
        read (30) mnode(ityp),nnode(ityp),
     *         ((node_bound((i-1)*nnode(ityp)+j+kelem),
     *             j=1,nnode(ityp)),i=1,mnode(ityp))
c        WRITE(*,*) 'THE FOLLOWING IS BUNDARY ELEMENT OF TYPE ',ityp
c        WRITE(*,*) 'NUM =',MNODE(ityp),' nnode =',nnode(ityp)
c        WRITE(*,*) 'NODE ='
c        WRITE(*,*) ((NODE_BOUND((I-1)*nnode(ityp)+J+kelem),
c     *             J=1,nnode(ityp)),I=1,mnode(ityp))
        kelem=kelem+mnode(ityp)*nnode(ityp)
c        read boundary material from elem file
        if(material.eq.'y'.or.material.eq.'Y') then
          read (30) mmate(ityp),nmate(ityp),
     *         ((emate_bound((i-1)*nmate(ityp)+j+kemate),
     *             j=1,nmate(ityp)),i=1,mmate(ityp))
c          WRITE(*,*) 'MMATE =',MMATE(ityp),' NMATE =',NMATE(ityp)
c          WRITE (*,*) 'EMATE ='
c          WRITE (*,*) ((EMATE_BOUND((I-1)*NMATE(ityp)+J+kemate),
c     *          J=1,NMATE(ityp)),I=1,MMATE(ityp))
          kemate=kemate+mmate(ityp)*nmate(ityp)
        endif
      enddo
      close(30)
c
c     get the boundary elem number from local to global in iblk subdomain
c
      do iblk=1,numblk
         num_iblk(iblk)=0
         do ityp=1,numtyp_body
           num_iblk(iblk)=num_iblk(iblk)+mnodei(iblk,ityp)
         enddo
      enddo
      kelem=0
      do ityp=numtyp_body+1,numtyp
        num=mnode(ityp)
        nnod=nnode(ityp)
        nne=nnod
        if(material.eq.'y'.or.material.eq.'Y') then
           nne=nnod-1
        endif
        do ne=1,num
          do j = 1,nne
             nb_elm(j)=0
          enddo
          nblk_ne=0
          do j = 1,nne
            nodi = NODE_bound(kelem+(ne-1)*NNod+J)
            iblk=idnode(nodi)+1
            do k=1,nblk_ne
               if(iblk.eq.nb_elm(k)) goto 100
            enddo
            nblk_ne=nblk_ne+1
            nb_elm(nblk_ne)=iblk
100         continue
          enddo
          do k=1,nblk_ne
            iblk=nb_elm(k)
            mnodei(iblk,ityp)=mnodei(iblk,ityp)+1
            mlgelm_no(iblk,num_iblk(iblk)+mnodei(iblk,ityp))=ne
          end do
        end do
        do iblk=1,numblk
          num_iblk(iblk)=num_iblk(iblk)+mnodei(iblk,ityp)
          nnodei(iblk,ityp)=nnod
        enddo
        kelem=kelem+num*nnod
      enddo
c
C.....OPEN AND READ ID FILE
      open(32,file=filename(4),form='unformatted',status='old')
      READ (32) NSKP,NSKP,((ID(I,J),J=1,KDGOF),I=1,KNODE)
      close(32)
c      call getpart_3(kdgof,knode,ID)
c
c     get global equational number ..........
c
c
      do j=1,kdgof
        do i=1,knode
          nodvar(i,j)=id(i,j)
        enddo
      enddo
 
      neq=0
      do i=1,knode
        inod=nod_ord(i)
        do j=1,kdgof
          if(nodvar(inod,j).gt.0) then
            neq=neq+1
            nodvar(inod,j)=neq
          endif
        enddo
      enddo
 
      do i=1,knode
        inod=nod_ord(i)
        do j=1,kdgof
          if( nodvar(inod,j) .lt. -1) then
            N = -nodvar(inod,j)-1
            nodvar(inod,j) = nodvar(N,j)
          endif
        enddo
      enddo
C.......Start Data partition.....................................
c
c     get element ,equation number,id0 information and then send them
c
      do iblk=1,numblk
        npartition = iblk
        print *,'Begin ',iblk,'  subdomain data partition...'
c
c       get local to global node number array of iblk domain
c
        nod_iblk=knode_iblk(iblk)
        do i=1,nod_iblk
          inode(i)=mlgnod(iblk,i)
        enddo
c
c       get the number of all element and material at iblk subdomain firstly
c
        kelem_iblk=0
        kemate_iblk=0
        do ityp=1,numtyp_body
          kelem_iblk=kelem_iblk
     *        +nnodei(iblk,ityp)*mnodei(iblk,ityp)
          kemate_iblk=kemate_iblk
     *        +mmate(ityp)*nmate(ityp)
        enddo
c
c     get the local enode array of boundary elements at iblk subdomain
c
        ne1=0
        do ityp=1,numtyp_body
          ne1=ne1+mnodei(iblk,ityp)
        enddo
        nebd1=0
        kelem1=0
        do ityp=numtyp_body+1,numtyp
          num=mnodei(iblk,ityp)
          nnod=nnodei(iblk,ityp)
          nne=nnod
          if(material .eq. 'y' .or. material .eq. 'Y') nne=nne-1
          do ne=1,num
            ne_g=mlgelm_no(iblk,ne1+ne)
            do inod=1,nne
              nod_g=node_bound(kelem1+(ne_g-1)*nnod+inod)
              if((idnode(nod_g)+1) .eq.iblk) then
                jnod=mglnod(nod_g)
                goto 400
              endif
              do jnod=1,knode_iblk(iblk)
                if(nod_g.eq.mlgnod(iblk,jnod)) goto 400
              enddo
              print *,'Fatal Error for getting local node array!!!'
400           continue
              nodei_bound(nebd1+(ne-1)*nnod+inod)
     *            =jnod
            enddo
            if(material.eq.'y'.or.material.eq.'Y') then
               nodei_bound(nebd1+ne*nnod)
     *                    =node_bound(kelem1+ne_g*nnod)
            endif
          enddo
          ne1=ne1+num
          nebd1=nebd1+num*nnod
          kelem1=kelem1+mnode(ityp)*nnode(ityp)
          kelem_iblk=kelem_iblk
     *        +nnodei(iblk,ityp)*mnodei(iblk,ityp)
          kemate_iblk=kemate_iblk
     *        +mmate(ityp)*nmate(ityp)
        enddo
cc
cc       send some constants to iblk subdomain for malloc memory
cc
        call sendint(iblk,nsource,kelem_iblk)
        call sendint(iblk,nsource,kemate_iblk)
c
c       get local to global nv file for each subdomain
C       ipool  : the local equation number map to global
C               ( array maplg in sgetpart ,file lgnv in slave process)
C       iapool : the index of inside or external equation number, 1: inside ,0: external
C               (imaplg in sgetpart, file idxlgnv in slave process)
c
        neqi = 0
        n_update =0
        do i=1,nod_iblk
          nodg = inode(i)
          if(iorder(iblk,i).eq.1) then
            do j=1,kdgof
              if(id(nodg,j).gt.0) then
                neqi=neqi+1
                n_update=n_update+1
                ipool(neqi)=nodvar(nodg,j)
                iapool(neqi)=1
              elseif(id(nodg,j) .lt. -1) then
                N = -id(nodg,j)-1
                do k=1,nod_iblk
                  if( inode(k) .eq. N) exit
                enddo
                if(iblk .eq. 0) then
                  do k=1,knode_iblk(iblk)
                    if(mlgnod(iblk,k) .eq. 73211) then
                      print *,'11111111111111111111111'
                    endif
                  enddo
                endif
                if( inode(k) .ne. N) then
                  write(*,'(2a)') 'Error!! Slave point in block ',
     +                            'diffrent with Master point!!!'
                  write(*,'(a,12i7)') 'iblk, K, nod_iblk=',iblk,k,nod_iblk
                  write(*,'(a,12i7)') 'i,N,nodg=',i, N,nodg
                  write(*,'(a,12i7)')'idnode(N),idnode(nodg)=',idnode(N), idnode(nodg)
                  stop 123
                endif
              endif
            end do
          else
            do j=1,kdgof
              if(id(nodg,j).gt.0) then
                neqi=neqi+1
                ipool(neqi)=nodvar(nodg,j)
                iapool(neqi)=0
              elseif(id(nodg,j) .lt. -1) then
                N = -id(nodg,j)-1
                do k=1,nod_iblk
                  if( inode(k) .eq. N) exit
                enddo
                if( inode(k) .ne. N) then
                  neqi=neqi+1
                  ipool(neqi)=nodvar(nodg,j)
                  iapool(neqi)=0
                endif
              endif
            end do
          endif
        enddo
        print *,'getpart: iblk, n_update =', iblk, n_update
c
c       send lgnv ,idxlgnv and nupdate to iblk domain
c
        call sendint(iblk,nsource,neqi)
        call sendai(iblk,nsource,ipool,neqi)
        call sendai(iblk,nsource,iapool,neqi)
        call sendint(iblk,nsource,n_update)
c
c       get id0 and send it in iblk subdomain
c
        call getpart_1(iblk,nsource,kdgof,knode,nod_iblk,
     *                 inode,id,ipool)
c        WRITE(*,*) 'KNODE,KDGOF,IBLK = ',nod_iblk,kdgof,iblk
c
c       deal with element and mate information for each subdomain and multiplier
c
c       get the elements and materials of every type at iblk subdomain
c       and then send them
c
        kelem=0
        kemate=1
        do ityp=1,numtyp_body
cc         get the body element of ityp type at iblk domain
          mi_ni=mnodei(iblk,ityp)*nnodei(iblk,ityp)
c          if(mi_ni .gt. maxip) then
c            print *,"no enough space of elements ipool",ityp,iblk
c            stop
c          endif
          nnod=nnodei(iblk,ityp)
          do i=1,mnodei(iblk,ityp)
            do j=1,nnod
              ipool((I-1)*nnod+J)=node_iblk(iblk,(I-1)*nnod+J+kelem)
            enddo
          enddo
c         send the element
          call sendint(iblk,nsource,mnodei(iblk,ityp))
          call sendint(iblk,nsource,nnod)
          call sendai(iblk,nsource,ipool,mi_ni)
          kelem=kelem+mi_ni
c
cc         get the material of itype type at iblk domain
c
          if(material .eq. 'y' .or. material .eq. 'Y') then
            mi_ni=mmate(ityp)*nmate(ityp)
c         send the matireal
            call sendint(iblk,nsource,mmate(ityp))
            call sendint(iblk,nsource,nmate(ityp))
            call sendar(iblk,nsource,emate_body(kemate),mi_ni)
            kemate=kemate+mi_ni
          endif
        enddo
c
        kelem=1
        kemate=1
        do ityp=numtyp_body+1,numtyp
cc         get the boundary element of ityp type at iblk domain
          mi_ni=mnodei(iblk,ityp)*nnodei(iblk,ityp)
          nnod=nnodei(iblk,ityp)
          nne=nnod
          if(material.eq.'y' .or. material.eq.'Y') nne=nnod-1
c         send the element
          call sendint(iblk,nsource,mnodei(iblk,ityp))
          call sendint(iblk,nsource,nnod)
          call sendai(iblk,nsource,nodei_bound(kelem),mi_ni)
          kelem=kelem+mi_ni
c
cc         get the material of itype type at iblk domain
c
          if(material .eq. 'y' .or. material .eq. 'Y') then
            mi_ni=mmate(ityp)*nmate(ityp)
c         send the material
            call sendint(iblk,nsource,mmate(ityp))
            call sendint(iblk,nsource,nmate(ityp))
            call sendar(iblk,nsource,emate_bound(kemate),mi_ni)
            kemate=kemate+mi_ni
          endif
        enddo
C       end do of ityp
 
      enddo    !enddo of iblk=1,numblk
c
c       deal with disp0 disp1 disp2 disp3 file for each subdomain cc
c
C.....OPEN AND READ DISP0 FILE
      rewind(33)
      READ (33) NUMNOD,NODDOF,((U012(I,J),J=1,NODDOF),I=1,NUMNOD)
      close(33)
c
c       get disp0 data and send them in iblk subdomain
c
      do iblk=1,numblk
        call getpart_2(iblk,nsource,kdgof,knode,knode_iblk(iblk),numblk,
     *                      max_lknode,mlgnod,u012,apool)
      enddo    !enddo of iblk=1,numblk
c
C.....OPEN AND WRITE UNOD FILE
c
      open(33,file=filename(5),form='unformatted',status='unknown')
      write(33) ((u012(i,j),i=1,numnod),j=1,noddof)
      close(33)
c
C.....OPEN AND READ DISP1 FILE
c
      inquire(file=filename(6),exist=filflg)
      if(filflg) then
        open(34,file=filename(6),form='unformatted',status='old')
        READ (34) NUMNOD,NODDOF,((U012(I,J),J=1,NODDOF),I=1,NUMNOD)
        close(34)
c
c       get disp1 data and send them in iblk subdomain
c
        do iblk=1,numblk
          call getpart_2(iblk,nsource,kdgof,knode,knode_iblk(iblk),
     &                   numblk,max_lknode,mlgnod,u012,apool)
        enddo    !enddo of iblk=1,numblk
      endif
c
C.....OPEN AND READ DISP2 FILE
c
      inquire(file=filename(7),exist=filflgdisp(1))
      if (filflgdisp(1)) then
        ndisp=ndisp+1
        open(35,file=filename(7),form='unformatted',status='old')
        READ (35) NUMNOD,NODDOF,((U012(I,J),J=1,NODDOF),I=1,NUMNOD)
        close(35)
c
c       get disp2 data and send them in iblk subdomain
c
        do iblk=1,numblk
          call getpart_2(iblk,nsource,kdgof,knode,knode_iblk(iblk),
     &                   numblk,max_lknode,mlgnod,u012,apool)
        enddo    !enddo of iblk=1,numblk
 
      endif
c
C.....OPEN AND READ DISP3 FILE
c
      inquire(file=filename(8),exist=filflgdisp(2))
      if (filflgdisp(2)) then
        open(36,file=filename(8),form='unformatted',status='old')
        READ (36) NUMNOD,NODDOF,((U012(I,J),J=1,NODDOF),I=1,NUMNOD)
        close(36)
c
c       get disp3 data and send them in iblk subdomain
c
        do iblk=1,numblk
        call getpart_2(iblk,nsource,kdgof,knode,knode_iblk(iblk),
     &                 numblk,max_lknode,mlgnod,u012,apool)
        enddo    !enddo of iblk=1,numblk
      endif
 
1000  format(i10,3e15.6)
1100  format(10i10)
1200  format(i10,e15.5)
1201  format(i10,10e15.5)
 
      return
      end
c
      subroutine getpart_1(iblk,nsource,kdgof,knode,nod_iblk,
     *                      inode,id,ipool)
      implicit real*8 (a-h,o-z)
      dimension ipool(1),id(knode,kdgof),inode(1)
c
        do i=1,nod_iblk
          nodg=inode(i)
          do j=1,kdgof
            N = id(nodg,j)
            if( N .lt. -1) then
              jgnod = -N-1
              do k=1,nod_iblk
                if(inode(k) .eq. jgnod) exit
              enddo
              if( inode(k) .eq. jgnod) then
                ipool(i+(j-1)*nod_iblk)=-k-1
              else
                ipool(i+(j-1)*nod_iblk)=1
              endif
            else
              ipool(i+(j-1)*nod_iblk) = N
            endif
          end do
        end do
cc
cc       send id to iblk domain
cc
c        WRITE(*,*) 'KNODE,KDGOF,IBLK = ',nod_iblk,kdgof,iblk
        call sendai(iblk,nsource,ipool,nod_iblk*kdgof)
c
      return
      end
c
c
      subroutine getpart_2(iblk,nsource,kdgof,knode,nod_iblk,numblk,
     *                      max_lknode,mlgnod,u,apool)
      implicit real*8 (a-h,o-z)
      dimension apool(1),u(knode,kdgof)
      dimension mlgnod(numblk,max_lknode)
c
        do j=1,kdgof
          do i=1,nod_iblk
            apool(i+(j-1)*nod_iblk)=u(mlgnod(iblk,i),j)
          end do
        end do
cc
cc       send id to iblk domain
cc
c        WRITE(*,*) 'KNODE,KDGOF,IBLK = ',nod_iblk,kdgof,iblk
        call sendar(iblk,nsource,apool,nod_iblk*kdgof)
c
      return
      end
c
      subroutine getpart_3(kdgof,knode,id)
      implicit real*8 (a-h,o-z)
      dimension id(knode,kdgof),ms(1000),is(1000)
      do 1001 k=1,kdgof
      m = 0
      do 801 n=1,knode
      if (id(n,k).le.1) goto 801
       j=id(n,k)
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
        id(n,k)=1
       else
        id(n,k)=-j0-1
       endif
801   continue
1001  continue
      return
      end
 
