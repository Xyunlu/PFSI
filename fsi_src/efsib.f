      subroutine efsib(iorder,iorder1,iorder2,
     & coor,eud,evd,ewd,
     & eus,evs,ews,node,
     & emate,coor0,ur,ubc,
     & nodvar,u,euf,eu,
     & emass,f,sml,apool,
     & iblk,nsource,knode,kdgof,knodeg,kdgofg,knode2,kdgof2,
     & neq,kvar,kvarg,kelem,kemate,kcoor,time,dt,
     & it,itn,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      dimension nodvar(knode,kdgof),u(knode,kdgof),coor(knode,kcoor),
     + euf(knode2,kdgof2),coor0(knodeg,kcoor),
     + ur(knodeg,kdgof),iorder(knodeg),iorder1(knode),iorder2(knode2),
     + ubc(knode2,kdgof2),
     * eu(knode,kdgof),eud(knode),evd(knode),ewd(knode),
     * eus(knode),evs(knode),ews(knode),Emass(kvar),
     & F(KVAR),EMATE(kemate),NODE(KELEM),SML(100000)
      dimension apool(kvarg)
      double precision,dimension(:),allocatable :: force, rotiner
 
      dimension mnode(100),nnode(100),mmate(100),nmate(100)
 
      common /rdata/ xo,yo,zo
      common /idata/ id4updown
C.....
      open (41,file='init.dat',form='formatted',status='old')
      read (41,*) xmax,ymax,zmax,xo,yo,zo
      read (41,*) skip
      read (41,*) skip
      read (41,*) skip
      read (41,*) skip
      read (41,*) skip
      read (41,*) skip
      read (41,*) id4updown
      close(41)
C .................................................................
6     FORMAT (1X,26I3)
7     FORMAT (1X,6E12.3)
1001  FORMAT(1X,9I7)
 
 
      do i=1,knode
        do j=1,kdgof
          nodvar(i,j) = 1
        enddo
      enddo
 
C.......COMPUTE NODVAR
      NEQ = 0
      DO 100 I=1,KNODE
      DO 100 J=1,KDGOF
      IF (NODVAR(I,J).GE.-1) then
      NEQ = NEQ+1
      NODVAR(I,J) = NEQ
      else
      N = -NODVAR(I,J)-1
      NODVAR(I,J) = NODVAR(N,J)
      endif
      U(I,J) = 0.0
100   CONTINUE
C      WRITE(*,*) 'KDGOF =',KDGOF,' KNODE =',KNODE
C      WRITE (*,*) 'NODVAR ='
C      WRITE (*,6) ((NODVAR(I,J),J=1,KDGOF),I=1,KNODE)
 
C.......OPEN COOR file
 
      numtyp = 1
C      IF (IT.EQ.0) THEN
C#OPENFILE.sub unknown
C      ELSE
C#OPENFILE.sub OLD
C      ENDIF
 
      DO 110 I=1,NEQ
      Emass(i) = 0.0
110   CONTINUE
C.......OPEN ELEM0 file
C.......NUMTYP number of element types
      NUMEL=0
      JNODE=1
      JMATE=1
      DO 2000 ITYP=1,NUMTYP
C.......INPUT ENODE
      if (MNODE(ITYP).EQ.0) GOTO 2001
      NUM = MNODE(ITYP)
      NNE = NNODE(ITYP)
      nne = nne-1
      K=0
      DO 115 I=1,NNE
      INOD = NODE(I)
       do idof=1,kdgof
       IF (nodvar(inod,idof).NE.0) K=K+1
       enddo
115     CONTINUE
c      WRITE(*,*) 'K =',K
      kk=k*k
      k0=1
      k1=k0+k*k
      k2=k1+k
      k3=k2+k
      k4=k3+k*k
      k5=k4+k*k
      k6=k5+k
      k7=k6+k
c        WRITE(*,*) 'MMATE =',MMATE,' NMATE =',NMATE
c        WRITE (*,*) 'EMATE ='
c        WRITE (*,18) ((EMATE((I-1)*NMATE+J),J=1,NMATE),
c     *  I=1,MMATE)
 
c      ec = 0.D0
c      sml(k2) = 0.D0
      allocate(force(num), rotiner(num))
 
      CALL fsibSUB(ITYP,IT,TIME,DT,K,KK,NUM,NNE,
     & KNODE,KDGOF,NODVAR,
     & KCOOR,COOR,
     & MNODE(ITYP),NNODE(ITYP),NODE(JNODE),
     & MMATE(ITYP),NMATE(ITYP),EMATE(JMATE),
     *sml(k0),sml(k1),sml(k2),sml(k3),sml(k4),sml(k5),sml(k6),
     &eud,evd,ewd,eus,evs,ews,eu,
     &Emass,
     & U,rotiner,force)
      JNODE=JNODE+MNODE(ITYP)*NNODE(ITYP)
2001  JMATE=JMATE+MMATE(ITYP)*NMATE(ITYP)
 
c      rotiner = ec
c      force = sml(k2)
      call sendar(nsource,iblk,rotiner,num)
      call sendar(nsource,iblk,force,num)
 
2000    CONTINUE
 
cc#CLOSFILE.sub
 
C....
C....  For inertia compute
C....
 
C...  Send coor0 to master processor
      ntemp = 0
      do inod=1,knodeg
        if(iorder(inod) .eq. 1) then
          do ic=1,kcoor
            ntemp = ntemp+1
            apool(ntemp) = coor0(inod,ic)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
C...  Send ur to master processor
      ntemp = 0
      do inod=1,knodeg
        if(iorder(inod) .eq. 1) then
          do idgof=1,kdgof
            ntemp = ntemp+1
            apool(ntemp) = ur(inod,idgof)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
C.... Read subdomain file and Compute ur at master processor
C...  Send ubc to master processor
      ntemp = 0
      do inod=1,knode2
        if(iorder2(inod) .eq. 1) then
          do idgof=1,kdgof2
            ntemp = ntemp+1
            apool(ntemp) = ubc(inod,idgof)
          enddo
        endif
      enddo
      call sendar(nsource,iblk,apool,ntemp)
 
C... Update ubc
      call recvar(iblk,nsource,apool,knode2*kdgof2)
      ntemp = 0
      do inod=1,knode2
        do jgof=1,kdgof2
          ntemp = ntemp + 1
          ubc(inod,jgof) = apool(ntemp)
        enddo
      enddo
 
C... Update ur
      call recvar(iblk,nsource,apool,knodeg*kdgof)
      ntemp = 0
      do inod=1,knodeg
        do jgof=1,kdgof
          ntemp = ntemp + 1
          ur(inod,jgof) = apool(ntemp)
        enddo
      enddo
 
      END
 
      SUBROUTINE fsibSUB(ITYP,IT,TIME,DT,K,KK,NUM,NNE,
     & KNODE,KDGOF,NODVAR,
     & NCOOR,COOR,
     & MNODE,NNODE,NODE,MMATE,NMATE,EMATE,
     *es,em,ef,Estifn,Estifv,Emassn,Emassv,
     *eud,evd,ewd,eus,evs,ews,eu,Emass,
     & U,rotiner,force)
      implicit real*8 (a-h,o-z)
      DIMENSION NODVAR(KNODE,KDGOF),COOR(KNODE,NCOOR),NODE(*),
     & U(KNODE,KDGOF),EMATE(*),rotiner(num),force(num),
     *es(k,k),em(k),ef(k),eud(knode),
     *evd(knode),ewd(knode),eus(knode),evs(knode),
     *ews(knode),eu(knode,kdgof),Estifn(k,k),Estifv(kk),
     *Emassn(k),Emassv(k),Emass(1),
     & R(500),PRMT(500),COEF(500),LM(500)
17      FORMAT (1X,15I5)
18      FORMAT (1X,8e9.2)
 
      DO 1000 NE=1,NUM
        ec = 0.D0
        ef(1) = 0.D0
        force(ne) = 0.D0
        rotiner(ne) = 0.D0
 
        NR=0
        DO 130 J=1,NNE
          JNOD = NODE((NE-1)*NNODE+J)
          IF (JNOD.LT.0) JNOD = -JNOD
          coef(j+0*nne)=eud(jnod)
          coef(j+1*nne)=evd(jnod)
          coef(j+2*nne)=ewd(jnod)
          coef(j+3*nne)=eus(jnod)
          coef(j+4*nne)=evs(jnod)
          coef(j+5*nne)=ews(jnod)
          DO 120 I=1,NCOOR
           NR=NR+1
120       R(NR) = COOR(JNOD,I)
130     CONTINUE
        IMATE = NODE(NNODE*NE)
        if (imate.eq.4) then
           l=0
           m=0
           goto 999
        endif
        DO 140 J=1,NMATE
140      PRMT(J) = EMATE((IMATE-1)*NMATE+J)
        PRMT(NMATE+1)=TIME
        PRMT(NMATE+2)=DT
        PRMT(NMATE+3)=IMATE
        prmt(NMATE+4)=NE
        prmt(NMATE+5)=NUM
        prmt(nmate+6)=it
        prmt(nmate+7)=nmate
 
        goto 1
1       call inertia(r,coef,prmt,es,em,ec,ef,ne)
        goto 2
2       continue
 
999     CONTINUE
 
        rotiner(ne) = ec
        force(ne) = ef(1)
 
1000  CONTINUE
 
      RETURN
      END
 
