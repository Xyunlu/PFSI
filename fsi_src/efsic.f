      subroutine efsic(nodvar,coor,bfu,
     & jdiag,jdiagaz,na,node,
     & emate,a,f,u,
     & eu,eu1,sml,u1,
     & lm,uu,iblk,kdgof,knode,maxa,
     & neq,neq1,kvar,kelem,kemate,kcoor,time,dt,
     & it,itn,cc,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      dimension nodvar(knode,kdgof),u(knode,kdgof),coor(knode,kcoor),
     *eu(knode,kdgof),eu1(knode,kdgof),
     &jdiag(kvar),jdiagaz(kvar),na(maxa),A(MAXA),F(kvar),EMATE(kemate),
     &NODE(KELEM),SML(100000),u1(kvar),lm(1000),
     *uu(kvar),bfu(knode,kdgof)
 
      DIMENSION MNODE(100),NNODE(100),MMATE(100),NMATE(100)
      CHARACTER*1 MATERIAL
      open(88,file='control.dat',status='old')
      read(88,*) emax,numstep,emax,numstep
      close(88)
 
c@sbegin
6     format (1x,10i5)
7     format (1x,6e12.5)
 
      ksml = 100000
C.... open time file
c ....read nv file
c       write (*,*) 'nodvar ========================'
c       write (*,6) ((nodvar(i,j),i=1,kdgof),j=1,knode)
c ....read coor file
c       write(*,*) 'coor =========================='
c       write(*,7) ((coor(j,i),j=1,ncoor),i=1,knode)
c ....read bfd file
c      write (*,*) 'bf ==========================='
c      write(*,7) ((u(i,j),j=1,kdgof),i=1,knode)
      do i=1,knode
         do j=1,kdgof
            u(i,j)=bfu(i,j)
         enddo
      enddo
      numtyp = 1
 
c.....  read diag file
c        read(29) (jdiag(i),i=1,neq)
c        read(29) (jdiagaz(i),i=1,neq1)
c       write (*,*) 'jdiag jdiagaz ================'
c        write(*,*) (jdiag(i),i=1,neq)
c        write(*,*) (jdiagaz(i),i=1,neq1)
c.....read na file
c        read(31) (na(i),i=1,maxa)
c       write (*,*) 'na============================'
c        do i = 1,neq
c        nl = jdiagaz(i)
c        nu = jdiagaz(i+1)-1
c        write(*,*) (na(j),j=nl,nu)
c        end do
C.......OPEN ELEM0 file
 
      DO 111 I=1,MAXA
      A(I) = 0.0
111   CONTINUE
 
      DO 2300 I=1,NEQ
2300    CONTINUE
      NUMEL=0
      JNODE=1
      JMATE=1
C.....OPEN EMATE+ENODE+ELOAD file
      do 2000 ityp=1,numtyp
C.......INPUT ENODE
        if (MNODE(ITYP).EQ.0) GOTO 2001
        num=mnode(ityp)
        nne = nnode(ityp)
        knum=num*nne
        if (knum.gt.kelem) then
        write(*,*) 'it is error: knum gt kelem'
        stop 0000
        endif
c       write(*,*) 'num =',num,' nnode =',nne
c       write(*,*) 'node ='
c       write(*,6) ((node(jnode+(i-1)*nne+j-1),j=1,nne),i=1,num)
      nne = nne-1
        k=0
        do 115 I=1,nne
          INOD = node(I)
       do idof=1,kdgof
          if (nodvar(inod,idof).ne.0) k=k+1
       enddo
115     continue
c        write(*,*) 'k =',k
      kk=k*k
      k0=1
      k1=k0+k*k
      k2=k1+k
      k3=k2+k
      k4=k3+k*k
      k5=k4+k*k
        if (k5.ge.ksml) then
        write(*,*) 'it is error: k5 gt ksml'
        stop 0000
        endif
c        WRITE(*,*) 'MMATE =',MMATE,' NMATE =',NMATE
c        WRITE (*,*) 'EMATE ='
c        WRITE (*,18) ((EMATE(jmate+(I-1)*NMATE(ityp)+J-1),
c     *  J=1,NMATE(ityp)), I=1,MMATE(ityp))
      CALL fsicSUB(iblk,ITYP,IT,TIME,DT,K,KK,NUM,NNE,
     & KNODE,KDGOF,NODVAR,
     & KCOOR,COOR,
     & MNODE(ITYP),NNODE(ITYP),NODE(JNODE),
     & MMATE(ITYP),NMATE(ITYP),EMATE(JMATE),
     & neq,neq1,MAXA,na,A,
     & jdiag,jdiagaz,
     *sml(k0),sml(k1),sml(k2),sml(k3),sml(k4),
     &eu,eu1,
     & U)
      JNODE=JNODE+MNODE(ITYP)*NNODE(ITYP)
2001  JMATE=JMATE+MMATE(ITYP)*NMATE(ITYP)
2000  continue
 
      DO 2050 IJ=1,NEQ
      F(IJ)=0.0D0
2050    CONTINUE
      do inod=1,knode
       do idof=1,kdgof
         ij=nodvar(inod,idof)
         if (ij.gt.0) then
           if (ij.gt.neq) neq = ij
           f(ij)=f(ij)+u(inod,idof)
         endif
       enddo
      enddo
c       write(*,*) 'a(ii) ======='
c       write(*,7) (a(jdiag(ii)),ii=1,neq)
c       write(*,*) 'b(ii) ======='
c       write(*,7) (b(jdiag(ii)),ii=1,neq)
c       write(*,*) '******** call redu ***********'
c        call reduin(a,f,jdiag,neq,maxa,3)
cc        write(*,*) 'a =========='
cc        write(*,17) (a(jdiag(ii)),ii=1,neq)
cc        write(*,*) 'f =========='
cc        write(*,17) (f(ii),ii=1,neq)
        return
        end
 
      SUBROUTINE fsicSUB(iblk,ITYP,IT,TIME,DT,K,KK,NUM,NNE,
     & KNODE,KDGOF,NODVAR,
     & NCOOR,COOR,
     & MNODE,NNODE,NODE,MMATE,NMATE,EMATE,
     & neq,neq1,MAXA,na,A,
     & jdiag,jdiagaz,
     *es,em,ef,Estifn,Estifv,eu,eu1,
     & U)
      implicit real*8 (a-h,o-z)
      DIMENSION NODVAR(KNODE,KDGOF),COOR(KNODE,NCOOR),NODE(*),
     & U(KNODE,KDGOF),EMATE(*),
     & A(MAXa),na(MAXa),
     *es(k,k),em(k),ef(k),eu(knode,kdgof),
     *eu1(knode,kdgof),Estifn(k,k),Estifv(kk),
     & jdiag(neq),jdiagaz(neq1),
     & R(500),PRMT(500),COEF(500),LM(500)
      common /area/ ar
      integer, dimension(:),allocatable :: ielmlg
      double precision,dimension(:),allocatable :: area
17      FORMAT (1X,15I5)
18      FORMAT (1X,8e9.2)
 
      open(22,file='mlgelm_no_2',form='formatted',status='old')
      read(22,*) numblk, numtyp
      do ib=1, iblk-1
        read(22,*) numi
        read(22,*) (ntmp, i=1,numi)
      enddo
      read(22,*) numi
      if( numi .ne. num) then
        write(*,*) 'Error! Read NUM from mlgelm_no Error!!, numi, num=',numi,num
        call endjob(ierr)
      endif
      allocate(ielmlg(numi))
      read(22,*) (ielmlg(i),i=1,numi)
      close(22)
      open (41,file='fcellarea',form='unformatted',status='old')
      read (41) numg
      allocate(area(numg))
      rewind(41)
      read (41) nnn,(area(i),i=1,nnn)
      close(41)
 
      DO 1000 NE=1,NUM
      NR=0
      DO 130 J=1,NNE
      JNOD = NODE((NE-1)*NNODE+J)
        IF (JNOD.LT.0) JNOD = -JNOD
      DO 120 I=1,NCOOR
        NR=NR+1
120     R(NR) = COOR(JNOD,I)
130      CONTINUE
      IMATE = NODE(NNODE*NE)
      if (imate.eq.4) then
        l=0
        m=0
        goto 1000
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
 
      neg = ielmlg(ne)
      ar = area(neg)
 
      goto 1
1     call laplace(r,coef,prmt,es,em,ec,ef,ne)
      goto 2
2     continue
cc#rwsml.sub
 
c       write(*,*) 'es ef ='
c       do 555 i=1,k
c555    write(*,18) (es(i,j),j=1,k)
c       write(*,18) (em(i),i=1,k)
c       write(*,18) (ef(i),i=1,k)
 
cc    if (it.gt.0) then
      do 201 i=1,k
      do 201 j=1,k
      Estifn(i,j)=0.0
201   continue
      do 202 i=1,k
      Estifn(i,i)=Estifn(i,i)
      do 202 j=1,k
      Estifn(i,j)=Estifn(i,j)+es(i,j)
202   continue
 
      l=0
      m=0
      i=0
      do 700 inod=1,nne
      nodi=node((ne-1)*nnode+inod)
      do 600 idgf=1,kdgof
      inv=nodvar(nodi,idgf)
      if (inv.eq.0) goto 600
      i=i+1
      if (inv.lt.0) goto 305
      l=l+1
      lm(l)=inv
      u(nodi,idgf)=u(nodi,idgf)
     *+ef(i)
305   j=0
      do 500 jnod=1,nne
      nodj=node((ne-1)*nnode+jnod)
      do 400 jdgf=1,kdgof
      jnv=nodvar(nodj,jdgf)
      if (jnv.eq.0) goto 400
      j=j+1
 
      if (jnv.lt.0) goto 400
      if (inv.lt.0) goto 310
      m=m+1
      Estifv(m)=Estifn(i,j)
310   continue
 
 
        if (inv.lt.0)
     *  U(NODJ,JDGF)=U(NODJ,JDGF)-ESTIFN(I,J)*U(NODI,IDGF)
400   continue
500   continue
600   continue
700   continue
      call nzaddmbs(a,na,jdiag,jdiagaz,l,lm,estifv,neq1,maxa)
      do i=1,l
        j=lm(i)
      enddo
1000  continue
      deallocate(ielmlg, area)
 
      return
      end
 
 
