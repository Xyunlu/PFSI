      subroutine starta(id,nodvar,coor,
     & dsp0,bfu,u0,u1,
     & u2,u3,um,vm,
     & wm,umn,vmn,wmn,
     & evlx,evly,evlz,evlxn,
     & evlyn,evlzn,ur,ua,
     & node,emate,jdiag,jdiagaz,
     & na,r,lm,jdiagr,
     & nsource,iblk,knode,kdgof,maxa,maxb,neq,neq1,
     & kvar,kelem,kemate,kcoor,numetyp,numelem,time,material,
     & itn,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      DIMENSION  NODVAR(KNODE,KDGOF),coor(knode,kcoor),
     &  jdiag(kvar),jdiagaz(kvar),na(maxa),id(knode,kdgof),
     &  U0(KNODE,KDGOF),u1(KNODE,KDGOF),u2(knode,kdgof),
     &  u3(knode,kdgof),BFU(KNODE,KDGOF),dsp0(knode,kdgof),
     &  node(kelem),emate(kemate),r(3),
     *  lm(1000),jdiagr(kvar)
      Dimension um(knode),vm(knode),wm(knode),
     &          umn(knode),vmn(knode),wmn(knode),
     &          evlx(knode),evly(knode),evlz(knode),
     &          evlxn(knode),evlyn(knode),evlzn(knode)
      dimension ur(knode,kcoor), ua(knode,kcoor)
 
      dimension mnode(100),nnode(100),mmate(100),nmate(100)
      CHARACTER*1 MATERIAL
 
6       FORMAT (1X,10I5)
7       FORMAT (1X,6E12.5)
 
      itn = 1
C...... READ ID0 FILE
      NEQ=0
      do inod=1,knode
       do idof=1,kdgof
        IF (id(inod,idof).LT.1) then
          nodvar(inod,idof)=id(inod,idof)
        else
          NEQ = NEQ + 1
          nodvar(inod,idof) = NEQ
        endif
       enddo
      enddo
      neq1=neq+1
      do inod=1,knode
       do idof=1,kdgof
         IF (nodvar(inod,idof).LT.-1) then
           N = -nodvar(inod,idof)-1
27         CONTINUE
           IF (nodvar(inod,idof).LT.-1) THEN
             N=-nodvar(inod,idof)-1
             GOTO 27
           ENDIF
           nodvar(inod,idof) = NODVAR(n,idof)
         endif
       enddo
      enddo
C.......WRITE NV FILE
C.......open coor0 file
 
C.... read disp0 file
 
      do inod=1,knode
       do idof=1,kdgof
        bfu(inod,idof)=dsp0(inod,idof)
       enddo
      enddo
 
C....  write bfd file
 
      zo=0.d0
c      do 321 n=1,knode
c       do 100 j=1,kcoor
c100    r(j)=coor(n,j)
c       do 200 j=1,kdgof
c        u0(n,j)=bound(r,zo,j)
c        u1(n,j)=bound1(r,zo,j)
c        u2(n,j)=bound2(r,zo,j)
c200    continue
c321   continue
 
      zero = 0.D0
      do i=1,knode
        do j=1,kdgof
          u1(i,j) = u0(i,j)
          u2(i,j) = zero
        enddo
      enddo
 
      do i=1,knode
        um(i) = 0.D0
        vm(i) = 0.D0
        wm(i) = 0.D0
        umn(i) = 0.D0
        vmn(i) = 0.D0
        wmn(i) = 0.D0
        evlx(i) = 0.D0
        evly(i) = 0.D0
        evlz(i) = 0.D0
        evlxn(i) = 0.D0
        evlyn(i) = 0.D0
        evlzn(i) = 0.D0
        do j=1,kcoor
          ur(i,j) = 0.D0
          ua(i,j) = 0.D0
        enddo
      enddo
 
C.... read and update initial value
 
C.... write unodd
 
C.... write unodr
 
C.... write unodac
 
      DO I=1,NEQ
         JDIAG(I)=0
         JDIAGR(I)=0
         JDIAGAZ(I)=0
      ENDDO
      jdiagaz(neq1)=0
      DO I=1,MAXA
         NA(I)=0
      ENDDO
c
      NUMELEM=0
      jnode=1
      jmate=1
      NDMAX=0
      nlast = 0
      nnai = 0
C...  Open elem0 file
      DO 2000 ITYP=1,NUMETYP
C..... READ ELEM0 FILE
      IF(mnode(ITYP).EQ.0) GOTO 2000
c      WRITE(*,*) 'NUM =',Mnode(ityp),' NNODE =',NNODE(ityp)
c      WRITE(*,*) 'NODE ='
c      WRITE(*,6) ((NODE(jnode-1+(I-1)*NNODE(ityp)+J),
c     * J=1,NNODE(ityp)),I=1,Mnode(ityp))
      num = mnode(ityp)
      nne = nnode(ityp)
      IF (MATERIAL.EQ.'Y' .OR. MATERIAL.EQ.'y') THEN
         NNE = NNE-1
      ENDIF
 
      DO 3000 NE=1,NUM
        L=0
        DO 3600 INOD=1,NNE
           NODI=NODE(jnode-1+(NE-1)*NNODE(ityp)+INOD)
           DO 3700 IDGF=1,KDGOF
              INV=NODVAR(NODI,IDGF)
              IF (INV.LE.0) GOTO 3700
              L=L+1
              LM(L)=INV
3700       CONTINUE
3600    CONTINUE
        NUMELEM=NUMELEM+1
c        WRITE (*,*) 'L,LM =',L
c        WRITE (*,'(1X,15I5)') (LM(I),I=1,L)
        NDMAX=MAX(NDMAX,L)
c----------------------------------------------------------
c
c       setup a link table which get the connection message of all unknowns
c
        do i = 1,L
          nunsi = lm(i)
c
          if(jdiagr(nunsi).eq.0) then
            jdiag(nunsi) = nlast + 1
            do j = 1,L
              nunsj = lm(j)
              nlast = nlast + 1
              if(nlast.ge.maxa) print *,'no enough memory for NA'
              na(nlast) = nunsj
              nlast = nlast + 1
              if(nlast.ge.maxa) print *,'no enough memory for NA'
              na(nlast) = nlast + 1
            end do
            jdiagr(nunsi) = l
          end if
c
          if(jdiagr(nunsi).gt.0) then
            nbtemp = jdiag(nunsi)
            do j = 1,L
              nunsj = lm(j)
              noyes = 0
              nbtemp = jdiag(nunsi)
              do k = 1,jdiagr(nunsi)
                kun = na(nbtemp)
                nbtemp0 = nbtemp
                nbtemp = na(nbtemp+1)
                if(kun.eq.nunsj) noyes = 1
              end do
c
              if(noyes.eq.0) then
                jdiagr(nunsi) = jdiagr(nunsi) + 1
                nlast = nlast + 1
                if(nlast.ge.maxa) print *,'no enough memory for NA'
                na(nbtemp0+1) = nlast
                na(nlast) = nunsj
                nlast = nlast + 1
                if(nlast.ge.maxa) print *,'no enough memory for NA'
                na(nlast) = nlast + 1
              end if
            end do
          end if
        end do
c
c----------------------------------------------------------
3000  CONTINUE
      jnode=jnode+mnode(ityp)*nnode(ityp)
      jmate=jmate+mmate(ityp)*nmate(ityp)
2000  CONTINUE
c------------------------------------------------------------
c
c       get the who Na memory structure and add again!
c
      jdiagaz(1) = 1
      do i = 1,neq
         jdiagaz(i+1) = jdiagaz(i) + jdiagr(i)
      end do
      maxa = jdiagaz(neq+1) -1
      maxb = maxa + 10
c      print *,'neq,neq1,maxa == ',neq,neq1,maxa
c
c      print *,(jdiag(i),i=1,neq)
c      print *,(jdiagr(i),i=1,neq)
c      print *,(jdiagaz(i),i=1,neq1)
c      do i = 1,neq
c         np = jdiag(i)
c         do j = 1,jdiagr(i)
c            lm(j) = na(np)
c            np = na(np+1)
c         end do
c         print *,(lm(j),j=1,jdiagr(i))
c      end do
c      print *,(na(i),i=1,maxa*2+2)
c
c      Then get the who matrix information again
c      Get Na dimension contains
c
      do i = 1,maxa
         na(i)=0
      end do
c
      do i=1,neq
         jdiag(i)=0
      end do
c
c      do i=1,neq1
c         jdiagr(i)=0
c      end do
c
      NUMELEM=0
      NDMAX=0
      jnode=1
      DO 5000 ITYP=1,NUMETYP
         NUM=Mnode(ITYP)
         NNe=Nnode(ITYP)
         if(num.eq.0) goto 5000
         IF (MATERIAL.EQ.'Y' .OR. MATERIAL.EQ.'y') THEN
            NNE = NNE-1
         ENDIF
cc        WRITE(*,*) 'NODE ='
cc        W
 
cc     * J=1,NNODE(ityp)),I=1,NUM)
         DO 6000 NE=1,NUM
            L=0
            DO 6600 INOD=1,NNE
               NODI=NODE(jnode-1+(NE-1)*NNODE(ityp)+INOD)
               DO 6700 IDGF=1,KDGOF
                  INV=NODVAR(NODI,idgf)
                  IF (INV.LE.0) GOTO 6700
                  L=L+1
                  LM(L)=INV
6700           CONTINUE
6600        CONTINUE
            NUMELEM=NUMELEM+1
C            WRITE (*,*) 'L,LM =',L
C            WRITE (*,'(1X,15I5)') (LM(I),I=1,L)
c
            call nzdmbsadd(neq,neq1,maxa,na,jdiag,jdiagr,
     *             jdiagaz,L,lm)
c
6000     CONTINUE
         jnode=jnode+mnode(ityp)*nnode(ityp)
         jmate=jmate+mmate(ityp)*nmate(ityp)
5000  CONTINUE
c
c       order Na dimension
c
      call nzdmbsorder(neq,neq1,maxa,na,jdiag,jdiagr,jdiagaz,L,lm)
c
c       get position of diag elements, error check
c
      call dmbsgetdiag(neq,neq1,maxa,na,jdiag,jdiagr,jdiagaz,L,lm)
c
c
c.....WRITE SYS FILE
c      WRITE(31) NUMELEM,NEQ,NUMETYP,MAXA
c.....WRITE JDIAG FILE
c      WRITE(35) (JDIAG(I),I=1,NEQ)
c      WRITE(35) (JDIAGAZ(I),I=1,NEQ1)
c      print *,iblk,'jdiag,jdiagaz =============='
c      print *,(jdiag(i),i=1,neq)
c      print *,(jdiagaz(i),i=1,neq1)
c      print *,'na ==== '
c      print *,(na(i),i=1,maxa)
c
c     write out Na file
c
c      write(36)(na(i),i=1,maxa)
c      do i = 1,maxa
c         na1(i)=na(i)
c      end do
      RETURN
      END
 
