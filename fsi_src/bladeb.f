      subroutine bladeb(coorr,coefr,
     & prmt,estif,emass,edamp,eload,num)
c .... coorr ---- nodal coordinate value
c .... coefr ---- nodal coef value
      implicit real*8 (a-h,o-z)
      dimension estif(12,12),elump(12),emass(12),
     & eload(12)
      dimension prmt(*),coef(17),coefr(3,17),
     & efuna(12),coorr(2,3),coor(2)
      common /rbladeb/ru(3,21),rv(3,21),rw(3,21),
     & rp(3,21),
     & cu(3,3),cv(3,3),cw(3,3),cp(3,3)
c .... store shape functions and their partial derivatives
c .... for all integral points
      common /vbladeb/rctr(2,2),crtr(2,2),coefd(17,5),coefc(17,5)
      common /dbladeb/ refc(2,7),gaus(7),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(4),kdord(4),kvord(12,4)
c .... nnode ---- the number of nodes
c .... nrefc ---- the number of numerical integral points
c .... ndisp ---- the number of unknown functions
c .... nrefc ---- the number of reference coordinates
c .... nvar ---- the number of unknown varibles var
c .... refc ---- reference coordinates at integral points
c .... gaus ---- weight number at integral points
c .... nvard ---- the number of var for each unknown
c .... kdord ---- the highest differential order for each unknown
c .... kvord ---- var number at integral points for each unknown
      stiffness=prmt(1)
      time=prmt(2)
      dt=prmt(3)
      imate=prmt(4)+0.5
      ielem=prmt(5)+0.5
      nelem=prmt(6)+0.5
      it=prmt(7)+0.5
      nmate=prmt(8)+0.5
      itime=prmt(9)+0.5
      ityp=prmt(10)+0.5
      onx=0.d0
      ony=0.d0
      onz=1.d0
      if (num.eq.1) call bladebi
c .... initialize the basic data
      do 10 i=1,nvar
      emass(i)=0.0
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      do 999 igaus=1,ngaus
      call bladebt(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
c .... coordinate transfer from reference to original system
c .... rctr ---- Jacobi's matrix
c .... crtr ---- inverse matrix of Jacobi's matrix
      x=coor(1)
      y=coor(2)
      r=refc(1,igaus)
      s=refc(2,igaus)
      call ebladeb(refc(1,igaus),coef,coorr,coefr,coefd)
c .... compute coef functions and their partial derivatives
      iu=(igaus-1)*3+1
      iv=(igaus-1)*3+1
      iw=(igaus-1)*3+1
      ip=(igaus-1)*3+1
      if (num.gt.1) goto 2
c .... the following is the shape function caculation
      call bladeb1(refc(1,igaus),ru(1,iu),rctr,crtr)
      call bladeb2(refc(1,igaus),rv(1,iv),rctr,crtr)
      call bladeb3(refc(1,igaus),rw(1,iw),rctr,crtr)
      call bladeb4(refc(1,igaus),rp(1,ip),rctr,crtr)
2     continue
c .... the following is the shape function transformation
c .... from reference coordinates to original coordinates
      call shapn(nrefc,ncoor,3,ru(1,iu),cu,crtr,1,3,3)
      call shapn(nrefc,ncoor,3,rv(1,iv),cv,crtr,1,3,3)
      call shapn(nrefc,ncoor,3,rw(1,iw),cw,crtr,1,3,3)
      call shapn(nrefc,ncoor,3,rp(1,ip),cp,crtr,1,3,3)
c .... the coef function transformation
c .... from reference coordinates to original coordinates
      call shapc(nrefc,ncoor,17,coefd,coefc,crtr,2,5,5)
      un=coef(1)
      vn=coef(2)
      wn=coef(3)
      pn=coef(4)
      unn=coef(5)
      vnn=coef(6)
      wnn=coef(7)
      pnn=coef(8)
      umn=coef(9)
      vmn=coef(10)
      wmn=coef(11)
      vlx=coef(12)
      vly=coef(13)
      vlz=coef(14)
      gx=coef(15)
      gy=coef(16)
      gz=coef(17)
      weigh=det*gaus(igaus)
      do 100 i=1,12
      efuna(i) = 0.0
100   continue
      do 101 i=1,3
      iv=kvord(i,1)
      stif=+cu(i,1)*onx
      efuna(iv)=efuna(iv)+stif
101   continue
      do 102 i=1,3
      iv=kvord(i,2)
      stif=+cv(i,1)*ony
      efuna(iv)=efuna(iv)+stif
102   continue
      do 103 i=1,3
      iv=kvord(i,3)
      stif=+cw(i,1)*onz
      efuna(iv)=efuna(iv)+stif
103   continue
c .... the following is the stiffness computation
      do 202 iv=1,12
      do 201 jv=1,12
      stif=+efuna(iv)*efuna(jv)*stiffness
      estif(iv,jv)=estif(iv,jv)+stif*weigh
201    continue
202    continue
c .... the following is the mass matrix computation
      stif=0.d0
      elump(1)=stif*weigh
      stif=0.d0
      elump(5)=stif*weigh
      stif=0.d0
      elump(9)=stif*weigh
      stif=0.d0
      elump(2)=stif*weigh
      stif=0.d0
      elump(6)=stif*weigh
      stif=0.d0
      elump(10)=stif*weigh
      stif=0.d0
      elump(3)=stif*weigh
      stif=0.d0
      elump(7)=stif*weigh
      stif=0.d0
      elump(11)=stif*weigh
      stif=0.d0
      elump(4)=stif*weigh
      stif=0.d0
      elump(8)=stif*weigh
      stif=0.d0
      elump(12)=stif*weigh
      do 301 i=1,nvard(1)
      iv = kvord(i,1)
      emass(iv)=emass(iv)+elump(iv)*cu(i,1)
301   continue
      do 302 i=1,nvard(2)
      iv = kvord(i,2)
      emass(iv)=emass(iv)+elump(iv)*cv(i,1)
302   continue
      do 303 i=1,nvard(3)
      iv = kvord(i,3)
      emass(iv)=emass(iv)+elump(iv)*cw(i,1)
303   continue
      do 304 i=1,nvard(4)
      iv = kvord(i,4)
      emass(iv)=emass(iv)+elump(iv)*cp(i,1)
304   continue
c .... the following is the load vector computation
      do 501 i=1,3
      iv=kvord(i,1)
      stif=+cu(i,1)*0.d0
      eload(iv)=eload(iv)+stif*weigh
501   continue
999   continue
998   continue
      return
      end

      subroutine bladebi
      implicit real*8 (a-h,o-z)
      common /dbladeb/ refc(2,7),gaus(7),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(4),kdord(4),kvord(12,4)
c .... initial data
c .... refc ---- reference coordinates at integral points
c .... gaus ---- weight number at integral points
c .... nvard ---- the number of var for each unknown
c .... kdord ---- the highest differential order for each unknown
c .... kvord ---- var number at integral points for each unknown
      ngaus=  7
      ndisp=  4
      nrefc=  2
      ncoor=  2
      nvar = 12
      nnode=  3
      kdord(1)=1
      nvard(1)=3
      kvord(1,1)=1
      kvord(2,1)=5
      kvord(3,1)=9
      kdord(2)=1
      nvard(2)=3
      kvord(1,2)=2
      kvord(2,2)=6
      kvord(3,2)=10
      kdord(3)=1
      nvard(3)=3
      kvord(1,3)=3
      kvord(2,3)=7
      kvord(3,3)=11
      kdord(4)=1
      nvard(4)=3
      kvord(1,4)=4
      kvord(2,4)=8
      kvord(3,4)=12
      refc(1,1)=1.
      refc(2,1)=0.
      gaus(1)=1./20.*0.5
      refc(1,2)=0.5
      refc(2,2)=0.5
      gaus(2)=2./15.*0.5
      refc(1,3)=0.
      refc(2,3)=1.
      gaus(3)=1./20.*0.5
      refc(1,4)=0.
      refc(2,4)=0.5
      gaus(4)=2./15.*0.5
      refc(1,5)=0.
      refc(2,5)=0.
      gaus(5)=1./20.*0.5
      refc(1,6)=0.5
      refc(2,6)=0.
      gaus(6)=2./15.*0.5
      refc(1,7)=1./3.
      refc(2,7)=1./3.
      gaus(7)=9./20.*0.5
      end


      subroutine bladebt(nnode,nrefc,ncoor,refc,coor,coorr,
     & rc,cr,det,coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     & coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call tbladeb(refc,coor,coorr,coefr,rc)
      n=nrefc
      m=n*2
      det = 1.0
      do 10 i=1,n
      do 10 j=1,n
      if (i.le.ncoor) a(i,j) = rc(i,j)
      if (i.gt.ncoor) a(i,j)=1.0
      a(i,n+j)=0.0
      if (i.eq.j) a(i,n+i) = 1.0
10    continue
c     write(*,*) 'a ='
c     do 21 i=1,n
c21   write(*,8) (a(i,j),j=1,m)
      do 400 i=1,n
      amax = 0.0
      l = 0
      do 50 j=i,n
      c = a(j,i)
      if (c.lt.0.0) c = -c
      if (c.le.amax) goto 50
      amax = c
      l = j
50    continue
      do 60 k=1,m
      c = a(l,k)
      a(l,k) = a(i,k)
      a(i,k) = c
60    continue
      c = a(i,i)
      det = c*det
      do 100 k=i+1,m
100   a(i,k) = a(i,k)/c
      do 300 j=1,n
      if (i.eq.j) goto 300
      do 200 k=i+1,m
200   a(j,k) = a(j,k)-a(i,k)*a(j,i)
c     write(*,*) 'i =',i,'  j =',j,'  a ='
c     do 11 ii=1,n
c11   write(*,8) (a(ii,jj),jj=1,m)
300   continue
400   continue
      do 500 i=1,nrefc
      do 500 j=1,ncoor
500   cr(i,j) = a(i,n+j)
c     write(*,*) 'a ='
c     do 22 i=1,n
c22   write(*,8) (a(i,j),j=1,m)
c     write(*,*) 'rc ='
c     do 24 i=1,ncoor
c24   write(*,8) (rc(i,j),j=1,nrefc)
c     write(*,*) 'cr ='
c     do 23 i=1,nrefc
c23   write(*,8) (cr(i,j),j=1,ncoor)
c     write(*,*) 'det =',det
      if (det.lt.0.0) det=-det
c     write(*,*) 'det =',det
8     format(1x,6f12.3)
      end
 

      subroutine bladeb1(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(2),shpr(3,3),rctr(2,2),crtr(2,2)
      external fbladeb1
      r=refc(1)
      s=refc(2)
      call dshap(fbladeb1,refc,shpr,2,3,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fbladeb1(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccbladeb/ xa(3),ya(3),una(3),vna(3),
     &wna(3),pna(3),unna(3),vnna(3),wnna(3),pnna(3),
     &umna(3),vmna(3),wmna(3),vlxa(3),vlya(3),vlza(3),
     &gxa(3),gya(3),gza(3)
      common /vbladeb/ rctr(2,2),crtr(2,2),coefd(17,5),coefc(17,5)
      dimension refc(2)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      r=refc(1)
      s=refc(2)
      goto (1,2,3) n
1     fbladeb1=+r 
      goto 1000
2     fbladeb1=+s 
      goto 1000
3     fbladeb1=+(+1.-r-s) 
      goto 1000
1000  return
      end

      subroutine bladeb2(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(2),shpr(3,3),rctr(2,2),crtr(2,2)
      external fbladeb2
      r=refc(1)
      s=refc(2)
      call dshap(fbladeb2,refc,shpr,2,3,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fbladeb2(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccbladeb/ xa(3),ya(3),una(3),vna(3),
     &wna(3),pna(3),unna(3),vnna(3),wnna(3),pnna(3),
     &umna(3),vmna(3),wmna(3),vlxa(3),vlya(3),vlza(3),
     &gxa(3),gya(3),gza(3)
      common /vbladeb/ rctr(2,2),crtr(2,2),coefd(17,5),coefc(17,5)
      dimension refc(2)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      r=refc(1)
      s=refc(2)
      goto (1,2,3) n
1     fbladeb2=+r 
      goto 1000
2     fbladeb2=+s 
      goto 1000
3     fbladeb2=+(+1.-r-s) 
      goto 1000
1000  return
      end

      subroutine bladeb3(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(2),shpr(3,3),rctr(2,2),crtr(2,2)
      external fbladeb3
      r=refc(1)
      s=refc(2)
      call dshap(fbladeb3,refc,shpr,2,3,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fbladeb3(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccbladeb/ xa(3),ya(3),una(3),vna(3),
     &wna(3),pna(3),unna(3),vnna(3),wnna(3),pnna(3),
     &umna(3),vmna(3),wmna(3),vlxa(3),vlya(3),vlza(3),
     &gxa(3),gya(3),gza(3)
      common /vbladeb/ rctr(2,2),crtr(2,2),coefd(17,5),coefc(17,5)
      dimension refc(2)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      r=refc(1)
      s=refc(2)
      goto (1,2,3) n
1     fbladeb3=+r 
      goto 1000
2     fbladeb3=+s 
      goto 1000
3     fbladeb3=+(+1.-r-s) 
      goto 1000
1000  return
      end

      subroutine bladeb4(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(2),shpr(3,3),rctr(2,2),crtr(2,2)
      external fbladeb4
      r=refc(1)
      s=refc(2)
      call dshap(fbladeb4,refc,shpr,2,3,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fbladeb4(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccbladeb/ xa(3),ya(3),una(3),vna(3),
     &wna(3),pna(3),unna(3),vnna(3),wnna(3),pnna(3),
     &umna(3),vmna(3),wmna(3),vlxa(3),vlya(3),vlza(3),
     &gxa(3),gya(3),gza(3)
      common /vbladeb/ rctr(2,2),crtr(2,2),coefd(17,5),coefc(17,5)
      dimension refc(2)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      r=refc(1)
      s=refc(2)
      goto (1,2,3) n
1     fbladeb4=+r 
      goto 1000
2     fbladeb4=+s 
      goto 1000
3     fbladeb4=+(+1.-r-s) 
      goto 1000
1000  return
      end

      subroutine ebladeb(refc,coef,coorr,coefr,coefd)
c .... compute coef value and their partial derivatives
c .... by reference coordinate value
      implicit real*8 (a-h,o-z)
      dimension refc(2),coef(17),coorr(2,3),coefr(3,17),coefd(17,2)
      external febladeb
      r=refc(1)
      s=refc(2)
      call dcoef(febladeb,refc,coef,coefd,2,17,2)
c .... coef value and their partial derivatives caculation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function febladeb(refc,n)
c .... coef function caculation
      implicit real*8 (a-h,o-z)
      dimension refc(2)
      common /ccbladeb/ xa(3),ya(3),un(3),vn(3),wn(3),pn(3),
     &unn(3),vnn(3),wnn(3),pnn(3),umn(3),vmn(3),
     &wmn(3),vlx(3),vly(3),vlz(3),gx(3),gy(3),
     &gz(3)
      common /vbladeb/ rctr(2,2),crtr(2,2),coefd(17,5),coefc(17,5)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      r=refc(1)
      s=refc(2)
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     & 16,17) n
1     febladeb=+(+r)*un(1)+(+s)*un(2)+(+(+1.-r-s))*un(3)
      goto 1000
2     febladeb=+(+r)*vn(1)+(+s)*vn(2)+(+(+1.-r-s))*vn(3)
      goto 1000
3     febladeb=+(+r)*wn(1)+(+s)*wn(2)+(+(+1.-r-s))*wn(3)
      goto 1000
4     febladeb=+(+r)*pn(1)+(+s)*pn(2)+(+(+1.-r-s))*pn(3)
      goto 1000
5     febladeb=+(+r)*unn(1)+(+s)*unn(2)+(+(+1.-r-s))*unn(3)
      goto 1000
6     febladeb=+(+r)*vnn(1)+(+s)*vnn(2)+(+(+1.-r-s))*vnn(3)
      goto 1000
7     febladeb=+(+r)*wnn(1)+(+s)*wnn(2)+(+(+1.-r-s))*wnn(3)
      goto 1000
8     febladeb=+(+r)*pnn(1)+(+s)*pnn(2)+(+(+1.-r-s))*pnn(3)
      goto 1000
9     febladeb=+(+r)*umn(1)+(+s)*umn(2)+(+(+1.-r-s))*umn(3)
      goto 1000
10     febladeb=+(+r)*vmn(1)+(+s)*vmn(2)+(+(+1.-r-s))*vmn(3)
      goto 1000
11     febladeb=+(+r)*wmn(1)+(+s)*wmn(2)+(+(+1.-r-s))*wmn(3)
      goto 1000
12     febladeb=+(+r)*vlx(1)+(+s)*vlx(2)+(+(+1.-r-s))*vlx(3)
      goto 1000
13     febladeb=+(+r)*vly(1)+(+s)*vly(2)+(+(+1.-r-s))*vly(3)
      goto 1000
14     febladeb=+(+r)*vlz(1)+(+s)*vlz(2)+(+(+1.-r-s))*vlz(3)
      goto 1000
15     febladeb=+(+r)*gx(1)+(+s)*gx(2)+(+(+1.-r-s))*gx(3)
      goto 1000
16     febladeb=+(+r)*gy(1)+(+s)*gy(2)+(+(+1.-r-s))*gy(3)
      goto 1000
17     febladeb=+(+r)*gz(1)+(+s)*gz(2)+(+(+1.-r-s))*gz(3)
      goto 1000
1000  return
      end

      subroutine tbladeb(refc,coor,coorr,coefr,rc)
c .... compute coordinate value and Jacobi's matrix rc
c .... by reference coordinate value
      implicit real*8 (a-h,o-z)
      dimension refc(2),coor(2),coorr(2,3),coefr(3,17),rc(2,2)
      common /ccbladeb/ x(3),y(3),un(3),vn(3),wn(3),pn(3),
     &unn(3),vnn(3),wnn(3),pnn(3),umn(3),vmn(3),
     &wmn(3),vlx(3),vly(3),vlz(3),gx(3),gy(3),
     &gz(3)
      external ftbladeb
      do 100 n=1,3
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
100   continue
      do 200 n=1,3
      un(n)=coefr(n,1)
      vn(n)=coefr(n,2)
      wn(n)=coefr(n,3)
      pn(n)=coefr(n,4)
      unn(n)=coefr(n,5)
      vnn(n)=coefr(n,6)
      wnn(n)=coefr(n,7)
      pnn(n)=coefr(n,8)
      umn(n)=coefr(n,9)
      vmn(n)=coefr(n,10)
      wmn(n)=coefr(n,11)
      vlx(n)=coefr(n,12)
      vly(n)=coefr(n,13)
      vlz(n)=coefr(n,14)
      gx(n)=coefr(n,15)
      gy(n)=coefr(n,16)
      gz(n)=coefr(n,17)
200   continue
      r=refc(1)
      s=refc(2)
      call dcoor(ftbladeb,refc,coor,rc,2,2,1)
c .... coordinate value and their partial derivatives caculation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function ftbladeb(refc,n)
c .... coordinate transfer function caculation
      implicit real*8 (a-h,o-z)
      dimension refc(2)
      common /ccbladeb/ x(3),y(3),un(3),vn(3),wn(3),pn(3),
     &unn(3),vnn(3),wnn(3),pnn(3),umn(3),vmn(3),
     &wmn(3),vlx(3),vly(3),vlz(3),gx(3),gy(3),
     &gz(3)
      common /vbladeb/ rctr(2,2),crtr(2,2),coefd(17,5),coefc(17,5)
      r=refc(1)
      s=refc(2)
      goto (1,2) n
1     ftbladeb=+(+r)*x(1)+(+s)*x(2)+(+(+1.-r-s))*x(3)
      goto 1000
2     ftbladeb=+(+r)*y(1)+(+s)*y(2)+(+(+1.-r-s))*y(3)
      goto 1000
1000  return
      end

