      subroutine laplace(coorr,coefr,
     & prmt,estif,emass,edamp,eload,num)
c .... coorr ---- nodal coordinate value
c .... coefr ---- nodal coef value
      implicit real*8 (a-h,o-z)
      dimension estif(12,12),elump(12),emass(12),
     & eload(12)
      dimension prmt(*),
     & efuna(12),efunb(12),efunc(12),efund(12),
     & efune(12),efunf(12),coorr(3,4),coor(3)
      common /rlaplace/ru(4,16),rv(4,16),rw(4,16),
     & cu(4,4),cv(4,4),cw(4,4)
c .... store shape functions and their partial derivatives
c .... for all integral points
      common /vlaplace/rctr(3,3),crtr(3,3)
      common /dlaplace/ refc(3,4),gaus(4),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(12,3)
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
      common /area/ ar
      tm1=prmt(1)
      tm2=prmt(2)
      tm3=prmt(3)
      tm4=prmt(4)
      tm5=prmt(5)
      tm6=prmt(6)
      tm7=prmt(7)
      tm8=prmt(8)
      time=prmt(9)
      dt=prmt(10)
      imate=prmt(11)+0.5
      ielem=prmt(12)+0.5
      nelem=prmt(13)+0.5
      it=prmt(14)+0.5
      nmate=prmt(15)+0.5
      itime=prmt(16)+0.5
      ityp=prmt(17)+0.5
      smu=1.d0+ar
      slambda=smu
      pe=smu*(3.d0*slambda+2.d0*smu)/(slambda+smu)
      pv=slambda/(slambda+smu)/2.d0
      fact = pe/(1.d0+pv)/(1.d0-pv*2.d0)
      if (num.eq.1) call laplacei
c .... initialize the basic data
      do 10 i=1,nvar
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      do 999 igaus=1,ngaus
      call laplacet(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det)
c .... coordinate transfer from reference to original system
c .... rctr ---- Jacobi's matrix
c .... crtr ---- inverse matrix of Jacobi's matrix
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      iu=(igaus-1)*4+1
      iv=(igaus-1)*4+1
      iw=(igaus-1)*4+1
      if (num.gt.1) goto 2
c .... the following is the shape function caculation
      call laplace1(refc(1,igaus),ru(1,iu),rctr,crtr)
      call laplace2(refc(1,igaus),rv(1,iv),rctr,crtr)
      call laplace3(refc(1,igaus),rw(1,iw),rctr,crtr)
2     continue
c .... the following is the shape function transformation
c .... from reference coordinates to original coordinates
      call shapn(nrefc,ncoor,4,ru(1,iu),cu,crtr,1,4,4)
      call shapn(nrefc,ncoor,4,rv(1,iv),cv,crtr,1,4,4)
      call shapn(nrefc,ncoor,4,rw(1,iw),cw,crtr,1,4,4)
      weigh=det*gaus(igaus)
      do 100 i=1,12
      efuna(i) = 0.0
      efunb(i) = 0.0
      efunc(i) = 0.0
      efund(i) = 0.0
      efune(i) = 0.0
      efunf(i) = 0.0
100   continue
      do 101 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,2) 
      efuna(iv)=efuna(iv)+stif
101   continue
      do 102 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,3) 
      efunb(iv)=efunb(iv)+stif
102   continue
      do 103 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,4) 
      efunc(iv)=efunc(iv)+stif
103   continue
      do 104 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,4) 
      efund(iv)=efund(iv)+stif
104   continue
      do 105 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,3) 
      efund(iv)=efund(iv)+stif
105   continue
      do 106 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,4) 
      efune(iv)=efune(iv)+stif
106   continue
      do 107 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,2) 
      efune(iv)=efune(iv)+stif
107   continue
      do 108 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,3) 
      efunf(iv)=efunf(iv)+stif
108   continue
      do 109 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,2) 
      efunf(iv)=efunf(iv)+stif
109   continue
c .... the following is the stiffness computation
      do 202 iv=1,12
      do 201 jv=1,12
      stif=+efuna(iv)*efuna(jv)*fact*(1.d0-pv)
     & +efuna(iv)*efunb(jv)*fact*(pv)
     & +efuna(iv)*efunc(jv)*fact*(pv)
     & +efunb(iv)*efuna(jv)*fact*(pv)
     & +efunb(iv)*efunb(jv)*fact*(1.d0-pv)
     & +efunb(iv)*efunc(jv)*fact*(pv)
     & +efunc(iv)*efuna(jv)*fact*(pv)
     & +efunc(iv)*efunb(jv)*fact*(pv)
     & +efunc(iv)*efunc(jv)*fact*(1.d0-pv)
     & +efund(iv)*efund(jv)*fact*(0.5d0-pv)
     & +efune(iv)*efune(jv)*fact*(0.5d0-pv)
     & +efunf(iv)*efunf(jv)*fact*(0.5d0-pv)
      estif(iv,jv)=estif(iv,jv)+stif*weigh
201    continue
202    continue
c .... the following is the load vector computation
      do 501 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,1)*0.0
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,1)*0.0
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,1)*0.0
      eload(iv)=eload(iv)+stif*weigh
503   continue
999   continue
998   continue
      return
      end

      subroutine laplacei
      implicit real*8 (a-h,o-z)
      common /dlaplace/ refc(3,4),gaus(4),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(12,3)
c .... initial data
c .... refc ---- reference coordinates at integral points
c .... gaus ---- weight number at integral points
c .... nvard ---- the number of var for each unknown
c .... kdord ---- the highest differential order for each unknown
c .... kvord ---- var number at integral points for each unknown
      ngaus=  4
      ndisp=  3
      nrefc=  3
      ncoor=  3
      nvar = 12
      nnode=  4
      kdord(1)=1
      nvard(1)=4
      kvord(1,1)=1
      kvord(2,1)=4
      kvord(3,1)=7
      kvord(4,1)=10
      kdord(2)=1
      nvard(2)=4
      kvord(1,2)=2
      kvord(2,2)=5
      kvord(3,2)=8
      kvord(4,2)=11
      kdord(3)=1
      nvard(3)=4
      kvord(1,3)=3
      kvord(2,3)=6
      kvord(3,3)=9
      kvord(4,3)=12
      refc(1,1)=0.5854102
      refc(2,1)=0.1381966
      refc(3,1)=0.1381966
      gaus(1)=0.25
      refc(1,2)=0.1381966
      refc(2,2)=0.5854102
      refc(3,2)=0.1381966
      gaus(2)=0.25
      refc(1,3)=0.1381966
      refc(2,3)=0.1381966
      refc(3,3)=0.5854102
      gaus(3)=0.25
      refc(1,4)=0.1381966
      refc(2,4)=0.1381966
      refc(3,4)=0.1381966
      gaus(4)=0.25
      end


      subroutine laplacet(nnode,nrefc,ncoor,refc,coor,coorr,
     & rc,cr,det)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     & coorr(ncoor,nnode),coor(ncoor)
      call tlaplace(refc,coor,coorr,rc)
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
 

      subroutine laplace1(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external flaplace1
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(flaplace1,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function flaplace1(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /cclaplace/ xa(4),ya(4),za(4)
      common /vlaplace/ rctr(3,3),crtr(3,3)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     flaplace1=+rx 
      goto 1000
2     flaplace1=+ry 
      goto 1000
3     flaplace1=+rz 
      goto 1000
4     flaplace1=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine laplace2(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external flaplace2
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(flaplace2,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function flaplace2(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /cclaplace/ xa(4),ya(4),za(4)
      common /vlaplace/ rctr(3,3),crtr(3,3)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     flaplace2=+rx 
      goto 1000
2     flaplace2=+ry 
      goto 1000
3     flaplace2=+rz 
      goto 1000
4     flaplace2=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine laplace3(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external flaplace3
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(flaplace3,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function flaplace3(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /cclaplace/ xa(4),ya(4),za(4)
      common /vlaplace/ rctr(3,3),crtr(3,3)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     flaplace3=+rx 
      goto 1000
2     flaplace3=+ry 
      goto 1000
3     flaplace3=+rz 
      goto 1000
4     flaplace3=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine tlaplace(refc,coor,coorr,rc)
c .... compute coordinate value and Jacobi's matrix rc
c .... by reference coordinate value
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,4),rc(3,3)
      common /cclaplace/ x(4),y(4),z(4)
      external ftlaplace
      do 100 n=1,4
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(ftlaplace,refc,coor,rc,3,3,1)
c .... coordinate value and their partial derivatives caculation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function ftlaplace(refc,n)
c .... coordinate transfer function caculation
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /cclaplace/ x(4),y(4),z(4)
      common /vlaplace/ rctr(3,3),crtr(3,3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     ftlaplace=+(+rx)*x(1)+(+ry)*x(2)+(+rz)*x(3)+(+(+1.-rx-
     & ry-rz))*x(4)
      goto 1000
2     ftlaplace=+(+rx)*y(1)+(+ry)*y(2)+(+rz)*y(3)+(+(+1.-rx-
     & ry-rz))*y(4)
      goto 1000
3     ftlaplace=+(+rx)*z(1)+(+ry)*z(2)+(+rz)*z(3)+(+(+1.-rx-
     & ry-rz))*z(4)
      goto 1000
1000  return
      end

