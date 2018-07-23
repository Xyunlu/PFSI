      subroutine inertia(coorr,coefr,
     & prmt,estif,emass,edamp,eload,num)
c .... coorr ---- nodal coordinate value
c .... coefr ---- nodal coef value
      implicit real*8 (a-h,o-z)
      dimension estif(12,12),elump(12),emass(12),
     & eload(12)
      dimension prmt(*),coef(6),coefr(4,6),
     & coorr(3,4),coor(3)
      common /rinertia/ru(4,16),rv(4,16),rw(4,16),
     & cu(4,4),cv(4,4),cw(4,4)
c .... store shape functions and their partial derivatives
c .... for all integral points
      common /vinertia/rctr(3,3),crtr(3,3),coefd(6,9),coefc(6,9)
      common /dinertia/ refc(3,4),gaus(4),
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
      common /rdata/ xo,yo,zo
      common /idata/ id4updown
      t1=prmt(1)
      t2=prmt(2)
      t3=prmt(3)
      t4=prmt(4)
      dens=prmt(5)
      t6=prmt(6)
      t7=prmt(7)
      t8=prmt(8)
      time=prmt(9)
      dt=prmt(10)
      imate=prmt(11)+0.5
      ielem=prmt(12)+0.5
      nelem=prmt(13)+0.5
      it=prmt(14)+0.5
      nmate=prmt(15)+0.5
      itime=prmt(16)+0.5
      ityp=prmt(17)+0.5
      if (num.eq.1) call inertiai
c .... initialize the basic data
      do 10 i=1,nvar
      emass(i)=0.0
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      do 999 igaus=1,ngaus
      call inertiat(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
c .... coordinate transfer from reference to original system
c .... rctr ---- Jacobi's matrix
c .... crtr ---- inverse matrix of Jacobi's matrix
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call einertia(refc(1,igaus),coef,coorr,coefr,coefd)
c .... compute coef functions and their partial derivatives
      iu=(igaus-1)*4+1
      iv=(igaus-1)*4+1
      iw=(igaus-1)*4+1
      if (num.gt.1) goto 2
c .... the following is the shape function caculation
      call inertia1(refc(1,igaus),ru(1,iu),rctr,crtr)
      call inertia2(refc(1,igaus),rv(1,iv),rctr,crtr)
      call inertia3(refc(1,igaus),rw(1,iw),rctr,crtr)
2     continue
c .... the following is the shape function transformation
c .... from reference coordinates to original coordinates
      call shapn(nrefc,ncoor,4,ru(1,iu),cu,crtr,1,4,4)
      call shapn(nrefc,ncoor,4,rv(1,iv),cv,crtr,1,4,4)
      call shapn(nrefc,ncoor,4,rw(1,iw),cw,crtr,1,4,4)
c .... the coef function transformation
c .... from reference coordinates to original coordinates
      call shapc(nrefc,ncoor,6,coefd,coefc,crtr,2,9,9)
      ud=coef(1)
      vd=coef(2)
      wd=coef(3)
      us=coef(4)
      vs=coef(5)
      ws=coef(6)
      weigh=det*gaus(igaus)
      xct=x-xo+ud
      yct=y-yo+vd
      zct=z-zo+wd
      detf=(1.d0+coefc(1,1))*(1.d0+coefc(2,2))
     * *(1.d0+coefc(3,3))
      detf=detf+coefc(2,1)*coefc(3,2)
     * *coefc(1,3)
      detf=detf+coefc(3,1)*coefc(1,2)
     * *coefc(2,3)
      detf=detf-(1.d0+coefc(1,1))*coefc(3,2)
     * *coefc(2,3)
      detf=detf-coefc(2,1)*coefc(1,2)
     * *(1.d0+coefc(3,3))
      detf=detf-coefc(3,1)*(1.d0+coefc(2,2))
     * *coefc(1,3)
      if (id4updown.eq.0) then
      sinertia=xct**2+zct**2
      tarque=zct*us-xct*ws
      elseif (id4updown.eq.1) then
      sinertia=xct**2+yct**2
      tarque=xct*vs-yct*us
      elseif (id4updown.eq.-1) then
      sinertia=yct**2+zct**2
      tarque=yct*ws-zct*vs
      endif
c .... the following is the stiffness computation
      do 202 i=1,4
      iv=kvord(i,1)
      do 201 j=1,4
      jv=kvord(j,1)
      stif=+cu(i,1)*cu(j,1)*0.0
      estif(iv,jv)=estif(iv,jv)+stif*weigh
201    continue
202    continue
c .... the following is the mass matrix computation
      stif=1.d0
      elump(1)=stif*weigh
      stif=1.d0
      elump(4)=stif*weigh
      stif=1.d0
      elump(7)=stif*weigh
      stif=1.d0
      elump(10)=stif*weigh
      stif=1.d0
      elump(2)=stif*weigh
      stif=1.d0
      elump(5)=stif*weigh
      stif=1.d0
      elump(8)=stif*weigh
      stif=1.d0
      elump(11)=stif*weigh
      stif=1.d0
      elump(3)=stif*weigh
      stif=1.d0
      elump(6)=stif*weigh
      stif=1.d0
      elump(9)=stif*weigh
      stif=1.d0
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
c .... the following is the load vector computation
      stif=sinertia*dens*detf
      edamp=edamp+stif*weigh

      stif=tarque*dens*detf
      eload(1)=eload(1)+stif*weigh

999   continue
998   continue
      return
      end

      subroutine inertiai
      implicit real*8 (a-h,o-z)
      common /dinertia/ refc(3,4),gaus(4),
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
      refc(1,1)=0.5854102d0
      refc(2,1)=0.1381966d0
      refc(3,1)=0.1381966d0
      gaus(1)=0.25d0
      refc(1,2)=0.1381966d0
      refc(2,2)=0.5854102d0
      refc(3,2)=0.1381966d0
      gaus(2)=0.25d0
      refc(1,3)=0.1381966d0
      refc(2,3)=0.1381966d0
      refc(3,3)=0.5854102d0
      gaus(3)=0.25d0
      refc(1,4)=0.1381966d0
      refc(2,4)=0.1381966d0
      refc(3,4)=0.1381966d0
      gaus(4)=0.25d0
      end


      subroutine inertiat(nnode,nrefc,ncoor,refc,coor,coorr,
     & rc,cr,det,coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     & coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call tinertia(refc,coor,coorr,coefr,rc)
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
 

      subroutine inertia1(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external finertia1
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(finertia1,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function finertia1(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccinertia/ xa(4),ya(4),za(4),uda(4),
     &vda(4),wda(4),usa(4),vsa(4),wsa(4)
      common /vinertia/ rctr(3,3),crtr(3,3),coefd(6,9),coefc(6,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     finertia1=+rx 
      goto 1000
2     finertia1=+ry 
      goto 1000
3     finertia1=+rz 
      goto 1000
4     finertia1=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine inertia2(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external finertia2
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(finertia2,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function finertia2(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccinertia/ xa(4),ya(4),za(4),uda(4),
     &vda(4),wda(4),usa(4),vsa(4),wsa(4)
      common /vinertia/ rctr(3,3),crtr(3,3),coefd(6,9),coefc(6,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     finertia2=+rx 
      goto 1000
2     finertia2=+ry 
      goto 1000
3     finertia2=+rz 
      goto 1000
4     finertia2=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine inertia3(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external finertia3
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(finertia3,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function finertia3(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccinertia/ xa(4),ya(4),za(4),uda(4),
     &vda(4),wda(4),usa(4),vsa(4),wsa(4)
      common /vinertia/ rctr(3,3),crtr(3,3),coefd(6,9),coefc(6,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     finertia3=+rx 
      goto 1000
2     finertia3=+ry 
      goto 1000
3     finertia3=+rz 
      goto 1000
4     finertia3=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine einertia(refc,coef,coorr,coefr,coefd)
c .... compute coef value and their partial derivatives
c .... by reference coordinate value
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(6),coorr(3,4),coefr(4,6),coefd(6,3)
      external feinertia
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(feinertia,refc,coef,coefd,3,6,2)
c .... coef value and their partial derivatives caculation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function feinertia(refc,n)
c .... coef function caculation
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccinertia/ xa(4),ya(4),za(4),ud(4),vd(4),wd(4),
     &us(4),vs(4),ws(4)
      common /vinertia/ rctr(3,3),crtr(3,3),coefd(6,9),coefc(6,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     feinertia=+(+rx)*ud(1)+(+ry)*ud(2)+(+rz)*ud(3)+(+(+1.-
     & rx-ry-rz))*ud(4)
      goto 1000
2     feinertia=+(+rx)*vd(1)+(+ry)*vd(2)+(+rz)*vd(3)+(+(+1.-
     & rx-ry-rz))*vd(4)
      goto 1000
3     feinertia=+(+rx)*wd(1)+(+ry)*wd(2)+(+rz)*wd(3)+(+(+1.-
     & rx-ry-rz))*wd(4)
      goto 1000
4     feinertia=+(+rx)*us(1)+(+ry)*us(2)+(+rz)*us(3)+(+(+1.-
     & rx-ry-rz))*us(4)
      goto 1000
5     feinertia=+(+rx)*vs(1)+(+ry)*vs(2)+(+rz)*vs(3)+(+(+1.-
     & rx-ry-rz))*vs(4)
      goto 1000
6     feinertia=+(+rx)*ws(1)+(+ry)*ws(2)+(+rz)*ws(3)+(+(+1.-
     & rx-ry-rz))*ws(4)
      goto 1000
1000  return
      end

      subroutine tinertia(refc,coor,coorr,coefr,rc)
c .... compute coordinate value and Jacobi's matrix rc
c .... by reference coordinate value
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,4),coefr(4,6),rc(3,3)
      common /ccinertia/ x(4),y(4),z(4),ud(4),vd(4),wd(4),
     &us(4),vs(4),ws(4)
      external ftinertia
      do 100 n=1,4
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,4
      ud(n)=coefr(n,1)
      vd(n)=coefr(n,2)
      wd(n)=coefr(n,3)
      us(n)=coefr(n,4)
      vs(n)=coefr(n,5)
      ws(n)=coefr(n,6)
200   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(ftinertia,refc,coor,rc,3,3,1)
c .... coordinate value and their partial derivatives caculation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function ftinertia(refc,n)
c .... coordinate transfer function caculation
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccinertia/ x(4),y(4),z(4),ud(4),vd(4),wd(4),
     &us(4),vs(4),ws(4)
      common /vinertia/ rctr(3,3),crtr(3,3),coefd(6,9),coefc(6,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     ftinertia=+(+rx)*x(1)+(+ry)*x(2)+(+rz)*x(3)+(+(+1.-rx-
     & ry-rz))*x(4)
      goto 1000
2     ftinertia=+(+rx)*y(1)+(+ry)*y(2)+(+rz)*y(3)+(+(+1.-rx-
     & ry-rz))*y(4)
      goto 1000
3     ftinertia=+(+rx)*z(1)+(+ry)*z(2)+(+rz)*z(3)+(+(+1.-rx-
     & ry-rz))*z(4)
      goto 1000
1000  return
      end

