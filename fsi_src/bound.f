      real*8 function curve(x,s,t,a)
      implicit real*8 (a-h,o-z)
c
cc Parabolic curve
c
      curve=-a*(x+(s+t)/4.d0)**2+a*(s-t)**2/16.d0
      RETURN
      END

      real*8 function bound(r,u,hx,hy,t,pi,j)
      implicit real*8 (a-h,o-z)
      DIMENSION R(3)
      goto (1,2,3,4),j
1     bound=0.d0
      goto 5
2     bound=1.5d0*u*(r(1)*(hx-r(1))/(hx/2.d0)**2)*
     &      (r(3)*(hy-r(3))/(hy/2.d0)**2)
      if (t.le.2.d0-1.d-10) 
     & bound=bound*(1.d0-dcos(pi/2.d0*t))/2.d0
      goto 5
3     bound=0.d0
      goto 5
4     bound=0.d0
5     continue
      RETURN
      END

      real*8 function boundpro(r,u,hx,hy,t,pi,j)
      implicit real*8 (a-h,o-z)
      DIMENSION R(3)
      goto (1,2,3,4),j
1     boundpro=0.d0
      goto 5
2     boundpro=1.5d0*u*(1.d0-(r(1)**2+r(3)**2)/hx**2)
      if (t.le.2.d0-1.d-10) 
     & boundpro=boundpro*(1.d0-dcos(pi/2.d0*t))/2.d0
      goto 5
3     boundpro=0.d0
      goto 5
4     boundpro=0.d0
5     continue
      RETURN
      END
