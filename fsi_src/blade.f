      subroutine blade(coorr,coefr,
     & prmt,estif,emass,edamp,eload,num,iblk)
c .... coorr ---- nodal coordinate value
c .... coefr ---- nodal coef value
      implicit real*8 (a-h,o-z)
      dimension estif(16,16),elump(16),emass(16),
     & edamp(16),eload(16)
      dimension prmt(*),coef(14),coefr(4,14),
     & efuna(16),efunb(16),efunc(16),efund(16),
     & efune(16),efunf(16),coorr(3,4),coor(3)
      common /rblade/ru(4,16),rv(4,16),rw(4,16),
     & rp(4,16),
     & cu(4,4),cv(4,4),cw(4,4),cp(4,4)
c .... store shape functions and their partial derivatives
c .... for all integral points
      common /vblade/rctr(3,3),crtr(3,3),coefd(14,9),coefc(14,9)
      common /dblade/ refc(3,4),gaus(4),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(4),kdord(4),kvord(16,4)
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
      common /rdata/ vangle,angle,penalty,pen4sd
      common /idata/ id4updown

c      if( iblk.eq.1 .and. num.eq.64296) then
c        print *,'94296-coefr:',coefr
c      endif
      visc_kinema=prmt(1)
      pe=prmt(2)
      pv=prmt(3)
      denf=prmt(4)
      dens=prmt(5)
      penalty=prmt(6)
      gravity=prmt(7)
      angle0=prmt(8)
      time=prmt(9)
      dt=prmt(10)
      imate=prmt(11)+0.5
      ielem=prmt(12)+0.5
      nelem=prmt(13)+0.5
      it=prmt(14)+0.5
      nmate=prmt(15)+0.5
      itime=prmt(16)+0.5
      ityp=prmt(17)+0.5
c......
      ifluid=0
      isolid=0
      if (imate.le.2) then
      ifluid=1
      den=denf
      else
      isolid=1
      den=dens
      endif
      visc=visc_kinema*denf
      fact = pe/(1.d0+pv)/(1.d0-pv*2.d0)*isolid
      fx=0.d0
      fy=0.d0
      fz=0.d0
      if (id4updown.eq.1) fz=-den*gravity
      if (num.eq.1) call bladei
c .... initialize the basic data
      do 10 i=1,nvar
      emass(i)=0.0
      edamp(i)=0.0
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
c      if( iblk.eq.1 .and. num.eq.64296) then
c        print *,'64296: coefr(9-11,1) =',(coefr(ii,1),ii=9,11)
c        print *,'64296: coefr(9-11,2) =',(coefr(ii,2),ii=9,11)
c        print *,'64296: coefr(9-11,3) =',(coefr(ii,3),ii=9,11)
c        print *,'64296: coefr(9-11,4) =',(coefr(ii,4),ii=9,11)
c      endif
      do 999 igaus=1,ngaus
      call bladet(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
c      if( iblk.eq.1 .and. num.eq.64296) then
c        print *,'64296: det, nnode,nrefc,ncoor =',det,nnode,nrefc,ncoor
c        print *,'64296:coor =', coor
c        print *,'64296:coorr=', coorr
c        print *,'64296:rctr=', rctr
c        print *,'64296:crtr=', crtr
c        print *,'64296:coefr =', coefr
c      endif
c .... coordinate transfer from reference to original system
c .... rctr ---- Jacobi's matrix
c .... crtr ---- inverse matrix of Jacobi's matrix
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call eblade(refc(1,igaus),coef,coorr,coefr,coefd)
c      if( iblk.eq.1 .and. num.eq.64296) then
c        print *,'64296:igaus,coef =', igaus,(coef(ii),ii=1,14)
c        print *,'64296:refc(*,igaus) =', (refc(ii,igaus),ii=1,3)
c        print *,'64296:coorr :', coorr
c        print *,'64296:coefd :', coefd
c      endif
c .... compute coef functions and their partial derivatives
      iu=(igaus-1)*4+1
      iv=(igaus-1)*4+1
      iw=(igaus-1)*4+1
      ip=(igaus-1)*4+1
      if (num.gt.1) goto 2
c .... the following is the shape function caculation
      call blade1(refc(1,igaus),ru(1,iu),rctr,crtr)
      call blade2(refc(1,igaus),rv(1,iv),rctr,crtr)
      call blade3(refc(1,igaus),rw(1,iw),rctr,crtr)
      call blade4(refc(1,igaus),rp(1,ip),rctr,crtr)
2     continue
c .... the following is the shape function transformation
c .... from reference coordinates to original coordinates
      call shapn(nrefc,ncoor,4,ru(1,iu),cu,crtr,1,4,4)
      call shapn(nrefc,ncoor,4,rv(1,iv),cv,crtr,1,4,4)
      call shapn(nrefc,ncoor,4,rw(1,iw),cw,crtr,1,4,4)
      call shapn(nrefc,ncoor,4,rp(1,ip),cp,crtr,1,4,4)
c .... the coef function transformation
c .... from reference coordinates to original coordinates
      call shapc(nrefc,ncoor,14,coefd,coefc,crtr,2,9,9)
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
      weigh=det*gaus(igaus)
      do 100 i=1,16
      efuna(i) = 0.0
      efunb(i) = 0.0
      efunc(i) = 0.0
      efund(i) = 0.0
      efune(i) = 0.0
      efunf(i) = 0.0
100   continue
c      if( iblk.eq.1 .and. num.eq.64296) then
c        print *,'64296:igaus,umn,vmn,wmn =', igaus,umn,vmn,wmn
c      endif
      if (isolid.gt.0) then
        if (id4updown.eq.0) then
          r11=dcos(angle-angle0)
          r12=0.d0
          r13=dsin(angle-angle0)
          r21=0.d0
          r22=1.d0
          r23=0.d0
          r31=-r13
          r32=0.d0
          r33=r11
        elseif (id4updown.eq.1) then
          r11=dcos(angle-angle0)
          r12=-dsin(angle-angle0)
          r13=0.d0
          r21=-r12
          r22=r11
          r23=0.d0
          r31=0.d0
          r32=0.d0
          r33=1.d0
        elseif (id4updown.eq.-1) then
          r11=1.d0
          r12=0.d0
          r13=0.d0
          r21=0.d0
          r22=dcos(angle-angle0)
          r23=-dsin(angle-angle0)
          r31=0.d0
          r32=-r23
          r33=r22
        endif
c      if( iblk.eq.1 .and. num.eq.64296) then
c        print *,'64296: coefc(5-7,1) =',(coefc(ii,1),ii=5,7)
c        print *,'64296: coefc(9-11,1) =',(coefc(ii,1),ii=9,11)
c      endif
        uvwnx = coefc(5,1)*r11+coefc(6,1)
     * *r21+coefc(7,1)*r31
        uvwny = coefc(5,2)*r12+coefc(6,2)
     * *r22+coefc(7,2)*r32
        uvwnz = coefc(5,3)*r13+coefc(6,3)
     * *r23+coefc(7,3)*r33
        uvwnzy= coefc(5,3)*r12+coefc(6,3)
     * *r22+coefc(7,3)*r32+coefc(5,2)
     * *r13+coefc(6,2)*r23+coefc(7,2)
     * *r33
        uvwnxz= coefc(5,3)*r11+coefc(6,3)
     * *r21+coefc(7,3)*r31+coefc(5,1)
     * *r13+coefc(6,1)*r23+coefc(7,1)
     * *r33
        uvwnyx= coefc(5,2)*r11+coefc(6,2)
     * *r21+coefc(7,2)*r31+coefc(5,1)
     * *r12+coefc(6,1)*r22+coefc(7,1)
     * *r32
        uvwmx = coefc(9,1)*r11+coefc(10,1)
     * *r21+coefc(11,1)*r31
        uvwmy = coefc(9,2)*r12+coefc(10,2)
     * *r22+coefc(11,2)*r32
        uvwmz = coefc(9,3)*r13+coefc(10,3)
     * *r23+coefc(11,3)*r33
        uvwmzy= coefc(9,3)*r12+coefc(10,3)
     * *r22+coefc(11,3)*r32+coefc(9,2)
     * *r13+coefc(10,2)*r23+coefc(11,2)
     * *r33
        uvwmxz= coefc(9,3)*r11+coefc(10,3)
     * *r21+coefc(11,3)*r31+coefc(9,1)
     * *r13+coefc(10,1)*r23+coefc(11,1)
     * *r33
        uvwmyx= coefc(9,2)*r11+coefc(10,2)
     * *r21+coefc(11,2)*r31+coefc(9,1)
     * *r12+coefc(10,1)*r22+coefc(11,1)
     * *r32
        rotux=1.d0-r11
        rotvy=1.d0-r22
        rotwz=1.d0-r33
        rotuvwzy=-r32-r23
        rotuvwxz=-r31-r13
        rotuvwyx=-r21-r12
        unx = uvwnx*dt/2.0+uvwmx-rotux
        vny = uvwny*dt/2.0+uvwmy-rotvy
        wnz = uvwnz*dt/2.0+uvwmz-rotwz
        uvwzy= uvwnzy*dt/2.0+uvwmzy-rotuvwzy
        uvwxz= uvwnxz*dt/2.0+uvwmxz-rotuvwxz
        uvwyx= uvwnyx*dt/2.0+uvwmyx-rotuvwyx
      else
        r11=0.d0
        r12=0.d0
        r13=0.d0
        r21=0.d0
        r22=0.d0
        r23=0.d0
        r31=0.d0
        r32=0.d0
        r33=0.d0
        unx = 0.d0
        vny = 0.d0
        wnz = 0.d0
        uvwzy= 0.d0
        uvwxz= 0.d0
        uvwyx= 0.d0
      endif
c      if( iblk.eq.1 .and. num.eq.64296) then
c        print *,'64296: unx,vny,wnz =', unx,vny,wnz
c        print *,'64296: uvwzy,uvwxz,uvwyz=',uvwzy,uvwxz,uvwyx
c      endif
      do 101 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,2)*r11
      efuna(iv)=efuna(iv)+stif
101   continue
      do 102 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,2)*r21
      efuna(iv)=efuna(iv)+stif
102   continue
      do 103 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,2)*r31
      efuna(iv)=efuna(iv)+stif
103   continue
      do 104 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,3)*r12
      efunb(iv)=efunb(iv)+stif
104   continue
      do 105 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,3)*r22
      efunb(iv)=efunb(iv)+stif
105   continue
      do 106 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,3)*r32
      efunb(iv)=efunb(iv)+stif
106   continue
      do 107 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,4)*r13
      efunc(iv)=efunc(iv)+stif
107   continue
      do 108 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,4)*r23
      efunc(iv)=efunc(iv)+stif
108   continue
      do 109 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,4)*r33
      efunc(iv)=efunc(iv)+stif
109   continue
      do 110 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,4)*r12
     & +cu(i,3)*r13
      efund(iv)=efund(iv)+stif
110   continue
      do 111 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,4)*r22
     & +cv(i,3)*r23
      efund(iv)=efund(iv)+stif
111   continue
      do 112 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,4)*r32
     & +cw(i,3)*r33
      efund(iv)=efund(iv)+stif
112   continue
      do 113 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,4)*r11
     & +cu(i,2)*r13
      efune(iv)=efune(iv)+stif
113   continue
      do 114 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,4)*r21
     & +cv(i,2)*r23
      efune(iv)=efune(iv)+stif
114   continue
      do 115 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,4)*r31
     & +cw(i,2)*r33
      efune(iv)=efune(iv)+stif
115   continue
      do 116 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,3)*r11
     & +cu(i,2)*r12
      efunf(iv)=efunf(iv)+stif
116   continue
      do 117 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,3)*r21
     & +cv(i,2)*r22
      efunf(iv)=efunf(iv)+stif
117   continue
      do 118 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,3)*r31
     & +cw(i,2)*r32
      efunf(iv)=efunf(iv)+stif
118   continue
c .... the following is the stiffness computation
      do 202 i=1,4
      iv=kvord(i,1)
      do 201 j=1,4
      jv=kvord(j,1)
      stif=+cu(i,2)*cu(j,2)*visc*ifluid
     & +cu(i,3)*cu(j,3)*visc*ifluid
     & +cu(i,4)*cu(j,4)*visc*ifluid
     & +cu(i,2)*cu(j,2)*visc*ifluid
     & +cu(i,2)*cu(j,1)*(un-vlx)*denf*ifluid
     & +cu(i,3)*cu(j,1)*(vn-vly)*denf*ifluid
     & +cu(i,4)*cu(j,1)*(wn-vlz)*denf*ifluid
     & +cu(i,1)*cu(j,1)*coefc(1,1)*denf*ifluid
     | +cu(i,2)*cu(j,2)*(un-vlx)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
     & +cu(i,3)*cu(j,3)*(vn-vly)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
     & +cu(i,4)*cu(j,4)*(wn-vlz)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
201    continue
202    continue
      do 204 i=1,4
      iv=kvord(i,1)
      do 203 j=1,4
      jv=kvord(j,2)
      stif=+cu(i,3)*cv(j,2)*visc*ifluid
     & +cu(i,1)*cv(j,1)*coefc(2,1)*denf*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
203    continue
204    continue
      do 206 i=1,4
      iv=kvord(i,1)
      do 205 j=1,4
      jv=kvord(j,3)
      stif=+cu(i,4)*cw(j,2)*visc*ifluid
     & +cu(i,1)*cw(j,1)*coefc(3,1)*denf*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
205    continue
206    continue
      do 208 i=1,4
      iv=kvord(i,1)
      do 207 j=1,4
      jv=kvord(j,4)
      stif=-cu(i,2)*cp(j,1)*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
207    continue
208    continue
      do 210 i=1,4
      iv=kvord(i,2)
      do 209 j=1,4
      jv=kvord(j,1)
      stif=+cv(i,2)*cu(j,3)*visc*ifluid
     & +cv(i,1)*cu(j,1)*coefc(1,2)*denf*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
209    continue
210    continue
      do 212 i=1,4
      iv=kvord(i,2)
      do 211 j=1,4
      jv=kvord(j,2)
      stif=+cv(i,2)*cv(j,2)*visc*ifluid
     & +cv(i,3)*cv(j,3)*visc*ifluid
     & +cv(i,4)*cv(j,4)*visc*ifluid
     & +cv(i,3)*cv(j,3)*visc*ifluid
     & +cv(i,2)*cv(j,1)*(un-vlx)*denf*ifluid
     & +cv(i,3)*cv(j,1)*(vn-vly)*denf*ifluid
     & +cv(i,4)*cv(j,1)*(wn-vlz)*denf*ifluid
     & +cv(i,1)*cv(j,1)*coefc(2,2)*denf*ifluid
     | +cv(i,2)*cv(j,2)*(un-vlx)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
     & +cv(i,3)*cv(j,3)*(vn-vly)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
     & +cv(i,4)*cv(j,4)*(wn-vlz)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
211    continue
212    continue
      do 214 i=1,4
      iv=kvord(i,2)
      do 213 j=1,4
      jv=kvord(j,3)
      stif=+cv(i,4)*cw(j,3)*visc*ifluid
     & +cv(i,1)*cw(j,1)*coefc(3,2)*denf*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
213    continue
214    continue
      do 216 i=1,4
      iv=kvord(i,2)
      do 215 j=1,4
      jv=kvord(j,4)
      stif=-cv(i,3)*cp(j,1)*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
215    continue
216    continue
      do 218 i=1,4
      iv=kvord(i,3)
      do 217 j=1,4
      jv=kvord(j,1)
      stif=+cw(i,2)*cu(j,4)*visc*ifluid
     & +cw(i,1)*cu(j,1)*coefc(1,3)*denf*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
217    continue
218    continue
      do 220 i=1,4
      iv=kvord(i,3)
      do 219 j=1,4
      jv=kvord(j,2)
      stif=+cw(i,3)*cv(j,4)*visc*ifluid
     & +cw(i,1)*cv(j,1)*coefc(2,3)*denf*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
219    continue
220    continue
      do 222 i=1,4
      iv=kvord(i,3)
      do 221 j=1,4
      jv=kvord(j,3)
      stif=+cw(i,2)*cw(j,2)*visc*ifluid
     & +cw(i,3)*cw(j,3)*visc*ifluid
     & +cw(i,4)*cw(j,4)*visc*ifluid
     & +cw(i,4)*cw(j,4)*visc*ifluid
     & +cw(i,2)*cw(j,1)*(un-vlx)*denf*ifluid
     & +cw(i,3)*cw(j,1)*(vn-vly)*denf*ifluid
     & +cw(i,4)*cw(j,1)*(wn-vlz)*denf*ifluid
     & +cw(i,1)*cw(j,1)*coefc(3,3)*denf*ifluid
     | +cw(i,2)*cw(j,2)*(un-vlx)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
     & +cw(i,3)*cw(j,3)*(vn-vly)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
     & +cw(i,4)*cw(j,4)*(wn-vlz)**2*denf**2*det**(2.d0/
     | 3.d0)/visc*pen4sd*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
221    continue
222    continue
      do 224 i=1,4
      iv=kvord(i,3)
      do 223 j=1,4
      jv=kvord(j,4)
      stif=-cw(i,4)*cp(j,1)*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
223    continue
224    continue
      do 226 i=1,4
      iv=kvord(i,4)
      do 225 j=1,4
      jv=kvord(j,1)
      stif=-cp(i,1)*cu(j,2)*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
225    continue
226    continue
      do 228 i=1,4
      iv=kvord(i,4)
      do 227 j=1,4
      jv=kvord(j,2)
      stif=-cp(i,1)*cv(j,3)*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
227    continue
228    continue
      do 230 i=1,4
      iv=kvord(i,4)
      do 229 j=1,4
      jv=kvord(j,3)
      stif=-cp(i,1)*cw(j,4)*ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
229    continue
230    continue
      do 232 i=1,4
      iv=kvord(i,4)
      do 231 j=1,4
      jv=kvord(j,4)
      stif=-cp(i,2)*cp(j,2)*det**(2.d0/3.d0)/visc*penalty*
     | ifluid
     & -cp(i,3)*cp(j,3)*det**(2.d0/3.d0)/visc*penalty*
     | ifluid
     & -cp(i,4)*cp(j,4)*det**(2.d0/3.d0)/visc*penalty*
     | ifluid
      estif(iv,jv)=estif(iv,jv)+stif*weigh
231    continue
232    continue
      do 234 iv=1,16
      do 233 jv=1,16
      stif=+efuna(iv)*efuna(jv)*fact*(1.d0-pv)*dt/2.d0
     & +efuna(iv)*efunb(jv)*fact*(pv)*dt/2.d0
     & +efuna(iv)*efunc(jv)*fact*(pv)*dt/2.d0
     & +efunb(iv)*efuna(jv)*fact*(pv)*dt/2.d0
     & +efunb(iv)*efunb(jv)*fact*(1.d0-pv)*dt/2.d0
     & +efunb(iv)*efunc(jv)*fact*(pv)*dt/2.d0
     & +efunc(iv)*efuna(jv)*fact*(pv)*dt/2.d0
     & +efunc(iv)*efunb(jv)*fact*(pv)*dt/2.d0
     & +efunc(iv)*efunc(jv)*fact*(1.d0-pv)*dt/2.d0
     & +efund(iv)*efund(jv)*fact*(0.5d0-pv)*dt/2.d0
     & +efune(iv)*efune(jv)*fact*(0.5d0-pv)*dt/2.d0
     & +efunf(iv)*efunf(jv)*fact*(0.5d0-pv)*dt/2.d0
      estif(iv,jv)=estif(iv,jv)+stif*weigh
233    continue
234    continue
c .... the following is the mass matrix computation
      stif=den
      elump(1)=stif*weigh
      stif=den
      elump(5)=stif*weigh
      stif=den
      elump(9)=stif*weigh
      stif=den
      elump(13)=stif*weigh
      stif=den
      elump(2)=stif*weigh
      stif=den
      elump(6)=stif*weigh
      stif=den
      elump(10)=stif*weigh
      stif=den
      elump(14)=stif*weigh
      stif=den
      elump(3)=stif*weigh
      stif=den
      elump(7)=stif*weigh
      stif=den
      elump(11)=stif*weigh
      stif=den
      elump(15)=stif*weigh
      stif=0.d0
      elump(4)=stif*weigh
      stif=0.d0
      elump(8)=stif*weigh
      stif=0.d0
      elump(12)=stif*weigh
      stif=0.d0
      elump(16)=stif*weigh
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
c .... the following is the damping matrix computation
      stif=0.d0
      elump(1)=stif*weigh
      stif=0.d0
      elump(5)=stif*weigh
      stif=0.d0
      elump(9)=stif*weigh
      stif=0.d0
      elump(13)=stif*weigh
      stif=0.d0
      elump(2)=stif*weigh
      stif=0.d0
      elump(6)=stif*weigh
      stif=0.d0
      elump(10)=stif*weigh
      stif=0.d0
      elump(14)=stif*weigh
      stif=0.d0
      elump(3)=stif*weigh
      stif=0.d0
      elump(7)=stif*weigh
      stif=0.d0
      elump(11)=stif*weigh
      stif=0.d0
      elump(15)=stif*weigh
      stif=ifluid
      elump(4)=stif*weigh
      stif=ifluid
      elump(8)=stif*weigh
      stif=ifluid
      elump(12)=stif*weigh
      stif=ifluid
      elump(16)=stif*weigh
      do 401 i=1,nvard(1)
      iv = kvord(i,1)
      edamp(iv)=edamp(iv)+elump(iv)*cu(i,1)
401   continue
      do 402 i=1,nvard(2)
      iv = kvord(i,2)
      edamp(iv)=edamp(iv)+elump(iv)*cv(i,1)
402   continue
      do 403 i=1,nvard(3)
      iv = kvord(i,3)
      edamp(iv)=edamp(iv)+elump(iv)*cw(i,1)
403   continue
      do 404 i=1,nvard(4)
      iv = kvord(i,4)
      edamp(iv)=edamp(iv)+elump(iv)*cp(i,1)
404   continue
c .... the following is the load vector computation
      do 501 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,1)*(un*coefc(1,1)+vn*coefc(1,2)+
     ! wn*coefc(1,3))*denf*ifluid
     & +cu(i,1)*fx
      eload(iv)=eload(iv)+stif*weigh
c      if(iblk.eq.1 .and. num.eq.64296 .and. iv.eq.1) then
c        print *,'1--igaus,iv,stif,weigh,eload=',iguas,iv,stif,weigh,eload(iv)
c      endif
501   continue
      do 502 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,1)*(un*coefc(2,1)+vn*coefc(2,2)+
     ! wn*coefc(2,3))*denf*ifluid
     & +cv(i,1)*fy
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,1)*(un*coefc(3,1)+vn*coefc(3,2)+
     ! wn*coefc(3,3))*denf*ifluid
     & +cw(i,1)*fz
      eload(iv)=eload(iv)+stif*weigh
503   continue
      do 504 iv=1,16
      stif=-efuna(iv)*unx*fact*(1.d0-pv)
     & -efunb(iv)*unx*fact*(pv)
     & -efunc(iv)*unx*fact*(pv)
     & -efuna(iv)*vny*fact*(pv)
     & -efunb(iv)*vny*fact*(1.d0-pv)
     & -efunc(iv)*vny*fact*(pv)
     & -efuna(iv)*wnz*fact*(pv)
     & -efunb(iv)*wnz*fact*(pv)
     & -efunc(iv)*wnz*fact*(1.d0-pv)
     & -efund(iv)*uvwzy*fact*(0.5d0-pv)
     & -efune(iv)*uvwxz*fact*(0.5d0-pv)
     & -efunf(iv)*uvwyx*fact*(0.5d0-pv)
      eload(iv)=eload(iv)+stif*weigh
c      if(iblk.eq.1 .and. num.eq.64296 .and. iv.eq.1) then
c        print *,'2--igaus,iv,stif,weigh,eload=',iguas,iv,stif,weigh,eload(iv)
c        print *,'efuna,efunb,efunc =', efuna(iv), efunb(iv),efunc(iv)
c        print *,'unx,vny,wnz =', unx,vny,wnz
c        print *,'efund,efune,efunf =', efund(iv),efune(iv),efunf(iv)
c        print *,'unwzy,uvwxz,uvwyx=',unwzy,uvwxz,uvwyx
c        print *,'fact,pv =', fact, pv
c      endif
504   continue
999   continue
998   continue
      return
      end

      subroutine bladei
      implicit real*8 (a-h,o-z)
      common /dblade/ refc(3,4),gaus(4),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(4),kdord(4),kvord(16,4)
c .... initial data
c .... refc ---- reference coordinates at integral points
c .... gaus ---- weight number at integral points
c .... nvard ---- the number of var for each unknown
c .... kdord ---- the highest differential order for each unknown
c .... kvord ---- var number at integral points for each unknown
      ngaus=  4
      ndisp=  4
      nrefc=  3
      ncoor=  3
      nvar = 16
      nnode=  4
      kdord(1)=1
      nvard(1)=4
      kvord(1,1)=1
      kvord(2,1)=5
      kvord(3,1)=9
      kvord(4,1)=13
      kdord(2)=1
      nvard(2)=4
      kvord(1,2)=2
      kvord(2,2)=6
      kvord(3,2)=10
      kvord(4,2)=14
      kdord(3)=1
      nvard(3)=4
      kvord(1,3)=3
      kvord(2,3)=7
      kvord(3,3)=11
      kvord(4,3)=15
      kdord(4)=1
      nvard(4)=4
      kvord(1,4)=4
      kvord(2,4)=8
      kvord(3,4)=12
      kvord(4,4)=16
      refc(1,1)=0.5854102D0
      refc(2,1)=0.1381966D0
      refc(3,1)=0.1381966D0
      gaus(1)=0.25D0
      refc(1,2)=0.1381966D0
      refc(2,2)=0.5854102D0
      refc(3,2)=0.1381966D0
      gaus(2)=0.25D0
      refc(1,3)=0.1381966D0
      refc(2,3)=0.1381966D0
      refc(3,3)=0.5854102D0
      gaus(3)=0.25D0
      refc(1,4)=0.1381966D0
      refc(2,4)=0.1381966D0
      refc(3,4)=0.1381966D0
      gaus(4)=0.25D0
      end


      subroutine bladet(nnode,nrefc,ncoor,refc,coor,coorr,
     & rc,cr,det,coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     & coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call tblade(refc,coor,coorr,coefr,rc)
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
 

      subroutine blade1(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external fblade1
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fblade1,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fblade1(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccblade/ xa(4),ya(4),za(4),una(4),
     &vna(4),wna(4),pna(4),unna(4),vnna(4),wnna(4),
     &pnna(4),umna(4),vmna(4),wmna(4),vlxa(4),vlya(4),
     &vlza(4)
      common /vblade/ rctr(3,3),crtr(3,3),coefd(14,9),coefc(14,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     fblade1=+rx 
      goto 1000
2     fblade1=+ry 
      goto 1000
3     fblade1=+rz 
      goto 1000
4     fblade1=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine blade2(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external fblade2
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fblade2,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fblade2(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccblade/ xa(4),ya(4),za(4),una(4),
     &vna(4),wna(4),pna(4),unna(4),vnna(4),wnna(4),
     &pnna(4),umna(4),vmna(4),wmna(4),vlxa(4),vlya(4),
     &vlza(4)
      common /vblade/ rctr(3,3),crtr(3,3),coefd(14,9),coefc(14,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     fblade2=+rx 
      goto 1000
2     fblade2=+ry 
      goto 1000
3     fblade2=+rz 
      goto 1000
4     fblade2=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine blade3(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external fblade3
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fblade3,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fblade3(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccblade/ xa(4),ya(4),za(4),una(4),
     &vna(4),wna(4),pna(4),unna(4),vnna(4),wnna(4),
     &pnna(4),umna(4),vmna(4),wmna(4),vlxa(4),vlya(4),
     &vlza(4)
      common /vblade/ rctr(3,3),crtr(3,3),coefd(14,9),coefc(14,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     fblade3=+rx 
      goto 1000
2     fblade3=+ry 
      goto 1000
3     fblade3=+rz 
      goto 1000
4     fblade3=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine blade4(refc,shpr,rctr,crtr)
c .... compute shape functions and their partial derivatives
c .... shapr ---- store shape functions and their partial derivatives
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(4,4),rctr(3,3),crtr(3,3)
      external fblade4
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fblade4,refc,shpr,3,4,1)
c .... shape function and their derivatives computation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function fblade4(refc,n)
c .... shape function caculation
      implicit real*8 (a-h,o-z)
      common /ccblade/ xa(4),ya(4),za(4),una(4),
     &vna(4),wna(4),pna(4),unna(4),vnna(4),wnna(4),
     &pnna(4),umna(4),vmna(4),wmna(4),vlxa(4),vlya(4),
     &vlza(4)
      common /vblade/ rctr(3,3),crtr(3,3),coefd(14,9),coefc(14,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     fblade4=+rx 
      goto 1000
2     fblade4=+ry 
      goto 1000
3     fblade4=+rz 
      goto 1000
4     fblade4=+(+1.-rx-ry-rz) 
      goto 1000
1000  return
      end

      subroutine eblade(refc,coef,coorr,coefr,coefd)
c .... compute coef value and their partial derivatives
c .... by reference coordinate value
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(14),coorr(3,4),coefr(4,14),coefd(14,3)
      external feblade
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(feblade,refc,coef,coefd,3,14,2)
c .... coef value and their partial derivatives caculation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function feblade(refc,n)
c .... coef function caculation
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccblade/ xa(4),ya(4),za(4),un(4),vn(4),wn(4),
     &pn(4),unn(4),vnn(4),wnn(4),pnn(4),umn(4),
     &vmn(4),wmn(4),vlx(4),vly(4),vlz(4)
      common /vblade/ rctr(3,3),crtr(3,3),coefd(14,9),coefc(14,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14) n
1     feblade=+(+rx)*un(1)+(+ry)*un(2)+(+rz)*un(3)+(+(+1.-
     & rx-ry-rz))*un(4)
      goto 1000
2     feblade=+(+rx)*vn(1)+(+ry)*vn(2)+(+rz)*vn(3)+(+(+1.-
     & rx-ry-rz))*vn(4)
      goto 1000
3     feblade=+(+rx)*wn(1)+(+ry)*wn(2)+(+rz)*wn(3)+(+(+1.-
     & rx-ry-rz))*wn(4)
      goto 1000
4     feblade=+(+rx)*pn(1)+(+ry)*pn(2)+(+rz)*pn(3)+(+(+1.-
     & rx-ry-rz))*pn(4)
      goto 1000
5     feblade=+(+rx)*unn(1)+(+ry)*unn(2)+(+rz)*unn(3)+(+(+
     & 1.-rx-ry-rz))*unn(4)
      goto 1000
6     feblade=+(+rx)*vnn(1)+(+ry)*vnn(2)+(+rz)*vnn(3)+(+(+
     & 1.-rx-ry-rz))*vnn(4)
      goto 1000
7     feblade=+(+rx)*wnn(1)+(+ry)*wnn(2)+(+rz)*wnn(3)+(+(+
     & 1.-rx-ry-rz))*wnn(4)
      goto 1000
8     feblade=+(+rx)*pnn(1)+(+ry)*pnn(2)+(+rz)*pnn(3)+(+(+
     & 1.-rx-ry-rz))*pnn(4)
      goto 1000
9     feblade=+(+rx)*umn(1)+(+ry)*umn(2)+(+rz)*umn(3)+(+(+
     & 1.-rx-ry-rz))*umn(4)
      goto 1000
10     feblade=+(+rx)*vmn(1)+(+ry)*vmn(2)+(+rz)*vmn(3)+(+(+
     & 1.-rx-ry-rz))*vmn(4)
      goto 1000
11     feblade=+(+rx)*wmn(1)+(+ry)*wmn(2)+(+rz)*wmn(3)+(+(+
     & 1.-rx-ry-rz))*wmn(4)
      goto 1000
12     feblade=+(+rx)*vlx(1)+(+ry)*vlx(2)+(+rz)*vlx(3)+(+(+
     & 1.-rx-ry-rz))*vlx(4)
      goto 1000
13     feblade=+(+rx)*vly(1)+(+ry)*vly(2)+(+rz)*vly(3)+(+(+
     & 1.-rx-ry-rz))*vly(4)
      goto 1000
14     feblade=+(+rx)*vlz(1)+(+ry)*vlz(2)+(+rz)*vlz(3)+(+(+
     & 1.-rx-ry-rz))*vlz(4)
      goto 1000
1000  return
      end

      subroutine tblade(refc,coor,coorr,coefr,rc)
c .... compute coordinate value and Jacobi's matrix rc
c .... by reference coordinate value
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,4),coefr(4,14),rc(3,3)
      common /ccblade/ x(4),y(4),z(4),un(4),vn(4),wn(4),
     &pn(4),unn(4),vnn(4),wnn(4),pnn(4),umn(4),
     &vmn(4),wmn(4),vlx(4),vly(4),vlz(4)
      external ftblade
      do 100 n=1,4
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,4
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
200   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(ftblade,refc,coor,rc,3,3,1)
c .... coordinate value and their partial derivatives caculation
c .... compute partial derivatives by centered difference
c .... which is in the file ccshap.for of FEPG library
      return
      end

      real*8 function ftblade(refc,n)
c .... coordinate transfer function caculation
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccblade/ x(4),y(4),z(4),un(4),vn(4),wn(4),
     &pn(4),unn(4),vnn(4),wnn(4),pnn(4),umn(4),
     &vmn(4),wmn(4),vlx(4),vly(4),vlz(4)
      common /vblade/ rctr(3,3),crtr(3,3),coefd(14,9),coefc(14,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     ftblade=+(+rx)*x(1)+(+ry)*x(2)+(+rz)*x(3)+(+(+1.-rx-
     & ry-rz))*x(4)
      goto 1000
2     ftblade=+(+rx)*y(1)+(+ry)*y(2)+(+rz)*y(3)+(+(+1.-rx-
     & ry-rz))*y(4)
      goto 1000
3     ftblade=+(+rx)*z(1)+(+ry)*z(2)+(+rz)*z(3)+(+(+1.-rx-
     & ry-rz))*z(4)
      goto 1000
1000  return
      end

