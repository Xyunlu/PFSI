      subroutine cpunodd(um,vm,wm,
     & umn,vmn,wmn,vlx,
     & vly,vlz,vlxn,vlyn,
     & vlzn,um1,vm1,wm1,
     & umn1,vmn1,wmn1,vlx1,
     & vly1,vlz1,vlxn1,vlyn1,
     & vlzn1,knode,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      DIMENSION um(knode),vm(knode),wm(knode),
     +          umn(knode),vmn(knode),wmn(knode),
     +          vlx(knode),vly(knode),vlz(knode),
     +          vlxn(knode),vlyn(knode),vlzn(knode)
      DIMENSION um1(knode),vm1(knode),wm1(knode),
     +          umn1(knode),vmn1(knode),wmn1(knode),
     +          vlx1(knode),vly1(knode),vlz1(knode),
     +          vlxn1(knode),vlyn1(knode),vlzn1(knode)
 
      do i=1,knode
          Um1(i)=Um(i)
          vm1(i)=vm(i)
          wm1(i)=wm(i)
          Umn1(i)=Umn(i)
          vmn1(i)=vmn(i)
          wmn1(i)=wmn(i)
          vlx1(i)=vlx(i)
          vly1(i)=vly(i)
          vlz1(i)=vlz(i)
          vlxn1(i)=vlxn(i)
          vlyn1(i)=vlyn(i)
          vlzn1(i)=vlzn(i)
      enddo
 
      RETURN
      END
