      subroutine cpcoor0(u0,u1,knode,kcoor,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      DIMENSION U0(knode,kcoor),U1(knode,kcoor)
 
      do i=1,knode
        do j=1,kcoor
          U1(i,j)=U0(i,j)
        enddo
      enddo
 
      RETURN
      END
