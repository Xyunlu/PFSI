      subroutine cpunod(u0,u1,knode,kdgof,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      DIMENSION U0(knode,kdgof),U1(knode,kdgof)
 
      do i=1,knode
        do j=1,kdgof
          U1(i,j)=U0(i,j)
        enddo
      enddo
 
      RETURN
      END
