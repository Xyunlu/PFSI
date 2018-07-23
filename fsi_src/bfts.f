      subroutine bfts(bf,iblk,nsource,tmax,dt,
     & time,it,msstop,knode,kdgof,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      dimension bf(knode,kdgof)
 
7     format(5(1x,e15.5))
C.... recv msstop
      call recvint(iblk,nsource,msstop)
C.... read time file
      call recvr(iblk,nsource,tmax)
      call recvr(iblk,nsource,dt)
      call recvr(iblk,nsource,time)
      call recvint(iblk,nsource,it)
 
      do i=1,knode
         do j=1,kdgof
            bf(i,j)=0.
         enddo
      enddo
      call recvar(iblk,nsource,bf,knode*kdgof)
c      write(*,*) 'iblk,bf =====',iblk
c      write(*,7) bf
 
      end
 
