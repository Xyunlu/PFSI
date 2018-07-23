      subroutine spart(coor,iorder,tttt,
     & nsource,iblk,knode,kcoor,
     & mnodea,nnodea,mmatea,nmatea,
     & mnodeb,nnodeb,mmateb,nmateb,
     & mnodec,nnodec,mmatec,nmatec)
      implicit real*8 (a-h,o-z)
      dimension coor(knode,kcoor),iorder(knode),tttt(1)
 
      call recvar(iblk,nsource,coor,knode*kcoor)
      call recvai(iblk,nsource,iorder,knode)
 
      return
      end
 
