      subroutine sgetpart(maplg,imaplg,id,
     & node,emate,disp0,u0,
     & u1,u2,u3,tttt,
     & knode,kdgof,kelem,kemate,initno,numetyp,iblk,nsource,
     & neq,n_update,n_external,n_dataorg,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      dimension id(knode,kdgof),
     & disp0(knode,kdgof),u0(knode,kdgof),
     & u1(knode,kdgof),u2(knode,kdgof),
     & u3(knode,kdgof),
     & node(kelem),emate(kemate),
     & maplg(neq),imaplg(neq),
     & tttt(1)
 
      dimension mnode(100),nnode(100),mmate(100),nmate(100)
c.......   recieve lgnv file
      call recvai(iblk,nsource,maplg,neq)
c.......   recieve ilgnv file
      call recvai(iblk,nsource,imaplg,neq)
c......    recieve n_update constant and then compute some constants
      call recvint(iblk,nsource,n_update)
      n_external = neq-n_update
      n_dataorg = n_update +762
 
c.......   recieve id0 file
      call recvai(iblk,nsource,id,knode*kdgof)
 
c ............... element
      kelem = 1
      kemate = 1
      do ityp=1,numetyp
        call recvint(iblk,nsource,mnode(ityp))
        call recvint(iblk,nsource,nnode(ityp))
        call recvai(iblk,nsource,node(kelem),mnode(ityp)*nnode(ityp))
        kelem = kelem+nnode(ityp)*mnode(ityp)
        call recvint(iblk,nsource,mmate(ityp))
        call recvint(iblk,nsource,nmate(ityp))
        call recvar(iblk,nsource,emate(kemate),
     &   mmate(ityp)*nmate(ityp))
        kemate = kemate+nmate(ityp)*mmate(ityp)
      enddo
      kelem = kelem-1
      kemate = kemate-1
      print *,'knode,kelem and kemate are',knode,kelem,kemate,' at',iblk
 
c.......   recieve disp0 file
      call recvar(iblk,nsource,disp0,knode*kdgof)
 
c......    recieve initial value file disp1,disp2,disp3
      do i=1,kdgof
        do j=1,knode
          u0(j,i)=disp0(j,i)
        enddo
      enddo
      do i=1,kdgof
        do j=1,knode
          u1(j,i)=0.0d0
          u2(j,i)=0.0d0
          u3(j,i)=0.0d0
        enddo
      enddo
      if(initno.gt.0) then
        call recvar(iblk,nsource,u0,knode*kdgof)
      endif
      if(initno.gt.1) then
        call recvar(iblk,nsource,u1,knode*kdgof)
      endif
      if(initno.eq.3) then
        call recvar(iblk,nsource,u2,knode*kdgof)
      endif
 
      end
 
