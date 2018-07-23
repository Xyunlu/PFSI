      subroutine slavep(iblk)
      implicit real*8 (a-h,o-z)
      common /fsia/ kdgofa,kelema,kematea,numetypa,
     & numelema,neqa,kvara,
     & mnodea(100),nnodea(100),mmatea(100),nmatea(100)
      common /fsib/ kdgofb,kelemb,kemateb,numetypb,
     & numelemb,neqb,kvarb,
     & mnodeb(100),nnodeb(100),mmateb(100),nmateb(100)
      common /fsic/ kdgofc,kelemc,kematec,numetypc,
     & numelemc,neqc,kvarc,
     & mnodec(100),nnodec(100),mmatec(100),nmatec(100)
      character*1 material
      common /ia/ ia(150000000)
      common /ib/ ib(75000000)
      nsource = 0
      material = 'y'
      itn = 1
      cc = 1.0
      maxa = 75000000
      call recvint(iblk,nsource,nblk)
      call recvr(iblk,nsource,t0)
      call recvr(iblk,nsource,tmax)
      call recvr(iblk,nsource,dt)
      call recvint(iblk,nsource,knode)
      call recvint(iblk,nsource,kcoor)
      kf0 = 1
c ........ file coor0 ........
c ........ array coor ........
      kp1f1a1 = knode*kcoor
      kp1f1a1 = kp1f1a1*2
      kf1 = +kp1f1a1
      kf1 = kf1+kf0
      if (kf1.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf1 memory needed = ',kf1
      stop 55555
      endif
      kp1f1a1 = kf1-kp1f1a1
c ........ file order0 ........
c ........ array iorder ........
      kp1f2a1 = knode
      if (kp1f2a1/2*2 .lt. kp1f2a1) kp1f2a1 = kp1f2a1+1
      kf2 = +kp1f2a1
      kf2 = kf2+kf1
      if (kf2.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf2 memory needed = ',kf2
      stop 55555
      endif
      kp1f2a1 = kf2-kp1f2a1
c ............... temporary array  ...............
c ........ array tttt ........
      kp1b1 = 1
      kp1b1 = kp1b1*2
      kp1b0 = 1 
      kp1b1 = kp1b1+kp1b0
      call spart(ia(kp1f1a1),ia(kp1f2a1),ib(kp1b0),
     & nsource,iblk,knode,kcoor,
     & mnodea,nnodea,mmatea,nmatea,
     & mnodeb,nnodeb,mmateb,nmateb,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
c ...........................................
      maxaa = maxa
      numetypa = 2
      call recvint(iblk,nsource,kelema)
      call recvint(iblk,nsource,kematea)
      kdgofa = 4
      kvara = knode*kdgofa+1
      call recvint(iblk,nsource,neqa)
      initnoa = 0
c ........ file lgnva ........
c ........ array maplg ........
      kp2f1a1 = neqa
      if (kp2f1a1/2*2 .lt. kp2f1a1) kp2f1a1 = kp2f1a1+1
      kf3 = +kp2f1a1
      kf3 = kf3+kf2
      if (kf3.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf3 memory needed = ',kf3
      stop 55555
      endif
      kp2f1a1 = kf3-kp2f1a1
c ........ file idxlgnva ........
c ........ array imaplg ........
      kp2f2a1 = neqa
      if (kp2f2a1/2*2 .lt. kp2f2a1) kp2f2a1 = kp2f2a1+1
      kf4 = +kp2f2a1
      kf4 = kf4+kf3
      if (kf4.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf4 memory needed = ',kf4
      stop 55555
      endif
      kp2f2a1 = kf4-kp2f2a1
c ........ file ida0 ........
c ........ array id ........
      kp2f3a1 = knode*kdgofa
      if (kp2f3a1/2*2 .lt. kp2f3a1) kp2f3a1 = kp2f3a1+1
      kf5 = +kp2f3a1
      kf5 = kf5+kf4
      if (kf5.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf5 memory needed = ',kf5
      stop 55555
      endif
      kp2f3a1 = kf5-kp2f3a1
c ........ file elema0 ........
c ........ array node ........
      kp2f4a1 = kelema
      if (kp2f4a1/2*2 .lt. kp2f4a1) kp2f4a1 = kp2f4a1+1
c ........ array emate ........
      kp2f4a2 = kematea
      kp2f4a2 = kp2f4a2*2
      kf6 = +kp2f4a1+kp2f4a2
      kf6 = kf6+kf5
      if (kf6.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf6 memory needed = ',kf6
      stop 55555
      endif
      kp2f4a2 = kf6-kp2f4a2
      kp2f4a1 = kp2f4a2-kp2f4a1
c ........ file dispa0 ........
c ........ array disp0 ........
      kp2f5a1 = knode*kdgofa
      kp2f5a1 = kp2f5a1*2
      kf7 = +kp2f5a1
      kf7 = kf7+kf6
      if (kf7.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf7 memory needed = ',kf7
      stop 55555
      endif
      kp2f5a1 = kf7-kp2f5a1
c ........ file unoda ........
c ........ array u0 ........
      kp2f6a1 = knode*kdgofa
      kp2f6a1 = kp2f6a1*2
c ........ array u1 ........
      kp2f6a2 = knode*kdgofa
      kp2f6a2 = kp2f6a2*2
c ........ array u2 ........
      kp2f6a3 = knode*kdgofa
      kp2f6a3 = kp2f6a3*2
c ........ array u3 ........
      kp2f6a4 = knode*kdgofa
      kp2f6a4 = kp2f6a4*2
      kf8 = +kp2f6a1+kp2f6a2+kp2f6a3+kp2f6a4
      kf8 = kf8+kf7
      if (kf8.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf8 memory needed = ',kf8
      stop 55555
      endif
      kp2f6a4 = kf8-kp2f6a4
      kp2f6a3 = kp2f6a4-kp2f6a3
      kp2f6a2 = kp2f6a3-kp2f6a2
      kp2f6a1 = kp2f6a2-kp2f6a1
c ............... temporary array  ...............
c ........ array tttt ........
      kp2b1 = 1
      kp2b1 = kp2b1*2
      kp2b0 = 1 
      kp2b1 = kp2b1+kp2b0
      call sgetpart(ia(kp2f1a1),ia(kp2f2a1),ia(kp2f3a1),
     & ia(kp2f4a1),ia(kp2f4a2),ia(kp2f5a1),ia(kp2f6a1),
     & ia(kp2f6a2),ia(kp2f6a3),ia(kp2f6a4),ib(kp2b0),
     & knode,kdgofa,kelema,kematea,initnoa,numetypa,iblk,nsource,
     & neqa,n_updatea,n_externala,n_dataorga,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
      call recvint(iblk,nsource,nblk)
      call recvr(iblk,nsource,t0)
      call recvr(iblk,nsource,tmax)
      call recvr(iblk,nsource,dt)
      call recvint(iblk,nsource,knode1)
      call recvint(iblk,nsource,kcoor)
c ........ file coor1 ........
c ........ array coor ........
      kp3f1a1 = knode1*kcoor
      kp3f1a1 = kp3f1a1*2
      kf9 = +kp3f1a1
      kf9 = kf9+kf8
      if (kf9.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf9 memory needed = ',kf9
      stop 55555
      endif
      kp3f1a1 = kf9-kp3f1a1
c ........ file order1 ........
c ........ array iorder ........
      kp3f2a1 = knode1
      if (kp3f2a1/2*2 .lt. kp3f2a1) kp3f2a1 = kp3f2a1+1
      kf10 = +kp3f2a1
      kf10 = kf10+kf9
      if (kf10.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf10 memory needed = ',kf10
      stop 55555
      endif
      kp3f2a1 = kf10-kp3f2a1
c ............... temporary array  ...............
c ........ array tttt ........
      kp3b1 = 1
      kp3b1 = kp3b1*2
      kp3b0 = 1 
      kp3b1 = kp3b1+kp3b0
      call spart(ia(kp3f1a1),ia(kp3f2a1),ib(kp3b0),
     & nsource,iblk,knode1,kcoor,
     & mnodea,nnodea,mmatea,nmatea,
     & mnodeb,nnodeb,mmateb,nmateb,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
c ...........................................
      maxab = maxa
      numetypb = 1
      call recvint(iblk,nsource,kelemb)
      call recvint(iblk,nsource,kemateb)
      kdgofb = 3
      kvarb = knode1*kdgofb+1
      call recvint(iblk,nsource,neqb)
      initnob = 0
c ........ file lgnvb ........
c ........ array maplg ........
      kp4f1a1 = neqb
      if (kp4f1a1/2*2 .lt. kp4f1a1) kp4f1a1 = kp4f1a1+1
      kf11 = +kp4f1a1
      kf11 = kf11+kf10
      if (kf11.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf11 memory needed = ',kf11
      stop 55555
      endif
      kp4f1a1 = kf11-kp4f1a1
c ........ file idxlgnvb ........
c ........ array imaplg ........
      kp4f2a1 = neqb
      if (kp4f2a1/2*2 .lt. kp4f2a1) kp4f2a1 = kp4f2a1+1
      kf12 = +kp4f2a1
      kf12 = kf12+kf11
      if (kf12.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf12 memory needed = ',kf12
      stop 55555
      endif
      kp4f2a1 = kf12-kp4f2a1
c ........ file idb0 ........
c ........ array id ........
      kp4f3a1 = knode1*kdgofb
      if (kp4f3a1/2*2 .lt. kp4f3a1) kp4f3a1 = kp4f3a1+1
      kf13 = +kp4f3a1
      kf13 = kf13+kf12
      if (kf13.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf13 memory needed = ',kf13
      stop 55555
      endif
      kp4f3a1 = kf13-kp4f3a1
c ........ file elemb0 ........
c ........ array node ........
      kp4f4a1 = kelemb
      if (kp4f4a1/2*2 .lt. kp4f4a1) kp4f4a1 = kp4f4a1+1
c ........ array emate ........
      kp4f4a2 = kemateb
      kp4f4a2 = kp4f4a2*2
      kf14 = +kp4f4a1+kp4f4a2
      kf14 = kf14+kf13
      if (kf14.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf14 memory needed = ',kf14
      stop 55555
      endif
      kp4f4a2 = kf14-kp4f4a2
      kp4f4a1 = kp4f4a2-kp4f4a1
c ........ file dispb0 ........
c ........ array disp0 ........
      kp4f5a1 = knode1*kdgofb
      kp4f5a1 = kp4f5a1*2
      kf15 = +kp4f5a1
      kf15 = kf15+kf14
      if (kf15.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf15 memory needed = ',kf15
      stop 55555
      endif
      kp4f5a1 = kf15-kp4f5a1
c ........ file unodb ........
c ........ array u0 ........
      kp4f6a1 = knode1*kdgofb
      kp4f6a1 = kp4f6a1*2
c ........ array u1 ........
      kp4f6a2 = knode1*kdgofb
      kp4f6a2 = kp4f6a2*2
c ........ array u2 ........
      kp4f6a3 = knode1*kdgofb
      kp4f6a3 = kp4f6a3*2
c ........ array u3 ........
      kp4f6a4 = knode1*kdgofb
      kp4f6a4 = kp4f6a4*2
      kf16 = +kp4f6a1+kp4f6a2+kp4f6a3+kp4f6a4
      kf16 = kf16+kf15
      if (kf16.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf16 memory needed = ',kf16
      stop 55555
      endif
      kp4f6a4 = kf16-kp4f6a4
      kp4f6a3 = kp4f6a4-kp4f6a3
      kp4f6a2 = kp4f6a3-kp4f6a2
      kp4f6a1 = kp4f6a2-kp4f6a1
c ............... temporary array  ...............
c ........ array tttt ........
      kp4b1 = 1
      kp4b1 = kp4b1*2
      kp4b0 = 1 
      kp4b1 = kp4b1+kp4b0
      call sgetpart(ia(kp4f1a1),ia(kp4f2a1),ia(kp4f3a1),
     & ia(kp4f4a1),ia(kp4f4a2),ia(kp4f5a1),ia(kp4f6a1),
     & ia(kp4f6a2),ia(kp4f6a3),ia(kp4f6a4),ib(kp4b0),
     & knode1,kdgofb,kelemb,kemateb,initnob,numetypb,iblk,nsource,
     & neqb,n_updateb,n_externalb,n_dataorgb,
     & mnodeb,nnodeb,mmateb,nmateb)
c ..........................................................
      call recvint(iblk,nsource,nblk)
      call recvr(iblk,nsource,t0)
      call recvr(iblk,nsource,tmax)
      call recvr(iblk,nsource,dt)
      call recvint(iblk,nsource,knode2)
      call recvint(iblk,nsource,kcoor)
c ........ file coor2 ........
c ........ array coor ........
      kp5f1a1 = knode2*kcoor
      kp5f1a1 = kp5f1a1*2
      kf17 = +kp5f1a1
      kf17 = kf17+kf16
      if (kf17.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf17 memory needed = ',kf17
      stop 55555
      endif
      kp5f1a1 = kf17-kp5f1a1
c ........ file order2 ........
c ........ array iorder ........
      kp5f2a1 = knode2
      if (kp5f2a1/2*2 .lt. kp5f2a1) kp5f2a1 = kp5f2a1+1
      kf18 = +kp5f2a1
      kf18 = kf18+kf17
      if (kf18.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf18 memory needed = ',kf18
      stop 55555
      endif
      kp5f2a1 = kf18-kp5f2a1
c ............... temporary array  ...............
c ........ array tttt ........
      kp5b1 = 1
      kp5b1 = kp5b1*2
      kp5b0 = 1 
      kp5b1 = kp5b1+kp5b0
      call spart(ia(kp5f1a1),ia(kp5f2a1),ib(kp5b0),
     & nsource,iblk,knode2,kcoor,
     & mnodea,nnodea,mmatea,nmatea,
     & mnodeb,nnodeb,mmateb,nmateb,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
c ...........................................
      maxac = maxa
      numetypc = 1
      call recvint(iblk,nsource,kelemc)
      call recvint(iblk,nsource,kematec)
      kdgofc = 3
      kvarc = knode2*kdgofc+1
      call recvint(iblk,nsource,neqc)
      initnoc = 0
c ........ file lgnvc ........
c ........ array maplg ........
      kp6f1a1 = neqc
      if (kp6f1a1/2*2 .lt. kp6f1a1) kp6f1a1 = kp6f1a1+1
      kf19 = +kp6f1a1
      kf19 = kf19+kf18
      if (kf19.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf19 memory needed = ',kf19
      stop 55555
      endif
      kp6f1a1 = kf19-kp6f1a1
c ........ file idxlgnvc ........
c ........ array imaplg ........
      kp6f2a1 = neqc
      if (kp6f2a1/2*2 .lt. kp6f2a1) kp6f2a1 = kp6f2a1+1
      kf20 = +kp6f2a1
      kf20 = kf20+kf19
      if (kf20.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf20 memory needed = ',kf20
      stop 55555
      endif
      kp6f2a1 = kf20-kp6f2a1
c ........ file idc0 ........
c ........ array id ........
      kp6f3a1 = knode2*kdgofc
      if (kp6f3a1/2*2 .lt. kp6f3a1) kp6f3a1 = kp6f3a1+1
      kf21 = +kp6f3a1
      kf21 = kf21+kf20
      if (kf21.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf21 memory needed = ',kf21
      stop 55555
      endif
      kp6f3a1 = kf21-kp6f3a1
c ........ file elemc0 ........
c ........ array node ........
      kp6f4a1 = kelemc
      if (kp6f4a1/2*2 .lt. kp6f4a1) kp6f4a1 = kp6f4a1+1
c ........ array emate ........
      kp6f4a2 = kematec
      kp6f4a2 = kp6f4a2*2
      kf22 = +kp6f4a1+kp6f4a2
      kf22 = kf22+kf21
      if (kf22.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf22 memory needed = ',kf22
      stop 55555
      endif
      kp6f4a2 = kf22-kp6f4a2
      kp6f4a1 = kp6f4a2-kp6f4a1
c ........ file dispc0 ........
c ........ array disp0 ........
      kp6f5a1 = knode2*kdgofc
      kp6f5a1 = kp6f5a1*2
      kf23 = +kp6f5a1
      kf23 = kf23+kf22
      if (kf23.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf23 memory needed = ',kf23
      stop 55555
      endif
      kp6f5a1 = kf23-kp6f5a1
c ........ file unodc ........
c ........ array u0 ........
      kp6f6a1 = knode2*kdgofc
      kp6f6a1 = kp6f6a1*2
c ........ array u1 ........
      kp6f6a2 = knode2*kdgofc
      kp6f6a2 = kp6f6a2*2
c ........ array u2 ........
      kp6f6a3 = knode2*kdgofc
      kp6f6a3 = kp6f6a3*2
c ........ array u3 ........
      kp6f6a4 = knode2*kdgofc
      kp6f6a4 = kp6f6a4*2
      kf24 = +kp6f6a1+kp6f6a2+kp6f6a3+kp6f6a4
      kf24 = kf24+kf23
      if (kf24.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf24 memory needed = ',kf24
      stop 55555
      endif
      kp6f6a4 = kf24-kp6f6a4
      kp6f6a3 = kp6f6a4-kp6f6a3
      kp6f6a2 = kp6f6a3-kp6f6a2
      kp6f6a1 = kp6f6a2-kp6f6a1
c ............... temporary array  ...............
c ........ array tttt ........
      kp6b1 = 1
      kp6b1 = kp6b1*2
      kp6b0 = 1 
      kp6b1 = kp6b1+kp6b0
      call sgetpart(ia(kp6f1a1),ia(kp6f2a1),ia(kp6f3a1),
     & ia(kp6f4a1),ia(kp6f4a2),ia(kp6f5a1),ia(kp6f6a1),
     & ia(kp6f6a2),ia(kp6f6a3),ia(kp6f6a4),ib(kp6b0),
     & knode2,kdgofc,kelemc,kematec,initnoc,numetypc,iblk,nsource,
     & neqc,n_updatec,n_externalc,n_dataorgc,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
c ........ file ida0 ........
c ........ array id ........
      kp7f1a1 = knode*kdgofa
      if (kp7f1a1/2*2 .lt. kp7f1a1) kp7f1a1 = kp7f1a1+1
      kp7f1a1 = kf5-kp7f1a1

      ktemp = kf4-kp7f1a1
      if (ktemp.lt.0) then 
      kp7f1a1 = kp7f1a1+ktemp
      endif
c ........ file nva ........
c ........ array nodvar ........
      kp7f2a1 = knode*kdgofa
      if (kp7f2a1/2*2 .lt. kp7f2a1) kp7f2a1 = kp7f2a1+1
      kf25 = +kp7f2a1
      kf25 = kf25+kf24
      if (kf25.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf25 memory needed = ',kf25
      stop 55555
      endif
      kp7f2a1 = kf25-kp7f2a1
c ........ file coor0 ........
c ........ array coor ........
      kp7f3a1 = knode*kcoor
      kp7f3a1 = kp7f3a1*2
      kp7f3a1 = kf1-kp7f3a1

      ktemp = kf0-kp7f3a1
      if (ktemp.lt.0) then 
      kp7f3a1 = kp7f3a1+ktemp
      endif
c ........ file dispa0 ........
c ........ array dsp0 ........
      kp7f4a1 = knode*kdgofa
      kp7f4a1 = kp7f4a1*2
      kp7f4a1 = kf7-kp7f4a1

      ktemp = kf6-kp7f4a1
      if (ktemp.lt.0) then 
      kp7f4a1 = kp7f4a1+ktemp
      endif
c ........ file bfda ........
c ........ array bfu ........
      kp7f5a1 = knode*kdgofa
      kp7f5a1 = kp7f5a1*2
      kf26 = +kp7f5a1
      kf26 = kf26+kf25
      if (kf26.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf26 memory needed = ',kf26
      stop 55555
      endif
      kp7f5a1 = kf26-kp7f5a1
c ........ file unoda ........
c ........ array u0 ........
      kp7f6a1 = knode*kdgofa
      kp7f6a1 = kp7f6a1*2
c ........ array u1 ........
      kp7f6a2 = knode*kdgofa
      kp7f6a2 = kp7f6a2*2
c ........ array u2 ........
      kp7f6a3 = knode*kdgofa
      kp7f6a3 = kp7f6a3*2
c ........ array u3 ........
      kp7f6a4 = knode*kdgofa
      kp7f6a4 = kp7f6a4*2
      kp7f6a4 = kf8-kp7f6a4
      kp7f6a3 = kp7f6a4-kp7f6a3
      kp7f6a2 = kp7f6a3-kp7f6a2
      kp7f6a1 = kp7f6a2-kp7f6a1

      ktemp = kf7-kp7f6a1
      if (ktemp.lt.0) then 
      kp7f6a1 = kp7f6a1+ktemp
      kp7f6a2 = kp7f6a2+ktemp
      kp7f6a3 = kp7f6a3+ktemp
      kp7f6a4 = kp7f6a4+ktemp
      endif
c ........ file unodd ........
c ........ array um ........
      kp7f7a1 = knode
      kp7f7a1 = kp7f7a1*2
c ........ array vm ........
      kp7f7a2 = knode
      kp7f7a2 = kp7f7a2*2
c ........ array wm ........
      kp7f7a3 = knode
      kp7f7a3 = kp7f7a3*2
c ........ array umn ........
      kp7f7a4 = knode
      kp7f7a4 = kp7f7a4*2
c ........ array vmn ........
      kp7f7a5 = knode
      kp7f7a5 = kp7f7a5*2
c ........ array wmn ........
      kp7f7a6 = knode
      kp7f7a6 = kp7f7a6*2
c ........ array evlx ........
      kp7f7a7 = knode
      kp7f7a7 = kp7f7a7*2
c ........ array evly ........
      kp7f7a8 = knode
      kp7f7a8 = kp7f7a8*2
c ........ array evlz ........
      kp7f7a9 = knode
      kp7f7a9 = kp7f7a9*2
c ........ array evlxn ........
      kp7f7a10 = knode
      kp7f7a10 = kp7f7a10*2
c ........ array evlyn ........
      kp7f7a11 = knode
      kp7f7a11 = kp7f7a11*2
c ........ array evlzn ........
      kp7f7a12 = knode
      kp7f7a12 = kp7f7a12*2
      kf27 = +kp7f7a1+kp7f7a2+kp7f7a3+kp7f7a4+kp7f7a5
     & +kp7f7a6+kp7f7a7+kp7f7a8+kp7f7a9+kp7f7a10+kp7f7a11+kp7f7a12
      kf27 = kf27+kf26
      if (kf27.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf27 memory needed = ',kf27
      stop 55555
      endif
      kp7f7a12 = kf27-kp7f7a12
      kp7f7a11 = kp7f7a12-kp7f7a11
      kp7f7a10 = kp7f7a11-kp7f7a10
      kp7f7a9 = kp7f7a10-kp7f7a9
      kp7f7a8 = kp7f7a9-kp7f7a8
      kp7f7a7 = kp7f7a8-kp7f7a7
      kp7f7a6 = kp7f7a7-kp7f7a6
      kp7f7a5 = kp7f7a6-kp7f7a5
      kp7f7a4 = kp7f7a5-kp7f7a4
      kp7f7a3 = kp7f7a4-kp7f7a3
      kp7f7a2 = kp7f7a3-kp7f7a2
      kp7f7a1 = kp7f7a2-kp7f7a1
c ........ file unodr ........
c ........ array ur ........
      kp7f8a1 = knode*kcoor
      kp7f8a1 = kp7f8a1*2
      kf28 = +kp7f8a1
      kf28 = kf28+kf27
      if (kf28.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf28 memory needed = ',kf28
      stop 55555
      endif
      kp7f8a1 = kf28-kp7f8a1
c ........ file unodac ........
c ........ array ua ........
      kp7f9a1 = knode*kcoor
      kp7f9a1 = kp7f9a1*2
      kf29 = +kp7f9a1
      kf29 = kf29+kf28
      if (kf29.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf29 memory needed = ',kf29
      stop 55555
      endif
      kp7f9a1 = kf29-kp7f9a1
c ........ file elema0 ........
c ........ array node ........
      kp7f10a1 = kelema
      if (kp7f10a1/2*2 .lt. kp7f10a1) kp7f10a1 = kp7f10a1+1
c ........ array emate ........
      kp7f10a2 = kematea
      kp7f10a2 = kp7f10a2*2
      kp7f10a2 = kf6-kp7f10a2
      kp7f10a1 = kp7f10a2-kp7f10a1

      ktemp = kf5-kp7f10a1
      if (ktemp.lt.0) then 
      kp7f10a1 = kp7f10a1+ktemp
      kp7f10a2 = kp7f10a2+ktemp
      endif
c ........ file jdiaga ........
c ........ array jdiag ........
      kp7f11a1 = kvara
      if (kp7f11a1/2*2 .lt. kp7f11a1) kp7f11a1 = kp7f11a1+1
c ........ array jdiagaz ........
      kp7f11a2 = kvara
      if (kp7f11a2/2*2 .lt. kp7f11a2) kp7f11a2 = kp7f11a2+1
      kf30 = +kp7f11a1+kp7f11a2
      kf30 = kf30+kf29
      if (kf30.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf30 memory needed = ',kf30
      stop 55555
      endif
      kp7f11a2 = kf30-kp7f11a2
      kp7f11a1 = kp7f11a2-kp7f11a1
c ........ file naa ........
c ........ array na ........
      kp7f12a1 = maxaa
      if (kp7f12a1/2*2 .lt. kp7f12a1) kp7f12a1 = kp7f12a1+1
      kf31 = +kp7f12a1
      kf31 = kf31+kf30
      if (kf31.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf31 memory needed = ',kf31
      stop 55555
      endif
      kp7f12a1 = kf31-kp7f12a1
c ............... temporary array  ...............
c ........ array r ........
      kp7b1 = 3
      kp7b1 = kp7b1*2
c ........ array lm ........
      kp7b2 = 1000
      if (kp7b2/2*2 .lt. kp7b2) kp7b2 = kp7b2+1
c ........ array jdiagr ........
      kp7b3 = kvara
      if (kp7b3/2*2 .lt. kp7b3) kp7b3 = kp7b3+1
      kp7b0 = 1 
      kp7b1 = kp7b1+kp7b0
      kp7b2 = kp7b2+kp7b1
      kp7b3 = kp7b3+kp7b2
      call starta(ia(kp7f1a1),ia(kp7f2a1),ia(kp7f3a1),
     & ia(kp7f4a1),ia(kp7f5a1),ia(kp7f6a1),ia(kp7f6a2),
     & ia(kp7f6a3),ia(kp7f6a4),ia(kp7f7a1),ia(kp7f7a2),
     & ia(kp7f7a3),ia(kp7f7a4),ia(kp7f7a5),ia(kp7f7a6),
     & ia(kp7f7a7),ia(kp7f7a8),ia(kp7f7a9),ia(kp7f7a10),
     & ia(kp7f7a11),ia(kp7f7a12),ia(kp7f8a1),ia(kp7f9a1),
     & ia(kp7f10a1),ia(kp7f10a2),ia(kp7f11a1),ia(kp7f11a2),
     & ia(kp7f12a1),ib(kp7b0),ib(kp7b1),ib(kp7b2),
     & nsource,iblk,knode,kdgofa,maxaa,maxba,neqa,neq1a,
     & kvara,kelema,kematea,kcoor,numetypa,numelema,time,material,
     & itn,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
c -----------------------------------------------
c ........ file naa ........
c ........ array na ........
      kp7f12a1 = maxaa
      if (kp7f12a1/2*2 .lt. kp7f12a1) kp7f12a1 = kp7f12a1+1
      kf31 = +kp7f12a1
      kf31 = kf31+kf30
c -----------------------------------------------
c ........ file idc0 ........
c ........ array id ........
      kp8f1a1 = knode2*kdgofc
      if (kp8f1a1/2*2 .lt. kp8f1a1) kp8f1a1 = kp8f1a1+1
      kp8f1a1 = kf21-kp8f1a1

      ktemp = kf20-kp8f1a1
      if (ktemp.lt.0) then 
      kp8f1a1 = kp8f1a1+ktemp
      endif
c ........ file nvc ........
c ........ array nodvar ........
      kp8f2a1 = knode2*kdgofc
      if (kp8f2a1/2*2 .lt. kp8f2a1) kp8f2a1 = kp8f2a1+1
      kf32 = +kp8f2a1
      kf32 = kf32+kf31
      if (kf32.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf32 memory needed = ',kf32
      stop 55555
      endif
      kp8f2a1 = kf32-kp8f2a1
c ........ file coor2 ........
c ........ array coor ........
      kp8f3a1 = knode2*kcoor
      kp8f3a1 = kp8f3a1*2
      kp8f3a1 = kf17-kp8f3a1

      ktemp = kf16-kp8f3a1
      if (ktemp.lt.0) then 
      kp8f3a1 = kp8f3a1+ktemp
      endif
c ........ file dispc0 ........
c ........ array dsp0 ........
      kp8f4a1 = knode2*kdgofc
      kp8f4a1 = kp8f4a1*2
      kp8f4a1 = kf23-kp8f4a1

      ktemp = kf22-kp8f4a1
      if (ktemp.lt.0) then 
      kp8f4a1 = kp8f4a1+ktemp
      endif
c ........ file bfdc ........
c ........ array bfu ........
      kp8f5a1 = knode2*kdgofc
      kp8f5a1 = kp8f5a1*2
      kf33 = +kp8f5a1
      kf33 = kf33+kf32
      if (kf33.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf33 memory needed = ',kf33
      stop 55555
      endif
      kp8f5a1 = kf33-kp8f5a1
c ........ file unodc ........
c ........ array u0 ........
      kp8f6a1 = knode2*kdgofc
      kp8f6a1 = kp8f6a1*2
c ........ array u1 ........
      kp8f6a2 = knode2*kdgofc
      kp8f6a2 = kp8f6a2*2
c ........ array u2 ........
      kp8f6a3 = knode2*kdgofc
      kp8f6a3 = kp8f6a3*2
c ........ array u3 ........
      kp8f6a4 = knode2*kdgofc
      kp8f6a4 = kp8f6a4*2
      kp8f6a4 = kf24-kp8f6a4
      kp8f6a3 = kp8f6a4-kp8f6a3
      kp8f6a2 = kp8f6a3-kp8f6a2
      kp8f6a1 = kp8f6a2-kp8f6a1

      ktemp = kf23-kp8f6a1
      if (ktemp.lt.0) then 
      kp8f6a1 = kp8f6a1+ktemp
      kp8f6a2 = kp8f6a2+ktemp
      kp8f6a3 = kp8f6a3+ktemp
      kp8f6a4 = kp8f6a4+ktemp
      endif
c ........ file elemc0 ........
c ........ array node ........
      kp8f7a1 = kelemc
      if (kp8f7a1/2*2 .lt. kp8f7a1) kp8f7a1 = kp8f7a1+1
c ........ array emate ........
      kp8f7a2 = kematec
      kp8f7a2 = kp8f7a2*2
      kp8f7a2 = kf22-kp8f7a2
      kp8f7a1 = kp8f7a2-kp8f7a1

      ktemp = kf21-kp8f7a1
      if (ktemp.lt.0) then 
      kp8f7a1 = kp8f7a1+ktemp
      kp8f7a2 = kp8f7a2+ktemp
      endif
c ........ file jdiagc ........
c ........ array jdiag ........
      kp8f8a1 = kvarc
      if (kp8f8a1/2*2 .lt. kp8f8a1) kp8f8a1 = kp8f8a1+1
c ........ array jdiagaz ........
      kp8f8a2 = kvarc
      if (kp8f8a2/2*2 .lt. kp8f8a2) kp8f8a2 = kp8f8a2+1
      kf34 = +kp8f8a1+kp8f8a2
      kf34 = kf34+kf33
      if (kf34.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf34 memory needed = ',kf34
      stop 55555
      endif
      kp8f8a2 = kf34-kp8f8a2
      kp8f8a1 = kp8f8a2-kp8f8a1
c ........ file nac ........
c ........ array na ........
      kp8f9a1 = maxac
      if (kp8f9a1/2*2 .lt. kp8f9a1) kp8f9a1 = kp8f9a1+1
      kf35 = +kp8f9a1
      kf35 = kf35+kf34
      if (kf35.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf35 memory needed = ',kf35
      stop 55555
      endif
      kp8f9a1 = kf35-kp8f9a1
c ............... temporary array  ...............
c ........ array r ........
      kp8b1 = 3
      kp8b1 = kp8b1*2
c ........ array lm ........
      kp8b2 = 1000
      if (kp8b2/2*2 .lt. kp8b2) kp8b2 = kp8b2+1
c ........ array jdiagr ........
      kp8b3 = kvarc
      if (kp8b3/2*2 .lt. kp8b3) kp8b3 = kp8b3+1
      kp8b0 = 1 
      kp8b1 = kp8b1+kp8b0
      kp8b2 = kp8b2+kp8b1
      kp8b3 = kp8b3+kp8b2
      call startc(ia(kp8f1a1),ia(kp8f2a1),ia(kp8f3a1),
     & ia(kp8f4a1),ia(kp8f5a1),ia(kp8f6a1),ia(kp8f6a2),
     & ia(kp8f6a3),ia(kp8f6a4),ia(kp8f7a1),ia(kp8f7a2),
     & ia(kp8f8a1),ia(kp8f8a2),ia(kp8f9a1),ib(kp8b0),
     & ib(kp8b1),ib(kp8b2),nsource,iblk,knode2,kdgofc,
     & maxac,maxbc,neqc,neq1c,kvarc,kelemc,kematec,kcoor,
     & numetypc,numelemc,time,material,itn,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
c -----------------------------------------------
c ........ file nac ........
c ........ array na ........
      kp8f9a1 = maxac
      if (kp8f9a1/2*2 .lt. kp8f9a1) kp8f9a1 = kp8f9a1+1
      kf35 = +kp8f9a1
      kf35 = kf35+kf34
c -----------------------------------------------
c ........ file coor0 ........
c ........ array u0 ........
      kp9f1a1 = knode*kcoor
      kp9f1a1 = kp9f1a1*2
      kp9f1a1 = kf1-kp9f1a1

      ktemp = kf0-kp9f1a1
      if (ktemp.lt.0) then 
      kp9f1a1 = kp9f1a1+ktemp
      endif
c ........ file coor00 ........
c ........ array u1 ........
      kp9f2a1 = knode*kcoor
      kp9f2a1 = kp9f2a1*2
      kf36 = +kp9f2a1
      kf36 = kf36+kf35
      if (kf36.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf36 memory needed = ',kf36
      stop 55555
      endif
      kp9f2a1 = kf36-kp9f2a1
c ............... temporary array  ...............
      kp9b0 = 1 
      call cpcoor0(ia(kp9f1a1),ia(kp9f2a1),knode,kcoor,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
c ........ file coor0 ........
c ........ array u0 ........
      kp10f1a1 = knode*kcoor
      kp10f1a1 = kp10f1a1*2
      kp10f1a1 = kf1-kp10f1a1

      ktemp = kf0-kp10f1a1
      if (ktemp.lt.0) then 
      kp10f1a1 = kp10f1a1+ktemp
      endif
c ........ file cooreule ........
c ........ array u1 ........
      kp10f2a1 = knode*kcoor
      kp10f2a1 = kp10f2a1*2
      kf37 = +kp10f2a1
      kf37 = kf37+kf36
      if (kf37.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf37 memory needed = ',kf37
      stop 55555
      endif
      kp10f2a1 = kf37-kp10f2a1
c ............... temporary array  ...............
      kp10b0 = 1 
      call cpcoore(ia(kp10f1a1),ia(kp10f2a1),knode,kcoor,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
c ........ file unod ........
c ........ array u0 ........
      kp11f1a1 = knode*kdgofa
      kp11f1a1 = kp11f1a1*2
      kf38 = +kp11f1a1
      kf38 = kf38+kf37
      if (kf38.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf38 memory needed = ',kf38
      stop 55555
      endif
      kp11f1a1 = kf38-kp11f1a1
c ........ file unod_0 ........
c ........ array u1 ........
      kp11f2a1 = knode*kdgofa
      kp11f2a1 = kp11f2a1*2
      kf39 = +kp11f2a1
      kf39 = kf39+kf38
      if (kf39.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf39 memory needed = ',kf39
      stop 55555
      endif
      kp11f2a1 = kf39-kp11f2a1
c ............... temporary array  ...............
      kp11b0 = 1 
      call cpunod(ia(kp11f1a1),ia(kp11f2a1),knode,kdgofa,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
c ........ file unodd ........
c ........ array um ........
      kp12f1a1 = knode
      kp12f1a1 = kp12f1a1*2
c ........ array vm ........
      kp12f1a2 = knode
      kp12f1a2 = kp12f1a2*2
c ........ array wm ........
      kp12f1a3 = knode
      kp12f1a3 = kp12f1a3*2
c ........ array umn ........
      kp12f1a4 = knode
      kp12f1a4 = kp12f1a4*2
c ........ array vmn ........
      kp12f1a5 = knode
      kp12f1a5 = kp12f1a5*2
c ........ array wmn ........
      kp12f1a6 = knode
      kp12f1a6 = kp12f1a6*2
c ........ array vlx ........
      kp12f1a7 = knode
      kp12f1a7 = kp12f1a7*2
c ........ array vly ........
      kp12f1a8 = knode
      kp12f1a8 = kp12f1a8*2
c ........ array vlz ........
      kp12f1a9 = knode
      kp12f1a9 = kp12f1a9*2
c ........ array vlxn ........
      kp12f1a10 = knode
      kp12f1a10 = kp12f1a10*2
c ........ array vlyn ........
      kp12f1a11 = knode
      kp12f1a11 = kp12f1a11*2
c ........ array vlzn ........
      kp12f1a12 = knode
      kp12f1a12 = kp12f1a12*2
      kp12f1a12 = kf27-kp12f1a12
      kp12f1a11 = kp12f1a12-kp12f1a11
      kp12f1a10 = kp12f1a11-kp12f1a10
      kp12f1a9 = kp12f1a10-kp12f1a9
      kp12f1a8 = kp12f1a9-kp12f1a8
      kp12f1a7 = kp12f1a8-kp12f1a7
      kp12f1a6 = kp12f1a7-kp12f1a6
      kp12f1a5 = kp12f1a6-kp12f1a5
      kp12f1a4 = kp12f1a5-kp12f1a4
      kp12f1a3 = kp12f1a4-kp12f1a3
      kp12f1a2 = kp12f1a3-kp12f1a2
      kp12f1a1 = kp12f1a2-kp12f1a1

      ktemp = kf26-kp12f1a1
      if (ktemp.lt.0) then 
      kp12f1a1 = kp12f1a1+ktemp
      kp12f1a2 = kp12f1a2+ktemp
      kp12f1a3 = kp12f1a3+ktemp
      kp12f1a4 = kp12f1a4+ktemp
      kp12f1a5 = kp12f1a5+ktemp
      kp12f1a6 = kp12f1a6+ktemp
      kp12f1a7 = kp12f1a7+ktemp
      kp12f1a8 = kp12f1a8+ktemp
      kp12f1a9 = kp12f1a9+ktemp
      kp12f1a10 = kp12f1a10+ktemp
      kp12f1a11 = kp12f1a11+ktemp
      kp12f1a12 = kp12f1a12+ktemp
      endif
c ........ file unodd_0 ........
c ........ array um1 ........
      kp12f2a1 = knode
      kp12f2a1 = kp12f2a1*2
c ........ array vm1 ........
      kp12f2a2 = knode
      kp12f2a2 = kp12f2a2*2
c ........ array wm1 ........
      kp12f2a3 = knode
      kp12f2a3 = kp12f2a3*2
c ........ array umn1 ........
      kp12f2a4 = knode
      kp12f2a4 = kp12f2a4*2
c ........ array vmn1 ........
      kp12f2a5 = knode
      kp12f2a5 = kp12f2a5*2
c ........ array wmn1 ........
      kp12f2a6 = knode
      kp12f2a6 = kp12f2a6*2
c ........ array vlx1 ........
      kp12f2a7 = knode
      kp12f2a7 = kp12f2a7*2
c ........ array vly1 ........
      kp12f2a8 = knode
      kp12f2a8 = kp12f2a8*2
c ........ array vlz1 ........
      kp12f2a9 = knode
      kp12f2a9 = kp12f2a9*2
c ........ array vlxn1 ........
      kp12f2a10 = knode
      kp12f2a10 = kp12f2a10*2
c ........ array vlyn1 ........
      kp12f2a11 = knode
      kp12f2a11 = kp12f2a11*2
c ........ array vlzn1 ........
      kp12f2a12 = knode
      kp12f2a12 = kp12f2a12*2
      kf40 = +kp12f2a1+kp12f2a2+kp12f2a3+kp12f2a4+kp12f2a5
     & +kp12f2a6+kp12f2a7+kp12f2a8+kp12f2a9+kp12f2a10+kp12f2a11+kp12f2a12
      kf40 = kf40+kf39
      if (kf40.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf40 memory needed = ',kf40
      stop 55555
      endif
      kp12f2a12 = kf40-kp12f2a12
      kp12f2a11 = kp12f2a12-kp12f2a11
      kp12f2a10 = kp12f2a11-kp12f2a10
      kp12f2a9 = kp12f2a10-kp12f2a9
      kp12f2a8 = kp12f2a9-kp12f2a8
      kp12f2a7 = kp12f2a8-kp12f2a7
      kp12f2a6 = kp12f2a7-kp12f2a6
      kp12f2a5 = kp12f2a6-kp12f2a5
      kp12f2a4 = kp12f2a5-kp12f2a4
      kp12f2a3 = kp12f2a4-kp12f2a3
      kp12f2a2 = kp12f2a3-kp12f2a2
      kp12f2a1 = kp12f2a2-kp12f2a1
c ............... temporary array  ...............
      kp12b0 = 1 
      call cpunodd(ia(kp12f1a1),ia(kp12f1a2),ia(kp12f1a3),
     & ia(kp12f1a4),ia(kp12f1a5),ia(kp12f1a6),ia(kp12f1a7),
     & ia(kp12f1a8),ia(kp12f1a9),ia(kp12f1a10),ia(kp12f1a11),
     & ia(kp12f1a12),ia(kp12f2a1),ia(kp12f2a2),ia(kp12f2a3),
     & ia(kp12f2a4),ia(kp12f2a5),ia(kp12f2a6),ia(kp12f2a7),
     & ia(kp12f2a8),ia(kp12f2a9),ia(kp12f2a10),ia(kp12f2a11),
     & ia(kp12f2a12),knode,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
cpunodc c c
      MSstop = 0
1     continue
c ........ file bfda ........
c ........ array bf ........
      kp13f1a1 = knode *kdgofa
      kp13f1a1 = kp13f1a1*2
      kp13f1a1 = kf26-kp13f1a1

      ktemp = kf25-kp13f1a1
      if (ktemp.lt.0) then 
      kp13f1a1 = kp13f1a1+ktemp
      endif
c ............... temporary array  ...............
      kp13b0 = 1 
      call bfts(ia(kp13f1a1),iblk,nsource,tmax,dt,
     & time,it,msstop,knode ,kdgofa,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
      MSend_tol = 0
      itn_tol=1
2     continue
      MSend = 0
      itn=1
3     continue
c ........ file nva ........
c ........ array nodvar ........
      kp14f1a1 = knode*kdgofa
      if (kp14f1a1/2*2 .lt. kp14f1a1) kp14f1a1 = kp14f1a1+1
      kp14f1a1 = kf25-kp14f1a1

      ktemp = kf24-kp14f1a1
      if (ktemp.lt.0) then 
      kp14f1a1 = kp14f1a1+ktemp
      endif
c ........ file coor0 ........
c ........ array coor ........
      kp14f2a1 = knode*kcoor
      kp14f2a1 = kp14f2a1*2
      kp14f2a1 = kf1-kp14f2a1

      ktemp = kf0-kp14f2a1
      if (ktemp.lt.0) then 
      kp14f2a1 = kp14f2a1+ktemp
      endif
c ........ file bfda ........
c ........ array bfu ........
      kp14f3a1 = knode*kdgofa
      kp14f3a1 = kp14f3a1*2
      kp14f3a1 = kf26-kp14f3a1

      ktemp = kf25-kp14f3a1
      if (ktemp.lt.0) then 
      kp14f3a1 = kp14f3a1+ktemp
      endif
c ........ file unoda ........
c ........ array eu ........
      kp14f4a1 = knode*kdgofa
      kp14f4a1 = kp14f4a1*2
c ........ array eu1 ........
      kp14f4a2 = knode*kdgofa
      kp14f4a2 = kp14f4a2*2
      kp14f4a2 = kf8-kp14f4a2
      kp14f4a1 = kp14f4a2-kp14f4a1

      ktemp = kf7-kp14f4a1
      if (ktemp.lt.0) then 
      kp14f4a1 = kp14f4a1+ktemp
      kp14f4a2 = kp14f4a2+ktemp
      endif
c ........ file unodd ........
c ........ array eum ........
      kp14f5a1 = knode
      kp14f5a1 = kp14f5a1*2
c ........ array evm ........
      kp14f5a2 = knode
      kp14f5a2 = kp14f5a2*2
c ........ array ewm ........
      kp14f5a3 = knode
      kp14f5a3 = kp14f5a3*2
c ........ array eumn ........
      kp14f5a4 = knode
      kp14f5a4 = kp14f5a4*2
c ........ array evmn ........
      kp14f5a5 = knode
      kp14f5a5 = kp14f5a5*2
c ........ array ewmn ........
      kp14f5a6 = knode
      kp14f5a6 = kp14f5a6*2
c ........ array evlx ........
      kp14f5a7 = knode
      kp14f5a7 = kp14f5a7*2
c ........ array evly ........
      kp14f5a8 = knode
      kp14f5a8 = kp14f5a8*2
c ........ array evlz ........
      kp14f5a9 = knode
      kp14f5a9 = kp14f5a9*2
      kp14f5a9 = kf27-kp14f5a9
      kp14f5a8 = kp14f5a9-kp14f5a8
      kp14f5a7 = kp14f5a8-kp14f5a7
      kp14f5a6 = kp14f5a7-kp14f5a6
      kp14f5a5 = kp14f5a6-kp14f5a5
      kp14f5a4 = kp14f5a5-kp14f5a4
      kp14f5a3 = kp14f5a4-kp14f5a3
      kp14f5a2 = kp14f5a3-kp14f5a2
      kp14f5a1 = kp14f5a2-kp14f5a1

      ktemp = kf26-kp14f5a1
      if (ktemp.lt.0) then 
      kp14f5a1 = kp14f5a1+ktemp
      kp14f5a2 = kp14f5a2+ktemp
      kp14f5a3 = kp14f5a3+ktemp
      kp14f5a4 = kp14f5a4+ktemp
      kp14f5a5 = kp14f5a5+ktemp
      kp14f5a6 = kp14f5a6+ktemp
      kp14f5a7 = kp14f5a7+ktemp
      kp14f5a8 = kp14f5a8+ktemp
      kp14f5a9 = kp14f5a9+ktemp
      endif
c ........ file jdiaga ........
c ........ array jdiag ........
      kp14f6a1 = kvara
      if (kp14f6a1/2*2 .lt. kp14f6a1) kp14f6a1 = kp14f6a1+1
c ........ array jdiagaz ........
      kp14f6a2 = kvara
      if (kp14f6a2/2*2 .lt. kp14f6a2) kp14f6a2 = kp14f6a2+1
      kp14f6a2 = kf30-kp14f6a2
      kp14f6a1 = kp14f6a2-kp14f6a1

      ktemp = kf29-kp14f6a1
      if (ktemp.lt.0) then 
      kp14f6a1 = kp14f6a1+ktemp
      kp14f6a2 = kp14f6a2+ktemp
      endif
c ........ file naa ........
c ........ array na ........
      kp14f7a1 = maxaa
      if (kp14f7a1/2*2 .lt. kp14f7a1) kp14f7a1 = kp14f7a1+1
      kp14f7a1 = kf31-kp14f7a1

      ktemp = kf30-kp14f7a1
      if (ktemp.lt.0) then 
      kp14f7a1 = kp14f7a1+ktemp
      endif
c ........ file elema0 ........
c ........ array node ........
      kp14f8a1 = kelema
      if (kp14f8a1/2*2 .lt. kp14f8a1) kp14f8a1 = kp14f8a1+1
c ........ array emate ........
      kp14f8a2 = kematea
      kp14f8a2 = kp14f8a2*2
      kp14f8a2 = kf6-kp14f8a2
      kp14f8a1 = kp14f8a2-kp14f8a1

      ktemp = kf5-kp14f8a1
      if (ktemp.lt.0) then 
      kp14f8a1 = kp14f8a1+ktemp
      kp14f8a2 = kp14f8a2+ktemp
      endif
c ........ file jrafa ........
c ........ array a ........
      kp14f9a1 = maxaa
      kp14f9a1 = kp14f9a1*2
c ........ array f ........
      kp14f9a2 = kvara
      kp14f9a2 = kp14f9a2*2
      kf41 = +kp14f9a1+kp14f9a2
      kf41 = kf41+kf40
      if (kf41.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf41 memory needed = ',kf41
      stop 55555
      endif
      kp14f9a2 = kf41-kp14f9a2
      kp14f9a1 = kp14f9a2-kp14f9a1
c ............... temporary array  ...............
c ........ array u ........
      kp14b1 = knode*kdgofa
      kp14b1 = kp14b1*2
c ........ array ev ........
      kp14b2 = knode*kdgofa
      kp14b2 = kp14b2*2
c ........ array eue ........
      kp14b3 = knode*kdgofa
      kp14b3 = kp14b3*2
c ........ array edu ........
      kp14b4 = knode*kdgofa
      kp14b4 = kp14b4*2
c ........ array evf ........
      kp14b5 = knode*kdgofa
      kp14b5 = kp14b5*2
c ........ array evlxn ........
      kp14b6 = knode
      kp14b6 = kp14b6*2
c ........ array evlyn ........
      kp14b7 = knode
      kp14b7 = kp14b7*2
c ........ array evlzn ........
      kp14b8 = knode
      kp14b8 = kp14b8*2
c ........ array emass ........
      kp14b9 = kvara
      kp14b9 = kp14b9*2
c ........ array sml ........
      kp14b10 = 100000
      kp14b10 = kp14b10*2
c ........ array u1 ........
      kp14b11 = kvara
      kp14b11 = kp14b11*2
c ........ array lm ........
      kp14b12 = 1000
      if (kp14b12/2*2 .lt. kp14b12) kp14b12 = kp14b12+1
c ........ array uu ........
      kp14b13 = kvara
      kp14b13 = kp14b13*2
      kp14b0 = 1 
      kp14b1 = kp14b1+kp14b0
      kp14b2 = kp14b2+kp14b1
      kp14b3 = kp14b3+kp14b2
      kp14b4 = kp14b4+kp14b3
      kp14b5 = kp14b5+kp14b4
      kp14b6 = kp14b6+kp14b5
      kp14b7 = kp14b7+kp14b6
      kp14b8 = kp14b8+kp14b7
      kp14b9 = kp14b9+kp14b8
      kp14b10 = kp14b10+kp14b9
      kp14b11 = kp14b11+kp14b10
      kp14b12 = kp14b12+kp14b11
      kp14b13 = kp14b13+kp14b12
      call efsia(ia(kp14f1a1),ia(kp14f2a1),ia(kp14f3a1),
     & ia(kp14f4a1),ia(kp14f4a2),ia(kp14f5a1),ia(kp14f5a2),
     & ia(kp14f5a3),ia(kp14f5a4),ia(kp14f5a5),ia(kp14f5a6),
     & ia(kp14f5a7),ia(kp14f5a8),ia(kp14f5a9),ia(kp14f6a1),
     & ia(kp14f6a2),ia(kp14f7a1),ia(kp14f8a1),ia(kp14f8a2),
     & ia(kp14f9a1),ia(kp14f9a2),ib(kp14b0),ib(kp14b1),
     & ib(kp14b2),ib(kp14b3),ib(kp14b4),ib(kp14b5),
     & ib(kp14b6),ib(kp14b7),ib(kp14b8),ib(kp14b9),
     & ib(kp14b10),ib(kp14b11),ib(kp14b12),iblk,kdgofa,
     & knode,maxaa,neqa,neq1a,kvara,kelema,kematea,kcoor,
     & time,dt,it,itn,cc,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
c ........ file jdiaga ........
c ........ array jdiag1 ........
      kp15f1a1 = kvara
      if (kp15f1a1/2*2 .lt. kp15f1a1) kp15f1a1 = kp15f1a1+1
c ........ array jdiagaz1 ........
      kp15f1a2 = kvara
      if (kp15f1a2/2*2 .lt. kp15f1a2) kp15f1a2 = kp15f1a2+1
      kp15f1a2 = kf30-kp15f1a2
      kp15f1a1 = kp15f1a2-kp15f1a1

      ktemp = kf29-kp15f1a1
      if (ktemp.lt.0) then 
      kp15f1a1 = kp15f1a1+ktemp
      kp15f1a2 = kp15f1a2+ktemp
      endif
c ........ file jrafa ........
c ........ array val ........
      kp15f2a1 = maxaa
      kp15f2a1 = kp15f2a1*2
c ........ array f ........
      kp15f2a2 = kvara
      kp15f2a2 = kp15f2a2*2
      kp15f2a2 = kf41-kp15f2a2
      kp15f2a1 = kp15f2a2-kp15f2a1
c ........ file naa ........
c ........ array na ........
      kp15f3a1 = maxaa
      if (kp15f3a1/2*2 .lt. kp15f3a1) kp15f3a1 = kp15f3a1+1
      kp15f3a1 = kf31-kp15f3a1

      ktemp = kf30-kp15f3a1
      if (ktemp.lt.0) then 
      kp15f3a1 = kp15f3a1+ktemp
      endif
c ........ file lgnva ........
c ........ array maplg ........
      kp15f4a1 = neqa
      if (kp15f4a1/2*2 .lt. kp15f4a1) kp15f4a1 = kp15f4a1+1
      kp15f4a1 = kf3-kp15f4a1

      ktemp = kf2-kp15f4a1
      if (ktemp.lt.0) then 
      kp15f4a1 = kp15f4a1+ktemp
      endif
c ........ file idxlgnva ........
c ........ array imaplg ........
      kp15f5a1 = neqa
      if (kp15f5a1/2*2 .lt. kp15f5a1) kp15f5a1 = kp15f5a1+1
      kp15f5a1 = kf4-kp15f5a1

      ktemp = kf3-kp15f5a1
      if (ktemp.lt.0) then 
      kp15f5a1 = kp15f5a1+ktemp
      endif
c ........ file ua ........
c ........ array x ........
      kp15f6a1 = neqa
      kp15f6a1 = kp15f6a1*2
      kf42 = +kp15f6a1
      kf42 = kf42+kf41
      if (kf42.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf42 memory needed = ',kf42
      stop 55555
      endif
      kp15f6a1 = kf42-kp15f6a1
c ............... temporary array  ...............
c ........ array diag ........
      kp15b1 = kvara
      kp15b1 = kp15b1*2
c ........ array nbindex ........
      kp15b2 = maxba
      if (kp15b2/2*2 .lt. kp15b2) kp15b2 = kp15b2+1
c ........ array b ........
      kp15b3 = neqa
      kp15b3 = kp15b3*2
c ........ array nupdate ........
      kp15b4 = n_updatea
      if (kp15b4/2*2 .lt. kp15b4) kp15b4 = kp15b4+1
c ........ array ndata_org ........
      kp15b5 = n_dataorga
      if (kp15b5/2*2 .lt. kp15b5) kp15b5 = kp15b5+1
c ........ array nexternal ........
      kp15b6 = n_externala
      if (kp15b6/2*2 .lt. kp15b6) kp15b6 = kp15b6+1
c ........ array nextern_indx ........
      kp15b7 = n_externala
      if (kp15b7/2*2 .lt. kp15b7) kp15b7 = kp15b7+1
c ........ array nupdate_indx ........
      kp15b8 = n_updatea
      if (kp15b8/2*2 .lt. kp15b8) kp15b8 = kp15b8+1
c ........ array nodvar ........
      kp15b9 = knode*kdgofa
      if (kp15b9/2*2 .lt. kp15b9) kp15b9 = kp15b9+1
c ........ array jdiag ........
      kp15b10 = kvara
      if (kp15b10/2*2 .lt. kp15b10) kp15b10 = kp15b10+1
c ........ array jdiagaz ........
      kp15b11 = kvara
      if (kp15b11/2*2 .lt. kp15b11) kp15b11 = kp15b11+1
      kp15b0 = 1 
      kp15b1 = kp15b1+kp15b0
      kp15b2 = kp15b2+kp15b1
      kp15b3 = kp15b3+kp15b2
      kp15b4 = kp15b4+kp15b3
      kp15b5 = kp15b5+kp15b4
      kp15b6 = kp15b6+kp15b5
      kp15b7 = kp15b7+kp15b6
      kp15b8 = kp15b8+kp15b7
      kp15b9 = kp15b9+kp15b8
      kp15b10 = kp15b10+kp15b9
      kp15b11 = kp15b11+kp15b10
      call sfsia(ia(kp15f1a1),ia(kp15f1a2),ia(kp15f2a1),
     & ia(kp15f2a2),ia(kp15f3a1),ia(kp15f4a1),ia(kp15f5a1),
     & ia(kp15f6a1),ib(kp15b0),ib(kp15b1),ib(kp15b2),
     & ib(kp15b3),ib(kp15b4),ib(kp15b5),ib(kp15b6),
     & ib(kp15b7),ib(kp15b8),ib(kp15b9),ib(kp15b10),
     & knode,kdgofa,neqa,neq1a,n_updatea,n_externala,n_dataorga,kvara,
     & maxaa,maxba,iblk,nsource,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
c ........ file order0 ........
c ........ array iorder ........
      kp16f1a1 = knode
      if (kp16f1a1/2*2 .lt. kp16f1a1) kp16f1a1 = kp16f1a1+1
      kp16f1a1 = kf2-kp16f1a1

      ktemp = kf1-kp16f1a1
      if (ktemp.lt.0) then 
      kp16f1a1 = kp16f1a1+ktemp
      endif
c ........ file order1 ........
c ........ array iorder1 ........
      kp16f2a1 = knode1
      if (kp16f2a1/2*2 .lt. kp16f2a1) kp16f2a1 = kp16f2a1+1
      kp16f2a1 = kf10-kp16f2a1

      ktemp = kf9-kp16f2a1
      if (ktemp.lt.0) then 
      kp16f2a1 = kp16f2a1+ktemp
      endif
c ........ file order2 ........
c ........ array iorder2 ........
      kp16f3a1 = knode2
      if (kp16f3a1/2*2 .lt. kp16f3a1) kp16f3a1 = kp16f3a1+1
      kp16f3a1 = kf18-kp16f3a1

      ktemp = kf17-kp16f3a1
      if (ktemp.lt.0) then 
      kp16f3a1 = kp16f3a1+ktemp
      endif
c ........ file coor0 ........
c ........ array coor ........
      kp16f4a1 = knode*kcoor
      kp16f4a1 = kp16f4a1*2
      kp16f4a1 = kf1-kp16f4a1

      ktemp = kf0-kp16f4a1
      if (ktemp.lt.0) then 
      kp16f4a1 = kp16f4a1+ktemp
      endif
c ........ file bfda ........
c ........ array bf ........
      kp16f5a1 = knode*kdgofa
      kp16f5a1 = kp16f5a1*2
      kp16f5a1 = kf26-kp16f5a1

      ktemp = kf25-kp16f5a1
      if (ktemp.lt.0) then 
      kp16f5a1 = kp16f5a1+ktemp
      endif
c ........ file nva ........
c ........ array nodvar ........
      kp16f6a1 = knode*kdgofa
      if (kp16f6a1/2*2 .lt. kp16f6a1) kp16f6a1 = kp16f6a1+1
      kp16f6a1 = kf25-kp16f6a1

      ktemp = kf24-kp16f6a1
      if (ktemp.lt.0) then 
      kp16f6a1 = kp16f6a1+ktemp
      endif
c ........ file ua ........
c ........ array uvar ........
      kp16f7a1 = neqa
      kp16f7a1 = kp16f7a1*2
      kp16f7a1 = kf42-kp16f7a1
c ........ file unoda ........
c ........ array eu ........
      kp16f8a1 = knode*kdgofa
      kp16f8a1 = kp16f8a1*2
c ........ array eu1 ........
      kp16f8a2 = knode*kdgofa
      kp16f8a2 = kp16f8a2*2
c ........ array edu ........
      kp16f8a3 = knode*kdgofa
      kp16f8a3 = kp16f8a3*2
      kp16f8a3 = kf8-kp16f8a3
      kp16f8a2 = kp16f8a3-kp16f8a2
      kp16f8a1 = kp16f8a2-kp16f8a1

      ktemp = kf7-kp16f8a1
      if (ktemp.lt.0) then 
      kp16f8a1 = kp16f8a1+ktemp
      kp16f8a2 = kp16f8a2+ktemp
      kp16f8a3 = kp16f8a3+ktemp
      endif
c ........ file coor00 ........
c ........ array coor0 ........
      kp16f9a1 = knode*kcoor
      kp16f9a1 = kp16f9a1*2
      kp16f9a1 = kf36-kp16f9a1

      ktemp = kf35-kp16f9a1
      if (ktemp.lt.0) then 
      kp16f9a1 = kp16f9a1+ktemp
      endif
c ........ file unodd ........
c ........ array eum ........
      kp16f10a1 = knode
      kp16f10a1 = kp16f10a1*2
c ........ array evm ........
      kp16f10a2 = knode
      kp16f10a2 = kp16f10a2*2
c ........ array ewm ........
      kp16f10a3 = knode
      kp16f10a3 = kp16f10a3*2
c ........ array eumn ........
      kp16f10a4 = knode
      kp16f10a4 = kp16f10a4*2
c ........ array evmn ........
      kp16f10a5 = knode
      kp16f10a5 = kp16f10a5*2
c ........ array ewmn ........
      kp16f10a6 = knode
      kp16f10a6 = kp16f10a6*2
c ........ array evlx ........
      kp16f10a7 = knode
      kp16f10a7 = kp16f10a7*2
c ........ array evly ........
      kp16f10a8 = knode
      kp16f10a8 = kp16f10a8*2
c ........ array evlz ........
      kp16f10a9 = knode
      kp16f10a9 = kp16f10a9*2
c ........ array evlxn ........
      kp16f10a10 = knode
      kp16f10a10 = kp16f10a10*2
c ........ array evlyn ........
      kp16f10a11 = knode
      kp16f10a11 = kp16f10a11*2
c ........ array evlzn ........
      kp16f10a12 = knode
      kp16f10a12 = kp16f10a12*2
      kp16f10a12 = kf27-kp16f10a12
      kp16f10a11 = kp16f10a12-kp16f10a11
      kp16f10a10 = kp16f10a11-kp16f10a10
      kp16f10a9 = kp16f10a10-kp16f10a9
      kp16f10a8 = kp16f10a9-kp16f10a8
      kp16f10a7 = kp16f10a8-kp16f10a7
      kp16f10a6 = kp16f10a7-kp16f10a6
      kp16f10a5 = kp16f10a6-kp16f10a5
      kp16f10a4 = kp16f10a5-kp16f10a4
      kp16f10a3 = kp16f10a4-kp16f10a3
      kp16f10a2 = kp16f10a3-kp16f10a2
      kp16f10a1 = kp16f10a2-kp16f10a1

      ktemp = kf26-kp16f10a1
      if (ktemp.lt.0) then 
      kp16f10a1 = kp16f10a1+ktemp
      kp16f10a2 = kp16f10a2+ktemp
      kp16f10a3 = kp16f10a3+ktemp
      kp16f10a4 = kp16f10a4+ktemp
      kp16f10a5 = kp16f10a5+ktemp
      kp16f10a6 = kp16f10a6+ktemp
      kp16f10a7 = kp16f10a7+ktemp
      kp16f10a8 = kp16f10a8+ktemp
      kp16f10a9 = kp16f10a9+ktemp
      kp16f10a10 = kp16f10a10+ktemp
      kp16f10a11 = kp16f10a11+ktemp
      kp16f10a12 = kp16f10a12+ktemp
      endif
c ........ file bfdc ........
c ........ array ubc ........
      kp16f11a1 = knode2*kdgofc
      kp16f11a1 = kp16f11a1*2
      kp16f11a1 = kf33-kp16f11a1

      ktemp = kf32-kp16f11a1
      if (ktemp.lt.0) then 
      kp16f11a1 = kp16f11a1+ktemp
      endif
c ........ file unodac ........
c ........ array accel ........
      kp16f12a1 = knode1*kdgofb
      kp16f12a1 = kp16f12a1*2
      kp16f12a1 = kf29-kp16f12a1

      ktemp = kf28-kp16f12a1
      if (ktemp.lt.0) then 
      kp16f12a1 = kp16f12a1+ktemp
      endif
c ........ file unods ........
c ........ array ums ........
      kp16f13a1 = knode1
      kp16f13a1 = kp16f13a1*2
c ........ array vms ........
      kp16f13a2 = knode1
      kp16f13a2 = kp16f13a2*2
c ........ array wms ........
      kp16f13a3 = knode1
      kp16f13a3 = kp16f13a3*2
c ........ array vs ........
      kp16f13a4 = knode1*kdgofb
      kp16f13a4 = kp16f13a4*2
c ........ array acs ........
      kp16f13a5 = knode1*kdgofb
      kp16f13a5 = kp16f13a5*2
      kf43 = +kp16f13a1+kp16f13a2+kp16f13a3+kp16f13a4+kp16f13a5
      kf43 = kf43+kf42
      if (kf43.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf43 memory needed = ',kf43
      stop 55555
      endif
      kp16f13a5 = kf43-kp16f13a5
      kp16f13a4 = kp16f13a5-kp16f13a4
      kp16f13a3 = kp16f13a4-kp16f13a3
      kp16f13a2 = kp16f13a3-kp16f13a2
      kp16f13a1 = kp16f13a2-kp16f13a1
c ........ file unodf ........
c ........ array evf ........
      kp16f14a1 = knode*kdgofa
      kp16f14a1 = kp16f14a1*2
      kf44 = +kp16f14a1
      kf44 = kf44+kf43
      if (kf44.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf44 memory needed = ',kf44
      stop 55555
      endif
      kp16f14a1 = kf44-kp16f14a1
c ............... temporary array  ...............
c ........ array ev ........
      kp16b1 = knode*kdgofa
      kp16b1 = kp16b1*2
c ........ array eue ........
      kp16b2 = knode*kdgofa
      kp16b2 = kp16b2*2
c ........ array apool ........
      kp16b3 = kvara
      kp16b3 = kp16b3*2
      kp16b0 = 1 
      kp16b1 = kp16b1+kp16b0
      kp16b2 = kp16b2+kp16b1
      kp16b3 = kp16b3+kp16b2
      call ufsia(ia(kp16f1a1),ia(kp16f2a1),ia(kp16f3a1),
     & ia(kp16f4a1),ia(kp16f5a1),ia(kp16f6a1),ia(kp16f7a1),
     & ia(kp16f8a1),ia(kp16f8a2),ia(kp16f8a3),ia(kp16f9a1),
     & ia(kp16f10a1),ia(kp16f10a2),ia(kp16f10a3),ia(kp16f10a4),
     & ia(kp16f10a5),ia(kp16f10a6),ia(kp16f10a7),ia(kp16f10a8),
     & ia(kp16f10a9),ia(kp16f10a10),ia(kp16f10a11),ia(kp16f10a12),
     & ia(kp16f11a1),ia(kp16f12a1),ia(kp16f13a1),ia(kp16f13a2),
     & ia(kp16f13a3),ia(kp16f13a4),ia(kp16f13a5),ia(kp16f14a1),
     & ib(kp16b0),ib(kp16b1),ib(kp16b2),knode,knode1,
     & knode2,kdgofa,kdgofb,kdgofc,kvara,neqa,kcoor,time,
     & dt,it,itn,cc,nsource,iblk,msend,
     & mnodea,nnodea,mmatea,nmatea)
c ..........................................................
      if (MSend.eq.0) goto 3
      knodeg=knode
      kdgofg=kdgofa
      kdgof2=kdgofc
      kvarg=knodeg*kdgofg
c ........ file order0 ........
c ........ array iorder ........
      kp17f1a1 = knodeg
      if (kp17f1a1/2*2 .lt. kp17f1a1) kp17f1a1 = kp17f1a1+1
      kp17f1a1 = kf2-kp17f1a1

      ktemp = kf1-kp17f1a1
      if (ktemp.lt.0) then 
      kp17f1a1 = kp17f1a1+ktemp
      endif
c ........ file order1 ........
c ........ array iorder1 ........
      kp17f2a1 = knode1
      if (kp17f2a1/2*2 .lt. kp17f2a1) kp17f2a1 = kp17f2a1+1
      kp17f2a1 = kf10-kp17f2a1

      ktemp = kf9-kp17f2a1
      if (ktemp.lt.0) then 
      kp17f2a1 = kp17f2a1+ktemp
      endif
c ........ file order2 ........
c ........ array iorder2 ........
      kp17f3a1 = knode2
      if (kp17f3a1/2*2 .lt. kp17f3a1) kp17f3a1 = kp17f3a1+1
      kp17f3a1 = kf18-kp17f3a1

      ktemp = kf17-kp17f3a1
      if (ktemp.lt.0) then 
      kp17f3a1 = kp17f3a1+ktemp
      endif
c ........ file coor1 ........
c ........ array coor ........
      kp17f4a1 = knode1*kcoor
      kp17f4a1 = kp17f4a1*2
      kp17f4a1 = kf9-kp17f4a1

      ktemp = kf8-kp17f4a1
      if (ktemp.lt.0) then 
      kp17f4a1 = kp17f4a1+ktemp
      endif
c ........ file unods ........
c ........ array eud ........
      kp17f5a1 = knode1
      kp17f5a1 = kp17f5a1*2
c ........ array evd ........
      kp17f5a2 = knode1
      kp17f5a2 = kp17f5a2*2
c ........ array ewd ........
      kp17f5a3 = knode1
      kp17f5a3 = kp17f5a3*2
c ........ array eus ........
      kp17f5a4 = knode1
      kp17f5a4 = kp17f5a4*2
c ........ array evs ........
      kp17f5a5 = knode1
      kp17f5a5 = kp17f5a5*2
c ........ array ews ........
      kp17f5a6 = knode1
      kp17f5a6 = kp17f5a6*2
      kp17f5a6 = kf43-kp17f5a6
      kp17f5a5 = kp17f5a6-kp17f5a5
      kp17f5a4 = kp17f5a5-kp17f5a4
      kp17f5a3 = kp17f5a4-kp17f5a3
      kp17f5a2 = kp17f5a3-kp17f5a2
      kp17f5a1 = kp17f5a2-kp17f5a1

      ktemp = kf42-kp17f5a1
      if (ktemp.lt.0) then 
      kp17f5a1 = kp17f5a1+ktemp
      kp17f5a2 = kp17f5a2+ktemp
      kp17f5a3 = kp17f5a3+ktemp
      kp17f5a4 = kp17f5a4+ktemp
      kp17f5a5 = kp17f5a5+ktemp
      kp17f5a6 = kp17f5a6+ktemp
      endif
c ........ file elemb0 ........
c ........ array node ........
      kp17f6a1 = kelemb
      if (kp17f6a1/2*2 .lt. kp17f6a1) kp17f6a1 = kp17f6a1+1
c ........ array emate ........
      kp17f6a2 = kemateb
      kp17f6a2 = kp17f6a2*2
      kp17f6a2 = kf14-kp17f6a2
      kp17f6a1 = kp17f6a2-kp17f6a1

      ktemp = kf13-kp17f6a1
      if (ktemp.lt.0) then 
      kp17f6a1 = kp17f6a1+ktemp
      kp17f6a2 = kp17f6a2+ktemp
      endif
c ........ file coor00 ........
c ........ array coor0 ........
      kp17f7a1 = knodeg*kcoor
      kp17f7a1 = kp17f7a1*2
      kp17f7a1 = kf36-kp17f7a1

      ktemp = kf35-kp17f7a1
      if (ktemp.lt.0) then 
      kp17f7a1 = kp17f7a1+ktemp
      endif
c ........ file unodr ........
c ........ array ur ........
      kp17f8a1 = knodeg*kdgofb
      kp17f8a1 = kp17f8a1*2
      kp17f8a1 = kf28-kp17f8a1

      ktemp = kf27-kp17f8a1
      if (ktemp.lt.0) then 
      kp17f8a1 = kp17f8a1+ktemp
      endif
c ........ file bfdc ........
c ........ array ubc ........
      kp17f9a1 = knode2*kdgof2
      kp17f9a1 = kp17f9a1*2
      kp17f9a1 = kf33-kp17f9a1

      ktemp = kf32-kp17f9a1
      if (ktemp.lt.0) then 
      kp17f9a1 = kp17f9a1+ktemp
      endif
c ............... temporary array  ...............
c ........ array nodvar ........
      kp17b1 = knode1*kdgofb
      if (kp17b1/2*2 .lt. kp17b1) kp17b1 = kp17b1+1
c ........ array u ........
      kp17b2 = knode1*kdgofb
      kp17b2 = kp17b2*2
c ........ array euf ........
      kp17b3 = knode2*kdgof2
      kp17b3 = kp17b3*2
c ........ array eu ........
      kp17b4 = knode1*kdgofb
      kp17b4 = kp17b4*2
c ........ array emass ........
      kp17b5 = kvarb
      kp17b5 = kp17b5*2
c ........ array f ........
      kp17b6 = kvarb
      kp17b6 = kp17b6*2
c ........ array sml ........
      kp17b7 = 100000
      kp17b7 = kp17b7*2
c ........ array apool ........
      kp17b8 = kvarg
      kp17b8 = kp17b8*2
      kp17b0 = 1 
      kp17b1 = kp17b1+kp17b0
      kp17b2 = kp17b2+kp17b1
      kp17b3 = kp17b3+kp17b2
      kp17b4 = kp17b4+kp17b3
      kp17b5 = kp17b5+kp17b4
      kp17b6 = kp17b6+kp17b5
      kp17b7 = kp17b7+kp17b6
      kp17b8 = kp17b8+kp17b7
      call efsib(ia(kp17f1a1),ia(kp17f2a1),ia(kp17f3a1),
     & ia(kp17f4a1),ia(kp17f5a1),ia(kp17f5a2),ia(kp17f5a3),
     & ia(kp17f5a4),ia(kp17f5a5),ia(kp17f5a6),ia(kp17f6a1),
     & ia(kp17f6a2),ia(kp17f7a1),ia(kp17f8a1),ia(kp17f9a1),
     & ib(kp17b0),ib(kp17b1),ib(kp17b2),ib(kp17b3),
     & ib(kp17b4),ib(kp17b5),ib(kp17b6),ib(kp17b7),
     & iblk,nsource,knode1,kdgofb,knodeg,kdgofg,knode2,kdgof2,
     & neqb,kvarb,kvarg,kelemb,kemateb,kcoor,time,dt,
     & it,itn,
     & mnodeb,nnodeb,mmateb,nmateb)
c ..........................................................
c ........ file nvc ........
c ........ array nodvar ........
      kp18f1a1 = knode2*kdgofc
      if (kp18f1a1/2*2 .lt. kp18f1a1) kp18f1a1 = kp18f1a1+1
      kp18f1a1 = kf32-kp18f1a1

      ktemp = kf31-kp18f1a1
      if (ktemp.lt.0) then 
      kp18f1a1 = kp18f1a1+ktemp
      endif
c ........ file coor2 ........
c ........ array coor ........
      kp18f2a1 = knode2*kcoor
      kp18f2a1 = kp18f2a1*2
      kp18f2a1 = kf17-kp18f2a1

      ktemp = kf16-kp18f2a1
      if (ktemp.lt.0) then 
      kp18f2a1 = kp18f2a1+ktemp
      endif
c ........ file bfdc ........
c ........ array bfu ........
      kp18f3a1 = knode2*kdgofc
      kp18f3a1 = kp18f3a1*2
      kp18f3a1 = kf33-kp18f3a1

      ktemp = kf32-kp18f3a1
      if (ktemp.lt.0) then 
      kp18f3a1 = kp18f3a1+ktemp
      endif
c ........ file jdiagc ........
c ........ array jdiag ........
      kp18f4a1 = kvarc
      if (kp18f4a1/2*2 .lt. kp18f4a1) kp18f4a1 = kp18f4a1+1
c ........ array jdiagaz ........
      kp18f4a2 = kvarc
      if (kp18f4a2/2*2 .lt. kp18f4a2) kp18f4a2 = kp18f4a2+1
      kp18f4a2 = kf34-kp18f4a2
      kp18f4a1 = kp18f4a2-kp18f4a1

      ktemp = kf33-kp18f4a1
      if (ktemp.lt.0) then 
      kp18f4a1 = kp18f4a1+ktemp
      kp18f4a2 = kp18f4a2+ktemp
      endif
c ........ file nac ........
c ........ array na ........
      kp18f5a1 = maxac
      if (kp18f5a1/2*2 .lt. kp18f5a1) kp18f5a1 = kp18f5a1+1
      kp18f5a1 = kf35-kp18f5a1

      ktemp = kf34-kp18f5a1
      if (ktemp.lt.0) then 
      kp18f5a1 = kp18f5a1+ktemp
      endif
c ........ file elemc0 ........
c ........ array node ........
      kp18f6a1 = kelemc
      if (kp18f6a1/2*2 .lt. kp18f6a1) kp18f6a1 = kp18f6a1+1
c ........ array emate ........
      kp18f6a2 = kematec
      kp18f6a2 = kp18f6a2*2
      kp18f6a2 = kf22-kp18f6a2
      kp18f6a1 = kp18f6a2-kp18f6a1

      ktemp = kf21-kp18f6a1
      if (ktemp.lt.0) then 
      kp18f6a1 = kp18f6a1+ktemp
      kp18f6a2 = kp18f6a2+ktemp
      endif
c ........ file jrafc ........
c ........ array a ........
      kp18f7a1 = maxac
      kp18f7a1 = kp18f7a1*2
c ........ array f ........
      kp18f7a2 = kvarc
      kp18f7a2 = kp18f7a2*2
      kf45 = +kp18f7a1+kp18f7a2
      kf45 = kf45+kf44
      if (kf45.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf45 memory needed = ',kf45
      stop 55555
      endif
      kp18f7a2 = kf45-kp18f7a2
      kp18f7a1 = kp18f7a2-kp18f7a1
c ............... temporary array  ...............
c ........ array u ........
      kp18b1 = knode2*kdgofc
      kp18b1 = kp18b1*2
c ........ array eu ........
      kp18b2 = knode2*kdgofc
      kp18b2 = kp18b2*2
c ........ array eu1 ........
      kp18b3 = knode2*kdgofc
      kp18b3 = kp18b3*2
c ........ array sml ........
      kp18b4 = 100000
      kp18b4 = kp18b4*2
c ........ array u1 ........
      kp18b5 = kvarc
      kp18b5 = kp18b5*2
c ........ array lm ........
      kp18b6 = 1000
      if (kp18b6/2*2 .lt. kp18b6) kp18b6 = kp18b6+1
c ........ array uu ........
      kp18b7 = kvarc
      kp18b7 = kp18b7*2
      kp18b0 = 1 
      kp18b1 = kp18b1+kp18b0
      kp18b2 = kp18b2+kp18b1
      kp18b3 = kp18b3+kp18b2
      kp18b4 = kp18b4+kp18b3
      kp18b5 = kp18b5+kp18b4
      kp18b6 = kp18b6+kp18b5
      kp18b7 = kp18b7+kp18b6
      call efsic(ia(kp18f1a1),ia(kp18f2a1),ia(kp18f3a1),
     & ia(kp18f4a1),ia(kp18f4a2),ia(kp18f5a1),ia(kp18f6a1),
     & ia(kp18f6a2),ia(kp18f7a1),ia(kp18f7a2),ib(kp18b0),
     & ib(kp18b1),ib(kp18b2),ib(kp18b3),ib(kp18b4),
     & ib(kp18b5),ib(kp18b6),iblk,kdgofc,knode2,maxac,
     & neqc,neq1c,kvarc,kelemc,kematec,kcoor,time,dt,
     & it,itn,cc,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
c ........ file jdiagc ........
c ........ array jdiag1 ........
      kp19f1a1 = kvarc
      if (kp19f1a1/2*2 .lt. kp19f1a1) kp19f1a1 = kp19f1a1+1
c ........ array jdiagaz1 ........
      kp19f1a2 = kvarc
      if (kp19f1a2/2*2 .lt. kp19f1a2) kp19f1a2 = kp19f1a2+1
      kp19f1a2 = kf34-kp19f1a2
      kp19f1a1 = kp19f1a2-kp19f1a1

      ktemp = kf33-kp19f1a1
      if (ktemp.lt.0) then 
      kp19f1a1 = kp19f1a1+ktemp
      kp19f1a2 = kp19f1a2+ktemp
      endif
c ........ file jrafc ........
c ........ array val ........
      kp19f2a1 = maxac
      kp19f2a1 = kp19f2a1*2
c ........ array f ........
      kp19f2a2 = kvarc
      kp19f2a2 = kp19f2a2*2
      kp19f2a2 = kf45-kp19f2a2
      kp19f2a1 = kp19f2a2-kp19f2a1
c ........ file nac ........
c ........ array na ........
      kp19f3a1 = maxac
      if (kp19f3a1/2*2 .lt. kp19f3a1) kp19f3a1 = kp19f3a1+1
      kp19f3a1 = kf35-kp19f3a1

      ktemp = kf34-kp19f3a1
      if (ktemp.lt.0) then 
      kp19f3a1 = kp19f3a1+ktemp
      endif
c ........ file lgnvc ........
c ........ array maplg ........
      kp19f4a1 = neqc
      if (kp19f4a1/2*2 .lt. kp19f4a1) kp19f4a1 = kp19f4a1+1
      kp19f4a1 = kf19-kp19f4a1

      ktemp = kf18-kp19f4a1
      if (ktemp.lt.0) then 
      kp19f4a1 = kp19f4a1+ktemp
      endif
c ........ file idxlgnvc ........
c ........ array imaplg ........
      kp19f5a1 = neqc
      if (kp19f5a1/2*2 .lt. kp19f5a1) kp19f5a1 = kp19f5a1+1
      kp19f5a1 = kf20-kp19f5a1

      ktemp = kf19-kp19f5a1
      if (ktemp.lt.0) then 
      kp19f5a1 = kp19f5a1+ktemp
      endif
c ........ file uc ........
c ........ array x ........
      kp19f6a1 = neqc
      kp19f6a1 = kp19f6a1*2
      kf46 = +kp19f6a1
      kf46 = kf46+kf45
      if (kf46.gt.150000000) then
      write(*,*) 'exceed memory ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'kf46 memory needed = ',kf46
      stop 55555
      endif
      kp19f6a1 = kf46-kp19f6a1
c ............... temporary array  ...............
c ........ array diag ........
      kp19b1 = kvarc
      kp19b1 = kp19b1*2
c ........ array nbindex ........
      kp19b2 = maxbc
      if (kp19b2/2*2 .lt. kp19b2) kp19b2 = kp19b2+1
c ........ array b ........
      kp19b3 = neqc
      kp19b3 = kp19b3*2
c ........ array nupdate ........
      kp19b4 = n_updatec
      if (kp19b4/2*2 .lt. kp19b4) kp19b4 = kp19b4+1
c ........ array ndata_org ........
      kp19b5 = n_dataorgc
      if (kp19b5/2*2 .lt. kp19b5) kp19b5 = kp19b5+1
c ........ array nexternal ........
      kp19b6 = n_externalc
      if (kp19b6/2*2 .lt. kp19b6) kp19b6 = kp19b6+1
c ........ array nextern_indx ........
      kp19b7 = n_externalc
      if (kp19b7/2*2 .lt. kp19b7) kp19b7 = kp19b7+1
c ........ array nupdate_indx ........
      kp19b8 = n_updatec
      if (kp19b8/2*2 .lt. kp19b8) kp19b8 = kp19b8+1
c ........ array nodvar ........
      kp19b9 = knode2*kdgofc
      if (kp19b9/2*2 .lt. kp19b9) kp19b9 = kp19b9+1
c ........ array jdiag ........
      kp19b10 = kvarc
      if (kp19b10/2*2 .lt. kp19b10) kp19b10 = kp19b10+1
c ........ array jdiagaz ........
      kp19b11 = kvarc
      if (kp19b11/2*2 .lt. kp19b11) kp19b11 = kp19b11+1
      kp19b0 = 1 
      kp19b1 = kp19b1+kp19b0
      kp19b2 = kp19b2+kp19b1
      kp19b3 = kp19b3+kp19b2
      kp19b4 = kp19b4+kp19b3
      kp19b5 = kp19b5+kp19b4
      kp19b6 = kp19b6+kp19b5
      kp19b7 = kp19b7+kp19b6
      kp19b8 = kp19b8+kp19b7
      kp19b9 = kp19b9+kp19b8
      kp19b10 = kp19b10+kp19b9
      kp19b11 = kp19b11+kp19b10
      call sfsic(ia(kp19f1a1),ia(kp19f1a2),ia(kp19f2a1),
     & ia(kp19f2a2),ia(kp19f3a1),ia(kp19f4a1),ia(kp19f5a1),
     & ia(kp19f6a1),ib(kp19b0),ib(kp19b1),ib(kp19b2),
     & ib(kp19b3),ib(kp19b4),ib(kp19b5),ib(kp19b6),
     & ib(kp19b7),ib(kp19b8),ib(kp19b9),ib(kp19b10),
     & knode2,kdgofc,neqc,neq1c,n_updatec,n_externalc,n_dataorgc,kvarc,
     & maxac,maxbc,iblk,nsource,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
      knodeg=knode
      kdgofg=kdgofa
      kvarg=knodeg*kdgofg
c ........ file order2 ........
c ........ array iorder ........
      kp20f1a1 = knode2
      if (kp20f1a1/2*2 .lt. kp20f1a1) kp20f1a1 = kp20f1a1+1
      kp20f1a1 = kf18-kp20f1a1

      ktemp = kf17-kp20f1a1
      if (ktemp.lt.0) then 
      kp20f1a1 = kp20f1a1+ktemp
      endif
c ........ file order1 ........
c ........ array iorder1 ........
      kp20f2a1 = knode1
      if (kp20f2a1/2*2 .lt. kp20f2a1) kp20f2a1 = kp20f2a1+1
      kp20f2a1 = kf10-kp20f2a1

      ktemp = kf9-kp20f2a1
      if (ktemp.lt.0) then 
      kp20f2a1 = kp20f2a1+ktemp
      endif
c ........ file order0 ........
c ........ array iorder0 ........
      kp20f3a1 = knodeg
      if (kp20f3a1/2*2 .lt. kp20f3a1) kp20f3a1 = kp20f3a1+1
      kp20f3a1 = kf2-kp20f3a1

      ktemp = kf1-kp20f3a1
      if (ktemp.lt.0) then 
      kp20f3a1 = kp20f3a1+ktemp
      endif
c ........ file coor2 ........
c ........ array coor ........
      kp20f4a1 = knode2*kcoor
      kp20f4a1 = kp20f4a1*2
      kp20f4a1 = kf17-kp20f4a1

      ktemp = kf16-kp20f4a1
      if (ktemp.lt.0) then 
      kp20f4a1 = kp20f4a1+ktemp
      endif
c ........ file bfdc ........
c ........ array bf ........
      kp20f5a1 = knode2*kdgofc
      kp20f5a1 = kp20f5a1*2
      kp20f5a1 = kf33-kp20f5a1

      ktemp = kf32-kp20f5a1
      if (ktemp.lt.0) then 
      kp20f5a1 = kp20f5a1+ktemp
      endif
c ........ file nvc ........
c ........ array nodvar ........
      kp20f6a1 = knode2*kdgofc
      if (kp20f6a1/2*2 .lt. kp20f6a1) kp20f6a1 = kp20f6a1+1
      kp20f6a1 = kf32-kp20f6a1

      ktemp = kf31-kp20f6a1
      if (ktemp.lt.0) then 
      kp20f6a1 = kp20f6a1+ktemp
      endif
c ........ file uc ........
c ........ array uvar ........
      kp20f7a1 = neqc
      kp20f7a1 = kp20f7a1*2
      kp20f7a1 = kf46-kp20f7a1
c ........ file unodc ........
c ........ array eu1 ........
      kp20f8a1 = knode2*kdgofc
      kp20f8a1 = kp20f8a1*2
      kp20f8a1 = kf24-kp20f8a1

      ktemp = kf23-kp20f8a1
      if (ktemp.lt.0) then 
      kp20f8a1 = kp20f8a1+ktemp
      endif
c ........ file unodr ........
c ........ array ur ........
      kp20f9a1 = knodeg*kdgofc
      kp20f9a1 = kp20f9a1*2
      kp20f9a1 = kf28-kp20f9a1

      ktemp = kf27-kp20f9a1
      if (ktemp.lt.0) then 
      kp20f9a1 = kp20f9a1+ktemp
      endif
c ........ file unodd ........
c ........ array ulap ........
      kp20f10a1 = knodeg*kdgofc
      kp20f10a1 = kp20f10a1*2
c ........ array ulapn ........
      kp20f10a2 = knodeg*kdgofc
      kp20f10a2 = kp20f10a2*2
c ........ array vlap ........
      kp20f10a3 = knodeg*kdgofc
      kp20f10a3 = kp20f10a3*2
c ........ array vlapn ........
      kp20f10a4 = knodeg*kdgofc
      kp20f10a4 = kp20f10a4*2
      kp20f10a4 = kf27-kp20f10a4
      kp20f10a3 = kp20f10a4-kp20f10a3
      kp20f10a2 = kp20f10a3-kp20f10a2
      kp20f10a1 = kp20f10a2-kp20f10a1

      ktemp = kf26-kp20f10a1
      if (ktemp.lt.0) then 
      kp20f10a1 = kp20f10a1+ktemp
      kp20f10a2 = kp20f10a2+ktemp
      kp20f10a3 = kp20f10a3+ktemp
      kp20f10a4 = kp20f10a4+ktemp
      endif
c ........ file coor0 ........
c ........ array coor0 ........
      kp20f11a1 = knodeg*kcoor
      kp20f11a1 = kp20f11a1*2
      kp20f11a1 = kf1-kp20f11a1

      ktemp = kf0-kp20f11a1
      if (ktemp.lt.0) then 
      kp20f11a1 = kp20f11a1+ktemp
      endif
c ........ file cooreule ........
c ........ array coorlap ........
      kp20f12a1 = knodeg*kcoor
      kp20f12a1 = kp20f12a1*2
      kp20f12a1 = kf37-kp20f12a1

      ktemp = kf36-kp20f12a1
      if (ktemp.lt.0) then 
      kp20f12a1 = kp20f12a1+ktemp
      endif
c ........ file unoda ........
c ........ array ua ........
      kp20f13a1 = knodeg*kdgofg
      kp20f13a1 = kp20f13a1*2
c ........ array u1a ........
      kp20f13a2 = knodeg*kdgofg
      kp20f13a2 = kp20f13a2*2
c ........ array dua ........
      kp20f13a3 = knodeg*kdgofg
      kp20f13a3 = kp20f13a3*2
      kp20f13a3 = kf8-kp20f13a3
      kp20f13a2 = kp20f13a3-kp20f13a2
      kp20f13a1 = kp20f13a2-kp20f13a1

      ktemp = kf7-kp20f13a1
      if (ktemp.lt.0) then 
      kp20f13a1 = kp20f13a1+ktemp
      kp20f13a2 = kp20f13a2+ktemp
      kp20f13a3 = kp20f13a3+ktemp
      endif
c ........ file unods ........
c ........ array ums ........
      kp20f14a1 = knode1
      kp20f14a1 = kp20f14a1*2
c ........ array vms ........
      kp20f14a2 = knode1
      kp20f14a2 = kp20f14a2*2
c ........ array wms ........
      kp20f14a3 = knode1
      kp20f14a3 = kp20f14a3*2
c ........ array vs ........
      kp20f14a4 = knode1*kdgofb
      kp20f14a4 = kp20f14a4*2
c ........ array acs ........
      kp20f14a5 = knode1*kdgofb
      kp20f14a5 = kp20f14a5*2
      kp20f14a5 = kf43-kp20f14a5
      kp20f14a4 = kp20f14a5-kp20f14a4
      kp20f14a3 = kp20f14a4-kp20f14a3
      kp20f14a2 = kp20f14a3-kp20f14a2
      kp20f14a1 = kp20f14a2-kp20f14a1

      ktemp = kf42-kp20f14a1
      if (ktemp.lt.0) then 
      kp20f14a1 = kp20f14a1+ktemp
      kp20f14a2 = kp20f14a2+ktemp
      kp20f14a3 = kp20f14a3+ktemp
      kp20f14a4 = kp20f14a4+ktemp
      kp20f14a5 = kp20f14a5+ktemp
      endif
c ........ file unodac ........
c ........ array ac ........
      kp20f15a1 = knode1*kdgofb
      kp20f15a1 = kp20f15a1*2
      kp20f15a1 = kf29-kp20f15a1

      ktemp = kf28-kp20f15a1
      if (ktemp.lt.0) then 
      kp20f15a1 = kp20f15a1+ktemp
      endif
c ............... temporary array  ...............
c ........ array eu ........
      kp20b1 = knode2*kdgofc
      kp20b1 = kp20b1*2
c ........ array apool ........
      kp20b2 = kvarg
      kp20b2 = kp20b2*2
      kp20b0 = 1 
      kp20b1 = kp20b1+kp20b0
      kp20b2 = kp20b2+kp20b1
      call ufsic(ia(kp20f1a1),ia(kp20f2a1),ia(kp20f3a1),
     & ia(kp20f4a1),ia(kp20f5a1),ia(kp20f6a1),ia(kp20f7a1),
     & ia(kp20f8a1),ia(kp20f9a1),ia(kp20f10a1),ia(kp20f10a2),
     & ia(kp20f10a3),ia(kp20f10a4),ia(kp20f11a1),ia(kp20f12a1),
     & ia(kp20f13a1),ia(kp20f13a2),ia(kp20f13a3),ia(kp20f14a1),
     & ia(kp20f14a2),ia(kp20f14a3),ia(kp20f14a4),ia(kp20f14a5),
     & ia(kp20f15a1),ib(kp20b0),ib(kp20b1),knode2,kdgofc,
     & kvarc,knodeg,kdgofg,kvarg,knode1,kdgofb,neqc,kcoor,
     & time,dt,it,itn_tol,cc,nsource,iblk,msend,
     & msend_tol,
     & mnodec,nnodec,mmatec,nmatec)
c ..........................................................
4     continue
      if (MSend_tol.eq.0) goto 2
call post.bat
      if (MSstop.eq.0) goto 1
      end
