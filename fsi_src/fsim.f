      subroutine masterp
      implicit real*8(a-h,o-z)
      character*12 filename(20)
      logical flg,iflg
      common /stop/ msstop,msend,msend_tol

      open (21,file='init.dat',status='old')
      read (21,*) xmax,ymax,zmax,xo,yo,zo,radius,thick,radaxis,angle,
     &             radart_mult,ratiol
      read (21,*) nx,nz,ns,npt,nct,nsc,naxis,nyi,nyo
      read (21,*) ai,bi,ao,bo
      read (21,*) u_in,visc_f,pe,pv,denf,dens,gravity,omega
      read (21,*) idpoint,timespring
      close(21)
      pi=datan(1.d0)*4.d0
      angle=pi/180.d0*angle
C
      zero = 0.D0
      OPEN (21,FILE='vangle',FORM='FORMATTED',STATUS='unknown')
      if (idpoint.ne.1.and.dabs(omega).gt.1.d-20) then
        WRITE(21,*) omega,omega,angle,angle
      else
        WRITE(21,*) zero,zero,angle,angle
      endif
      CLOSE (21)

      nsource = 0
      filename(1)='fsia.io'
      filename(2)='elem0'
      filename(3)='coor0'
      filename(4)='id0'
      call Mpartition0(filename,0)
      filename(1)='fsia.io'
      filename(2)='elem0'
      filename(3)='disp0'
      filename(4)='id0'
      filename(5)='unod'
      filename(6)='disp1'
      filename(7)='disp2'
      filename(8)='disp3'
      call Mgetpart(filename,0)
      filename(1)='fsib.io'
      filename(2)='elemb0'
      filename(3)='coor1'
      filename(4)='idb0'
      call Mpartition(filename,1)
      filename(1)='fsib.io'
      filename(2)='elemb0'
      filename(3)='dispb0'
      filename(4)='idb0'
      filename(5)='unodb'
      filename(6)='dispb1'
      filename(7)='dispb2'
      filename(8)='dispb3'
      call Mgetpart(filename,1)
      filename(1)='fsic.io'
      filename(2)='elemc0'
      filename(3)='coor2'
      filename(4)='idc0'
      call Mpartition(filename,2)
      filename(1)='fsic.io'
      filename(2)='elemc0'
      filename(3)='dispc0'
      filename(4)='idc0'
      filename(5)='unodc'
      filename(6)='dispc1'
      filename(7)='dispc2'
      filename(8)='dispc3'
      call Mgetpart(filename,2)
c .... startm
c .... startm
      call system("cp coor0 coor0.0")
      call system("cp coor0 cooreule")
C.... cpunod
C.... cpunodd
C.... cpunodc
      MSstop = 0
1     continue
c .... bftm
      filename(1)='time'
      call Mbftm(filename,0)
      MSend_tol = 0
      itn_tol=1
2     continue
      MSend = 0
      itn=1
3     continue
c.... e@0
      call Mmfsia(filename,0)
      filename(1)='disp0'
      filename(2)='dispb0'
      filename(3)='dispc0'
      call Mmufsia(filename,0)
      if (MSend.eq.0) goto 3
      filename(1)='dispb0'
      filename(2)='disp0'
      filename(3)='dispc0'
      filename(4)='elemb0'
      call Mmufsib(filename,0)
c.....me@0
      call Mmfsic(filename,2)
      filename(1)='dispc0'
      filename(2)='disp0'
      filename(3)='dispb0'
      call Mmufsic(filename,0)
4     continue
      if (MSend_tol.eq.0) goto 2
call post.bat
      if (MSstop.eq.0) goto 1
      end
