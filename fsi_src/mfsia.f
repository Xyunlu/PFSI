      subroutine Mmfsia(filename,kkcom)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      common /jcom/ jcom(100)
      common /ia/ ia(150000000)

6     format (1x,10i5)
7     format (1x,6e12.5)
1001  format(1x,7i7)
      neqm=10

      kna1=neqm*1
      if (kna1/2*2 .lt. kna1) kna1=kna1+1
      kna5=764*1
      if (kna5/2*2 .lt. kna5) kna5=kna5+1
      kna3=neqm*1
      if (kna3/2*2 .lt. kna3) kna3=kna3+1
      kna4=neqm*1
      if (kna4/2*2 .lt. kna4) kna4=kna4+1
      kna9=neqm*1
      if (kna9/2*2 .lt. kna9) kna9=kna9+1
      kna10=neqm*1
      if (kna10/2*2 .lt. kna10) kna10=kna10+1
      kna2=neqm*2
      kna6=neqm*2
      kna7=neqm*2
      kna8=neqm*2
      jcom(1) = 1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
      kna7=kna7+kna6
      kna8=kna8+kna7
      kna9=kna9+kna8
      kna10=kna10+kna9
      if (kna10-1.gt.150000000) then
      write(*,*) 'exceed memory of array ia'
      write(*,*) 'memory of ia = 150000000'
      write(*,*) 'memory needed = ',kna10,' in prgram mfsia'
      stop 55555
      endif
      call mfsia(neqm,ia(kna0),ia(kna1),ia(kna2),
     *ia(kna3),ia(kna4),ia(kna5),ia(kna6),ia(kna7),
     *ia(kna8),ia(kna9),
     *filename)


      return
      end
      subroutine mfsia(neqm,nbindex,val,nupdate,
     *nexternal,ndata_org,bt,xt,xt1,nupdate_index,
     *nextern_index,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
       include "az_aztecf.h"
       include "mpif.h"
c
c....  define some array for call callaztec solver
c
       integer noptions(0:AZ_OPTIONS_SIZE)
       double precision params(0:AZ_PARAMS_SIZE)
      dimension nbindex(neqm),ndata_org(764),nupdate(neqm),
     *  nexternal(neqm),nupdate_index(neqm),nextern_index(neqm)
      dimension val(neqm)
      dimension bt(neqm),xt(neqm),xt1(neqm)

c@mbegin
      nsource = 0
      N_update = 1
      do i=1,1
        val(i)=1.d0
      enddo
      val(2) = 0.d0
c
      do i=1,1
        nupdate(i)=neqall+1
      enddo
      do i=1,1
        nbindex(i)=2
      enddo
      nbindex(2) = 2
      do i=1,1
        bt(i)=0.
      enddo
c
c      initialize AZTEC options
c
      call AZ_defaults(noptions, params)
      params(AZ_tol) = 1.0e-12
      params(AZ_tol) = 1.0e-8
c      noptions(AZ_solver) = AZ_cgs
c      noptions(AZ_solver) = AZ_gmres
c      noptions(AZ_solver) = AZ_bicgstab
c      noptions(AZ_precond) = AZ_sym_GS
c      noptions(AZ_precond) = AZ_dom_decomp
c      noptions(AZ_overlap) = AZ_diag
c      noptions(AZ_overlap) = 10
c      noptions(AZ_subdomain_solve) = AZ_ilut
c      noptions(AZ_scaling) = AZ_sym_diag
c      noptions(AZ_subdomain_solve) = AZ_ilu
c      params(AZ_ilut_fill) = 4
c      params(AZ_drop) = 1.0e-8
c      params(AZ_drop) = 0.0d0
cccccccccccccccccccccccccccccccccccccccccccccccc
c      noptions(AZ_max_iter)=100000
c      noptions(AZ_solver) = AZ_bicgstab
c      noptions(AZ_precond) = AZ_sym_GS
c      noptions(AZ_poly_ord) = 3
c      noptions(AZ_scaling) = AZ_sym_diag
      include "solverbld.opt"
cccccccccccccccccccccccccccccccccccccccccccccccc
c
      call callpsolver(0,0,0,763,0,nupdate,
     &                     nbindex,val,bt,xt,ndata_org,nexternal,
     &                     nupdate_index,nextern_index,noptions,params)
c
      return
      end
