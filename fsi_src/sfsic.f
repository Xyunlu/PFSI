      subroutine sfsic(jdiag1,jdiagaz1,val,
     & f,na,maplg,imaplg,
     & x,diag,nbindex,b,
     & nupdate,ndata_org,nexternal,nextern_indx,
     & nupdate_indx,nodvar,jdiag,jdiagaz,
     & knode,kdgof,neq,neq1,n_update,n_external,n_dataorg,kvar,
     & maxa,maxb,iblk,nsource,
     & mnode,nnode,mmate,nmate)
      implicit real*8 (a-h,o-z)
      include "az_aztecf.h"
      include "mpif.h"
c
      dimension jdiag1(kvar),jdiagaz1(kvar),na(maxa),diag(kvar),
     &     f(kvar),maplg(neq),imaplg(neq),
     &     nbindex(maxb),val(maxa),x(neq),b(neq),
     &     nupdate(n_update),ndata_org(n_dataorg),nexternal(n_external),
     &     nextern_indx(n_external),
     &     nupdate_indx(n_update),
     &     nodvar(knode,kdgof)
      dimension jdiag(kvar),jdiagaz(kvar)
c      dimension a1(maxb),na1(maxb),jdiag(kvar),jdiagaz(kvar)
      integer noptions(0:AZ_OPTIONS_SIZE)
      double precision params(0:AZ_PARAMS_SIZE)
c
c
6       format (1x,10i5)
7       format (1x,6e12.5)
c
C ...  READ JDIAG FILE


       do i=1,kvar
          jdiag(i)=jdiag1(i)
          jdiagaz(i)=jdiagaz1(i)
       enddo

C ...   TO READ JRAF FILE

c        WRITE(*,*) 'EMASS ======='
c        WRITE(*,7) (A(I),I=1,MAXA)
c        WRITE(*,*) 'F ========'
c        WRITE(*,7) (F(I),I=1,NEQ)
c
c...... READ Na FILE

c
c      WRITE(*,*) 'Einform ======='
c      WRITE(*,7) (NA(I),I=1,MAXA)
c
c...... READ lgNv FILE

c
c      WRITE(*,*) 'local to global nv  ======='
c      WRITE(*,7) ((maplg(I),I=1,NEQ)
c
c...... READ ilgNv FILE

c      print *,'jdiag == ',(jdiag(i),i=1,neq+1)
c      print *,'jdiagaz ==',(jdiagaz(i),i=1,neq+1)
c      print *,'1111111111na == '
c      do i=1,neq
c         nl=jdiagaz(i)
c         nu=jdiagaz(i+1)-1
c         print *,(na(j),j=nl,nu)
c      enddo
c      print *,'======================================'
c
c      WRITE(*,*) 'indicator of local to global nv  ======='
c      WRITE(*,7) ((imaplg(I),I=1,NEQ)
c
      do i=1,maxa
c         val(i)=a(i)
         nbindex(i)=na(i)
      enddo
c
c---------------    dmb2 ===> dmsr   ------------------
c      and delete extra boundary function matrix and information
c
c      DMSR----distributed modified sparse row
c      change the distrubuted gross stiff matrix to DMSR format
c
      k=0
      do i = 1,neq
        if(imaplg(i).eq.1) then
          k=k+1
          diag(k) = val(jdiag(i))
        endif
      end do
      if(k.ne.n_update) print *,'Error! k unequal n_update!!'
c
      nonzero = 0
      ka = 0
      n_update1=n_update+1
      k=0
      do i = 1,neq
        nl = jdiagaz(i)
        nu = jdiagaz(i+1)-1
        if(imaplg(i).eq.1) then
          k=k+1
          jdiagaz(k) = nonzero
          do j = nl,nu
            ka = ka + 1
            if(nbindex(j).eq.i) then
              if(dabs(diag(k)-val(j)).gt.1.0e-35) then
                 print *,'DMSR error'
c                 print *,'i,ka,j,nbindex ==',i,ka,j,nbindex(j)
c                 print *,'diag,val_j==',diag(k),val(j)
              endif
            end if
            if((nbindex(j).gt.0).and.(nbindex(j).ne.i)) then
              nonzero = nonzero + 1
              val(nonzero) = val(j)
c            nbindex(nonzero) = nbindex(ka)
              nbindex(nonzero) = maplg(nbindex(j))-1
            end if
          end do
c        else
c          ka=ka+nu-nl
        endif

      end do
      if(k.ne.n_update) print *,"error! k unequal n_update", k
      jdiagaz(n_update1)  = nonzero
c
c      print *,'nonzero == ',nonzero
      do i = nonzero,1,-1
        nbindex(i+n_update1) = nbindex(i)
        val(i+n_update1) = val(i)
      end do
c
c      nbindex(1)=n_update1
      do i = 1,n_update1
        nbindex(i) = jdiagaz(i) + n_update1
      end do
c
      do i = 1,n_update
        val(i) = diag(i)
      end do
      val(n_update1) = 0.0d0
      maxa = nbindex(n_update1)
c
c======     check the matrix information     =============
c      print *,'all of nbindex ===='
c      print *,(nbindex(i),i=1,nbindex(n_update1))
c      print *,'iblk,n_update,nbindex_diag == ',iblk,n_update
c      print *,(nbindex(i),i=1,n_update1)
c      print *,'nbindex ====='
c      do i=1,n_update
c         nl=nbindex(i)+1
c         nu=nbindex(i+1)
c         print *,(nbindex(j),j=nl,nu)
c      enddo
c      print *,'val ====='
c      do i=1,n_update
c         nl=nbindex(i)+1
c         nu=nbindex(i+1)
c         print *,(val(j),j=nl,nu)
c      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(maxa.gt.maxb) then
        print *,'Fatal error!!! maxa is greater than maxb!!!'
        print *,'maxa,maxb ===',maxa,maxb
        stop 1000
      endif

cc
cc  change local number to global number of NBINDEX array
cc
cc
c        do i = 1,n_update1
c        jdiagaz(i) = nbindex(i)
c        end do
cc
c        do i = 1,n_update1
c        diag(i) = val(i)
c        end do
cc
c        nnew = 0
c        k = 0
c        nbindex(1) = N_update1
c        do i = 1,n_update
c          nl = jdiagaz(i)
c          nu = jdiagaz(i+1)
c          nadd = 0
c          do j = nl,nu-1
c            nnew = nnew + 1
c            nadd = nadd + 1
c            nbindex(nnew+N_update+1) = maplg(nbindex(j))-1
c            val(nnew+N_update+1) = val(j)
c          end do
c          k = k + 1
c          nbindex(k+1) = nbindex(k)+nadd
c          val(k) = diag(i)
c        end do
c        val(k+1) = 0.0d0
c
c        maxa = nnew+N_update1
cc
cc---------------  dmbs to dmsr finished  -----------------
cc
c
c=========================================================
c       print out the dmsr matrix for error check
c
c      write(*,*) 'neq, new maxa=====',neq,maxa
c      write(*,*) 'new jdiag file ====='
c      write(*,*) (JDIAG(I),I=1,NEQ)
c      write(*,*) 'new jdiagaz file ====='
c      write(*,*) (JDIAGAZ(I),I=1,NEQ1)
c      write(*,*) 'new na file ====='
c      write(*,*) (nbindex(I),I=1,maxa)
c      write(*,*) 'new a file ====='
c      write(*,*) (val(I),I=1,maxa)
c      write(*,*) 'new f file ====='
c      write(*,*) (b(I),I=1,NEQ)
c      do i = 1,neq
c        nl = nbindex(i)
c        nu = nbindex(i+1)-1
c        write(*,*) i,(nbindex(j),j=nl,nu)
c        write(*,*) val(i),(val(j),j=nl,nu)
c      end do
c
c       get the right term of the algebra equation
c
      k=0
      do i = 1,neq
        if(imaplg(i).eq.1) then
          k=k+1
          b(k)=f(i)
        endif
      end do
c      write(*,*) 'new f file ====='
c      write(*,*) (b(I),I=1,NEQ)
c
c      give the initial of solution
c
      k=0
      do i = 1,neq
        if(imaplg(i).eq.1) then
          k=k+1
          x(k)=0.d0
        endif
      end do
c
c     get nupdate array for calling parallel solver
c
      k=0
      do i = 1,neq
        if(imaplg(i).eq.1) then
          k=k+1
          nupdate(k) = maplg(i)-1
        endif
      end do
c
c         data check!!!
c
        do i = 1,maxa
        if(nbindex(i).lt.0) then
          print *,'Negative row number, fatal error!!!',iblk
          return
        end if
        end do
c
c       call parallel solver and got the solution.
c
      maxa0 = maxa - 1
      n_update0=n_update-1
      n_external0=n_external-1
      n_dataorg0 = n_dataorg - 1
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
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c      noptions(AZ_max_iter)=100000
c      noptions(AZ_solver) = AZ_bicgstab
c      noptions(AZ_precond) = AZ_sym_GS
c      noptions(AZ_poly_ord) = 3
c      noptions(AZ_scaling) = AZ_sym_diag
      include "solverlap.opt"
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      call callpsolver(maxa0,n_update0,n_external0,n_dataorg0,
     &          n_update,nupdate,nbindex,val,b,x,ndata_org,nexternal,
     &          nupdate_indx,nextern_indx,noptions,params)
c
c       write out the solution and return
c       The solutions in iblk subdomain are saved to vector x and b.
c
c      do i = n_update1,neq
c        b(i) = 0.0d0
c      end do
c      print *,'iblk,x == ',iblk
c      print *,(x(i),i=1,neq)
c      print *,(b(i),i=1,neq)
      do i=1,neq
         x(i)=0.0d0
      enddo
      k = 0
      do i = 1,neq
        if(imaplg(i).eq.1) then
          k = k + 1
          x(i) = b(k)
        endif
      end do
c      print *,'sol=========='
c      print *,(x(i),i=1,n_update)
c
c      write(25) (b(i),i=1,neq)
c      write(25) (x(i),i=1,N_update)
c
c      do i = 1,neq
c        x(i) = 0.0d0
c      end do
c
c       send the solution at update points to server
c
c      call sendint(nsource,iblk,neq)
c      call sendai(nsource,iblk,maplg,neq)
c      print *,'iblk,n_update===',iblk,n_update
c      call sendint(nsource,iblk,n_update)
c      call sendar(nsource,iblk,b,n_update)
c      do k = 1,neq
c        x(k) = b(k)
c      end do
c      do k = 1,neq
c        b(k) = 0.0d0
c      end do
c
c      recieve the solution at external points from server
c
c      call sendint(nsource,iblk,neq)
c      call sendai(nsource,iblk,maplg,neq)
c      call recvint(iblk,nsource,n_external)
c      call recvar(iblk,nsource,b,n_external)
c
c      k = 0
c      do i = 1,neq
c        if(imaplg(i).eq.0) then
c          k = k + 1
c          x(i) = b(k)
c        endif
c      end do
c      if(k.ne.n_external) print *,"error! k unequal n_external ",k
c
c.......write out lsol file
c
c
      maxa = maxb - 10
c
1000  format(i5,10f25.16)
      return
      end

