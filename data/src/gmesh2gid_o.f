C....
C.... For FEPG 
C.... Change elems0 tp elem0   !! Global Mesh
C.... Change elems1 to elem1   !! Fluid Mesh for Dynamic grid calculation
C....
      implicit real*8 (a-h,o-z)
      integer,ALLOCATABLE::node(:,:),id(:,:),idm(:,:),masterslave(:,:),
     &mapg2l(:)
      real*8, ALLOCATABLE::coor(:,:),disp(:,:),vol(:)
      dimension r(3),dis(3,4)
      parameter(err=1.d-10)

12    format (1X,I10,15E21.10E4)
13    format (1X,I10,28I10)

      open (21,file='gmesh/init.dat',status='old')
      read (21,*) xmax,ymax,zmax,xo,yo,zo,radius,width,thick,radaxis,
     &             angle,radart_mult,ratiol,bladeang
      read (21,*) nx,nz,ns,npt,nct,nsc,naxis,nyi,nyo
      read (21,*) ai,bi,ao,bo
      read (21,*) u_in,visc_f,pe,pv,denf,dens,gravity,omega
      read (21,*) idpoint,timespring
      read (21,*) lshape,icur,amplify
      read (21,*) iarc,arcradout,narcout,arcradin,narcin,mspline
      read (21,*) id4updown
      read (21,*) index4pro
      close(21)
      pe=2.d0*pe*(1.d0+pv)
      pi=datan(1.d0)*4.d0
      angle=pi/180.d0*angle
      open (21,file='control.dat',status='old')
      read (21,*) skip
      read (21,*) skip
      read (21,*) skip
      read (21,*) penalty
      close(21)

C.... 解析整体流固耦合网格坐标数据
      open(21,file='gmesh/gmesh.cor',form='formatted',status='old')
      read(21,*) knode,kcoor
      allocate(coor(kcoor,knode))
      do j=1,knode
         read(21,*) (coor(i,j),i=1,kcoor)
         if (index4pro .gt. 0) then
            tmp=coor(3,j)
            coor(3,j)=coor(1,j)
            coor(1,j)=coor(2,j)
            coor(2,j)=tmp
            if (index4pro.eq.1) then
c
cc Reverse inlet and outlet
c
               coor(2,j)=-coor(2,j)
               coor(3,j)=-coor(3,j)
            endif
         endif
      enddo
      close(21)

C...  输出整体网格坐标数据到coor0
      write(*,'(a,i8)') 'Global Mesh : coor0, knode=', knode
      open(21,file='coor0',form='unformatted',status='unknown')
      write(21) knode,kcoor,((coor(i,j),i=1,kcoor),j=1,knode)
      close(21)

C.... 输出整体网格坐标信息到Gid后处理文件
      open(22,file='fsi.flavia.msh',form='formatted',status='unknown')
      write(22,'(a)') 'Mesh "fsi_w4_s Dimension ',
     +                 '3 Elemtype Tetrahedra Nnode 4'
      write(22,'(a)') 'Coordinates'
      do i=1,knode
        write(22,'(i6,3E12.5)') i,(coor(j,i),j=1,kcoor)
      enddo
      write(22,'(a)') 'End coordinates'

C...  解析整体网格单元数据
      open(21,file='gmesh/gmesh.elm',form='formatted',status='old')
      read(21,*) nelem,nnode
      allocate(node(nnode,nelem))
      do j=1,nelem
         read(21,*) (node(i,j),i=1,nnode)
      enddo
      close(21)

C...  输出整体网格单元数据到elem0
      write(*,'(a,i8)') 'Global Mesh : elem0, nelem =', nelem
      open(21,file='elem0',form='unformatted',status='unknown')
C...  Write w4 mesh data
      write(21) nelem,nnode,((node(i,j),i=1,nnode),j=1,nelem)
C...  Write Material Data
      mtij=5
      mtik=8
      write(21) mtij,mtik,
     &(visc_f,pe,pv,denf,dens,penalty,gravity,angle,mtii=1,mtij)
C.... Write the second Mesh Data
      mtij=0
      mtik=5
      write(21) mtij,mtik
      write(21) mtij,mtik
C.... Write the Second Mate data
      mtij=0
      mtik=1
      write(21) mtij,mtik
      write(21) mtij,mtik
C.... 关闭elem0文件
      close(21)
C.... 输出整体网格elem0到后处理文件
      write(22,'(a)') 'Elements'
      do i=1,nelem
        write(22,'(10I6)') i,(node(j,i),j=1,nnode)
      enddo
      write(22,'(a)') 'End elements'
      close(22)

C.... 解析整体网格的节点ID
      open(21,file='gmesh/gmesh.id',form='formatted',status='old')
      read(21,*) knode,kdgof
      allocate(id(kdgof,knode))
      allocate(idm(kdgof,knode))
c
c id= 1,idp= 1: unknown inside fluid
c id= 2,idp=-1: unknown inside solid
c id=-1,idp= 1: inlet surface
c id= 0,idp= 1: side surface
c id= 2,idp= 1: outlet surface
c id=-1,idp=-1: rotation of axis
c
      do j=1,knode
         read(21,*) (id(i,j),i=1,kdgof)
      enddo
      close(21)
c
      if (index4pro.eq.3) then
         do j=1,nelem
            if (node(nnode,j).eq.3) then
c
cc Fix the fluid nodes on the interace on the inlet as id = 1
c
               do i=1,nnode-1
                  nd=node(i,j)
                  if (id(1,nd).lt.-1 .and. coor(2,nd).lt.err) then
                     md = -id(1,nd)-1
                     do k=1,kdgof
                        id(k,md)=1
                     enddo
                  endif
               enddo
            endif
         enddo
      endif

c Record inlet id=0
      open(21,file='inlet.id',form='unformatted',status='unknown')
      write(21) knode,(id(1,j)+id(kdgof,j),j=1,knode)
      close(21)

C...  输出 inlet.id 到Gid后处理文件以便查看
      open(22,file='fsi.flavia.res',form='formatted',status='unknown')
      write(22,'(a)') 'GID Post Results File 1.0'
      write(22,*) 'Result "inlet" "Load Analysis" 1 Scalar OnNodes'
      write(22,*) 'ComponentNames "id"'
      write(22,*) 'Values'
      do i=1,knode
        write(22,12) i, dble(id(1,i)+id(kdgof,i))
      enddo
      write(22,*) 'End values'

C.... 解析流体网格和固体网格数据，为动网格计算准备数据
      open(31,file='gmesh/fluid.cor',form='formatted',status='old')
      read(31,*) knodef,kcoor
      open(41,file='gmesh/solid.cor',form='formatted',status='old')
      read(41,*) knodes,kcoor

      do j=1,knode
        do i=1,kcoor
          idm(i,j)=id(i,j)
ccc Outlet
          if (id(i,j) .eq. 2 .and. id(kdgof,j) .eq.1 ) id(i,j)=1
ccc Side
          if (id(i,j) .eq. 0) id(i,j)=-1
        enddo
        idm(kdgof,j)=0
        if (id4updown .eq. 0) then
          if (index4pro.eq.0) then
            if (coor(2,j).ge.yo-thick/2.d0-err.and.
     &          coor(2,j).le.yo+thick/2.d0+err.and.
     &          dabs(coor(1,j)-xo).lt.err.and.
     &          dabs(coor(3,j)-zo).lt.err) then
              do i=1,kdgof
                id(i,j)=-1
              enddo
            endif
          elseif (index4pro.eq.1) then
c
cc Fix id for outlet
c
            if (coor(2,j).gt.ymax-1.d-6) then
              radii=dsqrt(coor(1,j)**2+coor(3,j)**2)
c              print *,j,radii,radius,radii-radius
              if (dabs(radii-radius).lt.2.d-3) then
                do i=1,3
                  id(i,j)=-1
                enddo
              endif
            endif
          elseif (index4pro.gt.1) then
            do i=1,kcoor
              if (id(i,j).eq.-1.and.id(kdgof,j).eq.-1) then
                idm(i,j)=1
              endif
            enddo
            if (coor(2,j).gt.ymax-1.d-6) then
              do i=1,kcoor
                idm(i,j)=2
              enddo
            endif
            if (index4pro.eq.2) then
c
cc Fix id for central axis of solid
c
              do i=1,kcoor
                if (id(i,j).eq.-1.and.id(kdgof,j).eq.-1) then
                  id(i,j)=1
                endif
              enddo
              if (coor(2,j) .ge. yo-thick/2.d0-err .and.
     &            coor(2,j) .le. yo+thick/2.d0+err .and.
     &            dabs(coor(1,j)-xo) .lt. err .and.
     &            dabs(coor(3,j)-zo) .lt. err) then
                do i=1,kdgof
                  id(i,j)=-1
                enddo
              endif
            endif
          endif
        else
          if (coor(3,j) .ge. yo-thick/2.d0-err.and.
     &        coor(3,j) .le. yo+thick/2.d0+err.and.
     &        dabs(coor(1,j)-xo) .lt. err .and.
     &        dabs(coor(2,j)-zo) .lt. err) then
            do i=1,kdgof
              id(i,j)=-1
            enddo
          endif
        endif
        if ((j.le.knodef .or. j.gt.knodef+knodes) .and. id(kdgof,j).eq.-1)
     &      id(kdgof,j)=1
      enddo

      if (index4pro .eq. 3) then
        do j=1,nelem
          if (node(nnode,j).eq.4) then
c
cc Fix the inside base surface of the cylindrical axis as id = -1
c
            do i=1,nnode-1
              nd=node(i,j)
              if (id(1,nd).lt.-1) then
                md=-id(1,nd)-1
                do k=1,kcoor
                  id(k,nd)=-1
                  id(k,md)=-1
                enddo
                idm(kdgof,md)=-1
              endif
            enddo
          endif
          if (node(nnode,j).eq.3) then
c
cc Fix the fluid nodes on the interace on the inlet as id = 1
c
            do i=1,nnode-1
              nd=node(i,j)
              if (id(1,nd).lt.-1.and.coor(2,nd).lt.err) then
                md=-id(1,nd)-1
                do k=1,kdgof
                  id(k,md)=1
                  idm(k,md)=1
                enddo
              endif
            enddo
          endif
        enddo
      endif
C.... 整体计算网格ID
      open(21,file='id0',form='unformatted',status='unknown')
      write(21) knode,kdgof,((id(i,j),i=1,kdgof),j=1,knode)
      close(21)
      print *,'29013: id =', (id(i,29013),i=1,kdgof)

      write(22,*) 'Result "id0" "Load Analysis" 1 Vector OnNodes'
      write(22,*) 'ComponentNames "idu" "idv" "idw"'
      write(22,*) 'Values'
      do i=1,knode
        write(22,12) i, (dble(id(j,i)),j=1,3)
      enddo
      write(22,*) 'End values'
      write(22,*) 'Result "idp0" "Load Analysis" 1 Scalar OnNodes'
      write(22,*) 'ComponentNames "idp"'
      write(22,*) 'Values'
      do i=1,knode
        write(22,12) i, dble(id(kdgof,i))
      enddo
      write(22,*) 'End values'

      allocate(disp(kdgof,knode))
      OPEN(21,FILE='time0',FORM='FORMATTED')
      READ(21,*) T0
      close(21)

C...  计算初始条件和边界条件
      do j=1,knode
        do i=1,kdgof
          disp(i,j)=0.d0
        enddo
        r(1)=coor(1,j)
        if (id4updown .eq. 0) then
          if (idm(1,j).eq.-1 .and. id(kdgof,j).eq.1) then
            r(2)=coor(2,j)
            r(3)=coor(3,j)
            if (index4pro .eq. 0) then
              disp(2,j)=bound(r,u_in,xmax,zmax,t0,pi,2)
            else
              disp(2,j)=boundpro(r,u_in,xmax,zmax,t0,pi,2)
            endif
          endif
          if (idpoint.ne.1.and.dabs(omega).gt.err) then
c
cc Dirichlet BC for velocity on the rotation of axis TOP
c
            do i=1,kdgof
              if (id(i,j).ne.-1 .and. idm(kdgof,j).ne.-1) goto 111 
            enddo
            disp(1,j)= omega*(coor(3,j)-zo)
            disp(3,j)=-omega*(coor(1,j)-xo)
111         continue
c
cc Dirichlet BC for velocity on the rotation of axis BOT
c
          endif
        elseif (id4updown.eq.1) then
          if (idm(1,j).eq.-1.and.id(kdgof,j).eq.1) then
            r(2)=coor(3,j)
            r(3)=coor(2,j)
            if (index4pro.eq.0) then
              disp(3,j)=-bound(r,u_in,xmax,zmax,t0,pi,2)
            else
              disp(3,j)=-boundpro(r,u_in,xmax,zmax,t0,pi,2)
            endif
          endif
          if (idpoint.ne.1.and.dabs(omega).gt.err) then
c
cc Dirichlet BC for velocity on the rotation of axis TOP
c
            do i=1,kdgof
              if (id(i,j).ne.-1.and.idm(kdgof,j).ne.-1) goto 112
            enddo
            disp(1,j)=-omega*(coor(2,j)-zo)
            disp(2,j)= omega*(coor(1,j)-xo)
112         continue
c
cc Dirichlet BC for velocity on the rotation of axis BOT
c
          endif
        endif
      enddo
      
      write(*,'(a,i3)') 'kdgof = ', kdgof

C...  整体计算网格初边值条件      
      open(21,file='disp0',form='unformatted',status='unknown')
      write(21) knode,kdgof,((disp(i,j),i=1,kdgof),j=1,knode)
      close(21)

      write(22,*) 'Result "disp0" "Load Analysis" 1 Vector OnNodes'
      write(22,*) 'ComponentNames "u0" "v0" "w0"'
      write(22,*) 'Values'
      do i=1,knode
        write(22,12) i, (dble(disp(j,i)),j=1,3)
      enddo
      write(22,*) 'End values'
      write(22,*) 'Result "dispp0" "Load Analysis" 1 Scalar OnNodes'
      write(22,*) 'ComponentNames "p0"'
      write(22,*) 'Values'
      do i=1,knode
        write(22,12) i, dble(disp(kdgof,i))
      enddo
      write(22,*) 'End values'
      close(22)

      open(21,file='unod',form='unformatted',status='unknown')
      write(21) ((disp(i,j),j=1,knode),i=1,kdgof)
      close(21)

      do j=1,knodef
        read(31,*) (coor(i,j),i=1,kcoor)
        if (index4pro.gt.0) then
          tmp=coor(3,j)
          coor(3,j)=coor(1,j)
          coor(1,j)=coor(2,j)
          coor(2,j)=tmp
          if (index4pro.eq.1) then
c
cc Reverse inlet and outlet
c
            coor(2,j)=-coor(2,j)
            coor(3,j)=-coor(3,j)
          endif
        endif
      enddo
      close(31)

C..... 流体网格，动网格计算
      write(*,'(a,i8)') 'Fluid mesh for Dynamic grid: coor1, knodef=', knodef
      open(21,file='coor1',form='unformatted',status='unknown')
      write(21) knodef,kcoor,((coor(i,j),i=1,kcoor),j=1,knodef)
      close(21)

      open(22,file='fsi_DM.flavia.msh',form='formatted',status='unknown')
      write(22,'(a)') 'Mesh "fsi_w4_1 Dimension 3 ',
     +                'Elemtype Tetrahedra Nnode 4'
      write(22,'(a)') 'Coordinates'
      do i=1,knodef
        write(22,'(i6,3E12.5)') i,(coor(j,i),j=1,kcoor)
      enddo
      write(22,'(a)') 'End coordinates'

      open(21,file='gmesh/fluid.elm',form='formatted',status='old')
      read(21,*) nelemf,nnode
      do j=1,nelemf
         read(21,*) (node(i,j),i=1,nnode)
      enddo
      close(21)
      
C.... 流体网格，动网格计算-laplace      
      write(*,'(a,i8)') 'Fluid mesh for Dynamic grid: elem1, nelemf=', nelemf
c      open(21,file='elem1',form='unformatted',status='unknown')
      open(21,file='elemb0',form='unformatted',status='unknown')
C...  write w4 element data
      write(21) nelemf,nnode,((node(i,j),i=1,nnode),j=1,nelemf)
C.... Write Mate data
      mtij=5
      mtik=8
      write(21) mtij,mtik,
     &(visc_f,pe,pv,denf,dens,penalty,gravity,angle,mtii=1,mtij)
C.... Write the second element Data
      mtij=0
      mtik=5
      write(21) mtij,mtik
      write(21) mtij,mtik
C.... Write the second Mate Data
      mtij=0
      mtik=1
      write(21) mtij,mtik
      write(21) mtij,mtik

      close(21)
C.... 输出elems到后处理文件
      write(22,'(a)') 'Elements'
      do i=1,nelemf
        write(22,'(10I6)') i,(node(j,i),j=1,nnode)
      enddo
      write(22,'(a)') 'End elements'
      close(22)
c
cc Heron-type formula for the volume of a tetrahedron
c
      allocate(vol(nelemf))
      vol_max=err
      vol_min=1.d0/err
      do i=1,nelemf
         do j=1,nnode-2
            jn=node(j,i)
            do k=j+1,nnode-1
               kn=node(k,i)
               dis(j,k)=0.d0
               do l=1,kcoor
                  dis(j,k)=dis(j,k)+(coor(l,jn)-coor(l,kn))**2
               enddo
               dis(j,k)=dsqrt(dis(j,k))
            enddo
         enddo
         xu=(dis(2,4)-dis(1,2)+dis(1,4))*(dis(1,2)+dis(1,4)+dis(2,4))
         xl=(dis(1,2)-dis(1,4)+dis(2,4))*(dis(1,4)-dis(2,4)+dis(1,2))
         yu=(dis(3,4)-dis(2,3)+dis(2,4))*(dis(2,3)+dis(2,4)+dis(3,4))
         yl=(dis(2,3)-dis(2,4)+dis(3,4))*(dis(2,4)-dis(3,4)+dis(2,3))
         zu=(dis(1,4)-dis(1,3)+dis(3,4))*(dis(1,3)+dis(3,4)+dis(1,4))
         zl=(dis(1,3)-dis(3,4)+dis(1,4))*(dis(3,4)-dis(1,4)+dis(1,3))
         a=dsqrt(xl*yu*zu)
         b=dsqrt(yl*zu*xu)
         c=dsqrt(zl*xu*yu)
         d=dsqrt(xl*yl*zl)
         vol(i)=(-a+b+c+d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d)
         denorm=dis(3,4)*dis(1,4)*dis(2,4)*192.d0
         vol(i)=dsqrt(vol(i))/denorm
         vol_max=dmax1(vol_max,vol(i))
         vol_min=dmin1(vol_min,vol(i))
      enddo
      do i=1,nelemf
         vol(i)=(vol_max-vol_min)/vol(i)
      enddo
      open(31,file='fcellarea',form='unformatted',status='unknown')
      write(31) nelemf,(vol(i),i=1,nelemf)
      close(31)

C....  由dispb0(dispm0)产生，此处以没用
c      open(1,file='unodm',form='unformatted',status='unknown')
c      write(1) ((0.d0,j=1,knodef),i=1,kcoor)
c      close(1)
c
      allocate(mapg2l(knode))
      open(1,file='gmesh/gmesh.map',form='formatted',status='old')
      read(1,*) knode
      do j=1,knode
        read(1,*) mapg2l(j)
      enddo
      close(1)

      do j=1,knodes
        read(41,*) (coor(i,j),i=1,kcoor)
        if (index4pro.gt.0) then
          tmp=coor(3,j)
          coor(3,j)=coor(1,j)
          coor(1,j)=coor(2,j)
          coor(2,j)=tmp
          if (index4pro.eq.1) then
c
cc Reverse inlet and outlet
c
             coor(2,j)=-coor(2,j)
             coor(3,j)=-coor(3,j)
           endif
        endif
      enddo
      close(41)
      open(1,file='coor2',form='unformatted',status='unknown')
      write(1) knodes,kcoor,((coor(i,j),i=1,kcoor),j=1,knodes)
      close(1)

      open(21,file='fsi_s.flavia.msh',form='formatted',status='unknown')
      write(21,'(a)') 'Mesh "fsi_w4_2 Dimension 3 ',
     +                'Elemtype Tetrahedra Nnode 4'
      write(21,'(a)') 'Coordinates'
      do i=1,knodes
        write(21,'(i6,3E12.5)') i,(coor(j,i),j=1,kcoor)
      enddo
      write(21,'(a)') 'End coordinates'

      open(1,file='gmesh/solid.elm',form='formatted',status='old')
      read(1,*) nelems,nnode
      do j=1,nelems
         read(1,*) (node(i,j),i=1,nnode)
      enddo
      close(1)

C..... 角动量计算网格
c      open(1,file='elems2',form='unformatted',status='unknown')
      open(1,file='elemc0',form='unformatted',status='unknown')
      write(1) nelems,nnode,((node(i,j),i=1,nnode),j=1,nelems)

      mtij=5
      mtik=8
      write(1) mtij,mtik,
     &(visc_f,pe,pv,denf,dens,penalty,gravity,angle,mtii=1,mtij)

      mtij=0
      mtik=5
      write(1) mtij,mtik
      write(1) mtij,mtik

      mtij=0
      mtik=1
      write(1) mtij,mtik
      write(1) mtij,mtik

      close(1)
C.... 输出elems到后处理文件
      write(21,'(a)') 'Elements'
      do i=1,nelems
        write(21,'(10I6)') i,(node(j,i),j=1,nnode)
      enddo
      write(21,'(a)') 'End elements'
      close(21)

      open (41,file='subdomain.idx',form='unformatted',status='unknown')
      write(41) knode,knodef,knodes,knode-knodef-knodes
      write(41) (mapg2l(i),i=1,knode)
      close(41)
      open (41,file='subdomain',form='unformatted',status='unknown')
      write(41) (mapg2l(i),i=1,knode)
      close(41)

      nms=2
      allocate(masterslave(nms,knodef))

      masterslaveff=0
      do j=1,knodef
         do i=1,kcoor
ccc Outlet
            if (idm(i,j).eq.2 .and. id(kdgof,j).eq.1) idm(i,j)=-1
ccc side
            if (idm(i,j).eq.0) idm(i,j)=-1
            if (idm(kdgof,j).eq.-1) idm(i,j)=-1
         enddo
         if (idm(1,j).lt.-1) then

            masterslaveff=masterslaveff+1
            masterslave(1,masterslaveff)=j
            masterslave(2,masterslaveff)=-idm(1,j)-1

            do i=1,kcoor
               idm(i,j)=-1
            enddo
         endif
      enddo
      do j=1,knodes
        if (id(1,knodef+j).lt.-1) then
          n=-id(1,knodef+j)-1
          m=mapg2l(n)
          do i=1,kcoor
            idm(i,m)=-1
          enddo
        endif
      enddo

      open(1,file='masterslave.fsi',form='unformatted',status='unknown')
      write(1) masterslaveff,nms,((masterslave(i,j),i=1,nms),
     &         j=1,masterslaveff)
      close(1)

C.... 动网格计算数据—laplace
c      open(1,file='idm0',form='unformatted',status='unknown')
      open(1,file='idb0',form='unformatted',status='unknown')
      write(1) knodef,kcoor,((idm(i,j),i=1,kcoor),j=1,knodef)
      close(1)
      open(22,file='fsi_f.flavia.res',form='formatted',status='unknown')
      write(22,'(a)') 'GID Post Results File 1.0'
      write(22,*) 'Result "idm0" "Load Analysis" 1 Vector OnNodes'
      write(22,*) 'ComponentNames "idu0" "idv0" "idw0"'
      write(22,*) 'Values'
      do i=1,knodef
        write(22,12) i, (dble(idm(j,i)),j=1,kcoor)
      enddo
      write(22,*) 'End values'
      close(22)

C.... 动网格计算数据—laplace
c      open(1,file='dispm0',form='unformatted',status='unknown')
      open(1,file='dispb0',form='unformatted',status='unknown')
      write(1) knodef,kcoor,((0.d0,i=1,kcoor),j=1,knodef)
      close(1)

      do j=1,knodes
         do i=1,kcoor
           idm(i,j)=id(i,knodef+j)
         enddo
         m=mapg2l(knodef+j)
         if (id(1,knodef+j).lt.-1) then
           do i=1,kcoor
             idm(i,m)=-1
           enddo
        endif
      enddo

C..  角动量计算数据-c场
c      open(1,file='idm1',form='unformatted',status='unknown')
      open(1,file='idc0',form='unformatted',status='unknown')
      write(1) knodes,kcoor,((idm(i,j),i=1,kcoor),j=1,knodes)
      close(1)
      open(22,file='fsi_s.flavia.res',form='formatted',status='unknown')
      write(22,'(a)') 'GID Post Results File 1.0'
      write(22,*) 'Result "idm1" "Load Analysis" 1 Vector OnNodes'
      write(22,*) 'ComponentNames "idu1" "idv1" "idw1"'
      write(22,*) 'Values'
      do i=1,knodes
        write(22,12) i, (dble(idm(j,i)),j=1,kcoor)
      enddo
      write(22,*) 'End values'
      close(22)

C..  角动量计算数据-c场
c      open(1,file='dispm1',form='unformatted',status='unknown')
      open(1,file='dispc0',form='unformatted',status='unknown')
      write(1) knodes,kcoor,((0.d0,i=1,kcoor),j=1,knodes)
      close(1)

C...  角动量场初值！！！！？？？？应该可由dispc0生成，或者直接读取dispc0
c      open(1,file='unods',form='unformatted',status='unknown')
c      write(1) ((0.d0,j=1,knodes),i=1,kcoor)
c      close(1)

      end
