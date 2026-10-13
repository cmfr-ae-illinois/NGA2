!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use mathtools,         only: Pi
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use vtk_class,         only: vtk
   use partmesh_class,    only: partmesh
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use tracer_class,      only: tracer
   use irl_fortran_interface
   implicit none
   private

   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs           !< DDADI solver for velocity
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   type(tracer),      public :: pt

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(partmesh) :: pmesh
   type(vtk)      :: vtk_out
   type(event)    :: vtk_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,curvfile
   character(len = 1000) :: MonitorFileName

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi

   !> Problem definition
   real(WP), dimension(3) :: center,vel
   real(WP) :: radius
   real(WP) :: pert_amp    !< perturbation amplitude
   integer  :: pert_mode
   integer  :: npart

   !> output monitor variables
   real(WP) :: r_max, r_min, time_cap, time_nd

contains

   subroutine compute_rp_radii(state)
      class(vfs), intent(in) :: state
      call compute_rp_max_radius(state)
      call compute_rp_min_radius(state)
   end subroutine compute_rp_radii

   subroutine compute_rp_max_radius(this)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_MAX, MPI_SUM, MPI_DOUBLE_PRECISION
      use mathtools, only: Pi
      implicit none

      class(vfs), intent(in) :: this
      integer :: i,j,k,ierr,nk
      real(WP) :: radius_k
      real(WP) :: max_radius_local, max_radius_global
      real(WP), allocatable, dimension(:) :: slice_volume_local, slice_volume_global

      allocate(slice_volume_local(this%cfg%kmin:this%cfg%kmax))
      allocate(slice_volume_global(this%cfg%kmin:this%cfg%kmax))

      slice_volume_local  = 0.0_WP
      slice_volume_global = 0.0_WP

      ! Accumulate local contribution to each global z-slice
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               slice_volume_local(k) = slice_volume_local(k) + &
                  this%VF(i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do

      nk = this%cfg%kmax - this%cfg%kmin + 1

      ! Sum slice volumes entrywise across all ranks
      call MPI_ALLREDUCE(slice_volume_local, slice_volume_global, nk, &
         MPI_DOUBLE_PRECISION, MPI_SUM, this%cfg%comm, ierr)

      max_radius_local = 0.0_WP

      do k=this%cfg%kmin,this%cfg%kmax
         if (slice_volume_global(k) .gt. 0.0_WP) then
            radius_k = sqrt(slice_volume_global(k)/(Pi*this%cfg%dz(k)))
            if (radius_k .gt. max_radius_local) max_radius_local = radius_k
         end if
      end do

      call MPI_ALLREDUCE(max_radius_local, max_radius_global, 1, &
         MPI_DOUBLE_PRECISION, MPI_MAX, this%cfg%comm, ierr)

      r_max = max_radius_global

      deallocate(slice_volume_local, slice_volume_global)
   end subroutine compute_rp_max_radius

   subroutine compute_rp_min_radius(this)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_MIN, MPI_SUM, MPI_DOUBLE_PRECISION
      use mathtools, only: Pi
      implicit none

      class(vfs), intent(in) :: this
      integer :: i,j,k,ierr,nk
      real(WP) :: radius_k
      real(WP) :: min_radius_local, min_radius_global
      real(WP), allocatable :: slice_volume_local(:), slice_volume_global(:)

      allocate(slice_volume_local(this%cfg%kmin:this%cfg%kmax))
      allocate(slice_volume_global(this%cfg%kmin:this%cfg%kmax))

      slice_volume_local  = 0.0_WP
      slice_volume_global = 0.0_WP

      ! Accumulate local contribution to each global z-slice
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               slice_volume_local(k) = slice_volume_local(k) + &
                  this%VF(i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do

      nk = this%cfg%kmax - this%cfg%kmin + 1

      ! Sum slice volumes entrywise across all ranks
      call MPI_ALLREDUCE(slice_volume_local, slice_volume_global, nk, &
         MPI_DOUBLE_PRECISION, MPI_SUM, this%cfg%comm, ierr)

      min_radius_local = huge(1.0_WP)

      do k=this%cfg%kmin,this%cfg%kmax
         if (slice_volume_global(k) .gt. 0.0_WP) then
            radius_k = sqrt(slice_volume_global(k)/(Pi*this%cfg%dz(k)))
            if (radius_k .lt. min_radius_local) min_radius_local = radius_k
         end if
      end do

      call MPI_ALLREDUCE(min_radius_local, min_radius_global, 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, this%cfg%comm, ierr)

      r_min = min_radius_global

      deallocate(slice_volume_local, slice_volume_global)
   end subroutine compute_rp_min_radius


   !> modified bessel functions of first kind
   function bessel_i0(x) result(val)
      implicit none
      real(WP), intent(in) :: x
      real(WP) :: val
      real(WP) :: ax, y

      ax = abs(x)

      if (ax < 3.75_WP) then
         y = (x/3.75_WP)**2
         val = 1.0_WP + y*(3.5156229_WP + y*(3.0899424_WP + y*(1.2067492_WP + &
               y*(0.2659732_WP + y*(0.0360768_WP + y*0.0045813_WP)))))
      else
         y = 3.75_WP/ax
         val = (exp(ax)/sqrt(ax)) * (0.39894228_WP + y*(0.01328592_WP + &
               y*(0.00225319_WP + y*(-0.00157565_WP + y*(0.00916281_WP + &
               y*(-0.02057706_WP + y*(0.02635537_WP + y*(-0.01647633_WP + &
               y*0.00392377_WP))))))))
      end if
   end function bessel_i0

   function bessel_i1(x) result(val)
      implicit none
      real(WP), intent(in) :: x
      real(WP) :: val
      real(WP) :: ax, y

      ax = abs(x)

      if (ax < 3.75_WP) then
         y = (x/3.75_WP)**2
         val = x*(0.5_WP + y*(0.87890594_WP + y*(0.51498869_WP + &
               y*(0.15084934_WP + y*(0.02658733_WP + y*(0.00301532_WP + &
               y*0.00032411_WP))))))
      else
         y = 3.75_WP/ax
         val = (exp(ax)/sqrt(ax)) * (0.39894228_WP + y*(-0.03988024_WP + &
               y*(-0.00362018_WP + y*(0.00163801_WP + y*(-0.01031555_WP + &
               y*(0.02282967_WP + y*(-0.02895312_WP + y*(0.01787654_WP - &
               y*0.00420059_WP))))))))
         if (x < 0.0_WP) val = -val
      end if
   end function bessel_i1

   !> Level set function for sinusoidally perturbed ligament along z-axis
   function levelset_ligament(xyz,t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP) :: xrel,yrel,zrel,r,rloc,kz

      xrel = xyz(1) - center(1)
      yrel = xyz(2) - center(2)
      zrel = xyz(3) - center(3)

      r = sqrt(xrel*xrel + yrel*yrel)

      ! Axial perturbation
      kz   = 2.0_WP*Pi*real(pert_mode,WP)/(cfg%zL)
      ! kz   = 2.0_WP*Pi
      rloc = radius*(1.0_WP + pert_amp*cos(kz*zrel))

      ! Positive inside liquid
      G = rloc - r
   end function levelset_ligament


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class,only: VFhi,VFlo,lvira,elvira,mof,wmof,r2p,youngs,lvlset,plicnet,r2pnet, &
                            jibben,jibbenPU,jibbenSq,PUppic,PUplic,jibbenPLIC,cf,pucf
         use messager, only: die
         use mpi_f08,  only: MPI_WTIME
         use string,   only: str_medium,lowercase
         integer :: i,j,k,n,si,sj,sk,curvature_method,stencil_size,hf_backup_method
         integer :: reconstruction_method
         character(len=str_medium) :: read_reconstruction_method
         character(len=str_medium) :: read_curvature_method
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=5
         real(WP) :: start, finish
         ! Select interface reconstruction from the input file
         call param_read('Reconstruction method',read_reconstruction_method,default='plicnet')
         select case (trim(lowercase(read_reconstruction_method)))
         case ('lvira'); reconstruction_method=lvira
         case ('elvira'); reconstruction_method=elvira
         case ('mof'); reconstruction_method=mof
         case ('wmof'); reconstruction_method=wmof
         case ('r2p'); reconstruction_method=r2p
         case ('youngs'); reconstruction_method=youngs
         case ('lvlset'); reconstruction_method=lvlset
         case ('plicnet'); reconstruction_method=plicnet
         case ('r2pnet'); reconstruction_method=r2pnet
         case ('jibben'); reconstruction_method=jibben
         case ('jibbenpu'); reconstruction_method=jibbenPU
         case ('jibbensq'); reconstruction_method=jibbenSq
         case ('puppic'); reconstruction_method=PUppic
         case ('puplic'); reconstruction_method=PUplic
         case ('jibbenplic'); reconstruction_method=jibbenPLIC
         case ('cf'); reconstruction_method=cf
         case ('pucf'); reconstruction_method=pucf
         case default
            call die('[rayleigh_plateau] Unknown reconstruction method: '//trim(read_reconstruction_method))
         end select
         call vf%initialize(cfg=cfg,reconstruction_method=reconstruction_method,name='VOF')
         ! Reading ligament parameters
         call param_read('Ligament radius',radius)
         call param_read('Ligament center',center,default=[0.0_WP,0.0_WP,0.0_WP])
         call param_read('Perturbation amplitude',pert_amp)
         call param_read('Perturbation mode',pert_mode,default=1)      
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_ligament,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Perform PPIC reconstruction
         if (vf%ppic) call vf%build_quadratic_interface()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use mathtools, only: Pi
         integer :: i,j,k
         real(WP), dimension(3) :: xyz
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         !Capilary time scale
         time_cap = sqrt(fs%rho_l*radius**3/fs%sigma)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initial velocity field
         call param_read('Initial velocity',vel,default=[0.0_WP,0.0_WP,0.0_WP])
         fs%U=vel(1)
         fs%V=vel(2)
         fs%W=vel(3)
         initialize_rp_velocity: block
            real(WP) :: kz, kr0, I0kr0, I1kr0, I0kr, I1kr
            real(WP), parameter :: tiny_r = 1.0e-14_WP
            real(WP) :: xrel, yrel, zrel, r
            real(WP) :: B, c
            real(WP), dimension(3) :: xyz

            ! wave number
            kz = 2.0_WP*Pi*real(pert_mode,WP)/cfg%zL
            kr0   = kz*radius
            I0kr0 = bessel_i0(kr0)
            I1kr0 = bessel_i1(kr0)

            ! Linear inviscid growth-rate parameter
            c = sqrt((fs%sigma/(fs%rho_l*radius**3)) * (I1kr0/I0kr0) * kr0 * (1.0_WP - kr0**2))
            !c  = sqrt(max(c,0.0_WP))

            ! amplitude of velocity potential
            B = pert_amp*radius*c/(kz*I1kr0)

            ! U velocity at x-faces
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     xyz = [fs%cfg%x(i), fs%cfg%ym(j), fs%cfg%zm(k)]

                     ! if (levelset_ligament(xyz,0.0_WP) >= 0.0_WP) then
                        xrel = xyz(1) - center(1)
                        yrel = xyz(2) - center(2)
                        zrel = xyz(3) - center(3)
                        r    = sqrt(xrel*xrel + yrel*yrel)

                        if (r > tiny_r) then
                           I1kr = bessel_i1(kz*r)
                           fs%U(i,j,k) = B*kz*I1kr*(xrel/r)*cos(kz*zrel)
                        else
                           fs%U(i,j,k) = 0.0_WP
                        ! end if
                     end if
                  end do
               end do
            end do

            ! V velocity at y-faces
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     xyz = [fs%cfg%xm(i), fs%cfg%y(j), fs%cfg%zm(k)]

                     ! if (levelset_ligament(xyz,0.0_WP) >= 0.0_WP) then
                        xrel = xyz(1) - center(1)
                        yrel = xyz(2) - center(2)
                        zrel = xyz(3) - center(3)
                        r    = sqrt(xrel*xrel + yrel*yrel)

                        if (r > tiny_r) then
                           I1kr = bessel_i1(kz*r)
                           fs%V(i,j,k) = B*kz*I1kr*(yrel/r)*cos(kz*zrel)
                        else
                           fs%V(i,j,k) = 0.0_WP
                        ! end if
                     end if
                  end do
               end do
            end do

            ! W velocity at z-faces
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     xyz = [fs%cfg%xm(i), fs%cfg%ym(j), fs%cfg%z(k)]

                     ! if (levelset_ligament(xyz,0.0_WP) >= 0.0_WP) then
                        xrel = xyz(1) - center(1)
                        yrel = xyz(2) - center(2)
                        zrel = xyz(3) - center(3)
                        r    = sqrt(xrel*xrel + yrel*yrel)

                        I0kr = bessel_i0(kz*r)
                        fs%W(i,j,k) = -B*kz*I0kr*sin(kz*zrel)
                     ! end if
                  end do
               end do
            end do

         end block initialize_rp_velocity
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=2,name='plic')
         smesh%varname(1)='curv'
         smesh%varname(2)='interface_type'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         call add_surfgrid_variable(smesh,1,vf%curv)  
         call add_surfgrid_variable(smesh,1,vf%interface_type)    
      end block create_smesh

      ! Initialize our tracer particle tracker
      initialize_pt: block
         ! Get number of particles
         call param_read('Number of tracers',npart,default=0)
         ! Create solver
         pt=tracer(cfg=cfg,name='PT')
         ! Initialize with zero particles
         call pt%resize(0)
      end block initialize_pt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=0,nvec=2,name='tracers')
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='acceleration'
         call pt%update_partmesh(pmesh)
         do i=1,pt%np_
            pmesh%vec(:,1,i)=pt%p(i)%vel
            pmesh%vec(:,2,i)=pt%p(i)%acc
         end do
      end block create_pmesh

      ! Seed interface with tracers
      insert_tracers: block
         use mpi_f08, only: MPI_INTEGER8,MPI_MAX
         integer :: i,np0_,ierr
         integer(kind=8) :: maxid_,maxid  !< Keep track of maximum particle id
         real :: hk,thetak,phik,realn
         ! Initial number of particles
         np0_=pt%np_
         ! Determine id to assign to particle
         maxid_=0
         do i=1,pt%np_
            maxid_=max(maxid_,pt%p(i)%id)
         end do
         call MPI_ALLREDUCE(maxid_,maxid,1,MPI_INTEGER8,MPI_MAX,cfg%comm,ierr)

         ! Add new particles
         if (cfg%amRoot) then
            ! Create space for new particle
            call pt%resize(np0_+npart)
            ! Initialize parameters
            phik=0.0_WP;realn=REAL(npart,WP)
            do i=1,npart
               ! Helicoidal seeding based on Saff and Kuijlaars (1997)
               hk=-1.0_WP+2.0_WP*(REAL(i,WP)-1.0_WP)/(realn-1.0_WP); thetak=acos(hk)
               if (i.eq.1.or.i.eq.npart) then
                  phik=0.0_WP
               else
                  phik=phik+3.6_WP/sqrt(realn)/sqrt(1.0_WP-hk**2.0_WP)
               end if
               ! Set particle ID
               pt%p(np0_+i)%id=maxid+int(i,8)
               ! Seed particle on sphere
               pt%p(np0_+i)%pos(1)=center(1)+radius*sin(thetak)*cos(phik)
               pt%p(np0_+i)%pos(2)=center(2)+radius*sin(thetak)*sin(phik)
               pt%p(np0_+i)%pos(3)=center(3)+radius*cos(thetak)
               ! Localize the particle
               pt%p(np0_+i)%ind=cfg%get_ijk_global(pt%p(np0_+i)%pos,[cfg%imin,cfg%jmin,cfg%kmin])
               ! Make it an "official" particle
               pt%p(np0_+i)%flag=0
            end do
         end if
         ! Communicate particles
         call pt%sync()
         call pt%update_partmesh(pmesh)
         do i=1,pt%np_
            pmesh%vec(:,1,i)=pt%p(i)%vel
            pmesh%vec(:,2,i)=pt%p(i)%acc
         end do
      end block insert_tracers

      ! Add Ensight output
      create_vtk: block
         ! Create Ensight output from cfg
         vtk_out=vtk(cfg=cfg,name='RayleighPlateau')
         ! Create event for Ensight output
         vtk_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',vtk_evt%tper)
         ! Add variables to output
         call vtk_out%add_vector('velocity',Ui,Vi,Wi)
         call vtk_out%add_scalar('VOF',vf%VF)
         call vtk_out%add_scalar('curvature',vf%curv)
         call vtk_out%add_surface('plic',smesh) 
         call vtk_out%add_particle('tracers',pmesh)
         ! Output to vtk
         if (vtk_evt%occurs()) call vtk_out%write_data(time%t)
      end block create_vtk

      create_monitor: block
         ! Create simulation monitor
         call param_read('Monitor File Name',MonitorFileName)
         ! Initial nondimensional monitor values
         time_nd = 0.0_WP
         r_max = 0.0_WP
         r_min = 0.0_WP
         mfile=monitor(fs%cfg%amRoot,MonitorFileName)
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time_nd,'t/tc')
         call mfile%add_column(r_max,'Rmax')
         call mfile%add_column(r_min,'Rmin')

         call mfile%write()         
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old VOF
         vf%VFold=vf%VF

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced

         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,diagnostics_callback=compute_rp_radii)

         ! monitor variables for diagonistics
         time_nd = time%t / time_cap
         call mfile%write()

         ! Advance and project tracer particles
         call pt%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho_U,P=fs%P)
         call project_tracers(vf=vf,pt=pt)

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)

            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Output to vtk
         if (vtk_evt%occurs()) then
            ! Update surfmesh object
            update_smesh: block
               use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               call add_surfgrid_variable(smesh,1,vf%curv)   
               call add_surfgrid_variable(smesh,2,vf%interface_type)   
            end block update_smesh  
            ! Update partmesh object
            update_pmesh: block
              integer :: i
              call pt%update_partmesh(pmesh)
              do i=1,pt%np_
                  pmesh%vec(:,1,i)=pt%p(i)%vel
                  pmesh%vec(:,2,i)=pt%p(i)%acc
              end do
            end block update_pmesh
            call vtk_out%write_data(time%t)
         end if

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! call print_curvature_error()
      ! Get rid of all objects - need destructors
      ! monitor
      ! vtk
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)

   end subroutine simulation_final


   !> Make a surface scalar variable
   subroutine add_surfgrid_variable(smesh,var_index,A)
      use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
      use vfs_class,only: VFhi,VFlo
      implicit none
      class(surfmesh), intent(inout) :: smesh
      integer, intent(in) :: var_index
      real(WP), dimension(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_), intent(in) :: A 
      integer :: i,j,k,n,shape,np,nplane,nbt

      ! Fill out arrays
      if ((smesh%nPoly+smesh%nBezierTri).gt.0) then
         np=0; nbt = 0
         ! Start with quadratic surfaces
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  shape=getNumberOfTriangles(vf%interface_mixed_surface(i,j,k))
                  if (shape.gt.0) then
                     do n=1,shape
                        nbt=nbt+1
                        smesh%var(var_index,nbt)=A(i,j,k)
                     end do
                  end if
               end do
            end do
         end do
         ! Then do planes
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(vf%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        ! Set nplane variable
                        smesh%var(var_index,nbt+np)=A(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      else
         smesh%var(var_index,1)=1
      end if      

   end subroutine add_surfgrid_variable

   ! Project tracers on the interface
   subroutine project_tracers(vf,pt)
      use mathtools, only: normalize
      implicit none
      class(vfs), intent(inout) :: vf
      class(tracer), intent(inout) :: pt
      integer :: i, j, maxit
      real(WP) :: distx, disty, distz, dist, maxdist
      real(WP), dimension(3) :: oldvel, oldpos, normal
      ! Maximum projection iterations
      maxit=10
      ! Maximum distance from the interface wanted
      maxdist=1.0E-9_WP*vf%cfg%meshsize(0,0,0)
      ! Get gradient of distance function to the interface
      call fs%get_pgrad(vf%G,resU,resV,resW)
      ! Loop over local tracers
      do j=1,maxit
         do i=1,pt%np_
            ! Avoid particles with id=0
            if (pt%p(i)%id.eq.0) cycle
            oldpos=pt%p(i)%pos
            ! Interpolate distance function
            dist=cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=vf%G,bc='n')
            if (abs(dist).gt.maxdist) then
               ! Interpolate interface normal
               normal=normalize([cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resU,bc='n'),&
               &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resV,bc='n'),&
               &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resW,bc='n')])
               ! Move tracer along distance function gradient
               pt%p(i)%pos=pt%p(i)%pos-dist*normal
               ! Correct the position to take into account periodicity
               pt%p(i)%pos(1)=cfg%x(cfg%imin)+modulo(pt%p(i)%pos(1)-cfg%x(cfg%imin),cfg%xL)
               pt%p(i)%pos(2)=cfg%y(cfg%jmin)+modulo(pt%p(i)%pos(2)-cfg%y(cfg%jmin),cfg%yL)
               pt%p(i)%pos(3)=cfg%z(cfg%kmin)+modulo(pt%p(i)%pos(3)-cfg%z(cfg%kmin),cfg%zL)
               ! Relocalize the tracer
               pt%p(i)%ind=cfg%get_ijk_global(pt%p(i)%pos,pt%p(i)%ind)
            end if
         end do
         ! Communicate tracers across procs
         call pt%sync()
      end do
      do i=1,pt%np_
         ! Interpolate fluid quantities to particle location
         pt%p(i)%vel=cfg%get_velocity(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
      end do
   end subroutine project_tracers

end module simulation
 
