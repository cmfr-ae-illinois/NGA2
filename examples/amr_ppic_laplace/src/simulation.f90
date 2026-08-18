!> AMR Translating Laplace Equilibrium
!> Periodic in X/Y/Z
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrmpinc_class,    only: amrmpinc
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   use amrio_class,       only: amrio
   use string,            only: str_medium
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! Time integration
   type(timetracker) :: time

   ! Solver data
   type(amrmpinc), target :: fs
   type(amrdata) :: dQdt,Umag

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile,currentsfile
   real(WP) :: vel_rms,curv_std

   ! Restart data
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   ! Physical parameters
   real(WP) :: radius = 0.4_WP
   real(WP), dimension(3) :: drop_stretch,drop_vel
   real(WP) :: viscL_mol,viscG_mol
   real(WP) :: Tnu,Tsigma,TU,tdTnu,tdTsigma,tdTU,Usigma

   real(WP), parameter, public :: VFlo=1.0e-12_WP    ! Minimum VF value considered
   real(WP), parameter, public :: VFhi=1.0_WP-VFlo   ! Maximum VF value considered

contains

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! G=radius-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
      G=1.0_WP-sqrt(xyz(1)**2/(drop_stretch(1)*radius)**2+xyz(2)**2/(drop_stretch(2)*radius)**2+xyz(3)**2/(drop_stretch(2)*radius)**2)
      ! if (amr%nz.eq.1) G=radius-sqrt(xyz(1)**2+xyz(2)**2) ! Enable 2D case
      if (amr%nz.eq.1) G=1.0_WP-sqrt(xyz(1)**2/(drop_stretch(1)*radius)**2+xyz(2)**2/(drop_stretch(2)*radius)**2) ! Enable 2D case
   end function sphere_levelset

   !> Compute viscosity
   subroutine get_viscosity()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVisc
      real(WP), parameter :: myeps=1.0e-15_WP
      ! Loop over levels
      do lvl=0,amr%clvl()
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Use harmonic averaging
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(viscL_mol,myeps)+(1.0_WP-pVF(i,j,k,1))/max(viscG_mol,myeps))
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosity

   !> User-provided initialization for jet
   subroutine drop_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpinc_class, only: VFlo
      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF=>null(),pCL=>null(),pCG=>null(),pU=>null(),pV=>null(),pW=>null(),pQ=>null()
      real(WP), dimension(3) :: BL,BG  ! Dummy barycenters
      real(WP) :: dx,dy,dz,VF
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         pU=>solver%U%mf(lvl)%dataptr(mfi)
         pV=>solver%V%mf(lvl)%dataptr(mfi)
         pW=>solver%W%mf(lvl)%dataptr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Set initial velocity
            pU(i,j,k,1)=drop_vel(1); pV(i,j,k,1)=drop_vel(2); pW(i,j,k,1)=drop_vel(3)
            pQ(i,j,k,1)=drop_vel(1); pQ(i,j,k,2)=drop_vel(2); pQ(i,j,k,3)=drop_vel(3)
            ! Set interface
            ! Compute VF and barycenters from levelset
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=nref,VFlo=VFlo,VF=VF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=VF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine drop_init

   !> Compute turbulence stats
   subroutine compute_stats()
      use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy,amrex_box,amrex_mfiter
      use amrex_interface, only: amrmask_make_fine
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
      integer :: i,j,k,ierr,lvl
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_imultifab) :: mask
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pVF,pCurv
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask
      real(WP) :: curv_mean, curv_mean_weight, vol_total
      ! Update nondimensional times
      tdTnu=time%t/Tnu
      tdTsigma=time%t/Tsigma
      tdTU=time%t/TU
      ! Uses composite integration with fine masking to avoid double-counting
      vel_rms  =0.0_WP
      curv_std =0.0_WP
      vol_total=0.0_WP
      curv_mean=0.0_WP
      curv_mean_weight=0.0_WP
      do lvl=0,amr%clvl()
         ! Build fine mask for this level (if not finest)
         if (lvl.lt.amr%clvl()) then
            call amrex_imultifab_build(mask,amr%ba(lvl),amr%dm(lvl),1,0)
            call amrmask_make_fine(mask,amr%ba(lvl+1),[amr%rrefx(lvl),amr%rrefy(lvl),amr%rrefz(lvl)],0,1)
         end if
         ! Loop over all cells
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pU =>fs%U%mf(lvl)%dataptr(mfi)
            pV =>fs%V%mf(lvl)%dataptr(mfi)
            pW =>fs%W%mf(lvl)%dataptr(mfi)
            if (lvl.lt.amr%clvl()) pMask=>mask%dataptr(mfi)
            if (lvl.eq.amr%maxlvl) then
               pCurv=>fs%curv%dataptr(mfi)
               pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
            end if
            ! Loop over tile
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Skip cells covered by finer level
               if (lvl.lt.amr%clvl()) then
                  if (pMask(i,j,k,1).eq.0) cycle
               end if
               ! Accumulate vel rms and curv mean
               vel_rms=vel_rms+((pU(i,j,k,1)-drop_vel(1))**2+(pV(i,j,k,1)-drop_vel(2))**2+(pW(i,j,k,1)-drop_vel(3))**2)*amr%cell_vol(lvl)
               vol_total=vol_total+amr%cell_vol(lvl)
               if (lvl.eq.amr%maxlvl) then
                  if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) cycle 
                  curv_mean=curv_mean+pCurv(i,j,k,1)
                  curv_mean_weight=curv_mean_weight+1.0_WP
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
         if (lvl.lt.amr%clvl()) call amrex_imultifab_destroy(mask)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,vel_rms,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vol_total,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,curv_mean,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,curv_mean_weight,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      vel_rms=sqrt(vel_rms/vol_total)/Usigma
      curv_mean=curv_mean/curv_mean_weight
      ! Loop over all cells
      lvl=amr%maxlvl
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pCurv=>fs%curv%dataptr(mfi)
         pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Accumulate curv std
            if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) cycle 
            curv_std=curv_std+(pCurv(i,j,k,1)-curv_mean)**2
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE,curv_std,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      curv_std=sqrt(curv_std/curv_mean_weight)/curv_mean
   end subroutine compute_stats
   
   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none

      ! Create amrgrid
      create_amrgrid: block
         amr%name='laplace'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-10.0_WP; amr%xhi=+10.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.true.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         call param_read('Blocking factor',amr%nbloc)
         call param_read('Max grid size',amr%nmax)
         ! Handle 2D case
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves here
      handle_restart: block
         integer :: restart_step
         ! Initialize IO object
         call io%initialize(amr=amr,nfiles=1)
         ! Check if restarting
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         ! If restarting, read header
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time integration
      initialize_time: block
         ! Create time tracker and initialize
         time=timetracker(amRoot=amr%amRoot,name="Laplace Eq.")
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_time

      ! Create flow solver
      create_flow_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrdata_class,    only: interp_face_lin
         use amrmpinc_class,   only: BC_GAS,BC_USER
         use amrmg_class,      only: amrmg_outer_pcg_mlmg
         ! Create flow solver
         call fs%initialize(amr,name='laplace')
         ! Set initial conditions
         fs%user_init=>drop_init
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set densities
         fs%rhoG=1.0_WP; call param_read('Density ratio',fs%rhoL,default=1.0_WP); fs%rhoL=fs%rhoG*fs%rhoL
         ! Set molecular viscosities
         viscG_mol=1.0_WP; call param_read('Viscosity ratio',viscL_mol,default=1.0_WP); viscL_mol=viscG_mol*viscL_mol
         ! Set surface tension coefficient
         call param_read('Laplace number',fs%sigma); fs%sigma=fs%sigma*viscL_mol**2.0_WP/(fs%rhoL*2.0_WP*radius)
         Tnu=(2.0_WP*radius)**2.0_WP*fs%rhoL/viscL_mol
         Tsigma=sqrt(fs%rhoL*(2.0_WP*radius)**3/fs%sigma)
         Usigma=sqrt(fs%sigma/(fs%rhoL*(2.0_WP*radius)))
         ! Initial sphere stretching
         call param_read('Droplet stretching',drop_stretch,default=[1.0_WP,1.0_WP,1.0_WP])
         call param_read('Droplet velocity',drop_vel,default=[0.0_WP,0.0_WP,0.0_WP])
         TU=2.0_WP*radius/sqrt(dot_product(drop_vel,drop_vel))
         ! Set pressure convergence
         fs%psolver%outer_solver=amrmg_outer_pcg_mlmg
         fs%psolver%tol_rel=1.0e-5_WP
      end block create_flow_solver

      ! Create workspace array
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=3,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         ! fs%user_tagging=>my_tagger
         ! call param_read('Tagging Reynolds',Re_tag)
         ! Create initial grid from scratch or restore from checkpoint
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            ! Create initial grid
            call amr%init_from_scratch(time=time%t)
            ! Build PLIC
#ifdef USE_IRL
            call fs%build_ppic(time%t)
#else
            call fs%build_plic(time%t)
#endif
            ! call fs%build_plic(time%t)
            call fs%build_subVF()
            ! Initialize face velocities
            call fs%get_face_velocity()
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         end if
         ! Set viscosity: molecular + SGS
         call get_viscosity()
         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_visualization: block
         ! Create visualization object
         call viz%initialize(amr,'laplace',use_hdf5=.false.)
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(fs%Q,1,'U')
         call viz%add_scalar(fs%Q,2,'V')
         call viz%add_scalar(fs%Q,3,'W')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(fs%VF,1,'VF')
         call fs%smesh%write_as_vtu()
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_visualization

      ! Create monitor
      create_monitor: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(time%dt,time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(fs%CFL,'CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%psolver%res,'Pressure residual')
         call mfile%add_column(fs%psolver%niter,'Pressure iterations')
         call mfile%add_column(fs%divmax,'Divergence')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLst,'CFLst')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
         ! Create parasitic currents monitor
         currentsfile=monitor(amRoot=amr%amRoot,name='statistics')
         call currentsfile%add_column(time%n,'Timestep')
         call currentsfile%add_column(time%t,'Time')
         call currentsfile%add_column(tdTnu,'Time/T_nu')
         call currentsfile%add_column(tdTsigma,'Time/T_sigma')
         call currentsfile%add_column(tdTU,'Time/T_U')
         call currentsfile%add_column(curv_std,'NonDim Curv STD')
         call currentsfile%add_column(vel_rms, 'NonDim Vel RMS')
         call compute_stats()
         call currentsfile%write()
      end block create_monitor

   end subroutine simulation_init

   !> Run the simulation
   subroutine simulation_run()
      implicit none

      ! Time integration loop
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Store old interface and velocities
         call fs%store_old()

         ! Prepare for surface output
         if (viz_evt%occurs()) fs%update_smesh=.true.

         ! Sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity: U^{mid} = 0.5*(U + Uold)
            call fs%Q%lincomb(a=0.5_WP,src1=fs%Qold,b=0.5_WP,src2=fs%Q)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Increment velocity with advection+viscous terms
            call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
            call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
            call fs%Q%average_down(); call fs%Q%fill(time%t)

            ! Rebuild PLIC and sub-cell VF
#ifdef USE_IRL
            call fs%build_ppic(time%t)
#else
            call fs%build_plic(time%t)
#endif
            ! call fs%build_plic(time%t)
            call fs%build_subVF()

            ! Interpolate velocity to the faces
            call fs%get_face_velocity()

            ! Increment both velocities with current pressure term
            call fs%add_pressure(scale=time%dt,phi=fs%P)

            ! Add surface tension to both velocities
            call fs%add_surface_tension(scale=time%dt)

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Correct outflow for mass conservation
            ! call fs%correct_outflow()

            ! Prepare and solve pressure Poisson
            call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
            call fs%prepare_psolver()
            call fs%psolver%solve(rhs=fs%div)

            ! Correct both velocities with pressure increment
            call fs%add_pressure(scale=time%dt)

            ! Add pressure increment
            call fs%P%add(src=fs%psolver%sol)

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Update viscosity
         call get_viscosity()

         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)

         ! Monitor output
         call fs%get_info()
         call mfile%write()
         call cflfile%write()
         call compute_stats()
         call currentsfile%write()

         ! Visualization output
         if (viz_evt%occurs()) then
            call viz%write(time=time%t)
            fs%update_smesh=.false.
         end if

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/laplace_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if
         
      end do

   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
      call currentsfile%finalize()
   end subroutine simulation_final

end module simulation
