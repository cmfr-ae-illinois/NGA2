!> AMR integral surface tension class
module amrist_class
    use precision,        only: WP 
    use mathtools,        only: Pi
    use amrvof_class,     only: VFlo,VFhi,vol_eps,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER
    use amrmg_class,      only: amrmg
    use amrmpflow_class,  only: amrmpflow
    use amrmpinc_class,    only: amrmpinc
    use amrdata_class,    only: amrdata
    use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter
    use string,            only: str_medium
#ifdef USE_IRL
    use amrpic_class,      only: picmfab,pic_container
#endif
    use irl_fortran_interface
    implicit none
    private

    ! Expose Type
    public :: amrist

    !> AMR Integral Surface Tension Calculator Type
    type :: amrist
        ! Name 
        character(len=str_medium) :: name

        ! Flow Solver
        class(amrmpinc), pointer :: fsvf

        ! Options
        real(WP), dimension(3) :: MarangoniOption
        integer :: SurfaceTensionOption,CurvatureOption,SmoothingOption
        logical :: TwoD
        real(WP) :: PressureOption
        real(WP) :: PU_spread

        ! Stresses
        ! Stresses in x direction, then y direction, then z direction. 
        ! Face dependence is indexed on each stress.
        type(amrdata) :: ST_x_stresses,ST_y_stresses,ST_z_stresses
        type(amrdata) :: ST_x_force,ST_y_force,ST_z_force
        ! Flag to skip registration with amrgrid in case of inheritance
        logical :: skip_registration=.false.
    contains 
        procedure :: temp 
        ! Constructor/destructor
        procedure :: initialize ! Done
        procedure :: finalize ! Done

        ! Lifecycle callbacks
        procedure :: on_init
        procedure :: on_coarse
        procedure :: on_remake
        procedure :: on_clear

        ! Physics 
        procedure :: get_dQdt ! Done
        procedure :: add_surface_tension ! Done 
        ! Visualization

        ! ================== Private 
        ! Momentum Functions

        ! Surface Tension Functions
        procedure, private :: add_integral_surface_tension ! Done
        procedure, private :: add_csf_shift_integral_surface_tension_jump 

        ! Integral Surface Tension Methods
        procedure, private :: update_surface_tension_stresses !Done 
        procedure, private :: update_surface_tension_forces ! Done

        ! Ellipsoid Exact Values
        procedure, private :: update_surface_tension_stresses_ellipsoid
        procedure, private :: update_surface_tension_forces_ellipsoid

        ! Smoothing Functions

        ! Helpers 
        procedure, private :: add_interface_to_PU_neighborhood ! Done
    end type amrist 
contains
    subroutine temp(this)
        class(amrist), intent(inout) :: this
        print *, "This is a temporary subroutine in the amrist module."
    end subroutine temp

    ! ============================================================================
    ! INITIALIZATION / FINALIZATION
    ! ============================================================================

    subroutine initialize(this,fsvf,name)
        use amrdata_class, only: interp_const
        implicit none 
        class(amrist), target, intent(inout) :: this
        class(amrmpinc), target, intent(in) :: fsvf 
        character(len=*), intent(in), optional :: name
        integer :: lvl
        ! IRL REQUIRED
#ifdef USE_IRL
            print *, "Creating IST"
#else
            print *, "ERROR: amrist require IRL usage"
            STOP 
#endif
        if(present(name)) then 
            this%name = trim(name) 
        else 
            this%name ='UNNAMED_IST'
        endif
        ! Store flow solver
        ! Check if 3 or more ghost layers
        this%fsvf=>fsvf
        if(fsvf%nover .lt. 3) then  
            print *, "WARNING: amrinc object has less than 3 ghost layers. 3 required for AMRIST initializiation"
        endif
        ! Initialize AMRData Variables
        print *, "INITIALIZING STRESSES"
        call this%ST_x_stresses%initialize(this%fsvf%amr,name='ST_x_Stresses',ncomp=3,ng=this%fsvf%nover)
        call this%ST_y_stresses%initialize(this%fsvf%amr,name='ST_y_Stresses',ncomp=3,ng=this%fsvf%nover)
        call this%ST_z_stresses%initialize(this%fsvf%amr,name='ST_z_Stresses',ncomp=3,ng=this%fsvf%nover)

        call this%ST_x_force%initialize(this%fsvf%amr,name='ST_x_Force',ncomp=1,ng=this%fsvf%nover)
        call this%ST_y_force%initialize(this%fsvf%amr,name='ST_y_Force',ncomp=1,ng=this%fsvf%nover)
        call this%ST_z_force%initialize(this%fsvf%amr,name='ST_z_Force',ncomp=1,ng=this%fsvf%nover)
        print *, "ST_x initialized TEST TEST TEST"

        if(.not. this%skip_registration) then 
            select type(this) 
            type is(amrist)
                call this%fsvf%amr%add_on_init   (amrist_on_init,   c_loc(this))
                call this%fsvf%amr%add_on_coarse (amrist_on_coarse, c_loc(this))
                call this%fsvf%amr%add_on_remake (amrist_on_remake, c_loc(this))
                call this%fsvf%amr%add_on_clear  (amrist_on_clear,  c_loc(this))
            end select 
        endif
    end subroutine initialize 

    subroutine finalize(this)
        implicit none
        class(amrist), intent(inout) :: this

        ! Finalize AMRData variables
    end subroutine finalize 

    ! ============================================================================
    ! DISPATCHERS (module-level) - recover concrete amrist type
    ! ============================================================================

    !> Dispatch on_init: calls type-bound method then user callback
    subroutine amrist_on_init(ctx,lvl,time,ba,dm)
        implicit none
        type(c_ptr), intent(in) :: ctx
        integer, intent(in) :: lvl
        real(WP), intent(in) :: time
        type(amrex_boxarray), intent(in) :: ba
        type(amrex_distromap), intent(in) :: dm
        type(amrist), pointer :: this
        call c_f_pointer(ctx,this)
        call this%on_init(lvl,time,ba,dm)
    end subroutine amrist_on_init

    subroutine on_init(this,lvl,time,ba,dm) 
        implicit none
        class(amrist), intent(inout) :: this
        integer, intent(in) :: lvl
        real(WP), intent(in) :: time
        type(amrex_boxarray), intent(in) :: ba
        type(amrex_distromap), intent(in) :: dm

        ! Reset Stresses level layout and zero
        call this%ST_x_stresses%reset_level(lvl,ba,dm); call this%St_x_stresses%setval(val=0.0_WP,lvl=lvl)
        call this%ST_y_stresses%reset_level(lvl,ba,dm); call this%St_y_stresses%setval(val=0.0_WP,lvl=lvl)
        call this%ST_z_stresses%reset_level(lvl,ba,dm); call this%St_z_stresses%setval(val=0.0_WP,lvl=lvl)

        call this%ST_x_force%reset_level(lvl,ba,dm); call this%St_x_force%setval(val=0.0_WP,lvl=lvl)
        call this%ST_y_force%reset_level(lvl,ba,dm); call this%St_y_force%setval(val=0.0_WP,lvl=lvl)
        call this%ST_z_force%reset_level(lvl,ba,dm); call this%St_z_force%setval(val=0.0_WP,lvl=lvl)
    end subroutine on_init

    !> Dispatch on_coarse: calls type-bound method
    subroutine amrist_on_coarse(ctx,lvl,time,ba,dm)
        implicit none
        type(c_ptr), intent(in) :: ctx
        integer, intent(in) :: lvl
        real(WP), intent(in) :: time
        type(amrex_boxarray), intent(in) :: ba
        type(amrex_distromap), intent(in) :: dm
        type(amrist), pointer :: this
        call c_f_pointer(ctx,this)
        call this%on_coarse(lvl,time,ba,dm)
    end subroutine amrist_on_coarse

    subroutine on_coarse(this,lvl,time,ba,dm)
        implicit none
        class(amrist), intent(inout) :: this
        integer, intent(in) :: lvl
        real(WP), intent(in) :: time
        type(amrex_boxarray), intent(in) :: ba
        type(amrex_distromap), intent(in) :: dm

        ! Reset Stresses level layout and zero
        call this%ST_x_stresses%reset_level(lvl,ba,dm); call this%St_x_stresses%setval(val=0.0_WP,lvl=lvl)
        call this%ST_y_stresses%reset_level(lvl,ba,dm); call this%St_y_stresses%setval(val=0.0_WP,lvl=lvl)
        call this%ST_z_stresses%reset_level(lvl,ba,dm); call this%St_z_stresses%setval(val=0.0_WP,lvl=lvl)

        call this%ST_x_force%reset_level(lvl,ba,dm); call this%St_x_force%setval(val=0.0_WP,lvl=lvl)
        call this%ST_y_force%reset_level(lvl,ba,dm); call this%St_y_force%setval(val=0.0_WP,lvl=lvl)
        call this%ST_z_force%reset_level(lvl,ba,dm); call this%St_z_force%setval(val=0.0_WP,lvl=lvl)
    end subroutine on_coarse
    !> Dispatch on_remake: calls type-bound method
    subroutine amrist_on_remake(ctx,lvl,time,ba,dm)
        implicit none
        type(c_ptr), intent(in) :: ctx
        integer, intent(in) :: lvl
        real(WP), intent(in) :: time
        type(amrex_boxarray), intent(in) :: ba
        type(amrex_distromap), intent(in) :: dm
        type(amrist), pointer :: this
        call c_f_pointer(ctx,this)
        call this%on_remake(lvl,time,ba,dm)
    end subroutine amrist_on_remake

    subroutine on_remake(this,lvl,time,ba,dm)
        implicit none
        class(amrist), intent(inout) :: this
        integer, intent(in) :: lvl
        real(WP), intent(in) :: time
        type(amrex_boxarray), intent(in) :: ba
        type(amrex_distromap), intent(in) :: dm

        ! Reset Stresses level layout and zero
        call this%ST_x_stresses%reset_level(lvl,ba,dm); call this%St_x_stresses%setval(val=0.0_WP,lvl=lvl)
        call this%ST_y_stresses%reset_level(lvl,ba,dm); call this%St_y_stresses%setval(val=0.0_WP,lvl=lvl)
        call this%ST_z_stresses%reset_level(lvl,ba,dm); call this%St_z_stresses%setval(val=0.0_WP,lvl=lvl)

        call this%ST_x_force%reset_level(lvl,ba,dm); call this%St_x_force%setval(val=0.0_WP,lvl=lvl)
        call this%ST_y_force%reset_level(lvl,ba,dm); call this%St_y_force%setval(val=0.0_WP,lvl=lvl)
        call this%ST_z_force%reset_level(lvl,ba,dm); call this%St_z_force%setval(val=0.0_WP,lvl=lvl)
    end subroutine on_remake
    !> Dispatch on_clear: calls type-bound method
    subroutine amrist_on_clear(ctx,lvl)
        implicit none
        type(c_ptr), intent(in) :: ctx
        integer, intent(in) :: lvl
        type(amrist), pointer :: this
        call c_f_pointer(ctx,this)
        call this%on_clear(lvl)
    end subroutine amrist_on_clear

    subroutine on_clear(this,lvl)
        implicit none
        class(amrist), intent(inout) :: this
        integer, intent(in) :: lvl

        ! Reset Stresses level layout and zero
        call this%ST_x_stresses%clear_level(lvl)
        call this%ST_y_stresses%clear_level(lvl)
        call this%ST_z_stresses%clear_level(lvl)

        call this%ST_x_force%clear_level(lvl)
        call this%ST_y_force%clear_level(lvl)
        call this%ST_z_force%clear_level(lvl)

    end subroutine on_clear


    ! ============================================================================
    ! UTILITIES
    ! ============================================================================

    subroutine get_dQdt(this,dQdt,dt,time)
        class(amrist), intent(inout) :: this
        class(amrdata), intent(inout) :: dQdt                        ! Output: momentum RHS (cell-centered)
        real(WP), intent(in) :: dt,time
        call this%fsvf%get_dQdt(dQdt,dt,time)
    end subroutine get_dQdt 

    subroutine add_surface_tension(this,scale)
        class(amrist), intent(inout) :: this
        real(WP), intent(in) :: scale
        ! print *,"START add_surface_tension"
        SELECT CASE (this%SurfaceTensionOption)
            case (2)
                call this%add_integral_surface_tension(scale)
            case (3)
                call this%add_csf_shift_integral_surface_tension_jump(scale)
            CASE DEFAULT
                call this%fsvf%add_surface_tension(scale)
        END SELECT
        ! print *,"END add_surface_tension"
    end subroutine add_surface_tension

    subroutine add_interface_to_PU_neighborhood(this,neighborhood,i,j,k,mfi,solver)
        use f_PUNeigh_RectCub_class
        use f_SeparatorVariant_class
        use f_SeparatorUnion_class
        use amrpic_class, only: is_one_plane,is_full,is_empty,set_to_one_plane,is_paraboloid,get_plane,cut_rectcub_pic 
        use irl_fortran_interface
        class(amrist), intent(inout) :: this
        integer :: i,j,k,lvl
        type(PUNeigh_RectCub_type) :: neighborhood
        type(PUST_RectCub_type), optional :: solver
        real(WP) :: problo(3)
        real(WP) :: dx,dy,dz
        real(WP), dimension(1:3) :: cenInitial,projectedNormal,cen,planeNormal
        real(WP) :: weight_vf, weight_normal, alignment 
        type(amrex_mfiter) :: mfi 
        real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPLIC
        type(pic_container), dimension(:,:,:,:), contiguous, pointer :: pPIC
        type(pic_container) :: interface_component 
        type(SeparatorVariant_type) :: interface_component_variant
        type(RectCub_type) :: cell
        type(PlanarSep_type) :: planar_sep
        type(Poly_type) :: polygon
        real(WP), dimension(3) :: lo,hi

        ! print *, "START ADD INTERFACE"
        ! Solver Check 
        lvl    = this%fsvf%amr%maxlvl 
        dx = this%fsvf%amr%dx(lvl)
        dy = this%fsvf%amr%dy(lvl)
        dz = this%fsvf%amr%dz(lvl)
        

        ! print *, "GOT DX"
        call new(cell)
        call new(polygon)
        call new(planar_sep)
        ! print *, "NEWS"
        if(present(solver)) then 
            ! Compute Centroid of cell
            cenInitial = (/this%fsvf%amr%xlo+real(i  ,WP)*dx,this%fsvf%amr%ylo+real(j  ,WP)*dy,this%fsvf%amr%zlo+real(k  ,WP)*dz/)
            cenInitial = cenInitial + (/dx/2.0_WP,dy/2.0_WP,dz/2.0_WP/)

            ! Get PU Normal
            call getNormalPU(solver,cenInitial(1),cenInitial(2),cenInitial(3),this%PU_spread*this%fsvf%amr%dx(lvl),projectedNormal)
        endif
        ! Get VF and PIC Containers as arrays
        ! print *, " Getting Pointers"
        pVF => this%fsvf%VF%mf(lvl)%dataptr(mfi)
        pPIC => this%fsvf%PIC%dataptr(mfi)
        pPLIC=>this%fsvf%PLIC%dataptr(mfi)

        if(pVF(i,j,k,1) .gt. VFlo .and. pVF(i,j,k,1) .lt. VFhi) then ! Mixed Cell 
            ! First, get plane into 
            ! print *, "ADDING INTERFACE AT",i,j,k
            ! print *, "i,j,k =", i,j,k
            ! print *, "PLIC bounds:"
            ! print *, lbound(pPLIC,1), ubound(pPLIC,1)
            ! print *, lbound(pPLIC,2), ubound(pPLIC,2)
            ! print *, lbound(pPLIC,3), ubound(pPLIC,3)
            ! print *, lbound(pPLIC,4), ubound(pPLIC,4)
            if (is_full(pPLIC(i,j,k,:)).or.is_empty(pPLIC(i,j,k,:))) then 
                print *,"NO VALID PLANE"
                return 
            endif

            interface_component = pPIC(i,j,k,1)
            call new(interface_component_variant,interface_component)

            call setNumberOfPlanes(planar_sep,1)
            call setPlane(planar_sep,0,pPLIC(i,j,k,1:3),pPLIC(i,j,k,4))
            lo=[this%fsvf%amr%xlo+real(i  ,WP)*dx,this%fsvf%amr%ylo+real(j  ,WP)*dy,this%fsvf%amr%zlo+real(k  ,WP)*dz]
            hi=[this%fsvf%amr%xlo+real(i+1,WP)*dx,this%fsvf%amr%ylo+real(j+1,WP)*dy,this%fsvf%amr%zlo+real(k+1,WP)*dz]
            call construct_2pt(cell,lo,hi)
            call getPoly(cell,planar_sep,0,polygon)
            ! cen = calculateCentroid(polygon)
            cen = (/this%fsvf%amr%xlo+real(i  ,WP)*dx,this%fsvf%amr%ylo+real(j  ,WP)*dy,this%fsvf%amr%zlo+real(k  ,WP)*dz/)
            cen = cen + (/dx/2.0_WP,dy/2.0_WP,dz/2.0_WP/)
           
            
            if (getNumberOfVertices(polygon).lt.3) then
                ! print *,"WARNING: DEGENERATE POLYGON IN PU NEIGHBORHOOD" 
                ! print *, "i,j,k =", i,j,k
                ! print *, "xlo,ylo,zlo =", this%fsvf%amr%xlo, &
                !                         this%fsvf%amr%ylo, &
                !                         this%fsvf%amr%zlo
                ! print *, "dx,dy,dz =", dx,dy,dz
                ! print *, "lvl =", lvl
                ! print *, cen
            endif

           
            

            ! Compute VF Weight
            weight_vf = 1.0_WP 
            if(pVF(i,j,k,1) .lt. 0.1) then 
                weight_vf = 0.5_WP - 0.5_WP * COS(10.0_WP*Pi*pVF(i,j,k,1))
            endif 

            if(pVF(i,j,k,1) .gt. 0.9) then 
                weight_vf = 0.5_WP - 0.5_WP * COS(10.0_WP*Pi*(1.0_WP-pVF(i,j,k,1)))
            endif 

            ! Compute Normal Weight
            weight_normal = 1.0_WP
            if(present(solver)) then 
                if (getNumberOfVertices(polygon).lt.3) then
                    ! print *,"WARNING: DEGENERATE POLYGON IN PU NEIGHBORHOOD"
                endif
                alignment = sum(projectedNormal*planeNormal)
                weight_normal = max(alignment,0.0_WP)
            endif

            ! Add Member 
            call addMember(neighborhood,cen,weight_vf*weight_normal,interface_component_variant,0.0_WP)
        endif
    end subroutine add_interface_to_PU_neighborhood 

    ! ============================================================================
    ! PHYSICS METHODS
    ! ============================================================================

    subroutine add_integral_surface_tension(this,scale)
        use amrex_amr_module, only: amrex_multifab
        use amrex_interface, only: amrmfab_average_down_face
        implicit none
        class(amrist), intent(inout) :: this
        real(WP), intent(in) :: scale 

        type(amrex_multifab),dimension(:), allocatable :: STFx,STFy,STFz
        type(amrex_mfiter) :: mfi 
        type(amrex_box) :: bx
        real(WP), dimension(:,:,:,:), contiguous, pointer :: pSigma_x,pSigma_y,pSigma_z
        real(WP) :: dxi,dyi,dzi,VF_f,mysurf,mycurv
        integer :: lvl,i,j,k
        
        ! print *, "START add_integral_surface_tension"
        ! Guard: no surface tension or clvl<maxlvl
        if (this%fsvf%sigma.eq.0.0_WP.or.this%fsvf%amr%clvl().lt.this%fsvf%amr%maxlvl) return
        ! Build temp face flux mfabs
        allocate(STFx(0:this%fsvf%amr%clvl()),STFy(0:this%fsvf%amr%clvl()),STFz(0:this%fsvf%amr%clvl()))
        do lvl=0,this%fsvf%amr%clvl()
            call this%fsvf%amr%mfab_build(lvl,STFx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.]); call STFx(lvl)%setval(0.0_WP)
            call this%fsvf%amr%mfab_build(lvl,STFy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.]); call STFy(lvl)%setval(0.0_WP)
            call this%fsvf%amr%mfab_build(lvl,STFz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ]); call STFz(lvl)%setval(0.0_WP)
        end do

        ! Here we make the forces (The Hard Part)
        call this%update_surface_tension_forces(STFx,STFy,STFz)
        
        ! Average down face fluxes from finest to coarser levels
        do lvl=this%fsvf%amr%clvl()-1,0,-1
            call amrmfab_average_down_face(fmf=STFx(lvl+1),cmf=STFx(lvl),rr=[this%fsvf%amr%rrefx(lvl),this%fsvf%amr%rrefy(lvl),this%fsvf%amr%rrefz(lvl)],cgeom=this%fsvf%amr%geom(lvl))
            call amrmfab_average_down_face(fmf=STFy(lvl+1),cmf=STFy(lvl),rr=[this%fsvf%amr%rrefx(lvl),this%fsvf%amr%rrefy(lvl),this%fsvf%amr%rrefz(lvl)],cgeom=this%fsvf%amr%geom(lvl))
            call amrmfab_average_down_face(fmf=STFz(lvl+1),cmf=STFz(lvl),rr=[this%fsvf%amr%rrefx(lvl),this%fsvf%amr%rrefy(lvl),this%fsvf%amr%rrefz(lvl)],cgeom=this%fsvf%amr%geom(lvl))
        end do
        ! Apply fluxes
        call this%fsvf%apply_face_fluxes(scale,STFx,STFy,STFz)
        ! Destroy temps
        do lvl=0,this%fsvf%amr%clvl()
            call this%fsvf%amr%mfab_destroy(STFx(lvl))
            call this%fsvf%amr%mfab_destroy(STFy(lvl))
            call this%fsvf%amr%mfab_destroy(STFz(lvl))
        end do
        deallocate(STFx,STFy,STFz)
        ! print *, "END add_integral_surface_tension"
    end subroutine add_integral_surface_tension

    subroutine add_csf_shift_integral_surface_tension_jump(this,scale)
        use amrex_amr_module, only: amrex_multifab
        use amrex_interface, only: amrmfab_average_down_face
        implicit none 
        class(amrist), intent(inout) :: this
        type(amrex_multifab),dimension(:), allocatable :: STFx,STFy,STFz
        real(WP) :: scale
        print *, "This is the add_csf_shift_integral_surface_tension_jump subroutine in the amrist module."
    end subroutine add_csf_shift_integral_surface_tension_jump 




    subroutine update_surface_tension_stresses(this)
        use irl_fortran_interface
        use f_PUNeigh_RectCub_class
        class(amrist), intent(inout) :: this
        type(PUST_RectCub_type) :: solver
        real(WP) :: dx,dy,dz
        integer :: lvl,i,j,k 
        integer :: i_in,j_in, k_in 
        integer :: shift_first_index,shift_second_index
        real(WP), dimension(1:3) :: dvec,shift
        real(WP), dimension(1:3) :: pressure_cell_center,velocity_cell_center,face_center 
        real(WP), dimension(1:3) :: force
        real(WP), dimension(1:3) :: P0,P1,P2,P3 
        type(PUNeigh_RectCub_type) :: neighborhood,neighborhood_trimmed
        type(amrex_box) :: bx
        type(amrex_mfiter) :: mfi
        real(WP), dimension(:,:,:,:), contiguous, pointer :: pSigma_x,pSigma_y,pSigma_z

        ! print *, "START update_surface_tension_stresses"
        lvl = this%fsvf%amr%maxlvl 
        dx = this%fsvf%amr%dx(lvl)
        dy = this%fsvf%amr%dy(lvl)
        dz = this%fsvf%amr%dz(lvl)
        dvec = (/dx,dy,dz/)
        ! Create Neighborhood and Solver
        call new(neighborhood) 
        call new(neighborhood_trimmed) 
        call new(solver)

        ! Iterate Over Domain at Finest lvl 
        call this%fsvf%amr%mfiter_build(lvl,mfi)
        do while (mfi%next())
            ! print *,"GETTING POINTERS"
            pSigma_x => this%ST_x_stresses%mf(lvl)%dataptr(mfi)
            pSigma_y => this%ST_y_stresses%mf(lvl)%dataptr(mfi)
            pSigma_z => this%ST_z_stresses%mf(lvl)%dataptr(mfi)
            bx = mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)

                ! Here we make the Neighborhood 
                ! print *,"Neighborhood1"
                call emptyNeighborhood(neighborhood)
                do i_in = -3,3; do j_in = -3,3; do k_in = -3,3 
                    call this%add_interface_to_PU_neighborhood(neighborhood,i+i_in, j+j_in, k+k_in,mfi)
                end do; end do; end do;
                ! print *,"DONE NEIGHBORHOOD"
                call setNeighborhood(solver,neighborhood)
                call setKernelSize(solver, this%PU_spread*dx) 

                ! Remake Neighborhood with normal weights
                ! call emptyNeighborhood(neighborhood)
                ! do i_in = -3,3; do j_in = -3,3; do k_in = -3,3 
                !     call this%add_interface_to_PU_neighborhood(neighborhood,i+i_in, j+j_in, k+k_in,mfi,solver)
                ! end do; end do; end do;
                ! call setNeighborhood(solver,neighborhood)
                ! call setKernelSize(solver, this%PU_spread*dx) 

                ! Get Stresses 
                pressure_cell_center = (/this%fsvf%amr%xlo+real(i  ,WP)*dx,this%fsvf%amr%ylo+real(j  ,WP)*dy,this%fsvf%amr%zlo+real(k  ,WP)*dz/)
                pressure_cell_center = pressure_cell_center + dvec/2
                do i_in = 1,3 
                    ! Set up shift to velocity cell
                    shift = (/0.0_WP,0.0_WP,0.0_WP/)  
                    shift(i_in) = dvec(i_in)/2
                    ! Move center to velocity cell we are looking at
                    velocity_cell_center = pressure_cell_center - shift
                    do j_in = 1,3 
                        force = (/0.0_WP,0.0_WP,0.0_WP/)
                        ! Shift from velocity cell center ot the face center we are looking at:
                        shift = (/0.0_WP,0.0_WP,0.0_WP/)  
                        shift(j_in) = dvec(j_in)/2
                        face_center = velocity_cell_center + shift
                        ! Now that we have the face center and normal direction (all positive normals)
                        ! we can use the pattern below to get a CCW orientation:
                        ! -- =>  +- => ++ => -+
                        ! The directions we apply this two are the directions perpendicular to our normal, in CCW order. This is to say if our directional order is:
                        ! (x,y,z) then we use  (y,z) for x normal, (z,x) for y normal, and (x,y) for z normal
                        ! This is where we define those directions, in indices
                        shift_first_index = mod(j_in,3)+1
                        shift_second_index = mod(j_in+1,3)+1
                        ! Apply to get first corner
                        shift = dvec/2     
                        shift(j_in) = 0.0_WP
                        shift(shift_first_index) = -shift(shift_first_index)
                        shift(shift_second_index) = -shift(shift_second_index)
                        P0 = face_center + shift

                        ! Apply to get first corner
                        shift = dvec/2     
                        shift(j_in) = 0.0_WP
                        shift(shift_first_index) = shift(shift_first_index)
                        shift(shift_second_index) = -shift(shift_second_index)
                        P1 = face_center + shift

                        ! Apply to get first corner
                        shift = dvec/2     
                        shift(j_in) = 0.0_WP
                        shift(shift_first_index) = shift(shift_first_index)
                        shift(shift_second_index) = shift(shift_second_index)
                        P2 = face_center + shift

                        ! Apply to get first corner
                        shift = dvec/2     
                        shift(j_in) = 0.0_WP
                        shift(shift_first_index) = -shift(shift_first_index)
                        shift(shift_second_index) = shift(shift_second_index)
                        P3 = face_center + shift

                        ! Now that we have the face corners, run the code and get the force
                        ! if(this%vf%VF(i,j,k) .gt. 1e-12 .and. this%vf%VF(i,j,k) .lt. 1.0_WP - 1e-12) then 
                
                        call solveFace(solver,this%fsvf%sigma,P0,P1,P2,P3,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                        
                        ! Now that we have the force, store it properly 
                        if(i_in .eq. 1) then 
                            pSigma_x(i,j,k,j_in) = force(i_in) 
                        elseif(i_in .eq. 2) then 
                            pSigma_y(i,j,k,j_in) = force(i_in) 
                        elseif (i_in .eq. 3) then 
                            pSigma_z(i,j,k,j_in) = force(i_in) 
                        endif 
                        ! if(sum(force**2).gt. 0.0_WP) then 
                        !     print *,"FORCE: ",force
                        ! endif
                        
                        ! print *, "==================================="
                        ! print *, "i_in,j_in: ",i_in,j_in
                        ! print *, "pressure_cell_center: ", pressure_cell_center
                        ! print *, "velocity_cell_center: ", velocity_cell_center
                        ! print *, "face_center: ", face_center 
                        ! print *, "shift_first_index: ", shift_first_index 
                        ! print *, "shift_second_index: ", shift_second_index 
                        ! print *, "P0: ", P0
                        ! print *, "P1: ", P1
                        ! print *, "P2: ", P2
                        ! print *, "P3: ", P3
                        ! write(*,'(A)', advance='no') "For Desmos: [("
                        ! write(*,'(F10.5, A, F10.5, A)', advance='no') P0(shift_first_index), ",", P0(shift_second_index), "),("
                        ! write(*,'(F10.5, A, F10.5, A)', advance='no') P1(shift_first_index), ",", P1(shift_second_index), "),("
                        ! write(*,'(F10.5, A, F10.5, A)', advance='no') P2(shift_first_index), ",", P2(shift_second_index), "),("
                        ! write(*,'(F10.5, A, F10.5, A)', advance='yes') P3(shift_first_index), ",", P3(shift_second_index), ")]"
                        ! print *, "==================================="
                    enddo
                enddo
                ! print *, "OUT STRESS LOOP"
            end do; end do; end do
            ! print *,"OUT BOX LOOP"
        enddo 
        ! print *, "OUT MFI LOOP"
    end subroutine update_surface_tension_stresses 
    

    subroutine update_surface_tension_forces(this,STFx,STFy,STFz)
        use amrex_amr_module, only: amrex_multifab
        use amrex_interface, only: amrmfab_average_down_face
        implicit none 
        class(amrist), intent(inout) :: this
        type(amrex_multifab),dimension(:), allocatable, intent(inout) :: STFx,STFy,STFz
        type(amrex_box) :: bx
        type(amrex_mfiter) :: mfi
        integer :: lvl,i,j,k
        real(WP) :: dxi,dyi,dzi
        real(WP), dimension(:,:,:,:), contiguous, pointer :: pSigma_x,pSigma_y,pSigma_z
        real(WP), dimension(:,:,:,:), contiguous, pointer :: pSTFx,pSTFy,pSTFz
        real(WP), dimension(:,:,:,:), contiguous, pointer :: pSTFx_viz,pSTFy_viz,pSTFz_viz

        ! print *, "START update_surface_tension_forces"
        call this%update_surface_tension_stresses()
        
        lvl = this%fsvf%amr%maxlvl 
        dxi=1.0_WP/this%fsvf%amr%dx(lvl); dyi=1.0_WP/this%fsvf%amr%dy(lvl); dzi=1.0_WP/this%fsvf%amr%dz(lvl)
        call this%fsvf%amr%mfiter_build(lvl,mfi)
        do while (mfi%next())
            pSTFx =>STFx(lvl)%dataptr(mfi)
            pSTFy =>STFy(lvl)%dataptr(mfi)
            pSTFz =>STFz(lvl)%dataptr(mfi)

            pSTFx_viz =>this%ST_x_force%mf(lvl)%dataptr(mfi)
            pSTFy_viz =>this%ST_y_force%mf(lvl)%dataptr(mfi)
            pSTFz_viz =>this%ST_z_force%mf(lvl)%dataptr(mfi)

            pSigma_x => this%ST_x_stresses%mf(lvl)%dataptr(mfi)
            pSigma_y => this%ST_y_stresses%mf(lvl)%dataptr(mfi)
            pSigma_z => this%ST_z_stresses%mf(lvl)%dataptr(mfi)

            ! stresses are stored at cell centers
            ! X Faces
            bx=mfi%nodaltilebox(1) 
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                pSTFx(i,j,k,1)= (pSigma_x(i,j,k,1)-pSigma_x(i-1,j,k,1)) * dyi*dzi + &
                                (pSigma_x(i,j,k,2)-pSigma_x(i,j-1,k,2)) * dxi*dzi
                if(.not. this%TwoD) then 
                    pSTFx(i,j,k,1) = pSTFx(i,j,k,1) + &
                                (pSigma_x(i,j,k,3)-pSigma_x(i,j,k-1,3)) * dxi*dyi
                endif
                pSTFx_viz(i,j,k,1) = pSTFx(i,j,k,1)


            end do; end do; end do

            ! Y Faces
            bx=mfi%nodaltilebox(2) 
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                pSTFy(i,j,k,1)= (pSigma_y(i,j,k,1)-pSigma_y(i-1,j,k,1)) * dyi*dzi + &
                                (pSigma_y(i,j,k,2)-pSigma_y(i,j-1,k,2)) * dxi*dzi
                if(.not. this%TwoD) then 
                    pSTFy(i,j,k,1) = pSTFy(i,j,k,1) + &
                                (pSigma_y(i,j,k,3)-pSigma_y(i,j,k-1,3)) * dxi*dyi
                endif

                pSTFy_viz(i,j,k,1) = pSTFy(i,j,k,1)
            end do; end do; end do

            ! Z Faces
            bx=mfi%nodaltilebox(3) 
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                pSTFz(i,j,k,1)= (pSigma_z(i,j,k,1)-pSigma_z(i-1,j,k,1)) * dyi*dzi + &
                                (pSigma_z(i,j,k,2)-pSigma_z(i,j-1,k,2)) * dxi*dzi
                if(.not. this%TwoD) then 
                    pSTFz(i,j,k,1) = pSTFz(i,j,k,1) + &
                                (pSigma_z(i,j,k,3)-pSigma_z(i,j,k-1,3)) * dxi*dyi
                endif

                pSTFz_viz(i,j,k,1) = pSTFz(i,j,k,1)
            end do; end do; end do
        enddo 
        call this%fsvf%amr%mfiter_destroy(mfi)  
        ! print *, "END update_surface_tension_forces"       
    end subroutine update_surface_tension_forces


    ! ============================================================================
    ! ELLIPSOID METHODS
    ! ============================================================================
    subroutine update_surface_tension_stresses_ellipsoid(this)
        class(amrist), intent(inout) :: this
        print *, "This is the update_surface_tension_stresses_ellipsoid subroutine in the amrist module."
    end subroutine update_surface_tension_stresses_ellipsoid 

    subroutine update_surface_tension_forces_ellipsoid(this)
        class(amrist), intent(inout) :: this
        print *, "This is the update_surface_tension_forces_ellipsoid subroutine in the amrist module."
    end subroutine update_surface_tension_forces_ellipsoid 

    

end module amrist_class