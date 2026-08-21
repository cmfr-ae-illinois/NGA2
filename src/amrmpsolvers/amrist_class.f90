!> AMR integral surface tension class
module amrist_class
    use precision,        only: WP 
    use mathtools,        only: Pi
    use geometry,         only: cfg 
    use amrvof_class,     only: VFlo,VFhi,vol_eps,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER
    use amrmg_class,      only: amrmg
    use amrmpflow_class,  only: amrmpflow
    use amrdata_class,    only: amrdata
    use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter
#ifdef USE_IRL
    use amrpic_class,      only: picmfab,pic_container
#endif
    implicit none
    private

    ! Expose Type
    public :: amrist

    !> AMR Integral Surface Tension Calculator Type
    type :: amrist

        ! Options
        real(WP), dimension(3) :: MarangoniOption
        integer :: SurfaceTensionOption,CurvatureOption,SmoothingOption
        logical :: TwoD
        real(WP) :: PressureOption

    contains 
        procedure :: temp 
        ! Constructor/destructor
        procedure :: initialize 
        procedure :: finalize 

        ! Physics 
        procedure :: get_dmomdt
        procedure :: add_surface_tension_jump 
        ! Visualization

        ! ================== Private 
        ! Momentum Functions

        ! Surface Tension Functions
        procedure, private :: add_integral_surface_tension_jump 
        procedure, private :: add_csf_shift_integral_surface_tension_jump 

        ! Integral Surface Tension Methods
        procedure, private :: update_surface_tension_stresses
        procedure, private :: update_tension_forces_forces 

        ! Ellipsoid Exact Values
        procedure, private :: update_surface_tension_stresses_ellipsoid
        procedure, private :: update_tension_forces_forces_ellipsoid

        ! Smoothing Functions

        ! Helpers 
        procedure, private :: add_interface_to_PU_neighborhood
    end type amrist 
contains
    subroutine temp(this)
        class(amrist), intent(inout) :: this
        print *, "This is a temporary subroutine in the amrist module."
    end subroutine temp

    subroutine initialize(this)
        class(amrist), intent(inout) :: this
        print *, "This is the initialize subroutine in the amrist module."
    end subroutine initialize 

    subroutine finalize(this)
        class(amrist), intent(inout) :: this
        print *, "This is the finalize subroutine in the amrist module."
    end subroutine finalize 

    subroutine get_dmomdt(this)
        class(amrist), intent(inout) :: this
        print *, "This is the get_dmomdt subroutine in the amrist module."
    end subroutine get_dmomdt 

    subroutine add_surface_tension_jump(this)
        class(amrist), intent(inout) :: this
        print *, "This is the add_surface_tension_jump subroutine in the amrist module."
    end subroutine add_surface_tension_jump 

    subroutine add_integral_surface_tension_jump(this)
        class(amrist), intent(inout) :: this
        print *, "This is the add_integral_surface_tension_jump subroutine in the amrist module."
    end subroutine add_integral_surface_tension_jump 

    subroutine add_csf_shift_integral_surface_tension_jump(this)
        class(amrist), intent(inout) :: this
        print *, "This is the add_csf_shift_integral_surface_tension_jump subroutine in the amrist module."
    end subroutine add_csf_shift_integral_surface_tension_jump 

    subroutine update_surface_tension_stresses(this)
        class(amrist), intent(inout) :: this
        print *, "This is the update_surface_tension_stresses subroutine in the amrist module."
    end subroutine update_surface_tension_stresses 
    
    subroutine update_tension_forces_forces(this)
        class(amrist), intent(inout) :: this
        print *, "This is the update_tension_forces_forces subroutine in the amrist module."
    end subroutine update_tension_forces_forces 

    subroutine update_surface_tension_stresses_ellipsoid(this)
        class(amrist), intent(inout) :: this
        print *, "This is the update_surface_tension_stresses_ellipsoid subroutine in the amrist module."
    end subroutine update_surface_tension_stresses_ellipsoid 

    subroutine update_tension_forces_forces_ellipsoid(this)
        class(amrist), intent(inout) :: this
        print *, "This is the update_tension_forces_forces_ellipsoid subroutine in the amrist module."
    end subroutine update_tension_forces_forces_ellipsoid 

    subroutine add_interface_to_PU_neighborhood(this)
        class(amrist), intent(inout) :: this
        print *, "This is the add_interface_to_PU_neighborhood subroutine in the amrist module."
    end subroutine add_interface_to_PU_neighborhood 

end module amrist_class