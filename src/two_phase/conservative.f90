module conservative_st
    use precision,         only: WP
    use mathtools,         only: Pi
    use geometry,          only: cfg,cfg_x,cfg_y,cfg_z
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
    use sgrid_class,       only: sgrid
    use config_class,      only: config
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan,ieee_is_finite
    use irl_fortran_interface
    
    implicit none
    private

    public :: temp,conservative_st_type 

    ! Parameters
    real(WP), parameter :: VFlo=1.0e-10_WP                         !< Minimum VF value considered
    real(WP), parameter :: VFhi=1.0_WP-VFlo                        !< Maximum VF value considered
    integer, parameter, dimension(4,2) :: neighbors= reshape( [-1, 1, 0, 0, 0, 0,-1, 1], [4,2] )
    !> Conservative Surface Tension Type
    type :: conservative_st_type
        ! Things we store here
        ! Options
        real(WP), dimension(3) :: MarangoniOption
        integer :: SurfaceTensionOption,CurvatureOption,SmoothingOption
        integer :: LevelsetOption,MatchingOption,ParaboloidOption,NormalOption
        logical :: TwoD
        real(WP) :: PressureOption
        ! Solvers
        type(tpns),pointer,public :: fs => null()
        type(vfs),pointer,public :: vf => null()
        type(vfs),public :: vf_x,vf_y,vf_z
        type(timetracker),pointer,public :: time => null()
        ! Working and Display Arrays
        real(WP) :: PU_spread ! How far the PU spreads
        real(WP), dimension(:,:,:), allocatable :: sigma_xx,sigma_yy,sigma_xy,sigma_yx ! Used Stress Tensor components
        real(WP), dimension(:,:,:), allocatable :: sigma_xx_P,sigma_yy_P,sigma_xy_P,sigma_yx_P ! Pressure Component
        real(WP), dimension(:,:,:), allocatable :: sigma_xx_NoP,sigma_yy_NoP,sigma_xy_NoP,sigma_yx_NoP ! Tangent Component
        real(WP), dimension(:,:,:,:,:), allocatable :: sigma_3D,sigma_3D_Exact,sigma_3D_OLD,dSigma,sigma_3D_P,sigma_3D_NoP, sigma_3D_DP
        real(WP), dimension(:,:,:), allocatable :: force_potential_field,poisson_source ! Potential Field Tools
        real(WP), dimension(:,:,:), allocatable :: Fst_x,Fst_y,Fst_z,Fst_mag ! Forces sent to momentum equation
        real(WP), dimension(:,:,:), allocatable :: Fst_x_3D,Fst_y_3D,Fst_z_3D
        real(WP), dimension(:,:,:), allocatable :: Fst_x_3D_OLD,Fst_y_3D_OLD,Fst_z_3D_OLD
        real(WP), dimension(:,:,:), allocatable :: Fst_x_3D_diff,Fst_y_3D_diff,Fst_z_3D_diff
        real(WP), dimension(:,:,:), allocatable :: Fst_x_3D_Exact,Fst_y_3D_Exact,Fst_z_3D_Exact
        real(WP), dimension(:,:,:), allocatable :: PjxD,PjyD,PjzD,CSF_Mag ! CSF Forces
        real(WP), dimension(:,:,:), allocatable :: Pjx_ST,Pjy_ST,Pjz_ST ! Forces from UpdateSTForces
        real(WP), dimension(:,:,:), allocatable :: Pjx_Marangoni,Pjy_Marangoni,Pjz_Marangoni ! Storage
        real(WP), dimension(:,:,:), allocatable :: Pjx_NoMarangoni,Pjy_NoMarangoni,Pjz_NoMarangoni ! Storage
        real(WP), dimension(:,:,:), allocatable :: DPjx_Marangoni,DPjy_Marangoni,DPjz_Marangoni ! Storage
        real(WP), dimension(:,:,:), allocatable :: grad_vf_x, grad_vf_y,grad_vf_z ! Gradient of volume fraction, but absolute value
        integer, dimension(:,:,:), allocatable :: x_smoothing_tracker,y_smoothing_tracker,height_function_type,HF_dir ! Track where smoothing is applied
        real(WP), dimension(:,:,:), allocatable :: SurfaceTensionDiv ! Divergence of surface tension force
        real(WP), dimension(:,:,:), allocatable :: ST_Coeff,Hx,Hy,STx,STy,dStds
        real(WP), dimension(:,:,:), allocatable :: Gx,Gy,Gz,Pjx_Seric,Pjy_Seric,Pjz_Seric,Stp1,St0,Stm1,Hp1,H0,Hm1
        real(WP), dimension(:,:,:), allocatable :: mag_tx,mag_ty
        real(WP), dimension(:,:,:), allocatable :: Seric_Populated,Seric_Pop_Temp
        real(WP), dimension(:,:,:,:), allocatable :: Seric_Pop_Storage,Gx_Storage,Gy_Storage,curvature_storage
        ! Abu-Al-Saud Methods
        real(WP), dimension(:,:,:), allocatable :: phi,phi_x,phi_y,phi_xx,phi_xy,phi_yy
        real(WP), dimension(:,:,:), allocatable :: curvature_abu
        ! Initailize Boolean
        logical :: initialized = .false.
        ! Visualization
        real(WP), dimension(:,:,:), allocatable :: PartitionOfUnityValue,PUTangent_x,PUTangent_y,PUTangent_z,PUTangent_mag,PartitionOfUnityWeight,DistanceField
        real(WP), dimension(:,:,:), allocatable :: trueDistanceField 
    contains
        procedure :: temp 
        procedure :: init
        procedure :: get_dmomdt
        procedure :: add_surface_tension_jump
        ! procedure :: end ! this is to deallocate variables, clear pointers, etc.
        ! Visualization
        procedure, public :: updatePartitionOfUnity
        procedure, public :: getProjectedGeometry
        procedure, public :: getProjectedGeometryEllipsoid
        procedure, public :: updateStresses
        procedure, public :: updateForces
        procedure, public :: getValuePU
        !Private Functions
        procedure, private :: get_dmomdt_full_ST
        procedure, private :: get_dmomdt_no_ST

        procedure, private :: add_surface_tension_jump_full_ST
        procedure, private :: add_surface_tension_jump_no_ST
        procedure, private :: add_conservative_surface_tension_jump
        procedure, private :: add_CSF_Shift_surface_tension_jump
        procedure, private :: add_Seric_surface_tension_jump
        procedure, private :: add_3D_surface_tension_jump
        procedure, private :: add_3D_surface_tension_jump_face_paraboloids
        procedure, private :: add_3D_CSF_Shift_surface_tension_jump

        procedure, private :: updateSurfaceTensionStresses
        procedure, private :: updateSurfaceTensionStresses3D
        procedure, private :: updateSurfaceTensionStresses3DFaceParaboloids
        procedure, private :: updateSurfaceTensionStresses3DEllipsoid
        procedure, private :: updateSurfaceTensionForces
        procedure, private :: updateSurfaceTensionForces3D
        procedure, private :: updateSurfaceTensionForces3DEllipsoid


        procedure, private :: applyLaplacianSmoothing
        procedure, private :: applyGradientSmoothing
        procedure, private :: applyPoissonSmoothing
        procedure, private :: applyPeskinSmoothing

        procedure, private :: addCellToNeighborhood

        procedure, private :: get_surface_tension_coefficient
        procedure, private :: update_surface_tension_coefficients
        procedure, private :: update_columns
        procedure, private :: get_column_surface_tension_and_height
        procedure, private :: height_function_curvature
        procedure, private :: mixed_young_center_normal
        procedure, private :: young_normal_method
        procedure, private :: centered_columns_normal_method
        procedure, private :: PeskinDeltaBase
        procedure, private :: PeskinDelta
        
        ! Abu-Al-Saud Methods
        procedure, private :: populate_levelset_with_pu
        procedure, private :: update_levelset_derivatives
        procedure, private :: populate_diagonal
        procedure, private :: populate_off_diagonal
        procedure, private :: add_Saud_surface_tension_jump
        ! procedure, private :: updateSaudForces
        ! procedure, private :: height_function_surface_tension

        ! Momentum VF Methods
        procedure, private :: fill_momentum_vfs
        procedure, private :: compute_paraboloid_force

        
    end type conservative_st_type
contains 

subroutine temp(this)
    class(conservative_st_type), intent(inout) :: this
    print *, "This is a temporary subroutine in the conservative_st module."
end subroutine temp

subroutine init(this,fs_in,vf_in,time_in)
    use vfs_class,only: r2p,lvira,VFhi,VFlo,jibben,plicnet,jibbenPU,jibbenSq,PUppic,jibbenPLIC
    class(conservative_st_type), intent(inout) :: this
    class(tpns),target, intent(in) :: fs_in
    class(vfs),target, intent(in) :: vf_in
    class(timetracker),target, intent(in) :: time_in
    ! Assign fs_in to fs
    this%fs => fs_in
    this%vf => vf_in
    this%time => time_in
    ! Create Momentum VF options
    call this%vf_x%initialize(cfg=cfg_x,reconstruction_method = PUppic,name='X-VOF')
    call this%vf_y%initialize(cfg=cfg_y,reconstruction_method = PUppic,name='Y-VOF')
    call this%vf_z%initialize(cfg=cfg_z,reconstruction_method = PUppic,name='Z-VOF')
    call fill_momentum_vfs(this)
    ! Allocate arrays
    allocate(this%sigma_xx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_xy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%sigma_xx_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yy_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_xy_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yx_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%sigma_xx_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yy_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_xy_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%sigma_yx_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Fst_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_mag(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Fst_x_3D(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_y_3D(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_z_3D(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Fst_x_3D_OLD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_y_3D_OLD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_z_3D_OLD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Fst_x_3D_diff(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_y_3D_diff(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_z_3D_diff(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Fst_x_3D_Exact(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_y_3D_Exact(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Fst_z_3D_Exact(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%PjxD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PjyD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PjzD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%CSF_Mag(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Pjx_ST(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjy_ST(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjz_ST(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Pjx_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjy_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjz_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Pjx_NoMarangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjy_NoMarangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjz_NoMarangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%DPjx_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%DPjy_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%DPjz_Marangoni(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%SurfaceTensionDiv(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%grad_vf_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%grad_vf_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%grad_vf_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%x_smoothing_tracker(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%y_smoothing_tracker(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%force_potential_field(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%poisson_source(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%ST_Coeff(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Hx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Hy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%STx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%STy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%dStds(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%height_function_type(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%HF_dir(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Stp1(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%St0(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Stm1(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Hp1(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%H0(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Hm1(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Gx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Gy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Gz(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%mag_tx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%mag_ty(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Pjx_Seric(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjy_Seric(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Pjz_Seric(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%Seric_Populated(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Seric_Pop_Temp(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Seric_Pop_Storage(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,0:2))
    allocate(this%Gx_Storage(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,0:2))
    allocate(this%Gy_Storage(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,0:2))
    allocate(this%curvature_storage(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,0:2))
    allocate(this%sigma_3D(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,1:3,1:3))
    allocate(this%dSigma(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,1:3,1:3))
    allocate(this%sigma_3D_OLD(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,1:3,1:3))
    allocate(this%sigma_3D_Exact(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,1:3,1:3))
    allocate(this%sigma_3D_P(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,1:3,1:3))
    allocate(this%sigma_3D_NoP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,1:3,1:3))
    allocate(this%sigma_3D_DP(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_,1:3,1:3))
    
    allocate(this%phi(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%phi_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%phi_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%phi_xx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%phi_xy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%phi_yy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%curvature_abu(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%PartitionOfUnityValue(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PUTangent_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PUTangent_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PUTangent_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PUTangent_mag(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%PartitionOfUnityWeight(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%DistanceField(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%trueDistanceField(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
end subroutine init

subroutine get_dmomdt(this,resU,resV,resW)
    class(conservative_st_type), intent(inout) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)

    SELECT CASE (this%SurfaceTensionOption)
        CASE (2)
            ! write(*,'(A)') 'CASE 2 ================================================'
            call this%add_surface_tension_jump_no_ST(dt=this%time%dt,div=this%fs%div)
            call this%get_dmomdt_full_ST(resU,resV,resW)
        CASE (3)
            ! write(*,'(A)') 'CASE 3 ================================================'
            call this%get_dmomdt_no_ST(resU,resV,resW)
        CASE (4)
            ! Conservative Surface Tension Approach for get_dmomdt
            call this%fs%get_dmomdt(resU,resV,resW)
        CASE (5)
            ! CSF Correction Approach
            call this%fs%get_dmomdt(resU,resV,resW)
        CASE DEFAULT
            ! write(*,'(A)') 'CASE DEFAULT ================================================'
            call this%fs%get_dmomdt(resU,resV,resW)
    END SELECT
end subroutine get_dmomdt 


subroutine add_surface_tension_jump(this)
    class(conservative_st_type), intent(inout) :: this

    ! call this%update_surface_tension_coefficients() 
    ! call this%update_columns()
    
    SELECT CASE (this%SurfaceTensionOption)
        CASE (2)
            call this%add_surface_tension_jump_no_ST(dt=this%time%dt,div=this%fs%div)
        CASE (3)
            call this%add_surface_tension_jump_full_ST(dt=this%time%dt,div=this%fs%div)
        CASE (4) 
            call this%add_conservative_surface_tension_jump(dt=this%time%dt,div=this%fs%div)
        CASE (5)
            call this%add_CSF_Shift_surface_tension_jump(dt=this%time%dt,div=this%fs%div)
        case (6)
            call this%add_Seric_surface_tension_jump(dt=this%time%dt,div=this%fs%div)
        case (7)
            call this%add_Saud_surface_tension_jump(dt=this%time%dt,div=this%fs%div)
        case (8)
            call this%add_3D_surface_tension_jump(dt=this%time%dt,div=this%fs%div)
        case(9)
            call this%add_3D_surface_tension_jump_face_paraboloids(dt=this%time%dt,div=this%fs%div)
        case(10)
            call this%add_3D_CSF_Shift_surface_tension_jump(dt=this%time%dt,div=this%fs%div)
        CASE DEFAULT
            call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
    END SELECT
end subroutine add_surface_tension_jump


! PRIVATE FUNCTIONS
! ########################################################### MOMENTUM FUNCTIONS #############################################
subroutine get_dmomdt_full_ST(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(conservative_st_type), intent(inout) :: this
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,jj,kk

      ! Zero out drhoUVW/dt arrays
      drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
      ! Flux of rhoU
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FUX(i,j,k)=this%fs%FUX(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybu_x(:,i,j,k)*this%fs%rho_Uold(i:i+1,j,k))*sum(this%fs%hybu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))*sum(this%fs%itpu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))+sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FUY(i,j,k)=this%fs%FUY(i,j,k) &
               &              -sum(this%fs%hybu_y(:,i,j,k)*this%fs%rho_Uold(i,j-1:j,k))*sum(this%fs%hybu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))*sum(this%fs%itpv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))+sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FUZ(i,j,k)=this%fs%FUZ(i,j,k) &
               &              -sum(this%fs%hybu_z(:,i,j,k)*this%fs%rho_Uold(i,j,k-1:k))*sum(this%fs%hybu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))*sum(this%fs%itpw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))+sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoU
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoUdt(i,j,k)=sum(this%fs%divu_x(:,i,j,k)*this%fs%FUX(i-1:i,j,k))+&
               &              sum(this%fs%divu_y(:,i,j,k)*this%fs%FUY(i,j:j+1,k))+&
               &              sum(this%fs%divu_z(:,i,j,k)*this%fs%FUZ(i,j,k:k+1))+this%fs%Pjx(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoUdt)
      
      ! Flux of rhoV
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FVX(i,j,k)=this%fs%FVX(i,j,k) &
               &              -sum(this%fs%hybv_x(:,i,j,k)*this%fs%rho_Vold(i-1:i,j,k))*sum(this%fs%hybv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))*sum(this%fs%itpu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))+sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)))
               ! Fluxes on y-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FVY(i,j,k)=this%fs%FVY(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybv_y(:,i,j,k)*this%fs%rho_Vold(i,j:j+1,k))*sum(this%fs%hybv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))*sum(this%fs%itpv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))+sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FVZ(i,j,k)=this%fs%FVZ(i,j,k) &
               &              -sum(this%fs%hybv_z(:,i,j,k)*this%fs%rho_Vold(i,j,k-1:k))*sum(this%fs%hybv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))*sum(this%fs%itpw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))+sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoV
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoVdt(i,j,k)=sum(this%fs%divv_x(:,i,j,k)*this%fs%FVX(i:i+1,j,k))+&
               &              sum(this%fs%divv_y(:,i,j,k)*this%fs%FVY(i,j-1:j,k))+&
               &              sum(this%fs%divv_z(:,i,j,k)*this%fs%FVZ(i,j,k:k+1))+this%fs%Pjy(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoVdt)
      
      ! Flux of rhoW
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FWX(i,j,k)=this%fs%FWX(i,j,k) &
               &              -sum(this%fs%hybw_x(:,i,j,k)*this%fs%rho_Wold(i-1:i,j,k))*sum(this%fs%hybw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))*sum(this%fs%itpu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))+sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FWY(i,j,k)=this%fs%FWY(i,j,k) &
               &              -sum(this%fs%hybw_y(:,i,j,k)*this%fs%rho_Wold(i,j-1:j,k))*sum(this%fs%hybw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))*sum(this%fs%itpv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))+sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)))
               ! Fluxes on z-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FWZ(i,j,k)=this%fs%FWZ(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybw_z(:,i,j,k)*this%fs%rho_Wold(i,j,k:k+1))*sum(this%fs%hybw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))*sum(this%fs%itpw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))+sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)))
            end do
         end do
      end do
      ! Time derivative of rhoW
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoWdt(i,j,k)=sum(this%fs%divw_x(:,i,j,k)*this%fs%FWX(i:i+1,j,k))+&
               &              sum(this%fs%divw_y(:,i,j,k)*this%fs%FWY(i,j:j+1,k))+&
               &              sum(this%fs%divw_z(:,i,j,k)*this%fs%FWZ(i,j,k-1:k))+this%fs%Pjz(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoWdt)
   end subroutine get_dmomdt_full_ST

   subroutine get_dmomdt_no_ST(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(conservative_st_type), intent(inout) :: this
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,jj,kk

      
      ! Zero out drhoUVW/dt arrays
      drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
      ! Flux of rhoU
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FUX(i,j,k)=this%fs%FUX(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybu_x(:,i,j,k)*this%fs%rho_Uold(i:i+1,j,k))*sum(this%fs%hybu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))*sum(this%fs%itpu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))+sum(this%fs%grdu_x(:,i,j,k)*this%fs%U(i:i+1,j,k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FUY(i,j,k)=this%fs%FUY(i,j,k) &
               &              -sum(this%fs%hybu_y(:,i,j,k)*this%fs%rho_Uold(i,j-1:j,k))*sum(this%fs%hybu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))*sum(this%fs%itpv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k))+sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FUZ(i,j,k)=this%fs%FUZ(i,j,k) &
               &              -sum(this%fs%hybu_z(:,i,j,k)*this%fs%rho_Uold(i,j,k-1:k))*sum(this%fs%hybu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))*sum(this%fs%itpw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k))+sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoU
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoUdt(i,j,k)=sum(this%fs%divu_x(:,i,j,k)*this%fs%FUX(i-1:i,j,k))+&
               &              sum(this%fs%divu_y(:,i,j,k)*this%fs%FUY(i,j:j+1,k))+&
               &              sum(this%fs%divu_z(:,i,j,k)*this%fs%FUZ(i,j,k:k+1))+this%fs%Pjx(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoUdt)
      
      ! Flux of rhoV
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FVX(i,j,k)=this%fs%FVX(i,j,k) &
               &              -sum(this%fs%hybv_x(:,i,j,k)*this%fs%rho_Vold(i-1:i,j,k))*sum(this%fs%hybv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))*sum(this%fs%itpu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)) &
               &              +this%fs%visc_xy(i,j,k)*(sum(this%fs%grdv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))+sum(this%fs%grdu_y(:,i,j,k)*this%fs%U(i,j-1:j,k)))
               ! Fluxes on y-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FVY(i,j,k)=this%fs%FVY(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybv_y(:,i,j,k)*this%fs%rho_Vold(i,j:j+1,k))*sum(this%fs%hybv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))*sum(this%fs%itpv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k))+sum(this%fs%grdv_y(:,i,j,k)*this%fs%V(i,j:j+1,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               this%fs%FVZ(i,j,k)=this%fs%FVZ(i,j,k) &
               &              -sum(this%fs%hybv_z(:,i,j,k)*this%fs%rho_Vold(i,j,k-1:k))*sum(this%fs%hybv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))*sum(this%fs%itpw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k))+sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoV
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoVdt(i,j,k)=sum(this%fs%divv_x(:,i,j,k)*this%fs%FVX(i:i+1,j,k))+&
               &              sum(this%fs%divv_y(:,i,j,k)*this%fs%FVY(i,j-1:j,k))+&
               &              sum(this%fs%divv_z(:,i,j,k)*this%fs%FVZ(i,j,k:k+1))+this%fs%Pjy(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoVdt)
      
      ! Flux of rhoW
      do kk=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
         do jj=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do ii=this%fs%cfg%imin_,this%fs%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               this%fs%FWX(i,j,k)=this%fs%FWX(i,j,k) &
               &              -sum(this%fs%hybw_x(:,i,j,k)*this%fs%rho_Wold(i-1:i,j,k))*sum(this%fs%hybw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))*sum(this%fs%itpu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)) &
               &              +this%fs%visc_zx(i,j,k)*(sum(this%fs%grdw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))+sum(this%fs%grdu_z(:,i,j,k)*this%fs%U(i,j,k-1:k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               this%fs%FWY(i,j,k)=this%fs%FWY(i,j,k) &
               &              -sum(this%fs%hybw_y(:,i,j,k)*this%fs%rho_Wold(i,j-1:j,k))*sum(this%fs%hybw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))*sum(this%fs%itpv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)) &
               &              +this%fs%visc_yz(i,j,k)*(sum(this%fs%grdw_y(:,i,j,k)*this%fs%W(i,j-1:j,k))+sum(this%fs%grdv_z(:,i,j,k)*this%fs%V(i,j,k-1:k)))
               ! Fluxes on z-face
               i=ii-1; j=jj-1; k=kk-1
               this%fs%FWZ(i,j,k)=this%fs%FWZ(i,j,k)-this%fs%P(i,j,k) &
               &              -sum(this%fs%hybw_z(:,i,j,k)*this%fs%rho_Wold(i,j,k:k+1))*sum(this%fs%hybw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))*sum(this%fs%itpw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)) &
               &              +this%fs%visc   (i,j,k)*(sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1))+sum(this%fs%grdw_z(:,i,j,k)*this%fs%W(i,j,k:k+1)))
            end do
         end do
      end do
      ! Time derivative of rhoW
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               drhoWdt(i,j,k)=sum(this%fs%divw_x(:,i,j,k)*this%fs%FWX(i:i+1,j,k))+&
               &              sum(this%fs%divw_y(:,i,j,k)*this%fs%FWY(i,j:j+1,k))+&
               &              sum(this%fs%divw_z(:,i,j,k)*this%fs%FWZ(i,j,k-1:k)) ! +this%fs%Pjz(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%fs%cfg%sync(drhoWdt)
      
   end subroutine get_dmomdt_no_ST

! ############################################################### SURFACE TENSION JUMP SUBROUTINES ##########################################################
subroutine add_surface_tension_jump_full_ST(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Set
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do
    
    ! Add wall contact force to pressure jump
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%Pjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                        +sum(this%fs%divp_y(:,i,j,k)*this%fs%Pjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k))&
            &                        +sum(this%fs%divp_z(:,i,j,k)*this%fs%Pjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
        end do
        end do
    end do
end subroutine add_surface_tension_jump_full_ST

subroutine add_surface_tension_jump_no_ST(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    implicit none
    
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_,this%fs%cfg%jmino_,this%fs%cfg%kmino_) :: addToDiv
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf, magx,magy,magz,magdiv
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do
    
    ! Add wall contact force to pressure jump
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Calculate Magnitude of DPs
    magx = NORM2(this%fs%DPjx)
    magy = NORM2(this%fs%DPjy)
    magz = NORM2(this%fs%DPjz)
    magdiv = NORM2(this%fs%div)
    addToDiv = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                        +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k))&
            &                        +sum(this%fs%divp_z(:,i,j,k)*this%fs%DPjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
    
    magz = NORM2(addToDiv)
    write(*,'(A,4F35.10)') '   magx,magy,magz,magdiv: ', magx,magy,magz,magdiv
    if(magz.gt.0) then
        write(*,'(A,F10.5)')  '====================================MagToDiv: ',magz
    endif
end subroutine add_surface_tension_jump_no_ST

subroutine add_conservative_surface_tension_jump(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    ! Update Values
    call this%updateSurfaceTensionStresses()
    call this%updateSurfaceTensionStresses3D()
    call this%updateSurfaceTensionForces()
    call this%updateSurfaceTensionForces3D()
    
    this%fs%Pjx = this%Pjx_ST
    this%fs%Pjy = this%Pjy_ST
    this%fs%Pjz = this%Pjz_ST
    
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
        call this%applyLaplacianSmoothing()
        CASE(3)
        call this%applyGradientSmoothing()
        CASE(4)
        call this%applyPoissonSmoothing()
        case(5) ! This is using the Peskin Delta Function
        call this%applyPeskinSmoothing()
        CASE DEFAULT

    END SELECT

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                                  +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k)))
            ! &                           +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_conservative_surface_tension_jump

subroutine add_3D_surface_tension_jump(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    ! Update Values
    call this%updateSurfaceTensionStresses()
    call this%updateSurfaceTensionForces()
    call this%updateSurfaceTensionStresses3D()
    call this%updateSurfaceTensionForces3D()
    
    this%fs%Pjx = this%Fst_x_3D 
    this%fs%Pjy = this%Fst_y_3D
    this%fs%Pjz = this%Fst_z_3D

    ! print *,sqrt(sum(this%Fst_z_3D**2))
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
        call this%applyLaplacianSmoothing()
        CASE(3)
        call this%applyGradientSmoothing()
        CASE(4)
        call this%applyPoissonSmoothing()
        case(5) ! This is using the Peskin Delta Function 
        call this%applyPeskinSmoothing()
        CASE DEFAULT

    END SELECT

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                                  +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k))&
            &                                  +sum(this%fs%divp_z(:,i,j,k)*this%fs%DPjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_3D_surface_tension_jump

subroutine add_3D_surface_tension_jump_face_paraboloids(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz
    this%CSF_Mag = sqrt(this%PjxD**2 + this%PjyD**2 + this%PjzD**2)
    ! Update Values
    call this%updateSurfaceTensionStresses3DFaceParaboloids()
    call this%updateSurfaceTensionForces3D()
    
    this%fs%Pjx = this%Fst_x_3D 
    this%fs%Pjy = this%Fst_y_3D
    this%fs%Pjz = this%Fst_z_3D

    ! print *,sqrt(sum(this%Fst_z_3D**2))
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
        call this%applyLaplacianSmoothing()
        CASE(3)
        call this%applyGradientSmoothing()
        CASE(4)
        call this%applyPoissonSmoothing()
        case(5) ! This is using the Peskin Delta Function 
        call this%applyPeskinSmoothing()
        CASE DEFAULT

    END SELECT

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                                  +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k))&
            &                                  +sum(this%fs%divp_z(:,i,j,k)*this%fs%DPjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
            ! endif
            div(i,j,k)=div(i,j,k)+this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_3D_surface_tension_jump_face_paraboloids

   
subroutine add_CSF_Shift_surface_tension_jump(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    real(WP), dimension(3) :: OldMarangoniOption
    !! THIS IS WHERE WE GET THE CSF VALUES !!
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    !! THIS IS WHERE WE GET THE PCST VALUES, WITH THE CURRENT MARANGONI OPTION
    ! Update Values
    call updateSurfaceTensionStresses(this)
    call updateSurfaceTensionForces(this)
    call this%updateSurfaceTensionForces3D()
    ! Store Force Values
    this%Pjx_Marangoni = this%Pjx_ST
    this%Pjy_Marangoni = this%Pjx_ST
    this%Pjz_Marangoni = this%Pjx_ST

    !! THIS IS WHERE WE GET THE PCST VALUES WITH CONSTANT SURFACE TENSION COEFFICIENT
    OldMarangoniOption = this%MarangoniOption
    this%MarangoniOption = (/0.0_WP,0.0_WP,0.0_WP/)
    call updateSurfaceTensionStresses(this)
    call updateSurfaceTensionForces(this)
    call this%updateSurfaceTensionForces3D()
    this%Pjx_NoMarangoni = this%Pjx_ST
    this%Pjy_NoMarangoni = this%Pjx_ST
    this%Pjz_NoMarangoni = this%Pjx_ST

    this%MarangoniOption = OldMarangoniOption   

    ! Find the Difference caused by Marangoni
    this%DPjx_Marangoni= this%Pjx_Marangoni- this%Pjx_NoMarangoni
    this%DPjy_Marangoni = this%Pjy_Marangoni- this%Pjy_NoMarangoni
    this%DPjz_Marangoni = this%Pjz_Marangoni- this%Pjz_NoMarangoni

    ! Assign to Pjx so that we can smooth it
    this%fs%Pjx = this%DPjx_Marangoni
    this%fs%Pjy = this%DPjy_Marangoni
    this%fs%Pjz = this%DPjz_Marangoni
    ! We might need to smooth this so that the shape of the force field is the same as the CSF model
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
            call applyLaplacianSmoothing(this)
        CASE(3)
            call applyGradientSmoothing(this)
        CASE(4)
            call applyPoissonSmoothing(this)
        case(5) ! This is using the Peskin Delta Function
            call this%applyPeskinSmoothing()
        CASE DEFAULT

    END SELECT
    ! Repalce for visualization
    this%DPjx_Marangoni = this%fs%Pjx
    this%DPjy_Marangoni = this%fs%Pjy
    this%DPjz_Marangoni = this%fs%Pjz
    ! Now we have (F_PCST - F_PCST_CONST_ST). We now add this to the CSF values and use that
    this%fs%Pjx = this%DPjx_Marangoni + this%PjxD
    this%fs%Pjy = this%DPjy_Marangoni + this%PjyD
    this%fs%Pjz = this%DPjz_Marangoni + this%PjzD

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                             +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k)))
            ! &                           +sum(this%fs%divp_z(:,i,j,k)*this%fs%DPjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%fs%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_CSF_Shift_surface_tension_jump

subroutine add_3D_CSF_Shift_surface_tension_jump(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    real(WP), dimension(3) :: OldMarangoniOption
    !! THIS IS WHERE WE GET THE CSF VALUES !!
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    !! THIS IS WHERE WE GET THE PCST VALUES, WITH THE CURRENT MARANGONI OPTION
    ! Update Values
    call updateSurfaceTensionStresses3D(this)
    call updateSurfaceTensionForces3D(this)
    ! Store Force Values
    this%Pjx_Marangoni = this%Fst_x_3D
    this%Pjy_Marangoni = this%Fst_y_3D
    this%Pjz_Marangoni = this%Fst_z_3D

    !! THIS IS WHERE WE GET THE PCST VALUES WITH CONSTANT SURFACE TENSION COEFFICIENT
    OldMarangoniOption = this%MarangoniOption
    this%MarangoniOption = (/0.0_WP,0.0_WP,0.0_WP/)
    call updateSurfaceTensionStresses3D(this)
    call updateSurfaceTensionForces3D(this)
    this%Pjx_NoMarangoni = this%Fst_x_3D
    this%Pjy_NoMarangoni = this%Fst_y_3D
    this%Pjz_NoMarangoni = this%Fst_z_3D
    
    this%MarangoniOption = OldMarangoniOption   

    ! Find the Difference caused by Marangoni
    this%DPjx_Marangoni= this%Pjx_Marangoni- this%Pjx_NoMarangoni
    this%DPjy_Marangoni = this%Pjy_Marangoni- this%Pjy_NoMarangoni
    this%DPjz_Marangoni = this%Pjz_Marangoni- this%Pjz_NoMarangoni

    ! Assign to Pjx so that we can smooth it
    this%fs%Pjx = this%DPjx_Marangoni
    this%fs%Pjy = this%DPjy_Marangoni
    this%fs%Pjz = this%DPjz_Marangoni
    ! We might need to smooth this so that the shape of the force field is the same as the CSF model
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
            call applyLaplacianSmoothing(this)
        CASE(3)
            call applyGradientSmoothing(this)
        CASE(4)
            call applyPoissonSmoothing(this)
        case(5) ! This is using the Peskin Delta Function
            call this%applyPeskinSmoothing()
        CASE DEFAULT

    END SELECT
    ! Repalce for visualization
    this%DPjx_Marangoni= this%fs%Pjx
    this%DPjy_Marangoni = this%fs%Pjy
    this%DPjz_Marangoni = this%fs%Pjz
    ! Now we have (F_PCST - F_PCST_CONST_ST). We now add this to the CSF values and use that
    this%fs%Pjx = this%DPjx_Marangoni + this%PjxD
    this%fs%Pjy = this%DPjy_Marangoni + this%PjyD
    this%fs%Pjz = this%DPjz_Marangoni + this%PjzD

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                             +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k))&
            &                                  +sum(this%fs%divp_z(:,i,j,k)*this%fs%DPjz(i,j,k:k+1)/this%fs%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%fs%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_3D_CSF_Shift_surface_tension_jump

subroutine add_Seric_surface_tension_jump(this,dt,div,contact_model)
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k,dir,n,m
    real(WP) :: mycurv,mysurf, count,delta,myGradient
    real(WP), dimension(3) :: planeNormal, tangent,planeNormal1,planeNormal2

    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    ! Update Update Curvature and Surface Gradients
    do k=this%fs%cfg%kmin_-1,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_-1,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_-1,this%fs%cfg%imax_+1
            call this%mixed_young_center_normal(i,j,k,planeNormal)
            ! print *, "Plane Normal = ", planeNormal
            if(abs(planeNormal(1)) .gt. abs(planeNormal(2))) then ! Plane more vertical, X orient
                dir = int(sign(1.0_WP,planeNormal(1)))
            else
                dir = int(sign(2.0_WP,planeNormal(2)))
            endif
            if(this%vf%VF(i,j,k) .gt. VFlo .and. this%vf%VF(i,j,k) .lt. VFhi) then 
                this%HF_dir(i,j,k) = dir
                this%Seric_Populated(i,j,k) = 1.0_WP
            else
                this%HF_dir(i,j,k) = 0
                this%Seric_Populated(i,j,k) = -1.0_WP
            endif
            call this%height_function_curvature(i,j,k,dir,this%vf%curv(i,j,k),this%dStds(i,j,k),this%height_function_type(i,j,k))
            ! Calculate Tangent
            ! Above we have the plane normal and the direction. 
            ! If we have a positive y direction (2) then we want our tangent to have a positive x component (0,1) -> (1,0)
            ! If we ahve a negative y direction (-2) then we want our tangent to have a negative x component (0,-1) -> (-1,0)
            ! If we have a positive x direction (1) then we want our tangent to have a negative y component (1,0) -> (0,-1)
            ! If we have a negative x direction (-1) then we want our tangent to have a positive y component. (-1,0) -> (0,1)
            tangent = (/-planeNormal(2),planeNormal(1),planeNormal(3)/) ! 2D tangent construction
            ! Calculate G
            this%Gx(i,j,k) = this%dStds(i,j,k) * sign(1.0_WP,tangent(1))
            this%Gy(i,j,k) = this%dStds(i,j,k) * sign(1.0_WP,tangent(2))
            this%Gz(i,j,k) = 0.0_WP
        enddo
        enddo
    enddo
    call this%vf%cfg%sync(this%vf%curv)
    call this%vf%cfg%sync(this%dStds)
    call this%vf%cfg%sync(this%Gx)
    call this%vf%cfg%sync(this%Gy)
    call this%vf%cfg%sync(this%Gz)
    
    ! Apply Smoothing to G
    this%Seric_Pop_Storage(:,:,:,0) = this%Seric_Populated
    this%Gx_Storage(:,:,:,0) = this%Gx
    this%Gy_Storage(:,:,:,0) = this%Gy
    this%curvature_storage(:,:,:,0) = this%vf%curv

    
    ! Now that we have Gx and Gy, we need to calculat the forces
    this%Fst_x = 0.0_WP
    this%Fst_y = 0.0_WP
    this%Fst_z = 0.0_WP
    ! Calculate forces 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! These are the deltas
            this%grad_vf_x(i,j,k) = sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            this%grad_vf_y(i,j,k) = sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            delta = sqrt(this%grad_vf_y(i,j,k)*this%grad_vf_y(i,j,k) + this%grad_vf_x(i,j,k)*this%grad_vf_x(i,j,k))
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k)) ! Surface Area
            ! Next we want to smooth the surface tension 
            if (mysurf.gt.0.0_WP) then
                ! Curvature interp
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
                call this%mixed_young_center_normal(i-1,j,k,planeNormal1)
                call this%mixed_young_center_normal(i,j,k,planeNormal2)
                ! Tangent Interp
                planeNormal(1) = -(this%vf%SD(i-1,j,k)*planeNormal1(2)*this%fs%cfg%vol(i-1,j,k) + this%vf%SD(i,j,k)*planeNormal2(2)*this%fs%cfg%vol(i,j,k) )/mysurf
                ! Gradient Interp
                myGradient=sum(this%vf%SD(i-1:i,j,k)*this%dStds(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
                planeNormal(1) = 0.0_WP
                myGradient = 0.0_WP 
            end if
            this%Fst_x(i,j,k) = myGradient * planeNormal(1) * delta ! Tangent Component
            this%fs%Pjx(i,j,k) = this%ST_Coeff(i,j,k) * mycurv * this%grad_vf_x(i,j,k)! Normal Component 
            ! if(this%vf%VF(i,j,k) .gt. VFlo .and. this%vf%VF(i,j,k) .lt. VFhi) then 
            !     print *,"X Face",myGradient,this%dStds(i-1:i,j,k),this%vf%SD(i-1:i,j,k)
            ! endif

            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
                call this%mixed_young_center_normal(i,j-1,k,planeNormal1)
                call this%mixed_young_center_normal(i,j,k,planeNormal2)
                ! Tangent Interp
                planeNormal(2) = (this%vf%SD(i,j-1,k)*planeNormal1(1)*this%fs%cfg%vol(i,j-1,k) + this%vf%SD(i,j,k)*planeNormal2(1)*this%fs%cfg%vol(i,j,k) )/mysurf
                ! Gradient Interp
                myGradient=sum(this%vf%SD(i,j-1:j,k)*this%dStds(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
                planeNormal(2) = 0.0_WP
                myGradient = 0.0_WP 
            end if
            this%Fst_y(i,j,k) = myGradient * planeNormal(2) * delta ! Tangent Component
            this%fs%Pjy(i,j,k) = this%ST_Coeff(i,j,k) * mycurv * this%grad_vf_y(i,j,k)! Normal Component 
            ! Z face
            this%fs%Pjz(i,j,k)=0.0_WP
        end do
        end do
    end do
    
    this%fs%Pjx = this%Fst_x + this%fs%Pjx
    this%fs%Pjy = this%Fst_y + this%fs%Pjy
    this%fs%Pjz = this%Fst_z

    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
        call this%applyLaplacianSmoothing()
        CASE(3)
        call this%applyGradientSmoothing()
        CASE(4)
        call this%applyPoissonSmoothing()
        CASE DEFAULT

    END SELECT

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                                  +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k)))
            ! &                           +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_Seric_surface_tension_jump

! ############################################# PartitionOfUnity Functions
subroutine updateSurfaceTensionStresses(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver
    
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,startPoint,endPoint,Gcen,Lcen,force,pos
    real(WP), dimension(1:3) :: P0,P1,P2,P3
    integer, dimension(1:3) :: ind,indguess
    real(WP) :: dx,dy,C,xEval,yEval,zEval

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    C = 4.0_WP
    ! Create Neighborhood and solver
    call new(neighborhood) 
    call new(solver)
    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then

                ! Empty Neighborhood
                call emptyNeighborhood(neighborhood) 
                ! Add 7x7 (3 on each side) stencil, in plane
                ! write(*,'(A)') '=============== New: '
                do j_in = -3,3
                    do i_in = -3,3
                    if(this%vf%VF(i+i_in,j+j_in,k) .gt. this%vf%VFmin .and. this%vf%VF(i+i_in,j+j_in,k) .lt. this%vf%VFmax) then
                        ! For now, compute centroid as cell center
                        ! Gcen = vf%Gbary(:,i+i_in,j+j_in,k)
                        ! Lcen = vf%Lbary(:,i+i_in,j+j_in,k)
                        ! ! cen = vf%VF(i+i_in,j+j_in,k) * Lcen + (1-vf%VF(i+i_in,j+j_in,k)) * Gcen
                        ! ! cen = (Lcen+Gcen)/2
                        ! cen = (/vf%cfg%xm(i+i_in),vf%cfg%ym(j+j_in),vf%cfg%zm(k)/)
                        cen = calculateCentroid(this%vf%interface_polygon(1,i+i_in,j+j_in,k))
                        ! Get Plane
                        plane = this%vf%liquid_gas_interface(i+i_in,j+j_in,k)
                        ! if((vf%VF(i+i_in,j+j_in,k) .gt. VFlo .and. vf%VF(i+i_in,j+j_in,k) .lt. VFhi) .and. (i .eq. 16) .and. (j .eq. 24)) then
                        !    write(*,'(A,3I10.5)') 'Mid Cell: ', i,j,k
                        !    write(*,'(A,3I10.5)') 'Curr Cell: ', i+i_in,j+j_in,k 
                        !    write(*,'(A,3F10.5)') 'Cen: ', cen 
                        !    call printToScreen(plane)
                        call addMember(neighborhood,cen,1.0_WP,plane,0.0_WP)
                    endif
                    end do
                end do
                ! Set Neighborhood in solver
                call setNeighborhood(solver,neighborhood)
                ! ====== Get Stresses

                ! sigma_xx, we go on right of u cell
                ! The right of the u cell goes y(j) to y(j+1) at xm(i) 
                startPoint = (/this%fs%cfg%xm(i),this%fs%cfg%y(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%xm(i),this%fs%cfg%y(j+1),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the x component of this force, so take it and store it as sigma_xx 
                this%sigma_xx(i,j,k) = force(1)
                ! ====== For Debugging
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_xx_NoP(i,j,k) = force(1)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_xx_P(i,j,k) = force(1)
                this%sigma_xx_P(i,j,k) = this%sigma_xx_P(i,j,k) - this%sigma_xx_NoP(i,j,k)


                ! sigma_xy, we go on bottom of u cell
                ! The bottom of the u cell goes xm(i-1) to xm(i) at y(j)
                startPoint = (/this%fs%cfg%xm(i-1),this%fs%cfg%y(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%xm(i),this%fs%cfg%y(j),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the x component of this force, so take it and store it as sigma_xy 
                this%sigma_xy(i,j,k) = force(1)

                ! ====== For Debugging
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_xy_NoP(i,j,k) = force(1)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_xy_P(i,j,k) = force(1)
                this%sigma_xy_P(i,j,k) = this%sigma_xy_P(i,j,k) - this%sigma_xy_NoP(i,j,k)

                ! sigma_yx, we go on left of v cell
                ! The bottom of the u cell goes ym(j-1) to ym(j) at x(i)
                startPoint = (/this%fs%cfg%x(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%x(i),this%fs%cfg%ym(j-1),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the y component of this force, so take it and store it as sigma_yx 
                this%sigma_yx(i,j,k) = force(2)
                
                ! ====== For Debugging
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_yx_NoP(i,j,k) = force(2)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_yx_P(i,j,k) = force(2)
                this%sigma_yx_P(i,j,k) = this%sigma_yx_P(i,j,k) - this%sigma_yx_NoP(i,j,k)

                ! sigma_yy, we go on top of v cell
                ! The bottom of the u cell goes x(i) to x(i+1) at ym(j) 
                startPoint = (/this%fs%cfg%x(i+1),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                endPoint = (/this%fs%cfg%x(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                ! call solveEdge(solver,fs%sigma,startPoint,endPoint,radius,center,5.0_WP,force)
                ! We only want the y component of this force, so take it and store it as sigma_yy
                this%sigma_yy(i,j,k) = force(2)

                ! ====== For Debugging  
                ! No Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,0.0_WP,this%MarangoniOption,force)
                this%sigma_yy_NoP(i,j,k) = force(2)
                ! With Pressure
                call solveEdge(solver,this%fs%sigma,startPoint,endPoint,this%PU_spread*dx,1.0_WP,this%MarangoniOption,force)
                this%sigma_yy_P(i,j,k) = force(2)
                this%sigma_yy_P(i,j,k) = this%sigma_yy_P(i,j,k) - this%sigma_yy_NoP(i,j,k)
            ! endif
        end do
        end do
    end do
    ! Boundary Conditions
    call this%vf%cfg%sync(this%sigma_xx)
    call this%vf%cfg%sync(this%sigma_xy)
    call this%vf%cfg%sync(this%sigma_yx)
    call this%vf%cfg%sync(this%sigma_yy)

    call this%vf%cfg%sync(this%sigma_xx_NoP)
    call this%vf%cfg%sync(this%sigma_xy_NoP)
    call this%vf%cfg%sync(this%sigma_yx_NoP)
    call this%vf%cfg%sync(this%sigma_yy_NoP)

    call this%vf%cfg%sync(this%sigma_xx_P)
    call this%vf%cfg%sync(this%sigma_xy_P)
    call this%vf%cfg%sync(this%sigma_yx_P)
    call this%vf%cfg%sync(this%sigma_yy_P)
end subroutine updateSurfaceTensionStresses

subroutine updateSurfaceTensionStresses3D(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class, only: solveFace,solveFaceParaboloid

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in,k_in,shift_first_index,shift_second_index ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood,neighborhood_trimmed
    type(PUST_RectCub_type) :: solver
    integer, dimension(1:3) :: shift_index_normal,shift_index_momentum,index
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    type(SeparatorVariant_type) :: interface_component
    real(WP), dimension(1:3) :: cen,corner,dvec,shift,force,face_normal
    real(WP), dimension(1:3) :: force1,force2,force3,force4
    real(WP), dimension(1:3) :: P0,P1,P2,P3,pressure_cell_center,velocity_cell_center,face_center
    real(WP), dimension(1:3) :: R0,R1,R2,R3
    real(WP) :: dx,dy,dz
    real(WP) :: x_i,y_j,z_k

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    dz = this%vf%cfg%dz(1)
    dvec = (/dx,dy,dz/)
    shift = dvec
    ! Zero out values
    this%sigma_3D_OLD = 0.0_WP 
    this%sigma_3D = 0.0_WP
    ! Create Neighborhood and solver  
    call new(neighborhood) 
    call new(neighborhood_trimmed) 
    call new(solver)
    if(this%LevelsetOption .eq. 3) then 
         call this%fill_momentum_vfs()
    endif
    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then   

                ! Empty Neighborhood
                call emptyNeighborhood(neighborhood) 
                ! Add 7x7x7 (3 on each side) stencil, in plane  
                ! write(*,'(A)') '=============== New: ' 
                SELECT CASE(this%LevelsetOption)
                    CASE(3)
                        do j_in = -3,3
                            do i_in = -3,3
                                do k_in = -3,3 
                                    call this%addCellToNeighborhood(neighborhood,i+i_in,j+j_in,k+k_in)
                                enddo
                            end do
                        end do
                        ! Set Neighborhood in solver
                        call setNeighborhood(solver,neighborhood)
                        call setKernelSize(solver,this%PU_spread*this%vf%cfg%dx(1))

                        ! Empty Neighborhood
                        call emptyNeighborhood(neighborhood_trimmed) 
                        ! Add 7x7x7 (3 on each side) stencil, in plane  
                        ! write(*,'(A)') '=============== New: ' 
                        do j_in = -3,3
                            do i_in = -3,3
                                do k_in = -3,3 
                                    call this%addCellToNeighborhood(neighborhood_trimmed,i+i_in,j+j_in,k+k_in,solver)
                                enddo
                            end do
                        end do
                    CASE(2)
                        do j_in = -3,3
                            do i_in = -3,3
                                do k_in = -3,3 
                                    call this%addCellToNeighborhood(neighborhood,i+i_in,j+j_in,k+k_in)
                                enddo
                            end do
                        end do
                        ! Set Neighborhood in solver
                        call setNeighborhood(solver,neighborhood)
                        call setKernelSize(solver,this%PU_spread*this%vf%cfg%dx(1))

                        ! Empty Neighborhood
                        call emptyNeighborhood(neighborhood_trimmed) 
                        ! Add 7x7x7 (3 on each side) stencil, in plane  
                        ! write(*,'(A)') '=============== New: ' 
                        do j_in = -3,3
                            do i_in = -3,3
                                do k_in = -3,3 
                                    call this%addCellToNeighborhood(neighborhood_trimmed,i+i_in,j+j_in,k+k_in,solver)
                                enddo
                            end do
                        end do
                    CASE DEFAULT
                        do j_in = -3,3
                            do i_in = -3,3
                                do k_in = -3,3 
                                    call this%addCellToNeighborhood(neighborhood,i+i_in,j+j_in,k+k_in)
                                enddo
                            end do
                        end do
                        ! Set Neighborhood in solver
                        call setNeighborhood(solver,neighborhood)
                        call setKernelSize(solver,this%PU_spread*this%vf%cfg%dx(1))

                        ! Empty Neighborhood
                        call emptyNeighborhood(neighborhood_trimmed) 
                        ! Add 7x7x7 (3 on each side) stencil, in plane  
                        ! write(*,'(A)') '=============== New: ' 
                        do j_in = -3,3
                            do i_in = -3,3
                                do k_in = -3,3 
                                    call this%addCellToNeighborhood(neighborhood_trimmed,i+i_in,j+j_in,k+k_in,solver)
                                enddo
                            end do
                        end do
                END SELECT
                ! Set Neighborhood in solver
                call setNeighborhood(solver,neighborhood_trimmed)
                call setKernelSize(solver,this%PU_spread*this%vf%cfg%dx(1))
                ! ====== Get Stresses
                ! call printSolver(solver)
                ! We are going to loop over the stress components and calculate them
                x_i = this%fs%cfg%x(i)
                y_j = this%fs%cfg%y(j)
                z_k = this%fs%cfg%z(k)
                corner = (/x_i,y_j,z_k/)
                ! First, get the cell center of the cell we are in
                pressure_cell_center = (/this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                do i_in = 1,3 ! Cell Loop, This tells us the force direction
                    ! Set up shift to velocity cell
                    shift = (/0.0_WP,0.0_WP,0.0_WP/)  
                    shift(i_in) = dvec(i_in)/2
                    shift_index_momentum = (/0,0,0/)
                    shift_index_momentum(i_in) = 1
                    ! Move center to velocity cell we are looking at
                    velocity_cell_center = pressure_cell_center - shift
                    do j_in = 1,3 ! Face Loop, This tells us the normal direction
                        force = (/0.0_WP,0.0_WP,0.0_WP/)
                        ! Shift from velocity cell center ot the face center we are looking at:
                        shift = (/0.0_WP,0.0_WP,0.0_WP/)  
                        shift_index_normal = (/0,0,0/)
                        shift_index_normal(j_in) = 1
                        face_normal = (/0.0_WP,0.0_WP,0.0_WP/)
                        face_normal(j_in) = 1.0_WP
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
                        SELECT CASE (this%LevelsetOption)
                            CASE(3)
                                force = (/0.0_WP,0.0_WP,0.0_WP/)

                                ! Depending on the momentum cell, use the appropriate VF and then 
                                if(i_in .eq. 1) then 
                                    call compute_paraboloid_force(this,this%vf_x,i,j,k,j_in,P0,P1,P2,P3,force)
                                elseif(i_in .eq. 2) then 
                                    call compute_paraboloid_force(this,this%vf_y,i,j,k,j_in,P0,P1,P2,P3,force)
                                else
                                    call compute_paraboloid_force(this,this%vf_z,i,j,k,j_in,P0,P1,P2,P3,force)
                                endif
                            CASE(2)
                                force = (/0.0_WP,0.0_WP,0.0_WP/)
                                force1 = (/0.0_WP,0.0_WP,0.0_WP/)
                                force2 = (/0.0_WP,0.0_WP,0.0_WP/)
                                force3 = (/0.0_WP,0.0_WP,0.0_WP/)
                                force4 = (/0.0_WP,0.0_WP,0.0_WP/)
                               
                                if(i_in .eq. j_in) then 
                                    ! If on main diagonal, we don't need to do a multi-cell averaging. 
                                    if(this%VF%vf(i,j,k) .gt. 1e-12 .and. this%VF%vf(i,j,k) .lt. 1.0_WP - 1e-12) then 
                                        call solveFaceParaboloid(solver,this%fs%sigma,P0,P1,P2,P3,this%vf%liquid_gas_interface(i,j,k),this%PressureOption,this%MarangoniOption,force)
                                    endif
                                    if(sum(force*face_normal) .lt. -1e-12) then
                                        force = -force
                                    endif
                                else
                                    ! If off diagonal, then we are touching four cells. 
                                    ! The values in the normal direction (shift_index_normal) are seperated
                                    ! The values in the momentum direction (shift_index_momentum) are grouped
                                    ! Get Appropriate Half Face (Positive)
                                    call getHalfFace(P0,P1,P2,P3,i_in,j_in,2,R0,R1,R2,R3)
                                    ! 0 Normal, 0 Momentum
                                    index = (/i,j,k/)
                                    if(this%VF%vf(index(1),index(2),index(3)) .gt. 1e-12 .and. this%VF%vf(index(1),index(2),index(3)) .lt. 1.0_WP - 1e-12) then 
                                        call solveFaceParaboloid(solver,this%fs%sigma,R0,R1,R2,R3,this%vf%liquid_gas_interface(index(1),index(2),index(3)),this%PressureOption,this%MarangoniOption,force1)
                                    endif
                                    if(sum(force1*face_normal) .lt. -1e-12) then
                                        force1 = -force1
                                    endif
                                    ! 1 Normal, 0 Momentum
                                    index = (/i,j,k/)
                                    index = index + shift_index_normal
                                    if(this%VF%vf(index(1),index(2),index(3)) .gt. 1e-12 .and. this%VF%vf(index(1),index(2),index(3)) .lt. 1.0_WP - 1e-12) then 
                                        call solveFaceParaboloid(solver,this%fs%sigma,R0,R1,R2,R3,this%vf%liquid_gas_interface(index(1),index(2),index(3)),this%PressureOption,this%MarangoniOption,force3)
                                    endif
                                    if(sum(force3*face_normal) .lt. -1e-12) then
                                        force3 = -force3
                                    endif
                                    ! Get other half face (Negative)
                                     call getHalfFace(P0,P1,P2,P3,i_in,j_in,1,R0,R1,R2,R3)
                                    ! 0 Normal, 1 Momentum
                                    index = (/i,j,k/)
                                    index = index - shift_index_momentum
                                    if(this%VF%vf(index(1),index(2),index(3)) .gt. 1e-12 .and. this%VF%vf(index(1),index(2),index(3)) .lt. 1.0_WP - 1e-12) then 
                                        call solveFaceParaboloid(solver,this%fs%sigma,R0,R1,R2,R3,this%vf%liquid_gas_interface(index(1),index(2),index(3)),this%PressureOption,this%MarangoniOption,force2)
                                    endif
                                    if(sum(force2*face_normal) .lt. -1e-12) then
                                        force2 = -force2
                                    endif
                                    ! 1 Normal, 1 Momentum 
                                    index = (/i,j,k/)
                                    index = index + shift_index_normal - shift_index_momentum
                                    if(this%VF%vf(index(1),index(2),index(3)) .gt. 1e-12 .and. this%VF%vf(index(1),index(2),index(3)) .lt. 1.0_WP - 1e-12) then 
                                        call solveFaceParaboloid(solver,this%fs%sigma,R0,R1,R2,R3,this%vf%liquid_gas_interface(index(1),index(2),index(3)),this%PressureOption,this%MarangoniOption,force4)
                                    endif
                                    if(sum(force4*face_normal) .lt. -1e-12) then
                                        force4 = -force4
                                    endif

                                    ! Now we have all the forces, so sum them properly 
                                    ! Combine momentum direction
                                    force1 = force1+force2
                                    force3 = force3+force4
                                    ! Average Normal Direction
                                    ! force = (force1 + force3 + force2 + force4)/4.0_WP
                                    if(sum(force1**2) .gt. 1e-11 .and. sum(force3**2) .gt. 1e-11) then 
                                        force = (force1 + force3)/2.0_WP
                                    elseif(sum(force1**2) .gt. 1e-12 .or. sum(force3**2) .gt. 1e-12) then
                                        force = force1 + force3
                                    else 
                                        force = (/0.0_WP,0.0_WP,0.0_WP/)
                                    endif
                                endif
                                
                            CASE DEFAULT
                                call solveFace(solver,this%fs%sigma,P0,P1,P2,P3,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force)
                                if(sum(force*face_normal) .lt. -1e-12) then
                                    force = -force
                                endif
                                this%sigma_3D_OLD(i,j,k,i_in,j_in) = force(i_in)
                        END SELECT
                        ! endif
                        ! print *, force 
                        if(sum(force*face_normal) .lt. -1e-12) then
                            force = -force
                        endif
                        this%sigma_3D(i,j,k,i_in,j_in) = force(i_in)

                        force4 = (/0.0_WP,0.0_WP,0.0_WP/)
                        call solveFace(solver,this%fs%sigma,P0,P1,P2,P3,this%PU_spread*dx,this%PressureOption,this%MarangoniOption,force4)
                        this%sigma_3D_OLD(i,j,k,i_in,j_in) = force4(i_in)
                        this%dSigma(i,j,k,i_in,j_in) = this%sigma_3D(i,j,k,i_in,j_in)-this%sigma_3D_OLD(i,j,k,i_in,j_in)
                        ! Do Prints here 

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
                ! STOP
        end do
        end do
    end do
    ! Two Dimensional Adjustment
    if(this%TwoD) then 
        this%sigma_3D(:,:,:,:,3) = 0.0_WP
        this%sigma_3D(:,:,:,3,:) = 0.0_WP 
    endif
    ! Boundary Conditions 
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,1,1))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,1,2))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,1,3))

    call this%vf%cfg%sync(this%sigma_3D(:,:,:,2,1))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,2,2))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,2,3))

    call this%vf%cfg%sync(this%sigma_3D(:,:,:,3,1))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,3,2))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,3,3))

    ! Boundary Conditions 
    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,1,1))
    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,1,2))
    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,1,3))

    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,2,1))
    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,2,2))
    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,2,3))

    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,3,1))
    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,3,2))
    call this%vf%cfg%sync(this%sigma_3D_OLD(:,:,:,3,3))

    ! Boundary Conditions 
    call this%vf%cfg%sync(this%dSigma(:,:,:,1,1))
    call this%vf%cfg%sync(this%dSigma(:,:,:,1,2))
    call this%vf%cfg%sync(this%dSigma(:,:,:,1,3))

    call this%vf%cfg%sync(this%dSigma(:,:,:,2,1))
    call this%vf%cfg%sync(this%dSigma(:,:,:,2,2))
    call this%vf%cfg%sync(this%dSigma(:,:,:,2,3))

    call this%vf%cfg%sync(this%dSigma(:,:,:,3,1))
    call this%vf%cfg%sync(this%dSigma(:,:,:,3,2))
    call this%vf%cfg%sync(this%dSigma(:,:,:,3,3))

end subroutine updateSurfaceTensionStresses3D

subroutine updateSurfaceTensionStresses3DFaceParaboloids(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class, only: solveFace,solveFaceParaboloid
    use mathtools, only: normalize
    
    implicit none 
    class(conservative_st_type), intent(inout) :: this
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver
    type(SeparatorVariant_type) :: interface_component 
    ! Geometry and Forces
    integer, dimension(1:3) :: shift_index_normal,shift_index_momentum,index
    integer :: i,j,k,j_in,i_in,k_in
    integer :: shift_first_index,shift_second_index
    real(WP), dimension(1:3) :: cen, dvec, shift,face_normal
    real(WP), dimension(1:3) :: force
    real(WP), dimension(1:3) :: P0,P1,P2,P3,R0,R1,R2,R3 
    real(WP), dimension(1:3) :: pressure_cell_center, velocity_cell_center,face_center
    real(WP) :: dx,dy,dz
    ! Jibben Stuff
    integer :: ind,ii,jj,kk, icenter 
    integer, dimension(1:3) :: lower,upper,size_stencil
    real(WP), dimension(1:3) :: datum,normal,quarter_size
    real(WP) :: vf_sample,count,normal_mag,liquid_in_cell,cell_volume_total,vf_cell
    type(JibbenNeigh_type) :: jibben_neighborhood
    type(RectCub_type) :: cell
    logical :: found_interface 
    type(SepVM_type)   :: phase_moments
    real(WP) :: weight,distance,length_scale
    real(WP) :: cell_volume
    real(WP) :: liquid_volume
    real(WP) :: gas_volume
    real(WP) :: liquid_fraction
    real(WP) :: gas_fraction
    real(WP), dimension(1:3) :: N1,N2,N3,N4
    real(WP) :: W1,W2,W3,W4,VF1,VF2,VF3,VF4
    ! HF Paraboloid Stuff
    real(WP), dimension(-3:3,-3:3,-3:3) :: VF_Stencil
    real(WP), dimension(1:3) :: stencil_center 
    integer :: Nx,Ny,Nz,status

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    dz = this%vf%cfg%dz(1)
    dvec = (/dx,dy,dz/)
    shift = dvec
    ! Zero Out 
    this%sigma_3D = 0.0_WP 
    ! Create Neighborhood Empty
    call new(neighborhood) 
    call new(solver)
    call setNeighborhood(solver,neighborhood) 
    call new(interface_component)
    call new(phase_moments)
    ! print *, "Solver Neighbhorhood"
    ! Jibben Start
    call new(cell)
    call new(jibben_neighborhood)
    call reserve(jibben_neighborhood,48) ! 4x4x3 Stencil, with extra cells in offset direction
    ! print *, "Jibben Neighbhorhood"/
    ! Stress
    ! First, get the cell center of the cell we are in
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
    do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
    do i=this%fs%cfg%imin_,this%fs%cfg%imax_
        pressure_cell_center = (/this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
        do i_in = 1,3 
            ! Set up shift to velocity cell
            shift = (/0.0_WP,0.0_WP,0.0_WP/)  
            shift(i_in) = dvec(i_in)/2
            shift_index_momentum = (/0,0,0/)
            shift_index_momentum(i_in) = 1
            ! Move center to velocity cell we are looking at
            velocity_cell_center = pressure_cell_center - shift

            do j_in = 1,3 
                ! Shift from velocity cell center ot the face center we are looking at:
                shift = (/0.0_WP,0.0_WP,0.0_WP/)  
                shift_index_normal = (/0,0,0/)
                shift_index_normal(j_in) = 1
                face_normal = (/0.0_WP,0.0_WP,0.0_WP/)
                face_normal(j_in) = 1.0_WP
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
                ! print *, "Points"
                ! Now we need to get the paraboloid using the face centered formulation 
                ! Make Jibben Neighborhood 
                call setSize(jibben_neighborhood,0)
                ! print *, "SetSize"
                ind = 0 
                found_interface = .false. 
                ! Set neighborhood bounds
                lower = (/i-1, j-1, k-1/)
                upper = (/i+1, j+1, k+1/)
                if (i_in.ne.j_in) then
                    lower(i_in) = lower(i_in) - 1
                    upper(j_in) = upper(j_in) + 1
                endif
                do kk=lower(3),upper(3)
                do jj=lower(2),upper(2)
                do ii=lower(1),upper(1)
                    if(getNumberOfVertices(this%vf%interface_polygon(1,ii,jj,kk)) .gt. 0) then
                        vf_sample = this%VF%vf(ii,jj,kk)
                        ! VF Weight
                        weight = 1.0_WP 
                        if(vf_sample .lt. 0.1_WP) then 
                            weight = SIN(5.0_WP*Pi*vf_sample)**2
                        endif 
                        if(vf_sample .gt. 0.9_WP) then 
                            weight = SIN(5.0_WP*Pi*(1.0_WP-vf_sample))**2
                        endif 
                        call addMember(jibben_neighborhood,this%vf%interface_polygon(1,ii,jj,kk),weight)
                        if(ii.eq.i.and.jj.eq.j.and.kk.eq.k) then 
                            icenter=ind
                            call setCenterOfStencil(jibben_neighborhood,icenter)
                        endif
                        ind = ind+1
                    endif
                enddo
                enddo
                enddo
                ! print *, "Added Members" 
                ! Get Datum and Normal
                datum = face_center
                normal = (/0.0_WP,0.0_WP,0.0_WP/)
                found_interface = .false.
                count = 0.0_WP
                W1 = 0.0_WP 
                W2 = 0.0_WP 
                W3 = 0.0_WP 
                W4 = 0.0_WP 
                N1 = 0.0_WP 
                N2 = 0.0_WP
                N3 = 0.0_WP 
                N4 = 0.0_WP
                liquid_in_cell = 0.0_WP 
                cell_volume = 0.0_WP 
                vf_cell = 0.0_WP
                cell_volume_total = 0.0_WP
                ! Create a cell centered at the faceCenter 
                if(i_in .eq. j_in) then 
                    index = (/i,j,k/)
                    vf_sample = this%VF%vf(index(1),index(2),index(3))
                    vf_cell = vf_sample
                    if(vf_sample .gt. 1e-12 .and. vf_sample .lt. 1.0_WP - 1e-12) then 
                        normal = calculateNormal(this%vf%interface_polygon(1,i,j,k))
                        datum = calculateCentroid( &
                            this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        found_interface = .true. 
                    endif
                    count = 1.0_WP
                    W1 = 1.0_WP 
                    ! print *, "Eq Branch"
                else
                    ! ========================== First Pressure Cell
                    index = (/i,j,k/)
                    vf_sample = this%VF%vf(index(1),index(2),index(3))
                    VF1 = vf_sample
                    ! Construct quarter cell
                    ! Shift to quarter cell center
                    shift = 0.0_WP
                    shift(i_in) = 0.25_WP*dvec(i_in)
                    shift(j_in) = -0.25_WP*dvec(j_in)
                    ! Cell Size
                    quarter_size = dvec
                    quarter_size(i_in) = 0.5_WP*dvec(i_in)
                    quarter_size(j_in) = 0.5_WP*dvec(j_in)
                    ! Construct
                    call construct_2pt(cell, &
                        face_center + shift - 0.5_WP*quarter_size, &
                        face_center + shift + 0.5_WP*quarter_size)
                    ! Calculate total quarter-cell volume
                    cell_volume = calculateVolume(cell)
                    cell_volume_total = cell_volume_total + cell_volume
                    ! Other Math
                    if (vf_sample .ge. 1.0_WP - 1e-12_WP) then
                        ! Full liquid cell
                        liquid_in_cell = liquid_in_cell + cell_volume
                    else if (vf_sample .gt. 1e-12_WP) then
                        ! Weight
                        ! Surface Area Weight
                        ! weight = this%vf%SD(index(1),index(2),index(3))
                        ! VF Weight
                        weight = 1.0_WP 
                        if(vf_sample .lt. 0.1_WP) then 
                            weight = SIN(5.0_WP*Pi*vf_sample)**2
                        endif 
                        if(vf_sample .gt. 0.9_WP) then 
                            weight = SIN(5.0_WP*Pi*(1.0_WP-vf_sample))**2
                        endif 
                        !Distance Weights
                        ! Calculate interface centroid
                        ! cen = calculateCentroid( &
                        !     this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        ! ! Distance from momentum face center to interface centroid
                        ! distance = sqrt(sum((face_center - cen)**2))
                        ! ! Characteristic length scale
                        ! length_scale = 0.5_WP * minval(dvec)
                        ! ! Inverse-distance-squared weighting (bounded)
                        ! weight = 1.0_WP / (1.0_WP + (distance/length_scale)**2)
                        ! No Weight
                        ! weight = 1.0_WP
                        N1 = calculateNormal(this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        W1 = weight
                        normal = normal + weight*calculateNormal( &
                            this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        call getNormMoments(cell, &
                            this%vf%liquid_gas_interface(index(1),index(2),index(3)), &
                            phase_moments)
                        liquid_volume = getVolume(phase_moments, 0)
                        liquid_in_cell = liquid_in_cell + liquid_volume
                        count = count + weight
                        found_interface = .true.
                    endif
                    ! ========================== Second Pressure Cell
                    index = (/i,j,k/)
                    index = index+shift_index_normal
                    vf_sample = this%VF%vf(index(1),index(2),index(3))
                    VF2 = vf_sample
                    ! Construct quarter cell
                    ! Shift to quarter cell center
                    shift = 0.0_WP
                    shift(i_in) = 0.25_WP*dvec(i_in)
                    shift(j_in) = 0.25_WP*dvec(j_in)
                    ! Cell Size
                    quarter_size = dvec
                    quarter_size(i_in) = 0.5_WP*dvec(i_in)
                    quarter_size(j_in) = 0.5_WP*dvec(j_in)
                    ! Construct
                    call construct_2pt(cell, &
                        face_center + shift - 0.5_WP*quarter_size, &
                        face_center + shift + 0.5_WP*quarter_size)
                    ! Calculate total quarter-cell volume
                    cell_volume = calculateVolume(cell)
                    cell_volume_total = cell_volume_total + cell_volume
                    
                    ! Other Math
                    if (vf_sample .ge. 1.0_WP - 1e-12_WP) then
                        ! Full liquid cell
                        liquid_in_cell = liquid_in_cell + cell_volume
                    else if (vf_sample .gt. 1e-12_WP) then
                        ! Weight
                        ! Surface Area Weight
                        ! weight = this%vf%SD(index(1),index(2),index(3))
                        ! VF Weight
                        weight = 1.0_WP 
                        if(vf_sample .lt. 0.1_WP) then 
                            weight = SIN(5.0_WP*Pi*vf_sample)**2
                        endif 
                        if(vf_sample .gt. 0.9_WP) then 
                            weight = SIN(5.0_WP*Pi*(1.0_WP-vf_sample))**2
                        endif 
                        !Distance Weights
                        ! Calculate interface centroid
                        ! cen = calculateCentroid( &
                        !     this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        ! ! Distance from momentum face center to interface centroid
                        ! distance = sqrt(sum((face_center - cen)**2))
                        ! ! Characteristic length scale
                        ! length_scale = 0.5_WP * minval(dvec)
                        ! ! Inverse-distance-squared weighting (bounded)
                        ! weight = 1.0_WP / (1.0_WP + (distance/length_scale)**2)
                        ! No Weight
                        ! weight = 1.0_WP
                        N2 = calculateNormal(this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        W2 = weight
                        normal = normal + weight*calculateNormal( &
                            this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        call getNormMoments(cell, &
                            this%vf%liquid_gas_interface(index(1),index(2),index(3)), &
                            phase_moments)
                        liquid_volume = getVolume(phase_moments, 0)
                        liquid_in_cell = liquid_in_cell + liquid_volume
                        count = count + weight
                        found_interface = .true.
                    endif
                    ! ========================== Third Pressure Cell
                    index = (/i,j,k/)
                    index = index + shift_index_normal - shift_index_momentum 
                    vf_sample = this%VF%vf(index(1),index(2),index(3))
                    VF3 = vf_sample
                    ! Construct quarter cell
                    ! Shift to quarter cell center
                    shift = 0.0_WP
                    shift(i_in) = -0.25_WP*dvec(i_in)
                    shift(j_in) = 0.25_WP*dvec(j_in)
                    ! Cell Size
                    quarter_size = dvec
                    quarter_size(i_in) = 0.5_WP*dvec(i_in)
                    quarter_size(j_in) = 0.5_WP*dvec(j_in)
                    ! Construct
                    call construct_2pt(cell, &
                        face_center + shift - 0.5_WP*quarter_size, &
                        face_center + shift + 0.5_WP*quarter_size)
                    ! Calculate total quarter-cell volume
                    cell_volume = calculateVolume(cell)
                    cell_volume_total = cell_volume_total + cell_volume
                    
                    ! Other Math
                    if (vf_sample .ge. 1.0_WP - 1e-12_WP) then
                        ! Full liquid cell
                        liquid_in_cell = liquid_in_cell + cell_volume
                    else if (vf_sample .gt. 1e-12_WP) then
                        ! Weight
                        ! Surface Area Weight
                        ! weight = this%vf%SD(index(1),index(2),index(3))
                        ! VF Weight
                        weight = 1.0_WP 
                        if(vf_sample .lt. 0.1_WP) then 
                            weight = SIN(5.0_WP*Pi*vf_sample)**2
                        endif 
                        if(vf_sample .gt. 0.9_WP) then 
                            weight = SIN(5.0_WP*Pi*(1.0_WP-vf_sample))**2
                        endif 
                        !Distance Weights
                        ! Calculate interface centroid
                        ! cen = calculateCentroid( &
                        !     this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        ! ! Distance from momentum face center to interface centroid
                        ! distance = sqrt(sum((face_center - cen)**2))
                        ! ! Characteristic length scale
                        ! length_scale = 0.5_WP * minval(dvec)
                        ! ! Inverse-distance-squared weighting (bounded)
                        ! weight = 1.0_WP / (1.0_WP + (distance/length_scale)**2)
                        ! No Weight
                        ! weight = 1.0_WP
                        N3 = calculateNormal(this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        W3 = weight
                        normal = normal + weight*calculateNormal( &
                            this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        call getNormMoments(cell, &
                            this%vf%liquid_gas_interface(index(1),index(2),index(3)), &
                            phase_moments)
                        liquid_volume = getVolume(phase_moments, 0)
                        liquid_in_cell = liquid_in_cell + liquid_volume
                        count = count + weight
                        found_interface = .true.
                    endif
                    ! ========================== Fourth Pressure Cell
                    index = (/i,j,k/)
                    index = index - shift_index_momentum 
                    vf_sample = this%VF%vf(index(1),index(2),index(3))
                    VF4 = vf_sample
                    ! Construct quarter cell
                    ! Shift to quarter cell center
                    shift = 0.0_WP
                    shift(i_in) = -0.25_WP*dvec(i_in)
                    shift(j_in) = -0.25_WP*dvec(j_in)
                    ! Cell Size
                    quarter_size = dvec
                    quarter_size(i_in) = 0.5_WP*dvec(i_in)
                    quarter_size(j_in) = 0.5_WP*dvec(j_in)
                    ! Construct
                    call construct_2pt(cell, &
                        face_center + shift - 0.5_WP*quarter_size, &
                        face_center + shift + 0.5_WP*quarter_size)
                    ! Calculate total quarter-cell volume
                    cell_volume = calculateVolume(cell)
                    cell_volume_total = cell_volume_total + cell_volume
                    
                    ! Other Math
                    if (vf_sample .ge. 1.0_WP - 1e-12_WP) then
                        ! Full liquid cell
                        liquid_in_cell = liquid_in_cell + cell_volume
                    else if (vf_sample .gt. 1e-12_WP) then
                        ! Weight
                        ! Surface Area Weight
                        ! weight = this%vf%SD(index(1),index(2),index(3))
                        ! VF Weight
                        weight = 1.0_WP 
                        if(vf_sample .lt. 0.1_WP) then 
                            weight = SIN(5.0_WP*Pi*vf_sample)**2
                        endif 
                        if(vf_sample .gt. 0.9_WP) then 
                            weight = SIN(5.0_WP*Pi*(1.0_WP-vf_sample))**2
                        endif 
                        !Distance Weights
                        ! Calculate interface centroid
                        ! cen = calculateCentroid( & 
                        !     this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        ! ! Distance from momentum face center to interface centroid
                        ! distance = sqrt(sum((face_center - cen)**2))
                        ! ! Characteristic length scale
                        ! length_scale = 0.5_WP * minval(dvec)
                        ! ! Inverse-distance-squared weighting (bounded) 
                        ! weight = 1.0_WP / (1.0_WP + (distance/length_scale)**2)
                        ! No Weight
                        ! weight = 1.0_WP
                        N4 = calculateNormal(this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        W4 = weight
                        normal = normal + weight*calculateNormal( &
                            this%vf%interface_polygon(1,index(1),index(2),index(3)))
                        call getNormMoments(cell, &
                            this%vf%liquid_gas_interface(index(1),index(2),index(3)), &
                            phase_moments)
                        liquid_volume = getVolume(phase_moments, 0)
                        liquid_in_cell = liquid_in_cell + liquid_volume
                        count = count + weight
                        found_interface = .true.
                    endif

                    
                    if (cell_volume_total .gt. 0.0_WP) then
                        vf_cell = liquid_in_cell / cell_volume_total
                    else
                        vf_cell = 0.0_WP
                    endif
                    ! print *, "Not EQ Branch" 
                endif
                if (.not.ieee_is_finite(vf_cell)) then
                    print *, "ERROR: Nonfinite momentum VF"
                    print *, "Cell:",i,j,k
                    print *, "Directions:",i_in,j_in
                    error stop
                endif

                if (vf_cell.lt.-1e-8_WP.or.vf_cell.gt.1.0_WP+1e-8_WP) then
                    print *, "ERROR: Momentum VF outside [0,1]"
                    print *, "Cell:",i,j,k
                    print *, "Directions:",i_in,j_in
                    print *, "VF:",vf_cell
                    print *, "Liquid volume:",liquid_in_cell
                    print *, "Total volume:",cell_volume_total
                    error stop
                endif

                vf_cell = min(max(vf_cell,0.0_WP),1.0_WP)
                ! count = W1 + W2 + W3 + W4
                if(count .lt. 0.0) then 
                    count = -count
                endif
                if(found_interface) then 
                    normal = normal/count 
                    normal_mag = sqrt(sum(normal**2))
                    if(normal_mag .gt. 1e-12) then 
                        normal = normal/normal_mag 
                    else 
                        print *, "====================" 
                        print *, "SMALL NORMAL MAG"
                        print *, "N1,W1 = ",N1,W1,VF1
                        print *, "N2,W2 = ",N2,W2,VF2
                        print *, "N3,W3 = ",N3,W3,VF3
                        print *, "N4,W4 = ",N4,W4,VF4
                        print *, "Normal = ",normal
                        print *, "Count = ", count
                        print *, "===================="
                        normal = 0.0_WP
                    endif
                else
                    normal = (/0.0_WP,0.0_WP,0.0_WP/)
                endif
                if(this%NormalOption .eq. 2) then 
                    ! Use Height Function Normal
                    ! Update Upper and lower for the proper stencils 
                    lower = (/-3,-3,-3/)
                    upper = (/3,3,3/)
                    size_stencil = (/7,7,7/)
                    if(i_in .ne. j_in) then 
                        ! In the momentum direction, we take one off the top
                        size_stencil(i_in) = size_stencil(i_in)-1
                        ! In the Normal direction, take one off the bottom 
                        size_stencil(j_in) = size_stencil(j_in)-1
                    endif
                    ! Get VF
                    do kk=lower(3),upper(3)
                    do jj=lower(2),upper(2)
                    do ii=lower(1),upper(1)
                        VF_Stencil(ii,jj,kk) = this%vf%VF(i+ii,j+jj,k+kk)
                    enddo
                    enddo
                    enddo
                    ! Remake Paraboloid
                    Nx = size_stencil(1)
                    Ny = size_stencil(2)
                    Nz = size_stencil(3)
                endif
                ! Compute Force 
                ! print *,"computeForce"  
                force = (/0.0_WP,0.0_WP,0.0_WP/)
                call construct_2pt(cell, &
                        face_center - 0.5_WP*dvec, &
                        face_center + 0.5_WP*dvec)

                if(found_interface .and. vf_cell .gt. 1e-12 .and. vf_cell .lt. 1.0_WP - 1e-12) then 

                    SELECT CASE(this%ParaboloidOption)
                    CASE(1)
                        ! print *,"PRE NORMAL"
                        if(this%NormalOption .eq. 2) then 
                            call HF3D_LVIRA_Normal(VF_Stencil,i_in,j_in,i_in*j_in,normal)
                        endif
                        ! print *,"POST NORMAL",normal
                        ! Localize Jibben neighborhood
                        call setDelta(jibben_neighborhood, 2.5_WP*this%vf%cfg%meshsize(i,j,k))

                        call localize(jibben_neighborhood, datum, normal)

                        call reconstructJibben3D(jibben_neighborhood,interface_component)
                    CASE(2)
                        ! Use the Height Function Paraboloid

                        ! Update Upper and lower for the proper stencils 
                        lower = (/-3,-3,-3/)
                        upper = (/3,3,3/)
                        size_stencil = (/7,7,7/)
                        if(i_in .ne. j_in) then 
                            ! In the momentum direction, we take one off the top
                            size_stencil(i_in) = size_stencil(i_in)-1
                            ! In the Normal direction, take one off the bottom 
                            size_stencil(j_in) = size_stencil(j_in)-1
                        endif
                        ! Get VF
                        do kk=lower(3),upper(3)
                        do jj=lower(2),upper(2)
                        do ii=lower(1),upper(1)
                            VF_Stencil(ii,jj,kk) = this%vf%VF(i+ii,j+jj,k+kk)
                        enddo
                        enddo
                        enddo
                        ! Remake Paraboloid
                        Nx = size_stencil(1)
                        Ny = size_stencil(2)
                        Nz = size_stencil(3)
                        stencil_center = datum
                        status = -100 
                        call reconstructHF3D_Variant(VF_Stencil,stencil_center,dx,Nx,Ny,Nz,interface_component,status)
                    end select
                    ! We now have VF cell and
                    ! print *,"VF Cell=",vf_cell 
                    if(vf_cell .lt. 1e-12) then 
                        print *, "SMALL VF IN INTERFACE CELL"
                    endif
                    SELECT CASE(this%MatchingOption)
                        CASE(2)
                            ! print *,"MATCHING"
                            call matchVolumeFraction(cell,vf_cell,interface_component) 
                        CASE DEFAULT
                            ! print *,"NOT MATCHING"
                    END SELECT
    
                    call emptyNeighborhood(jibben_neighborhood)
                    ! Montiors
                    ! P
                    call solveFaceParaboloid(solver,this%fs%sigma,P0,P1,P2,P3, interface_component,1.0_WP,this%MarangoniOption,force)
                    this%sigma_3D_P(i,j,k,i_in,j_in) = force(i_in)
                    ! No P
                    call solveFaceParaboloid(solver,this%fs%sigma,P0,P1,P2,P3, interface_component,0.0_WP,this%MarangoniOption,force)
                    this%sigma_3D_NoP(i,j,k,i_in,j_in) = force(i_in)

                    this%sigma_3D_DP(i,j,k,i_in,j_in) = this%sigma_3D_P(i,j,k,i_in,j_in) - this%sigma_3D_NoP(i,j,k,i_in,j_in)
                    ! Real
                    call solveFaceParaboloid(solver,this%fs%sigma,P0,P1,P2,P3, interface_component,this%PressureOption,this%MarangoniOption,force)
                endif
                ! if(sum(force * face_normal) .lt. 0.0_WP) then 
                !     force = -force
                !     print *, "FLIPPING FORCE IN NGA"
                ! endif
                ! print *,"Force DOne"   
                this%sigma_3D(i,j,k,i_in,j_in) = force(i_in)
                if (ieee_is_nan(force(i_in))) then 
                    print *, "=========================="
                    print *, "VALUE IS NAN"
                    print *, "i,j,k = ", i,j,k
                    print *, "i_in,j_in = ", i_in,j_in 
                    print *, "VALUE = ", force(i_in)
                    print *, "vf_cell = ",vf_cell
                    print *, "normal = ", normal 
                    call printToScreen(interface_component)
                    call solveFaceParaboloid(solver,this%fs%sigma,P0,P1,P2,P3, interface_component,this%PressureOption,this%MarangoniOption,force,.true.)
                    print *, "=========================="
                endif
            enddo
        enddo
    enddo
    enddo
    enddo
    ! print *,"MethodDone"
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,1,1))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,1,2))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,1,3))

    call this%vf%cfg%sync(this%sigma_3D(:,:,:,2,1))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,2,2))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,2,3))

    call this%vf%cfg%sync(this%sigma_3D(:,:,:,3,1))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,3,2))
    call this%vf%cfg%sync(this%sigma_3D(:,:,:,3,3))

end subroutine updateSurfaceTensionStresses3DFaceParaboloids

subroutine updateSurfaceTensionStresses3DEllipsoid(this,A,v)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class, only: solveFace,solveFaceEllipsoid

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in,k_in,shift_first_index,shift_second_index ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver
    real(WP),dimension(3,3) :: A
    real(WP),dimension(3) :: v
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,force,corner,dvec,shift,forceB
    real(WP), dimension(1:3) :: P0,P1,P2,P3,pressure_cell_center,velocity_cell_center,face_center
    real(WP) :: dx,dy,dz
    real(WP) :: x_i,y_j,z_k

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    dz = this%vf%cfg%dz(1)
    dvec = (/dx,dy,dz/)
    shift = dvec
    ! Create Neighborhood and solver  
    call new(neighborhood) 
    call new(solver)
    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                
            ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then   
                ! Empty Neighborhood
                call emptyNeighborhood(neighborhood) 
                ! Set Neighborhood in solver
                call setNeighborhood(solver,neighborhood)
                ! ====== Get Stresses
                ! call printSolver(solver)
                ! We are going to loop over the stress components and calculate them
                x_i = this%fs%cfg%x(i)
                y_j = this%fs%cfg%y(j)
                z_k = this%fs%cfg%z(k)
                corner = (/x_i,y_j,z_k/)
                ! First, get the cell center of the cell we are in
                pressure_cell_center = (/this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
                do i_in = 1,3 ! Cell Loop, This tells us the force direction
                    ! Set up shift to velocity cell
                    shift = (/0.0_WP,0.0_WP,0.0_WP/)  
                    shift(i_in) = dvec(i_in)/2
                    ! Move center to velocity cell we are looking at
                    velocity_cell_center = pressure_cell_center - shift
                    do j_in = 1,3 ! Face Loop, This tells us the normal direction
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
                        
                        call solveFaceEllipsoid(solver,this%fs%sigma,P0,P1,P2,P3,A(:,1),A(:,2),A(:,3),v,this%PressureOption,this%MarangoniOption,force)
                       
                        ! endif
                        ! print *, force  
                        this%sigma_3D_Exact(i,j,k,i_in,j_in) = force(i_in)
                        ! Do Prints here 

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
                ! STOP
        end do
        end do
    end do
    if(this%TwoD) then 
        this%sigma_3D(:,:,:,:,3) = 0.0_WP
        this%sigma_3D(:,:,:,3,:) = 0.0_WP 
    endif
    ! Boundary Conditions 
    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,1,1))
    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,1,2))
    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,1,3))

    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,2,1))
    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,2,2))
    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,2,3))

    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,3,1))
    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,3,2))
    call this%vf%cfg%sync(this%sigma_3D_Exact(:,:,:,3,3))
    ! print *, "full sum stress",sqrt(sum(this%sigma_3D_Exact**2))  
end subroutine updateSurfaceTensionStresses3DEllipsoid

subroutine updateSurfaceTensionForces3D(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k 
    real(WP) :: dxi,dyi,dzi
    real(WP) :: Force_Mag
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            dxi = this%fs%cfg%dxi(i)
            dyi = this%fs%cfg%dyi(j)
            dzi = this%fs%cfg%dzi(k) 
            ! X Component
            this%Fst_x_3D(i,j,k) = (this%sigma_3D(i,j,k,1,1)-this%sigma_3D(i-1,j,k,1,1)) + &
                                (this%sigma_3D(i,j,k,1,2)-this%sigma_3D(i,j-1,k,1,2)) 
            if(.not. this%TwoD) then 
                this%Fst_x_3D(i,j,k) = this%Fst_x_3D(i,j,k) + &
                                (this%sigma_3D(i,j,k,1,3)-this%sigma_3D(i,j,k-1,1,3))
            endif

            this%Fst_x_3D(i,j,k) = this%Fst_x_3D(i,j,k) *dxi*dyi*dzi
            ! Y Component
            this%Fst_y_3D(i,j,k) = (this%sigma_3D(i,j,k,2,1)-this%sigma_3D(i-1,j,k,2,1)) + &
                                (this%sigma_3D(i,j,k,2,2)-this%sigma_3D(i,j-1,k,2,2))
            if(.not. this%TwoD) then 
                this%Fst_y_3D(i,j,k) = this%Fst_y_3D(i,j,k) + &
                                (this%sigma_3D(i,j,k,2,3)-this%sigma_3D(i,j,k-1,2,3))
            endif
                                
            this%Fst_y_3D(i,j,k) = this%Fst_y_3D(i,j,k)*dxi*dyi*dzi
            ! Z Component
            this%Fst_z_3D(i,j,k) = 0.0_WP 
            if(.not. this%TwoD) then 
                this%Fst_z_3D(i,j,k) = (this%sigma_3D(i,j,k,3,1)-this%sigma_3D(i-1,j,k,3,1)) + &
                                (this%sigma_3D(i,j,k,3,2)-this%sigma_3D(i,j-1,k,3,2)) + &
                                (this%sigma_3D(i,j,k,3,3)-this%sigma_3D(i,j,k-1,3,3))
            endif
            this%Fst_z_3D(i,j,k) = this%Fst_z_3D(i,j,k)*dxi*dyi*dzi

            Force_Mag = sqrt(this%Fst_x_3D(i,j,k)**2 + this%Fst_y_3D(i,j,k)**2 + this%Fst_z_3D(i,j,k)**2)
            this%Fst_mag(i,j,k) = Force_Mag
            ! if(Force_Mag .gt. 4e3) then 
            !     print *, "Large Force at i,j,k = ",i,j,k,Force_Mag
            ! endif
            ! OLD FOR VIZ
            ! X Component
            this%Fst_x_3D_OLD(i,j,k) = (this%sigma_3D_OLD(i,j,k,1,1)-this%sigma_3D_OLD(i-1,j,k,1,1)) + &
                                (this%sigma_3D_OLD(i,j,k,1,2)-this%sigma_3D_OLD(i,j-1,k,1,2)) 
            if(.not. this%TwoD) then 
                this%Fst_x_3D_OLD(i,j,k) = this%Fst_x_3D_OLD(i,j,k) + &
                                (this%sigma_3D_OLD(i,j,k,1,3)-this%sigma_3D_OLD(i,j,k-1,1,3))
            endif

            this%Fst_x_3D_OLD(i,j,k) = this%Fst_x_3D_OLD(i,j,k)*dxi*dyi*dzi
            ! Y Component
            this%Fst_y_3D_OLD(i,j,k) = (this%sigma_3D_OLD(i,j,k,2,1)-this%sigma_3D_OLD(i-1,j,k,2,1)) + &
                                (this%sigma_3D_OLD(i,j,k,2,2)-this%sigma_3D_OLD(i,j-1,k,2,2))
            if(.not. this%TwoD) then 
                this%Fst_y_3D_OLD(i,j,k) = this%Fst_y_3D_OLD(i,j,k) + &
                                (this%sigma_3D_OLD(i,j,k,2,3)-this%sigma_3D_OLD(i,j,k-1,2,3))
            endif
                                
            this%Fst_y_3D_OLD(i,j,k) = this%Fst_y_3D_OLD(i,j,k)*dxi*dyi*dzi
            ! Z Component
            this%Fst_z_3D_OLD(i,j,k) = 0.0_WP 
            if(.not. this%TwoD) then 
                this%Fst_z_3D_OLD(i,j,k) = (this%sigma_3D_OLD(i,j,k,3,1)-this%sigma_3D_OLD(i-1,j,k,3,1)) + &
                                (this%sigma_3D_OLD(i,j,k,3,2)-this%sigma_3D_OLD(i,j-1,k,3,2)) + &
                                (this%sigma_3D_OLD(i,j,k,3,3)-this%sigma_3D_OLD(i,j,k-1,3,3))
            endif
            this%Fst_z_3D_OLD(i,j,k) = this%Fst_z_3D_OLD(i,j,k)*dxi*dyi*dzi
        end do
        end do
    end do
    ! STOP
    ! Here, we would need to do assignment to Fst_xyz to actually use in the operators.   
    this%Fst_x_3D_diff = this%Fst_x_3D - this%Fst_x_3D_OLD
    this%Fst_y_3D_diff = this%Fst_y_3D - this%Fst_y_3D_OLD
    this%Fst_z_3D_diff = this%Fst_z_3D - this%Fst_z_3D_OLD

    call this%vf%cfg%sync(this%Fst_z_3D(:,:,:))
    call this%vf%cfg%sync(this%Fst_y_3D(:,:,:))
    call this%vf%cfg%sync(this%Fst_x_3D(:,:,:))

    call this%vf%cfg%sync(this%Fst_z_3D_OLD(:,:,:))
    call this%vf%cfg%sync(this%Fst_y_3D_OLD(:,:,:))
    call this%vf%cfg%sync(this%Fst_x_3D_OLD(:,:,:))

    call this%vf%cfg%sync(this%Fst_z_3D_diff(:,:,:))
    call this%vf%cfg%sync(this%Fst_y_3D_diff(:,:,:))
    call this%vf%cfg%sync(this%Fst_x_3D_diff(:,:,:))
end subroutine updateSurfaceTensionForces3D


subroutine updateSurfaceTensionForces3DEllipsoid(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k 
    real(WP) :: dxi,dyi,dzi
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            dxi = this%fs%cfg%dxi(i)
            dyi = this%fs%cfg%dyi(j)
            dzi = this%fs%cfg%dzi(k) 
            ! X Component
            this%Fst_x_3D_Exact(i,j,k) = (this%sigma_3D_Exact(i,j,k,1,1)-this%sigma_3D_Exact(i-1,j,k,1,1)) + &
                                (this%sigma_3D_Exact(i,j,k,1,2)-this%sigma_3D_Exact(i,j-1,k,1,2)) + &
                                (this%sigma_3D_Exact(i,j,k,1,3)-this%sigma_3D_Exact(i,j,k-1,1,3))
            this%Fst_x_3D_Exact(i,j,k) = this%Fst_x_3D_Exact(i,j,k) *dxi*dyi*dzi! This assumes uniform mesh
            ! Y Component
            this%Fst_y_3D_Exact(i,j,k) = (this%sigma_3D_Exact(i,j,k,2,1)-this%sigma_3D_Exact(i-1,j,k,2,1)) + &
                                (this%sigma_3D_Exact(i,j,k,2,2)-this%sigma_3D_Exact(i,j-1,k,2,2)) + &
                                (this%sigma_3D_Exact(i,j,k,2,3)-this%sigma_3D_Exact(i,j,k-1,2,3))
            this%Fst_y_3D_Exact(i,j,k) = this%Fst_y_3D_Exact(i,j,k) *dxi*dyi*dzi
            ! Z Component
            this%Fst_z_3D_Exact(i,j,k) = (this%sigma_3D_Exact(i,j,k,3,1)-this%sigma_3D_Exact(i-1,j,k,3,1)) + &
                                (this%sigma_3D_Exact(i,j,k,3,2)-this%sigma_3D_Exact(i,j-1,k,3,2)) + &
                                (this%sigma_3D_Exact(i,j,k,3,3)-this%sigma_3D_Exact(i,j,k-1,3,3))
            this%Fst_z_3D_Exact(i,j,k) = this%Fst_z_3D_Exact(i,j,k) *dxi*dyi*dzi
        end do
        end do
    end do
    ! Here, we would need to do assignment to Fst_xyz to actually use in the operators.  
    call this%vf%cfg%sync(this%Fst_z_3D_Exact(:,:,:))
    call this%vf%cfg%sync(this%Fst_y_3D_Exact(:,:,:))
    call this%vf%cfg%sync(this%Fst_x_3D_Exact(:,:,:))
end subroutine updateSurfaceTensionForces3DEllipsoid


subroutine updateSurfaceTensionForces(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k 

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            ! write(*,'(A,F10.5)') 'DX,i,j,k: ', cfg%dxmi(i)
            this%Fst_x(i,j,k) = (this%sigma_xx(i,j,k)-this%sigma_xx(i-1,j,k)) - (this%sigma_xy(i,j+1,k)-this%sigma_xy(i,j,k))
            this%Fst_x(i,j,k) = this%Fst_x(i,j,k) * this%fs%cfg%dxi(i)

            this%Fst_y(i,j,k) = (this%sigma_yy(i,j,k)-this%sigma_yy(i,j-1,k)) - (this%sigma_yx(i+1,j,k)-this%sigma_yx(i,j,k))
            this%Fst_y(i,j,k) = this%Fst_y(i,j,k) * this%fs%cfg%dyi(j)
            this%Fst_z(i,j,k) = 0.0_WP
        end do
        end do
    end do

    this%Pjx_ST = this%Fst_x
    this%Pjy_ST = this%Fst_y
    this%Pjz_ST = this%Fst_z
end subroutine updateSurfaceTensionForces

! ####################################################### SMOOTHING FUNCTIONS ################################################
subroutine applyLaplacianSmoothing(this)
    implicit none

    class(conservative_st_type), intent(inout) :: this

    integer :: i,j,k
    real(WP), parameter :: alpha = 0.5_WP

    ! Temporary arrays
    this%Fst_x = this%fs%Pjx
    this%Fst_y = this%fs%Pjy
    this%Fst_z = this%fs%Pjz

    ! ============================================================
    ! 3D Laplacian smoothing
    ! ============================================================

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_

                ! X component
                this%Fst_x(i,j,k) = &
                    (1.0_WP-alpha)*this%fs%Pjx(i,j,k) + &
                    (alpha/6.0_WP)*( &
                        this%fs%Pjx(i+1,j,k) + &
                        this%fs%Pjx(i-1,j,k) + &
                        this%fs%Pjx(i,j+1,k) + &
                        this%fs%Pjx(i,j-1,k) + &
                        this%fs%Pjx(i,j,k+1) + &
                        this%fs%Pjx(i,j,k-1))

                ! Y component
                this%Fst_y(i,j,k) = &
                    (1.0_WP-alpha)*this%fs%Pjy(i,j,k) + &
                    (alpha/6.0_WP)*( &
                        this%fs%Pjy(i+1,j,k) + &
                        this%fs%Pjy(i-1,j,k) + &
                        this%fs%Pjy(i,j+1,k) + &
                        this%fs%Pjy(i,j-1,k) + &
                        this%fs%Pjy(i,j,k+1) + &
                        this%fs%Pjy(i,j,k-1))

                ! Z component
                this%Fst_z(i,j,k) = &
                    (1.0_WP-alpha)*this%fs%Pjz(i,j,k) + &
                    (alpha/6.0_WP)*( &
                        this%fs%Pjz(i+1,j,k) + &
                        this%fs%Pjz(i-1,j,k) + &
                        this%fs%Pjz(i,j+1,k) + &
                        this%fs%Pjz(i,j-1,k) + &
                        this%fs%Pjz(i,j,k+1) + &
                        this%fs%Pjz(i,j,k-1))

            end do
        end do
    end do

    ! ============================================================
    ! Update surface tension arrays
    ! ============================================================

    this%fs%Pjx = this%Fst_x
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z

end subroutine applyLaplacianSmoothing

subroutine applyGradientSmoothing(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,ii,jj,kk,negativeIndex,positiveIndex
    real(WP) :: tot_grad_x,tot_grad_y,tot_force

    this%Fst_x = 0.0_WP;
    this%Fst_y = 0.0_WP;
    this%Fst_z = 0.0_WP;
    
    negativeIndex = -1
    positiveIndex = 1
    ! Calculate Gradients   
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%grad_vf_x(i,j,k) = abs((this%vf%VF(i  ,j  ,k) - this%vf%VF(i-1,j  ,k))  * this%vf%cfg%dxi(i))
            this%grad_vf_y(i,j,k) = abs((this%vf%VF(i  ,j  ,k) - this%vf%VF(i  ,j-1,k)) * this%vf%cfg%dyi(j))
            ! While we are looping everything, also check to see if everything has been processed
            if(abs(this%fs%Pjx(i,j,k)) .gt. 1e-10) then
                this%x_smoothing_tracker(i,j,k) = 1 ! Mark as needing smoothing
            endif

            if(abs(this%fs%Pjy(i,j,k)) .gt. 1e-10) then
                this%y_smoothing_tracker(i,j,k) = 1 ! Mark as needing smoothing
            endif

        end do
        end do
    end do
    ! Smoothing
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            tot_grad_x = 0.0_WP
            tot_grad_y = 0.0_WP
            ! X Direction
            if(abs(this%fs%Pjx(i,j,k)) .gt. 1e-12 .and. this%x_smoothing_tracker(i,j,k) == 1) then  ! This checks to see if we are in a mixed cell that has not yet been smoothed
                ! Current Cell
                tot_grad_x = this%grad_vf_x(i,j,k)
                tot_force = this%fs%Pjx(i,j,k)
                ! Going Left
                do while(abs(this%fs%Pjx(i+negativeIndex,j,k)) .gt. 1e-12 .and. (i+negativeIndex) .ge. this%fs%cfg%imin_)
                    tot_grad_x = tot_grad_x + this%grad_vf_x(i+negativeIndex,j,k)
                    tot_force = tot_force + this%fs%Pjx(i+negativeIndex,j,k)
                    negativeIndex = negativeIndex -1
                end do
                ! Going Right
                do while(abs(this%fs%Pjx(i+positiveIndex,j,k)) .gt. 1e-12 .and. (i+positiveIndex) .le. this%fs%cfg%imax_)
                    tot_grad_x = tot_grad_x + this%grad_vf_x(i+positiveIndex,j,k)
                    tot_force = tot_force + this%fs%Pjx(i+positiveIndex,j,k)
                    positiveIndex = positiveIndex +1
                end do
                ! Now Reassign Force Values in between 
                do ii = negativeIndex,positiveIndex
                    this%Fst_x(i+ii,j,k) = this%grad_vf_x(i+ii,j,k)*tot_force/(tot_grad_x + 1e-12_WP) ! Summing Components
                    
                    ! Mark as Smoothed
                    this%x_smoothing_tracker(i+ii,j,k) = 0
                end do
                ! Reset Indicies
                negativeIndex = -1
                positiveIndex = 1
                tot_force = 0.0_WP
            endif
            
            ! Y Direction  
            if(abs(this%fs%Pjy(i,j,k)) .gt. 1e-12 .and. this%y_smoothing_tracker(i,j,k) == 1) then  
                ! Current Cell
                tot_grad_y = this%grad_vf_y(i,j,k)
                tot_force = this%fs%Pjy(i,j,k)
                ! Going Down
                do while(abs(this%fs%Pjy(i,j+negativeIndex,k)) .gt. 1e-12 .and. (j+negativeIndex) .ge. this%fs%cfg%jmin_)
                    tot_grad_y = tot_grad_y + this%grad_vf_y(i,j+negativeIndex,k)
                    tot_force = tot_force + this%fs%Pjy(i,j+negativeIndex,k)
                    negativeIndex = negativeIndex -1
                end do
                ! Going Up
                do while(abs(this%fs%Pjy(i,j+positiveIndex,k)) .gt. 1e-12 .and. (j+positiveIndex) .le. this%fs%cfg%jmax_)
                    tot_grad_y = tot_grad_y + this%grad_vf_y(i,j+positiveIndex,k)
                    tot_force = tot_force + this%fs%Pjy(i,j+positiveIndex,k)
                    positiveIndex = positiveIndex +1
                end do
                ! Now Reassign Force Values in between 
                do jj = negativeIndex,positiveIndex
                    this%Fst_y(i,j+jj,k) = this%grad_vf_y(i,j+jj,k)*tot_force/(tot_grad_y + 1e-12_WP) ! Summing Components
                    ! Mark as Smoothed
                    this%y_smoothing_tracker(i,j+jj,k) = 0
                end do
                ! Reset Indicies
                negativeIndex = -1
                positiveIndex = 1
                tot_force = 0.0_WP
            endif 

        end do 
        end do
    end do
    ! write(*,'(A)') 'Gradient Smoothing Complete'
    this%fs%Pjx = this%Fst_x 
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z
    
end subroutine applyGradientSmoothing

subroutine applyPoissonSmoothing(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class
    use hypre_str_class, only: pcg_smg,gmres

    implicit none
    class(conservative_st_type), intent(inout) :: this
    type(hypre_str) :: poisson_solver ! Poisson Solver Object

    integer :: i,j,k,ii,jj,kk
    real(WP) :: tot_grad_x,tot_grad_y


    ! Initialize Poisson Solver
    poisson_solver=hypre_str(cfg=cfg,name='Smoothing',method=gmres,nst=7)
    poisson_solver%maxit=this%fs%psolv%maxit
    poisson_solver%rcvg=this%fs%psolv%rcvg
    ! Setup Solver
    ! Set 7-pt stencil map for the pressure solver
    poisson_solver%stc(1,:)=[ 0, 0, 0]
    poisson_solver%stc(2,:)=[+1, 0, 0]
    poisson_solver%stc(3,:)=[-1, 0, 0]
    poisson_solver%stc(4,:)=[ 0,+1, 0]
    poisson_solver%stc(5,:)=[ 0,-1, 0]
    poisson_solver%stc(6,:)=[ 0, 0,+1]
    poisson_solver%stc(7,:)=[ 0, 0,-1]
    ! Set Up Laplacian Operator
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            ! Set Laplacian
            poisson_solver%opr(1,i,j,k)=this%fs%divp_x(1,i,j,k)*this%fs%divu_x(-1,i+1,j,k)+&
            &                       this%fs%divp_x(0,i,j,k)*this%fs%divu_x( 0,i  ,j,k)+&
            &                       this%fs%divp_y(1,i,j,k)*this%fs%divv_y(-1,i,j+1,k)+&
            &                       this%fs%divp_y(0,i,j,k)*this%fs%divv_y( 0,i,j  ,k)+&
            &                       this%fs%divp_z(1,i,j,k)*this%fs%divw_z(-1,i,j,k+1)+&
            &                       this%fs%divp_z(0,i,j,k)*this%fs%divw_z( 0,i,j,k  )
            poisson_solver%opr(2,i,j,k)=this%fs%divp_x(1,i,j,k)*this%fs%divu_x( 0,i+1,j,k)
            poisson_solver%opr(3,i,j,k)=this%fs%divp_x(0,i,j,k)*this%fs%divu_x(-1,i  ,j,k)
            poisson_solver%opr(4,i,j,k)=this%fs%divp_y(1,i,j,k)*this%fs%divv_y( 0,i,j+1,k)
            poisson_solver%opr(5,i,j,k)=this%fs%divp_y(0,i,j,k)*this%fs%divv_y(-1,i,j  ,k)
            poisson_solver%opr(6,i,j,k)=this%fs%divp_z(1,i,j,k)*this%fs%divw_z( 0,i,j,k+1)
            poisson_solver%opr(7,i,j,k)=this%fs%divp_z(0,i,j,k)*this%fs%divw_z(-1,i,j,k  )
            ! Scale it by the cell volume
            ! poisson_solver%opr(:,i,j,k)=-poisson_solver%opr(:,i,j,k)*this%cfg%vol(i,j,k) I don't think I need this right now.
        end do
        end do
    end do
    ! Initialize the Poisson solver
    call poisson_solver%init()
    call poisson_solver%setup()
    poisson_solver%sol=this%fs%Pjx
    ! We now need to calculate the right hand side, which will be the divergence of the force field
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = (sum(this%fs%divp_x(:,i,j,k)*this%fs%Pjx(i:i+1,j,k))&
                                        +sum(this%fs%divp_y(:,i,j,k)*this%fs%Pjy(i,j:j+1,k))) ! No z direction right now
        end do
        end do
    end do
    
    poisson_solver%rhs = this%SurfaceTensionDiv
    this%poisson_source = this%SurfaceTensionDiv
    ! Solve the Poisson Equation
    call poisson_solver%solve()
    this%force_potential_field = poisson_solver%sol
    ! Now we need to calculate the gradient of the smoothed pressure field to get the new forces
    ! Calculate Forces
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%Fst_x(i,j,k) = (poisson_solver%sol(i,j,k) - poisson_solver%sol(i-1,j,k)) * this%fs%cfg%dxi(i)
            this%Fst_y(i,j,k) = (poisson_solver%sol(i,j,k) - poisson_solver%sol(i,j-1,k)) * this%fs%cfg%dyi(j)
            this%Fst_z(i,j,k) = 0.0_WP
        end do
        end do
    end do   
    ! Calculate Magnitude of Fst_x
    ! write(*,'(A,F10.5)') 'Max Fst_x before assign: ', maxval(abs(Fst_x))
    ! Assign
    this%fs%Pjx = this%Fst_x 
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z
    call poisson_solver%destroy()
    ! write(*,'(A)') 'Poisson Smoothing Complete'
end subroutine applyPoissonSmoothing

subroutine applyPeskinSmoothing(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k, ii,jj,kk
    real(WP) :: x_curr,y_curr,z_curr,delta_val,sum
    this%Fst_x = 0.0_WP;
    this%Fst_y = 0.0_WP;
    this%Fst_z = 0.0_WP;
    
    ! Smoothing
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            x_curr = this%fs%cfg%xm(i)
            y_curr = this%fs%cfg%ym(j)
            z_curr = this%fs%cfg%zm(k)
            ! Visit Adjacent Cels and add force value times delta to current force
            ! print *, "Applying Peskin",delta_val
            sum = 0.0_WP
            do ii = -2,2
                do jj = -2,2
                    do kk = -2,2
                        call this%PeskinDelta(i+ii,j+jj,k+kk,x_curr,y_curr,z_curr,delta_val)
                        ! print *, "Applying Peskin",delta_val
                        this%Fst_x(i,j,k) = this%Fst_x(i,j,k) + delta_val * this%fs%Pjx(i+ii,j+jj,k+kk)
                        this%Fst_y(i,j,k) = this%Fst_y(i,j,k) + delta_val * this%fs%Pjy(i+ii,j+jj,k+kk)
                        this%Fst_z(i,j,k) = this%Fst_z(i,j,k) + delta_val * this%fs%Pjz(i+ii,j+jj,k+kk)
                        sum = sum + delta_val
                    enddo
                enddo
            enddo
            ! print *, "Sum Peskin",sum
        end do
        end do
    end do

    this%fs%Pjx = this%Fst_x 
    this%fs%Pjy = this%Fst_y
    this%fs%Pjz = this%Fst_z
end subroutine applyPeskinSmoothing


subroutine PeskinDeltaBase(this,radius,val)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class
    use hypre_str_class, only: pcg_smg,gmres

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(out) :: val
    real(WP) :: radius
    
    if(abs(radius) < 2*this%vf%cfg%dx(1)) then
        val = 1 + cos(Pi * radius / (2*this%vf%cfg%dx(1)))
    else
        val = 0
    endif
    val = val/4 ! /(4*this%vf%cfg%dx(1))
    
end subroutine PeskinDeltaBase

subroutine PeskinDelta(this,i,j,k,x,y,z,val)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class
    use hypre_str_class, only: pcg_smg,gmres

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(out) :: val
    real(WP) :: x,y,z,temp 
    integer :: i,j,k
    ! print *, "================PeskinDelta============="
    val = 0.0_WP
    call this%PeskinDeltaBase(x-this%fs%cfg%xm(i),temp)
    ! print *, "R1,d1 = ",x-i*this%fs%cfg%dx(1),temp
    val = val + temp 
    call this%PeskinDeltaBase(y-this%fs%cfg%ym(j),temp)
    ! print *, "R2,d2 = ",y-i*this%fs%cfg%dx(1),temp
    val = val * temp 
    call this%PeskinDeltaBase(z-this%fs%cfg%zm(k),temp)
    ! print *, "R2,d2 = ",y-i*this%fs%cfg%dx(1),temp
    val = val * temp 
end subroutine PeskinDelta

subroutine get_surface_tension_coefficient(this,x,y,stc)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(out) :: stc
    real(WP), intent(in) :: x,y

    stc = this%fs%sigma + this%MarangoniOption(3) * (this%MarangoniOption(1) * x + this%MarangoniOption(2) * y)
end subroutine get_surface_tension_coefficient

subroutine update_surface_tension_coefficients(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k
    this%ST_Coeff = 0.0_WP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            if(this%vf%VF(i,j,k) .gt. VFlo .and. this%vf%VF(i,j,k) .lt. VFhi) then
                call get_surface_tension_coefficient(this,this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%ST_Coeff(i,j,k)) ! Surface Tension in interfacial cells
            else
                call get_surface_tension_coefficient(this,this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%ST_Coeff(i,j,k))
            endif
        end do
        end do
    end do

end subroutine update_surface_tension_coefficients

subroutine get_column_surface_tension_and_height(this,i0,j0,k0,dir,ret)
    use messager,  only: die
    use vfs_class, only: vfs
    implicit none

    class(conservative_st_type), intent(inout) :: this
    integer :: i0,j0,k0 !i0,j0,k0 tell us the original cell location
    integer :: dir  ! 1 = x oriented, 2 = y oriented, 3 = z oriented. Sign gives direction

    integer :: iplus,jplus,kplus ! Tells us how to add
    integer :: i,j,k,count ! Current Cell Location
    real(WP),dimension(5) :: ret ! [i,j,k,H,ST] where i,j,k are for final cell visited in bottom direction and H is the height. ST is the column wise surface tension coefficient.
    real(WP) :: H,ct,cb,ST,volSum
    LOGICAL :: bound,whileTest

    ! Initial Values
    i = i0; j = j0; k = k0
    ct = this%vf%VF(i0,j0,k0)
    H = this%vf%VF(i0,j0,k0)
    if((ct .gt. VFlo) .and. (ct .le. VFhi)) then 
        ST = this%ST_Coeff(i0,j0,k0) * ct
        volSum = ct
    else
        ST = 0
        volSum = 0
    endif
    count = 0

    ! Convert dir to Index Additions
    SELECT CASE (abs(dir))
        CASE (1) ! x oriented
            iplus = sign(1,dir)
            jplus = 0
            kplus = 0
        CASE (2) ! y oriented
            iplus = 0
            jplus = sign(1,dir)
            kplus = 0
        CASE (3) ! z oriented
            iplus = 0
            jplus = 0
            kplus = sign(1,dir)
    END SELECT

    ! Set Initial Bound
    if(ct .lt. VFhi) then 
        bound = .TRUE.
    else
        bound = .FALSE.
    endif

    ! Find Top
    do while(((.not. bound) .or. ( (ct.lt.VFhi) .and. (ct.gt. VFlo) )) ) ! If not at interface, or at interface exactly, then go to top neighbor, optional: .and. count .lt. 3
        ! Go To Top Neighbor
        i = i + iplus
        j = j + jplus
        k = k + kplus
        ! Update ct
        ct = this%vf%VF(i,j,k)
        ! Update H
        H = H + ct
        if((ct .gt. VFlo) .and. (ct .le. VFhi)) then
            bound = .true.
            volSum = volSum + ct
            ST = ST + ct * this%ST_Coeff(i,j,k)
        endif
        count = count + 1
    enddo

    ! Catch Inconsistent result
    if((ct .gt. VFlo)) then
        H = -1.0_WP
        ret = [-1.0_WP,-1.0_WP,-1.0_WP,-1.0_WP,-1.0_WP]
        RETURN
    endif

    ! Reset to center
    i = i0
    j = j0
    k = k0
    cb = this%vf%VF(i0,j0,k0)
    count = 0
    ! Set Initial Bound
    if(cb .gt. VFlo) then 
        bound = .TRUE.
    else
        bound = .FALSE.
    endif

    do while(((.not. bound) .or. ( (cb.lt.VFhi) .and. (cb.gt. VFlo) )) ) ! If not at interface, or at interface exactly, then go to Bottom neighbo, optional : .and. count .lt. 3
        ! Go To Bottom Neighbor
        i = i - iplus
        j = j - jplus
        k = k - kplus
        ! Update cb
        cb = this%vf%VF(i,j,k)
        ! Update H
        H = H + cb
        if((cb .gt. VFlo) .and. (cb .le. VFhi)) then
            bound = .true.
            volSum = volSum + cb
            ST = ST + cb * this%ST_Coeff(i,j,k)
        endif
        count = count + 1
    enddo

    ! Catch Inconsistent result
    if(cb .lt. VFhi) then
        ! print *, "Inconsistent 2"
        H = -2.0_WP
        ret = [-2.0_WP,-2.0_WP,-2.0_WP,-2.0_WP,-2.0_WP]
        RETURN
    endif

    !return N and H
    ret = [0.0_WP,0.0_WP,0.0_WP,0.0_WP,0.0_WP]
    ret(1) = i
    ret(2) = j
    ret(3) = k
    ret(4) = H
    ret(5) = ST/volSum
end subroutine get_column_surface_tension_and_height

subroutine update_columns(this) ! A visualization tool 
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), dimension(5) :: result
    real(WP), dimension(3) :: planeNormal
    integer :: i,j,k,dir

    this%STx = 0.0_WP
    this%STy = 0.0_WP
    this%Hx = 0.0_WP 
    this%Hy = 0.0_WP

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            planeNormal = calculateNormal(this%vf%interface_polygon(1,i,j,k))
            ! If the outward normal is in the positive X direction, the top dir is 1. Else -1. This is for when columns are expanding in the X direction (Hx)
            if(planeNormal(1) .gt. 0) then
                dir = 1
            elseif(planeNormal(1) .lt. 0) then
                dir = -1
            else
                dir = int(sign(1.0_WP,this%fs%cfg%xm(i)))
            endif
            call get_column_surface_tension_and_height(this,i,j,k,dir,result)
            ! Store
            this%STx(i,j,k) = result(5)
            this%Hx(i,j,k) = result(4)
            ! If the outward normal is in the positive Y direction, the top dir is 2. Else -2. This is for when columns are expanding in the Y direction (HY)
            if(planeNormal(2) .gt. 0) then
                dir = 2
            elseif(planeNormal(2) .lt. 0) then
                dir = -2
            else
                dir = int(sign(2.0_WP, this%fs%cfg%ym(j)))
            endif
            call get_column_surface_tension_and_height(this,i,j,k,dir,result)
            this%STy(i,j,k) = result(5)
            this%Hy(i,j,k) = result(4)
            ! Put Z direction here, if going into 3D.
        end do
        end do
    end do
end subroutine update_columns

subroutine height_function_curvature(this,i,j,k,dir,curv,dSTds,curv_type)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,dir,n,iplus,jplus,kplus
    real(WP), dimension(1) :: curv
    real(WP) :: h0,hm1,hp1,dx,hP,hPP,shift ! h0 is the center height, hm1 is the cell to the "left" and hp1 is the cell to the right.
    real(WP) :: st0,stm1,stp1,ds
    integer, dimension(3) :: cell0,cellm1,cellp1
    real(WP), dimension(5) :: cellAndHeight
    real(WP),intent(out),optional :: dSTds
    integer,intent(out),optional :: curv_type
    ! If Full/Empty, return 0
    if((this%vf%VF(i,j,k) .gt. VFhi) .or. (this%vf%VF(i,j,k) .lt. VFlo)) then
        curv = 0.0_WP
        dSTds = 0.0_WP
        RETURN
    endif
    ! Convert dir to index additions in perpendicular direction
    SELECT CASE (abs(dir))
        CASE (1) ! x oriented
        iplus = 0
        jplus = -(sign(int(1),dir))
        kplus = -(sign(int(1),dir))
        CASE (2) ! y oriented
        iplus = (sign(int(1),dir))
        jplus = 0
        kplus = (sign(int(1),dir))
        CASE (3) ! z oriented
        iplus = (sign(int(1),dir))
        jplus = (sign(int(1),dir))
        kplus = 0
    END SELECT

    ! Calculate H0,N0
    ! write (*,'(A)') 'COMPUTE H0 ================================================='
    call this%get_column_surface_tension_and_height(i,j,k,dir,cellAndHeight)
    cell0(1) = cellAndHeight(1)
    cell0(2) = cellAndHeight(2)
    cell0(3) = cellAndHeight(3)
    h0 = cellAndHeight(4)
    st0 = cellAndHeight(5)
    ! Right now, focus on 2D, so there are only 2
    ! Calculate hp1, Np1
    ! write (*,'(A)') 'COMPUTE Hp1 ================================================='
    ! write (*,'(A,3I6.2)') 'PLUS i,j,dir = ',iplus,jplus,dir
    call this%get_column_surface_tension_and_height(i+iplus,j+jplus,k,dir,cellAndHeight)
    cellp1(1) = cellAndHeight(1)
    cellp1(2) = cellAndHeight(2)
    cellp1(3) = cellAndHeight(3)
    hp1 = cellAndHeight(4)
    stp1 = cellAndHeight(5)
    ! Calculate hm1, Nm1
    ! write (*,'(A)') 'COMPUTE Hm1 ================================================='
    call this%get_column_surface_tension_and_height(i-iplus,j-jplus,k,dir,cellAndHeight)
    cellm1(1) = cellAndHeight(1)
    cellm1(2) = cellAndHeight(2)
    cellm1(3) = cellAndHeight(3)
    hm1 = cellAndHeight(4)
    stm1 = cellAndHeight(5)
    ! Consistency Check
    if((sign(1.0_WP,hp1) .ne. -1.0_WP) .and. (sign(1.0_WP,hm1) .ne. -1.0_WP)) then 
        ! consistent
        dx = (this%vf%cfg%xm(cell0(1)) - this%vf%cfg%xm(cell0(1)+iplus)) + (this%vf%cfg%ym(cell0(2)) - this%vf%cfg%ym(cell0(2)+jplus))
        dx = -dx
        
        h0 = h0*abs(dx)
        ! Move Heights to a common origin
        hm1 = hm1*abs(dx)
        shift = this%vf%cfg%xm(cellm1(1)) * (2-abs(dir))+ this%vf%cfg%ym(cellm1(2)) * (abs(dir)-1)- this%vf%cfg%xm(cell0(1)) * (2-abs(dir)) - this%vf%cfg%ym(cell0(2)) * (abs(dir)-1)
        hm1 = hm1 + shift*sign(1.0_WP,real(dir,WP))
    
        hp1 = hp1*abs(dx)
        shift = this%vf%cfg%xm(cellp1(1)) * (2-abs(dir))+ this%vf%cfg%ym(cellp1(2)) * (abs(dir)-1)- this%vf%cfg%xm(cell0(1)) * (2-abs(dir)) - this%vf%cfg%ym(cell0(2)) * (abs(dir)-1)
        hp1 = hp1 + shift*sign(1.0_WP,real(dir,WP)) 
        
        ! Now that everything is in the commmon frame, finite differences
        hP = (hp1-hm1)/(2*abs(dx))
        hPP = (hp1 - 2* h0 + hm1)/(dx*dx)

        ! Calculate Curvature
        curv = -hPP/((1.0_WP+hP*hP)**(1.5_WP))
        ! I have included a factor of -1 here. This is because all the curvatures are exactly negative what they should be.
        ! The current theory is that I should not multiply by abs(dx) for h0,hm1,hp1 but simply dx. This could fix that problem.
        ! However, doing this would require me to then reevaluate how the shift is done, and i think that would hurt me right now.
        ! As such, I will simply not be doing that right, and instead adding in the -1. I think in an ideal world I would fix this, 
        ! but we do not live in an ideal world. 

        ! In the case we want to surface gradient of surface tension
        ds = 2 * abs(dx) * sqrt(1+hP*hP)
        dSTds = -(stp1-stm1)/ds 
        ! I have learned here that the orientation of the surface, which in this implementation is CW, influences these derivatives.
        ! These negative signs in front of the curvature and the surface gradient are a result of this CW orientation and 
        curv_type = 1
      else 
        curv = this%vf%curv(i,j,k)
        curv_type = 2
      endif
      
      this%Stp1(i,j,k) = stp1 
      this%St0(i,j,k) = st0
      this%Stm1(i,j,k) = stm1

      this%Hp1(i,j,k) =  hp1 
      this%H0(i,j,k) =  h0
      this%Hm1(i,j,k) = hm1

end subroutine height_function_curvature


subroutine young_normal_method(this,i,j,k,YoungNormal)
    use messager,  only: die
    use vfs_class, only: vfs

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,n
    ! Young Method Variables
    real(WP) :: Cbar_i,Cbar_ip1
    real(WP), dimension(1:3,1:4) :: normals
    real(WP), dimension(1:3) :: YoungNormal

    ! Top Right Corner
    Cbar_i = (this%vf%VF(i,j,k) + this%vf%VF(i,j+1,k))/2
    Cbar_ip1 = (this%vf%VF(i+1,j,k) + this%vf%VF(i+1,j+1,k))/2 
    normals(1,1) = this%vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

    Cbar_i = (this%vf%VF(i,j,k) + this%vf%VF(i+1,j,k))/2
    Cbar_ip1 = (this%vf%VF(i,j+1,k) + this%vf%VF(i+1,j+1,k))/2 
    normals(2,1) = this%vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
    
    normals(3,1) = 0.0_WP

    ! Bottom Right Corner
    Cbar_i = (this%vf%VF(i,j,k) + this%vf%VF(i,j-1,k))/2
    Cbar_ip1 = (this%vf%VF(i+1,j,k) + this%vf%VF(i+1,j-1,k))/2 
    normals(1,2) = this%vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

    Cbar_i = (this%vf%VF(i,j-1,k) + this%vf%VF(i+1,j-1,k))/2 
    Cbar_ip1 = (this%vf%VF(i,j,k) + this%vf%VF(i+1,j,k))/2
    normals(2,2) = this%vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
    
    normals(3,2) = 0.0_WP

    ! Top Left Corner
    Cbar_i = (this%vf%VF(i-1,j,k) + this%vf%VF(i-1,j+1,k))/2
    Cbar_ip1 = (this%vf%VF(i,j,k) + this%vf%VF(i,j+1,k))/2 
    normals(1,3) = this%vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

    Cbar_i = (this%vf%VF(i,j,k) + this%vf%VF(i-1,j,k))/2 
    Cbar_ip1 = (this%vf%VF(i,j+1,k) + this%vf%VF(i-1,j+1,k))/2
    normals(2,3) = this%vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
    
    normals(3,3) = 0.0_WP

    ! Bottom Left Corner
    Cbar_i = (this%vf%VF(i-1,j,k) + this%vf%VF(i-1,j-1,k))/2
    Cbar_ip1 = (this%vf%VF(i,j,k) + this%vf%VF(i,j-1,k))/2 
    normals(1,4) = this%vf%cfg%dxmi(i)*(Cbar_i-Cbar_ip1)

    Cbar_i = (this%vf%VF(i,j-1,k) + this%vf%VF(i-1,j-1,k))/2 
    Cbar_ip1 = (this%vf%VF(i,j,k) + this%vf%VF(i-1,j,k))/2
    normals(2,4) = this%vf%cfg%dymi(j)*(Cbar_i-Cbar_ip1)
    
    normals(3,4) = 0.0_WP

    ! Combine into a single normal
    YoungNormal(:) = (normals(:,1)+normals(:,2)+normals(:,3)+normals(:,4))/4

    ! Normalize
    YoungNormal = YoungNormal/(sqrt(YoungNormal(1)*YoungNormal(1) + YoungNormal(2)*YoungNormal(2)))
   end subroutine young_normal_method

   subroutine centered_columns_normal_method(this,i,j,k,normal)
      use messager,  only: die
      use vfs_class, only: vfs

      implicit none
      class(conservative_st_type), intent(inout) :: this
      integer :: i,j,k,n
      ! Young Method Variables
      real(WP) :: Xmx,Xmy,Ymx,Ymy,norm
      real(WP), dimension(1:3,1:2) :: normals
      real(WP), dimension(1:3) :: normal

      ! X Major Orientation
      ! First get Sign of mx
      Xmx = sign(1.0_WP,(this%vf%VF(i+1,j,k)-this%vf%VF(i-1,j,k)))
      ! Now get value of my
      Xmy = 0.5 * (sum(this%vf%VF(i-1:i+1,j+1,k))-sum(this%vf%VF(i-1:i+1,j-1,k)))
      ! Normalize
      norm = abs(Xmx)+abs(Xmy)
      Xmx = Xmx/norm
      Xmy = Xmy/norm 
      ! Store
      normals(1,1) = Xmx
      normals(2,1) = Xmy
      normals(3,1) = 0.0_WP

      ! Y Major Orientation
      ! First get Sign of my
      Ymy = sign(1.0_WP,(this%vf%VF(i,j+1,k)-this%vf%VF(i,j-1,k)))
      ! Now get value of mx
      Ymx = 0.5 * (sum(this%vf%VF(i+1,j-1:j+1,k))-sum(this%vf%VF(i-1,j-1:j+1,k)))

      ! Normalize
      norm = abs(Ymx)+abs(Ymy)
      Ymx = Ymx/norm
      Ymy = Ymy/norm 
      ! Store
      normals(1,2) = Ymx
      normals(2,2) = Ymy
      normals(3,2) = 1.0_WP

      ! Pick the Better Normal
      if(abs(Ymy) .gt. abs(Xmx)) then
         normal = normals(:,2)
      else
         normal = normals(:,1)
      endif
      ! Normalize
      normal = -normal
      ! normal = -normal/(sqrt(normal(1)*normal(1) + normal(2)*normal(2)))
      ! ! STOP
   end subroutine centered_columns_normal_method
      
   subroutine mixed_young_center_normal(this,i,j,k,normal_return)
      use messager,  only: die
      use vfs_class, only: vfs

      implicit none
      class(conservative_st_type), intent(inout) :: this
      integer :: i,j,k,n
      ! Variables
      real(WP), dimension(1:3) :: YoungNormal,CCNormal
      real(WP) :: mYoung,mCC
      ! Return Value
      real(WP), dimension(1:3) :: normal_return

      ! ================================== Young's Method
      call young_normal_method(this,i,j,k,YoungNormal)
      ! if(abs(YoungNormal(1)) .gt. abs(YoungNormal(2))) then ! X Oriented
      !    ! Normalize to have Unit x
      !    YoungNormal = YoungNormal/abs(YoungNormal(1))
      !    mYoung = YoungNormal(2)
      ! else
      !    YoungNormal = YoungNormal/abs(YoungNormal(2))
      !    mYoung = YoungNormal(1)
      ! endif
      ! ================================== CC Method
      call centered_columns_normal_method(this,i,j,k,CCNormal)
      ! if(abs(CCNormal(1)) .gt. abs(CCNormal(2))) then ! X Oriented
      !    ! Normalize to have Unit x
      !    CCNormal = CCNormal/abs(CCNormal(1))
      !    mCC = CCNormal(2)
      ! else
      !    CCNormal = CCNormal/abs(CCNormal(2))
      !    mCC = CCNormal(1)
      ! endif

      ! ================================== Pick the one with smaller angular offset from aligned
      normal_return = YoungNormal ! Default
      if(CCNormal(3) .eq. 0.0_WP) then
         if(abs(CCNormal(1)) .le. abs(YoungNormal(1))) then
            normal_return = CCNormal
         endif
      else 
         if(abs(CCNormal(2)).le. abs(YoungNormal(1))) then
            normal_return = CCNormal
            normal_return(3) = 0.0_WP
         endif
      endif
      if(sqrt(normal_return(1)*normal_return(1)+normal_return(2)*normal_return(2)) .gt. 0.0_WP) then 
        normal_return = normal_return/sqrt(normal_return(1)*normal_return(1)+normal_return(2)*normal_return(2))
      else 
        normal_return =(/0.0_WP,0.0_WP,0.0_WP/)
      endif
   end subroutine mixed_young_center_normal

   ! Abu Al Saud Methods
subroutine populate_levelset_with_pu(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in ! Current Cell Location 
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver

    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,startPoint,endPoint,Gcen,Lcen,force,pos
    integer, dimension(1:3) :: ind,indguess
    real(WP) :: dx,dy,C,xEval,yEval,zEval

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    ! Create Neighborhood and solver
    call new(neighborhood) 
    call new(solver)
    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
        ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then

            ! Empty Neighborhood
            call emptyNeighborhood(neighborhood) 
            ! Add 7x7 (3 on each side) stencil, in plane
            ! write(*,'(A)') '=============== New: '
            do j_in = -3,3
                do i_in = -3,3
                if(this%vf%VF(i+i_in,j+j_in,k) .gt. this%vf%VFmin .and. this%vf%VF(i+i_in,j+j_in,k) .lt. this%vf%VFmax) then
                    
                    ! Get Plane
                    plane = this%vf%liquid_gas_interface(i+i_in,j+j_in,k)
                    cen = calculateCentroid(this%vf%interface_polygon(1,i+i_in,j+j_in,k))
                    call addMember(neighborhood,cen,1.0_WP,plane,0.0_WP)
                endif
                end do
            end do
            ! Set Neighborhood in solver
            call setNeighborhood(solver,neighborhood)
            ! Get values and populate phi
            xEval = this%fs%cfg%xm(i)
            yEval = this%fs%cfg%ym(j)
            zEval = this%fs%cfg%zm(k)
            call getValue(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),this%phi(i,j,k))   
        end do
        end do
    end do
    ! Boundary Conditions  
    call this%vf%cfg%sync(this%phi)
end subroutine populate_levelset_with_pu

subroutine update_levelset_derivatives(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver

    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,startPoint,endPoint,Gcen,Lcen,force,pos
    integer, dimension(1:3) :: ind,indguess
    real(WP) :: dx,dy,numer,denom

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)

    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%phi_x(i,j,k) = (this%phi(i+1,j,k) - this%phi(i-1,j,k))/(2*dx)
            this%phi_y(i,j,k) = (this%phi(i,j+1,k) - this%phi(i,j-1,k))/(2*dy)
            this%phi_xx(i,j,k) = (this%phi(i+1,j,k) - 2*this%phi(i,j,k) + this%phi(i-1,j,k))/(dx*dx)
            this%phi_yy(i,j,k) = (this%phi(i,j+1,k) - 2*this%phi(i,j,k) + this%phi(i,j-1,k))/(dy)
            this%phi_xy(i,j,k) = (this%phi(i+1,j+1,k) - this%phi(i-1,j+1,k) - this%phi(i+1,j-1,k)+this%phi(i-1,j-1,k))/(4*dx*dy)
            
            numer = this%phi_x(i,j,k)*this%phi_x(i,j,k)*this%phi_yy(i,j,k)
            numer = numer - 2*this%phi_x(i,j,k)*this%phi_y(i,j,k)*this%phi_xy(i,j,k)
            numer = numer + this%phi_y(i,j,k)*this%phi_y(i,j,k)*this%phi_xx(i,j,k)

            denom = this%phi_x(i,j,k)*this%phi_x(i,j,k)+this%phi_y(i,j,k)*this%phi_y(i,j,k)
            denom = denom ** 1.5

            this%curvature_abu(i,j,k) = numer/denom
        end do
        end do
    end do
    ! Boundary Conditions  
    call this%vf%cfg%sync(this%phi_x)
    call this%vf%cfg%sync(this%phi_x)
    call this%vf%cfg%sync(this%phi_xx)
    call this%vf%cfg%sync(this%phi_yy)
    call this%vf%cfg%sync(this%phi_xy)
end subroutine update_levelset_derivatives

subroutine populate_diagonal(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver

    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,startPoint,endPoint,Gcen,Lcen,force,pos
    integer, dimension(1:3) :: ind,indguess
    real(WP) :: dx,dy,numer,denom,xi,kappa,gamma
    real(WP) :: tx,ty

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)

    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%sigma_xx(i,j,k) = 0.0_WP
            this%sigma_yy(i,j,k) = 0.0_WP
            do i_in = -1,1
                ! Sigma_xx
                if(this%phi(i,j,k)*(this%phi(i,j,k)+this%phi(i,j+i_in,k)) .le. 0.0_WP) then 
                    xi = this%phi(i,j,k)/(this%phi(i,j,k)-this%phi(i,j+i_in,k))
                    tx = (this%phi(i,j+1,k) - this%phi(i,j-1,k))/2 + i_in * xi *(this%phi(i,j-1,k)-2*this%phi(i,j,k)+this%phi(i,j+1,k))
                    tx = tx/dx
                    kappa =this%curvature_abu(i,j,k) + xi*(this%curvature_abu(i,j+i_in,k)-this%curvature_abu(i,j,k))
                    gamma = this%fs%sigma
                    this%sigma_xx(i,j,k) = this%sigma_xx(i,j,k) + gamma * (abs(tx)/dx - sign(1.0_WP,this%phi(i,j,k))*kappa*(0.5-xi))
                endif
                ! Sigma_yy
                if(this%phi(i,j,k)*(this%phi(i,j,k)+this%phi(i+i_in,j,k)) .le. 0.0_WP) then 
                    xi = this%phi(i,j,k)/(this%phi(i,j,k)-this%phi(i+i_in,j,k))
                    tx = (this%phi(i+1,j,k) - this%phi(i-1,j,k))/2 + i_in * xi *(this%phi(i-1,j,k)-2*this%phi(i,j,k)+this%phi(i+1,j,k))
                    tx = tx/dx
                    kappa =this%curvature_abu(i,j,k) + xi*(this%curvature_abu(i+i_in,j,k)-this%curvature_abu(i,j,k))
                    gamma = this%fs%sigma
                    this%sigma_yy(i,j,k) = this%sigma_yy(i,j,k) + gamma * (abs(tx)/dx - sign(1.0_WP,this%phi(i,j,k))*kappa*(0.5-xi))
                endif
            enddo
        end do
        end do
    end do
    ! Boundary Conditions  
    call this%vf%cfg%sync(this%sigma_xx)
    call this%vf%cfg%sync(this%sigma_yy)
end subroutine populate_diagonal

subroutine populate_off_diagonal(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver

    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,startPoint,endPoint,Gcen,Lcen,force,pos
    integer, dimension(1:3) :: ind,indguess
    real(WP) :: dx,dy,numer,denom,xi,kappa,gamma
    real(WP) :: tx,ty

    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)

    ! Loop over real domain 
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%sigma_xy(i,j,k) = 0.0_WP
            this%sigma_yx(i,j,k) = 0.0_WP
            ! Sigma_xy
            if((this%phi(i-1,j,k) + this%phi(i-1,j-1,k))*(this%phi(i,j,k)+this%phi(i,j-1,k)) .gt. 0.0_WP) then 
                this%sigma_xy(i,j,k) = 0.0_WP
            else
                numer = this%phi(i-1,j,k)+this%phi(i-1,j-1,k)
                denom = this%phi(i-1,j,k)+this%phi(i-1,j-1,k)-this%phi(i,j,k)-this%phi(i,j-1,k)
                xi = numer/denom
                tx = this%phi(i-1,j,k)-this%phi(i-1,j-1,k)
                tx = tx + xi * (this%phi(i,j,k)-this%phi(i-1,j,k) + this%phi(i-1,j-1,k)-this%phi(i,j-1,k))
                tx = tx/dx 
                gamma = this%fs%sigma

                this%sigma_xy(i,j,k) = -gamma*sign(1.0_WP,this%phi(i,j,k)+this%phi(i,j-1,k))*tx/dx
            endif
            ! Sigma_yx
            if((this%phi(i,j-1,k) + this%phi(i-1,j-1,k))*(this%phi(i,j,k)+this%phi(i-1,j,k)) .gt. 0.0_WP) then 
                this%sigma_yx(i,j,k) = 0.0_WP
            else
                numer = this%phi(i,j-1,k)+this%phi(i-1,j-1,k)
                denom = this%phi(i,j-1,k)+this%phi(i-1,j-1,k)-this%phi(i,j,k)-this%phi(i-1,j,k)
                xi = numer/denom
                tx = this%phi(i,j-1,k)-this%phi(i-1,j-1,k)
                tx = tx + xi * (this%phi(i,j,k)-this%phi(i,j-1,k) + this%phi(i-1,j-1,k)-this%phi(i-1,j,k))
                tx = tx/dx 
                gamma = this%fs%sigma

                this%sigma_yx(i,j,k) = -gamma*sign(1.0_WP,this%phi(i,j,k)+this%phi(i-1,j,k))*tx/dx
            endif
        end do
        end do
    end do
    ! Boundary Conditions  
    call this%vf%cfg%sync(this%phi_x)
    call this%vf%cfg%sync(this%phi_x)
    call this%vf%cfg%sync(this%phi_xx)
    call this%vf%cfg%sync(this%phi_yy)
    call this%vf%cfg%sync(this%phi_xy)
end subroutine populate_off_diagonal

subroutine add_Saud_surface_tension_jump(this,dt,div,contact_model) ! Need to Update
    use messager,  only: die
    use vfs_class, only: vfs
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class 
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), intent(inout) :: dt     !< Timestep size over which to advance 
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer, intent(in), optional :: contact_model
    integer :: i,j,k
    real(WP) :: mycurv,mysurf
    
    ! Store old jump
    this%fs%DPjx=this%fs%Pjx
    this%fs%DPjy=this%fs%Pjy
    this%fs%DPjz=this%fs%Pjz
    ! Calculate pressure jump
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
            ! X face
            mysurf=sum(this%vf%SD(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i-1:i,j,k)*this%vf%curv(i-1:i,j,k)*this%fs%cfg%vol(i-1:i,j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjx(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divu_x(:,i,j,k)*this%vf%VF(i-1:i,j,k))
            ! Y face
            mysurf=sum(this%vf%SD(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j-1:j,k)*this%vf%curv(i,j-1:j,k)*this%fs%cfg%vol(i,j-1:j,k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjy(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divv_y(:,i,j,k)*this%vf%VF(i,j-1:j,k))
            ! Z face
            mysurf=sum(this%vf%SD(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))
            if (mysurf.gt.0.0_WP) then
                mycurv=sum(this%vf%SD(i,j,k-1:k)*this%vf%curv(i,j,k-1:k)*this%fs%cfg%vol(i,j,k-1:k))/mysurf
            else
                mycurv=0.0_WP
            end if
            this%fs%Pjz(i,j,k)=this%fs%sigma*mycurv*sum(this%fs%divw_z(:,i,j,k)*this%vf%VF(i,j,k-1:k))
        end do
        end do
    end do

    this%PjxD = this%fs%Pjx
    this%PjyD = this%fs%Pjy 
    this%PjzD = this%fs%Pjz

    ! Update Values
    call this%populate_levelset_with_pu()
    call this%update_levelset_derivatives()
    call this%populate_diagonal()
    call this%populate_off_diagonal
    call this%updateSurfaceTensionForces()
    call this%updateSurfaceTensionForces3D()

    this%fs%Pjx = this%Pjx_ST
    this%fs%Pjy = this%Pjy_ST
    this%fs%Pjz = this%Pjz_ST
    
    SELECT CASE (this%SmoothingOption)
        CASE(1)
        
        CASE(2)
        call this%applyLaplacianSmoothing()
        CASE(3)
        call this%applyGradientSmoothing()
        CASE(4)
        call this%applyPoissonSmoothing()
        case(5) ! This is using the Peskin Delta Function
        call this%applyPeskinSmoothing()
        CASE DEFAULT

    END SELECT

    ! Add wall contact force to pressure jump 
    if (present(contact_model)) then
        select case (contact_model)
        case (1)
        call this%fs%add_static_contact(vf=this%vf)
        case default
        call die('[tpns: add_surface_tension_jump] Unknown contact model!')
        end select
    end if
    
    ! Compute jump of DP
    this%fs%DPjx=this%fs%Pjx-this%fs%DPjx
    this%fs%DPjy=this%fs%Pjy-this%fs%DPjy
    this%fs%DPjz=this%fs%Pjz-this%fs%DPjz
    
    ! Add div(Pjump) to RP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
        do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            this%SurfaceTensionDiv(i,j,k) = dt*(sum(this%fs%divp_x(:,i,j,k)*this%fs%DPjx(i:i+1,j,k)/this%fs%rho_U(i:i+1,j,k))&
            &                                  +sum(this%fs%divp_y(:,i,j,k)*this%fs%DPjy(i,j:j+1,k)/this%fs%rho_V(i,j:j+1,k)))
            ! &                           +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
            ! if(abs(SurfaceTensionDiv(i,j,k)) .gt. 1e-12) then
            !    write(*,'(A,3F20.10)') "Surface Tension Div, rho_U,div: ", SurfaceTensionDiv(i,j,k),this%rho_U(i,j,k),div(i,j,k)
            ! endif
            this%fs%div(i,j,k)=this%fs%div(i,j,k)+this%SurfaceTensionDiv(i,j,k)
        end do
        end do
    end do
end subroutine add_Saud_surface_tension_jump

! Visualization Tools

subroutine updatePartitionOfUnity(this)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in,k_in,count ! Current Cell Location 
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver
    
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,Gcen,Lcen,force,pos,tangentHolder
    integer, dimension(1:3) :: ind,indguess
    real(WP) :: dx,dy,xEval,yEval,zEval,Scale
    ! write(*,'(A,3F35.10)') ' ============================================================================== Scale value: ', Scale 
    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dx(1)
    count = 0
    ! Create Neighborhood and solver
    call new(neighborhood) 
    call new(solver)
    ! Loop over real domain 
    do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
        do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
        do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            ! write(*,'(A,3I10.5)') ' ==== location: ', i,j,k 
            ! if(vf%VF(i,j,k) .gt. vf%VFmin .and. vf%VF(i,j,k) .lt. vf%VFmax) then  
                ! Get Points
                xEval = this%vf%cfg%xm(i)
                yEval = this%vf%cfg%ym(j)
                zEval = this%vf%cfg%zm(k)
                pos = (/xEval,yEval,zEval/)
                ! write(*,'(A,3F10.5)') ' ==== pos: ', pos
                ! Get Global Location
                indGuess = (/int(i),int(j),int(k)/)
                ind = this%vf%cfg%get_ijk_global(pos,indGuess)
                ! write(*,'(A,3I10.5)') ' ==== ind value: ', ind
                
                ! Now we have the global index in ind, so we can move around based on that.  
                ! Empty Neighborhood
                call emptyNeighborhood(neighborhood)
                ! Add 7x7 (3 on each side) stencil, in plane
                ! cen = (/vf%cfg%xm(ind(1)),vf%cfg%ym(ind(2)),vf%cfg%zm(ind(3))/)
                ! call addMember(neighborhood,cen,vf%liquid_gas_interface(ind(1),ind(2),ind(3)))
                do j_in = -3,3
                    do i_in = -3,3
                        do k_in = -3,3
                            call this%addCellToNeighborhood(neighborhood,i+i_in,j+j_in,k+k_in)
                        enddo
                    end do
                end do
                ! Set Neighborhood in solver 
                call setNeighborhood(solver,neighborhood)
                ! if(ind(1) .eq. 11 .and. ind(2) .eq. 22 .and. count .lt. 1) then 
                !    write(*,'(A,F10.8)') '> Time  =  ', time%t
                !    ! call printSolver(solver)
                !    count = count + 1
                ! endif
                ! ====== Get Value
                ! print *, "Get Value"
                call getValue(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),this%PartitionOfUnityValue(i,j,k))   
                ! write(*,'(A,F10.8)') '> Radius  =  ', radius
                ! call getValue(solver,xEval,yEval,zEval,0.25_WP,(/0.0_WP,0.0_WP,0.0_WP/),PartitionOfUnityValue(i,j,k))
                ! print *, "Get Normal"
                call getNormalPU(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),tangentHolder) 
                ! call getTangent(solver,xEval,yEval,zEval,radius,center,tangentHolder)
                this%PUTangent_x(i,j,k) = tangentHolder(1)
                this%PUTangent_y(i,j,k) = tangentHolder(2)
                this%PUTangent_z(i,j,k) = tangentHolder(3)
                this%PUTangent_mag(i,j,k) = sqrt(sum(tangentHolder**2))
                ! print *, "Get Weight"
                call getWeight(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),this%PartitionOfUnityWeight(i,j,k)) 
                this%distanceField(i,j,k) = this%PartitionOfUnityValue(i,j,k)/(this%PUTangent_mag(i,j,k)+1e-12)
                ! print *, "Get End"
                ! call getWeight(solver,xEval,yEval,zEval,radius,center,PartitionOfUnityWeight(i,j,k))
                ! write(*,'(A,F35.10)') ' ================= Partition of unity value: ', PartitionOfUnityValue(i,j,k)   
            ! endif 
        end do
        end do
    end do
end subroutine updatePartitionOfUnity

subroutine getProjectedGeometry(this,initialIndex,projectedPoint,normal,curvature,value,weight)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in,k_in,count ! Current Cell Location 
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver
    real(WP), dimension(1:3) :: initialPoint,projectedPoint,normal
    integer, dimension(1:3) :: initialIndex
    real(WP) :: curvature
    real(WP),optional :: value, weight
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,Gcen,Lcen,force,pos,tangentHolder,cenInitial
    real(WP) :: dx,dy,xEval,yEval,zEval,Scale
    ! write(*,'(A,3F35.10)') ' ============================================================================== Scale value: ', Scale 
    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    ! Create Neighborhood and solver
    call new(neighborhood) 
    call new(solver)
    call emptyNeighborhood(neighborhood)
    ! Add 7x7 (3 on each side) stencil, in plane 
    ! cen = (/vf%cfg%xm(ind(1)),vf%cfg%ym(ind(2)),vf%cfg%zm(ind(3))/)
    ! call addMember(neighborhood,cen,vf%liquid_gas_interface(ind(1),ind(2),ind(3)))
    i = initialIndex(1)
    j = initialIndex(2)
    k = initialIndex(3) 
    ! print *, " Making Neighborhood"
    do j_in = -3,3
        do i_in = -3,3
            do k_in = -3,3
                call this%addCellToNeighborhood(neighborhood,i+i_in,j+j_in,k+k_in)
            enddo
        end do
    end do
    ! Set Neighborhood in solver  
    ! print *,"Setting Neigh"
    call setNeighborhood(solver,neighborhood)
    call setKernelSize(solver,this%PU_spread*this%vf%cfg%dx(1))
    ! print *," Set Neighborhood" 
    ! Project to PU
    ! print *, "Centroid Getting"
    cenInitial = calculateCentroid(this%vf%interface_polygon(1,i,j,k))
    ! print *, "Centroid Got" 
    call projectToPU(solver,cenInitial,projectedPoint)
    ! print *,"Projected" 
    xEval = projectedPoint(1)
    yEval = projectedPoint(2)
    zEval = projectedPoint(3)

    ! if(abs(sqrt(sum((projectedPoint)**2))-0.5)/0.5 .gt. 1) then 
    !     print *, "Initial Point:", cenInitial
    !     print *, "Final Point: ", projectedPoint
    !     print *, "Distance: ", sqrt(sum((cenInitial-projectedPoint)**2))
    ! endif 
    ! ===== Get Value
    ! print *,"Value"
    if(present(value)) then 
        call getValue(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),value)   
    endif
    ! ===== Tangent Value 
    ! print *,"Normal"
    call getNormalPU(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),normal) 
    ! ===== Curvature
    ! print *,"Curvature"
    call getMeanCurvaturePU(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),curvature)
    ! ===== Optional Weight
    if(present(weight)) then 
        call getWeight(solver,xEval,yEval,zEval,this%PU_spread*this%vf%cfg%dx(1),weight) 
    endif
    
end subroutine getProjectedGeometry

subroutine getProjectedGeometryEllipsoid(this,initialIndex,A,v,projectedPoint,normal,curvature,evaluationPoint)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in,k_in,count ! Current Cell Location 
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver
    real(WP), dimension(1:3) :: initialPoint,projectedPoint,tangent
    real(WP), dimension(1:3), optional :: evaluationPoint,normal
    real(WP), optional :: curvature
    real(WP), dimension(3) :: v
    real(WP), dimension(3,3) :: A
    integer, dimension(1:3) :: initialIndex
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,Gcen,Lcen,force,pos,tangentHolder,cenInitial
    real(WP) :: dx,dy,xEval,yEval,zEval,Scale
    ! write(*,'(A,3F35.10)') ' ============================================================================== Scale value: ', Scale 
    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    ! Create Neighborhood and solver 
    call new(neighborhood) 
    call new(solver)
    call emptyNeighborhood(neighborhood)
    ! Add 7x7 (3 on each side) stencil, in plane 
    ! cen = (/vf%cfg%xm(ind(1)),vf%cfg%ym(ind(2)),vf%cfg%zm(ind(3))/)
    ! call addMember(neighborhood,cen,vf%liquid_gas_interface(ind(1),ind(2),ind(3)))
    i = initialIndex(1)
    j = initialIndex(2)
    k = initialIndex(3) 
    
    ! Set Neighborhood in solver 
    call setNeighborhood(solver,neighborhood)

    ! Project to PU
    cenInitial = calculateCentroid(this%vf%interface_polygon(1,i,j,k))
    call projectToEllipsoid(solver,cenInitial,A(:,1),A(:,2),A(:,3),v,projectedPoint)
    
    
    ! If the evaluation Point is present, use that to get normal and curvature of ellipsoide.
    ! Else use the Projected Point
    if(present(evaluationPoint)) then 
        xEval = evaluationPoint(1)
        yEval = evaluationPoint(2)
        zEval = evaluationPoint(3)
    else 
        xEval = projectedPoint(1)
        yEval = projectedPoint(2)
        zEval = projectedPoint(3)
    endif 
    ! Only evaluate if these are present 
    if(present(normal)) then 
        call getNormalEllipsoid(solver,xEval,yEval,zEval,A(:,1),A(:,2),A(:,3),v,normal)
    endif

    if(present(curvature)) then 
        call getMeanCurvatureEllipsoid(solver,xEval,yEval,zEval,A(:,1),A(:,2),A(:,3),v,curvature)
    endif


end subroutine getProjectedGeometryEllipsoid

subroutine getValuePU(this,pt,value)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in,k_in ! Current Cell Location 
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type) :: solver
    integer, dimension(1:3) :: initialIndex,indGuess
    real(WP) :: value
    ! Temp Items
    type(SeparatorVariant_type) :: plane
    real(WP), dimension(1:3) :: cen,pt
    real(WP) :: dx,dy
    ! write(*,'(A,3F35.10)') ' ============================================================================== Scale value: ', Scale 
    dx = this%vf%cfg%dx(1)
    dy = this%vf%cfg%dy(1)
    ! Create Neighborhood and solver
    call new(neighborhood) 
    call new(solver)
    call emptyNeighborhood(neighborhood)
    ! Add 7x7 (3 on each side) stencil, in plane 
    indGuess = (/1,1,1/)
    initialIndex = this%vf%cfg%get_ijk_global(pt,indGuess)
    ! print *,"Got Index"  
    i = initialIndex(1)
    j = initialIndex(2)
    k = initialIndex(3) 
    ! print *,"Assigned IJK",i,j,k
    do j_in = -3,3
        do i_in = -3,3
            do k_in = -3,3
                if(this%vf%VF(i+i_in,j+j_in,k+k_in) .gt. 1e-12 .and. this%vf%VF(i+i_in,j+j_in,k+k_in) .lt. 1.0_WP - 1e-12) then
                    if(i+i_in .gt. -2 .and. i+i_in .lt. this%vf%cfg%imax+3 .and. & 
                       j+j_in .gt. -2 .and. j+j_in .lt. this%vf%cfg%jmax+3 .and. & 
                       k+k_in .gt. -2 .and. k+k_in .lt. this%vf%cfg%kmax+3) then 
                        cen = calculateCentroid(this%vf%interface_polygon(1,i+i_in,j+j_in,k+k_in))
                        ! Get Plane
                        plane = this%vf%liquid_gas_interface(i+i_in,j+j_in,k+k_in)
                        call addMember(neighborhood,cen,1.0_WP,plane,0.0_WP)
                    endif
                endif
            enddo
        end do
    end do
    ! print *,"Set Neighborhood"  
    ! Set Neighborhood in solver  
    call setNeighborhood(solver,neighborhood)
    ! print *,"Done Solver"


    ! ===== Get Value
    call getValue(solver,pt(1),pt(2),pt(3),this%PU_spread*this%vf%cfg%dx(1),value) 
    ! print *,"Got Value"  
end subroutine getValuePU

subroutine updateStresses(this, A,v)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class
    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), dimension(3) :: v
    real(WP), dimension(3,3) :: A

    call this%updateSurfaceTensionStresses3D()
    
    call this%updateSurfaceTensionStresses3DEllipsoid(A,v)

end subroutine updateStresses

subroutine updateForces(this, A,v)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class
    implicit none
    class(conservative_st_type), intent(inout) :: this
    real(WP), dimension(3) :: v
    real(WP), dimension(3,3) :: A

    call this%updateSurfaceTensionForces3D()

    call this%updateSurfaceTensionForces3DEllipsoid()

end subroutine updateForces

subroutine addCellToNeighborhood(this,neighborhood,i,j,k,solver)
    use irl_fortran_interface
    use f_PUNeigh_RectCub_class
    use f_SeparatorVariant_class
    use f_PUST_RectCub_class, only: solveFace

    implicit none
    class(conservative_st_type), intent(inout) :: this
    integer :: i,j,k,j_in,i_in,k_in,shift_first_index,shift_second_index ! Current Cell Location
    type(PUNeigh_RectCub_type) :: neighborhood
    type(PUST_RectCub_type), optional :: solver
    ! Temp Items
    integer :: jmax
    type(SeparatorVariant_type) :: plane,facePlane
    real(WP), dimension(1:4) :: planeVector
    real(WP) :: signedDistance0,signedDistance1,signedDistance2,signedDistance3
    logical :: hasPositive, hasNegative
    real(WP), dimension(1:3) :: cen,corner,dvec,shift, cenInitial,projectedPoint,projectedNormal,planeNormal
    real(WP), dimension(1:3) :: P0,P1,P2,P3,pressure_cell_center,velocity_cell_center,face_center
    real(WP) :: dx,dy,dz
    real(WP) :: x_i,y_j,z_k
    real(WP) :: weight,weight_normal
    real(WP) :: alignment
    ! if(this%TwoD) then 
    !     jmax = 2 
    ! else 
    !     jmax = 3
    ! endif
    if(present(solver)) then 
        cenInitial = (/this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)/)
        ! print *, "Centroid Got" 
        ! call projectToPU(solver,cenInitial,projectedPoint)
        call getNormalPU(solver,cenInitial(1),cenInitial(2),cenInitial(3),this%PU_spread*this%vf%cfg%dx(1),projectedNormal)
        ! if(sqrt(sum(projectedNormal**2)) .gt. 1e-12) then  
        !     call projectToPU(solver,cenInitial,projectedPoint)
        !     call getNormalPU(solver,projectedPoint(1),projectedPoint(2),projectedPoint(3),this%PU_spread*this%vf%cfg%dx(1),projectedNormal)
        ! else  
        !     projectedNormal = (/0.0_WP,0.0_WP,0.0_WP/) 
        ! endif
    endif



    if(this%vf%VF(i,j,k) .gt. 1e-12 .and. this%vf%VF(i,j,k) .lt. 1.0_WP - 1e-12) then ! Mixed Cell
        ! First add the plane itself
        cen = calculateCentroid(this%vf%interface_polygon(1,i,j,k))
        plane = this%vf%liquid_gas_interface(i,j,k)
        ! Compute the weight associated with this
        weight = 1.0_WP 
        if(this%vf%VF(i,j,k) .lt. 0.1) then 
            weight = 0.5_WP - 0.5_WP * COS(10.0_WP*Pi*this%vf%VF(i,j,k))
        endif 

        if(this%vf%VF(i,j,k) .gt. 0.9) then 
            weight = 0.5_WP - 0.5_WP * COS(10.0_WP*Pi*(1.0_WP-this%vf%VF(i,j,k)))
        endif 

        ! Compute Weight Associated with Normal, if solver is present
        weight_normal = 1.0_WP
        if(present(solver)) then 
            planeNormal = calculateNormal(this%vf%interface_polygon(1,i,j,k))
            alignment = sum(projectedNormal*planeNormal)
            weight_normal = max(alignment,0.0_WP)
        endif

        call addMember(neighborhood,cen,weight*weight_normal,plane,0.0_WP)

        ! ! Now, find the faces that are hit by the plane.
        ! pressure_cell_center = (/this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)/)
        ! dx = this%vf%cfg%dx(1)
        ! dy = this%vf%cfg%dy(1)
        ! dz = this%vf%cfg%dz(1)
        ! dvec = (/dx,dy,dz/)
        ! shift = dvec
        ! do i_in = -1,1,2 ! Cell Loop, This tells us the force direction
        !     do j_in = 1,jmax ! Face Loop, This tells us the normal direction
        !         ! Shift from Pressure cell center ot the face center we are looking at:
        !         shift = (/0.0_WP,0.0_WP,0.0_WP/)  
        !         shift(j_in) = dvec(j_in)/2
        !         face_center = pressure_cell_center + real(i_in, WP)*shift
        !         ! Now that we have the face center and normal direction (all positive normals)
        !         ! we can use the pattern below to get a CCW orientation:
        !         shift_first_index = mod(j_in,3)+1
        !         shift_second_index = mod(j_in+1,3)+1
        !         ! Apply to get first corner
        !         shift = dvec/2     
        !         shift(j_in) = 0.0_WP
        !         shift(shift_first_index) = -shift(shift_first_index)
        !         shift(shift_second_index) = -shift(shift_second_index)
        !         P0 = face_center + shift

        !         ! Apply to get second corner
        !         shift = dvec/2     
        !         shift(j_in) = 0.0_WP
        !         shift(shift_first_index) = shift(shift_first_index)
        !         shift(shift_second_index) = -shift(shift_second_index)
        !         P1 = face_center + shift

        !         ! Apply to get third corner
        !         shift = dvec/2     
        !         shift(j_in) = 0.0_WP
        !         shift(shift_first_index) = shift(shift_first_index)
        !         shift(shift_second_index) = shift(shift_second_index)
        !         P2 = face_center + shift

        !         ! Apply to get fourth corner
        !         shift = dvec/2     
        !         shift(j_in) = 0.0_WP
        !         shift(shift_first_index) = -shift(shift_first_index)
        !         shift(shift_second_index) = shift(shift_second_index)
        !         P3 = face_center + shift

        !         ! Now that we have the face corners, evaluate to see if the plane is hit (different signs in at least one corner) 
        !         planeVector = getPlane(plane,0)
        !         ! Evaluate Signed Distance at each point
        !         signedDistance0 = sum(planeVector(1:3)*P0) - planeVector(4)
        !         signedDistance1 = sum(planeVector(1:3)*P1) - planeVector(4)
        !         signedDistance2 = sum(planeVector(1:3)*P2) - planeVector(4)
        !         signedDistance3 = sum(planeVector(1:3)*P3) - planeVector(4)
        !         ! Check if there is one positive and one negative signed distance
        !         hasPositive = (signedDistance0 .ge. 0.0_WP) .or. (signedDistance1 .ge. 0.0_WP) .or. &
        !                       (signedDistance2 .ge. 0.0_WP) .or. (signedDistance3 .ge. 0.0_WP)
        !         hasNegative = (signedDistance0 .lt. 0.0_WP) .or. (signedDistance1 .lt. 0.0_WP) .or. &
        !                       (signedDistance2 .lt. 0.0_WP) .or. (signedDistance3 .lt. 0.0_WP)
        !         ! If there is at least one positive and one negative, add that face
        !         if(hasPositive .and. hasNegative) then 
        !             cen = 0.25_WP*(P0+P1+P2+P3) ! Centroid at face center 
        !             call new(facePlane)
        !             planeVector(1:3) = 0.0_WP
        !             planeVector(j_in) = real(i_in, WP)
        !             planeVector(4) = dot_product(planeVector(1:3),cen) 

        !             call setNumberOfPlanes(facePlane,1)
        !             call setPlane(facePlane,0,planeVector(1:3),planeVector(4))
        !             call addMember(neighborhood,cen,1.0_WP,facePlane,0.0_WP)
        !         endif
        !     enddo
        ! enddo
    endif
end subroutine addCellToNeighborhood

subroutine fill_momentum_vfs(this) 
    integer :: i,j,k 
    real(WP) :: liq_vol, gas_vol,tot_vol
    class(conservative_st_type), intent(inout) :: this

    ! Initialize VFs
    ! Calculate rho_U/V/Wold using subcell phasic volumes
    this%vf_x%VF = 0.0_WP
    this%vf_y%VF = 0.0_WP
    this%vf_z%VF = 0.0_WP
    do k=this%vf%cfg%kmino_+1,this%vf%cfg%kmaxo_
        do j=this%vf%cfg%jmino_+1,this%vf%cfg%jmaxo_
            do i=this%vf%cfg%imino_+1,this%vf%cfg%imaxo_
                ! U-cell
                liq_vol=sum(this%vf%Lvol(0,:,:,i,j,k))+sum(this%vf%Lvol(1,:,:,i-1,j,k))
                gas_vol=sum(this%vf%Gvol(0,:,:,i,j,k))+sum(this%vf%Gvol(1,:,:,i-1,j,k))
                tot_vol=gas_vol+liq_vol
                this%vf_x%vf(i,j,k) = 0.0_WP
                if (tot_vol.gt.0.0_WP) this%vf_x%vf(i,j,k)=liq_vol/tot_vol
                ! V-cell
                liq_vol=sum(this%vf%Lvol(:,0,:,i,j,k))+sum(this%vf%Lvol(:,1,:,i,j-1,k))
                gas_vol=sum(this%vf%Gvol(:,0,:,i,j,k))+sum(this%vf%Gvol(:,1,:,i,j-1,k))
                tot_vol=gas_vol+liq_vol
                this%vf_y%vf(i,j,k) = 0.0_WP
                if (tot_vol.gt.0.0_WP) this%vf_y%vf(i,j,k)=liq_vol/tot_vol
                ! W-cell
                liq_vol=sum(this%vf%Lvol(:,:,0,i,j,k))+sum(this%vf%Lvol(:,:,1,i,j,k-1))
                gas_vol=sum(this%vf%Gvol(:,:,0,i,j,k))+sum(this%vf%Gvol(:,:,1,i,j,k-1))
                tot_vol=gas_vol+liq_vol
                this%vf_z%vf(i,j,k) = 0.0_WP
                if (tot_vol.gt.0.0_WP) this%vf_z%vf(i,j,k)=liq_vol/tot_vol
            end do
        end do
    end do
    ! Update the band
    call this%vf_x%update_band()
    call this%vf_y%update_band()
    call this%vf_z%update_band()
    ! Perform interface reconstruction from VOF field
    call this%vf_x%build_interface()
    call this%vf_y%build_interface()
    call this%vf_z%build_interface()
    ! Create discontinuous polygon mesh from IRL interface
    call this%vf_x%polygonalize_interface()
    call this%vf_y%polygonalize_interface()
    call this%vf_z%polygonalize_interface()
    ! Calculate distance from polygons
    call this%vf_x%distance_from_polygon()
    call this%vf_y%distance_from_polygon()
    call this%vf_z%distance_from_polygon()
    ! Calculate subcell phasic volumes
    call this%vf_x%subcell_vol()
    call this%vf_y%subcell_vol()
    call this%vf_z%subcell_vol()
    ! Calculate curvature
    call this%vf_x%get_curvature()
    call this%vf_y%get_curvature()
    call this%vf_z%get_curvature()
    ! Perform PPIC reconstruction
    if (this%vf_x%ppic) call this%vf_x%build_quadratic_interface()
    if (this%vf_y%ppic) call this%vf_y%build_quadratic_interface()
    if (this%vf_z%ppic) call this%vf_z%build_quadratic_interface()
    ! Reset moments to guarantee compatibility with interface reconstruction
    call this%vf_x%reset_volume_moments()
    call this%vf_y%reset_volume_moments()
    call this%vf_z%reset_volume_moments()
end subroutine fill_momentum_vfs

subroutine compute_paraboloid_force(this,vf_input,i,j,k,j_in,P0,P1,P2,P3,force)
    use f_PUST_RectCub_class
    class(conservative_st_type), intent(inout) :: this
    type(vfs) :: vf_input
    integer :: i,j,k,j_in
    real(WP), dimension(1:3) :: P0,P1,P2,P3 
    real(WP), dimension(1:3),intent(out) :: force

    integer,dimension(1:3) :: shift_index_normal
    real(WP),dimension(1:3) :: face_normal 
    real(WP),dimension(1:3) :: force1,force2
    real(WP) :: vf
    type(PUST_RectCub_type) :: solver
    
    ! Setup
    shift_index_normal = (/0,0,0/)
    shift_index_normal(j_in) = 1
    face_normal = (/0.0_WP,0.0_WP,0.0_WP/)
    face_normal(j_in) = 1.0_WP

    force1 = (/0.0_WP,0.0_WP,0.0_WP/)
    force2 = (/0.0_WP,0.0_WP,0.0_WP/)
    force  = (/0.0_WP,0.0_WP,0.0_WP/)

    call new(solver)
    ! Current Cell
    vf = vf_input%VF(i,j,k)
    if(vf .gt. 1e-12 .and. vf .lt. 1.0_WP - 1e-12) then 
        call solveFaceParaboloid(solver,this%fs%sigma,P0,P1,P2,P3,vf_input%liquid_gas_interface(i,j,k),this%PressureOption,this%MarangoniOption,force1)
        if(sum(force1*face_normal) .lt. -1e-12) then
            force1 = -force1
        endif
    endif

    ! Adjacent Cell Cell
    vf = vf_input%VF(i+shift_index_normal(1),j+shift_index_normal(2),k+shift_index_normal(3))
    if(vf .gt. 1e-12 .and. vf .lt. 1.0_WP - 1e-12) then 
        call solveFaceParaboloid(solver,this%fs%sigma,P0,P1,P2,P3,&
                                vf_input%liquid_gas_interface(i+shift_index_normal(1),j+shift_index_normal(2),k+shift_index_normal(3)),&
                                this%PressureOption,this%MarangoniOption,force2)
        if(sum(force2*face_normal) .lt. -1e-12) then
            force2 = -force2
        endif
    endif
    ! Combine face forces
    force = (force1 + force2)/2.0_WP
    ! if(sum(force1**2) .gt. 1e-11 .and. sum(force2**2) .gt. 1e-11) then 
    !     force = (force1 + force2)/2.0_WP
    ! elseif(sum(force1**2) .gt. 1e-12 .or. sum(force2**2) .gt. 1e-12) then
    !     force = force1 + force3
    ! else 
    !     force = (force1 + force2)/2.0_WP
    ! endif
    
end subroutine compute_paraboloid_force

subroutine getHalfFace(P0,P1,P2,P3,i_in,j_in,half,R0,R1,R2,R3)
    implicit none

    real(WP), intent(in) :: P0(3),P1(3),P2(3),P3(3)
    integer, intent(in) :: i_in,j_in,half
    real(WP), intent(out) :: R0(3),R1(3),R2(3),R3(3)

    real(WP) :: M01(3),M12(3),M23(3),M30(3)
    integer :: shift_first_index,shift_second_index

    ! Determine the two directions perpendicular to the face normal
    shift_first_index  = mod(j_in,3)+1
    shift_second_index = mod(j_in+1,3)+1

    ! Calculate edge midpoints
    M01 = 0.5_WP*(P0+P1)
    M12 = 0.5_WP*(P1+P2)
    M23 = 0.5_WP*(P2+P3)
    M30 = 0.5_WP*(P3+P0)

    ! Split along momentum direction
    if (i_in == shift_first_index) then

        if (half == 1) then
            ! Negative momentum half
            R0 = P0
            R1 = M01
            R2 = M23
            R3 = P3
        else
            ! Positive momentum half
            R0 = M01
            R1 = P1
            R2 = P2
            R3 = M23
        endif

    elseif (i_in == shift_second_index) then

        if (half == 1) then
            ! Negative momentum half
            R0 = P0
            R1 = P1
            R2 = M12
            R3 = M30
        else
            ! Positive momentum half
            R0 = M30
            R1 = M12
            R2 = P2
            R3 = P3
        endif

    else
        ! Diagonal component: no splitting necessary
        R0 = P0
        R1 = P1
        R2 = P2
        R3 = P3
    endif

end subroutine getHalfFace

end module conservative_st