module temp_transport
    use precision,         only: WP
    use mathtools,         only: Pi
    use geometry,          only: cfg
    use hypre_str_class,   only: hypre_str
    use ddadi_class,       only: ddadi
    use tpns_class,        only: tpns
    use vfs_class,         only: vfs
    use string,            only: str_medium
    use timetracker_class, only: timetracker
    use vtk_class,         only: vtk
    use partmesh_class,    only: partmesh
    use surfmesh_class,    only: surfmesh
    use event_class,       only: event
    use monitor_class,     only: monitor
    use sgrid_class,       only: sgrid
    use config_class,      only: config
    use iterator_class,    only: iterator
    use linsol_class,   only: linsol
    use irl_fortran_interface
    implicit none
    private


    public :: tads
    ! Parameters

    ! Boundary Conditoins for the temperature solver
    type :: bcond
        type(bcond), pointer :: next
        character(len=str_medium) :: name = 'UNNAMED_BCOND'
        integer :: type 
        type(iterator) :: itr
        character(len=1) :: face 
        integer :: dir 
        real(WP) :: rdir
    end type bcond
    ! Type Temperature Advection Diffusion Solver
    type :: tads 
        ! Things we store here
        ! Options
        integer :: interpolationOption = 2 ! 1 for upwind, 2 for central (default)
        ! Upwinding - Currently default values with no way to change them. 
        integer :: nst = 1
        integer :: stp1 = 1
        integer :: stp2 = 1
        integer :: stm1 = 0
        integer :: stm2 = 0
        
        ! Solvers
        type(tpns),pointer,public :: fs => null()
        type(vfs),pointer,public :: vf => null()
        type(timetracker),pointer,public :: time => null()
        ! Work Arrays
        real(WP), dimension(:,:,:,:), allocatable :: hybp_x,hybp_y,hybp_z   !< Hybrid interpolation for P cell 
        real(WP), dimension(:,:,:,:), allocatable :: itp_x,itp_y,itp_z   !< Central Interpolation for P cell
        real(WP), dimension(:,:,:,:), allocatable :: grd_x ,grd_y ,grd_z    !< Gradient Operator
        real(WP), dimension(:,:,:), allocatable :: T                        !< Temperature array
        real(WP), dimension(:,:,:), allocatable :: Told                     !< Old Temperature Array
        real(WP), dimension(:,:,:), allocatable :: H                        !< Temperature array
        real(WP), dimension(:,:,:), allocatable :: Hold                     !< Old Temperature Array
        real(WP), dimension(:,:,:), allocatable :: diff                     ! Holds Conduction values
        real(WP), dimension(:,:,:), allocatable :: cp                   !< Holds value of rho_cp    
        real(WP), dimension(:,:,:), allocatable :: rho                   !< Holds value of rho_cp
        real(WP), dimension(:,:,:), allocatable :: rho_cp                   !< Holds value of rho_cp
        ! Palmore Arrays
        real(WP), dimension(:,:,:), allocatable :: TG,TGold,TGExtrap !< Holds Gas Temperature Field
        real(WP), dimension(:,:,:), allocatable :: TL,TLold,TLExtrap !< Holds Liquid Temperature Field
        real(WP), dimension(:,:,:), allocatable :: TPmix,Tinterface !< Holds Liquid Temperature Field
        real(WP), dimension(:,:,:), allocatable :: uG !< Holds Gas Velocity Field
        real(WP), dimension(:,:,:), allocatable :: uL !< Holds Liquid Velocity Field
        ! Fluid Properties
        real(WP) :: rhoL, rhoG
        real(WP) :: cpL, cpG
        real(WP) :: kL, kG
        ! Implicit Scalar solver
        class(linsol), pointer :: implicit   !< Iterative linear solver object for an implicit prediction of the scalar residual
        ! Initialized Boolean
        logical :: initialized = .false.
    contains
        ! Public Methods
        procedure :: temp
        procedure :: init
        procedure :: add_bcond
        procedure :: apply_bcond
        procedure :: get_dHdt
        procedure :: get_dHdt_SL
        procedure :: solve_implicit
        procedure :: populate_temperature
        procedure :: populate_enthalpy
        ! Private Methods
        procedure :: populate_palmore_arrays
        procedure :: extrapolate_fields_palmore
        procedure :: step_temperature_palmore
        procedure :: mix_temperature_palmore
        procedure :: compute_Aslam_RHS
        procedure :: compute_liquid_face_fraction
        procedure :: compute_interface_temperature
    end type tads
contains
! Method Implementations here

subroutine temp(this)
    class(tads), intent(inout) :: this
    print *, "This is a temporary subroutine in the temp_transport module."
end subroutine temp

subroutine init(this,fs_in,vf_in,time_in)
    class(tads), intent(inout) :: this
    class(tpns),target, intent(in) :: fs_in
    class(vfs),target, intent(in) :: vf_in
    class(timetracker),target, intent(in) :: time_in
    integer :: i,j,k
    ! Assign fs_in to fs
    this%fs => fs_in
    this%vf => vf_in
    this%time => time_in
    
    ! Allocate Arrays
    ! == METRICS ==
    allocate(this%hybp_x( 0:+1,this%fs%cfg%imino_  :this%fs%cfg%imaxo_,this%fs%cfg%jmino_  :this%fs%cfg%jmaxo_,this%fs%cfg%kmino_  :this%fs%cfg%kmaxo_)); this%hybp_x=0.0_WP
    allocate(this%hybp_y( 0:+1,this%fs%cfg%imino_  :this%fs%cfg%imaxo_,this%fs%cfg%jmino_  :this%fs%cfg%jmaxo_,this%fs%cfg%kmino_  :this%fs%cfg%kmaxo_)); this%hybp_y=0.0_WP
    allocate(this%hybp_z( 0:+1,this%fs%cfg%imino_  :this%fs%cfg%imaxo_,this%fs%cfg%jmino_  :this%fs%cfg%jmaxo_,this%fs%cfg%kmino_  :this%fs%cfg%kmaxo_)); this%hybp_z=0.0_WP
    
    ! Allocate finite difference diffusivity interpolation coefficients
    allocate(this%itp_x(-1:0,this%fs%cfg%imin_:this%fs%cfg%imax_+1,this%fs%cfg%jmin_:this%fs%cfg%jmax_+1,this%fs%cfg%kmin_:this%fs%cfg%kmax_+1)) !< X-face-centered
    allocate(this%itp_y(-1:0,this%fs%cfg%imin_:this%fs%cfg%imax_+1,this%fs%cfg%jmin_:this%fs%cfg%jmax_+1,this%fs%cfg%kmin_:this%fs%cfg%kmax_+1)) !< Y-face-centered
    allocate(this%itp_z(-1:0,this%fs%cfg%imin_:this%fs%cfg%imax_+1,this%fs%cfg%jmin_:this%fs%cfg%jmax_+1,this%fs%cfg%kmin_:this%fs%cfg%kmax_+1)) !< Z-face-centered
    ! Create diffusivity interpolation coefficients to cell face
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                this%itp_x(:,i,j,k)=this%fs%cfg%dxmi(i)*[this%fs%cfg%xm(i)-this%fs%cfg%x(i),this%fs%cfg%x(i)-this%fs%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
                this%itp_y(:,i,j,k)=this%fs%cfg%dymi(j)*[this%fs%cfg%ym(j)-this%fs%cfg%y(j),this%fs%cfg%y(j)-this%fs%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
                this%itp_z(:,i,j,k)=this%fs%cfg%dzmi(k)*[this%fs%cfg%zm(k)-this%fs%cfg%z(k),this%fs%cfg%z(k)-this%fs%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
        end do
    end do

    ! Allocate finite difference velocity gradient operators
    allocate(this%grd_x(-1:0,this%fs%cfg%imin_:this%fs%cfg%imax_+1,this%fs%cfg%jmin_:this%fs%cfg%jmax_+1,this%fs%cfg%kmin_:this%fs%cfg%kmax_+1)) !< X-face-centered
    allocate(this%grd_y(-1:0,this%fs%cfg%imin_:this%fs%cfg%imax_+1,this%fs%cfg%jmin_:this%fs%cfg%jmax_+1,this%fs%cfg%kmin_:this%fs%cfg%kmax_+1)) !< Y-face-centered
    allocate(this%grd_z(-1:0,this%fs%cfg%imin_:this%fs%cfg%imax_+1,this%fs%cfg%jmin_:this%fs%cfg%jmax_+1,this%fs%cfg%kmin_:this%fs%cfg%kmax_+1)) !< Z-face-centered
    ! Create gradient coefficients to cell faces
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                this%grd_x(:,i,j,k)=this%fs%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in x from [xm,ym,zm] to [x,ym,zm]
                this%grd_y(:,i,j,k)=this%fs%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in y from [xm,ym,zm] to [xm,y,zm]
                this%grd_z(:,i,j,k)=this%fs%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in z from [xm,ym,zm] to [xm,ym,z]
            end do
        end do
    end do

    ! == Values == 
    allocate(this%T(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Told(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%H(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Hold(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%diff(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%rho_cp(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%cp(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%rho(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(this%TG(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%TGold(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%TGExtrap(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%TL(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%TLold(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%TLExtrap(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%uG(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%uL(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%TPmix(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(this%Tinterface(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    ! Update Init
    this%initialized = .true.
end subroutine init 

subroutine add_bcond(this,name,type,locator,face,dir)
    use string, only: lowercase
    use messager, only: die
    use iterator_class, only: locator_ftype
    implicit none
    class(tads), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    procedure(locator_ftype) :: locator
    character(len=1), intent(in) :: face
    integer, intent(in) :: dir
    type(bcond), pointer :: new_bc
    integer :: i,j,k,n

    ! 
    print *, "ADD_BCOND NOT IMPLEMENTED YET"
end subroutine add_bcond

subroutine apply_bcond(this,t,dt)
    use messager, only: die
    implicit none
    class(tads), intent(inout) :: this
    real(WP), intent(in) :: t,dt
    integer :: i,j,k,n,stag
    type(bcond), pointer :: my_bc
    
    ! Just do brute force on the Hold for the Dirichlet Top and Bottom
    call this%fs%cfg%sync(this%H)
    call this%fs%cfg%sync(this%T)
    
end subroutine apply_bcond

subroutine get_dHdt(this,dHdt ,U,V,W)
    implicit none
    class(tads), intent(inout) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: dHdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
    real(WP) :: H_upx,H_upy,H_upz,diff_x,diff_y,diff_z,U_upx,U_upy,U_upz,Ux,Uy,Uz
    integer :: i,j,k
    ! Updated diff
    ! Testing Quick Scheme
    this%nst=1
    this%stp1=-(this%nst+1)/2; this%stp2=this%nst+this%stp1-1
    this%stm1=-(this%nst-1)/2; this%stm2=this%nst+this%stm1-1
    do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_+1
        do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_+1
            do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_+1
                ! this%diff(i,j,k) = this%vf%VF(i,j,k)*this%kL + (1-this%vf%VF(i,j,k)) *this%kG ! Linear on VF
                this%diff(i,j,k) = 1.0_WP / (this%vf%VF(i,j,k)/this%kL + (1.0_WP-this%vf%VF(i,j,k))/this%kG) ! Harmoinc on VF
                this%cp(i,j,k) = (this%vf%VF(i,j,k)*this%cpL + (1.0_WP-this%vf%VF(i,j,k))*this%cpG)
                this%rho(i,j,k) = (this%vf%VF(i,j,k)*this%fs%rho_l + (1.0_WP-this%vf%VF(i,j,k))*this%fs%rho_g)
                this%rho_cp(i,j,k) = this%cp(i,j,k) * this%rho(i,j,k)                
            end do
        end do
    end do

    call this%fs%cfg%sync(this%diff)
    call this%fs%cfg%sync(this%cp)

    allocate(FX(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FY(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FZ(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! Calculate Upwind Values
                ! Linearly Interpolate U value 
                H_upx=this%H(i,j,k)
                H_upy=this%H(i,j,k)
                H_upz=this%H(i,j,k)
                diff_x=sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k))

                diff_y=sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k))

                diff_z=sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k))
                
                ! Fluxes on x-face
                FX(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*sum(this%H(i+this%stp1:i+this%stp2,j,k)) &
                &         -0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*sum(this%H(i+this%stm1:i+this%stm2,j,k)) &
                &         +diff_x*sum(this%grd_x(:,i,j,k)*(this%H(i-1:i,j,k)/this%cp(i-1:i,j,k)))
                ! Fluxes on y-face
                FY(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*sum(this%H(i,j+this%stp1:j+this%stp2,k)) &
                &         -0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*sum(this%H(i,j+this%stm1:j+this%stm2,k)) &
                &         +diff_y*sum(this%grd_y(:,i,j,k)*(this%H(i,j-1:j,k)/this%cp(i,j-1:j,k)))
                ! Fluxes on z-face
                FZ(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*sum(this%H(i,j,k+this%stp1:k+this%stp2)) &
                &         -0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*sum(this%H(i,j,k+this%stm1:k+this%stm2)) &
                &         +diff_z*sum(this%grd_z(:,i,j,k)*(this%H(i,j,k-1:k)/this%rho_cp(i,j,k-1:k)))

                ! Semi Lagrangian Attempts


            end do
        end do
    end do

    call this%fs%cfg%sync(FX)
    call this%fs%cfg%sync(FY)
    call this%fs%cfg%sync(FZ)

    ! Time derivative of rhoSC
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                dHdt(i,j,k)=sum(this%fs%divp_x(:,i,j,k)*FX(i:i+1,j,k))+&
                &           sum(this%fs%divp_y(:,i,j,k)*FY(i,j:j+1,k))+&
                &           sum(this%fs%divp_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
        end do
    end do
end subroutine get_dHdt

subroutine get_dHdt_SL(this,dHdt ,U,V,W,detailed_face_flux,dt)
    use irl_fortran_interface
    implicit none
    class(tads), intent(inout) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: dHdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    type(TagAccVM_SepVM_type), dimension(1:,this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:) :: detailed_face_flux !< Needs to be (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
    real(WP) :: H_upx,H_upy,H_upz,diff_x,diff_y,diff_z,U_upx,U_upy,U_upz,Ux,Uy,Uz
    real(WP) :: my_vol,my_vol1,my_vol2
    type(SepVM_type) :: my_SepVM
    integer, dimension(3) :: ind
    integer :: i,j,k,n
    real(WP), intent(in) :: dt  !< This is the time step size that was used to generate the detailed_face_flux geometric data
    
    ! print*, "Start"
    ! detailed_face_flux = this%vf%detailed_face_flux
    ! Testing Quick Scheme
    this%nst=1
    this%stp1=-(this%nst+1)/2; this%stp2=this%nst+this%stp1-1
    this%stm1=-(this%nst-1)/2; this%stm2=this%nst+this%stm1-1
    do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_+1
        do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_+1
            do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_+1
                ! this%diff(i,j,k) = this%vf%VF(i,j,k)*this%kL + (1-this%vf%VF(i,j,k)) *this%kG ! Linear on VF
                this%diff(i,j,k) = 1.0_WP / (this%vf%VF(i,j,k)/this%kL + (1.0_WP-this%vf%VF(i,j,k))/this%kG) ! Harmoinc on VF
                this%cp(i,j,k) = (this%vf%VF(i,j,k)*this%cpL + (1.0_WP-this%vf%VF(i,j,k))*this%cpG)
                this%rho(i,j,k) = (this%vf%VF(i,j,k)*this%fs%rho_l + (1.0_WP-this%vf%VF(i,j,k))*this%fs%rho_g)
                this%rho_cp(i,j,k) = this%cp(i,j,k) * this%rho(i,j,k)
            end do
        end do
    end do
    ! print *, "Props"

    ! Calculate minmod-limited gradient of SC everywhere
    ! do k=this%cfg%kmino_+1,this%cfg%kmaxo_-1
    !     do j=this%cfg%jmino_+1,this%cfg%jmaxo_-1
    !         do i=this%cfg%imino_+1,this%cfg%imaxo_-1
    !             ! No need to calculate gradient inside of wall cell
    !             if (this%mask(i,j,k).eq.1) cycle
    !             ! Get gradient
    !             grad(1,i,j,k)=minmod((this%H(i+1,j,k,nsc)-this%H(i,j,k,nsc))*this%cfg%dxmi(i+1),(this%H(i,j,k,nsc)-this%H(i-1,j,k,nsc))*this%cfg%dxmi(i))
    !             grad(2,i,j,k)=minmod((this%H(i,j+1,k,nsc)-this%H(i,j,k,nsc))*this%cfg%dymi(j+1),(this%H(i,j,k,nsc)-this%H(i,j-1,k,nsc))*this%cfg%dymi(j))
    !             grad(3,i,j,k)=minmod((this%H(i,j,k+1,nsc)-this%H(i,j,k,nsc))*this%cfg%dzmi(k+1),(this%H(i,j,k,nsc)-this%H(i,j,k-1,nsc))*this%cfg%dzmi(k))
    !         end do
    !     end do
    ! end do
        
    call this%fs%cfg%sync(this%diff)
    call this%fs%cfg%sync(this%cp)

    allocate(FX(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FY(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FZ(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    ! print*, "Entering Advection Loop"
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! Semi Lagrangian Fluxes
                ! Advective Fluxes\
                ! Flux on x-face
                ! print*, "If X"
                if (getSize(detailed_face_flux(1,i,j,k)).gt.0) then
                    ! print *, "size gotten"
                    ! Detailed geometric flux is available, use geometric fluxing
                    do n=0,getSize(detailed_face_flux(1,i,j,k))-1
                        ! Get cell index for nth object
                        ! print*, "ind X"
                        ind=this%fs%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux(1,i,j,k),n))
                        ! Get SepVM for nth object
                        ! print*, "SepVM X"
                        call getSepVMAtIndex(detailed_face_flux(1,i,j,k),n,my_SepVM)

                        ! Extract volume for first phase (0 is liquid, 1 is gas)\
                        ! print*, "vol1 X"
                        my_vol1=getVolume(my_SepVM,0)
                        my_vol2=getVolume(my_SepVM,1)
                        my_vol = my_vol1+my_vol2
                        ! Increment flux with first order estimate
                        ! print*, "Fx1 X"
                        FX(i,j,k)=FX(i,j,k)-my_vol1*this%fs%rho_l*this%cpL*this%Told(ind(1),ind(2),ind(3))/(my_vol+1e-6)

                        ! Increment flux with first order estimate
                        ! print*, "Fx2 X"
                        FX(i,j,k)=FX(i,j,k)-my_vol2*this%fs%rho_g*this%cpG*this%Told(ind(1),ind(2),ind(3))/(my_vol+1e-6)

                        ! Second order correction
                        !my_bar=getCentroid(my_SepVM,this%phase(nsc))
                        !FX(i,j,k)=FX(i,j,k)-my_vol*(sum(grad(:,ii,jj,kk)*my_bar(:)-my_barold(:)))
                    end do
                    ! Scale by cell face area and time step size
                    FX(i,j,k)=FX(i,j,k)/(dt*this%fs%cfg%dy(j)*this%fs%cfg%dz(k))
                else
                    ! No detailed geometric flux is available, use upwind flux
                    ! print*, "Upwind X" 
                    FX(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*sum(this%H(i+this%stp1:i+this%stp2,j,k)) &
                    &         -0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*sum(this%H(i+this%stm1:i+this%stm2,j,k))

                    ! SCm=0.0_WP; if (VFold(i-1,j,k).ne.real(this%phase(nsc),WP)) SCm=this%SC(i-1,j,k,nsc)+0.5_WP*grad(1,i-1,j,k)*this%cfg%dx(i-1)
                    !  SCp=0.0_WP; if (VFold(i  ,j,k).ne.real(this%phase(nsc),WP)) SCp=this%SC(i  ,j,k,nsc)-0.5_WP*grad(1,i  ,j,k)*this%cfg%dx(i  )
                    !  FX(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*SCm-0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*SCp
                end if
                ! print*, "IF Y"
                ! Flux on y-face
                if (getSize(detailed_face_flux(2,i,j,k)).gt.0) then
                    ! Detailed geometric flux is available, use geometric fluxing
                    do n=0,getSize(detailed_face_flux(2,i,j,k))-1
                        ! Get cell index for nth object
                        ! print*, "ind Y"
                        ind=this%fs%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux(2,i,j,k),n))
                        ! Get SepVM for nth object
                        ! print*, "SepVm Y"
                        call getSepVMAtIndex(detailed_face_flux(2,i,j,k),n,my_SepVM)

                        ! Extract volume for relevant phase
                        ! print*, "Vol1 Y"
                        my_vol1=getVolume(my_SepVM,0)
                        my_vol2=getVolume(my_SepVM,1)
                        my_vol = my_vol1+my_vol2
                        ! Increment flux with first order estimate
                        ! print*, "Fy1 Y"
                        FY(i,j,k)=FY(i,j,k)-my_vol1*this%fs%rho_l*this%cpL*this%Told(ind(1),ind(2),ind(3))/(my_vol+1e-6)
                        ! Increment flux with first order estimate
                        ! print*, "Fy2 Y"
                        FY(i,j,k)=FY(i,j,k)-my_vol2*this%fs%rho_g*this%cpG*this%Told(ind(1),ind(2),ind(3))/(my_vol+1e-6)
                        ! Second order correction
                        !my_bar=getCentroid(my_SepVM,this%phase(nsc))
                        !FY(i,j,k)=FY(i,j,k)-my_vol*(sum(grad(:,ii,jj,kk)*my_bar(:)-my_barold(:)))
                    end do
                    ! Scale by cell face area and time step size
                    FY(i,j,k)=FY(i,j,k)/(dt*this%fs%cfg%dx(i)*this%fs%cfg%dz(k))
                else
                    ! No detailed geometric flux is available, use upwind flux
                    ! print*, "Upwind Y"
                    FY(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*sum(this%H(i,j+this%stp1:j+this%stp2,k)) &
                    &         -0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*sum(this%H(i,j+this%stm1:j+this%stm2,k))
                    
                    ! SCm=0.0_WP; if (VFold(i,j-1,k).ne.real(this%phase(nsc),WP)) SCm=this%SC(i,j-1,k,nsc)+0.5_WP*grad(2,i,j-1,k)*this%cfg%dy(j-1)
                    !  SCp=0.0_WP; if (VFold(i,j  ,k).ne.real(this%phase(nsc),WP)) SCp=this%SC(i,j  ,k,nsc)-0.5_WP*grad(2,i,j  ,k)*this%cfg%dy(j  )
                    !  FY(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*SCm-0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*SCp
                end if

                ! print*, "IF Z"
                ! Flux on z-face
                if (getSize(detailed_face_flux(3,i,j,k)).gt.0) then
                    ! Detailed geometric flux is available, use geometric fluxing
                    do n=0,getSize(detailed_face_flux(3,i,j,k))-1
                        ! Get cell index for nth object
                        ! print*, "ind Z"
                        ind=this%fs%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux(3,i,j,k),n))
                        ! Get SepVM for nth object
                        ! print*, "SepVm Z"
                        call getSepVMAtIndex(detailed_face_flux(3,i,j,k),n,my_SepVM)

                        ! Extract volume for relevant phase
                        ! print*, "Vol1  Z"
                        my_vol1=getVolume(my_SepVM,0)
                        my_vol2=getVolume(my_SepVM,1)
                        my_vol = my_vol1+my_vol2
                        ! Increment flux with first order estimate
                        ! print*, "FZ1"
                        FZ(i,j,k)=FZ(i,j,k)-my_vol1*this%fs%rho_l*this%cpL*this%Told(ind(1),ind(2),ind(3))/(my_vol+1e-6)

                        ! Increment flux with first order estimate
                        ! print*, "FZ2"
                        FZ(i,j,k)=FZ(i,j,k)-my_vol2*this%fs%rho_g*this%cpG*this%Told(ind(1),ind(2),ind(3))/(my_vol+1e-6)
                        
                        ! Second order correction, only one phase
                        !my_bar=getCentroid(my_SepVM,this%phase(nsc))
                        !FZ(i,j,k)=FZ(i,j,k)-my_vol*(sum(grad(:,ii,jj,kk)*my_bar(:)-my_barold(:)))
                    end do
                    ! Scale by cell face area and time step size
                    FZ(i,j,k)=FZ(i,j,k)/(dt*this%fs%cfg%dx(i)*this%fs%cfg%dy(j))
                else
                    ! No detailed geometric flux is available, use upwind flux
                    FZ(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*sum(this%H(i,j,k+this%stp1:k+this%stp2)) &
                    &         -0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*sum(this%H(i,j,k+this%stm1:k+this%stm2))

                    ! SCm=0.0_WP; if (VFold(i,j,k-1).ne.real(this%phase(nsc),WP)) SCm=this%SC(i,j,k-1,nsc)+0.5_WP*grad(3,i,j,k-1)*this%cfg%dz(k-1)
                    !  SCp=0.0_WP; if (VFold(i,j,k  ).ne.real(this%phase(nsc),WP)) SCp=this%SC(i,j,k  ,nsc)-0.5_WP*grad(3,i,j,k  )*this%cfg%dz(k  )
                    !  FZ(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*SCm-0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*SCp
                end if
                ! print*, "endif"
                ! Diffusive Fluxes

                ! FX(i,j,k) = FX(i,j,k) + diff_x*sum(this%grd_x(:,i,j,k)*(this%H(i-1:i,j,k)/this%cp(i-1:i,j,k)))
                ! FY(i,j,k) = FY(i,j,k) + diff_y*sum(this%grd_y(:,i,j,k)*(this%H(i,j-1:j,k)/this%cp(i,j-1:j,k)))
                ! FZ(i,j,k) = FZ(i,j,k) + diff_z*sum(this%grd_z(:,i,j,k)*(this%H(i,j,k-1:k)/this%rho_cp(i,j,k-1:k)))
                
            end do
        end do
    end do
    ! print*, "Exit Advection Loop"
    ! print*, "Sync"
    call this%fs%cfg%sync(FX)
    call this%fs%cfg%sync(FY)
    call this%fs%cfg%sync(FZ)
    ! print*, "Divp"
    ! Time derivative of rhoSC
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                dHdt(i,j,k)=sum(this%fs%divp_x(:,i,j,k)*FX(i:i+1,j,k))+&
                &           sum(this%fs%divp_y(:,i,j,k)*FY(i,j:j+1,k))+&
                &           sum(this%fs%divp_z(:,i,j,k)*FZ(i,j,k:k+1))

                ! dHdt(i,j,k) = -FX(i+1,j,k)+FX(i,j,k)&
                !              &-FY(i,j+1,k)+FY(i,j,k)&
                !              &-FZ(i,j,k+1)+FZ(i,j,k)
            end do
        end do
    end do
    ! print*, "End"
    contains
      
      !> Minmod gradient
      function minmod(g1,g2) result(g)
         implicit none
         real(WP), intent(in) :: g1,g2
         real(WP) :: g
         if (g1*g2.le.0.0_WP) then
            g=0.0_WP
         else
            if (abs(g1).lt.abs(g2)) then
               g=g1
            else
               g=g2
            end if
         end if
      end function minmod

end subroutine get_dHdt_SL



subroutine solve_implicit(this,dt,resH,rhoU,rhoV,rhoW)
    implicit none
    class(tads), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: resH !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)    :: rhoU  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)    :: rhoV  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)    :: rhoW  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    integer :: i,j,k,sti,std
    
    ! If no implicit solver available, just divide by density and return
    if (.not.associated(this%implicit)) then
        resH=resH
        call this%fs%cfg%sync(resH)
        return
    end if

    ! Prepare convective operator
    this%implicit%opr(1,:,:,:)=1.0_WP; this%implicit%opr(2:,:,:,:)=0.0_WP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                ! Loop over divergence stencil
                do std=0,1
                    ! Loop over plus interpolation stencil
                    do sti=this%stp1,this%stp2
                        this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%fs%divp_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)+abs(rhoU(i+std,j,k)))
                        this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%fs%divp_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)+abs(rhoV(i,j+std,k)))
                        this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%fs%divp_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)+abs(rhoW(i,j,k+std)))
                    end do
                    ! Loop over minus interpolation stencil
                    do sti=this%stm1,this%stm2
                        this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%fs%divp_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)-abs(rhoU(i+std,j,k)))
                        this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%fs%divp_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)-abs(rhoV(i,j+std,k)))
                        this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%fs%divp_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)-abs(rhoW(i,j,k+std)))
                    end do
                end do
            end do
        end do
    end do
    
    ! Prepare diffusive operator
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%fs%divp_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grd_x(-1,i+1,j,k)+&
                &                                                                this%fs%divp_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grd_x( 0,i  ,j,k)+&
                &                                                                this%fs%divp_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grd_y(-1,i,j+1,k)+&
                &                                                                this%fs%divp_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grd_y( 0,i,j  ,k)+&
                &                                                                this%fs%divp_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grd_z(-1,i,j,k+1)+&
                &                                                                this%fs%divp_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grd_z( 0,i,j,k  ))
                this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%fs%divp_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grd_x( 0,i+1,j,k))
                this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%fs%divp_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grd_x(-1,i  ,j,k))
                this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%fs%divp_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grd_y( 0,i,j+1,k))
                this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%fs%divp_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grd_y(-1,i,j  ,k))
                this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%fs%divp_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grd_z( 0,i,j,k+1))
                this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%fs%divp_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grd_z(-1,i,j,k  ))
            end do
        end do
    end do
    
    ! Solve the linear system
    call this%implicit%setup()
    this%implicit%rhs=resH
    this%implicit%sol=0.0_WP
    call this%implicit%solve()
    resH=this%implicit%sol

    print *, "solve_implicit NOT IMPLEMENTED YET"
end subroutine solve_implicit

subroutine populate_temperature(this)
    implicit none
    class(tads) :: this
    integer :: i,j,k

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! Update rho*cp
                this%cp(i,j,k) = (this%vf%VF(i,j,k)*this%cpL + (1.0_WP-this%vf%VF(i,j,k))*this%cpG)
                this%rho(i,j,k) = (this%vf%VF(i,j,k)*this%fs%rho_l + (1.0_WP-this%vf%VF(i,j,k))*this%fs%rho_g)
                this%rho_cp(i,j,k) = this%cp(i,j,k) * this%rho(i,j,k)
                ! Get T
                this%T(i,j,k) = this%H(i,j,k) / this%rho_cp(i,j,k)
            end do
        end do
    end do
end subroutine populate_temperature

subroutine populate_enthalpy(this)
    implicit none
    class(tads) :: this
    
    integer :: i,j,k

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! Update rho_cp
                this%cp(i,j,k) = (this%vf%VF(i,j,k)*this%cpL + (1.0_WP-this%vf%VF(i,j,k))*this%cpG)
                this%rho(i,j,k) = (this%vf%VF(i,j,k)*this%fs%rho_l + (1.0_WP-this%vf%VF(i,j,k))*this%fs%rho_g)
                this%rho_cp(i,j,k) = this%cp(i,j,k) * this%rho(i,j,k)
                ! Get T
                this%H(i,j,k) = this%T(i,j,k) * this%rho_cp(i,j,k)
            end do
        end do
    end do
end subroutine populate_enthalpy



! Palmore Functions
subroutine populate_palmore_arrays(this)
    implicit none
    class(tads) :: this
    
    integer :: i,j,k

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! Update rho_cp
                this%cp(i,j,k) = (this%vf%VF(i,j,k)*this%cpL + (1.0_WP-this%vf%VF(i,j,k))*this%cpG)
                this%rho(i,j,k) = (this%vf%VF(i,j,k)*this%fs%rho_l + (1.0_WP-this%vf%VF(i,j,k))*this%fs%rho_g)
                this%rho_cp(i,j,k) = this%cp(i,j,k) * this%rho(i,j,k)

            end do
        end do
    end do
end subroutine populate_palmore_arrays

subroutine step_temperature_palmore(this,dHGdt,dHLdt ,U,V,W,dt)
    implicit none
    class(tads), intent(inout) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: dHGdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: dHLdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), intent(in) :: dt
    real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)  :: VFG,VFL
    real(WP), dimension(:,:,:), allocatable :: FX_G,FY_G,FZ_G,FX_L,FY_L,FZ_L,QI_G,QI_L
    real(WP) :: H_upx,H_upy,H_upz,diff_x,diff_y,diff_z,U_upx,U_upy,U_upz,Ux,Uy,Uz,diff_coeff
    real(WP) :: indicator_x,indicator_y,indicator_z,beta,dist,interface_temp
    real(WP), dimension(3) :: GBC_m,LBC_m,GBC_p,LBC_p,dBC,plicCenter
    integer :: i,j,k

    ! Testing Quick Scheme
    this%nst=1
    this%stp1=-(this%nst+1)/2; this%stp2=this%nst+this%stp1-1 ! stp1 = -1, stp2 = -1
    this%stm1=-(this%nst-1)/2; this%stm2=this%nst+this%stm1-1 ! stm1 = 0, stm2 = 0
    do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_+1
        do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_+1
            do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_+1
                ! this%diff(i,j,k) = this%vf%VF(i,j,k)*this%kL + (1-this%vf%VF(i,j,k)) *this%kG ! Linear on VF
                this%diff(i,j,k) = 1.0_WP / (this%vf%VF(i,j,k)/this%kL + (1.0_WP-this%vf%VF(i,j,k))/this%kG) ! Harmoinc on VF
                this%cp(i,j,k) = (this%vf%VF(i,j,k)*this%cpL + (1.0_WP-this%vf%VF(i,j,k))*this%cpG)
                this%rho(i,j,k) = (this%vf%VF(i,j,k)*this%rhoL + (1.0_WP-this%vf%VF(i,j,k))*this%rhoG)
                this%rho_cp(i,j,k) = this%cp(i,j,k) * this%rho(i,j,k)                
            end do
        end do
    end do

    call this%fs%cfg%sync(this%diff)
    call this%fs%cfg%sync(this%cp)

    VFG = 1.0_WP-this%vf%VFold
    VFL = this%vf%VFold

    allocate(FX_G(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FY_G(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FZ_G(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    allocate(FX_L(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FY_L(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FZ_L(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    ! Advection
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! Calculate Upwind Fluxes
                
                ! Gas Temperature

                ! Fluxes on x-face
                ! FX_G(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*sum(this%TG(i+this%stp1:i+this%stp2,j,k)*this%rhoG*this%cpG*VFG(i+this%stp1:i+this%stp2,j,k)) &
                ! &           -0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*sum(this%TG(i+this%stm1:i+this%stm2,j,k)*this%rhoG*this%cpG*VFG(i+this%stm1:i+this%stm2,j,k)) 
                ! ! Fluxes on y-face
                ! FY_G(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*sum(this%TG(i,j+this%stp1:j+this%stp2,k)*this%rhoG*this%cpG*VFG(i,j+this%stp1:j+this%stp2,k)) &
                ! &           -0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*sum(this%TG(i,j+this%stm1:j+this%stm2,k)*this%rhoG*this%cpG*VFG(i,j+this%stm1:j+this%stm2,k)) 
                ! ! Fluxes on z-face
                ! FZ_G(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*sum(this%TG(i,j,k+this%stp1:k+this%stp2)*this%rhoG*this%cpG*VFG(i,j,k+this%stp1:k+this%stp2)) &
                ! &           -0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*sum(this%TG(i,j,k+this%stm1:k+this%stm2)*this%rhoG*this%cpG*VFG(i,j,k+this%stm1:k+this%stm2)) 

                ! ! Liquid Temperature

                ! ! Fluxes on x-face
                ! FX_L(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*sum(this%TL(i+this%stp1:i+this%stp2,j,k)*this%rhoL*this%cpL*VFL(i+this%stp1:i+this%stp2,j,k)) &
                ! &           -0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*sum(this%TL(i+this%stm1:i+this%stm2,j,k)*this%rhoL*this%cpL*VFL(i+this%stm1:i+this%stm2,j,k)) 
                ! ! Fluxes on y-face
                ! FY_L(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*sum(this%TL(i,j+this%stp1:j+this%stp2,k)*this%rhoL*this%cpL*VFL(i,j+this%stp1:j+this%stp2,k)) &
                ! &           -0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*sum(this%TL(i,j+this%stm1:j+this%stm2,k)*this%rhoL*this%cpL*VFL(i,j+this%stm1:j+this%stm2,k)) 
                ! ! Fluxes on z-face
                ! FZ_L(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*sum(this%TL(i,j,k+this%stp1:k+this%stp2)*this%rhoL*this%cpL*VFL(i,j,k+this%stp1:k+this%stp2)) &
                ! &           -0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*sum(this%TL(i,j,k+this%stm1:k+this%stm2)*this%rhoL*this%cpL*VFL(i,j,k+this%stm1:k+this%stm2)) 

                ! WITHOUT THAT VOLUME FRACTION TERM

                FX_G(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*sum(this%TG(i+this%stp1:i+this%stp2,j,k)*this%rhoG*this%cpG) &
                &           -0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*sum(this%TG(i+this%stm1:i+this%stm2,j,k)*this%rhoG*this%cpG) 
                ! Fluxes on y-face
                FY_G(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*sum(this%TG(i,j+this%stp1:j+this%stp2,k)*this%rhoG*this%cpG) &
                &           -0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*sum(this%TG(i,j+this%stm1:j+this%stm2,k)*this%rhoG*this%cpG) 
                ! Fluxes on z-face
                FZ_G(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*sum(this%TG(i,j,k+this%stp1:k+this%stp2)*this%rhoG*this%cpG) &
                &           -0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*sum(this%TG(i,j,k+this%stm1:k+this%stm2)*this%rhoG*this%cpG) 

                ! Liquid Temperature

                ! Fluxes on x-face
                FX_L(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*sum(this%TL(i+this%stp1:i+this%stp2,j,k)*this%rhoL*this%cpL) &
                &           -0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*sum(this%TL(i+this%stm1:i+this%stm2,j,k)*this%rhoL*this%cpL) 
                ! Fluxes on y-face
                FY_L(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*sum(this%TL(i,j+this%stp1:j+this%stp2,k)*this%rhoL*this%cpL) &
                &           -0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*sum(this%TL(i,j+this%stm1:j+this%stm2,k)*this%rhoL*this%cpL) 
                ! Fluxes on z-face
                FZ_L(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*sum(this%TL(i,j,k+this%stp1:k+this%stp2)*this%rhoL*this%cpL) &
                &           -0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*sum(this%TL(i,j,k+this%stm1:k+this%stm2)*this%rhoL*this%cpL) 

                ! Diffusion Terms added here
                ! Slightly different treatment of mixed and full cells
                indicator_x = this%vf%VF(i,j,k) + this%vf%VF(i-1,j,k)
                indicator_y = this%vf%VF(i,j,k) + this%vf%VF(i,j-1,k)
                indicator_z = this%vf%VF(i,j,k) + this%vf%VF(i,j,k-1)
                ! An indicator will always be between 0 and 2. 0 is full gas, 2 is full liquid, between is mixed. Treat each case differencely. 

                ! Note I may be doing this wrong because right now both liquid and gas components updated flux in cell. 
                ! But maybe I should check to see if the current cell is mixed or full, then based on that only update the flux associated with the phases present. 
                
                ! X Flux
                
                if(indicator_x .gt. 1e-12 .and. indicator_x .lt. 2.0_WP - 1e-12) then ! Mixed X
                    ! Beta Calculation ***** 
                    beta = 1.0_WP
                    ! Get Barycenters
                    GBC_m = this%vf%Gbary(:,i-1,j,k)
                    LBC_m = this%vf%Lbary(:,i-1,j,k)

                    GBC_p = this%vf%Gbary(:,i,j,k)
                    LBC_p = this%vf%Lbary(:,i,j,k)
                    

                    ! Liquid
                    call compute_liquid_face_fraction(this,(/i,j,k/),(/i-1,j,k/),beta)
                    dBC = LBC_p-LBC_m
                    diff_coeff = dBC(1)/sqrt(sum(dBC**2))
                    FX_L(i,j,k) = FX_L(i,j,k) + beta * this%kL * diff_coeff * (this%TL(i,j,k)-this%TL(i-1,j,k))

                    ! Gas 
                    beta = 1.0_WP - beta
                    dBC = GBC_p-GBC_m
                    diff_coeff = dBC(1)/sqrt(sum(dBC**2)) ! dx/distance

                    FX_G(i,j,k) = FX_G(i,j,k) + beta * this%kG * diff_coeff * (this%TG(i,j,k)-this%TG(i-1,j,k))

                    
                    
                else if(indicator_x .lt. 1e-12) then ! Full Gas
                    FX_G(i,j,k) = FX_G(i,j,k) + this%kG * (this%TG(i,j,k) - this%TG(i-1,j,k))/this%fs%cfg%dx(i)
                    ! Do nothing to liquid flux

                else ! Full Liquid
                    FX_L(i,j,k) = FX_L(i,j,k) + this%kL * (this%TL(i,j,k) - this%TL(i-1,j,k))/this%fs%cfg%dx(i)
                    ! Do nothing to gas flux
                endif

                ! y Flux
                if(indicator_y .gt. 1e-12 .and. indicator_y .lt. 2.0_WP - 1e-12) then ! Mixed y
                    ! Get Barycenters
                    GBC_m = this%vf%Gbary(:,i,j-1,k)
                    LBC_m = this%vf%Lbary(:,i,j-1,k)

                    GBC_p = this%vf%Gbary(:,i,j,k)
                    LBC_p = this%vf%Lbary(:,i,j,k)
                    

                    ! Liquid
                    call compute_liquid_face_fraction(this,(/i,j,k/),(/i,j-1,k/),beta)
                    dBC = LBC_p-LBC_m
                    diff_coeff = dBC(2)/sqrt(sum(dBC**2))
                    FY_L(i,j,k) = FY_L(i,j,k) + beta * this%kL * diff_coeff * (this%TL(i,j,k)-this%TL(i,j-1,k))

                    ! Gas 
                    beta = 1.0_WP - beta
                    dBC = GBC_p-GBC_m
                    diff_coeff = dBC(2)/sqrt(sum(dBC**2)) ! dx/distance

                    FY_G(i,j,k) = FY_G(i,j,k) + beta * this%kG * diff_coeff * (this%TG(i,j,k)-this%TG(i,j-1,k))

                    

                else if(indicator_y .lt. 1e-12) then ! Full Gas
                    FY_G(i,j,k) = FY_G(i,j,k) + this%kG * (this%TG(i,j,k) - this%TG(i,j-1,k))/this%fs%cfg%dy(j)
                    ! Do nothing to liquid flux

                else ! Full Liquid
                    FY_L(i,j,k) = FY_L(i,j,k) + this%kL * (this%TL(i,j,k) - this%TL(i,j-1,k))/this%fs%cfg%dy(j)
                    ! Do nothing to gas flux
                endif

                !z Flux
                if(indicator_z .gt. 1e-12 .and. indicator_z .lt. 2.0_WP - 1e-12) then ! Mixed z                    
                    ! Get Barycenters
                    GBC_m = this%vf%Gbary(:,i,j,k-1)
                    LBC_m = this%vf%Lbary(:,i,j,k-1)

                    GBC_p = this%vf%Gbary(:,i,j,k)
                    LBC_p = this%vf%Lbary(:,i,j,k)

                    ! Liquid
                    call compute_liquid_face_fraction(this,(/i,j,k/),(/i,j,k-1/),beta)
                    dBC = LBC_p-LBC_m
                    diff_coeff = dBC(3)/sqrt(sum(dBC**2))
                    FZ_L(i,j,k) = FZ_L(i,j,k) + beta * this%kL * diff_coeff * (this%TL(i,j,k)-this%TL(i,j,k-1))

                    ! Gas 
                    beta = 1.0_WP - beta
                    dBC = GBC_p-GBC_m
                    diff_coeff = dBC(3)/sqrt(sum(dBC**2)) ! dx/distance

                    FZ_G(i,j,k) = FZ_G(i,j,k) + beta * this%kG * diff_coeff * (this%TG(i,j,k)-this%TG(i,j,k-1))

                    

                else if(indicator_z .lt. 1e-12) then ! Full Gas
                    FZ_G(i,j,k) = FZ_G(i,j,k) + this%kG * (this%TG(i,j,k) - this%TG(i,j,k-1))/this%fs%cfg%dx(i)
                    ! Do nothing to liquid flux

                else ! Full Liquid
                    FZ_L(i,j,k) = FZ_L(i,j,k) + this%kL * (this%TL(i,j,k) - this%TL(i,j,k-1))/this%fs%cfg%dx(i)
                    ! Do nothing to gas flux
                endif
                
            end do
        end do
    end do

    call this%fs%cfg%sync(FX_G)
    call this%fs%cfg%sync(FY_G)
    call this%fs%cfg%sync(FZ_G)

    call this%fs%cfg%sync(FX_L)
    call this%fs%cfg%sync(FY_L)
    call this%fs%cfg%sync(FZ_L)

    ! Time derivative of rhoSC
    
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                dHGdt(i,j,k)=sum(this%fs%divp_x(:,i,j,k)*FX_G(i:i+1,j,k))+&
                &            sum(this%fs%divp_y(:,i,j,k)*FY_G(i,j:j+1,k))+&
                &            sum(this%fs%divp_z(:,i,j,k)*FZ_G(i,j,k:k+1))

                dHLdt(i,j,k)=sum(this%fs%divp_x(:,i,j,k)*FX_L(i:i+1,j,k))+&
                &            sum(this%fs%divp_y(:,i,j,k)*FY_L(i,j:j+1,k))+&
                &            sum(this%fs%divp_z(:,i,j,k)*FZ_L(i,j,k:k+1))

                ! Add Interface Fluxes
                this%Tinterface(i,j,k) = 0.0_WP
                if(this%vf%VF(i,j,k) .gt. 1e-12 .and. this%vf%VF(i,j,k) .lt. 1.0_WP - 1e-12) then
                    call this%compute_interface_temperature((/i,j,k/),interface_temp,plicCenter)
                    ! Gas
                    GBC_p = this%vf%Gbary(:,i,j,k)
                    dBC = GBC_p-plicCenter
                    diff_coeff = sqrt(sum(dBC**2))
                    dHLdt(i,j,k) = dHLdt(i,j,k) + (this%KG * (this%TG(i,j,k) - interface_temp)/diff_coeff)*this%vf%SD(i,j,k)

                    ! Liquid
                    LBC_p = this%vf%Lbary(:,i,j,k)
                    dBC = LBC_p-plicCenter
                    diff_coeff = sqrt(sum(dBC**2))
                    dHLdt(i,j,k) = dHLdt(i,j,k) + (this%KG * (interface_temp - this%TL(i,j,k))/diff_coeff)*this%vf%SD(i,j,k)

                    ! Storage
                    this%Tinterface(i,j,k) = interface_temp
                endif
            end do
        end do
    end do

    ! Update Temperature
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                if(1.0_WP-this%vf%VF(i,j,k) .gt. 1e-12_WP) then 
                    ! this%TG(i,j,k) = (dHGdt(i,j,k) * dt + this%rhoG*VFG(i,j,k)*this%cpG*this%TGold(i,j,k))/(this%rhoG*(1.0_WP-this%vf%VF(i,j,k))*this%cpG)
                    this%TG(i,j,k) = (dHGdt(i,j,k) * dt + this%rhoG*this%cpG*this%TGold(i,j,k))/(this%rhoG*this%cpG)
                else
                    this%TG(i,j,k) = (dHGdt(i,j,k) * dt + this%rhoG*this%cpG*this%TGold(i,j,k))/(this%rhoG*this%cpG)
                    this%TG(i,j,k) = 0.0_WP
                endif

                if(this%vf%VF(i,j,k) .gt. 1e-12_WP) then
                    ! this%TL(i,j,k) = (dHLdt(i,j,k) * dt + this%rhoL*VFL(i,j,k)*this%cpL*this%TL(i,j,k))/(this%rhoL*this%vf%VF(i,j,k)*this%cpL)
                    this%TL(i,j,k) = (dHLdt(i,j,k) * dt + this%rhoL*this%cpL*this%TL(i,j,k))/(this%rhoL*this%cpL)
                else
                    this%TL(i,j,k) = (dHLdt(i,j,k) * dt + this%rhoL*this%cpL*this%TL(i,j,k))/(this%rhoL*this%cpL)
                    this%TL(i,j,k) = 0.0_WP
                endif
            end do
        end do
    end do

    call this%fs%cfg%sync(this%TG)
    call this%fs%cfg%sync(this%TL)

end subroutine step_temperature_palmore

subroutine mix_temperature_palmore(this)
    implicit none
    class(tads) :: this

    this%TPmix = this%TL*this%vf%VF + this%TG*(1.0_WP - this%vf%VF)

end subroutine mix_temperature_palmore

subroutine extrapolate_fields_palmore(this,field,on_value,out_field,dt)
    implicit none
    class(tads) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: on_value     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(inout) :: field
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out)  :: out_field   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)  :: dPhidt,k1,k2,k3,k4,temp
    real(WP), intent(in) :: dt
    integer :: i,j,k
    real(WP) :: mag,mag0,time

    ! Out-field Initial condition is just input field
    out_field = field
    temp = field
    call compute_Aslam_RHS(this,temp,on_value,dPhidt)
    mag = sum(dPhidt**2)/real(size(dPhidt),WP)
    mag0 = mag
    time = 0.0_WP
    do while((mag .gt. 1e-6 .and. mag .lt. 10*mag0) .and. time .lt. 2.0_WP)
        
        ! RK4 Time integration
        call compute_Aslam_RHS(this,temp            , on_value,k1)
        call compute_Aslam_RHS(this,temp + k1 * dt/2.0_WP, on_value,k2)
        call compute_Aslam_RHS(this,temp + k2 * dt/2.0_WP, on_value,k3)
        call compute_Aslam_RHS(this,temp + k3 * dt  , on_value,k4)

        ! Update Magnitude
        dPhidt = (k1+2.0_WP*k2+2.0_WP*k3+k4)/6.0_WP
        mag = sum(dPhidt**2)/real(size(dPhidt),WP)
        ! print *, time,mag
        ! Update Outfield
        out_field = temp + dPhidt * dt
        temp = out_field
        time = time + dt
    enddo

end subroutine extrapolate_fields_palmore

subroutine compute_Aslam_RHS(this,field,on_value,dPhidt)
    class(tads) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: field,on_value     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out)  :: dPhidt
    real(WP) :: grad_x, grad_y, grad_z,A,mag
    real(WP) :: normal_x,normal_y,normal_z
    real(WP), dimension(3) :: normal 
    integer :: i,j,k
    real(WP) :: vals(6), eps
    integer  :: cnt

    eps = 1e-14_WP
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                ! Calculate Activation 
                A = 1.0_WP
                if(on_value(i,j,k) .gt. (1.0_WP - 1e-12)) then ! Turn off in the phase where we have data, and just use the data. 
                    A = 0.0_WP
                endif
                ! Get Normal Value
                ! normal_x = (this%vf%VF(i+1,j,k)-this%vf%VF(i-1,j,k))/(2*this%fs%cfg%dx(i))
                ! normal_y = (this%vf%VF(i,j+1,k)-this%vf%VF(i,j-1,k))/(2*this%fs%cfg%dy(j))
                ! normal_z = (this%vf%VF(i,j,k+1)-this%vf%VF(i,j,k-1))/(2*this%fs%cfg%dz(k))

                normal_x = -(on_value(i+1,j,k)-on_value(i-1,j,k))/(2*this%fs%cfg%dx(i))
                normal_y = -(on_value(i,j+1,k)-on_value(i,j-1,k))/(2*this%fs%cfg%dy(j))
                normal_z = -(on_value(i,j,k+1)-on_value(i,j,k-1))/(2*this%fs%cfg%dz(k))

                mag = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z)
                if(mag .gt. 1e-12) then
                    normal_x = normal_x/(mag+1e-12)
                    normal_y = normal_y/(mag+1e-12)
                    normal_z = normal_z/(mag+1e-12)
                endif
                ! Calculate Gradient using Upwinding on normal direction
                if(normal_x .lt. 0.0) then 
                    grad_x = (field(i+1,j,k)-field(i,j,k))/(this%fs%cfg%dx(i))
                else
                    grad_x = (field(i,j,k)-field(i-1,j,k))/(this%fs%cfg%dx(i))
                endif

                if(normal_y .lt. 0.0) then 
                    grad_y = (field(i,j+1,k)-field(i,j,k))/(this%fs%cfg%dy(j))
                else
                    grad_y = (field(i,j,k)-field(i,j-1,k))/(this%fs%cfg%dy(j))
                endif

                if(normal_z .lt. 0.0) then 
                    grad_z = (field(i,j,k+1)-field(i,j,k))/(this%fs%cfg%dz(k))
                else
                    grad_z = (field(i,j,k)-field(i,j,k-1))/(this%fs%cfg%dz(k))
                endif
                

                ! Calculate Value
                dPhidt(i,j,k) = - A * (normal_x*grad_x + normal_y*grad_y + normal_z*grad_z)
            end do
        end do
    end do

end subroutine compute_Aslam_RHS

subroutine compute_liquid_face_fraction(this,index_cell,index_next,fraction)
    class(tads) :: this
    integer, dimension(3), intent(in) :: index_cell,index_next  
    real(WP), intent(out) :: fraction
    real(WP),dimension(1:3,1:8) :: cube_pts_cell,cube_pts_next
    integer :: i,j,k,changing_index
    real(WP), dimension(3) ::  n_cell,n_cell_proj,n_next,n_next_proj,center_cell,center_next,center_face ! Stores Normals
    real(WP), dimension(4) :: plane
    real(WP) :: mixedflag_cell, mixedflag_next ! Tells us if a cell is mixed or not - 0 for empty, -1 for mixed, 1 for full. +-
    real(WP) :: d_cell,d_cell_proj,d_next,d_next_proj ! Distances
    real(WP) :: total_face_area
    real(WP) :: liquid_volume,gas_volume,face_fraction_cell,face_fraction_next
    ! Cube Intersection
    type(RectCub_type) :: cube_cell,cube_next
    type(PlanarSep_type) :: plane_sep_cell,plane_sep_next
    type(SepVM_type) :: phase_moments_cell,phase_moments_next
    ! Since the fractions are not unique, we need to compute them for both and take the smaller one

    ! Establish the first cell 
    i = index_cell(1); j = index_cell(2); k = index_cell(3)
    cube_pts_cell(:,1)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j  ),this%fs%cfg%z(k+1)]; 
    cube_pts_cell(:,2)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j  ),this%fs%cfg%z(k  )];
    cube_pts_cell(:,3)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j+1),this%fs%cfg%z(k  )];
    cube_pts_cell(:,4)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j+1),this%fs%cfg%z(k+1)];
    cube_pts_cell(:,5)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j  ),this%fs%cfg%z(k+1)];
    cube_pts_cell(:,6)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j  ),this%fs%cfg%z(k  )];
    cube_pts_cell(:,7)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j+1),this%fs%cfg%z(k  )];
    cube_pts_cell(:,8)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j+1),this%fs%cfg%z(k+1)];
    call new(cube_cell)
    call construct(cube_cell,cube_pts_cell)

    if(this%vf%VF(i,j,k) .gt. 1e-12 .and. this%vf%VF(i,j,k) .lt. 1.0_WP - 1e-12) then 
        plane = getPlane(this%vf%liquid_gas_interface(i,j,k),0) 
        n_cell=plane(1:3)
        n_cell_proj = n_cell
        d_cell = plane(4)
        ! n_cell=calculateNormal(this%interface_polygon(1,i,j,k))
        mixedflag_cell = -1.0_WP
        center_cell = [this%fs%cfg%xm(i  ),this%fs%cfg%ym(j  ),this%fs%cfg%zm(k  )];
    endif
    ! Second Cell
    i = index_next(1); j = index_next(2); k = index_next(3)
    cube_pts_next(:,1)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j  ),this%fs%cfg%z(k+1)]; 
    cube_pts_next(:,2)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j  ),this%fs%cfg%z(k  )];
    cube_pts_next(:,3)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j+1),this%fs%cfg%z(k  )];
    cube_pts_next(:,4)=[this%fs%cfg%x(i+1),this%fs%cfg%y(j+1),this%fs%cfg%z(k+1)];
    cube_pts_next(:,5)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j  ),this%fs%cfg%z(k+1)];
    cube_pts_next(:,6)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j  ),this%fs%cfg%z(k  )];
    cube_pts_next(:,7)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j+1),this%fs%cfg%z(k  )];
    cube_pts_next(:,8)=[this%fs%cfg%x(i  ),this%fs%cfg%y(j+1),this%fs%cfg%z(k+1)];
    call new(cube_next)
    call construct(cube_next,cube_pts_next)

    if(this%vf%VF(i,j,k) .gt. 1e-12 .and. this%vf%VF(i,j,k) .lt. 1.0_WP - 1e-12) then 
        plane = getPlane(this%vf%liquid_gas_interface(i,j,k),0) 
        n_next=plane(1:3)
        n_next_proj = n_next
        d_next = plane(4)
        mixedflag_next = -1.0_WP
        center_next = [this%fs%cfg%xm(i  ),this%fs%cfg%ym(j  ),this%fs%cfg%zm(k  )];
    endif
    
    ! Get face ceneter
    center_face = (center_cell + center_next) * 0.5_WP 

    ! To get direction of the face, compare coordinates of face center and cell center. The coordinate that changes is the normal direction. 
    
    if(abs(center_face(1) - center_cell(1)) .gt. 1e-12) then 
        changing_index = 1
        total_face_area = this%fs%cfg%dy(j)*this%fs%cfg%dz(k)
    endif

    if(abs(center_face(2) - center_cell(2)) .gt. 1e-12) then 
        changing_index = 2
        total_face_area = this%fs%cfg%dx(i)*this%fs%cfg%dz(k)
    endif

    if(abs(center_face(3) - center_cell(3)) .gt. 1e-12) then 
        changing_index = 3
        total_face_area = this%fs%cfg%dy(j)*this%fs%cfg%dx(i)
    endif

    ! Now that we have the direction, project the normals and calculate new distances, and make new planes and compute new area fractions
    
    if(mixedflag_cell .lt. -0.5_WP) then 
        ! Projection into plane
        n_cell_proj(changing_index) = 0.0_WP
        ! Distance Update
        d_cell_proj = d_cell - n_cell(changing_index) * center_face(changing_index)
        ! Make New Plane
        call new(plane_sep_cell)
        ! Add plane
        call addPlane(plane_sep_cell,n_cell_proj,d_cell_proj)
        ! New Phase Moments
        call new(phase_moments_cell)
        call getNormMoments(cube_cell,plane_sep_cell,phase_moments_cell)
        ! Now we want to get the liquid area
        liquid_volume = getVolume(phase_moments_cell,0)
        gas_volume = getVolume(phase_moments_cell, 1)
        ! Now get face fraction, which is now equal to volume fraction due to projection
        face_fraction_cell = liquid_volume/(liquid_volume+gas_volume) 
    else
        ! If it is full or empty, the above variables are undefine, so we don't want to do this. We want to just ignore it. 
        face_fraction_cell = 1.0_WP ! We can do this since face fraction is [0,1] always. 
    endif
    
    if(mixedflag_next .lt. -0.5_WP) then 
        ! Projection into plane
        n_next_proj(changing_index) = 0.0_WP
        ! Distance Update
        d_next_proj = d_next - n_next(changing_index) * center_face(changing_index)
        ! Make New Plane
        call new(plane_sep_next)
        ! Add plane
        call addPlane(plane_sep_next,n_next_proj,d_next_proj)
        ! New Phase Moments
        call new(phase_moments_next)
        call getNormMoments(cube_next,plane_sep_next,phase_moments_next)
        ! Now we want to get the liquid area
        liquid_volume = getVolume(phase_moments_next,0)
        gas_volume = getVolume(phase_moments_next, 1)
        ! Now get face fraction, which is now equal to volume fraction due to projection
        face_fraction_next = liquid_volume/(liquid_volume+gas_volume) 
    else
        ! If it is full or empty, the above variables are undefine, so we don't want to do this. We want to just ignore it.
        face_fraction_next = 1.0_WP ! We can do this since face fraction is [0,1] always. 
    endif

    ! Return minimum of the two
    fraction = min(face_fraction_cell,face_fraction_next)

end subroutine compute_liquid_face_fraction 

subroutine compute_interface_temperature(this,index,tInterface,xPlic)
    class(tads) :: this
    integer, dimension(3),intent(in) :: index
    real(WP) ,intent(out) :: tInterface
    real(WP) :: tG,tL,dG,dL,numer,denom 
    real(WP), dimension(3), intent(out),optional :: xPlic
    real(WP),dimension(3) :: xL,xG

    tG = this%TG(index(1),index(2),index(3))
    tL = this%TL(index(1),index(2),index(3))

    if(this%vf%VF(index(1),index(2),index(3)) .gt. 1e-12 .and. this%vf%VF(index(1),index(2),index(3)) .lt. 1.0_WP - 1e-12) then 
        ! Mixed

        ! Centers
        xG = this%vf%Gbary(:,index(1),index(2),index(3))
        xL = this%vf%Lbary(:,index(1),index(2),index(3))
        xPlic = calculateCentroid(this%vf%interface_polygon(1,index(1),index(2),index(3)))
        ! Compute Distances
        dG = sqrt(sum((xG-xPlic)**2))
        dL = sqrt(sum((xL-xPlic)**2))

        ! Evaluate
        numer = this%kG *tG/dG + this%kL*tL/dL
        denom = this%kG/dG + this%kL/dL 
        tInterface = numer/denom
    else 
        ! Full
        tInterface = this%vf%VF(index(1),index(2),index(3)) * tL + (1.0_WP-this%vf%VF(index(1),index(2),index(3)))*tG
    endif
end subroutine compute_interface_temperature


end module temp_transport