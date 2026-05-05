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
        real(WP), dimension(:,:,:), allocatable :: rho_cp                   !< Holds value of rho_cp
        ! Fluid Properties
        real(WP) :: rho1, rho2
        real(WP) :: cp1, cp2
        real(WP) :: k1, k2
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
        procedure :: solve_implicit
        procedure :: populate_temperature
        procedure :: populate_enthalpy
        ! Private Methods
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

subroutine get_dHdt(this,dHdt ,rhoU,rhoV,rhoW)
    implicit none
    class(tads), intent(inout) :: this
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(out) :: dHdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: rhoU     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: rhoV     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in)  :: rhoW     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
    real(WP) :: H_upx,H_upy,H_upz,diff_x,diff_y,diff_z
    integer :: i,j,k
    ! Updated diff
    ! Testing Quick Scheme
    this%nst=1
    this%stp1=-(this%nst+1)/2; this%stp2=this%nst+this%stp1-1
    this%stm1=-(this%nst-1)/2; this%stm2=this%nst+this%stm1-1
    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! this%diff(i,j,k) = this%vf%VF(i,j,k)*this%k1 + (1-this%vf%VF(i,j,k)) *this%k2 ! Linear on VF
                this%diff(i,j,k) = 1.0_WP / (this%vf%VF(i,j,k)/this%k1 + (1.0_WP-this%vf%VF(i,j,k))/this%k2) ! Harmoinc on VF
                this%rho_cp(i,j,k) = (this%vf%VF(i,j,k)*this%rho1 + (1.0_WP-this%vf%VF(i,j,k))*this%rho2)* &
                &(this%vf%VF(i,j,k)*this%cp1 + (1.0_WP-this%vf%VF(i,j,k))*this%cp2)
            end do
        end do
    end do


    allocate(FX(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FY(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))
    allocate(FZ(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_))

    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_+1
        do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_+1
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_+1
                ! Calculate Upwind Values
                H_upx=this%H(i,j,k)
                H_upy=this%H(i,j,k)
                H_upz=this%H(i,j,k)
                diff_x=0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%diff(i+this%stp1:i+this%stp2,j,k)) &
                &     +0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%diff(i+this%stm1:i+this%stm2,j,k));

                diff_y=0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%diff(i,j+this%stp1:j+this%stp2,k)) &
                &     +0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%diff(i,j+this%stm1:j+this%stm2,k)) ;

                diff_z=0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%diff(i,j,k+this%stp1:k+this%stp2)) &
                &     +0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%diff(i,j,k+this%stm1:k+this%stm2));

                ! Fluxes on x-face
                FX(i,j,k)=-0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%H(i+this%stp1:i+this%stp2,j,k)) &
                &         -0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%H(i+this%stm1:i+this%stm2,j,k)) &
                &         +diff_x*sum(this%grd_x(:,i,j,k)*this%H(i-1:i,j,k))
                ! Fluxes on y-face
                FY(i,j,k)=-0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%H(i,j+this%stp1:j+this%stp2,k)) &
                &         -0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%H(i,j+this%stm1:j+this%stm2,k)) &
                &         +diff_y*sum(this%grd_y(:,i,j,k)*this%H(i,j-1:j,k))
                ! Fluxes on z-face
                FZ(i,j,k)=-0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%H(i,j,k+this%stp1:k+this%stp2)) &
                &         -0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%H(i,j,k+this%stm1:k+this%stm2)) &
                &         +diff_z*sum(this%grd_z(:,i,j,k)*this%H(i,j,k-1:k))
            end do
        end do
    end do
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
                this%rho_cp(i,j,k) = (this%vf%VF(i,j,k)*this%cp1 + (1.0_WP-this%vf%VF(i,j,k))*this%cp2)
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
                ! Update rho*cp ! (this%vf%VF(i,j,k)*this%rho1 + (1-this%vf%VF(i,j,k))*this%rho2)* &
                this%rho_cp(i,j,k) = (this%vf%VF(i,j,k)*this%cp1 + (1.0_WP-this%vf%VF(i,j,k))*this%cp2)
                ! Get T
                this%H(i,j,k) = this%T(i,j,k) * this%rho_cp(i,j,k)
            end do
        end do
    end do
end subroutine populate_enthalpy


end module temp_transport