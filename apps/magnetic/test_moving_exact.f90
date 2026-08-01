module problem_moving_exact
    use settings
    implicit none

    real(8) :: mu = 4*m_pi
    real(8) :: sigma0_ = 10d0
    
contains

! right hand side function
subroutine f_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)
    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*t**2/mu - sigma0_*t*r**2 + sigma0_*(r**3*(1d0-r)*t**2)
    end select
end subroutine f_func

! velocity function
subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)
    if (.not. allocated(f)) allocate(f(2))


    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = -r*(1d0-r)*x(1)
        f(2) = -r*(1d0-r)*x(2)
    end select
end subroutine u_func
 
! magnetic field on boundary
subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = t**2*r
    end select
end subroutine g_func

! exact solution of B
subroutine B_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),   dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = -t**2*x(2)
        f(2) = t**2*x(1)
    end select
end subroutine B_func

! exact solution of J (J = curl(B/mu))
subroutine J_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),   dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*t**2/mu
    end select
end subroutine J_func

subroutine Lorentz_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),   dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t

    real(8), allocatable :: val_J(:), val_f(:), val_B(:)

    call J_func(x, val_J, DERIV_NONE)
    call f_func(x, val_f, DERIV_NONE)
    call B_func(x, val_B, DERIV_NONE)

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = -(val_J(1) - val_f(1))*val_B(2)
        f(2) = (val_J(1) - val_f(1))*val_B(1)
    end select
end subroutine Lorentz_func


! initial solution of B
subroutine B0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
        f(2) = 0d0
    end select
end subroutine B0_func

! exact solution of E (E+uxB)
subroutine E_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8), dimension(:), allocatable :: val_u, val_B

    real(8) :: t
    common /global/ t

    if (.not. allocated(f)) allocate(f(1))

    call u_func(x, val_u, DERIV_NONE)
    call B_func(x, val_B, DERIV_NONE)

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = t*(x(1)**2+x(2)**2) + val_u(1)*val_B(2) - val_u(2)*val_B(1)
    end select
end subroutine E_func

! initial solution of E
subroutine E0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine E0_func

! function of sigma
subroutine sigma_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sigma0_
    end select
end subroutine sigma_func
    
end module problem_moving_exact

program test_moving_exact
    use settings
    use problem_moving_exact
    use Hydrodynamic
    use magnetic
    use mesh_generator
    use readmeshfile
    use matvec
    implicit none

    character(len=100) :: arg

    ! global variables
    real(8) :: t
    common /global/ t
    real(8) :: dt
    real(8) :: Tend = 1d0

    ! mesh
    integer :: M,Nlist(1)
    real(8),parameter :: Center(2) = [0d0,0d0], Rlist(1) = [1d0]
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! velocity vector
    type(vector)  :: u

    ! Lorentz force
    type(vector)  :: F

    ! error
    real(8) :: L2_err, relL2_err

    ! parameters of Gauss quadrature
    integer, parameter :: Gauss_type = QuadPt16
    integer, parameter :: Gauss_type_bdry = LinePt4
    
    logical :: first_step = .true.

    if(iargc() .ne. 2) then
        write(*,*) "Usage: ./test_moving_exact M dt"
        stop
    end if

    M = 40
    call getarg(1,arg)
    read(arg,*) M
    Nlist(1) = M/2*3
    ! generate mesh
    call QuadMeshOnSector(Center,Rlist,M,Nlist,elems,nodes)

    call HydrodynamicInitialize(elems,nodes)

    call magnetic_init(elems,nodes,B0_func,E0_func,sigma_func,1d0)

    ! time step
    dt = 0.025d0
    call getarg(2,arg)
    read(arg,*) dt
    t = 0d0
    do while(t<Tend-1d-8)
        
        ! Move mesh
        call getVelocity(u_func, u)
        call MoveMesh(u, 0.5d0*dt, nodes)
        call magnetic_update(nodes)
        
        ! Diffusion
        call magnetic_diffusion(t, dt, g_func, f_func, first_step)
        first_step = .false.
        
        ! Move mesh
        t = t+0.5d0*dt
        call getVelocity(u_func, u)
        call MoveMesh(u, 0.5d0*dt, nodes)
        call magnetic_update(nodes)
        
        t = t+0.5d0*dt
    end do

    ! compute error
    call magnetic_errorB(B_func, L2_err, relL2_err)
    write(*,*) "L2 error of B = ", L2_err
    write(*,*) "Relative L2 error of B = ", relL2_err
    t = t-0d0*dt
    call magnetic_errorE(E_func, L2_err)
    t = t+0d0*dt
    write(*,*) "L2 error of E = ", L2_err
    
    call magnetic_lorentz(F)
        
    call magnetic_errorLorentz(F, Lorentz_func, L2_err)
    write(*,*) "L2 error of Lorentz force = ", L2_err

    call magnetic_plot("output/B.vtk")
    

end program test_moving_exact