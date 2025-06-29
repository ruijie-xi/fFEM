module problem_static
    use settings
    implicit none

    real(8) :: mu = 4*m_pi
    
contains

! right hand side function
subroutine f_func(x,f,deriv_type)
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
end subroutine f_func

! velocity function
subroutine u_func(x,f,deriv_type)
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
end subroutine u_func

! current function
function current(t) result(S)
    implicit none
    real(8), intent(in) :: t
    real(8) :: S
    real(8),parameter :: S0 = 1.2d0,b = 0.5d0

    S = S0*(b/t)**(5d0/2d0)*(2d0*b/t-3d0)*exp(-b/t)
end function current
 
! magnetic field on boundary
subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*current(t)/(sqrt(x(1)**2+x(2)**2))
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

    real(8) :: r,D,tmp1,tmp2,BB

    real(8), parameter :: sigma0 = 1d0
    real(8), parameter :: S0 = 1.2d0
    real(8), parameter :: b = 0.5d0
    real(8), parameter :: Rw = 1d0

    r = sqrt(x(1)**2+x(2)**2)
    D = 4*m_pi*sigma0
    tmp1 = b/t
    tmp2 = 1d0 + 0.5d0*sqrt(D/b)*log(Rw/r)

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        if (r<1d-8) then
            f(1) = 0d0
            f(2) = 0d0
        else
            BB = 4d0*S0/r*tmp1**2.5d0*tmp2*(tmp1*tmp2**2-1.5d0)*exp(-tmp1*tmp2**2)
            f(1) = -BB*x(2)/r
            f(2) = BB*x(1)/r
        end if
    end select
end subroutine B_func


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
    real(8) :: sigma0 = 1d0
    real(8) :: r

    if (.not. allocated(f)) allocate(f(1))
    r = sqrt(x(1)**2+x(2)**2)
    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sigma0/r**2
    end select
end subroutine sigma_func
    
end module problem_static


program test_static
    use settings
    use problem_static
    use Hydrodynamic
    use magnetic
    use mesh_generator
    implicit none

    character :: arg*10

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
    type(vector) :: u

    ! Lorentz force
    real(8), dimension(:), allocatable :: F

    ! error
    real(8) :: L2_err, relL2err

    ! parameters of Gauss quadrature
    integer, parameter :: Gauss_type = QuadPt4
    integer, parameter :: Gauss_type_bdry = LinePt3

    if(iargc() .ne. 2) then
        write(*,*) "Usage: ./test_static M dt"
        stop
    end if

    M = 10
    call getarg(1, arg)
    read(arg,*) M
    Nlist(1) = M/2*3
    ! generate mesh
    call QuadMeshOnSector(Center,Rlist,M,Nlist,elems,nodes)
    
    write(*,*) "Number of elements = ", size(elems,2)
    write(*,*) "Number of nodes = ", size(nodes,2)

    call HydrodynamicInitialize(elems,nodes)

    call magnetic_init(elems, nodes, B0_func, E0_func, sigma_func)

    ! time step
    dt = 0.1d0
    call getarg(2, arg)
    read(arg,*) dt
    do while(t<Tend-1d-8)
        ! Step 1: Get Lorentz force
        call magnetic_lorentz(F)

        ! Step 2: Calculate velocity at t+dt
        call getVelocity(u_func, u)

        ! Step 3: Solve diffusion equation, advance time
        call magnetic_diffusion(dt, g_func, f_func)

        ! Step 4: Move mesh
        call MoveMesh(u, dt, nodes)
    end do

    ! compute error
    call magnetic_errorB(B_func, L2_err, relL2err)
    write(*,*) "L2 error of B = ", L2_err
    write(*,*) "Relative L2 error of B = ", relL2err
    
    call magnetic_plot("output/static/B.vtk")
    

end program test_static