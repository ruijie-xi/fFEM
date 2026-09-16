PROGRAM test_interpolation
    use settings
    use mesh, only: MeshInit, MeshFree
    use fe, only: fespaceInit, Interpolate, fespaceFree
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    
    implicit none

    type(mesh2D) :: Th
    type(fespace) :: Vh
    type(vector) :: u 
    procedure(func) :: test_func

    real(8) :: value
    real(8), dimension(:), allocatable :: values

    integer :: Nx,Ny
    integer :: i_loop, n_levels
    character(32) :: arg

    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    integer, parameter :: mesh_type = MESH_QUAD
    integer :: DOF_type = DOF_QuadRT1
    integer, parameter :: Gauss_type = QuadPt9
    integer, parameter :: fe_dim = 2

    Nx = 10
    Ny = 10

    n_levels = 3
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg,*) n_levels
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg,*) DOF_type
    end if
    if (n_levels < 1 .or. n_levels > 6) error stop 'Expected 1 to 6 refinement levels'
    if (DOF_type /= DOF_QuadRT1 .and. DOF_type /= DOF_QuadRT2) error stop 'Expected RT1 (121) or RT2 (122)'

    do i_loop = 1,n_levels
        write(*,*) "Nx = ", Nx, "Ny = ", Ny

        if (mesh_type == MESH_TRIANGLE) then
            call TriangleMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
        elseif (mesh_type == MESH_QUAD) then
            call SquareMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
        end if
        call MeshInit(elems,nodes,Th)
    
        call fespaceInit(Vh, Th, DOF_type, fe_dim)
        call Interpolate(u, test_func, Th, Vh)
        call ComputeIntegral(u, Th, Vh, Gauss_type, values)
        write(*,*) "Integral", values
        call ComputeNorm(u, Th, Vh, NORM_L2, Gauss_type, value)
        write(*,*) "Norm L2 ", value
        call ComputeNorm(u, Th, Vh, NORM_INF, Gauss_type, value)
        write(*,*) "Norm Linf ", value
        call ComputeError(test_func, u, Th, Vh, NORM_L2, Gauss_type, value)
        write(*,*) "Error L2", value
        call ComputeError(test_func, u, Th, Vh, NORM_HDIV, Gauss_type, value)
        write(*,*) "Error Hdiv", value

        call fespaceFree(Vh)
        call MeshFree(Th)

        Nx = Nx*2
        Ny = Ny*2
    end do

END PROGRAM test_interpolation

subroutine test_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    select case(deriv_type)
    case(DERIV_NONE)
        if (.not. allocated(f)) allocate(f(2))
        f(1) = x(1) + x(2)*x(2)
        f(2) = exp(x(1) + x(2)*x(2))
    case(DERIV_DX)
        if (.not. allocated(f)) allocate(f(2))
        f(1) = 1d0
        f(2) = exp(x(1) + x(2)*x(2))
    case(DERIV_DY)
        if (.not. allocated(f)) allocate(f(2))
        f(1) = 2d0*x(2)
        f(2) = exp(x(1) + x(2)*x(2))*2d0*x(2)
    case(DERIV_DIV)
        if (.not. allocated(f)) allocate(f(1))
        f(1) = 1d0 + exp(x(1) + x(2)*x(2))*2d0*x(2)

    end select
    ! f(1) = exp(x(1) + x(2)**2d0)
end subroutine test_func