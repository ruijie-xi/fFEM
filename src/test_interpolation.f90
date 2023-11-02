PROGRAM test_interpolation
    use settings
    use mesh, only: MeshInit, PrintMesh
    use fe, only: fespace_init, Interpolate
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    
    implicit none

    type(mesh2D) :: Th
    type(fespace) :: Vh
    real(8), allocatable, dimension(:) :: u
    procedure(func) :: test_func

    real :: t_test
    real(8) :: value
    real(8), dimension(:), allocatable :: values

    integer :: Nx,Ny

    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    character(len=100) :: arg
    integer, parameter :: Gauss_type = QuadPt9

    call getarg(1, arg)
    read(arg,*) Nx
    call getarg(2, arg)
    read(arg,*) Ny

    call timer_start()
    call SquareMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
    call MeshInit(elems,nodes,Th)
    call timer_end(t_test)
    write(*,*) "Time taken = ", t_test

    call timer_start()
    call fespace_init(Vh, Th, DOF_Q1, 2)
    call Interpolate(u, test_func, Th,Vh)
    call ComputeIntegral(u, Th, Vh, Gauss_type, values)
    write(*,*) "Integral", values
    call ComputeNorm(u, Th, Vh, NORM_L2, Gauss_type, value)
    write(*,*) "Norm L2 ", value
    call ComputeError(test_func, u, Th, Vh, NORM_L2, Gauss_type, value)
    write(*,*) "Error L2", value
    call ComputeError(test_func, u, Th, Vh, NORM_H1, Gauss_type, value)
    write(*,*) "Error H1", value
    call timer_end(t_test)

    write(*,*) "Time taken = ", t_test

    call PlotFunction(u, Th, Vh, "test.vtk")

END PROGRAM test_interpolation

subroutine test_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = x(1) + x(2)*x(2)
        f(2) = exp(x(1) + x(2)*x(2))
    case(DERIV_DX)
        f(1) = 1d0
        f(2) = exp(x(1) + x(2)*x(2))
    case(DERIV_DY)
        f(1) = 2d0*x(2)
        f(2) = exp(x(1) + x(2)*x(2))*2d0*x(2)
    end select
    ! f(1) = exp(x(1) + x(2)**2d0)
end subroutine test_func