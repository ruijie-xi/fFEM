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
    real(8), allocatable, dimension(:) :: u
    procedure(func) :: test_func

    real(8) :: value
    real(8), dimension(:), allocatable :: values

    integer :: Nx,Ny
    integer :: i_loop

    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    integer, parameter :: mesh_type = MESH_TRIANGLE
    integer, parameter :: DOF_type = DOF_P1
    integer, parameter :: Gauss_type = TrianglePt9
    integer, parameter :: fe_dim = 2

    Nx = 10
    Ny = 10

    do i_loop = 1,8
        write(*,*) "Nx = ", Nx, "Ny = ", Ny

        if (mesh_type == MESH_TRIANGLE) then
            call TriangleMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
        elseif (mesh_type == MESH_QUAD) then
            call SquareMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
        end if
        call MeshInit(elems,nodes,Th)
    
        call fespaceInit(Vh, Th, DOF_type, fe_dim)
        call Interpolate(u, test_func, Th,Vh)
        call ComputeIntegral(u, Th, Vh, Gauss_type, values)
        write(*,*) "Integral", values
        call ComputeNorm(u, Th, Vh, NORM_L2, Gauss_type, value)
        write(*,*) "Norm L2 ", value
        call ComputeError(test_func, u, Th, Vh, NORM_L2, Gauss_type, value)
        write(*,*) "Error L2", value
        call ComputeError(test_func, u, Th, Vh, NORM_H1, Gauss_type, value)
        write(*,*) "Error H1", value

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