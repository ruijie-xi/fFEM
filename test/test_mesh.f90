PROGRAM test_mesh
    use settings
    use mesh
    use fe, only: fespaceInit, Interpolate, fespaceFree
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    
    implicit none

    type(mesh2D) :: Th
    type(fespace) :: Vh
    type(VECTOR) :: u
    procedure(func) :: test_func

    integer :: Nx,Ny

    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    real(8) :: point(2)
    integer :: i_elem, i_point

    integer, parameter :: mesh_type = MESH_QUAD
    integer, parameter :: DOF_type = DOF_Q1
    integer, parameter :: Gauss_type = QuadPt9
    integer, parameter :: fe_dim = 2

    Nx = 10
    Ny = 10

    if (mesh_type == MESH_TRIANGLE) then
        call TriangleMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
    elseif (mesh_type == MESH_QUAD) then
        call SquareMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
    end if
    call MeshInit(elems,nodes,Th)

    call fespaceInit(Vh, Th, DOF_type, fe_dim)
    call Interpolate(u, test_func, Th,Vh)
    
    point = [0.5d0,0.5d0]
    call GetPointElement(Th, point, i_elem)

    print *, "Point: ", point
    print *, "Element: ", i_elem

    print *, "The ", i_elem, "th element is: "
    do i_point = 1, Th%N_le
        print *, "Node ", i_point, ": ", Th%NodeCoord(1, Th%ElemNodeConn(i_point, i_elem)), &
        Th%NodeCoord(2, Th%ElemNodeConn(i_point, i_elem))
    end do
    
    call fespaceFree(Vh)
    call MeshFree(Th)

END PROGRAM test_mesh

subroutine test_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
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