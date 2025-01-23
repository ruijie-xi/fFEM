module my_module
    implicit none
    
contains

subroutine DirichletBC(A, b, Th, Vh, bndy_func)
    use settings
    use matvec
    use fe
    type(MATRIX_TRIPLET) :: A
    type(vector) :: b
    type(mesh2D) :: Th
    type(fespace) :: Vh
    procedure(func) :: bndy_func

    integer :: i_edge,i_dof,i_bdry
    real(8) :: DbndyVal(1)
    real(8), parameter :: large = 1d10
    integer, dimension(:), allocatable :: dof_index

    real(8), dimension(1) :: x
    integer, dimension(1) :: ind

    do i_bdry = 1, Th%N_bdryedge
        i_edge = Th%BdryEdge(i_bdry)
        if (Th%BdryMarker(i_bdry) == 0) then
            call getEdgeDofIndex(Th, Vh, i_edge, dof_index)
            do i_dof = 1,size(dof_index)
                DbndyVal(1) = ComputeDof(bndy_func, Th, Vh, dof_index(i_dof))
                ind = dof_index(i_dof)
                x(1) = large
                call MatrixTripletAddValues(A, ind, ind, x)
                
                x(1) = DbndyVal(1)*large
                b%data(ind(1)) = b%data(ind(1))+x(1)
            end do
        end if
    end do

end subroutine DirichletBC

end module my_module

program test_possion
    use settings
    use mesh, only: MeshInit, PrintMesh, AddBdryMarker
    use fe, only: fespaceInit, Interpolate
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    use matvec
    use assembler
    use memory_usage
    use my_module
    use solver_umfpack2
    
    implicit none

    ! mesh
    type(MESH2D) :: Th
    integer :: Nx,Ny
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! fe space
    type(fespace) :: Vh
    integer, parameter :: mesh_type = MESH_QUAD
    integer, parameter :: DOF_type = DOF_Q2
    integer, parameter :: Gauss_type = QuadPt9
    integer, parameter :: Gauss_type_bdry = LinePt2
    
    ! functions
    procedure(func) :: u_func, rhs_func, one_func, NeumannBdry_func, g_func

    ! timer
    real :: t_test

    ! command line arguments
    character(len=100) :: arg

    ! assemble matrices
    integer, dimension(:,:), allocatable :: assemble_info
    type(MATRIX_TRIPLET) :: A
    type(vector) :: x
    type(vector) :: b
    real(8) :: L2error, H1error
    integer :: RSS

    ! get command line arguments
    call getarg(1, arg)
    read(arg,*) Nx
    call getarg(2, arg)
    read(arg,*) Ny

    ! generate mesh
    call timer_start()
    if (mesh_type == MESH_TRIANGLE) then
        call TriangleMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
    elseif (mesh_type == MESH_QUAD) then
        call SquareMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
    end if
    call MeshInit(elems,nodes,Th)
    call AddBdryMarker(Th, NeumannBdry_func, 2)
    call timer_end(t_test)
    write(*,*) "Mesh generation done. Time taken = ", t_test

    ! initialize fe space
    call timer_start()
    call fespaceInit(Vh, Th, DOF_type, 2)
    call timer_end(t_test)
    write(*,*) "Initialize Done. Time taken = ", t_test
    write(*,*) "Number of DOF = ", Vh%N_DOF

    ! assemble rhs and matrix
    call MatrixTripletInit(A, Vh%N_DOF, Vh%N_DOF, 5*Th%N_elem*Vh%N_local_basis*Vh%N_local_basis)
    call timer_start()
    allocate(assemble_info(4,5))
    assemble_info(1,:) = (/1, 1, DERIV_DX, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, DERIV_DY, 1, DERIV_DY/)
    assemble_info(3,:) = (/1, 2, DERIV_DX, 2, DERIV_DX/)
    assemble_info(4,:) = (/1, 2, DERIV_DY, 2, DERIV_DY/)
    call AssembleMatrixElement(one_func, [1d0,1d0,1d0,1d0], Th, Vh, 0, Vh, 0, assemble_info, Gauss_type, A)
    deallocate(assemble_info)
    
    call VectorInit(b, Vh%N_DOF)
    allocate(assemble_info(2,3))
    assemble_info(1,:) = (/1, 1, DERIV_NONE/)
    assemble_info(2,:) = (/2, 2, DERIV_NONE/)
    call AssembleVectorElement(rhs_func, [1d0,1d0], Th, Vh, 0, assemble_info, Gauss_type, b)
    deallocate(assemble_info)
    allocate(assemble_info(2,3))
    assemble_info(1,:) = (/1, 1, DERIV_NONE/)
    assemble_info(2,:) = (/2, 2, DERIV_NONE/)
    call AssembleVectorBdry(2, g_func, [1d0,1d0], Th, Vh, 0, assemble_info, Gauss_type_bdry, b)
    call timer_end(t_test)
    write(*,*) "Assemble Done. Time taken = ", t_test

    ! get memory usage
    call system_mem_usage(RSS)
    write(*,*) "Memory usage = ", RSS

    ! Dirichlet boundary condition
    call timer_start()
    call DirichletBC(A, b, Th, Vh, u_func)
    call timer_end(t_test)
    write(*,*) "Dirichlet Done. Time taken = ", t_test

    ! solve the linear system
    call VectorInit(x, Vh%N_DOF)
    call timer_start()
    block
        type(MATRIX_COLUMN) :: A_col
        call MatrixTriplet2Column(A, A_col)
        call SolverSolveUMFPACK2(A_col, b, x)
    end block
    call timer_end(t_test)
    write(*,*) "Solve Done. Time taken = ", t_test
    
    ! write solution to file
    call PlotFunction(x, Th, Vh, "output/u.vtk")
    write(*,*) "Solution written to output/u.vtk"

    ! compute error
    call ComputeError(u_func, x, Th, Vh, NORM_L2, Gauss_type, L2error)
    write(*,*) "L2 error = ", L2error
    call ComputeError(u_func, x, Th, Vh, NORM_H1, Gauss_type, H1error)
    write(*,*) "H1 error = ", H1error

end program test_possion

subroutine rhs_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 5d0*sin(x(1))*sin(2d0*x(2))
        f(2) = - (3d0*exp(x(1)+x(2)*x(2)) + 4d0*x(2)*x(2)*exp(x(1)+x(2)*x(2)))
    end select
end subroutine rhs_func

subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sin(x(1))*sin(2d0*x(2))
        f(2) = exp(x(1)+x(2)*x(2))
    case(DERIV_DX)
        f(1) = cos(x(1))*sin(2d0*x(2))
        f(2) = exp(x(1)+x(2)*x(2))
    case(DERIV_DY)
        f(1) = 2d0*sin(x(1))*cos(2d0*x(2))
        f(2) = 2d0*x(2)*exp(x(1)+x(2)*x(2))
    end select
end subroutine u_func

subroutine one_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 1d0
    end select
end subroutine one_func

subroutine NeumannBdry_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = x(1)-1d0
    end select
end subroutine NeumannBdry_func

subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = cos(x(1))*sin(2d0*x(2))
        f(2) = exp(x(1)+x(2)*x(2))
    end select
end subroutine g_func