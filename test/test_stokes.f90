! This is a simple test program to solve the Stokes equation
! using the P2-P1 element. The domain is the unit square
! and the exact solution is given by
! u = (2y*exp(x+y^2), -exp(x+y^2))
! p = sin(2*pi*x)*sin(2*pi*y)
! f = -div(u) + grad(p)
! g = 0
! TODO: The treatment of the boundary condition is slow. Need to optimize it.

module mymodule
    use settings
    use fe, only: getEdgeDofIndex, ComputeDof 
    use matvec
    implicit none
    
contains

subroutine DirichletBC(A, b, Th, Vh, Ph, bndy_func, p_func)
    implicit none
    type(MATRIX_COLUMN) :: A
    real(8),dimension(:) :: b
    type(mesh2D) :: Th
    type(fespace) :: Vh, Ph
    procedure(func) :: bndy_func, p_func

    integer :: i_edge,i_dof,i_bdry
    real(8) :: DbndyVal
    integer, dimension(:), allocatable :: dof_index

    do i_bdry = 1, Th%N_bdryedge
        i_edge = Th%BdryEdge(i_bdry)
        if (Th%BdryMarker(i_bdry) == 0) then
            call getEdgeDofIndex(Th, Vh, i_edge, dof_index)
            do i_dof = 1,size(dof_index)
                ! clear the row
                call MatrixColumnClearRow(A, dof_index(i_dof))
                call ComputeDof(DbndyVal, bndy_func, Th, Vh, dof_index(i_dof))
                call MatrixColumnSet(A, dof_index(i_dof), dof_index(i_dof), 1d0)
                b(dof_index(i_dof)) = DbndyVal
            end do
        end if
    end do
    
    ! clear the row
    i_dof = 1
    call MatrixColumnClearRow(A, i_dof + Vh%N_DOF)
    call ComputeDof(DbndyVal, p_func, Th, Ph, i_dof)
    call MatrixColumnSet(A, i_dof + Vh%N_DOF, i_dof + Vh%N_DOF, 1d0)
    b(Vh%N_DOF+i_dof) = DbndyVal

end subroutine DirichletBC
    
end module mymodule

program test_stokes
    use settings
    use mesh, only: MeshInit, AddBdryMarker
    use fe, only: fespaceInit
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    use solver_umfpack2
    use solver_petsc
    use memory_usage
    use assembler
    use matvec
    use mymodule
        
    implicit none

    ! mesh
    type(mesh2D) :: Th
    integer :: Nx,Ny
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! fe space
    type(fespace) :: Vh, Ph
    integer, parameter :: mesh_type = MESH_TRIANGLE
    integer, parameter :: DOF_type_u = DOF_P2, DOF_type_p = DOF_P1
    integer, parameter :: Gauss_type = TrianglePt9
    
    ! functions
    real(8), allocatable, dimension(:) :: x,u,p
    procedure(func) :: u_func, p_func, one_func, f_func, g_func

    ! linear system
    type(MATRIX_TRIPLET) :: A
    type(MATRIX_COLUMN) :: A_column
    real(8), dimension(:), allocatable :: b
    integer :: N_nz

    ! timer
    real :: t_test

    ! command line arguments
    character(len=100) :: arg

    ! assemble matrices
    integer, dimension(:,:), allocatable :: assemble_info
    real(8) :: L2error, H1error
    integer :: RSS

    real(8), dimension(:), allocatable :: p_integral

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
    call timer_end(t_test)
    write(*,*) "Mesh generation done. Time taken = ", t_test

    ! initialize fe space
    call timer_start()
    call fespaceInit(Vh, Th, DOF_type_u, 2)
    call fespaceInit(Ph, Th, DOF_type_p, 1)
    call timer_end(t_test)
    write(*,*) "Initialize Done. Time taken = ", t_test
    write(*,*) "Number of DOF for velocity = ", Vh%N_DOF
    write(*,*) "Number of DOF for pressure = ", Ph%N_DOF

    ! assemble rhs and matrix
    N_nz = Th%N_elem*(Vh%N_local_basis*Vh%N_local_basis*4+&
                    Vh%N_local_basis*Ph%N_local_basis*2+&
                    Ph%N_local_basis*Ph%N_local_basis)+&
                    Th%N_bdryedge*5
    call MatrixTripletInit(A, Vh%N_DOF+Ph%N_DOF, Vh%N_DOF+Ph%N_DOF, 2*N_nz)
    call timer_start()
    allocate(assemble_info(4,5))
    assemble_info(1,:) = (/1, 1, DERIV_DX, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, DERIV_DY, 1, DERIV_DY/)
    assemble_info(3,:) = (/1, 2, DERIV_DX, 2, DERIV_DX/)
    assemble_info(4,:) = (/1, 2, DERIV_DY, 2, DERIV_DY/)
    call AssembleMatrixElement(one_func, [1d0,1d0,1d0,1d0], Th, Vh, 0, Vh, 0, assemble_info, Gauss_type, A)
    deallocate(assemble_info)
    allocate(assemble_info(2,5))
    assemble_info(1,:) = (/1, 1, DERIV_NONE, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, DERIV_NONE, 2, DERIV_DY/)
    call AssembleMatrixElement(one_func, [-1d0,-1d0], Th, Ph, Vh%N_DOF, Vh, 0, assemble_info, Gauss_type, A)
    deallocate(assemble_info)
    allocate(assemble_info(2,5))
    assemble_info(1,:) = (/1, 1, DERIV_DX, 1, DERIV_NONE/)
    assemble_info(2,:) = (/1, 2, DERIV_DY, 1, DERIV_NONE/)
    call AssembleMatrixElement(one_func, [-1d0,-1d0], Th, Vh, 0, Ph, Vh%N_DOF, assemble_info, Gauss_type, A)
    deallocate(assemble_info)
    allocate(assemble_info(1,5))
    assemble_info(1,:) = (/1, 1, DERIV_NONE, 1, DERIV_NONE/)
    call AssembleMatrixElement(one_func, [1d-6], Th, Ph, Vh%N_DOF, Ph, Vh%N_DOF, assemble_info, Gauss_type, A)
    deallocate(assemble_info)

    call VectorInit(b, Vh%N_DOF+Ph%N_DOF)
    allocate(assemble_info(2,3))
    assemble_info(1,:) = (/1, 1, DERIV_NONE/)
    assemble_info(2,:) = (/2, 2, DERIV_NONE/)
    call AssembleVectorElement(f_func, [1d0,1d0], Th, Vh, 0, assemble_info, Gauss_type, b)
    deallocate(assemble_info)
    allocate(assemble_info(1,3))
    assemble_info(1,:) = (/1, 1, DERIV_NONE/)
    call AssembleVectorElement(g_func, [1d0], Th, Ph, Vh%N_DOF, assemble_info, Gauss_type, b)
    call timer_end(t_test)
    write(*,*) "Assemble Done. Time taken = ", t_test

    ! get memory usage
    call system_mem_usage(RSS)
    write(*,*) "Memory usage = ", RSS

    ! Assemble
    call MatrixTriplet2Column(A, A_column)
    call MatrixTripletFree(A)

    ! Dirichlet boundary condition
    call timer_start()
    call DirichletBC(A_column, b, Th, Vh, Ph, u_func, p_func)
    call timer_end(t_test)
    write(*,*) "Dirichlet Done. Time taken = ", t_test

    ! call MatrixColumnTrim(A_column)
    call MatrixColumn2Triplet(A_column, A)

    ! solve the linear system
    call VectorInit(x, Vh%N_DOF+Ph%N_DOF)
    call timer_start()
    call SolverSolvePETSC(A, b, x)
    call timer_end(t_test)
    write(*,*) "Solve Done. Time taken = ", t_test

    allocate(u(Vh%N_DOF), p(Ph%N_DOF))
    u = x(1:Vh%N_DOF)
    p = x(Vh%N_DOF+1:Vh%N_DOF+Ph%N_DOF)

    ! compute integral
    call ComputeIntegral(p, Th, Ph, Gauss_type, p_integral)
    p = p - p_integral(1)
    call ComputeIntegral(p, Th, Ph, Gauss_type, p_integral)
    write(*,*) "Integral of p = ", p_integral
    
    ! write solution to file
    call PlotFunction(u, Th, Vh, "output/u.vtk")
    write(*,*) "Solution written to output/u.vtk"

    ! compute error
    call ComputeError(u_func, u, Th, Vh, NORM_L2, Gauss_type, L2error)
    write(*,*) "L2 error of u = ", L2error
    call ComputeError(u_func, u, Th, Vh, NORM_H1, Gauss_type, H1error)
    write(*,*) "H1 error of u = ", H1error
    call ComputeError(p_func, p, Th, Ph, NORM_L2, Gauss_type, L2error)
    write(*,*) "L2 error of p = ", L2error



end program test_stokes


subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*x(2)*exp(x(1)+x(2)*x(2))
        f(2) = -exp(x(1)+x(2)*x(2))
    case(DERIV_DX)
        f(1) = 2d0*x(2)*exp(x(1)+x(2)*x(2))
        f(2) = -exp(x(1)+x(2)*x(2))
    case(DERIV_DY)
        f(1) = (2d0+4d0*x(2)*x(2))*exp(x(1)+x(2)*x(2))
        f(2) = -2d0*x(2)*exp(x(1)+x(2)*x(2))
    case(DERIV_DXX)
        f(1) = 2d0*x(2)*exp(x(1)+x(2)*x(2))
        f(2) = -exp(x(1)+x(2)*x(2))
    case(DERIV_DYY)
        f(1) = (12d0*x(2)+8d0*x(2)**3)*exp(x(1)+x(2)*x(2))
        f(2) = -(2d0+4d0*x(2)*x(2))*exp(x(1)+x(2)*x(2))
    end select
end subroutine u_func

subroutine p_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sin(2d0*pi*x(1))*sin(2d0*pi*x(2))
    case(DERIV_DX)
        f(1) = 2d0*pi*cos(2d0*pi*x(1))*sin(2d0*pi*x(2))
    case(DERIV_DY)
        f(1) = 2d0*pi*sin(2d0*pi*x(1))*cos(2d0*pi*x(2))
    end select
end subroutine

subroutine f_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0
    real(8), dimension(:), allocatable :: val

    procedure(func) :: u_func, p_func

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        call u_func(x,val,DERIV_DXX)
        f(1) = - val(1)
        f(2) = - val(2)
        deallocate(val)
        call u_func(x,val,DERIV_DYY)
        f(1) = f(1) - val(1)
        f(2) = f(2) - val(2)
        deallocate(val)
        call p_func(x,val,DERIV_DX)
        f(1) = f(1) + val(1)
        deallocate(val)
        call p_func(x,val,DERIV_DY)
        f(2) = f(2) + val(1)
        deallocate(val)
    end select
end subroutine f_func



subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine g_func

subroutine one_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 1d0
    end select
end subroutine one_func