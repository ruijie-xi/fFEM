module my_module
    use settings
    use timer
    implicit none
    
    real(8) :: coeff = 1d0
    
contains

subroutine one_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable, intent(out) :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 1d0
    end select
end subroutine one_func


subroutine zero_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable, intent(out) :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine zero_func

subroutine rhs_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable, intent(out) :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = ((2d0*pi)**2+(2d0*pi)**2)*sin(2d0*pi*x(1))*sin(2d0*pi*x(2))
    end select
end subroutine rhs_func

subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable, intent(out) :: f
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
end subroutine u_func


subroutine DirichletBC(A, x, b, Th, Vh, bndy_func)
    use ffem_quicksort, only: unique
    use settings
    use matvec
    use fe
    type(MATRIX_TRIPLET) :: A, Ae
    real(8), dimension(:), intent(inout) :: b,x
    type(mesh2D) :: Th
    type(fespace) :: Vh
    procedure(func) :: bndy_func
    
    real(8),parameter :: diag=1d0
    
    real :: t_test

    integer :: i_edge,i_dof,i_bdry
    integer, dimension(:), allocatable :: dof_index

    integer :: ind
    
    integer, dimension(:), allocatable :: dof_list 
    
    allocate(dof_list(Vh%N_DOF))
    dof_list = 0
    
    i_dof = 0

    do i_bdry = 1, Th%N_bdryedge
        i_edge = Th%BdryEdge(i_bdry)
        if (Th%BdryMarker(i_bdry) == 0) then
            call getEdgeDofIndex(Th, Vh, i_edge, dof_index)
            
            dof_list(i_dof+1:i_dof+size(dof_index)) = dof_index
            
            i_dof = i_dof + size(dof_index)
        end if
    end do
        
    dof_list = unique(dof_list(1:i_dof))
    
    x = 0d0
    do i_dof = 1,size(dof_list)
        ind = dof_list(i_dof)
        x(ind) = ComputeDof(bndy_func, Th, Vh, ind)
    end do
    
    ! clear row
    call timer_start()
    call MatrixTripletClearRow(A, dof_list)
    call timer_end(t_test)
    write(*,*) "Clear Row Done. Time taken = ", t_test
    
    ! clear column 
    call timer_start()
    Ae = MatrixTripletEliminateColumn(A, dof_list, x, b)
    call timer_end(t_test)
    write(*,*) "Eliminate Column Done. Time taken = ", t_test
    
    call timer_start()
    ! set diagonal
    do i_dof = 1,size(dof_list)
        ind = dof_list(i_dof)
        call MatrixTripletAddValue(A, ind, ind, diag)
        b(ind) = x(ind)*diag
    end do
    call timer_end(t_test)
    write(*,*) "Set Diagonal Done. Time taken = ", t_test
    
    

end subroutine DirichletBC

end module my_module
    

module Operator
    use settings
    use module_prolongation
    implicit none
    
    private
    
    public :: OpR, OpP, SetSpace
    
    type(MESH2D), pointer :: Th_fine, Th_coarse
    type(FESPACE), pointer :: Vh_fine, Vh_coarse
    
contains

subroutine SetSpace(Th_fine_, Vh_fine_, Th_coarse_, Vh_coarse_)
    type(MESH2D), target :: Th_fine_, Th_coarse_
    type(FESPACE), target :: Vh_fine_, Vh_coarse_
    
    Th_fine => Th_fine_
    Th_coarse => Th_coarse_
    Vh_fine => Vh_fine_
    Vh_coarse => Vh_coarse_
    
end subroutine

subroutine OpR(x, y)
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: y
    
    call OperatorTransfer(x, Th_fine, Vh_fine, Th_coarse, Vh_coarse, y)
end subroutine

subroutine OpP(x, y)
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: y
    
    call OperatorTransfer(x, Th_coarse, Vh_coarse, Th_fine, Vh_fine, y)
end subroutine

end module

program test_twogrid
    use settings
    use mesh, only: MeshInit, PrintMesh, AddBdryMarker
    use fe, only: fespaceInit, Interpolate
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use matvec
    use assembler
    use my_module
    use solver_umfpack2
    use module_smoother
    use module_prolongation
    use Operator
    use module_twogrid
    
    implicit none

    ! mesh
    type(MESH2D) :: Th_fine, Th_coarse
    integer :: N_fine, N_coarse
    integer, dimension(:,:), allocatable :: elems_coarse, elems_fine
    real(8), dimension(:,:), allocatable :: nodes_coarse, nodes_fine

    ! fe space
    type(fespace) :: Vh_fine, Vh_coarse
    integer, parameter :: DOF_type = DOF_P1
    integer, parameter :: Gauss_type = TrianglePt4
    
    ! timer
    real :: t_test

    ! command line arguments
    character(len=100) :: arg

    ! assemble matrices
    integer, dimension(:,:), allocatable :: assemble_info
    type(MATRIX_TRIPLET) :: A, Ac
    real(8), dimension(:), allocatable :: b,bc
    
    real(8), dimension(:), allocatable :: x_fine, x_coarse
    
    integer :: i
    
    if(command_argument_count() < 2) then
        write(*,*) "Usage: ./test_twogrid N_coarse N_fine"
        stop
    end if
    call getarg(1,arg)
    read(arg,*) N_coarse
    call getarg(2,arg)
    read(arg,*) N_fine
    
    write(*,*) "N_coarse = ", N_coarse
    write(*,*) "N_fine = ", N_fine

    ! generate mesh
    call timer_start()
    call TriangleMesh(0d0,1d0,0d0,1d0,N_coarse,N_coarse,elems_coarse,nodes_coarse)
    call TriangleMesh(0d0,1d0,0d0,1d0,N_fine,N_fine,elems_fine,nodes_fine)
    call MeshInit(elems_coarse, nodes_coarse, Th_coarse)
    call MeshInit(elems_fine, nodes_fine, Th_fine)
    call timer_end(t_test)
    write(*,*) "Mesh generation done. Time taken = ", t_test

    ! initialize fe space
    call timer_start()
    call fespaceInit(Vh_coarse, Th_coarse, DOF_type, 1)
    call fespaceInit(Vh_fine, Th_fine, DOF_type, 1)
    call timer_end(t_test)
    write(*,*) "Initialize Done. Time taken = ", t_test
    write(*,*) "Number of DOF = ", Vh_fine%N_DOF

    ! assemble rhs and matrix
    call timer_start()

    call MatrixTripletInit(A, Vh_fine%N_DOF, Vh_fine%N_DOF, 5*Th_fine%N_elem*Vh_fine%N_local_basis*Vh_fine%N_local_basis)
    
    call MatrixTripletInit(Ac, Vh_coarse%N_DOF, Vh_coarse%N_DOF, 5*Th_coarse%N_elem*Vh_coarse%N_local_basis*Vh_coarse%N_local_basis)
    
    allocate(assemble_info(2,5))
    assemble_info(1,:) = (/1, 1, DERIV_DX, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, DERIV_DY, 1, DERIV_DY/)
    call AssembleMatrixElement(one_func, [1d0,1d0], Th_fine, Vh_fine, 0, Vh_fine, 0, assemble_info, Gauss_type, A)
    
    call AssembleMatrixElement(one_func, [1d0,1d0], Th_coarse, Vh_coarse, 0, Vh_coarse, 0, assemble_info, Gauss_type, Ac)
    
    deallocate(assemble_info)
    allocate(assemble_info(1,3))
    assemble_info(1,:) = (/1, 1, DERIV_NONE/)
    
    allocate(b(Vh_fine%N_DOF))
    b=0d0
    call AssembleVectorElement(rhs_func, [1d0], Th_fine, Vh_fine, 0, assemble_info, Gauss_type, b)
    
    call timer_end(t_test)
    write(*,*) "Assemble Done. Time taken = ", t_test
    
    ! Dirichlet BC
    allocate(x_fine(Vh_fine%N_DOF), x_coarse(Vh_coarse%N_DOF))
    call DirichletBC(A, x_fine, b, Th_fine, Vh_fine, u_func)
    
    allocate(bc(Vh_coarse%N_DOF))
    bc = 0d0
    call DirichletBC(Ac, x_coarse, bc, Th_coarse, Vh_coarse, zero_func)
    
    ! solve
    solve: block
    
    type(MATRIX_COLUMN) :: A_col, Ac_col
    real(8) :: H1error, L2error 
    
    call MatrixTriplet2Column(A, A_col)
    call MatrixTriplet2Column(Ac, Ac_col)
    
    call SetSpace(Th_fine, Vh_fine, Th_coarse, Vh_coarse)
    
    call solver_twogrid(A_col, Ac_col, OpR, OpP, b, x_fine, 1d-8, 100, SMOOTHER_JACOBI, 10, 10, .true.)
    
    call ComputeError(u_func, x_fine, Th_fine, Vh_fine, NORM_H1, Gauss_type, H1error)
    call ComputeError(u_func, x_fine, Th_fine, Vh_fine, NORM_L2, Gauss_type, L2error)
    write(*,*) "H1 error = ", H1error
    write(*,*) "L2 error = ", L2error
    
    end block solve

end program test_twogrid



