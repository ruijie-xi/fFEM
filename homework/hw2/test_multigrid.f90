module my_module
    use settings
    use timer
    implicit none
        
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
        ! f(1) = ((coeff*pi)**2+(coeff*pi)**2)*sin(coeff*pi*x(1))*sin(coeff*pi*x(2))
        f(1) = -2d0*((1d0-6d0*x(1)**2)*x(2)**2*(1d0-x(2)**2) +&
         (1d0-6d0*x(2)**2)*x(1)**2*(1d0-x(1)**2))
        ! f(1) = 0d0
        ! f(1) = 2.0d0*(2.0d0-x(1)*x(1)-x(2)*x(2))
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
        ! f(1) = sin(coeff*pi*x(1))*sin(coeff*pi*x(2))
        f(1) = x(1)**2*(1d0-x(1)**2)*x(2)**2*(1d0-x(2)**2)
        ! f(1) = 0d0
        ! f(1) = (x(1)*x(1)-1.0d0)*(x(2)*x(2)-1.0d0)
    case(DERIV_DX)
        ! f(1) = coeff*pi*cos(coeff*pi*x(1))*sin(coeff*pi*x(2))
        f(1) = (2d0*x(1)-4d0*x(1)**3)*(1d0-x(2)**2)*x(2)**2
        ! f(1) = 0d0
        ! f(1) = 2.0d0*x(1)*(x(2)*x(2)-1.0d0)
    case(DERIV_DY)
        ! f(1) = coeff*pi*sin(coeff*pi*x(1))*cos(coeff*pi*x(2))
        f(1) = (2d0*x(2)-4d0*x(2)**3)*(1d0-x(1)**2)*x(1)**2
        ! f(1) = 0d0
        ! f(1) = 2.0d0*x(2)*(x(1)*x(1)-1.0d0)
    end select
end subroutine u_func


subroutine DirichletBC(A, x, b, Th, Vh, bndy_func)
    use ffem_quicksort, only: unique
    use settings
    use matvec
    use fe
    type(MATRIX_TRIPLET) :: A, Ae
    type(VECTOR) :: x,b
    type(mesh2D) :: Th
    type(fespace) :: Vh
    procedure(func) :: bndy_func
    
    real(8),parameter :: diag=1d0
    
    real :: t_test

    integer :: i_edge,i_dof,i_bdry
    integer, dimension(:), allocatable :: dof_index

    integer :: ind
    
    integer, dimension(:), allocatable :: dof_list 
    
    allocate(dof_list(2*Vh%N_DOF))
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
    
    x%data = 0d0
    do i_dof = 1,size(dof_list)
        ind = dof_list(i_dof)
        x%data(ind) = ComputeDof(bndy_func, Th, Vh, ind)
    end do
    
    ! clear row
    call timer_start()
    call MatrixTripletClearRow(A, dof_list)
    call timer_end(t_test)
    
    ! clear column 
    call timer_start()
    Ae = MatrixTripletEliminateColumn(A, dof_list, x, b)
    call timer_end(t_test)
    
    call timer_start()
    ! set diagonal
    do i_dof = 1,size(dof_list)
        ind = dof_list(i_dof)
        call MatrixTripletAddValue(A, ind, ind, diag)
        b%data(ind) = x%data(ind)*diag
    end do
    call timer_end(t_test)
    
end subroutine DirichletBC


end module my_module

program test_multigrid
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
    use module_multigrid
    use visualize
    
    implicit none

    ! mesh
    type(MESH2D), dimension(:), allocatable :: mesh_list
    integer, dimension(:), allocatable :: N_list 
    integer :: n_level 

    ! fe space
    type(fespace), dimension(:), allocatable :: fes_list 
    integer, parameter :: DOF_type = DOF_P1
    integer, parameter :: Gauss_type = TrianglePt9
    
    ! command line arguments
    character(len=100) :: arg

    ! assemble matrices
    integer, dimension(:,:), allocatable :: assemble_info
    type(MATRIX_TRIPLET), dimension(:), allocatable :: Mat_list 
    type(MATRIX_TRIPLET), dimension(:), allocatable :: Mass_list
    type(VECTOR) :: b
    type(VECTOR) :: x
    
    integer :: i, n_arg, n_pre, n_post, smoother
    character(len=100) :: filename
    
    ! command line argument: ./test_multigrid nlevel N1, ..., Nn, pre, post, smoother
    if(command_argument_count() < 6) then
        write(*,*) "Usage: ./test_twogrid N0 N1 ... Nn filename"
        stop
    end if
    
    n_arg = command_argument_count()
    
    call getarg(1,arg)
    read(arg,*) n_level
    allocate(N_list(n_level))
    do i = 1, n_level
        call getarg(i+1,arg)
        read(arg,*) N_list(i)
        write(*,*) "N(",i,") = ", N_list(i)
    end do
    
    call getarg(n_level+2,arg)
    read(arg,*) n_pre
    call getarg(n_level+3,arg)
    read(arg,*) n_post
    call getarg(n_level+4,arg)
    read(arg,*) smoother
    call getarg(n_level+5,arg)
    read(arg,"(A)") filename
    
    write(*,*) "n_pre = ", n_pre
    write(*,*) "n_post = ", n_post
    write(*,*) "smoother = ", smoother
    write(*,*) "filename = ", filename
    
    
    ! generate mesh
    mesh_generation: block 
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes
        
    allocate(mesh_list(n_level))
    do i = 1, n_level
        call TriangleMesh(0d0,1d0,0d0,1d0,N_list(i),N_list(i),elems,nodes)
        call MeshInit(elems, nodes, mesh_list(i))
    end do

    ! initialize fe space
    allocate(fes_list(n_level))
    do i = 1, n_level
        call fespaceInit(fes_list(i), mesh_list(i), DOF_type, 1)
    end do
    write(*,*) "Number of DOF = ", fes_list(n_level)%N_DOF
    
    end block mesh_generation
    
    ! assemble rhs and matrix
    
    allocate(Mat_list(n_level), Mass_list(n_level))
    
    allocate(assemble_info(2,5))
    assemble_info(1,:) = (/1, 1, DERIV_DX, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, DERIV_DY, 1, DERIV_DY/)
    do i = 1, n_level
        call MatrixTripletInit(Mat_list(i), fes_list(i)%N_DOF, fes_list(i)%N_DOF, &
        5*mesh_list(i)%N_elem*fes_list(i)%N_local_basis*fes_list(i)%N_local_basis)
        
        call AssembleMatrixElement(one_func, [1d0,1d0], mesh_list(i), fes_list(i), 0, &
        fes_list(i), 0, assemble_info, Gauss_type, Mat_list(i))
    end do
    
    deallocate(assemble_info)
    allocate(assemble_info(1,5))
    assemble_info(1,:) = (/1, 1, DERIV_NONE, 1, DERIV_NONE/)
    do i = 1, n_level
        call MatrixTripletInit(Mass_list(i), fes_list(i)%N_DOF, fes_list(i)%N_DOF, &
        2*mesh_list(i)%N_elem*fes_list(i)%N_local_basis*fes_list(i)%N_local_basis)
        
        call AssembleMatrixElement(one_func, [1d0], mesh_list(i), fes_list(i), 0, &
        fes_list(i), 0, assemble_info, Gauss_type, Mass_list(i))
    end do
    
    
    deallocate(assemble_info)
    allocate(assemble_info(1,3))
    assemble_info(1,:) = (/1, 1, DERIV_NONE/)
    
    call b%Init(fes_list(n_level)%N_DOF)
    call AssembleVectorElement(rhs_func, [1d0], mesh_list(n_level), fes_list(n_level), 0, assemble_info, Gauss_type, b)
    
    ! Dirichlet BC
    call x%Init(fes_list(n_level)%N_DOF)
    call DirichletBC(Mat_list(n_level), x, b, mesh_list(n_level), fes_list(n_level), u_func)
        
    block
        type(VECTOR) :: bc 
        call bc%Init(fes_list(1)%N_DOF)
        do i = 1, n_level
            call bc%Reset(fes_list(i)%N_DOF)
            if(i<n_level) then
                call DirichletBC(Mat_list(i), bc, bc, mesh_list(i), fes_list(i), zero_func)
            end if
            call DirichletBC(Mass_list(i), bc, bc, mesh_list(i), fes_list(i), zero_func)
        end do
    end block
            
    ! solve
    solve: block
    type(MATRIX_COLUMN), dimension(:), allocatable :: Mat_col_list
    type(MATRIX_COLUMN), dimension(:), allocatable :: Mass_col_list
    real(8) :: H1error, L2error 
    type(GridTransfer), dimension(:), allocatable :: Tr_list 
    
    allocate(Mat_col_list(n_level), Mass_col_list(n_level))
    allocate(Tr_list(n_level-1))
    
    do i = 1, n_level
        call MatrixTriplet2Column(Mat_list(i), Mat_col_list(i))
        call MatrixTriplet2Column(Mass_list(i), Mass_col_list(i))
    end do
    
    do i = 1, n_level-1
        call Tr_list(i)%SetSpace(mesh_list(i+1), fes_list(i+1), Mass_col_list(i+1), &
         mesh_list(i), fes_list(i), Mass_col_list(i))
    end do
    
    open(unit=10, file=filename, status="unknown")
    call solver_multigrid(mesh_list, fes_list, Mat_col_list, Tr_list, b, x, 1d-8, 100, &
    smoother, n_pre, n_post, .true., 10)
    close(10)
    
    call ComputeError(u_func, x, mesh_list(n_level), fes_list(n_level), NORM_H1, Gauss_type, H1error)
    call ComputeError(u_func, x, mesh_list(n_level), fes_list(n_level), NORM_L2, Gauss_type, L2error)
    write(*,*) "H1 error = ", H1error
    write(*,*) "L2 error = ", L2error
    
    end block solve
    
end program test_multigrid



