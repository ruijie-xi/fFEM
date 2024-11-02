module my_module
    use settings
    use timer
    implicit none
    
    real(8) :: coeff = 1d0
    
contains

subroutine x0_func(x,f,deriv_type)
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable , intent(out):: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sin(coeff*pi*x(1))*sin(coeff*pi*x(2))
    end select
end subroutine x0_func

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
    
    real(8),parameter :: diag=1d-1
    
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

program test_jacobi_GS
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
    
    implicit none

    ! mesh
    type(MESH2D) :: Th
    integer :: Nx,Ny
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! fe space
    type(fespace) :: Vh
    integer, parameter :: mesh_type = MESH_TRIANGLE
    integer, parameter :: DOF_type = DOF_P2
    integer, parameter :: Gauss_type = TrianglePt9
    
    ! timer
    real :: t_test
    
    ! initial value
    real(8) :: coeff_input = 1d0
    
    ! max iteration
    integer :: max_steps = 10

    ! command line arguments
    character(len=100) :: arg
    
    ! outputfile name
    character(len=100) :: outputdir
    
    ! smoothing type
    ! 1 for Jacobi, 2 for GS
    integer :: smooth_type = 1

    ! assemble matrices
    integer, dimension(:,:), allocatable :: assemble_info
    type(MATRIX_TRIPLET) :: A
    real(8), dimension(:), allocatable :: b,x

    ! get command line arguments
    if(iargc() /= 6) then
        write(*,*) "Usage: test_jacobi_GS Nx Ny coeff maxsteps outputfile smooth_type"
        stop
    end if
    
    call getarg(1, arg)
    read(arg,*) Nx
    call getarg(2, arg)
    read(arg,*) Ny
    call getarg(3, arg)
    read(arg,*) coeff_input
    call getarg(4, arg)
    read(arg,*) max_steps
    call getarg(5, arg)
    read(arg,"(A)") outputdir
    call getarg(6, arg)
    read(arg,*) smooth_type
    
    ! print input
    write(*,*) "Nx = ", Nx
    write(*,*) "Ny = ", Ny
    write(*,*) "coeff = ", coeff_input
    write(*,*) "max_steps = ", max_steps
    write(*,*) "outputdir = ", outputdir
    write(*,*) "smooth_type = ", smooth_type
    

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
    call fespaceInit(Vh, Th, DOF_type, 1)
    call timer_end(t_test)
    write(*,*) "Initialize Done. Time taken = ", t_test
    write(*,*) "Number of DOF = ", Vh%N_DOF

    ! assemble rhs and matrix
    call MatrixTripletInit(A, Vh%N_DOF, Vh%N_DOF, 5*Th%N_elem*Vh%N_local_basis*Vh%N_local_basis)
    allocate(b(Vh%N_DOF))
    b=0d0
    call timer_start()
    allocate(assemble_info(2,5))
    assemble_info(1,:) = (/1, 1, DERIV_DX, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, DERIV_DY, 1, DERIV_DY/)
    call AssembleMatrixElement(one_func, [1d0,1d0], Th, Vh, 0, Vh, 0, assemble_info, Gauss_type, A)
    call timer_end(t_test)
    write(*,*) "Assemble Done. Time taken = ", t_test
    
    b = 0d0
    
    allocate(x(Vh%N_DOF))
    x = 0d0
    
    ! Dirichlet BC
    call DirichletBC(A, x, b, Th, Vh, zero_func)

    ! solve the linear system
    block 
        type(MATRIX_COLUMN) :: A_col
        real(8), dimension(:), allocatable :: x0
        allocate(x0(Vh%N_DOF))
        
        x0 = 0d0
        
        ! select initial value 
        coeff = coeff_input
        call Interpolate(x0, x0_func, Th, Vh)
        
        ! set outputfile
        open(unit=10, file=trim(outputdir)//"/residue.dat", status='replace')
        
        call MatrixTriplet2Column(A, A_col)
        
        ! begin smoothing
        call timer_start()
        if (smooth_type == 1) then
            call Jacobi_smooth(A_col, b, x0, x, max_steps, 10)
        elseif (smooth_type == 2) then
            call GS_smooth(A_col, b, x0, x, max_steps, 10)
        end if
        call timer_end(t_test)
        write(*,*) "Solve Done. Time taken = ", t_test
        close(10)
    end block

end program test_jacobi_GS



