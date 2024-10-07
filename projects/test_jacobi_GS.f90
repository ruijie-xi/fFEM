module my_module
    use settings
    implicit none
    
    real(8) :: coeff = 1d0
    
contains

subroutine DirichletBC(A, b, Th, Vh, bndy_func)
    use matvec
    use fe
    type(MATRIX_TRIPLET) :: A
    real(8), dimension(:) :: b
    type(mesh2D) :: Th
    type(fespace) :: Vh
    procedure(func) :: bndy_func

    integer :: i_edge,i_dof,i_bdry
    real(8) :: DbndyVal(1)
    real(8), parameter :: large = 1d0
    integer, dimension(:), allocatable :: dof_index

    real(8), dimension(1) :: x
    integer, dimension(1) :: ind

    do i_bdry = 1, Th%N_bdryedge
        i_edge = Th%BdryEdge(i_bdry)
        if (Th%BdryMarker(i_bdry) == 0) then
            call getEdgeDofIndex(Th, Vh, i_edge, dof_index)
            do i_dof = 1,size(dof_index)
                call ComputeDof(DbndyVal(1), bndy_func, Th, Vh, dof_index(i_dof))
                ind = dof_index(i_dof)
                x(1) = large
                call MatrixTripletAddValues(A, ind, ind, x)
                
                x(1) = DbndyVal(1)*large
                b(ind(1)) = b(ind(1))+x(1)
            end do
        end if
    end do
    
end subroutine DirichletBC

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

end module my_module

program test_jacobi_GS
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
    use solver_petsc
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
    
    ! functions
    procedure(func) :: u_func, rhs_func, one_func, NeumannBdry_func, g_func

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
    integer :: RSS

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
    deallocate(assemble_info)
    allocate(assemble_info(1,3))
    assemble_info(1,:) = (/1, 1, DERIV_NONE/)
    call AssembleVectorElement(rhs_func, [1d0], Th, Vh, 0, assemble_info, Gauss_type, b)
    call timer_end(t_test)
    write(*,*) "Assemble Done. Time taken = ", t_test

    ! get memory usage
    call system_mem_usage(RSS)
    write(*,*) "Memory usage = ", RSS

    ! Dirichlet boundary condition
    call timer_start()
    ! call DirichletBC(A, b, Th, Vh, u_func)
    call timer_end(t_test)
    write(*,*) "Dirichlet Done. Time taken = ", t_test

    ! solve the linear system
    block 
        type(MATRIX_COLUMN) :: A_col
        real(8), dimension(:), allocatable :: x0,r
        integer :: i_step
        real(8) :: res, res_old, rate
        
        allocate(x(Vh%N_DOF), x0(Vh%N_DOF), r(Vh%N_DOF))
        x = 0d0
        x0 = 0d0
        r = 0d0
        b = 0d0
        
        ! select initial value 
        coeff = coeff_input
        call Interpolate(x0, x0_func, Th, Vh)
        
        ! set outputfile
        open(unit=10, file=trim(outputdir)//"/residue.dat", status='replace')
        
        call MatrixTriplet2Column(A, A_col)

        ! write initial residue
        r = b
        call AddMultMV(-1d0, A_col, x0, r)
        ! res = sqrt(sum(r**2))
        res = maxval(abs(r))
        write(10,*) 0, res
        write(*,*) "Step = 0, Residual = ", res
        
        ! begin smoothing
        call timer_start()
        do i_step = 1, max_steps
            
            ! smooth
            if(smooth_type == 1) then
                call jacobi_smooth(A_col, b, x0, x)
            elseif(smooth_type == 2) then
                call GS_smooth(A_col, b, x0, x)
            end if
            
            ! compute residue
            res_old = res
            r = b
            call AddMultMV(-1d0, A_col, x, r)
            ! res = sqrt(sum(r**2))
            res = maxval(abs(r))
            
            ! output residue
            rate = res/res_old
            write(*,*) "Step = ", i_step, "Residual = ", res, "Rate = ", rate
            write(10,*) i_step, res
            
            ! save old solution
            x0 = x
            
        end do
        call timer_end(t_test)
        write(*,*) "Solve Done. Time taken = ", t_test
        close(10)
    end block
    
    
    ! write solution to file
    call PlotFunction(x, Th, Vh, trim(outputdir)//"/u.vtk")
    write(*,*) "Solution written to file ", trim(outputdir)//"/u.vtk"

end program test_jacobi_GS

subroutine rhs_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 5d0*sin(x(1))*sin(2d0*x(2))
    end select
end subroutine rhs_func

subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sin(x(1))*sin(2d0*x(2))
    case(DERIV_DX)
        f(1) = cos(x(1))*sin(2d0*x(2))
    case(DERIV_DY)
        f(1) = 2d0*sin(x(1))*cos(2d0*x(2))
    end select
end subroutine u_func

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

