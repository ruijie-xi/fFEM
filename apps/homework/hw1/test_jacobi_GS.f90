module my_module
    use settings
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
    type(VECTOR) :: b,x

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
    call timer_start()
    allocate(assemble_info(2,5))
    assemble_info(1,:) = (/1, 1, DERIV_DX, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, DERIV_DY, 1, DERIV_DY/)
    call AssembleMatrixElement(one_func, [1d0,1d0], Th, Vh, 0, Vh, 0, assemble_info, Gauss_type, A)
    call timer_end(t_test)
    write(*,*) "Assemble Done. Time taken = ", t_test

    ! solve the linear system
    block 
        type(MATRIX_COLUMN) :: A_col
        type(VECTOR) :: x0,r
        integer :: i_step
        real(8) :: res, res_old, rate
        
        call x%Init(Vh%N_DOF)
        call x0%Init(Vh%N_DOF)
        call r%Init(Vh%N_DOF)
        call b%Init(Vh%N_DOF)
        
        ! select initial value 
        coeff = coeff_input
        call Interpolate(x0%data, x0_func, Th, Vh)
        
        ! set outputfile
        open(unit=10, file=trim(outputdir)//"/residue.dat", status='replace')
        
        call MatrixTriplet2Column(A, A_col)

        ! write initial residue
        r = b
        call AddMultMV(-1d0, A_col, x0, 1d0, r)
        res = r%Norm()
        write(10,*) 0, res
        write(*,*) "Step = 0, Residual = ", res
        
        ! begin smoothing
        call timer_start()
        do i_step = 1, max_steps
            
            ! smooth
            if(smooth_type == 1) then
                call jacobi_smooth(A_col, b%data, x0%data, x%data)
            elseif(smooth_type == 2) then
                call GS_smooth(A_col, b%data, x0%data, x%data)
            end if
            
            ! compute residue
            res_old = res
            r = b
            call AddMultMV(-1d0, A_col, x, 1d0, r)
            ! res = sqrt(sum(r**2))
            res = r%Norm()
            
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

end program test_jacobi_GS



