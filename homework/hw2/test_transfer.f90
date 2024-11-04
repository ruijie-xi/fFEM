module my_module
    use settings
    use timer
    implicit none
    
    real(8) :: coeff = 8d0
    
contains


subroutine x0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable, intent(out) :: f
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
        ! f(1) = -2d0*((1d0-6d0*x(1)**2)*x(2)**2*(1d0-x(2)**2) +&
        !  (1d0-6d0*x(2)**2)*x(1)**2*(1d0-x(1)**2))
        f(1) = 0d0
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
        ! f(1) = x(1)**2*(1d0-x(1)**2)*x(2)**2*(1d0-x(2)**2)
        f(1) = 0d0
    case(DERIV_DX)
        ! f(1) = coeff*pi*cos(coeff*pi*x(1))*sin(coeff*pi*x(2))
        ! f(1) = (2d0*x(1)-4d0*x(1)**3)*(1d0-x(2)**2)*x(2)**2
        f(1) = 0d0
    case(DERIV_DY)
        ! f(1) = coeff*pi*sin(coeff*pi*x(1))*cos(coeff*pi*x(2))
        ! f(1) = (2d0*x(2)-4d0*x(2)**3)*(1d0-x(1)**2)*x(1)**2
        f(1) = 0d0
    end select
end subroutine u_func


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
    use module_twogrid
    use visualize
    
    implicit none

    ! mesh
    type(MESH2D), dimension(:), allocatable :: mesh_list
    integer, dimension(:), allocatable :: N_list 
    integer :: n_level 

    ! fe space
    type(fespace), dimension(:), allocatable :: fes_list 
    integer, parameter :: DOF_type = DOF_P1
    integer, parameter :: Gauss_type = TrianglePt4
    
    type(GridTransfer), dimension(:), allocatable :: Trs
    
    ! command line arguments
    character(len=100) :: arg

    type(VECTOR),dimension(:), allocatable :: xs
    
    integer :: i
    
    if(command_argument_count() < 2) then
        write(*,*) "Usage: ./test_twogrid N0 N1 ... Nn"
        stop
    end if
    
    n_level = command_argument_count()
    allocate(N_list(n_level))
    
    do i = 1, n_level
        call getarg(i,arg)
        read(arg,*) N_list(i)
        write(*,*) "N(",i,") = ", N_list(i)
    end do

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
    
    allocate(xs(n_level))
    do i = 1,n_level
        call xs(i)%Init(fes_list(i)%N_DOF)
    end do
    
    call Interpolate(xs(n_level)%data, x0_func, mesh_list(n_level), fes_list(n_level));
    
    allocate(Trs(n_level-1))
    do i = 1, n_level-1
        call Trs(i)%SetSpace(mesh_list(i+1), fes_list(i+1), mesh_list(i), fes_list(i))
    end do
    
    do i = n_level-1, 1, -1
        call Trs(i)%FineToCoarse(xs(i+1), xs(i))
    end do
    
    block
    character(len=100) :: level_name
    do i = 1, n_level
        write(level_name,*) i
        call PlotFunction(xs(i), mesh_list(i), fes_list(i), "level_"//trim(adjustl(level_name))//"solution.vtk")
        print *, xs(i)%data
    end do
    end block
    
    

end program test_multigrid



