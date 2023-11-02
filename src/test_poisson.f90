program test_possion
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscpcdef.h>
    use petscvec
    use petscmat
    use petscksp
    use petscpc
    use settings
    use mesh, only: MeshInit, PrintMesh
    use fe, only: fespace_init, Interpolate
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    use solver
    use memory_usage
    use assembler
    
    implicit none

    ! mesh
    type(mesh2D) :: Th
    integer :: Nx,Ny
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! fe space
    type(fespace) :: Vh
    integer, parameter :: DOF_type = DOF_P1
    integer, parameter :: Gauss_type = TrianglePt9
    
    ! functions
    real(8), allocatable, dimension(:) :: u
    procedure(func) :: u_func, rhs_func, one_func

    ! timer
    real :: t_test

    ! command line arguments
    character(len=100) :: arg

    ! petsc
    integer :: ierr
    Mat :: A
    Vec :: b, x
    KSP :: ksp
    PC :: pc
    PetscScalar, pointer :: xx_v(:)

    ! assemble matrices
    integer, dimension(:,:), allocatable :: assemble_info
    real(8) :: L2error, H1error
    integer :: RSS

    ! get command line arguments
    call getarg(1, arg)
    read(arg,*) Nx
    call getarg(2, arg)
    read(arg,*) Ny

    ! generate mesh
    call timer_start()
    call TriangleMesh(0d0,1d0,0d0,1d0,Nx,Ny,elems,nodes)
    call MeshInit(elems,nodes,Th)
    call timer_end(t_test)
    write(*,*) "Mesh generation done. Time taken = ", t_test

    ! initialize fe space and petsc
    call timer_start()
    call fespace_init(Vh, Th, DOF_type, 2)
    call PetscInitialize("petsc_options.dat", ierr)
    call CreateMat(A, Vh%N_DOF, Vh%N_DOF, ierr)
    call CreateVec(b, Vh%N_DOF, ierr)
    call CreateVec(x, Vh%N_DOF, ierr)
    call CreateSolver(ksp, pc, A, ierr)
    call timer_end(t_test)
    write(*,*) "Initialize Done. Time taken = ", t_test


    ! assemble rhs and matrix
    call timer_start()
    allocate(assemble_info(4,6))
    assemble_info(1,:) = (/1, 1, 1, DERIV_DX, 1, DERIV_DX/)
    assemble_info(2,:) = (/1, 1, 1, DERIV_DY, 1, DERIV_DY/)
    assemble_info(3,:) = (/1, 1, 2, DERIV_DX, 2, DERIV_DX/)
    assemble_info(4,:) = (/1, 1, 2, DERIV_DY, 2, DERIV_DY/)
    call AssembleMatrixElement(one_func, Th, Vh, Vh, assemble_info, Gauss_type, ierr, A)
    deallocate(assemble_info)
    allocate(assemble_info(2,4))
    assemble_info(1,:) = (/1, 1, 1, DERIV_NONE/)
    assemble_info(2,:) = (/1, 2, 2, DERIV_NONE/)
    call AssembleVectorElement(rhs_func, Th, Vh, assemble_info, Gauss_type, ierr, b)
    deallocate(assemble_info)
    call timer_end(t_test)
    write(*,*) "Assemble Done. Time taken = ", t_test

    ! get memory usage
    call system_mem_usage(RSS)
    write(*,*) "Memory usage = ", RSS

    ! Dirichlet boundary condition
    call timer_start()
    call DirichletBC(A, b, Th, Vh, u_func, ierr)
    call timer_end(t_test)
    write(*,*) "Dirichlet Done. Time taken = ", t_test

    ! solve the linear system
    call timer_start()
    call Solve(A, b, ksp, pc, x, ierr)
    call timer_end(t_test)
    write(*,*) "Solve Done. Time taken = ", t_test
    
    ! write solution to file
    allocate(u(Vh%N_DOF))
    call VecGetArrayReadF90(x, xx_v, ierr)
    call PlotFunction(xx_v, Th, Vh, "output/u.vtk")
    write(*,*) "Solution written to output/u.vtk"

    ! compute error
    call ComputeError(u_func, xx_v, Th, Vh, NORM_L2, Gauss_type, L2error)
    write(*,*) "L2 error = ", L2error
    call ComputeError(u_func, xx_v, Th, Vh, NORM_H1, Gauss_type, H1error)
    write(*,*) "H1 error = ", H1error

    call PetscFinalize(ierr)

end program test_possion

subroutine rhs_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 5d0*pi*pi*sin(pi*x(1))*sin(2d0*pi*x(2))
        f(2) = - (3d0*exp(x(1)+x(2)*x(2)) + 4d0*x(2)*x(2)*exp(x(1)+x(2)*x(2)))
    end select
end subroutine rhs_func

subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8), parameter :: pi = 3.14159265358979323846264d0

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sin(pi*x(1))*sin(2d0*pi*x(2))
        f(2) = exp(x(1)+x(2)*x(2))
    case(DERIV_DX)
        f(1) = pi*cos(pi*x(1))*sin(2d0*pi*x(2))
        f(2) = exp(x(1)+x(2)*x(2))
    case(DERIV_DY)
        f(1) = 2d0*pi*sin(pi*x(1))*cos(2d0*pi*x(2))
        f(2) = 2d0*x(2)*exp(x(1)+x(2)*x(2))
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

subroutine DirichletBC(A, b, Th, Vh, bndy_func, ierr)
    use settings
    use solver
    use fe
    Mat :: A
    Vec :: b
    type(mesh2D) :: Th
    type(fespace) :: Vh
    procedure(func) :: bndy_func
    integer :: ierr

    integer :: i_edge,i_dof
    real(8) :: DbndyVal(1)
    real(8), parameter :: large = 1d10
    integer, dimension(:), allocatable :: dof_index

    real(8), dimension(1) :: x
    integer, dimension(1) :: ind

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    do i_edge = 1, Th%N_edge
        if (Th%EdgeMarker(i_edge) == 1) then
            call getEdgeDofIndex(Th, Vh, i_edge, dof_index)
            do i_dof = 1,size(dof_index)
                call ComputeDof(DbndyVal(1), bndy_func, Th, Vh, dof_index(i_dof))
                ind = dof_index(i_dof)-1
                x(1) = large
                call MatSetValues(A, 1, ind, 1, ind, x, ADD_VALUES, ierr)
                x(1) = DbndyVal(1)*large
                call VecSetValues(b, 1, ind, x, ADD_VALUES, ierr)
            end do
        end if
    end do

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

end subroutine DirichletBC
