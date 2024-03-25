#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscpcdef.h>


module solver_procs
    use settings
    implicit none
    
contains

! scalar one
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

subroutine CurlOperator(Th,E,curlE)
    use settings
    use mesh, only: getEdgeLength
    implicit none
    type(mesh2D),intent(in) :: Th
    real(8), intent(in), dimension(:) :: E
    real(8), intent(out), dimension(:), allocatable :: curlE

    integer :: i_edge,pt1,pt2
    real(8) :: length

    if(allocated(curlE)) deallocate(curlE)
    allocate(curlE(Th%N_Edge))

    do i_edge = 1,Th%N_Edge
        pt1 = Th%EdgeNodeConn(1,i_edge)
        pt2 = Th%EdgeNodeConn(2,i_edge)
        call getEdgeLength(Th, i_edge, length)
        curlE(i_edge) = (E(pt2)-E(pt1))/length
    end do

end subroutine CurlOperator


! boundary function - straight line
subroutine StraightBdry_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = min(abs(x(1)),abs(x(2)))
    end select
end subroutine StraightBdry_func

! boundary function - circle
subroutine CircleBdry_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = abs(sqrt(x(1)**2+x(2)**2) - 1d0)
    end select
end subroutine CircleBdry_func


subroutine UpdateMeshAndSolution(u,dt,Th,B,Bhat)
    use mesh, only: getMeshhmax, getEdgeLength
    implicit none
    type(mesh2D),intent(inout) :: Th
    real(8), intent(in), dimension(:) :: u
    real(8), intent(in) :: dt
    real(8), intent(inout), dimension(:) :: B,Bhat

    real(8), dimension(:), allocatable :: arc_length
    integer :: i_edge

    allocate(arc_length(Th%N_Edge))
    do i_edge = 1,Th%N_Edge
        call getEdgeLength(Th, i_edge, arc_length(i_edge))
    end do

    B = B * arc_length
    Bhat = Bhat * arc_length

    Th%NodeCoord(1,:) = Th%NodeCoord(1,:) + u(1:Th%N_Node)*dt
    Th%NodeCoord(2,:) = Th%NodeCoord(2,:) + u(Th%N_Node+1:2*Th%N_Node)*dt

    call getMeshhmax(Th)

    do i_edge = 1,Th%N_Edge
        call getEdgeLength(Th, i_edge, arc_length(i_edge))
    end do

    B = B / arc_length
    Bhat = Bhat / arc_length

end subroutine UpdateMeshAndSolution

    
end module solver_procs

module problem_static
    use settings
    implicit none

    real(8) :: mu = 4*m_pi
    real(8) :: Center(2) = [0d0,0d0]
    real(8) :: Rlist(1) = [1d0]
    integer :: M = 10
    integer :: Nlist(1) = [15]
    
contains

! right hand side function
subroutine f_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine f_func

! velocity function
subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
        f(2) = 0d0
    end select
end subroutine u_func

! current function
function current(t) result(S)
    implicit none
    real(8), intent(in) :: t
    real(8) :: S
    real(8),parameter :: S0 = 1.2d0,b = 0.5d0

    S = S0*(b/t)**(5d0/2d0)*(2d0*b/t-3d0)*exp(-b/t)
end function current
 
! magnetic field on boundary
subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*current(t)/(sqrt(x(1)**2+x(2)**2))
    end select
end subroutine g_func

! exact solution of B
subroutine B_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),   dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t

    real(8) :: r,D,tmp1,tmp2,BB

    real(8), parameter :: sigma0 = 1d0
    real(8), parameter :: S0 = 1.2d0
    real(8), parameter :: b = 0.5d0
    real(8), parameter :: Rw = 1d0

    r = sqrt(x(1)**2+x(2)**2)
    D = 4*m_pi*sigma0
    tmp1 = b/t
    tmp2 = 1d0 + 0.5d0*sqrt(D/b)*log(Rw/r)

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        if (r<1d-8) then
            f(1) = 0d0
            f(2) = 0d0
        else
            BB = 4d0*S0/r*tmp1**2.5d0*tmp2*(tmp1*tmp2**2-1.5d0)*exp(-tmp1*tmp2**2)
            f(1) = -BB*x(2)/r
            f(2) = BB*x(1)/r
        end if
    end select
end subroutine B_func


! initial solution of B
subroutine B0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
        f(2) = 0d0
    end select
end subroutine B0_func

! initial solution of E
subroutine E0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine E0_func

! function of sigma
subroutine sigma_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type
    real(8) :: sigma0 = 1d0
    real(8) :: r

    if (.not. allocated(f)) allocate(f(1))
    r = sqrt(x(1)**2+x(2)**2)
    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sigma0/r**2
    end select
end subroutine sigma_func
    
end module problem_static


module problem_moving_exact
    use settings
    implicit none

    real(8) :: mu = 4*m_pi
    real(8) :: sigma0_ = 10d0
    real(8) :: Center(2) = [0d0,0d0]
    real(8) :: Rlist(1) = [1d0]
    integer :: M = 100
    integer :: Nlist(1) = [150]
    
contains

! right hand side function
subroutine f_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)
    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*t**2/mu - sigma0_*t*r**2 + sigma0_*(r**3*(1d0-r)*t**2)
    end select
end subroutine f_func

! velocity function
subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)
    if (.not. allocated(f)) allocate(f(2))


    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = -r*(1d0-r)*x(1)
        f(2) = -r*(1d0-r)*x(2)
    end select
end subroutine u_func
 
! magnetic field on boundary
subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = t**2*r
    end select
end subroutine g_func

! exact solution of B
subroutine B_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),   dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = -t**2*x(2)
        f(2) = t**2*x(1)
    end select
end subroutine B_func


! initial solution of B
subroutine B0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
        f(2) = 0d0
    end select
end subroutine B0_func

! initial solution of E
subroutine E0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine E0_func

! function of sigma
subroutine sigma_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sigma0_
    end select
end subroutine sigma_func
    
end module problem_moving_exact


module problem_moving
    use settings
    implicit none

    real(8) :: mu = 4*m_pi
    real(8) :: Center(2) = [0d0,0d0]
    real(8) :: Rlist(3) = [4d-1,6d-1,1d0]
    integer :: M = 40
    integer :: Nlist(3) = [16,40,60]
    
contains

! right hand side function
subroutine f_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine f_func

! velocity function
subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)
    if (.not. allocated(f)) allocate(f(2))


    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 10d0*(1d0-r)*(-x(1))*exp(-10d0*t)
        f(2) = 10d0*(1d0-r)*(-x(2))*exp(-10d0*t)
    end select
end subroutine u_func
 
! magnetic field on boundary
subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*(1d0+exp(-10d0*t))/r
    end select
end subroutine g_func


! initial solution of B
subroutine B0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
        f(2) = 0d0
    end select
end subroutine B0_func

! initial solution of E
subroutine E0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine E0_func

! function of sigma
subroutine sigma_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: r
    r = sqrt(x(1)**2+x(2)**2)

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        if(r<4d-1) then
            f(1) = 1d-6
        elseif(r<6d-1) then
            f(1) = 1d0
        else
            f(1) = 1d-6
        end if
    end select
end subroutine sigma_func
    
end module problem_moving


program test_Magnetic
    use problem_moving
    use solver_procs
    use petscvec
    use petscmat
    use petscksp
    use petscpc
    use settings
    use mesh, only: MeshInit, PrintMesh, AddBdryMarker
    use fe, only: fespaceInit, Interpolate
    use mesh_generator, only:QuadMeshOnSector
    use timer
    use fe_utils
    use visualize
    use solver
    use memory_usage
    use assembler
    use tools
    
    implicit none

    ! global variables
    real(8) :: t
    common /global/ t

    ! mesh
    type(mesh2D) :: Th
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! fe space for E ,B, sigma and u
    type(fespace) :: Eh,Bh,Qh,Uh
    integer, parameter :: DOF_type_Q = DOF_Q0
    integer, parameter :: DOF_type_E = DOF_Q1
    integer, parameter :: DOF_type_B = DOF_QuadRT1
    integer, parameter :: DOF_type_U = DOF_Q1
    
    ! parameters of Gauss quadrature
    integer, parameter :: Gauss_type = QuadPt9
    integer, parameter :: Gauss_type_bdry = LinePt3

    ! theta-scheme
    real(8) :: t_end = 1d0
    real(8) :: dt = 0.1d0
    real(8) :: t_old
    integer :: i_t,N_t
    real(8) :: theta = 0.5d0

    ! functions
    real(8), dimension(:), allocatable :: B,Bhat,E,sigma,tmp,u

    ! timer
    real :: t_test

    ! memory
    integer :: RSS

    ! petsc
    integer :: errpetsc
    PetscScalar, pointer :: vec_ptr(:)
    
    ! matrices and vectors
    integer, dimension(:,:), allocatable :: assemble_info
    Mat :: matM,matM1,matTotal
    Vec :: vec_B,vec_g,vec_f,vecTotal,vec_sol

    ! solver
    KSP :: ksp
    PC :: pc

    ! error
    real(8) :: L2_err

    ! generate mesh
    call timer_start()
    call QuadMeshOnSector(Center,Rlist,M,Nlist,elems,nodes)
    call MeshInit(elems,nodes,Th)
    call AddBdryMarker(Th, StraightBdry_func, 1)
    call AddBdryMarker(Th, CircleBdry_func, 2)
    call timer_end(t_test)
    write(*,*) "Mesh generation done. Time taken = ", t_test

    ! initialize fe space and petsc
    call timer_start()
    call fespaceInit(Qh, Th, DOF_type_Q, 1)
    call fespaceInit(Eh, Th, DOF_type_E, 1)
    call fespaceInit(Bh, Th, DOF_type_B, 2)
    call fespaceInit(Uh, Th, DOF_type_U, 2)
    call PetscInitialize("input/petsc_options.dat", errpetsc)
    call timer_end(t_test)
    write(*,*) "Initialize Done. Time taken = ", t_test

    ! initialzation of solutions
    call Interpolate(B, B0_func, Th, Bh)
    call Interpolate(E, E0_func, Th, Eh)
    allocate(Bhat(Bh%N_DOF))
    call CurlOperator(Th, E, tmp)
    Bhat = B - (1-theta)*dt*tmp
    call Interpolate(sigma, sigma_func, Th, Qh)
    call Interpolate(u, u_func, Th, Uh)

    dt = Th%hmax/2d0
    N_t = nint(t_end/dt)
    dt = t_end/N_t
    t = 0d0
    t_old = 0d0
    
    do i_t = 1,N_t

        t_old = t
        t = t + dt

        print *, "t = ", t
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! assemble matrices and vectors !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! assemble matrix M
        call CreateMat(matM, Eh%N_DOF, Eh%N_DOF, errpetsc)
        allocate(assemble_info(1,7))
        assemble_info(1,:) = [1, 1, DERIV_NONE, 1, DERIV_NONE, 1, DERIV_NONE]
        call AssembleMatrixElementFE(one_func, [1d0], Th, Eh, 0, Eh, 0, sigma, Qh,&
         assemble_info, Gauss_type, errpetsc, matM)
        deallocate(assemble_info)
        call AssembleMat(matM, errpetsc)

        ! assemble matrix M1
        call CreateMat(matM1, Eh%N_DOF, Eh%N_DOF, errpetsc)
        allocate(assemble_info(2,5))
        assemble_info(1,:) = [1, 1, DERIV_DX, 1, DERIV_DX]
        assemble_info(2,:) = [1, 1, DERIV_DY, 1, DERIV_DY]
        call AssembleMatrixElement(one_func, [1d0/mu,1d0/mu], Th, Eh, 0, Eh, 0,&
         assemble_info, Gauss_type, errpetsc, matM1)
        deallocate(assemble_info)
        call AssembleMat(matM1, errpetsc)

        ! assemble vector vec_B
        call CreateVec(vec_B, Eh%N_DOF, errpetsc)
        allocate(assemble_info(2,5))
        assemble_info(1,:) = [1, 1, DERIV_DY, 1, DERIV_NONE]
        assemble_info(2,:) = [1, 1, DERIV_DX, 2, DERIV_NONE]
        call AssembleVectorElementFE(one_func, [1d0/mu, -1d0/mu], Th, Eh, 0, Bhat, Bh, &
        assemble_info, Gauss_type, errpetsc, vec_B)
        deallocate(assemble_info)
        call AssembleVec(vec_B, errpetsc)

        ! assemble vector vec_g
        call CreateVec(vec_g, Eh%N_DOF, errpetsc)
        allocate(assemble_info(1,3))
        assemble_info(1,:) = [1, 1, DERIV_NONE]
        call AssembleVectorBdry(2, g_func, [1d0/mu], Th, Eh, 0, assemble_info, Gauss_type_bdry, &
        errpetsc, vec_g)
        deallocate(assemble_info)
        call AssembleVec(vec_g, errpetsc)

        ! assemble vector vec_f
        call CreateVec(vec_f, Eh%N_DOF, errpetsc)
        allocate(assemble_info(1,3))
        assemble_info(1,:) = [1, 1, DERIV_NONE]
        call AssembleVectorElement(f_func, [1d0], Th, Eh, 0, assemble_info, Gauss_type, &
        errpetsc, vec_f)
        deallocate(assemble_info)
        call AssembleVec(vec_f, errpetsc)

        ! form the linear system
        call CreateMat(matTotal, Eh%N_DOF, Eh%N_DOF, errpetsc)
        call MatCopy(matM, matTotal, DIFFERENT_NONZERO_PATTERN, errpetsc)
        ! call MatAXPY(matTotal, 1d0, matM, DIFFERENT_NONZERO_PATTERN, errpetsc)
        call MatAXPY(matTotal, theta*dt, matM1, DIFFERENT_NONZERO_PATTERN, errpetsc)
        call CreateVec(vecTotal, Eh%N_DOF, errpetsc)
        call VecAXPY(vecTotal, 1d0, vec_B, errpetsc)
        call VecAXPY(vecTotal, 1d0, vec_g, errpetsc)
        call VecAXPY(vecTotal, -1d0, vec_f, errpetsc)

        ! call MatView(matTotal, PETSC_VIEWER_STDOUT_WORLD, errpetsc)
        ! call VecView(vec_g, PETSC_VIEWER_STDOUT_WORLD, errpetsc)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! solve the linear system !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call CreateVec(vec_sol, Eh%N_DOF, errpetsc)
        call CreateSolver(ksp, matTotal, matTotal, errpetsc)
        call Solve(matTotal, vecTotal, ksp, pc, vec_sol, errpetsc)

        if(allocated(tmp)) deallocate(tmp)
        call VecGetArrayReadF90(vec_sol, vec_ptr, errpetsc)
        call CurlOperator(Th, vec_ptr, tmp)
        B = Bhat - theta*dt*tmp
        Bhat = B - (1d0-theta)*dt*tmp

        !!!!!!!!!!!!!!!!!!!!!
        !! update the mesh !!
        !!!!!!!!!!!!!!!!!!!!!
        call UpdateMeshAndSolution(u,dt,Th,B,Bhat)
        call Interpolate(u, u_func, Th, Uh)


        !!!!!!!!!!!!!!!!!
        !! free memory !!
        !!!!!!!!!!!!!!!!!

        call MatDestroy(matM, errpetsc)
        call MatDestroy(matM1, errpetsc)
        call MatDestroy(matTotal, errpetsc)
        call VecDestroy(vec_B, errpetsc)
        call VecDestroy(vec_g, errpetsc)
        call VecDestroy(vec_f, errpetsc)
        call VecDestroy(vecTotal, errpetsc)
        call VecDestroy(vec_sol, errpetsc)
        call KSPDestroy(ksp, errpetsc)
        call PCDestroy(pc, errpetsc)


        call system_mem_usage(RSS)
        write(*,*) "Memory usage = ", RSS

        
    end do


    !!!!!!!!!!!!!!!!!!!
    !! compute error !!
    !!!!!!!!!!!!!!!!!!!

    ! call ComputeError(B_func, B, Th, Bh, NORM_L2, Gauss_type, L2_err)
    ! write(*,*) "L2 error of B = ", L2_err

    !!!!!!!!!!!!!!!!!
    !! save result !!
    !!!!!!!!!!!!!!!!!

    call PlotFunction(B,Th, Bh, "output/B.vtk")

    call PetscFinalize(errpetsc)


end program test_Magnetic