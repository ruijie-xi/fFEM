module mymod
    use settings
    use quadrature
    use fe_utils
    use mesh
    implicit none
    
contains

    !(f(u), \nabla v)
    subroutine LocalVectorflux(i_elem, coe_fun, coe_fun_dim, Th, Vh_test, i_dim_test, deriv_type_test, &
        u, Vh_u, i_dim_u, deriv_type_u, flux_fun, flux_fun_dim, Gauss_type, localvec)
        integer, intent(in) :: i_elem, i_dim_test, deriv_type_test, Gauss_type, coe_fun_dim,i_dim_u,deriv_type_u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) ::  Vh_test, Vh_u
        procedure(func) :: coe_fun, flux_fun
        integer, intent(in) :: flux_fun_dim
        real(8), intent(out), dimension(:), allocatable :: localvec
        real(8), dimension(:), intent(in) :: u

        ! Gauss quadrature
        real(8), dimension(:,:), allocatable :: x,x_ref
        real(8), dimension(:), allocatable :: w,w_ref

        ! Basis functions
        real(8), dimension(:,:,:), allocatable :: basis_test

        ! fe function value
        real(8), dimension(:,:), allocatable :: fe_value
        real(8), dimension(:),allocatable :: flux_value

        ! coefficient function
        real(8), dimension(:), allocatable :: coe_value,tmp

        integer :: i_pt

        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, x, w)
        call getGaussRefElement(Gauss_type, x_ref, w_ref)
        
        call BasisLocal2D(x_ref, Th, Vh_test, i_elem, deriv_type_test, basis_test)

        call FEfunctionQuadValue(u, Th, Vh_u, i_elem, deriv_type_u, Gauss_type, fe_value)
        
        ! coefficient function
        allocate(coe_value(size(x,2)))
        do i_pt = 1, size(x, 2)
            call coe_fun(x(:,i_pt), tmp, DERIV_NONE)
            coe_value(i_pt) = tmp(coe_fun_dim)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*coe_value(i_pt)
            call flux_fun(fe_value(:,i_pt), flux_value, DERIV_NONE)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*flux_value(flux_fun_dim)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*w(i_pt)
        end do

        ! local vector
        allocate(localvec(Vh_test%N_local_basis))
        localvec = sum(basis_test(i_dim_test,:,:), 2)

    end subroutine LocalVectorflux

    ! (f(u)\cdot n, v) on \partial K
    subroutine LocalVectorElemBdry(i_elem, coe_fun, coe_fun_dim, Th, Vh_test, i_dim_test, deriv_type_test,&
         Gauss_type, flux, localvec)
        integer, intent(in) :: i_elem, i_dim_test, deriv_type_test, Gauss_type, coe_fun_dim
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) ::  Vh_test
        procedure(func) :: coe_fun
        real(8), intent(out), dimension(:), allocatable :: localvec
        real(8), dimension(:,:,:), intent(in) :: flux

        ! Basis functions
        real(8), dimension(:,:,:), allocatable :: basis_test

        ! coefficient function
        real(8), dimension(:), allocatable :: coe_value,tmp

        integer :: i_pt
        integer :: i_localedge

        ! Gauss quadrature
        real(8),dimension(:,:), allocatable :: vertices
        real(8), dimension(:,:), allocatable :: x_ref,x
        real(8), dimension(:), allocatable :: w_ref,w

        allocate(localvec(Vh_test%N_local_basis))
        localvec = 0d0

        do i_localedge = 1, Th%N_le
            ! get Gauss quadrature points and weights (x_ref and w)
            call getRefLinePts(Th%mesh_type, i_localedge, vertices)
            call getGaussQuadAnyLine(vertices, Gauss_type, x_ref, w_ref)
            deallocate(vertices)
            call getAnyLinePts(Th, i_elem, i_localedge, vertices)
            call getGaussQuadAnyLine(vertices, Gauss_type, x, w)
            
            ! Basis functions
            call BasisLocal2D(x_ref, Th, Vh_test, i_elem, deriv_type_test, basis_test)
            
            ! coefficient function
            allocate(coe_value(size(x,2)))
            do i_pt = 1, size(x, 2)
                call coe_fun(x(:,i_pt), tmp, DERIV_NONE)
                coe_value(i_pt) = tmp(coe_fun_dim)
                basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*coe_value(i_pt)
                basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*w(i_pt)
                basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*flux(i_pt,i_localedge,i_elem)
            end do

            ! local vector
            localvec = localvec + sum(basis_test(i_dim_test,:,:), 2)

            ! deallocate
            deallocate(coe_value)
            deallocate(basis_test)
            deallocate(x_ref)
            deallocate(x)
            deallocate(w_ref)
            deallocate(w)
            deallocate(vertices)
        end do

    end subroutine LocalVectorElemBdry

    ! calculate f(u)\cdot n on \partial K
    subroutine calculate_flux(u, Th, Vh, flux_fun, flux, Gauss_type_bdry)
        implicit none
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        procedure(func) :: flux_fun
        real(8), dimension(:,:,:), allocatable, intent(out) :: flux
        integer, intent(in) :: Gauss_type_bdry
    
        integer :: i_elem, i_le, i_elem_neigh, i_le_neigh, i_edge, i_pt
        integer :: num_pts
    
        real(8),dimension(:),allocatable :: w_ref
        real(8),dimension(:,:),allocatable :: x_ref
        real(8), dimension(:,:),allocatable :: fe_value, fe_value_neigh
        real(8),dimension(:), allocatable :: flux_value, flux_value_neigh
    
        real(8),dimension(2) :: normal
        real(8),parameter :: alpha = 1d0
    
        call getGaussQuadRefLine(Gauss_type_bdry, x_ref, w_ref)
        num_pts = size(w_ref)
    
        allocate(flux(num_pts, Th%N_le, Th%N_elem))
        flux = 0d0
    
        do i_elem = 1,Th%N_elem
            do i_le = 1,Th%N_le
    
                ! skip boundary
                i_edge = abs(Th%ElemEdgeConn(i_le,i_elem))
                if (Th%EdgeElemConn(1,i_edge)==0) cycle
    
                ! get normal
                call getEdgeNormal(Th, abs(Th%ElemEdgeConn(i_le,i_elem)), normal)
                if (Th%ElemEdgeConn(i_le,i_elem) < 0) normal = -normal
    
                ! find the neighboring element
                i_elem_neigh = abs(Th%EdgeElemConn(1,i_edge))
                i_le_neigh = Th%EdgeIdxInElem(1,i_edge)
                if (i_elem_neigh == i_elem) then
                    i_elem_neigh = abs(Th%EdgeElemConn(2,i_edge))
                    i_le_neigh = Th%EdgeIdxInElem(2,i_edge)
                end if
    
                call FEfunctionQuadValueLine(u, Th, Vh, i_elem, i_le, &
                DERIV_NONE, Gauss_type_bdry, fe_value)
                call FEfunctionQuadValueLine(u, Th, Vh, i_elem_neigh, i_le_neigh, &
                DERIV_NONE, Gauss_type_bdry, fe_value_neigh)
    
                do i_pt = 1,num_pts
                    call flux_fun(fe_value(:,i_pt), flux_value, DERIV_NONE)
                    call flux_fun(fe_value_neigh(:,i_pt), flux_value_neigh, DERIV_NONE)
                    flux(i_pt,i_le,i_elem) = 0.5d0*(flux_value(1)+flux_value_neigh(1))*normal(1) + &
                    0.5d0*(flux_value(2)+flux_value_neigh(2))*normal(2) &
                    - 0.5d0 * alpha * (fe_value_neigh(1,i_pt) - fe_value(1,i_pt))
                end do
            
            end do
        end do
    
    end subroutine calculate_flux
    
end module mymod


program test_burgers
    use settings
    use mesh, only: MeshInit, PrintMesh, AddBdryMarker
    use fe, only: fespaceInit, Interpolate
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    use mymod
    
    
    implicit none

    ! mesh
    type(mesh2D) :: Th
    integer :: Nx,Ny
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! fe space
    type(fespace) :: Vh
    integer, parameter :: mesh_type = MESH_TRIANGLE
    integer, parameter :: DOF_type = DOF_P0
    integer, parameter :: Gauss_type = TrianglePt4
    integer, parameter :: Gauss_type_bdry = LinePt2

    ! functions
    procedure(func) :: u0_func,one_func,flux_func,left_bdry_func,right_bdry_func
    real(8), dimension(:), allocatable :: u,u_old

    ! timer
    real :: t_test

    ! command line arguments
    character(len=100) :: arg

    ! flux value : (i_pt, i_le, i_elem)
    real(8),dimension(:,:,:),allocatable :: flux

    ! loop index
    integer :: i_elem,i_le,i_edge

    ! local matrices and vectors
    real(8), dimension(:,:), allocatable :: local_mat_tmp, local_mat
    real(8), dimension(:), allocatable :: local_vec_tmp, local_vec, local_sol

    ! time advance
    real(8) :: dt = 0.0025d0
    real(8) :: t = 0d0
    real(8) :: t_end = 0.8d0
    integer :: i_time = 0, n_time_steps
    n_time_steps = int(t_end/dt)

    ! get command line arguments
    call getarg(1, arg)
    read(arg,*) Nx
    call getarg(2, arg)
    read(arg,*) Ny

    ! generate mesh
    call timer_start()
    call TriangleMesh(-1d0,1d0,-1d0,1d0,Nx,Ny,elems,nodes)
    call MeshInit(elems,nodes,Th)
    call AddBdryMarker(Th, left_bdry_func, 1) ! left boundary
    call AddBdryMarker(Th, right_bdry_func, 2) ! right boundary
    call timer_end(t_test)
    write(*,*) "Mesh generation done. Time taken = ", t_test

    call fespaceInit(Vh, Th, DOF_type, 1)

    call Interpolate(u, u0_func, Th, Vh)
    allocate(u_old(size(u)))

    do i_time = 1, n_time_steps
        u_old = u
        t = t + dt
        write(*,*) "Time step ", i_time, " of ", n_time_steps, " at time ", t

        ! calculate flux
        call calculate_flux(u_old, Th, Vh, flux_func, flux, Gauss_type_bdry)
        
        do i_elem = 1,Th%N_elem
            allocate(local_mat(Vh%N_local_basis,Vh%N_local_basis))
            local_mat = 0d0
            call LocalMatrix(i_elem, one_func, 1, Th, Vh, Vh, 1, DERIV_NONE, 1, DERIV_NONE, Gauss_type, local_mat_tmp)
            local_mat = local_mat + local_mat_tmp/dt
            deallocate(local_mat_tmp)

            allocate(local_vec(Vh%N_local_basis))
            local_vec = 0d0
            call LocalVectorFE(i_elem, one_func, 1, Th, Vh, 1, DERIV_NONE, u_old, Vh, 1, DERIV_NONE, Gauss_type, local_vec_tmp)
            local_vec = local_vec + local_vec_tmp/dt
            deallocate(local_vec_tmp)

            call LocalVectorflux(i_elem, one_func, 1, Th, Vh, 1, DERIV_DX, u_old, Vh, 1, DERIV_NONE,&
             flux_func, 1, Gauss_type, local_vec_tmp)
            local_vec = local_vec + local_vec_tmp
            deallocate(local_vec_tmp)
            call LocalVectorflux(i_elem, one_func, 1, Th, Vh, 1, DERIV_DY, u_old, Vh, 1, DERIV_NONE,&
             flux_func, 2, Gauss_type, local_vec_tmp)
            local_vec = local_vec + local_vec_tmp

            call LocalVectorElemBdry(i_elem, one_func, 1, Th, Vh, 1, DERIV_NONE, Gauss_type_bdry, flux, local_vec_tmp)
            local_vec = local_vec - local_vec_tmp
            deallocate(local_vec_tmp)

            ! Dirichlet boundary condition
            do i_le = 1,Th%N_le
                i_edge = abs(Th%ElemEdgeConn(i_le,i_elem))
                if (Th%Edge2Bdry(i_edge)==0) cycle
                if (Th%BdryMarker(Th%Edge2Bdry(i_edge))==1) then
                    local_mat(:,1) = 0d0
                    local_mat(1,:) = 0d0
                    local_mat(1,1) = 1d0
                    local_vec(1) = -0.5d0
                end if
                if (Th%BdryMarker(Th%Edge2Bdry(i_edge))==2) then
                    local_mat(:,1) = 0d0
                    local_mat(1,:) = 0d0
                    local_mat(1,1) = 1d0
                    local_vec(1) = -0.5d0
                end if
            end do

            ! solve
            allocate(local_sol(Vh%N_local_basis))
            local_sol(1) = local_vec(1)/local_mat(1,1)
            u(Vh%ElemDOF(:,i_elem)) = local_sol

            ! deallocate
            deallocate(local_mat)
            deallocate(local_vec)
            deallocate(local_sol)
        end do


    end do

    ! write solution to file
    call PlotFunction(u, Th, Vh, "output/u.vtk")
    write(*,*) "Solution written to output/u.vtk"

end program test_burgers

subroutine u0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        if(x(1)<-0.25d0) then
            f(1) = -0.5d0
        elseif (x(1)>=0.25d0) then
            f(1) = -0.5d0
        else
            f(1) = 1d0
        end if
    end select
end subroutine u0_func

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

subroutine left_bdry_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = x(1) + 1d0
    end select
end subroutine left_bdry_func

subroutine right_bdry_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = x(1) - 1d0
    end select
end subroutine right_bdry_func

subroutine flux_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0.5d0*x(1)**2
        f(2) = 0d0
    end select
end subroutine flux_func