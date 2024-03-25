module mymod
    use settings
    use quadrature
    use fe_utils
    use mesh
    implicit none
    
contains

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


    subroutine calculate_flux(u, Th, Vh, beta, flux, Gauss_type_bdry)
        implicit none
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        real(8), dimension(:,:,:), allocatable, intent(out) :: flux
        integer, intent(in) :: Gauss_type_bdry
        real(8),dimension(2),intent(in) :: beta
    
        integer :: i_elem, i_le, i_elem_tmp, i_le_tmp, i_edge
        integer :: num_pts
    
        real(8),dimension(:),allocatable :: w_ref
        real(8),dimension(:,:),allocatable :: x_ref
        real(8), dimension(:,:),allocatable :: fe_value
    
        real(8),dimension(2) :: normal
    
        call getGaussQuadRefLine(Gauss_type_bdry, x_ref, w_ref)
        num_pts = size(w_ref)
    
        allocate(flux(num_pts, Th%N_le, Th%N_elem))
        flux = 0d0
    
        do i_elem = 1,Th%N_elem
            do i_le = 1,Th%N_le
    
                ! skip boundary
                i_edge = abs(Th%ElemEdgeConn(i_le,i_elem))
                if (Th%EdgeElemConn(1,i_edge)==0) cycle
    
                call getEdgeNormal(Th, abs(Th%ElemEdgeConn(i_le,i_elem)), normal)
                if (Th%ElemEdgeConn(i_le,i_elem) < 0) normal = -normal
    
                if (beta(1)*normal(1) + beta(2)*normal(2) > 0d0) then
                    i_elem_tmp = i_elem
                    i_le_tmp = i_le
                else
                    i_elem_tmp = abs(Th%EdgeElemConn(1,i_edge))
                    i_le_tmp = Th%EdgeIdxInElem(1,i_edge)
                    if (i_elem_tmp == i_elem) then
                        i_elem_tmp = abs(Th%EdgeElemConn(2,i_edge))
                        i_le_tmp = Th%EdgeIdxInElem(2,i_edge)
                    end if
                end if
    
                call FEfunctionQuadValueLine(u, Th, Vh, i_elem_tmp, i_le_tmp, &
                DERIV_NONE, Gauss_type_bdry, fe_value)
    
                flux(:,i_le,i_elem) = (beta(1)*normal(1) + beta(2)*normal(2))*fe_value(1,:)
            
            end do
        end do
    
    end subroutine calculate_flux
    
end module mymod


program test_DG
    use settings
    use mesh, only: MeshInit, PrintMesh, AddBdryMarker
    use fe, only: fespaceInit, Interpolate
    use mesh_generator, only:TriangleMesh,SquareMesh
    use timer
    use fe_utils
    use visualize
    use memory_usage
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
    integer, parameter :: Gauss_type = TrianglePt9
    integer, parameter :: Gauss_type_bdry = LinePt3

    ! functions
    procedure(func) :: u0_func,one_func,beta_func
    real(8), dimension(:), allocatable :: u,u_old

    ! timer
    real :: t_test

    ! command line arguments
    character(len=100) :: arg

    ! velocity
    real(8) :: beta(2) = (/1d0, 0d0/)

    ! flux value : (i_pt, i_le, i_elem)
    real(8),dimension(:,:,:),allocatable :: flux

    ! loop index
    integer :: i_elem,i_le

    ! local matrices and vectors
    real(8), dimension(:,:), allocatable :: local_mat_tmp, local_mat
    real(8), dimension(:), allocatable :: local_vec_tmp, local_vec, local_sol

    ! time advance
    real(8) :: dt = 0.005d0
    real(8) :: t = 0d0
    real(8) :: t_end = 0.2d0
    integer :: i_time = 0, n_time_steps
    n_time_steps = int(t_end/dt)

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

    call fespaceInit(Vh, Th, DOF_type, 1)

    call Interpolate(u, u0_func, Th, Vh)
    allocate(u_old(size(u)))

    do i_time = 1, n_time_steps
        u_old = u
        t = t + dt
        write(*,*) "Time step ", i_time, " of ", n_time_steps, " at time ", t

        ! calculate flux
        call calculate_flux(u_old, Th, Vh, beta, flux, Gauss_type_bdry)
        
        do i_elem = 1,Th%N_elem
            allocate(local_mat(Vh%N_local_basis,Vh%N_local_basis))
            local_mat = 0d0
            call LocalMatrix(i_elem, one_func, 1, Th, Vh, Vh, 1, DERIV_NONE, 1, DERIV_NONE, Gauss_type, local_mat_tmp)
            local_mat = local_mat + local_mat_tmp/dt
            deallocate(local_mat_tmp)

            call LocalMatrix(i_elem, beta_func, 1, Th, Vh, Vh, 1, DERIV_NONE, 1, DERIV_DX, Gauss_type, local_mat_tmp)
            local_mat = local_mat - local_mat_tmp
            deallocate(local_mat_tmp)
            call LocalMatrix(i_elem, beta_func, 2, Th, Vh, Vh, 1, DERIV_NONE, 1, DERIV_DY, Gauss_type, local_mat_tmp)
            local_mat = local_mat - local_mat_tmp
            deallocate(local_mat_tmp)


            allocate(local_vec(Vh%N_local_basis))
            local_vec = 0d0
            call LocalVectorFE(i_elem, one_func, 1, Th, Vh, 1, DERIV_NONE, u_old, Vh, 1, DERIV_NONE, Gauss_type, local_vec_tmp)
            local_vec = local_vec + local_vec_tmp/dt
            deallocate(local_vec_tmp)

            call LocalVectorElemBdry(i_elem, one_func, 1, Th, Vh, 1, DERIV_NONE, Gauss_type_bdry, flux, local_vec_tmp)
            local_vec = local_vec - local_vec_tmp
            deallocate(local_vec_tmp)

            ! Dirichlet boundary condition
            do i_le = 1,Th%N_le
                if (Th%EdgeElemConn(1,abs(Th%ElemEdgeConn(i_le,i_elem)))==0) then
                    local_mat(:,1) = 0d0
                    local_mat(1,:) = 0d0
                    local_mat(1,1) = 1d0
                    local_vec(1) = 0d0
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

end program test_DG


subroutine u0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        if(x(1)>0.4d0 .and. x(1)<0.6d0 .and. x(2)>0.4d0 .and. x(2)<0.6d0) then
            f(1) = 1d0
        else
            f(1) = 0d0
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

subroutine beta_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 1d0
        f(2) = 0d0
    end select
end subroutine beta_func


