module fe_utils
    use settings
    use quadrature
    use fe, only: getLocalDof, BasisLocal2D
    use mesh, only: getAnyLinePts, getRefLinePts
    implicit none
    
contains

    ! Compute the integral of u over the domain
    subroutine ComputeIntegral(u, Th, Vh, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: Gauss_type
        real(8), intent(out), dimension(:), allocatable :: result

        integer :: i_elem
        real(8), dimension(:), allocatable :: val
        
        allocate(result(Vh%dim))
        allocate(val(Vh%dim))

        result = 0d0;

        do i_elem = 1, Th%N_elem
            call QuadIntegral(u, Th, Vh, i_elem, Gauss_type, val)
            result = result + val
        end do
    end subroutine ComputeIntegral

    ! Compute different norms of u
    subroutine ComputeNorm(u, Th, Vh, norm_type, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: norm_type, Gauss_type
        real(8), intent(out) :: result

        integer :: i_elem

        real(8) :: val

        select case (norm_type)
        case(NORM_L2)
            result = 0d0;
            do i_elem = 1, Th%N_elem
                call QuadNorm(u, Th, Vh, i_elem, DERIV_NONE, Gauss_type, val)
                result = result + val
            end do
            result = sqrt(result)

        case default
            print *, "Error: norm type not implemented"
            stop

        end select
    end subroutine ComputeNorm

    ! Compute different Errors of u and fun
    subroutine ComputeError(fun, u, Th, Vh, norm_type, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: norm_type, Gauss_type
        real(8), intent(out) :: result
        procedure(func) :: fun

        integer :: i_elem

        real(8) :: val

        select case (norm_type)
        case(NORM_L2)
            result = 0d0;
            do i_elem = 1, Th%N_elem
                call QuadError(fun, u, Th, Vh, i_elem, DERIV_NONE, Gauss_type, val)
                result = result + val
            end do
            result = sqrt(result)
        case(NORM_H1)
            result = 0d0;
            do i_elem = 1, Th%N_elem
                call QuadError(fun, u, Th, Vh, i_elem, DERIV_NONE, Gauss_type, val)
                result = result + val
                call QuadError(fun, u, Th, Vh, i_elem, DERIV_DX, Gauss_type, val)
                result = result + val
                call QuadError(fun, u, Th, Vh, i_elem, DERIV_DY, Gauss_type, val)
                result = result + val
            end do
            result = sqrt(result)

        case default
            print *, "Error: norm type not implemented"
            stop

        end select
    end subroutine ComputeError

    ! On the i_elem th element, compute \|u\|^2_{L^2} or \|dxu\|^2_{L^2}
    subroutine QuadNorm(u, Th, Vh, i_elem, deriv_type, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type, Gauss_type
        real(8), intent(out) :: result

        real(8), dimension(:,:), allocatable :: x
        real(8), dimension(:), allocatable :: w
        real(8), dimension(:,:),allocatable :: fe_value

        result = 0d0
        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, x, w)
        
        call FEfunctionQuadValue(u, Th, Vh, i_elem, deriv_type, Gauss_type, fe_value)

        result = sum(spread(w,1,Vh%dim)*fe_value**2)
        
    end subroutine QuadNorm

    subroutine QuadError(fun, u, Th, Vh, i_elem, deriv_type, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type, Gauss_type
        real(8), intent(out) :: result
        procedure(func) :: fun

        real(8), dimension(:,:), allocatable :: x
        real(8), dimension(:), allocatable :: w
        real(8), dimension(:,:),allocatable :: fe_value,exact_value,err_value
        real(8), dimension(:), allocatable :: tmp
        integer :: numpts,i_pt

        result = 0d0
        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, x, w)
        numpts = size(x, 2)
        allocate(exact_value(Vh%dim,numpts))
        allocate(err_value(Vh%dim,numpts))

        ! exact value
        do i_pt = 1, numpts
            call fun(x(:,i_pt), tmp, deriv_type)
            exact_value(:,i_pt) = tmp
        end do

        ! FE value
        call FEfunctionQuadValue(u, Th, Vh, i_elem, deriv_type, Gauss_type, fe_value)

        err_value = (fe_value - exact_value)**2
        result = sum(spread(w,1,Vh%dim)*err_value)
        
        
    end subroutine QuadError

    ! On the i_elem th element, compute \int u dx
    subroutine QuadIntegral(u, Th, Vh, i_elem, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, Gauss_type
        real(8), intent(out), dimension(:) :: result

        real(8), dimension(:,:), allocatable :: x
        real(8), dimension(:), allocatable :: w
        real(8), dimension(:,:),allocatable :: fe_value

        result = 0d0
        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, x, w)

        call FEfunctionQuadValue(u, Th, Vh, i_elem, DERIV_NONE, Gauss_type, fe_value)

        result = sum(spread(w,1,Vh%dim)*fe_value, 2)

    end subroutine QuadIntegral
    
    
    ! Compute the value of the FE function u at points (barycentric) of element i_elem
    function FEfunctionGetValue(u, Th, Vh, i_elem, x_ref, deriv_type) result(result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type
        real(8), dimension(:,:), allocatable :: result
        real(8), dimension(:,:), intent(in) :: x_ref

        real(8), dimension(:), allocatable :: local_u
        
        real(8), dimension(:,:,:), allocatable :: basis_values
        integer :: numpts,i_dim
        
        call BasisLocal2D(x_ref, Th, Vh, i_elem, deriv_type, basis_values)

        numpts = size(x_ref, 2)
        if(allocated(result)) deallocate(result)
        allocate(result(Vh%dim,numpts))
        do i_dim = 1,Vh%dim
            call getLocalDof(u, Vh, i_elem, i_dim, local_u)
            result(i_dim,:) = matmul(transpose(basis_values(i_dim,:,:)), local_u)
        end do 
    end function FEfunctionGetValue
    

    ! Compute the value of the FE function u at the quadrature points of element i_elem
    subroutine FEfunctionQuadValue(u, Th, Vh, i_elem, deriv_type, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type, Gauss_type
        real(8), dimension(:,:), intent(out), allocatable :: result

        real(8), dimension(:), allocatable :: local_u,w_ref
        real(8), dimension(:,:), allocatable :: x_ref
        real(8), dimension(:,:,:), allocatable :: basis_values
        integer :: numpts,i_dim
        
        
        call getGaussRefElement(Gauss_type, x_ref, w_ref)
        call BasisLocal2D(x_ref, Th, Vh, i_elem, deriv_type, basis_values)

        numpts = size(x_ref, 2)
        if(allocated(result)) deallocate(result)
        allocate(result(Vh%dim,numpts))
        do i_dim = 1,Vh%dim
            call getLocalDof(u, Vh, i_elem, i_dim, local_u)
            result(i_dim,:) = matmul(transpose(basis_values(i_dim,:,:)), local_u)
        end do 
    end subroutine FEfunctionQuadValue

    ! Compute the value of the FE function u at the quadrature points of i_le edge of element i_elem
    subroutine FEfunctionQuadValueLine(u, Th, Vh, i_elem, i_le, deriv_type, Gauss_type, result)
        real(8), dimension(:), intent(in) :: u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type, Gauss_type, i_le
        real(8), dimension(:,:), intent(out), allocatable :: result

        real(8), dimension(:), allocatable :: local_u,w
        real(8), dimension(:,:), allocatable :: x_ref,vertices
        real(8), dimension(:,:,:), allocatable :: basis_values
        integer :: numpts,i_dim
        
        call getRefLinePts(Th%mesh_type, i_le, vertices)
        call getGaussQuadAnyLine(vertices, Gauss_type, x_ref, w)
        call BasisLocal2D(x_ref, Th, Vh, i_elem, deriv_type, basis_values)

        numpts = size(x_ref, 2)
        allocate(result(Vh%dim,numpts))
        do i_dim = 1,Vh%dim
            call getLocalDof(u, Vh, i_elem, i_dim, local_u)
            result(i_dim,:) = matmul(transpose(basis_values(i_dim,:,:)), local_u)
        end do 
    end subroutine FEfunctionQuadValueLine


    ! Assemblers
    ! Local Matrix Assembler
    subroutine LocalMatrix(i_elem, coe_fun, coe_fun_dim, Th, Vh_trial, Vh_test, i_dim_trial, &
    deriv_type_trial, i_dim_test, deriv_type_test, Gauss_type, localmat)
        integer, intent(in) :: i_elem, i_dim_trial, deriv_type_trial, i_dim_test, deriv_type_test, Gauss_type, coe_fun_dim
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh_trial, Vh_test
        procedure(func) :: coe_fun
        real(8), intent(out), dimension(:,:), allocatable :: localmat

        ! Gauss quadrature
        real(8), dimension(:,:), allocatable :: x,x_ref
        real(8), dimension(:), allocatable :: w,w_ref

        ! Basis functions
        real(8), dimension(:,:,:), allocatable :: basis_trial, basis_test

        ! coefficient function
        real(8), dimension(:), allocatable :: coe_value,tmp

        integer :: i_pt

        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, x, w)
        call getGaussRefElement(Gauss_type, x_ref, w_ref)
        
        call BasisLocal2D(x_ref, Th, Vh_trial, i_elem, deriv_type_trial, basis_trial)
        call BasisLocal2D(x_ref, Th, Vh_test, i_elem, deriv_type_test, basis_test)
        
        ! coefficient function
        allocate(coe_value(size(x,2)))
        do i_pt = 1, size(x, 2)
            call coe_fun(x(:,i_pt), tmp, DERIV_NONE)
            coe_value(i_pt) = tmp(coe_fun_dim)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*coe_value(i_pt)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*w(i_pt)
        end do

        allocate(localmat(Vh_test%N_local_basis, Vh_trial%N_local_basis))
        localmat = matmul(basis_test(i_dim_test,:,:), transpose(basis_trial(i_dim_trial,:,:)))

    end subroutine LocalMatrix

    subroutine LocalMatrixFE(i_elem, coe_fun, coe_fun_dim, Th, Vh_trial, Vh_test, i_dim_trial, &
        deriv_type_trial, i_dim_test, deriv_type_test, u, Vh_u, i_dim_u, deriv_type_u, Gauss_type, localmat)
            integer, intent(in) :: i_elem, i_dim_trial, deriv_type_trial, i_dim_test, deriv_type_test, Gauss_type, coe_fun_dim
            type(mesh2D), intent(in) :: Th
            type(fespace), intent(in) :: Vh_trial, Vh_test
            procedure(func) :: coe_fun
            real(8), intent(out), dimension(:,:), allocatable :: localmat
            real(8), intent(in), dimension(:) :: u
            type(fespace) :: Vh_u
            integer :: i_dim_u
            integer :: deriv_type_u
    
            ! Gauss quadrature
            real(8), dimension(:,:), allocatable :: x,x_ref
            real(8), dimension(:), allocatable :: w,w_ref
    
            ! Basis functions
            real(8), dimension(:,:,:), allocatable :: basis_trial, basis_test
    
            ! coefficient function
            real(8), dimension(:), allocatable :: coe_value,tmp

            ! FE function value
            real(8), dimension(:,:), allocatable :: fe_value
    
            integer :: i_pt
    
            call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, x, w)
            call getGaussRefElement(Gauss_type, x_ref, w_ref)
            
            call BasisLocal2D(x_ref, Th, Vh_trial, i_elem, deriv_type_trial, basis_trial)
            call BasisLocal2D(x_ref, Th, Vh_test, i_elem, deriv_type_test, basis_test)

            call FEfunctionQuadValue(u, Th, Vh_u, i_elem, deriv_type_u, Gauss_type, fe_value)
            
            ! coefficient function
            allocate(coe_value(size(x,2)))
            do i_pt = 1, size(x, 2)
                call coe_fun(x(:,i_pt), tmp, DERIV_NONE)
                coe_value(i_pt) = tmp(coe_fun_dim)
                
                basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*coe_value(i_pt)*fe_value(i_dim_u,i_pt)
                basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*w(i_pt)
            end do
    
            allocate(localmat(Vh_test%N_local_basis, Vh_trial%N_local_basis))
            localmat = matmul(basis_test(i_dim_test,:,:), transpose(basis_trial(i_dim_trial,:,:)))
    
        end subroutine LocalMatrixFE

    ! Local Vector Assembler
    subroutine LocalVector(i_elem, coe_fun, coe_fun_dim, Th, Vh_test, i_dim_test, deriv_type_test, Gauss_type, localvec)
        integer, intent(in) :: i_elem, i_dim_test, deriv_type_test, Gauss_type, coe_fun_dim
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) ::  Vh_test
        procedure(func) :: coe_fun
        real(8), intent(out), dimension(:), allocatable :: localvec

        ! Gauss quadrature
        real(8), dimension(:,:), allocatable :: x,x_ref
        real(8), dimension(:), allocatable :: w,w_ref

        ! Basis functions
        real(8), dimension(:,:,:), allocatable :: basis_test

        ! coefficient function
        real(8), dimension(:), allocatable :: coe_value,tmp

        integer :: i_pt

        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, x, w)
        call getGaussRefElement(Gauss_type, x_ref, w_ref)
        
        call BasisLocal2D(x_ref, Th, Vh_test, i_elem, deriv_type_test, basis_test)
        
        ! coefficient function
        allocate(coe_value(size(x,2)))
        do i_pt = 1, size(x, 2)
            call coe_fun(x(:,i_pt), tmp, DERIV_NONE)
            coe_value(i_pt) = tmp(coe_fun_dim)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*coe_value(i_pt)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*w(i_pt)
        end do

        ! local vector
        allocate(localvec(Vh_test%N_local_basis))
        localvec = sum(basis_test(i_dim_test,:,:), 2)

    end subroutine LocalVector


    ! Local Vector Assembler with fe function
    subroutine LocalVectorFE(i_elem, coe_fun, coe_fun_dim, Th, Vh_test, i_dim_test, deriv_type_test, &
        u, Vh_u, i_dim_u, deriv_type_u, Gauss_type, localvec)
        integer, intent(in) :: i_elem, i_dim_test, deriv_type_test, Gauss_type, coe_fun_dim,i_dim_u,deriv_type_u
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) ::  Vh_test, Vh_u
        procedure(func) :: coe_fun
        real(8), intent(out), dimension(:), allocatable :: localvec
        real(8), dimension(:), intent(in) :: u

        ! Gauss quadrature
        real(8), dimension(:,:), allocatable :: x,x_ref
        real(8), dimension(:), allocatable :: w,w_ref

        ! Basis functions
        real(8), dimension(:,:,:), allocatable :: basis_test

        ! fe function value
        real(8), dimension(:,:), allocatable :: fe_value

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
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*fe_value(i_dim_u,i_pt)
            basis_test(:,:,i_pt) = basis_test(:,:,i_pt)*w(i_pt)
        end do

        ! local vector
        allocate(localvec(Vh_test%N_local_basis))
        localvec = sum(basis_test(i_dim_test,:,:), 2)

    end subroutine LocalVectorFE


    ! Local Vector Assembler (integral on line)
    subroutine LocalVectorLine(i_edge, i_elem, coe_fun, coe_fun_dim, Th, Vh_test, i_dim_test, deriv_type_test, Gauss_type, localvec)
        integer, intent(in) :: i_edge, i_elem, i_dim_test, deriv_type_test, Gauss_type, coe_fun_dim
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) ::  Vh_test
        procedure(func) :: coe_fun
        real(8), intent(out), dimension(:), allocatable :: localvec

        ! Basis functions
        real(8), dimension(:,:,:), allocatable :: basis_test

        ! coefficient function
        real(8), dimension(:), allocatable :: coe_value,tmp

        integer :: i_pt
        integer :: i_edge2elemlocal,i_localedge

        ! Gauss quadrature
        real(8),dimension(:,:), allocatable :: vertices
        real(8), dimension(:,:), allocatable :: x_ref,x
        real(8), dimension(:), allocatable :: w_ref,w
    
        ! Find the local edge index
        if(abs(Th%EdgeElemConn(1,i_edge)) == i_elem) then
            i_edge2elemlocal = 1
        elseif (abs(Th%EdgeElemConn(2,i_edge)) == i_elem) then
            i_edge2elemlocal = 2
        else
            print *, "Error: i_edge not on i_elem"
            stop
        end if
        i_localedge = Th%EdgeIdxInElem(i_edge2elemlocal,i_edge)
        
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
        end do

        ! local vector
        allocate(localvec(Vh_test%N_local_basis))
        localvec = sum(basis_test(i_dim_test,:,:), 2)

    end subroutine LocalVectorLine


end module fe_utils