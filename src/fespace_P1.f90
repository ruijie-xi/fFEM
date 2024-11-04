module fespace_P1
    use settings
    implicit none
    
contains

subroutine fespaceInit_P1(Vh,Th,dim)
    type(fespace), intent(out) :: Vh
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: dim

    integer :: i_dim

    Vh%N_local_basis = 3
    Vh%dim = dim
    Vh%basis_type = DOF_P1
    Vh%N_DOF = Th%N_node * Vh%dim
    Vh%isStack = 1
    allocate(Vh%ElemDOF(Vh%N_local_basis*Vh%dim, Th%N_elem))
    Vh%ElemDOF = 0
    do i_dim = 1, Vh%dim
        Vh%ElemDOF(3*(i_dim-1)+1:3*i_dim, :) = Th%ElemNodeConn + (i_dim-1)*Th%N_node
    end do
end subroutine fespaceInit_P1


! Given i_dof, find the coordinates of the corresponding node and the dimension
! Lagrangian P1 element
subroutine FindDofLocation_P1(Th, i_dof, pt, i_dim)
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: i_dof
    real(8), dimension(DIM__), intent(out) :: pt
    integer, intent(out) :: i_dim
    
    integer ::  i_node
    
    i_dim = i_dof / Th%N_node + 1
    i_node = mod(i_dof, Th%N_node)
    if (i_node == 0) then
        i_node = Th%N_node
        i_dim = i_dim - 1
    end if
    
    pt = Th%NodeCoord(:,i_node)
    
    
end subroutine

subroutine ComputeDof_P1(result,fun,Th,i_dof)
    type(mesh2D) :: Th
    procedure(func) :: fun
    real(8), intent(out) :: result
    integer,intent(in) :: i_dof
    
    integer :: i_dim
    real(8), dimension(DIM__) :: pt 
    real(8), dimension(:), allocatable :: val
    
    call FindDofLocation_P1(Th, i_dof, pt, i_dim)
    call fun(pt, val, DERIV_NONE)
    result = val(i_dim)

end subroutine ComputeDof_P1

subroutine getEdgeDofIndex_P1(Th, Vh, i_edge, dof_index)
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_edge
    integer, dimension(:), intent(out), allocatable :: dof_index

    integer :: i_dim
    integer :: i_node1, i_node2

    i_node1 = Th%EdgeNodeConn(1, i_edge)
    i_node2 = Th%EdgeNodeConn(2, i_edge)
    allocate(dof_index(2*Vh%dim))
    do i_dim = 1, Vh%dim
        dof_index(2*(i_dim-1)+1) = (i_dim-1)*Th%N_node + i_node1
        dof_index(2*(i_dim-1)+2) = (i_dim-1)*Th%N_node + i_node2
    end do

end subroutine getEdgeDofIndex_P1


subroutine BasisReferenceP1(refpts, deriv_type, result)
    real(8), dimension(3,3), parameter :: A = reshape((/ 1d0, 0d0, 0d0, -1d0, 1d0, 0d0, -1d0, 0d0, 1d0 /), (/3,3/))
    real(8), intent(in), dimension(:,:) :: refpts
    integer, intent(in) :: deriv_type
    real(8), dimension(:,:), allocatable :: b
    real(8), intent(out), dimension(:,:) :: result
    
    allocate(b(3, size(refpts, 2)))

    select case (deriv_type)
    case (DERIV_NONE)
        b(1,:) = 1d0
        b(2,:) = refpts(1,:)
        b(3,:) = refpts(2,:)
    case (DERIV_DX)
        b(1,:) = 0d0
        b(2,:) = 1d0
        b(3,:) = 0d0
    case (DERIV_DY)
        b(1,:) = 0d0
        b(2,:) = 0d0
        b(3,:) = 1d0
    case default
        print *, 'BasisP1: Unknown derivative type'
        stop
    end select

    result = matmul(A, b)
end subroutine BasisReferenceP1

subroutine BasisLocalP1(refpts, Th, Vh, i_elem, deriv_type, result)
    real(8), intent(in), dimension(:,:) :: refpts
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_elem, deriv_type
    real(8), intent(out), dimension(:,:,:), allocatable :: result

    real(8), dimension(:,:), allocatable :: basis_dxh, basis_dyh

    integer :: num_pts, i_dim

    real(8) :: x1,x2,x3,y1,y2,y3, detJ

    num_pts = size(refpts, 2)
    allocate(result(Vh%dim, Vh%N_local_basis, num_pts))

    select case (deriv_type)
    case (DERIV_NONE)
        call BasisReferenceP1(refpts, deriv_type, result(1,:,:))
    case (DERIV_DX)

        ! get coordinates and jacobian
        x1 = Th%NodeCoord(1, Th%ElemNodeConn(1, i_elem))
        x2 = Th%NodeCoord(1, Th%ElemNodeConn(2, i_elem))
        x3 = Th%NodeCoord(1, Th%ElemNodeConn(3, i_elem))
        y1 = Th%NodeCoord(2, Th%ElemNodeConn(1, i_elem))
        y2 = Th%NodeCoord(2, Th%ElemNodeConn(2, i_elem))
        y3 = Th%NodeCoord(2, Th%ElemNodeConn(3, i_elem))

        detJ = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

        allocate(basis_dxh(Vh%N_local_basis, num_pts))
        allocate(basis_dyh(Vh%N_local_basis, num_pts))
        call BasisReferenceP1(refpts, DERIV_DX, basis_dxh)
        call BasisReferenceP1(refpts, DERIV_DY, basis_dyh)

        result(1,:,:) = ((y3-y1)*basis_dxh + (y1-y2)*basis_dyh) / detJ
    
    case (DERIV_DY)
        ! get coordinates and jacobian
        x1 = Th%NodeCoord(1, Th%ElemNodeConn(1, i_elem))
        x2 = Th%NodeCoord(1, Th%ElemNodeConn(2, i_elem))
        x3 = Th%NodeCoord(1, Th%ElemNodeConn(3, i_elem))
        y1 = Th%NodeCoord(2, Th%ElemNodeConn(1, i_elem))
        y2 = Th%NodeCoord(2, Th%ElemNodeConn(2, i_elem))
        y3 = Th%NodeCoord(2, Th%ElemNodeConn(3, i_elem))
        detJ = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

        allocate(basis_dxh(Vh%N_local_basis, num_pts))
        allocate(basis_dyh(Vh%N_local_basis, num_pts))
        call BasisReferenceP1(refpts, DERIV_DX, basis_dxh)
        call BasisReferenceP1(refpts, DERIV_DY, basis_dyh)

        result(1,:,:) = ((x1-x3)*basis_dxh + (x2-x1)*basis_dyh) / detJ

    case default
        print *, 'BasisLocalP1: Unknown derivative type'
        stop

    end select

    do i_dim = 1, Vh%dim
        result(i_dim,:,:) = result(1,:,:);
    end do
end subroutine BasisLocalP1

end module fespace_P1