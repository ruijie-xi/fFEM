module fespace_P2
    use settings
    implicit none
    
contains

subroutine fespaceInit_P2(Vh,Th,dim)
    type(fespace), intent(out) :: Vh
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: dim

    integer :: i_dim

    Vh%N_local_basis = 6
    Vh%dim = dim
    Vh%basis_type = DOF_P2
    Vh%N_DOF = (Th%N_node + Th%N_edge) * Vh%dim
    allocate(Vh%ElemDOF(Vh%N_local_basis*Vh%dim, Th%N_elem))
    Vh%ElemDOF = 0
    do i_dim = 1, Vh%dim
        Vh%ElemDOF(6*(i_dim-1)+1:6*(i_dim-1)+3, :) = Th%ElemNodeConn + (i_dim-1)*(Th%N_node + Th%N_edge)
        Vh%ElemDOF(6*(i_dim-1)+4:6*i_dim, :) = abs(Th%ElemEdgeConn) + Th%N_node + (i_dim-1)*(Th%N_node + Th%N_edge)
    end do
end subroutine fespaceInit_P2

subroutine ComputeDof_P2(result,fun,Th,i_dof)
    type(mesh2D) :: Th
    procedure(func) :: fun
    real(8), intent(out) :: result
    integer,intent(in) :: i_dof

    integer :: idx, i_dim, i_node, i_edge
    integer :: i_node1, i_node2
    real(8), dimension(DIM__) :: coord
    real(8), dimension(:), allocatable :: val

    i_dim = i_dof / (Th%N_node+Th%N_edge) + 1
    idx = mod(i_dof, Th%N_node+Th%N_edge)
    if (idx == 0) then
        idx = Th%N_node+Th%N_edge
        i_dim = i_dim - 1
    end if

    if (idx<=Th%N_node) then
        i_node = idx
        call fun(Th%NodeCoord(:,i_node), val, DERIV_NONE)
        result = val(i_dim)
        return
    elseif (idx>Th%N_node) then
        i_edge = idx - Th%N_node
        i_node1 = Th%EdgeNodeConn(1, i_edge)
        i_node2 = Th%EdgeNodeConn(2, i_edge)
        coord = (Th%NodeCoord(:,i_node1) + Th%NodeCoord(:,i_node2)) / 2d0
        call fun(coord, val, DERIV_NONE)
        result = val(i_dim)
        return
    end if

end subroutine ComputeDof_P2

subroutine getEdgeDofIndex_P2(Th, Vh, i_edge, dof_index)
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_edge
    integer, dimension(:), intent(out), allocatable :: dof_index

    integer :: i_dim
    integer :: i_node1, i_node2

    i_node1 = Th%EdgeNodeConn(1, i_edge)
    i_node2 = Th%EdgeNodeConn(2, i_edge)
    allocate(dof_index(3*Vh%dim))
    do i_dim = 1, Vh%dim
        dof_index(3*(i_dim-1)+1) = (i_dim-1)*(Th%N_node+Th%N_edge) + i_node1
        dof_index(3*(i_dim-1)+2) = (i_dim-1)*(Th%N_node+Th%N_edge) + i_node2
        dof_index(3*(i_dim-1)+3) = Th%N_node + (i_dim-1)*(Th%N_node+Th%N_edge) + i_edge
    end do

end subroutine getEdgeDofIndex_P2


subroutine BasisReferenceP2(refpts, deriv_type, result)
    real(8), dimension(6,6), parameter :: A = &
    reshape((/ 1d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
                -3d0, -1d0, 0d0, 4d0, 0d0, 0d0,&
                -3d0, 0d0, -1d0, 0d0, 0d0, 4d0, &
                2d0, 2d0, 0d0, -4d0, 0d0, 0d0, &
                4d0, 0d0, 0d0, -4d0, 4d0, -4d0, &
                2d0, 0d0, 2d0, 0d0, 0d0, -4d0/),(/6,6/))
    real(8), intent(in), dimension(:,:) :: refpts
    integer, intent(in) :: deriv_type
    real(8), dimension(:,:), allocatable :: b
    real(8), intent(out), dimension(:,:) :: result
    
    allocate(b(6, size(refpts, 2)))

    select case (deriv_type)
    case (DERIV_NONE)
        b(1,:) = 1d0
        b(2,:) = refpts(1,:)
        b(3,:) = refpts(2,:)
        b(4,:) = refpts(1,:)**2
        b(5,:) = refpts(1,:)*refpts(2,:)
        b(6,:) = refpts(2,:)**2
    case (DERIV_DX)
        b(1,:) = 0d0
        b(2,:) = 1d0
        b(3,:) = 0d0
        b(4,:) = 2d0*refpts(1,:)
        b(5,:) = refpts(2,:)
        b(6,:) = 0d0
    case (DERIV_DY)
        b(1,:) = 0d0
        b(2,:) = 0d0
        b(3,:) = 1d0
        b(4,:) = 0d0
        b(5,:) = refpts(1,:)
        b(6,:) = 2d0*refpts(2,:)
    case default
        print *, 'BasisP2: Unknown derivative type'
        stop
    end select

    result = matmul(A, b)
end subroutine BasisReferenceP2

subroutine BasisLocalP2(refpts, Th, Vh, i_elem, deriv_type, result)
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
        call BasisReferenceP2(refpts, deriv_type, result(1,:,:))
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
        call BasisReferenceP2(refpts, DERIV_DX, basis_dxh)
        call BasisReferenceP2(refpts, DERIV_DY, basis_dyh)

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
        call BasisReferenceP2(refpts, DERIV_DX, basis_dxh)
        call BasisReferenceP2(refpts, DERIV_DY, basis_dyh)

        result(1,:,:) = ((x1-x3)*basis_dxh + (x2-x1)*basis_dyh) / detJ

    case default
        print *, 'BasisLocalP2: Unknown derivative type'
        stop

    end select

    do i_dim = 1, Vh%dim
        result(i_dim,:,:) = result(1,:,:);
    end do
end subroutine BasisLocalP2
    
end module fespace_P2