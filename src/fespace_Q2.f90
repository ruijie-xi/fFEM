module fespace_Q2
    use settings
    use geometry
    implicit none
    
contains

subroutine fespaceInit_Q2(Vh,Th,dim)
    type(fespace), intent(out) :: Vh
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: dim

    integer :: i_dim, i

    Vh%N_local_basis = 9
    Vh%dim = dim
    Vh%basis_type = DOF_Q2
    Vh%N_DOF = Vh%dim * (Th%N_node + Th%N_edge + Th%N_elem)
    Vh%isStack = 1
    allocate(Vh%ElemDOF(Vh%N_local_basis*Vh%dim, Th%N_elem))
    Vh%ElemDOF = 0
    do i_dim = 1, Vh%dim
        Vh%ElemDOF(9*(i_dim-1)+1:9*(i_dim-1)+4, :) = &
        Th%ElemNodeConn + (i_dim-1)*(Th%N_node + Th%N_edge + Th%N_elem)
        Vh%ElemDOF(9*(i_dim-1)+5:9*(i_dim-1)+8, :) = &
        abs(Th%ElemEdgeConn) + Th%N_node + (i_dim-1)*(Th%N_node + Th%N_edge + Th%N_elem)
        Vh%ElemDOF(9*i_dim, :) = Th%N_node + Th%N_edge + &
         (i_dim-1)*(Th%N_node + Th%N_edge + Th%N_elem) + (/(i,i=1,Th%N_elem)/)
    end do
    
end subroutine fespaceInit_Q2

subroutine ComputeDof_Q2(result,fun,Th,i_dof)
    type(mesh2D) :: Th
    procedure(func) :: fun
    real(8), intent(out) :: result
    integer,intent(in) :: i_dof

    integer :: i_dof0, i_dim
    integer :: i_node, i_edge, i_elem
    real(8), dimension(:), allocatable :: val

    i_dim = i_dof / (Th%N_node + Th%N_edge + Th%N_elem) + 1
    i_dof0 = mod(i_dof, Th%N_node + Th%N_edge + Th%N_elem)
    if (i_dof0 == 0) then
        i_dof0 = Th%N_node + Th%N_edge + Th%N_elem
        i_dim = i_dim - 1
    end if
    
    ! find the location of idof0
    if (i_dof0 <= Th%N_node) then
        i_node = i_dof0
        call fun(Th%NodeCoord(:,i_node), val, DERIV_NONE)
        result = val(i_dim)
        return
    else if (i_dof0 <= Th%N_node + Th%N_edge) then
        i_edge = i_dof0 - Th%N_node
        block 
            real(8), dimension(2) :: x1, x2, x 
            x1 = Th%NodeCoord(:,Th%EdgeNodeConn(1,i_edge))
            x2 = Th%NodeCoord(:,Th%EdgeNodeConn(2,i_edge))
            x = (x1 + x2) / 2d0
            call fun(x, val, DERIV_NONE)
            result = val(i_dim)
        end block
    else
        i_elem = i_dof0 - Th%N_node - Th%N_edge
        block 
            real(8), dimension(2) :: x1, x2, x3, x4, x
            x1 = Th%NodeCoord(:,Th%ElemNodeConn(1,i_elem))
            x2 = Th%NodeCoord(:,Th%ElemNodeConn(2,i_elem))
            x3 = Th%NodeCoord(:,Th%ElemNodeConn(3,i_elem))
            x4 = Th%NodeCoord(:,Th%ElemNodeConn(4,i_elem))
            x = (x1 + x2 + x3 + x4) / 4d0
            call fun(x, val, DERIV_NONE)
            result = val(i_dim)
        end block
    end if

end subroutine ComputeDof_Q2

subroutine getEdgeDofIndex_Q2(Th, Vh, i_edge, dof_index)
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
        dof_index(3*(i_dim-1)+1) = (i_dim-1)*(Th%N_node + Th%N_edge + Th%N_elem) + i_node1
        dof_index(3*(i_dim-1)+2) = (i_dim-1)*(Th%N_node + Th%N_edge + Th%N_elem) + i_node2
        dof_index(3*i_dim) = Th%N_node + i_edge + (i_dim-1)*(Th%N_node + Th%N_edge + Th%N_elem)
    end do

end subroutine getEdgeDofIndex_Q2


subroutine BasisReferenceQ2(refpts, deriv_type, result)
    real(8), dimension(9,9):: A
    real(8), intent(in), dimension(:,:) :: refpts
    integer, intent(in) :: deriv_type
    real(8), dimension(:,:), allocatable :: b
    real(8), intent(out), dimension(:,:) :: result
    A(:, 1) = (/ 1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
    A(:, 2) = (/ -3d0, -1d0, 0d0, 0d0, 4d0, 0d0, 0d0, 0d0, 0d0 /)
    A(:, 3) = (/ -3d0, 0d0, 0d0, -1d0, 0d0, 0d0, 0d0, 4d0, 0d0 /)
    A(:, 4) = (/ 2d0, 2d0, 0d0, 0d0, -4d0, 0d0, 0d0, 0d0, 0d0 /)
    A(:, 5) = (/ 9d0, 3d0, 1d0, 3d0, -12d0, -4d0, -4d0, -12d0, 16d0 /)
    A(:, 6) = (/ 2d0, 0d0, 0d0, 2d0, 0d0, 0d0, 0d0, -4d0, 0d0 /)
    A(:, 7) = (/ -6d0, -6d0, -2d0, -2d0, 12d0, 8d0, 4d0, 8d0, -16d0 /)
    A(:, 8) = (/ -6d0, -2d0, -2d0, -6d0, 8d0, 4d0, 8d0, 12d0, -16d0 /)
    A(:, 9) = (/ 4d0, 4d0, 4d0, 4d0, -8d0, -8d0, -8d0, -8d0, 16d0 /)
    
    allocate(b(9, size(refpts, 2)))

    select case (deriv_type)
    case (DERIV_NONE)
        b(1,:) = 1d0
        b(2,:) = refpts(1,:)
        b(3,:) = refpts(2,:)
        b(4,:) = refpts(1,:)**2
        b(5,:) = refpts(1,:)*refpts(2,:)
        b(6,:) = refpts(2,:)**2
        b(7,:) = refpts(1,:)**2*refpts(2,:)
        b(8,:) = refpts(1,:)*refpts(2,:)**2
        b(9,:) = refpts(1,:)**2*refpts(2,:)**2
    case (DERIV_DX)
        b(1,:) = 0d0
        b(2,:) = 1d0
        b(3,:) = 0d0
        b(4,:) = 2d0*refpts(1,:)
        b(5,:) = refpts(2,:)
        b(6,:) = 0d0
        b(7,:) = 2d0*refpts(1,:)*refpts(2,:)
        b(8,:) = refpts(2,:)**2
        b(9,:) = 2d0*refpts(1,:)*refpts(2,:)**2
    case (DERIV_DY)
        b(1,:) = 0d0
        b(2,:) = 0d0
        b(3,:) = 1d0
        b(4,:) = 0d0
        b(5,:) = refpts(1,:)
        b(6,:) = 2d0*refpts(2,:)
        b(7,:) = refpts(1,:)**2
        b(8,:) = 2d0*refpts(1,:)*refpts(2,:)
        b(9,:) = 2d0*refpts(1,:)**2*refpts(2,:)
    case default
        print *, 'BasisQ2: Unknown derivative type'
        stop
    end select

    result = matmul(A, b)
end subroutine BasisReferenceQ2

subroutine BasisLocalQ2(refpts, Th, Vh, i_elem, deriv_type, result)
    real(8), intent(in), dimension(:,:) :: refpts
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_elem, deriv_type
    real(8), intent(out), dimension(:,:,:), allocatable :: result

    real(8), dimension(2,4) :: vertices
    real(8), dimension(:,:), allocatable :: basis_dxh, basis_dyh
    real(8), dimension(:), allocatable :: J11,J12,J21,J22,J

    integer :: num_pts, i_dim

    num_pts = size(refpts, 2)
    allocate(result(Vh%dim, Vh%N_local_basis, num_pts))

    vertices = Th%NodeCoord(:, Th%ElemNodeConn(:, i_elem))
    call JacobianCoeffQuad(vertices, refpts, J11,J12,J21,J22)
    allocate(J(size(refpts,2)))
    J = J11*J22 - J12*J21

    select case (deriv_type)
    case (DERIV_NONE)
        call BasisReferenceQ2(refpts, deriv_type, result(1,:,:))
    case (DERIV_DX)
        allocate(basis_dxh(Vh%N_local_basis, num_pts))
        allocate(basis_dyh(Vh%N_local_basis, num_pts))
        call BasisReferenceQ2(refpts, DERIV_DX, basis_dxh)
        call BasisReferenceQ2(refpts, DERIV_DY, basis_dyh)
        result(1,:,:) = (spread(J22,1,size(basis_dxh,1))*basis_dxh &
                        - spread(J21,1,size(basis_dyh,1))*basis_dyh) &
                        /spread(J,1,size(basis_dxh,1))
    case (DERIV_DY)
        allocate(basis_dxh(Vh%N_local_basis, num_pts))
        allocate(basis_dyh(Vh%N_local_basis, num_pts))
        call BasisReferenceQ2(refpts, DERIV_DX, basis_dxh)
        call BasisReferenceQ2(refpts, DERIV_DY, basis_dyh)
        result(1,:,:) = (-spread(J12,1,size(basis_dxh,1))*basis_dxh &
                        + spread(J11,1,size(basis_dyh,1))*basis_dyh)&
                        /spread(J,1,size(basis_dxh,1))
    case default
        print *, 'BasisLocalQ2: Unknown derivative type'
        stop

    end select

    do i_dim = 1, Vh%dim
        result(i_dim,:,:) = result(1,:,:);
    end do
end subroutine BasisLocalQ2
    
end module fespace_Q2