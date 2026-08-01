module fespace_QuadRT1
    use settings
    use tools, only: assert
    use mesh, only: getEdgeTangent,getEdgeLength,getEdgeNormal
    use geometry, only: JacobianCoeffQuad,  ref_pts_edge_to_elem, ref_edge_normal
    use quadrature
    implicit none
    
contains

subroutine fespaceInit_QuadRT1(Vh,Th,dim)
    type(fespace), intent(out) :: Vh
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: dim

    Vh%N_local_basis = 4
    call assert(dim.eq.2, "fespaceInit_QuadRT1: dim must be 2")
    Vh%dim = dim
    Vh%basis_type = DOF_QuadRT1
    Vh%N_DOF = Th%N_edge
    Vh%isStack = 0
    allocate(Vh%ElemDOF(Vh%N_local_basis, Th%N_elem))
    Vh%ElemDOF = abs(Th%ElemEdgeConn)
    
end subroutine fespaceInit_QuadRT1

subroutine ComputeDof_QuadRT1(result,fun,Th,i_dof)
    type(mesh2D) :: Th
    procedure(func) :: fun
    real(8), intent(out) :: result
    integer,intent(in) :: i_dof

    integer :: i_edge,i_edge_in_elem, i_elem, i_pt
    integer :: orientation
    
    integer, parameter :: Gauss_type_1d = LinePt4
    integer, parameter :: Gauss_type_2d = QuadPt16
    real(8), dimension(:), allocatable :: w_line_h
    real(8), dimension(:,:), allocatable :: pts_line_h
    real(8), dimension(:), allocatable :: w_line 
    real(8), dimension(:,:), allocatable :: x_line
    real(8), dimension(:,:), allocatable :: pts_elem_h
    real(8), dimension(:), allocatable :: J11,J12,J21,J22
    real(8), dimension(:), allocatable :: normal_hat
    real(8), dimension(:,:), allocatable :: normal
    real(8) :: length
    
    real(8), dimension(:,:), allocatable :: val
    real(8), dimension(:), allocatable :: val_tmp
        
    i_edge = i_dof
    i_elem = abs(Th%EdgeElemConn(2,i_edge))
    i_edge_in_elem = Th%EdgeIdxInElem(2,i_edge)
    
    orientation = 1
    if(Th%ElemEdgeConn(i_edge_in_elem,i_elem)<0) orientation = -1
    
    call ref_edge_normal(Th%mesh_type, i_edge_in_elem, orientation, normal_hat)
    
    call getGaussQuadRefLine(Gauss_type_1d, pts_line_h, w_line_h) 
    pts_line_h = pts_line_h*0.5d0 + 0.5d0
    
    call ref_pts_edge_to_elem(Th%mesh_type, i_edge_in_elem, pts_line_h, orientation, pts_elem_h)
    
    call JacobianCoeffQuad(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), pts_elem_h, J11,J12,J21,J22)
    
    allocate(normal(2,size(pts_elem_h,2)))
    normal(1,:) = J22*normal_hat(1) - J21*normal_hat(2)
    normal(2,:) = J11*normal_hat(2) - J12*normal_hat(1)
    
    call getEdgeLength(Th, i_edge, length)
    
    normal = normal/length 
    
    call getGaussQuadAnyLine(Th%NodeCoord(:,Th%EdgeNodeConn(:,i_edge)), Gauss_type_1d, x_line, w_line)
    
    allocate(val(2,size(x_line,2)))
    do i_pt = 1, size(x_line,2)
        call fun(x_line(:,i_pt), val_tmp, DERIV_NONE)
        val(:,i_pt) = val_tmp
    end do
    
    result = sum((val(1,:)*normal(1,:) + val(2,:)*normal(2,:))*w_line)
    
end subroutine ComputeDof_QuadRT1

subroutine getEdgeDofIndex_QuadRT1(Th, Vh, i_edge, dof_index)
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_edge
    integer, dimension(:), intent(out), allocatable :: dof_index

    allocate(dof_index(1))
    dof_index(1) = i_edge

end subroutine getEdgeDofIndex_QuadRT1


subroutine BasisReferenceQuadRT1(refpts, deriv_type, result)
    real(8), dimension(4,3):: A1,A2
    real(8), intent(in), dimension(:,:) :: refpts
    integer, intent(in) :: deriv_type
    real(8), dimension(:,:), allocatable :: b
    real(8), intent(out), dimension(:,:,:) :: result

    A1(:, 1) = (/ 0d0, 0d0, 0d0, -1d0 /)
    A1(:, 2) = (/ 0d0, 1d0, 0d0, 1d0 /)
    A1(:, 3) = (/ 0d0, 0d0, 0d0, 0d0 /)

    A2(:, 1) = (/ -1d0, 0d0, 0d0, 0d0 /)
    A2(:, 2) = (/ 0d0, 0d0, 0d0, 0d0 /)
    A2(:, 3) = (/ 1d0, 0d0, 1d0, 0d0 /)
    
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
        print *, 'BasisQuadRT1: Unknown derivative type'
        stop
    end select

    result(1,:,:) = matmul(A1, b)
    result(2,:,:) = matmul(A2, b)

end subroutine BasisReferenceQuadRT1

subroutine BasisLocalQuadRT1(refpts, Th, Vh, i_elem, deriv_type, result)
    real(8), intent(in), dimension(:,:) :: refpts
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_elem, deriv_type
    real(8), intent(out), dimension(:,:,:), allocatable :: result

    real(8), dimension(:), allocatable :: J11,J12,J21,J22,J
    integer :: num_pts, num_basis, i

    num_basis = Vh%N_local_basis
    num_pts = size(refpts, 2) ! number of points in reference element
    
    ! Jacobian coefficients
    call JacobianCoeffQuad(Th%NodeCoord(:, Th%ElemNodeConn(:, i_elem)), refpts, J11,J12,J21,J22)
    allocate(J(size(refpts,2)))
    J = J11*J22 - J12*J21

    select case (deriv_type)
    case (DERIV_NONE)
        
        block
            real(8), dimension(:,:,:),allocatable :: basis_ref
            
            allocate(result(Vh%dim, Vh%N_local_basis, num_pts))
        
            allocate(basis_ref(Vh%dim, Vh%N_local_basis, num_pts))
            call BasisReferenceQuadRT1(refpts, DERIV_NONE, basis_ref)
            
            result(1,:,:) = spread(J11,1,num_basis)*basis_ref(1,:,:)/spread(J,1,num_basis) &
            + spread(J12,1,num_basis)*basis_ref(2,:,:)/spread(J,1,num_basis)
            result(2,:,:) = spread(J21,1,num_basis)*basis_ref(1,:,:)/spread(J,1,num_basis) &
            + spread(J22,1,num_basis)*basis_ref(2,:,:)/spread(J,1,num_basis)
        
        end block
        
    case (DERIV_DIV)
        
        block
            real(8), dimension(:,:,:),allocatable :: basis_ref_dxh,basis_ref_dyh,basis_ref_divh
            
            allocate(result(1,Vh%N_local_basis,num_pts))
            allocate(basis_ref_divh(1,Vh%N_local_basis,num_pts))
            allocate(basis_ref_dxh(2,Vh%N_local_basis,num_pts))
            allocate(basis_ref_dyh(2,Vh%N_local_basis,num_pts))
            
            call BasisReferenceQuadRT1(refpts, DERIV_DX, basis_ref_dxh)
            call BasisReferenceQuadRT1(refpts, DERIV_DY, basis_ref_dyh)
            
            basis_ref_divh(1,:,:) = basis_ref_dxh(1,:,:) + basis_ref_dyh(2,:,:)
            
            result(1,:,:) = basis_ref_divh(1,:,:)/spread(J,1,num_basis)
            
        end block
        
    case default
        print *, 'BasisLocalQuadRT1: Unknown derivative type'
        stop

    end select
    
    ! Transform to physical element to make the orientation match the mesh
    do i = 1,4
        if(Th%ElemEdgeConn(i,i_elem)<0) then
            result(:,i,:) = -result(:,i,:)
        end if
    end do
    
end subroutine BasisLocalQuadRT1

end module fespace_QuadRT1