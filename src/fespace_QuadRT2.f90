module fespace_QuadRT2
    use settings
    use tools, only: assert
    use mesh, only: getEdgeTangent,getEdgeLength,getEdgeNormal
    use geometry, only: JacobianCoeffQuad,  ref_pts_edge_to_elem, ref_edge_normal
    use quadrature
    implicit none
    private

    integer,parameter :: nDOFonNode = 0
    integer,parameter :: nDOFonEdge = 2
    integer,parameter :: nDOFonElem = 4
    
    public :: fespaceInit_QuadRT2
    public :: ComputeDof_QuadRT2
    public :: getEdgeDofIndex_QuadRT2
    public :: BasisReferenceQuadRT2
    public :: BasisLocalQuadRT2
    
contains

subroutine fespaceInit_QuadRT2(Vh,Th,dim)
    type(fespace), intent(out) :: Vh
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: dim
    
    
    integer :: i, i_dof

    Vh%N_local_basis = 12
    call assert(dim.eq.2, "fespaceInit_QuadRT2: dim must be 2")
    Vh%dim = dim
    Vh%basis_type = DOF_QuadRT2
    Vh%N_DOF = Th%N_edge*nDOFonEdge + Th%N_elem*nDOFonElem + Th%N_node*nDOFonNode
    Vh%isStack = 0
    allocate(Vh%ElemDOF(Vh%N_local_basis, Th%N_elem))
    
    ! numbering: edge dofs, elem dofs
    Vh%ElemDOF([1,3,5,7],:) = (abs(Th%ElemEdgeConn)-1)*nDOFonEdge + 1
    Vh%ElemDOF([2,4,6,8],:) = (abs(Th%ElemEdgeConn)-1)*nDOFonEdge + 2
    
    do i_dof = 1, nDOFonElem
        Vh%ElemDOF(8+i_dof,:) = Th%N_edge*nDOFonEdge + [(i,i=0,Th%N_elem)]*nDOFonElem + i_dof
    end do
    
    
end subroutine fespaceInit_QuadRT2

subroutine ComputeDof_QuadRT2(result,fun,Th,i_dof)
    type(mesh2D) :: Th
    procedure(func) :: fun
    real(8), intent(out) :: result
    integer,intent(in) :: i_dof

    integer :: i_edge,i_edge_in_elem, i_elem, i_in_edge, i_pt
    integer :: orientation
    
    integer, parameter :: Gauss_type_1d = LinePt4
    integer, parameter :: Gauss_type_2d = QuadPt16
    real(8), dimension(:), allocatable :: w_line_h
    real(8), dimension(:,:), allocatable :: pts_line_h
    real(8), dimension(:), allocatable :: w_line 
    real(8), dimension(:,:), allocatable :: x_line
    real(8), dimension(:,:), allocatable :: pts_elem_h
    real(8), dimension(:), allocatable :: w_elem_h
    real(8), dimension(:,:), allocatable :: pts_elem
    real(8), dimension(:), allocatable :: w_elem
    real(8), dimension(2,2) :: vertices
    real(8), dimension(:), allocatable :: J11,J12,J21,J22,detJ
    real(8), dimension(:), allocatable :: normal_hat
    real(8), dimension(:,:), allocatable :: normal
    real(8) :: length
    
    real(8), dimension(:,:), allocatable :: val, val_h
    real(8), dimension(:), allocatable :: val_tmp
    
    integer :: i_in_elem
    
    if(i_dof<=Th%N_edge*nDOFonEdge) then 
        i_edge = (i_dof-1)/nDOFonEdge + 1
        i_elem = abs(Th%EdgeElemConn(2,i_edge))
        i_edge_in_elem = Th%EdgeIdxInElem(2,i_edge)
        i_in_edge = i_dof-(i_edge-1)*nDOFonEdge
        
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
        
        if(i_in_edge==1) then
            result = sum((val(1,:)*normal(1,:) + val(2,:)*normal(2,:))*w_line*(1d0-pts_line_h(1,:)))
        else if(i_in_edge==2) then
            result = sum((val(1,:)*normal(1,:) + val(2,:)*normal(2,:))*w_line*pts_line_h(1,:))
        else
            error stop "ComputeDof_QuadRT2: i_in_edge must be 1 or 2"
        end if
    else
        
        i_elem = (i_dof-Th%N_edge*nDOFonEdge-1)/nDOFonElem + 1
        i_in_elem = i_dof - Th%N_edge*nDOFonEdge - (i_elem-1)*nDOFonElem
        
        call getGaussRefElement(Gauss_type_2d, pts_elem_h, w_elem_h)
        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type_2d, pts_elem, w_elem)
        
        call JacobianCoeffQuad(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), pts_elem_h, J11,J12,J21,J22)
        
        allocate(detJ(size(pts_elem_h,2)))
        detJ = J11*J22 - J12*J21
        
        allocate(val(2,size(pts_elem,2)))
        do i_pt = 1, size(pts_elem,2)
            call fun(pts_elem(:,i_pt), val_tmp, DERIV_NONE)
            val(:,i_pt) = val_tmp
        end do
        
        allocate(val_h(2,size(pts_elem,2)))
        val_h(1,:) = (J22*val(1,:) - J12*val(2,:))/abs(detJ)
        val_h(2,:) = (J11*val(2,:) - J21*val(1,:))/abs(detJ)
        
        if(i_in_elem==1) then
            result = sum(val_h(1,:)*(1d0-pts_elem_h(2,:))*w_elem)
        else if(i_in_elem==2) then
            result = sum(val_h(2,:)*pts_elem_h(1,:)*w_elem)
        else if(i_in_elem==3) then
            result = -sum(val_h(1,:)*(pts_elem_h(2,:))*w_elem)
        else if(i_in_elem==4) then
            result = -sum(val_h(2,:)*(1d0-pts_elem_h(1,:))*w_elem)
        else
            error stop "ComputeDof_QuadRT2: i_in_elem must be 1, 2, 3, or 4"
        end if
    end if
    

end subroutine ComputeDof_QuadRT2

subroutine getEdgeDofIndex_QuadRT2(Th, Vh, i_edge, dof_index)
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_edge
    integer, dimension(:), intent(out), allocatable :: dof_index

    allocate(dof_index(2))
    dof_index = [(i_edge-1)*2 + 1, i_edge*2]

end subroutine getEdgeDofIndex_QuadRT2


subroutine BasisReferenceQuadRT2(refpts, deriv_type, result)
    real(8), dimension(12,8):: A1,A2
    real(8), intent(in), dimension(:,:) :: refpts
    integer, intent(in) :: deriv_type
    real(8), dimension(:,:), allocatable :: b
    real(8), intent(out), dimension(:,:,:) :: result

    A1(1, :) = 0d0
    A1(2, :) = 0d0
    A1(3, :) = [0d0, -8d0, 0d0, 12d0, 12d0, 0d0, -18d0, 0d0]
    A1(4, :) = [0d0, 4d0, 0d0, -6d0, -12d0, 0d0, 18d0, 0d0]
    A1(5, :) = 0d0
    A1(6, :) = 0d0
    A1(7, :) = [2d0, -8d0, -6d0, 6d0, 24d0, 0d0, -18d0, 0d0]
    A1(8, :) = [-4d0, 16d0, 6d0, -12d0, -24d0, 0d0, 18d0, 0d0]
    A1(9, :) = [0d0, 24d0, 0d0, -24d0, -36d0, 0d0, 36d0, 0d0]
    A1(10, :) = 0d0
    A1(11, :) = [0d0, 12d0, 0d0, -12d0, -36d0, 0d0, 36d0, 0d0]
    A1(12, :) = 0d0
    
    A2(1,:) = [-4d0,6d0,16d0,0d0,-24d0,-12d0,0d0,18d0]
    A2(2,:) = [2d0,-6d0,-8d0,0d0,24d0,6d0,0d0,-18d0]
    A2(3,:) = 0d0
    A2(4,:) = 0d0
    A2(5,:) = [0d0,0d0,4d0,0d0,-12d0,-6d0,0d0,18d0]
    A2(6,:) = [0d0,0d0,-8d0,0d0,12d0,12d0,0d0,-18d0]
    A2(7,:) = 0d0
    A2(8,:) = 0d0
    A2(9,:) = 0d0
    A2(10,:) = [0d0,0d0,-12d0,0d0,36d0,12d0,0d0,-36d0]
    A2(11,:) = 0d0
    A2(12,:) = [0d0,0d0,-24d0,0d0,36d0,24d0,0d0,-36d0]
    
    allocate(b(8, size(refpts, 2)))

    select case (deriv_type)
    case (DERIV_NONE)
        b(1,:) = 1d0
        b(2,:) = refpts(1,:)
        b(3,:) = refpts(2,:)
        b(4,:) = refpts(1,:)*refpts(1,:)
        b(5,:) = refpts(1,:)*refpts(2,:)
        b(6,:) = refpts(2,:)*refpts(2,:)
        b(7,:) = refpts(1,:)*refpts(1,:)*refpts(2,:)
        b(8,:) = refpts(1,:)*refpts(2,:)*refpts(2,:)
    case (DERIV_DX)
        b(1,:) = 0d0
        b(2,:) = 1d0
        b(3,:) = 0d0
        b(4,:) = 2d0*refpts(1,:)
        b(5,:) = refpts(2,:)
        b(6,:) = 0d0
        b(7,:) = 2d0*refpts(1,:)*refpts(2,:)
        b(8,:) = refpts(2,:)*refpts(2,:)
    case (DERIV_DY)
        b(1,:) = 0d0
        b(2,:) = 0d0
        b(3,:) = 1d0
        b(4,:) = 0d0
        b(5,:) = refpts(1,:)
        b(6,:) = 2d0*refpts(2,:)
        b(7,:) = refpts(1,:)*refpts(1,:)
        b(8,:) = 2d0*refpts(1,:)*refpts(2,:)
    case default
        print *, 'BasisQuadRT2: Unknown derivative type'
        stop
    end select

    result(1,:,:) = matmul(A1, b)
    result(2,:,:) = matmul(A2, b)

end subroutine BasisReferenceQuadRT2

subroutine BasisLocalQuadRT2(refpts, Th, Vh, i_elem, deriv_type, result)
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
            call BasisReferenceQuadRT2(refpts, DERIV_NONE, basis_ref)

            ! Piola transform
            ! do i = 1, num_pts
            !     basis(:,:,i) = reshape([ &
            !     J11(i)*basis_ref(1,:,i) + J12(i)*basis_ref(2,:,i), & 
            !     J21(i)*basis_ref(1,:,i) + J22(i)*basis_ref(2,:,i) &
            !     ], [Vh%dim, Vh%N_local_basis]) / J(i)
            ! end do
            
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
            
            call BasisReferenceQuadRT2(refpts, DERIV_DX, basis_ref_dxh)
            call BasisReferenceQuadRT2(refpts, DERIV_DY, basis_ref_dyh)
            
            basis_ref_divh(1,:,:) = basis_ref_dxh(1,:,:) + basis_ref_dyh(2,:,:)
            
            result(1,:,:) = basis_ref_divh(1,:,:)/spread(J,1,num_basis)
            
        end block
        
    case default
        print *, 'BasisLocalQuadRT2: Unknown derivative type'
        stop

    end select
    
    ! Transform to physical element to make the orientation match the mesh
    do i = 1,4
        if(Th%ElemEdgeConn(i,i_elem)>0) then
            result(:,[2*i-1,2*i],:) = result(:,[2*i-1,2*i],:)
        else
            result(:,[2*i,2*i-1],:) = -result(:,[2*i-1,2*i],:)
        end if
    end do
    
end subroutine BasisLocalQuadRT2
    
end module fespace_QuadRT2