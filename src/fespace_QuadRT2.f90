module fespace_QuadRT2
    use settings
    use tools, only: assert
    use mesh, only: getEdgeTangent,getEdgeLength,getEdgeNormal
    use quadrature
    implicit none
    
contains

subroutine fespaceInit_QuadRT2(Vh,Th,dim)
    type(fespace), intent(out) :: Vh
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: dim
    
    integer,parameter :: nDOFonEdge = 2
    integer,parameter :: nDOFonElem = 4
    integer :: i, i_dof

    Vh%N_local_basis = 12
    call assert(dim.eq.2, "fespaceInit_QuadRT2: dim must be 2")
    Vh%dim = dim
    Vh%basis_type = DOF_QuadRT2
    Vh%N_DOF = Th%N_edge*nDOFonEdge + Th%N_elem*nDOFonElem
    Vh%isStack = 0
    allocate(Vh%ElemDOF(Vh%N_local_basis, Th%N_elem))
    
    ! numbering: edge dofs, elem dofs
    Vh%ElemDOF([1,3,5,7],:) = (abs(Th%ElemEdgeConn)-1)*nDOFonEdge + 1
    Vh%ElemDOF([2,4,6,8],:) = (abs(Th%ElemEdgeConn)-1)*nDOFonEdge + 2
    
    do i_dof = 1, nDOFonElem
        Vh%ElemDOF(8+i_dof,:) = Th%N_edge*nDOFonEdge + [(i,i=0,Th%N_elem)]*nDOFonElem + i_dof
    end do
    
    
end subroutine fespaceInit_QuadRT2

subroutine ComputeDof_QuadRT1(result,fun,Th,i_dof)
    type(mesh2D) :: Th
    procedure(func) :: fun
    real(8), intent(out) :: result
    integer,intent(in) :: i_dof

    integer :: i_edge,i_pt
    real(8), dimension(2,2) :: vertices
    real(8), dimension(:), allocatable :: w
    real(8), dimension(:,:), allocatable :: x
    real(8), dimension(:,:), allocatable :: val
    real(8), dimension(:), allocatable :: val_tmp
    real(8), dimension(2) :: normal
    real(8) :: length

    i_edge = i_dof
    vertices = Th%NodeCoord(:, Th%EdgeNodeConn(:, i_edge))

    call getGaussQuadAnyLine(vertices, LinePt3, x, w)
    allocate(val(2,size(x,2)))
    do i_pt = 1, size(x,2)
        call fun(x(:,i_pt), val_tmp, DERIV_NONE)
        val(:,i_pt) = val_tmp
    end do
    call getEdgeNormal(Th, i_edge, normal)
    call getEdgeLength(Th, i_edge, length)

    result = sum(w*val(1,:)*normal(1) + w*val(2,:)*normal(2))/length
    

end subroutine ComputeDof_QuadRT1

subroutine getEdgeDofIndex_QuadRT2(Th, Vh, i_edge, dof_index)
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_edge
    integer, dimension(:), intent(out), allocatable :: dof_index

    allocate(dof_index(2))
    dof_index = [(i_edge-1)*2 + 1, i_edge*2]

end subroutine getEdgeDofIndex_QuadRT2


subroutine BasisReferenceQuadRT1(refpts, deriv_type, result)
    real(8), dimension(12,8):: A1,A2
    real(8), intent(in), dimension(:,:) :: refpts
    integer, intent(in) :: deriv_type
    real(8), dimension(:,:), allocatable :: b
    real(8), intent(out), dimension(:,:,:) :: result

    A1(1, :) = 0d0
    A1(2, :) = 0d0
    A1(3, :) = [0d0, -8d0, 0d0, 12d0, 12d0, 0d0, 0d0, -18d0]
    A1(4, :) = [0d0, 4d0, 0d0, -6d0, -12d0, 0d0, 0d0, 18d0]
    A1(5, :) = 0d0
    A1(6, :) = 0d0
    A1(7, :) = [2d0, -8d0, -6d0, 6d0, 24d0, 0d0, 0d0, -18d0]
    A1(8, :) = [-4d0, 16d0, 6d0, -12d0, -24d0, 0d0, 0d0, 18d0]
    A1(9, :) = 0d0
    A1(10, :) = 0d0
    A1(11, :) = [0d0, 24d0, 0d0, -24d0, -36d0, 0d0, 0d0, 36d0]
    A1(12, :) = [0d0, -12d0, 0d0, 12d0, 36d0, 0d0, 0d0, -36d0]
    
    allocate(b(8, size(refpts, 2)))

    select case (deriv_type)
    case (DERIV_NONE)
        b(1,:) = 1d0
        b(2,:) = refpts(1,:)
        b(3,:) = refpts(2,:)
        b(4,:) = refpts(1,:)*refpts(1,:)
        b(5,:) = refpts(1,:)*refpts(2,:)
        b(6,:) = refpts(2,:)*refpts(2,:)
        b(7,:) = refpts(1,:)*refpts(2,:)*refpts(2,:)
        b(8,:) = refpts(1,:)*refpts(1,:)*refpts(2,:)
    case (DERIV_DX)
        b(1,:) = 0d0
        b(2,:) = 1d0
        b(3,:) = 0d0
        b(4,:) = 2d0*refpts(1,:)
        b(5,:) = refpts(2,:)
        b(6,:) = 0d0
        b(7,:) = refpts(2,:)*refpts(2,:)
        b(8,:) = 2d0*refpts(1,:)*refpts(2,:)
    case (DERIV_DY)
        b(1,:) = 0d0
        b(2,:) = 0d0
        b(3,:) = 1d0
        b(4,:) = 0d0
        b(5,:) = refpts(1,:)
        b(6,:) = 2d0*refpts(2,:)
        b(7,:) = refpts(1,:)*2d0*refpts(2,:)
        b(8,:) = refpts(1,:)*refpts(1,:)
    case default
        print *, 'BasisQuadRT2: Unknown derivative type'
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

    real(8), dimension(:,:,:),allocatable :: basis,basis_ref,basis_ref_dx,basis_ref_dy
    real(8), dimension(:,:,:),allocatable :: basis_dx,basis_dy,basis_ref_dxh,basis_ref_dyh

    real(8), dimension(2,4) :: vertices
    real(8), dimension(:), allocatable :: J11,J12,J21,J22,J

    real(8), dimension(4) :: signedlength
    real(8) :: arclength

    integer :: num_pts, num_basis, i

    num_basis = Vh%N_local_basis
    num_pts = size(refpts, 2) ! number of points in reference element
    allocate(result(Vh%dim, Vh%N_local_basis, num_pts)) ! result is a 3D array
    
    ! Jacobian coefficients
    vertices = Th%NodeCoord(:, Th%ElemNodeConn(:, i_elem))
    call JacobianCoeffQuad(vertices, refpts, J11,J12,J21,J22)
    allocate(J(size(refpts,2)))
    J = J11*J22 - J12*J21

    ! signed length
    do i = 1,4
        call getEdgeLength(Th, abs(Th%ElemEdgeConn(i,i_elem)), arclength)
        if (Th%ElemEdgeConn(i,i_elem) < 0) then
            signedlength(i) = -arclength
        else
            signedlength(i) = arclength
        end if
    end do

    select case (deriv_type)
    case (DERIV_NONE)
        allocate(basis_ref(Vh%dim, Vh%N_local_basis, num_pts))
        allocate(basis(Vh%dim, Vh%N_local_basis, num_pts))
        call BasisReferenceQuadRT1(refpts, DERIV_NONE, basis_ref)

        ! Piola transform
        basis(1,:,:) = spread(J11,1,num_basis)*basis_ref(1,:,:)/spread(J,1,num_basis) &
        + spread(J12,1,num_basis)*basis_ref(2,:,:)/spread(J,1,num_basis)
        basis(2,:,:) = spread(J21,1,num_basis)*basis_ref(1,:,:)/spread(J,1,num_basis) &
        + spread(J22,1,num_basis)*basis_ref(2,:,:)/spread(J,1,num_basis)

        ! Transform to physical element
        do i = 1,num_basis
            result(:,i,:) = signedlength(i)*basis(:,i,:)
        end do
        
    case (DERIV_DX)
        allocate(basis_ref_dxh(Vh%dim, Vh%N_local_basis, num_pts))
        allocate(basis_ref_dyh(Vh%dim, Vh%N_local_basis, num_pts))
        allocate(basis_ref_dx(Vh%dim, Vh%N_local_basis, num_pts))
        allocate(basis_dx(Vh%dim, Vh%N_local_basis, num_pts))

        call BasisReferenceQuadRT1(refpts, DERIV_DX, basis_ref_dxh)
        call BasisReferenceQuadRT1(refpts, DERIV_DY, basis_ref_dyh)
        basis_ref_dx(1,:,:) = spread(J22,1,num_basis)*basis_ref_dxh(1,:,:)/spread(J,1,num_basis) &
        - spread(J21,1,num_basis)*basis_ref_dyh(1,:,:)/spread(J,1,num_basis)
        basis_ref_dx(2,:,:) = spread(J22,1,num_basis)*basis_ref_dxh(1,:,:)/spread(J,1,num_basis) &
        - spread(J21,1,num_basis)*basis_ref_dyh(1,:,:)/spread(J,1,num_basis)

        ! Piola transform
        basis_dx(1,:,:) = spread(J11,1,num_basis)*basis_ref_dx(1,:,:)/spread(J,1,num_basis) &
        + spread(J12,1,num_basis)*basis_ref_dx(2,:,:)/spread(J,1,num_basis)
        basis_dx(2,:,:) = spread(J21,1,num_basis)*basis_ref_dx(1,:,:)/spread(J,1,num_basis) &
        + spread(J22,1,num_basis)*basis_ref_dx(2,:,:)/spread(J,1,num_basis)

        ! Transform to physical element
        do i = 1,num_basis
            result(:,i,:) = signedlength(i)*basis_dx(:,i,:)
        end do

    case (DERIV_DY)
        allocate(basis_ref_dxh(Vh%dim, Vh%N_local_basis, num_pts))
        allocate(basis_ref_dyh(Vh%dim, Vh%N_local_basis, num_pts))
        allocate(basis_ref_dy(Vh%dim, Vh%N_local_basis, num_pts))
        allocate(basis_dy(Vh%dim, Vh%N_local_basis, num_pts))

        call BasisReferenceQuadRT1(refpts, DERIV_DX, basis_ref_dxh)
        call BasisReferenceQuadRT1(refpts, DERIV_DY, basis_ref_dyh)
        basis_ref_dy(1,:,:) = -spread(J12,1,num_basis)*basis_ref_dxh(1,:,:)/spread(J,1,num_basis) &
        + spread(J11,1,num_basis)*basis_ref_dyh(1,:,:)/spread(J,1,num_basis)
        basis_ref_dy(2,:,:) = -spread(J12,1,num_basis)*basis_ref_dxh(2,:,:)/spread(J,1,num_basis) &
        + spread(J11,1,num_basis)*basis_ref_dyh(2,:,:)/spread(J,1,num_basis)

        ! Piola transform
        basis_dy(1,:,:) = spread(J11,1,num_basis)*basis_ref_dy(1,:,:)/spread(J,1,num_basis) &
        + spread(J12,1,num_basis)*basis_ref_dy(2,:,:)/spread(J,1,num_basis)
        basis_dy(2,:,:) = spread(J21,1,num_basis)*basis_ref_dy(1,:,:)/spread(J,1,num_basis) &
        + spread(J22,1,num_basis)*basis_ref_dy(2,:,:)/spread(J,1,num_basis)

        ! Transform to physical element
        do i = 1,num_basis
            result(:,i,:) = signedlength(i)*basis_dy(:,i,:)
        end do
        
    case default
        print *, 'BasisLocalQuadRT1: Unknown derivative type'
        stop

    end select
    
end subroutine BasisLocalQuadRT1


subroutine JacobianCoeffQuad(vertices, refpts, J11,J12,J21,J22)
    real(8), dimension(:,:), intent(in) :: vertices
    real(8), dimension(:,:), intent(in) :: refpts
    real(8), dimension(:), intent(out), allocatable :: J11,J12,J21,J22

    real(8) :: x1,x2,x3,x4,y1,y2,y3,y4
    real(8) :: cx1,cx2,cx3,cx4,cy1,cy2,cy3,cy4

    x1 = vertices(1,1)
    x2 = vertices(1,2)
    x3 = vertices(1,3)
    x4 = vertices(1,4)
    y1 = vertices(2,1)
    y2 = vertices(2,2)
    y3 = vertices(2,3)
    y4 = vertices(2,4)

    cx1 = x1 - x2 + x3 - x4
    cx2 = -x1 + x2
    cx3 = -x1 + x4
    cx4 = x1

    cy1 = y1 - y2 + y3 - y4
    cy2 = -y1 + y2
    cy3 = -y1 + y4
    cy4 = y1

    allocate(J11(size(refpts,2)))
    allocate(J12(size(refpts,2)))
    allocate(J21(size(refpts,2)))
    allocate(J22(size(refpts,2)))

    J11 = cx1*refpts(2,:) + cx2
    J12 = cx1*refpts(1,:) + cx3
    J21 = cy1*refpts(2,:) + cy2
    J22 = cy1*refpts(1,:) + cy3
end subroutine JacobianCoeffQuad
    
end module fespace_QuadRT2