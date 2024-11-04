module module_prolongation
    use settings
    use solver_umfpack2
    use module_smoother
    use matvec, only: AddMultMV
    use fe, only: FindDofLocation
    implicit none
    
    
    type GridTransfer
        type(MESH2D), pointer :: Th_fine, Th_coarse
        type(FESPACE), pointer :: Vh_fine, Vh_coarse
        
        type(MATRIX_COLUMN), pointer :: M_fine
        type(MATRIX_COLUMN), pointer :: M_coarse
        
    contains
    
        procedure :: SetSpace
        procedure :: FineToCoarse
        procedure :: CoarseToFine
        
    end type GridTransfer
    
contains

! Check if a point is in an element
! orientation of element vertices must be counter-clockwise
recursive function IsPointInElement(pt, vertices, elem_type) result(flag)
    real(8), dimension(:), intent(in) :: pt
    real(8), dimension(:,:), intent(in) :: vertices
    integer, intent(in) :: elem_type
    logical :: flag
    
    real(8), parameter :: eps = 1d-8
    real(8) :: x1, x2, x3, y1, y2, y3, x, y
    
    real(8) :: val1, val2, val3
    
    flag = .false.
    
    if(elem_type==MESH_TRIANGLE) then
        x1 = vertices(1,1)
        y1 = vertices(2,1)
        x2 = vertices(1,2)
        y2 = vertices(2,2)
        x3 = vertices(1,3)
        y3 = vertices(2,3)
        x = pt(1)
        y = pt(2)
        
        val1 = (x-x1)*(y2-y1)-(x2-x1)*(y-y1)
        val2 = (x-x2)*(y3-y2)-(x3-x2)*(y-y2)
        val3 = (x-x3)*(y1-y3)-(x1-x3)*(y-y3)
        
        flag = val1<eps .and. val2<eps .and. val3<eps
        
    else if(elem_type==MESH_QUAD) then
        
        flag = IsPointInElement(pt, vertices(:,1:3), MESH_TRIANGLE) .and. IsPointInElement(pt, vertices(:,[3,4,1]), MESH_TRIANGLE)
        
    else 
        error stop "IsPointInElement: Unsupported element type"
    end if    
end function

function FindPointElement(Th, pt) result(elem)
    use tools
    type(MESH2D), intent(in) :: Th
    real(8), dimension(:), intent(in) :: pt
    integer :: elem
    
    integer :: i_box, j_box
    integer :: Nx, Ny
    real(8) :: hx, hy
    integer :: i, i_elem
    
    real(8), parameter :: eps = 1d-12
    
    call assert(pt(1)>Th%xlim(1)-eps .and. pt(1)<Th%xlim(2)+eps &
     .and. pt(2)>Th%ylim(1)-eps .and. pt(2)<Th%ylim(2)+eps, 'Point is out of domain')
    
    ! Find the box containing the point
    Nx = size(Th%ElemInBoxes,1)
    Ny = size(Th%ElemInBoxes,2)
    hx = (Th%xlim(2)-Th%xlim(1))/Nx
    hy = (Th%ylim(2)-Th%ylim(1))/Ny
    i_box = floor((pt(1)-Th%xlim(1))/hx) + 1
    j_box = floor((pt(2)-Th%ylim(1))/hy) + 1
    
    i_box = min(Nx, max(1, i_box))
    j_box = min(Ny, max(1, j_box))
    
    do i = 1, Th%ElemInBoxes(i_box,j_box)%N_elem
        i_elem = Th%ElemInBoxes(i_box,j_box)%elem_idx(i)
        if(IsPointInElement(pt, Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Th%mesh_type)) then
            elem = i_elem
            return
        end if
    end do
    
    elem = -1
    write(*,*) "Point", pt, "not found in any element"
    error stop "FindPointElement: Point not found in any element"
end function

! Given a point pt and element vertices, return the reference point refpt
! Only works for triangles
! x = x1 + J x_hat => x_hat = J^-1 (x-x1)
! J = [x2-x1, x3-x1;
!      y2-y1, y3-y1]
function GetRefpt(pt, vertices) result(refpt)
    use tools, only: assert
    real(8), dimension(:), intent(in) :: pt
    real(8), dimension(:,:), intent(in) :: vertices
    real(8), dimension(2) :: refpt
    
    real(8) :: x1, x2, x3, y1, y2, y3, x, y
    real(8) :: J(2,2), Jinv(2,2), detJ
    
    call assert(size(vertices,2)==3, 'GetRefpt Only works for triangles')
    
    x1 = vertices(1,1)
    y1 = vertices(2,1)
    x2 = vertices(1,2)
    y2 = vertices(2,2)
    x3 = vertices(1,3)
    y3 = vertices(2,3)
    x = pt(1)
    y = pt(2)
    
    J(1,1) = x2-x1
    J(1,2) = x3-x1
    J(2,1) = y2-y1
    J(2,2) = y3-y1
    
    ! Inverse of J
    Jinv(1,1) = J(2,2)
    Jinv(1,2) = -J(1,2)
    Jinv(2,1) = -J(2,1)
    Jinv(2,2) = J(1,1)
    
    detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
    
    Jinv = Jinv/detJ
    
    refpt(1) = Jinv(1,1)*(x-x1) + Jinv(1,2)*(y-y1)
    refpt(2) = Jinv(2,1)*(x-x1) + Jinv(2,2)*(y-y1)
    
end function

! Given grid function x, interpolate it to the target grid function x_target
! only works for linear elements
subroutine OperatorTransfer(x, Th, Vh, Th_target, Vh_target, x_target)
    use tools, only: assert
    use fe_utils, only: FEfunctionGetValue
    type(MESH2D), intent(in) :: Th
    type(FESPACE), intent(in) :: Vh
    type(MESH2D), intent(in) :: Th_target
    type(FESPACE), intent(in) :: Vh_target
    type(VECTOR), intent(in) :: x
    type(VECTOR), intent(out) :: x_target
    
    integer :: i_elem, i_dof_tg, i_dim
    
    real(8), dimension(DIM__) :: tg_pt
    
    real(8), dimension(2) :: refpt
    
    real(8), dimension(2,1) :: refpts
    real(8), dimension(:,:), allocatable :: results
    
    
    call x_target%Init(Vh_target%N_DOF)
    allocate(results(1,1))
    
    do i_dof_tg = 1, Vh_target%N_DOF
        
        call FindDofLocation(Th_target, Vh_target, i_dof_tg, tg_pt, i_dim)
        
        i_elem = FindPointElement(Th, tg_pt)
        
        refpt = GetRefpt(tg_pt, Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)))
                        
        refpts = reshape(refpt, [2,1])
                        
        results = FEfunctionGetValue(x, Th, Vh, i_elem, refpts, DERIV_NONE)
        
        x_target%data(i_dof_tg) = results(i_dim,1)
        
    end do
    
end subroutine OperatorTransfer

subroutine SetSpace(self, Th_fine_, Vh_fine_, Mass_fine_, Th_coarse_, Vh_coarse_, Mass_coarse_)
    class(GridTransfer) :: self
    type(MESH2D), target :: Th_fine_, Th_coarse_
    type(FESPACE), target :: Vh_fine_, Vh_coarse_
    type(MATRIX_COLUMN), target :: Mass_fine_, Mass_coarse_
    
    self%Th_fine => Th_fine_
    self%Th_coarse => Th_coarse_
    self%Vh_fine => Vh_fine_
    self%Vh_coarse => Vh_coarse_    
    
    self%M_fine => Mass_fine_
    self%M_coarse => Mass_coarse_
    
end subroutine

subroutine FineToCoarse(self, x, y)
    class(GridTransfer) :: self
    type(VECTOR), intent(in) :: x
    type(VECTOR) :: y
    
    type(VECTOR) :: x_temp
    
    call x_temp%Init(x%size)
        
    call SolverSolveUMFPACK2(self%M_fine, x, x_temp)
        
    call OperatorTransfer(x_temp, self%Th_fine, self%Vh_fine, self%Th_coarse, self%Vh_coarse, y)
    
    call AddMultMV(1d0, self%M_coarse, y, 0d0, y)
    
end subroutine

subroutine CoarseToFine(self, x, y)
    class(GridTransfer) :: self
    type(VECTOR), intent(in) :: x
    type(VECTOR) :: y
    
    call OperatorTransfer(x, self%Th_coarse, self%Vh_coarse, self%Th_fine, self%Vh_fine, y)
end subroutine

end module module_prolongation