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
                
        type(MATRIX_COLUMN) :: P ! Prolongation matrix
        
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

function FindPointAllElements(Th, pt) result(elems)
    use tools
    type(MESH2D), intent(in) :: Th
    real(8), dimension(:), intent(in) :: pt
    integer, dimension(:),allocatable :: elems
    
    integer, dimension(:), allocatable :: elems_temp
    
    integer :: i_box, j_box
    integer :: Nx, Ny
    real(8) :: hx, hy
    integer :: i, i_elem
    
    real(8), parameter :: eps = 1d-12
    
    call assert(pt(1)>Th%xlim(1)-eps .and. pt(1)<Th%xlim(2)+eps &
     .and. pt(2)>Th%ylim(1)-eps .and. pt(2)<Th%ylim(2)+eps, 'Point is out of domain')
     
    allocate(elems_temp(Th%N_elem))
    elems_temp = 0
    
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
            elems_temp(i_elem) = 1
        end if
    end do
    
    if(allocated(elems)) deallocate(elems)
    allocate(elems(count(elems_temp==1)))
    
    elems = pack([(i, i=1, size(elems_temp))], elems_temp==1)
    
end function

! Given pts, find the elements containing the box
function FindBoxElements(Th, pts) result(elems)
    use tools
    type(MESH2D), intent(in) :: Th
    real(8), dimension(:,:), intent(in) :: pts
    integer, allocatable :: elems(:)
    
    integer, allocatable :: elems_temp(:)
    integer :: i_box_min, j_box_min, i_box_max, j_box_max
    integer :: i_box, j_box
    integer :: Nx, Ny
    real(8) :: hx, hy
    integer :: i, i_elem
    
    real(8) :: box_xmin, box_xmax, box_ymin, box_ymax
    
    real(8), parameter :: eps = 1d-12
    
    ! Get box from points
    box_xmin = minval(pts(1,:))
    box_xmax = maxval(pts(1,:))
    box_ymin = minval(pts(2,:))
    box_ymax = maxval(pts(2,:))
    
    ! Find the boxes containing the box
    Nx = size(Th%ElemInBoxes,1)
    Ny = size(Th%ElemInBoxes,2)
    hx = (Th%xlim(2)-Th%xlim(1))/Nx
    hy = (Th%ylim(2)-Th%ylim(1))/Ny
    i_box_min = floor((box_xmin-Th%xlim(1))/hx) + 1
    i_box_max = ceiling((box_xmax-Th%xlim(1))/hx) + 1
    j_box_min = floor((box_ymin-Th%ylim(1))/hy) + 1
    j_box_max = ceiling((box_ymax-Th%ylim(1))/hy) + 1
    
    i_box_min = min(Nx, max(1, i_box_min))
    i_box_max = min(Nx, max(1, i_box_max))
    j_box_min = min(Ny, max(1, j_box_min))
    j_box_max = min(Ny, max(1, j_box_max))
    
    allocate(elems_temp(Th%N_elem))
    elems_temp = 0
    
    do i_box = i_box_min, i_box_max
        do j_box = j_box_min, j_box_max
            do i = 1, Th%ElemInBoxes(i_box,j_box)%N_elem
                i_elem = Th%ElemInBoxes(i_box,j_box)%elem_idx(i)
                elems_temp(i_elem) = 1
            end do
        end do
    end do
    
    if(allocated(elems)) deallocate(elems)
    allocate(elems(count(elems_temp==1)))
    
    elems = pack([(i, i=1, size(elems_temp))], elems_temp==1)
end function

! Given elems, return unique dofs
function GetElemsDofs(Vh, elems) result(dofs)
    use ffem_quicksort, only: unique
    type(FESPACE), intent(in) :: Vh 
    integer, dimension(:), intent(in) :: elems 
    integer, allocatable, dimension(:) :: dofs
    
    integer :: i_elem,i
    
    allocate(dofs(size(elems)*size(Vh%ElemDOF,1)))
    
    do i = 1, size(elems)
        i_elem = elems(i)
        dofs(((i-1)*size(Vh%ElemDOF,1)+1):(i*size(Vh%ElemDOF,1))) = &
        Vh%ElemDOF(:,i_elem)
    end do
    dofs = unique(dofs)
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

subroutine SetSpace(self, Th_fine_, Vh_fine_, Th_coarse_, Vh_coarse_)
    class(GridTransfer) :: self
    type(MESH2D), target :: Th_fine_, Th_coarse_
    type(FESPACE), target :: Vh_fine_, Vh_coarse_
    
    self%Th_fine => Th_fine_
    self%Th_coarse => Th_coarse_
    self%Vh_fine => Vh_fine_
    self%Vh_coarse => Vh_coarse_    
        
    self%P = GetProlongationMatrix(self%Th_fine, self%Vh_fine, self%Th_coarse, self%Vh_coarse)
        
end subroutine

subroutine FineToCoarse(self, x, y)
    class(GridTransfer) :: self
    type(VECTOR), intent(inout) :: x
    type(VECTOR) :: y
    
    type(VECTOR) :: x_temp
    
    ! call x_temp%Init(x%size)
        
    ! call SolverSolveUMFPACK2(self%M_fine, x, x_temp)
    
    call AddMultTransposeMV(1d0, self%P, x, 0d0, y)
        
    ! call OperatorTransfer(x_temp, self%Th_fine, self%Vh_fine, self%Th_coarse, self%Vh_coarse, y)
    
    ! call AddMultMV(1d0, self%M_coarse, y, 0d0, y)
    
end subroutine

subroutine CoarseToFine(self, x, y)
    class(GridTransfer) :: self
    type(VECTOR), intent(inout) :: x
    type(VECTOR) :: y
    
    call AddMultMV(1d0, self%P, x, 0d0, y)
    
    ! call OperatorTransfer(x, self%Th_coarse, self%Vh_coarse, self%Th_fine, self%Vh_fine, y)
end subroutine

function GetProlongationMatrix(Th_fine, Vh_fine, Th_coarse, Vh_coarse) result(P)
    use fe_utils, only: FEfunctionGetValue
    type(MESH2D), intent(in) :: Th_fine, Th_coarse
    type(FESPACE), intent(in) :: Vh_fine, Vh_coarse
    type(MATRIX_COLUMN) :: P
    
    type(MATRIX_TRIPLET) :: P_triplet
    
    integer :: i_elem_c, i_dof_c, i_dof_f, i_dim
    integer :: i,j 
    integer, dimension(:), allocatable :: dof_fs
    integer, dimension(:), allocatable :: elems_f
    
    real(8), dimension(2,3) :: vertices
    real(8), dimension(2) :: pt,refpt
    real(8), dimension(:,:), allocatable :: results
    Type(VECTOR) :: x_coarse
    
    call assert(Th_coarse%mesh_type==MESH_TRIANGLE .and. Th_fine%mesh_type==MESH_TRIANGLE, 'Only works for triangles')
    
    call x_coarse%Init(Vh_coarse%N_DOF)
    
    call MatrixTripletInit(P_triplet, Vh_fine%N_DOF, Vh_coarse%N_DOF, 50*Vh_fine%N_DOF)
    
    do i_dof_c = 1,Vh_coarse%N_DOF
        call FindDofLocation(Th_coarse, Vh_coarse, i_dof_c, pt, i_dim) 
        elems_f = FindPointAllElements(Th_fine, pt)
        x_coarse%data = 0d0
        x_coarse%data(i_dof_c) = 1d0
        dof_fs = GetElemsDofs(Vh_fine, elems_f)
        do j = 1, size(dof_fs)
            i_dof_f = dof_fs(j)
            call FindDofLocation(Th_fine, Vh_fine, i_dof_f, pt, i_dim)
            i_elem_c = FindPointElement(Th_coarse, pt)
            vertices = Th_coarse%NodeCoord(:,Th_coarse%ElemNodeConn(:,i_elem_c))
            refpt = GetRefpt(pt, vertices)
            results = FEfunctionGetValue(x_coarse, Th_coarse, Vh_coarse, i_elem_c, reshape(refpt, [2,1]), DERIV_NONE)
            call MatrixTripletAddValue(P_triplet, i_dof_f, i_dof_c, results(i_dim,1))
        end do
    end do
    
    call MatrixTriplet2Column(P_triplet, P)
end function GetProlongationMatrix

end module module_prolongation