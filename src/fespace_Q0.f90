module fespace_Q0
    use settings
    use quadrature, only: getGaussAnyElement
    use mesh, only: getElementArea
    implicit none
    
contains

subroutine fespaceInit_Q0(Vh,Th,dim)
    type(fespace), intent(out) :: Vh
    type(mesh2D), intent(in) :: Th
    integer, intent(in) :: dim

    integer :: i_dim,i_elem

    Vh%N_local_basis = 1
    Vh%dim = dim
    Vh%basis_type = DOF_Q0
    Vh%N_DOF = Th%N_elem * Vh%dim
    Vh%isStack = 1
    allocate(Vh%ElemDOF(Vh%N_local_basis*Vh%dim, Th%N_elem))
    Vh%ElemDOF = 0
    do i_dim = 1, Vh%dim
        Vh%ElemDOF(i_dim, :) = [(i_elem, i_elem = 1, Th%N_elem)] + (i_dim-1)*Th%N_elem
    end do
    
end subroutine fespaceInit_Q0

subroutine ComputeDof_Q0(result,fun,Th, Vh,i_dof)
    type(mesh2D) :: Th
    type(fespace) :: Vh
    procedure(func) :: fun
    real(8), intent(out) :: result
    integer,intent(in) :: i_dof

    integer :: i_elem, i_dim, i_pt
    real(8), dimension(:,:), allocatable :: func_val
    real(8), dimension(:), allocatable :: val
    real(8), dimension(:,:), allocatable :: x
    real(8), dimension(:),allocatable :: w
    real(8) :: elem_area
    
    i_dim = i_dof / Th%N_elem + 1
    i_elem = mod(i_dof, Th%N_elem)
    if (i_elem == 0) then
        i_elem = Th%N_elem
        i_dim = i_dim - 1
    end if
    
    result = 0d0
    call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)),&
    QuadPt9,x,w)

    allocate(func_val(Vh%dim, size(x,2)))
    do i_pt = 1,size(x,2)
        call fun(x(:,i_pt), val, DERIV_NONE)
        func_val(:,i_pt) = val
        deallocate(val)
    end do

    call getElementArea(Th, i_elem, elem_area)
    result = sum(w*func_val(i_dim,:))/elem_area

end subroutine ComputeDof_Q0

subroutine BasisLocalQ0(refpts, Th, Vh, i_elem, deriv_type, result)
    real(8), intent(in), dimension(:,:) :: refpts
    type(mesh2D), intent(in) :: Th
    type(fespace), intent(in) :: Vh
    integer, intent(in) :: i_elem, deriv_type
    real(8), intent(out), dimension(:,:,:), allocatable :: result


    integer :: num_pts, i_dim

    num_pts = size(refpts, 2)
    allocate(result(Vh%dim, Vh%N_local_basis, num_pts))

    select case (deriv_type)
    case (DERIV_NONE)
        result = 1d0
    case (DERIV_DX)
        result = 0d0
    case (DERIV_DY)
        result = 0d0
    case default
        print *, 'BasisLocalQ0: Unknown derivative type'
        stop

    end select

    do i_dim = 1, Vh%dim
        result(i_dim,:,:) = result(1,:,:);
    end do
end subroutine BasisLocalQ0
    
end module fespace_Q0