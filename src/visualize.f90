module visualize
    use settings
    use WriterVTK
    use fe_utils,only: FEfunctionQuadValue
    implicit none
    
contains


    subroutine PlotFunction(u, Th, Vh, filename)
        implicit none
        type(VECTOR), intent(in) :: u
        type(mesh2D),intent(in) :: Th
        type(fespace),intent(in) :: Vh
        character(len=*), intent(in) :: filename

        real(8), dimension(:,:), allocatable :: node_value,tmp ! value of the function at each node
        integer,dimension(:), allocatable :: node_countvalue ! number of values for each node
        real(8), dimension(:), allocatable :: node_value_1d ! value of the function at each node
        integer :: i_elem,i_dim, Gauss_type
        integer,dimension(:), allocatable :: elemprocid

        allocate(node_countvalue(Th%N_node))
        allocate(node_value(Vh%dim, Th%N_node))
        allocate(tmp(Vh%dim, Th%N_le))
        node_countvalue = 0

        if (Th%N_le.eq.3) then
            Gauss_type = TriangleNodeAverage
        else if (Th%N_le.eq.4) then
            Gauss_type = QuadNodeAverage
        else
            stop 'PlotFunction: not implemented for this element type'
        end if

        do i_elem = 1, Th%N_elem
            call FEfunctionQuadValue(u, Th, Vh, i_elem, DERIV_NONE, Gauss_type, tmp)
            node_countvalue(Th%ElemNodeConn(:,i_elem)) = node_countvalue(Th%ElemNodeConn(:,i_elem)) + 1
            node_value(:,Th%ElemNodeConn(:,i_elem)) = node_value(:,Th%ElemNodeConn(:,i_elem)) + tmp
        end do

        do i_dim = 1, Vh%dim
            node_value(i_dim,:) = node_value(i_dim,:)/node_countvalue
        end do

        allocate(node_value_1d(Th%N_node*Vh%dim))
        node_value_1d = reshape(node_value, shape(node_value_1d))

        allocate(elemprocid(Th%N_elem))
        elemprocid = 0

        call writeoutputvtk(DIM__, Th%N_elem, Th%N_node, Th%N_le, Vh%dim, Th%NodeCoord, Th%ElemNodeConn, &
        elemprocid, node_value_1d, filename)

    end subroutine


    
end module visualize