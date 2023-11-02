module mesh
    use readmeshfile
    use settings
    use solver, only: CreateMat
    implicit none
contains

    subroutine MeshInit(Elems,Nodes,Th)
        type(mesh2D), intent(out) :: Th
        integer, dimension(:,:), intent(inout) :: Elems
        real(8), dimension(:,:), intent(in) :: Nodes

        integer :: N_node, N_elem, N_le

        N_node = size(Nodes,2)
        N_elem = size(Elems,2)
        N_le = size(Elems,1)

        ! check the orientation of the elements (triangle only)
        if (N_le .eq. 3) then
            call CheckOrientation(Elems,Nodes)
        else
            write(*,*) "Warning: Orientation check only support triangle mesh. not checked."
        end if

        ! initialze Th
        Th%N_node = N_node
        Th%N_elem = N_elem
        Th%N_le = N_le
        allocate(Th%ElemNodeConn(N_le,N_elem))
        Th%ElemNodeConn = Elems
        allocate(Th%NodeCoord(DIM__,N_node))
        Th%NodeCoord = Nodes

        ! generate topology information, including EdgeNodeConn, ElemEdgeConn, EdgeElemConn
        call GenerateTopo2D(Th)

        ! get hmax
        call GetMeshhmax(Th)
    end subroutine

    subroutine PrintMesh(Th)
        type(mesh2D), intent(in) :: Th
        integer :: i_elem,i_node,i_edge
        write(*,*) " "
        write(*,*) " Mesh statistics ..... "
        write(*,*) " nElem    ", Th%N_elem
        write(*,*) " nNode    ", Th%N_node
        write(*,*) " nEdge    ", Th%N_edge
        write(*,*) " hmax     ", Th%hmax

        write(*,*) " "
        write(*,*) " Element connectivity "
        write(*,*) " "
        write(*,*) " Elem#   Node1   Node2   Node3"
        write(*,*) " -----   -----   -----   -----"
        do i_elem=1,Th%N_elem
            write(*,*) i_elem, Th%ElemNodeConn(1,i_elem), Th%ElemNodeConn(2,i_elem), Th%ElemNodeConn(3,i_elem)
        end do

        write(*,*) " "
        write(*,*) " Node coordinates "
        write(*,*) " "
        write(*,*) " Node#   x-coord   y-coord"
        write(*,*) " -----   -------   -------"
        do i_node=1,Th%N_node
            write(*,*) i_node, Th%NodeCoord(1,i_node), Th%NodeCoord(2,i_node)
        end do

        write(*,*) " "
        write(*,*) " Edge connectivity "
        write(*,*) " "
        write(*,*) " Edge#   Node1   Node2"
        write(*,*) " -----   -----   -----"
        do i_edge=1,Th%N_edge
            write(*,*) i_edge, Th%EdgeNodeConn(1,i_edge), Th%EdgeNodeConn(2,i_edge)
        end do

        write(*,*) " "
        write(*,*) " Element-Edge connectivity "
        write(*,*) " "
        write(*,*) " Elem#   Edge1   Edge2   Edge3"
        write(*,*) " -----   -----   -----   -----"
        do i_elem=1,Th%N_elem
            write(*,*) i_elem, Th%ElemEdgeConn(1,i_elem), Th%ElemEdgeConn(2,i_elem), Th%ElemEdgeConn(3,i_elem)
        end do

        write(*,*) " "
        write(*,*) " Edge-Element connectivity "
        write(*,*) " "
        write(*,*) " Edge#   Elem1   Elem2"
        write(*,*) " -----   -----   -----"
        do i_edge=1,Th%N_edge
            write(*,*) i_edge, Th%EdgeElemConn(1,i_edge), Th%EdgeElemConn(2,i_edge)
        end do

        write(*,*) " "
        write(*,*) " Edge markers "
        write(*,*) " "
        write(*,*) " Edge#   Marker"
        write(*,*) " -----   -----"
        do i_edge=1,Th%N_edge
            write(*,*) i_edge, Th%EdgeMarker(i_edge)
        end do
        


    end subroutine PrintMesh

    subroutine GenerateTopo2D(Th)
        type(mesh2D), intent(inout) :: Th
        integer :: i_elem,i_le,i_edge,i_node1,i_node2,i_node_tmp
        integer :: factor
        integer :: count_edges
        integer, dimension(:,:), allocatable :: EdgeNodeConn
        integer, dimension(:,:), allocatable :: EdgeElemConn
        integer, dimension(:,:), allocatable :: Edge_map

        ! allocate memory
        allocate(Th%ElemEdgeConn(Th%N_le,Th%N_elem))
        allocate(EdgeElemConn(2,Th%N_elem*Th%N_le))
        allocate(EdgeNodeConn(2,Th%N_elem*Th%N_le))
        allocate(Edge_map(Th%N_node,Th%N_node))
        EdgeElemConn = -1
        Edge_map = 0

        count_edges = 0
        do i_elem=1,Th%N_elem
            do i_le = 1, Th%N_le
                factor = 1
                i_node1 = Th%ElemNodeConn(i_le,i_elem)
                i_node2 = Th%ElemNodeConn(mod(i_le,Th%N_le)+1,i_elem)
                
                ! sort the nodes of the edge
                if (i_node1 > i_node2) then
                    i_node_tmp = i_node1
                    i_node1 = i_node2
                    i_node2 = i_node_tmp
                    factor = -1
                end if

                ! search for the edge
                if (Edge_map(i_node1,i_node2).eq.0) then
                    count_edges = count_edges + 1
                    EdgeNodeConn(:,count_edges) = (/i_node1,i_node2/)
                    Edge_map(i_node1,i_node2) = count_edges
                end if

                ! Update ElemEdgeConn
                i_edge = Edge_map(i_node1,i_node2)
                Th%ElemEdgeConn(i_le,i_elem) = factor*Edge_map(i_node1,i_node2)

                ! Update EdgeElemConn
                if (EdgeElemConn(2, i_edge) > 0) then
                    EdgeElemConn(1, i_edge) = i_elem
                else
                    EdgeElemConn(2, i_edge) = i_elem
                end if
            end do
        end do

        ! Update EdgeNodeConn
        allocate(Th%EdgeNodeConn(2,count_edges))
        Th%EdgeNodeConn = EdgeNodeConn(:,1:count_edges)
        allocate(Th%EdgeElemConn(2,count_edges))
        Th%EdgeElemConn = EdgeElemConn(:,1:count_edges)
        Th%N_edge = count_edges

        ! allocate the edge markers
        allocate(Th%EdgeMarker(Th%N_edge))
        Th%EdgeMarker = 0

        where (Th%EdgeElemConn(1,:) .eq. -1)
            Th%EdgeMarker = 1
        end where

    end subroutine GenerateTopo2D

    subroutine GetEdgeLength(node1,node2,edge_length)
        real(8), dimension(:), intent(in) :: node1, node2
        real(8), intent(out) :: edge_length
        integer :: N_dim,i_dim

        N_dim = size(node1)
        edge_length = 0.0d0
        do i_dim=1,N_dim
            edge_length = edge_length + (node1(i_dim)-node2(i_dim))**2
        end do
        edge_length = sqrt(edge_length)
    end subroutine GetEdgeLength

    subroutine GetMeshhmax(Th)
        type(mesh2D), intent(inout) :: Th
        integer :: i_le,i_elem
        real(8) :: hmax,h

        if (Th%N_le .eq. 3) then
            hmax = 0.0d0
            do i_elem=1,Th%N_elem
                do i_le=1,Th%N_le
                    call GetEdgeLength(Th%NodeCoord(:,Th%ElemNodeConn(i_le,i_elem)), &
                        Th%NodeCoord(:,Th%ElemNodeConn(mod(i_le,Th%N_le)+1,i_elem)), h)
                    if (h > hmax) then
                        hmax = h
                    end if
                end do
            end do
            Th%hmax = hmax
        elseif (Th%N_le .eq. 4) then
            write(*,*) "GetMeshhmax: Not implemented yet."
        else
            write(*,*) "GetMeshhmax: Only support triangle or quadrilateral mesh."
            stop
        end if

        
    end subroutine GetMeshhmax

    !Check all the triangles has positive area
    subroutine CheckOrientation(ElemNodeConn,NodeCoord)
        integer, dimension(:,:), intent(inout) :: ElemNodeConn
        real(8), dimension(:,:), intent(in) :: NodeCoord
        real(8) :: area
        integer :: i_elem,count,N_elem,tmp
        real(8) :: x1,x2,x3,y1,y2,y3

        if (size(ElemNodeConn,1) .ne. 3) then
            write(*,*) "CheckOrientation: Only support triangle mesh."
            stop
        end if

        count = 0
        N_elem = size(ElemNodeConn,2)
        do i_elem = 1,N_elem
            x1 = NodeCoord(1, ElemNodeConn(1, i_elem))
            x2 = NodeCoord(1, ElemNodeConn(2, i_elem))
            x3 = NodeCoord(1, ElemNodeConn(3, i_elem))
            y1 = NodeCoord(2, ElemNodeConn(1, i_elem))
            y2 = NodeCoord(2, ElemNodeConn(2, i_elem))
            y3 = NodeCoord(2, ElemNodeConn(3, i_elem))

            area = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
            if (area < 0) then
                tmp = ElemNodeConn(1,i_elem)
                ElemNodeConn(1,i_elem) = ElemNodeConn(2,i_elem)
                ElemNodeConn(2,i_elem) = tmp
                count = count + 1
            end if
        end do

        write(*,*) "Orientation Checked. ", count, "elements fliped."

    end subroutine CheckOrientation

    
end module mesh