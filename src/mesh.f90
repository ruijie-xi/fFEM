module mesh
    use tools
    use matvec
    use readmeshfile
    use settings
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

        if (N_le .eq. 3) then
            write(*,*) "MeshInit: Triangle mesh detected."
            Th%mesh_type = MESH_TRIANGLE
        elseif (N_le .eq. 4) then
            write(*,*) "MeshInit: Quadrilateral mesh detected."
            Th%mesh_type = MESH_QUAD
        else
            write(*,*) "MeshInit: Only support triangle or quadrilateral mesh."
            stop
        end if

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
        call GenerateTopo2D_new(Th)

        ! get hmax
        call getMeshhmax(Th)
    end subroutine

    subroutine MeshFree(Th)
        type(mesh2D), intent(inout) :: Th
        deallocate(Th%ElemNodeConn)
        deallocate(Th%EdgeNodeConn)
        deallocate(Th%ElemEdgeConn)
        deallocate(Th%EdgeElemConn)
        deallocate(Th%EdgeIdxInElem)
        deallocate(Th%NodeCoord)
        deallocate(Th%BdryEdge)
        deallocate(Th%BdryMarker)
        deallocate(Th%Edge2Bdry)
    end subroutine

    subroutine PrintMesh(Th)
        type(mesh2D), intent(in) :: Th
        integer :: i_elem,i_node,i_edge,i_bdry
        write(*,*) " "
        write(*,*) " Mesh statistics ..... "
        write(*,*) " nElem    ", Th%N_elem
        write(*,*) " nNode    ", Th%N_node
        write(*,*) " nEdge    ", Th%N_edge
        write(*,*) " hmax     ", Th%hmax

        write(*,*) " "
        write(*,*) " Element connectivity "
        write(*,*) " "
        write(*,*) " Elem#   Node1   Node2   Node3   Node4"
        write(*,*) " -----   -----   -----   -----   -----"
        if (Th%N_le .eq. 3) then
            do i_elem=1,Th%N_elem
                write(*,*) i_elem, Th%ElemNodeConn(1,i_elem), Th%ElemNodeConn(2,i_elem), Th%ElemNodeConn(3,i_elem)
            end do
        elseif (Th%N_le .eq. 4) then
            do i_elem=1,Th%N_elem
                write(*,*) i_elem, Th%ElemNodeConn(1,i_elem), Th%ElemNodeConn(2,i_elem), Th%ElemNodeConn(3,i_elem),&
                 Th%ElemNodeConn(4,i_elem)
            end do
        else
            write(*,*) "PrintMesh: Only support triangle or quadrilateral mesh."
            stop
        end if

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
        write(*,*) " Edge index in elements "
        write(*,*) " "
        write(*,*) " Edge#   Elem1   Elem2"
        write(*,*) " -----   -----   -----"
        do i_edge=1,Th%N_edge
            write(*,*) i_edge, Th%EdgeIdxInElem(1,i_edge), Th%EdgeIdxInElem(2,i_edge)
        end do
        

        ! boundary information
        write(*,*) " "
        write(*,*) " Boundary Edges"
        write(*,*) " "
        write(*,*) " Bdry#   Edge   Marker"
        write(*,*) " -----   -----   -----"
        do i_bdry = 1, Th%N_bdryedge
            write(*,*) i_bdry, Th%BdryEdge(i_bdry), Th%BdryMarker(i_bdry)
        end do

    end subroutine PrintMesh

    subroutine GenerateTopo2D(Th)
        type(mesh2D), intent(inout) :: Th
        integer :: i_elem,i_le,i_edge,i_node1,i_node2,i_node_tmp,i_bdry
        integer :: factor
        integer :: count_edges
        integer, dimension(:,:), allocatable :: EdgeNodeConn
        integer, dimension(:,:), allocatable :: EdgeElemConn
        integer, dimension(:,:), allocatable :: EdgeIdxInElem
        integer, dimension(:,:), allocatable :: Edge_map

        integer, dimension(:), allocatable :: one2N_edge

        ! allocate memory
        allocate(Th%ElemEdgeConn(Th%N_le,Th%N_elem))
        allocate(EdgeElemConn(2,Th%N_elem*Th%N_le))
        allocate(EdgeIdxInElem(2,Th%N_elem*Th%N_le))
        allocate(EdgeNodeConn(2,Th%N_elem*Th%N_le))
        allocate(Edge_map(Th%N_node,Th%N_node))
        EdgeElemConn = 0
        EdgeIdxInElem = -1
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
                if (EdgeElemConn(2, i_edge) .ne. 0) then
                    EdgeElemConn(1, i_edge) = factor*i_elem
                    EdgeIdxInElem(1, i_edge) = i_le
                else
                    EdgeElemConn(2, i_edge) = factor*i_elem
                    EdgeIdxInElem(2, i_edge) = i_le
                end if
            end do
        end do

        ! Update EdgeNodeConn
        allocate(Th%EdgeNodeConn(2,count_edges))
        Th%EdgeNodeConn = EdgeNodeConn(:,1:count_edges)
        allocate(Th%EdgeElemConn(2,count_edges))
        Th%EdgeElemConn = EdgeElemConn(:,1:count_edges)
        allocate(Th%EdgeIdxInElem(2,count_edges))
        Th%EdgeIdxInElem = EdgeIdxInElem(:,1:count_edges)
        Th%N_edge = count_edges

        ! allocate the boundary edges and markers
        Th%N_bdryedge = count(Th%EdgeElemConn(1,:)==0)
        allocate(Th%BdryEdge(Th%N_bdryedge))
        allocate(Th%BdryMarker(Th%N_bdryedge))

        ! get the boundary edges
        one2N_edge = (/ (i_edge,i_edge=1,Th%N_edge) /)
        Th%BdryEdge = pack(one2N_edge,Th%EdgeElemConn(1,:)==0)
        Th%BdryMarker = 0

        allocate(Th%Edge2Bdry(Th%N_edge))
        Th%Edge2Bdry = 0
        do i_bdry = 1,Th%N_bdryedge
            Th%Edge2Bdry(Th%BdryEdge(i_bdry)) = i_bdry
        end do
        

    end subroutine GenerateTopo2D

    ! use sparse matrix to generate topology information
    subroutine GenerateTopo2D_new(Th)
        type(mesh2D), intent(inout) :: Th
        integer :: i_elem,i_le,i_edge,i_node1,i_node2,i_node_tmp,i_bdry,i_nz
        integer :: factor
        real(8) :: i_edge_tmp

        integer, dimension(:), allocatable :: one2N_edge

        type(MATRIX_TRIPLET) :: Edge_map_triplet
        type(MATRIX_COLUMN) :: Edge_map

        call MatrixTripletInit(Edge_map_triplet, Th%N_node, Th%N_node, Th%N_elem*Th%N_le)

        ! find all edges
        do i_elem=1,Th%N_elem
            do i_le = 1, Th%N_le
                i_node1 = Th%ElemNodeConn(i_le,i_elem)
                i_node2 = Th%ElemNodeConn(mod(i_le,Th%N_le)+1,i_elem)
                
                ! sort the nodes of the edge
                if (i_node1 > i_node2) then
                    i_node_tmp = i_node1
                    i_node1 = i_node2
                    i_node2 = i_node_tmp
                end if

                call MatrixTripletAddValues(Edge_map_triplet, [i_node1], [i_node2], [1d0])

            end do
        end do

        call MatrixTriplet2Column(Edge_map_triplet, Edge_map)
        call MatrixTripletFree(Edge_map_triplet)
        
        Th%N_edge = Edge_map%N_nz
        Edge_map%val = [(i_nz, i_nz=1,Th%N_edge)]

        call MatrixColumn2Triplet(Edge_map, Edge_map_triplet)

        ! generate topology information
        allocate(Th%ElemEdgeConn(Th%N_le,Th%N_elem))
        allocate(Th%EdgeNodeConn(2,Th%N_edge))
        allocate(Th%EdgeElemConn(2,Th%N_edge))
        allocate(Th%EdgeIdxInElem(2,Th%N_edge))
        Th%EdgeElemConn = 0
        Th%EdgeIdxInElem = -1

        do i_edge = 1,Th%N_edge
            Th%EdgeNodeConn(:,i_edge) = [Edge_map_triplet%row_idx(i_edge), Edge_map_triplet%col_idx(i_edge)]
        end do

        do i_elem = 1,Th%N_elem
            do i_le = 1,Th%N_le
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

                call MatrixColumnGet(Edge_map, i_node1, i_node2, i_edge_tmp)
                i_edge = int(i_edge_tmp)
                Th%ElemEdgeConn(i_le,i_elem) = factor*i_edge

                ! Update EdgeElemConn
                if (Th%EdgeElemConn(2, i_edge) .ne. 0) then
                    Th%EdgeElemConn(1, i_edge) = factor*i_elem
                    Th%EdgeIdxInElem(1, i_edge) = i_le
                else
                    Th%EdgeElemConn(2, i_edge) = factor*i_elem
                    Th%EdgeIdxInElem(2, i_edge) = i_le
                end if
            end do
        end do

        ! allocate the boundary edges and markers
        Th%N_bdryedge = count(Th%EdgeElemConn(1,:)==0)
        allocate(Th%BdryEdge(Th%N_bdryedge))
        allocate(Th%BdryMarker(Th%N_bdryedge))

        ! get the boundary edges
        one2N_edge = (/ (i_edge,i_edge=1,Th%N_edge) /)
        Th%BdryEdge = pack(one2N_edge,Th%EdgeElemConn(1,:)==0)
        Th%BdryMarker = 0

        allocate(Th%Edge2Bdry(Th%N_edge))
        Th%Edge2Bdry = 0
        do i_bdry = 1,Th%N_bdryedge
            Th%Edge2Bdry(Th%BdryEdge(i_bdry)) = i_bdry
        end do

    end subroutine GenerateTopo2D_new

    subroutine GetLength(node1,node2,edge_length)
        real(8), dimension(:), intent(in) :: node1, node2
        real(8), intent(out) :: edge_length
        integer :: N_dim,i_dim

        N_dim = size(node1)
        edge_length = 0.0d0
        do i_dim=1,N_dim
            edge_length = edge_length + (node1(i_dim)-node2(i_dim))**2
        end do
        edge_length = sqrt(edge_length)
    end subroutine GetLength

    subroutine getMeshhmax(Th)
        type(mesh2D), intent(inout) :: Th
        integer :: i_le,i_elem,j
        real(8) :: hmax,h

        if (Th%N_le .eq. 3) then
            hmax = 0.0d0
            do i_elem=1,Th%N_elem
                do i_le=1,Th%N_le
                    call GetLength(Th%NodeCoord(:,Th%ElemNodeConn(i_le,i_elem)), &
                        Th%NodeCoord(:,Th%ElemNodeConn(mod(i_le,Th%N_le)+1,i_elem)), h)
                    if (h > hmax) hmax = h
                end do
            end do
            Th%hmax = hmax
        elseif (Th%N_le .eq. 4) then
            hmax = 0.0d0
            do i_elem=1,Th%N_elem
                do i_le=1,Th%N_le
                    do j = i_le, Th%N_le
                        call GetLength(Th%NodeCoord(:,Th%ElemNodeConn(i_le,i_elem)), &
                            Th%NodeCoord(:,Th%ElemNodeConn(j,i_elem)), h)
                        if (h > hmax) hmax = h
                    end do
                end do
            end do
            Th%hmax = hmax
        else
            write(*,*) "getMeshhmax: Only support triangle or quadrilateral mesh."
            stop
        end if

        
    end subroutine getMeshhmax

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


    subroutine AddBdryMarker(Th, fun, marker)
        type(mesh2D), intent(inout) :: Th
        procedure(func) :: fun
        integer,intent(in) :: marker

        integer :: i_bdry, i_edge, i_node1, i_node2
        real(8),dimension(:),allocatable :: val1, val2
        real(8),parameter :: eps = 1d-8

        do i_bdry = 1,Th%N_bdryedge
            i_edge = Th%BdryEdge(i_bdry)
            i_node1 = Th%EdgeNodeConn(1,i_edge)
            i_node2 = Th%EdgeNodeConn(2,i_edge)
            
            call fun(Th%NodeCoord(:,i_node1), val1, DERIV_NONE)
            call fun(Th%NodeCoord(:,i_node2), val2, DERIV_NONE)

            call assert(size(val1)==1 .and. size(val2)==1, "val1 and val2 are scalars")
            if (abs(val1(1))<eps .and. abs(val2(1))<eps) then
                Th%BdryMarker(i_bdry) = marker
            end if
        end do
    end subroutine AddBdryMarker

    ! get the end points of a local line in reference element
    subroutine getRefLinePts(mesh_type, i_local_edge, vertices)
        integer, intent(in) :: mesh_type, i_local_edge
        real(8), dimension(:,:), intent(out), allocatable :: vertices

        allocate(vertices(2,2))

        if(mesh_type==MESH_TRIANGLE) then
            select case(i_local_edge)
            case(1)
                vertices(:,1) = (/0d0, 0d0/)
                vertices(:,2) = (/1d0, 0d0/)
            case(2)
                vertices(:,1) = (/1d0, 0d0/)
                vertices(:,2) = (/0d0, 1d0/)
            case(3)
                vertices(:,1) = (/0d0, 1d0/)
                vertices(:,2) = (/0d0, 0d0/)
            case default
                write(*,*) "getRefLinePts: i_local_edge not supported. "
                write(*,*) "i_local_edge = ",i_local_edge
                stop
            end select
        elseif (mesh_type==MESH_QUAD) then
            select case(i_local_edge)
            case(1)
                vertices(:,1) = (/0d0, 0d0/)
                vertices(:,2) = (/1d0, 0d0/)
            case(2)
                vertices(:,1) = (/1d0, 0d0/)
                vertices(:,2) = (/1d0, 1d0/)
            case(3)
                vertices(:,1) = (/1d0, 1d0/)
                vertices(:,2) = (/0d0, 1d0/)
            case(4)
                vertices(:,1) = (/0d0, 1d0/)
                vertices(:,2) = (/0d0, 0d0/)
            case default
                write(*,*) "getRefLinePts: i_local_edge not supported. "
                write(*,*) "i_local_edge = ",i_local_edge
                stop
            end select
        else
            write(*,*) "getRefLinePts: mesh_type not supported. "
        end if
    end subroutine getRefLinePts

    ! get the end points of a local line in any element
    subroutine getAnyLinePts(Th, i_elem, i_local_edge, vertices)
        type(mesh2D) :: Th
        integer, intent(in) :: i_local_edge,i_elem
        real(8), dimension(:,:), intent(out), allocatable :: vertices

        allocate(vertices(2,2))

        vertices = Th%NodeCoord(:,Th%ElemNodeConn([i_local_edge,mod(i_local_edge,Th%N_le)+1],i_elem))

    end subroutine getAnyLinePts

    subroutine getElementArea(Th, i_elem, area)
        type(mesh2D), intent(in) :: Th
        integer, intent(in) :: i_elem
        real(8), intent(out) :: area

        real(8) :: x1,x2,x3,x4,y1,y2,y3,y4
        
        select case(Th%mesh_type)
        case(MESH_TRIANGLE)
            x1 = Th%NodeCoord(1, Th%ElemNodeConn(1, i_elem))
            x2 = Th%NodeCoord(1, Th%ElemNodeConn(2, i_elem))
            x3 = Th%NodeCoord(1, Th%ElemNodeConn(3, i_elem))
            y1 = Th%NodeCoord(2, Th%ElemNodeConn(1, i_elem))
            y2 = Th%NodeCoord(2, Th%ElemNodeConn(2, i_elem))
            y3 = Th%NodeCoord(2, Th%ElemNodeConn(3, i_elem))
            area = 5d-1*abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        case(MESH_QUAD)
            x1 = Th%NodeCoord(1, Th%ElemNodeConn(1, i_elem))
            x2 = Th%NodeCoord(1, Th%ElemNodeConn(2, i_elem))
            x3 = Th%NodeCoord(1, Th%ElemNodeConn(3, i_elem))
            x4 = Th%NodeCoord(1, Th%ElemNodeConn(4, i_elem))
            y1 = Th%NodeCoord(2, Th%ElemNodeConn(1, i_elem))
            y2 = Th%NodeCoord(2, Th%ElemNodeConn(2, i_elem))
            y3 = Th%NodeCoord(2, Th%ElemNodeConn(3, i_elem))
            y4 = Th%NodeCoord(2, Th%ElemNodeConn(4, i_elem))
            area = 5d-1*abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) &
            + 5d-1*abs((x3-x1)*(y4-y1) - (x4-x1)*(y3-y1))
        case default
            write(*,*) "getElementArea: mesh_type not supported. "
        end select

    end subroutine getElementArea

    subroutine getEdgeNormal(Th, i_edge, normal)
        type(mesh2D), intent(in) :: Th
        integer, intent(in) :: i_edge
        real(8), dimension(2), intent(out) :: normal

        real(8) :: x1,x2,y1,y2

        x1 = Th%NodeCoord(1,Th%EdgeNodeConn(1,i_edge))
        x2 = Th%NodeCoord(1,Th%EdgeNodeConn(2,i_edge))
        y1 = Th%NodeCoord(2,Th%EdgeNodeConn(1,i_edge))
        y2 = Th%NodeCoord(2,Th%EdgeNodeConn(2,i_edge))

        normal(1) = y2-y1
        normal(2) = x1-x2
        normal = normal/sqrt(normal(1)**2+normal(2)**2)

    end subroutine getEdgeNormal

    subroutine getEdgeTangent(Th, i_edge, tangent)
        type(mesh2D), intent(in) :: Th
        integer, intent(in) :: i_edge
        real(8), dimension(2), intent(out) :: tangent

        real(8) :: x1,x2,y1,y2

        x1 = Th%NodeCoord(1,Th%EdgeNodeConn(1,i_edge))
        x2 = Th%NodeCoord(1,Th%EdgeNodeConn(2,i_edge))
        y1 = Th%NodeCoord(2,Th%EdgeNodeConn(1,i_edge))
        y2 = Th%NodeCoord(2,Th%EdgeNodeConn(2,i_edge))

        tangent(1) = x2-x1
        tangent(2) = y2-y1
        tangent = tangent/sqrt(tangent(1)**2+tangent(2)**2)

    end subroutine getEdgeTangent

    subroutine getEdgeLength(Th, i_edge, length)
        type(mesh2D), intent(in) :: Th
        integer, intent(in) :: i_edge
        real(8), intent(out) :: length

        real(8) :: x1,x2,y1,y2

        x1 = Th%NodeCoord(1,Th%EdgeNodeConn(1,i_edge))
        x2 = Th%NodeCoord(1,Th%EdgeNodeConn(2,i_edge))
        y1 = Th%NodeCoord(2,Th%EdgeNodeConn(1,i_edge))
        y2 = Th%NodeCoord(2,Th%EdgeNodeConn(2,i_edge))

        length = sqrt((x2-x1)**2+(y2-y1)**2)

    end subroutine getEdgeLength

    
end module mesh