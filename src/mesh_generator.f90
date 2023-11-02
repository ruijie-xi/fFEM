module mesh_generator
    implicit none
    
contains

subroutine SquareMesh(left,right,bottom,top,Nx,Ny,elems,nodes)

    integer, intent(in) :: Nx, Ny
    real(8), intent(in) :: left, right, bottom, top
    integer, intent(out), allocatable :: elems(:,:)
    real(8), intent(out), allocatable :: nodes(:,:)
    integer :: row,col

    allocate(elems(4,Nx*Ny))
    allocate(nodes(2,(Nx+1)*(Ny+1)))

    do row = 1,Ny+1
        do col = 1,Nx+1
            nodes(1,map_node(row,col)) = left + (col-1)*(right-left)/Nx
            nodes(2,map_node(row,col)) = bottom + (row-1)*(top-bottom)/Ny
        end do
    end do

    do row = 1,Ny
        do col = 1,Nx
            elems(1,map_elem(row,col)) = map_node(row,col)
            elems(2,map_elem(row,col)) = map_node(row,col+1)
            elems(3,map_elem(row,col)) = map_node(row+1,col+1)
            elems(4,map_elem(row,col)) = map_node(row+1,col)
        end do
    end do

    contains
    integer function map_elem(r,c) result(elem)
        integer, intent(in) :: r, c
        elem = (r-1)*Nx + c
    end function map_elem

    integer function map_node(r,c) result(node)
        integer, intent(in) :: r, c
        node = (r-1)*(Nx+1) + c
    end function map_node

end subroutine SquareMesh

subroutine PrintElems(elems)
    integer, intent(in) :: elems(:,:)
    integer :: i_elem
    do i_elem = 1,size(elems,2)
        print *, 'elem', i_elem, ':', &
            'nodes = ', elems(:,i_elem)
    end do
end subroutine PrintElems

subroutine PrintNodes(nodes)
    real(8), intent(in) :: nodes(:,:)
    integer :: i_node
    do i_node = 1,size(nodes,2)
        print *, 'node', i_node, ':', &
            'x = ', nodes(1,i_node), &
            'y = ', nodes(2,i_node)
    end do
end subroutine PrintNodes

subroutine PrintEdges(edges)
    integer, intent(in) :: edges(:,:)
    integer :: i_edge
    do i_edge = 1,size(edges,2)
        print *, 'edge', i_edge, ':', &
            'nodes = ', edges(:,i_edge)
    end do
end subroutine PrintEdges

subroutine TriangleMesh(left,right,bottom,top,Nx,Ny,elems,nodes)
    integer, intent(in) :: Nx, Ny
    real(8), intent(in) :: left, right, bottom, top
    integer, intent(out), allocatable :: elems(:,:)
    real(8), intent(out), allocatable :: nodes(:,:)
    integer :: row,col

    allocate(elems(3,2*Nx*Ny))
    allocate(nodes(2,(Nx+1)*(Ny+1)))

    do row = 1,Ny+1
        do col = 1,Nx+1
            nodes(1,map_node(row,col)) = left + (col-1)*(right-left)/Nx
            nodes(2,map_node(row,col)) = bottom + (row-1)*(top-bottom)/Ny
        end do
    end do

    do row = 1,Ny
        do col = 1,Nx
            elems(1,map_elem(row,col,1)) = map_node(row,col)
            elems(2,map_elem(row,col,1)) = map_node(row,col+1)
            elems(3,map_elem(row,col,1)) = map_node(row+1,col)
            elems(1,map_elem(row,col,2)) = map_node(row+1,col+1)
            elems(2,map_elem(row,col,2)) = map_node(row+1,col)
            elems(3,map_elem(row,col,2)) = map_node(row,col+1)
        end do
    end do

    contains
    integer function map_elem(r,c,idx) result(elem)
        integer, intent(in) :: r, c, idx
        elem = ((r-1)*Nx+c)*2-2+idx;
    end function map_elem

    integer function map_node(r,c) result(node)
        integer, intent(in) :: r, c
        node = (r-1)*(Nx+1) + c
    end function map_node

end subroutine TriangleMesh


subroutine FindAllEdges(elems,edges)
    integer, dimension(:,:), intent(in) :: elems
    integer, dimension(:,:), intent(out), allocatable :: edges
    integer, dimension(:,:), allocatable :: edges_tmp
    integer(8), dimension(:),allocatable :: edge_list
    integer(8) :: edge_info
    integer :: i_elem, i_localedge, i_node1, i_node2, i_node, count_edges
    logical :: edge_exists, orientation

    integer, parameter :: max_num_node = 1000000000
    integer :: num_elem,num_local_edge
    
    num_elem = size(elems,2)
    num_local_edge = size(elems,1)

    allocate(edge_list(num_elem*num_local_edge))
    allocate(edges_tmp(2,num_elem*num_local_edge))

    count_edges = 0
    do i_elem = 1,num_elem
        do i_localedge = 1,num_local_edge
            ! Find global indices of nodes that define the edge
            i_node1 = elems(i_localedge, i_elem)
            i_node2 = elems(mod(i_localedge,num_local_edge)+1, i_elem)

            ! Sort global indices in ascending order
            orientation = i_node1 < i_node2
            if (i_node1 > i_node2) then
                i_node = i_node1
                i_node1 = i_node2
                i_node2 = i_node
            end if
            
            edge_info = i_node1*max_num_node + i_node2
            
            ! Check if edge already exists in edges array
            edge_exists = .false.
            if (count_edges > 0) then
                edge_exists = any(edge_list(1:count_edges)==edge_info)
            end if

            ! If edge does not exist, add it to edges array
            if (.not.edge_exists) then
                count_edges = count_edges + 1
                edge_list(count_edges) = edge_info
                edges_tmp(:,count_edges) = [i_node1,i_node2]
            end if
        end do
    end do

    allocate(edges(2,count_edges))
    edges = edges_tmp(:,1:count_edges)

end subroutine FindAllEdges



    
end module mesh_generator