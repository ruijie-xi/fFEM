module mesh_generator
    use tools
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

subroutine QuadMeshOnSector(center,Rlist,M,Nlist,elems,nodes)
    real(8),intent(in) :: center(2) ! center of sector
    real(8),intent(in) :: Rlist(:)  ! list of radii of layers
    integer, intent(in) :: M        ! number of points along the circumference
    integer, intent(in) :: Nlist(:) ! number of points along the radius
    integer, intent(out), allocatable :: elems(:,:)
    real(8), intent(out), allocatable :: nodes(:,:)

    real(8), parameter :: pi = 3.141592653589793238462643383279d0

    ! loop variables
    integer :: i,j,k

    ! internal square
    real(8) :: R ! radius of the whole sector
    integer :: n ! number of layers
    real(8) :: S ! internal square side length
    real(8) :: P1(2), P2(2), P3(2), P4(2) ! coordinates of the square corners
    integer :: NS ! number of points along the square side 
    real(8),dimension(:,:),allocatable :: L1,L2,X1,Y1
    real(8), dimension(:), allocatable :: t
    integer, dimension(:,:),allocatable :: T1

    ! region between the square and boundary
    integer :: Nbndy
    real(8),dimension(:),allocatable :: Xbndy, Ybndy, Theta
    real(8),dimension(:),allocatable :: Xs,Ys,Xsold,Ysold
    real(8),dimension(:,:),allocatable :: X2,Y2
    integer,dimension(:,:),allocatable :: T2

    integer :: N_node, N_elem, NS_elem
    integer :: N_node_current, N_elem_current

    n = size(Rlist)
    R = Rlist(n)

    ! Square at the center
    S = R/4d0
    P1 = center + [0d0,S]
    P2 = center + [S,S]
    P3 = center
    P4 = center + [S,0d0]
    NS = nint(M/2d0)

    call linspace(0d0,1d0,NS,t)

    allocate(L1(2,NS),L2(2,NS))
    do i=1,NS
        L1(:,i) = P1 + (P2-P1)*t(i)
        L2(:,i) = P3 + (P4-P3)*t(i)
    end do

    allocate(X1(NS,NS),Y1(NS,NS))
    do i = 1,NS
        do j = 1,NS
            X1(i,j) = L1(1,i)*(1-t(j)) + L2(1,i)*t(j)
            Y1(i,j) = L1(2,i)*(1-t(j)) + L2(2,i)*t(j)
        end do
    end do

    allocate(T1(4,(NS-1)*(NS-1)))
    do i = 1,NS-1
        do j = 1,NS-1
            T1(:,(i-1)*(NS-1)+j) = [(i-1)*NS+j,i*NS+j,i*NS+j+1,(i-1)*NS+j+1]
        end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mesh region between square and circle boundary !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get boundary of square
    Nbndy = NS+NS-1
    allocate(Xbndy(Nbndy),Ybndy(Nbndy))
    do i = 1,NS
        Xbndy(i) = X1(NS,NS+1-i)
        Ybndy(i) = Y1(NS,NS+1-i)
    end do
    do i = 1,NS-1
        Xbndy(NS+i) = X1(NS-i,1)
        Ybndy(NS+i) = Y1(NS-i,1)
    end do
    
    call linspace(0d0,pi/2d0,Nbndy,Theta)

    N_node = NS*NS
    N_elem = (NS-1)*(NS-1)
    NS_elem = 2*(NS-1)
    N_node_current = N_node
    N_elem_current = N_elem
    do i = 1,n
        N_node = N_node + Nbndy*(Nlist(i)-1)
        N_elem = N_elem + NS_elem*(Nlist(i)-1)
    end do
    allocate(nodes(2,N_node),elems(4,N_elem))
    nodes(1,1:NS*NS) = reshape(X1,[NS*NS])
    nodes(2,1:NS*NS) = reshape(Y1,[NS*NS])
    elems(:,1:(NS-1)*(NS-1)) = T1

    allocate(Xs(Nbndy),Ys(Nbndy),Xsold(Nbndy),Ysold(Nbndy))
    
    deallocate(t)

    ! add layer
    do i = 1,n

        if(allocated(X2))then 
            deallocate(X2)
        end if
        if(allocated(Y2))then 
            deallocate(Y2)
        end if
        allocate(X2(Nbndy,Nlist(i)-1),Y2(Nbndy,Nlist(i)-1))

        if(i.eq.1)then
            Xsold = Xbndy
            Ysold = Ybndy
            Xs = center(1) + Rlist(i)*cos(Theta)
            Ys = center(2) + Rlist(i)*sin(Theta)
        else
            Xsold = Xs
            Ysold = Ys
            Xs = center(1) + Rlist(i)*cos(Theta)
            Ys = center(2) + Rlist(i)*sin(Theta)
        end if

        call linspace(0d0,1d0,Nlist(i),t)

        do j=1,Nbndy
            do k=1,Nlist(i)-1
                X2(j,k) = Xsold(j)*(1-t(k+1))+Xs(j)*t(k+1)
                Y2(j,k) = Ysold(j)*(1-t(k+1))+Ys(j)*t(k+1)
            end do
        end do

        if(allocated(T2))then
            deallocate(T2)
        end if
        allocate(T2(4,NS_elem*(Nlist(i)-1)))

        if (i.eq.1) then
            do j = 1,NS-1
                T2(:,j) = [NS*(NS-j),NS*(NS-j+1),NS*NS+j,NS*NS+j+1]
                T2(:,NS-1+j) = [NS-j,NS-j+1,NS*NS+j+NS-1,NS*NS+j+NS];
            end do
            do j = 2,Nlist(i)-1
                do k = 1,NS_elem
                    T2(:,(j-1)*NS_elem+k) = [NS*NS+(j-2)*(2*NS-1)+k+1,&
                    NS*NS+(j-2)*(2*NS-1)+k,NS*NS+(j-1)*(2*NS-1)+k,NS*NS+(j-1)*(2*NS-1)+k+1]
                end do
            end do
        else
            do j = 1,Nlist(i)-1
                do k = 1,NS_elem
                    T2(:,(j-1)*NS_elem+k) = [N_node_current+(j-2)*(2*NS-1)+k+1,&
                    N_node_current+(j-2)*(2*NS-1)+k,N_node_current+(j-1)*(2*NS-1)+k,&
                    N_node_current+(j-1)*(2*NS-1)+k+1]
                end do
            end do
        end if

        elems(:,N_elem_current+1:N_elem_current+NS_elem*(Nlist(i)-1)) = T2
        nodes(1,N_node_current+1:N_node_current+(Nlist(i)-1)*Nbndy) = reshape(X2,[(Nlist(i)-1)*Nbndy])
        nodes(2,N_node_current+1:N_node_current+(Nlist(i)-1)*Nbndy) = reshape(Y2,[(Nlist(i)-1)*Nbndy])

        N_node_current = N_node_current + (Nlist(i)-1)*Nbndy
        N_elem_current = N_elem_current + NS_elem*(Nlist(i)-1)

    end do

end subroutine QuadMeshOnSector



    
end module mesh_generator