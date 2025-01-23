module geometry
    use settings
    implicit none
    
contains 

! return the normal vector to the edge of the reference element
!       3
!   <--------
!   |       ^
! 4 |       | 2
!   |       |
!   -------->
!       1
subroutine ref_edge_normal(mesh_type, i_local_edge, orientation, normal)
integer, intent(in) :: mesh_type, i_local_edge, orientation
real(8), dimension(:), intent(out), allocatable :: normal

allocate(normal(2))
if(mesh_type==MESH_QUAD) then
    
    if(i_local_edge==1) then
        normal = [0.0d0, -1.0d0]
    else if(i_local_edge==2) then
        normal = [1.0d0, 0.0d0]
    else if(i_local_edge==3) then
        normal = [0.0d0, 1.0d0]
    else if(i_local_edge==4) then
        normal = [-1.0d0, 0.0d0]
    else 
        error stop "ref_edge_normal: i_local_edge must be 1,2,3,4"
    end if
    
    if(orientation==-1) normal = -normal
    
else 
    error stop "ref_edge_normal: mesh_type not supported yet"
end if
        
end subroutine ref_edge_normal

! given the refpts on the edge, return the refpts on the element
!       3
!   <--------
!   |       ^
! 4 |       | 2
!   |       |
!   -------->
!       1
! edge_pts_h: 1 x num_pts
! elem_pts_h: 2 x num_pts
subroutine ref_pts_edge_to_elem(mesh_type, i_local_edge, edge_pts_h, orientation, elem_pts_h)
integer, intent(in) :: mesh_type, i_local_edge, orientation
real(8), dimension(:,:), intent(in) :: edge_pts_h
real(8), dimension(:,:), intent(out), allocatable :: elem_pts_h

real(8), dimension(:,:), allocatable :: edge_pts 

integer :: num_pts

num_pts = size(edge_pts_h)

allocate(edge_pts(1, num_pts))
edge_pts = edge_pts_h

if(orientation<0) edge_pts = 1.0d0-edge_pts

allocate(elem_pts_h(2,num_pts))

if(mesh_type==MESH_QUAD) then

    if(i_local_edge==1) then
        elem_pts_h(1,:) = edge_pts(1,:)
        elem_pts_h(2,:) = 0.0d0
    else if(i_local_edge==2) then
        elem_pts_h(1,:) = 1.0d0
        elem_pts_h(2,:) = edge_pts(1,:)
    else if(i_local_edge==3) then
        elem_pts_h(1,:) = 1.0d0-edge_pts(1,:)
        elem_pts_h(2,:) = 1.0d0
    else if(i_local_edge==4) then
        elem_pts_h(1,:) = 0.0d0
        elem_pts_h(2,:) = 1.0d0-edge_pts(1,:)
    else 
        error stop "ref_pts_edge_to_elem: i_local_edge must be 1,2,3,4"
    end if
else
    error stop "ref_pts_edge_to_elem: mesh_type not supported yet"
end if

end subroutine
    

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
    
end module geometry