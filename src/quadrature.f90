module quadrature
    use settings
    implicit none
    
contains

    subroutine getGaussRefElement(Gauss_type, x, w)
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        if (Gauss_type > 100 .and. Gauss_type < 200) then
            call getGaussQuadRefTriangle(Gauss_type, x, w)
        else if (Gauss_type > 200) then
            call getGaussQuadRefQuad(Gauss_type, x, w)
        else
            print *, "Gauss quadrature type not implemented"
            stop
        end if
    end subroutine getGaussRefElement

    subroutine getGaussAnyElement(vertices, Gauss_type, x, w)
        real(8),intent(in), dimension(:,:) :: vertices
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        if (Gauss_type > 100 .and. Gauss_type < 200) then
            call getGaussQuadAnyTriangle(vertices, Gauss_type, x, w)
        elseif (Gauss_type > 200) then
            call getGaussQuadAnyQuad(vertices, Gauss_type, x, w)
        else
            print *, "Gauss quadrature type not implemented"
            stop
        end if
    end subroutine getGaussAnyElement


    subroutine getGaussQuadRefTriangle(Gauss_type, x, w)
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        integer :: N_pts

        select case (Gauss_type)
        case (TriangleNodeAverage)
            N_pts = 3
            allocate(x(2,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ 0d0, 0d0 /)
            x(:,2) = (/ 1d0, 0d0 /)
            x(:,3) = (/ 0d0, 1d0 /)
            w(1) = 1d0/3d0
            w(2) = 1d0/3d0
            w(3) = 1d0/3d0
        case (TrianglePt1)
            N_pts = 1
            allocate(x(2,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ 1d0/3d0, 1d0/3d0 /)
            w(1) = 1d0/2d0
        case (TrianglePt4)
            N_pts = 4
            allocate(x(2,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ (1d0/sqrt(3d0)+1d0)/2d0,(1d0-1d0/sqrt(3d0))*(1d0+1d0/sqrt(3d0))/4d0 /)
            x(:,2) = (/ (1d0/sqrt(3d0)+1d0)/2d0,(1d0-1d0/sqrt(3d0))*(1d0-1d0/sqrt(3d0))/4d0 /)
            x(:,3) = (/ (1d0-1d0/sqrt(3d0))/2d0,(1d0+1d0/sqrt(3d0))*(1d0+1d0/sqrt(3d0))/4d0 /)
            x(:,4) = (/ (1d0-1d0/sqrt(3d0))/2d0,(1d0+1d0/sqrt(3d0))*(1d0-1d0/sqrt(3d0))/4d0 /)
            w(1) = (1d0-1d0/sqrt(3d0))/8d0
            w(2) = (1d0-1d0/sqrt(3d0))/8d0
            w(3) = (1d0+1d0/sqrt(3d0))/8d0
            w(4) = (1d0+1d0/sqrt(3d0))/8d0
        case (TrianglePt9)
            N_pts = 9
            allocate(x(2,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ 1d0/2,1d0/4 /)
            x(:,2) = (/ (1d0+sqrt(3d0/5d0))/2d0,(1d0-sqrt(3d0/5d0))*(1d0+sqrt(3d0/5d0))/4d0 /)
            x(:,3) = (/ (1d0+sqrt(3d0/5d0))/2d0,(1d0-sqrt(3d0/5d0))*(1d0-sqrt(3d0/5d0))/4d0 /)
            x(:,4) = (/ (1d0-sqrt(3d0/5d0))/2d0,(1d0+sqrt(3d0/5d0))*(1d0+sqrt(3d0/5d0))/4d0 /)
            x(:,5) = (/ (1d0-sqrt(3d0/5d0))/2d0,(1d0+sqrt(3d0/5d0))*(1d0-sqrt(3d0/5d0))/4d0 /)
            x(:,6) = (/ 1d0/2d0,1d0*(1d0+sqrt(3d0/5d0))/4d0 /)
            x(:,7) = (/ 1d0/2d0,1d0*(1d0-sqrt(3d0/5d0))/4d0 /)
            x(:,8) = (/ (1d0+sqrt(3d0/5d0))/2d0,(1d0-sqrt(3d0/5d0))/4d0 /)
            x(:,9) = (/ (1d0-sqrt(3d0/5d0))/2d0,(1d0+sqrt(3d0/5d0))/4d0 /)
            w(1) = 64d0/81d0/8d0
            w(2) = 100d0/324d0*(1d0-sqrt(3d0/5d0))/8d0
            w(3) = 100d0/324d0*(1d0-sqrt(3d0/5d0))/8d0
            w(4) = 100d0/324d0*(1d0+sqrt(3d0/5d0))/8d0
            w(5) = 100d0/324d0*(1d0+sqrt(3d0/5d0))/8d0
            w(6) = 40d0/81d0/8d0
            w(7) = 40d0/81d0/8d0
            w(8) = 40d0/81d0*(1d0-sqrt(3d0/5d0))/8d0
            w(9) = 40d0/81d0*(1d0+sqrt(3d0/5d0))/8d0
        case default
            print *, "Gauss quadrature type not implemented"
            stop
        end select

    end subroutine getGaussQuadRefTriangle

    subroutine getGaussQuadAnyTriangle(vertices,Gauss_type, x, w)
        real(8),intent(in) :: vertices(2,3)
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        real(8), dimension(:,:), allocatable :: x_ref
        real(8), dimension(:), allocatable :: w_ref

        real(8) :: x1, x2, x3, y1, y2, y3, J
        integer :: N_pts
        
        x1 = vertices(1,1)
        x2 = vertices(1,2)
        x3 = vertices(1,3)
        y1 = vertices(2,1)
        y2 = vertices(2,2)
        y3 = vertices(2,3)
        J = abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))

        call getGaussQuadRefTriangle(Gauss_type, x_ref, w_ref)
        N_pts = size(w_ref)

        allocate(x(2,N_pts))
        allocate(w(N_pts))
        w = w_ref*J
        x(1,:) = x1 + (x2-x1) * x_ref(1,:) + (x3-x1) * x_ref(2,:)
        x(2,:) = y1 + (y2-y1) * x_ref(1,:) + (y3-y1) * x_ref(2,:)

    end subroutine getGaussQuadAnyTriangle

    subroutine getGaussQuadRefLine(Gauss_type, x, w)
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        integer :: N_pts

        select case (Gauss_type)
        case (LinePt1)
            N_pts = 1
            allocate(x(1,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ 0d0 /)
            w(1) = 2d0
        case (LinePt2)
            N_pts = 2
            allocate(x(1,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ -1d0/sqrt(3d0) /)
            x(:,2) = (/ 1d0/sqrt(3d0) /)
            w(1) = 1d0
            w(2) = 1d0
        case (LinePt3)
            N_pts = 3
            allocate(x(1,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ -sqrt(3d0/5d0) /)
            x(:,2) = (/ 0d0 /)
            x(:,3) = (/ sqrt(3d0/5d0) /)
            w(1) = 5d0/9d0
            w(2) = 8d0/9d0
            w(3) = 5d0/9d0
        case (LinePt4)
            N_pts = 4
            allocate(x(1,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ -sqrt(525d0-70d0*sqrt(30d0))/35d0 /)
            x(:,2) = (/ sqrt(525d0-70d0*sqrt(30d0))/35d0 /)
            x(:,3) = (/ -sqrt(525d0+70d0*sqrt(30d0))/35d0 /)
            x(:,4) = (/ sqrt(525d0+70d0*sqrt(30d0))/35d0 /)
            w(1) = (18d0+sqrt(30d0))/36d0
            w(2) = (18d0+sqrt(30d0))/36d0
            w(3) = (18d0-sqrt(30d0))/36d0
            w(4) = (18d0-sqrt(30d0))/36d0
        case default
            print *, "Gauss quadrature type not implemented"
            stop
        end select

    end subroutine getGaussQuadRefLine

    subroutine getGaussQuadAnyLine(vertices,Gauss_type, x, w)
        real(8),intent(in),dimension(2,2) :: vertices
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        real(8), dimension(:,:), allocatable :: x_ref
        real(8), dimension(:), allocatable :: w_ref
        real(8), dimension(2) :: pt1,pt2,mid
        real(8) :: arclength

        pt1 = vertices(:,1)
        pt2 = vertices(:,2)
        mid = (pt1+pt2)/2d0
        arclength = sqrt((pt2(1)-pt1(1))**2+(pt2(2)-pt1(2))**2)

        call getGaussQuadRefLine(Gauss_type, x_ref, w_ref)

        allocate(x(2,size(w_ref)))
        allocate(w(size(w_ref)))
        x(1,:) = mid(1) + x_ref(1,:) * (pt2(1)-pt1(1)) / 2d0
        x(2,:) = mid(2) + x_ref(1,:) * (pt2(2)-pt1(2)) / 2d0
        w = w_ref * arclength / 2d0


    end subroutine getGaussQuadAnyLine

    subroutine getGaussQuadRefQuad(Gauss_type, x, w)
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        integer :: N_pts, Gauss_type_line

        real(8), dimension(:,:), allocatable :: x_line
        real(8), dimension(:), allocatable :: w_line
        real(8), dimension(:,:), allocatable :: w_line_t
        real(8), dimension(:,:),allocatable :: w_line_2d

        real(8), dimension(:,:), allocatable :: one

        if (Gauss_type .eq. QuadNodeAverage) then
            N_pts = 4
            allocate(x(2,N_pts))
            allocate(w(N_pts))
            x(:,1) = (/ 0d0, 0d0 /)
            x(:,2) = (/ 1d0, 0d0 /)
            x(:,3) = (/ 1d0, 1d0 /)
            x(:,4) = (/ 0d0, 1d0 /)
            w(1) = 1d0/4d0
            w(2) = 1d0/4d0
            w(3) = 1d0/4d0
            w(4) = 1d0/4d0
            return
        end if
        
        Gauss_type_line = Gauss_type - 200

        call getGaussQuadRefLine(Gauss_type_line, x_line, w_line)
        w_line = w_line / 2d0
        x_line = x_line / 2d0 + 1d0/2d0

        N_pts = size(w_line)**2

        allocate(one(1,size(w_line)))
        one = 1d0

        allocate(x(2,N_pts))
        allocate(w(N_pts))

        x(1,:) = reshape(matmul(reshape(one,(/size(w_line), 1/)),x_line),(/N_pts/))
        x(2,:) = reshape(matmul(reshape(x_line,(/size(w_line), 1/)),one),(/N_pts/))

        allocate(w_line_t(size(w_line),1))
        w_line_t = reshape(w_line,(/size(w_line), 1/))
        allocate(w_line_2d(1,size(w_line)))
        w_line_2d(1,:) = w_line
        w = reshape(matmul(w_line_t,w_line_2d),(/N_pts/))
    end subroutine getGaussQuadRefQuad

    subroutine getGaussQuadAnyQuad(vertices,Gauss_type, x, w)
        real(8),intent(in) :: vertices(2,4)
        real(8), dimension(:,:), allocatable, intent(out) :: x
        real(8), dimension(:), allocatable, intent(out) :: w
        integer, intent(in) :: Gauss_type

        real(8), dimension(:,:), allocatable :: x_ref
        real(8), dimension(:), allocatable :: w_ref

        real(8) :: x1, x2, x3, x4, y1, y2, y3, y4
        real(8), dimension(:), allocatable :: J
        real(8) :: cx1, cx2, cx3, cy1, cy2, cy3
        integer :: N_pts
        
        x1 = vertices(1,1)
        x2 = vertices(1,2)
        x3 = vertices(1,3)
        x4 = vertices(1,4)
        y1 = vertices(2,1)
        y2 = vertices(2,2)
        y3 = vertices(2,3)
        y4 = vertices(2,4)

        call getGaussQuadRefQuad(Gauss_type, x_ref, w_ref)

        N_pts = size(w_ref)
        allocate(x(2,N_pts))
        allocate(w(N_pts))

        x(1,:) = x1*(x_ref(1,:)-1d0)*(x_ref(2,:)-1d0) &
                + x2*x_ref(1,:)*(1-x_ref(2,:)) &
                + x3*x_ref(1,:)*x_ref(2,:) &
                + x4*(1-x_ref(1,:))*x_ref(2,:)
        x(2,:) = y1*(x_ref(1,:)-1d0)*(x_ref(2,:)-1d0) &
                + y2*x_ref(1,:)*(1-x_ref(2,:)) &
                + y3*x_ref(1,:)*x_ref(2,:) &
                + y4*(1-x_ref(1,:))*x_ref(2,:)

        cx1 = x1-x2+x3-x4
        cx2 = x2-x1
        cx3 = x4-x1
        cy1 = y1-y2+y3-y4
        cy2 = y2-y1
        cy3 = y4-y1

        allocate(J(N_pts))
        J = (cx2*cy1-cx1*cy2)*x_ref(1,:) + (cx1*cy3-cx3*cy1)*x_ref(2,:) +cx2*cy3-cx3*cy2

        w = w_ref*abs(J)
        
    end subroutine getGaussQuadAnyQuad
    
end module quadrature