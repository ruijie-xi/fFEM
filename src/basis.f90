module basis
    use settings
    implicit none
    
contains

    subroutine BasisLocal2D(refpts, Th, Vh, i_elem, deriv_type, result)
        real(8), intent(in), dimension(:,:) :: refpts
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type
        real(8), intent(out), dimension(:,:,:),allocatable :: result

        integer :: num_pts

        num_pts = size(refpts, 2)
        allocate(result(Vh%dim, Vh%N_local_basis, num_pts))

        select case (Vh%basis_type)
        case (DOF_P1)
            call BasisLocalP1(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_Q1)
            call BasisLocalQ1(refpts, Th, Vh, i_elem, deriv_type, result)
        case default
            print *, 'BasisLocal2D: Unknown basis type'
            stop
        end select
    end subroutine BasisLocal2D

    subroutine BasisReference2D(refpts, basis_type, deriv_type, result)
        real(8), intent(in), dimension(:,:) :: refpts
        integer, intent(in) :: basis_type, deriv_type
        integer :: num_pts
        real(8), intent(out), dimension(:,:,:) :: result
        ! result is a 3d array: (dimension, basis_no, point_no)

        num_pts = size(refpts, 2)

        select case (basis_type)
            case (DOF_P1)
                call BasisReferenceP1(refpts, deriv_type, result(1,:,:))
            case default
                print *, 'BasisReference2D: Unknown basis type'
                stop
        end select

    end subroutine BasisReference2D

    subroutine BasisReferenceP1(refpts, deriv_type, result)
        real(8), dimension(3,3), parameter :: A = reshape((/ 1d0, 0d0, 0d0, -1d0, 1d0, 0d0, -1d0, 0d0, 1d0 /), (/3,3/))
        real(8), intent(in), dimension(:,:) :: refpts
        integer, intent(in) :: deriv_type
        real(8), dimension(:,:), allocatable :: b
        real(8), intent(out), dimension(:,:) :: result
        
        allocate(b(3, size(refpts, 2)))

        select case (deriv_type)
        case (DERIV_NONE)
            b(1,:) = 1d0
            b(2,:) = refpts(1,:)
            b(3,:) = refpts(2,:)
        case (DERIV_DX)
            b(1,:) = 0d0
            b(2,:) = 1d0
            b(3,:) = 0d0
        case (DERIV_DY)
            b(1,:) = 0d0
            b(2,:) = 0d0
            b(3,:) = 1d0
        case default
            print *, 'BasisP1: Unknown derivative type'
            stop
        end select

        result = matmul(A, b)
    end subroutine BasisReferenceP1

    subroutine BasisLocalP1(refpts, Th, Vh, i_elem, deriv_type, result)
        real(8), intent(in), dimension(:,:) :: refpts
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type
        real(8), intent(out), dimension(:,:,:), allocatable :: result

        real(8), dimension(:,:), allocatable :: basis_dxh, basis_dyh

        integer :: num_pts, i_dim

        real(8) :: x1,x2,x3,y1,y2,y3, detJ

        num_pts = size(refpts, 2)
        allocate(result(Vh%dim, Vh%N_local_basis, num_pts))

        select case (deriv_type)
        case (DERIV_NONE)
            call BasisReferenceP1(refpts, deriv_type, result(1,:,:))
        case (DERIV_DX)

            ! get coordinates and jacobian
            x1 = Th%NodeCoord(1, Th%ElemNodeConn(1, i_elem))
            x2 = Th%NodeCoord(1, Th%ElemNodeConn(2, i_elem))
            x3 = Th%NodeCoord(1, Th%ElemNodeConn(3, i_elem))
            y1 = Th%NodeCoord(2, Th%ElemNodeConn(1, i_elem))
            y2 = Th%NodeCoord(2, Th%ElemNodeConn(2, i_elem))
            y3 = Th%NodeCoord(2, Th%ElemNodeConn(3, i_elem))

            detJ = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

            allocate(basis_dxh(Vh%N_local_basis, num_pts))
            allocate(basis_dyh(Vh%N_local_basis, num_pts))
            call BasisReferenceP1(refpts, DERIV_DX, basis_dxh)
            call BasisReferenceP1(refpts, DERIV_DY, basis_dyh)

            result(1,:,:) = ((y3-y1)*basis_dxh + (y1-y2)*basis_dyh) / detJ
        
        case (DERIV_DY)
            ! get coordinates and jacobian
            x1 = Th%NodeCoord(1, Th%ElemNodeConn(1, i_elem))
            x2 = Th%NodeCoord(1, Th%ElemNodeConn(2, i_elem))
            x3 = Th%NodeCoord(1, Th%ElemNodeConn(3, i_elem))
            y1 = Th%NodeCoord(2, Th%ElemNodeConn(1, i_elem))
            y2 = Th%NodeCoord(2, Th%ElemNodeConn(2, i_elem))
            y3 = Th%NodeCoord(2, Th%ElemNodeConn(3, i_elem))
            detJ = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)

            allocate(basis_dxh(Vh%N_local_basis, num_pts))
            allocate(basis_dyh(Vh%N_local_basis, num_pts))
            call BasisReferenceP1(refpts, DERIV_DX, basis_dxh)
            call BasisReferenceP1(refpts, DERIV_DY, basis_dyh)

            result(1,:,:) = ((x1-x3)*basis_dxh + (x2-x1)*basis_dyh) / detJ

        case default
            print *, 'BasisLocalP1: Unknown derivative type'
            stop

        end select

        do i_dim = 1, Vh%dim
            result(i_dim,:,:) = result(1,:,:);
        end do
    end subroutine BasisLocalP1

    subroutine BasisReferenceQ1(refpts, deriv_type, result)
        real(8), dimension(4,4):: A
        real(8), intent(in), dimension(:,:) :: refpts
        integer, intent(in) :: deriv_type
        real(8), dimension(:,:), allocatable :: b
        real(8), intent(out), dimension(:,:) :: result

        A(:, 1) = (/ 1d0, 0d0, 0d0, 0d0 /)
        A(:, 2) = (/ -1d0, 1d0, 0d0, 0d0 /)
        A(:, 3) = (/ -1d0, 0d0, 0d0, 1d0 /)
        A(:, 4) = (/ 1d0, -1d0, 1d0, -1d0 /)
        
        allocate(b(4, size(refpts, 2)))

        select case (deriv_type)
        case (DERIV_NONE)
            b(1,:) = 1d0
            b(2,:) = refpts(1,:)
            b(3,:) = refpts(2,:)
            b(4,:) = refpts(1,:)*refpts(2,:)
        case (DERIV_DX)
            b(1,:) = 0d0
            b(2,:) = 1d0
            b(3,:) = 0d0
            b(4,:) = refpts(2,:)
        case (DERIV_DY)
            b(1,:) = 0d0
            b(2,:) = 0d0
            b(3,:) = 1d0
            b(4,:) = refpts(1,:)
        case default
            print *, 'BasisQ1: Unknown derivative type'
            stop
        end select

        result = matmul(A, b)
    end subroutine BasisReferenceQ1

    subroutine BasisLocalQ1(refpts, Th, Vh, i_elem, deriv_type, result)
        real(8), intent(in), dimension(:,:) :: refpts
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type
        real(8), intent(out), dimension(:,:,:), allocatable :: result

        real(8), dimension(2,4) :: vertices
        real(8), dimension(:,:), allocatable :: basis_dxh, basis_dyh
        real(8), dimension(:), allocatable :: J11,J12,J21,J22,J

        integer :: num_pts, i_dim

        num_pts = size(refpts, 2)
        allocate(result(Vh%dim, Vh%N_local_basis, num_pts))

        vertices = Th%NodeCoord(:, Th%ElemNodeConn(:, i_elem))
        call JacobianCoeffQuad(vertices, refpts, J11,J12,J21,J22)
        allocate(J(size(refpts,2)))
        J = J11*J22 - J12*J21

        select case (deriv_type)
        case (DERIV_NONE)
            call BasisReferenceQ1(refpts, deriv_type, result(1,:,:))
        case (DERIV_DX)
            allocate(basis_dxh(Vh%N_local_basis, num_pts))
            allocate(basis_dyh(Vh%N_local_basis, num_pts))
            call BasisReferenceQ1(refpts, DERIV_DX, basis_dxh)
            call BasisReferenceQ1(refpts, DERIV_DY, basis_dyh)
            result(1,:,:) = (spread(J22,1,size(basis_dxh,1))*basis_dxh &
                            - spread(J21,1,size(basis_dyh,1))*basis_dyh) &
                            /spread(J,1,size(basis_dxh,1))
        case (DERIV_DY)
            allocate(basis_dxh(Vh%N_local_basis, num_pts))
            allocate(basis_dyh(Vh%N_local_basis, num_pts))
            call BasisReferenceQ1(refpts, DERIV_DX, basis_dxh)
            call BasisReferenceQ1(refpts, DERIV_DY, basis_dyh)
            result(1,:,:) = (-spread(J12,1,size(basis_dxh,1))*basis_dxh &
                            + spread(J11,1,size(basis_dyh,1))*basis_dyh)&
                            /spread(J,1,size(basis_dxh,1))
        case default
            print *, 'BasisLocalQ1: Unknown derivative type'
            stop

        end select

        do i_dim = 1, Vh%dim
            result(i_dim,:,:) = result(1,:,:);
        end do
    end subroutine BasisLocalQ1


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
    
end module basis