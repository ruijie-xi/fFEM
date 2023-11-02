module settings
    implicit none

    ! Dimension : only support 2D
    integer, parameter :: DIM__ = 2

    ! Basis type
    integer, parameter :: DOF_P1 = 1
    ! integer, parameter :: DOF_P2 = 2
    integer, parameter :: DOF_Q1 = 101

    ! Derivative type
    integer, parameter :: DERIV_NONE = 0
    integer, parameter :: DERIV_DX = 1
    integer, parameter :: DERIV_DY = 2

    ! Norm type
    integer, parameter :: NORM_L2 = 0
    integer, parameter :: NORM_H1 = 1

    ! Gauss type
    integer, parameter :: LinePt1 = 1
    integer, parameter :: LinePt2 = 2
    integer, parameter :: LinePt3 = 3
    integer, parameter :: LinePt4 = 4
    integer, parameter :: TrianglePt1 = 101
    integer, parameter :: TrianglePt4 = 102
    integer, parameter :: TrianglePt9 = 103
    integer, parameter :: TriangleNodeAverage = 104
    integer, parameter :: QuadPt1 = 201
    integer, parameter :: QuadPt4 = 202
    integer, parameter :: QuadPt9 = 203
    integer, parameter :: QuadPt16 = 204
    integer, parameter :: QuadNodeAverage = 205

    interface
        subroutine func(x, f, deriv_type)
            real(8), intent(in), dimension(:) :: x
            real(8), dimension(:), intent(out), allocatable :: f
            integer, intent(in) :: deriv_type
        end subroutine
    end interface

    type :: mesh2D
        integer :: N_node, N_elem, N_edge, N_le
        integer, dimension(:,:), allocatable :: ElemNodeConn
        integer, dimension(:,:), allocatable :: EdgeNodeConn
        integer, dimension(:,:), allocatable :: ElemEdgeConn
        integer, dimension(:,:), allocatable :: EdgeElemConn
        real(8), dimension(:,:), allocatable :: NodeCoord
        integer, dimension(:), allocatable :: EdgeMarker
        real(8) :: hmax
    end type mesh2D

    type :: fespace
        integer :: dim
        integer :: N_local_basis
        integer :: N_DOF
        integer :: basis_type
        integer, dimension(:,:),allocatable :: ElemDOF
    end type
    
contains
    
end module settings