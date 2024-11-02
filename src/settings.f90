module settings
    implicit none

    ! Dimension : only support 2D
    integer, parameter :: DIM__ = 2

    ! Mesh type
    integer, parameter :: MESH_TRIANGLE = 3
    integer, parameter :: MESH_QUAD = 4

    ! Basis type
    integer, parameter :: DOF_P0 = 0
    integer, parameter :: DOF_P1 = 1
    integer, parameter :: DOF_P2 = 2
    integer, parameter :: DOF_DG1 = 11
    integer, parameter :: DOF_Q0 = 100
    integer, parameter :: DOF_Q1 = 101
    ! integer, parameter :: DOF_Q2 = 102
    integer, parameter :: DOF_QuadNedelec1 = 111
    integer, parameter :: DOF_QuadRT1 = 121
    integer, parameter :: DOF_QuadRT2 = 122

    ! Derivative type
    integer, parameter :: DERIV_NONE = 0
    integer, parameter :: DERIV_DX = 1
    integer, parameter :: DERIV_DY = 2
    integer, parameter :: DERIV_DXX = 3
    integer, parameter :: DERIV_DXY = 4
    integer, parameter :: DERIV_DYY = 5

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

    abstract interface
        subroutine func(x, f, deriv_type)
            real(8), intent(in), dimension(:) :: x
            real(8), dimension(:), intent(out), allocatable :: f
            integer, intent(in) :: deriv_type
        end subroutine
    end interface
    
    type :: element_list
        integer :: N_elem
        integer, dimension(:), allocatable :: elem_idx
    end type element_list

    type :: MESH2D
        integer :: mesh_type
        integer :: N_node, N_elem, N_edge, N_le, N_bdryedge, N_nodeinelem
        integer, dimension(:,:), allocatable :: ElemNodeConn
        integer, dimension(:,:), allocatable :: EdgeNodeConn
        integer, dimension(:,:), allocatable :: ElemEdgeConn
        integer, dimension(:,:), allocatable :: EdgeElemConn
        integer, dimension(:,:), allocatable :: EdgeIdxInElem
        real(8), dimension(:,:), allocatable :: NodeCoord
        integer, dimension(:), allocatable :: BdryEdge
        integer, dimension(:), allocatable :: BdryMarker
        integer, dimension(:), allocatable :: Edge2Bdry
        real(8) :: hmax
        real(8),dimension(2) :: xlim, ylim
        type(element_list), dimension(:,:), allocatable :: ElemInBoxes
    end type MESH2D

    type :: FESPACE
        integer :: dim
        integer :: N_local_basis
        integer :: N_DOF
        integer :: basis_type
        integer, dimension(:,:),allocatable :: ElemDOF
        integer :: isStack
    end type FESPACE

    ! triplet form of sparse matrix (COO format)
    ! allow duplicate nonzeros
    type :: MATRIX_TRIPLET
        integer :: N_row, N_col, N_nz
        integer, dimension(:), allocatable :: row_idx, col_idx
        real(8), dimension(:), allocatable :: val
        integer :: actual_nnz ! number of added non-zero elements
    end type MATRIX_TRIPLET

    ! column-oriented form of sparse matrix
    ! col_ptr(N_col+1): the indices of nonzero elements in the ith column is col_ptr(i)->col_ptr(i+1)-1
    ! row_idx(N_nz): the row indices of each nonzero elements
    ! val(N_nz): value of nonzero elements
    type :: MATRIX_COLUMN
        integer :: N_row, N_col, N_nz
        integer, dimension(:), allocatable :: col_ptr ! column pointer, N_col+1
        integer, dimension(:), allocatable :: row_idx ! row index of non-zero elements, N_nz
        real(8), dimension(:), allocatable :: val ! non-zero elements, N_nz
    end type MATRIX_COLUMN

    ! row-oriented form of sparse matrix
    ! row_ptr(N_row+1): the indices of nonzero elements in the ith row is row_ptr(i)->row_ptr(i+1)-1
    ! col_idx(N_nz): the column indices of each nonzero elements
    ! val(N_nz): value of nonzero elements
    type :: MATRIX_ROW
        integer :: N_row, N_col, N_nz
        integer, dimension(:), allocatable :: row_ptr ! row pointer, N_row+1
        integer, dimension(:), allocatable :: col_idx ! column index of non-zero elements, N_nz
        real(8), dimension(:), allocatable :: val ! non-zero elements, N_nz
    end type MATRIX_ROW

    ! some useful constants
    real(8), parameter :: m_pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286
    
contains
    
end module settings