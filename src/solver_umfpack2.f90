module solver_umfpack2
    use settings
    use matvec
    implicit none
    
interface SolverSolveUMFPACK2
    module procedure SolverSolveUMFPACK2_Column
    module procedure SolverSolveUMFPACK2_Triplet
end interface SolverSolveUMFPACK2

private :: SolverSolveUMFPACK2_Column, SolverSolveUMFPACK2_Triplet
    
contains

subroutine SolverSolveUMFPACK2_Column(A,b,x)
    type(MATRIX_COLUMN) :: A
    type(VECTOR) :: b
    type(VECTOR) :: x
    
    integer :: i_nz, i_col

    real(8), dimension(:), allocatable :: w

    integer :: icntl(20),keep(20),info(40)
    real(8) :: cntl(10),rinfo(20)

    real(8), dimension(:), allocatable :: value
    integer,dimension(:), allocatable :: index
    integer :: lvalue,lindex

    call ums2in(icntl, cntl, keep)

    lvalue = 100*A%N_nz
    lindex = 100*A%N_nz
    
    allocate(index(lindex),value(lvalue))
        
    do i_col = 1, A%N_col
        do i_nz = A%col_ptr(i_col), A%col_ptr(i_col+1)-1
            index(i_nz) = A%row_idx(i_nz)
            index(A%N_nz+i_nz) = i_col
            value(i_nz) = A%val(i_nz)
        end do
    end do

    call ums2fa(A%N_row, A%N_nz, 0, .false., lvalue, lindex, value, index, &
        keep, cntl, icntl, info, rinfo)
    
    allocate(w(4*A%N_row))
    call ums2so(A%N_row, 0, .false., lvalue, lindex, value, index, &
        keep, b%data, x%data, w, cntl, icntl, info, rinfo)
    

end subroutine SolverSolveUMFPACK2_Column

subroutine SolverSolveUMFPACK2_Triplet(A,b,x)
    type(MATRIX_TRIPLET) :: A
    type(vector) :: b
    type(vector) :: x
    
    real(8), dimension(:), allocatable :: w

    integer :: icntl(20),keep(20),info(40)
    real(8) :: cntl(10),rinfo(20)

    real(8), dimension(:), allocatable :: value
    integer,dimension(:), allocatable :: index
    integer :: lvalue,lindex
    
    call VectorReset(x, b%size)

    call ums2in(icntl, cntl, keep)
    icntl(3) = 1

    lvalue = 100*A%actual_nnz
    lindex = 100*A%actual_nnz
    allocate(index(lindex),value(lvalue))

    index(1:A%actual_nnz) = A%row_idx(1:A%actual_nnz)
    index(A%actual_nnz+1:2*A%actual_nnz) = A%col_idx(1:A%actual_nnz)
    value(1:A%actual_nnz) = A%val(1:A%actual_nnz)

    call ums2fa(A%N_row, A%actual_nnz, 0, .false., lvalue, lindex, value, index, &
        keep, cntl, icntl, info, rinfo)
    
    allocate(w(4*A%N_row))
    call ums2so(A%N_row, 0, .false., lvalue, lindex, value, index, &
        keep, b%data, x%data, w, cntl, icntl, info, rinfo)

end subroutine SolverSolveUMFPACK2_Triplet




    
end module solver_umfpack2