module solver_umfpack2
    use settings
    implicit none
contains


! TODO: dynamically allocate the size of value and index
subroutine SolverSolveUMFPACK2(A,b,x)
    type(MATRIX_TRIPLET) :: A
    real(8),dimension(:) :: b
    real(8),dimension(:) :: x

    real(8), dimension(:), allocatable :: w

    integer :: icntl(20),keep(20),info(40)
    real(8) :: cntl(10),rinfo(20)

    real(8), dimension(:), allocatable :: value
    integer,dimension(:), allocatable :: index
    integer :: lvalue,lindex

    call ums2in(icntl, cntl, keep)

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
        keep, b, x, w, cntl, icntl, info, rinfo)
    

end subroutine SolverSolveUMFPACK2



    
end module solver_umfpack2