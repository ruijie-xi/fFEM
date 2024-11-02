module solver_umfpack2
    use settings
    implicit none
contains

subroutine SolverSolveUMFPACK2(A,b,x)
    type(MATRIX_COLUMN) :: A
    real(8),dimension(:) :: b
    real(8),dimension(:) :: x
    
    integer :: i_nz, i_col

    real(8), dimension(:), allocatable :: w

    integer :: icntl(20),keep(20),info(40)
    real(8) :: cntl(10),rinfo(20)

    real(8), dimension(:), allocatable :: value
    integer,dimension(:), allocatable :: index
    integer :: lvalue,lindex

    call ums2in(icntl, cntl, keep)

    lvalue = 10*A%N_nz
    lindex = 10*A%N_nz
    
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
        keep, b, x, w, cntl, icntl, info, rinfo)
    

end subroutine SolverSolveUMFPACK2



    
end module solver_umfpack2