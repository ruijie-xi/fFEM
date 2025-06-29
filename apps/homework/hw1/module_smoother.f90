module module_smoother
    use settings
    use tools
    use matvec 
    implicit none
    
contains

subroutine Jacobi_smooth(A, b, x0, x)
    type(MATRIX_COLUMN), intent(in) :: A 
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(in) :: x0
    real(8), dimension(:), intent(out), allocatable :: x 
    
    integer :: i_nz, i_row
    integer :: row, col
    real(8), dimension(:), allocatable :: diag
    real(8), parameter :: omega = 2d0/3d0
    
    call assert(A%N_row==A%N_col, "A must be square matrix!")
    call assert(A%N_col==size(b), "A must be compatible with b!")
    call assert(A%N_col==size(x0), "A must be compatible with x!")
    call assert(A%N_col>0, "A must be non-empty matrix!")
    
    allocate(diag(A%N_col))
    allocate(x(A%N_col))
    diag = 0.d0
    x = b
    
    ! Jacobi smoothing
    i_nz = 0
    do col = 1, A%N_col
        do i_nz = A%col_ptr(col), A%col_ptr(col+1)-1
            row = A%row_idx(i_nz)
            
            x(row) = x(row) - A%val(i_nz)*x0(col)
            if (row==col) then
                diag(col) = A%val(i_nz)
            end if
        end do
        if (abs(diag(col))<1d-8) then
            write(*,*) "Diagonal element of column ", col, " is too small!"
            stop
        end if
    end do
    
    do i_row = 1, A%N_row
        x(i_row) = x(i_row)/diag(i_row)
    end do
    
    x = x0 + omega*x
    
    deallocate(diag)
    
end subroutine Jacobi_smooth


subroutine GS_smooth(A, b, x0, x)
    type(MATRIX_COLUMN), intent(in) :: A 
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(in) :: x0
    real(8), dimension(:), intent(out), allocatable :: x 
    
    integer :: i_nz
    integer :: row, col
    real(8), dimension(:), allocatable :: diag
    
    type(MATRIX_COLUMN) :: A_transpose
    
    call assert(A%N_row==A%N_col, "A must be square matrix!")
    call assert(A%N_col==size(b), "A must be compatible with b!")
    call assert(A%N_col==size(x0), "A must be compatible with x!")
    call assert(A%N_col>0, "A must be non-empty matrix!")
    
    ! transpose A
    A_transpose = MatrixColumnTranspose(A)
    
    allocate(diag(A%N_col))
    allocate(x(A%N_col))
    diag = 0.d0
    x = 0.0d0
    
    ! Gauss-Seidel smoothing
    i_nz = 0
    do col = 1, A_transpose%N_col
        x(col) = b(col)
        do i_nz = A_transpose%col_ptr(col), A_transpose%col_ptr(col+1)-1
            row = A_transpose%row_idx(i_nz)
            if (row==col) then
                diag(col) = A_transpose%val(i_nz)
            else if (row<col) then
                x(col) = x(col) - A_transpose%val(i_nz)*x(row)
            else
                x(col) = x(col) - A_transpose%val(i_nz)*x0(row)
            end if
        end do
        if (abs(diag(col))<1d-8) then
            write(*,*) "Diagonal element of column ", col, " is too small!"
            stop
        end if
        x(col) = x(col)/diag(col)
    end do
    
    deallocate(diag)
    
end subroutine GS_smooth
    
end module module_smoother