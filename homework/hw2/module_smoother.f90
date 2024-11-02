module module_smoother
    use settings
    use tools
    use matvec 
    implicit none
    
contains

subroutine Jacobi_smooth(A, b, x, rtol, max_steps, writeunit)
    type(MATRIX_COLUMN), intent(in) :: A 
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(inout) :: x
    real(8), intent(in) :: rtol
    integer, intent(in) :: max_steps
    integer, intent(in) :: writeunit
    
    integer :: i_nz, i_row
    integer :: row, col
    real(8), dimension(:), allocatable :: diag
    real(8), parameter :: omega = 2d0/3d0
    
    ! residue
    real(8), dimension(:), allocatable :: r,x_old
    real(8) :: res, res_old, rate
    
    integer :: i_step
    
    call assert(A%N_row==A%N_col, "A must be square matrix!")
    call assert(A%N_col==size(b), "A must be compatible with b!")
    call assert(A%N_col==size(x), "A must be compatible with x!")
    call assert(A%N_col>0, "A must be non-empty matrix!")
    
    allocate(diag(A%N_col))
    diag = 0.d0
    x = 0d0
    
    ! initial residue
    allocate(r(size(b)))
    r = b
    call AddMultMV(-1d0, A, x, r)
    res = maxval(abs(r))
    write(*,*) "Step = ", 0, "Residual = ", res
    write(writeunit,*) 0, res
    
    allocate(x_old(size(x)))
    x_old = x
    
    do i_step = 1, max_steps
        
        ! Jacobi smoothing
        x = b
        i_nz = 0
        do col = 1, A%N_col
            do i_nz = A%col_ptr(col), A%col_ptr(col+1)-1
                row = A%row_idx(i_nz)
                
                x(row) = x(row) - A%val(i_nz)*x_old(col)
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
        
        x = x_old + omega*x
        
        
        ! compute residue
        res_old = res
        r = b
        call AddMultMV(-1d0, A, x, r)
        res = maxval(abs(r))
        
        ! output residue
        rate = res/res_old
        write(*,*) "Step = ", i_step, "Residual = ", res, "Rate = ", rate
        write(writeunit,*) i_step, res
        
        ! save old solution
        x_old = x
        
        if(res<rtol) exit
        
    end do
    
end subroutine Jacobi_smooth


subroutine GS_smooth(A, b, x, rtol, max_steps, writeunit)
    type(MATRIX_COLUMN), intent(in) :: A 
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(inout) :: x 
    real(8), intent(in) :: rtol
    integer, intent(in) :: max_steps
    integer, intent(in) :: writeunit
    
    integer :: i_nz
    integer :: row, col
    real(8), dimension(:), allocatable :: diag
    
    type(MATRIX_COLUMN) :: A_transpose
    
    ! residue
    real(8), dimension(:), allocatable :: r, x_old
    real(8) :: res, res_old, rate
    integer :: i_step
    
    call assert(A%N_row==A%N_col, "A must be square matrix!")
    call assert(A%N_col==size(b), "A must be compatible with b!")
    call assert(A%N_col==size(x), "A must be compatible with x!")
    call assert(A%N_col>0, "A must be non-empty matrix!")
    
    ! transpose A
    A_transpose = MatrixColumnTranspose(A)
    
    allocate(diag(A%N_col))
    diag = 0.d0
    x = 0.0d0
    
    ! initial residue
    allocate(r(size(b)))
    r = b
    call AddMultMV(-1d0, A, x, r)
    res = maxval(abs(r))
    write(*,*) "Step = ", 0, "Residual = ", res
    write(writeunit,*) 0, res
    
    allocate(x_old(size(x)))
    x_old = x
    
    do i_step = 1, max_steps
    
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
                    x(col) = x(col) - A_transpose%val(i_nz)*x_old(row)
                end if
            end do
            if (abs(diag(col))<1d-8) then
                write(*,*) "Diagonal element of column ", col, " is too small!"
                stop
            end if
            x(col) = x(col)/diag(col)
        end do
        
        ! compute residue
        
        res_old = res
        r = b
        call AddMultMV(-1d0, A, x, r)
        res = maxval(abs(r))
        
        ! output residue
        rate = res/res_old
        write(*,*) "Step = ", i_step, "Residual = ", res, "Rate = ", rate
        write(writeunit,*) i_step, res
        
        ! save old solution
        x_old = x
        
        if(res<rtol) exit
        
    end do
        
end subroutine GS_smooth
    
end module module_smoother