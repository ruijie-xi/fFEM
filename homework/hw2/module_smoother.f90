module module_smoother
    use settings
    use tools
    use matvec 
    implicit none
    
    integer, parameter :: SMOOTHER_JACOBI = 1
    integer, parameter :: SMOOTHER_GS = 2
    
contains

subroutine Smooth(A, b, x, rtol, max_steps, smoother, writeunit)
    type(MATRIX_COLUMN), intent(in) :: A 
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(inout) :: x
    real(8), intent(in) :: rtol
    integer, intent(in) :: max_steps
    integer, intent(in) :: smoother
    integer, intent(in), optional :: writeunit
    
    type(MATRIX_COLUMN) :: At
    
    ! residue
    real(8), dimension(:), allocatable :: r
    real(8) :: res, res_old, rate
    
    integer :: i_step
    
    call assert(A%N_row==A%N_col, "A must be square matrix!")
    call assert(A%N_col==size(b), "A must be compatible with b!")
    call assert(A%N_col==size(x), "A must be compatible with x!")
    call assert(A%N_col>0, "A must be non-empty matrix!")
    
    if(smoother==SMOOTHER_GS) then
        At = MatrixColumnTranspose(A)
    end if
        
    ! initial residue
    allocate(r(size(b)))
    r = b
    call AddMultMV(-1d0, A, x, r)
    res = maxval(abs(r))
    
    if(present(writeunit)) then
        write(*,*) "Step = ", 0, "Residual = ", res
        write(writeunit,*) 0, res
    end if
    
    do i_step = 1, max_steps
        
        ! smoothing
        if(smoother==SMOOTHER_JACOBI) then
            call Jacobi_smooth(A, b, x)
        elseif(smoother==SMOOTHER_GS) then
            call GS_smooth(At, b, x)
        else
            write(*,*) "Unknown smoother!"
            stop
        end if
        
        ! compute residue
        res_old = res
        r = b
        call AddMultMV(-1d0, A, x, r)
        res = maxval(abs(r))
        
        ! output residue
        rate = res/res_old
        
        if(present(writeunit)) then
            write(*,*) "Step = ", i_step, "Residual = ", res, "Rate = ", rate
            write(writeunit,*) i_step, res
        end if
                
        if(res<rtol) exit
        
    end do
    
end subroutine Smooth

subroutine Jacobi_smooth(A, b, x)
    type(MATRIX_COLUMN), intent(in) :: A 
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(inout) :: x
    
    real(8), dimension(:), allocatable :: x_old
    
    integer :: i_nz, i_row
    integer :: row, col
    real(8), dimension(:), allocatable :: diag
    real(8), parameter :: omega = 2d0/3d0
    
    allocate(diag(A%N_col))
    diag = 0d0
    
    allocate(x_old(size(x)))
    x_old = x
    
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
    
end subroutine Jacobi_smooth

! Gauss-Seidel smoothing
! Pass the transpose of A
subroutine GS_smooth(At, b, x)
    type(MATRIX_COLUMN), intent(in) :: At
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(inout) :: x 
    real(8), dimension(:), allocatable :: x_old
    real(8), dimension(:), allocatable :: diag
    
    integer :: i_nz
    integer :: row, col
    
    allocate(x_old(size(x)))
    x_old = x
    
    allocate(diag(At%N_col))
    diag = 0d0
    
    ! Gauss-Seidel smoothing
    i_nz = 0
    do col = 1, At%N_col
        x(col) = b(col)
        do i_nz = At%col_ptr(col), At%col_ptr(col+1)-1
            row = At%row_idx(i_nz)
            if (row==col) then
                diag(col) = At%val(i_nz)
            else if (row<col) then
                x(col) = x(col) - At%val(i_nz)*x(row)
            else
                x(col) = x(col) - At%val(i_nz)*x_old(row)
            end if
        end do
        if (abs(diag(col))<1d-8) then
            write(*,*) "Diagonal element of column ", col, " is too small!"
            stop
        end if
        x(col) = x(col)/diag(col)
    end do
        
end subroutine GS_smooth
    
end module module_smoother