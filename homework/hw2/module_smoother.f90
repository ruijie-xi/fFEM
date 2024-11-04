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
    type(VECTOR), intent(in) :: b
    type(VECTOR), intent(inout) :: x
    real(8), intent(in) :: rtol
    integer, intent(in) :: max_steps
    integer, intent(in) :: smoother
    integer, intent(in), optional :: writeunit
    
    type(MATRIX_COLUMN) :: At
    
    ! residue
    type(VECTOR) :: r
    real(8) :: res, res_old, rate
    
    integer :: i_step
    
    call assert(A%N_row==A%N_col, "A must be square matrix!")
    call assert(A%N_col==b%size, "A must be compatible with b!")
    call assert(A%N_col==x%size, "A must be compatible with x!")
    call assert(A%N_col>0, "A must be non-empty matrix!")
    
    if(smoother==SMOOTHER_GS) then
        At = MatrixColumnTranspose(A)
    end if
        
    ! initial residue
    call r%Init(b%size)
    r = b
    call AddMultMV(-1d0, A, x, r)
    res = r%Norm()
    
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
        res = r%Norm()
        
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
    type(VECTOR), intent(in) :: b
    type(VECTOR), intent(inout) :: x
    
    type(VECTOR) :: x_old, diag
    
    integer :: i_nz, i_row
    integer :: row, col
    real(8), parameter :: omega = 1d0
    
    call diag%Init(A%N_col)
    call x_old%Init(x%size)
    
    x_old = x
    
    ! Jacobi smoothing
    x = b
    i_nz = 0
    do col = 1, A%N_col
        do i_nz = A%col_ptr(col), A%col_ptr(col+1)-1
            row = A%row_idx(i_nz)
            x%data(row) = x%data(row) - A%val(i_nz)*x_old%data(col)
            if (row==col) then
                diag%data(col) = A%val(i_nz)
            end if
        end do
        if (abs(diag%data(col))<1d-8) then
            write(*,*) "Diagonal element of column ", col, " is too small!"
            stop
        end if
    end do

    do i_row = 1, A%N_row
        x%data(i_row) = x%data(i_row)/diag%data(i_row)
    end do
    
    x%data = x_old%data + omega*x%data
    
end subroutine Jacobi_smooth

! Gauss-Seidel smoothing
! Pass the transpose of A
subroutine GS_smooth(At, b, x)
    type(MATRIX_COLUMN), intent(in) :: At
    type(VECTOR), intent(in) :: b
    type(VECTOR), intent(inout) :: x 
    type(VECTOR) :: x_old, diag
    
    
    integer :: i_nz
    integer :: row, col
    
    call x_old%Init(x%size)
    x_old = x
    
    call diag%Init(At%N_col)
    diag%data = 0d0
    
    ! Gauss-Seidel smoothing
    i_nz = 0
    do col = 1, At%N_col
        x%data(col) = b%data(col)
        do i_nz = At%col_ptr(col), At%col_ptr(col+1)-1
            row = At%row_idx(i_nz)
            if (row==col) then
                diag%data(col) = At%val(i_nz)
            else if (row<col) then
                x%data(col) = x%data(col) - At%val(i_nz)*x%data(row)
            else
                x%data(col) = x%data(col) - At%val(i_nz)*x_old%data(row)
            end if
        end do
        if (abs(diag%data(col))<1d-8) then
            write(*,*) "Diagonal element of column ", col, " is too small!"
            stop
        end if
        x%data(col) = x%data(col)/diag%data(col)
    end do
        
end subroutine GS_smooth
    
end module module_smoother