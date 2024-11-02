module module_twogrid
    use settings
    use module_smoother
    use solver_umfpack2
    implicit none
    
    abstract interface
        subroutine Operator(x, y)
            real(8), intent(in), dimension(:) :: x
            real(8), intent(out), dimension(:), allocatable :: y
        end subroutine
    end interface
    
contains

subroutine solver_twogrid(A, Ac, OpR, OpP, b, x, rtol, max_steps, smoother, n_pre, n_post, print_screen, writeunit)
    type(MATRIX_COLUMN), intent(in) :: A
    type(MATRIX_COLUMN), intent(in) :: Ac
    procedure(Operator) :: OpR, OpP
    real(8), dimension(:), intent(in) :: b
    real(8), dimension(:), intent(inout) :: x
    real(8), intent(in) :: rtol
    integer, intent(in) :: max_steps
    integer, intent(in) :: smoother
    integer, intent(in) :: n_pre, n_post
    logical, optional :: print_screen
    integer, intent(in), optional :: writeunit
    
    integer :: i_smooth
    
    
    type(MATRIX_COLUMN) :: At
    
    ! residue
    real(8), dimension(:), allocatable :: r, r_coarse, e_coarse, e
    real(8) :: res, res_old, rate
    
    integer :: i_step
    
    call assert(A%N_row==A%N_col, "A must be square matrix!")
    call assert(A%N_col==size(b), "A must be compatible with b!")
    call assert(A%N_col==size(x), "A must be compatible with x!")
    call assert(A%N_col>0, "A must be non-empty matrix!")
    
    if(smoother==SMOOTHER_GS) then
        At = MatrixColumnTranspose(A)
    end if
    
    allocate(r_coarse(Ac%N_row), e_coarse(Ac%N_row))
    r_coarse = 0d0
    e_coarse = 0d0
        
    ! initial residue
    allocate(r(size(b)), e(size(b)))
    r = b
    e = 0d0
    call AddMultMV(-1d0, A, x, r)
    res = sum(sqrt(r**2))/size(r)
    
    
    if(.not. present(print_screen)) print_screen = .true.
    if(print_screen) then
        write(*,*) "Step = ", 0, "Residual = ", res
    end if
    
    if(present(writeunit)) then
        write(writeunit,*) 0, res
    end if
    
    do i_step = 1, max_steps
        
        ! pre-smoothing
        do i_smooth = 1, n_pre
            if(smoother==SMOOTHER_JACOBI) then
                call Jacobi_smooth(A, b, x)
            elseif(smoother==SMOOTHER_GS) then
                call GS_smooth(At, b, x)
            else
                write(*,*) "Unknown smoother!"
                stop
            end if
        end do
        
        ! compute residue
        r = b
        call AddMultMV(-1d0, A, x, r)
        res = sum(sqrt(r**2))/size(r)
        if(print_screen)write(*,*) "Pre-smoothing: Step = ", n_pre, "Residual = ", res
        
        ! restrict to coarse grid
        call OpR(r, r_coarse)
        
        ! solve on coarse grid
        call SolverSolveUMFPACK2(Ac, r_coarse, e_coarse)
        
        ! prolongate to fine grid
        call OpP(e_coarse, e)
                
        ! update x
        x = x + e
        
        r = b
        call AddMultMV(-1d0, A, x, r)
        res = sum(sqrt(r**2))/size(r)
        if(print_screen)write(*,*) "Correction: Residual = ", res
        
        ! post-smoothing
        do i_smooth = 1, n_post
            if(smoother==SMOOTHER_JACOBI) then
                call Jacobi_smooth(A, b, x)
            elseif(smoother==SMOOTHER_GS) then
                call GS_smooth(At, b, x)
            else
                write(*,*) "Unknown smoother!"
                stop
            end if
        end do
    
        ! output residue
        res_old = res
        r = b
        call AddMultMV(-1d0, A, x, r)
        res = sum(sqrt(r**2))/size(r)
        rate = res/res_old
        
        if(print_screen) then
            write(*,*) "Step = ", i_step, "Residual = ", res, "Rate = ", rate
        end if
        
        if(present(writeunit)) then
            write(writeunit,*) i_step, res
        end if
                
        if(res<rtol) exit
        
    end do
    
    
    
end subroutine solver_twogrid

    
end module module_twogrid


