module module_twogrid
    use settings
    use module_smoother
    use solver_umfpack2
    use module_prolongation
    implicit none
    
contains

subroutine solver_multigrid(Mats, Trs, b, x, rtol, max_steps, smoother, n_pre, n_post, print_screen, writeunit)
    type(MATRIX_COLUMN), dimension(:), intent(in) :: Mats
    type(GridTransfer), dimension(:), intent(in) :: Trs
    type(VECTOR), intent(in) :: b
    type(VECTOR), intent(inout) :: x
    real(8), intent(in) :: rtol
    integer, intent(in) :: max_steps
    integer, intent(in) :: smoother
    integer, intent(in) :: n_pre, n_post
    logical, optional :: print_screen
    integer, intent(in), optional :: writeunit
    
    integer :: i_smooth
    integer :: n_level
    
    
    type(MATRIX_COLUMN), dimension(:), allocatable :: At
    
    ! residue
    type(VECTOR) :: r, r_coarse, e_coarse, e
    real(8) :: res, res_old, rate
    
    integer :: i_step, i
    
    n_level = size(Mats)
    
    ! initial residue
    r = b
    call AddMultMV(-1d0, Mats(n_level), x, r)
    res = sum(sqrt(r%data**2))/r%size
    
    
    if(.not. present(print_screen)) print_screen = .true.
    if(print_screen) then
        write(*,*) "Step = ", 0, "Residual = ", res
    end if
    
    if(present(writeunit)) then
        write(writeunit,*) 0, res
    end if
    
    do i_step = 1, max_steps
        
        call cycle()
    
        ! output residue
        res_old = res
        r = b
        call AddMultMV(-1d0, A, x, r)
        res = r%Norm()
        rate = res/res_old
        
        if(print_screen) then
            write(*,*) "Step = ", i_step, "Residual = ", res, "Rate = ", rate
        end if
        
        if(present(writeunit)) then
            write(writeunit,*) i_step, res
        end if
                
        if(res<rtol) exit
        
    end do
    
    
    
end subroutine solver_multigrid

recursive subroutine cycle(i_level, Mats, Mat_ts, Trs, b, x, smoother, n_pre, n_post)
    integer, intent(in) :: i_level
    type(MATRIX_COLUMN), dimension(:), intent(in) :: Mats
    type(MATRIX_COLUMN), dimension(:), intent(in) :: Mat_ts
    type(GridTransfer), dimension(:), intent(in) :: Trs
    type(VECTOR), intent(in) :: b
    type(VECTOR), intent(inout) :: x
    integer, intent(in) :: smoother
    integer, intent(in) :: n_pre, n_post
    
    integer :: i_smooth
    
    type(VECTOR) :: r, r_coarse, e_coarse, e
    
    if(i_level==1) then
        call SolverSolveUMFPACK2(Mats(i_level), b, x)
        return
    end if
    
    call r%Init(b%size)
    call e%Init(b%size)
    call r_coarse%Init(Mats(i_level-1)%N_row)
    call e_coarse%Init(Mats(i_level-1)%N_row)
    
    
    ! pre-smoothing
    do i_smooth = 1, n_pre
        if(smoother==SMOOTHER_JACOBI) then
            call Jacobi_smooth(Mats(i_level), b, x)
        elseif(smoother==SMOOTHER_GS) then
            call GS_smooth(Mat_ts(i_level), b, x)
        else
            write(*,*) "Unknown smoother!"
            stop
        end if
    end do
    
    ! compute residue
    r = b
    call AddMultMV(-1d0, Mats(i_level), x, r)
    
    ! restrict to coarse grid
    call Trs(i_level-1)%FineToCoarse(r, r_coarse)
    
    ! recursive call
    call cycle(i_level-1, Mats, Mat_ts, Trs, r_coarse, e_coarse, smoother, n_pre, n_post)
    
    call Trs(i_level-1)%CoarseToFine(e_coarse, e)
                
    ! update x
    call x%AddVector(e)
    
    ! post-smoothing
    do i_smooth = 1, n_post
        if(smoother==SMOOTHER_JACOBI) then
            call Jacobi_smooth(Mats(i_level), b, x)
        elseif(smoother==SMOOTHER_GS) then
            call GS_smooth(Mat_ts(i_level), b, x)
        else
            write(*,*) "Unknown smoother!"
            stop
        end if
    end do
    
end subroutine cycle

    
end module module_twogrid


