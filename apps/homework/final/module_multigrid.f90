module module_multigrid
    use settings
    use module_smoother
    use solver_umfpack2
    use module_prolongation
    use visualize
    implicit none
    
contains

subroutine solver_multigrid(mesh_list, fes_list, Mats, Trs, b, x, rtol, max_steps, smoother, n_pre, n_post, print_screen, writeunit)
    type(MATRIX_COLUMN), dimension(:), intent(in) :: Mats
    type(MESH2D), dimension(:), intent(in) :: mesh_list
    type(FESPACE), dimension(:), intent(in) :: fes_list
    type(GridTransfer), dimension(:), intent(in) :: Trs
    type(VECTOR), intent(in) :: b
    type(VECTOR), intent(inout) :: x
    real(8), intent(in) :: rtol
    integer, intent(in) :: max_steps
    integer, intent(in) :: smoother
    integer, intent(in) :: n_pre, n_post
    logical, optional :: print_screen
    integer, intent(in), optional :: writeunit
    
    integer :: n_level
    
    type(MATRIX_COLUMN), dimension(:), allocatable :: Mat_ts
    
    ! residue
    type(VECTOR) :: r
    real(8) :: res, res_old, rate, res0
    
    integer :: i_step, i
    
    n_level = size(Mats)
    
    ! initial residue
    r = b
    call AddMultMV(-1d0, Mats(n_level), x, 1d0, r)
    res = r%Norm()
    res0 = res
    
    if(.not. present(print_screen)) print_screen = .true.
    if(print_screen) then
        write(*,*) "Step = ", 0, "Residual = ", res
    end if
    
    if(present(writeunit)) then
        write(writeunit,*) 0, res
    end if
    
    if(smoother==SMOOTHER_GS) then 
        allocate(Mat_ts(n_level))
        do i = 1, n_level
            Mat_ts(i) = MatrixColumnTranspose(Mats(i))
        end do
    end if
    
    do i_step = 1, max_steps
        
        call cycle(n_level, Mats, Mat_ts, mesh_list, fes_list, Trs, b, x, smoother, n_pre, n_post)
    
        ! output residue
        res_old = res
        r = b
        call AddMultMV(-1d0, Mats(n_level), x, 1d0, r)
        res = r%Norm()
        rate = res/res_old
        
        if(print_screen) then
            write(*,*) "Step = ", i_step, "Residual = ", res, "Rate = ", rate
        end if
        
        if(present(writeunit)) then
            write(writeunit,*) i_step, res
        end if
                
        if(res/res0<rtol) exit
        
    end do
end subroutine solver_multigrid

recursive subroutine cycle(i_level, Mats, Mat_ts, mesh_list, &
fes_list, Trs, b, x, smoother, n_pre, n_post)
    integer, intent(in) :: i_level
    type(MATRIX_COLUMN), dimension(:), intent(in) :: Mats
    type(MATRIX_COLUMN), dimension(:), intent(in) :: Mat_ts
    type(MESH2D), dimension(:), intent(in) :: mesh_list
    type(FESPACE), dimension(:), intent(in) :: fes_list
    type(GridTransfer), dimension(:), intent(in) :: Trs
    type(VECTOR), intent(in) :: b
    type(VECTOR), intent(inout) :: x
    integer, intent(in) :: smoother
    integer, intent(in) :: n_pre, n_post
    
    integer, dimension(:, :), allocatable :: box_list
    
    integer :: i_smooth, i
    
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
        elseif(smoother==SMOOTHER_BOX) then
            box_list = reshape([(i, i=1, fes_list(i_level)%N_DOF)], [fes_list(i_level)%N_DOF/2, 2])
            call Box_smooth(Mats(i_level), b, x, box_list)
        else
            write(*,*) "Unknown smoother!"
            stop
        end if
    end do

    ! compute residue
    r = b
    call AddMultMV(-1d0, Mats(i_level), x, 1d0, r)
    
    write(*,*) "Level = ", i_level, "Residual = ", r%Norm()
            
    ! restrict to coarse grid
    call Trs(i_level-1)%FineToCoarse(r, r_coarse)
    
    write(*,*) "Level = ", i_level, "Residual coarse = ", r_coarse%Norm()
            
    ! recursive call
    call cycle(i_level-1, Mats, Mat_ts, mesh_list, fes_list, Trs, r_coarse, e_coarse, smoother, n_pre, n_post)
    
    write(*,*) "Level = ", i_level, "Error coarse = ", e_coarse%Norm()
    
    call Trs(i_level-1)%CoarseToFine(e_coarse, e)
    
    write(*,*) "Level = ", i_level, "Error = ", e%Norm()
                
    ! update x
    call x%AddVector(e)
    
    ! post-smoothing
    do i_smooth = 1, n_post
        if(smoother==SMOOTHER_JACOBI) then
            call Jacobi_smooth(Mats(i_level), b, x)
        elseif(smoother==SMOOTHER_GS) then
            call GS_smooth(Mat_ts(i_level), b, x)
        elseif(smoother==SMOOTHER_BOX) then
            box_list = reshape([(i, i=1, fes_list(i_level)%N_DOF)], [fes_list(i_level)%N_DOF/2, 2])
            call Box_smooth(Mats(i_level), b, x, box_list)
        else
            write(*,*) "Unknown smoother!"
            stop
        end if
    end do
    
end subroutine cycle

    
end module module_multigrid

