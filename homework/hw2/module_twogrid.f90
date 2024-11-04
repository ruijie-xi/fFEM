module module_twogrid
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
    real(8) :: res, res_old, rate
    
    integer :: i_step, i
    
    n_level = size(Mats)
    
    ! initial residue
    r = b
    call AddMultMV(-1d0, Mats(n_level), x, r)
    res = r%Norm()    
    
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
        call AddMultMV(-1d0, Mats(n_level), x, r)
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
    
    integer :: i_smooth
    
    type(VECTOR) :: r, r_coarse, e_coarse, e
    
    if(i_level==1) then
        do i_smooth = 1, 100
            call Jacobi_smooth(Mats(i_level), b, x)
        end do
        ! call SolverSolveUMFPACK2(Mats(i_level), b, x)
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
    
    call PlotFunction(x, mesh_list(i_level), fes_list(i_level), "pre_solution.vtk")
    call PlotFunction(r, mesh_list(i_level), fes_list(i_level), "pre_residual.vtk")
        
    ! restrict to coarse grid
    call Trs(i_level-1)%FineToCoarse(r, r_coarse)
            
    ! recursive call
    call cycle(i_level-1, Mats, Mat_ts, mesh_list, fes_list, Trs, r_coarse, e_coarse, smoother, n_pre, n_post)
    
    call Trs(i_level-1)%CoarseToFine(e_coarse, e)
                
    ! update x
    call x%AddVector(e)
    
    call PlotFunction(x, mesh_list(i_level), fes_list(i_level), "cgc_solution.vtk")
    
    
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


