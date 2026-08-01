module magnetic
    use settings
    use timer
    use fe, only: fespaceInit, Interpolate, BasisLocal2D
    use fe_utils
    use solver_umfpack2
    use assembler
    use tools
    use matvec
    use mesh
    use visualize
    implicit none

    private
    public :: magnetic_diffusion, magnetic_init, magnetic_lorentz, magnetic_update
    public :: magnetic_plot
    public :: magnetic_errorB, magnetic_errorE, magnetic_errorLorentz
    public :: magnetic_Bmax
    
    real(8) :: t_glb 
    common /global/ t_glb

    ! Gaussian quadrature type
    integer :: Gauss_type = QuadPt16, Gauss_type_bdry=LinePt4

    ! mesh and fespace
    type(MESH2D) :: Th
    type(FESPACE) :: Bh,Eh,Qh,Uh

    type(vector) :: B, E, sigma, E_int
    real(8) :: mu = 4d0*m_pi

    ! theta-scheme
    real(8) :: theta = 0.5d0
    
    real(8) :: outer_radius = 1d0
        
contains

subroutine magnetic_init(elems, nodes, B0_func, E0_func, sigma0_func, radius)
    implicit none
    integer, dimension(:,:), intent(inout) :: elems ! 单元编号
    real(8), dimension(:,:), intent(in) :: nodes ! 节点坐标
    procedure(func) :: B0_func, E0_func, sigma0_func ! 磁场，电场初始值，电导率函数
    real(8), intent(in) :: radius ! 边界半径

    integer, parameter :: DOF_type_Q = DOF_Q0
    integer, parameter :: DOF_type_E = DOF_Q2
    integer, parameter :: DOF_type_B = DOF_QuadRT2
    integer, parameter :: DOF_type_U = DOF_Q1

    real :: t_test

    ! generate mesh topological structure
    call MeshInit(elems,nodes,Th)
    outer_radius = radius
    call AddBdryMarker(Th, CircleBdry_func, 2, 1d-8)
    write(*,*) "Mesh generation done. "
    
    ! initialize fe spaces
    call timer_start()
    call fespaceInit(Qh, Th, DOF_type_Q, 1)
    call fespaceInit(Eh, Th, DOF_type_E, 1)
    call fespaceInit(Bh, Th, DOF_type_B, 2)
    call fespaceInit(Uh, Th, DOF_type_U, 2)
    call timer_end(t_test)
    write(*,*) "Initialize Done. Time taken = ", t_test

    ! initialize functions
    call Interpolate(B, B0_func, Th, Bh)
    call Interpolate(E, E0_func, Th, Eh)
    call Interpolate(E_int, E0_func, Th, Eh)
    call Interpolate(sigma, sigma0_func, Th, Qh)
    
end subroutine magnetic_init

subroutine magnetic_diffusion(to, dt, g_func, f_func, firststep)
    implicit none
    real(8), intent(in) :: to, dt ! 时间，时间步长
    procedure(func) :: g_func, f_func ! 边界函数，源项函数
    logical, intent(in) :: firststep ! 是否第一个时间步
    
    real(8) :: t_tmp ! temporary saving variable for t_glb
    ! functions
    type(vector) :: Bold, curlE

    ! matrices and vectors
    integer, dimension(:,:), allocatable :: assemble_info
    type(MATRIX_TRIPLET) :: matM,matM1,matTotal
    type(vector) :: vec_B,vec_g,vec_f,vecTotal
    
    ! save the old B
    call VectorCopy(Bold, B)
    
    print *, "Magnetic diffusion: advancing from ", to, " to ", to+dt

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! assemble matrices and vectors !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! assemble matrix M
    call MatrixTripletInit(matM, Eh%N_DOF, Eh%N_DOF, Th%N_elem*Eh%N_local_basis*Eh%N_local_basis)
    allocate(assemble_info(1,7))
    assemble_info(1,:) = [1, 1, DERIV_NONE, 1, DERIV_NONE, 1, DERIV_NONE]
    call AssembleMatrixElementFE(one_func, [1d0], Th, Eh, 0, Eh, 0, sigma, Qh,&
        assemble_info, Gauss_type, matM)
    deallocate(assemble_info)

    ! assemble matrix M1
    call MatrixTripletInit(matM1, Eh%N_DOF, Eh%N_DOF, 2*Th%N_elem*Eh%N_local_basis*Eh%N_local_basis)
    allocate(assemble_info(2,5))
    assemble_info(1,:) = [1, 1, DERIV_DX, 1, DERIV_DX]
    assemble_info(2,:) = [1, 1, DERIV_DY, 1, DERIV_DY]
    call AssembleMatrixElement(one_func, [1d0/mu,1d0/mu], Th, Eh, 0, Eh, 0,&
        assemble_info, Gauss_type, matM1)
    deallocate(assemble_info)

    ! assemble vector vec_B
    call VectorInit(vec_B, Eh%N_DOF)
    allocate(assemble_info(2,5))
    assemble_info(1,:) = [1, 1, DERIV_DY, 1, DERIV_NONE]
    assemble_info(2,:) = [1, 1, DERIV_DX, 2, DERIV_NONE]
    call AssembleVectorElementFE(one_func, [1d0/mu, -1d0/mu], Th, Eh, 0, Bold, Bh, &
    assemble_info, Gauss_type, vec_B)
    deallocate(assemble_info)

    ! assemble vector vec_g
    call VectorInit(vec_g, Eh%N_DOF)
    allocate(assemble_info(1,3))
    assemble_info(1,:) = [1, 1, DERIV_NONE]
    t_tmp = t_glb
    t_glb = to + 0.5d0*dt
    call AssembleVectorBdry(2, g_func, [1d0/mu], Th, Eh, 0, assemble_info, Gauss_type_bdry, vec_g)
    deallocate(assemble_info)
    t_glb = t_tmp

    ! assemble vector vec_f
    call VectorInit(vec_f, Eh%N_DOF)
    allocate(assemble_info(1,3))
    assemble_info(1,:) = [1, 1, DERIV_NONE]
    t_tmp = t_glb
    t_glb = to + 0.5d0*dt
    call AssembleVectorElement(f_func, [1d0], Th, Eh, 0, assemble_info, Gauss_type, vec_f)
    deallocate(assemble_info)
    t_glb = t_tmp

    ! form the linear system
    ! matTotal = matM + theta*dt*matM1
    ! vecTotal = vec_B + vec_g - vec_f
    call MatrixTripletAdd(1d0, matM, theta*dt, matM1, matTotal)
    call VectorInit(vecTotal, Eh%N_DOF)
    call VectorCopy(vecTotal, vec_B)
    call VectorAddVector(vecTotal, 1d0, vec_g, 1d0)
    call VectorAddVector(vecTotal, 1d0, vec_f, -1d0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! solve the linear system !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    block
        type(MATRIX_COLUMN) :: matTotal_col
        call MatrixTriplet2Column(matTotal, matTotal_col)
        call SolverSolveUMFPACK2(matTotal_col, vecTotal, E)
    end block
    call CurlOperator(E, curlE)
    
    call VectorAddVector(B, 1d0, curlE, -dt)
    
    ! extrapolate E_int
    ! call VectorCopy(E_int, E)
    ! if(abs(to)<1d-6) then
    !     call VectorAddVector(E_int, -1d0, E, 2d0)
    ! else 
    !     call VectorAddVector(E_int, -1d0/2d0, E, 3d0/2d0)
    ! end if
    call VectorAddVector(E_int, -1d0, E, 2d0)
    
    if (firststep) then
    call VectorAddVector(B, 0.5d0, Bold, 0.5d0)
    end if
    
end subroutine magnetic_diffusion

! scalar one
subroutine one_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 1d0
    end select
end subroutine one_func

! boundary function - circle
subroutine CircleBdry_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = abs(sqrt(x(1)**2+x(2)**2) - outer_radius)
    end select
end subroutine CircleBdry_func

subroutine CurlOperator(f,curlf)
    implicit none
    type(vector), intent(in) :: f
    type(vector), intent(out) :: curlf
    
    integer, dimension(:,:), allocatable :: assemble_info
    type(MATRIX_TRIPLET) :: M
    type(vector) :: rhs
    
    call VectorInit(curlf, Bh%N_DOF)
    call VectorInit(rhs, Bh%N_DOF)
    call MatrixTripletInit(M, Bh%N_DOF, Bh%N_DOF, Th%N_elem*Bh%N_local_basis*Bh%N_local_basis*2)
    
    allocate(assemble_info(2,5))
    assemble_info(1,:) = [1, 1, DERIV_NONE, 1, DERIV_NONE] 
    assemble_info(2,:) = [1, 2, DERIV_NONE, 2, DERIV_NONE]
    call AssembleMatrixElement(one_func, [1d0, 1d0], Th, Bh, 0, Bh, 0, assemble_info, Gauss_type, M)
    deallocate(assemble_info)
    
    allocate(assemble_info(2,5))
    assemble_info(1,:) = [1, 1, DERIV_NONE, 1, DERIV_DY]
    assemble_info(2,:) = [1, 2, DERIV_NONE, 1, DERIV_DX]
    call AssembleVectorElementFE(one_func, [1d0, -1d0], Th, Bh, 0, f, Eh, assemble_info, Gauss_type, rhs)
    deallocate(assemble_info)

    call SolverSolveUMFPACK2(M, rhs, curlf)
    
end subroutine CurlOperator

subroutine magnetic_errorB(B_func, L2err, relL2err)
    implicit none
    procedure(func) :: B_func
    real(8), intent(out) :: L2err, relL2err
    
    type(vector) :: zero
    real(8) :: L2norm_B

    call ComputeError(B_func, B, Th, Bh, NORM_L2, Gauss_type, L2err)
    
    call VectorInit(zero, Bh%N_DOF)
    call ComputeError(B_func, zero, Th, Bh, NORM_L2, Gauss_type, L2norm_B)
    relL2err = L2err / L2norm_B
    
end subroutine

subroutine magnetic_errorE(E_func, L2err)
    implicit none
    procedure(func) :: E_func
    real(8), intent(out) :: L2err

    call ComputeError(E_func, E_int, Th, Eh, NORM_L2, Gauss_type, L2err)
end subroutine

subroutine magnetic_errorLorentz(F, Lorentz_func, L2err)
    implicit none
    procedure(func) :: Lorentz_func
    type(vector), intent(in) :: F
    real(8), intent(out) :: L2err

    call ComputeError(Lorentz_func, F, Th, Uh, NORM_L2, Gauss_type, L2err)
end subroutine 

! JxB = \sigma E x B = (-sigma E By, sigma E Bx)
! patch recovery to Q1
subroutine magnetic_lorentz(F)
    implicit none
    type(vector), intent(out) :: F ! 输出的洛伦兹力

    real(8), dimension(:,:), allocatable :: val_sigma, val_E, val_B
    integer :: i_elem,i_dim,i_basis
    
    real(8), dimension(:,:), allocatable :: ref_pt_h
    real(8), dimension(:), allocatable :: w_h
    real(8), dimension(:,:), allocatable :: pts
    real(8), dimension(:), allocatable :: w
    
    real(8), dimension(:,:,:), allocatable :: basis_uh
    
    integer, dimension(:), allocatable :: idx_local_dof
    real(8), dimension(:,:), allocatable :: lorentz_value
    real(8), dimension(:,:), allocatable :: lorentz_to_basis

    call VectorInit(F, Uh%N_DOF)
    allocate(val_E(Eh%dim,Th%N_le), val_B(Bh%dim,Th%N_le), val_sigma(Qh%dim,Th%N_le))
    
    call getGaussRefElement(Gauss_type, ref_pt_h, w_h)
    allocate(lorentz_value(2,size(ref_pt_h,2)))
    allocate(lorentz_to_basis(2, Uh%N_local_basis))

    do i_elem = 1, Th%N_elem
    
        ! evaluate lorentz force at quadrature points
        call FEfunctionQuadValue(sigma, Th, Qh, i_elem, DERIV_NONE, Gauss_type, val_sigma)
        call FEfunctionQuadValue(E, Th, Eh, i_elem, DERIV_NONE, Gauss_type, val_E)
        call FEfunctionQuadValue(B, Th, Bh, i_elem, DERIV_NONE, Gauss_type, val_B)
        lorentz_value(1,:) = -val_sigma(1,:)*val_E(1,:)*val_B(2,:)
        lorentz_value(2,:) = val_sigma(1,:)*val_E(1,:)*val_B(1,:)
        
        call getGaussAnyElement(Th%NodeCoord(:,Th%ElemNodeConn(:,i_elem)), Gauss_type, pts, w)
        call BasisLocal2D(ref_pt_h, Th, Uh, i_elem, DERIV_NONE, basis_uh)
        
        do i_basis = 1, Uh%N_local_basis
            lorentz_to_basis(1,i_basis) = sum(lorentz_value(1,:)*basis_uh(1, i_basis, :)*w)
            lorentz_to_basis(2,i_basis) = sum(lorentz_value(2,:)*basis_uh(2, i_basis, :)*w)
        end do
        
        do i_dim = 1,DIM__
            call getLocalDofIndex(Uh, i_elem, i_dim, idx_local_dof)
            F%data(idx_local_dof) = F%data(idx_local_dof) + lorentz_to_basis(i_dim,:)
        end do
    end do
    
    block
        type(MATRIX_TRIPLET) :: M
        type(MATRIX_COLUMN) :: M_col
        type(vector) :: F_b
        integer, dimension(:,:), allocatable :: assemble_info
        
        call MatrixTripletInit(M, Uh%N_DOF, Uh%N_DOF, 2*Th%N_elem*Uh%N_local_basis*Uh%N_local_basis)
        allocate(assemble_info(2,5))
        assemble_info(1,:) = [1, 1, DERIV_NONE, 1, DERIV_NONE]
        assemble_info(2,:) = [1, 2, DERIV_NONE, 2, DERIV_NONE]
        call AssembleMatrixElement(one_func, [1d0, 1d0], Th, Uh, 0, Uh, 0, assemble_info, Gauss_type, M)
        
        call MatrixTriplet2Column(M, M_col)
        call VectorCopy(F_b, F)
        call SolverSolveUMFPACK2(M_col, F_b, F)
    end block

end subroutine magnetic_lorentz

subroutine magnetic_plot(filename)
    implicit none
    character(len=*) :: filename ! 输出vtk文件名

    call PlotFunction(B, Th, Bh, filename)
end subroutine

subroutine magnetic_update(nodes)
    use mesh, only: getMeshhmax
    implicit none
    real(8), dimension(:,:), intent(in) :: nodes ! 新节点坐标
    
    ! update mesh
    Th%NodeCoord = nodes
    call getMeshhmax(Th)

end subroutine magnetic_update

subroutine magnetic_Bmax(i_elem, Bmax)
    implicit none
    integer, intent(in) :: i_elem
    real(8), intent(out) :: Bmax 
    
    call QuadLinfNorm(B, Th, Bh, i_elem, DERIV_NONE, Gauss_type, Bmax)
    
end subroutine magnetic_Bmax
    
end module magnetic