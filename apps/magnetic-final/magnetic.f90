module magnetic
    use settings
    use timer
    use fe, only: fespaceInit, Interpolate
    use fe_utils
    use solver_umfpack2
    use assembler
    use tools
    use matvec
    use mesh
    use visualize
    implicit none

    private
    public :: magnetic_diffusion, magnetic_init, magnetic_update, magnetic_lorentz, magnetic_finalize
    public :: magnetic_plot

    common /global/ t
    real(8) :: t

    ! Gaussian quadrature type
    integer :: Gauss_type = QuadPt16, Gauss_type_bdry=LinePt4

    ! mesh and fespace
    type(MESH2D) :: Th
    type(FESPACE) :: Bh,Eh,Qh,Uh

    type(vector) :: B, E, sigma, sigma_real
    real(8),parameter :: mu = 4d0*m_pi
    real(8), dimension(:,:), allocatable :: current

    ! theta-scheme
    real(8) :: theta = 5d-1
    
    ! visualization
    character(len=100) :: vtkdir
    character(len=100) :: timeseriesfile = "time.vtk.series"
    integer :: timeseriesio
    
    contains
    
subroutine read_current_file(current_file)
    implicit none
    character(len=100),intent(in) :: current_file
    
    integer :: io,reason
    integer :: i,N
    LOGICAL :: flag
    
    real(8) :: tmp1,tmp2
    
    open(newunit=io, file=current_file, status="old", action="read")
    
    i=0
    flag = .false.
    do while(.NOT.flag)
        read(io, *, iostat=reason) tmp1, tmp2
        if(reason /= 0) then 
            flag=.true.
        else
            i = i+1
        end if
    end do
    
    N=i
    close(io)
    write(*,*) "read current file: detected ",N," current data."
    allocate(current(2,N))
    open(newunit=io, file=current_file, status="old", action="read")
    do i=1,N
        read(io, *) current(1,i), current(2,i)
    end do

    close(io)

    write(*,*) "read current file: done."
    
end subroutine 

subroutine get_current_value(time, current_value)
    implicit none
    real(8), intent(in) :: time
    real(8), intent(out) :: current_value

    integer :: i

    if(time<0) then
        write(*,*) "get_current_value: time should be non-negative."
        stop
    end if

    if(time<current(1,1)) then
        current_value = current(2,1)*time/current(1,1)
        return
    else
        do i=1,size(current,2)-1
            if(current(1,i)<=time .and. time<current(1,i+1)) then
                current_value = current(2,i) + (current(2,i+1)-current(2,i)) &
                * (time-current(1,i)) / (current(1,i+1)-current(1,i))
                return
            end if
        end do

        write(*,*) "get_current_value: time exceeds the range of current data."
        write(*,*) "time = ", time
        write(*,*) "set current value to the last value."
        current_value = current(2,size(current,2))
    end if

end subroutine

subroutine magnetic_init(N_elem, elems_in, N_node, x, y, code, N_material)
    implicit none
    integer, intent(in) :: elems_in(6,*)
    real(8), intent(in) :: x(*),y(*),code(*)
    integer, intent(in) :: N_elem, N_node, N_material
    character(len=100) :: current_file
    character(len=100) :: resist_file

    integer, parameter :: DOF_type_Q = DOF_Q0
    integer, parameter :: DOF_type_E = DOF_Q2
    integer, parameter :: DOF_type_B = DOF_QuadRT2
    integer, parameter :: DOF_type_U = DOF_Q1

    real :: t_test
    
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes
    real(8), dimension(:), allocatable :: resist,resist_real
    real(8), parameter :: resist_max = 1d6
    integer :: i, io, reason
    
    call timer_start()
    
    ! read in all the data filenames
    call read_config("mag_datafile", resist_file, current_file, vtkdir)
    
    allocate(elems(4,N_elem),nodes(2,N_node))
    
    call VectorInit(sigma, N_elem)
    call VectorInit(sigma_real, N_elem)
    
    allocate(resist(N_material),resist_real(N_material))
    ! read in resistivity
    open(newunit=io, file=resist_file, status="old", action="read")
    do i = 1, N_material
        read(io, *, iostat=reason) resist_real(i)
        if(reason/=0) then
            write(*,*) "error reading resistivity of material ",i
            stop
        else
            write(*,*) "resistivity of material ", i, ": ",resist_real(i)
            resist(i) = resist_real(i)
            if (resist_real(i)>resist_max) then
                write(*,*) "we will set it to ", resist_max, "in magnetic computation."
                resist(i) = resist_max
            end if
        end if
    end do
    close(io)
    
    do i=1,N_elem
          sigma%data(i) = 1d0/resist(elems_in(5,i))
          sigma_real%data(i) = 1d0/resist_real(elems_in(5,i))
          elems(:,i) = elems_in(1:4,i)
    end do
    
    do i = 1,N_node
        nodes(1,i) = x(i)
        nodes(2,i) = y(i)
    end do

    ! generate mesh topological structure
    call MeshInit(elems,nodes,Th)
    call MarkOuterBdry(2, code)
    write(*,*) "Mesh generation done. "
    !call PrintMesh(Th)

    ! initialize fe spaces
    call fespaceInit(Qh, Th, DOF_type_Q, 1)
    call fespaceInit(Eh, Th, DOF_type_E, 1)
    call fespaceInit(Bh, Th, DOF_type_B, 2)
    call fespaceInit(Uh, Th, DOF_type_U, 2)

    ! initialize functions
    call VectorInit(B, Bh%N_dof)
    call VectorInit(E, Eh%N_dof)
    
    ! read current data
    call read_current_file(current_file)
    
    ! time series file
    open(newunit=timeseriesio, file=trim(vtkdir)//"/"//timeseriesfile, status="replace")
    write(timeseriesio, "(A)") "{"
    write(timeseriesio, "(A)") '"file-series-version" : "1.0",'
    write(timeseriesio, "(A)") '"files" : ['
    close(timeseriesio)
    
    call timer_end(t_test)
    write(*,*) "Initialization of magnetic module done. Time taken = ", t_test
    
end subroutine magnetic_init


subroutine read_config(config_name, resist_file, current_file, vtk_dir)
    integer :: io
    character(len=*), intent(in) :: config_name
    character(len=100) :: data_file
    character(len=100) :: resist_file, current_file, vtk_dir
    
    open(newunit=io, file=config_name, status="old", action="read")
    read(io, '(A)') data_file
    close(io)
    
    open(newunit=io, file=data_file, status="old", action="read")
    call read_skip_comment(io, resist_file)
    write(*,*) "resistivity file: ", resist_file
    call read_skip_comment(io, current_file)
    write(*,*) "current file: ", current_file
    call read_skip_comment(io, vtk_dir)
    write(*,*) "vtr directory: ", vtk_dir
    close(io)
    
end subroutine read_config

subroutine read_skip_comment(io, line)
    integer, intent(in) :: io
    character(len=100), intent(out) :: line
    
    logical :: flag
    
    flag=.false.
    
    do while(.not.flag)
        read(io, '(A)') line
        if(line(1:1)/="#")then
            flag = .true.
        end if
    end do
    
end subroutine read_skip_comment


subroutine magnetic_finalize()
    
    open(newunit=timeseriesio, file=trim(vtkdir)//"/"//timeseriesfile, action='write',position='append')
    write(timeseriesio, "(A)") "  ]"
    write(timeseriesio, "(A)") "}"
    close(timeseriesio)

end subroutine

! advance the magnetic field from to -> to+dt 
subroutine magnetic_diffusion(to, dt)
    implicit none
    real(8), intent(in) :: to, dt
    
    ! functions
    type(vector) :: Bhat,curlE,Bold

    ! matrices and vectors
    integer, dimension(:,:), allocatable :: assemble_info
    type(MATRIX_TRIPLET) :: matM,matM1,matTotal
    type(vector) :: vec_B,vec_g,vecTotal

    ! current data
    real(8) :: current_value

    call CurlOperator(E, curlE)
    call VectorCopy(Bhat, B)
    call VectorAddVector(Bhat, 1d0, curlE, -(1-theta)*dt)
    call VectorCopy(Bold, B)
    
    t = to + dt
    print *, "Magnetic diffusion: advancing from ", to, " to ", t

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
    call AssembleVectorElementFE(one_func, [1d0/mu, -1d0/mu], Th, Eh, 0, Bhat, Bh, &
    assemble_info, Gauss_type, vec_B)
    deallocate(assemble_info)

    ! assemble vector vec_g
    call get_current_value(t, current_value)
    current_value = current_value * 2d0 ! B=2S(t)/r
    call VectorInit(vec_g, Eh%N_DOF)
    allocate(assemble_info(1,3))
    assemble_info(1,:) = [1, 1, DERIV_NONE]
    call AssembleVectorBdry(2, one_on_r_func, [current_value/mu], Th, Eh, 0, assemble_info, Gauss_type_bdry, vec_g)
    deallocate(assemble_info)

    ! form the linear system
    ! matTotal = matM + theta*dt*matM1
    ! vecTotal = vec_B + vec_g
    call MatrixTripletAdd(1d0, matM, theta*dt, matM1, matTotal)
    call VectorInit(vecTotal, Eh%N_DOF)
    call VectorCopy(vecTotal, vec_B)
    call VectorAddVector(vecTotal, 1d0, vec_g, 1d0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! solve the linear system !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call SolverSolveUMFPACK2(matTotal, vecTotal, E)

    call CurlOperator(E, curlE)
    call VectorCopy(B, Bhat)
    call VectorAddVector(B, 1d0, curlE, -theta*dt)
    
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

! one on r
subroutine one_on_r_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 1d0/sqrt(x(1)**2+x(2)**2)
    end select
end subroutine one_on_r_func

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

subroutine magnetic_update(x,y)
    use mesh, only: getMeshhmax
    implicit none
    real(8), intent(in) :: x(*),y(*)
    
    ! update mesh
    Th%NodeCoord(1,1:Th%N_node) = x(1:Th%N_node)
    Th%NodeCoord(2,1:Th%N_node) = y(1:Th%N_node)
    call getMeshhmax(Th)

end subroutine magnetic_update

subroutine magnetic_errorB(B_func, L2err)
    implicit none
    procedure(func) :: B_func
    real(8), intent(out) :: L2err

    call ComputeError(B_func, B, Th, Bh, NORM_L2, Gauss_type, L2err)
end subroutine

subroutine magnetic_errorE(E_func, L2err)
    implicit none
    procedure(func) :: E_func
    real(8), intent(out) :: L2err

    call ComputeError(E_func, E, Th, Eh, NORM_L2, Gauss_type, L2err)
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
subroutine magnetic_lorentz(F_out)
    implicit none
    real(8), dimension(2,*) :: F_out
    !real(8), dimension(:), allocatable :: F
    real(8), dimension(:,:), allocatable :: val_sigma, val_E, val_B
    real(8), dimension(:,:), allocatable :: node_value ! value of the function at each node
    real(8) :: elem_area
    integer,dimension(:), allocatable :: node_countvalue ! number of values for each node
    integer :: i_elem

    !allocate(F(Uh%N_DOF))

    allocate(node_countvalue(Th%N_node))
    allocate(node_value(Uh%dim,Th%N_node))
    allocate(val_E(Eh%dim,Th%N_le), val_B(Bh%dim,Th%N_le), val_sigma(Qh%dim,Th%N_le))

    node_countvalue = 0
    node_value = 0d0

    if (Th%N_le.eq.3) then
        Gauss_type = TriangleNodeAverage
    else if (Th%N_le.eq.4) then
        Gauss_type = QuadNodeAverage
    else
        stop 'GetLorentzForce : not implemented for this element type'
    end if

    do i_elem = 1, Th%N_elem
        call FEfunctionQuadValue(sigma_real, Th, Qh, i_elem, DERIV_NONE, Gauss_type, val_sigma)
        call FEfunctionQuadValue(E, Th, Eh, i_elem, DERIV_NONE, Gauss_type, val_E)
        call FEfunctionQuadValue(B, Th, Bh, i_elem, DERIV_NONE, Gauss_type, val_B)
        call getElementArea(Th, i_elem, elem_area)
        !node_countvalue(Th%ElemNodeConn(:,i_elem)) = node_countvalue(Th%ElemNodeConn(:,i_elem)) + 1
        node_value(1,Th%ElemNodeConn(:,i_elem)) = node_value(1,Th%ElemNodeConn(:,i_elem)) &
        - val_sigma(1,:)*val_E(1,:)*val_B(2,:)*elem_area/4d0
        node_value(2,Th%ElemNodeConn(:,i_elem)) = node_value(2,Th%ElemNodeConn(:,i_elem)) &
        + val_sigma(1,:)*val_E(1,:)*val_B(1,:)*elem_area/4d0
    end do
    
    F_out(1,1:Th%N_node) = node_value(1,:)
    F_out(2,1:Th%N_node) = node_value(2,:)

end subroutine magnetic_lorentz

subroutine magnetic_plot(itts, time)
    implicit none
    character(len=30) :: itstr, tstr
    character(len=30) :: Bfilename, Efilename
    integer,intent(in) :: itts
    real(8),intent(in) :: time
    
    write(itstr,"(I0)") itts
    write(tstr,"(E11.4)") time
    write(Bfilename, "(A)") "B_"//trim(itstr)//".vtk"
    write(Efilename, "(A)") "E_"//trim(itstr)//".vtk"


    call PlotFunction(B, Th, Bh, trim(vtkdir)//"/"//trim(Bfilename))
    call PlotFunction(E, Th, Eh, trim(vtkdir)//"/"//trim(Efilename))
    
    open(newunit=timeseriesio, file=trim(vtkdir)//"/"//trim(timeseriesfile), Access = 'append',Status='old')
    write(timeseriesio, "(A)") '    { "name" : "'//trim(Bfilename)//'", "time" : '//trim(tstr)//' },'
    close(timeseriesio)

end subroutine

subroutine MarkOuterBdry(marker, code)
    integer,intent(in) :: marker
    real(8) :: code(*)

    integer :: i_bdry, i_edge, i_node1, i_node2
    real(8),parameter :: eps = 5d-4

    do i_bdry = 1,Th%N_bdryedge
        i_edge = Th%BdryEdge(i_bdry)
        i_node1 = Th%EdgeNodeConn(1,i_edge)
        i_node2 = Th%EdgeNodeConn(2,i_edge)
            
        if (abs(code(i_node1)-3d0)<eps .and. abs(code(i_node2)-3d0)<eps) then
            Th%BdryMarker(i_bdry) = marker
        end if
    end do
end subroutine MarkOuterBdry
    
end module magnetic