module problem_mag
    use settings
    implicit none

    real(8) :: mu = 4*m_pi
    real(8) :: sigma0_ = 10d0
    
    real(8), dimension(:,:), allocatable :: current
    
contains

! right hand side function
subroutine f_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r

    r = sqrt(x(1)**2+x(2)**2)
    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine f_func

! velocity function
subroutine u_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))


    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
        f(2) = 0d0
    end select
end subroutine u_func
 
! magnetic field on boundary
subroutine g_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out),dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    real(8) :: t
    common /global/ t
    real(8) :: r
    
    ! current data
    real(8) :: current_value    

    r = sqrt(x(1)**2+x(2)**2)

    if (.not. allocated(f)) allocate(f(1))
    
    call get_current_value(t, current_value)

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 2d0*current_value/r
    end select
end subroutine g_func

! initial solution of B
subroutine B0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(2))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
        f(2) = 0d0
    end select
end subroutine B0_func

! initial solution of E
subroutine E0_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = 0d0
    end select
end subroutine E0_func

! function of sigma
subroutine sigma_func(x,f,deriv_type)
    use settings
    implicit none
    real(8), intent(in), dimension(:) :: x
    real(8), intent(out), dimension(:), allocatable :: f
    integer, intent(in) :: deriv_type

    if (.not. allocated(f)) allocate(f(1))

    select case(deriv_type)
    case(DERIV_NONE)
        f(1) = sigma0_
    end select
end subroutine sigma_func


subroutine read_current_file(current_file)
    implicit none
    character(len=*),intent(in) :: current_file
    
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
    
end module problem_mag

program test_mag
    use settings
    use problem_mag
    use Hydrodynamic
    use magnetic
    use mesh_generator
    use readmeshfile
    implicit none

    character(len=100) :: arg

    ! global variables
    real(8) :: t
    common /global/ t
    real(8) :: dt
    real(8) :: Tend = 1d0

    ! mesh
    integer :: M,Nlist(1)
    real(8),parameter :: Center(2) = [0d0,0d0], Rlist(1) = [1d0]
    integer, dimension(:,:), allocatable :: elems
    real(8), dimension(:,:), allocatable :: nodes

    ! velocity vector
    type(vector)  :: u

    ! Lorentz force
    type(vector)  :: F

    ! error
    real(8) :: L2_err, relL2_err

    ! parameters of Gauss quadrature
    integer, parameter :: Gauss_type = QuadPt16
    integer, parameter :: Gauss_type_bdry = LinePt4

    ! read mesh
    read_mesh: block
        call ReadNodeFile("nodes.txt", nodes)
        call ReadElemFile("elems.txt", elems)
    end block read_mesh

    call HydrodynamicInitialize(elems,nodes)
    
    ! read current data
    call read_current_file("current.dat")

    call magnetic_init(elems,nodes,B0_func,E0_func,sigma_func, 3d0)

    ! time step
    dt = 0.025d0
    t = 0d0
    do while(t<Tend-1d-8)
        ! Step 1: Get Lorentz force
        ! call magnetic_lorentz(F)

        ! Step 3: Solve diffusion equation, advance time
        call magnetic_diffusion(t, 0.5d0*dt, g_func, f_func, .false.)

        ! Step 4: Move mesh
        t = t+0.5d0*dt
        ! Calculate velocity at t+0.5*dt
        ! Move the mesh from t to t+dt (leapfrog)
        call getVelocity(u_func, u)
        call MoveMesh(u, dt, nodes)
        call magnetic_update(nodes)
        
        call magnetic_diffusion(t, 0.5d0*dt, g_func, f_func, .false.)
        
        t = t+0.5d0*dt
        
    end do

    call magnetic_plot("output/B.vtk")
    

end program test_mag