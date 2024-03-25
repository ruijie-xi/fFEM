program test_UMFPACK
    implicit none

    integer :: icntl(20),keep(20),info(40)
    real :: cntl(10),rinfo(20)

    integer, parameter :: n = 100000, ne = 2*100000
    integer,dimension(:),allocatable :: ia, ja
    real,dimension(:),allocatable :: xa,b,y,x, w

    integer,parameter :: lvalue=5000000, lindex=5000000
    real,dimension(:),allocatable :: value
    integer,dimension(:),allocatable :: index

    integer :: i

    ! initial settings
    call ums2in(icntl, cntl, keep)
    ! icntl(3) = 4
    
    ! triplet form of the matrix
    allocate(ia(ne),ja(ne),xa(ne))
    do i = 1,n
        ia(2*i-1) = i
        ja(2*i-1) = i
        xa(2*i-1) = 2.0*i
        ia(2*i) = i
        ja(2*i) = mod(i+3,n)+1
        xa(2*i) = -1.5*(i-3)
    end do

    ! rhs
    allocate(y(n),b(n),x(n),w(4*n))
    do i = 1,n
        y(i) = 3*i
    end do
    do i = 1,n
        b(i) = xa(2*i-1)*y(ja(2*i-1)) + xa(2*i)*y(ja(2*i))
    end do

    ! factorize
    allocate(value(lvalue),index(lindex))
    do i=1,ne
        index(i) = ia(i)
        index(ne+i) = ja(i)
        value(i) = xa(i)
    end do

    call ums2fa(n,ne,0,.false.,lvalue, lindex, value, index, &
             keep, cntl, icntl, info, rinfo)

    call ums2so(n,0,.false.,lvalue,lindex,value,index, &
            keep, b, x, w, cntl, icntl, info, rinfo)

    do i=1,n
        print *, x(i)
    end do



end program test_UMFPACK