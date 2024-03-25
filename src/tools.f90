module tools
    implicit none
    
contains

    subroutine assert(x,str)
        logical,intent(in) :: x
        character(len=*),intent(in) :: str

        if (.not. x) then
            write(*,*) "assertion failed: ",str
            stop
        end if
    end subroutine assert

    subroutine linspace(left,right,N,array)
        real(8),intent(in) :: left,right
        integer,intent(in) :: N
        real(8),dimension(:),allocatable,intent(out) :: array

        integer :: i
        
        call assert(.not.allocated(array),"array not allocated in linspace")
        allocate(array(N))
        do i = 1,N
            array(i) = left + (right-left)*(i-1d0)/(N-1d0)
        end do
        
    end subroutine linspace

    subroutine ViewArray(array)
        real(8),dimension(:),intent(in) :: array
        integer :: i

        do i = 1,size(array)
            write(*,*) i,array(i)
        end do
    end subroutine ViewArray

end module tools