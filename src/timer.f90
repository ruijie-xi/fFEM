! timer module
! author: Ruijie Xi
! 
! This module provides a timer to measure the time elapsed between two points.
!
! Usage:
!   use timer
!   call timer_start()
!   ! do something
!   call timer_end(time)
!   print *, time
module timer
    implicit none
    integer :: count_start,count_end,count_max
    real :: time_start,time_end,elapsed_time,count_rate
    public :: timer_start,timer_end
    private
contains

    subroutine timer_start()
        call system_clock(count_start,count_rate,count_max)
        time_start = real(count_start)/count_rate
    end subroutine timer_start

    subroutine timer_end(time)
        real, intent(out) :: time
        call system_clock(count_end,count_rate,count_max)
        time_end = real(count_end)/count_rate
        elapsed_time = time_end - time_start
        time = elapsed_time
    end subroutine timer_end
    
end module timer