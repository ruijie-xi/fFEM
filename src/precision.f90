! This is a module that defines the precision of the real numbers.
module precision
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: sp = selected_real_kind(6, 37)

    public :: dp, sp, print_precision

contains
        subroutine print_precision()
            write(*,*) "double precision: ", dp
            write(*,*) "single precision: ", sp
        end subroutine print_precision
    
end module precision