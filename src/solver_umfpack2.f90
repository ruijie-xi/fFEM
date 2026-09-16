module solver_umfpack2
    use settings
    use matvec, only: MatrixColumn2Triplet, VectorReset
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    private
    public :: SolverSolveUMFPACK2

    interface SolverSolveUMFPACK2
        module procedure SolverSolveUMFPACK2_Column
        module procedure SolverSolveUMFPACK2_Triplet
    end interface

contains

    subroutine SolverSolveUMFPACK2_Column(A, b, x, status)
        type(MATRIX_COLUMN), intent(in) :: A
        type(VECTOR), intent(in) :: b
        type(VECTOR), intent(inout) :: x
        integer, intent(out), optional :: status
        type(MATRIX_TRIPLET) :: triplet

        call MatrixColumn2Triplet(A, triplet)
        call SolverSolveUMFPACK2_Triplet(triplet, b, x, status)
    end subroutine

    subroutine SolverSolveUMFPACK2_Triplet(A, b, x, status)
        type(MATRIX_TRIPLET), intent(in) :: A
        type(VECTOR), intent(in) :: b
        type(VECTOR), intent(inout) :: x
        ! Zero means success; nonzero is UMFPACK's status or -1000 for invalid input.
        ! Without status, failure terminates rather than returning an invalid solution.
        integer, intent(out), optional :: status
        real(8), allocatable :: w(:), value(:)
        integer, allocatable :: index(:)
        integer :: icntl(20), keep(20), info(40), lvalue, lindex, nnz
        real(8) :: cntl(10), rinfo(20)

        if (present(status)) status = 0
        call VectorReset(x, max(0, A%N_col))
        if (A%N_row /= A%N_col .or. A%N_row <= 0 .or. b%size /= A%N_row) then
            call solver_failure(-1000, x, status)
            return
        end if
        if (.not. allocated(b%data)) then
            call solver_failure(-1000, x, status)
            return
        end if
        if (size(b%data) /= b%size .or. .not. all(ieee_is_finite(b%data))) then
            call solver_failure(-1000, x, status)
            return
        end if
        nnz = A%actual_nnz
        if (nnz <= 0) then
            call solver_failure(-1000, x, status)
            return
        end if
        if (any(A%row_idx(:nnz) < 1) .or. any(A%row_idx(:nnz) > A%N_row) .or. &
            any(A%col_idx(:nnz) < 1) .or. any(A%col_idx(:nnz) > A%N_col) .or. &
            .not. all(ieee_is_finite(A%val(:nnz)))) then
            call solver_failure(-1000, x, status)
            return
        end if

        call ums2in(icntl, cntl, keep)
        lvalue = max(100*nnz, 4*A%N_row)
        lindex = max(100*nnz, 3*nnz + 64*A%N_row + 16)
        allocate(index(lindex), value(lvalue), w(4*A%N_row))
        index(:nnz) = A%row_idx(:nnz)
        index(nnz+1:2*nnz) = A%col_idx(:nnz)
        value(:nnz) = A%val(:nnz)

        call ums2fa(A%N_row, nnz, 0, .false., lvalue, lindex, value, index, &
            keep, cntl, icntl, info, rinfo)
        ! UMFPACK status 2 means duplicate triplets were summed, which is valid FEM input.
        if (info(1) /= 0 .and. info(1) /= 2) then
            call solver_failure(info(1), x, status)
            return
        end if
        call ums2so(A%N_row, 0, .false., lvalue, lindex, value, index, &
            keep, b%data, x%data, w, cntl, icntl, info, rinfo)
        if (info(1) /= 0) then
            call solver_failure(info(1), x, status)
        else if (.not. all(ieee_is_finite(x%data))) then
            call solver_failure(-1000, x, status)
        end if
    end subroutine

    subroutine solver_failure(code, x, status)
        integer, intent(in) :: code
        type(VECTOR), intent(inout) :: x
        integer, intent(out), optional :: status
        x%data = ieee_value(0d0, ieee_quiet_nan)
        if (present(status)) then
            status = code
        else
            write(error_unit, *) 'SolverSolveUMFPACK2 failed, status = ', code
            error stop 1
        end if
    end subroutine
end module solver_umfpack2
