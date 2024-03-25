module matvec
    use settings
    use quicksort_module
    implicit none
    
contains

    ! Initialize a matrix in triplet format
    subroutine MatrixTripletInit(A, nrows, ncols, nnz)
        implicit none
        integer, intent(in) :: nrows, ncols, nnz
        type(MATRIX_TRIPLET), intent(out) :: A
        integer :: i
        
        A%N_row = nrows
        A%N_col = ncols
        A%N_nz = nnz
        
        allocate(A%row_idx(nnz), A%col_idx(nnz), A%val(nnz))

        do i = 1, nnz
            A%row_idx(i) = 0
            A%col_idx(i) = 0
            A%val(i) = 0.0
        end do

        A%actual_nnz = 0

    end subroutine

    subroutine MatrixTripletFree(A)
        implicit none
        type(MATRIX_TRIPLET), intent(inout) :: A

        deallocate(A%row_idx, A%col_idx, A%val)

        A%N_row = 0
        A%N_col = 0
        A%N_nz = 0
        A%actual_nnz = 0

    end subroutine

    subroutine MatrixTripletAddValues(A, idx_row, idx_col, val)
        implicit none
        type(MATRIX_TRIPLET), intent(inout) :: A
        integer, intent(in),dimension(:) :: idx_row, idx_col
        real(8), intent(in),dimension(:) :: val

        integer :: num_add_nnz

        if(size(idx_row) /= size(idx_col) .or. size(idx_row) /= size(val) .or. size(idx_col) /= size(val)) then
            print *, "Error: size of idx_row, idx_col, and val must be the same"
            stop
        end if
        num_add_nnz = size(idx_row)

        if(A%actual_nnz + num_add_nnz > A%N_nz) then
            print *, "Error: Adding too many non-zero elements"
            stop
        end if

        A%row_idx(A%actual_nnz+1:A%actual_nnz+num_add_nnz) = idx_row
        A%col_idx(A%actual_nnz+1:A%actual_nnz+num_add_nnz) = idx_col
        A%val(A%actual_nnz+1:A%actual_nnz+num_add_nnz) = val

        A%actual_nnz = A%actual_nnz + num_add_nnz

    end subroutine

    ! Add matrix-like values
    subroutine MatrixTripletAddMatrixValues(A, ind_row, ind_col, values)
        implicit none
        type(MATRIX_TRIPLET), intent(inout) :: A
        integer, intent(in), dimension(:) :: ind_row, ind_col
        real(8), intent(in), dimension(:,:) :: values

        integer :: i, j, k

        integer :: nrows, ncols, nnz

        nrows = size(ind_row)
        ncols = size(ind_col)
        nnz = nrows*ncols

        if(size(values,1) /= nrows .or. size(values,2) /= ncols) then
            print *, "Error: size of values must be the same as ind_row and ind_col"
            stop
        end if

        if(A%actual_nnz + nnz > A%N_nz) then
            print *, "Error: Adding too many non-zero elements"
            stop
        end if

        k = A%actual_nnz
        do i = 1, nrows
            do j = 1, ncols
                k = k + 1
                A%row_idx(k) = ind_row(i)
                A%col_idx(k) = ind_col(j)
                A%val(k) = values(i,j)
            end do
        end do

        A%actual_nnz = A%actual_nnz + nnz

    end subroutine

    subroutine MatrixTripletAdd(coeff1, mat1, coeff2, mat2, result)
        implicit none
        real(8), intent(in) :: coeff1, coeff2
        type(MATRIX_TRIPLET), intent(in) :: mat1, mat2
        type(MATRIX_TRIPLET), intent(out) :: result

        if(mat1%N_row /= mat2%N_row .or. mat1%N_col /= mat2%N_col) then
            print *, "Error: Matrices must have the same size"
            stop
        end if

        call MatrixTripletInit(result, mat1%N_row, mat1%N_col, mat1%N_nz + mat2%N_nz)

        call MatrixTripletAddValues(result, mat1%row_idx, mat1%col_idx, coeff1*mat1%val)
        call MatrixTripletAddValues(result, mat2%row_idx, mat2%col_idx, coeff2*mat2%val)
    end subroutine

    ! convert triplet format to column format
    ! TODO: matrix with general dimensions (m*n)
    subroutine MatrixTriplet2Column(A_triplet, A_column)
        implicit none
        type(MATRIX_TRIPLET), intent(in) :: A_triplet
        type(MATRIX_COLUMN), intent(out) :: A_column

        real(8), dimension(:), allocatable :: xx
        integer, dimension(:), allocatable :: ii,w,wp

        integer :: N_nz

        integer :: icntl(20),keep(20),info(40)
        real(8) :: cntl(10)

        integer :: i, i_col,nnz_col
        integer,dimension(:),allocatable :: index

        if(A_triplet%N_row /= A_triplet%N_col) then
            print *, "Error: Matrix must be square"
            stop
        end if

        if(A_triplet%actual_nnz /= A_triplet%N_nz) then
            print *, "Warning: Matrix is not fully filled, and will be converted to column format anyway."
        end if

        ! use UMFPACK to convert triplet to column format (square matrix only)

        call ums2in(icntl, cntl, keep)

        allocate(xx(2*A_triplet%actual_nnz))
        xx(1:A_triplet%actual_nnz) = A_triplet%val(1:A_triplet%actual_nnz)
        allocate(ii(A_triplet%actual_nnz+max(2*A_triplet%actual_nnz, A_triplet%N_row+1)))
        ii(1:A_triplet%actual_nnz) = A_triplet%row_idx(1:A_triplet%actual_nnz)
        ii(A_triplet%actual_nnz+1:2*A_triplet%actual_nnz) = A_triplet%col_idx(1:A_triplet%actual_nnz)

        allocate(w(A_triplet%N_row),wp(A_triplet%N_row+1))

        N_nz = A_triplet%actual_nnz

        call ums2co(A_triplet%N_row, N_nz, .false., &
        xx, size(xx), info, icntl, ii, size(ii), w, wp, 1)

        A_column%N_row = A_triplet%N_row
        A_column%N_col = A_triplet%N_col
        A_column%N_nz = N_nz

        allocate(A_column%row_idx(N_nz), A_column%col_ptr(A_triplet%N_row+1), A_column%val(N_nz))
        A_column%col_ptr = ii(1:A_column%N_col+1)
        A_column%row_idx = ii(A_column%N_col+2:A_column%N_col+1+N_nz)
        A_column%val = xx(1:N_nz)

        deallocate(xx, ii, w, wp)

        ! sort row indices
        do i_col = 1, A_column%N_col
            nnz_col = A_column%col_ptr(i_col+1) - A_column%col_ptr(i_col)
            if (nnz_col == 0) cycle
            allocate(index(nnz_col))
            index = [(i, i=1, nnz_col)]
            call quicksort(A_column%row_idx(A_column%col_ptr(i_col):A_column%col_ptr(i_col+1)-1), index)
            A_column%val(A_column%col_ptr(i_col):A_column%col_ptr(i_col+1)-1) = A_column%val(A_column%col_ptr(i_col)+index-1)
            deallocate(index)
        end do

    end subroutine

    subroutine MatrixColumn2Triplet(A_column, A_triplet)
        implicit none
        type(MATRIX_COLUMN), intent(in) :: A_column
        type(MATRIX_TRIPLET), intent(out) :: A_triplet

        integer :: i, j, k

        A_triplet%N_row = A_column%N_row
        A_triplet%N_col = A_column%N_col
        A_triplet%N_nz = A_column%N_nz

        allocate(A_triplet%row_idx(A_column%N_nz), A_triplet%col_idx(A_column%N_nz), A_triplet%val(A_column%N_nz))

        k = 0
        do j = 1, A_column%N_col
            do i = A_column%col_ptr(j), A_column%col_ptr(j+1)-1
                k = k + 1
                A_triplet%row_idx(k) = A_column%row_idx(i)
                A_triplet%col_idx(k) = j
                A_triplet%val(k) = A_column%val(i)
            end do
        end do

        A_triplet%actual_nnz = A_column%N_nz

    end subroutine

    ! get the value of the matrix at (i,j)
    subroutine MatrixColumnGet(A,i,j,value)
        implicit none
        type(MATRIX_COLUMN), intent(in) :: A
        integer, intent(in) :: i,j
        real(8), intent(out) :: value

        integer :: k

        value = 0d0
        do k = A%col_ptr(j), A%col_ptr(j+1)-1
            if(A%row_idx(k) == i) then
                value = A%val(k)
                return
            end if
        end do
    end subroutine

    subroutine MatrixColumnPrint(A)
        implicit none
        type(MATRIX_COLUMN), intent(in) :: A

        integer :: i

        print *, "N_row: ",A%N_row, "N_col: ",A%N_col, "N_nz: ",A%N_nz

        ! show column pointers
        print *, "column pointers"
        do i = 1, A%N_col+1
            print *, i, A%col_ptr(i)
        end do

        ! show row indices of non-zero elements
        print *, "row indices"
        do i = 1, A%N_nz
            print *, i, A%row_idx(i)
        end do

        ! show value
        print *, "values"
        do i = 1, A%N_nz
            print *, i, A%val(i)
        end do

    end subroutine

    subroutine MatrixTripletPrint(A)
        implicit none
        type(MATRIX_TRIPLET), intent(in) :: A

        integer :: i

        print *, "N_row: ",A%N_row, "N_col: ",A%N_col, "N_nz: ",A%N_nz

        ! show row indices of non-zero elements
        print *, "row indices"
        do i = 1, A%N_nz
            print *, i, A%row_idx(i)
        end do

        ! show column indices of non-zero elements
        print *, "column indices"
        do i = 1, A%N_nz
            print *, i, A%col_idx(i)
        end do

        ! show value
        print *, "values"
        do i = 1, A%N_nz
            print *, i, A%val(i)
        end do

    end subroutine

    subroutine VectorInit(vec, n)
        implicit none
        real(8), dimension(:), allocatable, intent(out) :: vec
        integer, intent(in) :: n

        allocate(vec(n))
        vec = 0d0

    end subroutine

    subroutine VectorPrint(vec)
        implicit none
        real(8), dimension(:), intent(in) :: vec

        integer :: i

        do i = 1, size(vec)
            print *, i, vec(i)
        end do

    end subroutine

    subroutine VectorFree(vec)
        implicit none
        real(8), dimension(:), allocatable, intent(inout) :: vec

        deallocate(vec)

    end subroutine


    
end module matvec