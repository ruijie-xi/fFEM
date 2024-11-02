! TODO: make it faster
module matvec
    use settings
    use ffem_quicksort
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
    
    subroutine MatrixTripletAddValue(A, row, col, val)
        implicit none
        type(MATRIX_TRIPLET), intent(inout) :: A
        integer, intent(in) :: row, col
        real(8), intent(in) :: val
        
        if(A%actual_nnz >= A%N_nz) then
            print *, "Error: Adding too many non-zero elements"
            stop
        end if
        
        A%actual_nnz = A%actual_nnz + 1
        A%row_idx(A%actual_nnz) = row
        A%col_idx(A%actual_nnz) = col
        A%val(A%actual_nnz) = val
    
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
    subroutine MatrixTriplet2Column(A_triplet, A_column)
        implicit none
        type(MATRIX_TRIPLET), intent(inout) :: A_triplet
        type(MATRIX_COLUMN), intent(out) :: A_column

        integer :: N_invalid,N_duplicate,N_nz
        integer :: i,row,col
        integer, dimension(:), allocatable :: w_col, w_row

        integer,dimension(:),allocatable :: col_ptr, row_idx
        real(8),dimension(:),allocatable :: val

        integer :: pdest, pcol

        integer, dimension(:), allocatable :: index
        integer :: nnz_col

        ! make sure the matrix is fully assembled
        if(A_triplet%actual_nnz /= A_triplet%N_nz) then
            print *, "Warning: Matrix is not fully filled. "
        end if
        
        N_nz = A_triplet%actual_nnz

        ! check invalid elements and count the number of non-zero elements in each column
        N_invalid = 0
        allocate(w_col(A_triplet%N_col))
        w_col = 0
        do i = 1, A_triplet%actual_nnz
            row = A_triplet%row_idx(i)
            col = A_triplet%col_idx(i)
            if(row < 1 .or. row > A_triplet%N_row .or. col < 1 .or. col > A_triplet%N_col) then
                ! invalid element
                N_invalid = N_invalid + 1
                print *, "Warning: invalid element at (",row,",",col,")"
            else
                w_col(col) = w_col(col) + 1
            end if
        end do

        if(N_invalid > 0) then
            print *, "Error: Found ",N_invalid," invalid elements"
            call exit(1)
        end if

        ! convert triplet to column format
        A_column%N_row = A_triplet%N_row
        A_column%N_col = A_triplet%N_col
        allocate(col_ptr(A_triplet%N_col+1))
        allocate(row_idx(N_nz))
        allocate(val(N_nz))

        ! construct column pointers
        col_ptr(1) = 1
        do i = 1, A_column%N_col
            col_ptr(i+1) = col_ptr(i) + w_col(i)
        end do
        w_col = col_ptr(1:A_column%N_col) ! save the last pointer position of each column

        ! construct row indices and values
        do i = 1, N_nz
            row = A_triplet%row_idx(i)
            col = A_triplet%col_idx(i)
            row_idx(w_col(col)) = row
            val(w_col(col)) = A_triplet%val(i)
            w_col(col) = w_col(col) + 1
        end do
        
        ! remove duplicate elements
        ! pcol: column pointer of the current column
        ! w_row: the last nz pointer position of each row
        ! pdest: nz pointer
        N_duplicate = 0
        pdest = 1
        allocate(w_row(A_triplet%N_row))
        w_row = 0
        do col = 1,A_column%N_col
            pcol = pdest
            do i = col_ptr(col),col_ptr(col+1)-1
                row = row_idx(i)
                if(w_row(row) >= pcol) then
                    ! duplicate element
                    val(w_row(row)) = val(w_row(row)) + val(i)
                    N_duplicate = N_duplicate + 1
                else
                    w_row(row) = pdest
                    if(pdest /= i) then
                        ! move element
                        row_idx(pdest) = row
                        val(pdest) = val(i)
                    end if
                    pdest = pdest + 1
                end if
            end do
            col_ptr(col) = pcol
        end do

        col_ptr(A_column%N_col+1) = pdest
        N_nz = pdest - 1
        if(N_duplicate > 0) then
            print *, "Warning: Removed ",N_duplicate," duplicate elements"
        end if

        ! construct A_column
        A_column%N_nz = N_nz
        allocate(A_column%row_idx(N_nz), A_column%val(N_nz))
        A_column%col_ptr = col_ptr
        A_column%row_idx = row_idx(1:N_nz)
        A_column%val = val(1:N_nz)

        deallocate(w_col, w_row)
        deallocate(col_ptr, row_idx, val)
        
        ! sort row indices
        do col = 1, A_column%N_col
            nnz_col = A_column%col_ptr(col+1) - A_column%col_ptr(col)
            if (nnz_col == 0) cycle
            allocate(index(nnz_col))
            index = [(i, i=1, nnz_col)]
            call quicksort(A_column%row_idx(A_column%col_ptr(col):A_column%col_ptr(col+1)-1), index)
            A_column%val(A_column%col_ptr(col):A_column%col_ptr(col+1)-1) = A_column%val(A_column%col_ptr(col)+index-1)
            deallocate(index)
        end do

    end subroutine

    ! deprecated
    ! convert triplet format to column format
    ! using UMFPACK2.0
    ! square matrix only
    subroutine MatrixTriplet2Column_(A_triplet, A_column)
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

    ! remove duplicate entries in a matrix in triplet format
    ! the space remains the same
    subroutine MatrixTripletAssemble(A)
        implicit none
        type(MATRIX_TRIPLET), intent(inout) :: A
        type(MATRIX_COLUMN) :: A_col
        integer :: i, j, k

        call MatrixTriplet2Column(A, A_col)

        A%row_idx = 0
        A%col_idx = 0
        A%val = 0.0

        k = 0
        do j = 1, A_col%N_col
            do i = A_col%col_ptr(j), A_col%col_ptr(j+1)-1
                k = k + 1
                A%row_idx(k) = A_col%row_idx(i)
                A%col_idx(k) = j
                A%val(k) = A_col%val(i)
            end do
        end do

        A%actual_nnz = A_col%N_nz

    end subroutine

    ! remove elements in the i_row-th row where i_row is in row_list
    subroutine MatrixTripletClearRow(A, row_list)
        use tools, only: assert
        type(MATRIX_TRIPLET), intent(inout) :: A
        integer, dimension(:), intent(in) :: row_list
        
        integer :: i_nz
        
        call assert(all(row_list > 0 .and. row_list <= A%N_row), "Invalid row index")
        
        do i_nz = 1, A%actual_nnz
            if (any(A%row_idx(i_nz) == row_list)) then
                A%val(i_nz) = 0.0
            end if
        end do

    end subroutine MatrixTripletClearRow
    
    ! remove nondiagonal elements in colth column where col is in col_list
    ! subtract A(row, col)*x(col) from b(row) 
    ! store the removed elements in Ae
    function MatrixTripletEliminateColumn(A, col_list, x, b) result(Ae)
        type(MATRIX_TRIPLET), intent(inout) :: A
        integer, dimension(:), intent(in) :: col_list
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(:), intent(inout) :: b
        type(MATRIX_TRIPLET) :: Ae
        
        integer :: n_col, col, row
        
        integer :: i
        
        
        n_col = size(col_list)
        
        
        call MatrixTripletInit(Ae, A%N_row, A%N_col, A%N_row*n_col)
        
        do i = 1, A%actual_nnz
            if (any(A%col_idx(i) == col_list)) then
                col = A%col_idx(i)
                row = A%row_idx(i)
                if(col/=row) then
                    call MatrixTripletAddValue(Ae, row, col, A%val(i))
                    b(row) = b(row) - A%val(i)*x(col)
                    A%val(i) = 0.0
                end if
            end if
        end do
        
        
    end function

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

    ! set the value of the matrix at (i,j)
    ! the non-zero element must exist (TODO: allow adding new elements)
    subroutine MatrixColumnSet(A,i,j,value)
        implicit none
        type(MATRIX_COLUMN), intent(inout) :: A
        integer, intent(in) :: i,j
        real(8), intent(in) :: value

        integer :: k

        logical :: changed = .false.

        do k = A%col_ptr(j), A%col_ptr(j+1)-1
            if(A%row_idx(k) == i) then
                A%val(k) = value
                changed = .true.
                return
            end if
        end do

        if(.not. changed) then
            print *, "Error: Element (",i,",",j,") does not exist"
            stop
        end if
    end subroutine MatrixColumnSet

    subroutine MatrixColumnPrint(A)
        implicit none
        type(MATRIX_COLUMN), intent(in) :: A

        integer :: i_col
        integer :: i_nz,i

        print *, "N_row: ",A%N_row, "N_col: ",A%N_col, "N_nz: ",A%N_nz

        i_nz = 0

        do i_col = 1, A%N_col
            do i = A%col_ptr(i_col), A%col_ptr(i_col+1)-1
                i_nz = i_nz + 1
                write(*,"(i9,' (',i9,i9,'): ',E11.5)") i_nz, A%row_idx(i), i_col, A%val(i)
            end do
        end do

    end subroutine

    ! remove zero elements in a matrix in column format
    subroutine MatrixColumnTrim(A)
        type(MATRIX_COLUMN), intent(inout) :: A
        
        real(8),dimension(:),allocatable :: val_tmp
        integer,dimension(:),allocatable :: row_idx_tmp,col_ptr_tmp

        real(8) :: tol = 1.0e-12

        logical, dimension(:), allocatable :: mask
        integer :: nnz_new, i_col, nnz_col

        allocate(val_tmp(A%N_nz),row_idx_tmp(A%N_nz),col_ptr_tmp(A%N_col+1))
        val_tmp = A%val
        row_idx_tmp = A%row_idx
        col_ptr_tmp = A%col_ptr
        allocate(mask(A%N_nz))
        mask = abs(A%val) > tol
        nnz_new = count(mask)

        deallocate(A%val,A%row_idx,A%col_ptr)
        allocate(A%val(nnz_new),A%row_idx(nnz_new),A%col_ptr(A%N_col+1))

        A%val = pack(val_tmp,mask)
        A%row_idx = pack(row_idx_tmp,mask)

        A%col_ptr(1) = 1
        do i_col = 1, A%N_col
            nnz_col = count(mask(col_ptr_tmp(i_col):col_ptr_tmp(i_col+1)-1))
            A%col_ptr(i_col+1) = A%col_ptr(i_col) + nnz_col
        end do

        A%N_nz = nnz_new

        deallocate(val_tmp,row_idx_tmp,col_ptr_tmp,mask)

    end subroutine MatrixColumnTrim

    ! remove elements in the i_row-th row
    subroutine MatrixColumnClearRow(A, i_row)
        type(MATRIX_COLUMN), intent(inout) :: A
        integer, intent(in) :: i_row

        where(A%row_idx == i_row)
            A%val = 0d0
        end where

    end subroutine MatrixColumnClearRow
    
    function MatrixColumnTranspose(A) result(A_t)
        type(MATRIX_COLUMN), intent(in) :: A
        type(MATRIX_COLUMN) :: A_t
        
        type(MATRIX_TRIPLET) :: A_triplet
        integer, dimension(:), allocatable :: tmp_row_idx
                
        call MatrixColumn2Triplet(A, A_triplet)
        
        allocate(tmp_row_idx(A_triplet%N_nz))
        tmp_row_idx = A_triplet%row_idx
        A_triplet%row_idx = A_triplet%col_idx
        A_triplet%col_idx = tmp_row_idx
        
        call MatrixTriplet2Column(A_triplet, A_t)
        
    end function

    subroutine MatrixColumnFree(A)
        implicit none
        type(MATRIX_COLUMN), intent(inout) :: A

        deallocate(A%row_idx, A%col_ptr, A%val)

        A%N_row = 0
        A%N_col = 0
        A%N_nz = 0

    end subroutine

    subroutine MatrixTripletPrint(A)
        implicit none
        type(MATRIX_TRIPLET), intent(in) :: A

        integer :: i

        print *, "N_row: ",A%N_row, "N_col: ",A%N_col, "N_nz: ",A%N_nz, "actual_nnz: ",A%actual_nnz

        do i = 1, A%N_nz
            write(*,"(i9,' (',i9,i9,'): ',E11.5)") i, A%row_idx(i), A%col_idx(i), A%val(i)
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
    
    ! y += alpha*A*x
    subroutine AddMultMV(alpha, A, x, y) 
        implicit none
        real(8), intent(in) :: alpha
        type(MATRIX_COLUMN), intent(in) :: A
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(:), intent(inout) :: y

        integer :: i, j

        do j = 1, A%N_col
            do i = A%col_ptr(j), A%col_ptr(j+1)-1
                y(A%row_idx(i)) = y(A%row_idx(i)) + alpha*A%val(i)*x(j)
            end do
        end do

    end subroutine


    
end module matvec