program test_matvec
    use settings
    use matvec
    implicit none
    
    ! matrix:
    ! [1, 0, 0, 2, 0]   [0, 3, 0, 0, 4]
    ! [0, 3, 0, 0, 4]   [0, 0, 0, 0, 0]
    ! [0, 0, 5, 0, 0] + [0, 0, 6, 0, 0]
    ! [6, 0, 0, 7, 0]   [0, 0, 0, 0, 9]
    ! [0, 8, 0, 0, 9]   [0, 0, 0, 0, 0]
    
    type(MATRIX_TRIPLET) :: A1, A2, A3, A3_triplet
    type(MATRIX_COLUMN) :: A3_column

    real(8) :: val

    call MatrixTripletInit(A1, 5, 5, 9)
    call MatrixTripletInit(A2, 5, 5, 4)

    A1%val = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    A1%row_idx = [1, 1, 2, 2, 3, 4, 4, 5, 5]
    A1%col_idx = [1, 4, 2, 5, 3, 1, 4, 2, 5]

    A2%val = [3, 4, 6, 9]
    A2%row_idx = [1, 1, 3, 4]
    A2%col_idx = [2, 5, 3, 5]

    call MatrixTripletAdd(1d0, A1, 1d0, A2, A3)

    call MatrixTriplet2Column(A3, A3_column)

    call MatrixColumnPrint(A3_column)

    call MatrixColumnGet(A3_column, 1, 1, val)
    print *, "the value at (1, 1) is ", val
    call MatrixColumnGet(A3_column, 3, 3, val)
    print *, "the value at (3, 3) is ", val
    call MatrixColumnGet(A3_column, 3, 5, val)
    print *, "the value at (3, 5) is ", val

    call MatrixColumn2Triplet(A3_column, A3_triplet)
    call MatrixTripletPrint(A3_triplet)

end program test_matvec