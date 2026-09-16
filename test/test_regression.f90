program test_regression
    use settings
    use matvec
    use quadrature
    use mesh, only: MeshInit, AddBdryMarker
    use mesh_generator, only: TriangleMesh, SquareMesh
    use fe, only: fespaceInit, Interpolate
    use fe_utils, only: ComputeError, FEfunctionGetValue, FEfunctionQuadValueLine
    use solver_umfpack2
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    type(MATRIX_TRIPLET) :: a, b, c
    type(MATRIX_COLUMN) :: ac, at, att
    type(VECTOR) :: rhs, sol, copied, u
    type(MESH2D) :: th
    type(FESPACE) :: vh
    real(8), allocatable :: pts(:,:), weights(:), nodes(:,:), value(:,:)
    integer, allocatable :: elems(:,:)
    real(8) :: scalar, err
    integer :: i, typ, rule, status
    integer, parameter :: types(10) = [DOF_P0, DOF_P1, DOF_P2, DOF_DG1, DOF_Q0, &
        DOF_Q1, DOF_Q2, DOF_QuadNedelec1, DOF_QuadRT1, DOF_QuadRT2]
    integer, parameter :: tri_rules(4) = [TriangleNodeAverage, TrianglePt1, TrianglePt4, TrianglePt9]
    character(32) :: mode

    call get_command_argument(1, mode)
    if (trim(mode) == 'singular-fatal') then
        call singular_matrix()
        call SolverSolveUMFPACK2(ac, rhs, sol)
        ! The Python runner checks both nonzero exit and the solver diagnostic.
        error stop 'Unexpected return from failed solver'
    end if

    ! Conversion preserves small entries, sums duplicates, and removes exact zeros.
    call MatrixTripletInit(a, 2, 3, 4)
    call MatrixTripletAddValues(a, [1,2,2,1], [1,3,3,2], [1d-14,3d0,4d0,0d0])
    call MatrixTriplet2Column(a, ac)
    call check(ac%N_nz == 2, 'conversion nonzero count')
    call MatrixColumnGet(ac, 1, 1, scalar)
    call check(scalar == 1d-14, 'small coefficients preserved')
    at = MatrixColumnTranspose(ac)
    call check(at%N_row == 3 .and. at%N_col == 2, 'rectangular transpose dimensions')
    call MatrixColumnGet(at, 3, 2, scalar)
    call check(scalar == 7d0, 'transpose value and duplicate summation')
    att = MatrixColumnTranspose(at)
    call check(all(att%col_ptr == ac%col_ptr), 'double transpose column pointers')
    call check(all(att%row_idx == ac%row_idx) .and. all(att%val == ac%val), 'double transpose entries')
    call MatrixColumnTrim(att, 1d-12)
    call check(att%N_nz == 1, 'explicit numerical drop tolerance')

    ! Unused triplet capacity must never become a matrix entry.
    call MatrixTripletInit(a, 2, 2, 4)
    call MatrixTripletAddValue(a, 1, 1, 1d0)
    call MatrixTripletInit(b, 2, 2, 4)
    call MatrixTripletAddValue(b, 2, 2, 2d0)
    call MatrixTripletAdd(2d0, a, 3d0, b, c)
    call check(c%actual_nnz == 2, 'triplet addition actual count')
    call MatrixTriplet2Column(c, ac)
    call MatrixColumnGet(ac, 2, 2, scalar)
    call check(scalar == 6d0, 'triplet addition value')
    call MatrixTripletInit(a, 2, 3, 0)
    call MatrixTriplet2Column(a, ac)
    at = MatrixColumnTranspose(ac)
    call check(at%N_nz == 0 .and. all(at%col_ptr == 1), 'empty transpose')

    ! Both public vector APIs use the same initialization and norm semantics.
    call rhs%Init(2)
    rhs%data = [1d0,2d0]
    call VectorCopy(copied, rhs)
    call check(copied%size == 2 .and. all(copied%data == rhs%data), 'copy to fresh vector')
    call copied%Reset(3)
    call check(copied%size == 3 .and. all(copied%data == 0d0), 'legacy vector reset')
    call VectorCopy(copied, rhs)
    call copied%AddVector(rhs, 2d0)
    call check(copied%Norm() == 6d0, 'legacy scaled addition and norm')
    call check(VectorNormLinf(copied) == copied%Norm(), 'norm API agreement')

    ! A nonsingular solve succeeds via both storage interfaces with uninitialized x.
    call MatrixTripletInit(a, 2, 2, 4)
    call MatrixTripletAddValues(a, [1,2,1,2], [1,1,2,2], [2d0,-1d0,-1d0,2d0])
    rhs%data = [0d0,3d0]
    call SolverSolveUMFPACK2(a, rhs, sol, status)
    call check(status == 0, 'triplet solve status')
    call check(maxval(abs(sol%data-[1d0,2d0])) < 1d-12, 'triplet solve residual')
    call MatrixTriplet2Column(a, ac)
    call VectorFree(sol)
    call SolverSolveUMFPACK2(ac, rhs, sol, status)
    call check(status == 0 .and. sol%size == 2, 'CSC solve initializes output')
    call check(maxval(abs(sol%data-[1d0,2d0])) < 1d-12, 'CSC solve residual')
    call MatrixTripletInit(a, 2, 2, 5)
    call MatrixTripletAddValues(a, [1,1,2,1,2], [1,1,1,2,2], [1d0,1d0,-1d0,-1d0,2d0])
    call SolverSolveUMFPACK2(a, rhs, sol, status)
    call check(status == 0, 'duplicate triplets are a valid solver input')
    call check(maxval(abs(sol%data-[1d0,2d0])) < 1d-12, 'duplicate triplet solve')
    call MatrixTripletInit(a, 1, 1, 1)
    call MatrixTripletAddValue(a, 1, 1, 1d-14)
    call VectorReset(rhs, 1)
    rhs%data = 2d-14
    call SolverSolveUMFPACK2(a, rhs, sol, status)
    call check(status == 0, 'one-by-one solver workspace')
    call check(abs(sol%data(1)-2d0) < 1d-12, 'scaled one-by-one solve')
    ! Singleton-only BTF blocks must not read uninitialized UMFPACK workspace flags.
    call MatrixTripletInit(a, 9, 9, 9)
    do i = 1, 9
        call MatrixTripletAddValue(a, i, i, 1d0)
    end do
    call VectorReset(rhs, 9)
    rhs%data = 1d0
    call SolverSolveUMFPACK2(a, rhs, sol, status)
    call check(status == 0, 'diagonal BTF factorization')
    call check(maxval(abs(sol%data-1d0)) < 1d-12, 'diagonal BTF solve')
    call singular_matrix()
    call SolverSolveUMFPACK2(ac, rhs, sol, status)
    call check(status /= 0, 'singular CSC failure status')
    call check(.not. any(ieee_is_finite(sol%data)), 'failed solution is invalidated')
    call SolverSolveUMFPACK2(a, rhs, sol, status)
    call check(status /= 0, 'singular triplet failure status')

    do i = 1, size(tri_rules)
        call getGaussRefElement(tri_rules(i), pts, weights)
        call check(abs(sum(weights)-0.5d0) < 1d-14, 'triangle constant integral')
        call check(abs(sum(weights*pts(1,:))-1d0/6d0) < 1d-14, 'triangle linear integral')
    end do

    ! Constant reproduction includes nonaffine quadrilaterals and shared-edge orientations.
    do i = 1, size(types)
        typ = types(i)
        if (i <= 4) then
            call TriangleMesh(0d0,1d0,0d0,1d0,2,2,elems,nodes)
            rule = TrianglePt9
        else
            call SquareMesh(0d0,1d0,0d0,1d0,2,2,elems,nodes)
            rule = QuadPt16
            where (abs(nodes(1,:)-0.5d0)<1d-10 .and. abs(nodes(2,:)-0.5d0)<1d-10)
                nodes(1,:) = nodes(1,:) + 0.07d0
                nodes(2,:) = nodes(2,:) - 0.04d0
            end where
        end if
        call MeshInit(elems,nodes,th)
        call AddBdryMarker(th, left_boundary, 7)
        call check(count(th%BdryMarker == 7) == 2, 'default boundary tolerance')
        call AddBdryMarker(th, left_boundary, 8, 1d-10)
        call check(count(th%BdryMarker == 8) == 2, 'explicit boundary tolerance')
        call fespaceInit(vh, th, typ, 2)
        call Interpolate(u, constant, th, vh)
        call ComputeError(constant,u,th,vh,NORM_L2,rule,err)
        call check(err < 1d-12, 'constant field reproduction')
        if (typ == DOF_QuadRT1 .or. typ == DOF_QuadRT2) then
            value = FEfunctionGetValue(u,th,vh,1,reshape([0.3d0,0.4d0],[2,1]),DERIV_DIV)
            call check(size(value,1) == 1, 'point divergence range dimension')
            call check(maxval(abs(value)) < 1d-11, 'constant point divergence')
            call FEfunctionQuadValueLine(u,th,vh,1,1,DERIV_DIV,LinePt4,value)
            call check(size(value,1) == 1, 'edge divergence range dimension')
            call check(maxval(abs(value)) < 1d-11, 'constant edge divergence')
        end if
    end do
    print *, 'PASS: algebra, solver failures, quadrature, vector API, boundary markers, and FE invariants'

contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine

    subroutine singular_matrix()
        call MatrixTripletInit(a, 2, 2, 4)
        call MatrixTripletAddValues(a, [1,2,1,2], [1,1,2,2], [1d0,1d0,1d0,1d0])
        call MatrixTriplet2Column(a, ac)
        call VectorInit(rhs, 2)
        rhs%data = [1d0,2d0]
    end subroutine

    subroutine constant(x, f, deriv_type)
        real(8), intent(in) :: x(:)
        real(8), allocatable, intent(out) :: f(:)
        integer, intent(in) :: deriv_type
        allocate(f(2))
        f = [1d0,2d0]
    end subroutine

    subroutine left_boundary(x, f, deriv_type)
        real(8), intent(in) :: x(:)
        real(8), allocatable, intent(out) :: f(:)
        integer, intent(in) :: deriv_type
        allocate(f(1))
        f(1) = x(1)
    end subroutine
end program
