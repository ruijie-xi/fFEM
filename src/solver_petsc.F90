module solver_petsc
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscpcdef.h>
      use petscvec
      use petscmat
      use petscksp
      use petscpc
      use matvec
    implicit none
    
contains

    subroutine SolverSolvePETSC(A, b, x)
        type(MATRIX_TRIPLET) :: A
        real(8), dimension(:) :: b
        real(8), dimension(:) :: x

        type(MATRIX_COLUMN) :: A_column

        integer(4) :: ierr
        Mat :: matA
        Vec :: vecb, vecx
        KSP :: ksp
        PC :: pc
        PetscScalar, pointer :: xx_v(:)

        call PetscInitialize("input/petsc_options.dat", ierr)
        call CreateMat(matA, A%N_row, A%N_col, ierr)
        call CreateVec(vecb, size(b), ierr)
        call CreateVec(vecx, size(x), ierr)
        call CreateSolver(ksp, matA, matA, ierr)

        call MatrixTriplet2Column(A, A_column)
        call MatrixTripletFree(A)
        call MatrixColumn2Triplet(A_column, A)

        call MatrixTriplet2PETSC(A, matA, ierr)
        call Vector2PETSC(b, vecb, ierr)

        call Solve(matA, vecb, ksp, pc, vecx, ierr)

        call VecGetArrayF90(vecx, xx_v, ierr)

        x = xx_v


    end subroutine  SolverSolvePETSC

    subroutine MatrixTriplet2PETSC(A, matA, ierr)
        implicit none
        type(MATRIX_TRIPLET) :: A
        Mat :: matA
        integer(4) :: ierr

        integer :: i_nz

        do i_nz = 1, A%N_nz
            call MatSetValue(matA, A%row_idx(i_nz)-1, A%col_idx(i_nz)-1, A%val(i_nz), INSERT_VALUES, ierr)
        end do
    end subroutine MatrixTriplet2PETSC

    subroutine Vector2PETSC(b, vecb, ierr)
        implicit none
        real(8), dimension(:) :: b
        Vec :: vecb
        integer(4) :: ierr

        PetscScalar, pointer :: b_v(:)

        call VecGetArrayF90(vecb, b_v, ierr)
        b_v = b
        call VecRestoreArrayF90(vecb, b_v, ierr)
    end subroutine Vector2PETSC


    subroutine CreateVec(vec, size_global, errpetsc)
        Vec :: vec
        integer :: size_global, errpetsc
        call VecCreate(PETSC_COMM_WORLD, vec, errpetsc)
        CHKERRQ(errpetsc)
        call VecSetSizes(vec, PETSC_DECIDE, size_global, errpetsc)
        CHKERRQ(errpetsc)
        call VecSetFromOptions(vec, errpetsc)
        CHKERRQ(errpetsc)
        call VecSetOption(vec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, errpetsc)

    end subroutine CreateVec

    subroutine CreateMat(mat, size_row, size_col, errpetsc)
        Mat :: mat
        integer :: size_row, size_col, errpetsc
        integer :: dummy = 50
        integer,dimension(:),allocatable :: nnzVec

        allocate(nnzVec(size_row))
        nnzVec = 50

        call MatCreate(PETSC_COMM_WORLD, mat, errpetsc)
        CHKERRQ(errpetsc)

        call MatSetSizes(mat, size_row, size_col, size_row, size_col, errpetsc)
        CHKERRQ(errpetsc)

        call MatSetFromOptions(mat, errpetsc)
        CHKERRQ(errpetsc)

        call MatMPIAIJSetPreallocation(mat, dummy, nnzVec, dummy, nnzVec, errpetsc)
        CHKERRQ(errpetsc)

        call MatSeqAIJSetPreallocation(mat, dummy, nnzVec, errpetsc)
        CHKERRQ(errpetsc)

        call MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, errpetsc)
        CHKERRQ(errpetsc)

        call MatSetOption(mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, errpetsc)
        CHKERRQ(errpetsc)

        call MatSetOption(mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, errpetsc)
        CHKERRQ(errpetsc)


    end subroutine CreateMat

    subroutine CreateSolver(solver, mat, Pmat, errpetsc)
        KSP :: solver
        PC :: pc
        Mat :: mat
        Mat :: Pmat
        integer :: errpetsc

      ! Create the KSP context
      call KSPCreate(PETSC_COMM_WORLD, solver, errpetsc)
      CHKERRQ(errpetsc)
      ! Set the operators for the KSP context
      call KSPSetOperators(solver, mat, Pmat, errpetsc)
      CHKERRQ(errpetsc)
      ! Set the KSP type
      call KSPSetType(solver, KSPCG, errpetsc)
      ! call KSPSetType(this%ksp, KSPBCGS, errpetsc)
      CHKERRQ(errpetsc)

      ! ! Set whether to use non-zero initial guess or not
      ! !KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
      ! !KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

      ! Set KSP options from the input file
      ! This is convenient as it allows to choose different options
      ! from the input files instead of recompiling the code
      call KSPSetFromOptions(solver, errpetsc)
      CHKERRQ(errpetsc)

      ! Get the PC context
      call KSPGetPC(solver, pc, errpetsc)
      CHKERRQ(errpetsc)
      ! Set the PC context
      call PCSetType(pc, PCBJACOBI, errpetsc)
      CHKERRQ(errpetsc)
      ! Set PC options from the input file
      call PCSetFromOptions(pc, errpetsc)
      CHKERRQ(errpetsc)
    end subroutine CreateSolver

    subroutine SetZeroMat(mat, errpetsc)
        Mat :: mat
        integer :: errpetsc
        call MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY, errpetsc)
        CHKERRQ(errpetsc)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY, errpetsc)
        CHKERRQ(errpetsc)
        call MatZeroEntries(mat, errpetsc)
        CHKERRQ(errpetsc)
    end subroutine SetZeroMat

    subroutine AssembleMat(mat, errpetsc)
        Mat :: mat
        integer :: errpetsc
        call MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY, errpetsc)
        CHKERRQ(errpetsc)
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY, errpetsc)
        CHKERRQ(errpetsc)
    end subroutine AssembleMat

    subroutine AssembleVec(vec, errpetsc)
        Vec :: vec
        integer :: errpetsc
        call VecAssemblyBegin(vec, errpetsc)
        CHKERRQ(errpetsc)
        call VecAssemblyEnd(vec, errpetsc)
        CHKERRQ(errpetsc)
    end subroutine AssembleVec

    subroutine SetZeroVec(vec, errpetsc)
        Vec :: vec
        integer :: errpetsc
        call VecAssemblyBegin(vec, errpetsc)
        CHKERRQ(errpetsc)
        call VecAssemblyEnd(vec, errpetsc)
        CHKERRQ(errpetsc)
        call VecZeroEntries(vec, errpetsc)
        CHKERRQ(errpetsc)
    end subroutine SetZeroVec

    subroutine Solve(mat,vec,solver,pc,sol,errpetsc)
        Mat :: mat
        Vec :: vec, sol
        KSP :: solver
        PC :: pc
        integer :: errpetsc

        CHARACTER (LEN=100) :: charTemp

        KSPConvergedReason reason;
        PetscInt its;

        call MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY, errpetsc)
        CHKERRQ(errpetsc)
        call MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY, errpetsc)
        CHKERRQ(errpetsc)

        ! Assemble the solnVec
        call VecAssemblyBegin(sol, errpetsc)
        CHKERRQ(errpetsc)
        call VecAssemblyEnd(sol, errpetsc)
        CHKERRQ(errpetsc)
        call VecZeroEntries(sol, errpetsc)
        CHKERRQ(errpetsc)

        ! Assemble the rhsVec 
        call VecAssemblyBegin(vec, errpetsc)
        CHKERRQ(errpetsc)
        call VecAssemblyEnd(vec, errpetsc)
        CHKERRQ(errpetsc)

        call PetscPrintf(PETSC_COMM_WORLD, " Solving the matrix system \n", errpetsc)

        ! Solve the matrix system
        call KSPSolve(solver, vec, sol, errpetsc)
        CHKERRQ( errpetsc)

        call KSPGetConvergedReason(solver, reason, errpetsc)

        IF(reason < 0) THEN
            call PetscPrintf(PETSC_COMM_WORLD, "Divergence.\n", errpetsc)
        ELSE
            call KSPGetIterationNumber(solver, its, errpetsc)

        WRITE(charTemp,*) "Convergence in", its, " iterations.", "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
        END IF

    end subroutine Solve

end module solver_petsc