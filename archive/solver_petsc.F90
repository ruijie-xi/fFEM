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
    implicit none
    
contains
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