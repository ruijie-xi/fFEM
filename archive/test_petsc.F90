program test_petsc
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
    use petscvec
    use petscmat
    use solver
    implicit none
  
    PetscInt ierr
    Vec vec
    Mat A
    real(8),dimension(2,2) :: local_mat

    call PetscInitialize(ierr)

    call CreateMat(A,100,100,ierr)

    local_mat(1,:) = [1,2]
    local_mat(2,:) = [3,4]

    call MatSetValues(A, 2, [6,8], 2 , [3,4] , transpose(local_mat), INSERT_VALUES, ierr)

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)
    call PetscFinalize(ierr)
  
  end program test_petsc