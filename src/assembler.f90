module assembler
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
    use petscvec
    use petscmat
    use settings
    use fe_utils
    
    implicit none
    
contains

    subroutine AssembleMatrixElement(coe_fun, Th, Vh_trial, start_trial,  Vh_test, start_test,&
    assemble_info, Gauss_type, errpetsc, A)
        procedure(func) :: coe_fun
        type(mesh2D),intent(in) :: Th
        type(fespace),intent(in) :: Vh_trial, Vh_test
        integer, intent(in) :: start_trial, start_test !start position of trial and test basis
        integer, intent(in) :: Gauss_type
        integer,dimension(:,:),intent(in) :: assemble_info
        ! asssemble_info: coe_num, coe_fun_dim, trial_dim, trial_derive_type, test_dim, test_derive_type
        integer, intent(out) :: errpetsc
        Mat :: A

        integer :: i_elem, i_assemble
        integer :: coe_num, coe_fun_dim, trial_dim, trial_derive_type, test_dim, test_derive_type

        integer,dimension(:),allocatable :: ind_trial, ind_test
        real(8),dimension(:,:),allocatable :: local_mat

        allocate(ind_trial(Vh_trial%N_local_basis),ind_test(Vh_test%N_local_basis))

        ! matrix pattern
        do i_assemble = 1, size(assemble_info,1)
            coe_num = assemble_info(i_assemble,1)
            coe_fun_dim = assemble_info(i_assemble,2)
            trial_dim = assemble_info(i_assemble,3)
            trial_derive_type = assemble_info(i_assemble,4)
            test_dim = assemble_info(i_assemble,5)
            test_derive_type = assemble_info(i_assemble,6)

            do i_elem = 1, Th%N_Elem
                ind_trial = Vh_trial%ElemDOF((trial_dim-1)*Vh_trial%N_local_basis+1:trial_dim*Vh_trial%N_local_basis,i_elem)-1
                ind_test = Vh_test%ElemDOF((test_dim-1)*Vh_test%N_local_basis+1:test_dim*Vh_test%N_local_basis,i_elem)-1
                ind_trial = ind_trial + start_trial
                ind_test = ind_test + start_test
        
                allocate(local_mat(Vh_test%N_local_basis,Vh_trial%N_local_basis))
                local_mat = 0.0d0
                local_mat = transpose(local_mat)
                call MatSetValues(A, Vh_test%N_local_basis, ind_test, Vh_trial%N_local_basis,&
                 ind_trial, local_mat, INSERT_VALUES, errpetsc)
                deallocate(local_mat)
            end do
        end do

        ! assemble matrix
        do i_assemble = 1, size(assemble_info,1)
            coe_num = assemble_info(i_assemble,1)
            coe_fun_dim = assemble_info(i_assemble,2)
            trial_dim = assemble_info(i_assemble,3)
            trial_derive_type = assemble_info(i_assemble,4)
            test_dim = assemble_info(i_assemble,5)
            test_derive_type = assemble_info(i_assemble,6)

            do i_elem = 1, Th%N_Elem
                ind_trial = Vh_trial%ElemDOF((trial_dim-1)*Vh_trial%N_local_basis+1:trial_dim*Vh_trial%N_local_basis,i_elem)-1
                ind_test = Vh_test%ElemDOF((test_dim-1)*Vh_test%N_local_basis+1:test_dim*Vh_test%N_local_basis,i_elem)-1
                ind_trial = ind_trial + start_trial
                ind_test = ind_test + start_test
        
                call LocalMatrix(i_elem, coe_fun, coe_fun_dim, Th, Vh_trial, Vh_test,&
                 trial_dim, trial_derive_type, test_dim, test_derive_type, Gauss_type, local_mat)
                local_mat = transpose(local_mat)
                call MatSetValues(A, Vh_test%N_local_basis, ind_test, Vh_trial%N_local_basis,&
                 ind_trial, local_mat, ADD_VALUES, errpetsc)
                deallocate(local_mat)
            end do

        end do

        call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, errpetsc)
        call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, errpetsc)

    end subroutine AssembleMatrixElement

    subroutine AssembleVectorElement(coe_fun, Th, Vh_test, start_test, assemble_info, Gauss_type, errpetsc, b)
        procedure(func) :: coe_fun
        type(mesh2D),intent(in) :: Th
        type(fespace),intent(in) :: Vh_test
        integer, intent(in) :: start_test !start position of test basis
        integer, intent(in) :: Gauss_type
        integer,dimension(:,:),intent(in) :: assemble_info
        ! asssemble_info: coe_num, coe_fun_dim,  test_dim, test_derive_type
        integer, intent(out) :: errpetsc
        Vec :: b

        integer :: i_elem, i_assemble
        integer :: coe_num, coe_fun_dim, test_dim, test_derive_type

        integer,dimension(:),allocatable :: ind_test
        real(8),dimension(:),allocatable :: local_vec

        allocate(ind_test(Vh_test%N_local_basis))

        do i_assemble = 1, size(assemble_info,1)
            coe_num = assemble_info(i_assemble,1)
            coe_fun_dim = assemble_info(i_assemble,2)
            test_dim = assemble_info(i_assemble,3)
            test_derive_type = assemble_info(i_assemble,4)

            do i_elem = 1, Th%N_Elem
                ind_test = Vh_test%ElemDOF((test_dim-1)*Vh_test%N_local_basis+1:test_dim*Vh_test%N_local_basis,i_elem)-1
                ind_test = ind_test + start_test
                
                call LocalVector(i_elem, coe_fun, coe_fun_dim, Th, Vh_test,&
                 test_dim, test_derive_type, Gauss_type, local_vec)
                call VecSetValues(b, Vh_test%N_local_basis, ind_test, local_vec, ADD_VALUES, errpetsc)
                deallocate(local_vec)
            end do

        end do

        call VecAssemblyBegin(b, errpetsc)
        call VecAssemblyEnd(b, errpetsc)

    end subroutine AssembleVectorElement

    subroutine AssembleVectorBdry(bdry_marker, coe_fun, Th, Vh_test, start_test, assemble_info, Gauss_type, errpetsc, b)
        procedure(func) :: coe_fun
        type(mesh2D),intent(in) :: Th
        type(fespace),intent(in) :: Vh_test
        integer, intent(in) :: bdry_marker !bdry marker
        integer, intent(in) :: start_test !start position of test basis
        integer, intent(in) :: Gauss_type
        integer,dimension(:,:),intent(in) :: assemble_info
        ! asssemble_info: coe_num, coe_fun_dim,  test_dim, test_derive_type
        integer, intent(out) :: errpetsc
        Vec :: b

        integer :: i_elem, i_assemble, i_bdry, i_edge
        integer :: coe_num, coe_fun_dim, test_dim, test_derive_type

        integer,dimension(:),allocatable :: ind_test
        real(8),dimension(:),allocatable :: local_vec

        allocate(ind_test(Vh_test%N_local_basis))

        do i_assemble = 1, size(assemble_info,1)
            coe_num = assemble_info(i_assemble,1)
            coe_fun_dim = assemble_info(i_assemble,2)
            test_dim = assemble_info(i_assemble,3)
            test_derive_type = assemble_info(i_assemble,4)

            do i_bdry = 1, Th%N_bdryedge
                i_edge = Th%BdryEdge(i_bdry)
                i_elem = abs(Th%EdgeElemConn(2,i_edge))

                ind_test = Vh_test%ElemDOF((test_dim-1)*Vh_test%N_local_basis+1:test_dim*Vh_test%N_local_basis,i_elem)-1
                ind_test = ind_test + start_test
                
                call LocalVectorLine(i_edge, i_elem, coe_fun, coe_fun_dim, Th, Vh_test,&
                 test_dim, test_derive_type, Gauss_type, local_vec)
                call VecSetValues(b, Vh_test%N_local_basis, ind_test, local_vec, ADD_VALUES, errpetsc)
                deallocate(local_vec)
            end do

        end do

        call VecAssemblyBegin(b, errpetsc)
        call VecAssemblyEnd(b, errpetsc)

    end subroutine AssembleVectorBdry
    
end module assembler