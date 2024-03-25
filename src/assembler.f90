module assembler
        use settings
        use fe_utils
        use matvec
        use fe, only: getLocalDofIndex
        implicit none
        
    contains
    
        subroutine AssembleMatrixElement(coe_fun, coe_num_list, Th, Vh_trial, start_trial, Vh_test, start_test,&
        assemble_info, Gauss_type, A)
            procedure(func) :: coe_fun
            real(8), dimension(:), intent(in) :: coe_num_list
            type(MESH2D),intent(in) :: Th
            type(FESPACE),intent(in) :: Vh_trial, Vh_test
            integer, intent(in) :: start_trial, start_test !start position of trial and test basis
            integer, intent(in) :: Gauss_type
            integer,dimension(:,:),intent(in) :: assemble_info
            ! asssemble_info: coe_fun_dim, trial_dim, trial_derive_type, test_dim, test_derive_type
            type(MATRIX_TRIPLET) :: A
    
            integer :: i_elem, i_assemble
            integer :: coe_fun_dim, trial_dim, trial_derive_type, test_dim, test_derive_type
            real(8) :: coe_num
    
            integer,dimension(:),allocatable :: ind_trial, ind_test
            real(8),dimension(:,:),allocatable :: local_mat

            allocate(ind_trial(Vh_trial%N_local_basis),ind_test(Vh_test%N_local_basis))
    
            ! assemble matrix
            do i_assemble = 1, size(assemble_info,1)
                coe_num = coe_num_list(i_assemble)
                coe_fun_dim = assemble_info(i_assemble,1)
                trial_dim = assemble_info(i_assemble,2)
                trial_derive_type = assemble_info(i_assemble,3)
                test_dim = assemble_info(i_assemble,4)
                test_derive_type = assemble_info(i_assemble,5)
    
                do i_elem = 1, Th%N_Elem
                    call getLocalDofIndex(Vh_trial, i_elem, trial_dim, ind_trial)
                    call getLocalDofIndex(Vh_test, i_elem, test_dim, ind_test)
    
                    ind_trial = ind_trial + start_trial
                    ind_test = ind_test + start_test
            
                    call LocalMatrix(i_elem, coe_fun, coe_fun_dim, Th, Vh_trial, Vh_test,&
                     trial_dim, trial_derive_type, test_dim, test_derive_type, Gauss_type, local_mat)

                    call MatrixTripletAddMatrixValues(A, ind_test, ind_trial, coe_num*local_mat)

                    deallocate(local_mat)
                end do
    
            end do
    
        end subroutine AssembleMatrixElement


        subroutine AssembleMatrixElementFE(coe_fun, coe_num_list, Th, Vh_trial, start_trial, Vh_test, start_test,&
        u, Vh_u, assemble_info, Gauss_type, A)
            procedure(func) :: coe_fun
            real(8), dimension(:), intent(in) :: coe_num_list
            type(mesh2D),intent(in) :: Th
            type(fespace),intent(in) :: Vh_trial, Vh_test
            integer, intent(in) :: start_trial, start_test !start position of trial and test basis
            integer, intent(in) :: Gauss_type
            integer,dimension(:,:),intent(in) :: assemble_info
            ! asssemble_info: coe_fun_dim, trial_dim, trial_derive_type, test_dim, test_derive_type,
            !                 u_dim, u_derive_type
            real(8), dimension(:), intent(in) :: u
            type(fespace),intent(in) :: Vh_u
            type(MATRIX_TRIPLET) :: A
    
            integer :: i_elem, i_assemble
            integer :: coe_fun_dim, trial_dim, trial_derive_type, test_dim, test_derive_type
            integer :: u_dim, u_derive_type
            real(8) :: coe_num
    
    
            integer,dimension(:),allocatable :: ind_trial, ind_test
            real(8),dimension(:,:),allocatable :: local_mat
    
            allocate(ind_trial(Vh_trial%N_local_basis),ind_test(Vh_test%N_local_basis))
    
            ! assemble matrix
            do i_assemble = 1, size(assemble_info,1)
                coe_num = coe_num_list(i_assemble)
                coe_fun_dim = assemble_info(i_assemble,1)
                trial_dim = assemble_info(i_assemble,2)
                trial_derive_type = assemble_info(i_assemble,3)
                test_dim = assemble_info(i_assemble,4)
                test_derive_type = assemble_info(i_assemble,5)
                u_dim = assemble_info(i_assemble,6)
                u_derive_type = assemble_info(i_assemble,7)
    
                do i_elem = 1, Th%N_Elem
                    call getLocalDofIndex(Vh_trial, i_elem, trial_dim, ind_trial)
                    call getLocalDofIndex(Vh_test, i_elem, test_dim, ind_test)
    
                    ind_trial = ind_trial + start_trial
                    ind_test = ind_test + start_test
            
                    call LocalMatrixFE(i_elem, coe_fun, coe_fun_dim, Th, Vh_trial, Vh_test,&
                        trial_dim, trial_derive_type, test_dim, test_derive_type, &
                        u, Vh_u, u_dim, u_derive_type, Gauss_type, local_mat)
                    call MatrixTripletAddMatrixValues(A, ind_test, ind_trial, coe_num*local_mat)
                    deallocate(local_mat)
                end do
    
            end do
    
        end subroutine AssembleMatrixElementFE
    
        subroutine AssembleVectorElement(coe_fun, coe_num_list, Th, Vh_test, start_test,&
        assemble_info, Gauss_type, b)
            procedure(func) :: coe_fun
            real(8), dimension(:), intent(in) :: coe_num_list
            type(mesh2D),intent(in) :: Th
            type(fespace),intent(in) :: Vh_test
            integer, intent(in) :: start_test !start position of test basis
            integer, intent(in) :: Gauss_type
            integer,dimension(:,:),intent(in) :: assemble_info
            ! asssemble_info: coe_fun_dim,  test_dim, test_derive_type
            real(8), dimension(:) :: b
    
    
            integer :: i_elem, i_assemble
            integer :: coe_fun_dim, test_dim, test_derive_type
            real(8) :: coe_num
    
            integer,dimension(:),allocatable :: ind_test
            real(8),dimension(:),allocatable :: local_vec
    
            allocate(ind_test(Vh_test%N_local_basis))
    
            do i_assemble = 1, size(assemble_info,1)
                coe_num = coe_num_list(i_assemble)
                coe_fun_dim = assemble_info(i_assemble,1)
                test_dim = assemble_info(i_assemble,2)
                test_derive_type = assemble_info(i_assemble,3)
    
                do i_elem = 1, Th%N_Elem
                    call getLocalDofIndex(Vh_test, i_elem, test_dim, ind_test)
                    ind_test = ind_test + start_test
                    
                    call LocalVector(i_elem, coe_fun, coe_fun_dim, Th, Vh_test,&
                     test_dim, test_derive_type, Gauss_type, local_vec)

                    b(ind_test) = b(ind_test) + coe_num*local_vec
                    deallocate(local_vec)
                end do
    
            end do
    
        end subroutine AssembleVectorElement


        subroutine AssembleVectorElementFE(coe_fun, coe_num_list, Th, Vh_test, start_test,&
        u, Vh_u, assemble_info, Gauss_type, b)
            procedure(func) :: coe_fun
            real(8), dimension(:), intent(in) :: coe_num_list
            type(mesh2D),intent(in) :: Th
            type(fespace),intent(in) :: Vh_test
            integer, intent(in) :: start_test !start position of test basis
            integer, intent(in) :: Gauss_type
            integer,dimension(:,:),intent(in) :: assemble_info
            ! asssemble_info: coe_fun_dim,  test_dim, test_derive_type
            !                u_dim, u_derive_type
            real(8), dimension(:), intent(in) :: u
            type(fespace),intent(in) :: Vh_u
            real(8), dimension(:) :: b
    
    
            integer :: i_elem, i_assemble
            integer :: coe_fun_dim, test_dim, test_derive_type
            integer :: u_dim, u_derive_type
            real(8) :: coe_num
    
            integer,dimension(:),allocatable :: ind_test
            real(8),dimension(:),allocatable :: local_vec
    
            allocate(ind_test(Vh_test%N_local_basis))
    
            do i_assemble = 1, size(assemble_info,1)
                coe_num = coe_num_list(i_assemble)
                coe_fun_dim = assemble_info(i_assemble,1)
                test_dim = assemble_info(i_assemble,2)
                test_derive_type = assemble_info(i_assemble,3)
                u_dim = assemble_info(i_assemble,4)
                u_derive_type = assemble_info(i_assemble,5)
    
                do i_elem = 1, Th%N_Elem
                    call getLocalDofIndex(Vh_test, i_elem, test_dim, ind_test)
                    ind_test = ind_test + start_test
                    
                    call LocalVectorFE(i_elem, coe_fun, coe_fun_dim, Th, Vh_test,&
                        test_dim, test_derive_type, u, Vh_u, u_dim, u_derive_type, Gauss_type, local_vec)
                    b(ind_test) = b(ind_test) + coe_num*local_vec
                    deallocate(local_vec)
                end do
    
            end do
    
        end subroutine AssembleVectorElementFE

        subroutine AssembleVectorBdry(bdry_marker, coe_fun, coe_num_list, Th, Vh_test, start_test, assemble_info, &
        Gauss_type, b)
            procedure(func) :: coe_fun
            real(8), dimension(:), intent(in) :: coe_num_list
            type(mesh2D),intent(in) :: Th
            type(fespace),intent(in) :: Vh_test
            integer, intent(in) :: bdry_marker !bdry marker
            integer, intent(in) :: start_test !start position of test basis
            integer, intent(in) :: Gauss_type
            integer,dimension(:,:),intent(in) :: assemble_info
            ! asssemble_info: coe_fun_dim,  test_dim, test_derive_type
            real(8), dimension(:) :: b
    
            integer :: i_elem, i_assemble, i_bdry, i_edge
            integer :: coe_fun_dim, test_dim, test_derive_type
            real(8) :: coe_num
    
            integer,dimension(:),allocatable :: ind_test
            real(8),dimension(:),allocatable :: local_vec
    
            allocate(ind_test(Vh_test%N_local_basis))
    
            do i_assemble = 1, size(assemble_info,1)
                coe_num = coe_num_list(i_assemble)
                coe_fun_dim = assemble_info(i_assemble,1)
                test_dim = assemble_info(i_assemble,2)
                test_derive_type = assemble_info(i_assemble,3)
    
                do i_bdry = 1, Th%N_bdryedge
                    if (Th%BdryMarker(i_bdry) /= bdry_marker) cycle
                    i_edge = Th%BdryEdge(i_bdry)
                    i_elem = abs(Th%EdgeElemConn(2,i_edge))
    
                    call getLocalDofIndex(Vh_test, i_elem, test_dim, ind_test)
                    ind_test = ind_test + start_test
                    
                    call LocalVectorLine(i_edge, i_elem, coe_fun, coe_fun_dim, Th, Vh_test,&
                        test_dim, test_derive_type, Gauss_type, local_vec)

                    b(ind_test) = b(ind_test) + coe_num*local_vec

                    deallocate(local_vec)
                end do
    
            end do
    
        end subroutine AssembleVectorBdry
        
    end module assembler