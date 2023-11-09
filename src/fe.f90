module fe
    use mesh
    use settings
    use fespace_P1
    use fespace_Q1
    use fespace_P2
    implicit none
    
contains

    subroutine fespaceInit(Vh, Th, basis_type, dim)
        type(fespace), intent(out) :: Vh
        type(mesh2D), intent(in) :: Th
        integer, intent(in) :: basis_type, dim

        select case(basis_type)
            case(DOF_P1)
                call fespaceInit_P1(Vh,Th,dim)
            case(DOF_Q1)
                call fespaceInit_Q1(Vh,Th,dim)
            case(DOF_P2)
                call fespaceInit_P2(Vh,Th,dim)
        end select

    end subroutine fespaceInit

    subroutine fespaceFree(Vh)
        type(fespace), intent(inout) :: Vh

        deallocate(Vh%ElemDOF)

    end subroutine fespaceFree

    subroutine Interpolate(u, fun, Th, Vh)
        type(fespace), intent(in) :: Vh
        type(mesh2D), intent(in) :: Th
        real(8), intent(out), dimension(:), allocatable :: u
        procedure(func) :: fun

        integer :: i_dof
        
        allocate(u(Vh%N_DOF))
        u = 0.d0
        
        do i_dof = 1, Vh%N_DOF
            call ComputeDof(u(i_dof), fun, Th, Vh, i_dof)
        end do
    end subroutine Interpolate

    subroutine ComputeDof(result,fun,Th,Vh,i_dof)
        type(mesh2D) :: Th
        type(fespace) :: Vh
        procedure(func) :: fun
        real(8), intent(out) :: result
        integer,intent(in) :: i_dof

        select case(Vh%basis_type)
            case(DOF_P1)
                call ComputeDof_P1(result,fun,Th,i_dof)
            case(DOF_Q1)
                call ComputeDof_Q1(result,fun,Th,i_dof)
            case(DOF_P2)
                call ComputeDof_P2(result,fun,Th,i_dof)
            case default
                print *, "Error: Unknown basis type"
                stop
        end select
    end subroutine ComputeDof

    

    subroutine getLocalDof(u, Vh, i_elem, i_dim, local_dof)
        real(8), dimension(:), intent(in) :: u
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, i_dim
        real(8), dimension(:), intent(out), allocatable :: local_dof

        integer, dimension(:), allocatable :: idx_local_dof
        
        integer :: N_local_dof

        N_local_dof = Vh%N_local_basis
        allocate(local_dof(N_local_dof))
        allocate(idx_local_dof(N_local_dof))
        idx_local_dof = Vh%ElemDOF((i_dim-1)*Vh%N_local_basis + 1:(i_dim-1)*Vh%N_local_basis + Vh%N_local_basis,i_elem)
        local_dof = u(idx_local_dof)
    end subroutine

    subroutine getEdgeDofIndex(Th, Vh, i_edge, dof_index)
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_edge
        integer, dimension(:), intent(out), allocatable :: dof_index

        select case (Vh%basis_type)
        case (DOF_P1)
            call getEdgeDofIndex_P1(Th, Vh, i_edge, dof_index)
        case (DOF_Q1)
            call getEdgeDofIndex_Q1(Th, Vh, i_edge, dof_index)
        case (DOF_P2)
            call getEdgeDofIndex_P2(Th, Vh, i_edge, dof_index)
        case default
            print *, "Error: Unknown basis type"
            stop
        end select
    end subroutine getEdgeDofIndex

    subroutine BasisLocal2D(refpts, Th, Vh, i_elem, deriv_type, result)
        real(8), intent(in), dimension(:,:) :: refpts
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, deriv_type
        real(8), intent(out), dimension(:,:,:),allocatable :: result

        integer :: num_pts

        num_pts = size(refpts, 2)
        allocate(result(Vh%dim, Vh%N_local_basis, num_pts))

        select case (Vh%basis_type)
        case (DOF_P1)
            call BasisLocalP1(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_Q1)
            call BasisLocalQ1(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_P2)
            call BasisLocalP2(refpts, Th, Vh, i_elem, deriv_type, result)
        case default
            print *, 'BasisLocal2D: Unknown basis type'
            stop
        end select
    end subroutine BasisLocal2D

end module fe