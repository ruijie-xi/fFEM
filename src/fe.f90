module fe
    use mesh
    use settings
    use fespace_P0
    use fespace_P1
    use fespace_DG1
    use fespace_Q0
    use fespace_Q1
    use fespace_P2
    use fespace_QuadNedelec1
    use fespace_QuadRT1
    implicit none
    
contains

    subroutine fespaceInit(Vh, Th, basis_type, dim)
        type(fespace), intent(out) :: Vh
        type(mesh2D), intent(in) :: Th
        integer, intent(in) :: basis_type, dim

        select case(basis_type)
            case(DOF_P0)
                call fespaceInit_P0(Vh,Th,dim)
            case(DOF_P1)
                call fespaceInit_P1(Vh,Th,dim)
            case(DOF_Q0)
                call fespaceInit_Q0(Vh,Th,dim)
            case(DOF_Q1)
                call fespaceInit_Q1(Vh,Th,dim)
            case(DOF_P2)
                call fespaceInit_P2(Vh,Th,dim)
            case(DOF_DG1)
                call fespaceInit_DG1(Vh,Th,dim)
            case(DOF_QuadNedelec1)
                call fespaceInit_QuadNedelec1(Vh,Th,dim)
            case(DOF_QuadRT1)
                call fespaceInit_QuadRT1(Vh,Th,dim)
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
        
        if (allocated(u)) deallocate(u)
        allocate(u(Vh%N_DOF))
        u = 0.d0
        
        do i_dof = 1, Vh%N_DOF
            u(i_dof) = ComputeDof(fun, Th, Vh, i_dof)
        end do
    end subroutine Interpolate

    function ComputeDof(fun,Th,Vh,i_dof) result(result)
        type(mesh2D) :: Th
        type(fespace) :: Vh
        procedure(func) :: fun
        real(8) :: result
        integer,intent(in) :: i_dof

        select case(Vh%basis_type)
        case(DOF_P0)
            call ComputeDof_P0(result,fun,Th, Vh, i_dof)
        case(DOF_P1)
            call ComputeDof_P1(result,fun,Th,i_dof)
        case(DOF_Q0)
            call ComputeDof_Q0(result,fun,Th,Vh,i_dof)
        case(DOF_Q1)
            call ComputeDof_Q1(result,fun,Th,i_dof)
        case(DOF_P2)
            call ComputeDof_P2(result,fun,Th,i_dof)
        case(DOF_DG1)
            call ComputeDof_DG1(result,fun,Th,i_dof)
        case(DOF_QuadNedelec1)
            call ComputeDof_QuadNedelec1(result,fun,Th,i_dof)
        case(DOF_QuadRT1)
            call ComputeDof_QuadRT1(result,fun,Th,i_dof)
        case default
            print *, "Error: Unknown basis type"
            stop
        end select
    end function ComputeDof

    subroutine getLocalDofIndex(Vh, i_elem, i_dim, idx_local_dof)
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, i_dim

        integer, intent(out), dimension(:), allocatable :: idx_local_dof
        
        integer :: N_local_dof

        N_local_dof = Vh%N_local_basis
        if (allocated(idx_local_dof)) deallocate(idx_local_dof)
        allocate(idx_local_dof(N_local_dof))
        if (Vh%isStack.eq.1) then
            idx_local_dof = Vh%ElemDOF((i_dim-1)*Vh%N_local_basis + 1:(i_dim-1)*Vh%N_local_basis + Vh%N_local_basis,i_elem)
        elseif (Vh%isStack.eq.0) then
            idx_local_dof = Vh%ElemDOF(:,i_elem)
        end if
    end subroutine

    subroutine getLocalDof(u, Vh, i_elem, i_dim, local_dof)
        real(8), dimension(:), intent(in) :: u
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_elem, i_dim
        real(8), dimension(:), intent(out), allocatable :: local_dof

        integer, dimension(:), allocatable :: idx_local_dof
        
        integer :: N_local_dof

        N_local_dof = Vh%N_local_basis
        allocate(local_dof(N_local_dof))
        call getLocalDofIndex(Vh, i_elem, i_dim, idx_local_dof)
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
        case (DOF_DG1)
            call getEdgeDofIndex_DG1(Th, Vh, i_edge, dof_index)
        case (DOF_QuadNedelec1)
            call getEdgeDofIndex_QuadNedelec1(Th, Vh, i_edge, dof_index)
        case (DOF_QuadRT1)
            call getEdgeDofIndex_QuadRT1(Th, Vh, i_edge, dof_index)
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
        case (DOF_P0)
            call BasisLocalP0(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_P1)
            call BasisLocalP1(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_Q0)
            call BasisLocalQ0(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_Q1)
            call BasisLocalQ1(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_P2)
            call BasisLocalP2(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_DG1)
            call BasisLocalDG1(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_QuadNedelec1)
            call BasisLocalQuadNedelec1(refpts, Th, Vh, i_elem, deriv_type, result)
        case (DOF_QuadRT1)
            call BasisLocalQuadRT1(refpts, Th, Vh, i_elem, deriv_type, result)
        case default
            print *, 'BasisLocal2D: Unknown basis type'
            stop
        end select
    end subroutine BasisLocal2D
    
    ! Find the location of a degree of freedom
    ! Lagrangian bases are supported
    subroutine FindDofLocation(Th, Vh, i_dof, pt, dim)
        type(mesh2D), intent(in) :: Th
        type(fespace), intent(in) :: Vh
        integer, intent(in) :: i_dof
        real(8), dimension(DIM__), intent(out) :: pt
        integer, intent(out) :: dim
        
        select case (Vh%basis_type)
        case (DOF_P1)
            call FindDofLocation_P1(Th, i_dof, pt, dim)
        case (DOF_P2)
            call FindDofLocation_P2(Th, i_dof, pt, dim)
        case default
            print *, 'FindDofLocation: Only Lagrangian bases are supported'
            stop
        end select
        
    end subroutine

end module fe