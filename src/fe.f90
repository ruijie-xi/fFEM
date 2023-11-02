module fe
    use mesh
    use settings
    implicit none
    
contains

    subroutine fespace_init(Vh, Th, basis_type, dim)
        type(fespace), intent(out) :: Vh
        type(mesh2D), intent(in) :: Th
        integer, intent(in) :: basis_type, dim

        integer :: i_dim

        select case(basis_type)
            case(DOF_P1)
                Vh%N_local_basis = 3
                Vh%dim = dim
                Vh%basis_type = basis_type
                Vh%N_DOF = Th%N_node * Vh%dim
                allocate(Vh%ElemDOF(Vh%N_local_basis*Vh%dim, Th%N_elem))
                Vh%ElemDOF = 0
                do i_dim = 1, Vh%dim
                    Vh%ElemDOF(3*(i_dim-1)+1:3*i_dim, :) = Th%ElemNodeConn + (i_dim-1)*Th%N_node
                end do
            
            case(DOF_Q1)
                Vh%N_local_basis = 4
                Vh%dim = dim
                Vh%basis_type = basis_type
                Vh%N_DOF = Th%N_node * Vh%dim
                allocate(Vh%ElemDOF(Vh%N_local_basis*Vh%dim, Th%N_elem))
                Vh%ElemDOF = 0
                do i_dim = 1, Vh%dim
                    Vh%ElemDOF(4*(i_dim-1)+1:4*i_dim, :) = Th%ElemNodeConn + (i_dim-1)*Th%N_node
                end do
        end select

    end subroutine fespace_init

    subroutine fespace_free(Vh)
        type(fespace), intent(inout) :: Vh

        deallocate(Vh%ElemDOF)

    end subroutine fespace_free

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

        integer :: i_node, i_dim
        real(8), dimension(:), allocatable :: val

        select case(Vh%basis_type)
            case(DOF_P1, DOF_Q1)
                i_dim = i_dof / Th%N_node + 1
                i_node = mod(i_dof, Th%N_node)
                if (i_node == 0) then
                    i_node = Th%N_node
                    i_dim = i_dim - 1
                end if
                call fun(Th%NodeCoord(:,i_node), val, DERIV_NONE)
                result = val(i_dim)
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

        integer :: i_dim
        integer :: i_node1, i_node2

        select case (Vh%basis_type)
        case (DOF_P1 , DOF_Q1)
            i_node1 = Th%EdgeNodeConn(1, i_edge)
            i_node2 = Th%EdgeNodeConn(2, i_edge)
            allocate(dof_index(2*Vh%dim))
            do i_dim = 1, Vh%dim
                dof_index(2*(i_dim-1)+1) = (i_dim-1)*Th%N_node + i_node1
                dof_index(2*(i_dim-1)+2) = (i_dim-1)*Th%N_node + i_node2
            end do
        case default
            print *, "Error: Unknown basis type"
            stop
        end select
    end subroutine getEdgeDofIndex
    
end module fe