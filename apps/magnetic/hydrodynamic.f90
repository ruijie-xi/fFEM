module Hydrodynamic
    use settings
    use mesh
    use fe
    implicit none

    private
    public :: HydrodynamicInitialize, getVelocity, MoveMesh

    type(mesh2D) :: Th
    type(FESPACE) :: Uh

contains

subroutine HydrodynamicInitialize(elems,nodes)
    implicit none
    integer, dimension(:,:), intent(inout) :: elems
    real(8), dimension(:,:), intent(in) :: nodes

    call MeshInit(elems,nodes,Th)

    call fespaceInit(Uh, Th, DOF_Q1, 2)

end subroutine HydrodynamicInitialize

subroutine getVelocity(u_func, u)
    implicit none
    type(vector) :: u
    procedure(func) :: u_func
    
    call Interpolate(u, u_func, Th, Uh)

end subroutine getVelocity

subroutine MoveMesh(u, dt, nodes)
    implicit none
    type(vector) :: u
    real(8), dimension(:,:), intent(out), allocatable :: nodes
    real(8), intent(in) :: dt

    Th%NodeCoord(1,:) = Th%NodeCoord(1,:) + u%data(1:Th%N_Node)*dt
    Th%NodeCoord(2,:) = Th%NodeCoord(2,:) + u%data(Th%N_Node+1:2*Th%N_Node)*dt

    if (allocated(nodes)) deallocate(nodes)
    allocate(nodes(2,Th%N_Node))
    nodes = Th%NodeCoord

end subroutine MoveMesh
    
end module Hydrodynamic