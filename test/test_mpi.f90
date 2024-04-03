program test_mpi
    implicit none
    include 'mpif.h'
    
    
    integer :: ierr, myrank, nprocs

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    print *, 'Hello from process', myrank, 'of', nprocs

    call MPI_FINALIZE(ierr)

    
end program test_mpi