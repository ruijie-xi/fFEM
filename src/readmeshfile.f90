module readmeshfile
    implicit none
    
contains

    subroutine ReadNodeFile(infileNodes, nNode, ndim, coords)
        character(len=100),intent(in) :: infileNodes  ! input file name
        Integer, intent(out) :: nNode                 ! number of nodes
        double precision, dimension(:,:), allocatable, intent(out) :: coords ! nodal coordinates
        integer :: ii, io, nn
        integer :: ndim, n_attr, n_bdrymarkers
        
        logical :: FILEEXISTS

        ! check if the file exists
        INQUIRE(file=infileNodes, EXIST=FILEEXISTS)
        IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileNodes, "does not exist"
        call EXIT(1)
        END IF

        ! Open the file and count number of nodes first
        OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")
        READ(1,*, iostat=io) nNode, ndim, n_attr, n_bdrymarkers
        ALLOCATE(coords(ndim,nNode))

        DO ii=1,nNode
            READ(1,*, iostat=io) nn, coords(1,ii), coords(2,ii)
            END DO
        
        CLOSE(1)

        write(*,*) "Node file read. Number of nodes = ", nNode
        
        
    end subroutine ReadNodeFile

    subroutine ReadElemFile(infileElems, nElem, npElem, elemNodeConn)
        character(len=100),intent(in) :: infileElems  ! input file name
        Integer, intent(out) :: nElem                 ! number of elements
        integer, dimension(:,:), allocatable, intent(out) :: elemNodeConn ! element-node connectivity
        integer :: ii, io, nn, n1, n2, n3
        integer :: npElem, n_attr
        logical :: FILEEXISTS

        ! check if the file exists
        INQUIRE(file=infileElems, EXIST=FILEEXISTS)
        IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileElems, "does not exist"
        call EXIT(1)
        END IF

        ! Open the file and count number of elems
        nElem = 0
        OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
        read(2,*) nElem, npElem, n_attr
        if(npElem .NE. 3) then
            write(*,*) "This program is hardcoded for triangular elements"
            call EXIT(1)
        end if

        ALLOCATE(elemNodeConn(npElem,nElem))

        DO ii=1,nElem
            READ(2,*, iostat=io) nn, n1, n2, n3
            elemNodeConn(1,ii) = n1
            elemNodeConn(2,ii) = n2
            elemNodeConn(3,ii) = n3
        END DO 
        CLOSE(2)

        write(*,*) "Element file read. Number of elements = ", nElem
    end subroutine ReadElemFile

    subroutine ReadEdgeFile(infileEdges, nEdge, edgeNodeConn, bndy_markers)
        character(len=100),intent(in) :: infileEdges  ! input file name
        Integer, intent(out) :: nEdge                 ! number of edges
        integer, dimension(:,:), allocatable, intent(out) :: edgeNodeConn ! edge-node connectivity
        integer, dimension(:), allocatable, intent(out) :: bndy_markers ! boundary markers
        integer :: ii, io, nn, n1, n2
        integer :: n_bndymarkers
        logical :: FILEEXISTS

        ! check if the file exists
        INQUIRE(file=infileEdges, EXIST=FILEEXISTS)
        IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileEdges, "does not exist"
        call EXIT(1)
        END IF

        ! Open the file and count number of edges
        nEdge = 0
        OPEN(3, file=infileEdges,STATUS="OLD",ACTION="READ")
        read(3,*) nEdge, n_bndymarkers

        ALLOCATE(edgeNodeConn(2,nEdge))
        ALLOCATE(bndy_markers(nEdge))

        DO ii=1,nEdge
            READ(3,*, iostat=io) nn, n1, n2, bndy_markers(ii)
            edgeNodeConn(1,ii) = n1
            edgeNodeConn(2,ii) = n2
        END DO 
        CLOSE(3)

        write(*,*) "Edge file read. Number of edges = ", nEdge
    end subroutine ReadEdgeFile
    
end module ReadMeshFile