module readmeshfile
    implicit none
    
contains

    subroutine ReadNodeFile(infileNodes, coords)
        character(len=*),intent(in) :: infileNodes  ! input file name
        Integer :: nNode                 ! number of nodes
        double precision, dimension(:,:), allocatable, intent(out) :: coords ! nodal coordinates
        integer :: unit, ii, io, nn
        integer :: ndim, n_attr, n_bdrymarkers
        
        logical :: FILEEXISTS

        ! check if the file exists
        INQUIRE(file=infileNodes, EXIST=FILEEXISTS)
        IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileNodes, "does not exist"
        call EXIT(1)
        END IF

        ! Open the file and count number of nodes first
        OPEN(newunit=unit, file=infileNodes,STATUS="OLD",ACTION="READ")
        READ(unit,*, iostat=io) nNode, ndim, n_attr, n_bdrymarkers
        ALLOCATE(coords(ndim,nNode))

        DO ii=1,nNode
            READ(unit,*, iostat=io) nn, coords(1,ii), coords(2,ii)
            END DO
        
        CLOSE(unit)

        write(*,*) "Node file read. Number of nodes = ", nNode
        
        
    end subroutine ReadNodeFile

    subroutine ReadElemFile(infileElems, elemNodeConn)
        character(len=*),intent(in) :: infileElems  ! input file name
        Integer :: nElem   ! number of elements
        integer, dimension(:,:), allocatable, intent(out) :: elemNodeConn ! element-node connectivity
        integer :: unit, ii, io, nn, n1, n2, n3, n4
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
        OPEN(newunit=unit, file=infileElems,STATUS="OLD",ACTION="READ")
        read(unit,*) nElem, npElem, n_attr

        ALLOCATE(elemNodeConn(npElem,nElem))

        DO ii=1,nElem
            if (npElem == 3) then
                READ(unit,*, iostat=io) nn, n1, n2, n3
                elemNodeConn(1,ii) = n1
                elemNodeConn(2,ii) = n2
                elemNodeConn(3,ii) = n3
            else if (npElem == 4) then
                READ(unit,*, iostat=io) nn, n1, n2, n3, n4
                elemNodeConn(1,ii) = n1
                elemNodeConn(2,ii) = n2
                elemNodeConn(3,ii) = n3
                elemNodeConn(4,ii) = n4
            end if
        END DO 
        CLOSE(unit)

        write(*,*) "Element file read. Number of elements = ", nElem
    end subroutine ReadElemFile
    
end module ReadMeshFile