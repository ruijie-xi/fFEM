         subroutine ums2co (n, nz, transa, xx, xsize, info, icntl,
     $          ii, isize, w, wp, who)
c
cc UMS2CO converts matrix from triplet to column-oriented form.
c
        integer isize, ii (isize), n, nz, w (n), wp (n+1), info (40),
     $          icntl (20), xsize, who
        real
     $          xx (xsize)
        logical transa
 
c=== ums2co ============================================================
c
c  unsymmetric-pattern multifrontal package (umfpack). version 2.0s.
c  copyright (c) 1995, timothy a. davis, university of florida, usa.
c  joint work with iain s. duff, rutherford appleton laboratory, uk.
c  september 1995. work supported by the national science foundation
c  (dms-9223088 and dms-9504974) and the state of florida; and by cray
c  research inc. through the allocation of supercomputing resources.
 
c***********************************************************************
c* notice:  "the umfpack package may be used solely for educational,   *
c* research, and benchmarking purposes by non-profit organizations and *
c* the u.s. government.  commericial and other organizations may make  *
c* use of umfpack solely for benchmarking purposes only.  umfpack may  *
c* be modified by or on behalf of the user for such use but at no time *
c* shall umfpack or any such modified version of umfpack become the    *
c* property of the user.  umfpack is provided without warranty of any  *
c* kind, either expressed or implied.  neither the authors nor their   *
c* employers shall be liable for any direct or consequential loss or   *
c* damage whatsoever arising out of the use or misuse of umfpack by    *
c* the user.  umfpack must not be sold.  you may make copies of        *
c* umfpack, but this notice and the copyright notice must appear in    *
c* all copies.  any other use of umfpack requires written permission.  *
c* your use of umfpack is an implicit agreement to these conditions."  *
c*                                                                     *
c* the ma38 package in the harwell subroutine library (hsl) has        *
c* equivalent functionality (and identical calling interface) as       *
c* umfpack.  it is available for commercial use.   technical reports,  *
c* information on hsl, and matrices are available via the world wide   *
c* web at http://www.cis.rl.ac.uk/struct/arcd/num.html, or by          *
c* anonymous ftp at seamus.cc.rl.ac.uk/pub.  also contact john         *
c* harding, harwell subroutine library, b 552, aea technology,         *
c* harwell, didcot, oxon ox11 0ra, england.                            *
c* telephone (44) 1235 434573, fax (44) 1235 434340,                   *
c* email john.harding@aeat.co.uk, who will provide details of price    *
c* and conditions of use.                                              *
c***********************************************************************
 
c=======================================================================
c  not user-callable.
 
c=======================================================================
c  description:
c=======================================================================
c
c  convert input matrix (ii,xx,n,nz) into from triplet form to column-
c  oriented form.  remove invalid entries and duplicate entries.
 
c=======================================================================
c  input:
c=======================================================================
c
c       n:              size of the matrix
c       nz:             number of nonzeros in the input matrix
c       transa:         if true then transpose input matrix
c       xx (1..nz):     values of triplet form
c       xsize:          size of xx, must be >= 2*nz
c       ii (1..2*nz):   row and col indices of triplet form
c       isize:          size of ii, must be >= max (2*nz,n+1) + nz
c       icntl:          integer control parameters
c       who:            who called ums2co, 1: ums2fa, 2: ums2rf
c
c       ii must be at least of size (nz + max (2*nz, n+1))
c       xx must be at least of size (nz + max (  nz, n+1))
c
c       input triplet matrix:
c          if (transa) is false:
c               ii (p)          row index, for p = 1..nz
c               ii (nz+p)       col index
c               xx (p)          value
c          if (transa) is true:
c               ii (p)          col index, for p = 1..nz
c               ii (nz+p)       row index
c               xx (p)          value
 
c=======================================================================
c  workspace:
c=======================================================================
c
c       w (1..n)
 
c=======================================================================
c  output:
c=======================================================================
c
c       nz:             number of nonzeros in the output matrix,
c                       after removing invalid entries, and summing up
c                       duplicate entries
c       ii (n+2..nz+n+1): row indices in column-form
c       xx (1..nz):     values in column-form.
c       info (1):       error flag
c       info (3):       invalid entries
c       info (2):       duplicate entries
c       info (5):       remaining valid entries
c       info (6):       remaining valid entries
c       info (7):       0
c       wp (1..n+1)     column pointers for column form
c       ii (1..n+1)     column pointers for column form
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutines:  ums2fa, ums2rf
c       subroutines called:     ums2er, ums2p2
c       functions called:       max
        intrinsic max
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer row, col, pdest, p, nz1, pcol, ip, xp, io, prl, ninvld,
     $          ndupl, i
        logical pr3
        real
     $          zero
        parameter (zero = 0.0)
 
c  row:     row index
c  col:     column index
c  pdest:   location of an entry in the column-form, for dupl. removal
c  p:       pointer
c  nz1:     number of entries after removing invalid or duplic. entries
c  pcol:    column col starts here after duplicates removed
c  ip:      column-form copy of matrix placed in ii (ip...ip+nz-1)
c  xp:      column-form copy of matrix placed in xx (xp...xp+nz-1)
c  ninvld:  number of invalid entries
c  ndupl:   number of duplicate entries
c  i:       a row index if transa is true, a column index otherwise
c  io:      i/o unit for warning messages (for invalid or dupl. entries)
c  prl:     printing level
c  pr3:     true if printing invalid and duplicate entries
 
c=======================================================================
c  executable statements:
c=======================================================================
 
c-----------------------------------------------------------------------
c  get arguments and check memory sizes
c-----------------------------------------------------------------------
 
        io = icntl (2)
        prl = icntl (3)
        pr3 = prl .ge. 3 .and. io .ge. 0
 
c-----------------------------------------------------------------------
c  count nonzeros in columns and check for invalid entries
c-----------------------------------------------------------------------
 
        ninvld = 0
        ndupl = 0
        do col = 1, n
           w (col) = 0
        end do
        nz1 = nz
        do 20 p = nz, 1, -1
           row = ii (p)
           col = ii (nz+p)
           if (row.lt.1.or.row.gt.n.or.col.lt.1.or.col.gt.n) then
c             this entry is invalid - delete it
              if (pr3) then
c                print the offending entry on the diagnostic i/o unit
                 call ums2p2 (who, 99, row, col, xx(p), io)
              endif
              ii (p)    = ii (nz1)
              ii (nz+p) = ii (nz+nz1)
              xx (p)    = xx (nz1)
              nz1 = nz1 - 1
           else
              if (transa) then
c                factorizing a transpose
                 w (row) = w (row) + 1
              else
c                factorizing a
                 w (col) = w (col) + 1
              endif
           endif
20      continue
        ninvld = nz - nz1
        if (ninvld .ne. 0) then
c          invalid entries found - set warning flag and continue
           call ums2er (who, icntl, info, 1, ninvld)
        endif
 
c-----------------------------------------------------------------------
c  convert triplet form to column-form
c-----------------------------------------------------------------------
 
        wp (1) = 1
        do i = 1, n
           wp (i+1) = wp (i) + w (i)
        end do
        do i = 1, n
           w (i) = wp (i)
        end do
 
c       ----------------------------------------------------------------
c       construct column-form in ii (2*nz+1..3*nz) and xx (nz+1..2*nz)
c       ----------------------------------------------------------------
 
        ip = max (2*nz, n+1)
        xp = nz
        if (transa) then
           do p = 1, nz1
              row = ii (p)
              col = ii (nz+p)
              ii (ip + w (row)) = col
              xx (xp + w (row)) = xx (p)
              w (row) = w (row) + 1
           end do
        else
           do 60 p = 1, nz1
              row = ii (p)
              col = ii (nz+p)
              ii (ip + w (col)) = row
              xx (xp + w (col)) = xx (p)
              w (col) = w (col) + 1
60         continue
        endif
 
c       ----------------------------------------------------------------
c       shift the matrix back to ii (n+2..nz+n+1) and xx (n+2..nz+n+1)
c       ----------------------------------------------------------------
 
        nz = nz1
cfpp$ nodepchk l
        do 70 p = 1, nz
           ii (n+1+p) = ii (ip+p)
           xx (p) = xx (xp+p)
70      continue
 
c-----------------------------------------------------------------------
c  remove duplicate entries by adding them up
c-----------------------------------------------------------------------
 
        do 80 row = 1, n
           w (row) = 0
80      continue
        pdest = 1
        do 100 col = 1, n
           pcol = pdest
           do 90 p = wp (col), wp (col+1)-1
              row = ii (n+1+p)
              if (w (row) .ge. pcol) then
c                this is a duplicate entry
                 xx (w (row)) = xx (w (row)) + xx (p)
                 if (pr3) then
c                   print the duplicate entry on the diagnostic i/o
c                   unit.  the row and column indices printed reflect
c                   the input matrix.
                    if (transa) then
                       call ums2p2 (who, 98, col, row, xx (p), io)
                    else
                       call ums2p2 (who, 98, row, col, xx (p), io)
                    endif
                 endif
              else
c                this is a new entry, store and record where it is
                 w (row) = pdest
                 if (pdest .ne. p) then
                    ii (n+1+pdest) = row
                    xx (pdest) = xx (p)
                 endif
                 pdest = pdest + 1
              endif
90         continue
           wp (col) = pcol
100     continue
        wp (n+1) = pdest
        nz1 = pdest - 1
        ndupl = nz - nz1
        if (ndupl .ne. 0) then
c          duplicate entries found - set warning flag and continue
           call ums2er (who, icntl, info, 2, ndupl)
        endif
        nz = nz1
 
c-----------------------------------------------------------------------
c  save column pointers in ii (1..n+1)
c-----------------------------------------------------------------------
 
        do col = 1, n+1
           ii (col) = wp (col)
        end do
 
        info (2) = ndupl
        info (3) = ninvld
        info (5) = nz
        info (6) = nz
        info (7) = 0
        if (nz .eq. 0) then
c          set error flag if all entries are invalid
           call ums2er (who, icntl, info, -2, -1)
        endif
        return
        end
