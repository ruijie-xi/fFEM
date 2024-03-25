 
        subroutine ums2of (w, n, rperm, cperm, nzoff,
     $          offp, offi, offx, pr,
     $          icntl, mp, mi, mx, mn, mnz, presrv, nblks, blkp,
     $          onz, who, info, nbelow)
c
cc UMS2OF is a utility routine which permutes off-diagonal blocks.
c
        integer n, nzoff, w (n+1), rperm (n), cperm (n), onz,
     $          offp (n+1), offi (onz), pr (n), icntl (20), mn, mnz,
     $          mp (mn+1), mi (mnz), nblks, blkp (nblks+1), who, nbelow,
     $          info (40)
        logical presrv
        real
     $          offx (onz), mx (mnz)
 
c=== ums2of ============================================================
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
c  permute the off-diagonal blocks according to final pivot permutation.
c  this routine is called only if the block-triangular-form (btf) is
c  used.
 
c=======================================================================
c  input:
c=======================================================================
c
c       n:              order of matrix
c       rperm (1..n):   the final row permutations, including btf
c                       if i is the k-th pivot row, then rperm (k) = i
c       cperm (1..n):   the final column permutations, including btf
c                       if j is the k-th pivot col, then cperm (k) = j
c       icntl:          integer control parameters, see ums2in
c       info:           integer informational parameters
c       who:            who called (1: ums2fa, 2: ums2rf)
c
c       if presrv is true then
c           mn:                 order of preserved matrix
c           mnz:                number of entries in preserved matrix
c           mp (1..mn+1):       column pointers of preserved matrix
c           mi (1..mnz):        row indices of preserved matrix
c           mx (1..mnz):        values of preserved matrix
c           blkp (1..nblks+1):  the index range of the blocks
c           nblks:              the number of diagonal blocks
c       else
c           mn:                 0
c           mnz:                nzoff
c           mp:                 unaccessed
c           offp (1..n+1):      column pointers for off-diagonal entries
c                               in original order
c           mi (1..mnz):        the row indices of off-diagonal entries,
c                               in original order
c           mx (1..mnz):        the values of off-diagonal entries,
c                               in original order
c           nblks:              0
c           blkp (1..nblks+1):  unaccessed
c           nzoff:              number of entries in off-diagonal blocks
 
c=======================================================================
c  workspace:
c=======================================================================
c
c       w (1..n)
 
c=======================================================================
c  output:
c=======================================================================
c
c       offp (1..n+1):          row pointers for off-diagonal part
c       offi (1..nzoff):        column indices in off-diagonal part
c       offx (1..nzoff):        values in off-diagonal part
c       nzoff:                  number of entries in off-diagonal blocks
c       pr (1..n):              inverse row permutation
c       nbelow:                 entries that are below the diagonal
c                               blocks (can only occur if who = 2)
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2f0, ums2ra, ums2r0
c       subroutines called:     ums2p2
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer row, col, p, blk, k, k1, k2, io, prl
        logical pr3
 
c  row:     row index
c  col:     column index
c  p:       pointer
c  blk:     current diagonal block
c  k:       kth pivot
c  k1,k2:   current diaogonal block is a (k1..k2, k1..k2)
c  io:      i/o unit for diagnostic messages
c  prl:     printing level
c  pr3:     true if printing entries below diagonal blocks (ums2rf)
 
c=======================================================================
c  executable statments:
c=======================================================================
 
        io = icntl (2)
        prl = icntl (3)
        pr3 = prl .ge. 3 .and. io .ge. 0
 
c-----------------------------------------------------------------------
c  compute inverse row permutation
c-----------------------------------------------------------------------
 
c       if original row i is the kth pivot row, then
c               rperm (k) = i
c               pr (i) = k
c       if original col j is the kth pivot col, then
c               cperm (k) = j
cfpp$ nodepchk l
        do k = 1, n
           pr (rperm (k)) = k
        end do
 
c-----------------------------------------------------------------------
c  construct row-oriented pointers for permuted row-form
c-----------------------------------------------------------------------
 
        w (1) = 1
        do 20 row = 2, n
           w (row) = 0
20      continue
        nbelow = 0
        if (presrv) then
           do 50 blk = 1, nblks
              k1 = blkp (blk)
              k2 = blkp (blk+1) - 1
              do 40 col = k1, k2
cfpp$ nodepchk l
                 do 30 p = mp (cperm (col)), mp (cperm (col)+1)-1
                    row = pr (mi (p))
                    if (row .lt. k1) then
c                      offdiagonal entry
                       w (row) = w (row) + 1
                    else if (row .gt. k2 .and. who .eq. 2) then
c                      this entry is below the diagonal block - invalid.
c                      this can only occur if who = 2 (ums2rf).
                       if (pr3) then
c                         print the original row and column indices:
                          call ums2p2 (2, 96, mi(p), col,mx(p),io)
                       endif
                       nbelow = nbelow + 1
                    endif
30               continue
40            continue
50         continue
        else
           do 70 col = 1, n
cfpp$ nodepchk l
              do 60 p = offp (col), offp (col+1) - 1
                 row = pr (mi (p))
                 w (row) = w (row) + 1
60            continue
70         continue
        endif
        do 80 row = 2, n
           w (row) = w (row) + w (row-1)
80      continue
        w (n+1) = w (n)
c       w (row) now points just past end of row in offi/x
 
c-----------------------------------------------------------------------
c  construct the row-oriented form of the off-diagonal values,
c  in the final pivot order.  the column indices in each row
c  are placed in ascending order (the access of offi/offx later on
c  does not require this, but it makes access more efficient).
c-----------------------------------------------------------------------
 
        if (presrv) then
           do 110 blk = nblks, 1, -1
              k1 = blkp (blk)
              k2 = blkp (blk+1) - 1
              do 100 col = k2, k1, - 1
cfpp$ nodepchk l
                 do 90 p = mp (cperm (col)), mp (cperm (col)+1)-1
                    row = pr (mi (p))
                    if (row .lt. k1) then
c                      offdiagonal entry
                       w (row) = w (row) - 1
                       offi (w (row)) = col
                       offx (w (row)) = mx (p)
                    endif
90               continue
100           continue
110        continue
        else
           do 130 col = n, 1, -1
cfpp$ nodepchk l
              do 120 p = offp (cperm (col)), offp (cperm (col) + 1) - 1
                 row = pr (mi (p))
                 w (row) = w (row) - 1
                 offi (w (row)) = col
                 offx (w (row)) = mx (p)
120           continue
130        continue
        endif
 
c-----------------------------------------------------------------------
c  save the new row pointers
c-----------------------------------------------------------------------
 
        do 140 row = 1, n+1
           offp (row) = w (row)
140     continue
 
        nzoff = offp (n+1) - 1
 
        return
        end
