 
        subroutine ums2ra (presrv, n, nz, cperm, rperm, pr,
     $          w, nblks, arx, ari, nzoff, nzdia,
     $          icntl, mp, blkp, mi, mx, info, offp, on, nzblk,
     $          cblk, kn, nz2, nbelow)
c
cc UMS2RA is a utility which converts a column-oriented matrix to arrowhead form.
c
        integer n, nz, cperm (n), rperm (n), pr (n), kn, w (kn+1),
     $          nblks, nzblk, ari (nzblk), nzoff, nzdia, mp (n+1),
     $          mi (nz), on, icntl (20), blkp (nblks+1), nz2,
     $          info (40), offp (on+1), cblk, nbelow
        logical presrv
        real
     $          arx (nzblk), mx (nz)
 
c=== ums2ra ============================================================
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
c  convert a column-oriented matrix into an arrowhead format.
 
c=======================================================================
c  input:
c=======================================================================
c
c       n               size of entire matrix
c       mi (1..nz):     row indices of column form of entire matrix
c       mx (1..nz):     values of column form of entire matrix
c       mp (1..n+1)     column pointers for entire matrix
c       cperm (1..n):   column permutations
c       rperm (1..n):   row permutations
c
c       if nblks > 1 and presrv
c           cblk:               the block to convert
c           kn:                 the size of the block to convert
c       else
c           cblk:               0
c           kn                  n, size of input matrix
 
c=======================================================================
c  output:
c=======================================================================
c
c       if nblks = 1 and not presrv
c
c           nzoff               0
c           nzdia               nz - (entries below in diagonal blocks)
c           nz2                 nzdia
c
c           mi (1..nz2)         arrowheads for the diagonal block
c           mx (1..nz2)
c           ari, arx            used as workspace
c           w (1..n+1)          pointer to each arrowhead in mi/mx
c
c           offp                not accessed
c
c       if nblks = 1 and presrv
c
c           nzoff               0
c           nzdia               nz - (entries below in diagonal blocks)
c           nz2                 nzdia
c
c           mi, mx              not modified
c           ari (1..nz2)        arrowheads for the diagonal block
c           arx (1..nz2)
c           w (1..n+1)          pointer to each arrowhead in ari/arx
c
c           offp                not accessed
c
c       else if nblks > 1 and not presrv
c
c           nzoff               number of entries in off-diagonal part
c           nzdia               number of entries in diagonal blocks
c                               (nz = nzoff + nzdia + entries below
c                               diagonal blocks)
c           nz2                 nzoff + nzdia
c
c           mi (nzoff+1..nz2)   arrowheads for each diagonal block
c           mx (nzoff+1..nz2)
c           ari, arx            used as workspace
c           w (1..n+1)          pointer to each arrowhead in mi/mx
c
c           offp (1..n+1)       row pointers for off-diagonal part
c           mi (1..nzoff)       col indices for off-diagonal part
c           mx (1..nzoff)       values for off-diagonal part
c
c       else (nblks > 1 and presrv)
c
c           nzoff               0
c           nzdia               nonzeros in the diagonal block, cblk
c           nz2                 nzdia
c
c           mi, mx              not modified
c           ari (1..nz2)        arrowheads for the diagonal block, cblk
c           arx (1..nz2)
c           w (1..kn+1)         pointer to each arrowhead in ari/arx
c
c           offp                not accessed
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2r0
c       subroutines called:     ums2of
c       functions called:       min
        intrinsic min
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer i, p, row, col, blk, base, k1, k2, k, b1, b2, k0
 
c  i:       loop index, arrowhead index
c  p:       pointer into column-form input matrix
c  row:     row index
c  col:     column index
c  blk:     current diagonal block
c  base:    where to start the construction of the arrowhead form
c  k1,k2:   current diagonal block is a (k1..k2, k1..k2)
c  k:       loop index, kth pivot
c  b1,b2:   convert blocks b1...b2 from column-form to arrowhead form
c  k0:      convert a (k0+1..., k0+1...) to arrowhead form
 
c=======================================================================
c  executable statements:
c=======================================================================
 
c-----------------------------------------------------------------------
c  if entire matrix is to be converted, then create the off-diagonal
c  part in row-oriented form in ari (1..nzoff) and arx (1..nzoff) and
c  compute inverse row permutation.  otherwise, the inverse row
c  permutation has already been computed.
c-----------------------------------------------------------------------
 
        nzoff = 0
        nbelow = 0
        if (nblks .eq. 1) then
           do k = 1, n
              pr (rperm (k)) = k
           end do
        else if (nblks .gt. 1 .and. .not. presrv) then
           call ums2of (w, n, rperm, cperm, nzoff,
     $        offp, ari, arx, pr,
     $        icntl, mp, mi, mx, n, nz, .true., nblks, blkp,
     $        nz, 2, info, nbelow)
        endif
 
c-----------------------------------------------------------------------
c  construct the arrowhead form for the diagonal block(s)
c-----------------------------------------------------------------------
 
        do 20 i = 1, kn+1
           w (i) = 0
20      continue
 
        base = nzoff + 1
 
        if (cblk .ne. 0) then
c          convert just cblk
           k0 = blkp (cblk) - 1
           b1 = cblk
           b2 = cblk
        else
c          convert all the block(s)
           k0 = 0
           b1 = 1
           b2 = nblks
        endif
 
        do 80 blk = b1, b2
 
c          -------------------------------------------------------------
c          get the starting and ending indices of this diagonal block
c          -------------------------------------------------------------
 
           if (nblks .gt. 1) then
              k1 = blkp (blk)
              k2 = blkp (blk+1) - 1
           else
              k1 = 1
              k2 = n
           endif
 
c          -------------------------------------------------------------
c          count the number of entries in each arrowhead
c          -------------------------------------------------------------
 
           do 40 col = k1, k2
              do 30 p = mp (cperm (col)), mp (cperm (col) + 1) - 1
                 row = pr (mi (p))
                 if (row .ge. k1 .and. row .le. k2) then
c                   this is in a diagonal block, arrowhead i
                    i = min (row, col) - k0
                    w (i) = w (i) + 1
                 endif
30            continue
40         continue
 
c          -------------------------------------------------------------
c          set pointers to point just past end of each arrowhead
c          -------------------------------------------------------------
 
           w (k2-k0+1) = w (k2-k0) + base
           do 50 i = k2-k0, k1-k0+1, -1
              w (i) = w (i+1) + w (i-1)
50         continue
           w (k1-k0) = w (k1-k0+1)
c          w (i+1-k0) points just past end of arrowhead i in ari/arx
 
c          -------------------------------------------------------------
c          construct arrowhead form, leaving pointers in final state
c          -------------------------------------------------------------
 
           do 70 col = k1, k2
              do 60 p = mp (cperm (col)), mp (cperm (col) + 1) - 1
                 row = pr (mi (p))
                 if (row .ge. k1 .and. row .le. k2) then
                    if (row .ge. col) then
c                      diagonal, or lower triangular part
                       i = col - k0 + 1
                       w (i) = w (i) - 1
                       ari (w (i)) = row - k1 + 1
                       arx (w (i)) = mx (p)
                    else
c                      upper triangular part, flag by negating col
                       i = row - k0 + 1
                       w (i) = w (i) - 1
                       ari (w (i)) = -(col - k1 + 1)
                       arx (w (i)) = mx (p)
                    endif
                 endif
60            continue
70         continue
 
           base = w (k1-k0)
           w (k2-k0+1) = 0
80      continue
 
        w (kn+1) = nzoff + 1
        nzdia = base - nzoff - 1
        nz2 = nzoff + nzdia
 
c       ----------------------------------------------------------------
c       if cblk = 0, the entire matrix has been converted:
c
c          w (i) now points just past end of arrowhead i in ari/arx
c          arrowhead i is located in ari/arx (w (i+1) ... w (i)-1),
c          except for the k2-th arrowhead in each block.  those are
c          located in ari/arx (base ... w (k2) - 1), where base is
c          w (blkp (blk-1)) if blk>1 or w (n+1) = nzoff + 1 otherwise.
c
c       otherwise, just one block has been converted:
c
c          w (i) now points just past end of arrowhead i in ari/arx,
c          where i = 1 is the first arrowhead of this block (not the
c          first arrowhead of the entire matrix).  arrowhead i is
c          located in ari/arx (w (i+1) ... w (i)-1).
c          this option is used only if nblks>1 and presrv is true.
c       ----------------------------------------------------------------
 
c-----------------------------------------------------------------------
c  if not preserved, overwrite column-form with arrowhead form
c-----------------------------------------------------------------------
 
        if (.not. presrv) then
           do 90 i = 1, nz
              mi (i) = ari (i)
              mx (i) = arx (i)
90         continue
        endif
 
        return
        end
