 
        subroutine ums2fb (xx, xsize, ii, isize, n, nz, nzdia, nzoff,
     $          nblks, cp, cperm, rperm, pr, pc,
     $          w, zperm, bp, offp,
     $          presrv, icntl)
c
cc UMS2FB is a utility function which finds permutations to block triangular form.
c
        integer n, nz, isize, ii (isize), nzdia, nzoff, nblks, cp (n+1),
     $          cperm (n), rperm (n), pr (n), pc (n), w (n), zperm (n),
     $          bp (n+1), offp (n+1), icntl (20), xsize
        logical presrv
        real
     $          xx (xsize)
 
c=== ums2fb ============================================================
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
c  find permutations to block triangular form:
c       1) permute the matrix so that it has a zero-free diagonal.
c       2) find strongly-connected components of the corresponding
c          graph.  each diagonal block corresponds to exactly one
c          strongly-connected component.
c       3) convert the matrix to block triangular form, unless it is
c          to be preserved in its original form.
 
c  calls harwell ma28 routines mc21b and mc13e, which can be obtained
c  separately from netlib.  send email to netlib@ornl.gov with the
c  message:
c       send mc13e.f mc21b.f from harwell
 
c=======================================================================
c  installation note:
c=======================================================================
c
c  if the ma28 harwell subroutine library routines mc21b and mc13e
c  (which perform the permutation to block-triangular-form) are not
c  available, then you may comment out all executable code in this
c  routine, or place a "return" statement as the first executable
c  statement (see below).  if you do make this modification, please do
c  not delete any original code.  add a comment and date to your
c  modifications.
 
c  sept. 14, 1995:  changed non-ansi "60 enddo" statement to
c       "60 continue".
 
c=======================================================================
c  input:
c=======================================================================
c
c       presrv:         true if original matrix is to be preserved
c       n:              order of matrix
c       nz:             entries in matrix
c       isize:          size of ii
c       xsize:          size of xx
c       cp (1..n+1):    column pointers
c       xx (1..nz):     values
c       ii (1..nz):     row indices
c       icntl:          integer control arguments
c
c          input matrix in column form is in:
c          xx (1..nz), ii (1..nz), n, nz, cp (1..n+1), where
c               ii (cp(col) ... cp(col+1)-1): row indices
c               xx (cp(col) ... cp(col+1)-1): values
c          if presrv is false then xsize and isize must be >= 2*nz
c          otherwise, xsize and isize must be >= nz
 
c=======================================================================
c  workspace:
c=======================================================================
c
c       pr (1..n), pc (1..n), w (1..n), zperm (1..n)
 
c======================================================================
c  output:
c=======================================================================
c
c       nblks: number of blocks
c       if (nblks > 1):
c
c           cperm (1..n), rperm (1..n): permutation to block form:
c               rperm (newrow) = oldrow
c               cperm (newcol) = oldcol
c
c           bp (n-nblks+1...n+1) holds the start/end of blocks 1..nblks
c
c           if (presrv is false) then
c
c              input matrix is converted to block-upper-tri. form,
c              using ii/xx (nz+1..2*nz) as workspace.
c              nzdia: nonzeros in diagonal blocks
c              nzoff: nonzeros in off-diagonal blocks
c              (nz = nzdia + nzoff)
c
c              off-diagonal column-oriented form in xx/ii (1..nzoff)
c              col is located in
c              xx/ii (offp (col) ... offp (col+1)-1)
c
c              diagonal blocks now in xx/ii (nzoff+1 .. nzoff+nzdia)
c              col is located in
c              xx/ii (cp (col) ... cp (col+1)-1)
c
c       else, nblks=1: and no other output is generated.
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2f0
c       subroutines called:     mc21b, mc13e (in ma28 hsl package)
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer col, ndiag, i, po, pb, blk, p, row, k1, k
 
c  ndiag:   number of entries on the diagonal
c  po:      pointer into off-diagonal part
c  pb:      pointer into diagonal blocks
c  blk:     block number for current diagonal block
c  k1:      column col is in diagonal block a (k1.., k1...)
c  k:       kth row/col in btf form is rperm(k)/cperm(k) in input matrix
c  p:       pointer
c  row:     row index
c  col:     column index
c  i:       general loop index
 
c=======================================================================
c  executable statements:
c       if (mc21b and mc13e not available during installation) return
c=======================================================================
 
        nzdia = nz
        nzoff = 0
 
c-----------------------------------------------------------------------
c compute the length of each column
c-----------------------------------------------------------------------
 
        do col = 1, n
           w (col) = cp (col+1) - cp (col)
        end do
 
c-----------------------------------------------------------------------
c find a column permutation for a zero-free diagonal
c-----------------------------------------------------------------------
 
        call mc21b (n, ii, nz, cp, w, zperm, ndiag, offp, cperm, pr, pc)
c          mc21b calling interface:
c          input:       n, ii (1..nz), nz, cp (n), w (n):
c                       n-by-n matrix, col is of length w (col),
c                       and its pattern is located in
c                       ii (cp (col) ... cp (col)+w(col)-1)
c          output:      zperm (n), the permutation, such that
c                       colold = zperm (col), and ndiag (number of
c                       structural nonzeros on the diagonal.
c                       matrix is structurally singular if ndiag < n
c          workspace:   offp, cperm, pr, pc
 
c-----------------------------------------------------------------------
c  permute the columns of the temporary matrix to get zero-free diagonal
c-----------------------------------------------------------------------
 
        do 20 col = 1, n
           offp (col) = cp (zperm (col))
           w (col) = cp (zperm (col)+1) - cp (zperm (col))
20      continue
 
c-----------------------------------------------------------------------
c  find a symmetric permutation into upper block triangular form
c  (that is, find the strongly-connected components in the graph).
c-----------------------------------------------------------------------
 
        call mc13e (n, ii, nz, offp, w, rperm, bp, nblks, cperm, pr, pc)
c          mc13e calling interface:
c          input:       n, ii (1..nz), nz, offp (n), w (n)
c                       n-by-n matrix, col of length w(col),
c                       in ii (offp(col) ... offp(col)+w(col)-1), where
c                       this permuted matrix has a zero-free diagonal
c                       (unless the matrix is structurally singular).
c          output:      rperm (n), bp (n+1), nblks
c                       old = rperm (new) is the symmetric permutation,
c                       there are nblks diagonal blocks, bp (i) is
c                       the position in new order of the ith block.
c          workspace:   cperm, pr, pc
 
c-----------------------------------------------------------------------
c  if more than one block, get permutations and block pointers,
c  and convert to block-upper-triangular form (unless matrix preserved)
c-----------------------------------------------------------------------
 
        if (nblks .ne. 1) then
 
c          -------------------------------------------------------------
c          find the composite column permutation vector (cperm):
c          -------------------------------------------------------------
 
           do 30 col = 1, n
              cperm (col) = zperm (rperm (col))
30         continue
 
c          -------------------------------------------------------------
c          convert to block-upper-triangular form, if not preserved
c          -------------------------------------------------------------
 
           if (.not. presrv) then
 
c             ----------------------------------------------------------
c             find the inverse permutation vectors, pr and pc
c             ----------------------------------------------------------
 
              do 40 k = 1, n
                 pc (cperm (k)) = k
                 pr (rperm (k)) = k
40            continue
 
c             ----------------------------------------------------------
c             construct flag array to determine if entry in block or not
c             ----------------------------------------------------------
 
              bp (nblks+1) = n+1
              do 60 blk = 1, nblks
                 do 50 i = bp (blk), bp (blk+1)-1
                    w (i) = bp (blk)
50               continue
60            continue
 
c             ----------------------------------------------------------
c             construct block-diagonal form in xx/ii (nz+1..nz+nzdia)
c             ----------------------------------------------------------
 
c             these blocks are in a permuted order (according to rperm
c             and cperm).  the row indices in each block range from 1
c             to the size of the block.
 
              pb = nz + 1
              do 80 col = 1, n
                 zperm (col) = pb
                 k1 = w (col)
cfpp$ nodepchk l
                 do 70 p = cp (cperm (col)), cp (cperm (col)+1)-1
                    row = pr (ii (p))
                    if (w (row) .eq. k1) then
c                      entry is in the diagonal block:
                       ii (pb) = row - k1 + 1
                       xx (pb) = xx (p)
                       pb = pb + 1
                    endif
70               continue
80            continue
c             zperm (n+1) == pb  ( but zperm (n+1) does not exist )
              nzdia = pb - (nz + 1)
              nzoff = nz - nzdia
 
c             diagonal blocks now in xx/ii (nz+1..nz+nzdia)
c             col is located in xx/ii (zperm (col) ... zperm (col+1)-1)
 
c             ----------------------------------------------------------
c             compress original matrix to off-diagonal part, in place
c             ----------------------------------------------------------
 
c             the rows/cols of off-diagonal form correspond to rows/cols
c             in the original, unpermuted matrix.  they are permuted to
c             the final pivot order and stored in a row-oriented form,
c             after the factorization is complete (by ums2of).
 
              po = 1
              do 100 col = 1, n
                 offp (col) = po
                 k1 = w (pc (col))
cfpp$ nodepchk l
                 do 90 p = cp (col), cp (col+1)-1
                    row = pr (ii (p))
                    if (w (row) .ne. k1) then
c                      offdiagonal entry
                       ii (po) = ii (p)
                       xx (po) = xx (p)
                       po = po + 1
                    endif
90               continue
100           continue
              offp (n+1) = po
 
c             off-diagonal form now in xx/ii (1..nzoff)
c             col is located in xx/ii(offp(col)..offp(col+1)-1)
 
c             ----------------------------------------------------------
c             move block-diagonal part into place
c             ----------------------------------------------------------
 
              pb = nz + 1
cfpp$ nodepchk l
              do 110 i = 0, nzdia - 1
                 ii (po+i) = ii (pb+i)
                 xx (po+i) = xx (pb+i)
110           continue
              do 120 col = 1, n
                 cp (col) = zperm (col) - nzdia
120           continue
c             cp (n+1) == nz+1  ( this is unchanged )
 
c             diagonal blocks now in xx/ii (nzoff+1 .. nzoff+nzdia)
c             col is located in xx/ii (cp (col) ... cp (col+1)-1)
 
           endif
 
c          -------------------------------------------------------------
c          shift bp (1 .. nblks+1) down to bp (1+n-nblks .. n+1), which
c          then becomes the blkp (1 .. nblks+1) array.
c          -------------------------------------------------------------
 
           bp (nblks+1) = n+1
cfpp$ nodepchk l
           do 130 blk = nblks + 1, 1, -1
              bp (blk + (n-nblks)) = bp (blk)
130        continue
        endif
 
        return
        end
