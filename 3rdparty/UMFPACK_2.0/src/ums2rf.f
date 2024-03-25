 
        subroutine ums2rf (n, ne, job, transa, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo)
c
cc UMS2RF computes the LU factorization of a sparse matrix.
c
        integer n, ne, job, lvalue, lindex, index (lindex), keep (20),
     $          icntl (20), info (40)
        real
     $          value (lvalue), cntl (10), rinfo (20)
        logical transa
 
c=== ums2rf ============================================================
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
c  user-callable.
 
c=======================================================================
c  description:
c=======================================================================
c
c  given a sparse matrix a, and a sparsity-preserving and numerically-
c  acceptable pivot order and symbolic factorization, compute the lu
c  factors, paq = lu.  uses the sparsity pattern and permutations from
c  a prior factorization by ums2fa or ums2rf.  the matrix a should have
c  the same nonzero pattern as the matrix factorized by ums2fa or
c  ums2rf.  the matrix can have different numerical values.  no
c  variations are made in the pivot order computed by ums2fa.  if a
c  zero pivot is encountered, an error flag is set and the
c  factorization terminates.
c
c  this routine can actually handle any matrix a such that (paq)_ij can
c  be nonzero only if (lu)_ij is be nonzero, where l and u are the lu
c  factors of the matrix factorized by ums2fa.  if btf (block triangular
c  form) is used, entries above the diagonal blocks of (paq)_ij can have
c  an arbitrary sparsity pattern.  entries for which (lu)_ij is not
c  present, or those below the diagonal blocks are invalid and ignored
c  (a warning flag is set and the factorization proceeds without the
c  invalid entries).  a listing of the invalid entries can be printed.
c
c  this routine must be preceded by a call to ums2fa or ums2rf.
c  a call to ums2rf can be followed by any number of calls to ums2so,
c  which solves a linear system using the lu factors computed by this
c  routine or by ums2fa.  a call to ums2rf can also be followed by any
c  number of calls to ums2rf.
 
c=======================================================================
c  arguments:
c=======================================================================
 
c           ------------------------------------------------------------
c  n:       an integer variable.
c           must be set by caller on input (not modified).
c           order of the matrix.  must be identical to the value of n
c           in the last call to ums2fa.
 
c           ------------------------------------------------------------
c  ne:      an integer variable.
c           must be set by caller on input (not modified).
c           number of entries in input matrix.  normally not modified
c           since the last call to ums2fa.
c           restriction:  1 <= ne < (keep (4)) / 2
 
c           ------------------------------------------------------------
c  job:     an integer variable.
c           must be set by caller on input (not modified).
c           if job=1, then a column-oriented form of the input matrix
c           is preserved, otherwise, the input matrix is overwritten
c           with its lu factors.  if iterative refinement is to done
c           (icntl (8) > 0), then job must be set to 1.  can be
c           the same, or different, as the last call to ums2fa.
 
c           ------------------------------------------------------------
c  transa:  a logical variable.
c           must be set by caller on input (not modified).
c           if false then a is factorized: paq = lu.  otherwise, a
c           transpose is factorized:  pa'q = lu.  normally the same as
c           the last call to ums2fa.
 
c           ------------------------------------------------------------
c  lvalue:  an integer variable.
c           must be set by caller on input (not modified).
c           size of the value array.  restriction:  lvalue >= 2*ne,
c           although a larger will typically be required to complete
c           the factorization.  the exact value required is computed
c           by the last call to ums2fa or ums2rf (info (23)).
c           this value assumes that the ne, job, and transa parameters
c           are the same as the last call.  some garbage collection may
c           occur if lvalue is set to info (23), but usually not
c           much.  we recommend lvalue => 1.2 * info (23).  the
c           lvalue parameter is usually the same as in the last call to
c           ums2fa, however.
 
c           ------------------------------------------------------------
c  lindex:  an integer variable.
c           must be set by caller on input (not modified).
c           size of the index array.  restriction:
c           lindex >= 3*ne+2*n+1 + (keep (5) - keep (4) + 1),
c           although a larger will typically be required to complete
c           the factorization.  the exact value required is computed
c           by the last call to ums2fa or ums2rf (info (22)).
c           this value assumes that the ne, job, and transa parameters
c           are the same as the last call.  no garbage collection ever
c           occurs in the index array, since ums2rf does not create
c           external fragmentation in index.  the lindex parameter is
c           usually the same as in the last call to ums2fa, however.
c           note that lindex >= keep (5) is also required, since
c           the pattern of the prior lu factors reside in
c           index (keep (4) ... keep (5)).
 
c           ------------------------------------------------------------
c  value:   a real array of size lvalue.
c           must be set by caller on input (normally from the last call
c           to ums2fa or ums2rf).  modified on output.  on input,
c           value (1..ne) holds the original matrix in triplet form.
c           on output, value holds the lu factors, and (optionally) a
c           column-oriented form of the original matrix - otherwise
c           the input matrix is overwritten with the lu factors.
 
c           ------------------------------------------------------------
c  index:   an integer array of size lindex.
c           must be set by caller on input (normally from the last call
c           to ums2fa or ums2rf).  modified on output.  on input,
c           index (1..2*ne) holds the original matrix in triplet form,
c           and index (keep (4) ... keep (5)) holds the pattern
c           of the prior lu factors.  on output, index holds the lu
c           factors, and (optionally) a column-oriented form of the
c           original matrix - otherwise the input matrix is overwritten
c           with the lu factors.
c
c           on input the kth triplet (for k = 1...ne) is stored as:
c                       a (row,col) = value (k)
c                       row         = index (k)
c                       col         = index (k+ne)
c           if there is more than one entry for a particular position,
c           the values are accumulated, and the number of such duplicate
c           entries is returned in info (2), and a warning flag is
c           set.  however, applications such as finite element methods
c           naturally generate duplicate entries which are then
c           assembled (added) together.  if this is the case, then
c           ignore the warning message.
c
c           on input, and the pattern of the prior lu factors is in
c               index (keep (4) ... keep (5))
c
c           on output, the lu factors and the column-oriented form
c           of a (if preserved) are stored in:
c               value (keep (1)...keep (2))
c               index (keep (3)...keep (5))
c           where keep (2) = lvalue, and keep (5) = lindex.
 
c           ------------------------------------------------------------
c  keep:    an integer array of size 20.
c
c           keep (1 ... 3):  need not be set by caller on input.
c               modified on output.
c               keep (1): new lu factors start here in value
c               keep (2) = lvalue: new lu factors end here in value
c               keep (3): new lu factors start here in index
c
c           keep (4 ... 5): must be set by caller on input (normally
c               from the last call to ums2fa or ums2rf). modified on
c               output.
c               keep (4):  on input, the prior lu factors start here
c               in index, not including the prior (optionally) preserved
c               input matrix, nor the off-diagonal pattern (if btf was
c               used in the last call to ums2fa).  on output, the new
c               lu factors needed for ums2rf start here in index.
c               keep (5):  on input, the prior lu factors end here in
c               index.  on output, keep (5) is set to lindex, which
c               is where the new lu factors end in index
c
c           keep (6 ... 8):  unused.  these are used by ums2fa only.
c               future releases may make use of them, however.
c
c           keep (9 ... 20): unused.  reserved for future releases.
 
c           ------------------------------------------------------------
c  cntl:    a real array of size 10.
c           must be set by caller on input (not modified).
c           real control arguments, see ums2in for a
c           description, which sets the default values.  the current
c           version of ums2rf does not actually use cntl.  it is
c           included to make the argument list of ums2rf the same as
c           ums2fa.  ums2rf may use cntl in future releases.
 
c           ------------------------------------------------------------
c  icntl:   an integer array of size 20.
c           must be set by caller on input (not modified).
c           integer control arguments, see ums2in for a description,
c           which sets the default values.  ums2rf uses icntl (1),
c           icntl (2), icntl (3), and icntl (7).
 
c           ------------------------------------------------------------
c  info:    an integer array of size 40.
c           need not be set by caller on input.  modified on output.
c           it contains information about the execution of ums2rf.
c
c           info (1): zero if no error occurred, negative if
c               an error occurred and the factorization was not
c               completed, positive if a warning occurred (the
c               factorization was completed).
c
c               these errors cause the factorization to terminate:
c
c               error   description
c               -1      n < 1 or n > maximum value
c               -2      ne < 1 or ne > maximum value
c               -3      lindex too small
c               -4      lvalue too small
c               -5      both lindex and lvalue are too small
c               -6      prior pivot ordering no longer acceptable
c               -7      lu factors are uncomputed, or are corrupted
c
c               with these warnings the factorization was able to
c               complete:
c
c               error   description
c               1       invalid entries
c               2       duplicate entries
c               3       invalid and duplicate entries
c               4       singular matrix
c               5       invalid entries, singular matrix
c               6       duplicate entries, singular matrix
c               7       invalid and duplicate entries, singular matrix
c
c               subsequent calls to ums2rf and ums2so can only be made
c               if info (1) is zero or positive.  if info (1)
c               is negative, then some or all of the remaining
c               info and rinfo arrays may not be valid.
c
c           info (2): duplicate entries in a.  a warning is set
c               if info (2) > 0.  however, the duplicate entries
c               are summed and the factorization continues.  duplicate
c               entries are sometimes intentional - for finite element
c               codes, for example.
c
c           info (3): invalid entries in a, indices not in 1..n.
c               these entries are ignored and a warning is set in
c               info (1).
c
c           info (4): invalid entries in a, not in prior lu
c               factors.  these entries are ignored and a warning is
c               set in info (1).
c
c           info (5): entries in a after adding duplicates and
c               removing invalid entries.
c
c           info (6): entries in diagonal blocks of a.
c
c           info (7): entries in off-diagonal blocks of a.  zero
c               if info (9) = 1.
c
c           info (8): 1-by-1 diagonal blocks.
c
c           info (9): blocks in block-triangular form.
c
c           info (10): entries below diagonal in l.
c
c           info (11): entries below diagonal in u.
c
c           info (12): entries in l+u+offdiagonal part.
c
c           info (13): frontal matrices.
c
c           info (14): zero.  used by ums2fa only.
c
c           info (15): garbage collections performed on value.
c
c           info (16): diagonal pivots chosen.
c
c           info (17): numerically acceptable pivots found in a.
c               if less than n, then a is singular (or nearly so).
c               the factorization still proceeds, and ums2so can still
c               be called.  the zero-rank active submatrix of order
c               n - info (17) is replaced with the identity matrix
c               (assuming btf is not in use).  if btf is in use, then
c               one or more of the diagonal blocks are singular.
c               ums2rf can be called if the value of info (17)
c               returned by ums2fa was less than n, but the order
c               (n - info (17)) active submatrix is still replaced
c               with the identity matrix.  entries residing in this
c               submatrix are ignored, their number is included in
c               info (4), and a warning is set in info (1).
c
c           info (18): memory used in index.
c
c           info (19): memory needed in index (same as info (18)).
c
c           info (20): memory used in value.
c
c           info (21): minimum memory needed in value
c               (or minimum recommended).  if lvalue is set to
c               info (21) on a subsequent call, then a moderate
c               number of garbage collections (info (15)) will
c               occur.
c
c           info (22): memory needed in index for the next call to
c               ums2rf.
c
c           info (23): memory needed in value for the next call to
c               ums2rf.
c
c           info (24): zero.  used by ums2so only.
c
c           info (25 ... 40): reserved for future releases
 
c           ------------------------------------------------------------
c  rinfo:   a real array of size 20.
c           need not be set by caller on input.  modified on output.
c           it contains information about the execution of ums2rf.
c
c           rinfo (1): total flop count in the blas
c
c           rinfo (2): total assembly flop count
c
c           rinfo (3): zero.  used by ums2fa only.
c
c           rinfo (4): level-1 blas flops
c
c           rinfo (5): level-2 blas flops
c
c           rinfo (6): level-3 blas flops
c
c           rinfo (7): zero.  used by ums2so only.
c
c           rinfo (8): zero.  used by ums2so only.
c
c           rinfo (9 ... 20): reserved for future releases
 
c=======================================================================
c  to be preserved between calls to ums2fa, ums2rf, ums2so:
c=======================================================================
c
c  when calling ums2so to solve a linear system using the factors
c  computed by ums2fa or ums2rf, the following must be preserved:
c
c       n
c       value (keep (1)...keep (2))
c       index (keep (3)...keep (5))
c       keep (1 ... 20)
c
c  when calling ums2rf to factorize a subsequent matrix with a pattern
c  similar to that factorized by ums2fa, the following must be
c  preserved:
c
c       n
c       index (keep (4)...keep (5))
c       keep (4 ... 20)
c
c  note that the user may move the lu factors to a different position
c  in value and/or index, as long as keep (1 ... 5) are modified
c  correspondingly.
 
c## end of user documentation ##########################################
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   user routine
c       subroutines called:     ums2er, ums2p1, ums2co, ums2r0
c       functions called:       max, min
        intrinsic max, min
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer i, nz, lux1, lui1, iuse, xuse, n1, nz1, nblks,
     $          lind2, luir1, lusiz, lui2, rpermp, cpermp,
     $          offpp, lublpp, blkpp, on, nzoff, ip2, io, prl
        logical presrv, badlu
        real
     $          zero
        parameter (zero = 0.0)
 
c  printing control:
c  -----------------
c  io:      i/o unit for diagnostic messages
c  prl:     printing level
c
c  matrix to factorize:
c  --------------------
c  nz:      number of entries, after removing invalid/duplicate entries
c  presrv:  true if original matrix to be preserved
c
c  memory usage:
c  -------------
c  iuse:    current memory usage in index
c  xuse:    current memory usage in value
c  lind2:   allocatable part of index is (1..lind2)
c
c  location and status of lu factors:
c  ----------------------------------
c  lui1:    integer part of lu factors start in index (lui1...)
c  luir1:   index (luir1 ... lui2) is needed for this call to ums2rf
c  lusiz:   size of index (luir1..lui2), needed from prior lu factors
c  lui2:    integer part of lu factors end in index (..lui2)
c  lux1:    real part of lu factors in value (lux1...lvalue)
c  ip2:     pointer into trailing part of lu factors in index
c  badlu:   if true, then lu factors are corrupted or not computed
c
c  arrays and scalars allocated in lu factors (in order):
c  ------------------------------------------------------
c  ...      lu factors of each diagonal block located here
c  lublpp:  lublkp (1..nblks) array in index (lublpp..lublpp+nblks-1)
c  blkpp:   blkp (1..nblks+1) array loc. in index (blkpp...blkpp+nblks)
c  offpp:   offp (1..n+1) array located in index (offpp...offpp+n)
c  on:      size of offp array
c  cpermp:  cperm (1..n) array located in index (cpermp...cpermp+n-1)
c  rpermp:  rperm (1..n) array located in index (rpermp...rpermp+n-1)
c  nblks:   number of diagonal blocks
c  nz1:     number of entries when prior matrix factorize
c  n1:      n argument in ums2fa or ums2rf when prior matrix factorized
c
c  other:
c  ------
c  i:       loop index
 
c=======================================================================
c  executable statements:
c=======================================================================
 
        io = icntl (2)
        prl = icntl (3)
 
c-----------------------------------------------------------------------
c  clear informational output, and the unneeded part of the keep array
c-----------------------------------------------------------------------
 
        do i = 1, 40
           info (i) = 0
        end do

        do 20 i = 1, 20
           rinfo (i) = zero
20      continue
        keep (1) = 0
        keep (2) = 0
        keep (3) = 0
 
c-----------------------------------------------------------------------
c  print input arguments if requested
c-----------------------------------------------------------------------
 
        call ums2p1 (2, 1,
     $          n, ne, job, transa, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo,
     $          zero, zero, 1, zero, 1)
 
c-----------------------------------------------------------------------
c  check input arguments
c-----------------------------------------------------------------------
 
        iuse = 0
        xuse = 0
        info (5) = ne
        info (6) = ne
        if (n .lt. 1) then
c          n is too small
           call ums2er (2, icntl, info, -1, -1)
           go to 9000
        endif
        if (ne .lt. 1) then
c          ne is too small
           call ums2er (2, icntl, info, -2, -1)
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  get pointers to integer part of prior lu factors
c-----------------------------------------------------------------------
 
        luir1 = keep (4)
        lui2 = keep (5)
        lusiz = lui2 - luir1 + 1
 
        badlu = luir1 .le. 0 .or. lui2-6 .lt. luir1 .or. lui2.gt.lindex
        if (badlu) then
           call ums2er (2, icntl, info, -7, 0)
c          error return, lu factors are corrupted:
           go to 9000
        endif
        if (2*ne .gt. luir1) then
           call ums2er (2, icntl, info, -2, luir1/2)
c          error return, ne is too large:
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  shift the prior lu factors down to the end of index.  if keep and
c  lindex are unmodified from the prior call to ums2fa, then
c  keep (5) = lindex, and this shift is not performed.
c-----------------------------------------------------------------------
 
        if (lui2 .lt. lindex) then
           do 30 i = lindex, lindex - lusiz + 1, -1
              index (i) = index (i - lindex + lui2)
30         continue
           luir1 = lindex - lusiz + 1
           keep (5) = lindex
           keep (4) = luir1
        endif
 
c-----------------------------------------------------------------------
c  get seven scalars (transa, nzoff, nblks, presrv, nz, n, ne) from lu
c-----------------------------------------------------------------------
 
c       ne1 = index (lindex), not required for ums2rf
        n1 = index (lindex-1)
        nz1 = index (lindex-2)
c       presr1 = index (lindex-3) .ne. 0, not required for ums2rf
        nblks = index (lindex-4)
c       nzoff1 = index (lindex-5), not required for ums2rf
c       trans1 = index (lindex-6) .ne. 0, not required for ums2rf
 
c-----------------------------------------------------------------------
c  get pointers to permutation vectors
c-----------------------------------------------------------------------
 
        rpermp = (lindex-6) - n
        cpermp = rpermp - n
        ip2 = cpermp - 1
 
c-----------------------------------------------------------------------
c  get pointers to block-triangular information, if btf was used
c-----------------------------------------------------------------------
 
        if (nblks .gt. 1) then
 
c          -------------------------------------------------------------
c          get pointers to btf arrays
c          -------------------------------------------------------------
 
           offpp = cpermp - (n+1)
           blkpp = offpp - (nblks+1)
           lublpp = blkpp - (nblks)
           ip2 = lublpp - 1
           on = n
 
        else
 
c          -------------------------------------------------------------
c          matrix was factorized as a single block, pass dummy arg.
c          -------------------------------------------------------------
 
           offpp = 1
           blkpp = 1
           lublpp = 1
           on = 0
 
        endif
 
        badlu = n .ne. n1 .or. nz1 .le. 0 .or. luir1 .gt. ip2 .or.
     $          nblks .le. 0 .or. nblks .gt. n
        if (badlu) then
           call ums2er (2, icntl, info, -7, 0)
c          error return, lu factors are corrupted:
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  get memory for conversion to column form
c-----------------------------------------------------------------------
 
        nz = ne
        iuse = 2*n+1 + max (2*nz,n+1) + nz + lusiz
        xuse = 2*nz
        info (18) = iuse
        info (20) = xuse
        info (19) = iuse
        info (21) = xuse
        info (23) = xuse
        lind2 = luir1 - 1
        if (lindex .lt. iuse) then
c          set error flag if out of integer memory
           call ums2er (2, icntl, info, -3, iuse)
        endif
        if (lvalue .lt. xuse) then
c          set error flag if out of real memory
           call ums2er (2, icntl, info, -4, xuse)
        endif
        if (info (1) .lt. 0) then
c          error return, if not enough integer and/or real memory
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  convert to column-oriented form and remove duplicates
c-----------------------------------------------------------------------
 
        call ums2co (n, nz, transa, value, lvalue, info, icntl,
     $     index, lind2-(2*n+1), index (lind2-2*n), index (lind2-n), 2)
        if (info (1) .lt. 0) then
c          error return, if all entries are invalid (nz is now 0)
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  current memory usage:
c-----------------------------------------------------------------------
 
c       index (1..n+1): column pointers.  input matrix is now in
c       index (1..nz+n+1) and value (1..nz)
c       col pattern: index (n+1+ index (col) ... n+1+ index (col+1))
c       col values:  value (     index (col) ...      index (col+1))
c       at this point, nz <= ne (nz = ne if there are no invalid or
c       duplicate entries; nz < ne otherwise).
c       pattern of prior lu factors and btf arrays are in
c       index (keep (4) ... keep (5))
 
        iuse = nz + (n+1) + lusiz
        xuse = nz
 
c-----------------------------------------------------------------------
c  refactorize
c-----------------------------------------------------------------------
 
        presrv = job .eq. 1
        if (presrv) then
 
c          -------------------------------------------------------------
c          keep a copy of the original matrix in column-oriented form
c          -------------------------------------------------------------
 
c          copy column pointers (cp (1..n+1) = ap (1..n+1))
           iuse = iuse + (n+1)
cfpp$ nodepchk l
           do 40 i = 1, n+1
              index (nz+n+1+i) = index (i)
40         continue
 
           call ums2r0 (n, nz, index (nz+n+2),
     $          value (nz+1), lvalue-nz,
     $          index (nz+2*n+3), lind2-(nz+2*n+2),
     $          lux1, lui1, iuse, xuse, nzoff, nblks,
     $          icntl, cntl, info, rinfo,
     $          presrv, index, index (n+2), value, n, nz,
     $          index (luir1), ip2 - luir1 + 1,
     $          index (lublpp), index (blkpp), index (offpp), on,
     $          index (cpermp), index (rpermp), ne)
           if (info (1) .lt. 0) then
c             error return, if ums2r0 fails
              go to 9000
           endif
c          adjust pointers to reflect index/value, not ii/xx:
           lux1 = lux1 + nz
           lui1 = lui1 + (nz+2*n+2)
 
c          move preserved copy of a to permanent place
           lux1 = lux1 - (nz)
           lui1 = lui1 - (nz+n+1)
           do 50 i = nz+n+1, 1, -1
              index (lui1+i-1) = index (i)
50         continue
           do 60 i = nz, 1, -1
              value (lux1+i-1) = value (i)
60         continue
 
        else
 
c          -------------------------------------------------------------
c          do not preserve the original matrix
c          -------------------------------------------------------------
 
           call ums2r0 (n, nz, index,
     $          value, lvalue,
     $          index (n+2), lind2-(n+1),
     $          lux1, lui1, iuse, xuse, nzoff, nblks,
     $          icntl, cntl, info, rinfo,
     $          presrv, 1, 1, zero, 0, 1,
     $          index (luir1), ip2 - luir1 + 1,
     $          index (lublpp), index (blkpp), index (offpp), on,
     $          index (cpermp), index (rpermp), ne)
           if (info (1) .lt. 0) then
c             error return, if ums2r0 fails
              go to 9000
           endif
c          adjust pointers to reflect index/value, not ii/xx:
           lui1 = lui1 + (n+1)
 
        endif
 
c-----------------------------------------------------------------------
c  wrap-up
c-----------------------------------------------------------------------
 
        if (transa) then
           index (lindex-6) = 1
        else
           index (lindex-6) = 0
        endif
        index (lindex-5) = nzoff
        index (lindex-4) = nblks
        if (presrv) then
           index (lindex-3) = 1
        else
           index (lindex-3) = 0
        endif
        index (lindex-2) = nz
        index (lindex-1) = n
        index (lindex) = ne
 
c       save location of lu factors
        keep (1) = lux1
        keep (2) = lvalue
        keep (3) = lui1
        keep (4) = luir1
        keep (5) = lindex
 
c       update memory usage information
        iuse = lindex - lui1 + 1
        xuse = lvalue - lux1 + 1
 
c-----------------------------------------------------------------------
c  print the output arguments if requested, and return
c-----------------------------------------------------------------------
 
c       error return label:
9000    continue
        if (info (1) .lt. 0) then
           keep (1) = 0
           keep (2) = 0
           keep (3) = 0
           keep (4) = 0
           keep (5) = 0
        endif
 
        info (18) = min (lindex, max (info (18), iuse))
        info (19) = info (18)
        info (22) = info (19)
        info (20) = min (lvalue, max (info (20), xuse))
 
        call ums2p1 (2, 2,
     $          n, ne, job, transa, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo,
     $          zero, zero, 1, zero, 1)
        return
        end
