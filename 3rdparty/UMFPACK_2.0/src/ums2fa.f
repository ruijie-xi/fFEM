 
        subroutine ums2fa (n, ne, job, transa, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo)
c
cc UMS2FA computes the LU factors of a sparse matrix.
c
        integer n, ne, job, lvalue, lindex, index (lindex), keep (20),
     $          icntl (20), info (40)
        real
     $          value (lvalue), cntl (10), rinfo (20)
        logical transa
 
c=== ums2fa ============================================================
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
c  given a sparse matrix a, find a sparsity-preserving and numerically-
c  acceptable pivot order and compute the lu factors, paq = lu.  the
c  matrix is optionally preordered into a block upper triangular form
c  (btf).  pivoting is performed within each diagonal block to maintain
c  sparsity and numerical stability.  the method used to factorize the
c  matrix is an unsymmetric-pattern variant of the multifrontal method.
c  most of the floating-point work is done in the level-3 blas (dense
c  matrix multiply).  in addition, approximate degrees are used in the
c  markowitz-style pivot search to reduce the symbolic overhead.  for
c  best performance, be sure to use an optimized blas library.
c
c  this routine is normally preceded by a call to ums2in, to
c  initialize the default control parameters.  ums2in need only be
c  called once.  a call to ums2fa can be followed by any number of
c  calls to ums2so, which solves a linear system using the lu factors
c  computed by this routine.  a call to ums2fa can also be followed by
c  any number of calls to ums2rf, which factorizes another matrix with
c  the same nonzero pattern as the matrix factorized by ums2fa (but with
c  different numerical values).
c
c  for more information, see t. a. davis and i. s. duff, "an
c  unsymmetric-pattern multifrontal method for sparse lu factorization",
c  siam j. matrix analysis and applications (to appear), also
c  technical report tr-94-038, cise dept., univ. of florida,
c  p.o. box 116120, gainesville, fl 32611-6120, usa.  the method used
c  here is a modification of that method, described in t. a. davis,
c  "a combined unifrontal/multifrontal method for unsymmetric sparse
c  matrices," tr-94-005.  (technical reports are available via www at
c  http://www.cis.ufl.edu/).  the appoximate degree update algorithm
c  used here has been incorporated into an approximate minimum degree
c  ordering algorithm, desribed in p. amestoy, t. a. davis, and i. s.
c  duff, "an approximate minimum degree ordering algorithm", siam j.
c  matrix analysis and applications (to appear, also tr-94-039).  the
c  approximate minimum degree ordering algorithm is implemented as mc47
c  in the harwell subroutine library (mc47 is not called by umfpack).
 
c=======================================================================
c  installation note:
c=======================================================================
c
c  requires the blas (basic linear algebra subprograms) and two routines
c  from the harwell subroutine library.  ideally, you should use
c  vendor-optimized blas for your computer.  if you do not have them,
c  you may obtain the fortran blas from 1.  send email to
c  netlib@ornl.gov with the two-line message:
c               send index from blas
c               send blas.shar from blas
c
c  to obtain the two harwell subroutine library (hsl) routines, send
c  email to netlib@ornl.gov with the message:
c               send mc21b.f mc13e.f from harwell
c  these two routines hsl contain additional licensing restrictions.
c  if you want to run umfpack without them, see the "installation
c  note:" comment in ums2fb.
c
c  to permamently disable any diagnostic and/or error printing, see
c  the "installation note:" comments in ums2p1 and ums2p2.
c
c  to change the default control parameters, see the
c  "installation note:" comments in ums2in.
 
c=======================================================================
c  arguments:
c=======================================================================
 
c           ------------------------------------------------------------
c  n:       an integer variable.
c           must be set by caller on input (not modified).
c           order of the matrix.  restriction:  1 <= n <= (maxint-5)/3,
c           where maxint is the largest representable positive integer.
 
c           ------------------------------------------------------------
c  ne:      an integer variable.
c           must be set by caller on input (not modified).
c           number of entries in input matrix.  restriction:  ne => 1.
 
c           ------------------------------------------------------------
c  job:     an integer variable.
c           must be set by caller on input (not modified).
c           if job=1, then a column-oriented form of the input matrix
c           is preserved, otherwise, the input matrix is overwritten
c           with its lu factors.  if iterative refinement is to done
c           in ums2so, (icntl (8) > 0), then job must be set to 1.
 
c           ------------------------------------------------------------
c  transa:  an integer variable.
c           must be set by caller on input (not modified).
c           if false then a is factorized: paq = lu.  otherwise, a
c           transpose is factorized:  pa'q = lu.
 
c           ------------------------------------------------------------
c  lvalue:  an integer variable.
c           must be set by caller on input (not modified).
c           size of the value array.  restriction:  lvalue >= 2*ne
c           is required to convert the input form of the matrix into
c           the internal represenation.  lvalue >= ne + axcopy is
c           required to start the factorization, where axcopy = ne if
c           job = 1, or axcopy = 0 otherwise.  during factorization,
c           additional memory is required to hold the frontal matrices.
c           the internal representation of the matrix is overwritten
c           with the lu factors, of size (keep (2) - keep (1) + 1
c           + axcopy), on output.
 
c           ------------------------------------------------------------
c  lindex:  an integer variable.
c           must be set by caller on input (not modified).
c           size of the index array.  restriction: lindex >= 3*ne+2*n+1,
c           is required to convert the input form of the matrix into
c           its internal representation.  lindex >= wlen + alen + acopy
c           is required to start the factorization, where
c           wlen <= 11*n + 3*dn + 8 is the size of the workspaces,
c           dn <= n is the number of columns with more than d
c           entries (d = max (64, sqrt (n)) is the default),
c           alen <= 2*ne + 11*n + 11*dn + dne is the size of the
c           internal representation of the matrix, dne <= ne is the
c           number of entries in such columns with more than d entries,
c           and acopy = ne+n+1 if job = 1, or acopy = 0 otherwize.
c           during factorization, the internal representation of size
c           alen is overwritten with the lu factors, of size
c           luilen = (keep (5) - keep (3) + 1 - acopy) on output.
c           additional memory is also required to hold the unsymmetric
c           quotient graph, but this also overwrites the input matrix.
c           usually about 7*n additional space is adequate for this
c           purpose.  just prior to the end of factorization,
c           lindex >= wlen + luilen + acopy is required.
 
c           ------------------------------------------------------------
c  value:   a real array of size lvalue.
c           must be set by caller on input.  modified on output.  on
c           input, value (1..ne) holds the original matrix in triplet
c           form.  on output, value holds the lu factors, and
c           (optionally) a column-oriented form of the original matrix
c           - otherwise the input matrix is overwritten with the lu
c           factors.
 
c           ------------------------------------------------------------
c  index:   an integer array of size lindex.
c           must be set by caller on input.  modified on output.  on
c           input, index (1..2*ne) holds the original matrix in triplet
c           form.  on output, index holds the lu factors, and
c           (optionally) a column-oriented form of the original matrix
c           - otherwise the input matrix is overwritten with the lu
c           factors.
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
c           on output, the lu factors and the column-oriented form
c           of a (if preserved) are stored in:
c               value (keep (1)...keep (2))
c               index (keep (3)...keep (5))
c           where keep (2) = lvalue, and keep (5) = lindex.
 
c           ------------------------------------------------------------
c  keep:    an integer array of size 20.
c
c           keep (1 ... 5):  need not be set by caller on input.
c               modified on output.
c               keep (1): lu factors start here in value
c               keep (2) = lvalue: lu factors end here in value
c               keep (3): lu factors start here in index
c               keep (4): lu factors needed for ums2rf start here
c                             in index
c               keep (5) = lindex: lu factors end here in index
c
c           keep (6 ... 8):  must be set by caller on input (not
c               modified).
c               integer control arguments not normally modified by the
c               user.  see ums2in for details, which sets the defaults.
c               keep (6) is the largest representable positive
c               integer.  keep (7) and keep (8) determine the
c               size of d, where columns with more than d original
c               entries are treated as a priori frontal matrices.
c
c           keep (9 ... 20): unused.  reserved for future releases.
 
c           ------------------------------------------------------------
c  cntl:    a real array of size 10.
c           must be set by caller on input (not modified).
c           real control arguments, see ums2in for a description,
c           which sets the defaults. ums2fa uses cntl (1) and cntl (2).
 
c           ------------------------------------------------------------
c  icntl:   an integer array of size 20.
c           must be set by caller on input (not modified).
c           integer control arguments, see ums2in for a description,
c           which sets the defaults.  ums2fa uses icntl (1..7).
 
c           ------------------------------------------------------------
c  info:    an integer array of size 40.
c           need not be set by caller on input.  modified on output.
c           it contains information about the execution of ums2fa.
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
c               these entries are ignored and a warning is set
c               in info (1).
c
c           info (4): zero.  used by ums2rf only.
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
c           info (14): garbage collections performed on index, when
c               memory is exhausted.  garbage collections are performed
c               to remove external fragmentation.  if info (14) is
c               excessively high, performance can be degraded.  try
c               increasing lindex if that occurs.  note that external
c               fragmentation in *both* index and value is removed when
c               either is exhausted.
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
c
c           info (18): memory used in index.
c
c           info (19): minimum memory needed in index
c               (or minimum recommended).  if lindex is set to
c               info (19) on a subsequent call, then a moderate
c               number of garbage collections (info (14)) will
c               occur.
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
c           it contains information about the execution of ums2fa.
c
c           rinfo (1): total flop count in the blas
c
c           rinfo (2): total assembly flop count
c
c           rinfo (3): total flops during pivot search
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
c  coding conventions:
c=======================================================================
c
c  this package is written in ansi fortran 77.  to make the code more
c  understandable, the following coding conventions are followed for all
c  routines in this package:
c
c  1) large code blocks are delimited with [...] comments.
c
c  2) goto usage:
c       a) goto's used to return if an error condition is found are
c          written as "go to 9000" or "go to 9010".
c       b) goto's used to exit loops prematurely are written as "go to",
c          and have a target label of 2000 or less.
c       c) goto's used to jump to the next iteration of a do loop or
c          while loop (or to implement a while loop) are written as
c          "goto".
c       no other goto's are used in this package.
c
c  this package uses the following cray compiler directives to help
c  in the vectorization of loops.  each of them operate on the
c  do-loop immediately following the directive.  other compilers
c  normally treat these directives as ordinary comments.
c
c       cfpp$ nodepchk l        disables data dependency check, and
c                               asserts that no recursion exists.
c       cfpp$ nolstval l        disables the saving of last values of
c                               transformed scalars (indexes or promoted
c                               scalars, especially those in array
c                               subscripts).  asserts that values do not
c                               need to be the same as in the scalar
c                               version (for later use of the scalars).
c       cdir$ shortloop         asserts that the loop count is always
c                               64 or less.
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   user routine
c       subroutines called:     ums2er, ums2p1, ums2co, ums2f0
c       functions called:       max, min
        intrinsic max, min
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer i, nz, lux1, lui1, iuse, xuse, luir1, nzoff, nblks,
     $          maxint, nmax, io, prl
        logical presrv
        real
     $          zero
        parameter (zero = 0.0)
 
c  printing control:
c  -----------------
c  io:      i/o unit for diagnostic messages
c  prl:     printing level
c
c  location of lu factors:
c  -----------------------
c  lux1:    real part of lu factors placed in value (lux1 ... lvalue)
c  lui1:    integer part of lu factors placed in index (lui1 ... lindex)
c  luir1:   index (luir1 ... lindex) must be preserved for ums2rf
c
c  memory usage:
c  -------------
c  iuse:    current memory usage in index
c  xuse:    current memory usage in value
c
c  matrix to factorize:
c  --------------------
c  nblks:   number of diagonal blocks (1 if btf not used)
c  nzoff:   entries in off-diagonal part (0 if btf not used)
c  nz:      entries in matrix after removing invalid/duplicate entries
c
c  other:
c  ------
c  maxint:  largest representable positive integer
c  nmax:    largest permissible value of n
c  i:       general loop index
c  presrv:  true if original matrix to be preserved
 
c=======================================================================
c  executable statements:
c=======================================================================
 
        io = icntl (2)
        prl = icntl (3)
 
c-----------------------------------------------------------------------
c  clear informational output, and keep array (except keep (6..8)):
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
        keep (4) = 0
        keep (5) = 0
 
c-----------------------------------------------------------------------
c  print input arguments if requested
c-----------------------------------------------------------------------
 
        call ums2p1 (1, 1,
     $          n, ne, job, transa, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo,
     $          zero, zero, 1, zero, 1)
 
c-----------------------------------------------------------------------
c  initialize and check inputs
c-----------------------------------------------------------------------
 
        iuse = 0
        xuse = 0
        info (5) = ne
        info (6) = ne
        maxint = keep (6)
        nmax = (maxint - 2) / 3
        if (n .lt. 1) then
c          n is too small
           call ums2er (1, icntl, info, -1, -1)
           go to 9000
        endif
        if (n .gt. nmax) then
c          n is too big
           call ums2er (1, icntl, info, -1, nmax)
           go to 9000
        endif
        if (ne .lt. 1) then
c          ne is too small
           call ums2er (1, icntl, info, -2, -1)
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  get memory for conversion to column form
c-----------------------------------------------------------------------
 
        nz = ne
        iuse = 2*n+1 + max (2*nz, n+1) + nz
        xuse = 2*nz
        info (18) = iuse
        info (20) = xuse
        info (19) = iuse
        info (21) = xuse
        if (lindex .lt. iuse) then
c          set error flag if out of integer memory:
           call ums2er (1, icntl, info, -3, iuse)
        endif
        if (lvalue .lt. xuse) then
c          set error flag if out of real memory:
           call ums2er (1, icntl, info, -4, xuse)
        endif
        if (info (1) .lt. 0) then
c          error return, if not enough integer and/or real memory:
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  convert to column-oriented form and remove duplicates
c-----------------------------------------------------------------------
 
        call ums2co (n, nz, transa, value, lvalue, info, icntl,
     $     index, lindex-(2*n+1), index(lindex-2*n), index(lindex-n), 1)
        if (info (1) .lt. 0) then
c          error return, if all entries invalid (nz is now 0):
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
 
        iuse = nz + (n+1)
        xuse = nz
 
c-----------------------------------------------------------------------
c  factorize
c-----------------------------------------------------------------------
 
        presrv = job .eq. 1
        if (presrv) then
 
c          -------------------------------------------------------------
c          keep a copy of the original matrix in column-oriented form
c          -------------------------------------------------------------
 
c          copy column pointers (cp (1..n+1) = ap (1..n+1))
           iuse = iuse + (n+1)
cfpp$ nodepchk l
           do 30 i = 1, n+1
              index (nz+n+1+i) = index (i)
30         continue
 
           call ums2f0 (n, nz, index (nz+n+2),
     $          value (nz+1), lvalue-nz,
     $          index (nz+2*n+3), lindex-(nz+2*n+2),
     $          lux1, lui1, iuse, xuse, nzoff, nblks,
     $          icntl, cntl, info, rinfo,
     $          presrv, index, index (n+2), value, n, nz, keep, ne)
           if (info (1) .lt. 0) then
c             error return, if ums2f0 fails
              go to 9000
           endif
c          adjust pointers to reflect index/value, not ii/xx:
           lux1 = lux1 + nz
           lui1 = lui1 + (nz+2*n+2)
 
c          move preserved copy of a to permanent place
           lux1 = lux1 - nz
           lui1 = lui1 - (nz+n+1)
           do 40 i = nz+n+1, 1, -1
              index (lui1+i-1) = index (i)
40         continue
           do 50 i = nz, 1, -1
              value (lux1+i-1) = value (i)
50         continue
 
        else
 
c          -------------------------------------------------------------
c          do not preserve the original matrix
c          -------------------------------------------------------------
 
           call ums2f0 (n, nz, index,
     $          value, lvalue,
     $          index (n+2), lindex-(n+1),
     $          lux1, lui1, iuse, xuse, nzoff, nblks,
     $          icntl, cntl, info, rinfo,
     $          presrv, 1, 1, zero, 0, 1, keep, ne)
           if (info (1) .lt. 0) then
c             error return, if ums2f0 fails
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
 
c       do not need preserved matrix (n+1+nz), or off-diagonal entries
c       (nzoff) for ums2rf:
        luir1 = lui1
        if (presrv) then
c          do not need preserved matrix for ums2rf
           luir1 = luir1 + n+1 + nz
        endif
        if (nblks .gt. 1) then
c          do not need off-diagonal part for ums2rf
           luir1 = luir1 + nzoff
        endif
 
c       save location of lu factors
        keep (1) = lux1
        keep (2) = lvalue
        keep (3) = lui1
        keep (4) = luir1
        keep (5) = lindex
 
c       update memory usage information
        iuse = lindex - lui1 + 1
        xuse = lvalue - lux1 + 1
        info (22) = info (22) + (lindex - luir1 + 1)
 
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
        info (20) = min (lvalue, max (info (20), xuse))
 
        call ums2p1 (1, 2,
     $          n, ne, job, transa, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo,
     $          zero, zero, 1, zero, 1)
        return
        end
