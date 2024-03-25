 
        subroutine ums2in (icntl, cntl, keep)
c
cc UMS2IN initializes program parameters to default values.
c
        integer icntl (20), keep (20)
        real
     $          cntl (10)
 
c=== ums2in ============================================================
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
c  initialize user-controllable parameters to default values, and
c  non-user controllable parameters.  this routine is normally
c  called once prior to any call to ums2fa.
c
c  this routine sets the default control parameters.  we recommend
c  changing these defaults under certain circumstances:
c
c  (1) if you know that your matrix has nearly symmetric nonzero
c       pattern, then we recommend setting icntl (6) to 1 so that
c       diagonal pivoting is preferred.  this can have a significant
c       impact on the performance for matrices that are essentially
c       symmetric in pattern.
c
c   (2) if you know that your matrix is not reducible to block
c       triangular form, then we recommend setting icntl (4) to 0
c       so that umfpack does not try to permute the matrix to block
c       triangular form (it will not do any useful work and will
c       leave the matrix in its irreducible form).  the work saved
c       is typically small, however.
c
c   the other control parameters typically have less effect on overall
c   performance.
 
c=======================================================================
c  installation note:
c=======================================================================
c
c  this routine can be modified on installation to reflect the computer
c  or environment in which this package is installed (printing control,
c  maximum integer, block size, and machine epsilon in particular).  if
c  you, the installer, modify this routine, please comment out the
c  original code, and add comments (with date) to describe the
c  installation.  do not delete any original code.
 
c=======================================================================
c  arguments:
c=======================================================================
 
c               --------------------------------------------------------
c  icntl:       an integer array of size 20.  need not be set by
c               caller on input.  on output, it contains default
c               integer control parameters.
c
c  icntl (1):   fortran output unit for error messages.
c               default: 6
c
c  icntl (2):   fortran output unit for diagnostic messages.
c               default: 6
c
c  icntl (3):   printing-level.
c               0 or less: no output
c               1: error messages only
c               2: error messages and terse diagnostics
c               3: as 2, and print first few entries of all input and
c                       output arguments.  invalid and duplicate entries
c                       are printed.
c               4 or more: as 2, and print all entries of all input and
c                       output arguments.  invalid and duplicate entries
c                       are printed.
c               default: 2
c
c  icntl (4):   whether or not to attempt a permutation to block
c               triangular form.  if nonzero, then attempt the
c               permutation.  if you know the matrix is not reducible
c               to block triangular form, then setting icntl (4) to
c               zero can save a small amount of computing time.
c               default: 1 (attempt the permutation)
c
c  icntl (5):   the number of columns to examine during the global
c               pivot search.  a value less than one is treated as one.
c               default: 4
c
c  icntl (6):   if not equal to zero, then pivots from the diagonal
c               of a (or the diagonal of the block-triangular form) are
c               preferred.  if the nonzero pattern of the matrix is
c               basically symmetric, we recommend that you change this
c               default value to 1 so that pivots on the diagonal
c               are preferred.
c               default: 0 (do not prefer the diagonal)
c
c  icntl (7):   block size for the blas, controlling the tradeoff
c               between the level-2 and level-3 blas.  values less than
c               one are treated as one.
c               default: 16, which is suitable for the cray ymp.
c
c  icntl (8):   number of steps of iterative refinement to perform.
c               values less than zero are treated as zero.  the matrix
c               must be preserved for iterative refinement to be done
c               (job=1 in ums2fa or ums2rf).
c               default: 0  (no iterative refinement)
c
c  icntl (9 ... 20):  set to zero.  reserved for future releases.
 
c               --------------------------------------------------------
c  cntl:        a real array of size 10.
c               need not be set by caller on input.  on output, contains
c               default real control parameters.
c
c  cntl (1):    pivoting tradeoff between sparsity-preservation
c               and numerical stability.  an entry a(k,k) is numerically
c               acceptable if:
c                  abs (a(k,k)) >= cntl (1) * max (abs (a(*,k)))
c               values less than zero are treated as zero (no numerical
c               constraints).  values greater than one are treated as
c               one (partial pivoting with row interchanges).
c               default: 0.1
c
c  cntl (2):    amalgamation parameter.  if the first pivot in a
c               frontal matrix has row degree r and column degree c,
c               then a working array of size
c                  (cntl (2) * c) - by - (cntl (2) * r)
c               is allocated for the frontal matrix.  subsequent pivots
c               within the same frontal matrix must fit within this
c               working array, or they are not selected for this frontal
c               matrix.  values less than one are treated as one (no
c               fill-in due to amalgamation).  some fill-in due to
c               amalgamation is necessary for efficient use of the blas
c               and to reduce the assembly operations required.
c               default: 2.0
c
c  cntl (3):    normally not modified by the user.
c               defines the smallest positive number,
c               epsilon = cntl (3), such that fl (1.0 + epsilon)
c               is greater than 1.0 (fl (x) is the floating-point
c               representation of x).  if the floating-point mantissa
c               is binary, then cntl (3) is 2 ** (-b+1), where b
c               is the number of bits in the mantissa (including the
c               implied bit, if applicable).
c
c               typical defaults:
c               for ieee double precision, cntl (3) = 2 ** (-53+1)
c               for ieee single precision, cntl (3) = 2 ** (-24+1)
c               for cray double precision, cntl (3) = 2 ** (-96+1)
c               for cray single precision, cntl (3) = 2 ** (-48+1)
c
c               a value of cntl (3) less than or equal to zero
c               or greater than 2 ** (-15) is treated as 2 ** (-15),
c               which assumes that any floating point representation
c               has at least a 16-bit mantissa.  cntl (3) is only
c               used in ums2s2 to compute the sparse backward error
c               estimates, rinfo (7) and rinfo (8), when
c               icntl (8) > 0 (the default is icntl (8) = 0,
c               so by default, cntl (3) is not used).
c
c  cntl (4 ... 10):  set to zero.  reserved for future releases.
 
c               --------------------------------------------------------
c  keep:        an integer array of size 20.
c               need not be set by the caller.  on output, contains
c               integer control parameters that are (normally) non-user
c               controllable (but can of course be modified by the
c               "expert" user or library installer).
c
c  keep (1 ... 5):  unmodified (see ums2fa or ums2rf for a description).
c
c  keep (6):    largest representable positive integer.  set to
c               2^31 - 1 = 2147483647 for 32-bit machines with 2's
c               complement arithmetic (the usual case).
c               default: 2147483647
c
c  keep (7) and keep (8): a column is treated as "dense" if
c               it has more than
c               max (0, keep(7), keep(8)*int(sqrt(float(n))))
c               original entries.  "dense" columns are treated
c               differently that "sparse" rows and columns.  dense
c               columns are transformed into a priori contribution
c               blocks of dimension cdeg-by-1, where cdeg is the number
c               of original entries in the column.  modifying these two
c               parameters can change the pivot order.
c               default:  keep (7) = 64
c               default:  keep (8) = 1
c
c  keep (9 ... 20):  set to zero.  reserved for future releases.
 
c## end of user documentation ##########################################
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   user routine
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer i
        real
     $          zero, tenth, two
        parameter (tenth = 0.1, two = 2.0, zero = 0.0)
 
c  i:       loop index
 
c=======================================================================
c  executable statments:
c=======================================================================
 
c       ----------------------------------------------------------------
c       integer control parameters:
c       ----------------------------------------------------------------
 
        icntl (1) = 6
        icntl (2) = 6
        icntl (3) = 2
        icntl (4) = 1
        icntl (5) = 4
        icntl (6) = 0
        icntl (7) = 16
        icntl (8) = 0
 
c       icntl (9 ... 20) is reserved for future releases:
        do i = 9, 20
           icntl (i) = 0
        end do
 
c       ----------------------------------------------------------------
c       real control parameters:
c       ----------------------------------------------------------------
 
        cntl (1) = tenth
        cntl (2) = two
 
c       ieee single precision:  epsilon = 2 ** (-24)
        cntl (3) = two ** (-23)
 
c       cntl (4 ... 10) is reserved for future releases:
        do 30 i = 4, 10
           cntl (i) = zero
30      continue
 
c       ----------------------------------------------------------------
c       integer control parameters in keep:
c       ----------------------------------------------------------------
 
        keep (6) = 2147483647
        keep (7) = 64
        keep (8) = 1
 
c       keep (9 ... 20) is reserved for future releases:
        do 20 i = 9, 20
           keep (i) = 0
20      continue
 
        return
        end
