 
        subroutine ums2r2 (cp, nz, n, xtail, xx, xsize, xuse, ari,
     $          cperm, rperm, icntl, cntl, info, rinfo, mc, mr,
     $          wir, wic, wpr, wpc, wm, wj, frdimc, frxp, frnext,
     $          frprev, nlu, lup, lui, noutsd, xrmax)
c
cc UMS2R2 is a utility which refactors part of a matrix.
c
        integer xsize, icntl (20), info (40), cperm (n), rperm (n),
     $          xtail, nz, n, ari (nz), cp (n+1), mr, mc, noutsd,
     $          wir (n), wic (n), wpr (mr), xrmax, wpc (mc), wm (mc),
     $          nlu, frdimc (nlu+2), frxp (nlu+2), xuse, wj (mc),
     $          frnext (nlu+2), frprev (nlu+2), lup (nlu), lui (*)
        real
     $          xx (xsize), cntl (10), rinfo (20)
 
c=== ums2r2 ============================================================
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
c  ums2r2 refactorizes the n-by-n input matrix at the head of xx
c  (in arrowhead form) and places its lu factors at the tail of
c  xx.  the input matrix is overwritten.   no btf information is
c  used in this routine.
 
c=======================================================================
c  input:
c=======================================================================
c
c       cp (1..n+1):    column pointers of arrowhead form
c       n:              order of input matrix
c       nz:             entries in input matrix
c       xsize:          size of xx
c       icntl:          integer control parameters, see ums2in
c       cntl:           real control parameters, see ums2in
c
c       ari (1..nz):            arrowhead format of a
c       xx (1..nz):             arrowhead format of a, see below
c       xx (nz+1..xsize):       undefined on input, used as workspace
c
c       nlu:            number of lu arrowheads
c       lup (1..nlu):   pointers to lu arrowheads in lui
c       lui (1.. ):     lu arrowheads
c
c       xuse:           memory usage in value
c
c       noutsd:         entries not in prior lu pattern
c
c       cperm (1..n):   column permutation
c       rperm (1..n):   row permutation
 
c=======================================================================
c  workspace:
c=======================================================================
c
c       wir (1..n)
c       wic (1..n)
c
c       wpr (1.. max ludegr)
c       wpc (1.. max ludegc)
c       wm  (1.. max ludegc)
c       wj  (1.. max ludegc)
c
c       frdimc (1..nlu+2)
c       frxp   (1..nlu+2)
c       frnext (1..nlu+2)
c       frprev (1..nlu+2)
 
c=======================================================================
c  output:
c=======================================================================
c
c       lui (1..):              lu arrowheads, modified luxp pointers
c       xx (1..xtail-1):        undefined on output
c       xx (xtail..xsize):      lu factors of this matrix, see below
c
c       info:           integer informational output, see ums2fa
c       rinfo:          real informational output, see ums2fa
c
c       xuse:           memory usage in value
c
c       noutsd:         entries not in prior lu pattern, incremented
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2r0
c       subroutines called:     ums2er, ums2p2, ums2rg, sgemv,
c                               sgemm, strsv, strsm
c       functions called:       abs, max
        intrinsic abs, max
 
c=======================================================================
c  description of data structures:
c=======================================================================
 
c-----------------------------------------------------------------------
c  matrix being factorized:
c-----------------------------------------------------------------------
c
c  the input matrix is held in an arrowhead format.  for the kth pivot,
c  the nonzeros in the pivot row (a (k, k...n)) and pivot column
c  (a (k...n, k)) are stored in the kth arrowhead.  the kth arrowhead
c  is located in:
c       ari (cp (k+1) ... cp (k)-1):    pattern
c       xx  (cp (k+1) ... cp (k)-1):    values
c
c  suppose p is in the range cp (k+1) to cp (k)-1.  if ari (p) is
c  greater than zero, then the entry is in row ari (p), column k,
c  with value xx (p).  if ari (p) is less than zero, then the entry is
c  in row k, column -ari (p), with value xx (p).  the arrowheads are
c  stored in reverse order (arrowhead n, n-1, ... 2, 1) in ari and xx.
c  note that cp (n+1) = 1 unless btf is in use and the original matrix
c  is not preserved.   in all cases, the real part of the arrowhead
c  format (xx (cp (n+1) ... cp (1)-1)) is overwritten with the lu
c  factors.  the integer part (ari (cp (n+1) ... cp (1)-1)) is not
c  overwritten, since ums2r2 does not require dynamic allocation of
c  integer memory.
 
c-----------------------------------------------------------------------
c  frontal matrices
c-----------------------------------------------------------------------
c
c   each unassembled frontal matrix (element) is stored as follows:
c       total size: fscal integers, (fdimr*fdimc) reals
c
c       if e is an unassembled element, and not the current frontal
c       matrix:
c
c       fluip = lup (e) pointer to lu arrowhead in ii
c       fdimc = frdimc (e)      column dimension of contribution block
c       fxp   = frxp (e)        pointer to contribution block in xx
c       next  = frnext (e)      pointer to next block in xx
c       prev  = frprev (e)      pointer to previous block in xx
c       fdegr = abs (lui (fluip+2))
c       fdegc = abs (lui (fluip+2))
c       xx (fxp ... )
c               a 2-dimensional array, c (1..fdimc, 1..fdimr), where
c               fdimr = fdegr if the contribution block is compressed,
c               or fdimr = lui (fluip+5) if not.  note, however, that
c               fdimr is not needed.  the contribution block is stored
c               in c (1..fdegc, 1..fdegr) in the c (1..fdimc,...) array.
c
c               if memory is limited, garbage collection will occur.
c               in this case, the c (1..fdimc, 1..fdimr) array is
c               compressed to be just large enough to hold the
c               unassembled contribution block,
c               c (1..fdegc, 1..fdegr).
 
c-----------------------------------------------------------------------
c  current frontal matrix
c-----------------------------------------------------------------------
c
c  ffxp points to current frontal matrix (contribution block and lu
c  factors).  for example, if fflefc = 4, fflefr = 6, luk = 3,
c  ffdimc = 8, ffdimr = 12, then "x" is a term in the contribution
c  block, "l" in l1, "u" in u1, "l" in l2, "u" in u2, and "." is unused.
c  xx (fxp) is "x". the first 3 pivot values (diagonal entries in u1)
c  are labelled 1, 2, and 3.  the frontal matrix is ffdimc-by-ffdimr.
c
c                   |----------- col 1 of l1 and l2, etc.
c                   v
c       x x x x x x l l l . . .
c       x x x x x x l l l . . .
c       x x x x x x l l l . . .
c       x x x x x x l l l . . .
c       u u u u u u 3 l l . . .         <- row 3 of u1 and u2
c       u u u u u u u 2 l . . .         <- row 2 of u1 and u2
c       u u u u u u u u 1 . . .         <- row 1 of u1 and u2
c       . . . . . . . . . . . .
 
c-----------------------------------------------------------------------
c  lu factors
c-----------------------------------------------------------------------
c
c   the lu factors are placed at the tail of xx.  if this routine
c   is factorizing a single block, then this description is for the
c   factors of the single block:
c
c       lui (1..):      integer info. for lu factors
c       xx (xtail..xsize):      real values in lu factors
c
c   each lu arrowhead (or factorized element) is stored as follows:
c   ---------------------------------------------------------------
c
c       total size: (7 + ludegc + ludegr + lunson) integers,
c                   (luk**2 + ludegc*luk + luk*ludegc) reals
c
c       if e is an lu arrowhead, then luip = lup (e).
c
c       luxp   = lui (luip) pointer to numerical lu arrowhead
c       luk    = lui (luip+1) number of pivots in lu arrowhead
c       ludegr = lui (luip+2) degree of last row of u (excl. diag)
c       ludegc = lui (luip+3) degree of last col of l (excl. diag)
c       lunson = lui (luip+4) number of children in assembly dag
c       ffdimr = lui (luip+5)
c       ffdimc = lui (luip+6)
c                       max front size for this lu arrowhead is
c                       ffdimr-by-ffdimc, or zero if this lu arrowhead
c                       factorized within the frontal matrix of a prior
c                       lu arrowhead.
c       lucp   = (luip + 7)
c                       pointer to pattern of column of l
c       lurp   = lucp + ludegc
c                       pointer to patter of row of u
c       lusonp = lurp + ludegr
c                       pointer to list of sons in the assembly dag
c       lui (lucp ... lucp + ludegc - 1)
c                       row indices of column of l
c       lui (lurp ... lurp + ludegr - 1)
c                       column indices of row of u
c       lui (lusonp ... lusonp + lunson - 1)
c                       list of sons
c       xx (luxp...luxp + luk**2 + ludegc*luk + luk*ludegr - 1)
c                       pivot block (luk-by-luk) and the l block
c                       (ludegc-by-luk) in a single (luk+ludegc)-by-luk
c                       array, followed by the u block in a
c                       luk-by-ludegr array.
c
c   pivot column/row pattern (also columns/rows in contribution block):
c       if the column/row index is negated, the column/row has been
c       assembled out of the frontal matrix into a subsequent frontal
c       matrix.  after factorization, the negative flags are removed.
c
c   list of sons:
c       1 <= son <= n:           son an luson
c       n+1 <= son <= 2n:        son-n is an uson
c       2n+n <= son <= 3n:       son-2n is a lson
 
c-----------------------------------------------------------------------
c  workspaces:
c-----------------------------------------------------------------------
c
c  wpc (1..ludegr):     holds the pivot column pattern
c                       (excluding the pivot row indices)
c
c  wpr (1..ludegr):     holds the pivot row pattern
c                       (excluding the pivot column indices)
c
c  wir (row) >= 0 for each row in pivot column pattern.
c               offset into pattern is given by:
c               wir (row) == offset - 1
c               otherwise, wir (1..n) is < 0
c
c  wic (col) >= 0 for each col in pivot row pattern.
c               wic (col) == (offset - 1) * ffdimc
c               otherwise, wic (1..n) is < 0
c
c  wm (1..degc) or wm (1..fdegc):       a gathered copy of wir
c  wj (1..degc) or wj (1..fdegc):       offset in pattern of a son
 
c-----------------------------------------------------------------------
c  memory allocation in xx:
c-----------------------------------------------------------------------
c
c   xx (1..xhead):      values of original entries in arrowheads of
c                       matrix, values of contribution blocks, followed
c                       by the current frontal matrix.
c
c   mtail = nlu+2
c   mhead = nlu+1:      frnext (mhead) points to the first contribution
c                       block in the head of xx.  the frnext and frprev
c                       arrays form a doubly-linked list.  traversing
c                       the list from mhead to mtail gives the
c                       contribution blocks in ascending ordering of
c                       address (frxp).  a block is free if frdimc <= 0.
c                       the largest known free block in xx is pfree,
c                       located in
c                       xx (frxp (pfree) ... frxp (pfree) + xfree -1),
c                       unless pfree = 0, in which case no largest free
c                       block is known.
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer swpcol, swprow, fdimc, k0, colpos, rowpos, pivot, ffpp,
     $          p, i, j, ludegr, ludegc, kpos, sp, ffrp, ffcp, type,
     $          fxp, lurp, lucp, next, fflefr, prev, xhead, fdegr,
     $          fflefc, k, xcdp, xdp, xsp, s, fdegc, flurp, flucp,
     $          col, e, row, mhead, mtail, uxp, luk, io, fluip, lusonp,
     $          ffsize, ffxp, ffdimr, ffdimc, xrdp, npiv, nb, lunson,
     $          xneed, ldimr, ldimc, lxp, prl, xp, luip, pfree, xfree,
     $          xs, luxp, fsp, flp, fdp, degc, nzu, nzl, xruse
        logical pr3, allcol, allrow
        real
     $          one, zero, x
        parameter (one = 1.0, zero = 0.0)
 
c  printing control:
c  -----------------
c  prl:     invalid entries printed if prl >= 3
c  io:      i/o unit for warning messages (printing invalid entries)
c  pr3:     true if invalid entries are to be printed when found
c
c  current working array:
c  ----------------------
c  ffxp:    current working array is in xx (ffxp ... ffxp+ffsize-1)
c  ffsize:  size of current working array in xx
c  ffdimr:  row degree (number of columns) of current working array
c  ffdimc:  column degree (number of rows) of current working array
c  fflefr:  row degree (number of columns) of current contribution block
c  fflefc:  column degree (number of rows) of current contribution block
c  ffrp:     u2 block is in xx (ffrp ...)
c  ffcp:     l2 block is in xx (ffcp ...)
c  ffpp:     location in xx of the current pivot value
c
c  current element:
c  ----------------
c  s:       current element being factorized
c  luip:    current element is in lui (luip ...)
c  luk:     number of pivots in current element
c  ludegc:  degree of pivot column (excluding pivots themselves)
c  ludegr:  degree of pivot row (excluding pivots themselves)
c  ldimr:   row degree (number of columns) of current element
c  ldimc:   column degree (number of row) of current element
c  lucp:    pattern of col(s) of current element in lui (lucp...)
c  lurp:    pattern of row(s) of current element in lui (lurp...)
c  lusonp:  list of sons of current element is in lui (lusonp...)
c  lunson:  number of sons of current element
c  sp:      pointer into list of sons of current element
c  luxp:    numerical values of lu arrowhead stored in xx (luxp ...)
c  lxp:     l2 block is stored in xx (lxp ...) when computed
c  uxp:     u2 block is stored in xx (uxp ...) when computed
c  nzu:     nonzeros above diagonal in u in current lu arrowhead
c  nzl:     nonzeros below diagonal in l in current lu arrowhead
c  swpcol:  the non-pivotal column to be swapped with pivot column
c  swprow:  the non-pivotal row to be swapped with pivot row
c  colpos:  position in wpr of the pivot column
c  rowpos:  position in wpc of the pivot row
c  kpos:    position in c to place pivot row/column
c  k:       current pivot is kth pivot of current element, k=1..luk
c  k0:      contribution block, c, has been updated with pivots 1..k0
c  npiv:    number of pivots factorized so far, excl. current element
c  pivot:   current pivot entry is a (pivot, pivot)
c  xcdp:    current pivot column is in xx (xcdp ...)
c  xrdp:    current pivot row is in xx (xrdp ...)
c
c  son, or element other than current element:
c  -------------------------------------------
c  e:       an element other than s (a son of s, for example)
c  fluip:   lu arrowhead of e is in lui (fluip ...)
c  fxp:     contribution block of son is in xx (fxp ...)
c  fdimc:   leading dimension of contribution block of a son
c  fdegr:   row degree of contribution block of son (number of columns)
c  fdegc:   column degree of contribution block of son (number of rows)
c  allcol:  true if all columns are present in son
c  allrow:  true if all rows are present in son
c  flucp:   pattern of col(s) of son in lui (flucp...)
c  flurp:   pattern of row(s) of son in lui (flurp...)
c  type:    an luson (type = 1), uson (type = 2) or lson (type = 3)
c  degc:    compressed column offset vector of son is in wj/wm (1..degc)
c
c  memory allocation:
c  ------------------
c  mhead:   nlu+1, head pointer for contribution block link list
c  mtail:   nlu+2, tail pointer for contribution block link list
c  prev:    frprev (e) of the element e
c  next:    frnext (e) of the element e
c  pfree:   frxp (pfree) is the largest known free block in xx
c  xfree:   size of largest known free block in xx
c  xneed:   bare minimum memory currently needed in xx
c  xhead:   xx (1..xhead-1) is in use, xx (xhead ..) is free
c  xruse:   estimated memory needed in xx for next call to ums2rf,
c           assuming a modest number of garbage collections
c  xs:      size of a block of memory in xx
c
c  other:
c  ------
c  xdp:     destination pointer, into xx
c  xsp:     source pointer, into xx
c  xp:      a pointer into xx
c  fsp:     source pointer, into xx
c  fsp:     destination pointer, into xx
c  flp:     last row/column in current contribution is in xx (flp...)
c  col:     a column index
c  row:     a row index
c  nb:      block size for tradeoff between level-2 and level-3 blas
c  p, i, j, x:  various uses
 
c=======================================================================
c  executable statements:
c=======================================================================
 
c       ----------------------------------------------------------------
c       get control parameters and initialize various scalars
c       ----------------------------------------------------------------
 
        io = icntl (2)
        prl = icntl (3)
        nb = max (1, icntl (7))
        npiv = 0
        xhead = cp (1)
        xtail = xsize + 1
        xneed = xuse
        xruse = xuse
        xrmax = max (xrmax, xruse)
        mhead = nlu+1
        mtail = nlu+2
        xfree = -1
        pfree = 0
        pr3 = prl .ge. 3 .and. io .ge. 0
 
c       ----------------------------------------------------------------
c       initialize workspaces
c       ----------------------------------------------------------------
 
        do i = 1, n
           wir (i) = -1
           wic (i) = -1
        end do
 
        do 20 e = 1, nlu+2
           frdimc (e) = 0
           frxp (e) = 0
           frnext (e) = 0
           frprev (e) = 0
20      continue
        frnext (mhead) = mtail
        frprev (mtail) = mhead
        frxp (mhead) = xhead
        frxp (mtail) = xhead
 
c       count the numerical assembly of the original matrix
        rinfo (2) = rinfo (2) + nz
 
c       current working array is empty:
        fflefr = 0
        fflefc = 0
        ffsize = 0
        ffxp = xhead
 
c=======================================================================
c  factorization [
c=======================================================================
 
        do 600 s = 1, nlu
 
c=======================================================================
c  get the next element to factorize
c=======================================================================
 
           luip = lup (s)
           luk = lui (luip+1)
           ludegc = lui (luip+3)
           ludegr = lui (luip+2)
           lunson = lui (luip+4)
           lucp = (luip + 7)
           lurp = lucp + ludegc
           lusonp = lurp + ludegr
           ldimc = luk + ludegc
           ldimr = luk + ludegr
 
c=======================================================================
c  start new frontal matrix or merge with prior contribution block [
c=======================================================================
 
c          =============================================================
           if (lui (luip+6) .ne. 0) then
c          start new contribution block
c          =============================================================
 
c             ----------------------------------------------------------
c             clear the prior offsets
c             ----------------------------------------------------------
 
              do 30 i = 1, fflefr
                 wic (wpr (i)) = -1
30            continue
              do 40 i = 1, fflefc
                 wir (wpc (i)) = -1
40            continue
 
c             ----------------------------------------------------------
c             save prior contribution block (s-1), if it exists
c             ----------------------------------------------------------
 
              xs = fflefr * fflefc
              if (ffsize .ne. 0) then
c                one more frontal matrix is finished
                 xneed = xneed - (ffsize - xs)
                 xruse = xruse - (ffsize - xs)
                 info (13) = info (13) + 1
c             else
c                prior contribution block does not exist
              endif
 
              if (fflefr .le. 0 .or. fflefc .le. 0) then
 
c                -------------------------------------------------------
c                if prior contribution block nonexistent or empty
c                -------------------------------------------------------
 
                 xuse = xuse - (xhead - frxp (mtail))
                 xhead = frxp (mtail)
 
              else
 
c                -------------------------------------------------------
c                prepare the prior contribution block for later assembly
c                -------------------------------------------------------
 
                 e = s - 1
 
c                count the numerical assembly
                 rinfo (2) = rinfo (2) + xs
 
                 if (xs .le. xfree) then
 
c                   ----------------------------------------------------
c                   compress and store in a freed block
c                   ----------------------------------------------------
 
c                   place the new block in the list
                    xfree = xfree - xs
                    if (pfree .eq. mtail) then
c                      place the new block at start of tail block
                       prev = frprev (mtail)
                       next = mtail
                       xdp = frxp (mtail)
                       frxp (mtail) = xdp + xs
                    else
c                      place the new block at end of block
                       prev = pfree
                       next = frnext (pfree)
                       xdp = frxp (next) - xs
                       if (xfree .eq. 0 .and. pfree .ne. mhead) then
c                         delete the free block if its size is zero
                          prev = frprev (prev)
                          pfree = 0
                          xfree = -1
                       endif
                    endif
                    do 60 j = 0, fflefr - 1
cfpp$ nodepchk l
                       do 50 i = 0, fflefc - 1
                          xx (xdp+j*fflefc+i) = xx (ffxp+j*ffdimc+i)
50                     continue
60                  continue
                    xuse = xuse - (xhead - frxp (mtail))
                    xhead = frxp (mtail)
                    frxp (e) = xdp
                    frdimc (e) = fflefc
 
                 else
 
c                   ----------------------------------------------------
c                   deallocate part of unused portion of frontal matrix
c                   ----------------------------------------------------
 
c                   leave the contribution block c (1:fflefc, 1:fflefr)
c                   at head of xx, with column dimension of ffdimc and
c                   space of size (fflefr-1)*ffdimc for the first
c                   fflefr columns, and fflefc for the last column.
                    xs = ffsize - (fflefc + (fflefr-1)*ffdimc)
                    xhead = xhead - xs
                    xuse = xuse - xs
                    prev = frprev (mtail)
                    next = mtail
                    frxp (mtail) = xhead
                    frxp (e) = ffxp
                    frdimc (e) = ffdimc
                 endif
 
                 frnext (prev) = e
                 frprev (next) = e
                 frnext (e) = next
                 frprev (e) = prev
 
              endif
 
              if (pfree .eq. mtail) then
                 pfree = 0
                 xfree = -1
              endif
 
c             ----------------------------------------------------------
c             allocate a new ffdimr-by-ffdimc frontal matrix
c             ----------------------------------------------------------
 
              ffdimc = lui (luip+6)
              ffdimr = lui (luip+5)
              ffsize = ffdimr * ffdimc
              ffxp = 0
 
c             ----------------------------------------------------------
c             allocate and zero the space, garbage collection if needed
c             ----------------------------------------------------------
 
              if (ffsize .gt. xtail-xhead) then
                 info (15) = info (15) + 1
                 call ums2rg (xx, xsize, xhead, xtail, xuse,
     $              lui, frdimc, frxp, frnext, frprev, nlu, lup,
     $              icntl, ffxp, ffsize, pfree, xfree)
              endif
 
              ffxp = xhead
              xhead = xhead + ffsize
              xuse = xuse + ffsize
              xneed = xneed + ffsize
              xruse = xruse + ffsize
              xrmax = max (xrmax, xruse)
              info (20) = max (info (20), xuse)
              info (21) = max (info (21), xneed)
              if (xhead .gt. xtail) then
c                error return, if not enough real memory:
                 go to 9000
              endif
 
c             ----------------------------------------------------------
c             zero the frontal matrix
c             ----------------------------------------------------------
 
              do 70 p = ffxp, ffxp + ffsize - 1
                 xx (p) = zero
70            continue
 
c             ----------------------------------------------------------
c             place pivot rows and columns in correct position
c             ----------------------------------------------------------
 
              do 80 k = 1, luk
                 wic (npiv + k) = (ldimr - k) * ffdimc
                 wir (npiv + k) =  ldimc - k
80            continue
 
c             ----------------------------------------------------------
c             get the pivot row pattern of the new lu arrowhead
c             ----------------------------------------------------------
 
              do 90 i = 0, ludegr - 1
                 col = lui (lurp+i)
                 wic (col) = i * ffdimc
                 wpr (i+1) = col
90            continue
 
c             ----------------------------------------------------------
c             get the pivot column pattern of the new lu arrowhead
c             ----------------------------------------------------------
 
              do 100 i = 0, ludegc - 1
                 row = lui (lucp+i)
                 wir (row) = i
                 wpc (i+1) = row
100           continue
 
c          =============================================================
           else
c          merge with prior contribution block
c          =============================================================
 
c             ----------------------------------------------------------
c             prior block is located at xx (ffxp ... ffxp + ffsize - 1).
c             it holds a working array c (1..ffdimc, 1..ffdimr), with a
c             prior contribution block in c (1..fflefc, 1..fflefr).
c             the last pivot column pattern is wpc (1..fflefc), and
c             the last pivot row pattern is wpr (1..fflefr).  the
c             offsets wir and wic are:
c             wir (wpc (i)) = i-1, for i = 1..fflefc, and -1 otherwise.
c             wic (wpr (i)) = (i-1)*ffdimc, for i = 1..fflefr, else -1.
c             the prior lu arrowhead is an implicit luson of the current
c             element (and is implicitly assembled into the same
c             frontal matrix).
c             ----------------------------------------------------------
 
c             ----------------------------------------------------------
c             zero the newly extended frontal matrix
c             ----------------------------------------------------------
 
c             zero the new columns in the contribution and lu blocks
c             c (1..ldimc, fflefr+1..ldimr) = 0
              do 120 j = fflefr, ldimr - 1
                 do 110 i = 0, ldimc - 1
                    xx (ffxp + j*ffdimc + i) = zero
110              continue
120           continue
 
c             c (fflefc+1..ldimc, 1..fflefr) = 0
c             zero the new rows in the contribution and u blocks
              do 140 i = fflefc, ldimc - 1
cfpp$ nodepchk l
                 do 130 j = 0, fflefr - 1
                    xx (ffxp + j*ffdimc + i) = zero
130              continue
140           continue
 
c             ----------------------------------------------------------
c             move pivot rows and columns into correct position
c             ----------------------------------------------------------
 
              do 220 k = 1, luk
 
c                -------------------------------------------------------
c                kth pivot of frontal matrix, (npiv+k)th pivot of lu
c                -------------------------------------------------------
 
                 pivot = npiv + k
 
c                -------------------------------------------------------
c                move the kth pivot column into position
c                -------------------------------------------------------
 
                 xsp = wic (pivot)
                 kpos = ldimr - k + 1
                 xdp = (kpos - 1) * ffdimc
                 wic (pivot) = xdp
 
                 if (xsp .ge. 0) then
c                   pivot column is already in current frontal matrix,
c                   shift into proper position
                    colpos = (xsp / ffdimc) + 1
                    fsp = ffxp + xsp
                    fdp = ffxp + xdp
 
                    if (fflefr .lt. kpos) then
 
                       if (fflefr .eq. colpos) then
 
c                         ----------------------------------------------
c                         move c(:,colpos) => c (:,kpos)
c                         c (:,colpos) = 0
c                         ----------------------------------------------
cfpp$ nodepchk l
                          do 150 i = 0, ldimc - 1
                             xx (fdp+i) = xx (fsp+i)
                             xx (fsp+i) = zero
150                       continue
 
                       else
 
c                         ----------------------------------------------
c                         move c(:,colpos) => c (:,kpos)
c                         move c(:,fflefr) => c (:,colpos)
c                         c (:,fflefr) = 0
c                         ----------------------------------------------
 
                          flp = ffxp + (fflefr - 1) * ffdimc
cfpp$ nodepchk l
                          do 160 i = 0, ldimc - 1
                             xx (fdp+i) = xx (fsp+i)
                             xx (fsp+i) = xx (flp+i)
                             xx (flp+i) = zero
160                       continue
 
                          swpcol = wpr (fflefr)
                          wpr (colpos) = swpcol
                          wic (swpcol) = xsp
                       endif
 
                    else if (colpos .ne. kpos) then
 
c                      -------------------------------------------------
c                      swap c (:,colpos) <=> c (:,kpos)
c                      -------------------------------------------------
cfpp$ nodepchk l
                       do 180 i = 0, ldimc - 1
                          x = xx (fdp+i)
                          xx (fdp+i) = xx (fsp+i)
                          xx (fsp+i) = x
180                    continue
 
                       swpcol = wpr (kpos)
                       wpr (colpos) = swpcol
                       wic (swpcol) = xsp
                    endif
 
                    fflefr = fflefr - 1
                 endif
 
c                -------------------------------------------------------
c                move the kth pivot row into position
c                -------------------------------------------------------
 
                 xsp = wir (pivot)
                 kpos = ldimc - k + 1
                 xdp = (kpos - 1)
                 wir (pivot) = xdp
 
                 if (xsp .ge. 0) then
c                   pivot row is already in current frontal matrix,
c                   shift into proper position
                    rowpos = xsp + 1
                    fsp = ffxp + xsp
                    fdp = ffxp + xdp
 
                    if (fflefc .lt. kpos) then
 
                       if (fflefc .eq. rowpos) then
 
c                         ----------------------------------------------
c                         move c(rowpos,:) => c (kpos,:)
c                         c (rowpos,:) = 0
c                         ----------------------------------------------
cfpp$ nodepchk l
                          do 190 j = 0, (ldimr - 1) * ffdimc, ffdimc
                             xx (fdp+j) = xx (fsp+j)
                             xx (fsp+j) = zero
190                       continue
 
                       else
 
c                         ----------------------------------------------
c                         move c(rowpos,:) => c (kpos,:)
c                         move c(fflefc,:) => c (rowpos,:)
c                         c (fflefc,:) = 0
c                         ----------------------------------------------
 
                          flp = ffxp + (fflefc - 1)
cfpp$ nodepchk l
                          do 200 j = 0, (ldimr - 1) * ffdimc, ffdimc
                             xx (fdp+j) = xx (fsp+j)
                             xx (fsp+j) = xx (flp+j)
                             xx (flp+j) = zero
200                       continue
 
                          swprow = wpc (fflefc)
                          wpc (rowpos) = swprow
                          wir (swprow) = xsp
                       endif
 
                    else if (rowpos .ne. kpos) then
 
c                      -------------------------------------------------
c                      swap c (rowpos,:) <=> c (kpos,:)
c                      -------------------------------------------------
cfpp$ nodepchk l
                       do 210 j = 0, (ldimr - 1) * ffdimc, ffdimc
                          x = xx (fdp+j)
                          xx (fdp+j) = xx (fsp+j)
                          xx (fsp+j) = x
210                    continue
 
                       swprow = wpc (kpos)
                       wpc (rowpos) = swprow
                       wir (swprow) = xsp
                    endif
 
                    fflefc = fflefc - 1
                 endif
 
220           continue
 
c             ----------------------------------------------------------
c             merge with pivot row pattern of new lu arrowhead
c             ----------------------------------------------------------
 
              i = fflefr
              do 230 p = lurp, lurp + ludegr - 1
                 col = lui (p)
                 if (wic (col) .lt. 0) then
                    wic (col) = i * ffdimc
                    i = i + 1
                    wpr (i) = col
                 endif
230           continue
 
c             ----------------------------------------------------------
c             merge with pivot column pattern of new lu arrowhead
c             ----------------------------------------------------------
 
              i = fflefc
              do 240 p = lucp, lucp + ludegc - 1
                 row = lui (p)
                 if (wir (row) .lt. 0) then
                    wir (row) = i
                    i = i + 1
                    wpc (i) = row
                 endif
240           continue
 
           endif
 
c=======================================================================
c  done initializing frontal matrix ]
c=======================================================================
 
c=======================================================================
c  assemble original arrowheads into the frontal matrix, and deallocate
c=======================================================================
 
c          -------------------------------------------------------------
c          current workspace usage:
c          -------------------------------------------------------------
 
c          wpc (1..ludegr):     holds the pivot column pattern
c                               (excluding the pivot row indices)
c
c          wpr (1..ludegr):     holds the pivot row pattern
c                               (excluding the pivot column indices)
c
c          c (1..ffdimr, 1..ffdimc):  space for the frontal matrix,
c               in xx (ffxp ... ffxp + ffsize - 1)
c
c          c (i,j) is located at xx (ffxp+((i)-1)+((j)-1)*ffdimc)
c
c          c (1..ludegc, 1..ludegr):            contribution block
c          c (ludegc+1..ludegc+luk, 1..ludegr):             u2 block
c          c (1..ludegc, ludegr+1..ludegr+luk):             l2 block
c          c (ludegc+1..ludegc+luk, ludegr+1..ludegr+luk):  l1\u1 block
c
c          wir (row) >= 0 for each row in pivot column pattern.
c               offset into pattern is given by:
c               wir (row) == offset - 1
c               also, wir (npiv+1 ... npiv+luk) is
c               ludegc+luk-1 ... ludegc, the offsets of the pivot rows.
c
c               otherwise, wir (1..n) is < 0
c
c          wic (col) >= 0 for each col in pivot row pattern.
c               wic (col) == (offset - 1) * ffdimc
c               also, wic (npiv+1 ... npiv+luk) is
c               ludegr+luk-1 ... ludegr, the offsets of the pivot rows.
c
c               otherwise, wic (1..n) is < 0
 
           do 260 k = 1, luk
              i = npiv + k
              xcdp = ffxp + wic (i)
              xrdp = ffxp + wir (i)
              do 250 p = cp (i+1), cp (i) - 1
                 j = ari (p)
                 if (j .gt. 0) then
c                   a diagonal entry, or lower triangular entry
c                   row = j, col = i
                    xp = xcdp + wir (j)
                    if (xp .lt. xcdp) then
c                      invalid entry - not in prior lu pattern
                       noutsd = noutsd + 1
                       if (pr3) then
c                         get original row and column index and print it
                          row = rperm (j)
                          col = cperm (i)
                          call ums2p2 (2, 97, row, col, xx (p), io)
                       endif
                    else
                       xx (xp) = xx (xp) + xx (p)
                    endif
                 else
c                   an upper triangular entry
c                   row = i, col = -j
                    xp = xrdp + wic (-j)
                    if (xp .lt. xrdp) then
c                      invalid entry - not in prior lu pattern
                       noutsd = noutsd + 1
                       if (pr3) then
c                         get original row and column index and print it
                          row = rperm (i)
                          col = cperm (-j)
                          call ums2p2 (2, 97, row, col, xx (p), io)
                       endif
                    else
                       xx (xp) = xx (xp) + xx (p)
                    endif
                 endif
250           continue
260        continue
 
c          deallocate the original arrowheads
           p = cp (npiv + luk + 1)
           xs = cp (npiv + 1) - p
           frxp (mhead) = p
           xneed = xneed - xs
           if (xs .gt. xfree) then
              xfree = xs
              pfree = mhead
           endif
 
c=======================================================================
c  assemble lusons, usons, and lsons into the frontal matrix [
c=======================================================================
 
           do 480 sp = lusonp, lusonp + lunson - 1
 
c             ----------------------------------------------------------
c             get the son and determine its type (luson, uson, or lson)
c             ----------------------------------------------------------
 
              e = lui (sp)
              if (e .le. n) then
c                luson
                 type = 1
              else if (e .le. 2*n) then
c                uson
                 e = e - n
                 type = 2
              else
c                lson
                 e = e - 2*n
                 type = 3
              endif
 
c             ----------------------------------------------------------
c             if fdimc=0 this is the implicit luson (already assembled)
c             ----------------------------------------------------------
 
              fdimc = frdimc (e)
              if (fdimc .ne. 0) then
 
c                -------------------------------------------------------
c                get scalar info of the son (it needs assembling)
c                -------------------------------------------------------
 
                 fxp = frxp (e)
                 fluip = lup (e)
                 fdegr = lui (fluip+2)
                 fdegc = lui (fluip+3)
                 allcol = fdegr .gt. 0
                 allrow = fdegc .gt. 0
                 fdegr = abs (fdegr)
                 fdegc = abs (fdegc)
                 flucp = (fluip + 7)
                 flurp = flucp + fdegc
 
c                use wm (1..fdegc) for offsets:
 
c                -------------------------------------------------------
                 if (type .eq. 1) then
c                this is an luson - assemble an entire frontal matrix
c                -------------------------------------------------------
 
c                   ----------------------------------------------------
                    if (allrow) then
c                   no rows assembled out of this luson yet
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
                       do 270 i = 0, fdegc-1
                          row = lui (flucp+i)
                          wm (i+1) = wir (row)
270                    continue
 
c                      -------------------------------------------------
                       if (allcol) then
c                      no rows or cols assembled out of luson yet
c                      -------------------------------------------------
 
                          do 290 j = 0, fdegr-1
                             col = lui (flurp+j)
                             xdp = ffxp + wic (col)
cfpp$ nodepchk l
                             do 280 i = 0, fdegc-1
                                xx (xdp + wm (i+1)) =
     $                          xx (xdp + wm (i+1)) +
     $                          xx (fxp + j*fdimc + i)
280                          continue
290                       continue
 
c                      -------------------------------------------------
                       else
c                      some columns already assembled out of luson
c                      -------------------------------------------------
 
                          do 310 j = 0, fdegr-1
                             col = lui (flurp+j)
                             if (col .gt. 0) then
                                xdp = ffxp + wic (col)
cfpp$ nodepchk l
                                do 300 i = 0, fdegc-1
                                   xx (xdp + wm (i+1)) =
     $                             xx (xdp + wm (i+1)) +
     $                             xx (fxp + j*fdimc + i)
300                             continue
                             endif
310                       continue
 
                       endif
 
c                   ----------------------------------------------------
                    else
c                   some rows already assembled out of luson
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
                       degc = 0
                       do 320 i = 0, fdegc-1
                          row = lui (flucp+i)
                          if (row .gt. 0) then
                             degc = degc + 1
                             wj (degc) = i
                             wm (degc) = wir (row)
                          endif
320                    continue
 
c                      -------------------------------------------------
                       if (allcol) then
c                      some rows already assembled out of luson
c                      -------------------------------------------------
 
                          do 340 j = 0, fdegr-1
                             col = lui (flurp+j)
                             xdp = ffxp + wic (col)
cfpp$ nodepchk l
                             do 330 i = 1, degc
                                xx (xdp + wm (i)) =
     $                          xx (xdp + wm (i)) +
     $                          xx (fxp + j*fdimc + wj (i))
330                          continue
340                       continue
 
c                      -------------------------------------------------
                       else
c                      rows and columns already assembled out of luson
c                      -------------------------------------------------
 
                          do 360 j = 0, fdegr-1
                             col = lui (flurp+j)
                             if (col .gt. 0) then
                                xdp = ffxp + wic (col)
cfpp$ nodepchk l
                                do 350 i = 1, degc
                                   xx (xdp + wm (i)) =
     $                             xx (xdp + wm (i)) +
     $                             xx (fxp + j*fdimc + wj (i))
350                             continue
                             endif
360                       continue
 
                       endif
                    endif
 
c                   ----------------------------------------------------
c                   deallocate the luson frontal matrix
c                   ----------------------------------------------------
 
                    frdimc (e) = 0
                    prev = frprev (e)
                    next = frnext (e)
                    xneed = xneed - fdegr*fdegc
                    xruse = xruse - fdegr*fdegc
 
                    if (frdimc (prev) .le. 0) then
c                      previous block is free - delete this block
                       frnext (prev) = next
                       frprev (next) = prev
                       e = prev
                       prev = frprev (e)
                    endif
 
                    if (frdimc (next) .le. 0) then
c                      next block is free - delete this block
                       frxp (next) = frxp (e)
                       if (e .le. nlu) then
                          frnext (prev) = next
                          frprev (next) = prev
                       endif
                       e = next
                       next = frnext (e)
                       if (frnext (mhead) .eq. mtail) then
c                         no blocks left except mhead and mtail
                          frxp (mtail) = frxp (mhead)
                       endif
                    endif
 
c                   get the size of the freed block
                    if (next .eq. 0) then
c                      this is the mtail block
                       xs = ffxp - frxp (e)
                    else
                       xs = frxp (next) - frxp (e)
                    endif
                    if (xs .gt. xfree) then
c                      keep track of the largest free block
                       xfree = xs
                       pfree = e
                    endif
 
c                -------------------------------------------------------
                 else if (type .eq. 2) then
c                uson:  assemble all possible columns
c                -------------------------------------------------------
 
c                   ----------------------------------------------------
                    if (allrow) then
c                   no rows assembled out of this uson yet
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
                       do 370 i = 0, fdegc-1
                          row = lui (flucp+i)
                          wm (i+1) = wir (row)
370                    continue
 
                       do 390 j = 0, fdegr-1
                          col = lui (flurp+j)
                          if (col .gt. 0) then
                             if (wic (col) .ge. 0) then
                                xdp = ffxp + wic (col)
cfpp$ nodepchk l
                                do 380 i = 0, fdegc-1
                                   xx (xdp + wm (i+1)) =
     $                             xx (xdp + wm (i+1)) +
     $                             xx (fxp + j*fdimc + i)
380                             continue
c                               flag this column as assembled
                                lui (flurp+j) = -col
                             endif
                          endif
390                    continue
 
c                   ----------------------------------------------------
                    else
c                   some rows already assembled out of this uson
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
                       degc = 0
                       do 400 i = 0, fdegc-1
                          row = lui (flucp+i)
                          if (row .gt. 0) then
                             degc = degc + 1
                             wj (degc) = i
                             wm (degc) = wir (row)
                          endif
400                    continue
 
                       do 420 j = 0, fdegr-1
                          col = lui (flurp+j)
                          if (col .gt. 0) then
                             if (wic (col) .ge. 0) then
                                xdp = ffxp + wic (col)
cfpp$ nodepchk l
                                do 410 i = 1, degc
                                   xx (xdp + wm (i)) =
     $                             xx (xdp + wm (i)) +
     $                             xx (fxp + j*fdimc + wj (i))
410                             continue
c                               flag this column as assembled
                                lui (flurp+j) = -col
                             endif
                          endif
420                    continue
 
                    endif
 
c                   flag this element as missing some columns
                    lui (fluip+2) = -fdegr
 
c                -------------------------------------------------------
                 else
c                lson:  assemble all possible rows
c                -------------------------------------------------------
 
c                   compute the compressed column offset vector
                    degc = 0
                    do 430 i = 0, fdegc-1
                       row = lui (flucp+i)
                       if (row .gt. 0) then
                          if (wir (row) .ge. 0) then
c                            this row will be assembled in loop below
                             degc = degc + 1
                             wj (degc) = i
                             wm (degc) = wir (row)
c                            flag this row as assembled
                             lui (flucp+i) = -row
                          endif
                       endif
430                 continue
 
c                   ----------------------------------------------------
                    if (allcol) then
c                   no columns assembled out of this lson yet
c                   ----------------------------------------------------
 
                       do 450 j = 0, fdegr-1
                          col = lui (flurp+j)
                          xdp = ffxp + wic (col)
cfpp$ nodepchk l
                          do 440 i = 1, degc
                             xx (xdp + wm (i)) =
     $                       xx (xdp + wm (i)) +
     $                       xx (fxp + j*fdimc + wj (i))
440                       continue
450                    continue
 
c                   ----------------------------------------------------
                    else
c                   some columns already assembled out of this lson
c                   ----------------------------------------------------
 
                       do 470 j = 0, fdegr-1
                          col = lui (flurp+j)
                          if (col .gt. 0) then
                             xdp = ffxp + wic (col)
cfpp$ nodepchk l
                             do 460 i = 1, degc
                                xx (xdp + wm (i)) =
     $                          xx (xdp + wm (i)) +
     $                          xx (fxp + j*fdimc + wj (i))
460                          continue
                          endif
470                    continue
 
                    endif
 
c                   flag this element as missing some rows
                    lui (fluip+3) = -fdegc
 
                 endif
 
              endif
 
480        continue
 
c=======================================================================
c  done assemblying sons into the frontal matrix ]
c=======================================================================
 
c=======================================================================
c  factorize the frontal matrix [
c=======================================================================
 
           k0 = 0
           fflefr = ldimr
           fflefc = ldimc
           ffcp = ffxp + fflefr * ffdimc
           ffrp = ffxp + fflefc
           ffpp = ffxp + fflefc + fflefr * ffdimc
 
           do 500 k = 1, luk
 
c             ----------------------------------------------------------
c             compute kth column of u1, and update pivot column
c             ----------------------------------------------------------
 
              if (k-k0-2 .gt. 0) then
c                u1 = l1 \ u1.  note that l1 transpose is stored, and
c                that u1 is stored with rows in reverse order.
                 call strsv ('u', 'n', 'u', k-k0-1,
     $                         xx (ffpp         ), ffdimc,
     $                         xx (ffpp - ffdimc), 1)
                 rinfo (5) = rinfo (5) + (k-k0-2)*(k-k0-1)
              endif
              if (k-k0-1 .gt. 0) then
c                l1 = l1 - l2*u1
                 call sgemv ('n', fflefc, k-k0-1,
     $                   -one, xx (ffcp         ), ffdimc,
     $                         xx (ffpp - ffdimc), 1,
     $                    one, xx (ffcp - ffdimc), 1)
                 rinfo (5) = rinfo (5) + 2*fflefc*(k-k0-1)
              endif
 
              ffcp = ffcp - ffdimc
              ffrp = ffrp - 1
              ffpp = ffpp - ffdimc - 1
              fflefr = fflefr - 1
              fflefc = fflefc - 1
 
c             ----------------------------------------------------------
c             divide pivot column by pivot
c             ----------------------------------------------------------
 
c             k-th pivot in frontal matrix located in xx (ffpp)
              x = xx (ffpp)
              if (abs (x) .eq. zero) then
c                error return, if pivot order from ums2fa not acceptable
                 go to 9010
              endif
              x = one / x
              do 490 p = ffcp, ffcp + fflefc - 1
                 xx (p) = xx (p) * x
490           continue
c             count this as a call to the level-1 blas:
              rinfo (4) = rinfo (4) + fflefc
              info (17) = info (17) + 1
 
c             ----------------------------------------------------------
c             compute u1 (k0+1..k, k..ldimc) and
c             update contribution block: rank-nb, or if last pivot
c             ----------------------------------------------------------
 
              if (k-k0 .ge. nb .or. k .eq. luk) then
                 call strsm ('l', 'u', 'n', 'u', k-k0, fflefr, one,
     $                      xx (ffpp), ffdimc,
     $                      xx (ffrp), ffdimc)
                 call sgemm ('n', 'n', fflefc, fflefr, k-k0,
     $                -one, xx (ffcp ), ffdimc,
     $                      xx (ffrp ), ffdimc,
     $                 one, xx (ffxp), ffdimc)
                 rinfo (6) = rinfo (6) + fflefr*(k-k0-1)*(k-k0)
     $                                         + 2*fflefc*fflefr*(k-k0)
                 k0 = k
              endif
 
500        continue
 
c=======================================================================
c  done factorizing the frontal matrix ]
c=======================================================================
 
c=======================================================================
c  save the new lu arrowhead [
c=======================================================================
 
c          allocate permanent space for the lu arrowhead
           xs = luk*ludegc + luk*ludegr + luk*luk
 
           if (xs .gt. xtail-xhead) then
              info (15) = info (15) + 1
              call ums2rg (xx, xsize, xhead, xtail, xuse,
     $              lui, frdimc, frxp, frnext, frprev, nlu, lup,
     $              icntl, ffxp, ffsize, pfree, xfree)
           endif
 
           xtail = xtail - xs
           luxp = xtail
           xuse = xuse + xs
           xneed = xneed + xs
           xruse = xruse + xs
           xrmax = max (xrmax, xruse)
           info (20) = max (info (20), xuse)
           info (21) = max (info (21), xneed)
           if (xhead .gt. xtail) then
c             error return, if not enough real memory:
              go to 9000
           endif
 
c          save the scalar data of the lu arrowhead
           lui (luip) = luxp
 
c          save column pattern (it may have been rearranged)
           do 510 i = 0, ludegc-1
              lui (lucp+i) = wpc (i+1)
510        continue
 
c          save row pattern (it may have been rearranged)
           do 520 i = 0, ludegr-1
              lui (lurp+i) = wpr (i+1)
520        continue
 
c          move the l1,u1 matrix, compressing the dimension from
c          ffdimc to ldimc.  the lu arrowhead grows on top of stack.
           xp = ffxp + (ldimr-1)*ffdimc + ldimc-1
           do 540 j = 0, luk-1
cfpp$ nodepchk l
              do 530 i = 0, luk-1
                 xx (luxp + j*ldimc + i) = xx (xp - j*ffdimc - i)
530           continue
540        continue
 
c          move l2 matrix, compressing dimension from ffdimc to ldimc
           if (ludegc .ne. 0) then
              lxp = luxp + luk
              xp = ffxp + (ldimr-1)*ffdimc
              do 560 j = 0, luk-1
cfpp$ nodepchk l
                 do 550 i = 0, ludegc-1
                    xx (lxp + j*ldimc + i) = xx (xp - j*ffdimc + i)
550              continue
560           continue
           endif
 
c          move the u2 block.
           if (ludegr .ne. 0) then
              uxp = luxp + luk * ldimc
              xp = ffxp + ldimc-1
              do 580 j = 0, ludegr-1
cfpp$ nodepchk l
                 do 570 i = 0, luk-1
                    xx (uxp + j*luk + i) = xx (xp + j*ffdimc - i)
570              continue
580           continue
           endif
 
c          one more lu arrowhead has been refactorized
           nzu = (luk*(luk-1)/2) + luk*ludegc
           nzl = (luk*(luk-1)/2) + luk*ludegr
           info (10) = info (10) + nzl
           info (11) = info (11) + nzu
 
c          -------------------------------------------------------------
c          clear the pivot row and column offsets
c          -------------------------------------------------------------
 
           do 590 pivot = npiv + 1, npiv + luk
              wir (pivot) = -1
              wic (pivot) = -1
590        continue
           npiv = npiv + luk
 
c=======================================================================
c  done saving the new lu arrowhead ]
c=======================================================================
 
600     continue
 
c=======================================================================
c  factorization complete ]
c=======================================================================
 
c=======================================================================
c  wrap-up:  store lu factors in their final form
c=======================================================================
 
c       ----------------------------------------------------------------
c       flag remaining arrowheads as invalid entries, if prior matrix
c       was singular.  print them if requested.
c       ----------------------------------------------------------------
 
        if (npiv .lt. n) then
           if (pr3) then
              do 620 i = npiv+1, n
                 do 610 p = cp (i+1), cp (i) - 1
                    j = ari (p)
                    if (j .gt. 0) then
c                      a diagonal entry, or lower triangular entry
c                      get original row and column index
                       row = rperm (j)
                       col = cperm (i)
                    else
c                      an upper triangular entry
c                      get original row and column index
                       row = rperm (i)
                       col = cperm (-j)
                    endif
                    call ums2p2 (2, 95, row, col, xx(p), io)
610              continue
620           continue
           endif
           noutsd = noutsd + (cp (npiv+1) - cp (n+1))
        endif
 
c       ----------------------------------------------------------------
c       deallocate all remaining input arrowheads and frontal matrices
c       ----------------------------------------------------------------
 
        if (ffsize .ne. 0) then
           info (13) = info (13) + 1
        endif
        xuse = xuse - (xhead - cp (n+1))
        xneed = xuse
        xhead = cp (n+1)
 
        if (nlu .eq. 0) then
c          lu factors are completely empty (a = 0).
c          add one real, to simplify rest of code.
c          otherwise, some arrays in ums2rf or ums2so would have
c          zero size, which can cause an address fault.
           xtail = xsize
           xuse = xuse + 1
           xruse = xuse
           xneed = xuse
           info (20) = max (info (20), xuse)
           info (21) = max (info (21), xneed)
        endif
 
        if (xhead .le. xtail) then
 
c          -------------------------------------------------------------
c          sufficient memory to complete the factorization
c          -------------------------------------------------------------
 
           if (nlu .eq. 0) then
c             zero the dummy entry, although it won't be accessed:
              xx (xtail) = zero
           endif
 
c          -------------------------------------------------------------
c          update pointers in lu factors
c          -------------------------------------------------------------
 
           do 630 s = 1, nlu
              luip = lup (s)
              luxp = lui (luip)
              lui (luip) = luxp - xtail + 1
630        continue
 
c          -------------------------------------------------------------
c          get memory usage estimate for next call to ums2rf
c          -------------------------------------------------------------
 
           xruse = xuse
           xrmax = max (xrmax, xruse)
           return
 
        endif
 
c=======================================================================
c  error conditions
c=======================================================================
 
c       error return label:
9000    continue
c       out of real memory
        call ums2er (2, icntl, info, -4, info (21))
        return
 
c       error return label:
9010    continue
c       original pivot order computed by ums2fa is no longer acceptable
        call ums2er (2, icntl, info, -6, 0)
        return
        end
