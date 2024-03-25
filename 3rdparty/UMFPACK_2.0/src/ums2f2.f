 
        subroutine ums2f2 (cp, nz, n, pn, cperm, rperm, itail, xtail,
     $          xx, xsize, ii, isize, icntl, cntl, info, rinfo, pgiven,
     $          iuse, xuse, wir, wic, wpr, wpc, wm, head,
     $          wj, rp, wc, wr, dn, dsiz, keep,
     $          rmax, cmax, totnlu, xrmax, xruse)
c
cc UMS2F2 is a utility routine which factors part of the matrix.
c
        integer xsize, isize, icntl (20), info (40), pn,
     $          itail, xtail, nz, n, ii (isize), cp (n+1), dn, dsiz,
     $          rperm (pn), cperm (pn), wir (n), wic (n), wpr (n),
     $          wpc (n), wm (n), head (n), rp (n+dn), wc (n+dn),
     $          wr (n+dn), iuse, xuse, wj (n), keep (20),
     $          rmax, cmax, totnlu, xrmax, xruse
        logical pgiven
        real
     $          xx (xsize), cntl (10), rinfo (20)
 
c=== ums2f2 ============================================================
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
c  ums2f2 factorizes the n-by-n input matrix at the head of ii/xx
c  (in expanded column-form) and places its lu factors at the tail of
c  ii/xx.  the input matrix is overwritten.   no btf information is
c  used in this routine, except that the btf permutation arrays are
c  modified to include the final permutations.
 
c=======================================================================
c  input:
c=======================================================================
c
c       cp (1..n+1):    column pointers of expanded column-form,
c                       undefined on output
c       n:              order of input matrix
c       nz:             entries in input matrix
c       isize:          size of ii
c       xsize:          size of xx
c       iuse:           memory usage in index
c       xuse:           memory usage in value
c       icntl:          integer control parameters, see ums2in
c       cntl:           real control parameters, see ums2in
c       keep (6)        integer control parameter, see ums2in
c       dn:             number of dense columns
c       dsiz:           entries required for col to be treated as dense
c       rmax:           maximum ludegr seen so far (see below)
c       cmax:           maximum ludegc seen so far (see below)
c       totnlu:         total number of lu arrowheads constructed so far
c       xrmax:          maximum real memory usage for ums2rf
c       xruse:          current real memory usage for ums2rf
c
c       pgiven:         true if cperm and rperm are defined on input
c       if pgiven then:
c          cperm (1..pn):       col permutation to btf, n = pn
c          rperm (1..pn):       row permutation to btf
c       else
c          cperm (1..pn):       unaccessed pn = 1
c          rperm (1..pn):       unaccessed
c
c       ii (1..nz+cscal*n):             expanded column-form, see below
c       ii (nz+cscal*n+1..isize):       undefined on input
c       xx (1..nz):                     expanded column-form, see below
c       xx (nz+1..xsize):               undefined on input
 
c=======================================================================
c  workspace:
c=======================================================================
c
c       wir (1..n)
c       wic (1..n)
c       wpr (1..n)
c       wpc (1..n)
c       wm (1..n)
c       head (n)
c       rp (1..n+dn)
c       wr (1..n+dn)
c       wc (1..n+dn)
 
c=======================================================================
c  output:
c=======================================================================
c
c       ii (1..itail-1):        undefined on output
c       ii (itail..isize):      lu factors of this matrix, see below
c       xx (1..xtail-1):        undefined on output
c       xx (xtail..xsize):      lu factors of this matrix, see below
c
c       info:           integer informational output, see ums2fa
c       rinfo:          real informational output, see ums2fa
c       if pgiven:
c          cperm (1..n): the final col permutations, including btf
c          rperm (1..n): the final row permutations, including btf
c
c       wic (1..n):     row permutations, not including btf
c       wir (1..n):     column permutations, not including btf
c
c       iuse:           memory usage in index
c       xuse:           memory usage in value
c       rmax:           maximum ludegr seen so far (see below)
c       cmax:           maximum ludegc seen so far (see below)
c       totnlu:         total number of lu arrowheads constructed so far
c       xrmax:          maximum real memory usage for ums2rf
c       xruse:          current real memory usage for ums2rf
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2f1
c       subroutines called:     ums2er, ums2fg, sgemv, sgemm
c       functions called:       isamax, abs, max, min
        integer isamax
        intrinsic abs, max, min
 
c=======================================================================
c  description of data structures:
c=======================================================================
 
c-----------------------------------------------------------------------
c  column/element/arrowhead pointers:
c-----------------------------------------------------------------------
c
c  the cp (1..n) array contains information about non-pivotal columns
c
c       p = cp (j)
c       if (p = 0) then j is a pivotal column
c       else i is a non-pivotal column
c
c  the rp (1..n) array contains information about non-pivotal rows,
c  unassembled frontal matrices (elements), and the lu arrowheads
c
c       p = rp (i)
c       if (i > n) then
c          i is an artificial frontal matrix (a dense column)
c          if (p = 0) then i is assembled, else unassembled
c       else if (p = 0) then i is pivotal but not element/arrowhead
c       else if (wc (i) >= 0 and wc (i) <= n) then
c          i is a non-pivotal row
c       else if (wc (i) = -(n+dn+2)) then
c          i is a pivotal row, an assembled element, and an lu arrowhead
c       else i an unassembled element
 
c-----------------------------------------------------------------------
c  matrix being factorized:
c-----------------------------------------------------------------------
c
c    each column is stored in ii and xx:
c    -----------------------------------
c
c       if j is a non-pivotal column, pc = cp (j):
c
c       csiz = ii (pc) size of the integer data structure for col j,
c                        including the cscal scalars
c       cdeg = ii (pc+1) degree of column j
c       cxp  = ii (pc+2) pointer into xx for numerical values
c       next = ii (pc+3) pointer to next block of memory in xx
c       prev = ii (pc+4) pointer to previous block of memory in xx
c       celn = ii (pc+5) number of elements in column j element list
c       clen = ii (pc+6) number of original entries in column j
c       cnxt = ii (pc+7) next column with same degree as col j
c       cprv = ii (pc+8) previous column with same degree as col j
c       cep = (pc+9) pointer to start of the element list
c       ii (cep ... cep + 2*celn - 1)
c                       element list (e,f) for the column
c       ii (cep + 2*celn ... pc + csiz - clen - 1)
c                       empty
c       ii (pc + csiz - clen ... pc + csiz - 1)
c                       row indices of original nonzeros in the column
c       xx (xp ... xp + clen - 1)
c                       numerical values of original nonzeros in the col
c
c       if cdeg = ii (pc+1) = -(n+2), then this is a singular column
c       if cdeg = -1, then this column is deallocated
c
c    each row is stored in ii only:
c    ------------------------------
c
c       if i is a non-pivotal row, pr = rp (i)
c
c       rsiz = ii (pr) size of the integer data structure for row i,
c                        including the rscal scalars
c       rdeg = ii (pr+1) degree of row i
c       reln = wr (i) number of elements in row i element list
c       rlen = wc (i) number of original entries in row i
c       rep  = (pr+2) pointer to start of the element list
c       ii (rep ... rep + 2*reln - 1)
c                       element list (e,f) for the row
c       ii (rep + 2*reln ... pr + rsiz - rlen - 1)
c                       empty
c       ii (pr + rsiz - rlen ... pr + rsiz - 1)
c                       column indices of original nonzeros in the row
c
c       if rdeg = -1, then this row is deallocated
 
c-----------------------------------------------------------------------
c  frontal matrices
c-----------------------------------------------------------------------
c
c   each unassembled frontal matrix (element) is stored as follows:
c       total size: fscal integers, (fdimr*fdimc) reals
c
c       if e is an unassembled element, ep = rp (e), and e is also
c       the first pivot row in the frontal matrix.
c
c       fluip  = ii (ep)        pointer to lu arrowhead in ii
c       fdimc  = ii (ep+1)      column dimension of contribution block
c       fxp    = ii (ep+2)      pointer to contribution block in xx
c       next   = ii (ep+3)      pointer to next block in xx
c       prev   = ii (ep+4)      pointer to previous block in xx
c       fleftr = ii (ep+5)      number of unassembled rows
c       fleftc = ii (ep+6)      number of unassembled columns
c       fextr = wr (e) - w0     external row degree of the frontal mtx
c       fextc = wc (e) - w0     external col degree of the frontal mtx
c       xx (fxp ... )
c               a 2-dimensional array, c (1..fdimc, 1..fdimr).
c               note that fdimr is not kept (it is not needed,
c               except for the current frontal).  if this is not the
c               current frontal matrix, then luip points to the
c               corresponding lu arrowhead, and the contribution block
c               is stored in c (1..ludegc, 1..ludegr) in the
c               c (1..fdimc, ...) array.
c
c               if memory is limited, garbage collection will occur.
c               in this case, the c (1..fdimc, 1..fdimr) array is
c               compressed to be just large enough to hold the
c               unassembled contribution block,
c               c (1..ludegc, 1..ludegr).
 
c-----------------------------------------------------------------------
c  artificial frontal matrices
c-----------------------------------------------------------------------
c
c   an artificial frontal matrix is an original column that is treated
c   as a c-by-1 frontal matrix, where c is the number of original
c   nonzeros in the column.  dense columns (c > dsiz) are treated this
c   way.  an artificial frontal matrix is just the same as a frontal
c   matrix created by the elimination of one or more pivots, except
c   that there is no corresponding lu arrowhead.  the row and column
c   patterns are stored in:
c
c       ep = rp (e), where e = n+1 .. n+dn, where there are dn
c                    artificial frontal matrices.
c
c       lucp = (ep+9)   pointer to row pattern (just one column index)
c       lurp = (ep+8) pointer to column pattern (fdimc row indices)
 
c-----------------------------------------------------------------------
c  current frontal matrix
c-----------------------------------------------------------------------
c
c  ffxp points to current frontal matrix (contribution block and lu
c  factors).  for example, if fflefc = 4, fflefr = 6, k = 3, and
c  gro = 2.0, then "x" is a term in the contribution block, "l" in l1,
c  "u" in u1, "l" in l2, "u" in u2, and "." is unused.  xx (fxp) is "x".
c  the first 3 pivot values (diagonal entries in u1) are 1,2, and 3.
c  for this frontal matrix, ffdimr = 12 (the number of columns), and
c  ffdimc = 8 (the number of rows).  the frontal matrix is
c  ffdimc-by-ffdimr
c
c                             |----------- col 1 of l1 and l2, etc.
c                             v
c       x x x x x x . . . l l l
c       x x x x x x . . . l l l
c       x x x x x x . . . l l l
c       x x x x x x . . . l l l
c       . . . . . . . . . . . .
c       u u u u u u . . . 3 l l         <- row 3 of u1 and u2
c       u u u u u u . . . u 2 l         <- row 2 of u1 and u2
c       u u u u u u . . . u u 1         <- row 1 of u1 and u2
 
c-----------------------------------------------------------------------
c  lu factors
c-----------------------------------------------------------------------
c
c   the lu factors are placed at the tail of ii and xx.  if this routine
c   is factorizing a single block, then this decription is for the
c   factors of the single block:
c
c       ii (itail):             xtail = start of lu factors in xx
c       ii (itail+1):           nlu = number of lu arrowheads
c       ii (itail+2):           npiv = number of pivots
c       ii (itail+3):           maximum number of rows in any
c                               contribution block (max ludegc)
c       ii (itail+4):           maximum number of columns in any
c                               contribution block (max ludegr)
c       ii (itail+5..itail+nlu+4): lup (1..nlu) array, pointers to each
c                               lu arrowhead, in order of their
c                               factorization
c       ii (itail+nlu+5...isize):integer info. for lu factors
c       xx (xtail..xsize):      real values in lu factors
c
c   each lu arrowhead is stored as follows:
c   ---------------------------------------
c
c       total size: (7 + ludegc + ludegr + nsons) integers,
c                   (luk**2 + ludegc*luk + luk*ludegc) reals
c
c       if e is an lu arrowhead, then luip = rp (e), and luip >= itail.
c       when ums2f2 returns, then luip is given by luip =
c       ii (itail+s+1), where s = 1..nlu is the position of the lu
c       arrowhead in the lu factors (s=1,2,.. refers to the first,
c       second,.. lu arrowhead)
c
c       luxp   = ii (luip) pointer to numerical lu arrowhead
c       luk    = ii (luip+1) number of pivots in lu arrowhead
c       ludegr = ii (luip+2) degree of last row of u (excl. diag)
c       ludegc = ii (luip+3) degree of last col of l (excl. diag)
c       nsons  = ii (luip+4) number of children in assembly dag
c       ludimr = ii (luip+5)
c       ludimc = ii (luip+5) max front size is ludimr-by-ludimc,
c                       or zero if this lu arrowhead factorized within
c                       the frontal matrix of a prior lu arrowhead.
c       lucp   = (luip + 7)
c                       pointer to pattern of column of l
c       lurp   = lucp + ludegc
c                       pointer to patter of row of u
c       lusonp = lurp + ludegr
c                       pointer to list of sons in the assembly dag
c       ii (lucp ... lucp + ludegc - 1)
c                       row indices of column of l
c       ii (lurp ... lurp + ludegr - 1)
c                       column indices of row of u
c       ii (lusonp ... lusonp + nsons - 1)
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
c       matrix.  after factorization, the negative flags are removed,
c       and the row/col indices are replaced with their corresponding
c       index in the permuted lu factors.
c
c   list of sons:
c       1 <= son <= n:           son an luson
c       n+1 <= son <= 2n:        son-n is an uson
c       2n+n <= son <= 3n:       son-2n is a lson
c       during factorzation, a son is referred to by its first
c       pivot column.  after factorization, they are numbered according
c       to their order in the lu factors.
 
c-----------------------------------------------------------------------
c  workspaces:
c-----------------------------------------------------------------------
c
c   wir (e):  link list of sons of the current element
c       wir (e) = -1 means that e is not in the list.
c       wir (e) = next+n+2 means that "next" is the element after e.
c       the end of the list is marked with wir (e) = -(n+2).
c       sonlst points to the first element in the list, or 0 if
c       the sonlst is empty.
c
c   wir (row), wic (col):  used for pivot row/col offsets:
c
c       if wir (row) >= 0 then the row is in the current
c       column pattern.  similarly for wic (col).
c
c       if wir (row) is set to "empty" (<= -1), then
c       the row is not in the current pivot column pattern.
c
c       similarly, if wic (col) is set to -2, then the column is
c       not in the current pivot row pattern.
c
c       if wic (col) = -1 then col is pivotal
c
c       after factorization, wir/c holds the pivot permutations.
c
c   wpr/c (1..n):  the first part is used for the current frontal
c           matrix pattern.  during factorization, the last part holds
c           a stack of the row and column permutations (wpr/c (n-k+1)
c           is the k-th pivot row/column).
c
c   head (1..n):        degree lists for columns.  head (d) is the
c                       first column in list d with degree d.
c                       the cnxt and cprv pointers are stored in the
c                       column data structure itself.
c                       mindeg is the least non-empty list
c
c   wm (1..n):          various uses
c   wj (1..degc) or wj (1..fdegc):      offset in pattern of a son
 
c-----------------------------------------------------------------------
c  memory allocation in ii and xx:
c-----------------------------------------------------------------------
c
c   ii (1..ihead):      rows and columns of active submatrix, and
c                       integer information for frontal matrices.
c   xx (1..xhead):      values of original entries in columns of
c                       matrix, values of contribution blocks, followed
c                       by the current frontal matrix.
c
c   mhead:              a pointer to the first block in the head of
c                       xx.  each block (a column or frontal matrix)
c                       contains a next and prev pointer for this list.
c                       if the list is traversed starting at mhead,
c                       then the pointers to the reals (cxp or fxp)
c                       will appear in strictly increasing order.
c                       note that the next, prev, and real pointers
c                       are in ii.  next and prev point to the next
c                       and previous block in ii, and the real pointer
c                       points to the real part in xx.
c
c   mtail:              the end of the memory list.
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer swpcol, swprow, fdimc, k0, colpos, rowpos, row2, rdeg2,
     $          p, i, j, ffrow, pivrow, pivcol, ludegr, ludegc, e1,
     $          fxp, lurp, lucp, ip, next, fflefr, pc, mnext, mprev,
     $          fflefc, fedegr, fedegc, k, xudp, xdp, xsp, xlp, s, col2,
     $          bestco, col, e, row, cost, srched, pr, f1, rscan, rep,
     $          kleft1, ffsize, ffxp, w0, ffdimr, ffdimc, kleft, xldp,
     $          ep, scan1, scan2, scan3, scan4, nzl, nzu, degc, cep
        integer mindeg, nsrch, npiv, eson, luip1, dnz, iworst, wxp,
     $          nb, lupp, nlu, nsons, ineed, xneed, ldimc, lxp, rlen2,
     $          rsiz, lsons, sonlst, xhead, ihead, deln, dlen,
     $          slist, xp, luip, rdeg, cdeg1, pfree, xfree, cdeg2,
     $          f, cdeg, mtail, mhead, rsiz2, csiz2, ip2, maxdr, maxdc,
     $          xs, is, luxp, fsp, flp, fdp, jj, usons, ndn, p2,
     $          csiz, celn, clen, reln, rlen, uxp, pc2, pr2
        integer cnxt, cprv, cxp, fluip, lusonp, fleftr, fleftc, maxint,
     $          fmaxr, fmaxc, slots, limit, rscal, cscal, fscal, extra,
     $          fmax, w0big, minmem, dummy1, dummy2, dummy3, dummy4
        logical symsrc, pfound, movelu, okcol, okrow, better
        real
     $          toler, maxval, relpt, gro, one, zero, x
        parameter (one = 1.0, zero = 0.0,
     $          rscal = 2, cscal = 9, fscal = 7,
     $          minmem = 24)
 
c  current element and working array, c:
c  -------------------------------------
c  ffxp:    current working array is in xx (ffxp ... ffxp+ffsize-1)
c  ffsize:  size of current working array in xx
c  ffdimr:  row degree (number of columns) of current working array
c  ffdimc:  column degree (number of rows) of current working array
c  fflefr:  row degree (number of columns) of current contribution block
c  fflefc:  column degree (number of rows) of current contribution block
c  fmaxr:   max row degree (maximum front size is fmaxr-by-fmaxc)
c  fmaxc:   max col degree (maximum front size is fmaxr-by-fmaxc)
c  fedegr:  extended row degree
c  fedegc:  extended column degree
c  ffrow:   current element being factorized (a pivot row index)
c  pivrow:  current pivot row index
c  pivcol:  current pivot column index
c  e1:      first pivot row in the frontal matrix
c  gro:     frontal matrix amalgamation growth factor
c  usons:   pointer to a link list of usons, in wc, assembled this scan3
c  lsons:   pointer to a link list of lsons, in wr, assembled this scan4
c  sonlst:  pointer to a link list of sons, in wir, of current element
c  swpcol:  the non-pivotal column to be swapped with pivot column
c  swprow:  the non-pivotal row to be swapped with pivot row
c  colpos:  position in wpr of the pivot column
c  rowpos:  position in wpc of the pivot row
c  k:       current pivot is kth pivot of current element
c  k0:      contribution block, c, has been updated with pivots 1..k0
c
c  lu arrowhead (a factorized element):
c  ------------------------------------
c  movelu:  true if a new lu arrowhead is to be created
c  luip:    current element is in ii (luip ...)
c  luip1:   first element from current frontal matrix in ii (luip1...)
c  ludegc:  degree of pivot column (excluding pivots themselves)
c  ludegr:  degree of pivot row (excluding pivots themselves)
c  lucp:    pattern of col(s) of current element in ii (lucp...)
c  lurp:    pattern of row(s) of current element in ii (lurp...)
c  lusonp:  list of sons of current element is in ii (lusonp...)
c  nsons:   number of sons of current element
c  ldimc:   column dimension (number of rows) of [l1\u1 l2] block
c  luxp:    numerical values of lu arrowhead stored in xx (luxp ...)
c  lxp:     l2 block is stored in xx (lxp ...) when computed
c  uxp:     u2 block is stored in xx (uxp ...) when computed
c  nzu:     nonzeros above diagonal in u in current lu arrowhead
c  nzl:     nonzeros below diagonal in l in current lu arrowhead
c
c  son, or element other than current element:
c  -------------------------------------------
c  e:       an element
c  eson:    an element
c  s:       a renumbered element (1..nlu) for ums2so and ums2rf
c  ep:      frontal matrix integer data struct. in ii (ep...ep+fscal-1)
c  fscal:   = 7, size of frontal matrix data structure
c  fluip:   lu arrowhead of e is in ii (fluip ...)
c  fxp:     contribution block of son is in xx (fxp ...)
c  fdimc:   leading dimension of contribution block of e
c  lucp:    pattern of col(s) of e in ii (lucp...)
c  lurp:    pattern of row(s) of e in ii (lurp...)
c  ludegr:  row degree of contribution block of e
c  ludegr:  column degree of contribution block of e
c  maxdr:   maximum ludegr for any lu arrowhead, for ums2rf
c  maxdc:   maximum ludegc for any lu arrowhead, for ums2rf
c  degc:    compressed column offset vector of son is in wj/wm (1..degc)
c  fleftr:  remaining row degree (number of columns) of a contrib. block
c  fleftc:  remaining column degree (number of rows) of a contrib. block
c  xudp:    pointer to a column of a prior contribution block
c  xldp:    pointer to a row of a prior contribution block
c
c  memory allocation:
c  ------------------
c  mhead:   head pointer for link list of blocks in xx
c  mtail:   tail pointer for link list of blocks in xx
c  mprev:   previous block, ii (p+4), of the block located at p
c  mnext:   next block, ii (p+3), of the block located at p
c  pfree:   ii (pfree+2) is the largest known free block in xx
c  xfree:   size of largest known free block in xx
c  xhead:   xx (1..xhead-1) is in use, xx (xhead ..xtail-1) is free
c  xtail:   xx (xtail..xsize) is in use, xx (xhead ..xtail-1) is free
c  xneed:   bare minimum memory currently needed in xx
c  ihead:   ii (1..ihead-1) is in use, ii (ihead ..itail-1) is free
c  itail:   ii (itail..isize) is in use, ii (ihead ..itail-1) is free
c  ineed:   bare minimum memory currently needed in ii
c  iworst:  worst possible current integer memory required
c  xs:      size of a block of memory in xx
c  is:      size of a block of memory in ii
c  wxp:     pointer to a temporary workspace in xx (wxp ... )
c  slots:   number of slots added to element lists during garbage coll.
c  minmem:  smallest isize allowed
c
c  wr and wc flag arrays:
c  ----------------------
c  w0:      marker value for wr (1..n) and wc (1..n) arrays
c  w0big:   largest permissible value of w0 (w0+n must not overflow)
c  fmax:    largest row/col degree of an element seen so far
c
c  a column:
c  ---------
c  pc:      pointer to a column, in ii (pc...)
c  pc2:     pointer to a column, in ii (pc2...)
c  csiz:    size of integer data structure of a column
c  csiz2:   size of integer data structure of a column
c  cscal:   = 9, number of scalars in data structure of a column
c  cdeg:    degree of a column
c  cdeg1:   degree of a column
c  cdeg2:   degree of a column
c  celn:    number of elements in the element list of a column
c  clen:    number of original entries that remain in a column
c  cnxt:    next column with same degree as this column
c  cprv:    previous column with same degree as this column
c  cep:     pointer to the element list of a column
c  cxp:     pointer to the numerical values in a column
c  limit:   maximum size for row/col data structure (excl. scalars)
c
c  dense columns:
c  --------------
c  dnz:     number of original entries that reside in "dense" columns
c  dn:      number of "dense" columns
c  ndn:     n + dn
c  extra:   number of extra slots to add to reconstructed "dense" cols
c
c  a row:
c  ------
c  pr:      pointer to a row, in ii (pr...)
c  pr2:     pointer to a row, in ii (pr2...)
c  rsiz:    size of integer data structure of a row
c  rsiz2:   size of integer data structure of a row
c  rscal:   = 2, number of scalars in data structure of a row
c  rdeg:    degree of a row
c  rdeg2:   degree of a row
c  reln:    number of elements in the element list of a row
c  rlen:    number of original entries that remain in a row
c  rlen2:   number of original entries that remain in a row
c  rep:     pointer to the element list of a row
c
c  pivot search:
c  -------------
c  cost:    approximate markowitz-cost of the current candidate pivot
c  bestco:  best approximate markowitz-cost seen so far
c  srched:  number of non-singular candidates searched so far
c  mindeg:  minimum degree of columns in active submatrix
c  nsrch:   maximum number of columns to search
c  slist:   pointer to a link list of searched columns, in ii
c  symsrc:  true if attempting to preserve symmetry
c  pfound:  true if pivot found during local search
c  okcol:   true if candidate pivot column is acceptable, so far
c  okrow:   true if candidate pivot row is acceptable, so far
c  toler:   pivot tolerance; abs(pivot) must be >= toler
c  maxval:  maximum absolute value in a candidate pivot column
c  relpt:   relative pivot tolerance (cntl (1))
c  npiv:    number of pivots factorized so far, incl. current element
c  kleft:   number of rows/columns remaining in active submatrix
c  kleft1:  kleft - 1
c  better:  if true, then candidate is better than the prior candidate
c
c  assembly:
c  ---------
c  f1:      degree prior to assembly next item
c  f:       offset into an element
c  rscan:   skip row assembly if more than rscan original entries
c  scan1:   start scan1 at wpc (scan1 ... fflefc)
c  scan2:   start scan2 at wpr (scan2 ... fflefr)
c  scan3:   start scan3 at wpr (scan3 ... fflefr)
c  scan4:   start scan4 at wpc (scan4 ... fflefc)
c  deln:    number of (e,f) tuples to delete from an element list
c  dlen:    number of original entries to delete from a row/col
c
c  allocated arrays:
c  -----------------
c  lupp:    lup (1..nlu) array located in ii (lupp...lupp+nlu-1)
c  nlu:     number of lu arrowheads
c
c  other:
c  ------
c  xdp:     destination pointer, into xx
c  xsp:     source pointer, into xx
c  xlp:     pointer into xx of location of last row/col in c
c  xp:      pointer into xx
c  ip:      pointer into ii
c  ip2:     pointer into ii
c  p2:      pointer into ii
c  fsp:     source pointer, into xx
c  fsp:     destination pointer, into xx
c  flp:     last row/column in current contribution is in xx (flp...)
c  col,col2: a column index
c  row,row2: a row index
c  nb:      block size for tradeoff between level-2 and level-3 blas
c  p, i, j, k, x:  various uses
c  jj:      loop index
c  maxint:  largest representable positive integer
c  next:    next pointer, for a link list
c  dummy1:  dummy loop index for main factorization loop
c  dummy2:  dummy loop index for global pivot search loop
c  dummy3:  dummy loop index for outer frontal matrix factorization loop
c  dummy4:  dummy loop index for inner frontal matrix factorization loop
 
c=======================================================================
c  executable statements:
c=======================================================================
 
c       ----------------------------------------------------------------
c       get control parameters and initialize various scalars
c       ----------------------------------------------------------------
 
        nsrch = max (1, icntl (5))
        symsrc = icntl (6) .ne. 0
        nb = max (1, icntl (7))
        relpt = max (zero, min (cntl (1), one))
        gro = max (one, cntl (2))
        maxint = keep (6)
        ndn = n + dn
        w0big = maxint - n
        w0 = ndn + 2
c       currently: w0 = n+dn+2 < 2n+2 < w0big = maxint - n
c       2n+2 < maxint - n must hold, so n < (maxint - 2) / 3 is the
c       largest that n can be.  this condition is checked in ums2fa.
        kleft = n
        npiv = 0
        nlu = 0
        mindeg = 1
        fmax = 1
        ihead = nz + cscal*n + 1
        xhead = nz + 1
        itail = isize + 1
        xtail = xsize + 1
c       cp (1) must equal 1, the first block
        xfree = -1
        pfree = 0
c       make sure integer space is at least of size minmem (simplifies
c       link list management and memory management)
        info (19) = max (info (19), iuse+minmem)
        if (ihead.gt.itail.or.isize.lt.minmem.or.xhead.gt.xtail) then
c          error return, if not enough integer and/or real memory:
           go to 9000
        endif
        bestco = 0
        limit = n + 2*ndn
        lsons = ndn + 1
        usons = ndn + 1
 
c       ----------------------------------------------------------------
c       initialize workspaces
c       ----------------------------------------------------------------
 
        do i = 1, n
           wir (i) = -1
           wic (i) = -2
           head (i) = 0
           wc (i) = 0
           wr (i) = 0
        end do
 
c       ----------------------------------------------------------------
c       initialize the link list for keeping track of real memory usage
c       ----------------------------------------------------------------
 
        mhead = 0
        mtail = 0
        do 20 col = 1, n
           pc = cp (col)
           clen = ii (pc+6)
           if (clen .gt. 0) then
c             place the column in the link list of blocks in xx
              if (mhead .eq. 0) then
                 mhead = pc
              endif
              ii (pc+4) = mtail
              ii (pc+3) = 0
              if (mtail .ne. 0) then
                 ii (mtail+3) = pc
              endif
              mtail = pc
           else
              ii (pc+2) = 0
              ii (pc+4) = 0
              ii (pc+3) = 0
           endif
20      continue
 
c       ----------------------------------------------------------------
c       convert dense columns to a-priori contribution blocks and
c       get the count of nonzeros in each row
c       ----------------------------------------------------------------
 
        e = n
        dnz = 0
        do 50 col = 1, n
           pc = cp (col)
           clen = ii (pc+6)
           cep = (pc+9)
           if (clen .gt. dsiz) then
c             this is a dense column - add to element list length
              dnz = dnz + clen
              do 30 ip = cep, cep + clen - 1
                 row = ii (ip)
                 wr (row) = wr (row) + 1
30            continue
c             convert dense column (in place) into a frontal matrix
              e = e + 1
              ep = pc
              rp (e) = ep
              fdimc = clen
              fleftc = clen
              fleftr = 1
              ii (ep+1) = fdimc
              ii (ep+5) = fleftr
              ii (ep+6) = fleftc
              wr (e) = w0-1
              wc (e) = w0-1
              lurp = (ep+8)
              ii (lurp) = col
              fmax = max (fmax, fleftc)
           else
c             this is a sparse column - add to orig entry length
              do 40 ip = cep, cep + clen - 1
                 row = ii (ip)
                 wc (row) = wc (row) + 1
40            continue
           endif
50      continue
 
c       ----------------------------------------------------------------
c       get memory for row-oriented form, and dense column element lists
c       ----------------------------------------------------------------
 
        pr = ihead
        csiz = cscal + 2
        is = (nz + rscal*n + dnz) + (dn * csiz)
        ihead = ihead + is
        iuse = iuse + is
        ineed = iuse
        xneed = xuse
        info (18) = max (info (18), iuse)
        info (19) = max (info (19), ineed)
        if (ihead .gt. itail) then
c          error return, if not enough integer memory:
           go to 9000
        endif
 
c       ----------------------------------------------------------------
c       if memory is available, add up to dsiz+6 extra slots in the
c       reconstructed dense columns to allow for element list growth
c       ----------------------------------------------------------------
 
        if (dn .gt. 0) then
           extra = min ((itail - ihead) / dn, dsiz + 6)
           csiz = csiz + extra
           is = dn * extra
           ihead = ihead + is
           iuse = iuse + is
           info (18) = max (info (18), iuse)
        endif
 
c       ----------------------------------------------------------------
c       construct row pointers
c       ----------------------------------------------------------------
 
        do 60 row = 1, n
           rp (row) = pr
           rep  = (pr+2)
           reln = wr (row)
           rlen = wc (row)
           rsiz = 2*reln + rlen + rscal
           ii (pr) = rsiz
           rdeg = reln + rlen
           ii (pr+1) = rdeg
           wm (row) = rep
           pr = pr + rsiz
60      continue
 
c       ----------------------------------------------------------------
c       construct row element lists for dense columns
c       ----------------------------------------------------------------
 
        pc = pr
        do 80 e = n+1, n+dn
           ep = rp (e)
           lucp = (ep+9)
           fdimc = ii (ep+1)
cfpp$ nodepchk l
           do f = 0, fdimc - 1
              row = ii (lucp+f)
              ii (wm (row)    ) = e
              ii (wm (row) + 1) = f
              wm (row) = wm (row) + 2
           end do
c          re-construct dense columns as just an element list,
c          containing a single element tuple (e,f), where f = 0
           lurp = (ep+8)
           col = ii (lurp)
           cp (col) = pc
           ii (pc) = csiz
           cdeg = fdimc
           ii (pc+1) = cdeg
           ii (pc+2) = 0
           ii (pc+4) = 0
           ii (pc+3) = 0
           ii (pc+5) = 1
           ii (pc+6) = 0
           ii (pc+7) = 0
           ii (pc+8) = 0
c          store the (e,0) tuple:
           cep = (pc+9)
           ii (cep  ) = e
           ii (cep+1) = 0
           pc = pc + csiz
80      continue
 
c       ----------------------------------------------------------------
c       construct the nonzero pattern of the row-oriented form
c       ----------------------------------------------------------------
 
        do 100 col = 1, n
           pc = cp (col)
           cep = (pc+9)
           clen = ii (pc+6)
cfpp$ nodepchk l
           do 90 p = cep, cep + clen - 1
              row = ii (p)
              ii (wm (row)) = col
              wm (row) = wm (row) + 1
90         continue
100     continue
 
c       count the numerical assembly of the original matrix
        rinfo (2) = rinfo (2) + nz
 
c       ----------------------------------------------------------------
c       initialize the degree lists
c       ----------------------------------------------------------------
 
c       do so in reverse order to try to improve pivot tie-breaking
        do 110 col = n, 1, -1
           pc = cp (col)
           cdeg = ii (pc+1)
           if (cdeg .le. 0) then
c             empty column - remove from pivot search
              cdeg = -(n+2)
              ii (pc+1) = cdeg
           else
              cnxt = head (cdeg)
              ii (pc+7) = cnxt
              ii (pc+8) = 0
              if (cnxt .ne. 0) then
                 ii (cp (cnxt)+8) = col
              endif
              head (cdeg) = col
           endif
110     continue
 
c=======================================================================
c=======================================================================
c  main factorization loop [
c=======================================================================
c=======================================================================
 
        do 1540 dummy1 = 1, n
c       (this loop is not indented due to its length)
 
c       ----------------------------------------------------------------
c       factorization is done if n pivots have been found
c       ----------------------------------------------------------------
 
        if (npiv .ge. n) then
           go to 2000
        endif
 
c=======================================================================
c  global pivot search, and initialization of a new frontal matrix [
c=======================================================================
 
        if (mtail .ne. 0 .and. ii (mtail+6) .eq. 0) then
c          tail block is free, delete it
           xp = ii (mtail+2)
           xuse = xuse - (xhead - xp)
           xhead = xp
           if (mtail .eq. pfree) then
              pfree = 0
              xfree = -1
           endif
           mtail = ii (mtail+4)
           if (mtail .ne. 0) then
              ii (mtail+3) = 0
           else
c             singular matrix.  no columns or contribution blocks left.
              mhead = 0
           endif
        endif
 
c=======================================================================
c  global pivot search:  find pivot row and column
c=======================================================================
 
        nsons = 0
        sonlst = 0
        srched = 0
        pivcol = 0
        slist = 0
 
        do 255 dummy2 = 1, n
 
c          -------------------------------------------------------------
c          get col from column upper-bound degree list
c          -------------------------------------------------------------
 
           col = 0
           do 140 cdeg = mindeg, n
              col = head (cdeg)
              if (col .ne. 0) then
c                exit out of loop if column found:
                 go to 150
              endif
140        continue
           if (col .eq. 0) then
c             exit out of loop if column not found (singular matrix):
              go to 260
           endif
c          loop exit label:
150        continue
           pc = cp (col)
           cnxt = ii (pc+7)
           if (cnxt .ne. 0) then
              ii (cp (cnxt)+8) = 0
           endif
           head (cdeg) = cnxt
           mindeg = cdeg
 
c          -------------------------------------------------------------
c          construct candidate column in wm and xx (wxp..wxp+cdeg-1)
c          -------------------------------------------------------------
 
           xs = cdeg
c          use wm (1..cdeg) for pattern [
c          use xx (wxp..wxp+xs-1) as workspace for values [
 
           if (xs .gt. xtail-xhead) then
 
              info (15) = info (15) + 1
              call ums2fg (xx, xsize, xhead, xtail, xuse,
     $                     ii, isize, ihead, itail, iuse,
     $                     cp, rp, dn, n, icntl, wir, wic, wr, wc,
     $                     0, 0, 0, 0, .false.,
     $                     pfree, xfree, mhead, mtail, slots)
c             at this point, iuse = ineed and xuse = xneed
              pc = cp (col)
           endif
 
           wxp = xhead
           xhead = xhead + xs
           xuse = xuse + xs
           xneed = xneed + xs
           info (20) = max (info (20), xuse)
           info (21) = max (info (21), xneed)
           if (xhead .gt. xtail) then
c             error return, if not enough real memory:
              go to 9000
           endif
 
c          -------------------------------------------------------------
c          assemble the elements in the element list
c          -------------------------------------------------------------
 
           cdeg = 0
           cep = (pc+9)
           celn = ii (pc+5)
           do 190 ip = cep, cep + 2*celn - 2, 2
              e = ii (ip)
              f = ii (ip+1)
              ep = rp (e)
              fdimc = ii (ep+1)
              fxp = ii (ep+2)
              if (e .le. n) then
                 fluip = ii (ep)
                 ludegc = ii (fluip+3)
                 lucp = (fluip + 7)
              else
                 ludegc = fdimc
                 lucp = (ep+9)
              endif
              xp = fxp + f * fdimc
c             split into 3 loops so that they all vectorize on a cray
              cdeg1 = cdeg
              do 160 p = lucp, lucp + ludegc - 1
                 row = ii (p)
                 if (row .gt. 0) then
                    if (wir (row) .le. 0) then
                       cdeg = cdeg + 1
                       wm (cdeg) = row
                    endif
                 endif
160           continue
              do 170 i = cdeg1+1, cdeg
                 row = wm (i)
                 wir (row) = i
                 xx (wxp+i-1) = zero
170           continue
cfpp$ nodepchk l
              do 180 j = 0, ludegc - 1
                 row = ii (lucp+j)
                 if (row .gt. 0) then
                    xx (wxp + wir (row) - 1) =
     $              xx (wxp + wir (row) - 1) + xx (xp+j)
                 endif
180           continue
190        continue
 
c          -------------------------------------------------------------
c          assemble the original entries in the column
c          -------------------------------------------------------------
 
           cdeg1 = cdeg
           clen = ii (pc+6)
           csiz = ii (pc)
           ip = pc + csiz - clen
           cxp = ii (pc+2)
cfpp$ nodepchk l
           do i = 0, clen - 1
              row = ii (ip+i)
              wm (cdeg+1+i) = row
              xx (wxp+cdeg+i) = xx (cxp+i)
           end do
           cdeg = cdeg + clen
 
c          -------------------------------------------------------------
c          update the degree of this column (exact, not upper bound)
c          -------------------------------------------------------------
 
           ii (pc+1) = cdeg
 
c          wm (1..cdeg) holds the pattern of col being searched.
c          xx (wxp..wxp+cdeg-1) holds the numerical values of col being
c          searched.  wir (wm (1..cdeg1)) is 1..cdeg1.
 
c          -------------------------------------------------------------
c          find the maximum absolute value in the column
c          -------------------------------------------------------------
 
           maxval = abs (xx (wxp - 1 + isamax (cdeg, xx (wxp), 1)))
           rinfo (3) = rinfo (3) + cdeg
           toler = relpt * maxval
           rdeg = n+1
 
c          -------------------------------------------------------------
c          look for the best possible pivot row in this column
c          -------------------------------------------------------------
 
           if (cdeg .ne. 0 .and. maxval .gt. zero) then
              if (symsrc) then
c                prefer symmetric pivots, if numerically acceptable
                 row = col
                 rowpos = wir (row)
                 if (rowpos .le. 0) then
c                   diagonal may be in original entries
                    do 210 i = cdeg1 + 1, cdeg1 + clen
                       if (wm (i) .eq. row) then
                          rowpos = i
c                         exit out of loop if symmetric pivot found:
                          go to 220
                       endif
210                 continue
c                   loop exit label:
220                 continue
                 endif
                 if (rowpos .gt. 0) then
c                   diagonal entry exists in the column pattern
                    x = abs (xx (wxp-1+rowpos))
                    if (x .ge. toler .and. x .gt. zero) then
c                      diagonal entry is numerically acceptable
                       pr = rp (row)
                       rdeg = ii (pr+1)
                    endif
                 endif
              endif
              if (rdeg .eq. n+1) then
c                continue searching - no diagonal found or sought for.
c                minimize row degree subject to abs(value) constraints.
                 row = n+1
                 do 230 i = 1, cdeg
                    row2 = wm (i)
                    pr = rp (row2)
                    rdeg2 = ii (pr+1)
c                   among those numerically acceptable rows of least
c                   (upper bound) degree, select the row with the
c                   lowest row index
                    better = rdeg2 .lt. rdeg .or.
     $                      (rdeg2 .eq. rdeg .and. row2 .lt. row)
                    if (better) then
                       x = abs (xx (wxp-1+i))
                       if (x.ge.toler .and. x.gt.zero) then
                          row = row2
                          rdeg = rdeg2
                          rowpos = i
                       endif
                    endif
230              continue
              endif
           endif
 
c          -------------------------------------------------------------
c          deallocate workspace
c          -------------------------------------------------------------
 
           xhead = xhead - xs
           xuse = xuse - xs
           xneed = xneed - xs
c          done using xx (wxp...wxp+xs-1) ]
 
c          -------------------------------------------------------------
c          reset work vector
c          -------------------------------------------------------------
 
           do 240 i = 1, cdeg1
              wir (wm (i)) = -1
240        continue
 
c          -------------------------------------------------------------
c          check to see if a pivot column was found
c          -------------------------------------------------------------
 
           if (rdeg .eq. n+1) then
 
c             ----------------------------------------------------------
c             no pivot found, column is zero
c             ----------------------------------------------------------
 
c             remove this singular column from any further pivot search
              cdeg = -(n+2)
              ii (pc+1) = cdeg
 
           else
 
c             ----------------------------------------------------------
c             save a list of the columns searched (with nonzero degrees)
c             ----------------------------------------------------------
 
              srched = srched + 1
              ii (pc+7) = slist
              slist = col
 
c             ----------------------------------------------------------
c             check if this is the best pivot seen so far
c             ----------------------------------------------------------
 
c             compute the true markowitz cost without scanning the row
c             wm (1..cdeg) holds pivot column, including pivot row index
c             wm (rowpos) contains the candidate pivot row index
              cost = (cdeg - 1) * (rdeg - 1)
              if (pivcol .eq. 0 .or. cost .lt. bestco) then
                 fflefc = cdeg
                 do 250 i = 1, fflefc-1
                    wpc (i) = wm (i)
250              continue
c                remove the pivot row index from pivot column pattern
                 wpc (rowpos) = wm (fflefc)
                 pivcol = col
                 pivrow = row
                 bestco = cost
              endif
           endif
 
c          done using wm (1..cdeg) for pattern ]
c          wpc (1..fflefc-1) holds pivot column (excl. pivot row index)
 
c          -------------------------------------------------------------
c          exit global pivot search if nsrch pivots have been searched
c          -------------------------------------------------------------
 
           if (srched .ge. nsrch) then
              go to 260
           endif
 
255     continue
c       exit label for loop 255:
260     continue
 
c=======================================================================
c  quit early if no pivot found (singular matrix detected)
c=======================================================================
 
        if (pivcol .eq. 0) then
c          complete the column permutation vector in
c          wpc (n-npiv+1 ... n) in reverse order
           k = n - npiv + 1
           do 270 col = 1, n
              if (cp (col) .ne. 0) then
c                this is a non-pivotal column
                 k = k - 1
                 wpc (k) = col
                 cp (col) = 0
              endif
270        continue
c          complete the row permutation vector in
c          wpr (n-npiv+1 ... n) in reverse order
           k = n - npiv + 1
           do 280 row = 1, ndn
              if (row .gt. n) then
c                this is an artificial frontal matrix
                 e = row
                 rp (e) = 0
              else if (rp (row) .ne. 0) then
                 rlen = wc (row)
                 if (rlen .ge. 0 .and. rlen .le. n) then
c                   this is a non-pivotal row
                    k = k - 1
                    wpr (k) = row
                    rp (row) = 0
                 else if (rlen .ne. -(ndn+2)) then
c                   this is an unassembled element: convert to lu
                    e = row
                    ep = rp (row)
                    wr (e) = -(ndn+2)
                    wc (e) = -(ndn+2)
                    fluip = ii (ep)
                    rp (e) = fluip
                 endif
              endif
280        continue
c          factorization is done, exit the main factorization loop:
           go to 2000
        endif
 
c=======================================================================
c  place the non-pivotal columns searched back in degree lists
c=======================================================================
 
        do 300 i = 1, srched
           col = slist
           pc = cp (col)
           slist = ii (pc+7)
           if (col .ne. pivcol) then
              cdeg = ii (pc+1)
              cnxt = head (cdeg)
              ii (pc+7) = cnxt
              ii (pc+8) = 0
              if (cnxt .ne. 0) then
                 ii (cp (cnxt)+8) = col
              endif
              head (cdeg) = col
              mindeg = min (mindeg, cdeg)
           endif
300     continue
 
c=======================================================================
c  construct pivot row pattern
c=======================================================================
 
c       at this point, wir (1..n) = -1 and wic (1..n) is -2 for
c       nonpivotal columns and -1 for pivotal columns.
c       wic (wpr (1..fflefr+1)) is set to zero in the code below.  it
c       will be set to the proper offsets in do 775, once ffdimc is
c       known (offsets are dependent on ffdimc, which is dependent on
c       fflefr calculated below, and the memory allocation).
 
c       ----------------------------------------------------------------
c       assemble the elements in the element list
c       ----------------------------------------------------------------
 
        pr = rp (pivrow)
        fflefr = 0
        rep = (pr+2)
        reln = wr (pivrow)
        do 330 ip = rep, rep + 2*reln - 2, 2
           e = ii (ip)
           ep = rp (e)
           if (e .le. n) then
              fluip = ii (ep)
              lucp = (fluip + 7)
              ludegr = ii (fluip+2)
              ludegc = ii (fluip+3)
              lurp = lucp + ludegc
c             split into two loops so that they both vectorize on a cray
              f1 = fflefr
              do 310 p = lurp, lurp + ludegr - 1
                 col = ii (p)
                 if (col .gt. 0) then
                    if (wic (col) .eq. -2) then
                       fflefr = fflefr + 1
                       wpr (fflefr) = col
                    endif
                 endif
310           continue
              do 320 i = f1+1, fflefr
                 wic (wpr (i)) = 0
320           continue
           else
c             this is an artifical element (a dense column)
              lurp = (ep+8)
              col = ii (lurp)
              if (wic (col) .eq. -2) then
                 fflefr = fflefr + 1
                 wpr (fflefr) = col
                 wic (col) = 0
              endif
           endif
330     continue
 
c       ----------------------------------------------------------------
c       assemble the original entries in the pivot row
c       ----------------------------------------------------------------
 
        rsiz = ii (pr)
        rlen = wc (pivrow)
        do 340 p = pr + rsiz - rlen, pr + rsiz - 1
           col = ii (p)
           if (wic (col) .eq. -2) then
              fflefr = fflefr + 1
              wpr (fflefr) = col
           endif
340     continue
c       the exact degree of the pivot row is fflefr
 
c=======================================================================
c  initialize the new frontal matrix
c=======================================================================
 
c       ffrow is the name of current frontal matrix
        ffrow = pivrow
        e1 = pivrow
        k = 1
        k0 = 0
        ffdimr = min (kleft, int (gro * fflefr))
        ffdimc = min (kleft, int (gro * fflefc))
        fmaxr = fflefr
        fmaxc = fflefc
        ffsize = ffdimc * ffdimr
        rscan = max (dsiz, ffdimr)
 
c       ----------------------------------------------------------------
c       compute the offsets for rows in the pivot column
c       and the offsets for columns in the pivot row
c       ----------------------------------------------------------------
 
        do i = 1, fflefc - 1
           wir (wpc (i)) = i - 1
        end do
        do 360 i = 1, fflefr
           wic (wpr (i)) = (i - 1) * ffdimc
360     continue
 
c       ----------------------------------------------------------------
c       remove the pivot column index from the pivot row pattern
c       ----------------------------------------------------------------
 
        col = wpr (fflefr)
        colpos = (wic (pivcol)/ffdimc)+1
        wpr (colpos) = col
        wic (col) = wic (pivcol)
        wic (pivcol) = (ffdimr - 1) * ffdimc
        wir (pivrow) = ffdimc - 1
 
c       ----------------------------------------------------------------
c       remove the pivot row/col from the nonzero count
c       ----------------------------------------------------------------
 
        fflefr = fflefr - 1
        fflefc = fflefc - 1
 
c       ----------------------------------------------------------------
c       allocate the working array, doing garbage collection if needed
c       also allocate space for a work vector of size ffdimc
c       ----------------------------------------------------------------
 
        if (ffsize + ffdimc .gt. xtail-xhead) then
           info (15) = info (15) + 1
           call ums2fg (xx, xsize, xhead, xtail, xuse,
     $                  ii, isize, ihead, itail, iuse,
     $                  cp, rp, dn, n, icntl, wir, wic, wr, wc,
     $                  0, 0, 0, 0, .false.,
     $                  pfree, xfree, mhead, mtail, slots)
c          at this point, iuse = ineed and xuse = xneed
        endif
 
        ffxp = xhead
        xhead = xhead + ffsize
        wxp = xhead
        xhead = xhead + ffdimc
        xuse = xuse + ffsize + ffdimc
        xneed = xneed + ffsize + ffdimc
        info (20) = max (info (20), xuse)
        info (21) = max (info (21), xneed)
        if (xhead .gt. xtail) then
c          error return, if not enough real memory:
           go to 9000
        endif
 
c       ----------------------------------------------------------------
c       get memory usage for next call to ums2rf
c       ----------------------------------------------------------------
 
        xruse = xruse + ffsize
        xrmax = max (xrmax, xruse)
 
c       ----------------------------------------------------------------
c       zero the working array
c       ----------------------------------------------------------------
 
c       zero the contribution block:
        do 380 j = 0, fflefr - 1
           do 370 i = 0, fflefc - 1
              xx (ffxp + j*ffdimc + i) = zero
370        continue
380     continue
 
c       zero the pivot row:
        do 390 j = 0, fflefr - 1
           xx (ffxp + j*ffdimc + ffdimc-1) = zero
390     continue
 
c       zero the pivot column:
        do 400 i = 0, fflefc - 1
           xx (ffxp + (ffdimr-1)*ffdimc + i) = zero
400     continue
 
c       zero the pivot entry itself:
        xx (ffxp + (ffdimr-1)*ffdimc + ffdimc-1) = zero
 
c       ----------------------------------------------------------------
c       current workspace usage:
c       ----------------------------------------------------------------
 
c       wpc (1..fflefc):        holds the pivot column pattern
c                               (excluding the pivot row index)
c       wpc (fflefc+1 .. n-npiv):       unused
c       wpc (n-npiv+1 .. n):            pivot columns in reverse order
c
c       wpr (1..fflefr):        holds the pivot row pattern
c                               (excluding the pivot column index)
c       wpr (fflefr+1 .. n-npiv):       unused
c       wpr (n-npiv+1 .. n):            pivot rows in reverse order
c
c       c (1..ffdimr, 1..ffdimc):  space for the new frontal matrix.
c
c       c (i,j) is located at xx (ffxp+((i)-1)+((j)-1)*ffdimc)
c
c       wir (row) >= 0 for each row in pivot column pattern.
c               offset into pattern is given by:
c               wir (row) == offset - 1
c               also, wir (pivrow) is ffdimc-1, the offset in c of
c               the pivot row itself.
c               otherwise, wir (1..n) is -1
c
c       wic (col) >= 0 for each col in pivot row pattern.
c               wic (col) == (offset - 1) * ffdimc
c               also, wic (pivcol) is (ffdimr-1)*ffdimc,
c               the offset in c of the pivot column itself.
c               otherwise, wic (1..n) is -2 for nonpivotal columns,
c               and -1 for pivotal columns
 
c       ----------------------------------------------------------------
c       remove the columns affected by this element from degree lists
c       ----------------------------------------------------------------
 
        do 410 j = 1, fflefr
           pc = cp (wpr (j))
           cdeg = ii (pc+1)
           if (cdeg .gt. 0) then
              cnxt = ii (pc+7)
              cprv = ii (pc+8)
              if (cnxt .ne. 0) then
                 ii (cp (cnxt)+8) = cprv
              endif
              if (cprv .ne. 0) then
                 ii (cp (cprv)+7) = cnxt
              else
                 head (cdeg) = cnxt
              endif
           endif
410     continue
 
c=======================================================================
c  initialization of new frontal matrix is complete ]
c=======================================================================
 
c=======================================================================
c  assemble and factorize the current frontal matrix [
c=======================================================================
 
c       for first pivot in frontal matrix, do all scans
        scan1 = 0
        scan2 = 0
        scan3 = 0
        scan4 = 0
 
        do 1395 dummy3 = 1, n
c       (this loop is not indented due to its length)
 
c=======================================================================
c  degree update and numerical assembly [
c=======================================================================
 
        kleft1 = kleft - 1
 
c       ----------------------------------------------------------------
c       scan1:  scan the element lists of each row in the pivot col
c               and compute the external column degree for each frontal
c       ----------------------------------------------------------------
 
        row = pivrow
        do 440 j = scan1, fflefc
           if (j .ne. 0) then
c             get a row;  otherwise, scan the pivot row if j is zero.
              row = wpc (j)
           endif
           pr = rp (row)
           rep = (pr+2)
           reln = wr (row)
cfpp$ nodepchk l
           do 430 p = rep, rep + 2*reln - 2, 2
              e = ii (p)
              if (wc (e) .lt. w0) then
c                this is the first time seen in either scan 1 or 2:
                 ep = rp (e)
                 fleftr = ii (ep+5)
                 fleftc = ii (ep+6)
                 wr (e) = fleftr + w0
                 wc (e) = fleftc + w0
              endif
              wc (e) = wc (e) - 1
430        continue
440     continue
 
c       ----------------------------------------------------------------
c       scan2:  scan the element lists of each col in the pivot row
c               and compute the external row degree for each frontal
c       ----------------------------------------------------------------
 
        col = pivcol
        do 460 j = scan2, fflefr
           if (j .ne. 0) then
c             get a col;  otherwise, scan the pivot col if j is zero.
              col = wpr (j)
           endif
           pc = cp (col)
           celn = ii (pc+5)
           cep = (pc+9)
cfpp$ nodepchk l
           do 450 p = cep, cep + 2*celn - 2, 2
              e = ii (p)
              if (wr (e) .lt. w0) then
c                this is the first time seen in either scan 1 or 2:
                 ep = rp (e)
                 fleftr = ii (ep+5)
                 fleftc = ii (ep+6)
                 wr (e) = fleftr + w0
                 wc (e) = fleftc + w0
              endif
              wr (e) = wr (e) - 1
450        continue
460     continue
 
c       ----------------------------------------------------------------
c       scan3:  scan the element lists of each column in pivot row
c               do degree update for the columns
c               assemble effective usons and lu-sons
c       ----------------------------------------------------------------
 
c       flag usons in wc (e) as scanned (all now unflagged) [
c       uses wc (e) for the link list.  wc (e) <= 0
c       means that e is in the list, the external column
c       degree is zero, and -(wc (e)) is the next element in
c       the uson list.
 
        col = pivcol
        do 700 jj = scan3, fflefr
 
c          -------------------------------------------------------------
c          assemble and update the degree of a column
c          -------------------------------------------------------------
 
           if (jj .ne. 0) then
c             get a col;  otherwise, scan the pivot col if jj is zero
              col = wpr (jj)
           endif
 
c          -------------------------------------------------------------
c          compute the degree, and partition the element list into
c          two parts.  the first part are not lusons or usons, and
c          are not assembled.  the second part is assembled.
c          -------------------------------------------------------------
 
           cdeg = 0
           deln = 0
           pc = cp (col)
           cep = (pc+9)
           celn = ii (pc+5)
           ip2 = cep + 2*celn - 2
           xudp = ffxp + wic (col)
cfpp$ nodepchk l
           do 470 ip = cep, ip2, 2
              e = ii (ip)
              if (wc (e) .gt. w0) then
c                this element cannot be assembled
                    cdeg = cdeg + (wc (e) - w0)
              else
c                delete this tuple from the element list
                 deln = deln + 1
                 wm (deln) = ip
              endif
470       continue
 
          if (deln .ne. 0) then
 
c             ----------------------------------------------------------
c             move the deleted tuples to the end of the element list
c             ----------------------------------------------------------
 
              p2 = ip2
              do i = deln, 1, -1
                 e = ii (wm (i)  )
                 f = ii (wm (i)+1)
                 ii (wm (i)  ) = ii (p2  )
                 ii (wm (i)+1) = ii (p2+1)
                 ii (p2  ) = e
                 ii (p2+1) = f
                 p2 = p2 - 2
              end do
 
c             ----------------------------------------------------------
c             assemble from lusons and usons (the deleted tuples)
c             ----------------------------------------------------------
 
              do 670 ip = p2 + 2, ip2, 2
 
c                -------------------------------------------------------
c                this is an luson or uson.  if fextc < 0 then this has
c                already been assembled.
c                -------------------------------------------------------
 
                 e = ii (ip)
                 if (wc (e) .lt. w0) then
c                   go to next iteration if already assembled
                    goto 670
                 endif
 
c                -------------------------------------------------------
c                get scalar info, add son to list if not already there
c                -------------------------------------------------------
 
                 ep = rp (e)
                 fdimc = ii (ep+1)
                 fxp = ii (ep+2)
                 fleftr = ii (ep+5)
                 fleftc = ii (ep+6)
                 if (e .le. n) then
                    fluip = ii (ep)
                    ludegr = ii (fluip+2)
                    ludegc = ii (fluip+3)
                    lucp = (fluip + 7)
                    lurp = lucp + ludegc
                    if (wir (e) .eq. -1) then
                       wir (e) = sonlst - n - 2
                       sonlst = e
                       nsons = nsons + 1
                    endif
                 else
c                   an artificial frontal matrix
                    ludegr = 1
                    ludegc = fdimc
                    lucp = (ep+9)
                    lurp = (ep+8)
                 endif
 
c                -------------------------------------------------------
                 if (wr (e) .eq. w0) then
c                this is an luson - assemble an entire frontal matrix
c                -------------------------------------------------------
 
c                   ----------------------------------------------------
                    if (ludegc .eq. fleftc) then
c                   no rows assembled out of this frontal yet
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
c                      use wm (1..ludegc for offsets) [
                       do 490 i = 0, ludegc-1
                          row2 = ii (lucp+i)
                          wm (i+1) = wir (row2)
490                    continue
 
c                      -------------------------------------------------
                       if (ludegr .eq. fleftr) then
c                      no rows or cols assembled out of frontal yet
c                      -------------------------------------------------
 
                          do 510 j = 0, ludegr-1
                             col2 = ii (lurp+j)
                             xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                             do 500 i = 0, ludegc-1
                                xx (xdp + wm (i+1)) =
     $                          xx (xdp + wm (i+1)) +
     $                          xx (fxp + j*fdimc + i)
500                          continue
510                       continue
 
c                      -------------------------------------------------
                       else
c                      only cols have been assembled out of frontal
c                      -------------------------------------------------
 
                          do 530 j = 0, ludegr-1
                             col2 = ii (lurp+j)
                             if (col2 .gt. 0) then
                                xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                                do 520 i = 0, ludegc-1
                                   xx (xdp + wm (i+1)) =
     $                             xx (xdp + wm (i+1)) +
     $                             xx (fxp + j*fdimc + i)
520                             continue
                             endif
530                       continue
                       endif
c                      done using wm (1..ludegc for offsets) ]
 
c                   ----------------------------------------------------
                    else
c                   some rows have been assembled out of this frontal
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
c                      use wm (1..ludegc for offsets) [
                       degc = 0
                       do 540 i = 0, ludegc-1
                          row2 = ii (lucp+i)
                          if (row2 .gt. 0) then
                             degc = degc + 1
                             wj (degc) = i
                             wm (degc) = wir (row2)
                          endif
540                    continue
 
c                      -------------------------------------------------
                       if (ludegr .eq. fleftr) then
c                      only rows assembled out of this frontal
c                      -------------------------------------------------
 
                          do 560 j = 0, ludegr-1
                             col2 = ii (lurp+j)
                             xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                             do 550 i = 1, degc
                                xx (xdp + wm (i)) =
     $                          xx (xdp + wm (i)) +
     $                          xx (fxp + j*fdimc + wj (i))
550                          continue
560                       continue
 
c                      -------------------------------------------------
                       else
c                      both rows and columns assembled out of frontal
c                      -------------------------------------------------
 
                          do 580 j = 0, ludegr-1
                             col2 = ii (lurp+j)
                             if (col2 .gt. 0) then
                                xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                                do 570 i = 1, degc
                                   xx (xdp + wm (i)) =
     $                             xx (xdp + wm (i)) +
     $                             xx (fxp + j*fdimc + wj (i))
570                             continue
                             endif
580                       continue
                       endif
c                      done using wm (1..ludegc for offsets) ]
                    endif
 
c                   ----------------------------------------------------
c                   deallocate the luson frontal matrix
c                   ----------------------------------------------------
 
                    wr (e) = -(ndn+2)
                    wc (e) = -(ndn+2)
                    if (e .le. n) then
                       rp (e) = fluip
                       ii (ep) = fscal
                       ineed = ineed - fscal
                    else
                       rp (e) = 0
                       ii (ep) = fdimc + cscal
                       ineed = ineed - (fdimc + cscal)
                    endif
                    ii (ep+1) = -1
                    ii (ep+6) = 0
                    mprev = ii (ep+4)
                    mnext = ii (ep+3)
                    xneed = xneed - ludegr*ludegc
                    if (mnext .ne. 0 .and. ii (mnext+6) .eq. 0) then
c                      next block is free - delete it
                       mnext = ii (mnext+3)
                       ii (ep+3) = mnext
                       if (mnext .ne. 0) then
                          ii (mnext+4) = ep
                       else
                          mtail = ep
                       endif
                    endif
                    if (mprev .ne. 0 .and. ii (mprev+6) .eq. 0) then
c                      previous block is free - delete it
                       ii (ep+2) = ii (mprev+2)
                       mprev = ii (mprev+4)
                       ii (ep+4) = mprev
                       if (mprev .ne. 0) then
                          ii (mprev+3) = ep
                       else
                          mhead = ep
                       endif
                    endif
c                   get the size of the freed block
                    if (mnext .ne. 0) then
                       xs = ii (mnext+2) - ii (ep+2)
                    else
                       xs = ffxp - ii (ep+2)
                    endif
                    if (xs .gt. xfree) then
c                      keep track of the largest free block
                       xfree = xs
                       pfree = ep
                    endif
 
c                   ----------------------------------------------------
c                   get memory usage for next call to ums2rf
c                   ----------------------------------------------------
 
                    xruse = xruse - ludegr*ludegc
 
c                -------------------------------------------------------
                 else if (wr (e) - w0 .le. fleftr/2) then
c                this is a uson - assemble all possible columns
c                -------------------------------------------------------
 
c                   ----------------------------------------------------
c                   add to uson list - to be cleared just after scan 3
c                   ----------------------------------------------------
 
                    wc (e) = -usons
                    usons = e
 
c                   ----------------------------------------------------
                    if (ludegc .eq. fleftc) then
c                   no rows assembled out of this uson frontal yet
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
c                      use wm (1..ludegc for offsets)
                       do 590 i = 0, ludegc-1
                          row2 = ii (lucp+i)
                          wm (i+1) = wir (row2)
590                    continue
 
                       do 610 j = 0, ludegr-1
                          col2 = ii (lurp+j)
                          if (col2 .gt. 0) then
                             if (wic (col2) .ge. 0) then
                                xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                                do 600 i = 0, ludegc-1
                                   xx (xdp + wm (i+1)) =
     $                             xx (xdp + wm (i+1)) +
     $                             xx (fxp + j*fdimc + i)
600                             continue
c                               flag this column as assembled from uson
                                ii (lurp+j) = -col2
                             endif
                          endif
610                    continue
 
c                   ----------------------------------------------------
                    else
c                   some rows already assembled out of this uson frontal
c                   ----------------------------------------------------
 
c                      compute the compressed column offset vector
c                      use wm (1..ludegc for offsets)
                       degc = 0
                       do 620 i = 0, ludegc-1
                          row2 = ii (lucp+i)
                          if (row2 .gt. 0) then
                             degc = degc + 1
                             wj (degc) = i
                             wm (degc) = wir (row2)
                          endif
620                    continue
 
                       do 640 j = 0, ludegr-1
                          col2 = ii (lurp+j)
                          if (col2 .gt. 0) then
                             if (wic (col2) .ge. 0) then
                                xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                                do 630 i = 1, degc
                                   xx (xdp + wm (i)) =
     $                             xx (xdp + wm (i)) +
     $                             xx (fxp + j*fdimc + wj (i))
630                             continue
c                               flag this column as assembled from uson
                                ii (lurp+j) = -col2
                             endif
                          endif
640                    continue
 
                    endif
 
                    fleftr = wr (e) - w0
                    ii (ep+5) = fleftr
 
c                -------------------------------------------------------
                 else
c                this is a uson - assemble just one column
c                -------------------------------------------------------
 
c                   get the offset, f, from the (e,f) tuple
                    f = ii (ip+1)
 
c                   ----------------------------------------------------
                    if (ludegc .eq. fleftc) then
c                   no rows assembled out of this uson yet
c                   ----------------------------------------------------
 
cfpp$ nodepchk l
                       do 650 i = 0, ludegc-1
                          row2 = ii (lucp+i)
                          xx (xudp + wir (row2)) =
     $                    xx (xudp + wir (row2)) +
     $                    xx (fxp + f*fdimc + i)
650                    continue
 
c                   ----------------------------------------------------
                    else
c                   some rows already assembled out of this uson
c                   ----------------------------------------------------
 
cfpp$ nodepchk l
                       do 660 i = 0, ludegc-1
                          row2 = ii (lucp+i)
                          if (row2 .gt. 0) then
                             xx (xudp + wir (row2)) =
     $                       xx (xudp + wir (row2)) +
     $                       xx (fxp + f*fdimc + i)
                          endif
660                    continue
                    endif
 
c                   ----------------------------------------------------
c                   decrement count of unassembled cols in frontal
c                   ----------------------------------------------------
 
                    ii (ep+5) = fleftr - 1
c                   flag the column as assembled from the uson
                    ii (lurp+f) = -col
                 endif
 
670           continue
 
c             ----------------------------------------------------------
c             update the count of (e,f) tuples in the element list
c             ----------------------------------------------------------
 
              ii (pc+5) = ii (pc+5) - deln
              ineed = ineed - 2*deln
           endif
 
c          -------------------------------------------------------------
c          assemble the original column and update count of entries
c          -------------------------------------------------------------
 
           clen = ii (pc+6)
           if (clen .gt. 0) then
              csiz = ii (pc)
              ip = pc + csiz - clen
              dlen = 0
cfpp$ nodepchk l
              do 680 i = 0, clen - 1
                 row = ii (ip+i)
                 if (wir (row) .ge. 0) then
c                   this entry can be assembled and deleted
                    dlen = dlen + 1
                    wm (dlen) = i
                 endif
680           continue
              if (dlen .ne. 0) then
                 cxp = ii (pc+2)
                 do 690 j = 1, dlen
                    i = wm (j)
                    row = ii (ip+i)
c                   assemble the entry
                    xx (xudp + wir (row)) =
     $              xx (xudp + wir (row)) + xx (cxp+i)
c                   and delete the entry
                    ii (ip +i) = ii (ip +j-1)
                    xx (cxp+i) = xx (cxp+j-1)
690              continue
                 clen = clen - dlen
                 cxp = cxp + dlen
                 ineed = ineed - dlen
                 xneed = xneed - dlen
                 ii (pc+6) = clen
                 if (clen .ne. 0) then
                    ii (pc+2) = cxp
                 else
c                   deallocate the real portion of the column:
                    mprev = ii (pc+4)
                    mnext = ii (pc+3)
                    if (mnext .ne. 0 .and. ii (mnext+6) .eq. 0) then
c                      next block is free - delete it
                       mnext = ii (mnext+3)
                       ii (pc+3) = mnext
                       if (mnext .ne. 0) then
                          ii (mnext+4) = pc
                       else
                          mtail = pc
                       endif
                    endif
                    if (mprev .ne. 0 .and. ii (mprev+6) .eq. 0) then
c                      previous block is free - delete it
                       ii (pc+2) = ii (mprev+2)
                       mprev = ii (mprev+4)
                       ii (pc+4) = mprev
                       if (mprev .ne. 0) then
                          ii (mprev+3) = pc
                       else
                          mhead = pc
                       endif
                    endif
                    if (pc .eq. mhead) then
c                      adjust the start of the block if this is head
                       ii (pc+2) = 1
                    endif
c                   get the size of the freed block
                    if (mnext .ne. 0) then
                       xs = ii (mnext+2) - ii (pc+2)
                    else
                       xs = ffxp - ii (pc+2)
                    endif
                    if (xs .gt. xfree) then
c                      keep track of the largest free block
                       xfree = xs
                       pfree = pc
                    endif
                 endif
              endif
              cdeg = cdeg + clen
           endif
 
c          -------------------------------------------------------------
c          compute the upper bound degree - excluding current front
c          -------------------------------------------------------------
 
           cdeg2 = ii (pc+1)
           cdeg = min (kleft1 - fflefc, cdeg2, cdeg)
           ii (pc+1) = cdeg
 
700     continue
 
c       ----------------------------------------------------------------
c       scan-3 wrap-up:  remove flags from assembled usons
c       ----------------------------------------------------------------
 
c       while (usons .ne. ndn+1) do
710     continue
        if (usons .ne. ndn+1) then
           next = -wc (usons)
           wc (usons) = w0
           usons = next
c       end while:
        goto 710
        endif
c       done un-flagging usons, all now unflagged in wc (e) ]
 
c       ----------------------------------------------------------------
c       scan4:  scan element lists of each row in the pivot column
c               do degree update for the rows
c               assemble effective lsons
c       ----------------------------------------------------------------
 
c       flag lsons in wr (e) (all are now unflagged) [
c       uses wr (e) for the link list.  wr (e) <= 0 means
c       that e is in the list, the external row degree is zero, and
c       -(wr (e)) is the next element in the lson list.
 
        row = pivrow
        do 840 jj = scan4, fflefc
 
c          -------------------------------------------------------------
c          assemble and update the degree of a row
c          -------------------------------------------------------------
 
           if (jj .ne. 0) then
c             get a row;  otherwise, scan the pivot row if jj is zero
              row = wpc (jj)
           endif
 
c          -------------------------------------------------------------
c          compute the degree, and partition the element list into
c          two parts.  the first part are not lusons or lsons, and
c          are not assembled.  the second part is assembled.
c          -------------------------------------------------------------
 
           rdeg = 0
           deln = 0
           pr = rp (row)
           rep = (pr+2)
           reln = wr (row)
           ip2 = rep + 2*reln - 2
cfpp$ nodepchk l
           do 720 ip = rep, ip2, 2
              e = ii (ip)
              if (wr (e) .gt. w0) then
                 rdeg = rdeg + (wr (e) - w0)
              else
                 deln = deln + 1
                 wm (deln) = ip
              endif
720        continue
 
           if (deln .ne. 0) then
 
c             ----------------------------------------------------------
c             move the deleted tuples to the end of the element list
c             ----------------------------------------------------------
 
              p2 = ip2
              do 730 i = deln, 1, -1
                 e = ii (wm (i)  )
                 f = ii (wm (i)+1)
                 ii (wm (i)  ) = ii (p2  )
                 ii (wm (i)+1) = ii (p2+1)
                 ii (p2  ) = e
                 ii (p2+1) = f
                 p2 = p2 - 2
730           continue
 
c             ----------------------------------------------------------
c             assemble from lsons (the deleted tuples)
c             ----------------------------------------------------------
 
              do 810 ip = p2 + 2, ip2, 2
 
c                -------------------------------------------------------
c                this is an luson or lson.  if fextr < 0 then this has
c                already been assembled.  all lusons have already been
c                assembled (in scan3, above).
c                -------------------------------------------------------
 
                 e = ii (ip)
                 if (wr (e) .lt. w0) then
c                   go to next iteration if already assembled
                    goto 810
                 endif
 
c                -------------------------------------------------------
c                get scalar info, add to son list if not already there
c                -------------------------------------------------------
 
                 ep = rp (e)
                 fdimc = ii (ep+1)
                 fxp = ii (ep+2)
                 fleftr = ii (ep+5)
                 fleftc = ii (ep+6)
                 if (e .le. n) then
                    fluip = ii (ep)
                    ludegr = ii (fluip+2)
                    ludegc = ii (fluip+3)
                    lucp = (fluip + 7)
                    lurp = lucp + ludegc
                    if (wir (e) .eq. -1) then
                       wir (e) = sonlst - n - 2
                       sonlst = e
                       nsons = nsons + 1
                    endif
                 else
c                   an artificial frontal matrix
                    ludegr = 1
                    ludegc = fdimc
                    lucp = (ep+9)
                    lurp = (ep+8)
                 endif
 
c                -------------------------------------------------------
                 if (wc (e) - w0 .le. fleftc/2) then
c                this is an lson - assemble all possible rows
c                -------------------------------------------------------
 
c                   ----------------------------------------------------
c                   add to lson list - to be cleared just after scan 4
c                   ----------------------------------------------------
 
                    wr (e) = -lsons
                    lsons = e
 
c                   compute the compressed column offset vector
c                   use wm (1..ludegc for offsets) [
                    degc = 0
                    do 740 i = 0, ludegc-1
                       row2 = ii (lucp+i)
                       if (row2 .gt. 0) then
                          if (wir (row2) .ge. 0) then
c                            this row will be assembled in loop below
                             degc = degc + 1
                             wj (degc) = i
                             wm (degc) = wir (row2)
c                            flag the row as assembled from the lson
                             ii (lucp+i) = -row2
                          endif
                       endif
740                 continue
 
c                   ----------------------------------------------------
                    if (ludegr .eq. fleftr) then
c                   no columns assembled out this lson yet
c                   ----------------------------------------------------
 
                       do 760 j = 0, ludegr-1
                          col2 = ii (lurp+j)
                          xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                          do 750 i = 1, degc
                             xx (xdp + wm (i)) =
     $                       xx (xdp + wm (i)) +
     $                       xx (fxp + j*fdimc + wj (i))
750                       continue
760                    continue
 
c                   ----------------------------------------------------
                    else
c                   some columns already assembled out of this lson
c                   ----------------------------------------------------
 
                       do 780 j = 0, ludegr-1
                          col2 = ii (lurp+j)
                          if (col2 .gt. 0) then
                             xdp = ffxp + wic (col2)
cfpp$ nodepchk l
                             do 770 i = 1, degc
                                xx (xdp + wm (i)) =
     $                          xx (xdp + wm (i)) +
     $                          xx (fxp + j*fdimc + wj (i))
770                          continue
                          endif
780                    continue
                    endif
 
c                   done using wm (1..ludegc for offsets) ]
                    fleftc = wc (e) - w0
                    ii (ep+6) = fleftc
 
c                -------------------------------------------------------
                 else
c                this is an lson - assemble just one row
c                -------------------------------------------------------
 
                    xldp = ffxp + wir (row)
c                   get the offset, f, from the (e,f) tuple
                    f = ii (ip+1)
 
c                   ----------------------------------------------------
                    if (ludegr .eq. fleftr) then
c                   no columns assembled out this lson yet
c                   ----------------------------------------------------
 
cfpp$ nodepchk l
                       do 790 j = 0, ludegr-1
                          col2 = ii (lurp+j)
                          xx (xldp + wic (col2)) =
     $                    xx (xldp + wic (col2)) +
     $                    xx (fxp + j*fdimc + f)
790                    continue
 
c                   ----------------------------------------------------
                    else
c                   some columns already assembled out of this lson
c                   ----------------------------------------------------
 
cfpp$ nodepchk l
                       do 800 j = 0, ludegr-1
                          col2 = ii (lurp+j)
                          if (col2 .gt. 0) then
                             xx (xldp + wic (col2)) =
     $                       xx (xldp + wic (col2)) +
     $                       xx (fxp + j*fdimc + f)
                          endif
800                    continue
                    endif
 
                    ii (ep+6) = fleftc - 1
c                   flag the row as assembled from the lson
                    ii (lucp+f) = -row
                 endif
 
810           continue
 
c             ----------------------------------------------------------
c             update the count of (e,f) tuples in the element list
c             ----------------------------------------------------------
 
              wr (row) = wr (row) - deln
              ineed = ineed - 2*deln
           endif
 
c          -------------------------------------------------------------
c          assemble the original row and update count of entries
c          -------------------------------------------------------------
 
           rlen = wc (row)
           if (rlen .gt. 0) then
c             do not scan a very long row:
              if (rlen .le. rscan) then
                 rsiz = ii (pr)
                 ip = pr + rsiz - rlen
                 dlen = 0
cfpp$ nodepchk l
                 do 820 p = ip, ip + rlen - 1
                    col = ii (p)
                    if (wic (col) .ne. -2) then
c                      this entry can be assembled and deleted
c                      if wic (col) = -1, it is an older pivot col,
c                      otherwise (>=0) it is in the current element
                       dlen = dlen + 1
                       wm (dlen) = p
                    endif
820              continue
                 if (dlen .ne. 0) then
                    do 830 j = 1, dlen
c                      delete the entry
                       ii (wm (j)) = ii (ip+j-1)
830                 continue
                    rlen = rlen - dlen
                    ineed = ineed - dlen
                    wc (row) = rlen
                 endif
              endif
              rdeg = rdeg + rlen
           endif
 
c          -------------------------------------------------------------
c          compute the upper bound degree - excluding current front
c          -------------------------------------------------------------
 
           rdeg2 = ii (pr+1)
           rdeg = min (kleft1 - fflefr, rdeg2, rdeg)
           ii (pr+1) = rdeg
 
840     continue
 
c       ----------------------------------------------------------------
c       scan-4 wrap-up:  remove flags from assembled lsons
c       ----------------------------------------------------------------
 
c       while (lsons .ne. ndn+1) do
850     continue
        if (lsons .ne. ndn+1) then
           next = -wr (lsons)
           wr (lsons) = w0
           lsons = next
c       end while:
        goto 850
        endif
c       done un-flagging lsons, all now unflagged in wr (e) ]
 
c=======================================================================
c  degree update and numerical assemble is complete ]
c=======================================================================
 
c=======================================================================
c  factorize frontal matrix until next pivot extends it [
c=======================================================================
 
        do 1324 dummy4 = 1, n
c       (this loop is not indented due to its length)
 
c       ----------------------------------------------------------------
c       wc (e) = fextc+w0, where fextc is the external column
c               degree for each element (ep = rp (e)) appearing in
c               the element lists for each row in the pivot column.
c               if wc (e) < w0, then fextc is defined as ii (ep+6)
c
c       wr (e) = fextr+w0, where fextr is the external row
c               degree for each element (ep = rp (e)) appearing in
c               the element lists for each column in the pivot row
c               if wr (e) < w0, then fextr is defined as ii (ep+5)
c
c       wir (row) >= 0 for each row in pivot column pattern.
c               offset into pattern is given by:
c               wir (row) == offset - 1
c               wir (pivrow) is the offset of the latest pivot row
c
c       wic (col) >= 0 for each col in pivot row pattern.
c               wic (col) == (offset - 1) * ffdimc
c               wic (pivcol) is the offset of the latest pivot column
c
c       wpr (1..fflefr) is the pivot row pattern (excl pivot cols)
c       wpc (1..fflefc) is the pivot col pattern (excl pivot rows)
c       ----------------------------------------------------------------
 
c=======================================================================
c  divide pivot column by pivot
c=======================================================================
 
c       k-th pivot in frontal matrix located in c(ffdimc-k+1,ffdimr-k+1)
        xdp = ffxp + (ffdimr - k) * ffdimc
        x = xx (xdp + ffdimc - k)
 
c       divide c(1:fflefc,ffdimr-k+1) by pivot value
        x = one / x
        do 870 p = xdp, xdp + fflefc-1
           xx (p) = xx (p) * x
870     continue
c       count this as a call to the level-1 blas:
        rinfo (4) = rinfo (4) + fflefc
 
c=======================================================================
c  a pivot step is complete
c=======================================================================
 
        kleft = kleft - 1
        npiv = npiv + 1
        info (17) = info (17) + 1
 
c       ----------------------------------------------------------------
c       the pivot column is fully assembled and scaled, and is now the
c       (npiv)-th column of l. the pivot row is the (npiv)-th row of u.
c       ----------------------------------------------------------------
 
        wpr (n-npiv+1) = pivrow
        wpc (n-npiv+1) = pivcol
        wir (pivrow) = -1
        wic (pivcol) = -1
 
c       ----------------------------------------------------------------
c       deallocate the pivot row and pivot column
c       ----------------------------------------------------------------
 
        rlen = wc (pivrow)
        ineed = ineed - cscal - rscal - rlen
        pr = rp (pivrow)
        pc = cp (pivcol)
        ii (pr+1) = -1
        ii (pc+1) = -1
        rp (pivrow) = 0
        cp (pivcol) = 0
 
c=======================================================================
c  local search for next pivot within current frontal matrix [
c=======================================================================
 
        fedegc = fflefc
        fedegr = fflefr
        pfound = .false.
        okcol = fflefc .gt. 0
        okrow = .false.
 
c       ----------------------------------------------------------------
c       find column of minimum degree in current frontal row pattern
c       ----------------------------------------------------------------
 
c       among those columns of least (upper bound) degree, select the
c       column with the lowest column index
        if (okcol) then
           colpos = 0
           pivcol = n+1
           cdeg = n+1
c          can this be vectorized?  this is the most intensive
c          non-vector loop.
           do 880 j = 1, fflefr
              col = wpr (j)
              pc = cp (col)
              cdeg2 = ii (pc+1)
              better = cdeg2 .ge. 0 .and.
     $                (cdeg2 .lt. cdeg .or.
     $                (cdeg2 .eq. cdeg .and. col .lt. pivcol))
              if (better) then
                 cdeg = cdeg2
                 colpos = j
                 pivcol = col
              endif
880        continue
           okcol = colpos .ne. 0
        endif
 
c=======================================================================
c  assemble candidate pivot column in temporary workspace
c=======================================================================
 
        if (okcol) then
           pc = cp (pivcol)
           clen = ii (pc+6)
           okcol = fedegc + clen .le. ffdimc
        endif
 
        if (okcol) then
 
c          -------------------------------------------------------------
c          copy candidate column from current frontal matrix into
c          work vector xx (wxp ... wxp+ffdimc-1) [
c          -------------------------------------------------------------
 
           p = ffxp + (colpos - 1) * ffdimc - 1
cfpp$ nodepchk l
           do 890 i = 1, fflefc
              xx (wxp-1+i) = xx (p+i)
890        continue
 
c          -------------------------------------------------------------
c          update candidate column with previous pivots in this front
c          -------------------------------------------------------------
 
           if (k-k0 .gt. 0 .and. fflefc .ne. 0) then
              call sgemv ('n', fflefc, k-k0,
     $          -one, xx (ffxp + (ffdimr - k) * ffdimc)        ,ffdimc,
     $                xx (ffxp + (colpos - 1) * ffdimc + ffdimc - k), 1,
     $           one, xx (wxp)                                      , 1)
              rinfo (3) = rinfo (3) + 2*fflefc*(k-k0)
           endif
 
c          -------------------------------------------------------------
c          compute extended pivot column in xx (wxp..wxp-1+fedegc).
c          pattern of pivot column is placed in wpc (1..fedegc)
c          -------------------------------------------------------------
 
c          assemble the elements in the element list
           cep = (pc+9)
           celn = ii (pc+5)
           do 930 ip = cep, cep + 2*celn - 2, 2
              e = ii (ip)
              f = ii (ip+1)
              ep = rp (e)
              fleftc = ii (ep+6)
              fdimc = ii (ep+1)
              fxp = ii (ep+2)
              if (e .le. n) then
                 fluip = ii (ep)
                 lucp = (fluip + 7)
                 ludegc = ii (fluip+3)
              else
                 lucp = (ep+9)
                 ludegc = fdimc
              endif
              xp = fxp + f * fdimc
c             split into 3 loops so that they all vectorize on a cray
              f1 = fedegc
              do 900 p = lucp, lucp + ludegc - 1
                 row = ii (p)
                 if (row .gt. 0) then
                    if (wir (row) .lt. 0) then
                       f1 = f1 + 1
                       wpc (f1) = row
                    endif
                 endif
900           continue
              okcol = f1 + clen .le. ffdimc
              if (.not. okcol) then
c                exit out of loop if column too long:
                 go to 940
              endif
              do 910 i = fedegc+1, f1
                 row = wpc (i)
                 wir (row) = i - 1
                 xx (wxp-1+i) = zero
910           continue
              fedegc = f1
cfpp$ nodepchk l
              do 920 j = 0, ludegc - 1
                 row = ii (lucp+j)
                 if (row .gt. 0) then
                    xx (wxp + wir (row)) =
     $              xx (wxp + wir (row)) + xx (xp+j)
                 endif
920           continue
930        continue
c          loop exit label:
940        continue
        endif
 
c=======================================================================
c  find candidate pivot row - unless candidate pivot column is too long
c=======================================================================
 
        if (okcol) then
 
c          -------------------------------------------------------------
c          assemble the original entries in the column
c          -------------------------------------------------------------
 
           csiz = ii (pc)
           ip = pc + csiz - clen
           cxp = ii (pc+2)
cfpp$ nodepchk l
           do i = 0, clen - 1
              row = ii (ip+i)
              wir (row) = fedegc + i
              wpc (fedegc+1+i) = row
              xx  (wxp+fedegc+i) = xx (cxp+i)
           end do
           fedegc = fedegc + clen
 
c          -------------------------------------------------------------
c          update degree of candidate column - excluding current front
c          -------------------------------------------------------------
 
           cdeg = fedegc - fflefc
           ii (pc+1) = cdeg
 
c          -------------------------------------------------------------
c          find the maximum absolute value in the column
c          -------------------------------------------------------------
 
           maxval = abs (xx (wxp-1 + isamax (fedegc, xx (wxp), 1)))
           rinfo (3) = rinfo (3) + fedegc
           toler = relpt * maxval
           rdeg = n+1
 
c          -------------------------------------------------------------
c          look for the best possible pivot row in this column
c          -------------------------------------------------------------
 
           if (maxval .gt. zero) then
              if (symsrc) then
c                prefer symmetric pivots, if numerically acceptable
                 pivrow = pivcol
                 rowpos = wir (pivrow) + 1
                 if (rowpos .gt. 0 .and. rowpos .le. fflefc) then
c                   diagonal entry exists in the column pattern
c                   also within the current frontal matrix
                    x = abs (xx (wxp-1+rowpos))
                    if (x.ge.toler .and. x.gt.zero) then
c                      diagonal entry is numerically acceptable
                       pr = rp (pivrow)
                       rdeg = ii (pr+1)
                    endif
                 endif
              endif
              if (rdeg .eq. n+1) then
c                continue searching - no diagonal found or sought for.
c                minimize row degree subject to abs(value) constraints.
                 pivrow = n+1
                 do 960 i = 1, fflefc
                    row2 = wpc (i)
                    pr = rp (row2)
                    rdeg2 = ii (pr+1)
c                   among those numerically acceptable rows of least
c                   (upper bound) degree, select the row with the
c                   lowest row index
                    better = rdeg2 .lt. rdeg .or.
     $                      (rdeg2 .eq. rdeg .and. row2 .lt. pivrow)
                    if (better) then
                       x = abs (xx (wxp-1+i))
                       if (x.ge.toler .and. x.gt.zero) then
                          pivrow = row2
                          rdeg = rdeg2
                          rowpos = i
                       endif
                    endif
960              continue
              endif
           else
c             remove this column from any further pivot search
              cdeg = -(n+2)
              ii (pc+1) = cdeg
           endif
           okrow = rdeg .ne. n+1
        endif
 
c       done using xx (wxp...wxp+ffdimc-1) ]
 
c=======================================================================
c  if found, construct candidate pivot row pattern
c=======================================================================
 
        if (okrow) then
 
c          -------------------------------------------------------------
c          assemble the elements in the element list
c          -------------------------------------------------------------
 
           pr = rp (pivrow)
           rep = (pr+2)
           reln = wr (pivrow)
           do 990 ip = rep, rep + 2*reln - 2, 2
              e = ii (ip)
              ep = rp (e)
              if (e .le. n) then
                 fluip = ii (ep)
                 lucp = (fluip + 7)
                 ludegr = ii (fluip+2)
                 ludegc = ii (fluip+3)
                 lurp = lucp + ludegc
                 fleftr = ii (ep+5)
                 okrow = fleftr .le. ffdimr
                 if (.not. okrow) then
c                   exit out of loop if row too long:
                    go to 1000
                 endif
c                split into two loops so that both vectorize on a cray
                 f1 = fedegr
                 do 970 p = lurp, lurp + ludegr - 1
                    col = ii (p)
                    if (col .gt. 0) then
                       if (wic (col) .eq. -2) then
                          f1 = f1 + 1
                          wpr (f1) = col
                       endif
                    endif
970              continue
                 okrow = f1 .le. ffdimr
                 if (.not. okrow) then
c                   exit out of loop if row too long:
                    go to 1000
                 endif
                 do 980 i = fedegr+1, f1
                    wic (wpr (i)) = (i - 1) * ffdimc
980              continue
                 fedegr = f1
              else
c                this is an artificial element (a dense column)
                 lurp = (ep+8)
                 col = ii (lurp)
                 if (wic (col) .eq. -2) then
                    wic (col) = fedegr * ffdimc
                    fedegr = fedegr + 1
                    wpr (fedegr) = col
                    okrow = fedegr .le. ffdimr
                    if (.not. okrow) then
c                      exit out of loop if row too long:
                       go to 1000
                    endif
                 endif
              endif
990        continue
c          loop exit label:
1000       continue
        endif
 
        if (okrow) then
 
c          -------------------------------------------------------------
c          assemble the original entries in the row
c          -------------------------------------------------------------
 
           rlen = wc (pivrow)
           if (rlen .gt. 0) then
              f1 = fedegr
              rsiz = ii (pr)
              p2 = pr + rsiz
c             split into two loops so that they both vectorize on a cray
              do 1010 p = p2 - rlen, p2 - 1
                 col = ii (p)
                 if (wic (col) .eq. -2) then
c                   this entry cannot be assembled, do not delete
                    f1 = f1 + 1
                    wpr (f1) = col
                 endif
1010          continue
              rlen2 = f1 - fedegr
              if (rlen2 .lt. rlen) then
c                delete one or more entries in the row
                 do 1020 i = fedegr+1, f1
                    ii (p2 - f1 + i - 1) = wpr (i)
1020             continue
                 ineed = ineed - (rlen - rlen2)
                 wc (pivrow) = rlen2
              endif
 
c             ----------------------------------------------------------
c             update the candidate row degree - excluding current front
c             ----------------------------------------------------------
 
              rdeg = f1 - fflefr
              ii (pr+1) = rdeg
 
c             ----------------------------------------------------------
c             pivot is found if candidate pivot row is not too long
c             ----------------------------------------------------------
 
              okrow = f1 .le. ffdimr
              if (okrow) then
                 do 1030 i = fedegr+1, f1
                    wic (wpr (i)) = (i - 1) * ffdimc
1030             continue
                 fedegr = f1
              endif
 
           else
 
c             ----------------------------------------------------------
c             update the candidate row degree - excluding current front
c             ----------------------------------------------------------
 
              rdeg = fedegr - fflefr
              ii (pr+1) = rdeg
           endif
        endif
 
c       ----------------------------------------------------------------
c       if pivot not found: clear wir and wic
c       ----------------------------------------------------------------
 
        pfound = okrow .and. okcol
        if (.not. pfound) then
           movelu = k .gt. 0
           do 1040 i = fflefr+1, fedegr
              wic (wpr (i)) = -2
1040       continue
           fedegr = fflefr
           do 1050 i = fflefc+1, fedegc
              wir (wpc (i)) = -1
1050       continue
           fedegc = fflefc
        else
           movelu = fedegc .gt. ffdimc - k .or. fedegr .gt. ffdimr - k
        endif
 
c       ----------------------------------------------------------------
c       wpr (1..fflefr)                 unextended pivot row pattern
c       wpr (fflefr+1 .. fedegr)        extended pattern, if pfound
c       wpr (fedegr+1 .. n-npiv)        empty space
c       wpr (n-npiv+1 .. n)             pivot row order
c
c       wpc (1..fflefc)                 unextended pivot column pattern
c       wpc (fflefc+1 .. fedegc)        extended pattern, if pfound
c       wpc (fedegc+1 .. n-npiv)        empty space
c       wpc (n-npiv+1 .. n)             pivot column order
c       ----------------------------------------------------------------
 
c=======================================================================
c  local pivot search complete ]
c=======================================================================
 
c=======================================================================
c  update contribution block: rank-nb, or if lu arrowhead to be moved
c=======================================================================
 
        if (k-k0 .ge. nb .or. movelu) then
           call sgemm ('n', 'n', fflefc, fflefr, k-k0,
     $          -one, xx (ffxp + (ffdimr - k) * ffdimc), ffdimc,
     $                xx (ffxp +  ffdimc - k)          , ffdimc,
     $           one, xx (ffxp)                        , ffdimc)
           rinfo (6) = rinfo (6) + 2*fflefc*fflefr*(k-k0)
           k0 = k
        endif
 
c=======================================================================
c  move the lu arrowhead if no pivot found, or pivot needs room
c=======================================================================
 
        if (movelu) then
 
c          allocate permanent space for the lu arrowhead
           ludegr = fflefr
           ludegc = fflefc
           xs = k*ludegc + k*ludegr + k*k
           is = 7 + ludegc + ludegr + nsons
           if (is .gt. itail-ihead .or. xs .gt. xtail-xhead) then
              if (is .gt. itail-ihead) then
c                garbage collection because we ran out of integer mem
                 info (14) = info (14) + 1
              endif
              if (xs .gt. xtail-xhead) then
c                garbage collection because we ran out of real mem
                 info (15) = info (15) + 1
              endif
              call ums2fg (xx, xsize, xhead, xtail, xuse,
     $                     ii, isize, ihead, itail, iuse,
     $                     cp, rp, dn, n, icntl, wir, wic, wr, wc,
     $                     ffxp, ffsize, wxp, ffdimc, .false.,
     $                     pfree, xfree, mhead, mtail, slots)
c             at this point, iuse = ineed and xuse = xneed
           endif
 
           itail = itail - is
           luip = itail
           iuse = iuse + is
           ineed = ineed + is
           xtail = xtail - xs
           luxp = xtail
           xuse = xuse + xs
           xneed = xneed + xs
           info (18) = max (info (18), iuse)
           info (19) = max (info (19), ineed)
           info (20) = max (info (20), xuse)
           info (21) = max (info (21), xneed)
           if (ihead .gt. itail .or. xhead .gt. xtail) then
c             error return, if not enough integer and/or real memory:
              go to 9000
           endif
 
c          -------------------------------------------------------------
c          get memory usage for next call to ums2rf
c          -------------------------------------------------------------
 
           xruse = xruse + xs
           xrmax = max (xrmax, xruse)
 
c          -------------------------------------------------------------
c          save the new lu arrowhead
c          -------------------------------------------------------------
 
c          save the scalar data of the lu arrowhead
           ii (luip) = luxp
           ii (luip+1) = k
           ii (luip+2) = ludegr
           ii (luip+3) = ludegc
           ii (luip+4) = nsons
           ii (luip+5) = 0
           ii (luip+6) = 0
           e = ffrow
           if (e .eq. e1) then
c             this is the first lu arrowhead from this global pivot
              luip1 = luip
           endif
           wr (e) = -(ndn+2)
           wc (e) = -(ndn+2)
 
c          save column pattern
           lucp = (luip + 7)
           do 1060 i = 0, ludegc-1
              ii (lucp+i) = wpc (i+1)
1060       continue
 
c          save row pattern
           lurp = lucp + ludegc
           do 1070 i = 0, ludegr-1
              ii (lurp+i) = wpr (i+1)
1070       continue
 
c          add list of sons after the end of the frontal matrix pattern
c          this list of sons is for the refactorization (ums2rf) only.
           lusonp = lurp + ludegr
           ip = lusonp
           e = sonlst
c          while (e > 0) do
1080       continue
           if (e .gt. 0) then
              ep = rp (e)
              if (wc (e) .eq. -(ndn+2)) then
c                luson
                 ii (ip) = e
              else if (wc (e) .eq. w0) then
c                uson
                 ii (ip) = e + n
              else if (wr (e) .eq. w0) then
c                lson
                 ii (ip) = e + 2*n
              endif
              next = wir (e) + n + 2
              wir (e) = -1
              e = next
              ip = ip + 1
c          end while:
           goto 1080
           endif
           nsons = 0
           sonlst = 0
 
c          move the l1,u1 matrix, compressing the dimension from
c          ffdimc to ldimc.  the lu arrowhead grows on top of stack.
           ldimc = k + ludegc
           xp = ffxp + (ffdimr-1)*ffdimc + ffdimc-1
           do 1100 j = 0, k-1
cfpp$ nodepchk l
              do 1090 i = 0, k-1
                 xx (luxp + j*ldimc + i) = xx (xp - j*ffdimc - i)
1090          continue
1100       continue
 
c          move l2 matrix, compressing dimension from ffdimc to ludegc+k
           if (ludegc .ne. 0) then
              lxp = luxp + k
              xp = ffxp + (ffdimr-1)*ffdimc
              do 1120 j = 0, k-1
cfpp$ nodepchk l
                 do 1110 i = 0, ludegc-1
                    xx (lxp + j*ldimc + i) = xx (xp - j*ffdimc + i)
1110             continue
1120          continue
           endif
 
c          move the u2 block.
           if (ludegr .ne. 0) then
              uxp = luxp + k * ldimc
              xp = ffxp + ffdimc-1
              do 1140 j = 0, ludegr-1
cfpp$ nodepchk l
                 do 1130 i = 0, k-1
                    xx (uxp + j*k + i) = xx (xp + j*ffdimc - i)
1130             continue
1140          continue
           endif
 
c          one more lu arrowhead has been created
           nlu = nlu + 1
           nzu = (k*(k-1)/2) + k*ludegc
           nzl = (k*(k-1)/2) + k*ludegr
           info (10) = info (10) + nzl
           info (11) = info (11) + nzu
 
c          no more rows of u or columns of l in current frontal array
           k = 0
           k0 = 0
 
           if (pfound) then
 
c             ----------------------------------------------------------
c             place the old frontal matrix as the only item in the son
c             list, since the next "implied" frontal matrix will have
c             this as its son.
c             ----------------------------------------------------------
 
              nsons = 1
              e = ffrow
              wir (e) = - n - 2
              sonlst = e
 
c             ----------------------------------------------------------
c             the contribution block of the old frontal matrix is still
c             stored in the current frontal matrix, and continues (in a
c             unifrontal sense) as a "new" frontal matrix (same array
c             but with a new name, and the lu arrowhead is removed and
c             placed in the lu factors).  old name is "ffrow", new name
c             is "pivrow".
c             ----------------------------------------------------------
 
              rp (e) = luip
              ffrow = pivrow
           endif
        endif
 
c=======================================================================
c  stop the factorization of this frontal matrix if no pivot found
c=======================================================================
 
c       (this is the only way out of loop 1395)
        if (.not. pfound) then
c          exit out of loop 1395 if pivot not found:
           go to 1400
        endif
 
c=======================================================================
c  update the pivot column, and move into position as (k+1)-st col of l
c=======================================================================
 
        xsp = (colpos - 1) * ffdimc
        xdp = (ffdimr - k - 1) * ffdimc
        fsp = ffxp + xsp
        fdp = ffxp + xdp
 
        if (k-k0 .gt. 0 .and. fflefc .ne. 0) then
           call sgemv ('n', fflefc, k-k0,
     $          -one, xx (fdp + ffdimc    ), ffdimc,
     $                xx (fsp + ffdimc - k), 1,
     $           one, xx (fsp             ), 1)
           rinfo (5) = rinfo (5) + 2*fflefc*(k-k0)
        endif
 
        if (fflefr .lt. ffdimr - k) then
 
           xlp = (fflefr - 1) * ffdimc
           if (fflefr .eq. colpos) then
 
c             ----------------------------------------------------------
c             move c(:,colpos) => c(:,ffdimr-k)
c             ----------------------------------------------------------
 
              if (ffdimc .le. 64) then
c                copy the gap (useless work but quicker if gap is small)
cdir$ shortloop
cfpp$ nodepchk l
                 do 1150 i = 0, ffdimc - 1
                    xx (fdp+i) = xx (fsp+i)
1150             continue
              else
c                copy only what needs to be copied
c                column of the contribution block:
cfpp$ nodepchk l
                 do 1160 i = 0, fflefc - 1
                    xx (fdp+i) = xx (fsp+i)
1160             continue
c                column of the u2 block
cfpp$ nodepchk l
                 do 1170 i = ffdimc - k, ffdimc - 1
                    xx (fdp+i) = xx (fsp+i)
1170             continue
              endif
 
           else
 
c             ----------------------------------------------------------
c             move c(:,colpos) => c(:,ffdimr-k)
c             move c(:,fflefr) => c(:,colpos)
c             ----------------------------------------------------------
 
              flp = ffxp + xlp
              if (ffdimc .le. 64) then
c                copy the gap (useless work but quicker if gap is small)
cdir$ shortloop
cfpp$ nodepchk l
                 do 1180 i = 0, ffdimc - 1
                    xx (fdp+i) = xx (fsp+i)
                    xx (fsp+i) = xx (flp+i)
1180             continue
              else
 
c                copy only what needs to be copied
c                columns of the contribution block:
cfpp$ nodepchk l
                 do 1190 i = 0, fflefc - 1
                    xx (fdp+i) = xx (fsp+i)
                    xx (fsp+i) = xx (flp+i)
1190             continue
c                columns of the u2 block:
cfpp$ nodepchk l
                 do 1200 i = ffdimc - k, ffdimc - 1
                    xx (fdp+i) = xx (fsp+i)
                    xx (fsp+i) = xx (flp+i)
1200             continue
              endif
 
              swpcol = wpr (fflefr)
              wpr (colpos) = swpcol
              wic (swpcol) = xsp
           endif
 
           if (fedegr .ne. fflefr) then
c             move column fedegr to column fflefr (pattern only)
              swpcol = wpr (fedegr)
              wpr (fflefr) = swpcol
              wic (swpcol) = xlp
           endif
 
        else if (colpos .ne. ffdimr - k) then
 
c          -------------------------------------------------------------
c          swap c(:,colpos) <=> c (:,ffdimr-k)
c          -------------------------------------------------------------
 
           if (ffdimc .le. 64) then
c             swap the gap (useless work but quicker if gap is small)
cdir$ shortloop
cfpp$ nodepchk l
cfpp$ nolstval l
              do 1210 i = 0, ffdimc - 1
                 x = xx (fdp+i)
                 xx (fdp+i) = xx (fsp+i)
                 xx (fsp+i) = x
1210          continue
           else
c             swap only what needs to be swapped
c             columns of the contribution block:
cfpp$ nodepchk l
cfpp$ nolstval l
              do 1220 i = 0, fflefc - 1
                 x = xx (fdp+i)
                 xx (fdp+i) = xx (fsp+i)
                 xx (fsp+i) = x
1220          continue
c             columns of the u2 block:
cfpp$ nodepchk l
cfpp$ nolstval l
              do 1230 i = ffdimc - k, ffdimc - 1
                 x = xx (fdp+i)
                 xx (fdp+i) = xx (fsp+i)
                 xx (fsp+i) = x
1230          continue
           endif
           swpcol = wpr (ffdimr - k)
           wpr (colpos) = swpcol
           wic (swpcol) = xsp
        endif
 
        wic (pivcol) = xdp
        fedegr = fedegr - 1
        scan2 = fflefr
        fflefr = fflefr - 1
 
c=======================================================================
c  move pivot row into position as (k+1)-st row of u, and update
c=======================================================================
 
        xsp = rowpos - 1
        xdp = ffdimc - k - 1
        fsp = ffxp + xsp
        fdp = ffxp + xdp
 
        if (fflefc .lt. ffdimc - k) then
 
           xlp = fflefc - 1
           if (fflefc .eq. rowpos) then
 
c             ----------------------------------------------------------
c             move c(rowpos,:) => c(ffdimc-k,:)
c             ----------------------------------------------------------
 
              if (ffdimr .le. 64) then
c                copy the gap (useless work but quicker if gap is small)
cdir$ shortloop
cfpp$ nodepchk l
                 do 1240 j = 0, (ffdimr - 1) * ffdimc, ffdimc
                    xx (fdp+j) = xx (fsp+j)
1240             continue
              else
c                copy only what needs to be copied
c                row of the contribution block:
cfpp$ nodepchk l
                 do 1250 j = 0, (fflefr - 1) * ffdimc, ffdimc
                    xx (fdp+j) = xx (fsp+j)
1250             continue
c                row of the l2 block:
cfpp$ nodepchk l
                 do 1260 j = (ffdimr - k - 1) * ffdimc,
     $                       (ffdimr - 1) * ffdimc, ffdimc
                    xx (fdp+j) = xx (fsp+j)
1260             continue
              endif
 
           else
 
c             ----------------------------------------------------------
c             move c(rowpos,:) => c(ffdimc-k,:)
c             move c(fflefc,:) => c(rowpos,:)
c             ----------------------------------------------------------
 
              flp = ffxp + xlp
              if (ffdimr .le. 64) then
c                copy the gap (useless work but quicker if gap is small)
cdir$ shortloop
cfpp$ nodepchk l
                 do 1270 j = 0, (ffdimr - 1) * ffdimc, ffdimc
                    xx (fdp+j) = xx (fsp+j)
                    xx (fsp+j) = xx (flp+j)
1270             continue
              else
c                copy only what needs to be copied
c                rows of the contribution block:
cfpp$ nodepchk l
                 do 1280 j = 0, (fflefr - 1) * ffdimc, ffdimc
                    xx (fdp+j) = xx (fsp+j)
                    xx (fsp+j) = xx (flp+j)
1280             continue
c                rows of the l2 block:
cfpp$ nodepchk l
                 do 1290 j = (ffdimr - k - 1) * ffdimc,
     $                       (ffdimr - 1) * ffdimc, ffdimc
                    xx (fdp+j) = xx (fsp+j)
                    xx (fsp+j) = xx (flp+j)
1290             continue
              endif
              swprow = wpc (fflefc)
              wpc (rowpos) = swprow
              wir (swprow) = xsp
           endif
 
           if (fedegc .ne. fflefc) then
c             move row fedegc to row fflefc (pattern only)
              swprow = wpc (fedegc)
              wpc (fflefc) = swprow
              wir (swprow) = xlp
           endif
 
        else if (rowpos .ne. ffdimc - k) then
 
c          -------------------------------------------------------------
c          swap c(rowpos,:) <=> c (ffdimc-k,:)
c          -------------------------------------------------------------
 
           if (ffdimr .le. 64) then
c             swap the gap (useless work but quicker if gap is small)
cdir$ shortloop
cfpp$ nodepchk l
cfpp$ nolstval l
              do 1300 j = 0, (ffdimr - 1) * ffdimc, ffdimc
                 x = xx (fdp+j)
                 xx (fdp+j) = xx (fsp+j)
                 xx (fsp+j) = x
1300          continue
           else
c             swap only what needs to be swapped
c             rows of the contribution block:
cfpp$ nodepchk l
cfpp$ nolstval l
              do 1310 j = 0, (fflefr - 1) * ffdimc, ffdimc
                 x = xx (fdp+j)
                 xx (fdp+j) = xx (fsp+j)
                 xx (fsp+j) = x
1310          continue
c             rows of the l2 block:
cfpp$ nodepchk l
cfpp$ nolstval l
              do 1320 j = (ffdimr - k - 1) * ffdimc,
     $                    (ffdimr - 1) * ffdimc, ffdimc
                 x = xx (fdp+j)
                 xx (fdp+j) = xx (fsp+j)
                 xx (fsp+j) = x
 
1320          continue
           endif
           swprow = wpc (ffdimc - k)
           wpc (rowpos) = swprow
           wir (swprow) = xsp
        endif
 
        wir (pivrow) = xdp
        fedegc = fedegc - 1
        scan1 = fflefc
        fflefc = fflefc - 1
 
        if (k-k0 .gt. 0 .and. fflefr .gt. 0) then
           call sgemv ('t', k-k0, fflefr,
     $       -one, xx (fdp + 1)                    , ffdimc,
     $             xx (fdp + (ffdimr - k) * ffdimc), ffdimc,
     $        one, xx (fdp)                        , ffdimc)
           rinfo (5) = rinfo (5) + 2*(k-k0)*fflefr
        endif
 
c=======================================================================
c  prepare for degree update and next local pivot search
c=======================================================================
 
c       ----------------------------------------------------------------
c       if only column pattern has been extended:
c               scan1:  new rows only
c               scan2:  no columns scanned
c               scan3:  all columns
c               scan4:  new rows only
c
c       if only row pattern has been extended:
c               scan1:  no rows scanned
c               scan2:  new columns only
c               scan3:  new columns only
c               scan4:  all rows
c
c       if both row and column pattern have been extended:
c               scan1:  new rows only
c               scan2:  new columns only
c               scan3:  all columns
c               scan4:  all rows
c
c       if no patterns have been extended:
c               scan1-4: none
c       ----------------------------------------------------------------
 
        if (fedegc .eq. fflefc) then
c          column pattern has not been extended
           scan3 = fflefr + 1
        else
c          column pattern has been extended.
           scan3 = 0
        endif
 
        if (fedegr .eq. fflefr) then
c          row pattern has not been extended
           scan4 = fflefc + 1
        else
c          row pattern has been extended
           scan4 = 0
        endif
 
c=======================================================================
c  finished with step k (except for assembly and scaling of pivot col)
c=======================================================================
 
        k = k + 1
 
c       ----------------------------------------------------------------
c       exit loop if frontal matrix has been extended
c       ----------------------------------------------------------------
 
        if (fedegr .ne. fflefr .or. fedegc .ne. fflefc) then
           go to 1325
        endif
 
1324    continue
c       exit label for loop 1324:
1325    continue
 
c=======================================================================
c  finished factorizing while frontal matrix is not extended ]
c=======================================================================
 
c=======================================================================
c  extend the frontal matrix [
c=======================================================================
 
c       ----------------------------------------------------------------
c       zero the newly extended frontal matrix
c       ----------------------------------------------------------------
 
c       fill-in due to amalgamation caused by this step is
c       k*(fedegr-fflefr+fedegc-fflefc)
 
        do 1350 j = fflefr, fedegr - 1
c          zero the new columns in the contribution block:
           do 1330 i = 0, fedegc - 1
              xx (ffxp + j*ffdimc + i) = zero
1330       continue
c          zero the new columns in u block:
           do 1340 i = ffdimc - k, ffdimc - 1
              xx (ffxp + j*ffdimc + i) = zero
1340       continue
1350    continue
 
cfpp$ nodepchk l
        do 1380 i = fflefc, fedegc - 1
c          zero the new rows in the contribution block:
cfpp$ nodepchk l
           do 1360 j = 0, fflefr - 1
              xx (ffxp + j*ffdimc + i) = zero
1360       continue
c          zero the new rows in l block:
cfpp$ nodepchk l
           do 1370 j = ffdimr - k, ffdimr - 1
              xx (ffxp + j*ffdimc + i) = zero
1370       continue
1380    continue
 
c       ----------------------------------------------------------------
c       remove the new columns from the degree lists
c       ----------------------------------------------------------------
 
        do 1390 j = fflefr+1, fedegr
           pc = cp (wpr (j))
           cdeg = ii (pc+1)
           if (cdeg .gt. 0) then
              cnxt = ii (pc+7)
              cprv = ii (pc+8)
              if (cnxt .ne. 0) then
                 ii (cp (cnxt)+8) = cprv
              endif
              if (cprv .ne. 0) then
                 ii (cp (cprv)+7) = cnxt
              else
                 head (cdeg) = cnxt
              endif
           endif
1390    continue
 
c       ----------------------------------------------------------------
c       finalize extended row and column pattern of the frontal matrix
c       ----------------------------------------------------------------
 
        fflefc = fedegc
        fflefr = fedegr
        fmaxr = max (fmaxr, fflefr + k)
        fmaxc = max (fmaxc, fflefc + k)
 
c=======================================================================
c  done extending the current frontal matrix ]
c=======================================================================
 
1395    continue
c       exit label for loop 1395:
1400    continue
 
c=======================================================================
c  done assembling and factorizing the current frontal matrix ]
c=======================================================================
 
c=======================================================================
c  wrap-up:  complete the current frontal matrix [
c=======================================================================
 
c       ----------------------------------------------------------------
c       store the maximum front size in the first lu arrowhead
c       ----------------------------------------------------------------
 
        ii (luip1+5) = fmaxr
        ii (luip1+6) = fmaxc
 
c       one more frontal matrix is finished
        info (13) = info (13) + 1
 
c       ----------------------------------------------------------------
c       add the current frontal matrix to the degrees of each column,
c       and place the modified columns back in the degree lists
c       ----------------------------------------------------------------
 
c       do so in reverse order to try to improve pivot tie-breaking
        do 1410 j = fflefr, 1, -1
           col = wpr (j)
           pc = cp (col)
c          add the current frontal matrix to the degree
           cdeg = ii (pc+1)
           cdeg = min (kleft, cdeg + fflefc)
           if (cdeg .gt. 0) then
              ii (pc+1) = cdeg
              cnxt = head (cdeg)
              ii (pc+7) = cnxt
              ii (pc+8) = 0
              if (cnxt .ne. 0) then
                 ii (cp (cnxt)+8) = col
              endif
              head (cdeg) = col
              mindeg = min (mindeg, cdeg)
           endif
1410    continue
 
c       ----------------------------------------------------------------
c       add the current frontal matrix to the degrees of each row
c       ----------------------------------------------------------------
 
cfpp$ nodepchk l
        do 1420 i = 1, fflefc
           row = wpc (i)
           pr = rp (row)
           rdeg = ii (pr+1)
           rdeg = min (kleft, rdeg + fflefr)
           ii (pr+1) = rdeg
1420    continue
 
c       ----------------------------------------------------------------
c       reset w0 so that wr (1..n) < w0 and wc (1..n) < w0.
c       also ensure that w0 + n would not cause integer overflow
c       ----------------------------------------------------------------
 
        w0 = w0 + fmax + 1
        if (w0 .ge. w0big) then
           w0 = ndn+2
           do 1430 e = 1, n+dn
              if (wr (e) .gt. ndn) then
c                this is a frontal matrix
                 wr (e) = w0-1
                 wc (e) = w0-1
              endif
1430       continue
        endif
 
c       ----------------------------------------------------------------
c       deallocate work vector
c       ----------------------------------------------------------------
 
        xuse = xuse - ffdimc
        xneed = xneed - ffdimc
        xhead = xhead - ffdimc
 
c       ----------------------------------------------------------------
c       get the name of this new frontal matrix, and size of
c       contribution block
c       ----------------------------------------------------------------
 
        e = ffrow
        xs = fflefr * fflefc
        fmax = max (fmax, fflefr, fflefc)
 
c       ----------------------------------------------------------------
c       get memory usage for next call to ums2rf
c       ----------------------------------------------------------------
 
        xruse = xruse - ffsize + xs
 
c       ----------------------------------------------------------------
c       if contribution block empty, deallocate and continue next step
c       ----------------------------------------------------------------
 
        if (fflefr .le. 0 .or. fflefc .le. 0) then
           rp (e) = luip
           xuse = xuse - ffsize
           xneed = xneed - ffsize
           xhead = ffxp
           do 1440 i = 1, fflefr
              wic (wpr (i)) = -2
1440       continue
           do 1450 i = 1, fflefc
              wir (wpc (i)) = -1
1450       continue
c          next iteration of main factorization loop 1540:
           goto 1540
        endif
 
c       ----------------------------------------------------------------
c       prepare the contribution block for later assembly
c       ----------------------------------------------------------------
 
        if (fscal .gt. itail-ihead) then
           info (14) = info (14) + 1
           call ums2fg (xx, xsize, xhead, xtail, xuse,
     $                  ii, isize, ihead, itail, iuse,
     $                  cp, rp, dn, n, icntl, wir, wic, wr, wc,
     $                  ffxp, ffsize, 0, 0, .false.,
     $                  pfree, xfree, mhead, mtail, slots)
c          at this point, iuse = ineed and xuse = xneed
        endif
 
        ep = ihead
        ihead = ihead + fscal
        iuse = iuse + fscal
        ineed = ineed + fscal
        info (18) = max (info (18), iuse)
        info (19) = max (info (19), ineed)
        if (ihead .gt. itail) then
c          error return, if not enough integer memory:
c          (highly unlikely to run out of memory at this point)
           go to 9000
        endif
 
        rp (e) = ep
        ii (ep) = luip
        ii (ep+5) = fflefr
        ii (ep+6) = fflefc
        wr (e) = w0-1
        wc (e) = w0-1
 
c       count the numerical assembly
        rinfo (2) = rinfo (2) + xs
 
        if (xs .le. xfree) then
 
c          -------------------------------------------------------------
c          compress and store the contribution block in a freed block
c          -------------------------------------------------------------
 
c          place the new block in the list in front of the free block
           xdp = ii (pfree+2)
           ii (pfree+2) = ii (pfree+2) + xs
           xfree = xfree - xs
           mprev = ii (pfree+4)
           if (xfree .eq. 0) then
c             delete the free block if its size is zero
              mnext = ii (pfree+3)
              pfree = 0
              xfree = -1
           else
              mnext = pfree
           endif
           if (mnext .ne. 0) then
              ii (mnext+4) = ep
           else
              mtail = ep
           endif
           if (mprev .ne. 0) then
              ii (mprev+3) = ep
           else
              mhead = ep
           endif
           do 1470 j = 0, fflefr - 1
cfpp$ nodepchk l
              do 1460 i = 0, fflefc - 1
                 xx (xdp + j*fflefc + i) = xx (ffxp + j*ffdimc + i)
1460          continue
1470       continue
           xhead = ffxp
           xuse = xuse - ffsize
           xneed = xneed - ffsize + xs
           ffdimc = fflefc
           ii (ep+1) = ffdimc
           ii (ep+2) = xdp
           ii (ep+3) = mnext
           ii (ep+4) = mprev
 
        else
 
c          -------------------------------------------------------------
c          deallocate part of the unused portion of the frontal matrix
c          -------------------------------------------------------------
 
c          leave the contribution block c (1..fflefc, 1..fflefr) at the
c          head of xx, with column dimension of ffdimc and in space
c          of size (fflefr-1)*ffdimc for the first fflefr columns, and
c          fflefc for the last column.
           xneed = xneed - ffsize + xs
           xs = ffsize - (fflefc + (fflefr-1)*ffdimc)
           xhead = xhead - xs
           xuse = xuse - xs
           ii (ep+1) = ffdimc
           ii (ep+2) = ffxp
           ii (ep+3) = 0
           ii (ep+4) = mtail
           if (mtail .eq. 0) then
              mhead = ep
           else
              ii (mtail+3) = ep
           endif
           mtail = ep
        endif
 
c       ----------------------------------------------------------------
c       add tuples to the amount of integer space needed - and add
c       limit+cscal to maximum need to account for worst-case possible
c       reallocation of rows/columns.  required integer memory usage
c       is guaranteed not to exceed iworst during the placement of (e,f)
c       tuples in the two loops below.
c       ----------------------------------------------------------------
 
        ineed = ineed + 2*(fflefr+fflefc)
        iworst = ineed + limit + cscal
        info (19) = max (info (19), iworst)
        info (18) = max (info (18), iworst)
 
c       ----------------------------------------------------------------
c       place (e,f) in the element list of each column
c       ----------------------------------------------------------------
 
        do 1500 i = 1, fflefr
           col = wpr (i)
           pc = cp (col)
           celn = ii (pc+5)
           csiz = ii (pc)
           clen = ii (pc+6)
c          clear the column offset
           wic (col) = -2
 
c          -------------------------------------------------------------
c          make sure an empty slot exists - if not, create one
c          -------------------------------------------------------------
 
           if (2*(celn+1) + clen + cscal .gt. csiz) then
 
c             ----------------------------------------------------------
c             no room exists - reallocate elsewhere
c             ----------------------------------------------------------
 
c             at least this much space is needed:
              is = 2 * (celn + 1) + clen
c             add some slots for growth: at least 8 tuples,
c             or double the size - whichever is larger (but with a total
c             size not larger than limit+cscal)
              is = min (is + max (16, is), limit)
              csiz2 = is + cscal
 
c             ----------------------------------------------------------
c             make sure enough room exists: garbage collection if needed
c             ----------------------------------------------------------
 
              if (csiz2 .gt. itail-ihead) then
c                garbage collection:
                 info (14) = info (14) + 1
                 call ums2fg (xx, xsize, xhead, xtail, xuse,
     $                        ii, isize, ihead, itail, iuse,
     $                        cp, rp, dn, n, icntl, wir, wic, wr, wc,
     $                        0, 0, 0, 0, .true.,
     $                        pfree, xfree, mhead, mtail, slots)
c                at this point, iuse+csiz2 <= iworst and xuse = xneed
                 pc = cp (col)
                 csiz = ii (pc)
              endif
 
c             ----------------------------------------------------------
c             get space for the new copy
c             ----------------------------------------------------------
 
              pc2 = ihead
              ihead = ihead + csiz2
              iuse = iuse + csiz2
              info (18) = max (info (18), iuse)
              if (ihead .gt. itail) then
c                error return, if not enough integer memory:
                 go to 9000
              endif
 
c             ----------------------------------------------------------
c             make the copy, leaving hole in middle for element list
c             ----------------------------------------------------------
 
c             copy the cscal scalars, and the element list
cfpp$ nodepchk l
              do 1480 j = 0, cscal + 2*celn - 1
                 ii (pc2+j) = ii (pc+j)
1480          continue
 
c             copy column indices of original entries (xx is unchanged)
cfpp$ nodepchk l
              do 1490 j = 0, clen - 1
                 ii (pc2+csiz2-clen+j) = ii (pc+csiz-clen+j)
1490          continue
 
              if (clen .gt. 0) then
c                place the new block in the memory-list
                 mnext = ii (pc2+3)
                 mprev = ii (pc2+4)
                 if (mnext .ne. 0) then
                    ii (mnext+4) = pc2
                 else
                    mtail = pc2
                 endif
                 if (mprev .ne. 0) then
                    ii (mprev+3) = pc2
                 else
                    mhead = pc2
                 endif
              endif
 
              cp (col) = pc2
              ii (pc2) = csiz2
 
c             ----------------------------------------------------------
c             deallocate the old copy of the column in ii (not in xx)
c             ----------------------------------------------------------
 
              ii (pc+1) = -1
              ii (pc+6) = 0
              pc = pc2
           endif
 
c          -------------------------------------------------------------
c          place the new (e,f) tuple in the element list of the column
c          -------------------------------------------------------------
 
           cep = (pc+9)
           ii (cep + 2*celn  ) = e
           ii (cep + 2*celn+1) = i - 1
           ii (pc+5) = celn + 1
1500    continue
 
c       ----------------------------------------------------------------
c       place (e,f) in the element list of each row
c       ----------------------------------------------------------------
 
        do 1530 i = 1, fflefc
           row = wpc (i)
           pr = rp (row)
           rsiz = ii (pr)
           reln = wr (row)
           rlen = wc (row)
c          clear the row offset
           wir (row) = -1
 
c          -------------------------------------------------------------
c          make sure an empty slot exists - if not, create one
c          -------------------------------------------------------------
 
           if (2*(reln+1) + rlen + rscal .gt. rsiz) then
 
c             ----------------------------------------------------------
c             no room exists - reallocate elsewhere
c             ----------------------------------------------------------
 
c             at least this much space is needed:
              is = 2 * (reln + 1) + rlen
c             add some extra slots for growth - for at least 8
c             tuples, or double the size (but with a total size not
c             larger than limit+rscal)
              is = min (is + max (16, is), limit)
              rsiz2 = is + rscal
 
c             ----------------------------------------------------------
c             make sure enough room exists: garbage collection if needed
c             ----------------------------------------------------------
 
              if (rsiz2 .gt. itail-ihead) then
c                garbage collection:
                 info (14) = info (14) + 1
                 call ums2fg (xx, xsize, xhead, xtail, xuse,
     $                        ii, isize, ihead, itail, iuse,
     $                        cp, rp, dn, n, icntl, wir, wic, wr, wc,
     $                        0, 0, 0, 0, .true.,
     $                        pfree, xfree, mhead, mtail, slots)
c                at this point, iuse+rsiz2 <= iworst and xuse = xneed
                 pr = rp (row)
                 rsiz = ii (pr)
              endif
 
c             ----------------------------------------------------------
c             get space for the new copy
c             ----------------------------------------------------------
 
              pr2 = ihead
              ihead = ihead + rsiz2
              iuse = iuse + rsiz2
              info (18) = max (info (18), iuse)
              if (ihead .gt. itail) then
c                error return, if not enough integer memory:
                 go to 9000
              endif
 
c             ----------------------------------------------------------
c             make the copy, leaving hole in middle for element list
c             ----------------------------------------------------------
 
c             copy the rscal scalars, and the element list
cfpp$ nodepchk l
              do 1510 j = 0, rscal + 2*reln - 1
                 ii (pr2+j) = ii (pr+j)
1510          continue
 
c             copy the original entries
cfpp$ nodepchk l
              do 1520 j = 0, rlen - 1
                 ii (pr2+rsiz2-rlen+j) = ii (pr+rsiz-rlen+j)
1520          continue
 
              rp (row) = pr2
              ii (pr2) = rsiz2
 
c             ----------------------------------------------------------
c             deallocate the old copy of the row
c             ----------------------------------------------------------
 
              ii (pr+1) = -1
              pr = pr2
           endif
 
c          -------------------------------------------------------------
c          place the new (e,f) tuple in the element list of the row
c          -------------------------------------------------------------
 
           rep = (pr+2)
           ii (rep + 2*reln  ) = e
           ii (rep + 2*reln+1) = i - 1
           wr (row) = reln + 1
1530    continue
 
c=======================================================================
c  wrap-up of factorized frontal matrix is complete ]
c=======================================================================
 
1540    continue
c       exit label for loop 1540:
2000    continue
 
c=======================================================================
c=======================================================================
c  end of main factorization loop ]
c=======================================================================
c=======================================================================
 
c=======================================================================
c  wrap-up:  store lu factors in their final form [
c=======================================================================
 
c       ----------------------------------------------------------------
c       deallocate all remaining columns, rows, and frontal matrices
c       ----------------------------------------------------------------
 
        iuse = iuse - (ihead - 1)
        xuse = xuse - (xhead - 1)
        ineed = iuse
        xneed = xuse
        ihead = 1
        xhead = 1
 
        if (nlu .eq. 0) then
c          lu factors are completely empty (a = 0).
c          add one integer and one real, to simplify rest of code.
c          otherwise, some arrays in ums2rf or ums2so would have
c          zero size, which can cause an address fault.
           itail = isize
           xtail = xsize
           iuse = iuse + 1
           xuse = xuse + 1
           ineed = iuse
           xneed = xuse
           ip = itail
           xp = xtail
        endif
 
c       ----------------------------------------------------------------
c       compute permutation and inverse permutation vectors.
c       use wir/c for the row/col permutation, and wpr/c for the
c       inverse row/col permutation.
c       ----------------------------------------------------------------
 
        do 2010 k = 1, n
c          the kth pivot row and column:
           row = wpr (n-k+1)
           col = wpc (n-k+1)
           wir (k) = row
           wic (k) = col
2010    continue
c       replace wpr/c with the inversion permutations:
        do 2020 k = 1, n
           row = wir (k)
           col = wic (k)
           wpr (row) = k
           wpc (col) = k
2020    continue
 
        if (pgiven) then
c          the input matrix had been permuted from the original ordering
c          according to rperm and cperm.  combine the initial
c          permutations (now in rperm and cperm) and the pivoting
c          permutations, and place them back into rperm and cperm.
           do 2030 row = 1, n
              wm (wpr (row)) = rperm (row)
2030       continue
           do 2040 row = 1, n
              rperm (row) = wm (row)
2040       continue
           do 2050 col = 1, n
              wm (wpc (col)) = cperm (col)
2050       continue
           do 2060 col = 1, n
              cperm (col) = wm (col)
2060       continue
c       else
c          the input matrix was not permuted on input.  rperm and cperm
c          in ums2f1 have been passed to this routine as wir and wic,
c          which now contain the row and column permutations.  rperm and
c          cperm in this routine (ums2f2) are not defined.
        endif
 
c       ----------------------------------------------------------------
c       allocate nlu+3 integers for xtail, nlu, npiv and lup (1..nlu)
c       ----------------------------------------------------------------
 
        is = nlu + 5
        luip1 = itail
        itail = itail - is
        iuse = iuse + is
        ineed = iuse
        info (18) = max (info (18), iuse)
        info (19) = max (info (19), ineed)
        if (ihead .le. itail) then
 
c          -------------------------------------------------------------
c          sufficient memory exist to finish the factorization
c          -------------------------------------------------------------
 
           ii (itail+1) = nlu
           ii (itail+2) = npiv
           lupp = itail+5
           if (nlu .eq. 0) then
c             zero the dummy entries, if lu factors are empty
              ii (ip) = 0
              xx (xp) = zero
           endif
 
c          -------------------------------------------------------------
c          convert the lu factors into the new pivot order
c          -------------------------------------------------------------
 
           s = 0
           maxdr = 1
           maxdc = 1
           do 2100 k = 1, n
              e = wir (k)
              luip = rp (e)
              if (luip .gt. 0) then
c                this is an lu arrowhead - save a pointer in lup:
                 s = s + 1
c                update pointers to lu arrowhead relative to start of lu
                 ii (lupp+s-1) = luip - luip1 + 1
                 luxp = ii (luip)
                 ii (luip) = luxp - xtail + 1
c                convert the row and column indices to their final order
c                pattern of a column of l:
                 p = (luip + 7)
                 ludegc = ii (luip+3)
                 maxdc = max (maxdc, ludegc)
                 do 2070 j = 1, ludegc
                    ii (p) = wpr (abs (ii (p)))
                    p = p + 1
2070             continue
c                pattern of a row of u:
                 ludegr = ii (luip+2)
                 maxdr = max (maxdr, ludegr)
                 do 2080 j = 1, ludegr
                    ii (p) = wpc (abs (ii (p)))
                    p = p + 1
2080             continue
c                convert the lusons, usons, and lsons:
                 nsons = ii (luip+4)
                 do 2090 j = 1, nsons
                    eson = ii (p)
                    if (eson .le. n) then
c                      an luson
                       ii (p) = wm (eson)
                    else if (eson .le. 2*n) then
c                      a uson
                       ii (p) = wm (eson-n) + n
                    else
c                      an lson
                       ii (p) = wm (eson-2*n) + 2*n
                    endif
                    p = p + 1
2090             continue
c                renumber this lu arrowhead
                 wm (e) = s
              endif
2100       continue
 
           cmax = max (cmax, maxdc)
           rmax = max (rmax, maxdr)
           totnlu = totnlu + nlu
 
           ii (itail+3) = maxdc
           ii (itail+4) = maxdr
 
c          -------------------------------------------------------------
c          get memory usage for next call to ums2rf
c          -------------------------------------------------------------
 
           xruse = xruse - nz
           return
        endif
 
c=======================================================================
c  lu factors are now stored in their final form ]
c=======================================================================
 
c=======================================================================
c  error conditions
c=======================================================================
 
c       error return label:
9000    continue
        if (ihead .gt. itail .or. isize .lt. minmem) then
c          error return if out of integer memory
           call ums2er (1, icntl, info, -3, info (19))
        endif
        if (xhead .gt. xtail) then
c          error return if out of real memory
           call ums2er (1, icntl, info, -4, info (21))
        endif
        return
        end
