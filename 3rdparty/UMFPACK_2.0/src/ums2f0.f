 
        subroutine ums2f0 (n, nz, cp, xx, xsize, ii, isize, xtail,
     $          itail, iuse, xuse, nzoff, nblks, icntl, cntl, info,
     $          rinfo, presrv, ap, ai, ax, an, anz, keep, ne)
c
cc UMS2F0 factors an unsymmetric sparse matrix in column form.
c
        integer n, nz, isize, ii (isize), icntl (20), info (40),
     $          cp (n+1), xsize, xtail, itail, iuse, xuse, an, anz,
     $          ap (an+1), ai (anz), keep (20), nzoff, nblks, ne
        logical presrv
        real
     $          xx (xsize), cntl (10), rinfo (20), ax (anz)
 
c=== ums2f0 ============================================================
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
c  factorize an unsymmetric sparse matrix in column-form, optionally
c  permuting the matrix to upper block triangular form and factorizing
c  each diagonal block.
 
c=======================================================================
c  input:
c=======================================================================
c
c       n:              order of matrix
c       nz:             entries in matrix, after removing duplicates
c                       and invalid entries.
c       ne:             number of triplets, unchanged from ums2fa
c       cp (1..n+1):    column pointers of input matrix
c       presrv:         if true then preserve original matrix
c       xsize:          size of xx
c       isize:          size of ii
c       iuse:           memory usage in index on input
c       xuse:           memory usage in value on input
c       icntl:          integer control parameters, see ums2in
c       cntl:           real control parameters, see ums2in
c       keep (6..8):    integer control parameters, see ums2in
c
c       if presrv is true:
c           an:                 = n, order of preserved matrix
c           anz:                = anz, order of preserved matrix
c           ap (1..an+1):       column pointers of preserved matrix
c           ai (1..nz):         row indices of preserved matrix
c           ax (1..nz):         values of preserved matrix
c           ii:                 unused on input
c           xx:                 unused on input
c       else
c           an:                 1
c           anz:                1
c           ap:                 unused
c           ai:                 unused
c           ax:                 unused
c           ii (1..nz):         row indices of input matrix
c           xx (1..nz):         values of input matrix
 
c=======================================================================
c  output:
c=======================================================================
c
c       xx (xtail ... xsize), xtail:
c
c                       lu factors are located in xx (xtail ... xsize),
c                       including values in off-diagonal part if matrix
c                       was permuted to block triangular form.
c
c       ii (itail ... isize), itail:
c
c                       lu factors are located in ii (itail ... isize),
c                       including pattern, row and column permutations,
c                       block triangular information, etc.  see umf2fa
c                       for more information.
c
c       info:           integer informational output, see ums2fa
c       rinfo:          real informational output, see ums2fa
c
c       iuse:           memory usage in index on output
c       xuse:           memory usage in value on output
c
c       nzoff:          entries in off-diagonal part (0 if btf not used)
c       nblks:          number of diagonal blocks (1 if btf not used)
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2fa
c       subroutines called:     ums2er, ums2fb, ums2f1, ums2of
c       functions called:       max
        intrinsic max
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer kn, nzdia, blkpp, lublpp, p, offip, xhead, row,
     $          offxp, offpp, ihead, k1, k2, blk, prp, p2, cpermp,
     $          rpermp, nsgltn, npiv, mnz, nsym, k, col, rmax, cmax,
     $          totnlu, xrmax, xruse
        logical trybtf, iout, xout
        real
     $          zero, one, a
        parameter (zero = 0.0, one = 1.0)
 
c  allocated array pointers:
c  -------------------------
c  blkpp:   blkp (1..nblks+1) array located in ii (blkpp..blkp+nblks)
c  lublpp:  lublkp (1..nblks) array loc. in ii (lublpp..lublpp+nblks-1)
c  offip:   offi (1..nzoff) array located in ii (offip..offip+nzoff-1)
c  offxp:   offx (1..nzoff) array located in xx (offxp..offxp+nzoff-1)
c  offpp:   offp (1..n+1) array located in ii (offpp..offpp+n)
c  cpermp:  cperm (1..n) array located in ii (cpermp..cpermp+n-1)
c  rpermp:  rperm (1..n) array located in ii (rpermp..rpermp+n-1)
c  prp:     pr (1..n) work array located in ii (prp..prp+n-1)
c
c  btf information:
c  ----------------
c  k1:      starting index of diagonal block being factorized
c  k2:      ending index of diagonal block being factorized
c  kn:      the order of the diagonal block being factorized
c  blk:     block number of diagonal block being factorized
c  trybtf:  true if btf is to be attempted (= icntl (4) .eq. 1)
c  nzdia:   number of entries in diagonal blocks (= nz if btf not used)
c  nsgltn:  number of 1-by-1 diagonal blocks ("singletons")
c  npiv:    number of numerically valid singletons
c  a:       numerical value of a singleton
c  mnz:     nzoff
c
c  memory usage:
c  -------------
c  xhead:   xx (1..xhead-1) is in use, xx (xhead..xtail-1) is free
c  ihead:   ii (1..ihead-1) is in use, ii (ihead..itail-1) is free
c  iout:    true if ums2f1 ran out of integer memory, but did not
c           set error flag
c  xout:    true if ums2f2 ran out of integer memory, but did not
c           set error flag
c
c  estimated memory for ums2rf:
c  ----------------------------
c  rmax:    largest contribution block is cmax-by-rmax
c  cmax:       "         "         "    "   "   "  "
c  totnlu:  total number of lu arrowheads in all diagonal blocks
c  xrmax:   estimated maximum real memory usage for ums2rf
c  xruse:   estimated current real memory usage for ums2rf
c
c  other:
c  ------
c  k:       loop index (kth pivot)
c  row:     row index
c  col:     column index
c  p:       pointer
c  p2:      pointer
c  nsym:    number of symmetric pivots chosen
 
c=======================================================================
c  executable statements:
c=======================================================================
 
c-----------------------------------------------------------------------
c  get input parameters and initialize
c-----------------------------------------------------------------------
 
        nblks = 1
        nzoff = 0
        nzdia = nz
        nsgltn = 0
        npiv = 0
        rmax = 1
        cmax = 1
        totnlu = 0
        if (presrv) then
c          original matrix is not in cp/ii/xx, but in ap/ai/ax:
           ihead = 1
           xhead = 1
        else
           ihead = nz + 1
           xhead = nz + 1
           endif
        itail = isize + 1
        xtail = xsize + 1
 
c-----------------------------------------------------------------------
c  allocate permutation arrays: cperm (1..n) and rperm (1..n), and
c  seven scalars:  transa, nzoff, nblks, presrv, nz, n, ne
c  (in that order) at tail of ii (in lu factors)
c-----------------------------------------------------------------------
 
        itail = itail - (2*n+7)
        iuse = iuse + (2*n+7)
        info (18) = max (info (18), iuse)
        info (19) = info (18)
        cpermp = itail
        rpermp = cpermp + n
        if (ihead .gt. itail) then
c          error return, if not enough integer memory:
           go to 9000
           endif
 
c-----------------------------------------------------------------------
c  find permutations to block upper triangular form, if requested.
c-----------------------------------------------------------------------
 
        trybtf = icntl (4) .eq. 1
        if (trybtf) then
 
c          -------------------------------------------------------------
c          get workspace at tail of ii of size 6n+2
c          -------------------------------------------------------------
 
           itail = itail - (n+1)
           offpp = itail
           itail = itail - (5*n+1)
           p = itail
           iuse = iuse + (6*n+2)
           info (18) = max (info (18), iuse)
           info (19) = info (18)
 
c          -------------------------------------------------------------
           if (presrv) then
c          find permutation, but do not convert to btf form
c          -------------------------------------------------------------
 
              if (ihead .gt. itail) then
c                error return, if not enough integer memory:
                 go to 9000
                 endif
              call ums2fb (ax, anz, ai, anz, n, nz, nzdia, nzoff,
     $           nblks, cp, ii (cpermp), ii (rpermp), ii(p), ii(p+n),
     $           ii (p+2*n), ii (p+3*n), ii (p+4*n), ii (offpp),
     $           presrv, icntl)
 
c          -------------------------------------------------------------
           else
c          find permutation, convert to btf form, and discard original
c          -------------------------------------------------------------
 
c             use additional size nz temporary workspace in ii and xx
              ihead = ihead + nz
              xhead = xhead + nz
              iuse = iuse + nz
              xuse = xuse + nz
              info (18) = max (info (18), iuse)
              info (20) = max (info (20), xuse)
              info (19) = info (18)
              info (21) = info (20)
              if (ihead .gt. itail .or. xhead .gt. xtail) then
c                error return, if not enough integer and/or real memory:
                 go to 9000
                 endif
              call ums2fb (xx, 2*nz, ii, 2*nz, n, nz, nzdia, nzoff,
     $              nblks, cp, ii (cpermp), ii (rpermp), ii(p), ii(p+n),
     $              ii (p+2*n), ii (p+3*n), ii (p+4*n), ii (offpp),
     $              presrv, icntl)
c             deallocate extra workspace in ii and xx
              ihead = ihead - nz
              xhead = xhead - nz
              iuse = iuse - nz
              xuse = xuse - nz
              endif
 
c          -------------------------------------------------------------
c          deallocate workspace, and allocate btf arrays if required
c          -------------------------------------------------------------
 
           if (nblks .gt. 1) then
c             replace (6*n+2) workspace at tail of ii with
c             blkp (1..nblks+1) and lublkp (1..nblks), offp (1..n+1)
              blkpp = offpp - (nblks+1)
              lublpp = blkpp - (nblks)
              itail = lublpp
              iuse = iuse - (6*n+2) + (2*nblks+n+2)
           else
c             the matrix is irreducible.  there is only one block.
c             remove everything at tail of ii, except
c             for the 2*n permutation vectors and the 7 scalars.
c             (transa, nzoff, nblks, presrv, nz, n, ne).
              itail = (isize + 1) - (2*n+7)
              iuse = iuse - (6*n+2)
              endif
 
           endif
 
c-----------------------------------------------------------------------
c current memory usage:
c-----------------------------------------------------------------------
 
c       if .not. presrv then
c               input matrix is now in ii (1..nz) and xx (1..nz)
c               off-diagonal part: in ii/xx (1..nzoff)
c                       col pattern: ii (offp (col) ... offp (col+1))
c                       col values:  xx (offp (col) ... offp (col+1))
c               diagonal blocks: in ii/xx (nzoff+1..nz)
c                       col pattern: ii (cp (col) ... cp (col+1))
c                       col values:  xx (cp (col) ... cp (col+1))
c               total: nz+n+1 integers, nz reals
c       else
c               input matrix is now in ai (1..nz) and ax (1..nz),
c               in original (non-btf) order:
c                       col pattern: ai (ap (col) ... ap (col+1))
c                       col values:  ax (ap (col) ... ap (col+1))
c               cp is a size n+1 integer workspace
c               total: nz+2*(n+1) integers, nz reals
c
c       if (nblks > 1) then
c               at tail of ii (in order): 2*nblks+n+2
c                       lublkp (1..nblks)
c                       blkp (1..nblks+1)
c                       offp (1..n+1)
c               total: (2*nblks+n+2) integers
c
c       remainder at tail of ii:
c               cperm (1..n)
c               rperm (1..n)
c               seven scalars: transa, nzoff, nblks, presrv, nz, n, ne
c
c   grand total current memory usage (including ii,xx,cp,ai,ap,ax):
c
c       presrv  nblks>1         integers, iuse =
c       f       f               nz+  (n+1)+(2*n+7)
c       f       t               nz+  (n+1)+(2*n+7)+(2*nblks+n+2)
c       t       f               nz+2*(n+1)+(2*n+7)
c       t       t               nz+2*(n+1)+(2*n+7)+(2*nblks+n+2)
c
c   real usage is xuse = nz
 
c       ----------------------------------------------------------------
c       get memory usage for next call to ums2rf
c       ----------------------------------------------------------------
 
        xrmax = 2*ne
        xruse = nz
 
c-----------------------------------------------------------------------
c factorization
c-----------------------------------------------------------------------
 
        if (nblks .eq. 1) then
 
c          -------------------------------------------------------------
c          factorize the matrix as a single block
c          -------------------------------------------------------------
 
           call ums2f1 (cp, n, ii (cpermp), ii (rpermp), nzoff,
     $          itail, xtail, xx, xsize, xuse, ii, itail-1, iuse,
     $          icntl, cntl, info, rinfo, nblks,
     $          ap, ai, ax, presrv, 1, an, anz, ii, keep,
     $          rmax, cmax, totnlu, xrmax, xruse, iout, xout)
           if (iout .or. xout) then
c             error return, if not enough integer and/or real memory:
              go to 9000
              endif
           if (info (1) .lt. 0) then
c             error return, if error in ums2f2:
              go to 9010
              endif
c          original matrix has been deallocated
           ihead = 1
           xhead = 1
 
c          -------------------------------------------------------------
c          make the index of the block relative to start of lu factors
c          -------------------------------------------------------------
 
           ii (itail) = 1
 
        else
 
c          -------------------------------------------------------------
c          factorize the block-upper-triangular form of the matrix
c          -------------------------------------------------------------
 
           prp = offpp
           if (presrv) then
c             count the off-diagonal entries during factorization
              nzoff = 0
c             compute temp inverse permutation in ii (prp..prp+n-1)
cfpp$ nodepchk l
              do 10 k = 1, n
                 ii (prp + ii (rpermp+k-1) - 1) = k
10            continue
           endif
 
           do 30 blk = nblks, 1, -1
 
c             ----------------------------------------------------------
c             factorize the kn-by-kn block, a (k1..k2, k1..k2)
c             ----------------------------------------------------------
 
c             get k1 and k2, the start and end of this block
              k1 = ii (blkpp+blk-1)
              k2 = ii (blkpp+blk) - 1
              kn = k2-k1+1
              if (.not. presrv) then
                 p = cp (k1)
                 cp (k2+1) = ihead
              endif
 
              if (kn .gt. 1) then
 
c                -------------------------------------------------------
c                factor the block (the block is not a singleton)
c                -------------------------------------------------------
 
                 call ums2f1 (cp (k1), kn,
     $              ii (cpermp+k1-1), ii (rpermp+k1-1), nzoff,
     $              itail, xtail, xx, xtail-1, xuse, ii, itail-1,
     $              iuse, icntl, cntl, info, rinfo, nblks,
     $              ap, ai, ax, presrv, k1, an, anz, ii (prp), keep,
     $              rmax, cmax, totnlu, xrmax, xruse, iout, xout)
                 if (iout .or. xout) then
c                   error return, if not enough int. and/or real memory:
                    go to 9000
                 endif
                 if (info (1) .lt. 0) then
c                   error return, if error in ums2f2:
                    go to 9010
                 endif
                 if (presrv) then
                    ihead = 1
                    xhead = 1
                 else
                    ihead = p
                    xhead = p
                 endif
 
c                -------------------------------------------------------
c                save the location of the lu factors in lubkp (blk)
c                -------------------------------------------------------
 
                 ii (lublpp+blk-1) = itail
 
              else
 
c                -------------------------------------------------------
c                get the value of singleton at a (k1,k1), if it exists
c                -------------------------------------------------------
 
                 a = zero
                 if (presrv) then
c                   find the diagonal entry in the unpermuted matrix
                    col = ii (cpermp + k1 - 1)
                    do 20 p2 = ap (col), ap (col + 1) - 1
                       row = ii (prp + ai (p2) - 1)
                       if (row .lt. k1) then
c                         entry in off-diagonal blocks
                          nzoff = nzoff + 1
                       else
                          a = ax (p2)
                       endif
20                  continue
                    ihead = 1
                    xhead = 1
                 else if (p .ne. ihead) then
                    a = xx (p)
                    ihead = p
                    xhead = p
                    iuse = iuse - 1
                    xuse = xuse - 1
                    xruse = xruse - 1
                 endif
 
c                -------------------------------------------------------
c                store the 1-by-1 lu factors of a singleton
c                -------------------------------------------------------
 
                 nsgltn = nsgltn + 1
                 if (a .eq. zero) then
c                   the diagonal entry is either not present, or present
c                   but numerically zero.  this is a singular matrix,
c                   replace with 1-by-1 identity matrix.
                    a = one
                 else
c                   increment pivot count
                    npiv = npiv + 1
                 endif
                 xtail = xtail - 1
c                note: if the matrix is not preserved and nonsingular
c                then we will not run out of memory at this point.
                 xuse = xuse + 1
                 xruse = xruse + 1
                 xrmax = max (xrmax, xruse)
                 info (20) = max (info (20), xuse)
                 info (21) = max (info (21), xuse)
c                error return, if not enough real memory:
                 if (xhead .gt. xtail) then
                    go to 9000
                 endif
                 ii (lublpp+blk-1) = -xtail
                 xx (xtail) = a
 
              endif
 
30         continue
 
c          -------------------------------------------------------------
c          make the index of each block relative to start of lu factors
c          -------------------------------------------------------------
 
cfpp$ nodepchk l
           do 40 p = lublpp, lublpp + nblks - 1
              if (ii (p) .gt. 0) then
                 ii (ii (p)) = ii (ii (p)) - xtail + 1
                 ii (p) = ii (p) - itail + 1
              else
c                this is a singleton
                 ii (p) = (-ii (p)) - xtail + 1
              endif
40         continue
 
c          -------------------------------------------------------------
c          allocate temporary workspace for pr (1..n) at head of ii
c          -------------------------------------------------------------
 
           prp = ihead
           ihead = ihead + n
           iuse = iuse + n
 
c          -------------------------------------------------------------
c          allocate a single entry in case the lu factors are empty
c          -------------------------------------------------------------
 
           if (nblks .eq. n) then
c             otherwise, arrays in ums2rf and ums2so would have
c             zero size, which can cause an address fault later on
              itail = itail - 1
              iuse = iuse + 1
              p2 = itail
           endif
 
c          -------------------------------------------------------------
c          allocate permanent copy of off-diagonal blocks
c          -------------------------------------------------------------
 
           itail = itail - nzoff
           offip = itail
           xtail = xtail - nzoff
           offxp = xtail
           iuse = iuse + nzoff
           xuse = xuse + nzoff
           xruse = xruse + nzoff
           xrmax = max (xrmax, xruse)
           info (18) = max (info (18), iuse)
           info (19) = max (info (19), iuse)
           info (20) = max (info (20), xuse)
           info (21) = max (info (21), xuse)
           if (ihead .gt. itail .or. xhead .gt. xtail) then
c             error return, if not enough integer and/or real memory:
              go to 9000
           endif
 
c          -------------------------------------------------------------
c          re-order the off-diagonal blocks according to pivot perm
c          -------------------------------------------------------------
 
c          use cp as temporary work array:
           mnz = nzoff
           if (presrv) then
              call ums2of (cp, n, ii (rpermp), ii (cpermp), nzoff,
     $          ii (offpp), ii (offip), xx (offxp), ii (prp),
     $          icntl, ap, ai, ax, an, anz, presrv, nblks, ii (blkpp),
     $          mnz, 1, info, p)
           else
              call ums2of (cp, n, ii (rpermp), ii (cpermp), nzoff,
     $          ii (offpp), ii (offip), xx (offxp), ii (prp),
     $          icntl, 0, ii, xx, 0, mnz, presrv, 0, 0,
     $          mnz, 1, info, p)
           endif
           if (nblks .eq. n) then
c             zero the only entry in the integer part of the lu factors
              ii (p2) = 0
           endif
 
c          -------------------------------------------------------------
c          deallocate pr (1..n), and ii/xx (1..nzoff) if present
c          -------------------------------------------------------------
 
           ihead = 1
           xhead = 1
           iuse = iuse - n
           if (.not. presrv) then
              iuse = iuse - nzoff
              xuse = xuse - nzoff
           endif
 
        endif
 
c-----------------------------------------------------------------------
c  normal and error return
c-----------------------------------------------------------------------
 
c       error return label:
9000    continue
        if (iout .or. ihead .gt. itail) then
c          set error flag if not enough integer memory
           call ums2er (1, icntl, info, -3, info (19))
        endif
        if (xout .or. xhead .gt. xtail) then
c          set error flag if not enough real memory
           call ums2er (1, icntl, info, -4, info (21))
        endif
 
c       error return label, for error from ums2f2:
9010    continue
 
        info (4) = 0
        nzdia = nz - nzoff
        info (5) = nz
        info (6) = nzdia
        info (7) = nzoff
        info (8) = nsgltn
        info (9) = nblks
        info (12) = info (10) + info (11) + n + info (7)
 
c       count the number of symmetric pivots chosen.  note that some of
c       these may have been numerically unacceptable.
        nsym = 0
        if (info (1) .ge. 0) then
           do k = 1, n
              if (ii (cpermp+k-1) .eq. ii (rpermp+k-1)) then
c                this kth pivot came from the diagonal of a
                 nsym = nsym + 1
              endif
           end do
        endif
        info (16) = nsym
 
        info (17) = info (17) + npiv
        rinfo (1) = rinfo (4) + rinfo (5) + rinfo (6)
 
        if (info (1) .ge. 0 .and. info (17) .lt. n) then
c          set warning flag if matrix is singular
           call ums2er (1, icntl, info, 4, info (17))
        endif
 
c       ----------------------------------------------------------------
c       determine an upper bound on the amount of integer memory needed
c       (lindex) for a subsequent call to ums2rf.  if block-upper-
c       triangular-form is not in use (info (9) is 1), then
c       this bound is exact.  if ne is higher in the call to ums2rf
c       than in the call to ums2fa, then add 3 integers for each
c       additional entry (including the 2 integers required for the
c       row and column indices of the additional triplet itself).
c       this estimate assumes that job and transa are the same in
c       ums2fa and ums2rf.
c       ----------------------------------------------------------------
 
c       (keep (5) - keep (4) + 1), is added to info (22)
c       in ums2fa, to complete the computation of the estimate.
 
        if (presrv) then
           info (22) = max (3*ne+2*n+1, ne+3*n+2,
     $                           2*nz+4*n+10+rmax+3*cmax+4*totnlu)
        else
           info (22) = max (3*ne+2*n+1, ne+3*n+2, 2*nz+3*n+2,
     $                             nz+3*n+ 9+rmax+3*cmax+4*totnlu)
        endif
 
c       ----------------------------------------------------------------
c       approximate the amount of real memory needed (lvalue) for a
c       subsequent call to ums2rf.  the approximation is an upper bound
c       on the bare minimum amount needed.  some garbage collection may
c       occur, but ums2rf is guaranteed to finish if given an lvalue of
c       size info (23) and if the pattern is the same.  if ne is
c       higher in the call to ums2rf than in the call to ums2fa, then
c       add 2 reals for each additional entry (including the 1 real
c       required for the value of the additional triplet itself).
c       this estimate assumes that job and transa are the same in
c       ums2fa and ums2rf.
c       ----------------------------------------------------------------
 
        info (23) = xrmax
        return
        end
