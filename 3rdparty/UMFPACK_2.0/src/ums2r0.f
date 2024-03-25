 
        subroutine ums2r0 (n, nz, cp, xx, xsize, ii, isize, xtail,
     $          itail, iuse, xuse, nzoff, nblks, icntl, cntl, info,
     $          rinfo, presrv, ap, ai, ax, an, anz, lui, luisiz,
     $          lublkp, blkp, offp, on, cperm, rperm, ne)
c
cc UMS2R0 is a utility which refactors a sparse matrix.
c
        integer n, nz, isize, ii (isize), icntl (20), info (40),
     $          cp (n+1), xsize, xtail, itail, iuse, xuse, an, anz,
     $          ap (an+1), ai (anz), luisiz, lui (luisiz), nblks,
     $          lublkp (nblks), blkp (nblks+1), on, offp (on+1),
     $          cperm (n), rperm (n), nzoff, ne
        logical presrv
        real
     $          xx (xsize), cntl (10), rinfo (20), ax (anz)
 
c=== ums2r0 ============================================================
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
c  refactorize an unsymmetric sparse matrix in column-form, optionally
c  permuting the matrix to upper block triangular form and factorizing
c  each diagonal block.
 
c=======================================================================
c  input:
c=======================================================================
c
c       n:              order of matrix
c       nz:             entries in matrix
c       cp (1..n+1):    column pointers of input matrix
c       presrv:         if true then preserve original matrix
c       xsize:          size of xx
c       isize:          size of ii
c       iuse:           memory usage in index on input
c       xuse:           memory usage in value on input
c       icntl:          integer control parameters, see ums2in
c       cntl:           real control parameters, see ums2in
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
c           an:                 0
c           anz:                1
c           ap:                 unused
c           ai:                 unused
c           ax:                 unused
c           ii (1..nz):         row indices of input matrix
c           xx (1..nz):         values of input matrix
c
c       information from prior lu factorization:
c
c       luisiz:                 size of lui
c       lui (1..luisiz):        patterns of lu factors, excluding
c                               prior preserved matrix (if it existed)
c                               and prior off-diagonal part (if it
c                               existed)
c       cperm (1..n):           column permutations
c       rperm (1..n):           row permutations
c       nblks:                  number of diagonal blocks for btf
c       if nblks > 1:
c           lublkp (1..nblks):  pointers to each diagonal lu factors
c           blkp (1..nblks+1):  index range of diagonal blocks
 
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
c                       the off-diagonal nonzeros, if nblks > 1
c
c       offp (1..n+1):  row pointers for off-diagonal part, if nblks > 1
c       info:           integer informational output, see ums2fa
c       rinfo:          real informational output, see ums2fa
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2rf
c       subroutines called:     ums2er, ums2r2, ums2p2, ums2ra, ums2of
c       functions called:       max
        intrinsic max
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer i, nzdia, p, ihead, nsgltn, nsym, wp, arip, arxp, npiv,
     $          wrksiz, nlu, prp, mc, mr, dummy1, dummy2, nz2, k, blk,
     $          k1, k2, kn, nzblk, col, row, prl, io, luip, mnz, arnz,
     $          xhead, offip, offxp, noutsd, nbelow, nzorig, xrmax
        real
     $          zero, one, a
        parameter (zero = 0.0, one = 1.0)
 
c  printing control:
c  -----------------
c  io:      i/o unit for diagnostic messages
c  prl:     printing level
c
c  allocated array pointers:
c  -------------------------
c  wp:      w (1..n+1), or w (1..kn+1), work array located in ii (wp...)
c  prp:     pr (1..n) work array located in ii (prp..prp+n-1)
c  arip:    ari (1..nz) array located in ii (arip..arip+nz-1)
c  arxp:    arx (1..nz) array located in xx (arxp..arxp+nz-1)
c  offip:   offi (1..nzoff) array located in ii (offip..offip+nzoff-1)
c  offxp:   offx (1..nzoff) array located in xx (offxp..offip+nzoff-1)
c
c  arrowhead-form matrix:
c  ----------------------
c  nz2:     number of entries in arrowhead matrix
c  nzblk:   number of entries in arrowhead matrix of a single block
c  arnz:    arrowhead form of blocks located in ii/xx (1..arnz)
c
c  btf information:
c  ----------------
c  k1:      starting index of diagonal block being factorized
c  k2:      ending index of diagonal block being factorized
c  kn:      the order of the diagonal block being factorized
c  blk:     block number of diagonal block being factorized
c  nsgltn:  number of 1-by-1 diagonal blocks ("singletons")
c  a:       numerical value of a singleton
c  mnz:     nzoff
c  noutsd:  entries in diagonal blocks, but not in lu (invalid)
c  nbelow:  entries below diagonal blocks (invalid)
c  nzoff:   entries above diagonal blocks (valid)
c  nzdia:   entries in diagonal blocks (valid)
c  nzorig:  total number of original entries
c
c  memory usage:
c  -------------
c  xhead:   xx (1..xhead-1) is in use, xx (xhead..xtail-1) is free
c  ihead:   ii (1..ihead-1) is in use, ii (ihead..itail-1) is free
c  wrksiz:  total size of work arrays need in ii for call to ums2r2
c  xrmax:   memory needed in value for next call to ums2rf
c
c  symbolic information and pattern of prior lu factors:
c  -----------------------------------------------------
c  nlu:     number of elements in a diagonal block
c  luip:    integer part of lu factors located in lui (luip...)
c  mr,mc:   largest frontal matrix for this diagonal block is mc-by-mr
c
c  other:
c  ------
c  k:       loop index (kth pivot)
c  i:       loop index
c  row:     row index
c  col:     column index
c  p:       pointer
c  nsym:    number of symmetric pivots chosen
c  dummy1:  argument returned by ums2ra, but not needed
c  dummy2:  argument returned by ums2ra, but not needed
 
c=======================================================================
c  executable statements:
c=======================================================================
 
        io = icntl (2)
        prl = icntl (3)
        nzorig = nz
 
        if (presrv) then
c          original matrix is not in cp/ii/xx, but in ap/ai/ax:
           ihead = 1
           xhead = 1
        else
           ihead = nz + 1
           xhead = nz + 1
        endif
 
        nzoff = 0
        nzdia = nz
        nsgltn = 0
        npiv = 0
        noutsd = 0
        nbelow = 0
        itail = isize + 1
        xtail = xsize + 1
 
c-----------------------------------------------------------------------
c current memory usage:
c-----------------------------------------------------------------------
 
c       if .not. presrv then
c               input matrix is now in ii (1..nz) and xx (1..nz)
c                       col pattern: ii (cp (col) ... cp (col+1))
c                       col values:  xx (cp (col) ... cp (col+1))
c               total: nz+n+1 integers, nz reals
c       else
c               input matrix is now in ai (1..nz) and ax (1..nz)
c                       col pattern: ai (ap (col) ... ap (col+1))
c                       col values:  ax (ap (col) ... ap (col+1))
c               cp is a size n+1 integer workspace
c               total: nz+2*(n+1) integers, nz reals
c
c       if (nblks > 1) then
c                       lublkp (1..nblks)
c                       blkp (1..nblks+1)
c                       offp (1..n+1)
c               total: (2*nblks+n+2) integers
c
c       cperm (1..n) and rperm (1..n)
c               total: 2*n integers
c
c   grand total current memory usage (including ii,xx,cp,ai,ap,ax
c       and lui):
c
c       presrv  nblks>1 integers, iuse =
c       f       f       luisiz + nz+  (n+1)+(2*n+7)
c       f       t       luisiz + nz+  (n+1)+(2*n+7)+(2*nblks+n+2)
c       t       f       luisiz + nz+2*(n+1)+(2*n+7)
c       t       t       luisiz + nz+2*(n+1)+(2*n+7)+(2*nblks+n+2)
c
c   real usage is xuse = nz
 
c-----------------------------------------------------------------------
c  get memory usage estimate for next call to ums2rf
c-----------------------------------------------------------------------
 
        xrmax = 2*ne
 
c-----------------------------------------------------------------------
c  convert matrix into arrowhead format (unless btf and preserved)
c-----------------------------------------------------------------------
 
        if (nblks .gt. 1 .and. presrv) then
 
c          -------------------------------------------------------------
c          btf is to be used, and original matrix is to be preserved.
c          it is converted and factorized on a block-by-block basis,
c          using the inverse row permutation (computed and stored in
c          offp (1..n)
c          -------------------------------------------------------------
 
           do k = 1, n
              offp (rperm (k)) = k
           end do
 
        else
 
c          -------------------------------------------------------------
c          convert the entire input matrix to arrowhead form
c          -------------------------------------------------------------
 
c          -------------------------------------------------------------
c          allocate workspace: w (n+1), pr (n), ari (nz), arx (nz)
c          -------------------------------------------------------------
 
           itail = itail - (2*n+1)
           iuse = iuse + 2*n+1
           prp = itail
           wp = prp + n
           iuse = iuse + nz
           xuse = xuse + nz
           arxp = xhead
           arip = ihead
           ihead = ihead + nz
           xhead = xhead + nz
           info (18) = max (info (18), iuse)
           info (20) = max (info (20), xuse)
           info (21) = max (info (21), xuse)
           if (ihead .gt. itail .or. xhead .gt. xtail) then
c             error return, if not enough integer and/or real memory:
              go to 9000
           endif
 
c          -------------------------------------------------------------
c          convert
c          -------------------------------------------------------------
 
           if (nblks .eq. 1) then
              if (presrv) then
                 call ums2ra (presrv, n, nz, cperm, rperm, ii (prp),
     $              ii (wp), nblks, xx (arxp), ii (arip), nzoff, nzdia,
     $              icntl, ap, blkp, ai, ax, info, offp, on, nz,
     $              0, n, nz2, i)
              else
                 call ums2ra (presrv, n, nz, cperm, rperm, ii (prp),
     $              ii (wp), nblks, xx (arxp), ii (arip), nzoff, nzdia,
     $              icntl, cp, blkp, ii, xx, info, offp, on, nz,
     $              0, n, nz2, i)
              endif
           else
c             note that presrv is false in this case
              call ums2ra (presrv, n, nz, cperm, rperm, ii (prp),
     $           ii (wp), nblks, xx (arxp), ii (arip), nzoff, nzdia,
     $           icntl, cp, blkp, ii, xx, info, offp, on, nz,
     $           0, n, nz2, nbelow)
           endif
 
c          -------------------------------------------------------------
c          copy the arrowhead pointers from w (1..n+1) to cp (1..n+1)
c          -------------------------------------------------------------
 
           do 20 i = 1, n+1
              cp (i) = ii (wp+i-1)
20         continue
 
c          -------------------------------------------------------------
c          deallocate w and pr.  if not presrv deallocate ari and arx
c          -------------------------------------------------------------
 
           iuse = iuse - (2*n+1)
           if (.not. presrv) then
c             ari and arx have been deallocated.
              xuse = xuse - nz
              iuse = iuse - nz
           endif
           itail = isize + 1
           xtail = xsize + 1
           nz = nz2
           ihead = nz + 1
           xhead = nz + 1
 
        endif
 
        info (5) = nz
        info (6) = nzdia
        info (7) = nzoff
        info (4) = nbelow
 
c-----------------------------------------------------------------------
c  refactorization
c-----------------------------------------------------------------------
 
c       ----------------------------------------------------------------
c       if nblks=1
c          arrowhead form is now stored in ii (1..nz) and xx (1..nz)
c          in reverse pivotal order (arrowhead n, n-1, ..., 2, 1).
c          the arrowhead form will be overwritten.
c       else if not presrv
c          off-diagonal part is in ii (1..nzoff), xx (1..nzoff),
c          (with row pointers offp (1..n+1)) followed by each diagonal
c          block (block 1, 2, ... nblks) in ii/xx (nzoff+1...nz).
c          each diagonal block is in arrowhead form, and in
c          reverse pivotal order (arrowhead k2, k2-1, ..., k1-1, k1).
c          the arrowhead form will be overwritten.
c       else (nblks > 1 and presrv)
c          ii and xx are still empty.  original matrix is in ap, ai,
c          and ax.  inverse row permutation (pr) is in offp (1..n).
c          the arrowhead form is not yet computed.
c       ----------------------------------------------------------------
 
        if (nblks .eq. 1) then
 
c          -------------------------------------------------------------
c          refactorize the matrix as a single block
c          -------------------------------------------------------------
 
           nlu = lui (2)
           mc = lui (4)
           mr = lui (5)
           wrksiz = 2*n + mr + 3*mc + 4*(nlu+2)
           itail = itail - wrksiz
           iuse = iuse + wrksiz
           p = itail
           info (18) = max (info (18), iuse)
           if (ihead .gt. itail) then
c             error return, if not enough integer memory:
              go to 9000
           endif
 
           call ums2r2 (cp, nz, n, xtail,
     $          xx, xsize, xuse, ii, cperm, rperm,
     $          icntl, cntl, info, rinfo, mc, mr,
     $          ii (p), ii (p+n), ii (p+2*n), ii (p+2*n+mr),
     $          ii (p+2*n+mr+mc), ii (p+2*n+mr+2*mc),
     $          ii (p+2*n+mr+3*mc), ii (p+2*n+mr+3*mc+(nlu+2)),
     $          ii (p+2*n+mr+3*mc+2*(nlu+2)),
     $          ii (p+2*n+mr+3*mc+3*(nlu+2)),
     $          nlu, lui (6), lui (nlu+6), noutsd,
     $          xrmax)
 
           if (info (1) .lt. 0) then
c             error return, if not enough real memory or bad pivot found
              go to 9010
           endif
 
c          -------------------------------------------------------------
c          deallocate workspace and original matrix (reals already done)
c          -------------------------------------------------------------
 
           iuse = iuse - wrksiz - nz
           itail = itail + wrksiz
           lui (1) = 1
           ihead = 1
           xhead = 1
 
        else
 
c          -------------------------------------------------------------
c          refactorize the block-upper-triangular form of the matrix
c          -------------------------------------------------------------
 
           if (presrv) then
c             count the entries in off-diagonal part
              nzoff = 0
           endif
 
           do 70 blk = nblks, 1, -1
 
c             ----------------------------------------------------------
c             factorize the kn-by-kn block, a (k1..k2, k1..k2)
c             ----------------------------------------------------------
 
c             get k1 and k2, the start and end of this block
              k1 = blkp (blk)
              k2 = blkp (blk+1) - 1
              kn = k2-k1+1
              a = zero
 
c             ----------------------------------------------------------
c             get pointers to, or place the block in, arrowhead form
c             ----------------------------------------------------------
 
              if (presrv) then
 
                 if (kn .gt. 1) then
 
c                   ----------------------------------------------------
c                   convert a single block to arrowhead format, using
c                   the inverse row permutation stored in offp
c                   ----------------------------------------------------
 
c                   ----------------------------------------------------
c                   compute nzblk, allocate ii/xx (1..nzblk), w(1..kn+1)
c                   ----------------------------------------------------
 
                    nzblk = 0
                    do 40 k = k1, k2
                       col = cperm (k)
cfpp$ nodepchk l
                       do 30 p = ap (col), ap (col+1) - 1
                          row = offp (ai (p))
                          if (row .lt. k1) then
c                            entry in off-diagonal part
                             nzoff = nzoff + 1
                          else if (row .le. k2) then
c                            entry in diagonal block
                             nzblk = nzblk + 1
                          endif
30                     continue
40                  continue
 
                    itail = itail - (kn+1)
                    wp = itail
                    ihead = nzblk + 1
                    xhead = nzblk + 1
                    iuse = iuse + nzblk + kn+1
                    xuse = xuse + nzblk
                    xrmax = max (xrmax, xuse)
                    info (18) = max (info (18), iuse)
                    info (20) = max (info (20), xuse)
                    info (21) = max (info (21), xuse)
                    if (ihead .gt. itail .or. xhead .gt. xtail) then
c                      error return, if not enough integer
c                      and/or real memory:
                       go to 9000
                    endif
 
c                   ----------------------------------------------------
c                   convert blk from column-form in ai/ax to arrowhead
c                   form in ii/xx (1..nzblk)
c                   ----------------------------------------------------
 
                    call ums2ra (presrv, n, nz, cperm, rperm, offp,
     $                 ii (wp), nblks, xx, ii, dummy1, dummy2,
     $                 icntl, ap, blkp, ai, ax, info, 0, 0, nzblk,
     $                 blk, kn, nz2, i)
 
c                   ----------------------------------------------------
c                   copy the arrowhead pointers from w (1..kn+1)
c                   to cp (k1 ... k2+1)
c                   ----------------------------------------------------
 
                    do 50 i = 0, kn
                       cp (k1+i) = ii (wp+i)
50                  continue
 
c                   cp (k1) is nzblk + 1 and cp (k2+1) is 1
 
c                   ----------------------------------------------------
c                   deallocate w (1..kn+1)
c                   ----------------------------------------------------
 
                    iuse = iuse - (kn+1)
                    itail = itail + (kn+1)
 
                 else
 
c                   ----------------------------------------------------
c                   get the value of singleton at a (k1,k1) if it exists
c                   ----------------------------------------------------
 
c                   find the diagonal entry in the unpermuted matrix,
c                   and count the entries in the diagonal and
c                   off-diagonal blocks.
                    col = cperm (k1)
                    do 60 p = ap (col), ap (col + 1) - 1
c                      inverse row permutation is stored in offp
                       row = offp (ai (p))
                       if (row .lt. k1) then
                          nzoff = nzoff + 1
                       else if (row .eq. k1) then
                          a = ax (p)
c                      else
c                         this is an invalid entry, below the diagonal
c                         block.  it will be detected (and optionally
c                         printed) in the call to ums2of below.
                       endif
60                  continue
 
                    ihead = 1
                    xhead = 1
                 endif
 
              else
 
c                -------------------------------------------------------
c                the block is located in ii/xx (cp (k2+1) ... cp (k1)-1)
c                and has already been converted to arrowhead form
c                -------------------------------------------------------
 
                 if (blk .eq. 1) then
c                   this is the last block to factorize
                    cp (k2+1) = nzoff + 1
                 else
                    cp (k2+1) = cp (blkp (blk-1))
                 endif
 
                 ihead = cp (k1)
                 xhead = ihead
 
                 if (kn .eq. 1) then
c                   singleton block in ii/xx (cp (k1+1) ... cp (k1)-1)
                    if (cp (k1) .gt. cp (k1+1)) then
                       a = xx (cp (k1) - 1)
                       ihead = ihead - 1
                       xhead = xhead - 1
                       iuse = iuse - 1
                       xuse = xuse - 1
                    endif
                 endif
 
              endif
 
c             ----------------------------------------------------------
c             factor the block
c             ----------------------------------------------------------
 
              if (kn .gt. 1) then
 
c                -------------------------------------------------------
c                the a (k1..k2, k1..k2) block is not a singleton.
c                block is now in ii/xx (cp (k2+1) ... cp (k1)-1), in
c                arrowhead form, and is to be overwritten with lu
c                -------------------------------------------------------
 
                 arnz = cp (k1) - 1
 
c                if (presrv) then
c                   ii/xx (1..arnz) holds just the current block, blk
c                else
c                   ii/xx (1..arnz) holds the off-diagonal part, and
c                   blocks 1..blk, in that order.
c                endif
 
                 luip = lublkp (blk)
c                luxp = lui (luip), not needed for refactorization
                 nlu = lui (luip+1)
c                npiv = lui (luip+2), not needed for refactorization
                 mc = lui (luip+3)
                 mr = lui (luip+4)
                 wrksiz = 2*kn + mr + 3*mc + 4*(nlu+2)
                 itail = itail - wrksiz
                 iuse = iuse + wrksiz
                 p = itail
                 info (18) = max (info (18), iuse)
                 if (ihead .gt. itail) then
c                   error return, if not enough integer memory:
                    go to 9000
                 endif
 
                 call ums2r2 (cp (k1), arnz, kn, xtail,
     $                xx, xtail-1, xuse, ii, cperm (k1), rperm (k1),
     $                icntl, cntl, info, rinfo, mc, mr,
     $                ii (p), ii (p+kn), ii (p+2*kn), ii (p+2*kn+mr),
     $                ii (p+2*kn+mr+mc), ii (p+2*kn+mr+2*mc),
     $                ii (p+2*kn+mr+3*mc), ii (p+2*kn+mr+3*mc+(nlu+2)),
     $                ii (p+2*kn+mr+3*mc+2*(nlu+2)),
     $                ii (p+2*kn+mr+3*mc+3*(nlu+2)),
     $                nlu, lui (luip+5), lui (luip+nlu+5), noutsd,
     $                xrmax)
 
                 if (info (1) .lt. 0) then
c                   error return, if not enough real memory or bad pivot
                    go to 9010
                 endif
 
c                -------------------------------------------------------
c                deallocate workspace and original matrix (reals
c                already deallocated in ums2r2)
c                -------------------------------------------------------
 
                 iuse = iuse - wrksiz
                 itail = itail + wrksiz
                 lui (luip) = xtail
                 iuse = iuse - (ihead - cp (k2+1))
                 ihead = cp (k2+1)
                 xhead = ihead
 
              else
 
c                -------------------------------------------------------
c                factor the singleton a (k1,k1) block, in a
c                -------------------------------------------------------
 
                 nsgltn = nsgltn + 1
                 if (a .eq. zero) then
c                   this is a singular matrix, replace with 1-by-1
c                   identity matrix.
                    a = one
                 else
c                   increment pivot count
                    npiv = npiv + 1
                 endif
                 xtail = xtail - 1
                 xuse = xuse + 1
                 xrmax = max (xrmax, xuse)
                 info (20) = max (info (20), xuse)
                 info (21) = max (info (21), xuse)
c                note: if the matrix is not preserved and nonsingular
c                then we will not run out of memory
                 if (xhead .gt. xtail) then
c                   error return, if not enough real memory:
                    go to 9000
                 endif
 
c                -------------------------------------------------------
c                store the 1-by-1 lu factors
c                -------------------------------------------------------
 
                 xx (xtail) = a
                 lublkp (blk) = -xtail
 
              endif
70         continue
 
c          -------------------------------------------------------------
c          make the index of each block relative to start of lu factors
c          -------------------------------------------------------------
cfpp$ nodepchk l
           do 80 blk = 1, nblks
              if (lublkp (blk) .gt. 0) then
                 lui (lublkp (blk)) = lui (lublkp (blk)) - xtail + 1
              else
c                this is a singleton
                 lublkp (blk) = (-lublkp (blk)) - xtail + 1
              endif
80         continue
 
c          -------------------------------------------------------------
c          store the off-diagonal blocks
c          -------------------------------------------------------------
 
           if (presrv) then
 
c             ----------------------------------------------------------
c             allocate temporary workspace for pr (1..n) at head of ii
c             ----------------------------------------------------------
 
              prp = ihead
              ihead = ihead + n
              iuse = iuse + n
 
c             ----------------------------------------------------------
c             allocate permanent copy of off-diagonal blocks
c             ----------------------------------------------------------
 
              itail = itail - nzoff
              offip = itail
              xtail = xtail - nzoff
              offxp = xtail
              iuse = iuse + nzoff
              xuse = xuse + nzoff
              xrmax = max (xrmax, xuse)
              info (18) = max (info (18), iuse)
              info (20) = max (info (20), xuse)
              info (21) = max (info (21), xuse)
              if (ihead .gt. itail .or. xhead .gt. xtail) then
c                error return, if not enough integer and/or real memory:
                 go to 9000
              endif
 
c             ----------------------------------------------------------
c             re-order the off-diagonal blocks according to pivot perm
c             ----------------------------------------------------------
 
c             use cp as temporary work array:
              mnz = nzoff
              if (nzoff .eq. 0) then
c                offi and offx are not accessed in ums2of.  set offip
c                and offxp to 1 (since offip = itail = isize+1, which
c                can generate an address fault, otherwise).
                 offip = 1
                 offxp = 1
              endif
              call ums2of (cp, n, rperm, cperm, nzoff,
     $             offp, ii (offip), xx (offxp), ii (prp),
     $             icntl, ap, ai, ax, an, anz, presrv, nblks, blkp,
     $             mnz, 2, info, nbelow)
 
c             ----------------------------------------------------------
c             deallocate pr (1..n)
c             ----------------------------------------------------------
 
              ihead = 1
              xhead = 1
              iuse = iuse - n
 
           else
 
c             off-diagonal entries are in ii/xx (1..nzoff); shift down
c             to ii/xx ( ... itail/xtail).  no extra memory needed.
              do 90 i = nzoff, 1, -1
                 ii (itail+i-nzoff-1) = ii (i)
                 xx (xtail+i-nzoff-1) = xx (i)
90            continue
              ihead = 1
              xhead = 1
              itail = itail - nzoff
              xtail = xtail - nzoff
           endif
 
        endif
 
c       ----------------------------------------------------------------
c       clear the flags (negated row/col indices, and negated ludegr/c)
c       ----------------------------------------------------------------
 
        do 100 i = 1, luisiz
           lui (i) = abs (lui (i))
100     continue
 
c-----------------------------------------------------------------------
c  normal and error return
c-----------------------------------------------------------------------
 
c       error return label:
9000    continue
        if (ihead .gt. itail) then
c          set error flag if not enough integer memory
           call ums2er (2, icntl, info, -3, info (18))
        endif
        if (xhead .gt. xtail) then
c          set error flag if not enough real memory
           call ums2er (2, icntl, info, -4, info (21))
        endif
 
c       error return label, for error return from ums2r2:
9010    continue
 
c       cp can now be deallocated in ums2rf:
        iuse = iuse - (n+1)
 
        info (4) = noutsd + nbelow
        nzdia = nzorig - nzoff - noutsd - nbelow
        info (5) = nzoff + nzdia
        info (6) = nzdia
        info (7) = nzoff
        info (8) = nsgltn
        info (9) = nblks
        info (12) = info (10) + info (11) + n + info (7)
 
c       count the number of symmetric pivots chosen.  note that some
c       of these may have been numerically unacceptable.
        nsym = 0
        do 110 k = 1, n
           if (cperm (k) .eq. rperm (k)) then
c             this kth pivot came from the diagonal of a
              nsym = nsym + 1
           endif
110     continue
        info (16) = nsym
 
        info (17) = info (17) + npiv
        rinfo (1) = rinfo (4) + rinfo (5) + rinfo (6)
 
c       set warning flag if entries outside prior pattern are present
        if (info (4) .gt. 0) then
           call ums2er (2, icntl, info, 1, -info (4))
        endif
 
c       set warning flag if matrix is singular
        if (info (1) .ge. 0 .and. info (17) .lt. n) then
           call ums2er (2, icntl, info, 4, info (17))
        endif
 
c       ----------------------------------------------------------------
c       return memory usage estimate for next call to ums2rf
c       ----------------------------------------------------------------
 
        info (23) = xrmax
 
        return
        end
