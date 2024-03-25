 
        subroutine ums2fg (xx, xsize, xhead, xtail, xuse,
     $          ii, isize, ihead, itail, iuse,
     $          cp, rp, dn, n, icntl, wir, wic, wr, wc,
     $          ffxp, ffsize, wxp, ffdimc, doslot,
     $          pfree, xfree, mhead, mtail, slots)
c
cc UMS2FG is a utility routine which performs garbage collection.
c
        integer n, dn, isize, ii (isize), ihead, itail, rp (n+dn),
     $          cp (n+1), icntl (20), wir (n), wic (n), xsize, xuse,
     $          iuse, xhead, xtail, ffxp, ffsize, wxp,
     $          ffdimc, wr (n+dn), wc (n+dn), pfree, xfree, mhead,
     $          mtail, slots
        logical doslot
        real
     $          xx (xsize)
 
c=== ums2fg ============================================================
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
c  garbage collection for ums2f2.
 
c=======================================================================
c  input:
c=======================================================================
c
c       ii/xx:          integer/real workspace, containing matrix being
c                       factorized and partially-computed lu factors
c       isize:          size of ii
c       xsize:          size of xx
c       xhead:          xx (1..xhead) is in use (matrix, frontal mtc's)
c       xtail:          xx (xtail..xsize) is in use (lu factors)
c       xuse:           memory usage in value
c       ihead:          ii (1..ihead) is in use (matrix, frontal mtc's)
c       itail:          ii (itail..isize) is in use (lu factors)
c       iuse:           memory usage in index
c       cp (1..n+1):    pointers to columns
c       rp (1..n+dn):   pointers to rows, frontal matrices, and lu
c                       arrowheads
c       dn:             number of dense columns
c       n:              order of matrix
c       icntl:          integer control parameters, see ums2in
c       wr (1..n):      see ums2f2
c       wc (1..n):      see ums2f2
c       ffxp:           pointer to current contribution block
c       ffsize:         size of current contribution block
c       mhead:          pointer to first block in memory list
c       mtail:          pointer to last block in memory list
c       doslot:         true if adding slots
c       if doslot:
c           wir (1..n)  if wir (row) >= 0 then add (or keep) an extra
c                       slot in the row's element list
c           wic (1..n)  if wir (col) >= 0 then add (or keep) an extra
c                       slot in the col's element list
 
c=======================================================================
c  output:
c=======================================================================
c
c       ii/xx:          external fragmentation is removed at head
c       xhead:          xx (1..xhead) is in use, reduced in size
c       xuse:           memory usage in value, reduced
c       ihead:          ii (1..ihead) is in use, reduced in size
c       iuse:           memory usage in index, reduced
c       pfree:          pointer to free block in memory list, set to 0
c       xfree:          size of free block in xx, set to -1
c       mhead:          pointer to first block in memory list
c       mtail:          pointer to last block in memory list
c       ffxp            current working array has been shifted
c       wxp             current work vector has been shifted
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2f2
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer what, fsiz, row, col, p, idp, xdp, i, e, ep, fdimc,
     $          ludegr, ludegc, j, pc, celn, clen, reln, rlen,
     $          csiz1, csiz2, rsiz1, rsiz2, fluip, cxp, fxp, rdeg,
     $          cdeg, cscal, rscal, fscal
        parameter (cscal = 9, rscal = 2, fscal = 7)
        logical slot
 
c  compression:
c  ------------
c  what:    what this block of memory is (a row, column, etc.)
c  idp:     int. destination pointer, current block moved to ii (idp...)
c  xdp:     real destination pointer, current block moved to xx (xdp...)
c  slot:    true if adding, or keeping, a size-2 hole in an element list
c
c  columns:
c  --------
c  cscal:   = 9, the number of scalars in column data structure
c  celn:    number of (e,f) tuples in element list of a column
c  clen:    number of unassembled original entries in a column
c  cdeg:    degree of a column (number of entries, including fill-in)
c  cxp:     a column is in xx (cxp...) prior to compression
c  pc:      column is in ii (pc ...) prior to compression
c  csiz1:   size of a column in ii, prior to compression
c  csiz2:   size of a column in ii, after compression
c  col:     column index
c
c  rows:
c  -----
c  rscal:   = 2, the number of scalars in row data structure
c  reln:    number of (e,f) tuples in element list of a row
c  rlen:    number of unassembled original entries in a row
c  rsiz1:   size of a row in ii, prior to compression
c  rsiz2:   size of a row in ii, after compression
c  rdeg:    degree of a row (number of entries, including fill-in)
c  row:     row index
c
c  frontal matrices:
c  -----------------
c  fscal:   = 7, the number of scalars in element data structure
c  fluip:   element is in ii (fluip...) prior to compression
c  fxp:     a frontal matrix is in xx (fxp...) prior to compression
c  e:       an element
c  fdimc:   column dimension (number of rows) of a frontal matrix
c  ludegr:  row degree (number of columns) of a contribution block
c  ludegc:  column degree (number of rows) of a contribution block
c  fsiz:    size of an artificial frontal matrix
c  ep:      an artificial frontal matrix is in ii (ep ...) prior to comp
c
c  other:
c  ------
c  p:       pointer
c  i:       general loop index
c  j:       general loop index
 
c=======================================================================
c  executable statments:
c=======================================================================
 
        slots = 0
 
c-----------------------------------------------------------------------
c   prepare the non-pivotal rows/cols and unassembled elements
c-----------------------------------------------------------------------
 
c       place the size of each block of memory at the beginning,
c       and mark the 2nd entry in each block with what it is
 
c       ----------------------------------------------------------------
c       mark the columns
c       ----------------------------------------------------------------
 
cfpp$ nodepchk l
        do col = 1, n
           pc = cp (col)
           if (pc .ne. 0) then
c             this is a non-pivotal, non-null column
              cdeg = ii (pc+1)
              cp (col) = cdeg
              ii (pc+1) = col+n
           endif
        end do
 
c       ----------------------------------------------------------------
c       mark the rows and frontal matrices
c       ----------------------------------------------------------------
 
cfpp$ nodepchk l
        do 20 row = 1, n
           p = rp (row)
           rlen = wc (row)
           if (p .eq. 0) then
c             a pivotal row
              continue
           else if (rlen .ge. 0 .and. rlen .le. n) then
c             this is a non-pivotal, non-null row
              rdeg = ii (p+1)
              rp (row) = rdeg
              ii (p+1) = row+2*n
           else if (wr (row) .eq. -(n+dn+2)) then
c             a pivotal row, and an assembled element
              continue
           else
c             this is an unassembled frontal matrix
c             the size is implicitly fscal
              fdimc = ii (p+1)
              rp (row) = fdimc
              ii (p+1) = row
           endif
20      continue
 
c       ----------------------------------------------------------------
c       mark the artificial frontal matrices
c       ----------------------------------------------------------------
 
cfpp$ nodepchk l
        do 30 e = n+1, n+dn
           ep = rp (e)
           if (ep .ne. 0) then
c             this is an unassembled artificial frontal matrix
c             the size is ii (ep+1) + cscal
              fdimc = ii (ep+1)
              rp (e) = fdimc
              ii (ep+1) = e+2*n
           endif
30      continue
 
c-----------------------------------------------------------------------
c  scan the link list and compress the reals
c-----------------------------------------------------------------------
 
        xdp = 1
        p = mhead
c       while (p .ne. 0) do
40      continue
        if (p .ne. 0) then
 
           what = ii (p+1)
 
c          -------------------------------------------------------------
           if (what .gt. 3*n) then
c          -------------------------------------------------------------
 
c             this is an unassembled artificial frontal matrix
              e = what - 2*n
              fxp = ii (p+2)
              ii (p+2) = xdp
cfpp$ nodepchk l
              do 50 j = 0, rp (e) - 1
                 xx (xdp+j) = xx (fxp+j)
50            continue
              xdp = xdp + rp (e)
 
c          -------------------------------------------------------------
           else if (what .eq. -1 .or. ii (p+6) .eq. 0) then
c          -------------------------------------------------------------
 
c             this is a real hole - delete it from the link list
              if (ii (p+4) .ne. 0) then
                 ii (ii (p+4)+3) = ii (p+3)
              else
                 mhead = ii (p+3)
              endif
              if (ii (p+3) .ne. 0) then
                 ii (ii (p+3)+4) = ii (p+4)
              else
                 mtail = ii (p+4)
              endif
 
c          -------------------------------------------------------------
           else if (what .le. n) then
c          -------------------------------------------------------------
 
c             this is an unassembled frontal matrix
              e = what
              fxp = ii (p+2)
              ii (p+2) = xdp
              fluip = ii (p)
              ludegr = ii (fluip+2)
              ludegc = ii (fluip+3)
              fdimc = rp (e)
              if (fdimc .eq. ludegc) then
c                contribution block is already compressed
cfpp$ nodepchk l
                 do 60 i = 0, (ludegr * ludegc) - 1
                    xx (xdp+i) = xx (fxp+i)
60               continue
              else
c                contribution block is not compressed
c                compress xx (fxp..) to xx (xdp..xdp+(ludegr*ludegc)-1)
                 do 80 j = 0, ludegr - 1
cfpp$ nodepchk l
                    do 70 i = 0, ludegc - 1
                       xx (xdp + j*ludegc + i) = xx (fxp + j*fdimc + i)
70                  continue
80               continue
                 rp (e) = ludegc
              endif
              xdp = xdp + ludegr*ludegc
 
c          -------------------------------------------------------------
           else if (what .le. 2*n) then
c          -------------------------------------------------------------
 
c             this is a column
              cxp = ii (p+2)
              ii (p+2) = xdp
              clen = ii (p+6)
cfpp$ nodepchk l
              do 90 j = 0, clen - 1
                 xx (xdp+j) = xx (cxp+j)
90            continue
              xdp = xdp + clen
 
c          -------------------------------------------------------------
           endif
c          -------------------------------------------------------------
 
c          -------------------------------------------------------------
c          get the next item in the link list
c          -------------------------------------------------------------
 
           p = ii (p+3)
 
c       end while:
        goto 40
        endif
 
        pfree = 0
        xfree = -1
 
c       ----------------------------------------------------------------
c       shift the current working array (if it exists)
c       ----------------------------------------------------------------
 
        if (ffxp .ne. 0) then
cfpp$ nodepchk l
           do 100 i = 0, ffsize - 1
              xx (xdp+i) = xx (ffxp+i)
100        continue
           ffxp = xdp
           xdp = xdp + ffsize
        endif
 
c       ----------------------------------------------------------------
c       shift the current work vector (if it exists)
c       ----------------------------------------------------------------
 
        if (wxp .ne. 0) then
           wxp = xdp
           xdp = xdp + ffdimc
        endif
 
c-----------------------------------------------------------------------
c  scan from the top of integer memory (1) to bottom (ihead) and
c  compress the integers
c-----------------------------------------------------------------------
 
        p = 1
        idp = p
c       while (p .lt. ihead) do:
110     continue
        if (p .lt. ihead) then
 
           what = ii (p+1)
 
c          -------------------------------------------------------------
           if (what .gt. 3*n) then
c          -------------------------------------------------------------
 
c             this is an unassembled artificial frontal matrix
              e = what - 2*n
              fsiz = rp (e) + cscal
              ii (p+1) = rp (e)
              rp (e) = idp
cfpp$ nodepchk l
              do 120 i = 0, fsiz - 1
                 ii (idp+i) = ii (p+i)
120           continue
c             shift pointers in the link list
              if (ii (idp+4) .ne. 0) then
                 ii (ii (idp+4)+3) = idp
              else
                 mhead = idp
              endif
              if (ii (idp+3) .ne. 0) then
                 ii (ii (idp+3)+4) = idp
              else
                 mtail = idp
              endif
              p = p + fsiz
              idp = idp + fsiz
 
c          -------------------------------------------------------------
           else if (what .eq. -1) then
c          -------------------------------------------------------------
 
c             this is a integer hole
              p = p + ii (p)
 
c          -------------------------------------------------------------
           else if (what .ge. 1 .and. what .le. n) then
c          -------------------------------------------------------------
 
c             this is an unassembled frontal matrix (fscal integers)
              e = what
              fdimc = rp (e)
              ii (p+1) = fdimc
              rp (e) = idp
cfpp$ nodepchk l
              do 130 i = 0, fscal - 1
                 ii (idp+i) = ii (p+i)
130           continue
c             shift pointers in the link list
              if (ii (idp+4) .ne. 0) then
                 ii (ii (idp+4)+3) = idp
              else
                 mhead = idp
              endif
              if (ii (idp+3) .ne. 0) then
                 ii (ii (idp+3)+4) = idp
              else
                 mtail = idp
              endif
              p = p + fscal
              idp = idp + fscal
 
c          -------------------------------------------------------------
           else if (what .le. 2*n) then
c          -------------------------------------------------------------
 
c             this is a non-pivotal column
              csiz1 = ii (p)
              col = what - n
              celn = ii (p+5)
              clen = ii (p+6)
              csiz2 = 2*celn + clen + cscal
              slot = doslot .and. wic (col) .ge. 0 .and. p .ge. idp+2
              if (slot) then
c                keep (or make) one extra slot for element list growth
                 csiz2 = csiz2 + 2
                 slots = slots + 2
              endif
              cdeg = cp (col)
              ii (p+1) = cdeg
              cp (col) = idp
              ii (p) = csiz2
c             copy the cscal scalars and the celn (e,f) tuples
cfpp$ nodepchk l
              do 140 i = 0, cscal + 2*celn - 1
                 ii (idp+i) = ii (p+i)
140           continue
              if (clen .gt. 0) then
c                shift pointers in the link list
                 if (ii (idp+4) .ne. 0) then
                    ii (ii (idp+4)+3) = idp
                 else
                    mhead = idp
                 endif
                 if (ii (idp+3) .ne. 0) then
                    ii (ii (idp+3)+4) = idp
                 else
                    mtail = idp
                 endif
              endif
              p = p + csiz1 - clen
              idp = idp + cscal + 2*celn
              if (slot) then
c                skip past the slot
                 idp = idp + 2
              endif
c             copy the clen original row indices
cfpp$ nodepchk l
              do 150 i = 0, clen - 1
                 ii (idp+i) = ii (p+i)
150           continue
              p = p + clen
              idp = idp + clen
 
c          -------------------------------------------------------------
           else
c          -------------------------------------------------------------
 
c             this is a non-pivotal row
              rsiz1 = ii (p)
              row = what - 2*n
              reln = wr (row)
              rlen = wc (row)
              rsiz2 = 2*reln + rlen + rscal
              slot = doslot .and. wir (row) .ge. 0 .and. p .ge. idp+2
              if (slot) then
c                keep (or make) one extra slot for element list growth
                 rsiz2 = rsiz2 + 2
                 slots = slots + 2
              endif
              rdeg = rp (row)
              ii (p+1) = rdeg
              rp (row) = idp
              ii (p) = rsiz2
c             copy the rscal scalars, and the reln (e,f) tuples
cfpp$ nodepchk l
              do 160 i = 0, rscal + 2*reln - 1
                 ii (idp+i) = ii (p+i)
160           continue
              p = p + rsiz1 - rlen
              idp = idp + rscal + 2*reln
              if (slot) then
c                skip past the slot
                 idp = idp + 2
              endif
c             copy the rlen original column indices
cfpp$ nodepchk l
              do 170 i = 0, rlen - 1
                 ii (idp+i) = ii (p+i)
170           continue
              p = p + rlen
              idp = idp + rlen
 
c          -------------------------------------------------------------
           endif
c          -------------------------------------------------------------
 
c          -------------------------------------------------------------
c          move to the next block
c          -------------------------------------------------------------
 
c       end while:
        goto 110
        endif
 
c-----------------------------------------------------------------------
c  deallocate the unused space
c-----------------------------------------------------------------------
 
        iuse = iuse - (ihead - idp)
        ihead = idp
        xuse = xuse - (xhead - xdp)
        xhead = xdp
        return
        end
