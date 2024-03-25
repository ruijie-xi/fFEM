 
        subroutine ums2rg (xx, xsize, xhead, xtail, xuse,
     $          lui, frdimc, frxp, frnext, frprev, nlu, lup,
     $          icntl, ffxp, ffsize, pfree, xfree)
c
cc UMS2RG is a utility which performs garbage collection for UMS2R2.
c
        integer lui (*), nlu, frdimc (nlu+2), frxp (nlu+2),
     $          frnext (nlu+2), frprev (nlu+2), lup (nlu),
     $          icntl (20), xsize, xuse, xhead, xtail, ffxp, ffsize,
     $          pfree, xfree
        real
     $          xx (xsize)
 
c=== ums2rg ============================================================
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
c  garbage collection for ums2r2.
 
c=======================================================================
c  input:
c=======================================================================
c
c       xx:             real workspace, containing matrix being
c                       factorized and partially-computed lu factors
c       xsize:          size of xx
c       xhead:          xx (1..xhead) is in use (matrix, frontal mtc's)
c       xtail:          xx (xtail..xsize) is in use (lu factors)
c       xuse:           memory usage in value
c       icntl:          integer control parameters, see ums2in
c       ffxp:           pointer to current contribution block
c       ffsize:         size of current contribution block
c       nlu:            number of lu arrowheads
c
c       frdimc (1..nlu+2)       leading dimension of frontal matrices
c       frxp (1..nlu+2)         pointer to frontal matrices in xx
c       frnext (1..nlu+2)       pointer to next block in xx
c       frprev (1..nlu+2)       pointer to previous block in xx
c       lup (1..nlu)            pointer to lu arrowhead patters in lui
c       lui (*)                 pattern of lu factors
 
c=======================================================================
c  output:
c=======================================================================
c
c       xx:             external fragmentation is removed at head
c       xhead:          xx (1..xhead) is in use, reduced in size
c       xuse:           memory usage in value, reduced
c       pfree:          pointer to free block in memory list, set to 0
c       xfree:          size of free block in xx, set to -1
c       frdimc          arrays for frontal matrices are compressed
c       frxp            frontal matrices have been shifted
c       ffxp            current working array has been shifted
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2r2
c       functions called:       abs
        intrinsic abs
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer xdp, i, e, fdimc, ludegr, ludegc, j, fluip, fxp,
     $          mhead, mtail
 
c  xdp:     real destination pointer, current block moved to xx (xdp...)
c  e:       an element
c  fdimc:   column dimension (number of rows) of a frontal matrix
c  ludegr:  row degree (number of columns) of a contribution block
c  ludegc:  column degree (number of rows) of a contribution block
c  fluip:   element is in lui (fluip...)
c  fxp:     element is in xx (fxp...) prior to compression
c  mhead:   nlu+1, head pointer for contribution block link list
c  mtail:   nlu+2, tail pointer for contribution block link list
c  i:       general loop index
c  j:       general loop index
 
c=======================================================================
c  executable statments:
c=======================================================================
 
c-----------------------------------------------------------------------
c  scan the link list and compress the reals
c-----------------------------------------------------------------------
 
        mhead = nlu+1
        mtail = nlu+2
        xdp = frxp (mhead)
        e = frnext (mhead)
 
c       while (e .ne. mtail) do
10      continue
        if (e .ne. mtail) then
 
           fdimc = frdimc (e)
 
c          -------------------------------------------------------------
           if (fdimc .eq. 0) then
c          -------------------------------------------------------------
 
c             this is a real hole - delete it from the link list
 
              frnext (frprev (e)) = frnext (e)
              frprev (frnext (e)) = frprev (e)
 
c          -------------------------------------------------------------
           else
c          -------------------------------------------------------------
 
c             this is an unassembled frontal matrix
              fxp = frxp (e)
              frxp (e) = xdp
              fluip = lup (e)
              ludegr = abs (lui (fluip+2))
              ludegc = abs (lui (fluip+3))
              if (fdimc .eq. ludegc) then
c                contribution block is already compressed
cfpp$ nodepchk l
                 do i = 0, (ludegr * ludegc) - 1
                    xx (xdp+i) = xx (fxp+i)
                 end do
              else
c                contribution block is not compressed
c                compress xx (fxp..) to xx (xdp..xdp+(ludegr*ludegc)-1)
                 do 40 j = 0, ludegr - 1
cfpp$ nodepchk l
                    do 30 i = 0, ludegc - 1
                       xx (xdp + j*ludegc + i) = xx (fxp + j*fdimc + i)
30                  continue
40               continue
                 frdimc (e) = ludegc
              endif
              xdp = xdp + ludegr*ludegc
 
           endif
 
c          -------------------------------------------------------------
c          get the next item in the link list
c          -------------------------------------------------------------
 
           e = frnext (e)
 
c       end while:
        goto 10
        endif
 
        frxp (mtail) = xdp
        pfree = 0
        xfree = -1
 
c       ----------------------------------------------------------------
c       shift the current working array (if it exists)
c       ----------------------------------------------------------------
 
        if (ffxp .ne. 0) then
cfpp$ nodepchk l
           do 50 i = 0, ffsize - 1
              xx (xdp+i) = xx (ffxp+i)
50         continue
           ffxp = xdp
           xdp = xdp + ffsize
        endif
 
c-----------------------------------------------------------------------
c  deallocate the unused space
c-----------------------------------------------------------------------
 
        xuse = xuse - (xhead - xdp)
        xhead = xdp
        return
        end
