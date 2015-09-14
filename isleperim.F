      subroutine isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
c=======================================================================
c         Island and Island Perimeter Mapping Routines
c
c
c          The main computational subroutine, expand, uses a "floodfill"
c          algorithm to expand one previously unmarked land
c          point to its entire connected land mass and its perimeter
c          ocean points.   Diagonally adjacent land points are
c          considered connected.  Perimeter "collisions" (i.e.,
c          ocean points that are adjacent to two unconnected
c          land masses) are detected and error messages generated.
c
c          Perimeter collisions must be removed because the stream 
c          function cannot be unambiguously defined on such cells.
c          They are fixed in one of three ways:
c          -Dfill_perimeter_violations fills the cell to land to form a
c            land bridge connection the two land masses.  Only one such
c            bridge is built, regardless of the number of perimeter
c            collisions between two land masses.
c          -Dwiden_perimeter_violations changes some land cells in one
c            of the islands into ocean cells to widen the channel between
c            the land masses to two "t" cells width.
c          -Dprompt_perimeter_violations enables interactive editing of
c            topography. 
c          The subroutine isleperim will not return until all perimeter
c          collisions are removed.
c
c          The subroutine expand uses a queue of size maxq of
c          coordinate pairs of candidate points.  Suggested
c          size for maxq is 4*(imt+jmt).  Queue overflow stops
c          execution with a message to increase the size of maxq.
c          Similarly a map with more that maxipp island perimeter
c          points or more than mnisle land masses stops execution
c          with an appropriate error message.
c
c          Computes map of land masses and island perimeters
c
c          Input:
c                  kmt = array of depths.  (0 for land) (>0 for ocean)
c          Outputs:
c                  map = map of land masses and their ocean perimeters
c                           mid-ocean cells are labeled 0
c                           land masses are labeled 1, 2, 3, ...,
c                           their perimeter ocean cells -1, -2, -3, ...,
c                  iperm = i coordinates of perimeter points
c                  jperm = j coordinates of perimeter points
c                  iofs = offset of each land mass in iperm, jperm
c                  nippts = number of island perimeter points by isle
c                  nisle = number of land masses
c          Array size inputs:
c                  imt = east/west array extents
c                  jmt = north/south array extents
c                  mnisle = maximum number of land masses
c                  maxipp = maximum number of island perimeter points
c
c          Arguments used in editing kmt field
c                  kmt_opt= selectable options for kmt changes
c                  kmt_changes = changes to kmt field
c                  nchanges = number of changes to kmt field
c                  i_del_kmt = i/o unit for changes in kmt
c                  xu = longitude (degrees) at "u" points
c                  yu = latitude (degrees) at "u" points
c                  zw = depth at bottom of cells
c
c         author: Charles Goldberg  ==> chg@gfdl.gov
c         December 1993/revised February 1995
c=======================================================================
c
      common /qsize/ maxqsize
c
      dimension kmt(imt,jmt)
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
#include "topog.h"
      dimension xu(imt), yu(jmt), zw(km)
c
      parameter (maxq=10000)
c
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
      integer ocean
c
      parameter (land=1, ocean=0)
      parameter (kmt_land=0, kmt_ocean=1)
c
      print '(/,a,/)','Finding perimeters of all land masses'
#if defined fill_perimeter_violations || \
   !defined widen_perimeter_violations
      print '(a,/)','Default action is "-Dfill_perimeter_violations"'
#else
      print '(a,/)','Default action is "-Dwiden_perimeter_violations"'
#endif
c
c     initialize number of changes to kmt
c
      nchanges = 0
c
    1 continue
c
c-----------------------------------------------------------------------
c     copy kmt to map changing notation
c     initially, 0 means ocean and 1 means unassigned land in map
c     as land masses are found, they are labeled 2, 3, 4, ...,
c     and their perimeter ocean cells -2, -3, -4, ...,
c     when no land points remain unassigned, land mass numbers are
c     reduced by 1 and their perimeter ocean points relabelled accordingly
c-----------------------------------------------------------------------
c
      do i=1,imt
        do j=1,jmt
          if (kmt(i,j) .gt. 0) then
            map(i,j) = ocean
          else
            map(i,j) = land
          end if
        end do
      end do
c
c-----------------------------------------------------------------------
c     find unassigned land points and expand them to continents
c-----------------------------------------------------------------------
c
      maxqsize = 0
      call qinit (iq, jq, qfront, qback)
      label = 2
      iofs(label) = 0
      nippts(label) = 0
      nerror = 0
#ifdef symmetry
      jnorth = jmt-1
#else
      jnorth = jmt
#endif
#ifdef cyclic
      iwest = 2
      ieast = imt-1
#else
      iwest = 1
      ieast = imt
#endif
      do j=jnorth,1,-1
        do i=iwest,ieast
          if (map(i,j) .eq. land) then
            call qpush (i, j, iq, jq, qfront, qback)
            call expand (map, label, iq, jq, qfront, qback, nerror
     &,                  iperm, jperm, iofs, nippts
     &,                  imt, jmt, km, mnisle, maxipp, kmt
     &,                  xu, yu, zw)
            print '(a,i2,a,i4)',
     &        'number of island perimeter points: nippts(',label-1,')=',
     &         nippts(label)
            label = label + 1
            if (label .gt. mnisle) then
              print '(a,i3,a)','ERROR==> mnisle=',mnisle,' is too small'
              stop '==> expand'
            end if
            iofs(label) = iofs(label-1) + nippts(label-1)
            nippts(label) = 0
          end if
        end do
      end do
      nisle = label - 1
c
c-----------------------------------------------------------------------
c     relabel land masses and their ocean perimeters
c------------------------------------------------------------------------
c
      do i=iwest,ieast
        do j=1,jnorth
          if (map(i,j) .ne. 0) then
            map(i,j) = map(i,j) - sign(1, map(i,j))
          end if
        end do
      end do
      do isle=2,nisle
        iofs(isle-1) = iofs(isle)
        nippts(isle-1) = nippts(isle)
      end do
      nisle = nisle - 1

#ifdef symmetry
      do i=iwest,ieast
        map(i,jmt) = map(i,jmt-1)
      end do
#endif
#ifdef cyclic
      do j=1,jmt
        map(1,j) = map(imt-1,j)
        map(imt,j) = map(2,j)
      end do
#endif
c
#ifdef debug_island_perimeters
      call showmap (map, imt, jmt, linewidth)
#endif
      print '(/a)',
     &       'The following changes fix "PERIMETER VIOLATIONS"'
c
      call enter_kmt_changes()
c
      if (nerror .gt. 0) then
        print *,' Island perimeter statistics:'
        print *,'maximum queue size was ',maxqsize
        print *,'number of land masses is ', nisle
        print *,'number of island perimeter points is ',
     &          nippts(nisle) + iofs(nisle)
        print *, ' '
        print *
     &, '==>Remapping land masses to see if PERIMETER VIOLATIONS remain'
        print *, ' '
        goto 1
      end if
c
      print *,' Island perimeter statistics:'
      print *,'maximum queue size was ',maxqsize
      print *,'number of land masses is ', nisle
      print *,'number of island perimeter points is ',
     &        nippts(nisle) + iofs(nisle)
c
#if !defined drive_topog && !defined permit_runtime_topog_changes
      if (auto_kmt_changes .and. n_del_kmt .ne. 0) then
        write (*,'(/a/a/a/a/a/)') 
     & ' =>Error: Automatic changes have been made to "kmt" to correct'
     &,'          PERIMETER VIOLATIONS. This may indicate a problem'
     &,'          with the "kmt" field. Using option'
     &,'         -Dpermit_runtime_topog_changes will accept the changes'
     &,'          and allow the model run to continue.'
        stop 'isleperim'
      endif
#endif

      return
      end



      subroutine showmap (map, imt, jmt, linewidth)
      dimension map(imt,jmt)
#if !defined narrow_map
      linewidth = 125
#else
      linewidth = 70
#endif
      print '(/,132a)',(' ',l=1,5+min(linewidth,imt)/2-13)
     &,     'Land Masses and Perimeters'
      istart = 0
      iremain = imt
      do isweep=1,imt/linewidth + 1
        iline = min(iremain, linewidth)
        iremain = iremain - iline
        if (iline .gt. 0) then
          print *, ' '
          print '(t6,32i5)', (istart+i+4,i=1,iline,5)
          do j=jmt,1,-1
            print '(i4,t6,160i1)', j,(mmm(map(istart+i,j)),i=1,iline)
          end do
          print '(t6,32i5)', (istart+i+4,i=1,iline,5)
          istart = istart + iline
        end if
      end do
      print *, ' '
      return
      end



      function mmm(m)
      if (m .eq. 0) then
        mmm = 0
      else if (m .gt. 0) then
        mmm = mod(m,10)
      else
        mmm = m
      end if
      return
      end

      subroutine expand (map, label, iq, jq, qfront, qback, nerror
     &,                  iperm, jperm, iofs, nippts
     &,                  imt, jmt, km, mnisle, maxipp, kmt
     &,                  xu, yu, zw)
c
c-----------------------------------------------------------------------
c          The subroutine expand uses a "flood fill" algorithm
c          to expand one previously unmarked land
c          point to its entire connected land mass and its perimeter
c          ocean points.   Diagonally adjacent land points are
c          considered connected.  Perimeter "collisions" (i.e.,
c          ocean points that are adjacent to two unconnected
c          land masses) are detected and error messages generated.
c
c          The subroutine expand uses a queue of size maxq of
c          coordinate pairs of candidate points.  Suggested
c          size for maxq is 4*(imt+jmt).  Queue overflow stops
c          execution with a message to increase the size of maxq.
c          Similarly a map with more that maxipp island perimeter
c          points or more than mnisle land masses stops execution
c          with an appropriate error message.
c-----------------------------------------------------------------------
c
      dimension map(imt,jmt), kmt(imt,jmt)

      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)

#include "topog.h"
      dimension xu(imt), yu(jmt), zw(km)
      character * 32 problem

      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
      logical qempty

      integer offmap, ocean
      parameter (offmap = -1)
      parameter (land = 1, ocean = 0)

      parameter (mnisle2=100)
      logical bridge_to(1:mnisle2)
c
      print '(a,i3)', 'Exploring land mass ',label-1
c
      if (mnisle2 .lt. mnisle) then 
        print '(a,i4,a)', 
     & 'ERROR:  change parameter (mnisle2=',mnisle,') in isleperim.F'
        stop '==>isleperim'
      end if
c
      do isle=1,mnisle
        bridge_to(isle) = .false.
      end do
c
c-----------------------------------------------------------------------
c     main loop:
c        Pop a candidate point off the queue and process it.
c-----------------------------------------------------------------------
c
 1000 continue

      if (qempty (iq, jq, qfront, qback)) then
        call qinit (iq, jq, qfront, qback)
        return
      else
        call qpop (i, j, iq, jq, qfront, qback)
c
c       case: (i,j) is off the map
        if (i .eq. offmap .or. j .eq. offmap) then
          goto 1000
c
c       case: map(i,j) is already labeled for this land mass
        else if (map(i,j) .eq. label) then
          goto 1000
c
c       case: map(i,j) is an ocean perimeter point of this land mass
        else if (map(i,j) .eq. -label) then
          goto 1000
c
c       case: map(i,j) is an unassigned land point
        else if (map(i,j) .eq. land) then
          map(i,j) = label
c         print *, 'labeling ',i,j,' as ',label
          call qpush (i,         jn(j,jmt), iq, jq, qfront, qback)
          call qpush (ie(i,imt), jn(j,jmt), iq, jq, qfront, qback)
          call qpush (ie(i,imt), j,         iq, jq, qfront, qback)
          call qpush (ie(i,imt), js(j,jmt), iq, jq, qfront, qback)
          call qpush (i,         js(j,jmt), iq, jq, qfront, qback)
          call qpush (iw(i,imt), js(j,jmt), iq, jq, qfront, qback)
          call qpush (iw(i,imt), j,         iq, jq, qfront, qback)
          call qpush (iw(i,imt), jn(j,jmt), iq, jq, qfront, qback)
          goto 1000
c
c       case: map(i,j) is an ocean point adjacent to this land mass
        else if (map(i,j) .eq. ocean .or. map(i,j) .lt. 0) then

c         subcase: map(i,j) is a perimeter ocean point of another mass
          if (map(i,j) .lt. 0) then
            nerror = nerror + 1
            print '(a,a,i3,a,i3,a,a,i3,a,i3)',
     &            'PERIMETER VIOLATION==> ',
     &            'map(',i,',',j,') is in the perimeter of both ',
     &            'land masses ', -map(i,j)-1, ' and ', label-1
c           if we just quit processing this point here, problem points
c           will be flagged several times.
c           if we relabel them, then they are only flagged once, but
c           appear in both island perimeters, which causes problems in
c           island integrals.  current choice is quit processing.
c
c           only fill first common perimeter point detected.
c           after the first land bridge is built, subsequent collisions
c           are not problems.
c
            if (.not. bridge_to(-map(i,j)-1)) then
              call clear_kmt_options ()
c
c             OK to set kmt to 0, but don`t change map because that would
c             require complete relabeling of land masses
c
c             option 1: fill common perimeter point to make land bridge
c
              call kmt_option (1, i, j, kmt(i,j), 0, kmt)
c
c             option 2: user selected interactive changes
c
c             option 3: change nearby land points on older land mass to ocean.
c                       we do not want to change points on the newer land
c                       mass because it may be incompletely explored at the
c                       current moment and some conflicting land points may
c                       be missed.
c
              do i1=-1,1
                do j1=-1,1
                  if (map(i+i1,j+j1) .eq. -map(i,j) .and.
     &                kmt(i+i1,j+j1) .eq. 0) then
                    call kmt_option (3, i+i1, j+j1
     &,                              kmt(i+i1,j+j1), kmt(i,j), kmt)
                  end if
                end do
              end do
c
              problem = 'perim'
              call select_option ('perim'
     &,                            i, j, kmt, xu, yu, zw)
c
c             see if option selected builds a land bridge
c
              do n=1,nchanges
                if (kmt_changes(n,1) .eq.i .and.
     &              kmt_changes(n,2) .eq.j .and.
     &              kmt_changes(n,4) .eq.0) then
                  bridge_to(-map(i,j)-1) = .true.
                end if
              end do
            end if
c
            goto 1000
          end if

c         case: map(i,j) is a ocean point--label it for current mass
          map(i,j) = -label
          nippts(label) = nippts(label) + 1
c         print *, 'iofs(label)=',iofs(label)
c         print *, 'nippts(label)=',nippts(label)
          if (iofs(label) + nippts(label) .gt. maxipp) then
            print *, 'ERROR==>  maxipp=',maxipp,' is not large enough'
            stop '==>expand'
          end if
          iperm(iofs(label) + nippts(label)) = i
          jperm(iofs(label) + nippts(label)) = j
          goto 1000
c
c       case: map(i,j) is probably labeled for another land mass
c       ************* this case should not happen **************
        else
          nerror = nerror + 1
          print '(a,a,i3,a,i3,a,a,i3,a,i3)',
     &          'ERROR ==>  ',
     &          'map(',i,',',j,') is labeled for both ',
     &          'land masses ', map(i,j)-1,' and ',label-1
        end if
        goto 1000
c
      end if
      return
      end


      subroutine qinit (iq, jq, qfront, qback)
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
c
      qfront = 1
      qback = 0
c
c     fake assignments to iq and jq to avoid "flint" warning
c
      iq(qfront) = 0
      jq(qfront) = 0
c
      return
      end


      subroutine qpush (i, j, iq, jq, qfront, qback)
c
      common /qsize/ maxqsize
c
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
c
      qback = qback + 1
      if (qback .gt. maxq) then
        if (qfront .ne. 1) then
c         shift queue left to make room
          ishift = qfront - 1
          do ip=qfront,qback-1
            iq(ip-ishift) = iq(ip)
            jq(ip-ishift) = jq(ip)
          end do
          qfront = 1
          qback = qback - ishift
        else
          print *, 'queue fault in qpush'
          stop '==>qpush'
        end if
      end if
      iq(qback) = i
      jq(qback) = j
c
      maxqsize = max(maxqsize, (qback-qfront))
c
      return
      end



      subroutine qpop (i, j, iq, jq, qfront, qback)
c
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
c
      i = iq(qfront)
      j = jq(qfront)
      qfront = qfront + 1
      return
      end


      function qempty (iq, jq, qfront, qback)
c
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
      logical qempty
c
      qempty = (qfront .gt. qback)
c      
      return
      end



      function jn(j,jmt)

c     j coordinate to the north of j

      integer offmap
      parameter (offmap = -1)

#ifdef symmetry
      if (j .lt. jmt-1) then
        jn = j + 1
      else if (j .eq. jmt-1) then
        jn = jmt-2
      else
        jn = offmap
      end if
#else
      if (j .lt. jmt) then
        jn = j + 1
      else
        jn = offmap
      end if
#endif
      return
      end


      function js(j,jmt)

c     j coordinate to the south of j

      integer offmap
      parameter (offmap = -1)

      if (j .gt. 1) then
        js = j - 1
      else
        js = offmap
      end if

      return
      end



      function ie(i,imt)

c     i coordinate to the east of i

      integer offmap
      parameter (offmap = -1)

#ifdef cyclic
      if (i .lt. imt-1) then
        ie = i + 1
      else
        ie = (i+1) - imt + 2
      end if
#else
      if (i .lt. imt) then
        ie = i + 1
      else
        ie = offmap
      end if
#endif
      return
      end



      function iw(i,imt)

c     i coordinate to the west of i

      integer offmap
      parameter (offmap = -1)

#ifdef cyclic
      if (i .gt. 2) then
        iw = i - 1
      else
        iw = (i-1) + imt - 2
      end if
#else
      if (i .gt. 1) then
        iw = i - 1
      else
        iw = offmap
      end if
#endif
      return
      end


