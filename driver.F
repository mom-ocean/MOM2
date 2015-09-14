      program driver
c
c=======================================================================
c
c     This is the main driver for the Modular Ocean Model (MOM). The  
c     structure allows for surface boundary conditions (SBC) to
c     come from a dataset or an atmospheric model. Integration time
c     is divided into a number of equal time segments and SBC are held 
c     fixed for each time segment. When coupling to an atmospheric
c     model, SBC are supplied each time segment (the coupling period)
c     and held fixed for that period. When supplying SBC to MOM from a
c     dataset, the length of each time segment is one ocean time step.
c     
c     author:        r. c. pacanowski      e-mail=> rcp@gfdl.gov
c                    a. rosati             e-mail=>  ar@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "csbc.h"
#include "iounit.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
c
      write (stdout,'(/,48x,a,//)') '====  Driving MOM 2  ===='
c
c-----------------------------------------------------------------------
c     initialize i/o units for MOM
c-----------------------------------------------------------------------
c
      call ioinit
c
c-----------------------------------------------------------------------
c     initialize timing analysis for various parts of MOM
c-----------------------------------------------------------------------
c
      call inittimers
c
c-----------------------------------------------------------------------
c     Initialize S.B.C. names and related items.
c     
c     Dimensional units for the S.B.C. are for the model where they
c     will be used. eg: Heat flux will be used as a S.B.C. for the
c     ocean in cal/cm**2/sec.
#ifdef coupled
c     If the atmosphere calculates heat flux in watts/m**2 then the
c     conversion factor is 0.2389e-4 to go from watts/m**2 to
c     cal/cm**2/sec.
c     Note: the conversion factors are multiplicative except for
c     SST where the conversion factor is added.
c
c     default convergence criteria are for extrapolating S.B.C. into
c     land on the model grid where they were constructed. units are for
c     the model grid where the S.B.C. was constructed.
c     as a guide, the
c     criteria should only be chosen large enough to give reasonable 
c     values near the coastlines. whether the field is extrapolated
c     well into the interior of continents is presumed unimportant.
#endif
c-----------------------------------------------------------------------
c
      do n=1,maxsbc
        mapsbc(n)  = 0
	coabc(n)   = c0
	crits(n)   = c0
	sbcname(n) = 'unknown'
	dunits(n)  = 'unknown'
      enddo
c
      write (stdout,*) ' '
      write (stdout,*)
     & ' Surface Boundary Conditions are defined as follows:'
      write (stdout,*) ' '
      write (stdout,*) ' for the ocean (applied to MOM):'
c
      do n=1,numsbc
c
c       order S.B.C.'s such that S.B.C. "n" has index "m"=mapsbc(n)
c       within arrays "sbcocn" and "sbcatm"
c
        mapsbc(n) = n
        m = mapsbc(n)
	if (m .eq. 1) then
	  sbcname(m) = ' taux     '
	  dunits(m)  = ' dynes/cm**2   '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-3
	else if (m .eq. 2) then
	  sbcname(m) = ' tauy     '
	  dunits(m)  = ' dynes/cm**2   '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-3
	else if (m .eq. 3) then
	  sbcname(m) = ' heat flux'
	  dunits(m)  = ' cal/cm**2/sec '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-5
	else if (m .eq. 4) then
	  sbcname(m) = ' salt flux'
	  dunits(m)  = ' gm/cm**2/sec  '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-6
	else if (m .eq. 5) then
	  sbcname(m) = 'short wave'
	  dunits(m)  = ' cal/cm**2/sec '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-5
	else if (m .eq. 6) then
	  sbcname(m) = ' SST      '
	  dunits(m)  = ' deg C         '
	  coabc(m)   = 0.0
	  crits(m)   = 0.05
	else if (m .eq. 7) then
	  sbcname(m) = ' SSS      '
	  dunits(m)  = ' (S-35.0)/1000 '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-5
	else if (m .eq. 8) then
	  sbcname(m) = ' surface u'
	  dunits(m)  = ' cm/sec        '
	  coabc(m)   = 1.0
	  crits(m)   = 0.1
	else if (m .eq. 9) then
	  sbcname(m) = ' surface v'
	  dunits(m)  = '  cm/sec       '
	  coabc(m)   = 1.0
	  crits(m)   = 0.1
	endif
        if (n .eq. numosbc+1) then
          write (stdout,*) ' '
          write (stdout,*) ' for the atmosphere (from MOM):'
	endif
	if (n .le. numosbc) then
	  write (stdout,*) ' S.B.C. # ',m,' is ',sbcname(m)
     &,', atmos native units are ', dunits(m)
#ifdef coupled
     &,', conversion factor to MOM  is ',coabc(m)
     &,', extrap crit =',crits(m)
#endif
        else
	  write (stdout,*) ' S.B.C. # ',m,' is ',sbcname(m)
     &,', MOM native units are ', dunits(m)
#ifdef coupled
     &,', conversion factor to atmos is ',coabc(m)
     &,', extrap crit =',crits(m)
#endif
	endif
      enddo
c
c-----------------------------------------------------------------------
c     do the introductory ocean setup once per run
c-----------------------------------------------------------------------
c
      call setocn
c
c-----------------------------------------------------------------------
c     do the introductory atmosphere setup once per run
c-----------------------------------------------------------------------
c
      write (stdout,'(/a36/)') ' ==> Note: the atmos setup follows:'
c
#ifdef timing
      call tic ('driver', 'setatm')
#endif
#ifdef simple_sbc
      write (stdout,'(/,a,/)')
     &'==> Note: Simple  S.B.C. are implemented in setvbc.F'
      dtatm = 0.0
#else
c
c     "setatm" must do the following:
c       1) set up the atmospheric S.B.C. grid definition
c       2) define the atmosphere land/sea mask
c       3) set the atmosphere time step "dtatm" {seconds}
c
      call setatm (dtatm)
#endif
#ifdef timing
      call toc ('driver', 'setatm')
#endif
c
#if !defined coupled
c
c     when the MOM S.B.C. come from a dataset, force the segment time
c     and atmospheric time step to one MOM time step. This will force
c     the number of segments to one and the number of time steps per
c     segment to represent the length of the run in days. 
c
      dtatm = dtts
      segtim = dtts*secday
#endif
c
c-----------------------------------------------------------------------
c     compute the number of ocean time steps "numots" for this run and 
c     the number of ocean time steps per ocean segment "ntspos".
c     compute the number of atmos time steps "numats" for this run and
c     the number of atmos time steps per atmos segment "ntspas".
c     divide the integration time "days" into "numseg" segments.
c     each will be length "segtim" days. Surface boundary conditions
c     are supplied every "segtim" days.
c-----------------------------------------------------------------------
c
      numots = nint(rundays/(dtts*secday))
      ntspos = nint(segtim/(dtts*secday))
      numats = nint(rundays/(dtatm*secday))
      ntspas = nint(segtim/(dtatm*secday))
      numseg = numots/ntspos
c
#ifdef coupled
      write (stdout,8800) rundays, numseg, segtim, ntspos, ntspas, dtts
     &,                   dtatm
#else
      write (stdout,8800) rundays, numseg, segtim, ntspos, dtts
#endif
c
c-----------------------------------------------------------------------
c     prepare initial condition S.B.C. for the models
c-----------------------------------------------------------------------
c
#if !defined simple_sbc
      write (stdout,'(/,1x,a)')
     &    '==> Initializing all S.B.C. for MOM now.'
c
c     initialize S.B.C. for the atmosphere (on the ocean S.B.C. grid)
c
      do n=1,numsbc
        m = mapsbc(n)
        do jrow=1,jmt
	  do i=1,imt
	    sbcocn(i,jrow,m) = c0
          enddo
	enddo
      enddo
c
c     load the tracers (SST & SSS) for each row "j". (zero on land)
c     load from the MW if fully opened otherwise load from disk
c
      isst = mapsbc(6)
      isss = mapsbc(7)
      if (wide_open_mw) then
        do jrow=1,jmt
	  do i=1,imt
	    if (isst .ne. 0) sbcocn(i,jrow,isst) = t(i,1,jrow,1,taup1)        
	    if (isss .ne. 0) sbcocn(i,jrow,isss) = t(i,1,jrow,2,taup1)
	  enddo
	enddo
      else
        do jrow=1,jmt
	  if (isst .ne. 0) call getst (jrow, sbcocn(1,1,isst), 1)         
	  if (isss .ne. 0) call getst (jrow, sbcocn(1,1,isss), 2)
        enddo
      endif
#endif
c
#ifdef coupled
      if (init .and. numasbc .gt. 0) then
c
c       initialize S.B.C. for the ocean (on the atmos S.B.C. grid)
c
        do n=1,numsbc
	  m = mapsbc(n)
          do jrow=1,jma
            do i=1,imap2
              sbcatm(i,jrow,m) = c0
            enddo
	  enddo
	enddo
c
      else
c
c       when restarting ... read in the necessary S.B.C. data
c       to preserve values in land for extrapolations.
c
        call getunit (iosbc, 'restrt_sbc.dta', 'u s r')
        read (iosbc) stamp, iotext
        read (iosbc) iimt, ijmt, iima, ijma, inbc
        read (iosbc) stamp, iotext
        read (iosbc) sbcatm
        read (iosbc) stamp, iotext
        read (iosbc) sbcocn
        write (stdout,9100) stamp
	call relunit (iosbc)
      endif
#endif
c
c-----------------------------------------------------------------------
c     check for consistancy in the S.B.C. setup
c-----------------------------------------------------------------------
c
      call chkcpl (dtatm, dtts)
c
c
c
c
c-----------------------------------------------------------------------
c     S T A R T    S E G M E N T    L O O P
c-----------------------------------------------------------------------
c
      do n=1,numseg
#ifdef coupled
        write (stdout,9000) 'A T M O S', n, stamp
c
c-----------------------------------------------------------------------
c       get the atmospheric S.B.C. from MOM
c-----------------------------------------------------------------------
c
        call gasbc (n)
#endif
c
#if !defined simple_sbc && !defined minimize_sbc_memory
c
c-----------------------------------------------------------------------
c       call the atmospheric model once for each time step until one
c       segment of "segtim" days is complete. hold atmos S.B.C. fixed
c       during each segment and predict average S.B.C. for MOM
c-----------------------------------------------------------------------
c
        do loop=1,ntspas
# ifdef timing
          call tic ('driver', 'atmos model')
# endif
          asegs   = (loop .eq. 1)
          asege   = (loop .eq. ntspas)
          afirst  = (loop*n .eq. 1)
          alast   = (loop*n .eq. numats)
c
          call release_all
          call atmos
          call iomngr_resume
c
# ifdef timing
          call toc ('driver', 'atmos model')
# endif
        enddo
#endif
c
#ifdef coupled
        write (stdout,9000) 'O C E A N', n, stamp
c
c-----------------------------------------------------------------------
c       get MOM S.B.C.s from the atmosphere
c-----------------------------------------------------------------------
c
        call gosbc (n)
#endif
c
c-----------------------------------------------------------------------
c       call the ocean model once for each time step until one 
c       segment of "segtim" days is complete. hold MOM S.B.C. fixed
c       during each segment and predict average S.B.C. for atmos
c-----------------------------------------------------------------------
c
        do loop=1,ntspos
          call mom 
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     E N D    S E G M E N T    L O O P
c-----------------------------------------------------------------------
c
c
c
c
#ifdef coupled
c
c-----------------------------------------------------------------------
c     save all S.B.C. for restarting coupled models
c-----------------------------------------------------------------------
c
      call getunit (iosbc, 'restrt_sbc.dta', 'u s r')
c
      iotext =  ' read (iosbc) imt, jmt, ima, jma, numsbc'
      write (iosbc) stamp, iotext
      write (iosbc) imt, jmt, ima, jma, numsbc
c
      iotext =
     &  ' read (iosbc) (((sbcatm(i,j,n),i=1,imap2),j=1,jma),n=1,numsbc)'        
      write (iosbc) stamp, iotext
      write (iosbc) sbcatm
c
      iotext =
     &  ' read (iosbc) (((sbcocn(i,j,n),i=1,imt),j=1,jmt),n=1,numsbc)'        
      write (iosbc) stamp, iotext
      write (iosbc) sbcocn
      write (stdout,9200) stamp
c
      call relunit (iosbc)
#endif
c
      write (stdout,9300)
c
c     show files and close all units 
c
      call showfiles
      call release_all
c
c     show timing analysis
c
#ifdef timing
      print '(////a/a//)'
     &,' WARNING: timing analysis significantly increases cpu time'
     &,'          so turn it off when making long integrations'
      call showtimers
#endif
c
#ifdef coupled
8800  format (//,1x,'The model has been configured to run for'
     &,g14.7,' days in ',i4,' segments of ',g14.7,' days each.'
     &,/1x,'There will be ',i6,' ocean time steps per segment'
     &,' and ',i6,' atmosphere time steps per segment.'/
     &,/1x,' The ocean "dtts" =', g14.7, ' seconds'
     &,', and the atmosphere "dtatm" =', g14.7,' seconds'
     &,//)
9000  format (//,1x,'==== E N T E R I N G    ',a9,'    S E G M E N T  #'
     &,       i4, ',  time = ',a32,' ===='//)
9100  format (10x,'==> Finished reading MOM S.B.C. data from file'
     &,    ' restrt_sbc.dta at time = ',a32)
9200  format (10x,'==> Finished writing MOM S.B.C. data to file '
     &,    'restrt_sbc.dta at time = ',a32)
#else
8800  format (//,1x,'The model has been configured to run for'
     &,g14.7,' days in ',i4,' segments of ',g14.7,' days each.'
     &,/1x,'There will be ',i6,' ocean time steps per segment.'
     &,/1x,' The ocean "dtts" =', g14.7, ' seconds'
     &,//)
#endif
9300  format (/,10x,' ==>  M.O.M. integration is complete.')
      stop
      end


      subroutine chkcpl (dtatm, dtts)
      logical errorc
#include "param.h"
#include "csbc.h"
#include "switch.h"
c
c-----------------------------------------------------------------------
c     do consistency checks before allowing model to continue
c-----------------------------------------------------------------------
c
      errorc = .false.
      write (stdout,*) ' '
      write (stdout,*) '    (checking MOM S.B.C. setup)'
c
#ifdef coupled
      if (jma .gt. jmt) then
          write (stdout,9000)
     & '==> Warning: the atmospheric "jma" is greater than "jmt".      '
     &,'    Verify that the atmosphere resoultion is less than the     '
     &,'    ocean resoultion.                                          '
      endif
      if (imap2 .gt. imt) then
          write (stdout,9000)
     & '==> Warning: the atmospheric "imap2" is greater than "imt".    '
     &,'    Verify that the atmosphere resoultion is less than the     '
     &,'    ocean resoultion.                                          '
      endif
      if (bwidth .eq. c0) then
        write (stdout,9000)
     & '==> Note: no blending region was set. Verify that the ocean   '
     &,'    domain is global.                                         '
      else
        write (stdout,9000)
     & '==> Warning: blending region for SST set. verify the blending '
     &,'    zone is correctly set by inspecting the "bzone" printout  '
     &,'    this assumes a limited domain ocean and global atmosphere '
        write (stdout,'(1x,a27,g10.3,a8)') ' Note: "bwidth" is set to '
     &,                                  bwidth, ' degrees'
      endif
#endif
      if (dtatm .eq. c0) then
          write (stdout,9000)
     & '==> Error: the atmospheric time step must be set in "setatm"  '
          errorc = .true.
          dtatm = 1.e-6
      endif
      critv = 1.e-6
      if (segtim .ne. c0) then
        r1 = rundays/segtim
      else
        r1 = 0.5
      endif
      r2 = segtim/(dtts*secday)
      r3 = segtim/(dtatm*secday)
      if (segtim .eq. c0) then
          write (stdout,9000)
     & '==> Error: coupling period "segtim" must be specified when    '
     &,'   the "coupled" option is enabled                            '
          errorc = .true.
      elseif (abs(r1-nint(r1)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of segments       '
     &,'    "segtim"  within "rundays" (the length of the run)        '
          errorc = .true.
      elseif (abs(r2-nint(r2)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of density time   '
     &,'    steps "dtts"  within "segtim" (the segment time)          '
          errorc = .true.
      elseif (abs(r3-nint(r3)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of atmos time     '
     &,'    steps "dtatm"  within "segtim" (the segment time)         '
          errorc = .true.
      endif
#if defined coupled && defined restorst
        write (stdout,9000)
     & '==> Warning: restoring to surface tracers ("restorst" enabled)'
     &, '   when "coupled" is also enabled                            '
#endif
#if defined coupled && defined simple_sbc
        write (stdout,9000)
     & '==> Error: ifdef "simple_sbc" is not compatible               '
     &, '   with "coupled"                                            '
        errorc = .true.
#endif
#ifdef coupled
        if (numasbc .eq. 0) then
          write (stdout,9000) 
     & '==> Warning: "coupled" enabled but number of atmosphere       '
     &,'    boundary conditions "numasbc" is zero!                    '
        endif
#endif
      if (numsbc .gt. maxsbc) then
        write (stdout,9000) 
     & '==> Error: "numsbc" is greater than "maxsbc" in "csbc.h"      '
        errorc = .true.
      endif
      write (stdout,*) '    (End of MOM S.B.C. checks) '
      write (stdout,*) ' '
      if (errorc) stop '=>chkcpl'
c
9000  format (/,(1x,a80))
      return
      end


      subroutine getst (jrow, ocnout, ntabc)
c
c-----------------------------------------------------------------------
c     read surface tracers from disk row "jrow"
c-----------------------------------------------------------------------
c
#include "param.h"
#include "iounit.h"
#include "mw.h"
#include "tmngr.h"
c
      dimension ocnout(imt,jmt)
c
      call getrow (latdisk(taup1disk), nslab, jrow
     &,          u(1,1,jmw,1,taup1), t(1,1,jmw,1,taup1))
      do i=1,imt
        ocnout(i,jrow) = t(i,1,jmw,ntabc,taup1)
      enddo
      return
      end

