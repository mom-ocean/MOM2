c=======================================================================
c
c               T I M E      M A N A G E R   M O D U L E
c
c     The time manager does three basic things for the model
c           1) It keeps all times for the model in arrays in "tmngr.h"
c           2) It updates times using the chosen time step "dt" and
c              calendar type without any roundoff or drift for
c              arbitrarily long integrations, even when using
c              32 bit arithmetic.
c           3) It sets switches in "switch.h" that trigger periodically
c              recurring events in the model such as diagnostics
c              and end-of-run.
c
c     Times are kept internal to the time manager in the form of
c     integer days (in the array "iday") and integer milliseconds
c     in a fractional day (in the array "msday").  A unique 
c     subscript for each time accesses all the arrays.  (See "tmngr.h"
c     for a list of predeclared time subscripts.) 
c
c     Times that have a calendar reference are called "full times"
c     because in addition to the primary time fields "iday" and
c     "msday", they also have integer secondary time fields "year",
c     "month", "day", "hour", "minute", and "second", and a 32
c     character time stamp, "tstamp", in the form written and read
c     by the model. They also have tertiary fields "dayofyear",
c     "dayofweek", "daysinmon", and "daysinyear".  For such times,
c     the primary fields "iday" and "msday" represent days and
c     fractional days since the calendar base, (0,0) =
c     December 31, 1899 at the start of the day, 00:00:00.
c
c     For differences of two calendar times, these calendar fields
c     make no sense, so these are kept as "short times" with only
c     the two primary fields "iday" and "msday".
c
c     tmngr can keep track of time using either
c           a) a leap-year corrected Gregorian calendar
c           b) a constant 365-day year calendar
c           c) a constant month calendar, 12 months/year, with an
c              arbitrary number of days per month
c
c     Although the model does not currently use this feature, the
c     subroutine "tmngr" is designed to handle time increments, "dt",
c     varying at every time step.
c
c     Currently implemented switches include end-of-day, end-of-
c     week, end-of-two-weeks, end-of-month, end-of-year, end-of-
c     run, mid-month, and switches active at prespecified intervals
c     from either start of run, initial conditions, or any other
c     reference time of the user's choosing.  It is relatively
c     easy to add additional switches by following the models of
c     switches already provided.
c
c     Finally, this module provides a collection of utility 
c     subroutines to perform arithmetic on the internal representations
c     of time, to convert between primary and secondary internal
c     representations, and between internal representations and
c     external representations such a real seconds or real days.     
c
c     authors: d. j. goldberg
c              c. h. goldberg         chg@gfdl.gov
c              (an improved version of the time manager from MOM 1)
c           
c=======================================================================
c


#ifdef test_tmngr
# define snapshots
# define energy_analysis
# define xbts
# define time_step_monitor
# include "util.F"
c
      program driver
c
c=======================================================================
c     To test the time manager module, execute the "run_tmngr" script
c
c     author:  c. h. goldberg         chg@gfdl.gov          
c=======================================================================
c
#include "stdunits.h"
#include "calendar.h"
#include "tmngr.h"
#include "switch.h"
c
c-----------------------------------------------------------------------
c     initialize all quantities that the model would be supplying to
c     the time manager
c-----------------------------------------------------------------------
c     set type of year (eqyear = F => real calendar with leap years)
c                      (eqyear = T => idealized calendar)
c
c     if calendar is idealized then there are two choices:
c     eqmon = F => the usual number of days per month (31, 28, 31, ...)
c                  and the length of each year is 365 days
c                  (monlen is not used in this case)
c     eqmon = T => a fixed number of days per month set by "monlen"
c                  and the length of a year is 12*monlen
c-----------------------------------------------------------------------
c
      eqyear = .true.
      eqmon  = .false.
      monlen = 30
c
c-----------------------------------------------------------------------
c     initialize clock to time at initial conditions
c-----------------------------------------------------------------------
c
      year0  = 1
      month0 = 1
      day0   = 1
c
      hour0  = 0
      min0   = 0
      sec0   = 0
c
c-----------------------------------------------------------------------
c     decide on reference for calculating logical switches
c
c     refrun  = T ==> reference to the start of each run
c     refinit = T ==> reference to initial condition time given by:
c                     year0, month0, day0, hour0, min0, sec0
c     refuser = T ==> reference to user specified time given by:
c                     ryear, rmonth, rday, rhour, rmin, rsec
c
c     choose by setting one of the above to TRUE
c-----------------------------------------------------------------------
c
      refrun  = .true.
      refinit = .false.
      refuser = .false.
c
c     set the date and time for referencing switches (if refuser = .T.)
c
      if (refuser) then
        ryear  = 1900
        rmonth = 1
        rday   = 1
        rhour  = 0
        rmin   = 0
        rsec   = 0
      end if
c
c-----------------------------------------------------------------------
c     set the time step "dt" (sec) and how long to integrate "runlen".
c     (note: the tmngr can accomodate varying "dt")
c     "rununits" can either be 'days', 'years', or 'months'. within 
c     tmngr, "runlen" will be automatically converted to units of days.
c-----------------------------------------------------------------------
c
      dt       = 3600.0
      runlen   = 1.2*365.0
      rununits = 'days'
c
c     the time manager can also handle run length specifications of the
c     following sorts. for example, to run 2 years:
c
c     runlen   = 2
c     rununits = 'years'
c            or
c     runlen   = 2*12.0
c     rununits = 'months'
c
c
c-----------------------------------------------------------------------
c     for use in MOM, the integration time "runlen" is broken into a
c     number of segments each of length "segtim" days. updated surface
c     boundary conditions are applied to MOM every "segtim" days. this
c     is useful when coupling to atmospheric models in which case
c     both models exchange surface boundary conditions every "segtim"
c     days where "segtim" is 1/(coupling frequency). when MOM is without
c     an atmospheric model, "segtim" is set to the time step (in days)
c     as done here.     
c-----------------------------------------------------------------------
c
      segtim = dt/86400.0
c
c-----------------------------------------------------------------------
c     set how often output is desired (units are in days)
c     (the switches and corresponding logicals are in "switch.h")
c-----------------------------------------------------------------------
c
      tsiint    = 10.0
      glenint   = 365.0
      snapint   = 365.0*1.25
      xbtint    = 90.0
      nmix      = 17
c
      write (stdout,'(///,tr5,a,g14.7,a)')
     & 'Driver will simulate integrating', runlen, rununits
      print *, '     using timesteps of dt = ',dt,' (sec) '
      print '(/,5x,a,/)'
     &,'the following diagnostic switches will be tested:'
      print *, '     time step integrals every ',tsiint,' days'
      print *, '     global energetics every ',glenint,' days'
      print *, '     snapshots every ',snapint,' days'
      print *, '     time averaged xbt''s every ',xbtint, ' days'
c
c-----------------------------------------------------------------------
c     initialize the time step counter "itt" and 
c     set the corresponding model time counters in the
c     form (integer day,integer millisecond) = (idayrestart,msrestart)
c     in MOM, this information is read from the restart file. if
c     starting from initial conditions, set idayrestart=msrestart=0.
c-----------------------------------------------------------------------
c
      itt          = 0
      idayrestart  = 0
      msrestart    = 0
c
c-----------------------------------------------------------------------
c     initialize the time manager with specified initial condition time,
c     user reference time, and model time.
c-----------------------------------------------------------------------
c
      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,              ryear, rmonth, rday, rhour, rmin, rsec
     &,              idayrestart, msrestart
     &,              runlen, rununits, rundays, dt)
c
c-----------------------------------------------------------------------
c     integrate the time manager for the required time
c     and print when things will happen
c-----------------------------------------------------------------------
c
      do n=1,1000000
c
        itt = itt + 1
        call tmngr (dt)
c
c-----------------------------------------------------------------------
c       calculate days of model integration since I.C.
c-----------------------------------------------------------------------
c
        daymod = realdays(imodeltime)
c
c-----------------------------------------------------------------------
c       print when things happen
c-----------------------------------------------------------------------
c
        if (itt .lt. 10) then
          write (*,9000) itt, stamp, daymod
        endif
c
        if (eomon) print 9001,' ...end of month on itt=', itt, ' ',
     &                     stamp, dayname(dayofweek(itime))
        if (midmon) print 9001,' ... mid month on itt  =', itt, ' ',
     &                      stamp, dayname(dayofweek(itime))
        if (eoyear) print 9001,' ==> end of year on itt=',itt,' ',
     &                      stamp, dayname(dayofweek(itime))
c
        if (xbtts) print 9001
     &,' --->switch activated to save averaged xbts on itt= '
     &,      itt, ' ',stamp, dayname(dayofweek(itime))
c
        if (glents) print 9001
     &,' --->switch activated to do energy analysis on itt= '
     &,      itt, ' ',stamp, dayname(dayofweek(itime))
c
        if (snapts) print 9001
     &,' --->switch activated to take snapshot on itt= '
     &,      itt, ' ',stamp, dayname(dayofweek(itime))
c
        if (tsits) print 9001
     &,' --->switch activated to do time step integrals on itt= '
     &,      itt, ' ',stamp, dayname(dayofweek(itime))
c
        if (eorun)  then
	  print 9001,' ==> end of run on itt=',itt,' ', stamp
     &,               dayname(dayofweek(itime))
          print '(/,10x,a,/)','Can your Rolex do this?'
          stop
	endif
c
      enddo
      print *,' Increase loop limit on "n" in "driver.F"'
      stop 'driver'
9000  format (1x, ' itt=',i7,1x,a,' model days since I.C.=',f12.4)
9001  format (a, i8, a, a, 1x, a)
      end
#endif








      subroutine tmngri (icyear, icmonth, icday, ichour, icmin, icsec
     &,                     rfyear, rfmonth, rfday, rfhour, rfmin, rfsec
     &,                     idayrestart, msrestart
     &,                     runlen0, rununits0, rundays0, timestep)
c
c=======================================================================
c     initialize internal time variables and enter initial values
c     for externally defined times (initial conditions, user reference
c     and model time at start of run).
c=======================================================================
c
#include "stdunits.h"
#include "switch.h"
#include "tmngr.h"
#include "calendar.h"
      character*(*) rununits0
c
      integer icyear, icmonth, icday, ichour, icmin, icsec
      integer rfyear, rfmonth, rfday, rfhour, rfmin, rfsec
      integer idayrestart, msrestart
      logical error, timeless
c
      write (stdout,'(//,10x,a,/)') 'Time manager initialization'
c
      call inittime
      first = .true.
      error = .false.
      call calendari (eqyear, eqmon, monlen,
     &            yrlen, daypm, msum, dayname, monname, error)
      if (error) then
        stop 'badcal'
      end if
c
c-----------------------------------------------------------------------
c     enter and print initial date and time and check bounds.
c-----------------------------------------------------------------------
c
      call getfulltime (initial)
      call setfulltime (initial, icyear, icmonth, icday, ichour, icmin
     &,                 icsec)
      write (stdout,9000) 'Initial Conditions: ', tstamp(initial)
      write (stdout,*) ' '
      call ckdate (initial, error)
      if (error) then
        stop '=>tmngri'
      end if
c
c----------------------------------------------------------------------
c     set model time counter.  
c----------------------------------------------------------------------
c
      call gettime  (imodeltime)
      call settime2 (imodeltime, idayrestart, msrestart)
      write (stdout,'(a,1pg14.7,a/)') 
     & '  Time since initial conditions =',
     & realdays(imodeltime), ' days'
c
c----------------------------------------------------------------------
c     calculate y/m/d and h/m/s of start of run.
c----------------------------------------------------------------------
c
      call getfulltime (irunstart)
      call addtime (initial, imodeltime, irunstart)
      call expandtime2 (irunstart)
      call getfulltime (itime)
      call copyfulltime (irunstart, itime)
      stamp = tstamp(itime)
      pstamp = stamp
c
c-----------------------------------------------------------------------
c     calculate real output quantities:
c       relyr  = years of model time since initial conditions
c       dayoyr = days since start of current year
c-----------------------------------------------------------------------
c
      dayoyr = dayofyear(itime) - 1 + (msday(itime)/86400000.0)
      relyr  = year(itime) - year(initial) + dayoyr/daysinyear(itime) -
     &         (dayofyear(initial) - 1 + (msday(initial)/86400000.0))/
     &         daysinyear(initial)
      prelyr = relyr
c
#ifdef debug_tmngr
      print *, 'startms     = ',msday(irunstart)
      print *, ' initday    = ',iday(initial)
      print *, 'startday    = ',iday(irunstart)
      print *, 'idayrestart = ',iday(imodeltime)
#endif
      write (stdout,9000) 'Start of Run: ', tstamp(irunstart)
      write (stdout,'(a,i10/)') 'Corresponding to time step "itt" =',itt
c
c----------------------------------------------------------------------
c       select reference time for computing diagnostic switches
c----------------------------------------------------------------------
c
      irefs = 0
      call getfulltime (iref)
      if (refrun) then
        irefs = irefs + 1
        write (stdout,*)
     &  ' "refrun = .true." selected.'
        write (stdout,'(a,a)')
     &  '  intervals for diagnostic switches are referenced to ',
     &  '  the beginning of each run.'
        call copyfulltime (irunstart, iref)
      endif
      if (refinit) then
        irefs = irefs + 1
        write (stdout,*)
     &  ' "refinit = .true." selected.'
        write (stdout,'(a,a)')
     &  '  intervals for diagnostic switches are referenced to ',
     &  '  the initial conditions time.'
        call copyfulltime (initial, iref)
      endif
      if (refuser) then
        irefs = irefs + 1
        write (stdout,*)
     &  ' "refuser = .true." selected.'
        write (stdout,'(a,a)')
     &  '  intervals for diagnostic switches are referenced to ',
     &  '  user specified date and time.'
        call getfulltime (iuser)
        call setfulltime (iuser, rfyear, rfmonth, rfday, rfhour, rfmin
     &,                   rfsec)
        call ckdate (iuser, error)
        if (error) then
          stop '=>tmngri'
        endif
        write (stdout,9000) '    Reference time:     ', tstamp(iuser)
        write (stdout,*) ' '
        call copyfulltime (iuser, iref)
      endif
c
      if (irefs .ne. 1) then
        write (stdout, *) 'You must choose exactly one of the'
     &  // ' options: refrun, refinit, or refuser.'
        stop '=>tmngr'
      endif
c
      call gettime (idt)
      call gettime (idtd2)
      call gettime (iusertime)
      call gettime (ireftime)
      call gettime (iruntime)
      call gettime (ihalfstep)
      call gettime (itimeper)
      call getfulltime (itemptime)
      call getfulltime (itemptime2)
c
c-----------------------------------------------------------------------
c     set a reference to Sunday Jan 0, 1900, 0:00:00, the base date
c     or earier if necessary.
c-----------------------------------------------------------------------
c
      call gettime  (isunday)
      call settime2 (isunday, 0, 0)
      if (timeless (initial, isunday)) then
        iday(isunday) = 14*int(iday(initial)/14) - 14
      end if
c
c-----------------------------------------------------------------------
c     initialize idt as "timestep" rounded to the nearest millisecond
c-----------------------------------------------------------------------
c
      iddt = id (timestep)
      msdt = msec (timestep)
      call settime2 (idt, iddt, msdt)
c
      call set_eorun (runlen0, rununits0, rundays0)
c
      write (stdout,'(/,10x,a,//)') 'Initialization completed'
c
9000  format (a,a)
      return
      end





      subroutine tmngr (dt)
c
c=======================================================================
c
c               T I M E      M A N A G E R   M O D U L E
c
c     tmngr keeps track of model time using either
c           a) leap-year corrected Gregorian calendar
c           b) constant 365-day year calendar
c           c) constant month calendar, 12 months/year, arbitrary
c              number of days per month
c
c     date and time are accurate to within one millisecond for arbitrary
c     length time steps (even on computers with 32 bit word lengths!).
c
c               julian base = Jan 0, 1900 at 00:00:00
c
c       For accuracy, all fundamental times are kept
c       in the form:  integer days (with Jan 1, 1900=1)
c                     non-negative integer milliseconds within the day
c
c     input:
c
c       dt     = length of time step in seconds. (need not be constant)
c
c       a switch.h file with all switch intervals in units of days.
c
c     outputs:
c
c       dt     = dt rounded to nearest millisec if needed
c
c       updated time fields: year, month, day, hour, minute, second,
c                            tstamp, dayofyear, dayofweek, daysinmonth,
c                            daysinyear
c       times updated:
c         itime  = "absolute" (y/m/d...) time after adding dt
c         ihalfstep = dt/2 beyond itime
c         imodeltime = time since initial conditions
c
c       stamp   = 32 character time stamp (m/d/y h:m:s)
c       pstamp  = stamp returned in previous call to tmngr
c
c       relyr  = model time in years since initial conditions
c       dayoyr = days + fractional days since start of calendar year
c       plus a whole collection of useful logical switches in the
c       switch.h file.
c
c     authors: d. j. goldberg
c              c. h. goldberg        chg@gfdl.gov
c             (an improved version of the time manager
c              from MOM by r.c. pacanowski and k. dixon)
c           
c=======================================================================
c
      logical alarm, timeless, timeequal
c
#include "stdunits.h"
#include "tmngr.h"
#include "switch.h"
c
c=======================================================================
c     date and time calculations (done every time step)
c=======================================================================
c
c     set flag "first" if first iteration of a run
c
      first = timeequal (itime, irunstart)
c
c-----------------------------------------------------------------------
c     round dt to the nearest millisecond
c-----------------------------------------------------------------------
c
      iddt = id (dt)
      msdt = msec (dt)
      call settime2 (idt, iddt, msdt)
      dt   = realsecs(idt)
c
c-----------------------------------------------------------------------
c     calculate half time step
c-----------------------------------------------------------------------
c
      iddtd2 = iddt/2
      msdtd2 = msdt/2 + 43200000*modulo (iddt, 2)
      call settime2 (idtd2, iddtd2, msdtd2)
c
c-----------------------------------------------------------------------
c     save previous values of stamp and relyr
c-----------------------------------------------------------------------
c
      pstamp = tstamp(itime)
      prelyr = relyr
c
c-----------------------------------------------------------------------
c     increment time counters
c-----------------------------------------------------------------------
c
      call addtime (itime, idt, itime)
      call expandtime2 (itime)
c
c     set current time stamp for MOM
c
      stamp = tstamp(itime)
c
      call addtime (itime, idtd2, ihalfstep)
c
c-----------------------------------------------------------------------
c     calculate number of days since reference time and start of run
c     all times are of form: (integer days, fractional day in millisec)
c----------------------------------------------------------------------
c
      call subtime (itime, irunstart, iruntime)
      call subtime (itime, initial, imodeltime)
      if (refuser) then 
        call subtime (itime, iuser, iusertime)
      end if
c
c-----------------------------------------------------------------------
c     calculate real output quantities:
c       relyr  = years of model time since initial conditions
c       dayoyr = days since start of current year
c-----------------------------------------------------------------------
c
      dayoyr = dayofyear(itime) - 1 + (msday(itime)/86400000.0)
      relyr  = year(itime) - year(initial) + dayoyr/daysinyear(itime) -
     &         (dayofyear(initial) - 1 + (msday(initial)/86400000.0))/
     &         daysinyear(initial)
c
c=======================================================================
c     set all time dependent logical switches in this section.
c=======================================================================
c
c-----------------------------------------------------------------------
c      here are some examples of setting logicals based on the calendar
c      and clock. "dayoyr" is time in days since the start of a year.
c
c      alarms go off when current time + dt/2 is later then the alarm
c      time.  arguments of alarm routine:
c            1) index to the alarm setting time
c            2) index to the calculated time+dt/2
c            3) period of the alarm in real days
c            4) index to the reference time
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of the day ?
c-----------------------------------------------------------------------
c
      eoday = alarm (ieoday, ihalfstep, 1.0, isunday)
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of the week ?
c-----------------------------------------------------------------------
c
      eoweek = alarm (ieoweek, ihalfstep, 7.0, isunday)
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of two weeks ?
c-----------------------------------------------------------------------
c
      eo2wks = alarm (ieo2wks, ihalfstep, 14.0, isunday)
c
c-----------------------------------------------------------------------
c     set initial alarm times in the first time step for the
c     end of month, mid month, end of year, and end of run switches.
c-----------------------------------------------------------------------
c
      if (first) then
c
c-----------------------------------------------------------------------
c       initialize end of month alarm
c-----------------------------------------------------------------------
c
        call getfulltime (ieomon)
        call copyfulltime (itime, ieomon)
        month (ieomon) = month (ieomon) + 1
        if (month(ieomon) .gt. 12) then
          month(ieomon) = 1
          year (ieomon) = year(ieomon) + 1
        end if
        day   (ieomon) = 1
        hour  (ieomon) = 0
        minute(ieomon) = 0
        second(ieomon) = 0
        call expandtime (ieomon)
c
c-----------------------------------------------------------------------
c       initialize mid month alarm
c-----------------------------------------------------------------------
c
        call getfulltime (imidmon)
        call copyfulltime (itime, imidmon)
        day   (imidmon) = daysinmon(imidmon)/2 + 1
        hour  (imidmon) = 12*modulo(daysinmon(imidmon), 2)
        minute(imidmon) = 0
        second(imidmon) = 0
        call expandtime (imidmon)
        if (timeless (imidmon, itime)) then
          month (imidmon) = month(imidmon) + 1
          if (month(imidmon) .gt. 12) then
            month(imidmon) = 1
            year (imidmon) = year(imidmon) + 1
          end if
          day   (imidmon) = daysinmon(imidmon)/2 + 1
          hour  (imidmon) = 12*modulo(daysinmon(imidmon), 2)
          minute(imidmon) = 0
          second(imidmon) = 0
          call expandtime (imidmon)
        end if
c
c-----------------------------------------------------------------------
c       initialize end of year alarm
c-----------------------------------------------------------------------
c
        call getfulltime (ieoyear)
        call copyfulltime (itime, ieoyear)
        year  (ieoyear) = year(ieoyear) + 1
        month (ieoyear) = 1
        day   (ieoyear) = 1
        hour  (ieoyear) = 0
        minute(ieoyear) = 0
        second(ieoyear) = 0
        call expandtime (ieoyear)
c
      end if
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of the month ?
c-----------------------------------------------------------------------
c
      eomon = timeless (ieomon, ihalfstep)
      if (eomon) then
        iday(ieomon) = iday(ieomon) + daysinmon(ieomon)
        call expandtime2 (ieomon)
      end if
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the middle of the month ?
c-----------------------------------------------------------------------
c
      midmon = timeless (imidmon, ihalfstep)
      if (midmon) then
        day(imidmon) = 1
        month(imidmon) = month(imidmon) + 1
        if (month(imidmon) .gt. 12) then
          month(imidmon) = 1
          year(imidmon)  = year(imidmon) + 1
        end if
        call expandtime (imidmon)
        day   (imidmon) = daysinmon(imidmon)/2 + 1
        hour  (imidmon) = 12*modulo(daysinmon(imidmon), 2)
        minute(imidmon) = 0
        second(imidmon) = 0
        call expandtime (imidmon)
      end if
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of the year ?
c-----------------------------------------------------------------------
c
      eoyear = timeless (ieoyear, ihalfstep)
      if (eoyear) then
        iday(ieoyear) = iday(ieoyear) + daysinyear(ieoyear)
        call expandtime2 (ieoyear)
      endif
c
c-----------------------------------------------------------------------
c     is it the last time step of the run ?
c-----------------------------------------------------------------------
c
      eorun = timeless (ieorun, ihalfstep)
c
c-----------------------------------------------------------------------
c     is time mixing to be done now ?
c-----------------------------------------------------------------------
c
      if (nmix .eq. 0 .or. nmix .eq. 1) then
        leapfrog = .true.
      else
        if (mod(itt,nmix) .eq. 1) then
          leapfrog = .false.
        else
          leapfrog = .true.
        endif
      endif
c
c-----------------------------------------------------------------------
c     the following are all done when within 1/2 time step of the
c     requested time interval in days. an interval < 0 turns off
c     the switch. an interval of zero turns on the switch every time
c     step.
c-----------------------------------------------------------------------


c
c-----------------------------------------------------------------------
c     is it the first time step of an ocean segment?
c     is it the last time step of an ocean segment?
c-----------------------------------------------------------------------
c
      if (first) then
        osegs = .true.
      else
        osegs = osege
      endif
c
      osege = alarm (iosege, ihalfstep, segtim, irunstart)
c
#ifdef tracer_averages
c
c-----------------------------------------------------------------------
c     are regional tracer averages to be done at this time ?
c-----------------------------------------------------------------------
c
      tavgts = alarm (iitavgint, ihalfstep, tavgint, iref)
#endif
c
c
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     is saving of estimated mixing coefficients to be done at
c     this time ?
c-----------------------------------------------------------------------
c
      cmixts = alarm (icmixint, ihalfstep, cmixint, iref)
#endif
#ifdef save_convection
c
c-----------------------------------------------------------------------
c     is saving of temperature change due to explicit convection
c     to be done at this time ?
c-----------------------------------------------------------------------
c
      exconvts = alarm (iexconvint, ihalfstep, exconvint, iref)
#endif
#ifdef energy_analysis
c
c-----------------------------------------------------------------------
c     are global energetics to be done at this time ?
c-----------------------------------------------------------------------
c
      glents = alarm (iglenint, ihalfstep, glenint, iref)
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     are tracer & momentum term balances to be done at this time ?
c-----------------------------------------------------------------------
c
      trmbts = alarm (iitrmbint, ihalfstep, trmbint, iref)
#endif
c
#ifdef matrix_sections
c
c-----------------------------------------------------------------------
c     are x-z printouts to be done at this time ?
c-----------------------------------------------------------------------
c
      prxzts = alarm (iprxzint, ihalfstep, prxzint, iref)
#endif
c
#ifdef show_external_mode
c
c-----------------------------------------------------------------------
c     is the external mode printout to be done at this time ?
c-----------------------------------------------------------------------
c
      extts = alarm (iextint, ihalfstep, extint, iref)
#endif
c
#ifdef meridional_overturning
c
c-----------------------------------------------------------------------
c     is time to write the vertical_meridional stream function?
c-----------------------------------------------------------------------
c
      vmsfts = alarm (ivmsfint, ihalfstep, vmsfint, iref)
#endif
c
#ifdef gyre_components
c
c-----------------------------------------------------------------------
c     are the gyre components to be done at this time ?
c-----------------------------------------------------------------------
c
      gyrets = alarm (igyreint, ihalfstep, gyreint, iref)
#endif
c
#ifdef diagnostic_surf_height
c
c-----------------------------------------------------------------------
c     is the diagnostic surf height calculation to be done now ?
c-----------------------------------------------------------------------
c
      dspts = alarm (idspint, ihalfstep, dspint, iref)
#endif
c
#ifdef time_step_monitor
c
c-----------------------------------------------------------------------
c     is printing of time step integrals to be done at this time ?
c-----------------------------------------------------------------------
c
      tsits = alarm (itsiint, ihalfstep, tsiint, iref)
#endif
c
#ifdef snapshots
c
c-----------------------------------------------------------------------
c     is time to write a snapshot of the MOM data ?
c-----------------------------------------------------------------------
c
      snapts = alarm (isnapint, ihalfstep, snapint, iref)

#endif
c
#ifdef trajectories
c
c-----------------------------------------------------------------------
c     is it time to write particle trajectories?
c     they must be written at the end of the run for restarting.
c-----------------------------------------------------------------------
c
      trajts = alarm (itrajint, ihalfstep, trajint, iref)
      trajts = trajts .or. eorun
#endif
c
#ifdef time_averages
c
c-----------------------------------------------------------------------
c     is it time to write time mean data on the "averaging" grid?
c-----------------------------------------------------------------------
c
      timats = alarm (itimavgint, ihalfstep, timavgint, iref)
      timats = timats .or. eorun
#endif
c
#ifdef xbts
c
c-----------------------------------------------------------------------
c     is it time to write time mean data on the "averaging" grid?
c-----------------------------------------------------------------------
c
      xbtts = alarm (iixbtint, ihalfstep, xbtint, iref)
      xbtts = xbtts .or. eorun
#endif
c
#ifdef meridional_tracer_budget
c
c-----------------------------------------------------------------------
c     is it time to write the meridional tracer budget?
c-----------------------------------------------------------------------
c
      tmbts = alarm (iitmbint, ihalfstep, tmbint, iref)
      tmbts = tmbts .or. eorun
#endif
c
#ifdef stability_tests
c
c-----------------------------------------------------------------------
c     is it time to do a stability diagnosis ?
c-----------------------------------------------------------------------
c
      stabts = alarm (istabint, ihalfstep, stabint, iref)
#endif
c
#ifdef show_zonal_mean_of_sbc
c
c-----------------------------------------------------------------------
c      is it time to do a boundary condition diagnosis ?
c-----------------------------------------------------------------------
c
      zmbcts = alarm (izmbcint, ihalfstep, zmbcint, iref)
#endif


c
      return
      end



      subroutine set_eorun (runlen0, rununits0, rundays0)
c=======================================================================
c     initialize end of run alarm.
c     place dummy arguments in corresponding variables in "switch.h"
c
c     inputs:
c       runlen0   = length of run [see rununits for units]
c       rununits0 = time units ('days', 'months', or 'years' for runlen0
c     output:
c       rundays0  = length of run converted to (real) days
c
c     at present, arbitrary real run lengths in 'days' are handled.
c     run lengths in real 'months' or 'years' are rounded to the
c     nearest integral number of months or years before use.
c
c     author:  d. j. goldberg          chg@gfdl.gov
c=======================================================================
c
#include "switch.h"
#include "stdunits.h"
#include "tmngr.h"
      character *(*) rununits0
c
c     copy inputs to common in "switch.h"
c
      runlen = runlen0
      rununits = rununits0
c
c     calculate end of run based on rununits0
c
      call getfulltime (ieorun)
      call copyfulltime (irunstart, ieorun)
      if (rununits0 .eq. 'years') then
        year(ieorun) = year(ieorun) + nint(runlen0)
        call expandtime (ieorun)
        day  (ieorun) = min (day(irunstart), daysinmon(ieorun))
        call expandtime (ieorun)
        if (real(nint(runlen0)) .ne. runlen0) then
          print '(/,a,1pg14.7,a,i3,a,/)', 
     &          'WARNING:  run length ',runlen0,' years rounded to ',
     &          nint(runlen0), ' years'
          runlen0 = real(nint(runlen0))
          runlen = runlen0
        end if
      elseif (rununits0 .eq. 'months') then
        month(ieorun) = month(ieorun) + nint(runlen0) - 1
        year (ieorun) = year (ieorun) + month(ieorun)/12
        month(ieorun) = modulo (month(ieorun), 12) + 1
        day  (ieorun) = 1
        call expandtime (ieorun)
        if (real(nint(runlen0)) .ne. runlen0) then
          print '(/,a,1pg14.7,a,i3,a,/)', 
     &          'WARNING:  run length ',runlen0,' months rounded to ',
     &          nint(runlen0), ' months'
          runlen0 = real(nint(runlen0))
          runlen = runlen0
        end if
        day  (ieorun) = min (day(irunstart), daysinmon(ieorun))
        call expandtime (ieorun)
      elseif (rununits0 .eq. 'days') then
        call settime3 (itemptime, runlen0)
        call addtime (irunstart, itemptime, ieorun)
        call expandtime2 (ieorun)
      else
        write (stdout, *)
     &       'Warning: No units given for run length--days assumed'
        call settime3 (itemptime, runlen0)
        call addtime (irunstart, itemptime, ieorun)
        call expandtime2 (ieorun)
      end if
c
      call subtime (ieorun, irunstart, itemptime)
      rundays0 = realdays (itemptime)
      rundays  = rundays0
c
      return
      end



      subroutine calendari (eqyear, eqmon, monlen,
     &                   yrlen, daypm, msum, dayname, monname, error)
c
c=======================================================================
c     set up the calendar by choosing one of the following:
c
c        a) fully leap-year corrected calendar     (eqyear=F eqmon=F)
c        b) equal 365-day years (variable months)  (eqyear=T eqmon=F)
c        c) 12 equal months of "monlen" days each  (eqyear=T eqmon=T)
c
c     inputs:
c        eqyear, eqmon = logicals degined as in "a", "b", and "c"
c        monlen = days per month when using option "c"
c
c     outputs:
c        yrlen   = length of year in days
c        daypm   = days per month
c        msum    = accululated days per month
c        dayname = character day names
c        monname = character month names
c        error   = .t. if something went wrong
c
c     author:  d. j. goldberg          chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      parameter (nmonth = 12, nday = 7)
      logical eqyear, eqmon, error
      character*10 daynamei
      character*12 monnamei
      character*(*) dayname, monname
      dimension daynamei(nday), monnamei(nmonth)
      dimension dayname(nday), monname(nmonth)
      integer daypm(nmonth), daypmi(nmonth), msum(nmonth), yrlen
      data daypmi/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daynamei/'sunday', 'monday', 'tuesday', 'wednesday',
     &            'thursday', 'friday', 'saturday'/
c
      data monnamei /'january', 'febuary', 'march', 'april', 'may'
     &,           'june', 'july', 'august', 'september', 'october'
     &,           'november', 'december'/
c
c-----------------------------------------------------------------------
c     first check that a valid combination of eqyear and eqmon
c     has been specified
c-----------------------------------------------------------------------
c
      if (eqmon .and. (.not. eqyear)) then
        write (stdout,999) eqyear, eqmon
        stop 'badcal'
      endif
c
      if (.not. (eqyear .and. eqmon)) then
        do i = 1, nmonth
          daypm(i) = daypmi(i)
        enddo
      else
        error = error .or. (monlen .le. 0)
        do i = 1, nmonth
          daypm(i) = monlen
        enddo
      endif
c
      msum(1) = 0
      do i = 2,nmonth
        msum(i) = msum(i-1) + daypm(i-1)
      enddo
      yrlen = msum(nmonth) + daypm(nmonth)
c
c     initialize day and month names
c
      do n=1,nday
        dayname(n) = daynamei(n)
      enddo
c
      do n=1,nmonth
        monname(n) = monnamei(n)
      enddo
c
      write (stdout,'(a)') 'Calendar selection:'
      if (.not. eqyear) then
        print '(/,a)', '  leap year corrected calendar'
      else
        print '(/,a,i4,a)', '  equal years of ', yrlen, ' days'
      endif
c
      print '(a,12i3,/)', '  days per month =',daypm
c
      return
c
999   format(/' error => An inappropriate calendar type was selected  ',
     &   /'     eqyear was set to ',l1,'   eqmon was set to ',l1,
     &   /'     valid combinations are:',
     &   /'     fully leap-year corrected calendar  (eqyear=F eqmon=F)'
     &   /'     equal 365-day years                 (eqyear=T eqmon=F)'
     &   /'     12 equal months of monlen days each (eqyear=T eqmon=T)')
c
      end



      subroutine ckdate (i, error)
c
c=======================================================================
c     do bounds checking on clock parameters
c     year is not checked since all years are ok.
c     at present, one extra day is allowed in month 2, even if year is
c     not a leap year or if there are no leap years.
c
c     author:  d. j. goldberg          chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
#include "calendar.h"
      logical error
c
      error = .false.
c
      if (month(i) .lt. 1 .or. month(i) .gt. 12) then
        write (stdout,*) ' Error:  month is out of bounds'
        error = .true.
      endif
c 
      if ((day(i) .lt. 1 .or. day(i) .gt. daypm(month(i))) .and. .not.
     &    (month(i) .eq. 2 .and. day(i) .eq. daypm(2)+1)) then
        write (stdout,*) ' Error:  day is out of bounds'
        error = .true.
      endif
c
      if (hour(i) .lt. 0 .or. hour(i) .gt. 23) then
        write (stdout,*) ' Error:  hour is out of bounds'
        error = .true.
      endif
c
      if (minute(i) .lt. 0 .or. minute(i) .gt. 59) then
        write (stdout,*) ' Error:  minute is out of bounds'
        error = .true.
      endif
c
      if (second(i) .lt. 0 .or. second(i) .gt. 59) then
        write (stdout,*) ' Error:  second is out of bounds'
        error = .true.
      endif
c
      return
      end


      subroutine d2ymd (inday,
     &                  iyear, imonth, iday, idoy, idow, ndim, ndiy)
c
c======================================================================
c     d2ymd takes the number of days since Dec 31, 1899 and converts
c     to year, month, day, and other defining quantities.
c     ie. inday=1   yields  Jan 1, 1900.
c     if eqyear=.true. then calculation is done in subroutine d2ymdc
c
c     inday  = input days since Dec. 31, 1899
c     eqyear = .false. ==> use leap year corrected calendar
c     eqyear = .true.  ==> use constant length year calendar
c
c     output:
c       iyear  = output year
c       imonth = output month
c       iday   = output day of the month
c       idoy   = output day of the current year
c       idow   = output day of the week   1=sun - 7=sat
c       ndim   = number of days in the current month
c       ndiy   = number of days in the current year
c
c     author:  d. j. goldberg          chg@gfdl.gov
c======================================================================
c
#include "calendar.h"
c
      logical leap
c
c     define the number of days per month for a non leap year (daypmi)
c
      integer daypmi(12)
      data daypmi/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c
c-----------------------------------------------------------------------
c     equal-year calculations are handled by d2ymdc
c-----------------------------------------------------------------------
c
      if (eqyear) then
        call d2ymdc (inday,
     &               iyear, imonth, iday, idoy, idow, ndim, ndiy)
        return
      endif
c
c-----------------------------------------------------------------------
c     calculate day of the week (idow)
c     inday = 0, i.e.,  Dec 31, 1899 was a Sunday.
c-----------------------------------------------------------------------
c
      idow = modulo (inday, 7) + 1
c
c-----------------------------------------------------------------------
c     convert to days after january 1, 1601 using a virtual calendar
c     consistent with current leap-year systems to extend backwards.
c     (ie. 1700, 1800, 1900 are not considered leap-years but 2000
c     is a leap-year)
c-----------------------------------------------------------------------
c
      nday  = inday + 109206
      ibasyr = 1601
c
c-----------------------------------------------------------------------
c     make corrections for years before 1601 (assuming modern leap year
c     conventions)  note that 400 years doesn`t affect the day of the
c     week.   146097 days = 20871 weeks exactly.
c-----------------------------------------------------------------------
c
      if (nday .lt. 0) then
        nfh     = 1 + (-nday - 1)/146097
        ibasyr  = ibasyr - 400*nfh
        nday    = nday + 146097*nfh
      endif
c
c-----------------------------------------------------------------------
c     peel off year.  Acurate from before Jan 1, 1900 until the
c     julian day overflows machine limitations.
c----------------------------------------------------------------------
c
      nfh   = nday/146097
      nday  = modulo (nday, 146097)
      nhund = nday/36524
      if (nhund .gt. 3) then
        nhund = 3
        nday  = 36524
      else
        nday  = modulo (nday, 36524)
      endif
      nfour = nday/1461
      nday  = modulo (nday, 1461)
      nex   = nday/365
      if (nex .gt. 3) then
        nex   = 3
        nday  = 365
      else
        nday  = modulo (nday, 365)
      endif
      leap  = (nex .eq. 3) .and. ((nfour .ne. 24) .or. (nhund .eq. 3))
      if (leap) then
        ndiy = 366
      else
        ndiy = 365
      endif
      iyear = ibasyr + 400*nfh + 100*nhund + 4*nfour + nex
      idoy  = nday + 1
c
c-----------------------------------------------------------------------
c     peel off month and day.
c-----------------------------------------------------------------------
c
      iday = idoy
      do m=1,12
        if (leap .and. (m .eq. 2)) then
          if (iday .le. (daypmi(2)+1)) then
            imonth = 2
            ndim = daypmi(2) + 1
            goto 200
          endif
          iday = iday - daypmi(2) - 1
        else
          if (iday .le. daypmi(m))  then
            imonth = m
            ndim = daypmi(m)
            goto 200
          endif
          iday = iday - daypmi(m)
        endif
      enddo
c
 200  continue
c
      return
      end


      subroutine d2ymdc (inday,
     &                   iyear, imonth, iday, idoy, idow, ndim, ndiy)
c======================================================================
c     inverse of "ymd2dc"
c     d2ymdc takes the number of days since Dec 31, 1899 (or the last
c     day of December 1899 in case equal months of length < 31 days are
c     used) and converts to year, month, day, and other defining
c     quantities.
c     eg: inday=1   yields  Jan 1, 1900.
c
c     equal length years are assumed (no leap years here)
c
c     input:
c       inday  = days since Dec. 31, 1899
c
c     output:
c       iyear  = output year
c       imonth = output month
c       iday   = output day of the month
c       idoy   = day of the current year
c       idow   = day of the week   1=sun - 7=sat
c       ndim   = number of days in the current month
c       ndiy   = number of days in the current year
c
c     author: d. j. goldberg          chg@gfdl.gov
c======================================================================
#include "calendar.h"
c
c-----------------------------------------------------------------------
c     calculate day of the week (sunday=1)
c     inday = 0, i.e., Dec 31, 1899 was a Sunday.
c-----------------------------------------------------------------------
c
      idow = modulo (inday, 7) + 1
c
c-----------------------------------------------------------------------
c     peel off the year
c-----------------------------------------------------------------------
c
      if (inday .le. 0) then
        iyr = (inday / yrlen) - 1
      else
        iyr = (inday - 1) / yrlen
      end if
      iyear = iyr + 1900
      idoy  = inday - yrlen*iyr
c
c-----------------------------------------------------------------------
c     peel off the month and day.
c-----------------------------------------------------------------------
c
      iday = idoy
      do m=1,12
        if (iday .le. daypm(m))  then
          imonth = m
          ndim = daypm(m)
          goto 200
        endif
        iday = iday - daypm(m)
      enddo
c
 200  continue
      ndiy = yrlen
c
      return
      end


      subroutine ymd2d (iyear, imonth, iday,
     &                  nday, idoy, idow, ndim, ndiy)
c
c=======================================================================
c     inverse of "d2ymd"
c     ymd2d takes a date in the year, month, day form and converts it
c     to days since Dec 31. 1899 (nday).  It also needs the
c     cumulative sums of days from the beginning of the year to each
c     month "msum". It returns the number of days in the
c     specified year (ndiy).
c
c     equal year calculations are passed on to subroutine ymd2dc
c
c     input:
c       iyear  = input year
c       imonth = input month
c       iday   = input day of the month
c
c     output:
c       nday   = days since Dec. 31, 1899
c       idoy   = day of the current year
c       idow   = day of the week   1=sun - 7=sat
c       ndim   = number of days in the current month
c       ndiy   = number of days in the current year
c
c     author: d. j. goldberg          chg@gfdl.gov
c======================================================================
c
#include "calendar.h"
c
      logical leap
c
c-----------------------------------------------------------------------
c     have subroutine ymd2dc do equal-year calculations
c-----------------------------------------------------------------------
c
      if (eqyear) then
        call ymd2dc (iyear, imonth, iday, nday, idoy, idow, ndim, ndiy)
        return
      end if
c
      if (mod (iyear, 400) .eq. 0) then
        leap = .true.
      elseif ((mod(iyear,4) .eq. 0) .and. (mod(iyear,100) .ne. 0)) then
        leap = .true.
      else
        leap = .false.
      endif
      if (leap) then
        ndiy = 366
      else
        ndiy = 365
      endif
      nday = iday + msum(imonth)
      if (leap .and. (imonth .gt. 2)) then
        nday = nday + 1
      endif
      idoy = nday
      iyr = iyear - 1601
c
c-----------------------------------------------------------------------
c     make corrections for years before 1601.
c-----------------------------------------------------------------------
c
      if (iyr .lt. 0) then
        nfh  = 1 + (-iyr)/400
        nday = nday - 146097*nfh
        iyr  = iyr + 400*nfh
      endif
      nday = nday + 365*iyr
c
c-----------------------------------------------------------------------
c     correct for leap-years between 1601 and iyear-1.
c     A virtual calendar consistent with current leap-year systems
c     is used to extend backwards.  (ie. 1700, 1800, 1900 are not
c     considered leap-year but 2000 is a leap-year)
c-----------------------------------------------------------------------
c
      iyr = iyr/4
      nday = nday + iyr
      iyr = iyr/25
      nday = nday - iyr
      iyr = iyr/4
      nday = nday + iyr
c
c-----------------------------------------------------------------------
c     nday is now in days since Dec. 31, 1600.  Convert to days since
c     Dec 31, 1899.
c-----------------------------------------------------------------------
c
      nday = nday - 109207
      idow = modulo (nday, 7) + 1
      ndim = daypm(imonth)
      if (leap .and. (imonth .eq. 2)) then
        ndim = ndim + 1
      end if
c
      return
      end


      subroutine ymd2dc (iyear, imonth, iday,
     &                   nday, idoy, idow, ndim, ndiy)
c
c=======================================================================
c     inverse of "d2ymdc"
c     ymd2dc takes a date in the year, month, day form and converts
c     to days since Dec 31. 1899 (nday).  It also returns the number
c     of days in the specified year (ndiy).
c
c     ymd2dc uses a year of constant length with no leap years
c
c     input:
c       iyear  = input year
c       imonth = input month
c       iday   = input day of the month
c
c     output:
c       nday   = days since Dec. 31, 1899
c       idoy   = day of the current year
c       idow   = day of the week   1=sun - 7=sat
c       ndim   = number of days in the current month
c       ndiy   = number of days in the current year
c
c     author: d. j. goldberg          chg@gfdl.gov
c=======================================================================
c
#include "calendar.h"
c
      idoy = iday + msum(imonth)
      iyr = iyear - 1900
      nday = idoy + yrlen*iyr
      idow = modulo (nday, 7) + 1
      ndim = daypm(imonth)
      ndiy = yrlen
c
      return
      end


      function alarm (ialarm, ihalf, timeper, irefer)
c
c=======================================================================
c     inputs:
c
c     ialarm  = index to the absolute alarm time
c
c     ihalf   = index to the time one half time step ahead of current
c               time
c
c     timeper = specified time period of interest (in days. eg: 1.0
c               day, 36 hours (1.5 days), a week (7.0 days), two weeks
c               (14.0 days), 365.0 days ... etc)
c
c     irefer  = index to the time to which the alarm is referenced
c
c     outputs:
c
c     alarm:
c       if timeper > 0 then
c         alarm = true when the current time + dt/2 is later than the
c                 alarm time.  The alarm is then incremented
c       if timeper = 0 then
c         alarm = true (i.e., do this option every step)
c       if timeper < 0 then
c         alarm = false (i.e., this option is disabled)
c
c     if first=.T. then alarm sets the initial alarm time
c
c     authors: d.j.goldberg
c              c.h.goldberg          chg@gfdl.gov
c
c=======================================================================
c
      logical alarm, timeless
c
#include "tmngr.h"
#include "switch.h"
c
      logical ltemp
c
      if (first) then
c
c-----------------------------------------------------------------------
c       initialize the alarm on first iteration
c-----------------------------------------------------------------------
c
        call gettime (ialarm)
        if (timeper .gt. 0.0) then
          if (timeless (irunstart, irefer)) then
c
c-----------------------------------------------------------------------
c           if reference time is later than run start time then the
c           initial alarm is set to reference time.  (No alarms will
c           go off until the time reaches the reference time.)
c-----------------------------------------------------------------------
c
            call copytime (irefer, ialarm)
          else
c
c-----------------------------------------------------------------------
c           set alarm to first time of the form:  reftime + i * timeper
c           that is at least dt/2 after runstart.
c-----------------------------------------------------------------------
c
            call subtime (irunstart, irefer, ireftime)
            call addtime (ireftime, idtd2, ireftime)
            realintervals = realdays(ireftime)/timeper
c
c-----------------------------------------------------------------------
c           i = least integer greater than realintervals (ceiling)
c-----------------------------------------------------------------------
c
            i = int(realintervals) + 1
            i = i - int(i - realintervals)
            call settime3 (itimeper, timeper)
            call multime (i, itimeper, ialarm)
            call addtime (ialarm, irefer, ialarm)
          end if
        else 
          call addtime (itime, idt, ialarm)
        end if
      end if 
c
c-----------------------------------------------------------------------
c     check alarm
c-----------------------------------------------------------------------
c
      if (timeper .gt. 0.0) then
        ltemp = timeless (ialarm, ihalf)
        if (ltemp) then
c
c-----------------------------------------------------------------------
c         increment the alarm time
c-----------------------------------------------------------------------
c
          call settime3 (itimeper, timeper)
 100      continue
          call addtime (ialarm, itimeper, ialarm)
          if (timeless (ialarm, ihalf)) goto 100
        endif
      elseif (timeper .lt. 0.0) then
        ltemp = .false.
        call addtime (itime, idt, ialarm)
      else
        ltemp = .true.
        call addtime (itime, idt, ialarm)
      end if
      alarm = ltemp
c
      return
      end


      subroutine mkstmp (stamp, year, month, day, hour, min, sec)
c
c=======================================================================
c     make a 32 character time stamp from integer inputs
c
c     inputs: integer year, month, day, hour, min, sec
c     output: stamp
c     author:  d. j. goldberg          chg@gfdl.gov
c=======================================================================
c
      integer year, month, day, hour, min, sec
      character*32 stamp
c
      write (stamp,'(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)')
     &             'm/d/y=',month,'/',day,'/'
     &,             year,', h:m:s=', hour,':'
     &,             min,':', sec
      return
      end

      subroutine rdstmp (stmp, year, month, day, hour, min, sec)
c
c=======================================================================
c     convert 32 character time stamp into integer year month day etc
c
c     author:  d. j. goldberg          chg@gfdl.gov
c=======================================================================
c
      character *32 stmp
      character *6  junk6
      character *8  junk8
      character *1  junk1
      integer       year, month, day, hour, min, sec
c
      read  (stmp, '(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)')
     &       junk6, month, junk1, day, junk1, year, junk8, hour,
     &       junk1, min, junk1, sec
c
      return
      end



      subroutine incstamp (stamp1, days, stamp2)
c
c=======================================================================
c     increment a time stamp by a real number of days
c     input:
c      stamp1 =  character date and time stamp
c      days   =  days to increment stamp (may be negative)
c     output:
c      stamp2 = stamp1 incremented by days
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
      character*(*) stamp1, stamp2
      call settimefromstamp (stamp1, itemptime)
      call settime2 (itemptime2, id(days*86400.0), msec(days*86400.0))
      call addtime (itemptime, itemptime2, itemptime)
      call expandtime2 (itemptime)
      stamp2 = tstamp(itemptime)
      return
      end

      subroutine inctime (index1, days, index2)
c
c=======================================================================
c     increment a full time by a real number of days
c     input:
c      index1 =  subscript of a full time 
c      days   =  days to increment full time (may be negative)
c      index2 = subscript of the answer
c     output:
c      time arrays for index2 (year, month, day, etc) in tmngr.h 
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
      call settime2 (itemptime, id(days*86400.0), msec(days*86400.0))
      call addtime (index1, itemptime, index2)
      call expandtime2 (index2)
      return
      end

      subroutine settimefromstamp (stamp, index)
c
c=======================================================================
c     input:
c      stamp =  character dat and time stamp 
c      index = subscript of the answer
c     output:
c      time arrays for index (year, month, day, etc) in tmngr.h 
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
      character*(*) stamp
      call rdstmp (stamp, iyear, imonth, iday, ihour, imin, isec)
      call setfulltime (index, iyear, imonth, iday, ihour, imin, isec)
      return
      end



      subroutine addtime (index1, index2, index)
c
c=======================================================================
c     add two times given in (integer day, nonneg integer ms) form
c     input:
c      index1 =  subscript of the first time into time arrays in tmngr.h 
c      index2 =  subscript of the second time
c      index  = subscript of the answer
c     output:
c      iday(index) = integer day number of answer in tmngr.h
c      msday(index)= millisec fractional day of answer in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
c
      mst = msday(index1) + msday(index2)
      it  = iday(index1) + iday(index2)
      if (mst .ge. 86400000) then
        mst = mst - 86400000
        it  = it + 1
      end if
      iday(index) = it
      msday(index) = mst
c
      return
      end


      subroutine subtime (index1, index2, index)
c
c=======================================================================
c     subtract two times given in (integer day, nonneg integer ms) form
c     input:
c      index1 =  subscript of the first time into time arrays in tmngr.h 
c      index2 =  subscript of the second time
c      index  = subscript of the answer
c     output:
c      iday(index) = integer day number of answer in tmngr.h
c      msday(index)= millisec fractional day of answer in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
c
      mst = msday(index1) - msday(index2)
      it  = iday(index1) - iday(index2)
      if (mst .lt. 0) then
        mst = mst + 86400000
        it = it - 1
      end if
      iday(index)  = it
      msday(index) = mst
c
      return
      end


      subroutine multime (n, index1, index)
c
c=======================================================================
c     multiply time (integer day, nonneg integer ms) by an integer
c     input:
c      n      = integer multiple    
c      index1 = subscript of the time
c      index  = subscript of the answer
c     output:
c      iday(index) = integer day number of answer in tmngr.h
c      msday(index)= millisec fractional day of answer in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
c
      double precision dmst
c
      dmst = n*dble(msday(index1))
      it   = n*iday(index1)
      it   = it + int(dmst/86400000.0d0)
      mst  = nint(mod (dmst, 86400000.0d0))
      iday(index)  = it
      msday(index) = mst
c
      return
      end


      subroutine ms2hms (ms, hour, min, sec)
c
c=======================================================================
c     converts fraction of a day in milliseconds to hour/min/sec
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
      integer ms, hour, min, sec
c
      hour  = ms / 3600000
      msrem = ms - hour * 3600000
      min   = msrem / 60000
      sec   = (msrem - min * 60000) / 1000
c
      return
      end


      function hms2ms (hour, min, sec)
c
c=======================================================================
c     converts fraction of a day in hour/min/sec to milliseconds
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
      integer hms2ms, hour, min, sec
c
      hms2ms = hour*3600000 + min*60000 + sec*1000
c
      return
      end


      function id (realsec)
c
c=======================================================================
c     converts time in real seconds to the integer day part
c     see also "msec" 
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
      id = int(realsec/86400.0)
c
      return
      end


      function msec (realsec)
c
c=======================================================================
c     extracts integer milliseconds from time in real seconds
c     see also "id" 
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
      msec = nint (1000.0*(realsec-86400.0*id(realsec)))
c
      return
      end


      function realsecs (index)
c
c=======================================================================
c     converts integer days and milliseconds to real seconds
c     input:
c      index = index of any time in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
c
      realsecs = 86400.0*iday(index) + msday(index)/1000.0
c
      return
      end


      function realdays (index)
c
c=======================================================================
c     converts integer days and milliseconds to real days
c     input:
c      index = index of any time in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
c
      realdays = iday(index) + msday(index)/86400000.0
c
      return
      end


      function modulo (m, n)
c
c=======================================================================
c     Fortran 90 intrinsic
c     similar to mod, but remainder has sign of n.
c     always gives positive remainder when n .gt. 0, even if m .lt. 0
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
      if (m .ge. 0) then
        modulo = mod(m,n)
      else
        modulo = abs(n) - mod(-m-1,n) - 1
      end if
c
      return
      end


      subroutine inittime
c
c=======================================================================
c     initialize time indices for gettime and getfulltime.
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
c
      nextfulltime = 1
      nexttime = nfulltimes + 1
c
      return
      end


      subroutine getfulltime (index)
c
c=======================================================================
c     allocates and returns an index for a full time and increments the
c     next available full time index counter.
c     
c     output:
c      index = index of a full time in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      if (nextfulltime .gt. nfulltimes) then
        write (stdout, "(a,a,i4,a)") '     Too many full times.',
     &                    '  Increase nfulltimes = ', nfulltimes,
     &                    ' in tmngr.h'
        stop 'getfulltime'
      endif
      index = nextfulltime
      nextfulltime = nextfulltime + 1
c
      return
      end


      subroutine gettime (index)
c
c=======================================================================
c     returns an index for a time and increments the next available
c     time index counter.
c     output:
c      index = index of a short time (iday, msday only) in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      if (nexttime .gt. ntimes) then
        write (stdout, "(a,a,i4,a)") '     Too many times.',
     &                    '  Increase ntimes = ', ntimes,
     &                    ' in tmngr.h'
        stop 'gettime'
      endif
      index = nexttime
      nexttime = nexttime + 1
c
      return
      end


      subroutine copyfulltime (index1, index)
c
c=======================================================================
c     copy all fields of a full time structure from one index to
c     another.
c     input:
c      index1 = index of a full time in tmngr.h
c      index  = index of a full time in tmngr.h
c     output:
c      full time fields set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
c
      if (badindex (index1, 'full')) stop 'copyfulltime'
      if (badindex (index,  'full')) stop 'copyfulltime'
      iday  (index) = iday  (index1)
      msday (index) = msday (index1)
      year  (index) = year  (index1)
      month (index) = month (index1)
      day   (index) = day   (index1)
      hour  (index) = hour  (index1)
      minute(index) = minute(index1)
      second(index) = second(index1)
      tstamp (index) = tstamp (index1)
      dayofyear (index) = dayofyear (index1)
      dayofweek (index) = dayofweek (index1)
      daysinmon (index) = daysinmon (index1)
      daysinyear(index) = daysinyear(index1)
c
      return
      end


      subroutine copytime (index1, index)
c
c=======================================================================
c     copy both fields of a short time structure from one index to
c     another.
c     input:
c      index1 = index of a time in tmngr.h
c      index  = index of a short time in tmngr.h
c     output:
c      short time fields set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
c
      if (badindex (index1, 'short')) stop 'copytime'
      if (badindex (index,  'short')) stop 'copytime'
      iday (index) = iday (index1)
      msday(index) = msday(index1)
c
      return
      end


      subroutine settime2 (index, it, mst)
c
c=======================================================================
c     set a time by placing given day and millisecond information in
c     a time structure array indexed by index.
c     input:
c      index  = index of a short time in tmngr.h
c      it     = integer day part
c      mst    = fractional part of day in millisec 
c     output:
c      short time fields set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
c
      if (mst .ge. 86400000 .or. mst .lt. 0) then
        write (stdout, *) '     Invalid millisecond designation'
        stop 'settime2'
      end if
      if (badindex (index,  'short')) stop 'settime2'
      iday(index)  = it
      msday(index) = mst
c
      return
      end


      subroutine settime3 (index, rd)
c
c=======================================================================
c     set a time structure given a time in real days
c     input:
c      index  = index of a short time in tmngr.h
c      rd     = time in real days
c     output:
c      short time fields set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
c
      if (badindex (index,  'short')) stop 'settime3'
c
c-----------------------------------------------------------------------
c     greatest integer less than     (floor)
c-----------------------------------------------------------------------
c
      it  = int(rd) - 1
      it  = it + int(rd - it)
      mst = nint(86400000.0*(rd-it))
      iday(index)  = it
      msday(index) = mst
c
      return
      end


      subroutine setfulltime2 (index, it, mst)
c
c=======================================================================
c     set a full time by placing given day and millisecond information
c     in a time structure array indexed by index and expanding it
c     using the d2ymd conversion and making a stamp.
c     input:
c      index  = index of a full time in tmngr.h
c      it     = integer day part
c      mst    = fractional part of day in millisec 
c     output:
c      full time fields set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
c
      if (badindex (index, 'full')) stop 'setfulltime2'
      call settime2 (index, it, mst)
      call d2ymd (it, year(index), month(index), day(index),
     &            dayofyear(index), dayofweek(index), daysinmon(index),
     &            daysinyear(index))
      call ms2hms(mst, hour(index), minute(index), second(index))
      call mkstmp (tstamp(index), year(index), month(index), day(index),
     &             hour(index), minute(index), second(index))
c
      return
      end


      subroutine setfulltime (index, yr, mon, dy, hr, min, sec)
c
c=======================================================================
c     set a full time by specifying the y/m/d, h/m/s and expanding it
c     using the ymd2d conversion and making a stamp.
c     input:
c      index  = index of a full time in tmngr.h
c      yr     = year to set
c      mon    = month to set
c      dy     = day to set
c      hr     = hour to set
c      min    = minute to set
c      sec    = second to set
c     output:
c      full time fields set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
      integer yr, mon, dy, hr, min, sec, hms2ms
c
      if (badindex (index, 'full')) stop 'setfulltime'
      year  (index) = yr
      month (index) = mon
      day   (index) = dy
      hour  (index) = hr
      minute(index) = min
      second(index) = sec
      call ymd2d (yr, mon, dy, iday(index),
     &            dayofyear(index), dayofweek(index),
     &            daysinmon(index), daysinyear(index))     
      msday(index) = hms2ms(hr, min, sec)
      call mkstmp (tstamp(index), yr, mon, dy, hr, min, sec)
c
      return
      end


      subroutine expandtime (index)
c
c=======================================================================
c     expand a full time specified by the y/m/d, h/m/s
c     using the ymd2d conversion and making a stamp.
c     input:
c      index  = index of a full time in tmngr.h
c      y/m/d, h/m/s must be set for time index
c     output:
c      full time fields (iday,msday, tstamp) set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
      integer hms2ms
c
      if (badindex (index, 'full')) stop 'expandtime'
      call ymd2d (year(index), month(index), day(index),
     &            iday(index), dayofyear(index), dayofweek(index),
     &            daysinmon(index), daysinyear(index))     
      msday(index) = hms2ms(hour(index), minute(index), second(index))
      call mkstmp (tstamp(index), year(index), month(index), day(index),
     &             hour(index), minute(index), second(index))
c
      return
      end


      subroutine expandtime2 (index)
c
c=======================================================================
c     expand a full time specified by absolute day and millisecond
c     information using d2ymd conversion and making a stamp.
c     input:
c      index  = index of a full time in tmngr.h
c      iday, msday must be set for time index
c     output:
c      full time fields (year,month,day,hour,minute,second,tstamp)
c      set for index in tmngr.h
c
c     author:  c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
c
      if (badindex (index, 'full')) stop 'expandtime2'
      call d2ymd (iday(index), year(index), month(index),
     &            day(index), dayofyear(index), dayofweek(index),
     &            daysinmon(index), daysinyear(index))
      call ms2hms(msday(index), hour(index), minute(index),
     &            second(index))
      call mkstmp (tstamp(index), year(index), month(index), day(index),
     &             hour(index), minute(index), second(index))
c
      return
      end


      function timeless (index1, index2)
c
c=======================================================================
c     compare times referenced by index1 and index2.  timeless is
c     true if time1 is less than time2.
c
c     author:  d. j. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical timeless, badindex
c
      if (badindex (index1, 'short')) stop 'timeless'
      if (badindex (index2, 'short')) stop 'timeless'
      if (iday(index1) .lt. iday(index2)) then
        timeless = .true.
      elseif ((iday (index1) .eq. iday (index2)) .and.
     &        (msday(index1) .lt. msday(index2))) then
        timeless = .true.
      else
        timeless = .false.
      end if
c
      return
      end


      function timeequal (index1, index2)
c
c=======================================================================
c     compare times referenced by index1 and index2.  timeequal is
c     true if time1 is equal to  time2.
c
c     author:  d. j. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical timeequal, badindex
c
      if (badindex (index1, 'short')) stop 'timeequal'
      if (badindex (index2, 'short')) stop 'timeequal'
      if (iday(index1) .eq. iday(index2) .and. 
     &    msday(index1) .eq. msday(index2)) then
        timeequal = .true.
      else
        timeequal = .false.
      end if
c
      return
      end


      function timemore (index1, index2)
c
c=======================================================================
c     compare times referenced by index1 and index2.  timemore is
c     true if time1 is more than time2.
c
c     author:  d. j. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical timemore, badindex
c
      if (badindex (index1, 'short')) stop 'timemore'
      if (badindex (index2, 'short')) stop 'timemore'
      if (iday(index1) .gt. iday(index2)) then
        timemore = .true.
      elseif ((iday (index1) .eq. iday (index2)) .and.
     &        (msday(index1) .gt. msday(index2))) then
        timemore = .true.
      else
        timemore = .false.
      end if
c
      return
      end


      function badindex (index, timetype)

c
c=======================================================================
c     check to see if a given index is in bounds in tmngr.h
c     type = 'full'  means a full time
c     type = 'short' means a short time
c
c     author:  d. j. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "tmngr.h"
c
      logical badindex
      character *(*) timetype
c
      if (timetype .eq. 'full') then
        if (index .lt. 1 .or. index .gt. nfulltimes) then
          badindex = .true.
          print *, '     Invalid full time index = ',index
        else
          badindex = .false.
        end if
      else
        if (index .lt. 1 .or. index .gt. ntimes) then
          badindex = .true.
          print *, '     Invalid time index = ',index
        else
          badindex = .false.
        end if
      end if
c
      return
      end


