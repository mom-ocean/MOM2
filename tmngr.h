c====================== include file "tmngr.h" =========================
c
c                       time manager variables
c
c-----------------------------------------------------------------------
c     time manager inputs:
c-----------------------------------------------------------------------
c
c     how to choose a reference time:
c
c     refrun  = (true,false) to base calculation for diagnostic switches
c              on (the start of each job, other reference time)
c              example:
c              suppose each job submission integrates
c              for one month but the number of days per month changes.
c              setting "refrun" = true and setting
c              "timavgint" = (days in month)/3 will give 3 averaging
c              periods per month of approximately 10 days each. the
c              only restriction is that "timavgint"is an integral number
c              of time steps (if not then "timavgint" is reset to insure
c              this condition. other diagonstic switches do not have
c              this restriction).
c
c     refinit = (true, false) for basing calculation of logical switches
c              on (initial conditions, other reference time)
c              example: if term balances are desired every 20 days
c              (trmbint=20.0) and refinit = true, then they
c              will be done every 20 days starting from initial
c              condition time.
c
c     refuser = (true, false) to base calculations of logical switches
c              on (user-chosen reference time, other reference time)
c              if refuser = true, the user must also supply values for
c              ryear, rmonth, rday, rhour, rmin, rsec (integer)
c              example: if term balances are desired every 20 days
c              (trmbint=20.0) and refuser = true, then they will be done
c              every 20 days counting from reference time, ignoring the
c              initial condition time. for comparing diagnostics from
c              various experiments with different initial condition
c              times, refuser = true will be more appropriate. setting
c              refuser = true and choosing the reference time to be
c              the initial condition time is the same as refinit = true.
c
c     summary of how to choose the time for referencing calculations
c     of logical switches
c
c     refrun  = T ==>  referenced to the start of each run
c     refinit = T ==>  referenced to initial condition time given by:
c                     year0, month0, day0, hour0, min0, sec0
c     refuser = T ==>  referenced to user specified reference time so
c                     must set: ryear, rmonth, rday, rhour, rmin, rsec
c
c-----------------------------------------------------------------------
c
c     time variable arrays
c
c     arrays "iday" and "msday" contain the primary internal 
c     representation of all times within the time manager. they are
c     referenced by using a subscript to indicate which time.
c
c     iday    = integer days (since Dec 31, 1899 when specifying a date)
c     msday   = non-negative integer milliseconds
c
c     it is desirable to have time information expanded to include the
c     following secondary time fields:
c
c     year       = 
c     month      = 
c     day        = 
c     hour       = 
c     minute     = 
c     second     =
c     tstamp     = 32 character date and time stamp m/d/y h:m:s
c     dayofyear  = integer day of the year (1..yrlen)
c     dayofweek  = 1=sun - 7=sat
c     daysinmon  = days in the month
c     daysinyear = days in the year
c
c     those times for which primary and secondary information is 
c     maintained by the time manager are called "full times". those for
c     which only primary information is kept are called "short times"
c
c     indices to  "full times" (including year, month ,day, etc).
c
c     itime     = simulation time corresponding to "itt"
c     initial   = time of the initial conditions
c     irunstart = time of the start of the run
c     iuser     = user defined reference time
c     iref      = one of the three above selected by logicals
c                 (refinit, refrun, refuser)
c
c     indices to  "short times". ("iday", "msday" only)
c
c     isunday    = time of a sunday for week and two week switches
c     ihalfstep  = dt/2 beyond itime
c     imodeltime = time since initial conditions
c     iruntime   = time since run start
c     iusertime  = time since user specified reference time
c     idt        = integer days and milliseconds of dt
c     idtd2      = integer days and milliseconds of dt/2
c
c     ireftime   = time used locally in alarm function
c     itimeper   = locally use for decomposition of real timeper
c
c     for any time index (short or full) the internal representation
c     may be converted to either real days or real seconds using
c     the functions:
c                  realdays(index)
c                  realsecs(index)
c
c     dayoyr  = relative day number referenced to the beginning
c               of the current year.  (real)
c     relyr   = number of years (and fractional years) of model
c               integration (for time tau+1 {itt}) relative to 
c               initial condition
c     prelyr  = relyr for previous time step
c
c     stamp   = 32 character date and time for current model timestep
c     pstamp  = 32 character date and time for previous model timestep
c
c     itt     = current time step counter (from initial cond.)
c     itt0    = time step at start of current run
c
c               variables used for initialization
c
c     irstdy  = integer number of days at start of run
c     msrsdy  = fractional day in millisec at start of run
c
c     year0   = year of initial conditions
c     month0  = month of initial conditions
c     day0    = day of initial conditions
c     hour0   = hour of initial conditions
c     min0    = minute of initial conditions
c     sec0    = second of initial conditions  
c     
c     ryear   = year of user specified reference time
c     rmonth  = month of user specified reference time
c     rday    = day of user specified reference time
c     rhour   = hour of user specified reference time
c     rmin    = minute of user specified reference time
c     rsec    = second of user specified reference time  
c
c-----------------------------------------------------------------------
c
c
      logical refrun, refinit, refuser
c
      parameter (ntimes = 100, nfulltimes = 20)
      integer iday(ntimes), msday(ntimes)
      integer year(nfulltimes), month(nfulltimes), day(nfulltimes)
      integer hour(nfulltimes), minute(nfulltimes), second(nfulltimes)
      integer dayofyear(nfulltimes), dayofweek(nfulltimes)
      integer daysinmon(nfulltimes), daysinyear(nfulltimes)
      integer year0, month0, day0, hour0, min0, sec0
      integer ryear, rmonth, rday, rhour, rmin, rsec
      character *32 tstamp(nfulltimes), pstamp, stamp
c
      common /tmngrl/ refrun, refinit, refuser
c
      common /tmngrii/ nextfulltime, nexttime
      common /tmngrii/ initial, iref, irunstart, itime, iuser
      common /tmngrii/ iruntime, imodeltime, ireftime, iusertime
      common /tmngrii/ ihalfstep, itimeper, isunday
      common /tmngrii/ itemptime, itemptime2
      common /tmngrii/ idt, idtd2
      common /tmngrii/ iday, msday
      common /tmngrii/ year, month, day, hour, minute, second
      common /tmngrii/ dayofyear, dayofweek, daysinmon, daysinyear
      common /tmngrii/ itt0, itt, irstdy, msrsdy
c
      common /tmngrr/ dayoyr
      common /tmngrr/ relyr, prelyr
c
      common /tmngrc/ tstamp, stamp, pstamp
c
      common /tmngrii/ year0, month0, day0, hour0, min0, sec0
      common /tmngrii/ ryear, rmonth, rday, rhour, rmin, rsec
c

