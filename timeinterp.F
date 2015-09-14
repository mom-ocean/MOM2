#ifdef test_timeinterp
# include "tmngr.F"
# include "util.F"
c
      program driver
c
c=======================================================================
c     To test the time interploator module, use the "run_timeinterp"
c     run script
c
c     author: R. C. Pacanowski        e-mail=> rcp@gfdl.gov              
c=======================================================================
c
#include "stdunits.h"
#include "calendar.h"
#include "tmngr.h"
#include "switch.h"
      logical period, readit, leapyrdata
      character*12 recname
      parameter (nrec=12)
      character*32 dstamp(nrec)
      dimension tdrec(nrec), aprec(nrec), recname(nrec)
c
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
      nmix      = 17
c
      write (stdout,'(///,tr5,a,g14.7,a)')
     &'Driver will test the time interpolator by simulating integrating'
     &, runlen, rununits
      print *, '     using timesteps of dt = ',dt,' (sec) '
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
c     Setup to test interpolating data from climatological monthly 
c     datasets to current model time (at time step "itt")
c
c     choose either leap year corrected interpolation data using
c     28.2425 days in February, or non leap year correct data using
c     28 days
c
c     define the dataset as periodic and choose an interpolation
c     method.  see subroutine "timeinterp" included within this module
c     to explain the various methods of interpolation.
c-----------------------------------------------------------------------
c
      leapyrdata = .true.
      period     = .true.
      method     = 3
c
c     Set the starting time for the first record.
c     Build the time stamps and periods for each record. Note that the
c     time stamps mark the end of each record. (1st one is for january)
c
      isbcyear  = 0
      isbcmonth = 1
      isbcday   = 1
      isbchour  = 0
      isbcmin   = 0
      isbcsec   = 0
      isbcmon   = mod(isbcmonth,nrec) + 1
c      
      call mkstmp (dstamp(1), isbcyear, isbcmon, isbcday, isbchour
     &,            isbcmin, isbcsec)
      print *,' End of 1st record =',dstamp(1), ' month = ',isbcmonth
c
c     use "daypm" array from calendar.h to build record periods and
c     date stamps for end of records
c
      do m=1,nrec
        mm         = mod(m+isbcmonth+nrec-2,nrec) + 1
        if (leapyrdata .and. mm .eq. 2) then
          aprec(m) = daypm(mm) + 0.2425
        else
	  aprec(m)   = daypm(mm)
        endif
	recname(m) = monname(mm)
	if (m .gt. 1) call incstamp (dstamp(m-1), aprec(m-1), dstamp(m))
      enddo
c
c     initialize time interpolation arrays
c
      call timeinterpi (nrec, dstamp(1), aprec, tdrec, isbcstart
     &,                 period)
c
      do m=1,nrec
	print '(a,i3,3a,f10.4,a,f10.4)',' data record # ',m
     &, ' is for ',recname(m), ' tdrec=',tdrec(m), ' aprec=',aprec(m)
      enddo
c
c     express model time in days after start of S.B.C. by adding time
c     of I.C. to current model time then subtract time at start of
c     S.B.C.. Note that "itemptime" was allocated in tmngri and is
c     only needed as a temporary.
c
      call addtime (initial, imodeltime, itemptime)
      call subtime (itemptime, isbcstart, itemptime)
      daysbc = realdays(itemptime)
c
c     initialize pointers to memory buffers
c      iprevm will point to memory holding data before model timestep
c      inextm will point to memory holding data after model timestep     
c     also, set the number of datasets being interpolated to 1
c     
      iprevm = 1
      inextm = 2
      ntdrec = 1
c
      call timeinterp (daysbc, ntdrec, tdrec, aprec, nrec, period
     &,                method, inextd, iprevd, wprev, readit, inextm
     &,                iprevm)
c
c     read data from disk records "inextd" and "iprevd" into
c     memory buffers "inextm" and "iprevm" to initialize them
c
      write (*,'(/a,a,a,i2/a,a,a,i2/)')  
     &  ' read disk data record for ', recname(iprevd)
     &, ' into memory buffer ',iprevm 
     &, ' and disk data record for ', recname(inextd)
     &, ' into memory buffer ', inextm
          print '(a,i3,a,i3,a,f10.8,a,i2,a,i2)'
     &,   ' inextd=',inextd,' iprevd=',iprevd, ' wprev=',wprev
     &,   ' inextm=',inextm,' iprevm=',iprevm
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
c       test for calculating when to read monthly S.B.C. (surface
c       boundary conditions) for use in interpolating them to the   
c       current model time.
c
c       express model time in days after start of S.B.C. by adding time
c       of I.C. to model time since I.C. then subtract time at start of
c       S.B.C.. Note that "itemptime" was allocated in tmngri and is
c       only needed as a temporary.
c-----------------------------------------------------------------------
c
        call addtime (initial, imodeltime, itemptime)
        call subtime (itemptime, isbcstart, itemptime)
        daysbc = realdays(itemptime)
c
        call timeinterp (daysbc, ntdrec, tdrec, aprec, nrec, period
     &,                  method, inextd, iprevd, wprev, readit, inextm
     &,                  iprevm)
c
c-----------------------------------------------------------------------
c       print when things happen
c-----------------------------------------------------------------------
c
        if (itt .lt. 10) then
          write (*,9000) itt, stamp, daymod, daysbc
        endif
c
c         show when to update the memory buffers holding monthly S.B.C.
c         by reading data from disk record "inextd"  into
c         memory buffers "inextm"
c
          if (readit) then
            write (*,'(/a,i7,a,g14.7,a,a,a,i2,a,i1,a)')  
     &     ' ==> time to read data on itt=',itt, ', day= ',daymod
     &,    ' from ', recname(inextd), ' into buffer ', inextm
     &,   ' (previous data is now in ',iprevm,')'
c            write (*,'(1x,a,i5,a,i2,a,i2,a,e14.7,a,i2,a,i2)')
c     &    ' itt=',itt,' inextd=',inextd,' iprevd=',iprevd
c     &,   ' wprev=',wprev, ' inextm=',inextm,' iprevm=',iprevm
          endif
	  if (wprev .lt. 0.005 .or. wprev .gt. 0.995) then
	    print '(1x,a,i7,a,a,a,f12.8,a,f12.8)'
     &,     ' itt=',itt, ' ',stamp, ' weight for (prev,next) data ='
     &,      wprev,',',1.0-wprev
	  endif	
c
        if (eomon) print 9001,' ...end of month on itt=', itt, ' ',
     &                     stamp, dayname(dayofweek(itime))
        if (midmon) print 9001,' ... mid month on itt  =', itt, ' ',
     &                      stamp, dayname(dayofweek(itime))
        if (eoyear) print 9001,' ==> end of year on itt=',itt,' ',
     &                      stamp, dayname(dayofweek(itime))
c
        if (eorun)  then
	  print 9001,' ==> end of run on itt=',itt,' ', stamp
     &,               dayname(dayofweek(itime))
          stop
	endif
c
      enddo
      print *,' Increase loop limit on "n" in "driver.F"'
      stop 'driver'
9000  format (1x, ' itt=',i7,1x,a,' model days since I.C.=',f12.4
     &,      '. days since start of S.B.C.=', f12.4)
9001  format (a, i8, a, a, 1x, a)
      end
#endif



      subroutine timeinterpi (nrec, stamp1, aprec, tdrec, isbcstart
     &,                       period)
c              
c=======================================================================
c     initializes time center of each data record based on the time
c     stamp and average period.
c
c     author: c. h. goldberg        e-mail=> chg@gfdl.gov
c             r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "tmngr.h"
      dimension aprec(nrec), tdrec(nrec)
      character*(*) stamp1
      logical period
      data isbcend /0/
      save isbcend
c
c     define each climatological data record to be at the center of the
c     month starting with month "isbcmonth"
c
      sum = 0.0
      do m=1,nrec
        sum      = sum + aprec(m)
        tdrec(m) = sum - 0.5*aprec(m)
      enddo
c
c     calculate time at start of first record: "isbcstart"
c
      if (isbcend .eq. 0) call getfulltime (isbcend)
      call getfulltime (isbcstart)
c      
      call rdstmp (stamp1, isbcyear, isbcmon, isbcday, isbchour
     &,            isbcmin, isbcsec)
      call setfulltime (isbcend, isbcyear, isbcmon, isbcday, isbchour
     &,            isbcmin, isbcsec)
      call inctime (isbcend, -aprec(1), isbcstart)
c
c     check integrity of data record times. also when using datasets
c     as periodic, add 0.2425 days to febuary and adjust subsequent
c     months to account for this change when using the time manager 
c     with a leap year calendar.           
c
      if (period) call checkinterp (nrec, tdrec, aprec) 
      return
      end


      subroutine timeinterp (tm, n, tdrec, aprec, ndr, period, method
     &,                     ia, ib, wb, change, inext, iprev)
c
c=======================================================================
c
c     time interpolator ... constructs indices & weight needed for
c     linearly interpolating data defined at arbitrary time intervals
c     (midpoints of years, months, days or  random intervals) to
c     the time of the current model time step.
c
c     inputs:
c
c     tm     = the time at which the data is desired (units of "tdrec")
c
c     tdrec  = the times at which the data records in the dataset are
c              defined. times must be monotonically increasing and are
c              assumed to be at the centers of the averaging periods.
c              (eg: the centers of the months if using monthly averaged
c               climatology. units are arbitrary)
c
c     aprec  = array of averaging periods for the data records
c              (eg: the number of days per month)
c
c     ndr    = number of data records in the dataset. (eg: 12 if using
c              monthly climatology)
c
c     period = (true,false) if the dataset is to be treated as
c              (perodic, not periodic). if periodic, then the model
c               time is always mapped into the dataset. if not, then
c               record 1 is used for all model time before the 
c               beginning of the dataset and record "ndr" is used for
c               all model time after the end of the dataset.
c
c     method = interpolation scheme desired.  (0..3)
c                0 = no interpolation; the average value is used
c                    for all times in the entire averaging period.
c                    (preserves the integral over averaging periods,
c                    but is discontinuous at period boundaries.)
c                1 = linear interpolation between the middles of
c                    two adjacent averaging periods.
c                    (continuous but does not preserve integral for
c                    unequal periods.)
c                2 = equal linear interpolation.  Assumes that the
c                    value on the boundary between two adjacent
c                    averaging periods is the unweighted average of
c                    the two average values.  Linearly interpolates
c                    between the midperiod and period boundary.
c                    (continuous but does not preserve integral for
c                    unequal periods.)
c                3 = equal area (midperiod to midperiod) interpolation
c                    chooses a value for the boundary between two
c                    adjacent periods such that linear interpolation 
c                    between the two midperiods and this value will
c                    preserve the integral midperiod to midperiod.
c                Note that methods 1,2, and 3 are equivalent if
c                all periods lengths are equal.
c
c     n      = a number denoting which dataset is being interpolated
c              (each dataset should be referenced by a unique number
c               starting with 1 for the 1st, 2 for the 2nd, ...etc)
c
c
c     outputs:
c
c     ia     = index for pointing to the next data record which will be
c              reached by the model. (eg: ahead of the model. "ia" would
c              be 3 if "tm" was beyond the  middle of {but still within} 
c              february)
c     ib     = index for pointing to the data record which was just
c              passed by the model. (eg: behind the model. "ib" would
c              be 2 if "tm" was beyond the middle of {but still within}
c              february)
c     inext  = index to memory buffer containing data from "ia"
c     iprev  = index to memory buffer containing data from "ib"
c     wb     = interpolation weight for defining data at "tm"
c              schematically the interpolation is defined by:
c
c              data(iprev) <== disk data "ib"              
c              data(inext) <== disk data "ia"
c              data(tm) = wb*data(iprev) + (1-wb)*data(inext)
c
c     change = logical for sensing when "ia" and "ib" change.
c              when change = T then it is time to read the disk
c              and update "inext" and "iprev"
c     
c     author:        r. c. pacanowski      e-mail=> rcp@gfdl.gov
c     revised by:    d. j. goldberg        e-mail=> chg@gfdl.gov         
c=======================================================================
c
      logical change, period
      parameter (maxsets=10, iflag=-99999)
      dimension iaold(maxsets), tdrec(ndr), aprec(ndr)
      dimension imethod(maxsets)
      data iaold /maxsets*iflag/
      save iaold, imethod
c
c-----------------------------------------------------------------------
c     statement function
c-----------------------------------------------------------------------
c
      frac(r) = r - int(r)
c
      if (n .gt. maxsets) then
        write (*,'(a,i10,a,i10)') 'Error: n=', n, ' maxsets=',maxsets
        stop '=>timeinterp'
      end if
c
      if (iaold(n) .eq. iflag) then
        write (*,'(/1x,a,i2,a,i3/)')
     &      'Assigning interpolation method ',method, ' to dataset # ',n
        imethod(n) = method
      endif 
c
      if (method .ne. imethod(n)) then
        write (*,'(/a,i2,a,i3/a,i2,a/)')
     &   'Error: trying to use method ',method, ' on dataset # ',n
     &,  'originally, method ',imethod(n),' was used in timeinterp'
        stop 
      endif 
c
      if (period) then
c
c       define the position of the dataset in time
c      
        dstart = tdrec(1) - 0.5*aprec(1)
        dend   = tdrec(ndr) + 0.5*aprec(ndr)
        dlen   = dend - dstart
c
c       map the model time into the dataset assuming dataset periodicity
c
        if (tm .lt. dstart) then
          d = dstart - tm
          f = frac(d/dlen)
          time = dend - f*dlen
        elseif (tm .gt. dend) then
          d = tm - dend
          f = frac(d/dlen)
          time = dstart + f*dlen
        else
          time = tm
        end if
      else
c
c       define the position of the dataset in time. no periodicity
c      
        dstart = tdrec(1)
        dend   = tdrec(ndr)
        dlen   = dend - dstart
c
c       map the model time into the dataset. assume data is constant
c       before the beginning and after the end of the dataset
c
        if (tm .lt. dstart) then
          time = dstart
        elseif (tm .gt. dend) then
          time = dend 
        else
          time = tm
        end if
      end if
c
c     calculate record pointers and weighting for interpolation of
c     dataset records to the model time step.
c
      ib = indp (time, tdrec, ndr)
      if (tdrec(ib) .gt. time) ib = ib - 1
      if (period) then
        ia = mod(ib, ndr) + 1
        if (ib .lt. 1) ib = ndr
      else
        ia = ib + 1
        if (ia .gt. ndr) ia = ib
        if (ib .lt. 1)   ib = ia
      end if
c
c     find whether "time" is closer to midpoint of record "ia" or ib"
c     ic is the index of the closest midpoint
c     io is the index of the other midpoint
c
      startaft = tdrec(ia) - 0.5*aprec(ia)
      if (time .ge. startaft .and. time .le. tdrec(ia)) then
        ic = ia
        io = ib
      else
        ic = ib
        io = ia
      end if
c
c     dtmid = distance from "time" to midpoint of closer record
c     dtbnd = distance from "time" to boundary of closer record
c     dtomid = distance from "time" to midpoint of other record
c
      dtmid  = abs(time - tdrec(ic))
      dtbnd  = 0.5*aprec(ic) - dtmid
      dtomid = 0.5*aprec(io) + dtbnd
c
c-----------------------------------------------------------------------
c     3) equal area (midperiod to midperiod) interpolation formula
c-----------------------------------------------------------------------
c
      if (method .eq. 3) then
        wc = 2.0*dtbnd/aprec(ic) + 2.0*dtmid/(aprec(ic) + aprec(io))
c
c-----------------------------------------------------------------------
c     2) equal linear interpolation
c             value on period boundary assumed to be average of values
c             on the two adjacent periods.
c-----------------------------------------------------------------------
c
      elseif (method .eq. 2) then
        wc = (2.0*dtbnd + dtmid)/aprec(ic)
c
c-----------------------------------------------------------------------
c     1) linear interpolation
c-----------------------------------------------------------------------
c
      elseif (method .eq. 1) then
        wc = dtomid/(dtmid + dtomid)
c
c-----------------------------------------------------------------------
c     0) no interpolation
c-----------------------------------------------------------------------
c
      elseif (method .eq. 0) then
        wc = 1.0
      else
c
c-----------------------------------------------------------------------
c     anyting else is not allowed for (unless you want to add one!)
c-----------------------------------------------------------------------
c
        print *,'=>Error: method = ',method,' not allowed in timeinterp'
        stop 
      end if
c
      if (ib .eq. ic) then
        wb = wc
      else
        wb = 1.0 - wc
      end if
      if (wc .lt. 0.0 .or. wc .gt. 1.0) then
        print *,' ic=',ic,' io=',io, ' dtmid=',dtmid,' dtbnd=',dtbnd
     &,' dtomid=',dtomid, ' time=',time, ' ia=',ia,' ib=',ib
     &, ' wc=',wc
        print *,' =>Error: bad interpolation wieght in timeinterp'
        stop
      endif
c
c     refresh pointers to memory buffers when reading disk data
c
      if (iaold(n) .ne. ia) then
        change = .true.
        itemp = iprev
        iprev = inext
        inext = itemp
      else
        change = .false.
      end if
      iaold(n) = ia
c
#ifdef debug_timeinterp
      print *,'timeinterp: ic=',ic,' io=',io, ' dtmid=',dtmid,' dtbnd='
     &, dtbnd,' dtomid=',dtomid, ' time=',time, ' ia=',ia,' ib=',ib
     &, ' wc=',wc, ' change=',change, ' iprev=',iprev, ' inext=',inext
#endif
      return
      end



      subroutine checkinterp (ntdrec, tdrec, aprec)
c              
c=======================================================================
c     check for consistency between interpolation period centers "tdrec"
c     and period lengths "aprec".
c     adjust tdrec and aprec for leap years
c     check for and compensate for some mismatches between data and 
c     calendar
c
c     authors: c. h. goldberg           chg@gfdl.gov
c              r. c. pacanowski         rcp@gfdl.gov
c=======================================================================
c
#include "calendar.h"
#include "stdunits.h"
      dimension tdrec(ntdrec), aprec(ntdrec)
      logical febdone, monthly, error
c
c     test for consistency of tdrec and aprec times
c
      monthly = .true.
      error   = .false.
      sum     = 0.5*aprec(1)
      do m=2,ntdrec
        sum = sum + 0.5*(aprec(m) + aprec(m-1))
	if (abs(tdrec(m) - sum) .gt. 0.01*tdrec(m)) then
	  error = .true.
	  write (stdout,*) 'Error in time interpolation data'
	  write (stdout,*) 'Date for middle of record ',m
     &,                   ' is not centered'
	endif
	if (.not.(28.0 .le. aprec(m) .and. aprec(m) .le. 32)) then
	  monthly = .false.
	endif
      enddo
c
      if (error) then
        write (stdout,*) 'STOP in checkinterp'
c	stop
      endif
c
      dlen = (tdrec(ntdrec) + 0.5*aprec(ntdrec)) - 
     &       (tdrec(1) - 0.5*aprec(1))
c
c     if using leap years, add 1/4 day to feburary (or last record in
c     feb if data is other than monthly. eg: daily)
c
      if (.not.eqyear) then
        if (mod(dlen, real(yrlen)) .lt. 0.01) then
c
c         calendar has leap years but data does not, add 1/4 day to
c         feburary (or last record in feb if data is other than monthly.
c         eg: daily)
c     
          write (stdout, '(/,a,a)')
     &           'Checkinterp: Modifying equal year interpolation' 
     &,          ' data for use with leap year calendar'
          febdone = .false.
	  time    = 0.0
          do m=1,ntdrec
	    time = time + aprec(m)
	    if (time .ge. yrlen) then
	      time = time - yrlen
	      febdone = .false.
	    endif
	    if (time .ge. msum(3)) then
	      if (.not. febdone) then
	        aprec(m)   = aprec(m) + 0.2425
                write (stdout, '(a,i4)')
     &                'Checkinterp: Adding 0.2425 days to record ',m
	        febdone = .true.
	      endif
	    endif
	  enddo
          sum = tdrec(1) - 0.5*aprec(1)
          do m=1,ntdrec
            sum = sum + aprec(m)
            tdrec(m) = sum - 0.5*aprec(m)
c	    print *,' m=',m,' tdrec=',tdrec(m), ' aprec=',aprec(m)
          enddo
c	  
        elseif (mod(dlen, real(yrlen)) - 0.2425*nint(dlen/yrlen)
     &          .le. 0.01) then
c
c         calendar has leap years and data is leap year corrected by adding
c         0.2425 days per year.  interpolation data is consistent.
c
        else
c
c         calendar has leap years but data is neither leap year
c         compensated nor an exact number of years, it is not clear
c         what user wants.
c
          write (stdout,*) 'Problem in checkinterp'
	  write (stdout,*) 'Calendar uses leap years, but interpolation'
     &,	     ' data is neither an integer number of years or leap year'
     &,		   ' corrected by adding 0.2425 days per year.'
          stop
	endif
c
      else
        if (mod(dlen, real(yrlen)) .lt. 0.01) then
c
c         calendar uses equal years and data is an integral number of these
c         years.  interpolation data is consistent.
c       
        elseif (mod(dlen, real(yrlen)) - 0.2425*nint(dlen/yrlen)
     &          .le. 0.01) then
c
c         calendar uses equal years, but data is leap year corrected.
c         subtract 1/4 day from feburary (or last record in feb if data is other
c         than monthly. eg: daily)
c     
          write (stdout, '(/,a,a)')
     &                     'Checkinterp: Modifying leap year corrected'
     &,                  ' interpolation data for use with equal years'
          febdone = .false.
	  time    = 0.0
          do m=1,ntdrec
	    time = time + aprec(m)
	    if (time .ge. yrlen + 0.2425) then
	      time = time - yrlen - 0.2425
	      febdone = .false.
	    endif
	    if (time .ge. msum(3)+0.2425) then
	      if (.not. febdone) then
	        aprec(m)   = aprec(m) - 0.2425
                write (stdout, '(a,i4)')
     &           'Checkinterp: Subtracting 0.2425 days from record ',m
	        febdone = .true.
	      endif
	    endif
	  enddo
          sum = tdrec(1) - 0.5*aprec(1)
          do m=1,ntdrec
            sum = sum + aprec(m)
            tdrec(m) = sum - 0.5*aprec(m)
c	    print *,' m=',m,' tdrec=',tdrec(m), ' aprec=',aprec(m)
          enddo
c	  
        else
c
c         calendar has equal years but data is neither leap year
c         compensated nor an exact number of years, it is not clear
c         what user wants.
c
          write (stdout,*) 'Problem in checkinterp'
	  write (stdout,*) 'Calendar uses equal years, but interpolation'
     &,	     ' data is neither an integer number of years or leap year'
     &,		   ' corrected by adding 0.2425 days per year.'
          stop
c
        endif          
c
      endif
      return
      end


