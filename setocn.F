      subroutine setocn
c
c=======================================================================
c
c     set up everything which must be done only once per run
c
c=======================================================================
c
      logical error, vmixset, hmixset
      parameter (ifdmax=100)
      character*30 cifdef(ifdmax)
c
#include "param.h"
#include "accel.h"
#include "calendar.h"
#ifdef sponges
# include "sponge.h"
#else
      logical annlev
#endif
#include "coord.h"
#if defined fourfil || defined firfil
# include "cpolar.h"
#endif
#include "cprnts.h"
#include "cregin.h"
#include "csbc.h"
#ifdef shortwave
# include "cshort.h"
#endif
#include "csnap.h"
#include "ctmb.h"
#include "diag.h"
#include "docnam.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "index.h"
#include "iounit.h"
#ifdef isopycmix
# include "isopyc.h"
#endif
#include "levind.h"
#include "mw.h"
#ifdef trajectories
# include "ptraj.h"
#else
      logical initpt
#endif
#include "scalar.h"
#include "stab.h"
#include "state.h"
#include "switch.h"
#include "tmngr.h"
#include "topog.h"
#include "vmixc.h"
#if defined fourfil || defined firfil
      dimension kmz(imt,jmt)
#endif
      character*80 momver
c
      namelist /contrl/ init, runlen, rununits, segtim, restrt, initpt
      namelist /mbcin/  mapsbc, sbcname, dunits, coabc, crits
     &,                 numpas, bwidth
      namelist /tsteps/ dtts, dtuv, dtsf
      namelist /riglid/ mxscan, sor, tolrsf, tolrsp, tolrfs
      namelist /mixing/ am, ah, ambi, ahbi
     &,                 kappa_m, kappa_h, cdbot
     &,                 aq, aidif
     &,                 ncon, nmix, eb, acor, dampts, dampdz, annlev
     &,                 alph, gam, theta
      namelist /diagn/  tsiint, tavgint, itavg, tmbint, itmb, stabint
     &,                 zmbcint, glenint, trmbint, itrmb, vmsfint
     &,                 gyreint, extint, prxzint, trajint, exconvint
     &,                 dspint, snapint, snapls, snaple, snapde
     &,                 timavgint, cmixint 
     &,                 prlat, prslon, prelon, prsdpt, predpt
     &,                 slatxy, elatxy, slonxy, elonxy
     &,                 cflons, cflone, cflats, cflate, cfldps, cfldpe
     &,                 maxcfl, xbtint
      namelist /io/     expnam, iotavg 
     &,                 iotmb, iotrmb
     &,                 ioglen, iovmsf, iogyre
     &,                 ioprxz, ioext, iodsp, iotsi, iozmbc, iotraj
     &,                 ioxbt
      namelist /ictime/ year0, month0, day0, hour0, min0, sec0
     &,                 ryear, rmonth, rday, rhour, rmin, rsec
     &,                 refrun, refinit, refuser, eqyear, eqmon, monlen
c
c-----------------------------------------------------------------------
c                                INITIALIZATIONS
c
c     grep xxx *.h to find the description (where xxx is what you want)
c
c-----------------------------------------------------------------------
#ifdef timing
      call tic ('driver', 'setocn')
#endif
c
c     is the MW (memory window) opened all the way or not?
c
      if (jmw .eq. jmt) then
        wide_open_mw = .true.
      else
        wide_open_mw = .false.
      endif
c
c-----------------------------------------------------------------------
c     time at initial conditions.
c-----------------------------------------------------------------------
c
      year0  =1
      month0 =1
      day0   =1
      hour0  =0
      min0   =0
      sec0   =0
c
c-----------------------------------------------------------------------
c     set whether calculations for logical switches are referenced to
c
c     refrun  = T ==>  the start of each run
c     refinit = T ==>  initial condition time given by:
c                     year0, month0, day0, hour0, min0, sec0
c     refuser = T ==>  user specified reference time
c                     ryear, rmonth, rday, rhour, rmin, rsec
c
c     choose by setting one of the above to TRUE
c-----------------------------------------------------------------------
c
      refrun  = .false.
      refinit = .true.
      refuser = .false.
c
c     set the date and time for referencing switches
c     (only needed if rfuser = .T.)
c
      ryear  = 1900
      rmonth = 1
      rday   = 1
      rhour  = 0
      rmin   = 0
      rsec   = 0
c
c-----------------------------------------------------------------------
c     set type of year (eqyear = F => real calendar with leap years)
c                      (eqyear = T => idealized calendar)
c
c     if calendar is idealized then there are two choices:
c     eqmon = F => the usual number of days per month (31, 28, 31, ...)
c                  and the length of a year is 365 days
c                  (monlen is not used in this case)
c     eqmon = T => a fixed number of days per month set by "monlen"
c                  and the length of a year is 12*monlen
c-----------------------------------------------------------------------
c
      eqyear   = .false.
      eqmon    = .false.
      monlen   = 30
c
c     stability diagnostic
c
      call stabi
c
c     misc variables
c
      error    = .false.
      vmixset  = .false.
      hmixset  = .false.
      momver   = ' GFDL MOM 2.0'
      expnam   = ' MOM test case  '
      init     =.true.
      initpt   =.true.
      bwidth   = 0.0
      numpas   = 100
c
c     initialize ifdef list
c
      do n=1,ifdmax
        cifdef(n) = '                              '
      enddo
c
      nifdef = 0
      cifdef(1) = 'no ifdefs enabled (turned on)'
c
c     integration time related variables 
c
      runlen   = 5.0
      rununits = 'days'
      segtim   = 1.0
      dtts     = -3600.0
      dtuv     = -3600.0
      dtsf     = -3600.0
      do k=1,km
        dtxcel(k) = 1.0
      enddo
c
c     control for leapfrog/mixing time steps
c
      nmix     = 17
      eb       = .true.
c
c     tolerance settings for poisson solvers
c
      tolrsf   = 1.0e8
      tolrsp   = 5.0e4
      tolrfs   = 5.0e4
c
c     vertical and horizontal mixing + newtonain damping variables
c
      kappa_h    = 1.0
      kappa_m    = 20.0
      am             = 1.0e9
      ah             = 2.0e7
      ambi           = 1.0e23
      ahbi           = 5.0e22
      cdbot          = 0.0
      visc_cbu_limit = 1.0e6
      diff_cbt_limit = 1.0e6
c
c     damping time scales and depths for test case
c
      do n=1,nt
        dampts(n)    = 50.0
        dampdz(n)    = 26.575e2
      enddo
      annlev         = .false.
c
c     number of convective passes for standard explicit convection
c
      ncon     = 1
c
c     physical constants
c
      rho0     = 1.035
      rho0r    = c1/rho0
      grav     = 980.6
      radius   = 6370.0e5
c
c     time centering for equations
c
#if defined implicitvmix || defined isopycmix
      aidif    = 0.5
#else
      aidif    = 0.0
#endif
c
#ifdef rigid_lid_surface_pressure
      alph     = 1.0
      gam      = 0.0
      theta    = 1.0
#endif
c
#ifdef implicit_free_surface
      alph     = 0.3333333
      gam      = 0.3333333
      theta    = 0.5
#endif
c
#ifdef linearized_density
c
c     eliminate the pressure effect when using linearized density
c
      do k=1,km
        to(k) = 0.0
        so(k) = 0.0
      enddo
#endif
c
c
c
c-----------------------------------------------------------------------
c     write MOM version number
c-----------------------------------------------------------------------
c
      write (stdout,'(/1x,a14,a80/)') 'MOM version = ',momver
c
c-----------------------------------------------------------------------
c     provide for change in above presets using "namelist"
c-----------------------------------------------------------------------
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, contrl)
      write (stdout,contrl)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, mbcin)
      write (stdout,mbcin)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, tsteps)
      write (stdout,tsteps)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, riglid)
      write (stdout,riglid)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, mixing)
      write (stdout,mixing)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, diagn)
      write (stdout,diagn)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, io)
      write (stdout,io)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist', 'fsr')
      read  (ioun, ictime)
      write (stdout,ictime)
      call relunit (ioun)
c
      call getunit (iodoc, 'document.dta', 'fsr')
      write (iodoc, contrl)
      write (iodoc, mbcin)
      write (iodoc, tsteps)
      write (iodoc, riglid)
      write (iodoc, mixing)
      write (iodoc, diagn)
      write (iodoc, io)
      write (iodoc, ictime)
      write (iodoc,'(a,i6)') 'imt=', imt, 'jmt=', jmt, 'km=', km
     &,                       'nt=', nt, 'nvar=', nvar
c
c     user specified tracer names are placed into "trname" here.
c
      do m=1,nt
        trname(m) = '**unknown***'
      enddo
      trname(1) = 'potentl_temp'
      trname(2) = 'salinity    '
      write (iodoc,'(a,a)') 'trname(1)=', trname(1)
     &,                     'trname(2)=', trname(2)
      call relunit (iodoc)
c
c-----------------------------------------------------------------------
c     open i/o units needed for prognostic variables
c     nbuf = number of i/o buffers (currently not used) 
c-----------------------------------------------------------------------
c
#if defined ramdrive 
      call getunitnumber (kflds)
      call getunitnumber (latdisk(1))
      call getunitnumber (latdisk(2))
#endif
c
      nbuf = 2
      call ostart (kflds, 'kflds', nkflds*nwds, nwds, nbuf)
      call ostart (latdisk(1), 'latdisk1', jmt*nslab, nslab, nbuf)
      call ostart (latdisk(2), 'latdisk2', jmt*nslab, nslab, nbuf)
c
c-----------------------------------------------------------------------
c     set up the grids in x (longitude), y (latitude), and z (depth)
c     corresponding to Arakawa "b" gird system
c-----------------------------------------------------------------------
c
      call grids
c
c-----------------------------------------------------------------------
c     check that the density coefficients are consistent with the
c     depths generated by grids
c-----------------------------------------------------------------------
c
      cksum = checksum (zt, km, 1)
      if (abs(cksum - cksumzt) .gt. 1.e-6*cksumzt )then
        write (stdout,'(a/a)')
     & '=>Error: density coeffs were calculated for depth levels that'
     &,'         differ from those used in this model run'
        write (stdout,*) 'checksum from setocn   = ',cksum
	write (stdout,*) 'checksum from dncoef.h = ',cksumzt
        stop 'setocn'
      endif
c
#ifdef linearized_advection
c
c     set the initial idealized stratification as a function of
c     temperature only.
c
      tzero = 7.5
      tone  = 10.0
      z0    = 30.0e2
      bigl  = 80.0e2
      write (stdout,'(//a/)')
     & 'Temperature stratification for linearized advection'
      do k=1,km
        tbarz(k) = tzero*(1.0-tanh((zt(k)-bigl)/z0)) + 
     &              tone*(1.0-zt(k)/zt(km))
        write (stdout,'(1x,"k=",i3," tbarz=",e14.7)') k, tbarz(k) 
      enddo     
#endif
c
c-----------------------------------------------------------------------
c     set the earth`s rotation rate to computer precision
c-----------------------------------------------------------------------
c
      omega  = pi/43082.0
      call getunit (iodoc, 'document.dta', 'fsa')
      write (iodoc,'(a,e14.7)') 'omega=', omega
      call relunit (iodoc)
c
c-----------------------------------------------------------------------
c     if the MW is not fully opened, then set time level indicators in
c     the MW ("tau-1" "tau" "tau+1") to constant values.
c-----------------------------------------------------------------------
c
      if (.not. wide_open_mw) then
        taum1 = -1
        tau   =  0
        taup1 = +1
      endif
c
c-----------------------------------------------------------------------
c     set prognostic quantities to either initial conditions or restart
c-----------------------------------------------------------------------
c
      if (init) then
c
c       generate topographic levels "kmt" on "t" cells.
c
        call topog (kmt, map, xt, yt, zt, xu, yu, zw, imt, jmt, km)
c
c       construct depth arrays associated with "u" cells
c
        call depth_u (kmt, imt, jmt, zw, km, kmu, h, hr)
c
c       initialize two dimensional fields on disk
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c       initialize surface pressure fields in memory
c
        do jrow=1,jmt
          do i=1,imt
            ps(i,jrow,1)     = c0
	    ps(i,jrow,2)     = c0
	    pguess(i,jrow)   = c0
	    ubar(i,jrow,1)   = c0
	    ubar(i,jrow,2)   = c0
	    ubarm1(i,jrow,1) = c0
	    ubarm1(i,jrow,2) = c0
          enddo
	enddo
#endif
#ifdef stream_function
c
c       initialize stream function fields in memory
c
        do n=1,2
          do jrow=1,jmt
            do i=1,imt
              psi(i,jrow,n) = c0
	      ptd(i,jrow)   = c0
            enddo
          enddo
	enddo
c
c       initialize stream function guess fields on disk
c       block`s 1 & 2 are for the stream function guess field on disk
c
        do n=1,nkflds
          call oput (kflds, nwds, n, ptd(1,1))
	enddo
#endif
c
c       initialize all latitude rows
c
        call rowi
c
c       initialize time step counter = 0
c
        itt    = 0
	irstdy = 0
	msrsdy = 0
      else
c
c       initialize all prognostic quantities from the restart
c
        call read_restart
c
c       compute a topography checksum
c
        cksum = 0.0
        do jrow=1,jmt
          do i=1,imt
	    cksum = cksum + i*jrow*kmt(i,jrow)
	  enddo
        enddo
        write (stdout,*) ' "kmt" checksum = ', cksum
      endif
c
c-----------------------------------------------------------------------
c     initialize the time manager with specified initial conditions
c     time, user reference time, model time, and how long to integrate.
c-----------------------------------------------------------------------
c
      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,              ryear, rmonth, rday, rhour, rmin, rsec
     &,              irstdy, msrsdy
     &,              runlen, rununits, rundays, dtts)
c
#ifdef stability_tests
c
c-----------------------------------------------------------------------
c     convert starting and ending longitudes for the stability tests
c     to nearest model grid points.
c-----------------------------------------------------------------------
c
      if (stabint .ge. c0) then
        iscfl  = max(indp (cflons, xt, imt), 2)
        cflons = xt(iscfl)
        iecfl  = min(indp (cflone, xt, imt), imt-1)
        cflone = xt(iecfl)
        jscfl  = max(indp (cflats, yt, jmt), 2)
        cflats = yt(jscfl)
        jecfl  = min(indp (cflate, yt, jmt), jmt-1)
        cflate = yt(jecfl)
        kscfl  = indp (cfldps, zt, km)
        cfldps = zt(kscfl)
        kecfl  = indp (cfldpe, zt, km)
        cfldpe = zt(kecfl)
      endif
#endif
c
#ifdef restorst
c
c-----------------------------------------------------------------------
c      damp surface tracers back to data. schematically, the restoring
c      term will be = (dampdz/(dampts*86400)) * (data - t)
c      where dampdz is the thicness (cm) and dampts is the time
c      scale (days)
c-----------------------------------------------------------------------
c
      do n=1,nt
        write (stdout,'(/,1x,a,i2,a,a,/,1x,1pe14.7,a,1pe14.7,a,/)')
     &  ' Surface tracer #',n,' will be damped back to data using a'
     &, ' Newtonian restoring time scale of '
     &,  dampts(n),' days. and a level thickness =', dampdz(n),' cm.'
      enddo
      call getunit (iodoc, 'document.dta', 'fsa')
      write (iodoc,'(a,(1x,e14.7))') 'dampts=', dampts
      write (iodoc,'(a,(1x,e14.7))') 'dampdz=',dampdz
      call relunit (iodoc)
#endif
c
#if defined fourfil || defined firfil
c
c     compute sin and cos values for vector correction before filtering
c
      fx =  1.0e-10
      fxa = dxt(1)/radius
c
      do i=2,imtm1
        fxb = fxa*float(i-2)
        spsin(i) = sin(fxb)
        spcos(i) = cos(fxb)
        if (abs(spsin(i)) .lt. fx) spsin(i) = c0
        if (abs(spcos(i)) .lt. fx) spcos(i) = c0
      enddo
c
      spsin(1) = c0
      spcos(1) = c0
      spsin(imt) = c0
      spcos(imt) = c0
c
c     set up model indices for filtering high latitudes
c
      jfrst  = indp (rjfrst, yt, jmt)
      jft0   = indp (rjft0, yt, jmt)
      jft1   = indp (rjft1, yt, jmt)
      jft2   = indp (rjft2, yt, jmt)
      jfu0   = indp (rjfu0, yu, jmt)
      jfu1   = indp (rjfu1, yu, jmt)
      jfu2   = indp (rjfu2, yu, jmt)
      call getunit (iodoc, 'document.dta', 'fsa')
      write (iodoc,'(a,i6)') 'jfrst=', jfrst, 'jft0=',jft0
     &, 'jft1=',jft1, 'jft2=',jft2, 'jfu0=',jfu0, 'jfu1=',jfu1
     &, 'jfu2=',jfu2
      call relunit (iodoc)
      jskpt  = jft2-jft1
      jskpu  = jfu2-jfu1
      njtbft = (jft1-jfrst+1)+(jmtm1-jft2+1)
      njtbfu = (jfu1-jfrst+1)+(jmtm1-jfu2+1)
      if (njtbft .gt. jmtfil .or. njtbfu .gt. jmtfil) then
        write (stdout,9599) njtbft, njtbfu
        write (stderr,9599) njtbft, njtbfu
        stop '>ocean 1'
      endif
9551  format (/' ==== start and end indices for fourier filtering ====')
9552  format (' only 11 sets of indices fit across the page.',
     &       '  others will not be printed.'/)
9553  format (/,' == filtering indices for t,s ==')
9554  format (/,' == filtering indices for u,v ==')
9555  format (/,' == filtering indices for stream function ==')
9599  format (/,' error => jmtfil must be >= max(njtbft,njtbfu)',
     &        /,'          njtbft=',i8,' njtbfu=',i8)
c
#  ifdef firfil
c
c     set "numflt" and "numflu" to filter more at higher latitudes
c     using "jft0" and "jfu0" as reference latitude rows.
c
      numfmx = imt * p25
      refcos = cst(jft0)
      write(stdout,9501) refcos, yt(jft0)
      do jrow=jfrst,jmt
        if ((jrow .le. jft1 .or. jrow .ge. jft2) .and. jrow .ge. jfrst)
     &    then
          jj = jrow - jfrst + 1
          if (jrow .ge. jft2) jj = jj - jskpt + 1
          numflt(jj) = max(1,int(refcos/cst(jrow)))
          if (numflt(jj) .gt. numfmx) numflt(jj) = numfmx
          write(stdout,9502) jrow, jj, numflt(jj), yt(jrow)
	  if (jj .gt. jmtfil) then
	    write (stdout,'(1x,a,i4,a,i4,a)')
     &      'Error: jj exceeds jmtfil. jj=',jj, ' jmtfil=',jmtfil
     &,     '.  increase jmtfil'
            stop "=>setocn"
	  endif
        endif
      enddo
      refcos = csu(jfu0)
      write(stdout,9503) refcos, yu(jfu0)
      do jrow=jfrst,jmt
        if ((jrow .le. jfu1 .or. jrow .ge. jfu2) .and. jrow .ge. jfrst)
     &    then
          jj = jrow - jfrst + 1
          if (jrow .ge. jfu2) jj = jj - jskpu + 1
          numflu(jj) = max(1,int(refcos/csu(jrow)))
          if (numflu(jj) .gt. numfmx) numflu(jj) = numfmx
          write(stdout,9502) jrow, jj, numflu(jj), yu(jrow)
	  if (jj .gt. jmtfil) then
	    write (stdout,'(1x,a,i4,a,i4,a)')
     &      'Error: jj exceeds jmtfil. jj=',jj, ' jmtfil=',jmtfil
     &,     '.  increase jmtfil'
            stop "=>setocn"
	  endif
        endif
      enddo
c
9501  format(/' firfil reference cosine for tracers =',e12.6,' (',
     &       f7.3,' deg)'/'     jrow      jj numflt(jj)  latitude')
9502  format(1x,3i8,6x,f7.3)
9503  format(/' firfil reference cosine for velocities =',e12.6,' (',
     &       f7.3,' deg)'/'     jrow      jj numflu(jj)  latitude')
#  endif
#endif
c
#ifdef sponges
c
c-----------------------------------------------------------------------
c     set latitude functions for newtonian damping term in sponge layers
c     near artificial northern and southern boundaries in limited
c     domain models.
c     schematically: damping = - sponge(j) * (t(i,k,j) - levitus(i,k,j))
c     all related data is assumed to have been prepared using the
c     "sponge" routines in the MOM dataset.
c     disk resource factor of 13 is for 12 months + 1 annual mean
c-----------------------------------------------------------------------
c
      lrec = 4 + 4 + 2*jmt + 4*imt*km
      write (opt_sponge,'(a,1x,i8)') 'u d words =',lrec
      write (stdout,'(/a,1pg10.3,a)')
     & ' Sequential access disk resource for file "sponge.mom" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew1, 'sponge.mom', 'u s r ieee')
      write (stdout,'(/a,1pg10.3,a)')
     & ' Direct access disk resource for file "sponges" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew2, 'sponges', opt_sponge) 
c
      sum = 0.0
      cksum = 0.0
      do m=1,13
        read (ionew1)
        read (ionew1) sstamp, spdpm, im, kk, jm, j1, j2, mm
     &, spngs, spngn, spbuf1
        write (ionew2, rec=m) sstamp, spdpm, im, kk, jm, spngs, spngn
     &, spbuf1
        if (m .le. 12) then
          spgdpm(m) = spdpm
	  sum = sum + spdpm
	  tspng(m) = sum - 0.5*spdpm
	endif
	do n=1,4
	  cksum = cksum + checksum (spbuf1(1,1,n), imt, km)
	enddo
      enddo
      print *,' checksum for sponge data = ',cksum
      if (annlev) then
        write (stdout,'(/a)')
     &  ' => Annual mean levitus data will be used for Newtonain sponge'
        do m=1,12
	  write (ionew2, rec=m) sstamp, spdpm, im, kk, jm, spngs, spngn
     &,   spbuf1
	enddo
      else
        write (stdout,'(/a)')
     &  ' => Monthly levitus data will be used for Newtonain sponge'
      endif
      write (stdout,'(1a/)')  ' Newtonian damping zone setup:'
      do jrow=1,jmt
        if (spngs(jrow) .ne. c0) then
          write (stdout,'(a,i3,a,f7.2,a,e10.3,a,e10.3)') ' jrow='
     &,   jrow,', yt(jrow)=',yt(jrow), ', spngs(jrow) =',spngs(jrow)
     &,   ', Newtonain time scale (days)=',secday/spngs(jrow)
        endif
        if (spngn(jrow) .ne. c0) then
          write (stdout,'(a,i3,a,f7.2,a,e10.3,a,e10.3)') ' jrow='
     &,   jrow,', yt(jrow)=',yt(jrow), ', spngn(jrow) =',spngn(jrow)
     &,   ', Newtonain time scale (days)=',secday/spngn(jrow)
        endif
        if (spngn(jrow) .ne. c0 .and. spngs(jrow) .ne. c0) then
	  write (stdout,'(/a/)')
     &    ' Error: Overlapping north and south sponges not allowed'
          stop '=>setocn'
	endif
      enddo
      write (stdout,*) ' '
c
c     set the current model time in days and initialize interpolation
c     information
c
      begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
c
      if (.not.eqyear) then
        write (stdout,*) '=>Warning leap year being used with '
     &, ' climatological sponges?'
      endif
      iprev  = 1
      inext  = 2
      indxsp = 5
      method = 3
      call timeinterp(begtim, indxsp, tspng, spgdpm, 12, .true., method
     &,               inextd, iprevd, wprev, readsp, inext, iprev)
c
      read (ionew2, rec=iprevd) stprev, spdpmp, im, kk, jm, spngs
     &,    spngn, spbuf1
      do n=1,4
        do k=1,km
          do i=1,imt
            spbuf(i,k,n,iprev) = spbuf1(i,k,n)
          enddo
        enddo
      enddo
      read (ionew2, rec=inextd) stnext, spdpmn, im, kk, jm, spngs
     &,    spngn, spbuf1 
      do n=1,4
        do k=1,km
          do i=1,imt
            spbuf(i,k,n,inext) = spbuf1(i,k,n)
          enddo
        enddo
      enddo
c
      write (stdout,'(2(/a,i3,1x,a,a,i2)/a,g14.7,1x,a,a,g14.7/)')
     & ' reading sponge record ', iprevd, stprev, ' into buffer ', iprev
     &,' reading sponge record ', inextd, stnext, ' into buffer ', inext
     &,' for day =', begtim, stamp, ' weight =',wprev
      write (stdout,'(/a,i2,a,i2/)') 'sponge is dataset index ',indxsp
     &,' for time interpolation using method #',method
      call relunit (ionew1)
      call relunit (ionew2)
#endif
c
#ifdef shortwave
c
c-----------------------------------------------------------------------
c     Solar Shortwave energy penetrates below the ocean surface. Clear
c     water assumes energy partitions between two exponentials as 
c     follows:
c     
c     58% of the energy decays with a 35 cm e-folding scale
c     42% of the energy decays with a 23 m e-folding scale
c
c     if the thickness of the first ocean level "dzt(1)" is 50 meters,
c     then shortwave penetration wouldn't matter. however, for 
c     dzt(1) = 10 meters, the effect can be significant. this can be
c     particularly noticable in the summer hemisphere.
c
c     Paulson and Simpson (1977 Irradiance measurements in the upper
c                               ocean JPO 7, 952-956)
c     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
c                  A General Circulation Model for Upper Ocean
c                  Simulaton (Rosati and Miyakoda JPO vol 18,Nov 1988)
c-----------------------------------------------------------------------
c
      write (stdout,'(/a/)') 
     &' => Shortwave penetration is a double exponential as follows:'
      rpart  = 0.58
      efold1 = 35.0e0
      efold2 = 23.0e2
      rscl1  = 1.0/efold1
      rscl2  = 1.0/efold2
      call getunit (iodoc, 'document.dta', 'fsa')
      write (iodoc,'(a,e14.7)') 'rpart=', rpart, 'efold1=',efold1
     &, 'efold2=',efold2
      call relunit (iodoc)
c
c     note that pen(0) is set 0.0 instead of 1.0 to compensate for the
c     shortwave part of the total surface flux in "stf(i,1)"
c
      pen(0) = c0
c
      do k=1,km
        swarg1 = -min(zw(k)*rscl1,70.0)
        swarg2 = -min(zw(k)*rscl2,70.0)
        pen(k) = rpart*exp(swarg1) + (c1-rpart)*exp(swarg2)
	divpen(k) = (pen(k-1) - pen(k))/dzt(k)
	write (stdout,9200) k, zw(k)*1.e-2, pen(k), zt(k)*1.e-2
     &,                     divpen(k)
      enddo
      write (stdout,*) ' '
9200  format (1x,'k=',i3,' zw(k)=',f8.2,'(m)  pen(k)=',e14.7
     &,       ' zt(k)=',f8.2,'(m)  divpen(k)=',e14.7)
#endif
c
c-----------------------------------------------------------------------
c     compute surface area and volume of ocean ("t" cells and "u" cells)
c     (note that areas are defined at each level)
c-----------------------------------------------------------------------
c
      do k=1,km
        tcella(k) = c0
        ucella(k) = c0
      enddo
      ocnp   = c0
      tcellv = c0
      ucellv = c0
c
c     this comment directive turns off autotasking on the YMP
c     for the following loop
c
cfpp$ noconcur l
      do jrow=2,jmtm1
        do i=2,imtm1
          if (kmt(i,jrow) .gt. 0) then
	    do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) + cst(jrow)*dxt(i)*dyt(jrow)
	    enddo
            tcellv = tcellv + cst(jrow)*dxt(i)*dyt(jrow)*zw(kmt(i,jrow))
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
          if (kmu(i,jrow) .gt. 0) then
	    do k=1,kmu(i,jrow)
              ucella(k) = ucella(k) + csu(jrow)*dxu(i)*dyu(jrow)
	    enddo
            ucellv = ucellv + csu(jrow)*dxu(i)*dyu(jrow)*zw(kmu(i,jrow))
          endif
        enddo
      enddo
c
      write (stdout,9341) tcella(1), tcellv, ucella(1), ucellv
c
#if defined tracer_averages || defined term_balances || gyre_components
c
c---------------------------------------------------------------------
c     set the horizontal regional masks and names to be used when
c     computing averages on the "t" grid in subroutine "region.F".
c     also set the vertical regional masks and names for use (along
c     with the horizontal ones) in term balance calculations for
c     tracers & momentum. For term balance calculations the number
c     of masks is the product of horizonatl & vertical regions.
c---------------------------------------------------------------------
c
# ifdef readrmsk
c
c    read in horizontal & vertical regional masks ("mskhr" & "mskvr")
c    and names ("hregnm" & "vregnm") on unit iormsk
c
      call getunit (iormsk, 'region_masks', 'fsr')
      call reg1st (iormsk, .true., .false., .false., .false., .true.)
      call relunit (iormsk)
# else
c
c     set up the horizontal regions:
c     specify "mskhr" and "hregnm" values rather than reading them in
c     from a file (arbitrarily defaulted below for 5 zonal bands).
c     note: there must be "nhreg" calls ... one for each horizontal
c     region. The form is:
c     call sethr (region number, starting lon, ending lon, starting lat,
c                 ending lat)
c     where the lon and lat limits are for the edges of the "t" cells
c     "sethr" will fit the region using the nearest model grid points
c
      do jrow=1,jmt
        do i=1,imt
          mskhr(i,jrow) = 0
        enddo
      enddo
c
      call sethr (1, xu(1), xu(imtm1), -90.0, -60.1)
      call sethr (2, xu(1), xu(imtm1), -60.1, -24.9)
      call sethr (3, xu(1), xu(imtm1), -24.9,  24.9)
      call sethr (4, xu(1), xu(imtm1),  24.9,  59.9)
      call sethr (5, xu(1), xu(imtm1),  59.9,  90.0)
c
      do jrow=1,jmt
        mskhr(1,jrow)   = 0
        mskhr(imt,jrow) = 0
      enddo
c
c     set up the vertical regions:
c     specify "mskvr" and "vregnm" values. also the starting & ending
c     levels for the vertical regions (arbitrarily defaulted for two
c     vertical regions). Note: there must be "nvreg" calls... one for
c     each vertical region. 
c     The form is:
c     call setvr (region number, starting depth, ending depth)
c     where the depth limits are for the edges of the "t" cells
c     "setvr" will fit the region using the nearest model grid points
c
      n = 1
      call setvr (n, 0.0e2, 245.0e2)
      n = n + 1
      call setvr (n, 245.0e2, 1570.e2)
c
c     Note: additional vertical regions can be added as follows:
c
c      n = n + 1
c      call setvr (n, 1570.0e2, 2500.e2)
c
      if (n .gt. nvreg) then
        write (stdout,*) "=>Error: increase parameter 'nvreg' "
	write (stdout,*) "         stop in setocn"
	stop
      endif
c
      do k=1,km
        mskvr(k) = 0
      enddo
      do n=1,nvreg
        ks = llvreg(n,1)
        ke = llvreg(n,2)
        do k=ks,ke
          mskvr(k) = n
        enddo
      enddo
# endif
#ifdef tracer_averages
c
c     print out regional mask info for tracer average diagnostic
c
      if ((iotavg .ne. stdout .or. iotavg .lt. 0) .and. itavg) then
	call getunit (iou, 'tracer_avg.dta','u s a ieee')
        call reg1st (iou, .false., .false., .false., .true., .false.)
        call relunit (iou)
      endif
      if ((iotavg .eq. stdout .or. iotavg .lt. 0) .and. itavg) then
        call reg1st (stdout, .true., .false., .false., .true., .false.)
      endif
#endif
#ifdef term_balances
c
c     print out regional mask info for term balance diagnostic
c
      if ((iotrmb .ne. stdout .or. iotrmb .lt. 0) .and. itrmb) then
	call getunit (iou, 'term_bal.dta','u s a ieee')
        call reg1st (iou, .false., .false., .false., .true., .false.)
        call relunit (iou)
      endif
#endif
c
c-----------------------------------------------------------------------
c     compute the surface area and volume of the ocean regions. index 0
c     refers to the sum of all regions.
c     (values used in subroutine region are done in terms of meters,
c     rather than centimeters)
c-----------------------------------------------------------------------
c
      areag = c0
      volgt  = c0
c
      do k=1,km
        volgk(k) = c0
      enddo
c
      do n=0,numreg
        areat(n) = c0
        areau(n) = c0
        volt(n)  = c0
        volu(n)  = c0
      enddo
      do mask=1,nhreg
        areab(mask) = c0
        volbt(mask) = c0
        do k=1,km
          volbk(mask,k) = c0
        enddo
      enddo
c
      do jrow=2,jmtm1
        do i=2,imtm1
          mask = mskhr(i,jrow)
          kz   = kmt(i,jrow)
          if (kz .gt. 0 .and. mask .gt. 0) then
            boxat = cst(jrow) * dxt(i) * dyt(jrow)
            if (kmu(i,jrow) .ne. 0) then
              boxau = csu(jrow) * dxu(i) * dyu(jrow)
            else
              boxau = c0
            endif
            areat(0) = areat(0) +  boxat
            areau(0) = areau(0) +  boxau
	    areab(mask)  = areab(mask) + boxat * 1.e-4
	    do k=1,kz
              volbk(mask,k) = volbk(mask,k) + boxat * dzt(k) * 1.e-6
              dvolt   = boxat * dzt(k)
              if (kmu(i,jrow) .ge. k) then
                dvolu   = boxau * dzt(k)
              else
                dvolu   = c0
              endif
              n = nhreg*(mskvr(k)-1) + mask
	      if (n .gt. 0) then
                volt(0) = volt(0) +  dvolt
                volu(0) = volu(0) +  dvolu
                volt(n) = volt(n) +  dvolt
                volu(n) = volu(n) +  dvolu
                do nv=1,nvreg
                  ks = llvreg(nv,1)
                  if (k .eq. ks) then
                    areat(n) = areat(n) +  boxat
                    areau(n) = areau(n) +  boxau
                  endif
                enddo
	      endif
            enddo
          endif
        enddo
      enddo
c
      do mask=0,numreg
        if (areat(mask) .eq. c0) then
          rareat(mask) = c0
        else
          rareat(mask) = c1/areat(mask)
        endif
c
        if (areau(mask) .eq. c0) then
          rareau(mask) = c0
        else
          rareau(mask) = c1/areau(mask)
        endif
c
        if (volt(mask) .eq. c0) then
          rvolt(mask) = c0
        else
          rvolt(mask) = c1/volt(mask)
        endif
c
        if (volu(mask) .eq. c0) then
          rvolu(mask) = c0
        else
          rvolu(mask) = c1/volu(mask)
        endif
      enddo
c
      do mask=1,nhreg
        do k=1,km
          volbt(mask) = volbt(mask) + volbk(mask,k)
          volgk(k) = volgk(k) + volbk(mask,k)
        enddo
        areag = areag + areab(mask)
        volgt = volgt + volbt(mask)
      enddo
c
      if (iotavg .ne. stdout .or. iotavg .lt. 0) then
	call getunit (iou, 'tracer_avg.dta','u s a ieee')
        call reg1st (iou, .false., .true., .true., .false., .false.)
	call relunit (iou)
      endif
      if (iotavg .eq. stdout .or. iotavg .lt. 0) then
        call reg1st (stdout, .true., .true., .true., .false., .false.)
      endif
c
c     compute and print statistics for regions
c
      sum = c0
      do n=1,numreg
        sum = sum + volt(n)
      enddo
      sum    = 100.0*sum/tcellv
      pctocn = 100.0*ocnp/float((imt-2)*(jmt-2)*km)
      diffa  = 100.0 * (c1 - (tcella(1) - 10000.0*areag)/tcella(1))
c
      write (stdout,9342) diffa, numreg, sum, pctocn
9342  format ('  the horizontal regional masks cover',f8.3
     &, '% of the total ocean surface area.'/
     &, '  there are ', i6, ' regions over which tracer & '
     &, 'momentum balances will be computed',/,'  accounting for '
     &, f6.2,'% of the total ocean volume.'/
     &, 1x,f6.2,'% of the grid points lie within the ocean.'/)
#endif
c
c-----------------------------------------------------------------------
c     find all land mass perimeters for poisson solvers
c-----------------------------------------------------------------------
c
      auto_kmt_changes = .false.
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
c     set mask for land mass perimeters on which to perform calculations
c     imask(-n) = .false.  [no equations ever on dry land mass n]
c     imask(0)  = .true.   [equations at all mid ocean points]
c     imask(n)  = .true./.false [controls whether there will be
c                                equations on the ocean perimeter of
c                                land mass n]
c     note: land mass 1 is the northwest-most land mass
c           usually includes the "north pole", and at low resolutions,
c           the "main continent"
c     for the numbering of the other landmasses, see generated map(i,j)
c           by selecting option -Dshow_details
c
      do isle=-mnisle,mnisle
        if (isle .ge. 0 .and. isle .le. nisle) then
          imask(isle) = .true.
        else
          imask(isle) = .false.
        end if
      end do
#ifdef symmetry
c
c     do not perform island integrals for land masses whose perimeters
c     touch the equator in models symmetric about the equator
c
      do isle=1,nisle
        do n=1,nippts(isle)
          j = jperm(iofs(isle)+n)
          if (j .eq. jmt-1) then 
	    imask(isle) = .false.
	  endif
        end do
      end do
#endif
c
c     user-specified changes to island mask
c       imask(1) = .true.
c       imask(2) = .true.
c
c     there are problems if imask is set .true. for a nonexistent
c     island.
c
c     print diagnostic information
c
      do isle=-mnisle,mnisle
        if (imask(isle)) then
	  if (isle .eq. 0) then
            print '(a)','=> calculations enabled for mid ocean points'
	  else
            print '(2a,i3)','=> calculations enabled for ocean ',
     &                      'perimeter of land mass',isle
	  end if
	end if
      end do
      do isle=0,nisle
        if (.not. imask(isle)) then
            print '(2a,i3)','=> calculations disabled for ocean ',
     &                      'perimeter of land mass',isle
	end if
      end do
c
c     imain is the land mass on which dpsi is normalized to 0
c     if imain is 0, then dpsi is not normalized.
c     default value of imain is land mass with longest perimeter
c
      imain = min(2,nisle)
      do isle=1,nisle
        if (nippts(isle) .gt. nippts(imain)) then
          imain = isle
        end if
      end do
c
c     if any island perimeter is not calculated, imain must be one such
c
      do isle=1,nisle
        if (.not.imask(isle)) then
          imain = isle
        end if
      enddo
#ifdef symmetry
c
c     do not normilize dpsi when using symmetry about the equator
c
      imain = 0
#endif
c
      if (imain .gt. 0 .and. imain .le. nisle) then
        print '(a,i4)', 'dpsi normalized to zero on land mass',imain
      else if (imain .eq. 0) then
        print *, 'no normalization on dpsi'
      else
        print *, 'ERROR: illegal value for choice of normalization ',
     &           'land mass, imain =', imain
      end if
      print *,' (user may set "imain" to any valid land mass number)'
c
c---------------------------------------------------------------------
c     compute checksum of density coefficients
c---------------------------------------------------------------------
c
      print *,' '
      call print_checksum (c(1,1), km, 9
     &,                   ' density coefficient checksum = ')
c
#if defined stream_function
c
c-----------------------------------------------------------------------
c     checksum the starting stream function.
c-----------------------------------------------------------------------
c
      call print_checksum (psi(1,1,1), imt, jmt
     &, ' checksum for psi(,,1) = ')
      call print_checksum (psi(1,1,2), imt, jmt
     &, ' checksum for psi(,,2) = ')
#endif
c
#if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c     compute an array to indicate "interior" stream function grid cells
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        kmz(1,jrow)   = 0
        kmz(imt,jrow) = 0
      enddo
c
      do i=1,imt
        kmz(i,1)   = 0
        kmz(i,jmt) = 0
      enddo
c
      do jrow=2,jmtm1
        do i=2,imt
          kmz(i,jrow) = min(kmu(i-1,jrow-1), kmu(i,jrow-1)
     &,                     kmu(i-1,jrow), kmu(i,jrow))
        enddo
      enddo
#  ifdef cyclic
      do jrow=1,jmt
        kmz(1,jrow) = kmz(imtm1,jrow)
      enddo
#  endif
c
c
c---------------------------------------------------------------------
c     find and print start & end indices for filtering
c---------------------------------------------------------------------
c
      write (stdout,9551)
      if (lsegf.gt.11) write (stdout,9552)
      write (stdout,9553)
      call findex (kmt, jmtfil, km, jft1, jft2, imt, istf, ietf)
      write (stdout,9554)
      call findex (kmu, jmtfil, km, jfu1, jfu2, imt, isuf, ieuf)
      write (stdout,9555)
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
      call findex (kmu, jmtfil, 1, jfu1, jfu2, imt, iszf, iezf)
# endif
# if defined stream_function
      call findex (kmz, jmtfil, 1, jft1, jft2, imt, iszf, iezf)
# endif
#endif
c
c---------------------------------------------------------------------
c     print the timestep multipliers
c---------------------------------------------------------------------
c
      write (stdout,9601) (dtxcel(k),k=1,km)
c
c-----------------------------------------------------------------------
c     initialize various things
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
	  uhat(i,jrow,1) = c0
	  uhat(i,jrow,2) = c0
	  divf(i,jrow) = c0
#endif
#ifdef stream_function
          ztd(i,jrow) = c0
#endif
          zu(i,jrow,1)  = c0
          zu(i,jrow,2)  = c0
        enddo
      enddo
c
c     coriolis factors
c
#if defined beta_plane || defined f_plane
c
c     beta plane with f = f0 + beta*y where f0 is at yu(1)
c     if f_plane then beta = 0
c
      degtcm = radius/radian
# ifdef f_plane
      dfdy   = c0
# else
      dfdy   = c2*omega*csu(1)/radius
# endif
      do jrow=1,jmt
        y            = (yu(jrow)-yu(1))*degtcm
        cori(jrow,1) = c2*omega*sine(1)  + y*dfdy
        cori(jrow,2) = -c2*omega*sine(1) + y*dfdy
      enddo
# ifdef f_plane
      write (stdout,'(//,a,1pg14.7,a,g14.7,a,g14.7//)')
     & ' "f plane" approximation set: f0 = ', cori(1,1)
     &, ' at latitude (deg) =', yu(1), ' with beta =',dfdy
# else
      write (stdout,'(//,a,1pg14.7,a,g14.7,a,g14.7//)')
     & ' "beta plane" approximation set: f0 = ', cori(1,1)
     &, ' at latitude (deg) =', yu(1), ' with beta =',dfdy
# endif
#else
      do jrow=1,jmt
        cori(jrow,1) = c2*omega*sine(jrow)
        cori(jrow,2) = -c2*omega*sine(jrow)
      enddo
#endif
c
c     metric diffusion factors
c
#if defined consthmix
# ifdef biharmonic
      amix = sqrt(abs(ambi))
# else
      amix = am
# endif
      do jrow=1,jmt
# if defined beta_plane || defined f_plane
c
c       set metric factors = 0 on beta and f plane
c       (strictly, they should be non zero... but they are small)
c
	am3(jrow)   = c0
	am4(jrow,1) = c0
	am4(jrow,2) = c0
	write (stdout,*)
     & '=>Warning: "beta_plane" and "f_plane" approximation sets '
     &, ' lateral mixing and non-linear metric terms to zero.'
# else
	am3(jrow)   = amix*(c1-tng(jrow)*tng(jrow))/(rho0*radius**2)
	am4(jrow,1) = -amix*c2*sine(jrow)/(rho0*radius*csu(jrow)
     &                                     *csu(jrow))
	am4(jrow,2) = -am4(jrow,1)
# endif
      enddo
#endif
c
c     metric advection factors
c
      do jrow=1,jmt
# if defined beta_plane || defined f_plane
c
c       set metric factors = 0 on beta and f plane
c       (strictly, they should be non zero... but they are small)
c
        advmet(jrow,1) = c0
	advmet(jrow,2) = c0
# else
        advmet(jrow,1) = tng(jrow)/radius
	advmet(jrow,2) = -advmet(jrow,1)
# endif
      enddo
c
c     diffusive flux through bottom of cells
c
      do j=jsmw,jemw
        do k=0,km
	  do i=1,imt
	    diff_fb(i,k,j) = c0
	  enddo
	enddo
      enddo	
c
c-----------------------------------------------------------------------
c     initialize diagnostics
c-----------------------------------------------------------------------
c
#ifdef meridional_tracer_budget
c
c-----------------------------------------------------------------------
c     set basins and initialize arrays
c-----------------------------------------------------------------------
c
      if (tmbint .ge. c0) then
c
c       set all points to signify basin # 1
c
        do jrow=1,jmt
          do i=2,imtm1
            msktmb(i,jrow) = 1
          enddo
          msktmb(1,jrow)   = 0
          msktmb(imt,jrow) = 0
        enddo
c
c       divide "msktmb" into other basins (2 .. ntmbb) here if desired.
c
c
c       set all land points to basin # 0
c
        do jrow=1,jmt
          do i=1,imt
            if (kmt(i,jrow) .eq. 0) msktmb(i,jrow) = 0
          enddo
        enddo
c
        if (ntmbb .gt. 1) then
          write (stdout,*) 
     &  ' Basin arrangement for Meridional Tracer Budget diagnostic'
          call iplot (msktmb, imt, imt, jmt)
        endif
c
c       write out the meridional tracer budget basin mask 
c
        if ((iotmb .ne. stdout .or. iotmb .lt. 0) .and. itmb) then
	  call getunit (iu, 'tracer_bud.dta','u s a ieee')
          iotext =
     &  ' read (iotmb) imt, jmt, ((msktmb(i,j),i=1,imt),j=1,jmt)'
          write (iu) stamp, iotext, expnam
          write (iu) imt, jmt, ((msktmb(i,j),i=1,imt),j=1,jmt)
          call relunit (iu)
        endif
      endif
c
      numtmb = 0
      do mask=0,ntmbb
        do jrow=1,jmt
          smdvol(jrow,mask)  = c0
        enddo
      enddo
      do mask=0,ntmbb
        do n=1,nt
          do jrow=1,jmt
            tstor(jrow,n,mask) = c0
            tdiv(jrow,n,mask)  = c0
            tflux(jrow,n,mask) = c0
            tdif(jrow,n,mask)  = c0
            tsorc(jrow,n,mask) = c0
          enddo
        enddo
      enddo
#endif
c
#if defined bryan_lewis_vertical || defined bryan_lewis_horizontal
c
c-----------------------------------------------------------------------
c     initialize Bryan_Lewis tracer diffusion coefficients
c-----------------------------------------------------------------------
c
      call blmixi
#endif
c
#ifdef trajectories
c
c-----------------------------------------------------------------------
c     initialize particle trajectories
c-----------------------------------------------------------------------
c
      if (initpt) call ptraji
        
#endif
#ifdef time_averages
c
c-----------------------------------------------------------------------
c     initialize time mean "averaging" grid data
c-----------------------------------------------------------------------
c
      call avgi
c
#endif
#ifdef xbts
c
c-----------------------------------------------------------------------
c     initialize XBT locations and averaging arrays 
c-----------------------------------------------------------------------
c
      call xbti
c
#endif
#ifdef ppvmix
c
c-----------------------------------------------------------------------
c     initialize pacanowski-philander vertical mixing scheme
c-----------------------------------------------------------------------
c
      call ppmixi (error, cifdef, ifdmax, nifdef, vmixset)
#endif
#ifdef smagnlmix
c
c-----------------------------------------------------------------------
c     initialize Smagorinsky nonlinear horizontal mixing scheme
c-----------------------------------------------------------------------
c
      call smagnli (error, cifdef, ifdmax, nifdef, hmixset)
#endif
#ifdef isopycmix
c
c-----------------------------------------------------------------------
c    initialize the isopycnal mixing
c-----------------------------------------------------------------------
c
      call isopi (error, cifdef, ifdmax, nifdef, am, ah)
c
# ifdef held_larichev
c
c-----------------------------------------------------------------------
c     initialize held_larichev horizontal mixing coefficients to
c     replace constant isopycnal mixing coefficient
c-----------------------------------------------------------------------
c
      call hlmixi (error, cifdef, ifdmax, nifdef)
# endif
#endif
c
c-----------------------------------------------------------------------
c     do all consistency checks last
c-----------------------------------------------------------------------
c
      call checks (error, vmixset, hmixset)
c
c-----------------------------------------------------------------------
c     list options enabled for this model run
c-----------------------------------------------------------------------
c
      call options (cifdef, ifdmax, nifdef)
c
#ifdef timing
      call toc ('driver', 'setocn')
#endif
      return
c
9341  format (//,'  Global ocean statistics:'
     &,/,'  the total ocean surface area (t cells) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (t cells)       =',1pe15.8,'cm**3'
     &,/,'  the total ocean surface area (u cells) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (u cells)       =',1pe15.8,'cm**3')
9401  format(/t50,'number of levels on "u,v" grid')
9402  format(/t50,'number of levels on "t,s" grid')
9601  format(/,' "dtxcel(km)" tracer timestep multipliers:',/,10(f9.3))
      end





      subroutine depth_u (kmt, imt, jmt, zw, km, kmu, h, hr)
c
c=======================================================================
c     calculate depth arrays associated with "u" cells.
c
c     input:
c       kmt = number of oecan "t" cells from surface to bottom of ocean
c       imt = longitudinal dimension of arrays
c       jmt = latitudinal dimension of arrays
c       zw  = depth to bottom of "t" cells
c       km  = max number of depths
c
c     output:
c       kmu = number of ocean "u" cells from surface to bottom of ocean
c       h   = depth to ocean floor over "u" cells (cm)
c       hr  = reciprocal of "h"       
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      dimension kmt(imt,jmt), kmu(imt,jmt), h(imt,jmt), hr(imt,jmt)
      dimension zw(km)
c
c-----------------------------------------------------------------------
c     set some constants
c-----------------------------------------------------------------------
c
      c0 = 0.0
      c1 = 1.0
c
c-----------------------------------------------------------------------
c     compute number of vertical levels on the "u" grid
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        kmu(imt,jrow) = 0
      enddo
c
      do i=1,imt
        kmu(i,jmt) = 0
      enddo
c
      do jrow=1,jmt-1
        do i=1,imt-1
           kmu(i,jrow) = min (kmt(i,jrow), kmt(i+1,jrow), kmt(i,jrow+1)
     &                      , kmt(i+1,jrow+1))
        enddo
      enddo
#ifdef cyclic
      do jrow=1,jmt
        kmu(imt,jrow) = kmu(2,jrow)
      enddo
#endif
#ifdef symmetry
      do i=1,imt
        kmu(i,jmt) = kmu(i,jmt-2)
      enddo
#endif
c
c---------------------------------------------------------------------
c     compute depths and reciprocal depths over "u" cells
c---------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
          hr(i,jrow) = c0
          h(i,jrow)  = c0
          if (kmu(i,jrow) .ne. 0) then
            hr(i,jrow) = c1/zw(kmu(i,jrow))
            h (i,jrow) = zw(kmu(i,jrow))
          endif
        enddo
      enddo
c
      return
      end




      subroutine options (cifdef, ifdmax, nifdef)
      character*(*) cifdef(ifdmax)
#include "stdunits.h"
c
c-----------------------------------------------------------------------
c     list "ifdef" options that are enabled (turned on)
c-----------------------------------------------------------------------
c
      logical error
c
      error = .false.
c
        write (stdout,'(/,10x,a,/)')
     &  ' For available options ---- consult the manual'
c
        write (stdout,'(/,a)')
     &   ' THE FOLLOWING "IFDEF" OPTIONS HAVE BEEN ENABLED (TURNED ON):'
#ifdef coupled
        nifdef = nifdef + 1
        cifdef(nifdef) = 'coupled     '
#endif
#ifdef simple_sbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'simple_sbc'
#endif
#ifdef time_mean_sbc_data
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_mean_sbc_data '
#endif
#ifdef time_varying_sbc_data
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_varying_sbc_data '
#endif
#ifdef minimize_sbc_memory
        nifdef = nifdef + 1
        cifdef(nifdef) = 'minimize_sbc_memory '
#endif
#ifdef trace_coupled_fluxes
        nifdef = nifdef + 1
        cifdef(nifdef) = 'trace_coupled_fluxes'
#endif
#ifdef linearized_advection
        nifdef = nifdef + 1
        cifdef(nifdef) = 'linearized_advection'
#endif
#ifdef linearized_density
        nifdef = nifdef + 1
        cifdef(nifdef) = 'linearized_density'
#endif
#ifdef timing
        nifdef = nifdef + 1
        cifdef(nifdef) = 'timing      '
#endif
#ifdef free_slip
        nifdef = nifdef + 1
        cifdef(nifdef) = 'free_slip'
#endif
#ifdef f_plane
        nifdef = nifdef + 1
        cifdef(nifdef) = 'f_plane'
#endif
#ifdef beta_plane
        nifdef = nifdef + 1
        cifdef(nifdef) = 'beta_plane'
#endif
#ifdef rectangular_box
        nifdef = nifdef + 1
        cifdef(nifdef) = 'rectangular_box'
#endif
#ifdef flat_bottom
        nifdef = nifdef + 1
        cifdef(nifdef) = 'flat_bottom '
#endif
#ifdef scripps_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'scripps_kmt '
#endif
#ifdef read_my_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'read_my_kmt '
#endif
#ifdef idealized_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'idealized_kmt '
#endif
#ifdef generate_a_grid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'generate_a_grid'
#endif
#ifdef read_my_grid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'read_my_grid'
#endif
#ifdef equivalence_mw
        nifdef = nifdef + 1
        cifdef(nifdef) = 'equivalence_mw'
#endif
#ifdef cray_ymp
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cray_ymp'
#endif
#ifdef sgi
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sgi'
#endif
#ifdef solid_walls
        nifdef = nifdef + 1
        cifdef(nifdef) = 'solid_walls'
#endif
#ifdef cyclic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cyclic      '
#endif
#ifdef symmetry
        nifdef = nifdef + 1
        cifdef(nifdef) = 'symmetry    '
#endif
#ifdef source_term
        nifdef = nifdef + 1
        cifdef(nifdef) = 'source_term '
#endif
#ifdef firfil
        nifdef = nifdef + 1
        cifdef(nifdef) = 'firfil      '
#endif
#ifdef fourfil
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fourfil     '
#endif
#if defined implicit_free_surface
        nifdef = nifdef + 1
        cifdef(nifdef) = 'implicit_free_surface'
#endif
#if defined rigid_lid_surface_pressure
        nifdef = nifdef + 1
        cifdef(nifdef) = 'rigid_lid_surface_pressure'
#endif
#ifdef stream_function
        nifdef = nifdef + 1
        cifdef(nifdef) = 'stream_function'
#endif
#ifdef sf_5_point
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sf_5_point'
#endif
#ifdef sf_9_point
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sf_9_point'
#endif
#ifdef conjugate_gradient
        nifdef = nifdef + 1
        cifdef(nifdef) = 'conjugate_gradient'
#endif
#ifdef oldrelax
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldrelax    '
#endif
#ifdef hypergrid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'hypergrid   '
#endif
#ifdef implicitvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'implicitvmix'
#endif
#ifdef fullconvect
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fullconvect'
#endif
#ifdef isopycmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'isopycmix   '
#endif
#ifdef held_larichev
        nifdef = nifdef + 1
        cifdef(nifdef) = 'held_larichev'
#endif
#ifdef gent_mcwilliams
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gent_mcwilliams'
#endif
#ifdef save_conection
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_convection'
#endif
#ifdef save_mixing_coeff
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_mixing_coeff'
#endif
#ifdef isopycmixspatialvar
        nifdef = nifdef + 1
        cifdef(nifdef) = 'isopycmixspatialvar'
#endif
#ifdef consthmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'consthmix   '
#endif
#ifdef biharmonic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'biharmonic  '
#endif
#ifdef shortwave
        nifdef = nifdef + 1
        cifdef(nifdef) = 'shortwave   '
#endif
#ifdef tcvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tcvmix      '
        write (stdout,*) ' Error ==> tcvmix not implemented yet'
	error = .true.	
#endif
#ifdef leq
        nifdef = nifdef + 1
        cifdef(nifdef) = 'leq         '
        write (stdout,*) ' Error ==> leq not fully checked out yet'
	error = .true.	
#endif
#ifdef lalg
        nifdef = nifdef + 1
        cifdef(nifdef) = 'lalg        '
        write (stdout,*) ' Error ==> lalg not fully checked out yet'
	error = .true.	
#endif
#ifdef readrmsk
        nifdef = nifdef + 1
        cifdef(nifdef) = 'readrmsk  '
#endif
#ifdef restorst
        nifdef = nifdef + 1
        cifdef(nifdef) = 'restorst    '
#endif
#ifdef sponges
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sponges    '
#endif
#ifdef bryan_lewis_horizontal
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bryan_lewis_horizontal'
#endif
#ifdef bryan_lewis_vertical
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bryan_lewis_vertical'
#endif
#ifdef constvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'constvmix   '
#endif
#ifdef snapshots
        nifdef = nifdef + 1
        cifdef(nifdef) = 'snapshots   '
#endif
#ifdef trajectories
        nifdef = nifdef + 1
        cifdef(nifdef) = 'trajectories'
#endif
#ifdef lyapunov
        nifdef = nifdef + 1
        cifdef(nifdef) = 'lyapunov    '
#endif
#ifdef time_averages
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_averages   '
#endif
#ifdef time_averages_disk
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_averages_disk'
#endif
#ifdef meridional_overturning
        nifdef = nifdef + 1
        cifdef(nifdef) = 'meridional_overturning'
#endif
#ifdef gyre_components
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gyre_components   '
#endif
#ifdef tracer_averages
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_averages   '
#endif
#ifdef show_zonal_mean_of_sbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'show_zonal_mean_of_sbc'
#endif
#ifdef show_external_mode
        nifdef = nifdef + 1
        cifdef(nifdef) = 'show_external_mode'
#endif
#ifdef diagnostic_surf_height
        nifdef = nifdef + 1
        cifdef(nifdef) = 'diagnostic_surf_height'
#endif
#ifdef term_balances
        nifdef = nifdef + 1
        cifdef(nifdef) = 'term_balances'
#endif
#ifdef meridional_tracer_budget
        nifdef = nifdef + 1
        cifdef(nifdef) = 'meridional_tracer_budget'
#endif
#ifdef time_step_monitor
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_step_monitor   '
#endif
#ifdef energy_analysis
        nifdef = nifdef + 1
        cifdef(nifdef) = 'energy_analysis'
#endif
#ifdef stability_tests
        nifdef = nifdef + 1
        cifdef(nifdef) = 'stability_tests'
#endif
#ifdef matrix_sections
        nifdef = nifdef + 1
        cifdef(nifdef) = 'matrix_sections'
#endif
#ifdef xbts
        nifdef = nifdef + 1
        cifdef(nifdef) = 'xbts        '
#endif
#ifdef levitus_ic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'levitus_ic   '
#endif
#ifdef idealized_ic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'idealized_ic'
#endif
#ifdef ramdrive
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ramdrive    '
#endif
#ifdef  crayio
        nifdef = nifdef + 1
        cifdef(nifdef) = 'crayio      '
#endif
#ifdef fio
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fio         '
#endif
#ifdef show_details
        nifdef = nifdef + 1
        cifdef(nifdef) = 'show_details'
#endif
#ifdef netcdf
        nifdef = nifdef + 1
        cifdef(nifdef) = 'netcdf'
#endif
c
#ifdef coupled
        write (stdout,'(/a/a/)')
     &  ' Note ==> MOM S.B.C. will come from an atmospheric model'
     &, '          with spatial interpolations for all S.B.C.'
#else
# if !defined simple_sbc
        write (stdout,'(/a/a/)')
     & ' Note ==> MOM S.B.C. will come from a prepared dataset'
     &,'          (requiring no spatial interpolations)'
# else
        write (stdout,'(/a,a/)')
     & ' Note ==> MOM S.B.C. are internally generated idealized zonal'
     &,' means.'
# endif
#endif
      if (error) then
        write (stdout,'(a)')
     & ' Regrets: sorry but the above options are not implemented yet.'
	stop '>options'
      endif
      nacros = 4
      nn = nifdef/float(nacros) + 1
      do i=1,nn
        ns = 1+nacros*(i-1)
        ne = ns + nacros-1
        if (ne .gt. ifdmax) then
          ne = ifdmax
          write (stdout,'(/a/)')
     &'=>Error: option list is truncated. increase parameter "ifdmax"'
          stop 'options'
	endif
        if (ns .le. nifdef) write (stdout,9101) (cifdef(n),n=ns,ne)
      enddo
      write (stdout,'(/a/)')  ' END OF ENABLED OPTIONS LIST'
c
      return
9101  format(4(1x,a30))
      end


      subroutine rowi
c
c=======================================================================
c
c     initialize prognositc quantities at "tau-1" and "tau"
c
c     inputs:
c
c     kmt  = number of vertical levels on "t" cells
c     yt   = latitudes of "t" points
c     zt   = depths of "t" points  
c
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "iounit.h"
#include "levind.h"
#include "mw.h"
#ifdef levitus_ic
      character*32 stampl
      dimension xtlev(imt), ztlev(km)
      dimension tbuf(imt,km), sbuf(imt,km)
#endif
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data on disk.
c     for latitude rows they point to latdisk(1) or latdisk(2)
c     for 2D fields they point to records on kflds
c-----------------------------------------------------------------------
c
      itt   = 0
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data in the MW
c-----------------------------------------------------------------------
c
      if (wide_open_mw) then
c
c       rotate time levels instead of moving data
c
        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
      endif
c
      ucksum = 0.0
      vcksum = 0.0
      tcksum = 0.0
      scksum = 0.0
c
c-----------------------------------------------------------------------
c     initialize every latitude jrow either in the MW (when wide opened)
c     or on disk (when jmw < jmt)
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
c
        if (wide_open_mw) then
	  j = jrow
	else
	  j = jmw
	endif
c
c-----------------------------------------------------------------------
c       zero out all variables. velocities are internal modes only
c-----------------------------------------------------------------------
c
        do k=1,km
          do i=1,imt
            u(i,k,j,1,taup1) = c0
            u(i,k,j,2,taup1) = c0
          enddo
	enddo
	do n=1,nt
	  do k=1,km
	    do i=1,imt
	      t(i,k,j,n,taup1) = c0
	    enddo
	  enddo
	enddo
c
c-----------------------------------------------------------------------
c       set the Temperature and Salinity
c-----------------------------------------------------------------------
c
#ifdef levitus_ic
c
c       set T & S using climatological monthly Levitus data prepared
c       (by user) from the MOM database using subroutine "ic.F"
c       note: initial condition time in "tmngr.F" should be set to the
c       middle of the starting month where this data is defined. 
c
        if (jrow .eq. 1) then
          call getunit (iolevt, 'temp.mom.ic', 'u s r ieee')
          call getunit (iolevs, 'salt.mom.ic', 'u s r ieee')
        endif
c
        read (iolevt)
        read (iolevt) stampl, dpml, im, kk, jro, alat, mm, xtlev, ztlev
     &,               tbuf
        read (iolevs)
        read (iolevs) stampl, dpml, im, kk, jro, alat, mm, xtlev, ztlev
     &,               sbuf
c
        if (im .ne. imt .or. kk .ne. km .or. jro .ne. jrow) then
          write (stdout,'(/a,a,/,6(a,i5)/)')
     &     ' ==> Error reading levitus data for ', stampl
     &,    ' jrow=',jrow, ', jro=',jro, ', imt=', imt, ', im=',im
     &,    ', km=',km, ', kk=',kk
          stop '=>rowi'
        endif
c
c       transfer the data
c
        do k=1,km
	  do i=1,imt
	    t(i,k,j,1,taup1) = tbuf(i,k)
	    t(i,k,j,2,taup1) = sbuf(i,k)
	  enddo
	enddo
c
        if (jrow .eq. jmt) then
          call relunit (iolevt)
	  call relunit (iolevs)
	  write (stdout,'(/a,a,/a,g14.7,a/)')
     &    '=> Initialized T & S from Levitus data for ', stampl
     &,   '   averaged over previous ', dpml, ' days'
        endif
#endif
#ifdef idealized_ic
c
c       construct highly idealized initial density profile
c
        do i=1,imt
          do k=1,km
            t(i,k,j,1,taup1) = theta0 (yt(jrow), zt(k))
            t(i,k,j,2,taup1) = 0.0349 - 0.035
          enddo
        enddo
        if (jrow .eq. jmt) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S to internally generated idealized'
     &,                    ' values for a test case.'
        endif
#endif
#ifdef linearized_advection
c
c       initialize density profile with tbarz
c
        do i=1,imt
          do k=1,km
            t(i,k,j,1,taup1) = tbarz(k)
            t(i,k,j,2,taup1) = 0.0349 - 0.035
          enddo
        enddo
        if (jrow .eq. jmt) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S for linearized advection'
        endif
#endif
c
c-----------------------------------------------------------------------
c       zero out Temperature and Salinity in land points
c-----------------------------------------------------------------------
c
        do i=1,imt
          kz = kmt(i,jrow)
          do k=1,km
            if (k .gt. kz) then
              t(i,k,j,1,taup1) = c0
              t(i,k,j,2,taup1) = c0
            endif
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       checksum the initial conditions
c-----------------------------------------------------------------------
c
        ucksum = ucksum + checksum (u(1,1,j,1,taup1), imt, km)
        vcksum = vcksum + checksum (u(1,1,j,2,taup1), imt, km)
        tcksum = tcksum + checksum (t(1,1,j,1,taup1), imt, km)
        scksum = scksum + checksum (t(1,1,j,2,taup1), imt, km)
c
c-----------------------------------------------------------------------
c       initialize every latitude jrow either on disk (when jmw < jmt)
c       or in the MW (when the last jrow is complete and jmw = jmt)
c-----------------------------------------------------------------------
c
        if (wide_open_mw) then
	  if (jrow .eq. jmt) then
	    call copy_all_rows (taup1, tau)
	    call copy_all_rows (tau, taum1)
	  endif
	else
          call putrow (latdisk(taudisk),  nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
          call putrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
	endif
c
      enddo
      write (stdout,*) ' I.C. checksum for t =',tcksum
      write (stdout,*) ' I.C. checksum for s =',scksum
      write (stdout,*) ' I.C. checksum for u =',ucksum
      write (stdout,*) ' I.C. checksum for v =',vcksum
      return
      end


      subroutine sethr (nr, xstart, xend, ystart, yend)
c
c=======================================================================
c     discretizes the horizontal region to nearset model grid points
c
c     nr     = the horizontal region number
c     xstart = starting longitude at edge of "t" box region
c     xend   = ending longitude at edge of  "t" box region
c     ystart = starting latitude at edge of  "t" box region
c     yend   = ending latitude at edge of "t" box region
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "cregin.h"
#include "levind.h"
c
c     find the nearest "t" box indicies within the region
c
      jsr = min(indp (ystart, yu, jmt)+1, jmt)
      jer = indp (yend, yu, jmt)
c
      isr = min(indp (xstart, xu, imt)+1, imt)
      ier = indp (xend, xu, imt)
c
c     define "edges" of "t" box region
c
      if (isr .eq. 1) then
        xsrl = xu(1) - dxudeg(1)
      else
        xsrl = xu(isr-1)
      endif
      xerl = xu(ier)
      if (jsr .eq. 1) then
        ysrl = yu(1) - dyudeg(1)
      else
        ysrl = yu(jsr-1)
      endif
      yerl = yu(jer)
c
      write (hregnm(nr),9000) xsrl, xerl, ysrl, yerl
      write (stdout,*) ' Defining horizontal region # ',nr
     &, ' as "t" cells within ', hregnm(nr)
      if (isr .gt. ier) then
        write (stdout,*) ' Error: isr=',isr,' >  ier=',ier,' in sethr'
        stop '=>sethr'
      endif
      if (jsr .gt. jer) then
        write (stdout,*) ' Error: jsr=',jsr,' >  jer=',jer, 'in sethr'
        stop '=>sethr'
      endif
      do j=jsr,jer
        do i=isr,ier
          if (kmt(i,j) .gt. 0)  then
            mskhr(i,j) = nr
          endif
        enddo
      enddo
9000  format ('lon: ',f5.1,' => ',f5.1,'  lat: ',f5.1,' => ',f5.1)
      return
      end


      subroutine setvr (nr, zstart, zend)
c
c=======================================================================
c     discretizes the vertical region to nearset model grid points
c
c     nr     = the vertical region number
c     zstart = starting depth at edge of "t" box region in cm.
c     zend   = ending depth at edge of "t" box region in cm.
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "cregin.h"
c
c     find the nearest "t" box indicies within the region
c
      if (zstart .lt. p5*zw(1)) then
        ksr = 1
        ztopb = 0.0
      else
        ksr = min(indp (zstart, zw, km)+1, km)
        ztopb = zw(ksr-1)*0.01
      endif
      ker = indp (zend, zw, km)
      llvreg(nr,1) = ksr
      llvreg(nr,2) = ker
c
      write (vregnm(nr),9000) ztopb, zw(ker)*0.01
      write (stdout,*) ' Defining vertical region # ',nr
     &,  ' as "t" cells within ',  vregnm(nr)
      if (ksr .gt. ker) then
        write (stdout,*) ' Error: ksr=',ksr,' >  ker=',ker, 'in setvr'
        stop '=>setvr'
      endif
9000  format (' dpt:',f6.1, '=>',f6.1, 'm')
      return
      end

      function theta0 (ydeg, depth)
c
c=======================================================================
c     this subroutine returns estimates of global mean potential
c     temperature for model initialization as a function of depth.
c     it is used to produce a reference thermal stratification for the
c     upper 2000m of the MOM`s test case.  below 2000m, the
c     potential temperature returned is 2.0 degrees C.  surface
c     values are set slightly above 18.4 degrees C at the reference
c     latitude "reflat".
c     the estimates are produced from a 7th order ploynomial fit to
c     the annual mean world ocean potential temperature observations
c     of Levitus (1982).
c
c     input [units]:
c       a latitdue (ydeg): [degrees]
c       a zt value (depth): [centimeters]
c     output [units]:
c       potential temperature estimate (est): [degrees centigrade]
c
c     variables:
c       coeft     = coefficients for the polynomial fit of potential
c                   temperature vs. depth
c       reflat    = reference latitude at which observed surface
c                   temperatures approximately equal coeft(1)
c       factor    = the ratio of the cosine of the latitude requested
c                   ("ydeg") to the reference latitude ("reflat")
c                   used to scale the upper 2000 meters of the vertical
c                   temperature profile
c       tmin,tmax = the minumum and maximum potential temperatures
c                   allowed at the time of model initialization
c
c     reference:
c       Levitus, S., Climatological atlas of the world ocean, NOAA
c     Prof. Paper 13, US Gov`t printing Office, Washington, DC, 1982.
c
      parameter (ndeg=7)
      dimension coeft(ndeg+1)
      save coeft, tmin, tmax, reflat
      data coeft / 0.184231944E+02,-0.430306621E-01, 0.607121504E-04
     &           ,-0.523806281E-07, 0.272989082E-10,-0.833224666E-14
     &           , 0.136974583E-17,-0.935923382E-22/
      data tmin, tmax, reflat /2.0, 25.0, 34.0/
c
c=======================================================================
c
      c0 = 0.0
      pi = atan(1.0) * 4.0
      refcos = abs(cos(pi*reflat/180.))
c
      coslat = abs(cos(pi*ydeg/180.))
      factor = coslat/refcos
      z = depth * 0.01
c
      if (z .gt. 2000.) then
        est = 2.0
      else
        est = c0
        bb = 1.0
        do nn=1,ndeg+1
c          if (nn.gt.1) bb = z**(nn-1)
          est = est + coeft(nn)*bb
	  bb = bb*z
        enddo
        est = est * factor
      endif
c
      if (est .gt. tmax) est = tmax
      if (est .lt. tmin) est = tmin
c
      theta0 = est
      return
      end
