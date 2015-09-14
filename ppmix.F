#ifdef test_ppmix
# define ppvmix
# include "grids.F"
# include "iomngr.F"
# include "size_check.F"
# include "state.F"
# include "tmngr.F"
# include "util.F"
# define driver_only
      program driver
c
c=======================================================================
c
c           PACANOWSKI-PHILANDER VERTICAL MIXING MODULE
c
c     To test the pacanowski-philander vertical mixing scheme in a 
c     simple one dimensional model (at one latitude and longitude
c     within the MOM grid): 
c     
c        1) setup the grid. (see grids.F)
c
c        2 if the number of vertical levels is changed... run_denscoef       
c
c        3) compile and run this module using the script "run_ppmix"
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      logical error, vmixset, mixing
      parameter (ifdmax=100)
      character*30 cifdef(ifdmax)
# include "param.h"
# include "accel.h"
# include "calendar.h"
# include "coord.h"
# include "grdvar.h"
# include "iounit.h"
# include "mw.h"
# include "scalar.h"
# include "state.h"
# include "switch.h"
# include "tmngr.h"
# include "vmixc.h"
      dimension tt(km,2)
# include "dncoef.h"
# include "fdifm.h"
# include "fdift.h"
c
      print '(//,25x,a/)',
     &' T E S T I N G   A   1-D   V E R S I O N   O F   P P M I X'
c
c     initialize physical constants
c
      radius   = 6370.0e5
      grav     = 980.6
      rho0     = 1.035
c
c     initialize time step accelerators
c
      do k=1,km
        dtxcel(k) = 1.0
      enddo
c
c-----------------------------------------------------------------------
c     set up the grids in x (longitude), y (latitude), and z (depth)
c     corresponding to Arakawa "b" gird system
c-----------------------------------------------------------------------
c
      call grids
c
c-----------------------------------------------------------------------
c     initialize clock and calendar
c-----------------------------------------------------------------------
c
      year0  = 1994
      month0 = 5
      day0   = 1
c
      hour0  = 0
      min0   = 0
      sec0   = 0
      eqyear = .false.
      eqmon  = .false.
      monlen = 30
      refrun  = .false.
      refinit = .true.
      refuser = .false.
      if (refuser) then
        ryear  = 1900
        rmonth = 1
        rday   = 1
        rhour  = 0
        rmin   = 0
        rsec   = 0
      end if
      runlen   = 120.0
      rununits = 'days'
      idayrestart  = 0
      msrestart = 0
      dtts    = 3600.0
      dtuv    = 3600.0
      segtim  = dtts/86400.0
c
c     note: nmix is set but mixing timesteps are ignored if we
c     do a robert time filter every time step
c
      nmix    = 11
c#define robert_time_filter
#ifdef robert_time_filter
        write (stdout,*) ' Note: robert time filter is applied every ts'
#else
        write (stdout,*) ' Note: a forward ts is taken every nmix ts'
#endif
      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,            ryear, rmonth, rday, rhour, rmin, rsec
     &,            idayrestart, msrestart
     &,            runlen, rununits, rundays, dtts)
c
c-----------------------------------------------------------------------
c     prescribe some initial stratification
c-----------------------------------------------------------------------
c
      z0 = 30.0e2
      hh = 80.0e2
      zm = zt(km)
      t0 = 7.5
      t1 = 10.0
      print '(/,10x, a/)', 'Initial conditions'
      do k=1,km
        tt(k,1) = t0*(1.0 - tanh((zt(k)-hh)/z0)) + t1*(1.0-zt(k)/zm)
	tt(k,2) = 0.0349 - 0.035
	print *, 'k=',k,' zt(k)=',zt(k), ' temp=',tt(k,1)
     &, ' salt=',tt(k,2)
      enddo
c
c     set time levels
c
      taum1 = -1
      tau   =  0
      taup1 = +1
c
c     set I.C. for u,v,t,s  and set the land/sea masks to all sea (1.0)
c
      do j=1,jmw
        do i=1,imt
	  do k=1,km
	    umask(i,k,j) = 1.0
	    tmask(i,k,j) = 1.0
	    do n=1,2
	      u(i,k,j,n,tau)   = 0.0
	      u(i,k,j,n,taum1) = 0.0
	      u(i,k,j,n,taup1) = 0.0
	      t(i,k,j,n,tau)   = tt(k,n)
	      t(i,k,j,n,taum1) = tt(k,n)
	      t(i,k,j,n,taup1) = tt(k,n)
	    enddo
	  enddo
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     pick a latitude for the one dimensional model
c-----------------------------------------------------------------------
c
c      alat = 45.0
      alat = 0.0
      jrow = indp (alat, yu, jmt)
c
c-----------------------------------------------------------------------
c     set the coriolis factors for this latitude and choose
c     a centered coriolis term (explicit coriolis)
c-----------------------------------------------------------------------
c
      omega        = pi/43082.0
      cori(jrow,1) = c2*omega*sine(jrow)
      cori(jrow,2) = -c2*omega*sine(jrow)
      gcor         = 1.0
c
c-----------------------------------------------------------------------
c     initialize the vertical mixing scheme
c-----------------------------------------------------------------------
c
      error   = .false.
      nifdef  = 0
      vmixset = .false.
      call ioinit
c
      call ppmixi (error, cifdef, ifdmax, nifdef, vmixset)
      if (error) stop '=>driver'
c
c-----------------------------------------------------------------------
c     integrate equations for all k at one point (i,j)
c     set "joff=0" to prevent shifting MW northward
c     "j" is the row in the MW and "jrow" is the latitude row on disk
c      j=2 corresponds to the latitude of jrow
c-----------------------------------------------------------------------
c
      is   = imt/2
      ie   = imt/2
      joff = 0
c
      print '(/a,g14.7/)','  1-D model latitude = ', yu(jrow)
c
      do itt=1,100000
c
c       set all switches to control when things happens 
c
        call tmngr (dtts)
	mixing = .not. leapfrog
c
c       broadcast the central point to cover the memory window
c
        icent = is
        jcent = 2
        do n=1,2
	  do j=1,jmw
	    do i=is-1,ie+1
	      do k=1,km
	        u(i,k,j,n,tau) = u(icent,k,jcent,n,tau)
		u(i,k,j,n,taup1) = u(icent,k,jcent,n,taup1)
	        t(i,k,j,n,tau) = t(icent,k,jcent,n,tau)
		t(i,k,j,n,taup1)   = t(icent,k,jcent,n,taup1)
	      enddo
	    enddo
	  enddo
	enddo
#ifdef robert_time_filter
c
c       robert time filter for velocity
c
        smooth=.01
        do n=1,2
	  do j=1,jmw
	    do i=is-1,ie+1
	      do k=1,km
		u(i,k,j,n,tau)   = u(i,k,j,n,tau)
     &          + smooth*(0.5*(u(i,k,j,n,taup1) + u(i,k,j,n,taum1))
     &          - u(i,k,j,n,tau)) 
		t(i,k,j,n,tau)   = t(i,k,j,n,tau)
     &          + smooth*(0.5*(t(i,k,j,n,taup1) + t(i,k,j,n,taum1))
     &          - t(i,k,j,n,tau)) 
	      enddo
	    enddo
	  enddo
	enddo
#endif
c
c       move "tau"   variables ==> "tau-1" positions
c            "tau+1" variables ==> "tau"   positions
c
        do n=1,2
	  do j=1,jmw
	    do i=is-1,ie+1
	      do k=1,km
	        u(i,k,j,n,taum1) = u(i,k,j,n,tau)
		u(i,k,j,n,tau)   = u(i,k,j,n,taup1)
	        t(i,k,j,n,taum1) = t(i,k,j,n,tau)
		t(i,k,j,n,tau)   = t(i,k,j,n,taup1)
	      enddo
	    enddo
	  enddo
	enddo
c
#ifdef robert_time_filter
        dtu = 2.0*dtuv
	dtt = 2.0*dtts
#else
        if (mixing) then
c
c       move "tau"   variables ==> "tau-1" positions and set timestep
c
          do n=1,2
	    do j=1,jmw
	      do i=is-1,ie+1
	        do k=1,km
	          u(i,k,j,n,taum1) = u(i,k,j,n,tau)
	          t(i,k,j,n,taum1) = t(i,k,j,n,tau)
	        enddo
	      enddo
	    enddo
	  enddo
c
          dtu = dtuv
	  dtt = dtts
        else
          dtu = 2.0*dtuv
	  dtt = 2.0*dtts
        endif
#endif
c
c-----------------------------------------------------------------------
c       calculate mixing coefficients for i=is..ie and j=2..jmw-1
c       (variables are loaded for j=1,jmw. calculations from j=2,jmw-1)
c        output is:
c       "diff_cbt" = diffusive coeff at bottom of "t" cell
c       "visc_cbu" = diffusive coeff at bottom of "u" cell
c-----------------------------------------------------------------------
c
        call ppmix (joff, 1, jmw, is, ie)
c
c       prescribe some surface and bottom fluxes
c       wind components are cgs
c
	windx = -7.0e2
	windy = 3.0e2
	ro    = 1.0e-3
	cd    = 1.2e-3
	wind  = sqrt(windx**2 + windy**2)
	k     = 1
	airt  = tt(k,1) - 0.5
	damp  = 1.0/(30.0*86400.0)
        do j=jsmw,jemw
          do i=is,ie
            stf(i,j,1) = damp*(airt - t(i,k,j,1,taum1))/dzt(1)
	    stf(i,j,2) = 0.0
	    btf(i,j,1) = 0.0
	    btf(i,j,2) = 0.0
            smf(i,j,1) = ro*cd*wind*(windx - u(i,k,j,1,taum1))
	    smf(i,j,2) = ro*cd*wind*(windy - u(i,k,j,2,taum1))
	    bmf(i,j,1) = 0.
	    bmf(i,j,2) = 0.
	  enddo
        enddo
c
c-----------------------------------------------------------------------
c                   SOLVE THE 1-D TRACER EQUATIONS
c       construct diffusive flux for i=is..ie and and solve for j=2
c       "diff_fb" = diffusive flux at bottom of "t" cell 
c-----------------------------------------------------------------------
c
	do n=1,2
          do j=2,2
            do k=1,km-1
              do i=is,ie
                diff_fb(i,k,j) = diff_cbt(i,k,j)*(
     &                    t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))*dzwr(k)
              enddo
            enddo
            do i=is,ie
              diff_fb(i,0,j)  = stf(i,j,n)
              diff_fb(i,km,j) = btf(i,j,n)
            enddo
c
c           solve tracer eqns for each depth
c
	    i = is
	    do k=1,km
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + dtt*DIFF_Tz(i,k,j)
c
c             set lateral boundaries
c
              do jj=1,jmw
                do ii=is-1,ie+1
                  t(ii,k,jj,n,taup1) = t(i,k,j,n,taup1)
	        enddo
	      enddo
            enddo                   
          enddo                   
        enddo
c
c-----------------------------------------------------------------------
c                   SOLVE THE 1-D MOMENTUM EQUATIONS
c       construct diffusive flux for i=is..ie and and solve for j=2
c       "diff_fb" = diffusive flux at bottom of "u" cell 
c-----------------------------------------------------------------------
c
	do n=1,2
          do j=2,2
            do k=1,km-1
              do i=is,ie
                diff_fb(i,k,j) = visc_cbu(i,k,j)*
     &                          (u(i,k,j,n,taum1) - u(i,k+1,j,n,taum1))
     &                          *dzwr(k)
              enddo
            enddo
            do i=is,ie
              diff_fb(i,0,j)  = smf(i,j,n)
              diff_fb(i,km,j) = bmf(i,j,n)
            enddo
c
c           solve momentum eqns for each depth
c
	    i = is
	    do k=1,km
              u(i,k,j,n,taup1) = u(i,k,j,n,taum1) + dtu*( 
     &                           CORIOLIS(i,k,j,n) + DIFF_Uz(i,k,j))
c
c             set lateral boundaries
c
              do jj=1,jmw
                do ii=is-1,ie+1
                  u(ii,k,jj,n,taup1) = u(i,k,j,n,taup1)
	        enddo
	      enddo
            enddo                   
          enddo
	enddo
c
c       show some results
c
        if (eoday) then
          i = is
	  j = 2
          print '(/a,i6,a,a)',' End of day at itt=',itt, '  ',stamp
          print '(3(a,g11.4))'
     &,   ' taux=',smf(i,j,1), ' tauy=',smf(i,j,2), ' hflx=',stf(i,j,1)
          do k=1,6
	    print '(a,i2,8(a,g11.4))' 
     &,     'k=',k,' u=',u(i,k,j,1,taup1),' v=',u(i,k,j,2,taup1), ' t='
     &,     t(i,k,j,1,taup1), ' diff_cbt=',diff_cbt(i,k,j)
     &,     ' visc_cbu=', visc_cbu(i,k,j)
	  enddo
        endif
        if (eorun) stop
      enddo
	
      stop
      end
# ifdef timing
      subroutine tic (a, b)
      character*(*) a, b
      entry toc (a, b)
      return
      entry ticr (a, b)
      return
      end
# endif

#endif



#ifdef ppvmix
      subroutine ppmixi (error, cifdef, ifdmax, nifdef, vmixset)
c
      logical vmixset, error
      character*(*) cifdef(ifdmax)
c
# include "size.h"
# include "accel.h"
# include "coord.h"
# include "iounit.h"
# include "scalar.h"
# include "stdunits.h"
# include "vmixc.h"
c
c=======================================================================
c     Initialization for the Pacanowski/Philander vertical mixing scheme
c     Pacanowski & Philander (JPO vol 11, #11, 1981). 
c    
c     input:
c       dzt    = thickness of vertical levels (cm)
c       km     = number of vertical levels
c       yt     = latitude of grid points (deg)
c       jmt    = number of latitudes
c       dtxcel = time step accelerator as a function of level
c       dtts   = density time step (sec)
c       dtuv   = internal mode time step (sec)
c       error  = logical to signal problems
c       cifdef = array of character strings for listing enabled "ifdefs"
c       ifdmax = size of "cifdef"
c       nifdef = current number of enabled "ifdefs" 
c       vmixset= logical to determine if a vertical mixing scheme was
c                chosen
c
c     output:
c       wndmix = min value for mixing at surface to simulate high freq
c                wind mixing (if absent in forcing). (cm**2/sec)
c       fricmx = maximum mixing (cm**2/sec)
c       diff_cbt_back = background "diff_cbt" (cm**2/sec)
c       visc_cbu_back = background"visc_cbu" t(cm**2/sec)
c       diff_cbt_limit = largest "diff_cbt" (cm**2/sec)
c       visc_cbu_limit = largest "visc_cbu" (cm**2/sec)
c       cifdef = array of character strings for listing enabled "ifdefs"
c       nifdef = incremented by 1 if this routine is called 
c       error  = true if some inconsistancy was found
c       vmixset= true
c       
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
c
      namelist /ppmix/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
     &,                visc_cbu_limit, diff_cbt_limit
c
      write (stdout,'(/,20x,a,/)')
     & 'P P V M I X    I N I T I A L I Z A T I O N'
c
c-----------------------------------------------------------------------
c     initialize variables (all mixing units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      wndmix    = 10.0
      fricmx    = 50.0
      diff_cbt_back =  0.1
      visc_cbu_back =  1.0
      dzmin     = 1.e10
      p25       = 0.25
      c0        = 0.0
c
# ifdef implicitvmix
c
c     simulate convective adjustment with large mixing coefficient
c     limits
c
      visc_cbu_limit = 1.0e6
      diff_cbt_limit = 1.0e6
# else
c
c     in regions of gravitational instability set mixing limits to the
c     maximum consistant with the "cfl" criterion. convective adjustment
c     will also act on the instability.
c
      visc_cbu_limit = fricmx
      diff_cbt_limit = fricmx
# endif
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of above settings + documentation
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist', 'fsr')
      read (io,ppmix,end=100)
100   continue
c
c-----------------------------------------------------------------------
c     set no-flux condition on density difference across bottom level
c-----------------------------------------------------------------------
c
      do j=1,jmw
        do i=1,imt
	  rhom1z(i,km,j) = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     add character string to "ifdef option list" indicating that this
c     option is enabled
c-----------------------------------------------------------------------
c
      nifdef = nifdef + 1
      cifdef(nifdef) = 'ppvmix      '
c
c-----------------------------------------------------------------------
c     check for problems
c-----------------------------------------------------------------------
c
# if defined ppvmix && !defined implicitvmix && defined isopycmix
        write (stdout,'(/,(1x,a))')
     & '==> Error:  "ppvmix" must use "implicitvmix" when "isopycmix" '
     &,'            is also enabled. Also "aidif" should = 0.5        '
        error = .true.
# endif
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "ppvmix" cannot be enabled because another   '
     &,'            vertical mixing scheme has been enabled           '        
        error = .true.
      else
        vmixset = .true.
      endif
c
      do k=1,km
        dzmin = min(dzmin,dzt(k))
      enddo
      if (dzmin .ge. 25.e2) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "ppvmix" may not work well with coarse vertical  '
     &,'              resolution                                      '
      endif
c
      extlat = c0
      do jrow=1,jmt
        extlat = max(abs(yt(jrow)),extlat)
      enddo
      if (extlat .gt. 10.0) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "ppvmix" may not work well outside the tropics   '
     &,'              where vertical shear is small unless solar      '
     &,'              shortwave penetration into the ocean is         '
     &,'              accounted for by enabeling  "shortwave"         '
      endif
c
# if !defined implicitvmix
      do k=1,km
        if ((dtts*dtxcel(k)*fricmx)/dzt(k)**2 .ge. p25) then
          write (stdout,'(/,(1x,a))')
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'           "fricmx".  use a smaller "dtts", "dtxcel", and/or  '
     &,'           "fricmx" .... or enable "implicitvmix"             '
          write (stdout,'(a48,i3)') ' at level =',k
        endif
        if ((dtts*dtxcel(k)*diff_cbt_limit)/dzt(k)**2 .ge. p25) then
          write (stdout,'(/,(1x,a))')
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'            "diff_cbt_limit". use a smaller "dtts", "dtxcel"  '
     &,'         ,and/or  "diff_cbt_limit" ...or enable "implicitvmix"'
          write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
c
      if ((dtuv*fricmx)/dzmin**2 .ge. p25) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'            "fricmx". use a smaller "dtuv" and/or "fricmx"    '
     &,'             or enable "implicitvmix"                         '
      endif
c
      if ((dtuv*visc_cbu_limit)/dzmin**2 .ge. p25) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'            "visc_cbu_limit". use a smaller "dtuv" or         '
     &,'            "visc_cbu_limit" or enable "implicitvmix"         '
      endif
# else
        write (stdout,'(/,(1x,a))')
     & '==> Warning: enabeling "implicitvmix" with "ppvmix" uses      '
     &,'             variables defined at "tau" rather than at "tau-1"'
     &,'             as was done in MOM 1.x'
# endif
# ifdef bryan_lewis_vertical
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a))')
     & '==> Warning: "bryan_lewis_vertical" tracer diffus coefficients'
     &,'              will  be added to "ppvmix" diffus coefficients  '        
     &,'              Note that diff_cbt_back is being reset to zero  '
     &,'              while diff_cbu_back is unchanged                '
      diff_cbt_back =  0.0
# endif
c
c     write out namelist values
c
      write (stdout,ppmix)
      call relunit (io)
      call getunit (iodoc, 'document.dta', 'f s a')
      write (iodoc, ppmix)
      call relunit (iodoc)
      return
      end


      subroutine ppmix (joff, js, je, is, ie)
c
      integer tlev
# include "size.h"
# include "coord.h"
# include "grdvar.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "vmixc.h"
c
c=======================================================================
c     Compute vertical mixing coefficients based on...
c     Pacanowski & Philander (JPO vol 11, #11, 1981). 
c
c     Note: this parameterization was designed for equatorial models
c     and may not do a good job in mid or high latitudes. Simulations
c     in these regions (where vertical shear is small) are improved with
c     the addition of solar short wave penetration into the ocean which 
c     reduces buoyancy and enhances vertical mixing.
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c               "joff" > 0 moves variables
c      js     = starting row for loading variables to calculate
c               coefficients. calculations start at jstrt=max(js-1,jsmw)
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the 
c               longitude direction
c      ie     = ending index for calculating coefficients in the 
c               longitude direction
c      km     = number of vertical levels
c      grav   = gravity (cm/sec**2)
c      umask  = land/sea mask on "u" grid (land=0.0, sea=1.0)
c      tmask  = land/sea mask on "t" grid (land=0.0, sea=1.0)
c      fricmx = max viscosity (cm**2/sec)
c      wndmix = min viscosity at bottom of 1st level to simulate 
c               missing high frequency windstress components (cm**2/sec)
c      visc_cbu_back = background "visc_cbu" (cm**2/sec)
c      diff_cbt_back = background "diff_cbt" (cm**2/sec)
c      visc_cbu_limit = largest "visc_cbu" in regions of gravitational
c                      instability (cm**2/sec)
c      diff_cbt_limit = largest "diff_cbt" in regions of gravitational 
c                       instability (cm**2/sec)
c
c     outputs:
c
c      riu     = richardson number at bottom of "u" cells
c      rit     = richardson number at bottom of "t" cells
c      visc_cbu = viscosity coefficient at bottom of "u" cells (cm**2/s)
c      diff_cbt = diffusion coefficient at bottom of "t" cells (cm**2/s)
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      dimension ro(imt,km,1:jmw)
c
# ifdef timing
      call tic ('vmixc', 'ppmix')
# endif
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      c0    = 0.0
      c1    = 1.0
      c5    = 5.0
      p25   = 0.25
      epsln = 1.e-25
      fx    = -p25*grav
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     set time level 
c-----------------------------------------------------------------------
c
# ifdef implicitvmix
      tlev = tau
# else
      tlev = taum1
# endif
c
c-----------------------------------------------------------------------
c     set "ro" (density) at j=1 for 1st memory window otherwise ... move 
c     variables from top two rows to bottom two rows to eliminate
c     redundant calculation (three rows if using "bi-harmonic")
c-----------------------------------------------------------------------
c
      if (joff .eq. 0) then
        do k=1,km
	  do i=istrt-1,iend+1
	    ro(i,k,1) = c0
	  enddo
	enddo
      else
        call moveppmix (istrt-1, iend+1)
      endif
c
c-----------------------------------------------------------------------
c     compute density difference across bottom of "t" cells at tau-1
c     for rows js through je in the MW. Set density difference = zero
c     across bottom and in land areas 
c-----------------------------------------------------------------------
c
      do ks=1,2
        call statec (t(1,1,1,1,tlev), t(1,1,1,2,tlev), ro(1,1,jsmw)
     &,              max(js,jsmw), je, istrt-1, iend+1, ks)
        do j=js,je
          do k=ks,km-1,2
            do i=istrt-1,iend+1
              rhom1z(i,k,j) = (ro(i,k,j) - ro(i,k+1,j))*tmask(i,k+1,j)
            enddo
          enddo
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute richardson numbers on bottom of "u" cells
c     ("jsriu" is also appropriate for the bi-harmonic mixing case)
c     Pacanowski & Philander (JPO vol 11, #11, 1981). Eq: 3 
c-----------------------------------------------------------------------
c
      jsriu = max(js,jsmw)-1
      do j=jsriu,je-1
        do k=1,km-1
          t1 = fx*dzw(k)
          do i=istrt-1,iend
            riu(i,k,j) = t1*umask(i,k+1,j)*(
     &                   rhom1z(i,k,j+1) + rhom1z(i+1,k,j+1) +
     &                   rhom1z(i,k,j)   + rhom1z(i+1,k,j)) /
     &                   ((u(i,k,j,1,tlev) - u(i,k+1,j,1,tlev))**2 +
     &                    (u(i,k,j,2,tlev) - u(i,k+1,j,2,tlev))**2 + 
     &                    epsln)
          enddo
        enddo
	call setbcx (riu(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute richardson numbers on bottom of "t" cells as average
c     of four nearest richardson numbers on bottom of "u" cells
c     ("jstrt" is also appropriate for the bi-harmonic mixing case)
c     Pacanowski & Philander (JPO vol 11, #11, 1981). average of Eq 3 
c-----------------------------------------------------------------------
c
      jstrt = max(js-1,jsmw)
      do j=jstrt,je-1
        do k=1,km-1
          do i=istrt,iend
            rit(i,k,j) = p25*(riu(i,k,j) + riu(i-1,k,j)
     &                     + riu(i,k,j-1) + riu(i-1,k,j-1))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute vertical viscosity coeff on "u" cell bottoms
c     Pacanowski & Philander (JPO vol 11, #11, 1981). Eq: 1 
c-----------------------------------------------------------------------
c
      do j=jstrt,je-1
        do k=1,km-1
          do i=istrt,iend
            t1             = c1/(c1 + c5*riu(i,k,j))
            visc_cbu(i,k,j) = fricmx*t1**2 + visc_cbu_back
          enddo
        enddo
      enddo
c
c     in regions of gravitational instability, reset the vertical
c     mixing coefficient to its limit
c
      do j=jstrt,je-1
        do k=1,km-1
          do i=istrt,iend
            if (riu(i,k,j) .lt. c0) visc_cbu(i,k,j) = visc_cbu_limit
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute vertical diffusivity coeff on "t" cell bottoms
c     Pacanowski & Philander (JPO vol 11, #11, 1981). Eq: 2 
c-----------------------------------------------------------------------
c
      do j=jstrt,je-1
        do k=1,km-1
          do i=istrt,iend
            t1       = c1/(c1 + c5*rit(i,k,j))
            diff_cbt(i,k,j) = fricmx*t1**3 + diff_cbt_back
          enddo
        enddo
      enddo
c
c     in regions of gravitational instability, reset the vertical
c     mixing coefficient to the limit
c
      do j=jstrt,je-1
        do k=1,km-1
          do i=istrt,iend
            if (rit(i,k,j) .lt. c0) diff_cbt(i,k,j) = diff_cbt_limit
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     approximation for high freq wind mixing near the surface
c     set no flux through bottom of bottom level "km"
c-----------------------------------------------------------------------
c
      do j=jstrt,je-1
        do i=istrt,iend
          if (diff_cbt(i,1,j) .lt. wndmix) diff_cbt(i,1,j) = wndmix
          if (visc_cbu(i,1,j) .lt. wndmix) visc_cbu(i,1,j) = wndmix
	  diff_cbt(i,km,j) = c0
	  visc_cbu(i,km,j) = c0
        enddo
      enddo
c
#ifdef bryan_lewis_vertical
c
c-----------------------------------------------------------------------
c     add Bryan-Lewis mixing if wanted
c-----------------------------------------------------------------------
c
      do j=jstrt,je-1
        do k=1,km-1
          do i=istrt,iend
            diff_cbt(i,k,j) = diff_cbt(i,k,j) + Ahv(k)
          enddo
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c       set diffusion and viscosity coeffs to zero on land box bottoms
c-----------------------------------------------------------------------
c
      do j=jstrt,je-1
        do k=1,km-1
          do i=istrt,iend
            visc_cbu(i,k,j) = umask(i,k+1,j)*visc_cbu(i,k,j)
            diff_cbt(i,k,j) = tmask(i,k+1,j)*diff_cbt(i,k,j)
          enddo
        enddo
        call setbcx (visc_cbu(1,1,j), imt, km)
        call setbcx (diff_cbt(1,1,j), imt, km)
      enddo
c
# ifdef matrix_sections
      if (prxzts .and. eots) then
        call diagpp (joff, jstrt, je-1)
      endif
# endif
# ifdef trace_indices
      write (stdout,'(2x,7(a,i4))')
     & "=> In ppmix: js=",js," je=",je," joff=",joff
     &," jstrt=",jstrt," jsriu=",jsriu," jrow=",jstrt+joff
     &," to ",je-1+joff 
# endif
# ifdef timing
      call toc ('vmixc', 'ppmix')
# endif
c
      return
      end



# ifdef matrix_sections
      subroutine diagpp (joff, js, je)
#  include "param.h"
#  include "coord.h"
#  include "cprnts.h"
#  include "iounit.h"
#  include "switch.h"
#  include "tmngr.h"
#  include "vmixc.h"
c
# ifdef timing
      call tic ('diagnostic', 'matrix sections')
# endif
      do j=js,je
	jrow = j + joff
        reltim = relyr
        do jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt)) then
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2
c
            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cbt', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cbt(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
	      call getunit (io, 'sections.dta','u s a ieee')
              write (stdout,*) ' => diff_cbt ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cbt(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cbt(1,1,j), imt*km)
              call relunit (io)
            endif
c
            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_cbu', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_cbu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
	      call getunit (io, 'sections.dta','u s a ieee')
              write (stdout,*) ' => visc_cbu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_cbu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_cbu(1,1,j), imt*km)
              call relunit (io)
            endif
          endif
        enddo
      enddo
# ifdef timing
      call toc ('diagnostic', 'matrix sections')
# endif
      return
9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
      end
# endif

      
      subroutine moveppmix (is, ie)
c
c=======================================================================
c     as the MW moves northward, move data from the last two rows
c     into the first two rows. (last 3 rows if using biharmonic option)
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
# include "param.h"
# include "vmixc.h"
c
      nrows = jmw - ncrows
      do move=1,nrows
        jfrom = jmw - (nrows - move)
	jto   = move
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (1:jmw)
c-----------------------------------------------------------------------
c
        do k=1,km
	  do i=is,ie
	    rhom1z(i,k,jto)  = rhom1z(i,k,jfrom)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (1:jemw)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jemw) then
          do k=1,km
	    do i=is,ie
	      riu(i,k,jto)  = riu(i,k,jfrom)
            enddo
          enddo
        endif
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (jsmw:jmw)
c-----------------------------------------------------------------------
c
        if (jto .ge. jsmw) then
c         nothing to move
        endif
c
# ifdef biharmonic
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (jsme:jemw)
c-----------------------------------------------------------------------
c
        if (jto .eq. jsmw) then
          do k=1,km
	    do i=is,ie
	      rit(i,k,jto)  = rit(i,k,jfrom)
            enddo
          enddo
        endif
# endif
# ifdef trace_indices
          write (stdout,'(2x,2(a,i4))')
     & "=> In moveppmix: moving variables on row ",jfrom," to row ",jto
# endif
      enddo
      return
      end
#else
      subroutine ppmix
      return
      end
#endif
