      subroutine mom
c
c=======================================================================
c
c                        GFDL Modular Ocean Model
c
c
c             A primitive equation ocean model implemented by
c
c           R.C. Pacanowski   C. Goldberg   A. Rosati   K. Dixon
c
c                                at the
c
c              Geophysical Fluid Dynamics Laboratory /NOAA
c                         Princeton, NJ. 08542
c
c                     based on the pioneering work of
c
c         Kirk Bryan: A numerical method for the study of the  of the
c         circulation world ocean: 1969, J. Computat. Phys 4 347-376
c
c                                 and
c
c         the invaluable work of Mike Cox & Bert Semtner on earlier
c                       fortran implementations.
c
c
c        The GFDL Modular Ocean Model (acronym MOM) is a three
c        dimensional primitive equation ocean model intended  to be
c        a flexible tool useful for ocean and coupled air-sea modeling
c        applications over a wide range of space & time scales.
c        It is also intended to run efficiently on scalar, vector &
c        parallel architectures. The programming approach is
c        modular and additions to this model are encouraged to follow
c        this approach. Additional modules will be added with time and
c        new versions will be released as waranted. Please
c        refer to the read_me file for more details until a proper user
c        guide is available.
c
c        Requirements:
c
c        Standard fortran 77 is used (except for namelist, do enddo, and
c        variable names longer than 6 characters)
c        The preprocessor "cpp" (available on systems using "c" or UNIX)
c
c        Please refer questions or problems directly to the author(s)
c        who are responsible for their sections of code.
c
c
c        Disclaimer:
c
c        MOM is an ocean modeling research tool developed at GFDL.
c        Others may use it freely but we assume no responsibility
c        for problems or incorrect use of MOM. It is left to the user to
c        satisfy (him/her)self that a particular configuration is
c        working correctly. To this end, many of the included
c        diagnostics will be helpful.
c       
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      external blkdta
c
#include "param.h"
#include "emode.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
        save cps, itt1st
c
#ifdef timing
        call tic ('driver', 'ocean')
        call tic ('ocean', 'misc')
#endif
c
c-----------------------------------------------------------------------
c       integrate one time step
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c       update timestep counter, set corresponding model time, and set
c       time dependent logical switches which determine program flow.
c-----------------------------------------------------------------------
c
        itt = itt + 1
        call tmngr (dtts)
c
c-----------------------------------------------------------------------
c       read system clock to calculate compute_time/grid_point/timestep 
c       (neglecting overhead of "setocn")
c-----------------------------------------------------------------------
c
        if (first) then
	  itt1st = itt
          cps    = cpusec()
        endif
c
c-----------------------------------------------------------------------
c       initialize diagnostic variables
c-----------------------------------------------------------------------
c
        call diagi
c
c-----------------------------------------------------------------------
c       adjust various quantities for leapfrog/mixing timesteps
c
c       leapfrog----------> h(tau+1) = h(tau-1) + 2dt*F(tau)
c                          
c       forward-----------> tau-1 <= tau
c                           h(tau+1) = h(tau-1) + dt*F(tau)
c
c       euler backward:     tau-1 <= tau
c         euler1----------> h(tau` ) = h(tau-1) + dt*F(tau)
c         euler2----------> h(tau+1) = h(tau-1) + dt*F(tau')
c-----------------------------------------------------------------------
c
        if (leapfrog) then
c
c         normal leapfrog time step
c
          euler1  = .false.
	  euler2  = .false.
	  forward = .false.
	  eots    = .true.
c
          c2dtts  = c2*dtts
          c2dtuv  = c2*dtuv
          c2dtsf  = c2*dtsf
        else
c
c         mixing time step (forward step or euler backward step)
c
          if (eb) then
            euler1  = .true.
	    euler2  = .false.
	    forward = .false.
	    eots    = .false.
	  else
            euler1  = .false.
	    euler2  = .false.
	    forward = .true.
	    eots    = .true.
	  endif
c
          c2dtts = dtts
          c2dtuv = dtuv
          c2dtsf = dtsf
c
#ifdef stream_function
          do jrow=1,jmt
            do i=1,imt
              psi(i,jrow,2) = psi(i,jrow,1)
            enddo
	  enddo
#endif
        endif
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
        if (leapfrog) then
	  apgr = alph
        else
	  apgr = theta
          do jrow=1,jmt
            do i=1,imt
              ubarm1(i,jrow,1) = ubar(i,jrow,1)
	      ubarm1(i,jrow,2) = ubar(i,jrow,2)
            enddo
	  enddo
        endif
#endif
c
#ifdef implicit_free_surface
        if (leapfrog) then
          if (acor .ne. c0) acor = alph
        else
          if (acor .ne. c0) acor = theta
        endif
#endif
c
c-----------------------------------------------------------------------
c       set time centering "gcor" for coriolis term
c-----------------------------------------------------------------------
c
        if (acor .eq. c0) then
	  gcor = c1
        elseif (acor .ne. c0) then
	  gcor = c0
        endif
c
#if defined implicit_free_surface
c
c       set "gcor" for leapfrog or mixing time steps
c
        if (leapfrog) then
	  if (acor .eq. c0) then
	    gcor = c1
	  else
	    gcor = gam
	  endif
        else
          gcor = c0
        endif
#endif
c
c-----------------------------------------------------------------------
c       update pointers to tau-1, tau, & tau+1 data on disk.
c       for latitude rows they point to latdisk(1) or latdisk(2)
c       for 2D fields they point to records on kflds
c-----------------------------------------------------------------------
c
        taum1disk = mod(itt+1,2) + 1
        taudisk   = mod(itt  ,2) + 1
        taup1disk = taum1disk
c
c-----------------------------------------------------------------------
c       update pointers (indices) to tau-1, tau, & tau+1 data in the MW
c-----------------------------------------------------------------------
c
        if (wide_open_mw) then
c
c         rotate time levels instead of moving data
c
          taum1 = mod(itt+0,3) - 1
          tau   = mod(itt+1,3) - 1
          taup1 = mod(itt+2,3) - 1
        else
c
c         they are being held constant in time.
c
	endif
c
c=======================================================================
c
c               SOLVE THE BAROCLINIC AND TRACER EQUATIONS
c
c       Since all latitude rows may not fit into central memory, a
c       flexible MW (memory window) approach is used. The minimum MW
c       holds 3 latitude rows and the maximum MW holds "jmt" latitude
c       rows in central memory. Choose the size to fit into available
c       central memory. The MW is loaded with variables from disk as
c       many times as needed to solve latitude rows 2 through "jmt-2".
c
c       Example using a MW with 3 rows (jmw=3)
c
c       "loadmw" loads variables from the first 3 latitude jrows into
c       rows js=1 through je=3 in the 1st MW (mw=1). Equations are
c       computed for j=2 in the MW (corresponding to latitude jrow=2)
c       then written to disk. For the second MW (mw=2), "loadmw"
c       first moves variables from j=2 to j=1, then moves variables from
c       j=3 to j=2 in the MW, before loading latitude jrow=4 variables
c       into row js=je=3 in the MW. Equations are computed for j=2 in
c       the MW (corresponding to latitude jrow=3) then written to disk. The
c       process continues until latitude jrows 2 through jmt-1 are
c       computed.
c
c       Example using a MW with 5 rows (jmw=5)
c
c       "loadmw" loads variables from the first 5 latitude jrows into
c       rows js=1 through je=5 in the 1st MW (mw=1). Equations are
c       computed for j=2,3,4 in the MW (latitude jrows=2,3,4)
c       and written to disk. For the second MW (mw=2),
c       "loadmw" first moves variables from j=2 to j=1, then moves
c       variables from j=3 to j=2 in the MW, before loading jrow 6,7,8
c       variables into rows js=3 to je=5 in the MW. Equations
c       are computed for j=2,3,4 in the MW (latitude jrows=5,6,7)
c       then written to disk. The process continues until latitude
c       jrows 2 through jmt-1 are computed. Note that the last MW may
c       not always contain 3 computed rows (depending on "jmt").
c
c       Note: 
c
c       When the MW is fully opened (jmw=jmt), all latitude rows
c       reside in the MW (none on disk). Instead of reading/writing
c       from MW to disk, data is moved between time levels within the
c       MW.
c
c=======================================================================
c
1000    continue
#ifdef trace_indices
      if (leapfrog) write (stdout,'(/1x,a)')
     &   'Starting a --->Leapfrog timestep'
      if (forward) write (stdout,'(/1x,a)')
     &   'Starting a --->Forward timestep'
      if (euler1) write (stdout,'(/1x,a)')
     &   'Starting the first step of an ---> Euler backward timestep'
      if (euler2) write (stdout,'(/1x,a)')
     &   'Starting the second step of an --->Euler backward tstep'
      write (stdout,'(//1x,7(a,i4))')
     & "=> In mom: itt=",itt," disk indices: taum1disk =",taum1disk
     &," taudisk=",taudisk," taup1disk=",taup1disk
     &,", memory window indices: taum1=",taum1," tau=",tau
     &," taup1 ",taup1
#endif
c
c-----------------------------------------------------------------------
c       establish a queue for reading all "tau-1" and "tau" latitude
c       rows from disk into the MW. (only needed when i/o is slow to
c       overlap computation. otherwise this does nothing) 
c-----------------------------------------------------------------------
c
        if (wide_open_mw) then
c
c         no need to do anything
c
        else
          do jrow = 1,jmt
	    if (leapfrog .or. euler2) then
              call ofind (latdisk(taum1disk), nslab, jrow)
	    endif
            call ofind (latdisk(taudisk), nslab, jrow)
          enddo
	endif
c
c-----------------------------------------------------------------------
c       compute latitude rows within one MW, then move northward
c       and repeat until latitudes (2..jmt-1) are computed. 
c-----------------------------------------------------------------------
c
        do mw = 1,maxmw
c
c-----------------------------------------------------------------------
c         joff = offset relating row "j" in the MW to latitude "jrow"
c         js   = starting row within the MW for LOADING latitude row
c         je   = ending row within the MW for LOADING latitude row
c         is   = starting index longitude
c         ie   = ending index longitude
c
c         first MW     : load latitude rows into js=1      ... je=jmw
c         1 < MW < last: load latitude rows into js=jsmw+1 ... je=jmw
c         last MW      : load latitude rows into js=jsmw+1 ... je<=jmw
c-----------------------------------------------------------------------
c
	  joff = (mw-1)*ncrows
#ifdef biharmonic
	  js   = (jsmw+1)*(1-1/mw) + 1
#else
	  js   = jsmw*(1-1/mw) + 1
#endif
          je   = min(jmw,jmt-joff)
	  is   = 2
	  ie   = imt - 1
#ifdef trace_indices
      write (stdout,'(//1x,7(a,i4))')
     & "=> In mom: itt=",itt," mw =",mw," js=",js," je=",je
     &," joff=",joff," jrows=",js+joff," to ",je+joff
#endif
c
c-----------------------------------------------------------------------
c         load prognostic and related variables into the MW
c-----------------------------------------------------------------------
c
	  call loadmw (joff, js, je, is, ie, latdisk(taum1disk)
     &,                latdisk(taudisk))
c
c-----------------------------------------------------------------------
c         calculate advection velocities for momentum and tracers     
c-----------------------------------------------------------------------
c
          call adv_vel (joff, js, je, is, ie)
c
#ifdef isopycmix
c
c-----------------------------------------------------------------------
c         calculate isopycnal diffusion tensor components (and
c         gent_mcwilliams advective velocities) for use with tracers
c-----------------------------------------------------------------------
c
          call isopyc (joff, js, je, is, ie)
#endif
c
c-----------------------------------------------------------------------
c         set vertical mixing coefficients for momentum and tracers
c-----------------------------------------------------------------------
c
          call vmixc (joff, js, je, is, ie)
c
c-----------------------------------------------------------------------
c         set horizontal mixing coefficients for momentum and tracers
c-----------------------------------------------------------------------
c
          call hmixc (joff, js, je, is, ie)
c
c
c-----------------------------------------------------------------------
c         set which MW rows to calculate: jscalc through jecalc
c
c         1 <= MW < last: calculate rows jsmw...jemw
c                                   rows jsmw...jemw-1 for biharmonic
c         last MW       : calculate rows jsmw...je-1 (MW may not be full
c                                   rows jsmw...je-2 for biharmonic
c-----------------------------------------------------------------------
c
          jscalc = jsmw
#ifdef biharmonic
          jecalc = je - 2
	  if (mw .eq. maxmw) jecalc = je - 1
#else
	  jecalc = je-1
#endif
c
c-----------------------------------------------------------------------
c         set vertical boundary conditions for momentum and tracers
c-----------------------------------------------------------------------
c
          call setvbc (joff, jscalc, jecalc, is, ie)
c
c-----------------------------------------------------------------------
c         compute internal mode velocities and tracers
c-----------------------------------------------------------------------
c
          call clinic (joff, jscalc, jecalc, is, ie)
          call tracer (joff, jscalc, jecalc, is, ie)
c
c-----------------------------------------------------------------------
c         calculate diagnostics
c-----------------------------------------------------------------------
c
          call diag (joff, jscalc, jecalc, is, ie)
c
c-----------------------------------------------------------------------
c         write prognostic variables from the MW to disk "tau+1"
c-----------------------------------------------------------------------
c
          if (wide_open_mw) then
c
c           do nothing since variables are already in "tau+1" MW
c
	  else
            call putmw (joff, jscalc, jecalc, latdisk(taup1disk))
	  endif
c
        enddo
c
c=======================================================================
c
c               SOLVE THE BAROTROPIC EQUATION
c
c=======================================================================
c
        call tropic (c2dtsf, acor, cori(1,1), itt)
c
c-----------------------------------------------------------------------
c       if this is the 1st pass of an euler backward timestep, set the
c       disk pointers so the proper time levels are read on the 2nd pass 
c       and go back to do the 2nd pass.
c-----------------------------------------------------------------------
c
        if (euler1) then
          eots      = .true.
	  euler1    = .false.
	  euler2    = .true.
          ntaux     = taum1disk
          taum1disk = taudisk
          taudisk   = taup1disk
          taup1disk = ntaux
          go to 1000
        endif
	if (wide_open_mw .and. euler2) then
c
c         shuffle "tau" and "tau+1" after euler backward to
c         insure data is in the right place for the next timestep
c
	  call euler_shuffle
c
c         re-establish correct pointers for this timestep
c
          taum1 = mod(itt+0,3) - 1
          tau   = mod(itt+1,3) - 1
          taup1 = mod(itt+2,3) - 1
	endif
c
c-----------------------------------------------------------------------
c       output all remaining diagnostics
c-----------------------------------------------------------------------
c
        call diago
c
c-----------------------------------------------------------------------
c     if it`s the last timestep then clean things up otherwise return
c-----------------------------------------------------------------------
c
      if (eorun) then
c
c       calculate cp_time/grid_point/timestep (neglecting setocn)
c
        ncalls = itt - itt1st + 1
        cpe    = cpusec()
        write (stdout,'(/,1x, a, f9.3, a,1x,a,1pg14.7,a)')
     &    'Timing Summary: ', cpe-cps,  ' cp seconds,'
     &,    '(',(cpe-cps)/(ncalls*imt*jmt*km)
     &,    ' cp seconds per grid point per time step)'
c
c       write a restart
c
        if (restrt) call write_restart
c
        write (stdout,'(1x,a)') 'MOMdone'
        call oclose (kflds)
	call oclose (latdisk(1))
	call oclose (latdisk(2))
      endif
c
#ifdef timing
      call toc ('ocean', 'misc')
      call toc ('driver', 'ocean')
#endif
c
      return
      end


      function cpusec ()
c
c     user supplied function to measure cpu time in seconds
c
      cpusec = second()
      return
      end
      
