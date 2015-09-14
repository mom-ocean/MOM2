      subroutine checks (errorc, vmixset, hmixset)
      logical errorc, vmixset, hmixset
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "hmixc.h"
#include "iounit.h"
#include "levind.h"
#ifdef isopycmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
#ifdef xbts
# include "xbt.h"
#endif
c
c-----------------------------------------------------------------------
c     do consistency checks before allowing MOM to continue
c-----------------------------------------------------------------------
c
      write (stdout,'(/,20x,a,/)')
     &         'G E N E R A L    C O N S I S T E N C Y    C H E C K S'
c
      if (imt .lt. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  parameter "imt" less than 3 is not allowed'        
        errorc = .true.
      endif
c
      if (jmt .lt. 4) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  parameter "jmt" less than 4 is not allowed'        
        errorc = .true.
      endif
c
#ifdef linearized_advection
# if defined levitus_ic || defined idealized_ic
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "levitus_ic" or "idealized_ic" cannot be          '
     &,'            be enabled when using "linearized_advection".     '        
        errorc = .true.
# endif
# if !defined linearized_density
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "linearized_density" must be enabled when using   '
     &,'            "linearized_advection".                           '        
        errorc = .true.
# endif
#else
# if !defined levitus_ic && !defined idealized_ic
      write (stdout,'(/,(1x,a))')
     & '==> Error:  either "levitus_ic" or "idealized_ic" must be     '
     &,'            enabled to specify the type of initial conditions '        
        errorc = .true.
# endif
# if defined levitus_ic && defined idealized_ic
      write (stdout,'(/,(1x,a))')
     & '==> Error:  either "levitus_ic" or "idealized_ic" must be     '
     &,'            enabled... but not both.                          '        
        errorc = .true.
# endif
#endif
c
#if defined minimize_sbc_memory && defined time_mean_sbc_data
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "minimize_sbc_memory" and "time_mean_sbc_data"    '
     &,'            are incompatible. "minimize_sbc_memory" is        '        
     &,'            intended for option "time_varying_sbc_data".      '        
        errorc = .true.
#endif
#ifdef consthmix
      if (hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "consthmix"  cannot be enabled because another     '
     &,'            horizontal mixing scheme has been enabled         '        
        errorc = .true.
      else
        hmixset = .true.
      endif
#endif
      if (.not.hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal mixing scheme has been enabled  '
        errorc = .true.
      endif
c
#ifdef constvmix
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "constvmix"  cannot be enabled because another     '
     &,'            vertical mixing scheme has been enabled           '        
        errorc = .true.
      else
c
c       set vmixset = true for enabeling "constvmix"
c
        vmixset = .true.
      endif
#endif
#ifdef bryan_lewis_vertical
      if (.not. vmixset) then
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a/1x,a/1x,a))')
     & '==> Error: "bryan_lewis_vertical" tracer diffusion requires   '
     &,'           "constvmix", "ppvmix" or "tcvmix" to be enabled.   '
     &,'           If "constvmix" is enabled, "bryan_lewis_vertical"  '
     &,'           coefficients replace kappa_h. If "ppvmix" or       '
     &,'           "tcvmix" is enabled, "bryan_lewis_vertical"        '
     &,'           coefficients are added to the predicted coeffs.    '        
        errorc = .true.
      endif
#endif
#if defined bryan_lewis_horizontal
# if !defined consthmix
        write (stdout,'(/,(1x,a/1x,a))')
     & '==> Error: "bryan_lewis_horizontal" tracer diffusion requires '
     &,'            "consthmix" to be enabled.                        '
        errorc = .true.
# endif
#endif
c
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No vertical mixing scheme has been enabled  '
        errorc = .true.
      endif
#ifdef tcvmix
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tcvmix" cannot be enabled because another         '
     &,'            vertical mixing scheme has been enabled           '        
        errorc = .true.
      else
        vmixset = .true.
      endif
#endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: there is no vertical mixing scheme enabled         '
        errorc = .true.
      endif
#if defined rigid_lid_surface_pressure
        write (stdout,9000)
     & '==> Warning: when using the "rigid_lid_surface_pressure" method'
     &,'             a checkerboard may develop in "ps". Using option '
     &,'             "remove_ps_checkerboard" may help to remove it.'
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
# if defined hypergrid || defined oldrelax
        write (stdout,9000)
     & '==> Error: "hypergrid" or "oldrelax" is not allowed when using'
     &,'           "implicit_free_surface" or '
     &,'            "rigid_lid_surface_pressure" options.'
        errorc = .true.  
# endif
# if defined sf_5_point
        write (stdout,9000)
     & '==> Warning:  "sf_5_point" is ignored when using the'
     &,'              "implicit_free_surface" or '
     &,'              "rigid_lid_surface_pressure" option.'
     &,'              conjugate_gradient with 9 point numerics is used.'
# endif
#endif
c
#if !defined stream_function && !defined implicit_free_surface
# if !defined rigid_lid_surface_pressure
        write (stdout,9000)
     & '==> Error: choose "stream_function", "implicit_free_surface"'
     &,'           or the "rigid_lid_surface_pressure"  method'
        errorc = .true.  
# endif
#endif
c
#ifdef stream_function
# if !defined sf_5_point && !defined sf_9_point
        write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using'
     &,'           the stream_function method'
        errorc = .true.  
# endif
# if defined sf_5_point && defined sf_9_point
        write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using'
     &,'           the stream_function method (not both)'
        errorc = .true.  
# endif
#endif
c
c     look for potential load imbalance 
c
      n = mod((jmt-2),ncrows)
      if (n .ne. 0) then
        write (stdout
     &,'(1x,a,i4,/14x,a,i4,a,i4/14x,a,i4,/14x,a,i4,a,/14x,a,i4,a)')
     &'==> Warning: the number of rows (jmw) in the MW (memory window)='
     &, jmw, 'with calculations starting on MW row =',jsmw
     &,' and ending on MW row =',jemw
#ifdef biharmonic
     &                               -1
#endif
     &,'disk latitude rows to be calculated are from 2 through', jmt-1
     &,'to compute these rows, the MW is moved northward '
     &, maxmw,' times.','however, the last time there are only ',n
     &,' row(s) to calculate.'
        write (stdout,'(/14x,a,/,14x,a,/,14x,a)')
     &  'If you are using multiprocessors, this will result in a load'
     &, 'imbalance on the last MW. The following is a list of "jmw"'
     &, 'which will not cause a load imbalance:'
#ifdef biharmonic
        jbuf = 3
        jjs  = 4
#else
        jbuf = 2
        jjs  = 3
#endif
        do jj=jjs,jmt-2
          ncrow = jj-jbuf
          if (mod(jmt-2,ncrow) .eq. 0) then
	    write (stdout,'(14x,a,i4,a,i4)')
     &      'jmw = ', jj, ',  number of MW moves = ',(jmt-2)/ncrow
	  endif
        enddo
      endif
#ifdef biharmonic
      if (jmw .lt. 4) then
        write (stdout,9000)
     & '==> Error: the MW can not have fewer than 4 rows (jmw >3)     '
     &,'           when using the "biharmonic" mixing option'
        write (stdout,*)'          you have set jmw=',jmw
        errorc = .true.  
      endif
      if (jmw .gt. 4) then
        write (stdout,9000)
     & '==> Warning: "jmw" > 4 ("jmw"=4 will use the minimum memory) '
      endif
#else
      if (jmw .lt. 3) then
        write (stdout,9000)
     & '==> Error: the MW can not have fewer than 3 rows (jmw >2)     '
        write (stdout,*)'          you have set jmw=',jmw
        errorc = .true.  
      endif
      if (jmw .gt. 3) then
        write (stdout,9000)
     & '==> Warning: "jmw" > 3 ("jmw"=3 will use the minimum memory) '
      endif
#endif
      if (jmw .gt. jmt) then
        write (stdout,9000)
     & '==> Error: the MW can not have more rows than "jmt"           '
        write (stdout,*)'          you have set jmw=',jmw, ', jmt=',jmt
        errorc = .true.  
      endif
      if (jmw .eq. jmt) then
        write (stdout,9000)
     & '==> Warning: The MW is open all the way ("jmw" = "jmt") which '
     &,'             is meant for "mpp" configurations. Note that the '
     &,'             latitude rows are kept in the MW and not on disk!' 
      endif
c
      if (nkflds .lt. 2) then
        write (stdout,9000)
     & '==> Error: "nkflds" must be at least 2                        '
        write (stdout,*)'          nkflds is set = ',nkflds
        errorc = .true.  
      endif
c
#ifdef equivalence_mw
        write (stdout,9000)
     & '==> Warning: The MW is equivalenced to the poisson solver     '
     &,'           workspace to reduce memory requirements. Verify    ' 
     &,'           that results are correct by comparing to a case    ' 
     &,'           with the "equivalence_mw" option removed.          ' 
      if (jmw .eq. jmt) then
        write (stdout,9000)
     & '==> Error: The MW cannot be equivalenced to the poisson solver'
     &,'           workspace but when "jmw" = "jmt"'
     &,'           remove the "equivalence_mw" option!' 
        errorc = .true.
      endif
#endif

c
#ifdef restorst
      if (dampts(1) .eq. c0 .or. dampts(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping time scale "dampts" must be set        '
     &,'           when enabeling the "restost" otpion                '
        errorc = .true.
      endif
      if (dampdz(1) .eq. c0 .or. dampdz(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping thickness "dampdz" must be set nonzero '
     &,'           when enabeling the "restost" otpion                '
        errorc = .true.
      endif
#else
      if (dampts(1) .ne. c0 .or. dampts(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping time scale "dampts" is > zero but    '
     &,'           the "restost" otpion is not enabled                '
      endif
      if (dampdz(1) .ne. c0 .or. dampdz(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping thickness "dampdz" is > zero but     '
     &,'           the "restost" otpion is not enabled                '
      endif
#endif
#if !defined ramdrive && !defined crayio && !defined fio
        write (stdout,9000)
     & '==> Error: either "ramdrive", crayio" or "fio"                '
     &,'           option must be turned on.                          '
        errorc = .true.
#endif
#if defined ramdrive && defined crayio
        write (stdout,9000)
     & '==> Error:  "ramdrive" and "crayio" are incompatible          '
     &,'             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined ramdrive && defined fio
        write (stdout,9000)
     & '==> Error:  "ramdrive" and "fio" are incompatible             '
     & '             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined crayio && defined fio
        write (stdout,9000)
     & '==> Error:  "crayio" and "fio" are incompatible               '
     &,'             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined crayio && !defined cray_ymp
        write (stdout,9000)
     & '==> Error:  "cray_ymp" must be enabled for "crayio" to work   '
        errorc = .true.
#endif
c
#ifdef isopycmix
      if (aidif .ne. p5) then
         write (stdout,9000)
     & '==> Warning:  isopycmix should use aidif = 0.5 for accuracy   '
     &,'          (this is the setting for the Crank-Nicholson scheme)'
      endif
# if defined consthmix && !defined biharmonic
        write (stdout,9000)
     & '==> Note: consthmix will only affect mixing of momentum       '
     &,'    since isopycmix was specified for tracer diffusion.       '
     &,'    kappa_h and Ah will be used as backgroud mixing coeffs    '
# endif
# if defined biharmonic
        write (stdout,9000)
     & '==> Error: biharmonic is incompatable with isopycmix          '
     &,'           (use "smagnlmix" instead)                          '
        errorc = .true.
# endif
#else
# ifdef gent_mcwilliams
        write (stdout,9000)
     & '==> Error:  "isopycmix" must be enabled for "gent_mcwilliams" '
     &,'            to work.                                          '
        errorc = .true.
# endif
#endif
c
c
c
c
#if defined implicit_free_surface && defined stream_function
        write (stdout,9000)
     & '==> Error:  "implicit_free_surface" and "stream_function" are '
     &,'             incompatible. choose one or the other option.    '
        errorc = .true.
#endif
#if defined rigid_lid_surface_pressure && defined stream_function
        write (stdout,9000)
     & '==> Error: "rigid_lid_surface_pressure" and "stream_function" '
     &,'           are incompatible. choose one or the other option.  '
        errorc = .true.
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
#ifdef diagnostic_surf_height
        write (stdout,9000)
     & '==> Warning:  "diagnostic_surf_height" does not apply under'
     &,'              the "rigid_lid_surface_pressure" or'
     &,'              "implicit_free_surface" options.'
#endif
#endif
#ifdef stream_function
# if !defined conjugate_gradient && !defined oldrelax
#  if !defined hypergrid
        write (stdout,9000)
     & '==> Error: either "conjugate_gradient", "oldrelax", or   '
     &,'         "hypergrid" option must be used with stream_function.'
        errorc = .true.
#  endif
# endif
# if defined oldrelax && defined hypergrid
        write (stdout,9000)
     & '==> Error:  "oldrelax" and "hypergrid" are incompatible       '
     &,'             choose one or the other option.                  '
        errorc = .true.
# endif
# if defined oldrelax && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error:  "oldrelax" and "conjugate_gradient" are'
     &,'             incompatible. choose one or the other option.    '
        errorc = .true.
# endif
# if defined oldrelax && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error:  "oldrelax" and "conjugate_gradient" are'
     &,'             incompatible. choose one or the other option.    '
        errorc = .true.
# endif
# if defined hypergrid && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error:  "hypergrid" and "conjugate_gradient" are'
     &,'             incompatible. choose one or the other option.    '
        errorc = .true.
# endif
# if defined hypergrid && defined conjugate_gradient
        write (stdout,9000)
     & '==> Error:  "hypergrid" and "conjugate_gradient" are          '
     &,'             incompatible. choose one or the other option.    '
        errorc = .true.
# endif
# ifdef sf_9_point
        write (stdout,9000)
     & '==> Warning: although "sf_9_point" conserves total energy     '
     &,'             it has a 2 grid point null space.                '
# endif
# ifdef sf_5_point
        write (stdout,9000)
     & '==> Warning: although "sf_5_point" has no null space, it does '
     &,'             not conserve total energy.                       '
# endif
# if defined sf_9_point && defined oldrelax
        write (stdout,9000)
     & '==> Error:  "oldrelax" cannot be used with 9 point numerics'
     &,'             use sf_5_point with oldrelax.'
        errorc = .true.
# endif
#endif
c
#ifdef free_slip
        write (stdout,9000)
     & '==> Warning: "free-slip" only accounts for no transfer of '
     &,'              tangential viscous stress to lateral boundaries.'
     &,'              Non-linear effects have been neglected for now,'
     &,'              but a complete "free_slip" is being researched' 
     &,'              Note: do not use yet... still needs work' 
#endif
c
#if defined biharmonic && defined free_slip
        write (stdout,9000)
     & '==> Error: "biharmonic" does not account for "free_slip"      '
        errorc = .true.
#endif
#if defined biharmonic && !defined consthmix
        write (stdout,9000)
     & '==> Error: "consthmix" must be enabled when using "biharmonic"'
        errorc = .true.
#endif
#if defined isopycmix && defined consthmix
        if ((ah+ahisop) .gt. 1.e11) then
          write (stdout,9000)
     & '==> Error: "ahisop"+"ah" is too large for the                 '
     &,'             "isopycmix" mixing option                        '
          errorc = .true.
        endif
#endif
c
        if (dtsf .lt. c0) then
          write (stdout,9000)
     & '==> Error: need to set the external mode time step "dtsf"     '
          errorc = .true.
	endif
c
        if (dtuv .lt. c0) then
          write (stdout,9000)
     & '==> Error: need to set the internal mode time step "dtuv"     '
          errorc = .true.
	endif
c
        if (dtts .lt. c0) then
          write (stdout,9000)
     & '==> Error: need to set the density time step "dtts"           '
          errorc = .true.
	endif
c
        critv = 1.e-6
        if (mod(rundays,dtts*secday) .gt. critv) then
          t1 = nint(rundays/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: there must be an integral number of density time '
     &,'    steps within "rundays" (the integration time).            '
          write (stdout,*) '               (changed "rundays" from '
     &,   rundays,' days to ', t1,' days to insure this condition)'
          rundays = t1
        endif
c
        if (itmb) then
          write (stdout,9000)
     & '==> Warning: "itmb" is set to "true". set it to "false" in    '
     &,'            subsequent runs to prevent the time independent   '
     &,'            basin mask from being written more than once. This'
     &,'            reduces the size of the diagnostic file.          '
	endif
c
        if (itrmb) then
          write (stdout,9000)
     & '==> Warning: "itrmb" is set to "true". set it to "false" in   '
     &,'            subsequent runs to prevent the time independent   '
     &,'            region masks from being written more than once.   '
     &,'            This reduces the size of the diagnostic file.     '
	endif
c
        if (itavg) then
          write (stdout,9000)
     & '==> Warning: "itavg" is set to "true". set it to "false" in   '
     &,'            subsequent runs to prevent the time independent   '
     &,'            region masks from being written more than once.   '
     &,'            This reduces the size of the diagnostic file.     '
	endif
#ifdef trajectories
        if (trajint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: "trajectories" is set but "trajint"=0. this means'
     &,'            particle positions will be written every time step'
        endif
        if (trajint .lt. c0) then
          write (stdout,9000)
     & '==> Error: "trajectories" is set but "trajint" < 0. this means'
     &,'            particle positions will not be written.           '
          errorc = .true.
        endif
# ifndef lyapunov
          write (stdout,9000)
     & '==> Warning: "trajectories" is enabled but "lyapunov" is not. '
# endif
#else
# ifdef lyapunov
          write (stdout,9000)
     & '==> Warning: "trajectories" is not enabled but "lyapunov" is. '
# endif
#endif
#ifdef meridional_tracer_budget
        if (mod(tmbint,dtts*secday) .gt. critv .and. tmbint .gt. c0)
     &    then
          t1 = nint(tmbint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "tmbint" does not contain an integral number of  '
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "tmbint" from '
     &,   tmbint,' days to ', t1,' days to insure this condition)'
          tmbint = t1
        endif
        if (tmbint .eq. c0) then
          write (stdout,9000)
     & '==> Error: averaging period "tmbint"  must be specified when  '
     &,'   the "meridional_tracer_budget" option is enabled           '
          errorc = .true.
        endif
#else
      if (tmbint .gt. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "tmbint" is > zero but the  '
     &,'             the "meridional_tracer_budget" option is not on. '
      endif
#endif
#ifdef time_averages
        if (mod(timavgint,dtts*secday).gt.critv .and. timavgint .gt. c0)
     &    then
          t1 = nint(timavgint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "timavgint" does not contain an integral number of'
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "timavgint" from '
     &,   timavgint,' days to ', t1,' days to insure this condition)'
          timavgint = t1
        endif
        if (timavgint .eq. c0) then
          write (stdout,9000)
     & '==> Error: averaging period "timavgint" must be specified when'
     &,'   the "time_averages" option is enabled                      '
          errorc = .true.
        endif
# ifndef time_averages_disk
          write (stdout,9000)
     & '==> Warning: "option "time_averages" is enabled but           '
     &,      '"time_averages_disk" is not. extra memory will be used. '
     &,      'Note that"time_averages_disk" does not work with option '
     &,      ' "ramdrive"                                             '
# endif
#else
      if (timavgint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "timavgint"  is > zero but  '
     &,'           the "time_averages" option is not enabled          '
      endif
#endif
#if defined time_averages_disk && !defined time_averages
          write (stdout,9000)
     & '==> Warning:  "time_averages_disk" is enabled but '
     &,' "time_averages" is not!  '
#endif
#if defined time_averages_disk && defined ramdrive
          write (stdout,9000)
     & '==> Error: "time_averages_disk" is incompatible with "ramdrive"'
        errorc = .true.
#endif
#ifdef xbts
        if (mod(xbtint,dtts*secday) .gt. critv .and. xbtint .gt. c0)
     &    then
          t1 = nint(xbtint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "xbtint" does not contain an integral number of  '
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "xbtint" from '
     &,   xbtint,' days to ', t1,' days to insure this condition)'
          xbtint = t1
        endif
        if (xbtint .eq. c0) then
          write (stdout,9000)
     & '==> Error: averaging period "xbtint"  must be specified when  '
     &,'   the "xbts" option is enabled                               '
          errorc = .true.
        endif
        if (numxbt .eq. 0) then
          write (stdout,9000)
     &   '==> Error: the number of XBT stations detected is zero      '
     &,  '           at least one station must be specified           '
          errorc = .true.
        endif
#else
      if (xbtint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "xbtint"  is > zero but     '
     &,'             the "xbts" option is not enabled                 '
      endif
#endif
c
      if ((dtuv .ne. dtsf) .or. (dtuv .ne. dtts)) then
        write (stdout,9000)
     & '==> Warning: use of unequal time steps implies the transient  '
     &,'             response is unimportant and multiple equilibria  '
     &,'             do not exist.                                    '
      endif
c
c     check for mixing coefficients larger than stability permits
c
      dymin  = dyt(2)
      dxmin  = dxt(2)
      dzmin  = dzt(1)
      xlmax  = dtxcel(1)
#ifdef isopycmix
      fzimax = fzisop(1)
#endif
      do jrow=2,jmtm1
        dymin  = min(dymin,dyt(jrow))
      enddo
      do i=2,imtm1
        dxmin  = min(dxmin,dxt(i))
      enddo
      do k=2,km
        xlmax  = max(xlmax,dtxcel(k))
        dzmin  = min(dzmin,dzt(k))
#ifdef isopycmix
        fzimax = max(fzimax,fzisop(k))
#endif
      enddo
c
      if (xlmax .gt. c1) then
        write (stdout,9000)
     & '==> Warning: use of accelerated time steps implies the        '
     &,'             transient response is unimportant and multiple   '
     &,'             equilibria do not exist. stability tests will    '
     &,'             use "dtts" multiplied by the maximum "dtxcel"    '
      endif
c
c
c
#if defined consthmix
      dtxl = dtts*xlmax
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
# ifdef isopycmix
        if ((dtxl*(ah+ahisop*fzimax))/dxdymn .ge. p25) then
# else
        if ((dtxl*ah)/dxdymn .ge. p25) then
# endif
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
# ifdef isopycmix
     & '==> Warning: lateral diffusive criteria exceeded for "ah" +   '
     &,'             "ahisop". use a smaller "dtts", "dtxcel", and/or '
     &,'             "ah" + "ahisop"                                  '
# else
     & '==> Warning: lateral diffusive criteria exceeded for "ah".    '
     &,'              use a smaller "dtts", "dtxcel", and/or "ah"     '
# endif
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
        if ((dtuv*am)/dxdymn .ge. p25) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: lateral diffusive criteria exceeded for "am".    '
     &,'             use a smaller "dtuv" and/or "am"                 '
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
#endif
#ifdef sponges
# ifndef source_term
      write (stdout,9000)
     & '==> Error: ifdef "sponges" is enabled but "source_term"     '
     &,'             is not. "source_term" needs to be enabled  '
      errorc = .true.
# endif
#endif
#ifdef shortwave
# ifndef source_term
      write (stdout,9000)
     & '==> Error: ifdef "shortwave" is enabled but "source_term"     '
     &,'             is not. "source_term" needs to be enabled  '
      errorc = .true.
# endif
# ifndef coupled
      write (stdout,9000)
     & '==> Warning: ifdef "shortwave" is enabled but "coupled"       '
     &,'             is not. Is the shortwave component of radiation  '
     &,'             defined ?                                        '
# endif
#else
      if (dzt(1) .lt. 20.0e2) then
        write (stdout,9000)
     & '==> Warning: if shallow mixed layers develop, then enabling   '
     &,'             ifdef "shortwave" may help to deepen them. note  '
     &,'             that either you or the atmosphere must provide   '
     &,'             the solar short wave as a boundary condition.    '
      endif
#endif
#ifdef constvmix
# ifdef implicitvmix
        write (stdout,9000)
     & '==> Error: "implicitvmix" is not allowed with "constvmix"     '
     &,'            use option "fullconvect" instead of "implicitvmix"'
        errorc = .true.
# endif
# if defined implicitvmix && defined isopycmix
        write (stdout,9000)
     & '==> Error: "implicitvmix" is not allowed with "constvmix"     '
     &,'            and "isopycmix". "isopycmix" implies implicit     '
     &,'            vertical mixing for tracers but not for momentum. '
     &,'            Vertical mixing schemes other than "constvmix"    '
     &,'            would allow "implicitvmix" with "isopycmix"       '
        errorc = .true.
# endif
      do k=1,km
        if ((dtts*dtxcel(k)*kappa_h)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_h"'
     &,'             use a smaller "dtts", "dtxcel", and/or "kappa_h" '
         write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*kappa_m)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_m"'
     &,'           use a smaller "dtuv" and/or "kappa_m"              '
      endif
#endif
#ifdef implicitvmix
      write (stdout,9000)
     & '==> Warning: "implicitvmix" turns convective adjustment off.  '
     &,'              It will be handled by the vertical mixing       '
#else
# ifdef fullconvect
      write (stdout,9000)
     & '==> Warning: the full convective scheme is enabled.           '
     &,'             it will ignore "ncon" and remove all instability '
# else
      write (stdout,9000)
     & '==> Warning: the standard convective scheme is enabled. Verify'
     &,'              that "ncon" is set to what is wanted:           '
      write (stdout,'(32x,a,i3)') '"ncon" is currently set to ',ncon
# endif
#endif
#if defined tcvmix && !defined implicitvmix
      do k=1,km
        if ((dtts*dtxcel(k)*diff_cbt_limit)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'            "diff_cbt_limit".  use a smaller "dtts", "dtxcel" '
     &,'            ,and/or "diff_cbt_limit"                          '
          write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*visc_cbu_limit)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'             "visc_cbu_limit". use a smaller "dtuv" or        '
     &,'             "visc_cbu_limit"                                 '
      endif
#endif
c
c     check range of implicit factors
c
#if defined implicitvmix || defined isopycmix
      if (aidif .ne. p5) then
        write (stdout,9000)
     & '==> Warning: implicitvmix should use aidif = 0.5 for accuracy '
     &,'          (this is the setting for the Crank-Nicholson scheme)'
      endif
#else
      if (aidif .ne. c0) then
        write (stdout,9000)
     & '==> Warning: "aidif" was not = 0.0  but "implicivmix" is not  '
     &,'              defined, so vertical mixing  will be  explicit  '
      endif
#endif
      if (acor .lt. c0 .or. acor .gt. c1) then
        write (stdout,9000)
     & '==> Error: "acor" must be between 0.0 & 1.0                   '
          errorc = .true.
      endif
c
c     check for resolved inertial oscillation
c
      if (acor .eq. c0) then
        fmax = max(abs(c2*omega*sine(2)),epsln)
        do jrow=3,jmtm2
          fmax = max(fmax,abs(c2*omega*sine(jrow)))
        enddo
        if (dtuv .gt. p25*(c2*pi)/fmax) then
          write (stdout,9000)
     & '==> Warning: the inertial oscillation is not resolved.        '
     &,'             reduce "dtuv" or use "acor" > 0.0                '
        endif
      endif
#ifdef consthmix
c
c-----------------------------------------------------------------------
c     search for topographic instabilities (based  on the  work of
c     Peter Killworth  ...  eqn 11 from ocean modeling nov 1987)
c-----------------------------------------------------------------------
c
      num   = 50
      do j=2,jmtm1
        dysq = dyt(j)**2
        do i=2,imtm1
          if (kmu(i+1,j-1) .ne. 0 .and. kmu(i+1,j) .ne. 0) then
            dxsq = (dxt(i)*cst(j))**2
            clix = am*dtuv/dxsq
            h1   = zw(kmu(i+1,j-1))
            h2   = zw(kmu(i+1,j))
            hx   = (8.0*h1*h2/(h1+h2)**2 + dxsq/dysq)/(4.0 + dxsq/dysq)
            if (clix .ge. hx .and. num .ge. 0) then
              num = num - 1
              write (stdout,'(1x,a56,a31,i4,a1,i4,a1)')
     &        '==>   Warning: Killworth topographic roughness condition'
     &,       ' exceeded at location (i,j) = (',i+1,',',j,')'
              if (num .eq. 0) then
                write (stdout,9000)
     &         '==> Warning: msgs terminated after 50 cases were found '
              endif
            endif
          endif
        enddo
      enddo
#endif
c
c     verify that the domain boundary is valid
c
#if !defined cyclic && !defined solid_walls
        write (stdout,9000)
     & '==> Error: Either "cyclic" or "solid_walls" must be enabled'
        errorc = .true.
#endif
#ifdef solid_walls
      jw = 0
      je = 0
      do j=1,jmt
        if (kmt(1,j) .ne. 0) jw = j
	if (kmt(imt,j) .ne. 0) je = j
      enddo
      if (je .ne. 0 .or. jw .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the eastern or western boundary. This is not       '
     &,'           correct when enabeling option "solid_walls"        '
        if (jw .ne. 0) write (stdout,*) '           at i=1 and j=',jw
        if (je .ne. 0) write (stdout,*) '           at i=imt and j=',je
      endif 
#endif
      in = 0
      is = 0
      do i=1,imt
        if (kmt(i,1) .ne. 0) is = i
	if (kmt(i,jmt) .ne. 0) in = i
      enddo
      if (is .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the southern boundary.                             '
        write (stdout,*) '           at j=1 and i=',is
      endif 
#ifndef symmetry
      if (in .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the northern boundary.                             '
        write (stdout,*) '           at j=jmt and i=',in
      endif 
#endif
c
c     verify that each ocean point is at least 2 levels deep
c
      numk = 0
      do jrow=1,jmt
        do i=1,imt
	  if (kmt(i,jrow) .eq. 1) then
            numk = numk + 1
	    errorc = .true.
	    write (stdout,*)
     &	            ' Error: kmt(',i,',',jrow,') = 1 is not allowed'
	  endif
        enddo
      enddo
      if (numk .ne. 0) then	
        write (stdout,9000)
     & '==> Error: "kmt" must be at least 2 levels deep at all ocean  '
     &,'           points.                                            '
      endif
c
      write (stdout,'(/,20x,a,/)')
     &         ' E N D    C O N S I S T E N C Y    C H E C K S'
      if (errorc) stop '=>checks'
c
      return
9000  format (/,(1x,a))
      end
