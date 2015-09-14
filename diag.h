c====================== include file "diag.h" ==========================
c
c     variables used for computing diagnostics:
c
c     tcella  = "t" cell surface area cm**2 (entire ocean)
c     ucella  = "u" cell surface area cm**2 (entire ocean)
c     tcellv  = "t" cell volume cm**3 (entire ocean)
c     ucellv  = "u" cell volume cm**3 (entire ocean)
c
      common /cdiag/ tcellv, ucellv, tcella(km), ucella(km)
c
#ifdef time_step_monitor
c
c     ektot    = "total" kinetic energy per unit volume at "tau". units
c                ergs/cm**3 = dyn/cm**2 = g/cm/sec**2 = 10**-7 J/cm**3.
c                ektot is the "total" ke in the sense that it considers
c                both the internal and external modes summed over the
c                entire ocean volume. The contributions of
c                vertical motions are neglected on the basis of scaling
c                arguments (i.e., w**2 << (u**2 + v**2).
c     dtabs    = absolute value of rate of change of tracer per unit
c                volume centered at "tau"
c     tbar     = first moment of tracer at "tau"
c     travar   = variance = second moment of tracer about mean at "tau"
c
      common /cdiag/ ektot(0:km,jmt), dtabs(0:km,nt,jmt)
      common /cdiag/ travar(0:km,nt,jmt), tbar(0:km,nt,jmt)
#endif
#ifdef energy_analysis
c
c     engint   = volume averaged internal mode energy integral
c                components
c     engext   = volume averaged external mode energy integral
c                components
c     buoy     = volume averaged buoyancy
c
c     tcerr    = maximum "t" cell continuity error
c     ucerr    = maximum "u" cell continuity error
c     itcerr   = "i" index corresponding to "tcerr"
c     jtcerr   = "jrow" index corresponding to "tcerr"
c     ktcerr   = "k" index corresponding to "tcerr"
c     iucerr   = "i" index corresponding to "ucerr"
c     jucerr   = "jrow" index corresponding to "ucerr"
c     kucerr   = "k" index corresponding to "ucerr"
c
c     wtbot    = maximum "adv_vbt" error at ocean bottom
c     iwtbot   = "i" index corresponding to "wtbot"
c     jwtbot   = "jrow" index corresponding to "wtbot"
c     kwtbot   = "k" index corresponding to "wtbot"
c     wubot    = maximum "adv_vbu" at ocean bottom
c     iwubot   = "i" index corresponding to "wubot"
c     jwubot   = "jrow" index corresponding to "wubot"
c     kwubot   = "k" index corresponding to "wubot"
c
c     wtlev    = zonally integrated adv_vbt for each level
c     wulev    = zonally integrated adv_vbu for each level
c
      common /cdiag/ buoy(0:km,jmt), engint(0:km,8,jmt), engext(8,jmt)
      common /cdiag/ tcerr(jmt), ucerr(jmt)
      common /cdiagi/ itcerr(jmt), jtcerr(jmt), ktcerr(jmt) 
      common /cdiagi/ iucerr(jmt), jucerr(jmt), kucerr(jmt) 
      common /cdiag/ wtbot(jmt), wubot(jmt) 
      common /cdiagi/ iwtbot(jmt), jwtbot(jmt), kwtbot(jmt)
      common /cdiagi/ iwubot(jmt), jwubot(jmt), kwubot(jmt)
      common /cdiag/ wtlev(km,0:jmt), wulev(km,0:jmt)
#endif
#ifdef term_balances
# include "termbal.h"
#endif
#ifdef gyre_components
c
c     ttn      = northward transport of tracer components
c
c     ttn2     = northward transport of tracers for ocean basins
c                  (.,.,.,0)       Global
c                  (.,.,.,1:nhreg) Ocean basins
c                also,
c                  (6,.,.,.) total transport due to advection 
c                  (7,.,.,.) total transport due to diffusion
c                  (8,.,.,.) total transport
# if defined isopycmix && defined gent_mcwilliams
c                  (9,.,.,.) total transport due to isopycnal advection
# endif
      common /gyres/ ttn(8,jmt,ntmin2)
# if defined isopycmix && defined gent_mcwilliams
      common /gyres/ ttn2(6:9,jmt,nt,0:nhreg)
# else
      common /gyres/ ttn2(6:8,jmt,nt,0:nhreg)
# endif
#endif
#ifdef meridional_overturning
c
c     vsf      = vertical_meridional stream function
c
      common /cdiag/ vsf(jmt,km)
#endif
#ifdef show_zonal_mean_of_sbc
c
c     zmsmf    = zonal mean surface momentum flux
c     zmstf    = zonal mean surface tracer flux
c     zmsm     = zonal mean surface momentum
c     zmst     = zonal mean surface tracers
c     zmau     = surface area weighting for "u" latitudes
c     zmat     = surface area weighting for "t" latitudes
c
      common /cdiag/ zmsmf(jmt,2), zmstf(jmt,nt), zmau(jmt), zmat(jmt)
      common /cdiag/ zmsm(jmt,2), zmst(jmt,nt)
#endif
c
