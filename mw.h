c======================= include file "mw.h" ===========================
c
c                     M E M O R Y    W I N D O W
c
c                         Data layout on disk
c
c     On disk, a latitude consists of a row of prognostic quantities. 
c     Latitudes are ordered from south to north on two logical disk
c     units: one for time level "tau" and one for time level "tau-1".
c     Newly computed "tau+1" quantities are written to the "tau-1" disk
c     which serves double duty as a "tau+1" disk.
c     
c                         Data layout in memory
c
c     A memory window "MW" is constructed to hold a group of "jmw"
c     adjacent latitude rows of prognostic quantities from disks "tau-1"
c     and "tau".  Parameter "jmw" controls the size of the MW and can be
c     set anywhere from a miminum of three latitude rows to all "jmt"
c     latitudes.  In addition, the MW holds diagnostic quantities
c     (density, hydrostatic pressure gradients, and advective
c     velocities) along with work arrays for constructing intermediate
c     fluxes used in solving the tracer and momentum equations. A
c     latitude row in the MW is referred to by index "j" and corresponds
c     to the physical latitude row "jrow" on disk.
c
c                   Data flow between disk and memory
c
c     The MW is loaded with prognostic data from the first "jww"
c     latitude rows on disks "tau" and "tau-1". As the tracer and
c     momentum equations are solved for rows j=2 through jmw-1 in the
c     MW, the solutions are written to disk "tau+1". When finished, the
c     MW is moved northward by moving quantities from the northernmost
c     two rows into the southernmost two rows within the MW. The
c     remaining MW rows are then loaded with more latitude rows from
c     disk. The process continues until all latitude rows 2 through
c     jmt-1 on disk "tau+1" have been updated.
c
c                   Parallel processing considerations
c
c     When the MW is opened all the way ("jmw" = "jmt"), there is no
c     need for disk units since all three time levels for all prognostic
c     quantities at all latitudes reside in the MW. Therefore, there is
c     no need to read data, write data or move the MW northward (or
c     re-calculate the three dimensional land/sea mask arrays). After a
c     time step has finished, time level pointers on the next time step
c     are simply updated to reflect where the "tau-1", "tau" and "tau+1"
c     data is. In principle, MOM 2 can be viewed as a sucession of
c     triply nested do loops with parallelism exposed at the do loop
c     level (fine grained parallelism) as opposed to the latitude row
c     level (coarse grained parallelism) of MOM 1. This structure
c     significantly reduces the memory requirements compared to MOM 1
c     when executing on multiple processors. On one processor, the
c     memory requirement is comparable to MOM 1.
c
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
c     taum1 = tau-1 time level for variables in MW
c     tau   = tau   time level for variables in MW
c     taup1 = tau+1 time level for variables in MW
c
      integer taum1, tau, taup1
      common /mwi/ taum1, tau, taup1
c
c-----------------------------------------------------------------------
c     MW arrays for prognostic equations:
c-----------------------------------------------------------------------
c
c     u(i,k,j,n,tau) = total velocity where:
c      i   = index for longitude
c      k   = index for depth
c      j   = index for latitude row within MW
c      n   = component (1 = zonal, 2 = meridional)
c      tau = time level (tau-1, tau, tau+1)
c             
c      (only internal modes are on disk and at tau+1 in the MW)
c
c     t(i,k,j,n,tau) = tracer where:
c      i   = index for longitude
c      k   = index for depth
c      j   = index for latitude row within MW
c      n   = component (1 = temperature, 2 = salinity)
c            if nt > 2 then other tracers are allowed.
c      tau = time level (tau-1, tau, tau+1)
c
c     note: temperature is potential temperature in degrees Celsius and 
c           salinity is in "model units", the deviation from 0.035 grams
c           of salt/cm**3 of water, or, assuming a water density of
c           1 gram/cm**3, the deviation from 0.035 g of salt/g of water.
c           one can convert model units to the more common units of parts
c           per thousand (ppt) by adding 0.035 grams/cm**3 to the model
c           units and then multiplying by 1000. 
c     
      common /mw/ u(imt,km,jmw,2,-1:1), t(imt,km,jmw,nt,-1:1)
c
c-----------------------------------------------------------------------
c     MW arrays for diagnostic equations and workspace:
c-----------------------------------------------------------------------
c
c     diagnostic advective velocities are in units of cm/sec
c
c     adv_vet  = advective velocity on the eastern face of a "T" cell
c     adv_vnt  = advective velocity on the northern face of a "T" cell
c     adv_veu  = advective velocity on the eastern face of a "u" cell
c     adv_vnu  = advective velocity on the northern face of a "u" cell
c     adv_vbt  = advective velocity on the bottom face of a "T" cell
c     adv_vbu  = advective velocity on the bottom face of a "u" cell
c
c     rho      = density at center of a "T" cell in units of gm/cm**3
c                note: there is an arbitrary constant which is only a
c                function of depth in "rho". It is related to 
c                subtracting a reference level density for purposes of
c                accuracy.
c
c     grad_p   = hydrostatic pressure gradient for "u" cell. There are
c                two components: (1,2) is for (dp/dx, dp/dy)
c
      common /mw/ adv_vet(imt,km,jsmw:jmw), adv_vnt(imt,km,1:jmw)
      common /mw/ adv_veu(imt,km,jsmw:jemw)
      common /mw/ adv_vnu(imt,km,1:jemw)
      common /mw/ adv_vbt(imt,0:km,jsmw:jmw)
#ifdef trajectories
      common /mw/ adv_vbu(imt,0:km,1:jemw)
#else
      common /mw/ adv_vbu(imt,0:km,jsmw:jemw)
#endif
      common /mw/ rho(imt,km,jsmw:jmw)
      common /mw/ grad_p(imt,km,jsmw:jemw,2)
c
c     tmask = tracer cell land/sea mask   = (0.0, 1.0) on (land, sea)
c     umask = velocity cell land/sea mask = (0.0, 1.0) on (land, sea)
c
      common /mw/ tmask(imt,km,1:jmw), umask(imt,km,1:jmw)
c
#ifdef biharmonic
c
c     delsq = del**2 of prognostic variables
c
# ifdef tcvmix
      parameter (nvarbh=nt+2+1)
#  ifdef leq
      parameter (nvarbh=nt+2+1+1)
#  endif
# else
      parameter (nvarbh=nt+2)
# endif
      common /mw/ del2(imt,km,1:jmw,nvarbh)
#endif
c
c
c     these workspace arrays are recalculated for each component of the
c     equations so do not have to be moved as the MW moves northward. 
c
c     adv_fe   = advective flux across the eastern face of a cell
c     adv_fn   = advective flux across the northern face of a cell
c                (removed and put directly into the statement functions
c                 as a trial optimization. It saves space at the expense
c                 of redundant calculation when the MW is opened > 3)
c     adv_fb   = advective flux across the bottom face of a cell
c
c     diff_fe  = diffusive flux across the eastern face of a cell
c     diff_fn  = diffusive flux across the northern face of a cell
c     diff_fb  = diffusive flux across the bottom face of a cell
c     source   = source term
c
      common /mw/ adv_fe(imt,km,jsmw:jemw)
c     common /mw/ adv_fn(imt,km,1:jemw)
      common /mw/ adv_fb(imt,0:km,jsmw:jemw)
c
      common /mw/ diff_fe(imt,km,jsmw:jemw)
#if !defined consthmix || defined biharmonic || defined isopycmix
      common /mw/ diff_fn(imt,km,1:jemw)
#endif
      common /mw/ diff_fb(imt,0:km,jsmw:jemw)
#if defined isopycmix
      common /mw/ diff_fbiso(imt,0:km,jsmw:jemw)
#endif
c
#ifdef source_term
      common /mw/ source(imt,km,jsmw:jemw)
#endif
c
#if defined implicitvmix || defined isopycmix
      common /mw/ zzi(imt,km,jsmw:jemw)
#endif
c
c     these grid factors are for optimizations and are recalculated as
c     the MW moves northward so they do not have to be moved.
c     
      common /grdfact/ cstdxtr(imt,jsmw:jemw), cstdxur(imt,jsmw:jemw)
      common /grdfact/ cstdxt2r(imt,jsmw:jemw),ah_cstdxur(imt,jsmw:jemw)
      common /grdfact/ rho0csudxur(imt,jsmw:jemw)
      common /grdfact/ csudxtr(imt,jsmw:jemw)
      common /grdfact/ csudxu2r(imt,jsmw:jemw),am_csudxtr(imt,jsmw:jemw)
c
c     these variables are either constant or globally dimensiond by
c     "jmt", so they do not need to be moved as the MW moves northward
c     
c     advmet = coeff for metric advection.
c     cori   = coriolis parameter for velocity component "n"
c
      common /advec/ advmet(jmt,2)
c
      common /coriol/ cori(jmt,2)
#ifdef equivalence_mw
c
c     if the MW storage is > workspace for the poisson solvers, then
c     the poisson solver workspace may be equivalenced to the MW to
c     save space. this is not allowed when the MW is opened all the way
c     amount of savings is dependent on problem size
c
      dimension cf(imt,jmt,-1:1,-1:1)
      equivalence (cf,u(1,1,1,1,-1))
#else
      common /extwrk/ cf(imt,jmt,-1:1,-1:1)
#endif
c
#ifdef linearized_advection
c
c     initial temperature structure used for linearlized advection
c
      common /linadv/ tbarz(km)
#endif
#include "cvbc.h"
c
