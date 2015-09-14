c====================== include file "index.h" =========================
c
c    various starting & ending indices for controlling where quantities
c    are filtered. this removes the time step restriction due to
c    convergence of meridians.
c
#if defined fourfil || defined firfil
c    istf  = starting index for filtering "t" grid points
c    ietf  = ending index for filtering "t" grid points
c    isuf  = starting index for filtering "u,v" grid points
c    ieuf  = ending index for filtering "u,v" grid points
c    iszf  = starting index for filtering "vorticity" grid points
c            (also for divergence when using implicit free surface)
c    iezf  = ending index for filtering "vorticity" grid points
c            (also for divergence when using implicit free surface)
c
c
c     define latitudinal domain for filtering
c
c     filter t to yield equiv dx at lat rjft0 from lat rjfrst to rjft1
c       and lat rjft2 to yt(jmtm1)
c     filter u to yield equiv dx at lat rjfu0 from lat rjfrst to rjfu1
c       and lat rjfu2 to yt(jmtm2)
c
c     lsegf  = max number of longitudinal strips per filtering latitude 
c     jmtfil = max number of latitudes to be filtered
c     numflt = specifies the number of filter applications to use
c              for the tracers, vorticity, or divergence when using the
c               finite impulse response filter
c     numflu = specifies the number of filter applications to use
c              for velocities when using the finite impulse response
c               filter
c
c
c
      parameter (lsegf=7, jmtfil=21)
c
      common /indexr/ rjfrst, rjft0, rjft1, rjft2, rjfu0, rjfu1, rjfu2
      common /indexi/   jfrst,  jft0,  jft1,  jft2,  jfu0,  jfu1,  jfu2
      common /indexi/   jskpt, jskpu, njtbft, njtbfu
# ifdef firfil
      common /indexi/   numflt(jmtfil), numflu(jmtfil)
# endif
      common /indexi/ istf(jmtfil,lsegf,km), ietf(jmtfil,lsegf,km)
      common /indexi/ isuf(jmtfil,lsegf,km), ieuf(jmtfil,lsegf,km)
      common /indexi/ iszf(jmtfil,lsegf),    iezf(jmtfil,lsegf)
#endif
c
