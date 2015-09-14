c====================== include file "emode.h" =========================
c
c
c     variables for  external mode
c
c     psi   = stream function (,,1) is for tau; (,,2) is for tau-1
c     zu    = vertically averaged forcing from momentum equations
c             (,,1) is zonal and (,,2) is meridional component
c     ztd   = curl of "zu" for the stream function equation
c     ptd   = time change of stream function
c     h     = depth over "u" points
c     hr    = reciprocal depth over "u" points
c     res   = residual from elliptic solver
c
c     map   = land mass map distinguishing, ocean, land, and perimiters
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c     surface pressure/free surface time centering parameters:
c
c     alph,gam = parameters to center the coriolis and surface pressure 
c                gradients in time (leapfrog)
c
c     theta    =parameters to center the coriolis and surface pressure 
c                gradients in time (mixing ts )
c
c     apgr     = is = to alpha/theta leapfrog/mixing ts
c
# ifdef implicit_free_surface
c
c     to set free surface time centering parameters...
c     note: for proper time centering of coriolis
c     and pressure gradients alph, gam should
c     satisfy gam = 1 -2*alph.  for stability
c     alph should be > 1/4.  recommended values
c     are alph = gam = 1/3.
# endif
      common /emode/  alph, gam, theta, apgr
#endif
c
c     mxscan  = max number of allowable scans for poisson solvers
c     mscan   = actual number of scans taken by poisson solvers
c     sor     = sucessive over-relaxation constant
c     tolrsf  = tolerence for stream function calculation. 
c               the solution is halted when it is within "tolrsf"
c               of the "true" solution assuming geometric convergence.
c     tolrsp  = tolerence for surface pressure calculation
c               the solution is halted when it is within "tolrsp"
c               of the "true" solution assuming geometric convergence.
c     tolrfs  = tolerence for implicit free surface calculation
c               the solution is halted when it is within "tolrfs"
c               of the "true" solution assuming geometric convergence.
c     esterr  = estimated maximum error in elliptic solver assuming
c               geometric convergence
c    
      common /emode/  mxscan, mscan, tolrsf, tolrsp, tolrfs, sor
      common /emode/  esterr
c
c     nisle = number of land masses
c     nippts= number of land mass perimeter points
c     iperm = "i" coordinate for the land mass perimeter point
c     jperm = "j" coordinate for the land mass perimeter point
c     iofs  = offset for indexing into the land mass perimeter points
c     imask = controls whether calculations get done on perimeters
c     set mask for land mass perimeters on which to perform calculations
c     imask(-n) = .false.  [no equations ever on dry land mass n]
c     imask(0)  = .true.   [equations at all mid ocean points]
c     imask(n)  = .true./.false [controls whether there will be
c                                equations on the ocean perimeter of
c                                land mass n]
c     note: land mass 1 is the northwest-most land mass
c     for the numbering of the other landmasses, see generated map(i,j)
      logical imask
      common /emodei/ map(imt,jmt)
      common /emodei/ nippts(mnisle), iofs(mnisle), iperm(maxipp)
      common /emodei/ jperm(maxipp), nisle, imain
      common /emodel/ imask (-mnisle:mnisle)
c
      common /emode/ ptd(imt,jmt), res(imt,jmt), hr(imt,jmt), h(imt,jmt)
      common /emode/ zu(imt,jmt,2)
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c     ps    = surface pressure (,,1) is for tau; (,,2) is for tau-1
c     divf  = barotropic divergence of uncorrected ubar & vbar (rhs of 
c             surface pressure eqn)
c     ubar  = barotropic velocity defined on "u" point for "tau" 
c             (,,1) is zonal and (,,2) is meridional velocity
c     ubarm1= barotropic velocity defined on "u" point for "tau-1" 
c             (,,1) is zonal and (,,2) is meridional velocity
      common /emode/ uhat(imt,jmt,2),   pguess(imt,jmt)
      common /emode/ ps(imt,jmt,2)
      common /emode/ divf(imt,jmt),   ubar(imt,jmt,2), ubarm1(imt,jmt,2)
#endif
c
#ifdef stream_function
      common /emode/ psi(imt,jmt,2), ztd(imt,jmt)
#endif
      character*16 variable
      logical converged
      common /emodel/ converged
      common /emodec/ variable
c
