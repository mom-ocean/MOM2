c====================== include file "vmixc.h" =========================
c
c         vertical mixing coefficients and related variables
c
c     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
c     kappa_m = constant vertical viscosity coefficient (cm**2/sec)
c
c     visc_cbu  = viscosity coeff at bottom of "u" cell 
c     diff_cbt  = diffusion coeff at bottom of "T" cell
c     visc_cbu_limit = largest allowable "visc_cbu" 
c     diff_cbt_limit = largest allowable "diff_cbt" 
c     aidif = coefficient for implicit time differencing for
c             vertical diffusion. aidif=1 gives the fully implicit
c             case. aidif=0 gives the fully explicit case
c             recommendation: aidif=0.5 for accuracy
c             note: not used unless "implicitvmix" or "isopycmix"
c                   is enabled
c
c     author:   r.c.pacanowski     e-mail  rcp@gfdl.gov
c=======================================================================
c
      common /vmixr/ visc_cbu_limit, diff_cbt_limit, aidif
      real kappa_h,  kappa_m
      common /vmixr/ kappa_h, kappa_m
c
#ifdef constvmix
c
c     variables for constant vertical mixing of momentum
c
      real kappa_m_dzwr
      common /vmixr/ visc_cbu, kappa_m_dzwr(km)
c
c     variables for constant vertical mixing of tracers
c
# if defined bryan_lewis_vertical
      common /vmixr/ Ahv(km)
# endif
# if defined isopycmix
      common /vmixr/ diff_cbt(imt,km,jsmw:jemw)
# else
      real kappa_h_dzwr
      common /vmixr/ kappa_h_dzwr(km)
#  if defined bryan_lewis_vertical
      common /vmixr/ diff_cbt(km)
#  else
      common /vmixr/ diff_cbt
#  endif
# endif
#endif
c
#ifdef ppvmix
c
c     variables for pacanowski-philander vertical diffusion
c
c     rhom1z = d(rho)/dz at bottom of "T" cell at tau-1
c     riu    = richardson number at bottom of "u" cell
c     rit    = richardson number at bottom of "T" cell
c     fricmx = max vertical mixing coefficient
c     wndmix = min vertical mixing in level 1 to simulate wind mixing
c     diff_cbt_back = background "diff_cbt"
c     visc_cbu_back = background "visc_cbu"
c
      common /vmixr/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
      common /vmixr/ rhom1z(imt,km,1:jmw)
      common /vmixr/ riu(imt,km,1:jemw), rit(imt,km,jsmw:jemw)
      common /vmixr/ visc_cbu(imt,km,jsmw:jemw)
      common /vmixr/ diff_cbt(imt,km,jsmw:jemw)
# ifdef bryan_lewis_vertical
      common /vmixr/ Ahv(km)
# endif
#endif
c
#ifdef tcvmix
c
c     variables for mellor-yamada turbulence closure scheme
c
c     diff_cbt_back = background "diff_cbt"
c     visc_cbu_back = background "visc_cbu"
c     tke    = turbulent kinetic energy tendency
c     advq2  = advection of tke
c     hdq2   = horizontal diffusion of tke
c     sprod  = shear production
c     bprod  = buoyancy production
c     dissp  = dissipation of tke
c
# ifdef bryan_lewis_vertical
      common /vmixr/ Ahv(km)
# endif
      common /vmixr/ bprod(imt,km), sprod(imt,km), dissp(imt,km)
      common /vmixr/ tke(imt,km), bprol(imt,km)
      common /vmixr/ advq2(imt,km), fvsq(imt,km)
      common /vmixr/ hdq2(imt,km)
      common /vmixr/ el(imt,km,nslabs), elmax(imt,km)
      common /vmixr/ eeq(imt,km), ffq(imt,km)
      common /vmixr/ gmp1(imt,km), gh(imt,km)
      common /vmixr/ sm(imt,km), sh(imt,km)
#ifdef leq
      common /vmixr/ tkel(imt,km)
      common /vmixr/ advq2l(imt,km)
      common /vmixr/ hdq2l(imt,km)
#endif
      common /vmixr/ a1, b1, a2, b2, cc1, e1, e2, e3
      common /vmixr/ vdqlim, diff_cbt_back, visc_cbu_back
      common /vmixr/ alpha, vk, sq
      common /vmixr/ rhoh2o
      common /vmixr/ small, big, deps
      common /vmixr/ p032, p1, p154, p4, p48
      common /vmixr/ c15, c100
      common /vmixi/ jpr1,ipr1,jpr2,ipr2
c
#endif
