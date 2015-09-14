c======================= include file "hmixc.h" ========================
c
c                    horizontal mixing coefficients
c
c     visc_cnu = viscosity coeff for northern face of "u" cell
c     visc_ceu = viscosity coeff for eastern face of "u" cell
c     diff_cnt = diffusion coeff for northern face of "T" cell
c     diff_cet = diffusion coeff for eastern face of "T" cell
c     
c     am     = constant lateral viscosity coeff for momentum
c     ah     = constant lateral diffusion coeff for tracers
c     am3    = viscosity coeff for metric term on "u" cell
c     am4    = another viscosity coeff for metric term on "u" cell
c     ambi   = constant lateral biharmonic viscosity coeff for momentum
c     ahbi   = constant lateral biharmonic diffusion coeff for tracers
c
c     author:   r.c.pacanowski     e-mail  rcp@gfdl.gov
c=======================================================================
c
#if defined consthmix
      common /diffus/ am, ambi, am3(jmt), am4(jmt,2)
      common /diffus/ ah, ahbi
      common /diffus/ visc_ceu, visc_cnu
      common /diffus/ amc_north(jmt), amc_south(jmt)
c
# if defined bryan_lewis_horizontal
c
c     bryan_lewis mixing case
c
      common /diffus/ Ahh(km)
      common /diffus/ diff_cnt(km), diff_cet(km) 
      common /diffus/ ahc_north(jmt,km), ahc_south(jmt,km) 
# else
      common /diffus/ diff_cnt, diff_cet
      common /diffus/ ahc_north(jmt), ahc_south(jmt) 
# endif
#endif
c
#ifdef smagnlmix
c
c     non-linear horizontal viscosity after Smagorinsky 1963,
c     as described in Rosati & Miyakoda (jpo,vol 18,#11,1988)
c     see Smagorinsky 1963, Mon Wea Rev, 91, 99-164. 
c     Also see Deardorff 1973 J. Fluid Eng. Sep., 429-438.
c
c     strain = tension(1) and shearing(2) rates of strain
c     smag_metric  = metric term
c     diff_c_back = background diffusion coeff for "t" cell (cm**2/sec)
c     
      common /diffus/ strain(imt,km,1:jemw,2)
      common /diffus/ am_lambda(imt,km,1:jemw), am_phi(imt,km,1:jemw)
      common /diffus/ smag_metric(imt,km,jsmw:jemw)
      common /diffus/ diff_c_back
      common /diffus/ visc_ceu(imt,km,jsmw:jemw)
      common /diffus/ visc_cnu(imt,km,1:jemw)
      common /diffus/ diff_cet(imt,km,jsmw:jemw)
      common /diffus/ diff_cnt(imt,km,1:jemw)
#endif
c
#ifdef held_larichev
c
c     variables for held_larichev diffusion
c
c     hl_u     = Held/Larichev diffusion coefficient (cm**2/sec)
c                defined over "u" cells
c     hl_n     = avg_x(hl_u) defined at northern face of "T" cells
c     hl_e     = avg_y(hl_u) defined at eastern face of "T" cells
c     hl_b     = avg_x(avg_y(hl_u))
c     hl_depth = depth of integration (cm)
c     hl_back  = background mixing coeff (cm**2/sec)
c     hl_max   = max allowed mixing coeff (cm**2/sec)
c     droz     = vertical difference of rho
c
      common /diffus/ hl_depth, hl_back, hl_max
      common /diffus/ hl_u(imt,1:jemw)
      common /diffus/ hl_n(imt,1:jemw), hl_e(imt,jsmw:jemw)
      common /diffus/ hl_b(imt,jsmw:jemw)
      common /diffus/ droz(imt,km,jmw), rich_inv(imt,km,1:jemw)
#endif
