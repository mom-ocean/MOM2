c====================== include file "diaga.h" =========================
c
c     variables used for computing diagnostics:
c
#ifdef save_convection
c
c     convect0  = temperature before explicit convection
c     convect1  = time rate of change of temperature due to 
c                 explicit convection. set to epsln over land
c                 for use in identifying land cells
c
      common /exconv/ excnv0(imt,km, jsmw:jemw)
      common /exconv/ excnv1(imt,km,jsmw:jemw)
#endif
#ifdef save_mixing_coeff
c
c     ce = ceoff on east face of cell (1 is for momentum, 2 for tracers)
c     cn = ceoff on north face of cell (1 is for momentum, 2 for tracers)
c     cb = ceoff on bottom face of cell(1 is for momentum,2 for tracers)
c
      common /mixcoef/ ce(imt,km,jsmw:jemw,2)
      common /mixcoef/ cn(imt,km,jsmw:jemw,2)
      common /mixcoef/ cb(imt,km,jsmw:jemw,2)
#endif
