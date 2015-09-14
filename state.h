c====================== include file "state.h" =========================
c
c     variables for equation of state
c
c     to = reference temperture for level
c     to = reference salinity for level
c     c  = polynomial coefficients for equation of state
c     cksumzt = checksum of depths used to calculate density coeffs
c
      common /stater/ to(km), so(km), c(km,9), cksumzt
