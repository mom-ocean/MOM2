c====================== include file "cvbc.h" ==========================
c
c     vertical boundary condition variables:
c
c     smf = surface momentum flux
c          1 => zonal wind stress (dynes/cm**2)
c          2 => meridional wind stress (dynes/cm**2)
c     bmf = bottom momentum flux
c          1 => zonal bottom drag (dynes/cm**2)
c          2 => meridional bottom drag (dynes/cm**2)
c     stf = surface tracer flux
c          1 => surface heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
c                                 (assuming rho*cp = 1 cal/degC/cm**3)
c          2 => surface salt flux (grams of salt/cm**2/sec)
c                                 (assuming rho of water = 1 g/cm**3)
c     btf = bottom tracer flux (for consistency but normally zero!)
c          1 => bottom heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
c                                (assuming rho*cp = 1 cal/degC/cm**3)
c          2 => bottom salt flux (grams of salt/cm**2/sec)
c                                (assuming rho of water = 1 g/cm**3)
c
c     Note:  to convert "stf" and "bmf" salt fluxes to freshwater fluxes
c            one may use the following relationship:
c                salt flux  = -(P-E+R) * rho * Sref
c            where...
c            P-E+R represents a precipitation minus evaporation
c                  plus runoff rate in cm of water/sec
c            rho is the density if water taken to be 1 g/cm**3, and
c            Sref is a reference salinity in units of grams of salt per
c                  gram of water (in units of "parts per part" such
c                  as 0.035 ... not "parts per thousand")
c            Depending upon the application of interest, users may wish
c            to set Sref to either be a constant over the entire model
c            domain or for it to be the locally predicted salinity of
c            the uppermost model level (S(1)).  If one desires to ensure
c            that a global average P-E+R flux of zero translates to zero
c            trend in the salt content of the ocean, then a constant Sref
c            is needed.
c
      common /cvbc/ smf(imt,jsmw:jemw,2),  bmf(imt,jsmw:jemw,2)
      common /cvbc/ stf(imt,jsmw:jemw,nt), btf(imt,jsmw:jemw,nt)
c
