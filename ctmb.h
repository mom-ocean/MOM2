c====================== include file "ctmb.h" ==========================
c
c                      meridional tracer balance
c
c     to calculate the meridional tracer balance, each term in the
c     tracer equation is integrated in depth,longitude and time for each
c     latitude. the time integration is over "tmbint"(see switch.h) days
c     for the calculation to make sense, the longitudinal
c     integration must go around the world (in the case of an open
c     latitude {eg: 60 deg S}) or between land masses (in the case of
c     a closed latitude {eg: 30 deg S from south america to africa})
c     the domain may be divided into arbitrary basins (atlantic, indian
c     pacific) but for certain latitudes, results from basins may have
c     to be combined to satisfy the above conditions. 
c
c
c     ntmbb    = number of basins over which the tracer meridional
c                balance will be calculated. the test case assumes 1
c                basin as defined in "setocn.F"
c     tstor    = average latitudinal storage of tracer
c     tdiv     = average latitudinal divergence of tracer
c     tflux    = average latitudinal surface flux of tracer
c     tdif     = average latitudinal diffusion of tracer
c     tsorc    = average latitudinal additional source of tracer
c     smdvol   = ocean area at latitude
c     msktmb   = mask for basin numbers (1..ntmbb. 0 is reserved for
c                sum of all basins)
c     numtmb   = number of time steps over which the terms have been
c                accumulated
      parameter (ntmbb=1)
      common /tmbr/ tstor(jmt,nt,0:ntmbb), tdiv(jmt,nt,0:ntmbb)
      common /tmbr/ tflux(jmt,nt,0:ntmbb), tdif(jmt,nt,0:ntmbb)
      common /tmbr/ tsorc(jmt,nt,0:ntmbb), smdvol(jmt,0:ntmbb)
      common /tmbi/ msktmb(imt,jmt), numtmb
c
