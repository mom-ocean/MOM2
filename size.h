c====================== include file "size.h" ==========================
c                                                                       
c                        grid size parameters             
c                                                                       
c     imt    = number of grid points in longitudinal direction          
c     jmt    = number of grid points (rows) in latitudinal direction           
c     km     = number of grid points in the vertical direction          
c     jmw    = number of latitude rows in the MW (memory window).       
c              min jmw = 3 (4 with the biharmonic mixing option)
c              max jmw = jmt (all latitude rows in the MW. none on disk)
c     jsmw   = 1st calculated row within the MW     
c     jemw   = last calculated row within the MW       
c     nt     = number of tracers (temperature, salinity, ...)
c     mnisle = maximum number of islands (unconnected land masses)
c     maxipp = maximum number of island perimeter points
c
c                                                                       
      parameter (imt=  92, jmt=  61, km= 15)
#ifdef biharmonic
      parameter (jmw=4, jsmw=2, jemw=jmw-1)
#else
      parameter (jmw=3, jsmw=2, jemw=jmw-1)
#endif
      parameter (nt=2)                                           
      parameter (mnisle=50, maxipp=5000)
