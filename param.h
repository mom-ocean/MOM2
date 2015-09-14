c====================== include file "param.h" =========================
c
c     main parameter file which sets ocean characteristics:
c
c     ncrows = number of calculated rows within the MW. (the remaining 
c              rows are buffer rows).
c     maxmw  = maximum number of times data is loaded from disk into the
c              MW to solve latitude rows 2 through jmt-1
c     lseg   = maximum number of longitudinal stream function segments
c     numjpt = number of pointers for internal mode velocities
c              (also used for del**2 fields with biharmonic option)
c     nvar   = number of prognostic variables
c     nvarbh = number of prognostic variables using biharmonic mixing
c     nhreg = number of regions in the horizontal used for averaging
c             tracers.
c     nvreg = number of regions in the vertical used for term balance
c             calculations. note "nvreg" is not used for tracer averages
c     numreg= total number of regions ( = product of nhreg & nvreg)
c             used for term balance calculations
c     nlatpr = maximum number of latitudes for matrix printouts
c              on diagnostic time steps
c
c
#include "size.h"
#ifdef biharmonic
      parameter (ncrows = jmw - 3)
#else
      parameter (ncrows = jmw - 2)
#endif
      parameter (nmwlat=(jmt-2)/ncrows)                                  
      parameter (maxmw = nmwlat + (jmt-3)/(ncrows*nmwlat))               
c
#ifdef isopycmix
c
c     nrpl = number of reference presure levels used in isopycnal 
c            tracer mixing
c
      parameter (nrpl=5)
#endif
      parameter (lseg=5, nlatpr=10)
      parameter (nhreg=5, nvreg=2, numreg=nhreg*nvreg)
#ifdef tcvmix
      parameter (nvar=nt+2+4)
# ifdef leq
      parameter (nvar=nt+2+4+1)
# endif
#else
      parameter (nvar=nt+2)
#endif
      parameter (imtm1=imt-1, kmm1=km-1)                                
      parameter (imtp1=imt+1, imtm2=imt-2
     &,          jmtp1=jmt+1, jmtm1=jmt-1, jmtm2=jmt-2
#ifdef symmetry
     &,          jscan=jmtm2+1
#else
     &,          jscan=jmtm2
#endif
     &,          kmp1=km+1, kmp2=km+2
     &,          imtkm=imt*km, nwds=imt*jmt, nkflds=2
     &,          nslab=imt*nvar*km, ntmin2=nt+1/nt)
c
c     add parameter constants
c
#include "pconst.h"
#include "stdunits.h"
c
