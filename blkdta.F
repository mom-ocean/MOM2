      block data blkdta
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "cprnts.h"
#include "csbc.h"
#include "csnap.h"
#include "emode.h"
#include "index.h"
#include "iounit.h"
#ifdef isopycmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "stab.h"
#include "state.h"
#include "tmngr.h"
#include "calendar.h"
c
c     start data initializations
c
#include "dncoef.h"

c
c     initialize how often diagnostics are to be written
c
      data glenint /-365.25/
      data prxzint /-365.25/
      data extint  /-365.25/
      data vmsfint /-365.25/
      data gyreint /-365.25/
      data dspint  /-365.25/
      data stabint /-365.25/
      data zmbcint /-365.25/
      data trajint /-365.25/
c
      data tsiint    /1.0/
      data cmixint   /-365.25/
      data exconvint /-365.25/
      data tavgint   /-365.25/
      data itavg     /.true./
      data tmbint    /-365.25/
      data itmb      /.true./
      data trmbint   /-365.25/
      data itrmb     /.true./
      data snapint   /-365.25/
      data timavgint /-365.25/
      data xbtint    /-365.25/
c
c     default starting and ending latitudes (degrees) for snapshots.
c     also depth range from surface (centimeters). 
c
      data snapls /-90.0/
      data snaple /90.0/
      data snapde /5000.0e2/
c
c     default starting & ending coordinates for diagnostic matrix prints 
c     (lats and lons in degrees, depths in cm)
c     see "cprnts.h"
c
      data prlat /-60.0, 0.0, 27.0, 55.0, 6*100.0/
      data prslon /180.0, 180.0, 180.0, 180.0, 6*0.0/
      data prelon /250.0, 250.0, 250.0, 250.0, 6*0.0/
      data prsdpt /nlatpr*0.0/
      data predpt /nlatpr*6000.e2/
      data slatxy, elatxy, slonxy, elonxy /-90.0, 90.0, 3.0, 357.0/
c
c     default some control #`s for how diagnostics are to be written
c      (see "iounit.h" for explanation) 
c
      data restrt /.false./
      data iotavg  /-1/
      data iotmb   /-1/
      data ioglen  /-1/
      data iotrmb  /-1/
      data iogyre  /-1/
      data iovmsf  /-1/
      data ioprxz  /-1/
      data ioext   /-1/
      data iodsp   /-1/
      data iotsi   /-1/
      data iozmbc  /-1/
      data ioxbt   /-1/
c
c
#if defined firfil || defined fourfil
c
c     set latitudes used in filtering of tracer and velocity fields
c     see "index.h"
c
      data rjfrst /-81.0/
      data rjft0, rjft1, rjft2 /-51.0, -57.0, 57.0/
      data rjfu0, rjfu1, rjfu2 /-54.0, -60.0, 60.0/
c
#endif
#ifdef firfil
c     "numflt" specifies the number of filter applications to use
c              for the tracers
c     "numflu" specifies the number of filter applications to use
c              for velocities and vorticity
c
      data numflt /jmtfil*1/
      data numflu /jmtfil*1/
c
#endif
c
      end
