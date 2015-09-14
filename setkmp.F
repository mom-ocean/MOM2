      subroutine setkmp (alat1, slon1, elon1, alat2, slon2, elon2, num)
c
c-----------------------------------------------------------------------
c     set the topography mask "kmt(i,j)" = "num" within the area of the
c     parallelogram bounded by vertices:
c     (alat1,slon1), (alat1,elon1), (alat2,slon1), & (alat2,elon2)
c
c-----------------------------------------------------------------------
c
#include "param.h"
#include "coord.h"
#include "levind.h"
c
c     convert the four vertices into model indices
c     (js,is1), (js,ie1), (je,is2), (je,ie2)
c
      j1 = indp (alat1, yt, jmt)
      j2 = indp (alat2, yt, jmt)
      js = min (j1,j2)
      je = max (j1,j2)
c
      i1  = indp (slon1, xt, imt)
      i2  = indp (elon1, xt, imt)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)
c
      i1  = indp (slon2, xt, imt)
      i2  = indp (elon2, xt, imt)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)
c
      is = is1
      ie = ie1
c
c     fill in the area bounded by (js,is1), (js,ie1), (je,is2), (je,ie2)
c
      if (js .eq. je) then
        rdj = c1
      else
        rdj = c1/(je-js)
      endif
      do 100 j=js,je
        do 90 i=is,ie
          kmt(i,j) = num
90      continue
        is = nint(rdj*((j-js)*is2 + (je-j)*is1))
        ie = nint(rdj*((j-js)*ie2 + (je-j)*ie1))
100   continue
      return
      end
