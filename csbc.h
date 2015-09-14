c====================== include file "csbc.h" ==========================
c
c                    surface boundary conditions (S.B.C.)
c
c
c     numosbc = number of ocean S.B.C..
c
c     numasbc = number of atmospheric S.B.C..
c
c     maxsbc  = maximum number of allowable S.B.C. fields.
c
c     numsbc  = total number of S.B.C. this is the sum of the 
c               atmosphere S.B.C. and ocean S.B.C. fields. it must
c               be less than or equal to "maxsbc".
c
c     sbcocn  = S.B.C. fields on the ocean grid.
c               n=1..numosbc are for the incomming ocean S.B.C.
c               n=numosbc+1..numosbc+numasbc are for the outgoing S.B.C.
c
#ifdef restorst
c     dampts  = time scale for damping surface tracers (days) to data
c     dampdz  = ocean level thickness for converting newtonian damping
c               to a surface flux
#endif
c
c     ntspos  = number of time steps per ocean segment
c
c     sbcname  = S.B.C. character string names (eg: "SST", "taux", etc)
c               note that these character strings should be changed with
c               care since they are used to identify fields for special
c               treatment
c
c     dunits  = dimensional units for S.B.C. (eg: "deg C")
c
c     mapsbc  = indices for relating how S.B.C. fields are stored
c               within arrays "sbcocn" and "sbcatm"  
c               eg: m = mapsbc(n) where sbcocn(,,m) is S.B.C. "n"
c
c     restrictions on re-arranging S.B.C. fields and "mapsbc":
c
c     1) all ocean S.B.C. (taux, tauy, etc) should be grouped together. 
c        relative ordering within the group is unimportant but the
c        group must preceed the atmos S.B.C. group.
c
c     2) all atmos S.B.C. (sst, ...) should also be grouped together.
c        relative ordering within the group is unimportant.
c
c     the default arrangement is defined as:
c
c             (these five ocean S.B.C. come from the atmos )
c     mapsbc(1) references taux
c     mapsbc(2) references tauy
c     mapsbc(3) references heat flux
c     mapsbc(4) references salt flux
c      if the atmos supplied fresh water flux then it
c      would have to be converted to ...
c      salt flux  = -(P-E+R) * rho * Sref
c      where...
c      P-E+R represents a precipitation minus evaporation
c            plus runoff rate in cm of water/sec
c      rho is the density if water taken to be 1 g/cm**3, and
c      Sref is a reference salinity in units of grams of salt per
c            gram of water (in units of "parts per part" such
c            as 0.035 ... not "parts per thousand")
c      Depending upon the application of interest, users may wish
c      to set Sref to either be a constant over the entire model
c      domain or for it to be the locally predicted salinity of
c      the uppermost model level (S(1)).  If one desires to ensure
c      that a global average P-E+R flux of zero translates to zero
c      trend in the salt content of the ocean, then a constant Sref
c      is needed.      
c     mapsbc(5) references short wave
c
c             (these four atmos S.B.C. come from the ocean)
c     mapsbc(6) references SST
c     mapsbc(7) references SSS
c     mapsbc(8) references surface u
c     mapsbc(9) references surface v
c
c
c     when re-ordering S.B.C., search for "mapsbc"  
c     (eg: grep mapsbc *.F) to insure the new indexing is consistant
c     with the re-ordering. also make changes in "sbcname", "dunits", 
c     "coabc", and "crits" to maintain consistancy.
c
c     example 1: suppose you want only the first six S.B.C...
c                in "csbc.h", set  numosbc=5, numasbc=1
c                this reduces memory requirements.
c                no other changes are necessary
c
c     example 2: suppose you want to remove the short wave boundary
c                condition (to save more memory) from example 1 . 
c                a) in "csbc.h", set numosbc=4, numasbc=1
c                b) grep "mapsbc(5)" *.F to find and remove references
c                c) grep "mapsbc(6)" *.F to change references to 
c                                        mapsbc(5)
c
#ifdef coupled
c
c     ntspas  = number of time steps per atmosphere segment. if not
c               "coupled" then ntspas = ntspos
c
c     sbcatm  = S.B.C. fields on the atmosphere grid.
c               n=1..numosbc are for the outgoing ocean S.B.C.
c               n=numosbc+1..numosbc+numasbc are for the incomming S.B.C
c
c     abcgx   = longitudes (degs) for the boundary condition fields on
c               the atmosphere grid. (similar to "xt" in the ocean)
c
c     abcgy   = latitudes (degs) for the boundary condition fields on
c               the atmosphere grid. (similar to "yt" in the ocean)
c
c     abcgcs  = cosine of "abcgy"
c
c     abcgdx  = width (degs) of the boundary condition grid boxes
c               in the longitude direction on the atmosphere grid.
c
c     abcgdy  = height (degs) of the boundary condition grid boxes
c               in the latitude direction on the atmosphere grid.
c
c     aland   = mask of (0,1) for (land,ocean) on the atmospheric grid
c
c
c     asegs, asege, define the beginning and end of segments for
c     initialization and averaging purposes in the atmosphere.
c
c     asegs   = true on the 1st time step of an atmosphere segment.
c               otherwise false.
c
c     asege  =  true on the last time step of an atmosphere segment.
c               otherwise false.
c
c     afirst =  true on the 1st atmospheric time step of the run
c     alast  =  true on the last atmospheric time step of the run
c
c
c     isocn, ieocn, jsocn, jeocn locate the ocean grid domain in terms
c     of atmosphere grid indices
c
c     isocn   = index of the first atmos grid point (in longitude) that
c               falls within the ocean domain
c     ieocn   = index of the last atmos grid point (in longitude) that
c               falls within the ocean domain
c     jsocn   = index of the first atmos grid point (in latitude) that
c               falls within the ocean domain
c     jeocn   = index of the last atmos grid point (in latitude) that
c               falls within the ocean domain
c
c
c     iombc   = unit for reading/writing MOM boundary condition restart
c               data
c
c     coabc   = conversion factors for converting S.B.C.
c               from their dimensional units (in the model where they
c               were constructed) to the other model`s units
c               note: if not "coupled" then data is assumed to be 
c               in units required by MOM so... no conversion is done.
c
c     crits   = convergence criteria for forcing S.B.C. into land on the
c               model grid where they were constructed.
c
c     numpas  = number of passes used to extrapolate data into land
c               regions on the model grids where it was constructed.
c
c     sstpre  = prescribed value of "sst" outside of the ocean model
c               domain. this is used only when the ocean model is NOT
c               global in extent (the atmospheric model needs "sst"
c               globally). see "bwidth" below for further explanation.
c
c     bwidth  = when the ocean domain is NOT global, "bwidth" specifies
c               the width of a blending zone (in degrees) surrounding
c               the ocean domain. the purpose of this is to provide a
c               smooth transition between "sst" inside the ocean and
c               the prescribed sst "sstpre" outside of the ocean (since
c               the atmosphere always requires global "sst" as a b.c.).
c               for example:
c               when the ocean domain is a cyclic strip between two
c               latitudes, "sstpre" can be set to a polar temperature
c               (outside the ocean domain) and "bwidth" can be set to
c               extend to within one atmospheric grid point of the pole.
c               when the ocean domain is a basin, "sstpre" should
c               really be changed to a function of latitude and
c               and longitude.
c               setting "bwidth" = 0 implies a global ocean with no
c               blending zone.
c
c
c     bzone   = is defined by setting a non zero "bwidth". this assumes
c               the ocean is NOT global (see "bwidth"). "bzone"
c               defines buffer zones exterior to the ocean domain when
c               the ocean domain is less than global. "bzone" is
c               defined on the atmosphere b.c. grid. it defines
c               three regions as shown below:
c
c                 bzone=0 is the blending zone of width "bwidth"
c                 bzone=1 is the non global ocean domain
c                 bzone=2 is everything external to the blending zone
c
c               Note: these zones are automatically defined by the
c               grid definitions of atmos & ocean models. the only
c               freedom is to choose "bwidth" (described above)
c
c                x--------------------------------------------x
c                |      global atmos domain                   |
c                |        (bzone=2)                           |
c                |    x----------------------------------x    |
c                |    |          blending zone           |    |
c                |    |          (bzone=0)               |    |
c                |    |    x------------------------x    |    |
c                |    |    |                        |    |    |
c                |    |    |       ocean domain     |    |    |
c                |    |    |        (bzone=1)       |    |    |
c                |    |    |                        |    |    |
c                |    |    x------------------------x    |    |
c                |    |                                  |    |
c                |    |                                  |    |
c                |    x----------------------------------x    |
c                |                                            |
c                |                                            |
c                x--------------------------------------------x
c
#endif
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c
      logical asegs, asege, afirst, alast
      character*10 sbcname
      character*15 dunits
c
      parameter (maxsbc=9, numosbc=5, numasbc=2, numsbc=numosbc+numasbc) 
c
#if defined simple_sbc
c
c     there is no "sbcocn" array. its effect is implemented in setvbc.F
c
#else
      common /csbc/ sbcocn(imt,jmt,numsbc) 
#endif
      common /csbc/ dampts(nt), dampdz(nt)
      common /csbc/ coabc(maxsbc), crits(maxsbc)
      common /csbc1/ mapsbc(maxsbc), numpas, iombc, ntspos, ntspas   
      common /csbc2/ sbcname(maxsbc), dunits(maxsbc)
      common /csbc3/ asegs, asege, afirst, alast 
#ifdef coupled
# include "aparam.h"
      integer aland, bzone
      common /csbc4/ aland(imap2,jma), bzone(imap2,jma), isocn, ieocn
      common /csbc4/ jsocn, jeocn
      common /csbc5/ sbcatm(imap2,jma,numsbc)                      
      common /csbc5/ abcgx(imap2), abcgy(jma), abcgdx(imap2) 
      common /csbc5/ abcgdy(jma), abcgcs(jma), sstpre, bwidth 
#endif
c
