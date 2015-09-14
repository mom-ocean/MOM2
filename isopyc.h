c======================== include file "isopyc.h" ======================
#ifdef isopycmix
c
c     isopycnal diffusion variables:
c
c     ahisop = isopycnal tracer diffusivity (cm**2/sec) 
c     athkdf = isopycnal thickness diffusivity (cm**2/sec)
c     dptlim = depth limits for the reference pressure levels (in cm). 
c              the mid point of the two neighboring "dptlim" elements is
c              used as the reference pressure (depth) for that interval,
c              i.e.,
c
c             reference pressure level        reference depth
c              ----------------------- -----------------------------
c                          1           0.5*(dptlim(1)+dptlim(2))
c                          2           0.5*(dptlim(2)+dptlim(3))
c                          .           .
c                        nrpl          0.5*(dptlim(nrpl)+dptlim(nrpl+1))
c
c              REMARK: the first and the last elements of "dptlim" must
c                      be the depth at the top (0m) and the maximum bottom
c                      depth, respectively. also, the elements of "dptlim"
c                      must be in increasing order.
c     e      = scratch array 
c     K1     = (1,3) component of the isopycnal mixing tensor computed
c              at the center of the eastern face of the "t" grid cell
c     K2     = (2,3) component of the isopycnal mixing tensor computed
c              at the center of the northern face of the "t" grid cell
c     K3     = (3,.) components of the isopycnal mixing tensor 
c              computed at the top face of the "t" grid cell
c               (,,1) --> (3,1) component
c               (,,2) --> (3,2) component
c               (,,3) --> (3,3) component 
c              REMARK: (1,1) and (2,2) components of the tensor are
c                      assumed to be unity. also, (1,2) and (2,1)
c                      components are assumed to be zero.
c     fzisop = function containing the vertical variation of the isopycnal
c              diffusion coefficient. "fzisop" multiplies "ahisop".
c     kisrpl = isopycnal reference pressure levels for the "t" grid 
c              point levels computed based on the depth (pressure) at the
c              "t" grid points
c     krplin = indices indicating the location of the reference pressure 
c              depth in the 20-level table of polynomial expansion
c              variables
c     slmxr  = reciprocal of maximum slope of isopycnals allowed in mixing
c              scheme to prevent excessively large vertical mixing that
c              could create numerical instabilities. furthermore, the
c              form of the isopycnal diffusion tensor incorporates the
c              assumption that horizontal density gradients are much
c              smaller than vertical gradients. "slmxr" is also used to
c              satisfy this assumption. a value of 100 for "slmxr"
c              translates to a slope of 1:100.
c
c     rhoi    = potential density at "t" cell centers
c
# ifdef isopycmixspatialvar
c     dciso1 = isopycnal tracer diffusivity coefficients modified based
c              on the slopes of the isopycnal surfaces on the east face
c              of "T" cells.
c     dciso2 = isopycnal tracer diffusivity coefficients modified based
c              on the slopes of the isopycnal surfaces on the north face
c              of "T" cells.
c     dslope = half length of the interval in which "ahisop" changes
c              with a steep slope from about 0.9*"ahisop" to about
c              0.1*"ahisop"
c     slopec = slope at which "ahisop" is equal to half of its original
c              value
c
c              REMARK: 0 <= "slopec", "dslope" <= 1/"slmxr".
c              REMARK: because the vertical gradient of density must be
c                      less than zero for isopycnal mixing, 1/"slmxr" is
c                      actually a negative maximum slope. this fact is
c                      taken into account in "isopyc.F". consequently,
c                      "slmxr" must be a positive number in "blkdta.F".
# endif
# ifdef gent_mcwilliams
c     adv_vetiso = zonal isopycnal mixing velocity computed at the 
c                  center of the eastern face of the "t" cells
c     adv_vntiso = meridional isopycnal mixing velocity computed at
c                  the center of the northern face of the "t" cells
c                  (Note: this includes the cosine as in "adv_vnt")
c     adv_vbtiso = vertical isopycnal mixing velocity computed at the
c                  center of the top face of the "t" cells
c     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
c                  the bottom face of the "t" cells
# endif
c
      real K1, K2, K3
      common /cisop/ rhoi(imt,km,1:jmw,nrpl) 
      common /cisop/ e(imt,kmp1,1:jemw,3)
      common /cisop/ K1(imt,km,jsmw:jemw,3:3), K2(imt,km,1:jemw,3:3) 
      common /cisop/ K3(imt,km,jsmw:jemw,1:3)   
      common /cisop/ ahisop, athkdf, dptlim(nrpl+1), fzisop(km), slmxr 
      common /cisopi/ kisrpl(km), krplin(nrpl)
# ifdef gent_mcwilliams
      common /cisop/ adv_vetiso(imt,km,jsmw:jemw)
      common /cisop/ adv_vntiso(imt,km,1:jemw)
      common /cisop/ adv_vbtiso(imt,0:km,jsmw:jemw)   
      common /cisop/ adv_fbiso(imt,0:km,jsmw:jemw)   
# endif
# ifdef isopycmixspatialvar
      common /cisop/ dciso1(imt,km,jsmw:jemw), dciso2(imt,km,1:jemw)
      common /cisop/ dslope, slopec
# endif
c
c     author:    Gokhan Danabasoglu (previous MOM version)
c                                    gokhan@isis.cgd.ucar.edu
c                r.c.pacanowski  (modifications for MOM 2) rcp@gfdl.gov
c                A.Rosati        (modifications for MOM 2) ar@gfdl.gov
#endif
