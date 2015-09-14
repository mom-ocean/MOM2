      subroutine state (t, s, rho, js, je, is, ie)
c
c=======================================================================
c
c     state computes normalized densities by using a 3rd
c     order polynomial fit to the knudsen formula, for each level
c     subtract normalizing constants from temperature and salinity
c     and compute polynomial approximation of knudsen density.
c
c     note.. for precision purposes, there is a depth dependent
c     constant subtracted from the density returned by this routine.
c     so... this routine should be used only for horizontal gradients
c     of density.
c
c     inputs:
c
c     t  = the input row of temperatures (potential deg C)
c     s  = the input row of salinities (units: (ppt-35)/1000)
c     js = starting row for computing density within the MW
c     je = ending row for computing density within the MW
c     is = starting longitude index for computing density within the MW
c     ie = ending longitude index for computing density within the MW
c
c     output:
c
c     rho = normalized densities
c     These densities are in cgs units(g/cm3) and represent 
c     the in situ density at a level minus a depth dependent
c     normalization. The complete in situ density is given by 
c     rho_complete(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k) 
c                           + rho_norm(k)*10-3,
c     where rho_norm(k) are the depth dependent normalization densities  
c     [in sigma units (density-1)*1000] given at the bottom of dncoef.h
c
c=======================================================================
c
c
#include "param.h"
#include "state.h"
      dimension t(imt,km,jmw), s(imt,km,jmw), rho(imt,km,jsmw:jmw)
#include "dens.h"
c
#ifdef timing
      call tic ('state', 'state')
#endif
      do j=js,je
        do k=1,km
          do i=is,ie
            rho(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k)
          enddo
        enddo
      enddo
#ifdef timing
      call toc ('state', 'state')
#endif
c
      return
c
c
c
      entry statec (t, s, rho, js, je, is, ie, ind)
c
c
c=======================================================================
c
c     statec computes, for one row, the normalized densities by using
c     a 3rd order polynomial fit to the knudsen formula. For
c     purposes of checking vertical stability between adjacent
c     levels, the reference depth for pressure dependence in
c     the knudsen formula must be held constant. that level is
c     determined by "ind".
c
c     inputs:
c
c     t   = the input row of temperatures (potential deg C)
c     s   = the input row of salinities (units: (ppt-35)/1000)
c     js  = starting row for computing density within the MW
c     je  = ending row for computing density within the MW
c     is  = starting longitude index for computing density within the MW
c     ie  = ending longitude index for computing density within the MW
c     ind = 1 for comparing levels 1 to 2, 3 to 4, etc.
c           (coefficients for the lower of the 2 levels are used)
c           2 for comparing levels 2 to 3, 4 to 5, etc.
c           (coefficients for the lower of the 2 levels are used)
c
c     output:
c
c     rho = normalized densities
c     These densities are in cgs units(g/cm3) and represent 
c     the in situ density at a level minus a depth dependent
c     normalization. The complete in situ density is given by 
c     rho_complete(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k) 
c                           + rho_norm(k)*10-3,
c     where rho_norm(k) are the depth dependent normalization densities  
c     [in sigma units (density-1)*1000] given at the bottom of dncoef.h
c
c=======================================================================
c
c
#ifdef timing
      call tic ('state', 'statec')
#endif
      if (ind .lt. 1 .or. ind .gt. 2) then
        write (stderr,99) ind
        stop '=>statec'
      endif
c
      do j=js,je
        do l=1,km,2
          if (ind .eq. 1) then
            k = min(l+1,km)
          else
            k = l
          endif
          do i=is,ie
            rho(i,l,j) = dens (t(i,l,j)-to(k), s(i,l,j)-so(k), k)
          enddo
        enddo
      enddo
c
      do j=js,je
        do l=2,km,2
          if (ind .eq. 1) then
            k = l
          else
            k = min(l+1,km)
          endif
          do i=is,ie
            rho(i,l,j) = dens (t(i,l,j)-to(k), s(i,l,j)-so(k), k)
          enddo
        enddo
      enddo
#ifdef timing
      call toc ('state', 'statec')
#endif
      return
   99 format(/' error => bad "ind" in statec: ind =',i10)
      end
