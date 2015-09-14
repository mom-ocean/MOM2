      subroutine setvbc (joff, js, je, is, ie)
c
c=======================================================================
c     set momentum and tracer vertical boundary conditions 
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#include "mw.h"
#ifdef simple_sbc
      dimension star(4)
#endif
c
#ifdef timing
      call tic ('setvbc', 'setvbc')
#endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c----------------------------------------------------------------------
c       set no flux condition for all tracers at surface & bottom.
c----------------------------------------------------------------------
c
      do n=1,nt
        do j=js,je
          do i=istrt,iend
            stf(i,j,n) = c0
            btf(i,j,n) = c0
          enddo
        enddo
      enddo
c
#if defined minimize_sbc_memory
c
c----------------------------------------------------------------------
c     construct "sbcocn" for latitudes corresponding to memory window
c     rows js through je.
c----------------------------------------------------------------------
c
      call atmos (joff, js, je)
#endif
c
c----------------------------------------------------------------------
c       apply surface tracer and momentum fluxes from the atmosphere
c       code is for 2 tracer and 2 momentum fluxes.
c----------------------------------------------------------------------
c
      do n=1,2
        if (n .ne. 1 .and. n .ne. 2) then
          write (stdout,*) '=>Error: setvbc n =',n
	  stop '=>setvbc'
        endif
c
c       mapsbc(1) = taux
c       mapsbc(2) = tauy
c       mapsbc(3) = heat flux
c       mapsbc(4) = fresh water flux
c
        do j=js,je
          jrow = j + joff
# ifdef simple_sbc
c
          call bcest (yt(jrow), yu(jrow), star)
c
c         star(1) is taux in dynes/cm**2
c         star(2) is tauy in dynes/cm**2
c         star(3) is SST in deg C
c         star(4) is SSS in parts per part
c
c         convert surface salinity to model units
c
          star(4) = (star(4) - 35.0) * 0.001
# endif
          itflx = mapsbc(2+n)
          if (itflx .ne. 0) then
            do i=istrt,iend
# ifdef simple_sbc
#  ifdef restorst
              damp       = dampdz(n)*secday/dampts(n)
              stf(i,j,n) = damp*(star(n+2) - t(i,1,j,n,taum1))
#  else
              stf(i,j,n) = c0
#  endif
# else
              stf(i,j,n) = sbcocn(i,jrow,itflx)
# endif
            enddo
          endif
          itau = mapsbc(n)
          if (itau .ne. 0) then
            do i=istrt,iend
# ifdef simple_sbc
              smf(i,j,n) = star(n)
# else
              smf(i,j,n) = sbcocn(i,jrow,itau)
# endif
            enddo
          else
            do i=istrt,iend
              smf(i,j,n) = c0
            enddo
          endif
        enddo
      enddo
c
c----------------------------------------------------------------------
c       set bottom drag
c----------------------------------------------------------------------
c
      do n=1,2
        if (cdbot .eq. c0) then
          do j=js,je
            do i=istrt,iend
	      bmf(i,j,n) = c0
            enddo
          enddo
        else
          do j=js,je
            jrow = j + joff
            do i=istrt,iend
              kz = kmu(i,jrow)
              if (kz .ne. 0) then
                uvmag    = sqrt(u(i,kz,j,1,taum1)**2 + 
     &                          u(i,kz,j,2,taum1)**2)
                bmf(i,j,n) = cdbot*u(i,kz,j,n,taum1)*uvmag
              else
                bmf(i,j,n) = c0
              endif
            enddo
          enddo
        endif
c
c----------------------------------------------------------------------
c     apply zonal boundary conditions
c----------------------------------------------------------------------
c
        call setbcx (stf(1,jsmw,n), imt, jemw-jsmw+1)
        call setbcx (btf(1,jsmw,n), imt, jemw-jsmw+1)
        call setbcx (smf(1,jsmw,n), imt, jemw-jsmw+1)
        call setbcx (bmf(1,jsmw,n), imt, jemw-jsmw+1)
      enddo
c
# ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In setvbc: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
     & 
# endif
#ifdef timing
      call toc ('setvbc', 'setvbc')
# endif
c
      return
      end
