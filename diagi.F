      subroutine diagi
c
c=======================================================================
c     initialize diagnostics quantities
c=======================================================================
c
#include "param.h"
#include "stab.h"
#ifdef tracer_averages
# include "ctmb.h"
#endif
#include "coord.h"
#include "ctavg.h"
#include "diag.h"
#include "iounit.h"
#include "switch.h"
#include "tmngr.h"
c
#ifdef timing
      call tic ('diagnostic', 'initialization')
#endif
#ifdef stability_tests
      if (stabts) then
        numcfl = 0
        cflup  = c0
        cflvp  = c0
        cflwtp = c0
        cflwup = c0
        reynx  = c0
        reyny  = c0
        reynz  = c0
        peclx  = c0
        pecly  = c0
        peclz  = c0
        call getunit (iostab, 'iostab', 'fsr')
	rewind iostab
	endfile iostab
	rewind iostab
        call relunit (iostab)
      endif
#endif
c
#ifdef energy_analysis
      if (glents) then
        do jrow=1,jmt
	  do k=1,km
	    wtlev(k,jrow) = c0
	    wulev(k,jrow) = c0
	  enddo
          do k=0,km
            buoy(k,jrow) = c0
	  enddo
          do ll=1,8
            do k=0,km
              engint(k,ll,jrow) = c0
            enddo
	  enddo
          do ll=1,8
            engext(ll,jrow) = c0
	  enddo
	  tcerr(jrow) = c0
	  ucerr(jrow) = c0
	  itcerr(jrow) = 1
	  jtcerr(jrow) = 1
	  ktcerr(jrow) = 1
	  iucerr(jrow) = 1
	  jucerr(jrow) = 1
	  kucerr(jrow) = 1
	  wtbot(jrow)  = c0
	  wubot(jrow)  = c0
	  iwtbot(jrow) = 1
	  jwtbot(jrow) = 1
	  kwtbot(jrow) = 1
	  iwubot(jrow) = 1
	  jwubot(jrow) = 1
	  kwubot(jrow) = 1
        enddo
      endif
#endif
c
#ifdef term_balances
      if (trmbts) then
        do n=1,nt
          ustf(n,1) = ' unknown  '
          ustf(n,2) = ' unknown units '
          if (n .eq. 1) then
            ustf(n,1) = ' stf(1)   = '
            ustf(n,2) = ' cal/cm**2/sec '
          endif
          if (n .eq. 2) then
            ustf(n,1) = ' stf(2)   = '
            ustf(n,2) = '     cm/sec '
          endif
	enddo
c
        do nreg=0,numreg
          if (nreg .gt. 0) then
            avgw(nreg) = c0
            do ll=1,17
              do k=0,km
                termbm(k,ll,1,nreg) = c0
                termbm(k,ll,2,nreg) = c0
              enddo
            enddo
          endif
          do n=1,nt
            do ll=1,15
              do k=0,km
                termbt(k,ll,n,nreg) = c0
              enddo
            enddo
          enddo
	enddo
c
        do nreg=0,nhreg
          smflx(1,nreg) = c0
          smflx(2,nreg) = c0
          do n=1,nt
            stflx(n,nreg) = c0
            asst(n,nreg)  = c0
          enddo
        enddo
      endif
#endif
c
#ifdef gyre_components
      if (gyrets) then
        do jrow=1,jmt
          do m=1,ntmin2
            do ll=1,8
              ttn(ll,jrow,m) = c0
            enddo
          enddo
	enddo
c
        do n=0,nhreg
          do m=1,nt
            do jrow=1,jmt
#if defined isopycmix && defined gent_mcwilliams
              do ll=6,9
#else
              do ll=6,8
#endif
                ttn2(ll,jrow,m,n) = c0
              enddo
            enddo
          enddo
        enddo
      endif
#endif
c
#ifdef meridional_overturning
      if (vmsfts) then
        do jrow=1,jmt
          do k=1,km
            vsf(jrow,k) = c0
          enddo
	enddo
      endif
#endif
c
#ifdef tracer_averages
      if (tavgts) then
        do n=1,nt
          do mask=1,nhreg
            sumbf(mask,n) = c0
            do k=1,km
              sumbk(mask,k,n) = c0
              avgbk(mask,k,n) = c0
            enddo
          enddo
	enddo
      endif
#endif
c
#ifdef time_step_monitor
      if (tsits) then
        do jrow=1,jmt
          do k=0,km
            ektot(k,jrow) = c0
            do n=1,nt
              dtabs(k,n,jrow)  = c0
              tbar(k,n,jrow)   = c0
              travar(k,n,jrow) = c0
            enddo
          enddo
	enddo
      endif
#endif
c
#ifdef show_zonal_mean_of_sbc
      if (zmbcts) then
        do jrow=1,jmt
	  zmau(jrow) = c0
	  zmat(jrow) = c0
          do n=1,2
            zmsmf(jrow,n) = c0
            zmsm(jrow,n)  = c0
	  enddo
          do n=1,nt
            zmstf(jrow,n) = c0
            zmst(jrow,n)  = c0
          enddo
	enddo
      endif
#endif
#ifdef timing
      call toc ('diagnostic', 'initialization')
#endif
      return
      end
