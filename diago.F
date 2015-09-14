      subroutine diago
c
c=======================================================================
c
c     write out diagnostic output
c
c=======================================================================
c
#include "param.h"
#include "ctmb.h"
#include "coord.h"
#include "cprnts.h"
#include "cregin.h"
#ifdef tracer_averages
# include "ctavg.h"
#endif
#include "diag.h"
#include "docnam.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#include "stab.h"
#include "switch.h"
#include "tmngr.h"
#include "vmixc.h"
c#define debug_extmode
#if defined stream_function && defined debug_extmode
      dimension ext(imt,jmt,2)
#endif
#if defined diagnostic_surf_height && defined stream_function
# include "mw.h"
      character*8 bc_symm
c
c     psgrad= gradient of sea surface pressure (1,2) = (u,v) component
c     dsp   = diagnostic sea surface pressure converted to height
c     divf  = divergence of sea surface pressure gradients
c
      dimension psgrad(imt,jmt,2), divf(imt,jmt)
      dimension dsp(imt,jmt)
      save divf, dsp, numdsp
c
c-----------------------------------------------------------------------
c     reconstruct the diagnostic surface pressure
c-----------------------------------------------------------------------
c
# ifdef timing
      call tic ('diagnostic', 'diagnostic surf height')
# endif
c
c     initialize the surface pressure fields and averaging counter
c
      if (first) then
        numdsp = 0
        do jrow=1,jmt
          do i=1,imt
            psgrad(i,jrow,1)  = c0
            psgrad(i,jrow,2)  = c0
	    divf(i,jrow) = c0
	    dsp(i,jrow)  = c0
	    res(i,jrow)  = c0
          enddo
        enddo
      endif
c
c     increment the counter for averaging the divergence
c
      numdsp = numdsp + 1
c
c     construct surface pressure gradients
c     (note: uext & vext are returned but not needed)
c
      call calc_psgrad (psgrad, uext, vext, 1, jmt, 1, imt)
c
      do i=1,imt
        psgrad(i,jmt,1) = c0
        psgrad(i,jmt,2) = c0
      enddo
      call setbcx (psgrad(1,1,1), imt, jmt)
      call setbcx (psgrad(1,1,2), imt, jmt)
c
c     compute the divergence of the sea surface pressure gradients
c
      call spforc (psgrad, dxu, dyu, csu, h, res)
c
c     accumulate the divergence
c
      do jrow=1,jmt
        do i=1,imt
	  divf(i,jrow) = divf(i,jrow) + res(i,jrow)
	enddo
      enddo
      if (dspts) then
c
        write (stdout,'(///,40x,a,/)') 
     & 'D I A G O N O S T I C    S U R F A C E    H E I G H T'
c
c       average the divergence, zero the guess
c
        rnum = c1/numdsp
        do jrow=1,jmt
          do i=1,imt
	    divf(i,jrow) = rnum*divf(i,jrow)
	    dsp(i,jrow)  = c0
	  enddo
        enddo       
c
c       initialize coefficients for the conjugate gradient solver
c
        call spc9pt (dxu, dyu, csu, h, cf)
        npt    = 9
	type   = 'surfpres'
	bc_symm = 't even'
c
c       number of islands must be zero for surface pressure
c       dspcrt = tolerence of solution for dsp
c
	nislsp = 0
	dspcrt = tolrsp
c
# ifdef timing
        call tic ('congrad', 'diagnostic surf press (diago)')
# endif
c
        call congr (npt, type, bc_symm, dsp, dsp, divf, res 
     &,           cf
     &,           mxscan, mdscan, dspcrt
     &,           imask, iperm, jperm, iofs, nislsp, nippts
     &,           converged, esterr)
c
        write (stdout,'(/a,i5,a,e14.7)') 
     &  '=> diagnostic surface pressure scans=',mdscan,' estimated err='
     &, esterr
c
        if (.not.converged) then
	  write (stdout,'(/,a,/)')
     & '=> Warning: convergence not reached for diag surface pressure'
	  converged = .true.
	endif
c
# ifdef timing
        call toc ('congrad', 'diagnostic surf press (diago)')
# endif
c
c       subtract the null space from the surface pressure
c
        call checkerboard (dsp, map)
        call border (dsp, 't even')
c
c       remove mean from the surface pressure
c
	call zero_level (dsp, 'surf press', map, dxt, dyt, cst)
        call border (dsp, 't even')
c
c       convert surface pressure to surface elevation and write out
c
        rgrav = c1/grav
        do jrow=1,jmt
          do i=1,imt
	    dsp(i,jrow) = rgrav*dsp(i,jrow)
	  enddo
        enddo
c
        if (iodsp .eq. stdout .or. iodsp .lt. 0) then
          is  = indp (slonxy, xt, imt)
          ie  = indp (elonxy, xt, imt)
          js  = indp (slatxy, yt, jmt)
          je  = indp (elatxy, yt, jmt)
          scl = c1
	  write (stdout,'(/1x,a,i5,a)')
     &   'Diagnostic surf height is averaged over ',numdsp,' time steps'
          write (stdout,7100) 'Diagnostic surface height (cm)'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
c
          call matrix (dsp(1,1), imt, is, ie, -js, -je, scl)
c
7100      format(1x,a30,1x,'ts=',i7
     &,    ', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     &,    ', scaling=',1pg10.3)
c
        endif
        if (iodsp .ne. stdout .or. iodsp .lt. 0) then
c
	  call getunit (io, 'diag_surf.dta','u s a ieee')
c
          write (stdout,'(a,a,f8.3,a,a,i10,a)')
     &    ' => Diagnostic surface height (averaged '
     &,   'over ',extint,' days) written '
     &,  'unformatted to file diag_surf.dta on ts=', itt, ' ',stamp
          reltim = prelyr
c
          iotext = 'read (iodsp) imt, jmt, reltim, extint, xt, yt'
          write (io) pstamp, iotext, expnam
          write (io) imt, jmt, reltim, extint, xt, yt
c
          iotext = 'read (iodsp) ((dsp(i,j),i=1,imt),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, dsp, imt*jmt)
c
          call relunit (io)
        endif
c
c       reset counter for next average. zero the accumulated divergence
c
        numdsp = 0
        do jrow=1,jmt
          do i=1,imt
	    divf(i,jrow) = c0
	  enddo
        enddo
      endif
# ifdef timing
      call toc ('diagnostic', 'diagnostic surf height')
# endif
#endif
c
#ifdef timing
      call tic ('diagnostic', 'output (diago)')
#endif
c
#ifdef tracer_averages
c
c-----------------------------------------------------------------------
c     compute tracer averages under horizontal regions and write them
c     out
c-----------------------------------------------------------------------
c
      if (tavgts) then
c
c       initialize sums for horizontal regions & tracer averages
c
        do m=1,nt
          sumgf(m) = c0
          avggf(m) = c0
          sumgt(m) = c0
          avggt(m) = c0
          do mask=1,nhreg
            sumbt(mask,m) = c0
            avgbt(mask,m) = c0
            avgbf(mask,m) = c0
          enddo
          do k=1,km
            sumgk(k,m) = c0
            avggk(k,m) = c0
          enddo
        enddo
c
c       compute sums for tracer averages over horizontal regions
c
        do m=1,nt
          do mask=1,nhreg
            sumgf(m) = sumgf(m) + sumbf(mask,m)
            do k=1,km
              sumbt(mask,m) = sumbt(mask,m) + sumbk(mask,k,m)
              sumgk(k,m) = sumgk(k,m) + sumbk(mask,k,m)
            enddo
            sumgt(m) = sumgt(m) + sumbt(mask,m)
          enddo
        enddo
c
        do k=1,km
          if (volgk(k) .gt. c0) then
            rvolgk = c1 / volgk(k)
            do m=1,nt
              avggk(k,m) = sumgk(k,m) * rvolgk
              do mask=1,nhreg
                if (volbk(mask,k) .gt. c0) then
                  avgbk(mask,k,m) = sumbk(mask,k,m) / volbk(mask,k)
                endif
              enddo
            enddo
          endif
        enddo
c
        rvolgt = c1 / volgt
        rareag = c1 / areag
        do m=1,nt
          avggt(m) = sumgt(m) * rvolgt
          avggf(m) = sumgf(m) * rareag
          do mask=1,nhreg
            if (volbt(mask) .gt. c0) then
              avgbt(mask,m) = sumbt(mask,m) / volbt(mask)
            endif
            if (areab(mask) .gt. c0) then
              avgbf(mask,m) = sumbf(mask,m) / areab(mask)
            endif
          enddo
        enddo
c
c       write out regional tracer means
c
        if (iotavg .eq. stdout .or. iotavg .lt. 0) then
	  write (stdout,'(//,30x,a,/)') 'T R A C E R    A V E R A G E S'
          do m=1,nt
            write(stdout,9004) trname(m), itt, stamp
            write(stdout,9001) (mask,mask=1,nhreg)
            do k=1,km
              write(stdout,9002) k, avggk(k,m),
     &                             (avgbk(mask,k,m),mask=1,nhreg)
            enddo
            write(stdout,9003) avggt(m), (avgbt(mask,m),mask=1,nhreg)
            write(stdout,9014) m, avggf(m), (avgbf(msk,m),msk=1,nhreg)
          enddo
        endif
        if (iotavg .ne. stdout .or. iotavg .lt. 0) then
c
	  call getunit (io, 'tracer_avg.dta','u s a ieee')
c
          reltim = prelyr
          write (stdout,*) ' => Regional tracer averages written'
     &,      ' unformatted to file tracer_avg.dta on ts=',itt, ' ',stamp
c
          iotext = 'read (iotavg) reltim, nt, nhreg, km'
          write (io) pstamp, iotext, expnam
          write (io) reltim, nt, nhreg, km
c
          iotext = 'read (iotavg) ((avggk(k,n),k=1,km),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggk, km*nt)
c
          iotext =
     &    'read (iotavg) (((avgbk(l,k,n),l=1,nhreg),k=1,km),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbk, nhreg*km*nt)
c
          iotext = 'read (iotavg) (avggt(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggt, nt)
c
          iotext = 'read (iotavg) (avggf(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggf, nt)
c
          iotext = 'read (iotavg) ((avgbt(l,n),l=1,nhreg),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbt, nhreg*nt)
c
          iotext = 'read (iotavg) ((avgbf(l,n),l=1,nhreg),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbf, nhreg*nt)
c
	  call relunit (io)
        endif
      endif
9001  format('    k','  All Regions ',9(1x,i7,5x))
9002  format(1x,i4,10(1x,e12.6))
9003  format('  AVG',10(1x,e12.6))
9004  format(/' Volume Weighted Averages for ',a12,' on ts =',i10,a33)
9014  format(/' FLX',i1,10(1x,e12.6),/)
#endif
c
#ifdef time_step_monitor
c
c-----------------------------------------------------------------------
c     print integrals for monitoring the time step
c-----------------------------------------------------------------------
c
      if (tsits) then
        do jrow=2,jmtm1
          do k=km,1,-1
            ektot(0,1) = ektot(0,1) + ektot(k,jrow)
	  enddo
          do m=1,nt
            do k=km,1,-1
              tbar(0,m,1)   = tbar(0,m,1) + tbar(k,m,jrow)
              travar(0,m,1) = travar(0,m,1) + travar(k,m,jrow)
              dtabs(0,m,1)  = dtabs(0,m,1) + dtabs(k,m,jrow)
            enddo
	  enddo
	enddo
        ektot(0,1) = ektot(0,1)/ucellv
        do m=1,nt
          tbar(0,m,1)=tbar(0,m,1)/tcellv
          travar(0,m,1)=travar(0,m,1)/tcellv - tbar(0,m,1)**2
          dtabs(0,m,1)=dtabs(0,m,1)/tcellv
	enddo
        if (iotsi .eq. stdout .or. iotsi .lt. 0) then
          write (stdout,9601) itt, stamp, ektot(0,1), dtabs(0,1,1)
     &,                      dtabs(0,2,1), tbar(0,1,1), tbar(0,2,1)
     &,                      travar(0,1,1), travar(0,2,1), mscan
        endif
        if (iotsi .ne. stdout .or. iotsi .lt. 0) then
c
	  call getunit (io, 'ts_intgrls.dta','u s a ieee')
c
          write (stdout,*) ' => Time step integrals written'
     &,     ' unformatted to file ts_intgrls.dta on ts=',itt, ' ',stamp
          reltim = prelyr
          totke  = ektot(0,1)
          tbar1  = tbar(0,1,1)
          tbar2  = tbar(0,2,1)
          tvar1  = travar(0,1,1)
          tvar2  = travar(0,2,1)
          dtabs1 = dtabs(0,1,1)
          dtabs2 = dtabs(0,2,1)
          iotext = 
     &'read(io) reltim, tke, dtabs, dsabs, tbar, sbar, tvar, tvar, nscn'
          write (io) pstamp, iotext, expnam
          write (io) reltim, totke, dtabs1, dtabs2, tbar1, tbar2, tvar1
     &,              tvar2, mscan
c
          call relunit (io)
        endif
      endif
9601  format (1x,'ts=',i7, 1x, a32, ', ke=', 1pe13.6,' |dT|=',1pe13.6
     &,      ' |dS|=',1pe13.6,' Tbar=',1pe13.6,' Sbar=',1pe13.6
     &,      ' Tvar=',1pe13.6,' Svar=',1pe13.6, ' scans=',i4)
#endif
c
#ifdef stability_tests
c
c-----------------------------------------------------------------------
c     show stability and CFL conditions, reynolds and peclet numbers 
c-----------------------------------------------------------------------
c
      if (stabts) then
        write (stdout,'(///20x,a/15x,a,/)')
     & ' S T A B I L I T Y     A N A L Y S I S'
     &, '(The indicated locations are the most unstable ones)'
c
        write (stdout,*) ' longitudinal domain: ',cflons, ' to ',cflone
        write (stdout,*) ' latitudinal  domain: ',cflats, ' to ',cflate
        write (stdout,*) ' depth domain (m)   : ',cfldps*0.01
     &,                  ' to ',cfldpe*0.01
c
        write (stdout,'(/60x,a/)') ' CFL summary'
        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     &   ' Local U velocity (',cflum,') is ',cflup,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflu,',',jcflu,',',kcflu,')'
     &,  ' (lon,lat,dpt) = (',xu(icflu),',',yu(jcflu),','
     &,  0.01*zt(kcflu),')'
c
        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     &   ' Local V velocity (',cflvm,') is ',cflvp,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflv,',',jcflv,',',kcflv,')'
     &,  ' (lon,lat,dpt) = (',xu(icflv),',',yu(jcflv),','
     &,  0.01*zt(kcflv),')'
c
        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     & ' Local adv_vbu    (',cflwum,') is ',cflwup,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflwu,',',jcflwu,',',kcflwu,')'
     &,  ' (lon,lat,dpt) = (',xu(icflwu),',',yu(jcflwu),','
     &,  0.01*zw(kcflwu),')'
c
        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     & ' Local adv_vbt    (',cflwtm,') is ',cflwtp,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflwt,',',jcflwt,',',kcflwt,')'
     &,  ' (lon,lat,dpt) = (',xu(icflwt),',',yu(jcflwt),','
     &,  0.01*zw(kcflwt),')'
c
        fx = 100.0
        if (cflup .gt. fx .or. cflvp .gt. fx .or. cflwup .gt. fx .or.
     &      cflwtp .gt. fx) then
          write (stdout,*)
     &      ' => Warning. CFL exceeded... computational mode exists!'
        endif
c
        write (stdout,'(/60x,a24/)') ' Reynolds number summary'
c
        write (stdout,10300) reynx, ireynx, jreynx, kreynx
     &,                      xu(ireynx), yu(jreynx), 0.01*zt(kreynx)
        write (stdout,10310) reynu, reynmu
c
        write (stdout,10400) reyny, ireyny, jreyny, kreyny
     &,                      xu(ireyny), yu(jreyny), 0.01*zt(kreyny)
        write (stdout,10410) reynv, reynmv
c
        write (stdout,10500) reynz, ireynz, jreynz, kreynz
     &,                      xu(ireynz), yu(jreynz), 0.01*zt(kreynz)
        write (stdout,10510) reynw, reynmw
c
        write (stdout,'(/60x,a22/)') ' Peclet number summary'
c
        write (stdout,10600) peclx, ipeclx, jpeclx, kpeclx
     &,                      xu(ipeclx), yu(jpeclx), 0.01*zt(kpeclx)
        write (stdout,10610) peclu, peclmu
c
        write (stdout,10700) pecly, ipecly, jpecly, kpecly
     &,                      xu(ipecly), yu(jpecly), 0.01*zt(kpecly)
        write (stdout,10710) peclv, peclmv
c
        write (stdout,10800) peclz, ipeclz, jpeclz, kpeclz
     &,                      xu(ipeclz), yu(jpeclz), 0.01*zt(kpeclz)
        write (stdout,10810) peclw, peclmw
c
c       show ficticious tracer extremums
c
        call getunit (iostab, 'iostab', 'fsr')
        rewind iostab
	write (stdout,'(/,10x,a/,11x,a,1pe10.3,a/)')
     &  'Spurious creation of local tracer extremum (if any) follow:'
     &, '(where tracer exceeds local extremum by ',tdig,'*tracer)'
	do num=1,1000
          read (iostab,'(i4, i4, i4, i2, 3g14.7)', end=101, err=101)
     &           i, k, jrow, n, tnew, tsml, tbig
	  write (stdout,'(1x,a,i4,a,i4,a,i4,a,i2,3(a,g14.7))')
     &    't(i,k,jrow,n) = t(',i,',',k,',',jrow,',',n, ') = '
     &,   tnew,' : local min was ',tsml, ' : local max was ', tbig
          if (num .eq. 100) then
	    write (stdout,'(/a/)') 'Bailing out after 100 locations...'
	    go to 101
	  endif
	enddo
101     continue
        call relunit (iostab)
        write (stdout,'(/60x,a/)') ' End Stability Analysis'
      endif
10300 format (1x,'Maximim zonal Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10310 format (1x,' local U =',1pe9.2, ' and  mixing =',e9.2)
10400 format (1x,'Maximim meridional Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10410 format (1x,' local V =',1pe9.2, ' and  mixing =',e9.2)
10500 format (1x,'Maximim vertical Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10510 format (1x,' local Wu =',1pe9.2, ' and  mixing =',e9.2)
10600 format (1x,'Maximim zonal Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10610 format (1x,' local U =',1pe9.2, ' and  mixing =',e9.2)
10700 format (1x,'Maximim meridional Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10710 format (1x,' local V =',1pe9.2, ' and  mixing =',e9.2)
10800 format (1x,'Maximim vertical Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10810 format (1x,' local Wt =',1pe9.2, ' and  mixing =',e9.2)
#endif
c
#ifdef energy_analysis
c
c-----------------------------------------------------------------------
c     add external mode component of total work done 
c-----------------------------------------------------------------------
c
      if (glents) then
# ifdef timing
        call toc ('diagnostic', 'output (diago)')
# endif
        call ge3 (c2dtuv)
# ifdef timing
        call ticr ('diagnostic', 'output (diago)')
# endif
      endif
c
c-----------------------------------------------------------------------
c     integrate previously computed energy components vertically
c-----------------------------------------------------------------------
c
      if (glents) then
        jte  = 1
	jue  = 1
        jwte = 1
	jwue = 1
	do k=1,km
	  wtlev(k,0) = c0
	  wulev(k,0) = c0
	enddo
        do jrow=2,jmt
          do k=km,1,-1
            buoy(0,1) = buoy(0,1) + buoy(k,jrow)
	    wtlev(k,0) = wtlev(k,0) + wtlev(k,jrow)
	    wulev(k,0) = wulev(k,0) + wulev(k,jrow)
          enddo
          do ll=1,8
            do k=km,1,-1
              engint(0,ll,1) = engint(0,ll,1) + engint(k,ll,jrow)
            enddo
	    engext(ll,1) = engext(ll,1) + engext(ll,jrow)
          enddo
	  if (abs(tcerr(jrow)) .gt. abs(tcerr(jte))) jte = jrow
	  if (abs(ucerr(jrow)) .gt. abs(ucerr(jte))) jue = jrow
	  if (abs(wtbot(jrow)) .gt. abs(wtbot(jwte))) jwte = jrow
	  if (abs(wubot(jrow)) .gt. abs(wubot(jwue))) jwue = jrow
        enddo
        buoy(0,1) = buoy(0,1)/ucellv
        do ll=1,8
          engint(0,ll,1) = engint(0,ll,1)/ucellv
          engext(ll,1)   = engext(ll,1)/ucellv
        enddo
c
        plicin = engint(0,1,1) - engint(0,2,1) - engint(0,3,1)
     &           - engint(0,4,1) - engint(0,5,1) - engint(0,6,1)
        plicex = engext(1,1) - engext(2,1) - engext(3,1)
     &            - engext(4,1)  - engext(5,1) - engext(6,1)
        buoerr = buoy(0,1) - engint(0,6,1) - engext(6,1)
        enleak = engint(0,2,1) + engint(0,3,1) + engext(2,1)
     &         + engext(3,1)
c
        if (ioglen .eq. stdout .or. ioglen .lt. 0) then
	  write (stdout,'(///40x,a,/)') 'E N E R G Y    A N A L Y S I S'
          write (stdout,9100)
     &              'Globally averaged work done'
     &,                     ' for ts =', itt, stamp, ucellv, ucella(1)
          write (stdout,9101) ' time rate of change ',engint(0,1,1)
     &,                      engint(0,1,1), engext(1,1), engext(1,1)
          write (stdout,9101) ' horizontal advection', engint(0,2,1)
     &,                      engint(0,2,1), engext(2,1), engext(2,1)
          write (stdout,9101) ' vertical advection  ',engint(0,3,1)
     &,                      engint(0,3,1), engext(3,1), engext(3,1)
          write (stdout,9101) ' horizontal friction ',engint(0,4,1)
     &,                      engint(0,4,1), engext(4,1), engext(4,1)
          write (stdout,9101) ' vertical friction   ',engint(0,5,1)
     &,                      engint(0,5,1), engext(5,1), engext(5,1)
          write (stdout,9101) ' pressure forces     ',engint(0,6,1)
     &,                      engint(0,6,1), engext(6,1), engext(6,1)
          write (stdout,9101) ' ficticious sources  ',plicin
     &,                       plicin, plicex, plicex
          write (stdout,9101) ' work by wind        ',engint(0,7,1)
     &,                      engint(0,7,1), engext(7,1), engext(7,1)
          write (stdout,9101) ' bottom drag         ',engint(0,8,1)
     &,                      engint(0,8,1), engext(8,1), engext(8,1)
          write (stdout,9110) buoy(0,1), buoy(0,1), buoerr, buoerr
     &,                      enleak, enleak
          write (stdout,9111) tcerr(jte), itcerr(jte), jtcerr(jte)
     &,                       ktcerr(jte)
     &,                       ucerr(jue), iucerr(jue), jucerr(jue)
     &,                       kucerr(jue)
          write (stdout,9112) wtbot(jwte), iwtbot(jwte), jwtbot(jwte)
     &,                       kwtbot(jwte)
     &,                       wubot(jwue), iwubot(jwue), jwubot(jwue)
     &,                       kwubot(jwue)
c
          write (stdout,'(/a,a,//a,a,a,a,a,/1x,a,a)')
     &    'Average vertical velocity through bottom of "T" and "u" '  
     &,   'cells at each level','level', '   adv_vbt err '
     &, '  "T" cell area ', '    adv_vbu err ','  "U" cell area'
     &, '(Note: adv_vbu err only goes to zero when non-zero values on'
     &, ' boundary cells are taken into account)'
          do k=1,km
	    if (tcella(k) .ne. c0) then
	      wtlev(k,0) = wtlev(k,0)/tcella(k)
	    else
	      wtlev(k,0) = c0
	    endif
	    if (ucella(k) .ne. c0) then
	      wulev(k,0) = wulev(k,0)/ucella(k)
	    else
	      wulev(k,0) = c0
	    endif
	    write (stdout,'(i4,4(2x,e14.7))') 
     &      k, wtlev(k,0), tcella(k), wulev(k,0), ucella(k)
          enddo
        endif
c
        if (ioglen .ne. stdout .or. ioglen .lt. 0) then
c
	  call getunit (io, 'energy_int.dta','u s a ieee')
c
          write (stdout, *)
     &      ' ==> Global energy integrals written unformatted'
     &,     ' to file energy_int.dta on ts =', itt, stamp
          reltim = prelyr
          iotext = 'read(ioglen) reltim,(engint(i),engext(i),i=1,8)'
          iotext(46:) = ',plicin,plicex,buoy,buoerr,enleak'
          write (io) pstamp, iotext, expnam
          write (io)  reltim, (engint(0,i,1),engext(i,1),i=1,8)
     &,                 plicin, plicex, buoy(0,1), buoerr, enleak
          call relunit (io)
        endif
      endif
9100  format(///,1x,
     &/1x,a,a,i10,a/1x,'ocean volume =',e16.9,' cm**3'
     &, ', ocean surface area =',e16.9,' cm**2'//' work by:',14x
     &,     'internal mode                         external mode'/)
9101  format(a21,2(1pe15.6, ' (',z16,' hex)'))
9110  format(/' work by buoyancy forces   =',1pe14.6, ' (',z16,' hex)'/
     &,       ' energy conversion error   =',1pe14.6, ' (',z16,' hex)'/
     &,       ' nonlinear error           =',1pe14.6, ' (',z16,' hex)'/)
9111  format(/' max "t" cell continuity error =',1pe14.6, ' at location'
     &,       ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')'
     &,      /' max "u" cell continuity error =',1pe14.6, ' at location'
     &,       ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')')
9112  format(/' max bottom "adv_vbt" (error)       =',1pe14.6
     &, ' at location', ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')'
     &,      /' max bottom "adv_vbu" (slope vel)   =',1pe14.6
     &, ' at location',' (i,jrow,k) = ','(',i4,',',i4,',',i4,')') 
#endif
c
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     add the external mode part of d/dt into the momentum balance,
c     the external mode part of the implicit coriolis term, and the
c     surface pressure gradients into the specified volumes
c-----------------------------------------------------------------------
c
      if (trmbts) then
c
# ifdef timing
        call toc ('diagnostic', 'output (diago)')
# endif
        call utb3
# ifdef timing
        call ticr ('diagnostic', 'output (diago)')
# endif
c
c-----------------------------------------------------------------------
c     integrate previously computed term balance components vertically
c-----------------------------------------------------------------------
c
        do n=0,numreg
          if (n .gt. 0) then
            nv = (n-1)/nhreg + 1
            ks = llvreg(nv,1)
            ke = llvreg(nv,2)
            do ll=1,17
              do k=ke,ks,-1
                termbm(0,ll,1,n) = termbm(0,ll,1,n) + termbm(k,ll,1,n)
                termbm(0,ll,2,n) = termbm(0,ll,2,n) + termbm(k,ll,2,n)
	      enddo
	    enddo
          else
            ks = 1
            ke = km
          endif
          do m=1,nt
            do ll=1,15
              do k=ke,ks,-1
                termbt(0,ll,m,n) = termbt(0,ll,m,n) + termbt(k,ll,m,n)
              enddo
	    enddo
c
c           construct change due to convection and filtering
c
            dtconv = termbt(0,10,m,n) - termbt(0,9,m,n)
	    dtfilt = termbt(0,1,m,n) - termbt(0,10,m,n)
            termbt(0,9,m,n) = dtconv
	    termbt(0,10,m,n) = dtfilt
	  enddo
	enddo
c
c       normalize integrals by appropriate volume (or area)
c
        do n=0,numreg
          do m=1,nt
            if (n .le. nhreg)  then
              stflx(m,n) = stflx(m,n)*rareat(n)
              asst(m,n)  = asst(m,n)*rareat(n)
            endif
            do ll=1,15
              termbt(0,ll,m,n) = termbt(0,ll,m,n)*rvolt(n)
	    enddo
	  enddo
	enddo
c
        do n=0,numreg
          if (n .le. nhreg) then
	    smflx(1,0) = smflx(1,0) + smflx(1,n)
	    smflx(2,0) = smflx(1,0) + smflx(2,n)
            smflx(1,n) = smflx(1,n)*rareau(n)
            smflx(2,n) = smflx(2,n)*rareau(n)
          endif
          if (n .gt. 0) then
            avgw(n) = avgw(n)*rvolu(n)
            do ll=1,17
              termbm(0,ll,1,n) = termbm(0,ll,1,n)*rvolu(n)
              termbm(0,ll,2,n) = termbm(0,ll,2,n)*rvolu(n)
	    enddo
          endif
	enddo
	smflx(1,0) = smflx(1,0)*rareau(0)
	smflx(2,0) = smflx(2,0)*rareau(0)
c
        if (iotrmb .eq. stdout .or. iotrmb .lt. 0) then
c
          write (stdout,'(///,40x,a,/)') 'T E R M    B A L A N C E S'
c
          n = 0
          taux = smflx(1,n)
          tauy = smflx(2,n)
          write (stdout,10106)
     &                   'All regions added together for ts ='
     &,                     itt, stamp, volt(n), areat(n)
     &,                     volu(n), areau(n)
          write (stdout,10104)
          write (stdout,10098) n, ' smf(1)   = ', taux,' dynes/cm**2   '
          write (stdout,10098) n, ' smf(2)   = ', tauy,' dynes/cm**2   '
          do m=1,nt
            write (stdout,10098) n, ustf(m,1), stflx(m,n), ustf(m,2)
	  enddo
          write (stdout,10098) n, ' tot heat = '
     &,                           (termbt(0,15,1,n)*volt(n))
     &,                                        ' deg C * cm**3 '
          write (stdout,10098) n, ' sst      = ',asst(1,n)
     &,                                       ' deg C         '
c
          do n=0,numreg
            if (n .eq. 1) then
              write (stdout,10050)
     &  'Regional averaged Momentum & Tracer Term Balances for ts =   '
     &,       itt, stamp
            endif
            iv = 0
            if (n .gt. 0) then
              iv = (n-1)/nhreg + 1
              ih = n - (iv-1)*nhreg
              write (stdout,10100)
     &        'Momentum terms averaged over region #'
     &,          n, ': ', hregnm(ih), vregnm(iv), volu(n), areau(n)
              write (stdout,10104)
              write (stdout,10101) n,' dU/dt   = ', termbm(0,1,1,n)
     &,                              ' dV/dt   = ', termbm(0,1,2,n)
              write (stdout,10101) n,' -Px     = ', termbm(0,2,1,n)
     &,                              ' -Py     = ', termbm(0,2,2,n)
              write (stdout,10101) n,' -surf Px= ', termbm(0,12,1,n)
     &,                              ' -surf Py= ', termbm(0,12,2,n)
              write (stdout,10101) n,' -(UU)x  = ', termbm(0,3,1,n)
     &,                              ' -(UV)x  = ', termbm(0,3,2,n)
              write (stdout,10101) n,' -(VU)y  = ', termbm(0,4,1,n)
     &,                              ' -(VV)y  = ', termbm(0,4,2,n)
              write (stdout,10101) n,' -(WU)z  = ', termbm(0,5,1,n)
     &,                              ' -(WV)z  = ', termbm(0,5,2,n)
              write (stdout,10101) n,'ADV_Umet = ', termbm(0,13,1,n)
     &,                              '-ADV_Vmet= ', termbm(0,13,2,n)
              write (stdout,10101) n,'  DIFF_Ux= ', termbm(0,6,1,n)
     &,                              '  DIFF_Vx= ', termbm(0,6,2,n)
              write (stdout,10101) n,'  DIFF_Uy= ', termbm(0,7,1,n)
     &,                              '  DIFF_Vy= ', termbm(0,7,2,n)
              write (stdout,10101) n,'  DIFF_Uz= ', termbm(0,8,1,n)
     &,                              '  DIFF_Vz= ', termbm(0,8,2,n)
              write (stdout,10101) n,'DIFF_Umet= ', termbm(0,9,1,n)
     &,                              'DIFF_Vmet= ', termbm(0,9,2,n)
              write (stdout,10101) n,'  fV     = ', termbm(0,10,1,n)
     &,                              ' -fU     = ', termbm(0,10,2,n)
              write (stdout,10101) n,'  source = ', termbm(0,11,1,n)
     &,                              '  source = ', termbm(0,11,2,n)
              erru = c0
              errv = c0
              do lll=2,13
                erru = erru + termbm(0,lll,1,n)
                errv = errv + termbm(0,lll,2,n)
              enddo
              write (stdout,10101) n,'  error  = ', termbm(0,1,1,n)-erru
     &,                              '  error  = ', termbm(0,1,2,n)-errv
c
              write (stdout,*) ' '
              write (stdout,10101) n,' -U(U)x  = ', termbm(0,14,1,n)
     &,                              ' -U(V)x  = ', termbm(0,14,2,n)
              write (stdout,10101) n,' -V(U)y  = ', termbm(0,15,1,n)
     &,                              ' -V(V)y  = ', termbm(0,15,2,n)
              write (stdout,10101) n,' -W(U)z  = ', termbm(0,16,1,n)
     &,                              ' -W(V)z  = ', termbm(0,16,2,n)
c
c             mass conservation within volume: Ux + Vy + Wz
c
              contu = (termbm(0,3,1,n) + termbm(0,4,1,n) + 
     &                termbm(0,5,1,n)) - (termbm(0,14,1,n) +
     &                termbm(0,15,1,n) + termbm(0,16,1,n))
              write (stdout,10101) n,' mass err= ', contu
c
              write (stdout,10101) n,'  ubar   = ', termbm(0,17,1,n)
     &,                              '  vbar   = ', termbm(0,17,2,n)
     &,                              '  wbar   = ', avgw(n)
            endif
            if (iv .eq. 1) then
              write (stdout,10101) n,'  surf Uz= ', smflx(1,n)
     &,                              '  surf Vz= ', smflx(2,n)
            endif
c
            if (n .eq. 0) then
              write (stdout,10051)
     &   'Global averaged (entire ocean) Tracer Term Balances for ts = '
     &,       itt, stamp
            else
              write (stdout,10100)
     &        'Tracer terms averaged over region   #'
     &,        n, ': ', hregnm(ih), vregnm(iv), volt(n), areat(n)
            endif
c
            do m=1,nt
              dchg            = termbt(0,2,m,n) + termbt(0,3,m,n) +
     &                          termbt(0,4,m,n) + termbt(0,5,m,n) +
     &                          termbt(0,6,m,n) + termbt(0,7,m,n) +
     &                          termbt(0,8,m,n) + termbt(0,9,m,n) +
     &                          termbt(0,10,m,n)
              terr(m)        = termbt(0,1,m,n) - dchg
            enddo
            maxm = (nt-1)/7 + 1
            do mloop=1,maxm
              ms = (mloop-1)*7 + 1
              me = min(ms + 6,nt)
              write (stdout,10103) (trname(m),m=ms,me)
              write (stdout,10102) n,' dT/dt   = ', (termbt(0,1,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(UT)x  = ', (termbt(0,2,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(VT)y  = ', (termbt(0,3,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(WT)z  = ', (termbt(0,4,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Tx= ', (termbt(0,5,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Ty= ', (termbt(0,6,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Tz= ', (termbt(0,7,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  source = ', (termbt(0,8,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  convct = ', (termbt(0,9,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  filter = ', (termbt(0,10,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  error  = ', (terr(m)
     &,                                           m=ms,me)
c
              write (stdout,*) ' '
              write (stdout,10102) n,' -U(T)x  = ', (termbt(0,11,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -V(T)y  = ', (termbt(0,12,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -W(T)z  = ', (termbt(0,13,m,n)
     &,                                           m=ms,me)
c
c             mass conservation within volume: Ux + Vy + Wz
c
              cont = (termbt(0,2,1,n) + termbt(0,3,1,n) + 
     &               termbt(0,4,1,n)) - (termbt(0,11,1,n) +
     &               termbt(0,12,1,n) + termbt(0,13,1,n))
              write (stdout,10102) n,' mass err= ', cont
c
              write (stdout,10102) n,'  chg var= ', (termbt(0,14,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  tbar   = ', (termbt(0,15,m,n)
     &,                                           m=ms,me)
              if (iv .eq. 1) then
                write (stdout,10102) n,'  surflx = ', (stflx(m,n)
     &,                                           m=ms,me)
                taux = smflx(1,n)
                tauy = smflx(2,n)
                write (stdout,10105) ' Regionally averaged quantities:'
                write (stdout,'(1x)')
                write (stdout,10098) n
     &,            ' smf(1)   = ', taux,' dynes/cm**2   '
                write (stdout,10098) n
     &,            ' smf(2)   = ', tauy,' dynes/cm**2   '
                do m=1,nt
                  write (stdout,10098) n
     &,            ustf(m,1), stflx(m,n), ustf(m,2)
                enddo
                if (ms .eq. 1) then
                  write (stdout,10098) n, ' tot heat = '
     &,             (termbt(0,15,1,n)*volt(n)),' deg C * cm**3 '
                  write (stdout,10098) n
     &,            ' sst      = ', asst(1,n),' deg C         '
                endif
              endif
            enddo
          enddo
        endif
c
        if (iotrmb .ne. stdout .or. iotrmb .lt. 0) then
c
	  call getunit (io, 'term_bal.dta','u s a ieee')
c
          write (stdout, *)
     &     ' ==> Term balances written unformatted to file term_bal.dta'
     &,    ' for ts =', itt, ' ',stamp
c
          reltim = prelyr
c
          iotext =' read (iotrmb) reltim, nt, numreg, nhreg, km'
          iotext(45:) =', ((ustf*15(n,i),n=1,nt),i=1,2)'
          write (io) pstamp, iotext, expnam
          write (io)  reltim, nt, numreg, nhreg, km, ustf
c
          iotext ='read (iotrmb) (hregnm*40(n),n=1,nhreg)'
          iotext(39:)=', (vregnm*20(n),n=1,nvreg)'
          write (io) pstamp, iotext, expnam
          write (io)  hregnm, vregnm
c
          iotext ='read (iotrmb) (trname*12(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          write (io)  trname
c
          iotext ='read (iotrmb) (volu(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, volu, numreg+1)
c
          iotext ='read (iotrmb) (volt(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, volt, numreg+1)
c
          iotext ='read (iotrmb) (areau(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, areau, numreg+1)
c
          iotext ='read (iotrmb) (areat(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, areat, numreg+1)
c
          iotext ='read (iotrmb) ((smflx(i,l),i=1,2),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, smflx, 2*(nhreg+1))
c
          iotext ='read (iotrmb) ((stflx(n,l),n=1,nt),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, stflx, nt*(nhreg+1))
c
          iotext ='read (iotrmb) ((asst(n,l),n=1,nt),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, asst, nt*(nhreg+1))
c
          iotext ='read (iotrmb) (avgw(l),l=1,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgw, numreg)
c
          iotext = 'read (iotrmb) (((termbm'
          iotext(24:) = '(i,n,l),i=1,17),n=1,2),l=1,numreg)'
          write (io) pstamp, iotext, expnam
          write (io) (((termbm(0,i,n,l),i=1,17),n=1,2),l=1,numreg)
c
          iotext = 'read (iotrmb) (((termbt'
          iotext(24:) = '(i,n,l),i=1,15),n=1,nt),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          write (io) (((termbt(0,i,n,l),i=1,15),n=1,nt),l=0,numreg)
c
          call relunit (io)
        endif
      endif
10050 format(///,1x,/1x,a,i10,a/)
10051 format(///,1x,a,i10,a/)
10098 format(1x,'(',i4,')',a12,(1pe16.8,a15))
10100 format(/1x,a,i5,a,a,a/1x,'ocean volume =',e16.9,' cm**3'
     &, ', horizontal ocean area =',e16.9,' cm**2'/)
10101 format(1x,'(',i4,')',a11,1pe15.7, 2x, a11,1pe15.7, 2x, a11
     &,      1pe15.7)
10102 format(1x,'(',i4,')',a11,7(1pe15.7))
10103 format(' Region',11x,8a15)
10104 format(' Region')
10105 format(/8x,a32)
10106 format(///,1x,
     &/1x,a,i10,a/1x,'ocean "t" volume =',e16.9,' cm**3 '
     &, ', ocean "t" surface area =',e16.9,' cm**2 '
     &/ 1x,'ocean "u" volume =',e16.9,' cm**3, '
     &, ' ocean "u" surface area =',e16.9,' cm**2'/)
#endif
c
#ifdef gyre_components
c
c-----------------------------------------------------------------------
c     write out the gyre_components
c-----------------------------------------------------------------------
c
c     convert heat transport to petawatts,
c     salt transport to 10**10 cm**3/sec
c
      if (gyrets) then
        pwatts = 4.186e-15
        csalt  = 1.e-10
        do jrow=1,jmt
          do ll=1,8
            ttn(ll,jrow,1)=ttn(ll,jrow,1)*pwatts
            ttn(ll,jrow,2)=ttn(ll,jrow,2)*csalt
          enddo
	enddo
        if (iogyre .eq. stdout .or. iogyre .lt. 0) then
c
          write (stdout,'(///,40x,a,/)') 'G Y R E   C O M P O N E N T S'
c
          write (stdout,8195)
          do jrow=2,jmtm2
            l = jmt - jrow
            write (stdout,8196) l, (ttn(i,l,1),i=1,8)
     &,                            (ttn(i,l,2),i=1,8)
          enddo
c
          do m=0,nhreg
            do jrow=1,jmt
#if defined isopycmix && defined gent_mcwilliams
              do ll=6,9
#else
              do ll=6,8
#endif
                ttn2(ll,jrow,1,m) = ttn2(ll,jrow,1,m)*4.186e-15
                ttn2(ll,jrow,2,m) = ttn2(ll,jrow,2,m)*1.e-10
              enddo
            enddo
	  enddo
c
          write (stdout,82001)
          do m=0,nhreg
            if (m .ne. 0) then
              write (stdout,8201) hregnm(m)
            else
              write (stdout,8202)
            endif
            write (stdout,8203)
            do jrow=1,jmt
#if defined isopycmix && defined gent_mcwilliams
              write (stdout,8204) jrow,(ttn2(i,jrow,1,m),i=6,9)
#else
              write (stdout,8204) jrow,(ttn2(i,jrow,1,m),i=6,8)
#endif
            enddo
	  enddo
c
          write (stdout,8205)
          do m=0,nhreg
            if (m .ne. 0) then
              write (stdout,8201) hregnm(m)
            else
              write (stdout,8202)
            endif
            write (stdout,8203)
            do jrow=1,jmt
#if defined isopycmix && defined gent_mcwilliams
              write (stdout,8204) jrow,(ttn2(i,jrow,2,m),i=6,9)
#else
              write (stdout,8204) jrow,(ttn2(i,jrow,2,m),i=6,8)
#endif
            enddo
	  enddo
        endif
        if (iogyre .ne. stdout .or. iogyre .lt. 0) then
c
	  call getunit (io, 'gyre_comp.dta','u s a ieee')
c
          write (stdout,*) ' => Gyre components written'
     &,      ' unformatted to file gyre_comp.dta on ts=',itt, ' ',stamp
          reltim = prelyr
c
          iotext = 'read (iogyre) jmt, ntmin2, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, ntmin2, reltim
c
          iotext =
     &      'read (iogyre) (((ttn(l,j,n),l=1,8),j=1,jmt),n=1,ntmin2)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, ttn, 8*jmt*ntmin2)
c
#if defined isopycmix && defined gent_mcwilliams
          iotext =
     &    'read (iogyre) (((ttn2(6..9,1..jmt,1..2,0..nhreg)'
          write (io) pstamp, iotext, expnam
	  len = 4*jmt*nt*(1+nhreg)
#else
          iotext =
     &    'read (iogyre) (((ttn2(6..8,1..jmt,1..2,0..nhreg)'
          write (io) pstamp, iotext, expnam
	  len = 3*jmt*nt*(1+nhreg)
#endif
          call wrufio (io, ttn2, len)
c
          call relunit (io)
        endif
      endif
8195  format(//,30x,'Gyre Components'/
     &,6x,'northward transport of heat (x10**15 watts)'
     &       ,22x,'northward transport of salt (x10**10 cm**3/sec)',/,
     &       3x,2(3x,'x mean  x eddy  z mean  z eddy   ekman tot adv  ',
     &       'diffus   total'))
8196  format(i4,8f8.3,1x,8f8.3)
82001 format (/,6x,'northward transport of heat (x10**15 watts)',/)
8201  format (/,22x,a40,/)
8202  format (/,22x,' Global ',/)
#ifdef isopycmix
8203  format (8x,'total advection',2x,'total diffusion',2x,'total',
     >        12x,'isop. advection')
8204  format (2x,i3,3x,4(e12.6,5x))
#else
8203  format (8x,'total advection',2x,'total diffusion',2x,'total')
8204  format (2x,i3,3x,3(e12.6,5x))
#endif
8205  format (/,6x,'northward transport of salt (x10**10 cm**3/sec)',/)
#endif
c
#ifdef meridional_overturning
c
c-----------------------------------------------------------------------
c     write out the meridional overturning (mass)
c-----------------------------------------------------------------------
c
      if (vmsfts) then
        if (iovmsf .eq. stdout .or. iovmsf .lt. 0) then
c
          write (stdout,'(///,40x,a,/)')
     &    'M E R I D I O N A L    O V E R T U R N I N G'
c
          scl = 1.e12
          write (stdout,8194)
          js = indp (slatxy, yt, jmt)
          je = indp (elatxy, yt, jmt)
          call matrix (vsf, jmt, js, je, 1, km, scl)
        endif
        if (iovmsf .ne. stdout .or. iovmsf .lt. 0) then
c
	  call getunit (io, 'overturn.dta','u s a ieee')
c
          write (stdout,*) ' => Meridional overturning of mass'
     &,    ' written unformatted to file overturn.dta on ts=',itt
     &,    ' ',stamp
          reltim = prelyr
c
          iotext = 'read (iovmsf) jmt, km, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, km, reltim
c
          iotext = 'read (iovmsf)  (zw(k),k=1,km)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zw, km)
c
          iotext = 'read (iovmsf)  (yu(j),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, yu, jmt)
c
          iotext = 'read (iovmsf)  ((vsf(j,k),j=1,jmt),k=1,km)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, vsf, jmt*km)
c
          call relunit (io)
        endif
      endif
8194  format(/,' meridional overturning stream function (sverdrups)')
#endif
c
#ifdef show_external_mode
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c-----------------------------------------------------------------------
c     write out the external mode velocities and surface pressure
c-----------------------------------------------------------------------
c
      if (extts .and. eots) then
        if (ioext .eq. stdout .or. ioext .lt. 0) then
c
          write (stdout,'(///,40x,a,/)') 'E X T E R N A L    M O D E'
c
          is = indp (slonxy, xt, imt)
          ie = indp (elonxy, xt, imt)
          js = indp (slatxy, yt, jmt)
          je = indp (elatxy, yt, jmt)
          scl=1.0
c
          write (stdout,8000) ' ubar '
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ubar(1,1,1), imt, is, ie, -js, -je, scl)
          write (stdout,8000) ' vbar '
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ubar(1,1,2), imt, is, ie, -js, -je, scl)
          scl = 0.0
          write (stdout,8000) ' div ps'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (divf(1,1), imt, is, ie, -js, -je, scl)
          scl = grav
c          write (stdout,8000) ' surface press(gm/cm/sec**2)'
          write (stdout,8000)
     &   'sea surface height (surface_press/grav) in cm'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ps(1,1,1), imt, is, ie, -js, -je, scl)
        endif
        if (ioext .ne. stdout .or. ioext .lt. 0) then
c
#ifdef rigid_lid_surface_pressure 
	  call getunit (io, 'surf_press.dta','u s a ieee')
#endif
#ifdef implicit_free_surface 
	  call getunit (io, 'ifree_surf.dta','u s a ieee')
#endif
c
          reltim = relyr
c
          iotext = 'read (ioext) imt, jmt, reltim'
          write (io) stamp, iotext, expnam
          write (io) imt, jmt, reltim
c
          iotext = 'read (ioext)  (xt(i),i=1,imt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, xt, imt)
c
          iotext = 'read (ioext)  (yt(j),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, yt, jmt)
c
          iotext = 'read (ioext)  ((ps(i,j),i=1,imt),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, ps(1,1,1), imt*jmt)
          write (stdout,*) ' => Surface pressure written unformatted'
     &,              ' to unit ',io, ' on ts=',itt, stamp
#ifdef rigid_lid_surface_pressure 
     &,              ' to file surf_press.dta on ts=',itt, stamp
#endif
#ifdef implicit_free_surface 
     &,              ' to file ifree_surf.dta on ts=',itt, stamp
#endif
c
          call relunit (io)
        endif
      endif
# endif
# ifdef stream_function
c
c-----------------------------------------------------------------------
c     write out the stream function
c-----------------------------------------------------------------------
c
      if (extts .and. eots) then
        if (ioext .eq. stdout .or. ioext .lt. 0) then
c
          write (stdout,'(///,40x,a,/)') 'E X T E R N A L    M O D E'
c
          is = indp (slonxy, xt, imt)
          ie = indp (elonxy, xt, imt)
          js = indp (slatxy, yt, jmt)
          je = indp (elatxy, yt, jmt)
#  ifdef debug_extmode
c
c         construct ubar and vbar for comparison with
c         rigid_lid_surface_pressure
c
          scl=1.0
          do jrow=1,jmt-1
	    j = jrow
            do i=2,imt-1
              diag1       = psi(i+1,jrow+1,1) - psi(i  ,jrow,1)
              diag0       = psi(i  ,jrow+1,1) - psi(i+1,jrow,1)
	      ext(i,j,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
	      ext(i,j,2)  = (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
            enddo
	  enddo
	  call border (ext(1,1,1), 'u even')
	  call border (ext(1,1,2), 'u odd')
c
          write (stdout,8000) ' ubar '
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ext(1,1,1), imt, is, ie, -js, -je, scl)
          write (stdout,8000) ' vbar '
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ext(1,1,2), imt, is, ie, -js, -je, scl)
#  endif
c
          scl=1.e12
          write (stdout,8000) ' stream function (sverdrups)'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (psi(1,1,1), imt, is, ie, -js, -je, scl)
        endif
        if (ioext .ne. stdout .or. ioext .lt. 0) then
c
	  call getunit (io, 'psi.dta','u s a ieee')
c
          reltim = relyr
c
          iotext = 'read (ioext) imt, jmt, reltim'
          write (io) stamp, iotext, expnam
          write (io) imt, jmt, reltim
c
          iotext = 'read (ioext)  (xt(i),i=1,imt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, xt, imt)
c
          iotext = 'read (ioext)  (yt(j),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, yt, jmt)
c
          iotext = 'read (ioext)  ((psi(i,j),i=1,imt),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, psi(1,1,1), imt*jmt)
          write (stdout,*) ' => Stream function written unformatted'
     &,              ' to file psi.dta on ts=',itt, stamp
c
          call relunit (io)
        endif
      endif
# endif
8000  format(1x,a,1x,'ts=',i7
     &,', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     &,', scaling=',1pg10.3)
#endif
c
#ifdef meridional_tracer_budget
c
c-----------------------------------------------------------------------
c     construct and write out the averaged components for the
c     meridional tracer balance
c-----------------------------------------------------------------------
c
      if (tmbts) then
c
c       quantities have been weighted by volume at each cell and
c       summed over all ocean longitudes and depths for each latitude.
c       average quantities over "numtmb" days and do units conversion
c       to petawatts for temperature and giga cm**3/sec for salt
c
        tmbavg = numtmb * dtts	
        rnum = c1 / float(numtmb)
c
        do mask=1,ntmbb
          do jrow=1,jmt
            smdvol(jrow,mask) = rnum*smdvol(jrow,mask)
          enddo
	enddo
c
        do mask=1,ntmbb
          do n=1,nt
	    if (n .eq. 1) then
	      convrt = rnum * 4.186e-15
	    else if (n .eq. 2) then
	      convrt = rnum * 1.e-9
	    else
	      convrt = c1
	    endif
            do jrow=1,jmt
              tstor(jrow,n,mask) = convrt*tstor(jrow,n,mask)
              tdiv(jrow,n,mask)  = convrt*tdiv(jrow,n,mask)
              tflux(jrow,n,mask) = convrt*tflux(jrow,n,mask)
              tdif(jrow,n,mask)  = convrt*tdif(jrow,n,mask)
              tsorc(jrow,n,mask) = convrt*tsorc(jrow,n,mask)
            enddo
          enddo
	enddo
c
c       accumulate terms for all basins and store in basin # 0
c
        do jrow=1,jmt
          smdvol(jrow,0) = c0
          do mask=1,ntmbb
            smdvol(jrow,0) = smdvol(jrow,0) + smdvol(jrow,mask)
          enddo
	enddo
c
        do n=1,nt
          do jrow=1,jmt
            tstor(jrow,n,0) = c0
            tdiv(jrow,n,0)  = c0
            tflux(jrow,n,0) = c0
            tdif(jrow,n,0)  = c0
            tsorc(jrow,n,0) = c0
          enddo
	enddo
c
        do mask=1,ntmbb
          do n=1,nt
            do jrow=1,jmt
              tstor(jrow,n,0) = tstor(jrow,n,0) + tstor(jrow,n,mask)
              tdiv(jrow,n,0)  = tdiv(jrow,n,0)  + tdiv(jrow,n,mask)
              tflux(jrow,n,0) = tflux(jrow,n,0) + tflux(jrow,n,mask)
              tdif(jrow,n,0)  = tdif(jrow,n,0)  + tdif(jrow,n,mask)
              tsorc(jrow,n,0) = tsorc(jrow,n,0) + tsorc(jrow,n,mask)
            enddo
          enddo
	enddo
c
c       write out the results
c
        if (iotmb .eq. stdout .or. iotmb .lt. 0) then
c
          write (stdout,'(///,20x,a,/)')
     &	   'M E R I D I O N A L    T R A C E R    B U D G E T'
c
          do n=1,nt
c
c           construct sum in latitude for all basins
c
            sumsto = c0
            sumdiv = c0
            sumflx = c0
            sumdif = c0
            sumsor = c0
            sumvol = c0
            do jrow=1,jmt
              sumsto = sumsto + tstor(jrow,n,0)
              sumdiv = sumdiv + tdiv(jrow,n,0)
              sumflx = sumflx + tflux(jrow,n,0)
              sumdif = sumdif + tdif(jrow,n,0)
              sumsor = sumsor + tsorc(jrow,n,0)
	      sumvol = sumvol + smdvol(jrow,0)
            enddo
c
	    write (stdout,*) ' '
	    if (n .eq. 1) then
	      write (stdout,*) '(Tracer # 1 results are in petawatts)'
	    else if (n .eq. 2) then
	      write (stdout,*) 
     &       '(Tracer # 2 results are in giga cm**3/sec)'
	    else
	      write (stdout,*) ' '
	    endif
            write (stdout,8050) n, tmbavg*secday, stamp
c
            do jrow=jmt,1,-1
              terror = tstor(jrow,n,0) - (tdiv(jrow,n,0) + 
     &               tflux(jrow,n,0) + tdif(jrow,n,0) + tsorc(jrow,n,0))
              write (stdout,8100) jrow, yt(jrow), tstor(jrow,n,0)
     &,       tdiv(jrow,n,0), tflux(jrow,n,0), tdif(jrow,n,0)
     &,       tsorc(jrow,n,0), terror, smdvol(jrow,0)
            enddo
            serror = sumsto - (sumdiv + sumflx + sumdif + sumsor)
            write (stdout,8200)  sumsto, sumdiv, sumflx, sumdif
     &,                          sumsor, serror, sumvol
          enddo
        endif
        if (iotmb .ne. stdout .or. iotmb .lt. 0) then
c
	  call getunit (io, 'tracer_bud.dta','u s a ieee')
c
          write (stdout,*) ' => Meridional Tracer Budget written '
     &,   ' unformatted to file tracer_bud.dta on ts=',itt, ' ',stamp
          reltim = relyr
c
          iotext = 'read (iotmb) imt, jmt, nt, ntmbb, reltim, avgper'
          write (io) stamp, iotext, expnam
          write (io) imt, jmt, nt, ntmbb, reltim, tmbavg*secday
c
          iotext = 'read (iotmb) (yt(j),j=1,jmt), (dyt(j),j=1,jmt)'
          write (io) stamp, iotext, expnam
	  write (io) (yt(j),j=1,jmt), (dyt(j),j=1,jmt)
c
          iotext = 
     &   'read (iotmb) (((tstor(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tstor, jmt*nt*(ntmbb+1))
c
          iotext =
     &    'read (iotmb) (((tdiv(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tdiv, jmt*nt*(ntmbb+1))
c
          iotext =
     &    'read (iotmb) (((tflux(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tflux, jmt*nt*(ntmbb+1))
c
          iotext = 
     &    'read (iotmb) (((tdif(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tdif, jmt*nt*(ntmbb+1))
c
          iotext = 
     &    'read (iotmb) (((tsorc(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tsorc, jmt*nt*(ntmbb+1))
c
          iotext = 
     &    'read (iotmb) ((smdvol(j,l),j=1,jmt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, smdvol, jmt*(ntmbb+1))
c
          call relunit (io)
        endif
c
c       zero quantities in preparation for next average
c
        numtmb = 0
        do mask=0,ntmbb
          do jrow=1,jmt
            smdvol(jrow,mask)  = c0
          enddo
	enddo
        do mask=0,ntmbb
          do m=1,nt
            do jrow=1,jmt
              tstor(jrow,m,mask) = c0
              tdiv(jrow,m,mask)  = c0
              tflux(jrow,m,mask) = c0
              tdif(jrow,m,mask)  = c0
              tsorc(jrow,m,mask) = c0
            enddo
          enddo
	enddo
      endif
8050  format (/1x,'Meridional Tracer Balance for tracer #',i2
     &,' averaged over depth and longitude',/1x
     &,' and averaged over a period of ', f10.2
     &,' days ending on ', a32//
     &,' jrow    lat      (T)t      (VT)y   DIFF_Tz     DIFF_Ty'
     &,'      source       error    volume'/)
8100  format (1x,i4,2x,f6.2,1x,7(1pe10.3,1x))
8200  format (1x,'global sum =',1x,7(1pe10.3,1x))
#endif
c
#ifdef show_zonal_mean_of_sbc
c
c-----------------------------------------------------------------------
c     print the zonal mean boundary conditions and related items
c-----------------------------------------------------------------------
c
      if (zmbcts) then
        cmmday = 10.0 * 86400.0
        cwatts = 4.186e4
        do jrow=1,jmt
          if (zmau(jrow) .ne. c0) then
            zmsmf(jrow,1) = zmsmf(jrow,1) / zmau(jrow)
            zmsmf(jrow,2) = zmsmf(jrow,2) / zmau(jrow)
            zmsm(jrow,1)  = zmsm(jrow,1)  / zmau(jrow)
            zmsm(jrow,2)  = zmsm(jrow,2)  / zmau(jrow)
          endif
          if (zmat(jrow) .ne. c0) then
            zmstf(jrow,1) = cwatts*zmstf(jrow,1)
            zmstf(jrow,2) = cmmday*zmstf(jrow,2)
            zmst(jrow,2)  = 1.0e3*zmst(jrow,2)
            do m=1,nt
              zmstf(jrow,m) = zmstf(jrow,m) / zmat(jrow)
              zmst(jrow,m)  = zmst(jrow,m)  / zmat(jrow)
	    enddo
          endif
	enddo
c
        if (iozmbc .eq. stdout .or. iozmbc .lt. 0) then
c
          write (stdout,'(///,27x,a,/)')
     &  'Z O N A L   M E A N   O F    S U R F A C E   B.  C.'
c
          write (stdout,10200)
          do jj=1,jmt
	    jrow = jmt + 1 - jj
            write (stdout,10201) jrow, yt(jrow), zmsmf(jrow,1)
     &,     zmsmf(jrow,2), zmstf(jrow,1), zmstf(jrow,2)
     &,     zmsm(jrow,1), zmsm(jrow,2), zmst(jrow,1), zmst(jrow,2) 
	  enddo
        endif
        if (iozmbc .ne. stdout .or. iozmbc .lt. 0) then
c
	  call getunit (io, 'zmean_sbc.dta','u s a ieee')
c
          write (stdout,*) ' => Zonal mean surface bc written'
     &,   ' unformatted to file zmean_sbc.dta on ts=',itt, ' ',stamp
          reltim = prelyr
c
          iotext = 'read (iozmbc) jmt, nt, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, nt, reltim
c
          iotext = 'read (iozmbc) (yt(j),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, yt, jmt)
c
          iotext = 'read (iozmbc) (yu(j),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, yu, jmt)
c
          iotext = 'read (iozmbc) ((zmsmf(j,n),j=1,jmt),n=1,2)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmsmf, jmt*2)
c
          iotext = 'read (iozmbc) ((zmstf(j,n),j=1,jmt),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmstf, jmt*nt)
c
          iotext = 'read (iozmbc) ((zmsm(j,n),j=1,jmt),n=1,2)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmsm, jmt*2)
c
          iotext = 'read (iozmbc) ((zmst(j,n),j=1,jmt),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmst, jmt*nt)
c
          call relunit (io)
        endif
c
      endif
10200 format(//
     &,26x,'Zonal Mean of Surface Boundary Condition diagnostic'
     &, //,1x,' row',1x,'  lat ', '   Taux    ','   Tauy    '
     &,' Heat Flux ',' Prec-Evap ',' Surface u ',' Surface v '
     &,'   SST     ','   SSS    ',/,12x
     &,' dyn/cm**2 ',' dyn/cm**2 ',' watts/m**2','   mm/day  '
     &,'   cm/sec  ','   cm/sec  ','   deg C   ','   ppt-35  ',/)
10201 format (1x, i4, 1x, f6.2, 8(1pe11.3))
#endif
c
#ifdef time_averages
c
c-----------------------------------------------------------------------
c     save the time mean averages on the "averaging" grid
c-----------------------------------------------------------------------
c
      if ((timats .or. eorun) .and. timavgint .gt. c0) then
# ifdef timing
        call toc ('diagnostic', 'output (diago)')
# endif
        call avgout
# ifdef timing
        call ticr ('diagnostic', 'output (diago)')
# endif
      endif
#endif
c
#ifdef xbts
c
c-----------------------------------------------------------------------
c     add in surface pressure gradients and missing external mode
c     parts of d/dt and implicit coriolis term for all XBTs
c-----------------------------------------------------------------------
c
# ifdef timing
      call toc ('diagnostic', 'output (diago)')
# endif
      call uxbt3
# ifdef timing
      call ticr ('diagnostic', 'output (diago)')
# endif
c
c-----------------------------------------------------------------------
c     save all time averaged XBT data at end of averaging period
c-----------------------------------------------------------------------
c
      if ((xbtts .or. eorun) .and. xbtint .gt. c0)then
# ifdef timing
        call toc ('diagnostic', 'output (diago)')
# endif
        call xbto
# ifdef timing
        call ticr ('diagnostic', 'output (diago)')
# endif
      endif
#endif
#ifdef timing
      call toc ('diagnostic', 'output (diago)')
#endif
c
      return
      end
