#ifdef time_varying_sbc_data
# if defined minimize_sbc_memory
      subroutine atmos (joff, js, je)
# else
      subroutine atmos
# endif
c
c-----------------------------------------------------------------------
c     this "atmos" model is an EXAMPLE of how to supply MOM with surface
c     boundary conditions which are a function of latitude, longitude,
c     and time (monthly data).
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
#include "param.h"
#include "csbc.h"
#include "ctdbc.h"
#include "tmngr.h"
#include "switch.h"
c
#if defined minimize_sbc_memory
c
c-----------------------------------------------------------------------
c     Work only on rows js through je in the memory window
c-----------------------------------------------------------------------
c
#else
c
c-----------------------------------------------------------------------
c     Work on all rows assuming S.B.C. arrays are (imt,jmt)
c-----------------------------------------------------------------------
c
      joff = 0
      js   = 1
      je   = jmt
#endif
c
c-----------------------------------------------------------------------
c     determine the disk pointers, time weight interpolation factor,
c     and whether or not it is time to bring in new S.B.C. from disk
c     based on the time (days) in MOM since dec 31, 1899 midnight.
c
c     express model time in days after start of S.B.C. by adding time
c     of I.C. to current model time then subtract time at start of
c     S.B.C.. Note that "itemptime" was allocated in settmngr and is
c     only needed as a temporary.
#if !defined minimize_sbc_memory
c     need to add "dt" to the model time because the call to 
c     atmos precedes the time stepping loop which calls mom, so the
c     model time has not yet been incremented when atmos executes.
#endif
c-----------------------------------------------------------------------
c
      if (joff .eq. 0) then
        do n=1,ntdbc
          call addtime (initial, imodeltime, itemptime)
#if !defined minimize_sbc_memory
            call addtime (itemptime, idt, itemptime)
#endif
          call subtime (itemptime, isbcstart(n), itemptime)
          daysbc(n) = realdays(itemptime)
#ifdef debug_atmos
          print *,' tp:n=',n,' daysbc(n)=',daysbc(n)
#endif
        enddo
      endif
c
c-----------------------------------------------------------------------
c     determine the disk pointers, time weight interpolation factor,
c     and whether or not it is time to bring in new S.B.C. from disk
c     based on the time (days) in MOM since dec 31, 1899 midnight.
c-----------------------------------------------------------------------
c
      if (joff .eq. 0) then
        do n=1,ntdbc
c
c         initialize the pointers (indices) to memory buffers and choose
c         interpolation method
c
          call timeinterp (daysbc(n), n, tdrec(1,n), aprec(1,n)
     &,      ntdrec(n), period(n), method, inextd(n), iprevd(n)
     &,      wprev(n), rdtdbc(n), inextm(n), iprevm(n))
        enddo
      endif
c
c-----------------------------------------------------------------------
c     read in data for each S.B.C. when necessary
c     indices for "obc" which are used to relate how S.B.C. are
c     ordered within "sbcocn":
c
c     1 is for hellerman taux
c     2 is for hellerman tauy
c     3 is for levitus sst
c     4 is for levitus sss
c-----------------------------------------------------------------------
c
      do n=1,ntdbc
#if defined minimize_sbc_memory
c
c       read in the next data record and previous data record 
c       from disk for rows js through je
c
c        call getunit (iotdbc(n), nameout(n), opt_tdbc)
        call read_obc_data (iotdbc(n), iprevd(n), rstamp1
     &,                     obc(1,js,n,iprevm(n)), js, je, joff)
        call read_obc_data (iotdbc(n), inextd(n), rstamp2
     &,                     obc(1,js,n,inextm(n)), js, je, joff)
c        call relunit (iotdbc(n))
#else
        if (rdtdbc(n)) then
c
c         read in the next data record from disk
c
          write (stdout,'(/tr1,a,i2,a,i6,a,a,a,i3,a,a,a,a32)') 
     &      '==> Reading S.B.C. #',n, ' from rec ',inextd(n)
     &,     '[ending at ', rstamp2,'] into index ',inextm(n),'  for '
     &,     nameout(n), ' at MOM time',stamp
c
          call getunit (iotdbc(n), nameout(n), opt_tdbc)
          call read_obc_data (iotdbc(n), inextd(n), rstamp2
     &,                       obc(1,js,n,inextm(n)), js, je, joff)
          call relunit (iotdbc(n))
        endif
#endif
      enddo
c
c-----------------------------------------------------------------------
c     time dependent boundary condition section:
c
c     interpolate the monthly climatological values brought in
c     from disk to the time step level for MOM
c-----------------------------------------------------------------------
c
c     index 1 is for hellerman taux
c     index 2 is for hellerman tauy
c     index 3 is for surface tracer flux #1 (using sst)
c     index 4 is for surface tracer flux #2 (using sss) 
c     index 5 is for ocean short wave (not used here)
c     index 6 is ocean model SST
c     index 7 is ocean model SSS
c
c
c     indices given by "mapsbc" relate how the S.B.C. are ordered
c     within array "sbcocn"
c
      itaux = mapsbc(1)
      itauy = mapsbc(2)
      ihflx = mapsbc(3)
      ipme  = mapsbc(4)
c      ishort= mapsbc(5)
      isst  = mapsbc(6)
      isss  = mapsbc(7)
c
      n = 1
      wnext = c1-wprev(n)
      do j=js,je
        jrow = joff + j
        do i=1,imt
          sbcocn(i,jrow,itaux)  =  wprev(n)*obc(i,j,n,iprevm(n))
     &                             +  wnext*obc(i,j,n,inextm(n))
        enddo
      enddo
c
      n = 2
      wnext = c1-wprev(n)
      do j=js,je
        jrow = joff + j
        do i=1,imt
          sbcocn(i,jrow,itauy)  =  wprev(n)*obc(i,j,n,iprevm(n))
     &                            +   wnext*obc(i,j,n,inextm(n))
        enddo
      enddo
c
#ifdef restorst
c
c-----------------------------------------------------------------------
c     set flux for restoring surface tracers.
c     n=(3,4) restores model (SST,SSS) to levitus (SST,SSS) on
c     timescale "dampts" using level thickness "dampdz"
c     note: set surface tracer flux to zero if not restoring. 
c-----------------------------------------------------------------------
c
      n = 3
      wnext = c1 - wprev(n)
      do j=js,je
        jrow = joff + j
        gamma = dampdz(1)*secday/dampts(1)
        do i=1,imt
          datax           = wprev(n)*obc(i,j,n,iprevm(n))
     &                      +  wnext*obc(i,j,n,inextm(n))
          sbcocn(i,jrow,ihflx) = gamma*(datax - sbcocn(i,jrow,isst))
        enddo
      enddo
c
      n = 4
      wnext = c1 - wprev(n)
      do j=js,je
        jrow = joff + j
        gamma = dampdz(2)*secday/dampts(2)
        do i=1,imt
          datax             = wprev(n)*obc(i,j,n,iprevm(n))
     &                        +  wnext*obc(i,j,n,inextm(n))
          sbcocn(i,jrow,ipme) = gamma*(datax - sbcocn(i,jrow,isss))
        enddo
      enddo
#else
      do j=js,je
        jrow = joff + j
        do i=1,imt
          sbcocn(i,jrow,ihflx)  = c0
          sbcocn(i,jrow,ipme)   = c0
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     zero out SST and SSS accumulators for the next ocean segment
c-----------------------------------------------------------------------
c
      if (numasbc .ne. 0 .and. eots) then
        do n=1,numasbc
          m= mapsbc(numosbc + n)
          do j=js,je
            jrow = joff + j
            do i=1,imt
              sbcocn(i,jrow,m) = c0
            enddo
          enddo
        enddo
      endif
c
      return
      end
#else
      subroutine nodataatmos
      return
      end
#endif
