      subroutine setatm (dtatm)
c
c-----------------------------------------------------------------------
c     "setatm" must do the following:
c
c         1) set up any needed atmospheric initializations 
c         2) set up the atmospheric boundary condition grids in degrees
c         3) define the atmospheric land/sea mask
c         4) set the atmospheric time step in seconds
c     
c     coded by:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
#include "stdunits.h"
c
#ifdef simple_sbc
      write (stdout,'(/,a,a,/)')
     &'==> Error: this "setatm (ATMOS)" is incompatible with'
     &,' ifdef "simple_sbc"'
      stop "=>setatm"
#else
      write (stdout,'(/,a,/,a/)')
     & '==> MOM S.B.C. are comming from: ATMOS which simulates an'
     &,'    atmospheric model with resolution different than MOM.'
#endif
#ifndef coupled
c
c     verify that "coupled" is  enabled.
c
      write (stdout,*)
     &' ==> Error: "coupled" should be enabled for ATMOS '
      stop '=>setatm'
#endif
c
c-----------------------------------------------------------------------
c     set up the atmospheric b.c. grids
c-----------------------------------------------------------------------
c
      call agrids
c
c-----------------------------------------------------------------------
c     define the atmospheric land/sea mask
c-----------------------------------------------------------------------
c
      call amask
c
c-----------------------------------------------------------------------
c     set the atmospheric time step in seconds
c-----------------------------------------------------------------------
c
      dtatm = 200.0
c
      return
      end

