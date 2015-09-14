c======================= include file "calendar.h"======================
c
c                       calendar specification arrays
c
c-----------------------------------------------------------------------
c
c     eqyear = true to select a calendar in which each year
c              has the same number of days (i.e., no leap years)
c              false selects a julian calendar
c
c     eqmon  = true to force all months to have the same number of days
c              false => the usual 31, 28, 31, 30, ..., days per month. 
c              only used when eqyear = true 
c
c     dayname = character names of days
c
c     monname = character names of months
c
c     monlen = the length of each month (in days) when eqmon is true
c
c     yrlen  = the length of a typical (non-leap) year in days
c
c     daypm  = array of month lengths in days   (non-leap)
c
c     msum   = array of cumulative days preceding each month
c              (again, non-leap)
c
c-----------------------------------------------------------------------
      logical eqyear, eqmon
      integer daypm, msum, yrlen, monlen
c
      character*10 dayname
      character*12 monname
c
      common /calenc/ dayname(7), monname(12)
      common /calenl/ eqyear, eqmon
      common /caleni/ daypm(12), msum(12), yrlen, monlen

