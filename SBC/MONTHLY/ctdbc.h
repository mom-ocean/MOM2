c====================== include file "ctdbc.h" =========================
c
c           time dependent S.B.C. data for MOM 
c
#ifdef time_varying_sbc_data
c
c
c     ntdbc  = number of ocean S.B.C.
c               1 is for taux (dynes/cm**2)
c               2 is for tauy (dynes/cm**2)
c               3 is sst (converted to heat flux in "atmos.F")
c               4 is sss (converted to fresh water flux in "atmos.F")
c
c               for generality in data assimilation situations, each
c               S.B.C. is treated separately.  the default case assumes,
c               however, that all S.B.C. are monthly.
c
c
c     ntdrec = number of time dependent records on disk.
c               note: default case assumes 12 monthly S.B.C. records
c                     (defined at the centers of the months) from
c                     the MOM climatological database.
c
c     tdrec   = times at which the data records are defined.
c               times are assumed to define the centers of the records.
c               (ie: when using  monthly climatological data and the
c                first record is january, "tdrec(1)" is defined in days 
c                as day 15.5)
c               
c     aprec   = period lengths for the data records. (number of days
c               per month)
c
c     isbcstart= index to time at beginning of first S.B.C. data record
c
c     method   = interpolation method = (0,1,2,3) = (none, linear,
c                equal linear, equal area) see timeinterp for details. 
c
c     dstamp   = time stamps marking the end of each S.B.C. data record
c
c     period  = (true,false)  if the data records are to be treated as
c               (perodic, not periodic). if periodic, then the model
c               time is always mapped into the data records. if not,
c               then record 1 is used for all model time before the
c               beginning of the dataset and record "ntdrec" is used 
c               for all model time after the end of the dataset.
c
c     daysbc  = model time in days after start of S.B.C.
c
c     iprevd  = index for pointing to the nearest S.B.C. data on disk
c               which the current model time has passed.
c     inextd  = index for pointing to the nearest S.B.C. data on disk 
c               which the current model time has not yet reached.
c     iprevm  = index for pointing to the nearest S.B.C. data in memory
c               which the current model time has passed.
c     inextm  = index for pointing to the nearest S.B.C. data in memory 
c               which the current model time has not yet reached.
c     wprev   = interpolation factor for defining data at the current
c               time step (the weight attributed to the data from the 
c               "iprevd" time interval)
c               schematically the interpolation is defined at the time
c               step by:
c               data = wprev*data(iprevd) + (1-wprev)*data(inextd)
c     obc     = ocean S.B.C.
c     
c     iotdbc  = units for time dependent monthly boundary conditions
c               based on hellerman stress and levitus temperature and 
c               salt from the MOM database.
c     
c     coded by:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c 
      parameter (ntdbc = 4, maxbcr=12)
c
      character*32 rstamp1, rstamp2, dstamp
      character*40 namein, nameout
      character*64 opt_tdbc
      logical period, rdtdbc              
c
      common /ctdbci/ ntdrec(ntdbc), iprevd(ntdbc), inextd(ntdbc)
      common /ctdbci/ iprevm(ntdbc), inextm(ntdbc), iotdbc(ntdbc)
      common /ctdbci/ isbcstart(ntdbc), method
      common /ctdbcr/ wprev(ntdbc), aprec(maxbcr,ntdbc)
      common /ctdbcr/ tdrec(maxbcr,ntdbc), avgper(maxbcr), daysbc(ntdbc) 
      common /ctdbcl/ period(ntdbc), rdtdbc(ntdbc)
      common /ctdbcc/ rstamp1, rstamp2, dstamp(maxbcr,ntdbc)
      common /ctdbcc/ namein(ntdbc), nameout(ntdbc), opt_tdbc
# ifdef minimize_sbc_memory
      common /ctdbcr/ obc(imt,jmw,ntdbc,2)
# else
      common /ctdbcr/ obc(imt,jmt,ntdbc,2)
# endif
#endif
