c====================== include file "sponge.h" ========================
c
c     newtonian damping variables for sponge regions adjacent to
c     artificial southern and northern boundaries for use in limited 
c     domain basins. data must be prepared using the "sponge" routines 
c     included in PREP_DATA.
c
c     spngs  = coefficient for  damping T & S back to southern boundary
c     spngn  = coefficient for  damping T & S back to northern boundary
c
c     iprevd = pointer to previous month on disk
c              (the  month whose mid point has just been passed)
c     inextd = pointer to next month on disk
c              (the next month whose mid point hasn`t been reached yet)
c     iprev  = pointer to memory buffer for the previous month data
c              (the  month whose mid point has just been passed)
c     inext  = pointer to memory buffer for the next month data
c              (the next month whose mid point hasn`t been reached yet)
c     spbuf  = buffer for holding previous and next month disk data
c     spbuf  = sponge buffer data (imt,km,4,2) laid out as follows:
c              T(imt,km,1), S(imt,km,2) for the southern boundary
c              T(imt,km,3), S(imt,km,4) for the northern boundary
c     annlev = (t,f) = (replace seasonal data by annual means, use 
c                       seasonal data)
c     spgdpm = period in days for each monthly record
c     tspng  = time at midpoints of monthly records (days)
c     indxsp = index of dataset (needed for the interpolator)
c     readsp = (true,false) = (read, do not read) sponge data
c     wprev  = interpolation weight for previous month`s data
c
c     symbolically:
c     spdata(at time step) = (1-wprev)*spbuf(inext) + wprev*spbuf(iprev)
c
      common /cnewti/ inext, iprev, inextd, iprevd, indxsp  
      common /cnewt/ wprev, spgdpm(12), tspng(12)
      common /cnewt/ spbuf(imt,km,4,2)
      common /cnewt/ spngs(jmt), spngn(jmt)
c
c     array spbuf1 is only used locally to increase I/O efficiency
c     when reading and writing sponge data from disk.
c
      dimension spbuf1(imt,km,4)
c
      character*32 sstamp, stprev, stnext
      character*32 opt_sponge
      common /cnewtc/ sstamp, stprev, stnext, opt_sponge
      logical annlev, readsp
      common /cnewtl/ annlev, readsp
