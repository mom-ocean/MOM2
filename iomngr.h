c=====begin include "iomngr.h"==========================================
c
#include "stdunits.h"
c
c   arrays:
c     inuse        = does iomngr know unit number is currently in use
c     hide_file    = file which because of system quirks should not
c                    be closed by relunit
c     scratch_file = file should be deleted when released
c     unix_name    = file known to unix by name [exceptions: cray_ymp
c                    "word" and "sdsalloc" files]
c     ifile        = index into fname array corresponding to unit number
c     fname        = list of file names used since ioinit
c     iunit        = unit number corresponding to file name
c   scalars:
c     iohist       = unit number of io_history file
c     nfiles       = number of files names used so far
c
      parameter (maxunit = 99, maxfilenames = 200)
      logical inuse, hide_file, scratch_file, unix_name
      character * 64 fname
      common /iomngrl/ inuse(1:maxunit), hide_file(0:maxfilenames)
      common /iomngrl/ scratch_file(1:maxunit), unix_name(1:maxunit)
      common /iomngri/ ifile(1:maxunit)
      common /iomngrc/ fname(0:maxfilenames)
      common /iomngri/ iunit(0:maxfilenames), iohist, nfiles
c


