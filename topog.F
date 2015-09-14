#ifdef drive_topog
# define driver_only
# include "grids.F"
# include "util.F"
# include "isleperim.F"
# include "iomngr.F"
# include "size_check.F"
# ifdef sgi_display
#  include "glstuff.F"
# endif


      program driver
c
c=======================================================================
c
c              T O P O G R A P H Y   M O D U L E
c
c     To generate a geometry and topography, look at the options in
c     the "run_topog" script, change them and execute the script.
c
c     topics:
c                  "idealized_kmt", "scripps_kmt", "rectangular_box"
c                  "flat bottom", importing, and exporting kmt file
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c                  
c     topics:
c                  delta.kmt file, violations, warnings
c                  interractions, grpahics, etc...
c     author:      c. h. goldberg        e-mail=> chg@gfdl.gov
c=======================================================================
c
#include "size.h"
#include "stdunits.h"
#include "coord.h"
#ifdef sgi_display
# include "glstuff.h"
#endif
      dimension kmt(imt,jmt), map(imt,jmt)
c
      call ioinit
c
#ifdef sgi_display
      call initdisplay
#endif

      write (stdout,'(//1x,a,//)') 'Construct the MOM grid'
c
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      maxlen = max(imt, jmt, km)
      call gcoord (maxlen, imt2, jmt2, km2
     &,            dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c     check that returned grid sizes match those in file "size.h"
c
      call size_check (imt2, jmt2, km2, 'driver', 'stop')
#ifdef sgi_display
      call scale_map (imainmap,'scale map', xu, yu, zw, imt2, jmt2, km2)
c     call display_grid ('rescaled grid', xu, yu, zw, imt2, jmt2, km2)
c     call waitmouse
#endif
c
c-----------------------------------------------------------------------
c     generate the kmt field on the grid defined by grids.F
c-----------------------------------------------------------------------
c
      write (stdout,'(/,20x,a/)') ' Driving the topog module...'
      call topog (kmt, map, xt, yt, zt, xu, yu, zw, imt2, jmt2, km2)
c
      call release_all
      stop
      end
#endif



      subroutine topog (kmt, map, xt, yt, zt, xu, yu, zw, imt2,jmt2,km2)
c
c-----------------------------------------------------------------------
c     construct the topographic mask "kmt" which determines the geometry
c     and topography by defining the number of vertical levels on model
c     "t" points.
c
c     The kmt field is constructed in the following order:
c
c      (A) Construct a kmt field using one of the following:
c        -Dscripps_kmt : interpolate Scripps 1 degree topography to
c                        "kmt" on MOM grid
c        -Dread_my_kmt : read in your own "kmt" from a file "kmt.dta"
c        -Didealized_kmt : generate an "idealized" not very realistic
c                          topography without reading the database
c        -Drectangular_box : sets "kmt" to a flat bottom rectangular box
c
c       Options:
c        -Dcyclic        : sets cyclic conditons in zonal direction
c        -Dflat_bottom   : sets "kmt" over ocean to maximum levels (km)
c
c      (B) Apply user modifications to kmt from the file
c          "delta.kmt.###x###x###.h", where ###x###x### is the grid size
c           "imt x jmt x km".  All changes made in this run will be
c           appended to the file. This file may be examined with an
c           editor but changes should not be made because column
c           alignment is important. If no such file is supplied, an
c           empty file is created.
c
c      (C) Limit minimum depth based on user-defined parameter "kmt_min"
c          "kmt_min" is the minimum permitted number of levels at an ocean
c          "t" cell.  kmt_min must be at least 2.
c       Options:
c        -Dfill_shallow: makes all "t" cells for which
c                        kmt(i,jrow) < kmt_min into land.  (default)
c        -Ddeepen_shallow: deepens such cell to depth kmt_min.
c        -Dround_shallow: changes depth to the nearer of 0 and
c                         zw(kmt_min)
c
c      (D) Fix perimeter violations (cells in the ocean perimeter of
c                                    two different land masses)
c       Options:
c        -Dfill_perimeter_violations : build a land bridge between
c                                      these land masses (default)
c        -Dwiden_perimeter_violations: widen straits between two
c                                      different land masses to a
c                                      minimum of two "t" cells
c
c      (E) Detect and handle remaining nonfatal problem cells
c        -Dfill_nonadvective : "t" cells with land at all four corners
c                               are filled to make them land. (default)
c        -Dleave_nonadvective : "t" cells with land at all four corners
c                               are allowed to remain ocean
c
c        -Dfill_potholes : "t" cells that are deeper than all four
c                          nearest neighbors are filled to the maximum
c                          depth of a nearest neighbor (default)
c        -Dleave_potholes : "t" cells that are deeper than all four
c                           nearest neighbors remain unchanged
c
c      (F) Interactive options
c        <at present, only available for SGI graphic workstations>
c        -Dprompt_shallow: allows interactive editing of shallow cells
c        -Dprompt_perimeter_violations : allows interactive editing of
c                                        cells flagged as perimeter
c                                        violations
c        -Dprompt_nonadvective : allows interactive editing of
c                                nonadvective cells
c        -Dprompt_potholes : allows interactive editing of "potholes"
c        -Dinteractive_edit_topography : activates unrestricted (free)
c                                      interactive editing of topography
c
c      (G) Final check of minimum depth and perimeter violations
c
c-----------------------------------------------------------------------
c
c     Programming note: kmt options are selected in 3 stages to allow
c     for maximum flexibility, control, and ease of selection.
c       (1) subroutine set_kmt_option:  presets up to 3 options
c       (2) subroutine accept_option: selects an option
c       (3) subroutine enter_kmt_changes: writes a group of kmt changes
c           to the file "delta.kmt.###x###x###.h"
c
c-----------------------------------------------------------------------
c     input:
c
c     xt  = longitudes of model grid "t" points (degrees)
c     yt  = latitudes of model grid "t" points (degrees)
c     zt  = depths of box "u/v" points on model grid (cm)
c     xu  = longitudes of model grid "u" points (degrees)
c     yu  = latitudes of model grid "u" points (degrees)
c     zw  = depths of box bottoms on model grid "t" boxes (cm)
c     imt = number of model longitudes
c     jmt = number of model latitudes
c     km  = maximum number of levels
c
c     a file "delta.kmt.###x###x###.h", where ###x###x### is the grid
c            size "imt x jmt x km", of user-modifications to the base
c            kmt field selected in part (A) above.
c
c     output:
c
c     kmt = number of vertical levels on "t" points
c           0 is a land point.
c           n implies n boxes between the surface and bottom
c             where (n >= kmt_min) & (n <= km)
c     map = array indicating land masses and their perimeters
c           0 is a mid ocean point
c           n > 0 is a land mass
c           n < 0 is the ocean perimeter of land mass |n|
c     file "delta.kmt.###x###x###.h" with changes appended.
c
c-----------------------------------------------------------------------
c
#include "stdunits.h"
#include "size.h"
#include "topog.h"
#ifdef sgi_display
# include "glstuff.h"
#endif
c
c     kmt_min is the minimum permitted number of levels at an ocean
c     "t" cell.  kmt_min must be at least 2.
c
      parameter (kmt_min = 2)
c
      dimension kmt(imt,jmt), xt(imt), yt(jmt), zt(km)
      dimension xu(imt), yu(jmt), zw(km)
c
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)

      character *32 title, title2
      character *32 del_kmt_filename, sdate
      character *80 del_kmt_line
      character *4 si, sj, sk
      character *5 sn_del_kmt
      character *1 sc
      character*80 iotext
      character*60 expnam
      character*32 stamp
c
c     the variable auto_kmt_changes tells whether any changes have been
c     made to the kmt field as a result of define options or interactive
c     user actions.  nchanges counts how many changes in each group.
c
      auto_kmt_changes = .false.
      nchanges = 0
c
c-----------------------------------------------------------------------
c     check that grid sizes in argument list and file "size.h" agree
c-----------------------------------------------------------------------
c
      call size_check (imt2, jmt2, km2, 'topog', 'stop')
c
      ncase = 0
c
      write (stdout,'(//,t37,a,/)')
     & 'T O P O G R A P H Y    G E N E R A T I O N'
      write (stdout,'(a,i4,/)') 'kmt_min = ', kmt_min
c
c     show options
c
      call topog_options
c
#if defined sgi_display && !defined drive_topog
      call initdisplay
      print *, 'imainmap=',imainmap
      call scale_map (imainmap,'scale map', xu, yu, zw, imt2, jmt2, km2)
#endif
c
#ifdef scripps_kmt
c
c-----------------------------------------------------------------------
c     construct model "kmt" based on Scripps 1 deg topography data.
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a/)')
     & ' Interpolating Scripps 1 degree topography to the MOM grid:'
c
      call announce_phase('scripps topog')
      call scripp (kmt, xt, yt, zt, xu, yu, zw, imt, jmt, km)
# ifdef sgi_display
c     call display_kmt ('kmt/Scripps', xu, yu, zw, kmt, imt, jmt, km)
c     call display_grid ('kmt/Scripps with grid', xu, yu, zw,imt,jmt,km)
c     call waitmouse
# endif

#endif
#ifdef read_my_kmt
c
      ncase = ncase + 1
      write (stdout,'(/a/)')
     & ' =>Reading a "kmt" field from file "kmt.dta" '
      call announce_phase('read my kmt')
      call getunit (io, 'kmt.dta', 'u s r ieee')
c
c     skip header record and read dimensional info
c
      read (io) stamp, iotext, expnam 
      read (io) imax, jmax, kmax
c
c     check that grid sizes in argument list and file "size.h" agree
c
      call size_check (imax, jmax, kmax, 'read kmt.dta', 'stop')
c
c     skip header record and read kmt field
c
      read (io) stamp, iotext, expnam
      read (io) kmt
      call relunit (io)
c
      nerror = 0
      do i=1,imt
        do j=1,jmt
          if (kmt(i,j) .lt. 0 .or. kmt(i,j) .gt. km) then
            nerror = nerror + 1
          end if
        end do
      end do
      if (nerror .gt. 0) then
        write (stdout,'(a,i4,a,a,i4,/)')
     &    '==>ERROR: There are ', nerror, ' depth values in the file '
     &,   '"kmt.dta" that are outside the range 0 <= kmt(i,j) <= km = '
     &,   km
        stop '==>topog'
      end if
#endif
#ifdef idealized_kmt
c
c-----------------------------------------------------------------------
c     construct a highly "idealized" world ... piece by piece
c
c     note: the purpose of this geometry/topography is to automatically
c           map into arbitrary resolution as grid dimensions "imt" and
c           "jmt" are changed, thereby facilitating the implementation
c           and verification of the model on various computer platforms
c           without referencing the topographic data base.  Although it
c           somewhat resembles the real world, it is NOT realistic.
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a/a/)')
     &' Constructing an idealized test case "kmt" which only crudely'
     &,' approximates world geometry. Topography is unrealistic.'
c
      call announce_phase('idealized topog')
      call idealized (kmt, xt, yt, zw, imt, jmt, km)
#endif
#ifdef rectangular_box
      ncase = ncase + 1
      write (stdout,'(/a,i4,a/)')
     & 'Constructing a flat bottomed rectangular box with ',km,' levels'
c
      call announce_phase('flat bottom box')
      do jrow=2,jmt-1
        do i=2,imt-1
          kmt(i,jrow) = km
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     check that only one method of generating topography was chosen
c-----------------------------------------------------------------------
c
      if (ncase .eq. 0) then
        write (stdout,'(/a/a/)')
     &  '=>Error: One of the following options must be enabled:'
     &, '      idealized_kmt, scripps_kmt, read_my_kmt, rectangular_box'
        stop '=>topog'
      elseif (ncase .gt. 1) then
        write (stdout,'(/a/a/)')
     &  '=>Error: Only one of the following options may be enabled:'
     &, '      idealized_kmt, scripps_kmt, read_my_kmt, rectangular_box'
        stop '=>topog'
      endif
c
#ifdef flat_bottom
      call announce_phase('flat bottom')
      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .ne. 0) then
            kmt(i,jrow) = km
          endif
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     apply boundary conditions
c-----------------------------------------------------------------------
c
      call kmtbc (kmt, imt, jmt)
c
      call area_volume (kmt, xu, yu, zw)
c
c-----------------------------------------------------------------------
c     read "delta.kmt" file of user-selected changes to kmt
c     filename = "delta.kmt.imt x jmt x km.h"
c-----------------------------------------------------------------------
c
      call announce_phase('delta kmt')
      call i2a(imt, si)
      call i2a(jmt, sj)
      call i2a(km,  sk)
      write (del_kmt_filename, '(a,a,a,a,a,a,a)')
     &      'delta.kmt.', si(1:len_trim(si)), 'x',
     &       sj(1:len_trim(sj)), 'x', sk(1:len_trim(sk)), '.h'
c
      call getunit(io_del_kmt, del_kmt_filename, 'fsr')
      write (stdout, '(/,a,a,/)')
     &     'Reading user modifications to kmt from file ',
     &     del_kmt_filename(1:len_trim(del_kmt_filename))
c
c     verify that "delta.kmt" file used same topography as current run
c
#ifdef scripps_kmt
        title = 'scripps topography'
#endif
#ifdef read_my_kmt
        title = 'kmt.dta topography'
#endif
#ifdef idealized_kmt
        title = 'idealized topography'
#endif
#ifdef rectangular_box
        title = 'rectangular box topography'
#endif
c
c     check file header
c
      n_del_kmt = 0
# if defined read_my_kmt
c
c     when importing a finished kmt field, do not read the delta file
c
#else
      read (io_del_kmt, '(t1,a1,t7,a32)', end=666) sc, title2
      if (title2 .ne. title) then
        print '(a)', '==>ERROR topog'
        print '(a,a,a,a,/,a,a)','file '
     &,       del_kmt_filename(1:len_trim(del_kmt_filename))
     &, ' has header: ', title2
     &, 'but current run uses topographic option: ',title
        stop '==>topog'
      end if
      print '(t1,a6,t7,a)', 'c     ', title2
c
c     read modifications to kmt field from file, ignoring comment lines
c
      do n=1,100000
        read (io_del_kmt, '(a80)', end=666) del_kmt_line
        if (del_kmt_line(1:1) .ne. 'c') then
          read (del_kmt_line,'(t1,a1,t7,tr4,i4,tr1,i4,t25,i4)', end=666)
     &                 sc, i, jrow, kmt_change
          n_del_kmt = n_del_kmt + 1
          kmt(i,jrow) = kmt_change
# ifdef show_topog_details
          print '(t7,a4,i4,a1,i4,a,t25,i4)',
     &          'kmt(',i,',',jrow,') = ', kmt_change
        else
          print '(a)', del_kmt_line
# endif
        end if
      end do
  666 continue
      if (n_del_kmt .gt. 0) then
        write (sn_del_kmt, '(i5)') n_del_kmt
      else
        sn_del_kmt = '   no'
      end if
      print '(/,a5,a,/)', sn_del_kmt, ' user modifications read'
# ifndef show_topog_details
      if (n_del_kmt .gt. 0) then
        print '(a/)', '   Use -Dshow_topog_details to see modifications'
      endif
# endif
      print '(a,a,/)', 'End of user modifications to kmt from file ',
     &     del_kmt_filename(1:len_trim(del_kmt_filename))
      call relunit(io_del_kmt)
#endif
c
      call area_volume (kmt, xu, yu, zw)
c
c-----------------------------------------------------------------------
c     apply boundary conditions
c-----------------------------------------------------------------------
c
      call kmtbc (kmt, imt, jmt)
c
#ifdef sgi_display
      call display_kmt('kmt + delta.kmt', xu, yu, zw, kmt, imt, jmt, km)
c     call waitmouse
      call display_grid ('kmt+delta.kmt ', xu, yu, zw,imt,jmt,km)
c     call waitmouse
#endif
c-----------------------------------------------------------------------
c     open file "delta.kmt.###x###x###.h" for additional changes to kmt
c-----------------------------------------------------------------------
c
#if defined read_my_kmt
c
c     when importing a finished kmt field, open a scratch file to
c     prevent writing over any existing delta file.
c
      call getunit(io_del_kmt, 'scratch_chgs', 'f s append')
#else
      call getunit(io_del_kmt, del_kmt_filename, 'f s append')
#endif
c
c     if this is a new file, write a header saying which kind of
c     base kmt field this delta.kmt file modifies
c
      if (n_del_kmt .eq. 0) then
        write (io_del_kmt, '(t1,a6,t7,a)') 'c     ', title
      end if
c-----------------------------------------------------------------------
c     limit the minimum number of levels. kmt_min should be >= 2
c-----------------------------------------------------------------------
c
      call announce_phase('min depth')
      call min_depth (kmt, kmt_min, xu, yu, zw)
#ifdef sgi_display
      call show_iter (0,'exit')
      call display_kmt ('kmt/shallow', xu, yu, zw, kmt, imt,jmt,km)
c     call waitmouse
#endif
c
c-----------------------------------------------------------------------
c     test for "island perimeter violations"
c     channels between different land masses must be at least 2 ocean
c     cells wide to prevent conflicting values of the stream function
c     on these island perimeter points.
c-----------------------------------------------------------------------
c
      call announce_phase('isle perimeters')
      write (stdout,'(/a/)')
     & 'Searching for and correcting island PERIMETER VIOLATIONS...'
      write (io_del_kmt, '(t1,a6,t7,a)') 'c     ',
     &                                   'isle perimeter changes'
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)

      call area_volume (kmt, xu, yu, zw)
#ifdef sgi_display
# ifdef prompt_perimeter_violations
      call clear_edit_window
# endif
c     call display_kmt ('kmt/perim', xu, yu, zw, kmt, imt, jmt, km)
c     call display_grid ('kmt/perim with grid', xu, yu, zw, imt,jmt,km)
c     call waitmouse
#endif
c
c-----------------------------------------------------------------------
c     detect isolated bays:
c     These are situations where an ocean point on the "t" grid is
c     surrounded by 4 land points on the "u" grid (the 4 verticies of
c     the "t" cell) which prevents lateral advection of tracers (T & S).
c     These points sometimes lead to numerical instabilities and can be
c     eliminated by changing the ocean "t" cell to a land value.
c
c     Isolated bays will be converted into land unless the directive
c     -Dleave_nonadvective is selected.
c
c     The directive -Dprompt_nonadvective enables interactive editing
c     of kmt at and in a neighborhood of a nonadvective "isolated bay"
c-----------------------------------------------------------------------
c
      call announce_phase('nonadvective')
      call nonadvective(kmt,map,iperm,jperm,iofs,nippts,nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
c-----------------------------------------------------------------------
c     detect potholes:
c     These are situations where an ocean point on the "t" grid is
c     deeper than its 4 surrounding neighbor cells, which prevents
c     lateral and vertical advection of tracers at this depth.
c     These points usually cause no problems, but they also do not
c     contribute to ocean circulation.  They may be eliminated by
c     changing kmt at this "t" grid location to the maximum of kmt on
c     the 4 neighboring locations.
c
c     Potholes will be converted into land unless the directive
c     -Dleave_potholes is selected.
c
c     The directive -Dprompt_potholes enables interactive editing
c     of kmt at and in a neighborhood of "potholes"
c-----------------------------------------------------------------------
c
      call announce_phase('potholes')
      call potholes (kmt,map,iperm,jperm,iofs,nippts,nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
#ifdef interactive_edit_topography
c
c-----------------------------------------------------------------------
c     free interactive editing of topography
c-----------------------------------------------------------------------
c
      call announce_phase('free editing')
      write (stdout, '(t1,a6,t7,a)') 'c     ',
     &                                   'free interactive editing'
      write (io_del_kmt, '(t1,a6,t7,a)') 'c     ',
     &                                   'free interactive editing'
      call clear_kmt_options ()
      i_opt = 2
      call interactive ('edit', 2, i_opt
     &,           imt/2, jmt/2, kmt, imt/2, jmt/2
     &,           xu, yu, zw)
      call accept_option(i_opt, kmt)
c
      call enter_kmt_changes()
c
      call area_volume (kmt, xu, yu, zw)
#endif
c
c-----------------------------------------------------------------------
c     ALL USER MODIFICATIONS TO TOPOG MUST PRECEDE THESE FINAL CHECKS
c-----------------------------------------------------------------------
c
c     limit the minimum number of levels. kmt_min should be >= 2
c
      call announce_phase('min depth')
      call min_depth (kmt, kmt_min, xu, yu, zw)
#ifdef sgi_display
c     call display_kmt ('kmt/shallow', xu, yu, zw, kmt, imt,jmt,km)
c     call waitmouse
#endif
c
c-----------------------------------------------------------------------
c     apply boundary conditions
c-----------------------------------------------------------------------
c
      call kmtbc (kmt, imt, jmt)
c
c-----------------------------------------------------------------------
c     Finally, retest for "island perimeter violations" that might have
c     been introduced by define options or interactive editing.
c     Channels between different land masses must be at least 2 ocean
c     cells wide to prevent conflicting values of the stream function
c     on these island perimeter points.
c-----------------------------------------------------------------------
c
      print '(a,/)'
     &,     '==>Final recheck of island perimeters'
c
      call announce_phase('isle perimeters')
      write (io_del_kmt, '(t1,a6,t7,a)') 'c     ',
     &                                   'isle perimeter changes'
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
#  ifndef skip_island_map
      print '(/a/)'
     &,   'Use -Dskip_island_map to not show the following map'
      call showmap (map, imt, jmt, linewidth)
#  else
      print '(/a/)'
     &,   'Remove -Dskip_island_map to show the island perimeter map'
#  endif
      call area_volume (kmt, xu, yu, zw)
# ifdef sgi_display
      call clear_edit_window
      call display_kmt ('kmt/from topog', xu, yu, zw, kmt, imt, jmt, km)
      call display_depths ('depths', xu, yu, zw, imt, jmt, km)
c     call display_grid ('kmt/from topog', xu, yu, zw, imt,jmt,km)
c     call waitmouse
# endif
c
      call relunit(io_del_kmt)
c
c-----------------------------------------------------------------------
c     show the kmt field
c-----------------------------------------------------------------------
c
#  ifndef skip_kmt_map
      print '(/a/)'
     &,   'Use -Dskip_kmt_map to not show the following map'
      write (stdout,'(/,20x,a/)') ' The "kmt" field follows:'
      call iplot (kmt, imt, imt, jmt)
#  else
      print '(/a/)'
     &,   'Remove -Dskip_kmt_map to show the kmt map'
#  endif
c
c---------------------------------------------------------------------
c     compute a topography checksum
c---------------------------------------------------------------------
c
      cksum = 0.0
      do jrow=1,jmt
        do i=1,imt
          cksum = cksum + i*jrow*kmt(i,jrow)
        enddo
      enddo
      write (stdout,*) ' "kmt" checksum = ', cksum
#if defined read_my_kmt
      write (stdout,'(/a/)') '==> The kmt field has been imported.'
#else
      write (stdout,'(/a,a/a,a/)')
     & '==> CAUTION: If any changes have been made to the interpolated'
     &,' base kmt field.'
     &,'             they can be viewed in file: '
     &,  del_kmt_filename(1:len_trim(del_kmt_filename))
#endif
#if defined write_my_kmt
      write (stdout,'(/a/)')  '==> The kmt field has been exported.'
#endif
c
#ifdef write_my_kmt
c
c---------------------------------------------------------------------
c     export the finished kmt field
c---------------------------------------------------------------------
c
      write (stdout,'(/a/)')
     & ' =>Writing a "kmt" field to file "kmt.dta" '
      call getunit (io, 'kmt.dta', 'u s r ieee')
c
      stamp = ' no stamp'
      iotext = 'read (io) imax, jmax, kmax'
      expnam = '  '
      write (io) stamp, iotext, expnam
      write (io) imt, jmt, km
      iotext = 'read (io) ((kmt(i,j),i=1,imax),j=1,jmax)'
      write (io) stamp, iotext, expnam
      write (io) kmt
c
      call relunit (io)
#endif
#if !defined drive_topog && !defined permit_runtime_topog_changes
c
c     Automatic changes to kmt are always permitted when topog is run
c     in stand-alone mode, but are normally not permitted when topog
c     is called from the model if there is a nonempty "delta.kmt" file.
c
c     It is assumed that the "delta.kmt" file creates the user''s choice
c     of topography, and that any additional changes made during a model
c     run of MOM are probably errors.
c
      if (auto_kmt_changes .and. n_del_kmt .ne. 0) then
        write (stdout,'(/a/a/a/a/)')
     & ' =>Error: Automatic changes have been made to "kmt".'
     &,'          This may indicate a problem with the "kmt" field.'
     &,'          Using option -Dpermit_runtime_topog_changes will'
     &,'          accept the changes and allow the run to continue.'
        stop 'topog'
      endif
#endif
#if defined sgi_display && !defined drive_topog
      call make_extra_map_windows
      print *, 'tp: imap2 =', imap2
      call scale_map (imap2, 'map2', xu, yu, zw, imt2, jmt2, km2)
      call scale_map (imap3, 'map3', xu, yu, zw, imt2, jmt2, km2)
#endif
      return
      end



      subroutine topog_options
      parameter (maxifdefs=30)
      character*(30) ifdefs(maxifdefs)
#include "stdunits.h"
        write (stdout,'(/,10x,a,/)')
     &  ' A V A I L A B L E    O P T I O N S   F O R   T O P O G'
#ifdef show_available_options
        write (stdout,9000)
        write (stdout,9100)
        write (stdout,9200)
c
#else
        write (stdout,'(/a/)')
     & 'Use option "show_available_options" to see them'
#endif
        write (stdout,'(/,a)')
     &   ' THE FOLLOWING "IFDEF" OPTIONS HAVE BEEN ENABLED (TURNED ON):'
        n = 0
#ifdef drive_topog
        n = n + 1
        ifdefs(n) = 'drive_topog'
#endif
#ifdef scripps_kmt
        n = n + 1
        ifdefs(n) = 'scripps_kmt'
#endif
#ifdef read_my_kmt
        n = n + 1
        ifdefs(n) = 'read_my_kmt'
#endif
#ifdef idealized_kmt
        n = n + 1
        ifdefs(n) = 'idealized_kmt'
#endif
#ifdef rectangular_box
        n = n + 1
        ifdefs(n) = 'rectangular_box'
#endif
#ifdef flat_bottom
        n = n + 1
        ifdefs(n) = 'flat_bottom'
#endif
#ifdef cyclic
        n = n + 1
        ifdefs(n) = 'cyclic'
#endif
#ifdef fill_shallow
        n = n + 1
        ifdefs(n) = 'fill_shallow'
#endif
#ifdef deepen_shallow
        n = n + 1
        ifdefs(n) = 'deepen_shallow'
#endif
#ifdef round_shallow
        n = n + 1
        ifdefs(n) = 'round_shallow'
#endif
#ifdef prompt_shallow
        n = n + 1
        ifdefs(n) = 'prompt_shallow'
#endif
#ifdef fill_perimeter_violations
        n = n + 1
        ifdefs(n) = 'fill_perimeter_violations'
#endif
#ifdef widen_perimeter_violations
        n = n + 1
        ifdefs(n) = 'widen_perimeter_violations'
#endif
#ifdef prompt_perimeter_violations
        n = n + 1
        ifdefs(n) = 'prompt_perimeter_violations'
#endif
#ifdef fill_nonadvective
        n = n + 1
        ifdefs(n) = 'fill_nonadvective'
#endif
#ifdef leave_nonadvective
        n = n + 1
        ifdefs(n) = 'leave_nonadvective'
#endif
#ifdef prompt_nonadvective
        n = n + 1
        ifdefs(n) = 'prompt_nonadvective'
#endif
#ifdef fill_potholes
        n = n + 1
        ifdefs(n) = 'fill_potholes'
#endif
#ifdef leave_potholes
        n = n + 1
        ifdefs(n) = 'leave_potholes'
#endif
#ifdef prompt_potholes
        n = n + 1
        ifdefs(n) = 'prompt_potholes'
#endif
#ifdef interactive_edit_topography
        n = n + 1
        ifdefs(n) = 'interactive_edit_topography'
#endif
      if (n .gt. maxifdefs) then
        write (stdout,*) n
     &, '=>Error: increase "maxifdefs" to ',n,' in "topog_options"'
        stop
      endif
      nacros = 4
      nn = n/float(nacros) + 1
      do i=1,nn
        ns = 1+nacros*(i-1)
        ne = ns + nacros-1
        if (ne .gt. n) ne = n
        if (ns .le. n) write (stdout,9101) (ifdefs(n1),n1=ns,ne)
      enddo
      write (stdout,'(/a/)')  ' END OF ENABLED OPTIONS LIST'
      return
9101  format(4(1x,a30))
#ifdef show_available_options
9000  format(/
     & //' "kmt" defines the number of grid cells from the surface'
     &  /'  to bottom of ocean for each latitude and longitude point'
     & //' drive_topog  = constructs topography in a stand alone mode'
     &  /'               (without MOM). use run_topog script to execute'
     &  /' scripps_kmt  = interpolate Scripps 1 deg topography to "kmt"'
     &  /'                on grid defined by grids.F module'
     &  /' read_my_kmt  = read in your own "kmt" from a file "kmt.dta'
     &  /' idealized_kmt= generate a crude (not very realistic) "kmt"'
     &  /'                without using the DATABASE'
     &  /' rectangular_box = generate a flat bottomed "kmt". Enabeling'
     &  /'                "cyclic" will open the E-W sides to make a'
     &  /'                re-entrant channel'
     &  /' flat_bottom = sets "kmt" over ocean to maximum levels (km)'
     &  /' cyclic      = sets an E-W cyclic condition on "kmt"'
     &   )
9100  format(/
     &  /' Non-Interactive options'
     & //' A shallow place is one where kmt(i,j) < "kmt_min" '
     &  /' fill_shallow   = fills shallow regions to land (kmt(i,j)=0)'
     &  /' deepen_shallow = deepens shallow regions to "kmt_min"'
     &  /' round_shallow  = rounds shallow regions to nearer of zero'
     &  /'                  or zw(kmt_min)'
     &///' Perimeter violations occur when two distinct land masses'
     &  /' are separated by only one ocean cell'
     &  /' fill_perimeter_violations = build a land bridge between'
     &  /'                             these two land masses'
     &  /'                             (default)         '
     &  /' widen_perimeter_violations= widen straits between two land'
     &  /'                             masses to a minimum of two '
     &  /'                             "t" cells'
     &///' Non-advective cells have land at all four corners'
     &  /' and are not influenced by advection'
     &  /' fill_nonadvective = make them land. (default)'
     &  /' leave_nonadvective = leave them as ocean'
     &///' Potholes are "t" cells deeper than all 4 nearest neighbors'
     &  /' fill_potholes = fill to the maximum depth of a nearest'
     &  /'                  neighbor (default)'
     &  /' leave_potholes = leave them alone'
     &   )
9200  format(/
     &  /' Interactive options (only at GFDL for SGI workstations)'
     & //' prompt_shallow = allows interactive editing of regions'
     &  /'                  shallower than "kmt_min"'
     &  /' prompt_perimeter_violations = allows interactive editing of'
     &  /'                               cells flagged as perimeter'
     &  /'                               violations'
     &  /' prompt_nonadvective = allows interactive editing of'
     &  /'                       nonadvective cells'
     &  /' prompt_potholes = allows interactive editing of "potholes"'
     &  /' interactive_edit_topography = allows free interactive'
     &  /'                               editing of "kmt"'
     &   )
#endif
      end





      subroutine area_volume(kmt, xu, yu, zw)
c
c-----------------------------------------------------------------------
c     compute surface area and volume of ocean ("t" cells and "u" cells)
c     (note that areas are defined at each level)
c-----------------------------------------------------------------------
c
#include "size.h"
      dimension kmt(imt, jmt), xu(imt), yu(jmt), zw(km)
      dimension tcella(km)
      dimension cst(jmt)
      dimension dxt(imt), dyt(jmt)
c
      do k=1,km
        tcella(k) = 0
      enddo
      ocnp   = 0
      tcellv = 0
      pi = 4.0 * atan(1.0)
      degrad = pi / 180.0
      radius = 4.0e9 / (2.0 * pi)
      do i=2,imt-1
        dxt(i) = radius * degrad * (xu(i) - xu(i-1))
      end do
      do jrow = 2,jmt-1
        cst(jrow) = cos(0.5*(yu(jrow)+yu(jrow-1))*degrad)
        dyt(jrow) = radius * degrad * (yu(jrow) - yu(jrow-1))
      end do
c
      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) +
     &           cst(jrow)*dxt(i)*dyt(jrow)
            enddo
            tcellv = tcellv + cst(jrow)*dxt(i)*dyt(jrow)*zw(kmt(i,jrow))
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
        enddo
      enddo
c
      print '(/)'
      print '(a,f8.0)',   'number of ocean t cells =', ocnp
      print '(a,e15.6,a)','surface area (t cells)  =',tcella(1),' cm**2'
      print '(a,e15.6,a)', 'ocean volume (t cells)  =', tcellv,' cm**3'
      print '(/)'
c
      return
      end



      subroutine kmtbc (kmt, imt, jmt)
c
c-----------------------------------------------------------------------
c     set lateral boundary conditions on kmt
c     set cyclic, solid wall and symmetry conditions on "t" grid
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt)
c
#ifdef cyclic
c
c     set cyclic conditions on eastern and western boundary
c
      do jrow=1,jmt
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
      enddo
#else
c
c     set solid wall conditions on eastern and western boundary
c
      do jrow=1,jmt
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
      enddo
#endif
#ifdef symmetry
c
c     set symmetry conditions at "u" row jmt-1
c
      do i=1,imt
        kmt(i,jmt) = kmt(i,jmt-1)
      enddo
#else
c
c     close the domain at the northernmost row
c
      do i=1,imt
        kmt(i,jmt) = 0
      enddo
#endif
c
c     close the domain at the southernmost row
c
      do i=1,imt
        kmt(i,1)   = 0
      enddo
      return
      end



#ifdef scripps_kmt
      subroutine scripp (kmt, xt, yt, zt, xu, yu, zw
     &,                  imt2, jmt2, km2)
c
c-----------------------------------------------------------------------
c     construct "kmt" from scripps 1 deg topography.
c-----------------------------------------------------------------------
c
# include "stdunits.h"
      parameter (iw=362, jw=182)
      parameter (lenw=20*iw)
# include "size.h"
c
      dimension kmt(imt2,jmt2), xt(imt2), yt(jmt2), zt(km2)
      dimension xu(imt2), yu(jmt2), zw(km2)
      dimension a(imt,jmt)
c
c     file 'scripps.top' has real *4 ieee format data
c
      real*4 lons4, lats4, cx4, cy4, d4
      dimension d(iw,jw), cx(iw), cy(jw), dt(iw), cxt(iw), work(lenw)
      dimension d4(iw,jw), cx4(iw), cy4(jw)
c
c     check that grid sizes in argument list and file "size.h" agree
c
      call size_check (imt2, jmt2, km2, 'scripps', 'stop')
c
c-----------------------------------------------------------------------
c     read  Scripps 1 deg topography data
c-----------------------------------------------------------------------
c
      call getunit (io_scripps, 'scripps.top','usr ieee')
      read (io_scripps)
      read (io_scripps) lons4, lats4, cx4, cy4, d4
      call relunit (io_scripps)
      lons = lons4
      lats = lats4
      do i=1,iw
        cx(i) = cx4(i)
      end do
      do j=1,jw
        cy(j) = cy4(j)
      end do
      do i=1,iw
        do j=1,jw
          d(i,j) = d4(i,j)
        end do
      end do
      write (stdout,'(/,a,/)')
     & ' => Preparing model "kmt" from Scripps 1 deg Topography'
c
c-----------------------------------------------------------------------
c     translate Scripps grid longitudes and data so model grid lies
c     entirely within the Scripps grid.
c-----------------------------------------------------------------------
c
      call tranlon (d, iw, iw-2, jw, dt, cx, xt(2), imt-2, cxt)
c
#ifndef skip_translation_details
      write (stdout,'(/a/)')
     & 'Use -Dskip_translation_details to not show the following'
      write(stdout,'(///,50x, a30,/)') 'Scripps 1 deg topography grid'
      write (stdout,'(//,1x,a30,/)') '            latitudes:'
      write (stdout,'(1x,10f10.5)') cy
      write (stdout,'(//,1x,a30,/)') '  original longitudes:'
      write (stdout,'(1x,10f10.5)') cx
      write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
      write (stdout,'(1x,10f10.5)') cxt
#else
      write (stdout,'(/a/)')
     & 'Remove -Dskip_translation_details to show translation details'
#endif
c
c-----------------------------------------------------------------------
c     interpolate Scripps topography to MOM "t" grid points
c     intrp = 1 => model grid is coarser than 1 deg Scripps data
c     intrp = 2 => model grid is finer than 1 deg Scripps data
c-----------------------------------------------------------------------
c
c     USER INPUT: set "intrp" appropriately if this test fails
c
      avgdx = (xu(imt-1) - xu(1))/(imt-2)
      avgdy = (yu(jmt-1) - yu(1))/(jmt-2)
      if (avgdx .gt. 1.0 .and. avgdy .gt. 1.0) then
        intrp = 1
      else
        intrp = 2
      end if
c     END USER INPUT
      if (intrp .eq. 1) then
        print '(/a/a/)'
     &,' "intrp=1" assumes MOM resolution is coarser than Scripps 1 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ftc (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      else if (intrp .eq. 2) then
        print '(/a/a/)'
     &,' "intrp=2" assumes MOM resolution is finer than Scripps 1 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ctf (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      endif
c
      do jrow=1,jmt
# ifdef cyclic
        a(1,jrow)   = a(imt-1,jrow)
        a(imt,jrow) = a(2,jrow)
# else
        a(1,jrow)   = a(2,jrow)
        a(imt,jrow) = a(imt-1,jrow)
# endif
      enddo
c
c-----------------------------------------------------------------------
c     convert depths from meters to cm and change sign.
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
          if (a(i,jrow) .lt. 0.0) then
            a(i,jrow) = -a(i,jrow)*100.0
          else
            a(i,jrow) = 0.0
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     set topographic modifications on model grid here
c-----------------------------------------------------------------------
c
c     set artificial walls at northern & southern boundaries
c
      do i=1,imt
        a(i,1 )  = 0.0
        a(i,jmt) = 0.0
      enddo
c
c-----------------------------------------------------------------------
c     discretize topography "kmt" to nearest "zw" levels
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=2,imt-1
          kmt(i,jrow) = 0
          if (a(i,jrow) .ne. 0.0) then
            kmt(i,jrow) = indp (a(i,jrow), zw, km)
          endif
        enddo
# ifdef cyclic
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
# else
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
# endif
      enddo
c
      return
      end

#endif



#ifdef idealized_kmt
      subroutine idealized (kmt, xt, yt, zw, imt, jmt, km)
c
c-----------------------------------------------------------------------
c     construct a highly "idealized" world ... piece by piece
c
c     note: the purpose of this geometry/topography is to automatically
c           map into arbitrary resolution as grid dimensions "imt" and
c           "jmt" are changed, thereby facilitating the implementation
c           and verification of the model on various computer platforms
c           without referencing the topographic data base.  Although it
c           somewhat resembles the real world, it is NOT realistic.
c
c     author: r. c. pacanowski            e-mail==>rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt), xt(imt), yt(jmt), zw(km)
c
      do jrow=2,jmt-1
        do i=2,imt-1
          kmt(i,jrow) = km
        enddo
      enddo
c
      do jrow=1,jmt
        kmt(imt,jrow) = kmt(2,jrow)
        kmt(1,jrow)   = kmt(imt-1,jrow)
      enddo
      do i=1,imt
        kmt(i,jmt) = 0
        kmt(i,1)   = 0
      enddo
c
c     antarctica
c
      call setkmt (kmt, xt, yt, imt, jmt, -90.0, 0.0, 360.0, -80.0, 0.0
     &,            360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 360.0-25.0, 360.0
     &,            -70.0, 360.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 0.0, 360.0, -70.0, 0.0
     &,            170.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 360.0-135.0, 360.0-60.0
     &,            -68.0, 360.0-75.0, 360.0-60.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -70.0, 0.0, 155.0, -67.0, 50.0
     &,            145.0, 0)
c
c     australia
c
      call setkmt (kmt, xt, yt, imt, jmt, -35.0, 116.0, 120.0, -31.0
     &,            114.0, 130.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -38.0, 140.0, 151.0, -31.0
     &,            130.0, 151.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -31.0, 115.0, 153.0, -20.0
     &,            113.0, 149.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 113.0, 149.0, -11.0
     &,            131.0, 143.0, 0)
c
c     south america
c
      call setkmt (kmt, xt, yt, imt, jmt, -50.0, 360.0-74.0, 360.0-68.0
     &,            -40.0, 360.0-73.0, 360.0-62.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -40.0, 360.0-73.0, 360.0-62.0
     &,            -20.0, 360.0-70.0, 360.0-40.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 360.0-70.0, 360.0-40.0
     &,            -16.0, 360.0-81.0, 360.0-35.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -16.0, 360.0-81.0, 360.0-35.0
     &,            0.0, 360.0-80.0, 360.0-50.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 0.0, 360.0-80.0, 360.0-50.0
     &,            11.0, 360.0-75.0, 360.0-60.0, 0)
c
c     central america
c
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 360.0-78.0, 360.0-75.0
     &,            20.0, 360.0-105.0, 360.0-97.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 20.0, 360.0-105.0, 360.0-97.0
     &,            30.0, 360.0-115.0, 360.0-94.0, 0)
c
c     north america
c
      call setkmt (kmt, xt, yt, imt, jmt, 25.0, 360.0-82.0, 360.0-80.0
     &,            30.0, 360.0-85.0, 360.0-81.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 30.0, 360.0-115.0, 360.0-80.0
     &,            40.0, 360.0-124.0, 360.0-74.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 360.0-124.0, 360.0-74.0
     &,            50.0, 360.0-124.0, 360.0-57.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 50.0, 360.0-124.0, 360.0-57.0
     &,            60.0, 360.0-140.0, 360.0-64.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 360.0-165.0, 360.0-64.0
     &,            65.0, 360.0-140.0, 360.0-64.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 65.0, 360.0-140.0, 360.0-64.0
     &,            70.0, 360.0-162.0, 360.0-72.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 360.0-162.0, 360.0-140.0
     &,            72.0, 360.0-157.0, 360.0-157.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 360.0-130.0, 360.0-70.0
     &,            75.0, 360.0-120.0, 360.0-80.0, 0)
c
c     greenland
c
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 360.0-45.0, 360.0-45.0
     &,            75.0, 360.0-58.0, 360.0-19.0, 0)
c
c     africa
c
      call setkmt (kmt, xt, yt, imt, jmt, -35.0, 19.0, 28.0, 6.0, 8.0
     &,            50.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 0.0, 50.0, 18.0, 0.0
     &,            56.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 18.0, 0.0, 56.0, 26.0, 0.0
     &,            59.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 360.0-10.0, 360.0, 18.0
     &,            360.0-18.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 18.0, 360.0-18.0, 360.0, 26.0
     &,            360.0-15.0, 360.0, 0)
c
c     northern africa &  europe & asia
c
      call setkmt (kmt, xt, yt, imt, jmt, 26.0, 360.0-15.0, 360.0, 40.0
     &,            360.0-7.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 360.0-7.0, 360.0, 50.0
     &,            360.0, 360.0, 0)
c
      call setkmt (kmt, xt, yt, imt, jmt, 8.0, 77.0, 78.0, 26.0, 65.0
     &,            90.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 4.0, 99.0, 100.0, 26.0, 90.0
     &,            115.0, 0)
c
      call setkmt (kmt, xt, yt, imt, jmt, 26.0, 0.0, 126.0, 40.0, 0.0
     &,            122.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 0.0, 130.0, 50.0, 0.0
     &,            140.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 50.0, 0.0, 140.0, 60.0, 8.0
     &,            140.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 8.0, 163.0, 65.0, 13.0
     &,            180.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 65.0, 13.0, 188.0, 70.0, 20.0
     &,            180.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 70.0, 180.0, 75.0, 90.0
     &,            100.0, 0)
c
c     add an "idealized" undulating topography
c
      bot = zw(km)
      pi  = 4.0*atan(1.0)
      do jrow=2,jmt
        do i=2,imt-1
          if (kmt(i,jrow) .ne. 0) then
            arg = bot*(1-0.4*abs(cos((jrow*pi)/jmt)*sin((i*2*pi)/imt)))
            kmt(i,jrow) = indp (arg, zw, km)
          endif
        enddo
      enddo
c
c     add an "idealized" ridges
c
      level = indp (0.666*zw(km), zw, km)
c
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 360.0-20.0, 360.0-10.0
     &,            30.0, 360.0-45.0, 360.0-35.0, level)
      call setkmt (kmt, xt, yt, imt, jmt, 30.0, 360.0-45.0, 360.0-35.0
     &,            60.0, 360.0-20.0,  360.0-30.0, level)
      call setkmt (kmt, xt, yt, imt, jmt, -60.0,360.0-100.0, 360.0-130.0
     &,            40.0, 360.0-160.0, 180.0, level)
c
      level = indp (0.5*zw(km), zw, km)
c
      call setkmt (kmt, xt, yt, imt, jmt, -50.0, 360.0-120.0
     &,            360.0-120.0, 30.0, 190.0, 190.0, level)
      return
      end



      subroutine setkmt (kmt, xt, yt, imt, jmt, alat1, slon1, elon1
     &,                  alat2, slon2, elon2, num)
c
c-----------------------------------------------------------------------
c     set the topography mask "kmt(i,j)" = "num" within the area of
c     the parallelogram bounded by vertices:
c     (alat1,slon1), (alat1,elon1), (alat2,slon2), & (alat2,elon2)
c
c     inputs:
c
c     xt = longitudes of "t" points in degrees
c     yt = latitudes of "t" points in degrees
c     imt = number of model longitudes
c     jmt = number of model latitudes
c     alat1 = southern latitude of parallelogram (degrees)
c     slon1 = starting longitude of southern edge of parallelogram (deg)
c     elon1 = ending longitude of southern edge of parallelogram (deg)
c     alat2 = northern latitude of parallelogram (degrees)
c     slon2 = starting longitude of notthern edge of parallelogram (deg)
c     elon2 = ending longitude of notthern edge of parallelogram (deg)
c     num   = number of vertical levels
c
c     outputs:
c
c     kmt   = mask of vertical levels on model "t" points
c
c     author: r. c. pacanowski            e-mail==>rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt), xt(imt), yt(jmt)
c
c     convert the four vertices into model indices
c     (js,is1), (js,ie1), (je,is2), (je,ie2)
c
      j1 = indp (alat1, yt, jmt)
      j2 = indp (alat2, yt, jmt)
      js = min (j1,j2)
      je = max (j1,j2)
c
      i1  = indp (slon1, xt, imt)
      i2  = indp (elon1, xt, imt)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)
c
      i1  = indp (slon2, xt, imt)
      i2  = indp (elon2, xt, imt)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)
c
      is = is1
      ie = ie1
c
c     fill in the area bounded by (js,is1), (js,ie1), (je,is2), (je,ie2)
c     the nudging of 1.e-5 is to insure that the test case resolution
c     generates the same topography and geometry on various computers.
c
      c1 = 1.0
      if (js .eq. je) then
        rdj = c1
      else
        rdj = c1/(je-js)
      endif
      do jrow=js,je
        do i=is,ie
          kmt(i,jrow) = num
        enddo
        is = nint(rdj*((jrow-js)*is2 + (je-jrow)*is1) + 1.0e-5)
        ie = nint(rdj*((jrow-js)*ie2 + (je-jrow)*ie1) + 1.0e-5)
      enddo
      return
      end
#endif






      subroutine clear_kmt_options ()
c
c-----------------------------------------------------------------------
c     clear all potential options for kmt changes
c-----------------------------------------------------------------------
c
#include "topog.h"
c
      do i_opt=1,max_opt
        do j_opt=1,len_opt
          do k_opt=1,4
            kmt_opt(i_opt, j_opt, k_opt) = 0
          end do
        end do
      end do
      return
      end


      subroutine kmt_option (i_opt, i, jrow, kmt_old, kmt_new, kmt)
c
c-----------------------------------------------------------------------
c     add a change for kmt(i,jrow) from kmt_old to kmt_new to the
c     options list for option i_opt
c-----------------------------------------------------------------------
c
#include "topog.h"
#include "size.h"
      dimension kmt(imt,jmt)
c
      if (kmt_new .ne. kmt_old) then
c
c       if current option is full, accept it and clear it
c
        if (kmt_opt(i_opt,len_opt,1) .ne. 0) then
          call accept_option(i_opt, kmt)
          do n=1,len_opt
            do k_opt=1,4
              kmt_opt(i_opt,n,k_opt) = 0
            end do
          end do
        end if
c
c      enter the change in kmt_opt    
c
       j_opt = len_opt
       do nn=len_opt,1,-1
          if (kmt_opt(i_opt,nn,1) .eq. 0) then
            j_opt = nn
          end if
        end do
        kmt_opt(i_opt,j_opt,1) = i
        kmt_opt(i_opt,j_opt,2) = jrow
        kmt_opt(i_opt,j_opt,3) = kmt_old
        kmt_opt(i_opt,j_opt,4) = kmt_new
      end if
#ifdef debug_interactive
      print '(a,6i5)'
     &,     '    kmt_option:',i_opt,j_opt,i,jrow,kmt_old,kmt_new
#endif
      return
      end


      subroutine select_option (problem
     &,                         i, jrow, kmt, xu, yu, zw)
c
c-----------------------------------------------------------------------
c     use define options to select a default option for changing kmt,
c     then accept and enter the changes.
c     if prompting is enabled, enter interactive editing to accept or
c     modify the default option for each kind of editing of kmt.
c-----------------------------------------------------------------------
c
#include "size.h"
#include "topog.h"
      dimension kmt(imt, jmt)
      character *(*) problem
      dimension xu(imt), yu(jmt), zw(km)
c
c     select default option based on ifdef options
c     then if corresponding "prompt" option is selected, allow user to
c     accept or reject default, or make other modifications to kmt
c
      if (problem .eq. 'perim') then
#if defined fill_perimeter_violation || !defined widen_perimeter_violations
        idefault = 1
#else
        idefault = 3
#endif
        i_opt = idefault
#if defined prompt_perimeter_violations
        call interactive ('perim', idefault, i_opt
     &,                   i, jrow, kmt, i, jrow
     &,                   xu, yu, zw)
#else
        call accept_option(i_opt, kmt)
#endif
c
      else if (problem .eq. 'bay') then
#if defined fill_nonadvective || !defined leave_nonadvective
        idefault = 1
#else
        idefault = 3
#endif
#if !defined drive_topog && !defined permit_runtime_topog_changes
        if (n_del_kmt .ne. 0) then
          idefault = 3
        endif
#endif
        i_opt = idefault
#if defined prompt_nonadvective
        call interactive ('bay', idefault, i_opt
     &,                   i, jrow, kmt, i, jrow
     &,                   xu, yu, zw)
#else
        call accept_option(i_opt, kmt)
#endif
c
       else if (problem .eq. 'pothole') then
#if defined fill_potholes || !defined leave_potholes
        idefault = 1
#else
        idefault = 3
#endif
#if !defined drive_topog && !defined permit_runtime_topog_changes
        if (n_del_kmt .ne. 0) then
          idefault = 3
        endif
#endif
        i_opt = idefault
#if defined prompt_potholes
        call interactive ('pothole', idefault, i_opt
     &,                   i, jrow, kmt, i, jrow
     &,                   xu, yu, zw)
#else
        call accept_option(i_opt, kmt)
#endif
c
       else if (problem .eq. 'shallow') then
#if defined round_shallow ||  \
    !(defined fill_shallow || defined deepen_shallow)
        idefault = 2
#else
# if defined fill_shallow
        idefault = 1
# else
#  if defined deepen_shallow
        idefault = 3
#  endif
# endif
#endif
#if !defined drive_topog && !defined permit_runtime_topog_changes \
    && !defined round_shallow
        if (n_del_kmt .ne. 0) then
c
c         clear option 2
c
          i_opt = 2
          do j_opt=1,len_opt
            do k_opt=1,4
              kmt_opt(i_opt,j_opt,k_opt) = 0
            end do
          end do
          idefault = 2
        endif
#endif
        i_opt = idefault
#if defined prompt_shallow
        call interactive ('shallow', idefault, i_opt
     &,                   i, jrow, kmt, i, jrow
     &,                   xu, yu, zw)
#else
        call accept_option(i_opt, kmt)
#endif
      end if
c
      return
      end



      subroutine accept_option(i_opt,kmt)
c
c-----------------------------------------------------------------------
c     copy selected option to kmt and kmt_changes
c-----------------------------------------------------------------------
c
#include "size.h"
#include "topog.h"
      dimension kmt(imt, jmt)
c
#ifdef debug_interactive
      print '(/,a,)', 'accepting option'
#endif
      do j_opt=1,len_opt
        if (kmt_opt(i_opt, j_opt, 1) .ne. 0 .and.
     &      kmt_opt(i_opt, j_opt, 2) .ne. 0) then
          ii = kmt_opt(i_opt,j_opt,1)
          jj = kmt_opt(i_opt,j_opt,2)
          kk = kmt_opt(i_opt,j_opt,4)
          kmt(ii,jj) = kk
#ifdef debug_interactive
          print '(2i3,a,i4,a,i4,a,i4)'
     &,         i_opt, j_opt,'  kmt(',ii,',',jj,') = ', kk
#endif
          nn = nchanges + 1
          if (nn .gt. max_change) then
            print '(/,a,a,/,a,a,i6,a)'
     &,           'WARNING: kmt_changes buffer full.  '
     &,           'Changes written early to delta.kmt file.'
     &,           'To avoid this message, '
     &,            'increase max_change = ', max_change,' in topog.h'
            call enter_kmt_changes()
          end if
          nchanges = nchanges + 1
          do jch=1,4
            kmt_changes(nchanges,jch) = kmt_opt(i_opt, j_opt, jch)
          end do
          kmt(kmt_changes(nchanges,1), kmt_changes(nchanges,2)) =
     &      kmt_changes(nchanges,4)
        end if
      end do
c
      return
      end



      subroutine undo_option(i_opt,kmt)
c
c-----------------------------------------------------------------------
c     copy selected option to kmt and kmt_changes exchanging new and
c     old values for kmt.
c-----------------------------------------------------------------------
c
#include "size.h"
#include "topog.h"
      dimension kmt(imt, jmt)
c
      do j_opt=1,len_opt
        if (kmt_opt(i_opt, j_opt, 1) .ne. 0 .and.
     &      kmt_opt(i_opt, j_opt, 2) .ne. 0) then
          ii = kmt_opt(i_opt,j_opt,1)
          jj = kmt_opt(i_opt,j_opt,2)
          kk = kmt_opt(i_opt,j_opt,3)
          kmt(ii,jj) = kk
          nn = nchanges + 1
          if (nn .gt. max_change) then
            print '(/,a,a,/,a,a,i6,a)'
     &,           'WARNING: kmt_changes buffer full.  '
     &,           'Changes written early to delta.kmt file.'
     &,           'To avoid this message, '
     &,            'increase max_change = ', max_change,' in topog.h'
            call enter_kmt_changes()
          end if
          nchanges = nchanges + 1
          do jch=1,4
            if (jch .le. 2) then
              kmt_changes(nchanges,jch) = kmt_opt(i_opt, j_opt, jch)
            else
              kmt_changes(nchanges,jch) = kmt_opt(i_opt, j_opt, 7-jch)
            end if
          end do
        end if
      end do
c
      return
      end



      subroutine enter_kmt_changes()
c
c-----------------------------------------------------------------------
c     copy accumulated changes from the array kmt_changes
c     to kmt and to the file "delta.kmt.###x###x###.h"
c-----------------------------------------------------------------------
c
#include "size.h"
#include "topog.h"
      dimension kmt(imt,jmt)

      print '(/,a,/)', 'Modifications to kmt field'
      if (nchanges .eq. 0) then
        print '(t4,a,/)', '<none>'
      else
        auto_kmt_changes = .true.
      end if
      do n=1,nchanges
        kmt(kmt_changes(n,1), kmt_changes(n,2)) = kmt_changes(n,4)
        write (io_del_kmt, '(t7,a4,i4,a1,i4,a,t25,i4)')
     &      'kmt(',kmt_changes(n,1),',',kmt_changes(n,2),') = '
     &,      kmt_changes(n,4)
#ifdef show_topog_details
        print '(t7,a4,i4,a1,i4,a,a,i4,a,i4)'
     &,      'kmt(',kmt_changes(n,1),',',kmt_changes(n,2),') changed '
     &,      'from ',kmt_changes(n,3),' to ',kmt_changes(n,4)
#else
        if (n .eq. 1) then
          print '(a)',' Use -Dshow_topog_details to see modifications'
        endif
#endif
      end do
      print '(/)'
      nchanges = 0
      return
      end



#if defined prompt_potholes || defined prompt_nonadvective \
 || defined prompt_perimeter_violations || defined prompt_shallow \
 || defined interactive_edit_topography
c
      subroutine interactive (problem, idefault, i_opt
     &,           icen2, jcen2, kmt, ied2, jed2
     &,           xu, yu, zw)
c
c-----------------------------------------------------------------------
c     interactively selects an option to fix a particular "problem"
c     cell.  user may also at this time freely edit any number of
c     other "t" cells selected from either the edit window or the
c     main map window
c-----------------------------------------------------------------------
c
# ifdef sgi_display
c   Interactive options:
c       click left mouse on problem cell to select option 1
c                           (usually default option)
c       click right mouse on problem cell to select option 3
c                            (usually alternative option)
c       click middle mouse to select a cell in either the edit window or
c                  the main map window for free editing.
c   Free edit options:
c        left mouse on cell makes cell one level shallower
c        right mouse on cell deepens cell by one level
c        middle mouse chooses new cell without changes
c        same three options on "mouse button window" change level faster
c        any mouse button on color bar selects that depth
c        any mouse button on "end free edit" ends free editing
c
# endif
c     input:
c       problem : problem type: 'shallow', 'perim', 'bay', 'pothole'.
c                               or 'edit'
c       ied2, jed2 : coordinates of problem "t" cell
c       icen2, jcen2 : coordinates of center of edit window
c       idefault : default for this kind of problem selected by user
c                  define options
c       i_opt : tentative option selected so far
c       kmt_opt : collection of pre-defined options for fixing problem
c
c     output:
c       i_opt : user-selected option number
c       kmt_changes : changes made to kmt field
c       nchanges : number of such changes
c
c     author: c. h. goldberg              e-mail==>chg@gfdl.gov
c-----------------------------------------------------------------------
c
# include "size.h"
# include "topog.h"
      dimension kmt(imt, jmt)

      character *(*) problem
      character *(32) title
      dimension xu(imt), yu(jmt), zw(km)
      character choice
      character*(32) locstr
      integer iedit, jedit, ieditwidth, jeditwidth
      character*(20) formatstr

# if defined sgi_display
#  include "glstuff.h"
c
c     copy input arguments so that changes made in this subroutine do
c     not cause changes in the input arguments
c
      icen = icen2
      jcen = jcen2
      ied  = ied2
      jed  = jed2
c
c     initial interactive edit state
c
  100 continue
#  ifdef debug_interactive
      print '(a,a,a,2i4)',
     &   'interactive editing for ',problem(1:len_trim(problem)),
     &   ' at', ied, jed
#  endif
c
c     "Back door" test for early exit from interactive editing.
c     Putting the cursor in the iteration window leaves all nonadvective
c     cells and potholes.
c
      if (problem .ne. 'perim' .and. icurrentwin.eq.iwin(iterwin)) then
        i_opt = 3
        return
      end if
c
      if (problem .eq. 'shallow') then
        call prompt ('fill', 'free edit', 'deepen', 'ignore', 'exit')
      else if (problem .eq. 'perim') then
        call prompt ('fill', 'free edit', 'widen', 'ignore', ' ')
      else if (problem .eq. 'bay') then
        call prompt ('fill', 'free edit', 'leave', 'ignore', 'exit')
      else if (problem .eq. 'pothole') then
        call prompt ('fill', 'free edit', 'leave', 'ignore', 'exit')
      else if (problem .eq. 'edit') then
        call prompt (' ', 'free edit', ' ', 'end free edit', 'exit')
      else
        call prompt ('default', 'free edit', 'alternative',
     &               'ignore', 'exit')
      end if
      title = problem
      call display_zoom (title, xu, yu, zw, kmt
     &,                  icen,jcen,ied,jed
     &,                  imt, jmt, km)
      call waitmouse
c
  101 continue
      call wait_cursor(iwinnr, x, y,
     &                 iselect, jselect, xu, yu,
     &                 imt, jmt, km)
c
c     process user action in quit window
c
      if (iwinnr .eq. iquit) then
c
c       ignore problem and exit interactive editing
c
        i_opt = 2
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
        call clear_edit_window
        goto 999
c
c     start of free edit session unrelated to "problem" cell
c
      else if (problem .eq. 'edit') then
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
        if (iwinnr .eq. imainmap .or. iwinnr .eq. izoom) then
          goto 400
        else
          goto 300
        end if
c
c     clicking left mouse button on center box or left third of
c     "mouse button" window selects option 1
c
      else if (leftmouse .and.
     &    ((iselect .eq. ied) .and. (jselect .eq. jed) .or.
     &     (iwinnr .eq. iprompt) .and. (nint(x) .eq. 1))) then
        i_opt = 1
        call accept_option(i_opt, kmt)
        goto 210
c
c     clicking right mouse button on center box or right third of
c     "mouse button" window selects option 3
c
      else if (rightmouse .and.
     &         ((iselect .eq. ied) .and. (jselect .eq. jed) .or.
     &          (iwinnr .eq. iprompt) .and. (nint(x) .eq. 3))) then
        i_opt = 3
        call accept_option(i_opt, kmt)
        goto 230
c
c     middle mouse selects user-specified "free" editing
c     in this editing mode: leftmouse   reduces kmt
c                           rightmouse  increases kmt
c
      else if (middlemouse) then
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
        if ((iwinnr .eq. iprompt) .and. (nint(x) .eq. 2)) then
          goto 300
        else if (iwinnr .eq. imainmap .or. iwinnr .eq. izoom) then
          goto 400
        else
          call waitmouse
          goto 101
        end if
#  ifdef debug_interactive
c
c     test for early exit from interactive editing
c
      else if (problem.ne.'perim'.and.icurrentwin.eq.iwin(iterwin))then
        i_opt = 3
        return
#  endif
      else
c
c       cursor did not select valid alternative at main menu
c
        call waitmouse
        goto 101
      end if
c
c     The next statement should not be reached
c
      print *,'Error in subroutine interactive just above 210'
      stop '==>interactive'
c
c-----------------------------------------------------------------------
c
  210 continue
c
c     redisplay zoom window and confirm choice
c
      call prompt ('accept & next', 'accept & edit', 'undo', ' ', ' ')
      title = 'confirm'
      call display_zoom (title, xu, yu, zw, kmt
     &,                  icen,jcen,ied,jed
     &,                  imt, jmt, km)
      call waitmouse

  211 continue
      call wait_cursor(iwinnr, x, y,
     &                 iselect, jselect, xu, yu,
     &                 imt, jmt, km)
c
      if (leftmouse .and.
     &  ((iselect .eq. ied) .and. (jselect .eq. jed) .or.
     &   (iwinnr .eq. iprompt) .and. (nint(x) .eq. 1))) then
c
c       accept option 1 and exit
c
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' accepted'
#  endif
        goto 999
c
      else if (rightmouse .and.
     &  ((iselect .eq. ied) .and. (jselect .eq. jed) .or.
     &   (iwinnr .eq. iprompt) .and. (nint(x) .eq. 3))) then
c
c       undo option 1
c
        call undo_option(i_opt, kmt)
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' undone'
#  endif
        goto 100
c
      else if (middlemouse) then
c
c       accept option 1 and enter free editing
c
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' accepted--enter free edit'
#  endif
        if ((iwinnr .eq. iprompt) .and. (nint(x) .eq. 2)) then
          goto 300
        else if (iwinnr .eq. imainmap .or. iwinnr .eq. izoom) then
          goto 400
        end if
c
      else
c
c       cursor did not select a valid option at confirm option 1
c
        call waitmouse
        goto 211
      end if
c-----------------------------------------------------------------------
c
  230 continue
c
c     redisplay zoom window
c
      call prompt ('undo', 'accept & edit', 'accept & next', ' ', ' ')
      title = 'confirm'
      call display_zoom (title, xu, yu, zw, kmt
     &,                  icen,jcen,ied,jed
     &,                  imt, jmt, km)
      call waitmouse
c
  231 continue
      call wait_cursor(iwinnr, x, y,
     &                 iselect, jselect, xu, yu,
     &                 imt, jmt, km)
c
      if (rightmouse .and.
     &      ((iselect .eq. ied) .and. (jselect .eq. jed) .or.
     &        (iwinnr .eq. iprompt) .and. (nint(x) .eq. 3))) then
c
c       accept option 3 and exit
c
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                              icen, jcen, .false.)
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' accepted'
#  endif
        goto 999
c
      else if (leftmouse .and.
     &         ((iselect .eq. ied) .and. (jselect .eq. jed) .or.
     &          (iwinnr .eq. iprompt) .and. (nint(x) .eq. 1))) then
c
c       undo option 3
c
        call undo_option(i_opt, kmt)
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' undone'
#  endif
        goto 100
c
      else if (middlemouse) then
c
c       accept option 3 and enter free editing
c
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' accepted--enter free edit'
#  endif
        if ((iwinnr .eq. iprompt) .and. (nint(x) .eq. 2)) then
          goto 300
        else if (iwinnr .eq. imainmap .or. iwinnr .eq. izoom) then
          goto 400
        end if
c
      else
c
c       cursor did not select a valid option at confirm option 3
c
        call waitmouse
        goto 231
      end if
c
c-----------------------------------------------------------------------
c
  300 continue
c
c     "free edit" mode, but edit cell not yet selected
c
c     make sure last button push was in main map or edit/zoom window
c
      call prompt(' ','select edit location',' ','end free edit','exit')
c
  301 continue
c
c     process user action in quit window
c
      if (iwinnr .eq. iquit) then
c
c       ignore problem and exit interactive editing
c
        i_opt = 2
        call clear_edit_window
        goto 999
      end if
c
c     process user action in iteration window
c
      if (iwinnr .eq. iterwin) then
c
c       Immediate and final exit from topog!  delta.kmt file will have
c       all changes made so far, but topography is not guaranteed to be
c       free of violations.  After this exit, a new run of topog will
c       read all changes from the delta.kmt file and continue the
c       current run.  NO check will be made to insure that define
c       options have not changed between runs.
c
        i_opt = 2
        call clear_edit_window
        call enter_kmt_changes()
c
        call area_volume (kmt, xu, yu, zw)
        print '(/,a,/,a,a,)', 'WARNING:  USER TERMINATION OF TOPOG',
     &              '          Topography not guaranteed to fulfill',
     &              ' model requirements'
        stop ('==>topog')
      end if
c
      if (iwinnr .ne. imainmap .and. iwinnr .ne. izoom) then
        call waitmouse
        call find_cursor(iwinnr, x, y,
     &                   iselect, jselect, xu, yu,
     &                   imt, jmt, km)
        goto 301
      end if
      goto 400
c
c-----------------------------------------------------------------------
c
  400 continue
c
c     "free edit" mode
c     edit cell and edit window centering cell are implicitly selected
c     by last mouse click.  display edit window and get user action.
c
c     clear option 2
c
      i_opt = 2
      do j_opt=1,len_opt
        do k_opt=1,4
          kmt_opt(i_opt,j_opt,k_opt) = 0
        end do
      end do
c
  401 continue
c
c     change the centering of the edit window only if edit location is
c     selected from main map.
c
      if (iwinnr .eq. imainmap) then
        icen = iselect
        jcen = jselect
      end if
      iedit = iselect
      jedit = jselect
      kmt_orig = kmt(iedit,jedit)
#  ifdef debug_interactive
      print *, 'tp: user edit at (', iedit, jedit,')'
#  endif
      call prompt ('shallower','new edit cell','deeper',
     &             'end free edit',' ')
      title = 'edit:'
      call display_zoom (title, xu, yu, zw, kmt
     &,                  icen,jcen,iedit,jedit
     &,                  imt, jmt, km)
      call waitmouse
c
  402 continue
c
c     edit depth of kmt(iedit,jedit) directed by user mouse clicks
c
      call wait_cursor(iwinnr, x, y,
     &                 iselect, jselect, xu, yu,
     &                 imt, jmt, km)
c
c     mouse click in main map window
c
      if (iwinnr .eq. imainmap) then
        if (iselect .ne. iedit .or. jselect .ne. jedit) then
c
c         accept changes to previous edit cell and select new edit cell
c
          if (kmt(iedit,jedit) .ne. kmt_orig) then
            call kmt_option(i_opt,iedit,jedit,kmt_orig,kmt(iedit,jedit)
     &,                     kmt)
          end if
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.false.)
          call redisplay_main_zoom (xu, yu, zw, kmt,
     &                              icen, jcen, .false.)
          iedit = iselect
          jedit = jselect
          kmt_orig = kmt(iedit,jedit)
#  ifdef debug_interactive
          print *, 'tp: new main map user edit at (', iedit, jedit,')'
#  endif
        end if
        call accept_option(i_opt, kmt)
        goto 400
c
c     mouse click in edit window
c
      else if (iwinnr .eq. izoom) then
        if (iselect .ne. iedit .or. jselect .ne. jedit) then
c
c         accept changes to previous edit cell and select new edit cell
c
          if (kmt(iedit,jedit) .ne. kmt_orig) then
            call kmt_option(i_opt,iedit,jedit,kmt_orig,kmt(iedit,jedit)
     &,                     kmt)
          end if
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.false.)
          call redisplay_main_zoom (xu, yu, zw, kmt,
     &                              icen, jcen, .true.)
          iedit = iselect
          jedit = jselect
          kmt_orig = kmt(iedit,jedit)
#  ifdef debug_interactive
          print *, 'tp: new user edit at (', iedit, jedit,')'
#  endif
        end if
c
c       process user action in edit window
c
        if (leftmouse) then
c
c         decrease kmt(iedit,jedit) by 1
c
          kmt(iedit,jedit) = max (kmt(iedit,jedit)-1, 0)
          print '(a,2i4,a,i4)',
     &          'd) decreased kmt(',iedit,jedit,') to ',
     &          kmt(iedit,jedit)
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.true.)
          call waitmouse
          goto 402
c
        else if (rightmouse) then
c
c         increase kmt(iedit,jedit) by 1
c
          kmt(iedit,jedit) = min (kmt(iedit,jedit)+1, km)
          print '(a,2i4,a,i4)',
     &          'd) increased kmt(',iedit,jedit,') to ',
     &          kmt(iedit,jedit)
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.true.)
          call waitmouse
          goto 402
c
        else if (middlemouse) then
c
c         accept changes to edit cell
c
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.true.)
          if (kmt(iedit,jedit) .ne. kmt_orig) then
            call kmt_option(i_opt,iedit,jedit,kmt_orig,kmt(iedit,jedit)
     &,                     kmt)
          end if
          call redisplay_main_zoom (xu, yu, zw, kmt,
     &                              icen, jcen, .false.)
c
c         be sure not to erase option 2 which has accumulated changes
c         made in this edit window
c
          goto 401
c
        else
          goto 402
        end if
c
c     process user action in prompt window
c
      else if (iwinnr .eq. iprompt) then
        if (leftmouse .and. nint(x) .eq. 1) then
c
c         decrease kmt(iedit,jedit)
c
          kmt(iedit,jedit) = nint(0.5*kmt(iedit,jedit)-0.001)
          print '(a,2i4,a,i4)',
     &          'd) decreased kmt(',iedit,jedit,') to ',
     &          kmt(iedit,jedit)
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.true.)
          call waitmouse
          goto 402
        else if (rightmouse .and. nint(x) .eq. 3) then
c
c         increase kmt(iedit,jedit)
c
          kmt(iedit,jedit) = nint(0.5*(kmt(iedit,jedit)+km))
          print '(a,2i4,a,i4)',
     &          'd) increased kmt(',iedit,jedit,') to ',
     &          kmt(iedit,jedit)
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.true.)
          call waitmouse
          goto 402
        else if (middlemouse .and. nint(x) .eq. 2) then
c
c         accept changes to edit cell
c
          call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.true.)
          if (kmt(iedit,jedit) .ne. kmt_orig) then
            call kmt_option(i_opt,iedit,jedit,kmt_orig,kmt(iedit,jedit)
     &,                     kmt)
          end if
          call redisplay_main_zoom (xu, yu, zw, kmt,
     &                              icen, jcen, .false.)
          goto 300
        else
          call waitmouse
          goto 402
        end if
c
c     process user action in color bar window
c
      else if (iwinnr .eq. icolorbar) then
#  if defined kmt_colors || !defined depth_colors
        kmt(iedit,jedit) = nint((1.0 - x/maxscolor)*km)
#  else
        kmt(iedit,jedit) = 0
        z = (1.0 - x/maxscolor) * zw(km)
        diff = z
        do k=1,km
          if (abs(z-zw(k)) .lt. diff) then
            kmt(iedit,jedit) = k
            diff = abs(z-zw(k))
          end if
        end do
#  endif
#  ifdef debug_interactive
        print *,'tp:  colorbar:  x=',x,'  kmt-->', kmt(iedit,jedit)
#  endif
        print '(a,2i4,a,i4)',
     &        'd)   changed kmt(',iedit,jedit,') to ',
     &        kmt(iedit,jedit)
        call redisplay_edit_cell (kmt,iedit,jedit,xu,yu,zw,.true.)
        call waitmouse
        goto 402
c
c     process user action in quit window
c
      else if (iwinnr .eq. iquit) then
c
c       accept changes to current edit cell and end free edit
c
        if (kmt(iedit,jedit) .ne. kmt_orig) then
          call kmt_option(i_opt,iedit,jedit,kmt_orig,kmt(iedit,jedit)
     &,                   kmt)
        end if
        call redisplay_main_zoom (xu, yu, zw, kmt,
     &                            icen, jcen, .false.)
        call clear_edit_window
        goto 999
c
      else
c
c       invalid mouse click option selected: retry
c
        call waitmouse
        goto 402
      end if
c
c     The next statement should not be reached
c
      print *,'Error in interactive just above 999'
      stop '==>interactive'
c
c-----------------------------------------------------------------------
c
  999 continue
      call clear_edit_window
      call accept_option(i_opt, kmt)
#  ifdef debug_interactive
      print '(a,i1,a)','option ',i_opt,
     &      ' accepted. leaving interactive edit'
#  endif
c
# endif



#ifdef ascii_display
c
c-----------------------------------------------------------------------
c     Call the ASCII version of interactive
c-----------------------------------------------------------------------
c
      call ascii_interactive(problem,idefault,i_opt,icen2,jcen2,
     &                       kmt,ied2,jed2,xu,yu,zw)
#endif
      return
      end



      subroutine prompt(opt1, opt2, opt3, quit_opt, exit_opt)
c
c-----------------------------------------------------------------------
c     show user options in interactive editing of topography
c-----------------------------------------------------------------------
c
      character *(*) opt1, opt2, opt3, quit_opt, exit_opt

# if (defined debug_interactive) && !(defined ascii_display)
      print *,'tp: interactive options: (l)',opt1(1:len_trim(opt1)),
     &                              '   (m)',opt2(1:len_trim(opt2)),
     &                              '   (r)',opt3(1:len_trim(opt3))
# endif
# if defined sgi_display
      call gl_prompt(opt1, opt2, opt3, quit_opt, exit_opt)
# endif
#ifdef ascii_display
      call ascii_prompt(opt1,opt2,opt3,quit_opt,exit_opt)
#endif
      return
      end
#endif



      subroutine nonadvective(kmt,map,iperm,jperm,iofs,nippts,nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
c-----------------------------------------------------------------------
c     detect isolated bays:
c     These are situations where an ocean point on the "t" grid is
c     surrounded by 4 land points on the "u" grid (the 4 verticies of
c     the "t" cell) which prevents lateral advection of tracers (T & S).
c     These points sometimes lead to numerical instabilities and can be
c     eliminated by changing the ocean "t" cell to a land value.
c
c     Isolated bays will be converted into land unless the directive
c     -Dleave_nonadvective is selected.
c
c     The directive -Dprompt_nonadvective enables interactive editing
c     of kmt at and in a neighborhood of a nonadvective "isolated bay"
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt)
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
#include "topog.h"
      dimension xu(imt), yu(jmt), zw(km)
c
c     statement function : depth at u/v points
c
      ku(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
c
      print '(/,a,/)'
     &,'Looking for "isolated bays" having no possible advective flow.'
#if !defined leave_nonadvective
      print '(a,/)', 'Default action is -Dfill_nonadvective'
#else
      print '(a,/)', 'Default action is -Dleave_nonadvective'
#endif
      write (io_del_kmt, '(t1,a6,t7,a)') 'c     ',
     &                                   'modifying nonadvective cells'
c
2000  n = 0
      do i=2,imt-1
        do jrow=jmt-1,2,-1
          if (kmt(i,jrow) .ne. 0) then
            if (ku(i  ,jrow) .eq. 0 .and. ku(i  ,jrow-1) .eq. 0 .and.
     &          ku(i-1,jrow) .eq. 0 .and. ku(i-1,jrow-1) .eq. 0) then
              n = n + 1
              kmt_bay = kmt(i,jrow)
              call clear_kmt_options ()
c
c             option 1: fill cell to land
c
              call kmt_option (1, i, jrow, kmt_bay, 0, kmt)
c
c             option 2: interactive user changes
c             option 3: leave nonadvective bay (requires no kmt changes)
c
              call select_option ('bay'
     &,                           i, jrow, kmt, xu, yu, zw)
#ifdef debug_topog
c
c             print disposition of nonadvective bay
c
              if (kmt(i,jrow) .eq. kmt_bay) then
                print '(a,i4,a,i4,a)'
     &,         'b)  left bay at (',i,',',jrow,')'
              else if (kmt(i,jrow) .eq. 0) then
                print '(a,i4,a,i4,a)'
     &,         'b) made land at (',i,',',jrow,')'
              else
                print '(a,i4,a,i4,a,i4)'
     &,         'b) changed kmt(',i,',',jrow,') to ', kmt(i,jrow)
              end if
#endif
            end if
          endif
        enddo
      enddo
c
      call enter_kmt_changes()
c
      if (n .ne. 0) then
        print '(a,i4,a)'
     &,'=>Warning: detected ',n,' "bays" that had no advective flow'
        print '(a,/)'
     &,'=>Rechecking that changes give no island PERIMETER VIOLATIONS'
        call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
#if !defined drive_topog && !defined permit_runtime_topog_changes
        if (n_del_kmt .ne. 0) then
          go to 2001
        endif
#endif
#if !defined leave_nonadvective && !defined prompt_nonadvective
c
c       repeat unless leaving nonadvective bays is permitted
c
        go to 2000
#endif
      else
        print '(a)', 'No nonadvective bays detected'
      end if
2001  continue
c
      call area_volume (kmt, xu, yu, zw)
#ifdef sgi_display
# ifdef prompt_nonadvective
      call clear_edit_window
# endif
      call display_kmt ('kmt/bays fixed', xu, yu, zw, kmt, imt,jmt,km)
c     call display_grid ('kmt/bays fixed with grid',xu,yu,zw,imt,jmt,km)
c     call waitmouse
#endif
c
c-----------------------------------------------------------------------
c     apply boundary conditions
c-----------------------------------------------------------------------
c
      call kmtbc (kmt, imt, jmt)
c
      return
      end



      subroutine potholes (kmt,map,iperm,jperm,iofs,nippts,nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)
c
c-----------------------------------------------------------------------
c     detect potholes:
c     These are situations where an ocean point on the "t" grid is
c     deeper than its 4 surrounding neighbor cells, which prevents
c     lateral and vertical advection of tracers at this depth.
c     These points usually cause no problems, but they also do not
c     contribute to ocean circulation.  They may be eliminated by
c     changing kmt at this "t" grid location to the maximum of kmt on
c     the 4 neighboring locations.
c
c     Potholes will be converted into land unless the directive
c     -Dleave_potholes is selected.
c
c     The directive -Dprompt_potholes enables interactive editing
c     of kmt at and in a neighborhood of "potholes"
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt)
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
#include "topog.h"
      dimension xu(imt), yu(jmt), zw(km)
c
c
      print '(//,a,a,/,a/)'
     &,'=> Searching for "potholes", i.e., '
     &,'"t" cells deeper than all 4 neighbor cells'
     &,'    (usually not a problem but worth knowing where they are)'
#if !defined leave_potholes
      print '(a,/)', 'Default action is -Dfill_potholes'
#else
      print '(a,/)', 'Default action is -Dleave_potholes'
#endif
      write (io_del_kmt, '(t1,a6,t7,a)') 'c     ',
     &                                   'modifying pothole cells'
      n = 0
      do i=2,imt-1
        do jrow=jmt-1,2,-1
          if (kmt(i+1,jrow) .lt. kmt(i,jrow) .and.
     &        kmt(i-1,jrow) .lt. kmt(i,jrow) .and.
     &        kmt(i,jrow+1) .lt. kmt(i,jrow) .and.
     &        kmt(i,jrow-1) .lt. kmt(i,jrow)) then
            n = n + 1
            call clear_kmt_options ()
c
c           option 1: fill this cell to the max of its neighbors
c           note: this still is not enough to guarantee flow
c
            kmt_pot = kmt(i,jrow)
            kmt_fill = max(kmt(i+1,jrow), kmt(i-1,jrow),
     &                     kmt(i,jrow+1), kmt(i,jrow-1))
            call kmt_option (1, i, jrow, kmt_pot, kmt_fill, kmt)
c
c           option 2: interactive user changes
c           option 3: leave pothole (requires no kmt changes)
c
            call select_option ('pothole'
     &,                         i, jrow, kmt, xu, yu, zw)
#ifdef debug_topog
c
c           print disposition of pothole
c
            if (kmt(i,jrow) .eq. kmt_pot) then
              print '(a,i4,a,i4,a)'
     &,         'c)  left pothole at (',i,',',jrow,')'
            else if (kmt(i,jrow) .eq. 0) then
              print '(a,i4,a,i4,a)'
     &,         'c) filled pothole at (',i,',',jrow,')'
            else
              print '(a,i4,a,i4,a,i4)'
     &,         'c) changed kmt(',i,',',jrow,') to ', kmt(i,jrow)
            end if
#endif
          endif
        enddo
      enddo
c
      call enter_kmt_changes()
c
      call area_volume (kmt, xu, yu, zw)
#ifdef sgi_display
# ifdef prompt_potholes
      call clear_edit_window
# endif
      call display_kmt ('kmt/potholes', xu, yu, zw, kmt, imt,jmt,km)
c     call display_grid ('kmt fixed with grid',xu,yu,zw,imt,jmt,km)
c     call waitmouse
#endif
c
c-----------------------------------------------------------------------
c     apply boundary conditions
c-----------------------------------------------------------------------
c
      call kmtbc (kmt, imt, jmt)
c
      return
      end



      subroutine min_depth (kmt, kmt_min, xu, yu, zw)
c
c     limit the minimum number of levels. kmt_min should be >= 2
c
#include "size.h"
#include "stdunits.h"
#include "topog.h"
c
      dimension kmt(imt,jmt)
      dimension xu(imt), yu(jmt), zw(km)
c
      write (stdout,'(/a,i3,a/)')
     &  'Making every ocean "kmt(i,jrow)" at least kmt_min = ',kmt_min
     &, ' levels deep.'
#if defined round_shallow ||  \
    !(defined deepen_shallow || defined fill_shallow)
      print '(a,/)', 'Default action is -Dround_shallow'
#endif
#if defined fill_shallow
      print '(a,/)', 'Default action is -Dfill_shallow'
#endif
#if defined deepen_shallow
      print '(a,/)', 'Default action is -Ddeepen_shallow'
#endif
      write (io_del_kmt, '(t1,a6,t7,a)') 'c     ',
     &                                   'modifying shallow cells'
      n = 0
      do i=2,imt-1
        do jrow=jmt-1,2,-1
          if (kmt(i,jrow) .ne. 0 .and. kmt(i,jrow) .lt. kmt_min) then
            n = n + 1
            call clear_kmt_options ()
            kmt_shallow = kmt(i,jrow)
            kmt_fill = 0
c
c           option 1: fill this cell
c
            call kmt_option (1, i, jrow, kmt_shallow, kmt_fill, kmt)
c
#if defined prompt_shallow
c           option 2: interactive user changes
#else
c           option 2: round to nearer of 0 and zw(kmt_min)
c
            if (zw(kmt_shallow) .lt. 0.5*zw(kmt_min)) then
              call kmt_option (2, i, jrow, kmt_shallow, kmt_fill, kmt)
            else
              call kmt_option (2, i, jrow, kmt_shallow,
     &                         max(kmt(i,jrow), kmt_min, 2), kmt)
            end if
#endif
c
c           option 3: deepen cell to kmt_min [or at least 2]
c
            call kmt_option (3, i, jrow, kmt_shallow,
     &                       max(kmt(i,jrow), kmt_min, 2), kmt)
c
            call select_option ('shallow', i, jrow, kmt, xu, yu, zw)
c
c           print disposition of shallow cell
c
#ifdef debug_topog
            write (stdout,'(a,i4,a,i4,a,i4,a,i4)')
     &      'a) changed shallow kmt(',i,',',jrow,') from ', kmt_shallow
     &,     ' to ', kmt(i,jrow)
#endif
          endif
        enddo
      enddo
      write (stdout,'(/a/)')
     & 'End modifications for minimum number of levels'
c
      call enter_kmt_changes()
c
      call area_volume (kmt, xu, yu, zw)
#ifdef sgi_display
# ifdef prompt_shallow
      call clear_edit_window
# endif
#endif
c
      return
      end


      subroutine announce_phase(title)
      character *(*) title
c      print *, 'tp: starting phase ', title
#ifdef sgi_display
      call gl_announce_phase(title)
#endif
      return
      end



#ifdef ascii_display
c
c=======================================================================
c ASCII interface for topography editor
c
c Written by Hans Vahlenkamp    e-mail: hnv@gfdl.gov
c UCAR/GFDL
c Princeton University
c Last updated: 6/9/95
c=======================================================================
c
      subroutine ascii_interactive(problem,idefault,i_opt,icen2,jcen2,
     &                             kmt,ied2,jed2,xu,yu,zw)

#include "size.h"
#include "topog.h"
      dimension kmt(imt, jmt)

      character *(*) problem
      character *(32) title
      dimension xu(imt), yu(jmt), zw(km)
      character choice
      save choice
      character*(32) locstr
      integer iedit, jedit, ieditwidth, jeditwidth
      character*(20) formatstr

c     copy input arguments so that changes made in this subroutine do
c     not cause changes in the input arguments
c
      icen = icen2
      jcen = jcen2
      ied  = ied2
      jed  = jed2
c
c     initial interactive edit state
c
  100 continue
      print '(/,a,a,a,2i4)',
     &   'Interactive editing for ',problem(1:len_trim(problem)),
     &   ' at', ied, jed
      if (problem .ne. 'perim' .and. choice .eq. 'x') then
        i_opt=3
        return
      end if
      call ascii_display_zoom(xu,yu,zw,kmt,icen,jcen,imt,jmt,km)

  101 continue

      if (problem .eq. 'shallow') then
        call prompt('(f)ill','(e)dit anywhere','(d)eepen','(i)gnore',
     &              ' ')
        read "(a)",choice
        if (choice .eq. 'f') then
          i_opt=1
          call accept_option(i_opt,kmt)
          goto 210
        else if (choice .eq. 'e') then
          goto 300
        else if (choice .eq. 'd') then
          i_opt=3
          call accept_option(i_opt,kmt)
          goto 230
        else if (choice .eq. 'i') then
          i_opt=2
          goto 999
        else
          goto 101
        endif
      else if (problem .eq. 'perim') then
        call prompt('(f)ill','(e)dit anywhere','(w)iden','(i)gnore',' ')
        read "(a)",choice
        if (choice .eq. 'f') then
          i_opt=1
          call accept_option(i_opt,kmt)
          goto 210
        else if (choice .eq. 'e') then
          goto 300
        else if (choice .eq. 'w') then
          i_opt=3
          call accept_option(i_opt,kmt)
          goto 230
        else if (choice .eq. 'i') then
          i_opt=2
          goto 999
        else
          goto 101
        endif
      else if (problem .eq. 'bay') then
        call prompt('(f)ill','(e)dit anywhere','(l)eave','(i)gnore',
     &              '(x)leave all & exit')
        read "(a)",choice
        if (choice .eq. 'f') then
          i_opt=1
          call accept_option(i_opt,kmt)
          goto 210
        else if (choice .eq. 'e') then
          goto 300
        else if (choice .eq. 'l') then
          i_opt=3
          call accept_option(i_opt,kmt)
          goto 230
        else if (choice .eq. 'i') then
          i_opt=2
          goto 999
        else if (choice .eq. 'x') then
          i_opt=3
          goto 100
        else
          goto 101
        endif
      else if (problem .eq. 'pothole') then
        call prompt('(f)ill','(e)dit anywhere','(l)eave','(i)gnore',
     &              '(x)leave all & exit')
        read "(a)",choice
        if (choice .eq. 'f') then
          i_opt=1
          call accept_option(i_opt,kmt)
          goto 210
        else if (choice .eq. 'e') then
          goto 300
        else if (choice .eq. 'l') then
          i_opt=3
          call accept_option(i_opt,kmt)
          goto 230
        else if (choice .eq. 'i') then
          i_opt=2
          goto 999
        else
          goto 101
        endif
      else if (problem .eq. 'edit') then
        call prompt(' ','(e)dit anywhere',' ','(q)uit editing',' ')
        read "(a)",choice
        if (choice .eq. 'e') then
          goto 300
        else if (choice .eq. 'q') then
          goto 999
        else
          goto 101
        endif
      else
        call prompt('default','free edit','alternative','ignore',' ')
      endif
      read "(a)",choice

c-----------------------------------------------------------------------

c Confirm chosen option

  210 continue

      call ascii_display_zoom (xu,yu,zw,kmt,icen,jcen,imt,jmt,km)

  211 continue

      call prompt('(a)ccept & next (or return)','(e)dit & accept',
     &            '(u)ndo',' ',' ')
      read "(a)",choice

c "accept & next" - accept option 1 and exit
      if ((choice .eq. 'a') .or. (choice .eq. ' ')) then
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' accepted'
#  endif
        goto 999
c "undo" - undo option 1
      else if (choice .eq. 'u') then
        call undo_option(i_opt, kmt)
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' undone'
#  endif
        goto 100
c "edit & accept" - accept option 1 and enter free editing
      else if (choice .eq. 'e') then
#  ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' accepted--enter free edit'
#  endif
        goto 300
c Invalid option chosen
      else
        goto 211
      endif

c-----------------------------------------------------------------------

c Check if option should be accepted

  230 continue

      call ascii_display_zoom(xu,yu,zw,kmt,icen,jcen,imt,jmt,km)

  231 continue

      call prompt('(u)ndo','(e)dit & accept',
     &            '(a)ccept & next (or return)',' ',' ')
      read "(a)",choice

c "accept & next" - accept option 3 and exit
      if ((choice .eq. 'a') .or. (choice .eq. ' ')) then
#ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' accepted'
#endif
        goto 999
c "undo" - undo option 3
      else if (choice .eq. 'u') then
        call undo_option(i_opt, kmt)
#ifdef debug_interactive
        print '(a,i1,a)','option ',i_opt,' undone'
#endif
        goto 100
c "edit & accept" - accept option 3 and enter free editing
      else if (choice .eq. 'e') then
        goto 300
c Invalid option chosen
      else
        goto 231
      endif

c-----------------------------------------------------------------------

c "free edit" mode, but edit cell not yet selected

  300 continue

      Write (*,'(a,a,$)')
     & 'Enter edit cell: i j [depth] ... '
     &,'or "q"(quit free edit) or "exit"(stop/save) > '
      read '(a32)',locstr

      i=1
      ieditwidth=0
  301 continue
      if (((locstr(i:i) .lt. '0') .or. (locstr(i:i) .gt. '9')) .and.
     &    (i .le. len_trim(locstr))) then
        i=i+1
        ieditwidth=ieditwidth+1
        goto 301
      endif
  302 continue
      if (((locstr(i:i) .ge. '0') .and. (locstr(i:i) .le. '9')) .and.
     &    (i .le. len_trim(locstr))) then
        i=i+1
        ieditwidth=ieditwidth+1
        goto 302
      endif

      jeditwidth=0
  303 continue
      if (((locstr(i:i) .lt. '0') .or. (locstr(i:i) .gt. '9')) .and.
     &    (i .le. len_trim(locstr))) then
        i=i+1
        jeditwidth=jeditwidth+1
        goto 303
      endif
  304 continue
      if (((locstr(i:i) .ge. '0') .and. (locstr(i:i) .le. '9')) .and.
     &    (i .le. len_trim(locstr))) then
        i=i+1
        jeditwidth=jeditwidth+1
        goto 304
      endif

      keditwidth=0
  305 continue
      if (((locstr(i:i) .lt. '0') .or. (locstr(i:i) .gt. '9')) .and.
     &    (i .le. len_trim(locstr))) then
        i=i+1
        keditwidth=keditwidth+1
        goto 305
      endif
  306 continue
      if (((locstr(i:i) .ge. '0') .or. (locstr(i:i) .le. '9')) .and.
     &    (i .le. len_trim(locstr))) then
        i=i+1
        keditwidth=keditwidth+1
        goto 306
      endif

      if (locstr(1:1) .eq. 'q') then
        i_opt=2
        goto 999
      elseif (locstr(1:4) .eq. 'exit') then
c
c       Immediate and final exit from topog!  delta.kmt file will have
c       all changes made so far, but topography is not guaranteed to be
c       free of violations.  After this exit, a new run of topog will
c       read all changes from the delta.kmt file and continue the
c       current run.  NO check will be made to insure that define
c       options have not changed between runs.
c
        call enter_kmt_changes()
c
        call area_volume (kmt, xu, yu, zw)
        print '(/,a,/,a,a,)', 'WARNING:  USER TERMINATION OF TOPOG',
     &              '          Topography not guaranteed to fulfill',
     &              ' model requirements'
        stop ('==>topog')
      else
        nformats=0
        if (ieditwidth .gt. 0) then
          nformats=nformats+1
        endif
        if (jeditwidth .gt. 0) then
          nformats=nformats+1
        endif
        if (keditwidth .gt. 0) then
          nformats=nformats+1
        endif
        if (nformats .lt. 3) goto 307
        write (formatstr,'(a,i3,a,i3,a,i3,a)')
     &         '(i',ieditwidth,',i',jeditwidth,',i',keditwidth,')'
        read (locstr,formatstr,err=307) iedit,jedit,kdepth
        if ((kdepth .lt. 0) .or. (kdepth .gt. km)) then
          write (*,'(a,i3,a)') 'Invalid cell depth (k = 0 to ',km,')'
          goto 300
        endif
        goto 308

  307 continue

        write (formatstr,'(a,i3,a,i3,a)')
     &         '(i',ieditwidth,',i',jeditwidth,')'
        read (locstr,formatstr,err=300)
     &        iedit,jedit
        kdepth=-1

  308 continue

        if ((iedit .lt. 1) .or. (iedit .gt. imt) .or.
     &      (jedit .lt. 1) .or. (jedit .gt. jmt)) then
          write (*,'(a,i3,a,i3,a)')
     &      'Invalid cell location (i = 1 to ',imt,', j = 1 to ',jmt,')'
          goto 300
        endif

        goto 400

      endif

c-----------------------------------------------------------------------

  400 continue

c "free edit" mode

c clear option 2

      i_opt=2
      do j_opt=1,len_opt
        do k_opt=1,4
          kmt_opt(i_opt,j_opt,k_opt)=0
        enddo
      enddo

      kmt_orig=kmt(iedit,jedit)

#ifdef debug_interactive
      print *, 'tp: user edit at (', iedit, jedit,')'
#endif

      if (kdepth .ne. -1) then
        kmt(iedit,jedit)=kdepth
      endif

      call ascii_display_zoom(xu,yu,zw,kmt,iedit,jedit,imt,jmt,km)

  401 continue

      call prompt('(s)hallower','(n)ew edit cell','(d)eeper',
     &            '(q)uit editing',' ')
      read "(a)",choice

c "shallower" - decrease kmt(iedit,jedit) by 1
      if (choice .eq. 's') then
        kmt(iedit,jedit)=max(kmt(iedit,jedit)-1,0)
        print '(a,2i4,a,i4)',
     &        'd) decreased kmt(',iedit,jedit,') to ',
     &        kmt(iedit,jedit)
        call ascii_display_zoom(xu,yu,zw,kmt,iedit,jedit,imt,jmt,km)
        goto 401
c "deeper" - increase kmt(iedit,jedit) by 1
      else if (choice .eq. 'd') then
        kmt(iedit,jedit)=min(kmt(iedit,jedit)+1,km)
        print '(a,2i4,a,i4)',
     &        'd) increased kmt(',iedit,jedit,') to ',
     &        kmt(iedit,jedit)
        call ascii_display_zoom(xu,yu,zw,kmt,iedit,jedit,imt,jmt,km)
        goto 401
c "new edit cell" - accept changes to edit cell
      else if (choice .eq. 'n') then
        call ascii_display_zoom(xu,yu,zw,kmt,iedit,jedit,imt,jmt,km)
        if (kmt(iedit,jedit) .ne. kmt_orig) then
          call kmt_option(i_opt,iedit,jedit,kmt_orig,kmt(iedit,jedit),
     &                    kmt)
        endif
        goto 300
c "quit editing" - accept changes to current edit cell and end free edit
      else if (choice .eq. 'q') then
        if (kmt(iedit,jedit) .ne. kmt_orig) then
          call kmt_option(i_opt,iedit,jedit,kmt_orig,kmt(iedit,jedit)
     &,                   kmt)
        endif
        goto 999
c Invalid option chosen
      else
        goto 401
      endif

c-----------------------------------------------------------------------

  999 continue
      call accept_option(i_opt, kmt)
#  ifdef debug_interactive
      print '(a,i1,a)','option ',i_opt,
     &      ' accepted. leaving interactive edit'
#  endif

c-----------------------------------------------------------------------

      return
      end


      subroutine ascii_prompt(opt1, opt2, opt3, quit_opt, exit_opt)

      character *(*) opt1, opt2, opt3, quit_opt, exit_opt

      write (*,'(a,$)') 'Options: '
      if (opt1 .ne. ' ') then
        write (*,'(a,$)') opt1
      endif
      if (opt1 .eq. ' ') then
        write (*,'(a,$)') opt2
      else
        write (*,'(a,a,$)') ', ',opt2
      endif
      if (opt3 .ne. ' ') then
        write (*,'(a,a,$)') ', ',opt3
      endif
      if (quit_opt .ne. ' ') then
        write (*,'(a,a,$)') ', ',quit_opt
      endif
      if (exit_opt .ne. ' ') then
        write (*,'(a,a,$)') ', ',exit_opt
      endif
      write (*,'(a,$)') ' > '

      return
      end


c This common block contains a flag for printing the window size once

      block data iffirstsub
      common /firstflag/ iffirst
      data iffirst /1/
      end


      subroutine ascii_display_zoom(xu1,yu1,zw1,kmt,imid,jmid,imt2,
     &                              jmt2,km2)

#include "size.h"

      dimension xu1(imt), yu1(jmt)
      dimension zw1(km), kmt(imt,jmt)
      integer imid, jmid, imt2, jmt2, km2
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
      integer istart, iend, jstart, jend
      character*10 szwkmt
      integer ndashes
      common /firstflag/ iffirst

c Horizontal and vertical number of cells in the edit/zoom window

      parameter (iwidth=20, jheight=16)

c Verify/fix imt, jmt, km from "size.h" via "ijk_max.h"

      call size_check(imt2,jmt2,km2,'ascii_display_zoom','fix_size')

c Get ranges of grid points to find visible cells

      i1=imid-iwidth/2
      i2=imid+iwidth/2
      j1=jmid-jheight/2
      j2=jmid+jheight/2
      ndashes=(iwidth+1)*3+1
      istart=0
      iend=0
      jstart=0
      jend=0
      if (i1 .lt. 1) then
        istart=abs(i1-1)*3
        i1=1
      endif
      if (i2 .gt. imt2) then
        iend=abs(i2-imt2)*3
        i2=imt2
      endif
      if (j1 .lt. 1) then
        jstart=abs(j1)
        j1=1
      endif
      if (j2 .gt. jmt2) then
        jend=abs(j2-jmt2)
        j2=jmt2
      endif

c Print message once indicating location of width and height parameters

      if (iffirst .eq. 1) then
        write (*,'(/,a,/,a)')
     &     'The width and height of the window are parameters in the ',
     &     'ascii_display_zoom subroutine in the topog.F file.'
        write (*,'(a,i3,a,i3,a,/)')
     &         '(iwidth = ',iwidth,', jheight = ',jheight,')'
        iffirst=0
      endif

c Print longitude, latitude, and z value of selected cell

      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
      if (kmt(imid,jmid) .gt. 0) then
        write (szwkmt,'(f7.1)') zw(kmt(imid,jmid))/100
      else
        szwkmt='0'
      endif
      write (*,'(a,f5.1,a,f5.1,a,a,a)')
     &       '(lon = ', 0.5*(xu(imid)+xu(imid-1)),
     &       ', lat = ',0.5*(yu(jmid)+yu(jmid-1)),
     &       ')   z =',szwkmt(1:len_trim(szwkmt)),' m'

c Begin printing zoom window

      write (*,'(a,$)') ' '
      if (istart .ne. 0) then
        do i=1,istart
          write (*,'(a,$)') ' '
        enddo
      endif
      do i=i1,i2
        write (*,'(1x,i2.2,$)') i
      enddo
      write (*,'(/,a,$)') '+'
      do i=1,ndashes
        write (*,'(a,$)') '-'
      enddo
      write (*,'(a)') '+'

c Insert any leading blank lines if necessary

      if (jend .ne. 0) then
        do j=1,jend
          write (*,'(a,$)') '|'
          do i=1,ndashes
            write (*,'(a,$)') ' '
          enddo
          write (*,'(a)') '|'
        enddo
      endif

      do j=j2,j1,-1
      
c Insert any leading spaces on a line if necessary

        write (*,'(a,$)') '|'
        if (istart .ne. 0) then
          do i=1,istart
            write (*,'(a,$)') ' '
          enddo
        endif

c Print actual kmt values

        do i=i1,i2
          if ((i .eq. imid) .and. (j .eq. jmid)) then
            write (*,'(a1,$)') '>'
            if (kmt(i,j) .eq. 0) then
              write (*,'(a2,$)') '..'
            else
              write (*,'(i2.2,$)') kmt(i,j)
            endif
            write (*,'(a1,$)') '<'
            prevedit=1
            if (i .eq. i2) then
              lastedit=1
            endif
          else
            if (prevedit .eq. 1) then
              prevedit=0
              if (kmt(i,j) .eq. 0) then
                write (*,'(a2,$)') '..'
              else
                write (*,'(i2.2,$)') kmt(i,j)
              endif
            else
              if (kmt(i,j) .eq. 0) then
                write (*,'(1x,a2,$)') '..'
              else
                write (*,'(1x,i2.2,$)') kmt(i,j)
              endif
            endif
          endif
        enddo

c Insert any trailing spaces on a line if necessary

        if (iend .ne. 0) then
          do i=1,iend
            write (*,'(a,$)') ' '
          enddo
        endif
        if (lastedit .eq. 1) then
          write (*,'(a,i2.2)') '| ',j
          lastedit=0
        else
          write (*,'(a,i2.2)') ' | ',j
        endif

      enddo

c Insert any trailing blank lines if necessary

      if (jstart .ne. 0) then
        do j=1,jstart
          write (*,'(a,$)') '|'
          do i=1,ndashes
            write (*,'(a,$)') ' '
          enddo
          write (*,'(a)') '|'
        enddo
      endif
      write (*,'(a,$)') '+'
      do i=1,ndashes
        write (*,'(a,$)') '-'
      enddo
      write (*,'(a,/)') '+'

      return
      end


      subroutine make_xyz0 (xu1, yu1, zw1, xu0, yu0, zw0)
c
c-----------------------------------------------------------------------
c     extend longitude, latitude, and depth arrays to subscript 0
c-----------------------------------------------------------------------
c
#include "size.h"
      dimension xu1(1:imt), yu1(1:jmt), zw1(1:km)
      dimension xu0(0:imt), yu0(0:jmt), zw0(0:km)
c      
      call coord0(xu1, xu0, imt)
      call coord0(yu1, yu0, jmt)
      call coord0(zw1, zw0, km)
c     just to be on the safe side:
      zw0(0) = 0
      return
      end


      subroutine coord0 (xu1, xu0, imt)
c
c-----------------------------------------------------------------------
c     extend one coordinate axes to subscript 0
c-----------------------------------------------------------------------
c
      dimension xu1(1:imt), xu0(0:imt)
      do i=1,imt
        xu0(i) = xu1(i)
      end do
      xu0(0) = 2 * xu1(1) - xu1(2)
      return
      end

c-----------------------------------------------------------------------
c
c=======================================================================
c End of ASCII interface for topography editor
c=======================================================================
c
#endif
