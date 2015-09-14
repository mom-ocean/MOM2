#ifdef drive_grids
      program driver
c
c=======================================================================
c
c              G R I D   G E N E R A T O R   M O D U L E
c
c     To experiment with designing grids for MOM: 
c     
c        1) inputs: specify the grid in the USER INPUT section below
c                          
c        2) compile and run this module by choosing the appropriate 
c           options in the "run_grids" script and executing the script.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      parameter (maxlen=1500)
      dimension dxtdeg(maxlen), dytdeg(maxlen), dzt(maxlen)
      dimension dxudeg(maxlen), dyudeg(maxlen), dzw(0:maxlen)
      dimension xt(maxlen), xu(maxlen), yt(maxlen), yu(maxlen)
      dimension zt(maxlen), zw(maxlen)
c
#if defined read_my_grid
      call ioinit
#endif
c
c-----------------------------------------------------------------------
c     calculate resolution and coordinates for "t", "u", "w" grid cells
c-----------------------------------------------------------------------
c
      write (stdout,'(//,39x,a,/)') 'Driving the grid module:'
c
      call gcoord (maxlen, imt, jmt, km, dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
      if (imt .gt. maxlen .or. jmt .gt. maxlen .or. km .gt. maxlen) then
         write (stdout,*) ' => increase maxlen in the grids driver'
	 stop '=>driver'
      endif
c
c-----------------------------------------------------------------------
c     prompt user to change imt, jmt, and km in "size.h"
c-----------------------------------------------------------------------
c
      write (stdout,'(//1x,a/)') '==>GRID INSTALLATION DIRECTIONS'
      write (stdout,'(6x,a,/,6x,a,i4,a,i4,a,i3,a,/)')
     & 'To use this grid in MOM, change parameter in "size.h" to:'
     &,'parameter (imt=',imt,', jmt=',jmt,', km=',km,')'
c
      write (stdout,'(/a/)')
     &' (To change resolution, read the USER INPUT section in grids.F)'
      stop
      end

# define driver_only
# include "util.F"
# include "iomngr.F"

#endif



      subroutine gcoord (maxlen, imt, jmt, km, dxtdeg, dytdeg, dxudeg
     &,                  dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c=======================================================================
c
c                      G R I D    C O O R D I N A T E S
c
c     Construct grid point coordinates and resolution
c
c     input:
c
c     maxlen = maximum number of grid cells in latitude, longitude,
c              and depth
c
c              set grid specifications in USER INPUT section.
c
c     output:
c
c     imt    = number of longitudes
c     jmt    = number of latitudes
c     km     = number of depths
c     dxtdeg = width of "t" grid cells (degrees)
c     dytdeg = height of "t" grid cells (degrees)
c     dxudeg = width of "u" grid cells (degrees)
c     dyudeg = height of "u" grid cells (degrees)
c     dzt    = thickness of "t" grid cells (cm)
c     dzw    = thickness of "w" grid cells (cm)
c     xt     = longitude at centers of "t" grid cells (degrees)
c     xu     = longitude at centers of "u" grid cells (degrees)
c     yt     = latitude at centers of "t" grid cells (degrees)
c     yu     = latitude at centers of "u" grid cells (degrees)
c     zt     = depth at centers of "t" grid cells (centimeters)
c     zw     = depth at centers of "u" grid cells (centimeters)
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      parameter (maxbounds=11)
      dimension xt(maxlen), yt(maxlen), xu(maxlen), yu(maxlen)
      dimension zw(maxlen), zt(maxlen)
c
      dimension dxtdeg(maxlen), dytdeg(maxlen), dzt(maxlen)
      dimension dxudeg(maxlen), dyudeg(maxlen), dzw(0:maxlen)
c
      dimension x_lon(maxbounds),   dx_lon(maxbounds)
      dimension y_lat(maxbounds),   dy_lat(maxbounds)
      dimension z_depth(maxbounds), dz_depth(maxbounds)
c
c-----------------------------------------------------------------------
c     USER INPUT ===>  read this section and specify the grid.
c-----------------------------------------------------------------------
c
c     How to define a grid resolution:
c
c     The ocean model domain is composed of one or more regions in
c     latitude, longitude, and depth. Latitude and longitude are in
c     degrees east of Greenwich and depth is in centimeters. Each
c     region is defined by its bounds and resolution at those bounds.
c     Within each region, resolution may be constant or smoothly
c     varying but there must be an integral number of grid cells
c     contained within the region`s bounds.
c     
c     If the resolution at both bounds of a region is the same, then
c     resolution within the region is constant. If the bounding
c     resolutions differ, resolution varies smoothly across the region
c     to make the transition. The functional for the variation is
c     taken to be a cosine. Regions sharing a common boundary have the
c     same resolution at that boundary AND the first derivitive is
c     zero there to minimize effects of numerics. In the vertical, the
c     last region allows a stretching factor to provide a more drastic
c     fall off of resolution to the bottom if desired.
c 
c     Example:
c
c     Define a domain composed of two regions of longitude: The first
c     is bounded by longitudes 120E and 180E with constant 1.0 degree
c     resolution. The second extends east of 180E to 300E where the
c     resolution is desired to be 3 degrees (at 300E). To do this,
c     specify the following:
c     
c      parameter (nxlons=3)
c      data (x_lon(i), i=1,nxlons)  /120.0, 180.0, 300.0/
c      data (dx_lon(i),i=1,nxlons)  /  1.0,   1.0,   3.0/                                           
c
c      where...
c
c     "nxlons"  = number of bounding longitudes to define regions {3}
c     "x_lon"   = bounding longitudes {120.0E, 180.0E, 300.0E}
c     "dx_lon"  = resolution centered at "x_lon" {1.0, 1.0, 3.0}
c     
c     Note: Region bounds are not defined exactly as in MOM 1. Now they
c           are defined on "u" grid points (which are at the edges of
c           "T" cells) as described below.
c
c     The western edge of the 1st "T" cell will be at x_lon(1) and
c     the eastern edge of the last "T" cell within the first region
c     will be at x_lon(2). Since the resolution at both bounds is
c     1.0 degree, there will be a constant 1 degree resolution within
c     this first region.
c
c     The western edge of the 1st "T" cell within the second region is
c     at x_lon(2) and the eastern edge of the last "T" cell is at
c     x_lon(3). Since resolutions at the bounds of region 2 differ,
c     resolution will gradually taper from 1.0 degree at the western
c     bound to 3.0 degrees at the eastern bound. Note that the
c     resolution of the 1st "T" cell adjacent to x_lon(2) will
c     be slightly greater than 1.0 degree (it would be 1.0 degree if
c     the grid cell were centered exactly at x_lon(2) but it is not).
c
c     The number of "T" cells within each region will be equal to the
c     width of the region divided by the average resolution within the
c     region. The average resolution within the region is simply the
c     average of the bounding resolutions. In region #1, it will be
c     (1.0 + 1.0)/2 = 1.0 degree and in region #2 it will be
c     (1.0 + 3.0)/2 = 2.0 degree. Since region #1 is 60 degree wide,
c     the number of "T" cells will be 60.0/1.0 = 60 and for region #2
c     there will be 120.0/2.0 = 60 cells also.
c
c     How are "u" cells defined?     
c
c     Since MOM uses an Arakawa "B" grid, the "u" cells are offset from
c     the "T" cells. The "u" cell grid points are at the vertices of the
c     "T" cells and visa versa. The "u" cell with coordinates (i,j) is
c     on the northeast vertex of the "T" cell with coordinates (i,j). In
c     the vertical, "u" cells and "T" cells are at the same level.
c     However, advective vertical velocity is defined at the base of
c     "T" and "u" cells. Thus "w" cells can be thought of as vertical
c     advective velocity cells with the grid being staggared in the
c     vertical. The relationship between "T" (or "u") cells and "w"
c     cells is that the "w" cell with index (k) is at the base of the
c     "T" (or "u") cell with index k. 
c
c     What happens when the grid cells are stretched?
c
c     Grid stretching (variable resolution) is not done exactly as in
c     MOM 1. In MOM 1, the "u" cell resolution was the average of the
c     resolutions of the adjacent "T" cells. For instance: "dxt" was 
c     specified first... then dxu(i) = 0.5*(dxt(i+1) + dxt(i)).
c     In a stretched grid, the "T" grid points were always in the middle
c     of "T" cells and the "u" grid points were off center in "u" cells.
c
c     In MOM 2, it`s the other way around: the "T" cell resolution is
c     the average of the resolutions of the adjacent "u" cells. So,
c     "dxu" is specified first... then dxt(i) = 0.5*(dxu(i-1) + dxu(i))
c     Also, "u" grid points are always in the middle of "u" cells and
c     "T" grid points are off center in "T" cells.
c          
c     This allows for more accurate advection velocities to be defined
c     on the faces of the cells when the resolution varies. The
c     implication is that advection will be more accurate for tracers
c     in the horizontal and both tracers and momentum in the vertical.
c     Advection will not be more accurate for momentum in the horizontal
c     because "u" grid  points will always be in the center of the "u"
c     cells.
c
c     What if I want the resolution from MOM 1?
c
c     It is possible to closely approximate the horizontal resolution
c     from MOM 1 by enabeling the "-Dcentered_t" option. Whether this is
c     done or not, the "u" grid points are always on the edges of "T"
c     cells and the "T" grid points are always on the edges of "u" 
c     cells. In the vertical, MOM 1 used a gaussian for stretching the
c     resolution. MOM 2 uses the same cosine functional for varying
c     resolution in the vertical as in the horizontal. Additionally,
c     a stretching factor is used to give the functionality of enhanced
c     stretching towards the bottom of the ocean domain if desired.
c           
c     Another alternative is to read in whatever grid discretization you 
c     like by enabeling the -D"read_my_grid" option and being aware of
c     the potential consequences. The crucial point is that if you have
c     enough resolution to adequately resolve the scales of motion (but
c     who does?), it doesn`t much matter how the resolution is
c     distributed. However, keep in mind the above remarks about
c     advection.
c            
c-----------------------------------------------------------------------
c     The following specifies the global test grid resolution: 
c-----------------------------------------------------------------------
c
c     constant 4.0 degree resolution in longitude
c
      parameter (nxlons=2)
      data (x_lon(i), i=1,nxlons) / 0.0, 360.0/
      data (dx_lon(i),i=1,nxlons) / 4.0,   4.0/        
c
c     constant 3.0 degree resolution in latitude
c
      parameter (nylats=2)
      data (y_lat(j), j=1,nylats) / -87.0,  90.0/
      data (dy_lat(j),j=1,nylats) /  3.0,    3.0/        
c
c     two vertical regions: constant 25m in the upper 100m, and variable
c     resolution in the lower region. Note that by setting the stretch
c     factor "stretch_z" > 1.0, the last region in the vertical can
c     have the grid cells stretched further. To see the effect of
c     stretching on vertical cell thickness and number of vertical
c     cells, try using "stretch_z"=1.1 and gradually increase it.  
c
      parameter (nzdepths=3)
      data (z_depth(k), k=1,nzdepths) /  0.0e2,   100.0e2,  5600.0e2/
      data (dz_depth(k),k=1,nzdepths) / 25.0e2,    25.0e2,   975.0e2/        
      data stretch_z /1.0/
c
c     Some other examples:
c
c      4.0 degrees at the poles down to 2.0 deg at the equator
c
c      parameter (nylats=3)
c      data (y_lat(j), j=1,nylats) / -90.0,  0.0,  90.0/
c      data (dy_lat(j),j=1,nylats) /   4.0,  2.0,   4.0/        
c
c      set 1 deg at -30 deg lat down to 1/3 deg at -10 deg lat
c      constant 1/3 deg resolution from -10 deg to +10 deg lat and
c      1/3 deg at +10 deg lat to 1.0 deg at +30 deg lat 
c
c      parameter (nylats=4)
c      data (y_lat(j), j=1,nylats) / -30.0,  -10.0,    10.0,  30.0/
c      data (dy_lat(j),j=1,nylats) /   1.0,  0.3333,  0.3333,  1.0/        
c
c-----------------------------------------------------------------------
c     ==> end of USER INPUT
c-----------------------------------------------------------------------
c
c     set some constants
c
      p5 = 0.5
c
      write (stdout,'(/,a,/(/a))')
     &  'The following grid options have been selected:'
#ifdef generate_a_grid
     &,'-Dgenerate_a_grid'
#endif
#ifdef read_my_grid
     &,'-Dread_my_grid'
#endif
#ifdef centered_t
     &,'-Dcentered_t ("T" points centered in "T" cells as in MOM 1)'
#endif
c
      ncase = 0
#ifdef read_my_grid
c
      ncase = ncase + 1
      write (stdout,'(//a/a/)')
     & ' =>Reading the grid definition from file "grid.dta" '
     &,'   (this assumes it was created elsewhere and is available)'
      call getunit (io, 'grid.dta', 'u s r')
      read (io) imt, jmt, km
      read (io) (dxtdeg(i),i=1,imt)
     &,         (dytdeg(j),j=1,jmt)
     &,         (dxudeg(i),i=1,imt)
     &,         (dyudeg(j),j=1,jmt)
     &,         (dzt(k),k=1,km)
     &,         (dzw(k),k=0,km)
     &,         (xt(i),i=1,imt)
     &,         (xu(i),i=1,imt)
     &,         (yt(j),j=1,jmt)
     &,         (yu(j),j=1,jmt)
     &,         (zt(k),k=1,km)
     &,         (zw(k),k=1,km)
      call relunit (io)
# ifdef centered_t
      write (stdout,*)
     &'=>Error: option -Dcentered_t not compatible with -Dread_my_grid'
     &,'         ...remove one of them.'
# endif
#endif
#ifdef generate_a_grid
      ncase = ncase + 1
c
      write (stdout,'(//,36x,a,/)') 'G R I D   G E N E R A T I O N'
      write (stdout,'(/,43x,a,/)') 'Grid resolution:'
c
      if (nxlons .gt. maxbounds .or. nylats .gt. maxbounds .or.
     &    nzdepths .gt. maxbounds) then
        write (stdout,'(/a/)')
     &  ' Increase parameter "maxbounds" to contain all regions.'
        stop
      endif
c
c-----------------------------------------------------------------------
c     Calculate resolution in longitude. Add one boundary cell at the
c     start and end of the domain so that calculations are meaningful
c     for grid cells i=2,imt-1
c-----------------------------------------------------------------------
c
      nbpts = 2
      write (stdout,'(/a,i1,a)')
     & ' Generating the longitudinal resolution ( ', nbpts
     &, ' extra boundary points will be used).'
c
c     if a region contains the greenwich meridian, compensate by adding
c     360.0 degrees to all remaining regions
c
      do n=2,nxlons
        if (x_lon(n-1) .gt. x_lon(n)) then
	  x_lon(n) = x_lon(n) + 360.0
	  write (stdout,'(/,a,i1,a,a/)')
     &    ' Warning: adding 360.0 degrees to x_lon(',n,')'
     &,   ' to insure region boundaries increase monotonically'
	endif
      enddo
c
c     if width of domain exceeds 360.0 deg... limit width to 360.0
c
      if (x_lon(nxlons) - x_lon(1) .gt. 360.0) then
        write (stdout,'(/a,a,g14.7/)')
     &  '=>Warning: Domain width exceeds 360 deg. Restricting last'
     &, ' x_lon to ',x_lon(1)+360.0
          do m=1,nxlons
            write (stdout,'(i3,f10.5)') m, x_lon(m)
          end do
      endif
      if (x_lon(nxlons) - x_lon(1) .eq. 360.0) then
        if (dx_lon(nxlons) .ne. dx_lon(1)) then
          write (stdout,'(/a,a)')
     &    '=>Error: dx_lon(1) must equal dx_lon(last) when domain'
     &,   ' width = 360.0 degrees'
        endif
      endif

      stretch_x = 1.0
      call gcell (maxlen, nxlons, x_lon,  dx_lon, nbpts
     &,           imt, dxtdeg, dxudeg, stretch_x)
c
c-----------------------------------------------------------------------
c     Build the longitudinal grid points for a "B" grid
c     (account for an extra boundary point at the start)
c-----------------------------------------------------------------------
c
      xt(1) = x_lon(1) - p5*dx_lon(1)
      xu(1) = x_lon(1)
c
      do i=2,imt
        xu(i) = xu(i-1) +  dxtdeg(i)
        xt(i) = xt(i-1) +  dxudeg(i-1)
      enddo
c
c-----------------------------------------------------------------------
c     Calculate resolution in latitude. Add one boundary cell at the
c     start and end of the domain so that calculations are meaningful
c     for grid cells jrow=2,jmt-1
c-----------------------------------------------------------------------
c
      nbpts = 2
      write (stdout,'(/a,i1,a)')
     & ' Generating the latitudinal resolution ( ', nbpts
     &, ' extra boundary points will be used).'
c
c     if width of domain exceeds 180.0 deg... limit width to 180.0
c
      if (y_lat(nylats) - y_lat(1) .gt. 180.0) then
        write (stdout,'(/a,a,g14.7/)')
     &  '=>Warning: Latitudinal domain width exceeds 180 deg.'
     &, ' Restricting last y_lat to ',y_lat(1)+180.0
      endif
      do n=2,nylats
        if (y_lat(n-1) .gt. y_lat(n)) then
	  write (stdout,'(/,a,/a/)')
     &    ' =>Error: latitude boundaries must increase monotonically'
     &,   '          check the specifications in the USER INPUT section'
          do m=1,nylats
            write (stdout,'(i3,f10.5)') m, y_lat(m)
          end do
          stop '==>grids'
	endif
      enddo
      stretch_y = 1.0
      call gcell (maxlen, nylats, y_lat, dy_lat, nbpts
     &,           jmt, dytdeg, dyudeg, stretch_y)
c
c-----------------------------------------------------------------------
c     Build the latitudinal grid points on a "B" grid
c     (account for an extra boundary point at the start)
c-----------------------------------------------------------------------
c
      yt(1) = y_lat(1) - p5*dy_lat(1)
      yu(1) = y_lat(1)
      do jrow=2,jmt
        yu(jrow) = yu(jrow-1) +  dytdeg(jrow)
        yt(jrow) = yt(jrow-1) +  dyudeg(jrow-1)
      enddo
c
c-----------------------------------------------------------------------
c     Calculate resolution in depth. No boundary cells added here so
c     calculations are meaningful for k=1,km. Allow the bottom region 
c     to be stretched further if desired.
c-----------------------------------------------------------------------
c
      nbpts = 0
      write (stdout,'(/a,i1,a)')
     & ' Generating the vertical resolution ( ', nbpts
     &, ' extra boundary points will be used).'
      do n=2,nzdepths
        if (z_depth(n-1) .gt. z_depth(n)) then
	  write (stdout,'(/,a,/a/)')
     &    ' =>Error: depth boundaries must increase monotonically'
     &,   '          check the specifications in the USER INPUT section'
          stop '==>grids'
	endif
      enddo
      call gcell (maxlen, nzdepths, z_depth, dz_depth, nbpts
     &,           km, dzt, dzw(1), stretch_z)
c
c-----------------------------------------------------------------------
c     Build the vertical grid points on a "B" grid. The "T" and "u"
c     cells are staggared in the horizontal but at the same level in
c     the vertical. However, the "w" cells here refer to the vertical
c     advection velocities at the bottoms of the "u" and "T" cells.
c     (no extra boundary point at the start)
c-----------------------------------------------------------------------
c
      zt(1) = z_depth(1) + p5*dz_depth(1)
      zw(1) = z_depth(1) + dzt(1)
      do k=2,km
        zw(k) = zw(k-1) + dzt(k)
        zt(k) = zt(k-1) + dzw(k-1)
      enddo
c
c     set "w" cell thickness at surface and bottom to ocean part of cell
c
      dzw(0)  = zt(1)
      dzw(km) = zw(km) - zt(km)
c
c-----------------------------------------------------------------------
c     Print grid "t" cell resolution in longitude, latitude and depth
c     (also for "w" cells in depth)
c-----------------------------------------------------------------------
c
      write (stdout,9101) 
     &'Vertical resolution of cells "dzw(k)" k=0,', km, 'cm'
      write (stdout,9002) (dzw(k),k=0,km)
      write (stdout,9101) 
     &'Vertical resolution of cells "dzt(k)" k=1,', km, 'cm'
      write (stdout,9002) (dzt(k),k=1,km)
      write (stdout,9101) 
     &'Longitudinal resolution of cells "dxtdeg(i)" i=1,', imt, 'deg'
      write (stdout,9001) (dxtdeg(i),i=1,imt)
      write (stdout,9101) 
     &'Latitudinal resolution of cells "dytdeg(j)" jrow=1,', jmt, 'deg'
      write (stdout,9001) (dytdeg(jrow),jrow=1,jmt)
#endif
      if (ncase .eq. 0) then
        write (stdout,'(/a/a/)')
     & '=>Error: One of the following options must be enabled:'
     &,'         generate_a_grid    read_my_grid'
        stop '=>grid' 
      elseif (ncase .gt. 1) then
        write (stdout,'(/a/a/)')
     & '=>Error: Only one of the following options may be enabled:'
     &,'         generate_a_grid    read_my_grid'
        stop '=>grid' 
      endif
c
c-----------------------------------------------------------------------
c     Check if the "t" grid resolution is an average of the
c     "u" cell resolution. This insures more accurate advection of
c     tracers within a stretched grid.
c-----------------------------------------------------------------------
c
      num  = 0
      tolr = 1.e-5
      write (stdout,'(/)')
      do i=2,imt-1
        dxubar = p5*(dxudeg(i) + dxudeg(i-1))
        if (abs(dxubar-dxtdeg(i)) .gt. tolr) then
	  num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta x at i=',i
     &,   ' is not an average of adjacent "u" cell delta x`s'     
        endif     
      enddo
c
      do jrow=2,jmt-1
        dyubar = p5*(dyudeg(jrow) + dyudeg(jrow-1))
        if (abs(dyubar-dytdeg(jrow)) .gt. tolr) then
	  num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta y at jrow=',jrow
     &,   ' is not an average of adjacent "u" cell delta y`s' 
        endif     
      enddo
c
      tolr = 1.e0
      do k=2,km-1
        dzwbar = p5*(dzw(k) + dzw(k-1))
        if (abs(dzwbar-dzt(k)) .gt. tolr) then
	  num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta z at k=',k
     &,   ' is not an average of adjacent "w" cell delta z`s'     
        endif     
      enddo
c
      if (num .ne. 0) then
        write (stdout,'(/a/a/a/a//a,a/)')
     &  '==>Warning, At the above locations, advection of tracers is'
     &, 'not as accurate as it could be. If you are reading in your own'
     &, 'grid or constructing a grid as in MOM 1, we assume you want to'
     &, 'define the grid this way and we let you proceed from here...'
     &, 'Please read ALL the information in the USER INPUT section to '
     &, 'understand what this means' 
      endif
c
c-----------------------------------------------------------------------
c     Print all grid coordinates
c-----------------------------------------------------------------------
c
      write (stdout
     &,'(//,40x,a,//,a,g14.7,a,/a/,a,g14.7,a/a,/,a,g14.7,a)') 
     &  ' Grid Point Coordinate details: '
     &, ' The western edge of the 2nd "t" cell is at longitude:'
     &,   x_lon(1),' (deg)',' (the 1st "t" cell is a boundary cell)'
     &, ' The southern edge of the 2nd "t" cell is at latitude:'
     &,   y_lat(1),' (deg)',' (the 1st "t" cell is a boundary cell)'
     &,' The top edge of the 1st "t" cell is at z =',z_depth(1),' (cm)'
      write (stdout,'(/,a,g14.7,a/a/,a,g14.7,a/a/,a,g14.7,a/)') 
     &  ' The western edge of the 1st "u" cell is at longitude:', xt(1)
     &, ' (deg)',' (the 1st "u" cell is a boundary point)'
     &, ' The southern edge of the 1st "u" cell is at latitude:', yt(1)
     &, ' (deg)',' (the 1st "u" cell is a boundary point)'
     &, ' The top edge of the 1st "w" cell is at z =',zt(1),' (cm)'
      write (stdout,9103) km
      write (stdout,9002) (zt(k),k=1,km)
      write (stdout,9104) km
      write (stdout,9002) (zw(k),k=1,km)
      write (stdout,9105) jmt
      write (stdout,9001) (yt(jrow),jrow=1,jmt)
      write (stdout,9106) jmt
      write (stdout,9001) (yu(jrow),jrow=1,jmt)
      write (stdout,9107) imt
      write (stdout,9001) (xt(i),i=1,imt)
      write (stdout,9108) imt
      write (stdout,9001) (xu(i),i=1,imt)
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     insure that yu(jmt-1) = 0.0 (equator) when using symmetry     
c-----------------------------------------------------------------------
c
      if (yu(jmt-1) .ne. 0.0) then
        write (stdout,*) '=> Error: yu(jmt-1) must = 0.0 for symmetry'
	stop '=>gcoord'
      endif
#endif
c
c---------------------------------------------------------------------
c     compute a grid checksum
c---------------------------------------------------------------------
c
      cksum = 0.0
      cksum = cksum + checksum (xt, imt, 1)
      cksum = cksum + checksum (yt, jmt, 1)
      cksum = cksum + checksum (zt, km, 1)
      cksum = cksum + checksum (xu, imt, 1)
      cksum = cksum + checksum (yu, jmt, 1)
      cksum = cksum + checksum (zw, km, 1)
      cksum = cksum + checksum (dxtdeg, imt, 1)
      cksum = cksum + checksum (dytdeg, jmt, 1)
      cksum = cksum + checksum (dxudeg, imt, 1)
      cksum = cksum + checksum (dyudeg, jmt, 1)
      cksum = cksum + checksum (dzt, km, 1)
      cksum = cksum + checksum (dzw, km+1, 1)
      write (stdout,'(/)')
      write (stdout,*) 'Grid checksum = ',cksum
      write (stdout,'(/)')
      return
9001  format (1x,10f10.4)
9002  format (1x,10f10.2)
9101  format (/,  a,i4,' in units of ',a,' as follows:')
9103  format (/,' Depth to "t" & "u" grid points (cm): zt(k) k=1,',i3)
9104  format (/,' Depth to "w" grid points (cm): zw(k) k=1,',i3)
9105  format (/,' Latitude of "t" points (deg): yt(j) j=1,',i4)
9106  format (/,' Latitude of "u" points (deg): yu(j) j=1,',i4)
9107  format (/,' Longitude of "t" points (deg): xt(i) i=1,',i4)
9108  format (/,' Longitude of "u" points (deg): xu(i) i=1,',i4)
      end





      subroutine gcell (maxlen, n_bounds, bounds, d_bounds, nbpts
     &,                 num, deltat, deltau, stretch)
c
c=======================================================================
c
c              G R I D   C E L L   C O N S T R U C T I O N
c
c     A domain is composed of one or more regions:
c     Build "num" "t"  cells with resolution "deltat(n) n=1,num" 
c     within the domain composed of regions bounded by "bounds".
c     Also construct "num" "u"  cells of resolution "deltau(n) n=1,num"
c     with the relation between "t" and "u" cells given by:
#ifdef centered_t
c     deltau(n) = 0.5*(deltat(n+1) + deltat(n))
#else
c     deltat(n) = 0.5*(deltau(n-1) + deltau(n))
#endif
c     Resolution may be constant or smoothly varying within each
c     region AND there must be an integral number of grid cells within 
c     each region. The domain is the sum of all regions.
c
c     inputs:
c
c     maxlen   = maximum length of "deltat" and "deltau"
c     n_bounds = number of bounds needed to define the regions
c     bounds   = latitude, longitude, or depth at each bound
c     d_bounds = delta (resolution) at each of the "bounds"
c     nbpts    = number of extra boundary cells to add to the domain.
c                (usually one at the beginning and end)
c     stretch  = stretching factor for last region (should only be used
c                in the vertical to provide increased stretching of grid
c                points. "stretch" = 1.0 gives no increased stretching.
c                "stretch" = 1.2 gives increased stretching...etc
c             
c     outputs:
c
c     num    = total number of grid cells within the domain
c     deltau = resolution of "u" grid cells: n=1,num
c     deltat = resolution of "t" grid cells: n=1,num
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      dimension deltat(maxlen), deltau(maxlen)
      dimension d_bounds(n_bounds), bounds(n_bounds)
c
c     Set some constants
c
      p5 = 0.5
      pi = 4.0*atan(1.0)
c
c     Do all regions, one at a time, to construct the domain
c
      num  = 1
      do l=1,n_bounds-1
c
        write (stdout,'(2x,a,i2,a,g14.7,a,g14.7,a,g14.7,a,g14.7,a)')
     & ' region # ',l,'  going from ',bounds(l),' (res=',d_bounds(l)
     &,') to ',  bounds(l+1),' (res=',d_bounds(l+1),')' 
c
c       avg_res = average resolution of "t" cells within region
c       chg_res = change in resolution across the region
c       wid     = width of region
c       tol     = tolerance for fitting "t" cels within region width
c
c       provide for stretching last region if needed
c
        if (l .eq. n_bounds-1) then
          avg_res = p5*(d_bounds(l) + stretch*d_bounds(l+1))
	  chg_res = (stretch*d_bounds(l+1) - d_bounds(l))
	else
          avg_res = p5*(d_bounds(l) + d_bounds(l+1))
	  chg_res = (d_bounds(l+1) - d_bounds(l))
	endif
c
	tol = 1.e-5
	wid = abs(bounds(l+1) - bounds(l))
        an  = wid/avg_res
	m   = nint(an)
c
#ifdef centered_t
c
c       Calculate resolution of "T" cells: "deltat". Note that
c       "T" grid points are centered in these cells (as in MOM 1)
c       n = number of "t" cells fitting within the region boundaries
c
        sum = 0.0
	n   = 0
        do i = 1,100000
	  del = avg_res - p5*chg_res*cos((pi/m)*(i-0.5))
	  if (sum + del .le. wid*(1.0 + tol)) then
	    sum = sum + del
	    if (num+i-1 .gt. maxlen) then
	      write (stdout,*) "=>Error: maxlen exceeded in gcell. "
     &,                        " ...increase size of maxlen"
	      stop
	    endif
	    deltat(num+i-1) = del
	    n = n + 1
	  else
	    go to 100
	  endif
	enddo 
#else
c
c       Calculate resolution of "u" cells: "deltau"
c       "u" grid points will be centered in these cells
c       n = number of "t" cells fitting within the region boundaries
c       note: "sum" initially discounts half of the "u" cells widths
c       at the boundaries
c
        sum = 0.5*d_bounds(l) - 0.5*d_bounds(l+1)
	n   = 0
        do i = 1,100000
	  del = avg_res - p5*chg_res*cos((pi/m)*i)
	  if (sum + del .le. wid*(1.0 + tol)) then
	    sum = sum + del
	    if (num+i-1 .gt. maxlen) then
	      write (stdout,*) "=>Error: maxlen exceeded in gcell. "
     &,                        " ...increase size of maxlen"
	      stop
	    endif
	    deltau(num+i-1) = del
	    n = n + 1
	  else
	    go to 100
	  endif
	enddo 
#endif
100     continue
        if (l .eq. n_bounds-1 .and. stretch .ne. 1.0) then
          write (stdout,'(a,i3,a,f5.2)')
     & '    constructed ',n,' cells with a stretch factor of ', stretch
          write (stdout,'(/2(a,g14.7),/2(a,g14.7),/a,a/)')
     &    'Note: you specified the ocean bottom at ',bounds(l+1)
     &,   ' cm with a bottom cell thickness of ',d_bounds(l+1)
     &,   '      The stretch factor puts the bottom at ',bounds(l)+sum
#ifdef centered_t
     &,   ' cm with a bottom cell thickness of ',deltat(num+n-1)
#else
     &,   ' cm with a bottom cell thickness of '
     &,   0.5*(deltau(num+n-1) + deltau(num+n-2))
#endif
     &,   '      Adjust "stretch_z" in subroutine "gcoord" to get'
     &,   ' closer to the desired specifications if needed.'

	else
          write (stdout,'(a,g14.7,a)')
     &   '    constructed ',an,' grid cells for this region'
          if (abs(an-n) .gt. 0.01) then
            write (stdout, '(/,a,i2,/,a,g14.7/,a,g14.7,a//a/a)')
     & '==>Error: non integral number of cells in region #',l
     &,'          average resolution within region =',avg_res
     &,'          this implies ',an,' grid cells'
     &,'          Change grid specifications within USER INPUT section'
     &,'          Here is some help...'
            d_new = (2.0*wid)/(n-1) - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n-1,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            d_new = (2.0*wid)/n - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            d_new = (2.0*wid)/(n+1) - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n+1,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            stop '=>gcell'
          endif
        endif
        num = num + n
      enddo
c
c     adjust "num" to reflect the total number of cells contained in
c     all regions
c
      num = num - 1
c
      do i=1,num
#ifdef centered_t
c
c       build resolution for "u" cells: "deltau". Note that
c       variable resolution (stretched grid) implies "u" points are
c       off center as in MOM 1
c
        if (i .eq. num) then
          deltau(i) = d_bounds(n_bounds)
        else
	  deltau(i) = p5*(deltat(i+1) + deltat(i))
        endif
#else
c
c       build resolution for "T" cells: "deltat". Note that
c       variable resolution (stretched grid) implies "T" points are
c       off center
c
        if (i .eq. 1) then
          deltat(i) = p5*(d_bounds(1) + deltau(i))
        else
	  deltat(i) = p5*(deltau(i) + deltau(i-1))
        endif
#endif
      enddo     
c
c     add boundary points if needed
c
      if (nbpts .ne. 0) then
	do i=num,1,-1
	  deltat(i+1) = deltat(i) 
	  deltau(i+1) = deltau(i)
	enddo
	deltat(1)     = deltat(2)
	deltau(1)     = d_bounds(1)
	deltat(num+2) = deltat(num+1) 
	deltau(num+2) = deltau(num+1)
	num           = num + 2 
      endif
      return
      end



#if !defined driver_only
      subroutine grids
c
c=======================================================================
c     set up a staggered "B" grid for MOM and compute grid related
c     variables
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#ifndef implicitvmix
# include "accel.h"
#endif
#include "coord.h"
#include "grdvar.h"
#include "hmixc.h"
#include "scalar.h"
#include "vmixc.h"
c
c-----------------------------------------------------------------------
c     set some constants
c-----------------------------------------------------------------------
c
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      degtcm = radius/radian
c
c-----------------------------------------------------------------------
c     calculate coordinates for "t" and "u" grid cells.
c-----------------------------------------------------------------------
c
      maxlen = max(imt,jmt,km)
      call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c-----------------------------------------------------------------------
c     verify that the number of grid points match the number set in
c     the parameter statement in "size.h".
c-----------------------------------------------------------------------
c
      call size_check (imt2, jmt2, km2, 'sub grids', 'stop')
c
c-----------------------------------------------------------------------
c     convert grid resolution to cm
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        dyt(jrow) = dytdeg(jrow)*degtcm
        dyu(jrow) = dyudeg(jrow)*degtcm
      enddo
c
      do i=1,imt
        dxt(i) = dxtdeg(i)*degtcm
        dxu(i) = dxudeg(i)*degtcm
      enddo
#ifdef cyclic
      dxt(1)   = dxt(imt-1)
      dxt(imt) = dxt(2)
      dxu(1)   = dxu(imt-1)
      dxu(imt) = dxu(2)
#endif
c
c-----------------------------------------------------------------------
c     compute all quantities derived from the grid spacings
c-----------------------------------------------------------------------
c
      do k=1,km
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
      enddo
c
      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)
c
      do k=1,km
        dzwr(k-1)    = c1/dzw(k-1)
        dzw2r(k-1)   = p5/dzw(k-1)
      enddo
c
      do k=1,km
#if defined implicitvmix || defined isopycmix
        dztur(k) = c1/(dzw(k-1)*dzt(k))
        dztlr(k) = c1/(dzw(k)*dzt(k))
#endif
#ifdef tcvmix
        dzwur(k) = c1/(dzt(k)*dzw(k))
        if (k .lt. km) dzwlr(k) = c1/(dzt(k+1)*dzw(k))
#endif
        dztr(k)  = c1/dzt(k)
        rho0dztr(k) = c1/(rho0*dzt(k))
      enddo
c
      tiny = 1.e-20
      do jrow=1,jmt
        dytr(jrow)  = c1/dyt(jrow)
        dyt2r(jrow) = p5/dyt(jrow)
        dyt4r(jrow) = p25/dyt(jrow)
        dyur(jrow)  = c1/dyu(jrow)
        dyu2r(jrow) = p5/dyu(jrow)
        dyu4r(jrow) = p25/dyu(jrow)
#if defined beta_plane || defined f_plane
        phi(jrow)   = yu(1)/radian
	phit(jrow)  = yt(1)/radian
#else
        phi(jrow)   = yu(jrow)/radian
	phit(jrow)  = yt(jrow)/radian
#endif
        cst(jrow)   = cos(phit(jrow))
        csu(jrow)   = cos(phi (jrow))
        sine(jrow)  = sin(phi(jrow))
	if (cst(jrow)  .eq. 0.0) then
	  print '(/a,e14.7,a,i4,/a)'
     &  ,' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
	  cst(jrow)  = tiny
	endif
	if (csu(jrow)  .eq. 0.0) then
	  print '(/a,e14.7,a,i4,/a)'
     &  ,' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
	  csu(jrow)  = tiny
	endif
        cstr(jrow)     = c1/cst(jrow)
        csur(jrow)     = c1/csu(jrow)
        tng(jrow)      = sine(jrow)/csu(jrow)
	cstdytr(jrow)  = c1/(cst(jrow)*dyt(jrow))
	cstdyt2r(jrow) = cstdytr(jrow)*p5
        rho0csudyur(jrow)  = c1/(rho0*csu(jrow)*dyu(jrow))
        csudyu2r(jrow) = p5/(csu(jrow)*dyu(jrow))
        cst_dytr(jrow) = cst(jrow)/dyt(jrow)
        csu_dyur(jrow) = csu(jrow)/dyu(jrow)
c
      enddo
c
      do i=1,imt
        dxtr(i)  = c1/dxt(i)
        dxt2r(i) = p5/dxt(i)
        dxt4r(i) = p25/dxt(i)
        dxur(i)  = c1/dxu(i)
        dxu2r(i) = p5/dxu(i)
        dxu4r(i) = p25/dxu(i)
      enddo
c
      do i=2,imtm1
        dxmetr(i) = c1/(dxt(i) + dxt(i+1))
      enddo
c
      do i=1,imt
        duw(i) = (xu(i) - xt(i))*degtcm
      enddo
      do i=1,imtm1
	due(i) = (xt(i+1) - xu(i))*degtcm
      enddo
# ifdef cyclic
      due(imt) = due(2)
# else
      due(imt) = due(imtm1)
# endif
c
      do jrow=1,jmt
        dus(jrow) = (yu(jrow) - yt(jrow))*degtcm
      enddo
c
      do jrow=1,jmtm1
	dun(jrow) = (yt(jrow+1) - yu(jrow))*degtcm
      enddo
      dun(jmt) = dun(jmtm1)
c
# ifndef implicitvmix
c
c     for convection code, compute values needed to include effects
c     of tracer timestep acceleration on effective layer thicknesses
c
      do k=1,km
        dztxcl(k) = dzt(k)/dtxcel(k)
      enddo
c
      do k=1,kmm1
        dzwxcl(k) = c1/(dztxcl(k)+dztxcl(k+1))
      enddo
      dzwxcl(km) = c0
c
# endif
      return
      end
#endif


