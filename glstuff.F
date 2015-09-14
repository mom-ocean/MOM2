

#ifdef sgi_display
c-----------------------------------------------------------------------
c
c
c     sgi graphics stuff
c
c     author: c. h. goldberg          e.mail==>chg@gfdl.gov
c
c-----------------------------------------------------------------------
c
      subroutine initdisplay

c     colors uniformly distributed by depth: range = [-zw(km)..0.0]
#define depth_colors
c     colors uniformly distributed by kmt: range = [km..0]
c#define kmt_colors
c     lighten all colors by adding white [may print better]
c#define pastel_colors

#include "glstuff.h"
c#include "/usr/include/fdevice.h"
c#include "/usr/include/fgl.h"
c#include "/usr/include/fget.h"
#include "fdevice.h"
#include "fgl.h"
c#include "fget.h"
      integer *2 data
      integer*4 p(2)
#include "size.h"
      character *32 title
c
c-----------------------------------------------------------------------
c
c     initialize graphics display server
c
c      these next two initialization lines are currently unnecessary,
c      but worth remembering.  they open a pseudo remote DGL connection.
c
c      iserv = dglope ('chg', 3, DGLTSO)
c      print *, 'iserv=',iserv
c
c     If the next GL call is not made, the entire executable becomes
c     a background job at the first GL library call, and commands in
c     the run script after the execute command begin executing before
c     the run is completed.
c
      call foregr
c
      call initscolor
      nrwindows = 0
c
      imainmap = 1
      i = imainmap
      iwd(i) = 1280
      iht(i) = iwd(i)/2

c     make main map window (360 x 180) degrees

      ix1 = 0
      iy1 = 10
      ix2 = ix1 + iwd(i) - 1
      iy2 = iy1 + iht(i) - 1
      title = 'mainmap'
      call makewindow (imainmap, title, midgray, 'prefpo',
     &                 0.0, 360.0, -90.0, 90.0,
     &                 ix1, ix2, iy1, iy2)
      icurrentwin = iwin(imainmap)
c
c     make interation window
c
      iterwin = 2
      title = 'i'
      iwidth = 90
      iheight = 30
      call makewindow (iterwin, title, midgray, 'prefpo',
     &                 -0.5, 0.5, -0.5, 0.5,
     &                 ix2-100,ix2-100+iwidth, iy2+120,iy2+120+iheight)
c
c     make a color bar
c
      icolorbar = 3
      idy = 30
      iy3 = iy2+idy
      iy4 = iy3 + 50
      title = 'depths'
      call makewindow (icolorbar, title, midgray, 'prefpo',
     &                 -0.5, maxscolor+0.5, -0.5, 50.5,
     &                 ix1+iwd(imainmap)/3, ix2, iy2+idy, iy2+idy+50)
      do ix = 0, maxscolor
        call cpack (scolor(ix))
        call gl_rectf (ix-0.5, -0.5, ix+0.5, 50.0+0.5)
      end do

c     make edit/zoom map window

      izoom = 4
      i = izoom
      iwd(i) = iwd(imainmap)/3
      iht(i) = iwd(i)
      title = 'edit window'
      call makewindow (izoom, title, midgray, 'prefpo',
     &                 0.0, 360.0, -90.0, 90.0,
     &                 ix1, ix1+iwd(izoom), 1023-30-iht(izoom), 1023-30)

c
c     make mouse button prompt window

      iprompt = 5
      i = iprompt
      iwd(i) = iwd(imainmap)/3
      iht(i) = iwd(i)/4
      lz = 200
      title = 'mouse buttons'
      ix5 = ix1+iwd(izoom)+20
      iy5 = iy4 + 30
      call makewindow (iprompt, title, midgray, 'prefpo',
     &                 0.5, 3.5, -0.5, 0.5,
     &                 ix5, ix5+iwd(i), iy5, iy5+50)

c
c     make "leave interactive edit" button window

      iquit = 6
      i = iquit
      iwd(i) = iwd(iprompt)/3
      iht(i) = iht(iprompt)
      title = 'quit'
      ix6 = ix5+iwd(iprompt)/3
      iy6 = iy5 + 80
      call makewindow (iquit, title, midgray, 'prefpo',
     &                 -0.5, 0.5, -0.5, 0.5,
     &                 ix6, ix6+iwd(i), iy6, iy6+50)
c
c     define window numbers for two additional map windows
c     but do not make them yet
c
c
      call clear_edit_window
c
      call winset (iwin(imainmap))
c
c     enable input devices  [a winope call must precede these]
c
      call qdevic (KEYBD)
      call qdevic (LEFTMO)
      call qdevic (RIGHTM)
      call qdevic (MIDDLE)
      call qdevic (MOUSEX)
      call qdevic (MOUSEY)
      leftmouse   = .false.
      middlemouse = .false.
      rightmouse  = .false.
      mx = 0
      my = 0
c
#ifdef debug_glstuff
c     initialization complete, now test features
c
c     get window origin
      call getori (ixo, iyo)
      print *, 'window origin=', ixo, iyo
c     get number of bit planes
      nbp = getpla ()
      print *, 'number of bit planes = ',nbp
      call cmov2i (175, 87)
      call cpack (scolor(2*ngray))
      call charst ('hello', 5)
#endif

      return
      end
c-----------------------------------------------------------------------
c
      subroutine reldisplay

#include "glstuff.h"
c
      do i=imainmap,iquit
        call winclo (iwin(imainmap))
      end do

      return
      end


c-----------------------------------------------------------------------
c
      subroutine make_extra_map_windows

#include "glstuff.h"
      character *32 title
c
      call winclo(izoom)
      call winclo(iprompt)
      call winclo(iquit)

c     make second map window

      imap2 = 7
      i = imap2
      iwd(i) = iwd(imainmap)*45/100
      iht(i) = iwd(i)/2
c
      ix1 = 0
      iy1 = 10
      ix7 = ix1
      iy7 = iy1 + iht(imainmap) + 30
      ix8 = ix7 + iwd(i) - 1
      iy8 = iy7 + iht(i) - 1
      title = 'map2'
      call makewindow (imap2, title, midgray, 'prefpo',
     &                 0.0, 360.0, -90.0, 90.0,
     &                 ix7, ix8, iy7, iy8)

c     make third map window

      imap3 = 8
      i = imap3
      iwd(i) = iwd(imap2)
      iht(i) = iwd(i)/2
c
      ix9 = ix8 + 15
      iy9 = iy7
      ix10 = ix9 + iwd(i) - 1
      iy10 = iy9 + iht(i) - 1
      title = 'map3'
      call makewindow (imap3, title, midgray, 'prefpo',
     &                 0.0, 360.0, -90.0, 90.0,
     &                 ix9, ix10, iy9, iy10)

      return
      end
c-----------------------------------------------------------------------


      subroutine initscolor
c
c     make colors for a color bar
c
#include "glstuff.h"

      integer rcolor, gcolor, bcolor
      parameter (rcolor=1, gcolor=256, bcolor=256*256)
      integer r, g, b, maxc
      parameter (maxc=255)
      parameter (gamma=real(maxc)/ngray)
      parameter (ncolorsteps=20)

      do is = 0,ngray-1
        call rgbtopack (0,
     &                  nint(is * gamma),
     &                  nint(ngray * gamma),
     &                  scolor (is))
      end do
      do is = 0,ngray-1
        call rgbtopack (0,
     &                  nint(ngray * gamma),
     &                  nint((ngray-is) * gamma),
     &                  scolor (is+ngray))
      end do
      do is = 0,ngray-1
        call rgbtopack (nint(is * gamma),
     &                  nint(ngray * gamma),
     &                  0,
     &                  scolor (is+2*ngray))
      end do
      do is = 0,ngray-1
        call rgbtopack (nint(ngray * gamma),
     &                  nint((ngray-is) * gamma),
     &                  0,
     &                  scolor (is+3*ngray))
      end do
      call rgbtopack (nint(ngray * gamma),
     &                0,
     &                0,
     &                scolor (4*ngray))

c
c     create slightly darker contour bands
c
      do is = 0,maxscolor
        if (mod((is*ncolorsteps)/maxscolor,2) .eq. 0) then
          call packtorgb (scolor(is), r, g, b)
          call rgbtopack (r*9/10, g*9/10, b*9/10, scolor(is))
        end if
      end do

#ifdef pastel_colors
      do is = 0,maxscolor
        call packtorgb (scolor(is), r, g, b)
c       call rgbtopack ((r+maxc)/2, (g+maxc)/2, (b+maxc)/2, scolor(is))
        call rgbtopack ((2*r+maxc)/3, (2*g+maxc)/3, (2*b+maxc)/3,
     &                   scolor(is))
      end do
#endif

      return
      end

c-----------------------------------------------------------------------

      subroutine rgbtopack (r, g, b, rgb)

      integer r, g, b, rgb
      integer rcolor, gcolor, bcolor
      parameter (rcolor=1, gcolor=256, bcolor=256*256)

      rgb = r * rcolor + g * gcolor + b * bcolor

      return
      end

c-----------------------------------------------------------------------

      subroutine packtorgb (rgb, r, g, b)

      integer r, g, b, rgb
      integer rcolor, gcolor, bcolor
      parameter (rcolor=1, gcolor=256, bcolor=256*256)

      r = mod (rgb, 256)
      g = mod (rgb/256, 256)
      b = mod (rgb/(256*256), 256)

      return
      end

c-----------------------------------------------------------------------

      subroutine plotpoint (ix, iy, icolor)

      call cpack (icolor)
      call gl_rectf (ix-0.5, iy-0.5, ix+0.5, iy+0.5)

      return
      end

c-----------------------------------------------------------------------

      function ilcolor (v, vmin, vmax, cmin, cmax)
c
c     select a color based on log(abs(v))
c
      integer cmin, cmax
#include "glstuff.h"

        if (v .eq. 0) then
          call rgbtopack (64, 32, 0, ilcolor)
        else
          is = maxscolor * (log(abs(v))-log(vmin))/
     &                     (log(vmax)  -log(vmin))
          is = min (is, maxscolor)
c         is = max (is, 0)
          if (is .ge. 0) then
            ilcolor = scolor (mod(is,maxscolor+1))
          else
            is = max(is, -128)
            call rgbtopack (0,0,255+is, ilcolor)
c           ilcolor = 0
          end if
        end if

      return
      end
c-----------------------------------------------------------------------

      function iscolor (v, vmin, vmax, cmin, cmax)
c
c     select a color based on v
c
      integer r, g, b
      integer cmin, cmax

#include "glstuff.h"

        if (v .gt. vmax) then
          is = cmax
          iscolor = scolor(cmax)
c         call rgbtopack (255,192,192, iscolor)
        else if (v .lt. vmin) then
          is = cmin
          iscolor = scolor(cmin)
c         iscolor = 0
        else if (v .eq. 0.0) then
           iscolor = iwhite
        else
          is = cmin + (cmax-cmin) * (v-vmin)/(vmax-vmin)
          iscolor = scolor (mod(is,maxscolor+1))
        end if

      return
      end
c-----------------------------------------------------------------------

      subroutine set_kmt_color(k, km, zw)
c
c     set a GL color based on depth [either k or zw(k)]
c
#include "glstuff.h"
      dimension zw(km)
c
      if (k .ne. 0) then
#if defined kmt_colors || !defined depth_colors
        call cpack(iscolor(real(k), 0.0, real(km), maxscolor,0))
#else
        call cpack(iscolor(-zw(k), -zw(km), 0.0, 0, maxscolor))
#endif
       else
         call RGBcol(255,255,255)
        end if
c
       return
       end

c-----------------------------------------------------------------------

      subroutine waitmouse
c
c     wait for any mouse button click
c     left   mouse button sets logical "leftmouse"
c     middle mouse button sets logical "middlemouse"
c     right  mouse button sets logical "rightmouse"
c     the other mouse logicals are set to false
c     screen pixel coordinates of cursor are returned in (mx, my)
c     icurrentwin has sgi window number (not makewindow number)
c
#include "fdevice.h"
#include "fgl.h"
c#include "fget.h"

#include "glstuff.h"

      integer * 2 data

      call qreset
      leftmouse   = .false.
      middlemouse = .false.
      rightmouse  = .false.

 999  continue
      idev = qtest()
      if (idev .eq. 0) goto 999
 888  continue
      idev = qread (data)
#ifdef debug_mouse2
      print *, 'idev=',idev, ' data=',data
#endif
      if      (idev .eq. MOUSEX) then
        mx = data
#ifdef debug_mouse2
        print '(a,i4,a,i4,a)', 'MOUSE = (', mx,',',my,')'
#endif
      else if (idev .eq. MOUSEY) then
        my = data
#ifdef debug_mouse2
        print '(a,i4,a,i4,a)', 'MOUSE = (', mx,',',my,')'
#endif
      else if (idev .eq. LEFTMO .and. data .eq. 1) then
        leftmouse = .true.
#ifdef debug_mouse
        print *, 'LEFT MOUSE'
#endif
      else if (idev .eq. MIDDLE .and. data .eq. 1) then
        middlemouse = .true.
#ifdef debug_mouse
        print *, 'MIDDLE MOUSE'
#endif
      else if (idev .eq. RIGHTM .and. data .eq. 1) then
        rightmouse = .true.
#ifdef debug_mouse
        print *, 'RIGHT MOUSE'
#endif
      else if (idev .eq. INPTCH) then
        icurrentwin = data
#ifdef debug_mouse
        print '(a,i4)', 'Entered window ',icurrentwin
#endif
      end if
      if (.not. (leftmouse .or. middlemouse .or. rightmouse)) goto 999

      call qreset

      return
      end

c-----------------------------------------------------------------------

      subroutine sensemouse
c
c     Senses if one or more mouse buttons have been pressed.
c     Flips all mouse buttons which have been pressed.
c     Returns which mouse buttons have been pressed, but only after
c     all pressed buttons are released.
c
#include "fdevice.h"
#include "fgl.h"
c#include "fget.h"

#include "glstuff.h"

      integer * 2 data
      logical leftmouse_released, middlemouse_released
      logical rightmouse_released
      logical leftmouse_pressed, middlemouse_pressed
      logical rightmouse_pressed

      leftmouse_pressed   = .false.
      middlemouse_pressed = .false.
      rightmouse_pressed  = .false.
      leftmouse_released   = .false.
      middlemouse_released = .false.
      rightmouse_released  = .false.

c     *****  THIS SUBROUTINE HAS BEEN FIXED BUT NOT TESTED *****

#define stopmouse_OK
#ifdef stopmouse_OK

  999 continue
      idev = qtest ()

      if (idev .ne. 0) then
        idev = qread (data)
        if (idev .eq. RIGHTM) then
          if (data .eq. 1) then
            rightmouse_pressed = .true.
            rightmouse = .not. rightmouse
          else if (rightmouse) then
            rightmouse_released = .true.
          end if
        else if (idev .eq. MIDDLE) then
          if (data .eq. 1) then
            middlemouse_pressed = .true.
            middlemouse = .not. middlemouse
          else if (middlemouse) then
            middlemouse_released = .true.
          end if
        else if (idev .eq. LEFTMO) then
          if (data .eq. 1) then
            leftmouse_pressed = .true.
            leftmouse = .not. leftmouse
          else if (leftmouse) then
            leftmouse_released = .true.
          end if
        end if
c       goto 999
      else
c
c       there is nothing in the event queue
c       if a mouse button has been pushed, wait for it to be released
c

        if (leftmouse_pressed .and. .not. leftmouse_released) then
          goto 999
        end if
        if (middlemouse_pressed .and. .not. middlemouse_released) then
          goto 999
        end if
        if (rightmouse_pressed .and. .not. rightmouse_released) then
          goto 999
        end if
      end if

#endif
      return
      end

c-----------------------------------------------------------------------

      subroutine makewindow (iwinnr, title, icolor, pref,
     &                       x1,  x2,  y1,  y2,
     &                       ix1, ix2, iy1, iy2)

#include "size.h"
#include "glstuff.h"
      character * 32 title
      character * 6 pref
      integer*4 winope

      nrwindows = max(iwinnr, nrwindows)

      iwd(iwinnr) = ix2-ix1+1
      iht(iwinnr) = iy2-iy1+1
      if (pref .eq. 'prefpo') then
        call prefpo (ix1, ix2, iy1, iy2)
      else
        call prefsi (iwd(iwinnr), iht(iwinnr))
      end if

      win_name(iwinnr) = title
      iwin(iwinnr) = winope(title, len_trim(title))
      call wintit(title, len_trim(title))
      xl(iwinnr)=x1
      xr(iwinnr)=x2
      yl(iwinnr)=y1
      yr(iwinnr)=y2
      call ortho2 (xl(iwinnr), xr(iwinnr), yl(iwinnr), yr(iwinnr))
c     I need real*4 arguments to ortho--since xl, xr, yl, yr are
c     declared real*4, this trick keeps them as real*4, even when the
c     compiler directive -r8 is invoked
      call RGBmod
      call gconfi
      call cpack(icolor)
      call clear
c
#ifdef debug_glstuff
      print '(a,i2,a,i2,a,a)', 
     &      'tp: made window ',iwinnr,'-->',iwin(iwinnr),
     &      '  title: ',title
#endif
c
      return
      end
c-----------------------------------------------------------------------


      subroutine scale_map (iwinnr,title,xu1,yu1,zw1, imt2, jmt2, km2)
c
c     choose an appropriate orthogonal projection for a map window
c     based on mom grid.
c
#include "size.h"
#include "glstuff.h"
      character *(*) title
      dimension xu1(imt), yu1(jmt), zw1(km)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)

      call size_check (imt2, jmt2, km2, 'scale_main', 'fix_size')
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
c
#ifdef debug_glstuff
      print *, 'tp:  scale_map: iwinnr=',iwinnr
      print *, 'tp:             iwin(iwinnr)=', iwin(iwinnr)
#endif
c
      call winset (iwin(iwinnr))
      call wintit (title, len_trim(title))
c
c     Check to see if highest longitudes exceed 360 degrees.
c     Change boundaries of main map to compensate.
c
      if (xu(imt) .gt. 360.0) then
        xl(iwinnr)=xu(imt) - 360.0
        xr(iwinnr)=xu(imt)
        yl(iwinnr)=yl(imainmap)
        yr(iwinnr)=yr(imainmap)
      end if
c
c     Check to see if limited area simulation is less than half the
c     globe in latitude or longitude. If so, enlarge, preserving shape.
c
      xscale = (xu(imt) - xu(0)) / 360.0
      yscale = (yu(jmt) - yu(0)) / 180.0
      xyscale = 1.00 * max(xscale, yscale)
      if (xscale .lt. 0.5 .or. yscale .lt. 0.5 .or. .true.) then
        xc = 0.5 * (xu(0) + xu(imt))
        yc = 0.5 * (yu(0) + yu(jmt))
        xl(iwinnr) = xc - 180.0 * xyscale
        xr(iwinnr) = xc + 180.0 * xyscale
        yl(iwinnr) = yc -  90.0 * xyscale
        yr(iwinnr) = yc +  90.0 * xyscale
      end if
c
c     remap selected map window
c
      call ortho2(xl(iwinnr),xr(iwinnr),yl(iwinnr),yr(iwinnr))
      call RGBmod
      call gconfi
      call RGBcol (128, 128, 128)
      call clear
c
#ifdef debug_glstuff
c  tp:
c  
      print *, 'xl,xr,yl,yr=',
     &          xl(iwinnr),xr(iwinnr),yl(iwinnr),yr(iwinnr)
#endif
c
      return
      end


c-----------------------------------------------------------------------

      subroutine clear_edit_window
#include "glstuff.h"
      character *32 opt2
      real*4 xm, ym
c
c     clear edit window
c
      call winset(iwin(izoom))
      call wintit(' ', 1)
      call cpack(iwhite)
      call clear
      call cpack(ired)
      call cmov2((xl(izoom)+xr(izoom))/2, (yl(izoom)+yr(izoom))/2)
      call charst('done', 4)
c
c     clear mouse prompt window
c
      call winset(iwin(iprompt))
      call wintit(' ', 1)
      call cpack(iwhite)
      call clear
      scale_factor = 32.0
      call RGBcol (200, 0, 100)
      opt2 = 'click any button to continue'
      xm = 2-len_trim(opt2)/scale_factor
      ym = -0.2
      call cmov2(xm, ym)
      call charst(opt2, len_trim(opt2))
c
c     clear quit free edit window
c
      call winset(iwin(iquit))
      call wintit(' ', 1)
      call cpack(iwhite)
      call clear
c
      call winset(iwin(izoom))
      return
      end

c-----------------------------------------------------------------------


      subroutine gl_announce_phase(title)
#include "glstuff.h"
      character *(*) title
      character *32 title2
c
c     change title on main map window
c
      call winset(iwin(imainmap))
      title2 = title
      call wintit(title2, len_trim(title2))
      call winset(iwin(izoom))

      return
      end

c-----------------------------------------------------------------------


      subroutine show_iter (iter, title)

#include "glstuff.h"
#include "size.h"
      character *6 str_iter, title2
      character *(*) title
      real*4 xm, ym
      

c     call sensemouse

      call winset (iwin(iterwin))
      title2 = title
      call wintit(title2, len_trim(title2))
      call RGBcol (250, 250, 0)
      it = nint(real(iter)/imt * maxscolor)
      call cpack (scolor(min(it, maxscolor)))
      call clear
      call RGBcol (0, 0, 0)
      write (str_iter, '(i6)') iter

      call cmov2i (12,8)
      xm = -0.4
      ym = -0.2
      call cmov2 (xm, ym)
      call charst (str_iter, len_trim(str_iter))

      call winset (iwin(imainmap))
      return
      end

c-----------------------------------------------------------------------
      subroutine display_vector (vec,title,amax,window)
c
c     displays a 2-dimensional field of real values.
c     using a balanced scalar plot scaled from -amax <= v <= amax
c     window must be scaled for latitude and longitude using ortho2.
c
c     inputs:
c       vec = the field to display
c       title = title to write to window frame
c       amax  = maximum absolute value used for scaling colors
c       iwinnr= window number [assigned by makewindow]
#include "size.h"
#include "glstuff.h"
      character *(*) title, window
      dimension vec (imt,jmt)

      do iwinnr=1,min(nrwindows, maxwindows)
        if (window .eq. win_name(iwinnr)) then
          call display_vec(vec,title,1.0,amax,'b',iwinnr)
        end if
      end do      
      return
      end

c-----------------------------------------------------------------------
      subroutine display_vec (vec,title,amin,amax,mode,iwinnr)
c
c     displays a 2-dimensional field of real values.
c     window must be scaled for latitude and longitude using ortho2.
c
c     inputs:
c       vec = the field to display
c       title = title to write to window frame
c       amin  = minimum absolute value used for scaling colors
c       amax  = maximum absolute value used for scaling colors
c       mode  = 's' : scalar plot scaled from  amin <= v <= amax
c             = 'b' : scalar plot scaled from -amax <= v <= amax
c             = 'l' : log plot scaled from ln(amin) <= ln(v) <= ln(amax)
c             color reversed modes
c             = 'S' : scalar plot scaled from  amin <= v <= amax
c             = 'B' : scalar plot scaled from -amax <= v <= amax
c       iwinnr= window number [assigned by makewindow]
c       xu    = longitude of u/v points
c       yu    = latitude of u/v points

#include "size.h"
#include "coord.h"
#include "glstuff.h"
      character *(*) title
      character *(*) mode
      dimension vec (imt,jmt)
      dimension xu0(0:imt), yu0(0:jmt), zw0(0:km)
      real *4 x1, y1, x2, y2

      call sensemouse
      if (middlemouse) return

      call make_xyz0(xu, yu, zw, xu0, yu0, zw0)

      call winset (iwin(iwinnr))
      call wintit (title, len_trim(title))
      do i=1,imt
        do j=jmt,1,-1
          if (mode .eq. 's') then
            icolor = iscolor (vec(i,j), amin, amax, 0, maxscolor)
          else if (mode .eq. 'b') then
            icolor = iscolor (vec(i,j), -amax, amax, 0, maxscolor)
          else if (mode .eq. 'l') then
            icolor = ilcolor (vec(i,j), amin, amax, 0, maxscolor)
          else if (mode .eq. 'S') then
            icolor = iscolor (vec(i,j), amin, amax, maxscolor, 0)
          else if (mode .eq. 'B') then
            icolor = iscolor (vec(i,j), -amax, amax, maxscolor, 0)
          else
            icolor = iscolor (vec(i,j), amin, amax, 0, maxscolor)
          end if
          call cpack(icolor)
          call gl_rectf(xu0(i-1), yu0(j-1), xu0(i), yu0(j))
        end do
      end do
      call winset (iwin(imainmap))

      return
      end

c-----------------------------------------------------------------------

      subroutine gl_rectf(x1, y1, x2, y2)
c
c     filled rectangle to work with real*8 [-r8] or real*4 arguments
c
      real *4 x1_gl, y1_gl, x2_gl, y2_gl

      x1_gl = x1
      x2_gl = x2
      y1_gl = y1
      y2_gl = y2
      call rectf(x1_gl, y1_gl, x2_gl, y2_gl)

      return
      end



c-----------------------------------------------------------------------

      subroutine gl_rect(x1, y1, x2, y2)
c
c     rectangle to work with real*8 [-r8] or real*4 arguments
c
      real *4 x1_gl, y1_gl, x2_gl, y2_gl

      x1_gl = x1
      x2_gl = x2
      y1_gl = y1
      y2_gl = y2
      call rect(x1_gl, y1_gl, x2_gl, y2_gl)

      return
      end



c-----------------------------------------------------------------------
      subroutine display_vec_3x3 (vec,title,amax,iwinnr)
c
c     displays NE quadrants of a 2-dimensional field of 3x3 arrays
c     of real values.  [like cf(imt,jmt,-1:1,-1:1)]
c     window must be scaled for latitude and longitude using ortho2.
c
c     inputs:
c       vec = the field to display
c       title = title to write to window frame
c       amax  = maximum absolute value used for scaling colors
c       iwinnr= window number [assigned by makewindow]
c       xu    = longitude of u/v points
c       yu    = latitude of u/v points

#include "size.h"
#include "coord.h"
#include "glstuff.h"
      character *(*) title
      character *(*) mode
      dimension vec (imt,jmt,-1:1,-1:1)
      dimension xu0(0:imt), yu0(0:jmt), zw0(0:km)
      real *4 x1,x2,y1,y2

      call sensemouse
      if (middlemouse) return

      call make_xyz0(xu, yu, zw, xu0, yu0, zw0)

      call winset (iwin(iwinnr))
      call cpack(iwhite)
      call clear
      call wintit (title, len_trim(title))
      do i=1,imt
        do j=jmt,1,-1
          do i1=-1,1
            do j1=-1,1
              icolor = iscolor(vec(i,j,i1,j1),-amax,amax,0,maxscolor)
              call cpack(icolor)
              x1 = 0.5*(xu0(i-1)+xu0(i))+(i1  )*(xu0(i)-xu0(i-1)) / 4.0
              x2 = 0.5*(xu0(i-1)+xu0(i))+(i1+1)*(xu0(i)-xu0(i-1)) / 4.0
              y1 = 0.5*(yu0(j-1)+yu0(j))+(j1  )*(yu0(j)-yu0(j-1)) / 4.0
              y2 = 0.5*(yu0(j-1)+yu0(j))+(j1+1)*(yu0(j)-yu0(j-1)) / 4.0
              call rectf(x1, y1, x2, y2)
            end do
          end do
        end do
      end do
      call winset (iwin(imainmap))

      return
      end



c-----------------------------------------------------------------------
      subroutine display_grid (title, xu1, yu1, zw1, imt2, jmt2, km2)

#include "size.h"
#include "glstuff.h"
      character *(*) title
      character *32 title2
      dimension xu1(imt), yu1(jmt), zw1(km)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
      real*4 p(2), zmax
      character * 4 si

      call size_check (imt2, jmt2, km2, 'display_grid', 'fix_size')
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)

      call winset (iwin(imainmap))
      call wintit (title, len_trim(title))
      call RGBcol (200, 200, 200)
      call subpix(.true.)
      igrid_space=5
      do i=1,imt,igrid_space
        call bgnlin
        p(1)=xu(i)
        p(2)=yu(1)
        call v2f(p)
        p(1)=xu(i)
        p(2)=yu(jmt)
        call v2f(p)
        call endlin
      end do
      do j=1,jmt,igrid_space
        call bgnlin
        p(1)=xu(1)
        p(2)=yu(j)
        call v2f(p)
        p(1)=xu(imt)
        p(2)=yu(j)
        call v2f(p)
        call endlin
      end do
c
c     print longitudes to window
c
      call RGBcol (200, 0, 100)
      do i=1,imt,imt/18
        call cmov2i(nint(xu(i))-3, 0)
        write (si, '(i3)') nint(xu(i))
        call charst(si, len_trim(si))
      end do
c
c     display depths
c
      call display_depths (title, xu1, yu1, zw1, imt2, jmt2, km2)

      call winset (iwin(imainmap))

      return
      end

c-----------------------------------------------------------------------
      subroutine display_depths (title, xu1, yu1, zw1, imt2, jmt2, km2)

#include "size.h"
#include "glstuff.h"
      character *(*) title
      character *32 title2
      dimension xu1(imt), yu1(jmt), zw1(km)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
      real*4 p(2), zmax
      character * 4 si

      call size_check (imt2, jmt2, km2, 'display_depths', 'fix_size')
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
c
c     display depths
c
      call winset (iwin(izoom))
      title2 = 'depths'
      call wintit (title2, len_trim(title2))
      call subpix(.true.)
      xl(izoom)=0.0
      xr(izoom)=360.0
      yl(izoom)=-zw(km)
      yr(izoom)=0.0
      call ortho2 (xl(izoom), xr(izoom), yl(izoom), yr(izoom))
      call RGBmod
      call gconfi
      call RGBcol (0, 128, 128)
      call clear
c
      do k=1,km
        call set_kmt_color(k, km, zw)
        zl = -zw(k)
        zr = -zw(k-1)
        call rectf(xl(izoom),zl, xr(izoom), zr)
      end do
c
      igrid_space=5
      call RGBcol (255, 255, 200)
      do k=0,km,igrid_space
        call bgnlin
        p(1)=xl(izoom)
        p(2)=-zw(k)
        call v2f(p)
        p(1)=xr(izoom)
        p(2)=-zw(k)
        call v2f(p)
        call endlin
      end do

      call winset (iwin(imainmap))

      return
      end

c-----------------------------------------------------------------------
      subroutine display_kmt(title, xu1, yu1, zw1, kmt, imt2, jmt2, km2)

#include "size.h"
      character * 32 title
      dimension xu1(imt), yu1(jmt), zw1(km), kmt(imt,jmt)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
      real*4 p(2)
      character * 4 si

#include "glstuff.h"
c
c     verify/fix imt, jmt, km from "size.h" via "ijk_max.h"
c
      call size_check(imt2, jmt2, km2, 'display_kmt', 'fix_size')
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
c
c     select main map window
c
      call winset (iwin(imainmap))
      call wintit (title, len_trim(title))
      call cpack(midgray)
      call clear
      call subpix(.true.)
      do i=1,imt
        do j=1,jmt
          call set_kmt_color(kmt(i,j), km, zw)
          call gl_rectf(xu(i-1), yu(j-1), xu(i), yu(j))
        end do
      end do

      call RGBcol (200, 0, 100)
      do i=1,imt,imt/18
        call cmov2i(nint(xu(i))-3, 0)
        write (si, '(i3)') nint(xu(i))
        call charst(si, len_trim(si))
      end do

      return
      end

c-----------------------------------------------------------------------
      subroutine display_zoom (title, xu1, yu1, zw1, kmt
     &,                        imid, jmid, iedit, jedit
     &,                        imt2, jmt2, km2)

#include "size.h"
#include "glstuff.h"
      character * 32 title, full_title
      dimension xu1(imt), yu1(jmt), zw1(km), kmt(imt,jmt)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
      real*4 p(2)
      character * 4 s_kmt
      logical outline
c
c     verify/fix imt, jmt, km from "size.h" via "ijk_max.h"
c
      call size_check(imt2, jmt2, km2, 'display_zoom', 'fix_size')
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
c
c     select edit/zoom window
c
      call winset (iwin(izoom))
      write (full_title, '(a,a,i3,a,i3,a)')
     &       title(1:len_trim(title)),'(',imid,',',jmid,')'
      call wintit (full_title, len_trim(full_title))
      call subpix(.true.)
c
c     select world coordinates for edit window
c
      xwidth = 10*(xu(imid)-xu(imid-1))
      yheight = xwidth * (real(iht(izoom))/real(iwd(izoom)))
      xt = 0.5 * (xu(imid-1) + xu(imid))
      yt = 0.5 * (yu(jmid-1) + yu(jmid))
      xl(izoom)=xt-0.5*xwidth
      xr(izoom)=xt+0.5*xwidth
      yl(izoom)=yt-0.5*yheight
      yr(izoom)=yt+0.5*yheight
      call ortho2 (xl(izoom), xr(izoom), yl(izoom), yr(izoom))
      call RGBmod
      call gconfi
      call RGBcol (128, 128, 128)
      call clear
      char_ht=8.0/real(iht(izoom)) * yheight
c
c     display cells in edit window
c
c     find which cells will be visible
c
      call find_visible(izoom, imid, jmid, xu1, yu1, i1, i2, j1, j2)
c
      do i=i1,i2
        do j=j1,j2
          call set_kmt_color(kmt(i,j), km, zw)
          call gl_rectf(xu(i-1), yu(j-1), xu(i), yu(j))
c
c         print kmt(i,j) on display of cell
c
          if (kmt(i,j) .ne. 0) then
            call cmov2(xu(i-1), 0.5*(yu(j-1)+yu(j))-char_ht)
            call cpack (iblack)
            write (s_kmt, '(i3)') kmt(i,j)
            call charst (s_kmt, len_trim(s_kmt))
          end if
        end do
      end do
c
c     outline edit rectangle
c
      outline = .true.
      call redisplay_edit_cell (kmt,iedit,jedit,xu1,yu1,zw1,outline)
c
c     outline zoom area on main map
c
      call outline_main_zoom (xu1, yu1, i1, i2, j1, j2)
c
      call winset(iwin(izoom))
c
#ifdef debug_mouse
      print '(a,i4,a,i4,a)','display_zoom: MOUSE=(',mx,',',my,')'
#endif
      return
      end


c-----------------------------------------------------------------------
      subroutine redisplay_main_zoom (xu1, yu1, zw1, kmt
     &,                               imid, jmid, outline)
c
c     redisplay in main window those cells that appear in zoom window
c
#include "size.h"
#include "glstuff.h"
      dimension xu1(imt), yu1(jmt), zw1(km), kmt(imt,jmt)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
      logical outline
c
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
c
c     select main map window
c
      call winset (iwin(imainmap))
      call subpix(.true.)
c
c     find which cells would be visible in edit/zoom window
c
      call find_visible(izoom, imid, jmid, xu1, yu1, i1, i2, j1, j2)
c
      do i=i1,i2
        do j=j1,j2
          call set_kmt_color(kmt(i,j), km, zw)
          call gl_rectf(xu(i-1), yu(j-1), xu(i), yu(j))
        end do
      end do
c
      if (outline) then
        call outline_main_zoom (xu1, yu1, i1, i2, j1, j2)
      end if
c
      call winset (iwin(izoom))

      return
      end
c-----------------------------------------------------------------------

      subroutine outline_main_zoom (xu1, yu1, i1, i2, j1, j2)
c
c     outline zoom area on main map
c
#include "size.h"
#include "glstuff.h"
      dimension xu1(imt), yu1(jmt)
      dimension xu(0:imt), yu(0:jmt)
      real*4 e

      call coord0(xu1, xu, imt)
      call coord0(yu1, yu, jmt)
c
      call winset(iwin(imainmap))
      call cpack(midgray)
      e = 2*(xr(izoom) - xl(izoom))/40
      call gl_rect(max(xu(i1-1), xl(izoom))+e,
     &             max(yu(j1-1), yl(izoom))+e,
     &          min(xu(i2), xr(izoom))-e, min(yu(j2), yr(izoom))-e)
      call cpack(iwhite)
      e = (xr(izoom) - xl(izoom))/40
      call gl_rect(max(xu(i1-1), xl(izoom))+e,
     &             max(yu(j1-1), yl(izoom))+e,
     &          min(xu(i2), xr(izoom))-e, min(yu(j2), yr(izoom))-e)
      call RGBcol(255, 0, 0)
      call gl_rect(max(xu(i1-1), xl(izoom)), max(yu(j1-1), yl(izoom)),
     &          min(xu(i2), xr(izoom)), min(yu(j2), yr(izoom)))
      return
      end

c-----------------------------------------------------------------------

      subroutine find_visible(iwinnr,imid,jmid,xu1,yu1,i1,i2,j1,j2)
c
c     finds a range of subscripts for display by comparing longitudes xu
c     and latitudes yu to bounds of window iwinnr.
c
c     inputs:
c       iwinnr = window number [from makewindow]
c       imid, jmid = subscripts within the range to be displayed
c       xu, yu = longitudes and latitudes corresponding to subscripts
c
c     outputs:
c       i1 = smallest i to display
c       i2 = largest i to display
c       j1 = smallest j to display
c       j2 = largest j to display
c
#include "size.h"
#include "glstuff.h"
      dimension xu1(imt), yu1(jmt)
      dimension xu(0:imt), yu(0:jmt)
      call coord0(xu1, xu, imt)
      call coord0(yu1, yu, jmt)
      do i1=imid,2,-1
        if (xu(i1-1) .le. xl(iwinnr)) then
          goto 111
        end if
      end do
  111 continue
      do i2=imid,imt-1
        if (xu(i2) .ge. xr(iwinnr)) then
          goto 222
        end if
      end do
  222 continue
      do j1=jmid,2,-1
        if (yu(j1-1) .le. yl(iwinnr)) then
          goto 333
        end if
      end do
  333 continue
      do j2=jmid,jmt-1
        if (yu(j2) .ge. yr(iwinnr)) then
          goto 444
        end if
      end do
  444 continue

c      old way
c
c      i1 = max(imid-idelta,2)
c      i2 = min(imid+,imt)
c      j1 = max(jmid-jdelta,2)
c      j2 = min(jmid+jdelta,jmt)

      return
      end


      subroutine redisplay_edit_cell (kmt,iedit,jedit,
     &                                xu1,yu1,zw1,outline)
c
c     redisplay edit cell in edit window
c
#include "size.h"
#include "glstuff.h"
#include "fgl.h"
      dimension kmt(imt,jmt), xu1(imt), yu1(jmt), zw1(km)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
      character *3 s_kmt
      character *32 slld
      character *10 szwkmt
      real *4 x1, y1, x2, y2, x3, y3
      logical outline
c
      call winset(iwin(izoom))
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
      call set_kmt_color(kmt(iedit,jedit), km, zw)
      call gl_rectf(xu(iedit-1), yu(jedit-1),
     &           xu(iedit), yu(jedit))
c
c     write digits kmt(iedit,jedit) on cell
c
      xwidth = xr(izoom) - xl(izoom)
      char_ht=8.0/real(iht(izoom)) * xwidth
      call cmov2(xu(iedit-1), 0.5*(yu(jedit-1)+yu(jedit))-char_ht)
      call cpack (iblack)
      write (s_kmt, '(i3)') kmt(iedit,jedit)
      call charst (s_kmt, len_trim(s_kmt))
c
c     outline edit cell
c
      if (outline) then
        call RGBcol(200,0,200)
        call gl_rect(xu(iedit-1), yu(jedit-1),
     &               xu(iedit), yu(jedit))
      end if
c
c     display longitude, latitude, and depth
c
      scalefactor = (xr(izoom)-xl(izoom))/(iwd(izoom)/9.0)
      max_slld = 32
      if (kmt(iedit,jedit) .ne. 0) then
        write (szwkmt, '(f7.1)') zw(kmt(iedit,jedit))/1e2
      else
        szwkmt = ' 0'
      end if
      write (slld,'(a,f5.1,a,f5.1,a,a,a)')
     &      '(', 0.5*(xu(iedit)+xu(iedit-1)), ','
     &,          0.5*(yu(jedit)+yu(jedit-1)), ')  z ='
     &,          szwkmt(1:len_trim(szwkmt)), ' m'
      x1 = (xl(izoom)+xr(izoom))/2 - scalefactor*len_trim(slld)/2.0
      y1 = yr(izoom)-0.07 * (yr(izoom)-yl(izoom))
      x2 = (xl(izoom)+xr(izoom))/2 - scalefactor*max_slld/2
      y2 = y1 - scalefactor * 0.4
      x3 = x2 + scalefactor * max_slld
      y3 = y2 + scalefactor * 2
      call cpack (iwhite)
      call rectf(x2, y2, x3, y3)
      call cmov2(x1, y1)
      call RGBcol (200, 0, 100)
      call charst(slld, len_trim(slld))
c
      return
      end


c-----------------------------------------------------------------------
      subroutine redisplay_main_cell (kmt, iedit, jedit, xu1, yu1, zw1)
c
c     redraw rectangle on main window
c
#include "size.h"
#include "glstuff.h"
#include "fgl.h"
      dimension kmt(imt,jmt), xu1(imt), yu1(jmt), zw1(km)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
c
      call winset(iwin(imainmap))
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
      call set_kmt_color(kmt(iedit,jedit), km, zw)
      call gl_rectf(xu(iedit-1), yu(jedit-1),
     &           xu(iedit), yu(jedit))
      call winset(iwin(izoom))
c
      return
      end

c-----------------------------------------------------------------------
      subroutine redisplay_option(iwinnr, kmt_opt, i_opt, xu1, yu1, zw1)
c
c     redisplay an entire option in the chosen window
c
#include "size.h"
#include "glstuff.h"
#include "fgl.h"
      parameter (max_opt=3, len_opt=10)
      dimension kmt_opt(max_opt, len_opt, 4)
      dimension xu1(imt), yu1(jmt), zw1(km)
      dimension xu(0:imt), yu(0:jmt), zw(0:km)
c
      call make_xyz0(xu1, yu1, zw1, xu, yu, zw)
c
      if (iwinnr .ne. imainmap .and. iwinnr .ne. izoom) then
        print '(a,a,i2)'
     &,  'WARNING: program bug--subroutine redisplay_option called with'
     &,  ' invalid window: ',iwinnr
        return
      end if
c
      call winset(iwin(iwinnr))
      do j_opt=1,len_opt
        if (kmt_opt(i_opt,j_opt,1) .ne. 0) then
          ii = kmt_opt(i_opt,j_opt,1)
          jj = kmt_opt(i_opt,j_opt,2)
          kk = kmt_opt(i_opt,j_opt,4)
          call set_kmt_color(kk, km, zw)
          call gl_rectf(xu(ii-1), yu(jj-1), xu(ii), yu(jj))
        end if
      end do
      call winset(iwin(izoom))
c
      return
      end

c-----------------------------------------------------------------------
      subroutine wait_cursor (iwinnr, x, y
     &,                       iselect, jselect
     &,                       xu1, yu1
     &,                       imt2, jmt2, km2)
c
c     locates cursor position, and if not in window for which we have
c     world coordinates, waits for a mouse button click in such a window
c
      dimension xu1(imt2), yu1(jmt2)

  222 continue
      call find_cursor(iwinnr, x, y,
     &                 iselect, jselect, xu1, yu1,
     &                 imt2, jmt2, km2)
      if (iwinnr .eq. 0) then
        call waitmouse
        goto 222
      end if
      return
      end

c-----------------------------------------------------------------------

      subroutine find_cursor (iwinnr, x, y
     &,                       iselect, jselect
     &,                       xu1, yu1
     &,                       imt2, jmt2, km2)
c
c     if last mouse button press was with cursor in one of the
c     glstuff-defined windows, then
c
c     outputs:
c       iwinnr = window number [given by makewindow]
c       x,y    = world coordinates of cursor as defined in makewindow
c                or call to ortho2
c       iselect, jselect = i, j of cursor in main map or edit/zoom
c                          window
c       window iwinnr is selected as the active window
c
c     if last mouse button press was not in glstuff-defined window then
c        iwinnr = 0
c        x,y,iselect,jselect = [undefined]
c
c-----------------------------------------------------------------------
#include "size.h"
#include "glstuff.h"
      dimension xu1(imt), yu1(jmt), kmt(imt,jmt)
      dimension xu(0:imt), yu(0:jmt)
c
      call size_check(imt2, jmt2, km2, 'find_cursor', 'fix_size')
      call coord0(xu1, xu, imt)
      call coord0(yu1, yu, jmt)
c
c     check if cursor is in a window for which we have x,y information
c
      iwinnr = 0
      do i=imainmap, iquit
        if (icurrentwin .eq. iwin(i)) then
          iwinnr = i
        end if
      end do

      if (iwinnr .ne. 0) then
        call winset (iwin(iwinnr))
        call getori (ixo, iyo)
        x = xl(iwinnr)+(mx-ixo)*(xr(iwinnr)-xl(iwinnr))/(iwd(iwinnr)-1)
        y = yl(iwinnr)+(my-iyo)*(yr(iwinnr)-yl(iwinnr))/(iht(iwinnr)-1)
#ifdef debug_mouse
        print '(a,2f10.2,a,i2)', 'Selected (x,y)=', x, y,
     &             ' in window', iwinnr
#endif
        if (icurrentwin .ne. iwin(iwinnr)) then
          print '(a,i4)', 'ERROR: icurrentwin = ', icurrentwin
          print *, ' iwinnr=',iwinnr,'  iwin(iwinnr)=',iwin(iwinnr)
        end if
c
c       select i in range 1..imt
c
        iselect = 1
        do iselect = 2,imt2-1
          if (xu(iselect) .ge. x) goto 100
        end do
  100   continue
c
c       select j in range 1..jmt-1
c
        jselect = 1
        do jselect = 2, jmt2-2
          if (yu(jselect) .ge. y) goto 200
        end do
  200   continue

#ifdef debug_mouse
        print '(a,i4,a,i4,a)','Selected (',iselect,',',jselect,')'
c       print '(a,i4,a,i4,a)','MOUSE=(',mx,',',my,')'
#endif
      end if

      return
      end



c-----------------------------------------------------------------------
      subroutine gl_prompt(opt1, opt2, opt3, quit_opt, exit_opt)
c-----------------------------------------------------------------------
c     show user options in interactive editing of topography
c-----------------------------------------------------------------------
      character *(*) opt1, opt2, opt3, quit_opt, exit_opt
      real*4 xm, ym
#  include "glstuff.h"

      call winset(iwin(iprompt))
      call cpack(iwhite)
      call clear
      scale_factor = 32.0
      ym = -0.2
      call RGBcol (200, 0, 100)
      xm = 1-len_trim(opt1)/scale_factor
      call cmov2(xm, ym)
      call charst(opt1, len_trim(opt1))
      xm = 2-len_trim(opt2)/scale_factor
      call cmov2(xm, ym)
      call charst(opt2, len_trim(opt2))
      xm = 3-len_trim(opt3)/scale_factor
      call cmov2(xm, ym)
      call charst(opt3, len_trim(opt3))
c
      call winset(iwin(iquit))
      call cpack(iwhite)
      call clear
      call RGBcol (200, 0, 100)
      xm = -len_trim(quit_opt)/scale_factor
      call cmov2(xm, ym)
      call charst(quit_opt, len_trim(quit_opt))
c
      call winset(iwin(iterwin))
      call cpack(iwhite)
      call clear
      call RGBcol (200, 0, 100)
      xm = -len_trim(exit_opt)/scale_factor
      call cmov2(xm, ym)
      call charst(exit_opt, len_trim(exit_opt))
c
      call winset(iwin(izoom))

      return
      end

c-----------------------------------------------------------------------
      subroutine make_xyz0 (xu1, yu1, zw1, xu0, yu0, zw0)
c-----------------------------------------------------------------------
c     extend longitude, latitude, and depth arrays to subscript 0
c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------

      subroutine coord0 (xu1, xu0, imt)
c-----------------------------------------------------------------------
c     extend one coordinate axes to subscript 0
c-----------------------------------------------------------------------
      dimension xu1(1:imt), xu0(0:imt)
      do i=1,imt
        xu0(i) = xu1(i)
      end do
      xu0(0) = 2 * xu1(1) - xu1(2)
      return
      end

c-----------------------------------------------------------------------
#ifndef len_trim_exists
# define len_trim_exists
c=====begin file len_trim.F=============================================

      function len_trim (string)

      character * (*) string
      integer len_trim, i

      len_trim = 0
      do i = len (string), 1, -1
        if (string(i:i) .NE. ' ') then
          len_trim = i
          goto 999
        end if
      end do
  999 continue

      return
      end

c=====end file len_trim.F=============================================
#endif
#else
      subroutine glstuff
      return
      end
#endif



