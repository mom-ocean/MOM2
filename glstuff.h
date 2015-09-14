c-----------------------------------------------------------------------
c
c
c
c     graphics stuff
c
      parameter (ngray=255)
      parameter (maxscolor=4*ngray)
      parameter (maxwindows=24)
      parameter (igray=$010101)
      parameter (iblack=0, iwhite=256*256*256-1, midgray=128*igray)
      parameter (ired=255, igreen=255*256, iblue=255*256*256)
      integer scolor
      dimension scolor(0:maxscolor)
      logical rightmouse, leftmouse, middlemouse
      real *4 xl, xr, yl, yr, zl, zr
      character *32 win_name
      common /glstuff/scolor, nrwindows,
     &                iwin(1:maxwindows), win_name(1:maxwindows),
     &                iht(1:maxwindows), iwd(1:maxwindows), 
     &                xl(1:maxwindows), xr(1:maxwindows), 
     &                yl(1:maxwindows), yr(1:maxwindows),
     &                imainmap, izoom, icolorbar, iterwin, icurrentwin,
     &                iprompt, iquit, imap2, imap3
      common /mousestuff/ leftmouse, middlemouse, rightmouse,
     &                    mx, my
c
c
c
c
c-----------------------------------------------------------------------
