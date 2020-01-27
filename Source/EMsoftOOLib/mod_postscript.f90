! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list 
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this 
!        list of conditions and the following disclaimer in the documentation and/or 
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
!        of its contributors may be used to endorse or promote products derived from 
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

module mod_postscript
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! A collection of postscript output routines used to create a variety of graphics output
  !!
  !! This class contains a number of Postscript generating subroutines
  !! which can be used to generate PS-drawings from a Fortran program.
  !! The routines are based on the Pascal version written by J. Fransaer
  !! of the Catholic University of Leuven, Belgium.
  !! Translated into Fortran-90 by M. De Graef on 12/31/92.
  !!
  !! A typical program would start as follows\n
  !!
  !!       call self%openfile(scale)\n
  !!       call self%drawframe(...)\n
  !!       ....
  !!
  !! and end with \n
  !!
  !!       call self%closefile\n
  !!
  !! The code to draw spheres (sp) is taken from the appendix of 
  !! Earl J. Kirklands book on Advanced Computing in Electron Microscopy,
  !! and slightly modified to include color.  
  !! 
  !! All dimensions are in inches 
  !! pspage      = pagenumber for multi-page files
  !! psfigwidth  = width of drawing (pagewidth - 2 inches for margins = 6.5)
  !! psfigheight = height of drawing (pageheight - 2 = 9.0)
  !! psscale     = scale factor for overall file
  !! psname      = PostScript file name
  !! psdash      = used to store a dash pattern
  !! fonts       = array with standard PostScript fonts (more may be added)


use mod_kinds
use mod_global

IMPLICIT NONE

private

! the following postscript preamble is inspired on the one from the old EMS package
! (before it was converted to JEMS) written by Pierre Stadelmann
character(55),parameter :: PSpreamble(23) = (/ &
        "%!PS-Adobe-3.0                                         ", &
        "%%Creator:                                             ", &
        "%%Title:                                               ", &
        "%%Pages: (atend)                                       ", &
        "%%EndComments                                          ", &
        "/M {moveto} def /N {newpath} def /L {lineto} def       ", &
        "/S {stroke} def /T {translate} def /R {rotate} def     ", &
        "/F {fill} def /Cl {closepath} def                      ", &
        "/circle {N 0 0 1 0 360 arc Cl F} def                   ", &
        "/sp { gsave T scale 1.0 -0.04 0 { pop 3 array astore   ", &
        "{1.02 mul} forall 3 copy setrgbcolor -0.025 0.030 T    ", &
        "circle 0.93 0.93 scale } for                           ", &
        "pop pop pop grestore } def                             ", &
        "/frame {1.0 setgray N left rad add bottom M            ", &
        "right bottom right top rad arcto L right top left top  ", &
        "rad arcto L left top left bottom rad arcto L left      ", &
        "bottom right bottom rad arcto L Cl F 0.0 setgray N     ", &
        "left rad add bottom M right bottom right top rad       ", &
        "arcto L right top left top rad arcto L left top left   ", &
        "bottom rad arcto L left bottom right bottom rad arcto  ", &
        "L Cl S } def                                           ", &
        "%%EndProlog                                            ", &
        "72 dup scale                                           " /)
!DEC$ ATTRIBUTES DLLEXPORT :: PSpreamble


! font-related stuff
character(25),public,parameter :: PSlbl = "Written by MDG, 2001-2020"
character(20),public,parameter :: PSfonts(5) = (/"Symbol              ", &
                                          "Times-Bold          ", &
                                          "Times-BoldItalic    ", &
                                          "Times-Italic        ", &
                                          "Times-Roman         "/)
!DEC$ ATTRIBUTES DLLEXPORT :: PSlbl
!DEC$ ATTRIBUTES DLLEXPORT :: PSfonts


!> Shannon-Prewitt ionic radii in nanometer
real(kind=sgl), public, dimension(98) :: ATOM_SPradii=(/0.010,0.010,0.680,0.300,0.160,0.150,0.148,0.146,0.133,0.500, &
                                                0.098,0.065,0.450,0.380,0.340,0.190,0.181,0.500,0.133,0.940, &
                                                0.068,0.060,0.740,0.690,0.670,0.640,0.630,0.620,0.720,0.740, &
                                                0.113,0.073,0.580,0.202,0.196,0.500,0.148,0.110,0.880,0.770, &
                                                0.067,0.068,0.500,0.500,0.500,0.860,0.126,0.970,0.132,0.930, &
                                                0.076,0.222,0.219,0.500,0.167,0.129,0.104,0.111,0.500,0.108, &
                                                0.050,0.104,0.500,0.970,0.500,0.990,0.500,0.960,0.500,0.940, &
                                                0.050,0.050,0.680,0.600,0.520,0.500,0.500,0.500,0.137,0.112, &
                                                0.140,0.132,0.740,0.230,0.227,0.500,0.175,0.137,0.111,0.990, &
                                                0.090,0.083,0.500,0.108,0.500,0.500,0.500,0.500/)
!DEC$ ATTRIBUTES DLLEXPORT :: ATOM_SPradii

!> atomic (metallic) radii in nanometer (0.100 if not known/applicable)
real(kind=sgl), public, dimension(98) :: ATOM_MTradii=(/0.100,0.100,0.156,0.112,0.100,0.100,0.100,0.100,0.100,0.100, &
                                                0.191,0.160,0.142,0.100,0.100,0.100,0.100,0.100,0.238,0.196, &
                                                0.160,0.146,0.135,0.128,0.136,0.127,0.125,0.124,0.128,0.137, &
                                                0.135,0.139,0.125,0.116,0.100,0.100,0.253,0.215,0.181,0.160, &
                                                0.147,0.140,0.135,0.133,0.134,0.137,0.144,0.152,0.167,0.158, &
                                                0.161,0.143,0.100,0.100,0.270,0.224,0.187,0.182,0.182,0.181, &
                                                0.100,0.100,0.204,0.178,0.177,0.175,0.176,0.173,0.174,0.193, &
                                                0.173,0.158,0.147,0.141,0.137,0.135,0.135,0.138,0.144,0.155, &
                                                0.171,0.174,0.182,0.168,0.100,0.100,0.100,0.100,0.100,0.180, &
                                                0.163,0.154,0.150,0.164,0.100,0.100,0.100,0.100/)
!DEC$ ATTRIBUTES DLLEXPORT :: ATOM_MTradii

!> atom colors for PostScript drawings
character(3), public, dimension(98) :: ATOM_color=(/'blu','grn','blu','blu','red','bro','blu','red','grn','grn', &
                                            'blu','pnk','grn','blu','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','blu','grn','grn', &
                                            'blu','blu','grn','red','pnk','cyn','blu','grn'/)

real(kind=sgl), public, dimension(3,92) :: ATOM_colors = reshape( (/ &
                                              0.90000,0.90000,0.15000, &
                                              0.00000,0.90000,0.15000, &
                                              0.32311,0.53387,0.69078, &
                                              0.61572,0.99997,0.61050, &
                                              0.53341,0.53341,0.71707, &
                                              0.06577,0.02538,0.00287, &
                                              0.50660,0.68658,0.90000, &
                                              0.90000,0.00000,0.00000, &
                                              0.09603,0.80000,0.74127, &
                                              0.90000,0.12345,0.54321, &
                                              0.78946,0.77423,0.00002, &
                                              0.41999,0.44401,0.49998, &
                                              0.09751,0.67741,0.90000, &
                                              0.00000,0.00000,0.90000, &
                                              0.53486,0.51620,0.89000, &
                                              0.90000,0.98070,0.00000, &
                                              0.94043,0.96999,0.37829, &
                                              0.33333,0.00000,0.33333, &
                                              0.65547,0.58650,0.69000, &
                                              0.36245,0.61630,0.77632, &
                                              0.98804,0.41819,0.90000, &
                                              0.31588,0.53976,0.67494, &
                                              0.83998,0.08402,0.67625, &
                                              0.90000,0.00000,0.60000, &
                                              0.90000,0.00000,0.60000, &
                                              0.71051,0.44662,0.00136, &
                                              0.00000,0.00000,0.68666, &
                                              0.22939,0.61999,0.60693, &
                                              0.78996,0.54162,0.14220, &
                                              0.52998,0.49818,0.50561, &
                                              0.74037,0.90000,0.18003, &
                                              0.66998,0.44799,0.19431, &
                                              0.53341,0.53341,0.71707, &
                                              0.92998,0.44387,0.01862, &
                                              0.96999,0.53349,0.24250, &
                                              0.25000,0.75000,0.50000, &
                                              0.90000,0.00000,0.60000, &
                                              0.00000,0.90000,0.15256, &
                                              0.00000,0.00000,0.90000, &
                                              0.00000,0.90000,0.00000, &
                                              0.50660,0.68658,0.90000, &
                                              0.35003,0.52340,0.90000, &
                                              0.80150,0.69171,0.79129, &
                                              0.92998,0.79744,0.04651, &
                                              0.90000,0.06583,0.05002, &
                                              0.00002,0.00005,0.76999, &
                                              0.09751,0.67741,0.90000, &
                                              0.55711,0.50755,0.90000, &
                                              0.22321,0.72000,0.33079, &
                                              0.29000,0.26679,0.28962, &
                                              0.53999,0.42660,0.45495, &
                                              0.90000,0.43444,0.13001, &
                                              0.42605,0.14739,0.66998, &
                                              0.13001,0.90000,0.24593, &
                                              0.90000,0.00000,0.60000, &
                                              0.74342,0.39631,0.45338, &
                                              0.72000,0.24598,0.15122, &
                                              0.00000,0.00000,0.90000, &
                                              0.90000,0.44813,0.23003, &
                                              0.20000,0.90000,0.11111, &
                                              0.90000,0.00000,0.00000, &
                                              0.99042,0.02402,0.49194, &
                                              0.90000,0.00000,0.60000, &
                                              0.66998,0.44799,0.19431, &
                                              0.23165,0.09229,0.57934, &
                                              0.90000,0.87648,0.81001, &
                                              0.00000,0.20037,0.64677, &
                                              0.53332,0.53332,0.53332, &
                                              0.15903,0.79509,0.98584, &
                                              0.15322,0.99164,0.95836, &
                                              0.18293,0.79933,0.59489, &
                                              0.83000,0.09963,0.55012, &
                                              0.34002,0.36210,0.90000, &
                                              0.46000,0.19898,0.03679, &
                                              0.06270,0.56999,0.12186, &
                                              0.18003,0.24845,0.90000, &
                                              0.33753,0.30100,0.43000, &
                                              0.25924,0.25501,0.50999, &
                                              0.90000,0.79663,0.39001, &
                                              0.47999,0.47207,0.46557, &
                                              0.70549,0.83000,0.74490, &
                                              0.38000,0.32427,0.31919, &
                                              0.62942,0.42309,0.73683, &
                                              0.90000,0.00000,0.00000, &
                                              0.50000,0.33333,0.33333, &
                                              0.72727,0.12121,0.50000, &
                                              0.90000,0.00000,0.00000, &
                                              0.46310,0.90950,0.01669, &
                                              0.66667,0.66666,0.00000, &
                                              0.14893,0.99596,0.47105, &
                                              0.53332,0.53332,0.53332, &
                                              0.47773,0.63362,0.66714 /), (/3,92/))
!DEC$ ATTRIBUTES DLLEXPORT :: ATOM_colors


  type, public :: PostScript_T
   private 
     integer(kind=irg)      :: pspage
     integer(kind=irg)      :: psunit = 20
     integer(kind=irg)      :: imanum
     real(kind=sgl)         :: psdash(20)
     real(kind=sgl)         :: psfigwidth
     real(kind=sgl)         :: psfigheight
     real(kind=sgl)         :: psscale
     character(fnlen)       :: psname

   contains
    private 
      procedure, pass(self) :: openFile_
      procedure, pass(self) :: closeFile_
      procedure, pass(self) :: newpage_
      procedure, pass(self) :: cellinfo_
      procedure, pass(self) :: clippath_
      procedure, pass(self) :: translate_
      procedure, pass(self) :: move_
      procedure, pass(self) :: draw_
      procedure, pass(self) :: line_gray_
      procedure, pass(self) :: setlinewidth_
      procedure, pass(self) :: square_
      procedure, pass(self) :: filledsquare_
      procedure, pass(self) :: cross_
      procedure, pass(self) :: sphere_
      procedure, pass(self) :: arc_
      procedure, pass(self) :: circle_
      procedure, pass(self) :: filledcircle_
      procedure, pass(self) :: drawframe_
      procedure, pass(self) :: drawrect_
      procedure, pass(self) :: line_
      procedure, pass(self) :: setdash_
      procedure, pass(self) :: closepathS_
      procedure, pass(self) :: stroke_
      procedure, pass(self) :: gsave_
      procedure, pass(self) :: grestore_
      procedure, pass(self) :: closepath_
      procedure, pass(self) :: newpath_
      procedure, pass(self) :: text_
      procedure, pass(self) :: textv_
      procedure, pass(self) :: texttitle_
      procedure, pass(self) :: textvtitle_
      procedure, pass(self) :: textint_
      procedure, pass(self) :: textvar_
      procedure, pass(self) :: textvardbl_
      procedure, pass(self) :: textballoon_
      procedure, pass(self) :: balloon_
      procedure, pass(self) :: setfont_
      procedure, pass(self) :: Printhkl_
      procedure, pass(self) :: DumpIndices_ 
      procedure, pass(self) :: PrintIndices_
      procedure, pass(self) :: DumpImage_
      procedure, pass(self) :: DumpImageDistort_
      procedure, pass(self) :: DrawSPFrame_
      procedure, pass(self) :: DrawcellFrame_
      procedure, pass(self) :: getpsscale_
      procedure, pass(self) :: StereoProj_
      procedure, pass(self) :: DumpZAP_
      procedure, pass(self) :: DumpPP_


      generic, public :: openFile => openFile_ 
!DEC$ ATTRIBUTES DLLEXPORT :: openFile
      generic, public :: closeFile => closeFile_ 
!DEC$ ATTRIBUTES DLLEXPORT :: closeFile
      generic, public :: newpage => newpage_
!DEC$ ATTRIBUTES DLLEXPORT :: newpage
      generic, public :: cellinfo => cellinfo_
!DEC$ ATTRIBUTES DLLEXPORT :: cellinfo
      generic, public :: clippath => clippath_
!DEC$ ATTRIBUTES DLLEXPORT :: clippath
      generic, public :: translate => translate_
!DEC$ ATTRIBUTES DLLEXPORT :: translate
      generic, public :: move => move_
!DEC$ ATTRIBUTES DLLEXPORT :: move
      generic, public :: draw => draw_
!DEC$ ATTRIBUTES DLLEXPORT :: draw
      generic, public :: line_gray => line_gray_
!DEC$ ATTRIBUTES DLLEXPORT :: line_gray
      generic, public :: setlinewidth => setlinewidth_
!DEC$ ATTRIBUTES DLLEXPORT :: setlinewidth
      generic, public :: square => square_
!DEC$ ATTRIBUTES DLLEXPORT :: square
      generic, public :: filledsquare => filledsquare_
!DEC$ ATTRIBUTES DLLEXPORT :: filledsquare
      generic, public :: cross => cross_
!DEC$ ATTRIBUTES DLLEXPORT :: cross
      generic, public :: sphere => sphere_
!DEC$ ATTRIBUTES DLLEXPORT :: sphere
      generic, public :: arc => arc_
!DEC$ ATTRIBUTES DLLEXPORT :: arc 
      generic, public :: circle => circle_
!DEC$ ATTRIBUTES DLLEXPORT :: circle
      generic, public :: filledcircle => filledcircle_
!DEC$ ATTRIBUTES DLLEXPORT :: filledcircle
      generic, public :: drawframe => drawframe_
!DEC$ ATTRIBUTES DLLEXPORT :: drawframe
      generic, public :: drawrect => drawrect_
!DEC$ ATTRIBUTES DLLEXPORT :: drawrect
      generic, public :: line => line_
!DEC$ ATTRIBUTES DLLEXPORT :: line
      generic, public :: setdash => setdash_
!DEC$ ATTRIBUTES DLLEXPORT :: setdash
      generic, public :: closepathS => closepathS_
!DEC$ ATTRIBUTES DLLEXPORT :: closepathS
      generic, public :: stroke => stroke_
!DEC$ ATTRIBUTES DLLEXPORT :: stroke
      generic, public :: gsave => gsave_
!DEC$ ATTRIBUTES DLLEXPORT :: gsave
      generic, public :: grestore => grestore_
!DEC$ ATTRIBUTES DLLEXPORT :: grestore
      generic, public :: closepath => closepath_
!DEC$ ATTRIBUTES DLLEXPORT :: closepath
      generic, public :: newpath => newpath_
!DEC$ ATTRIBUTES DLLEXPORT :: newpath
      generic, public :: text => text_
!DEC$ ATTRIBUTES DLLEXPORT :: text
      generic, public :: textv => textv_
!DEC$ ATTRIBUTES DLLEXPORT :: textv
      generic, public :: texttitle => texttitle_
!DEC$ ATTRIBUTES DLLEXPORT :: texttitle
      generic, public :: textvtitle => textvtitle_
!DEC$ ATTRIBUTES DLLEXPORT :: textvtitle
      generic, public :: textint => textint_
!DEC$ ATTRIBUTES DLLEXPORT :: textint
      generic, public :: textvar => textvar_
!DEC$ ATTRIBUTES DLLEXPORT :: textvar
      generic, public :: textvardbl => textvardbl_
!DEC$ ATTRIBUTES DLLEXPORT :: textvardbl
      generic, public :: textballoon => textballoon_
!DEC$ ATTRIBUTES DLLEXPORT :: textballoon
      generic, public :: balloon => balloon_
!DEC$ ATTRIBUTES DLLEXPORT :: balloon
      generic, public :: setfont => setfont_
!DEC$ ATTRIBUTES DLLEXPORT :: setfont
      generic, public :: Printhkl => Printhkl_
!DEC$ ATTRIBUTES DLLEXPORT :: Printhkl
      generic, public :: DumpIndices => DumpIndices_
!DEC$ ATTRIBUTES DLLEXPORT :: DumpIndices
      generic, public :: PrintIndices => PrintIndices_
!DEC$ ATTRIBUTES DLLEXPORT :: PrintIndices
      generic, public :: DumpImage => DumpImage_
!DEC$ ATTRIBUTES DLLEXPORT :: DumpImage
      generic, public :: DumpImageDistort => DumpImageDistort_
!DEC$ ATTRIBUTES DLLEXPORT :: DumpImageDistort
      generic, public :: DrawSPFrame => DrawSPFrame_
!DEC$ ATTRIBUTES DLLEXPORT :: DrawSPFrame
      generic, public :: DrawcellFrame => DrawcellFrame_
!DEC$ ATTRIBUTES DLLEXPORT :: DrawcellFrame
      generic, public :: getpsscale => getpsscale_
!DEC$ ATTRIBUTES DLLEXPORT :: getpsscale
      generic, public :: StereoProj => StereoProj_
!DEC$ ATTRIBUTES DLLEXPORT :: StereoProj
      generic, public :: DumpZAP => DumpZAP_ 
!DEC$ ATTRIBUTES DLLEXPORT :: DumpZAP
      generic, public :: DumpPP => DumpPP_
!DEC$ ATTRIBUTES DLLEXPORT :: DumpPP

  end type PostScript_T

  ! the constructor routine for this class 
  interface PostScript_T
    module procedure :: PS_constructor
  end interface PostScript_T

contains 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! we begin with the functions/subroutines that are public in this class
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(PostScript_T) function PS_constructor( progdesc, EMsoft, imanum ) result(PS)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! constructor for the PostScript_T Class 

use mod_EMsoft 

IMPLICIT NONE

character(fnlen), INTENT(IN)    :: progdesc 
type(EMsoft_T), INTENT(INOUT)   :: EMsoft
integer(kind=irg), INTENT(IN)   :: imanum

  PS%imanum = imanum
  call PS%openfile(progdesc, EMsoft)

end function PS_constructor

!--------------------------------------------------------------------------
recursive subroutine openfile_(self, progdesc, EMsoft, dontask)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! open postscript file and dump the preamble to the file

use mod_EMsoft
use mod_io

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)     :: self
character(fnlen),INTENT(IN)           :: progdesc
type(EMsoft_T),INTENT(INOUT)          :: EMsoft
logical,INTENT(IN),optional           :: dontask		
 !! optional parameter to select file opening route

type(IO_T)                            :: Message 

real(kind=sgl)    		                :: fw, fh		!< page format parameters
integer(kind=irg) 		                :: i			!< loop counter
character(fnlen)                      :: gname

! define the writeable portion of the page (should be made more user-friendly by adding A4 format...)
 self%psfigwidth=6.5 
 self%psfigheight=9.0
 self%psscale=1.0

! open file and dump Prolog and Comments sections
 if (present(dontask)) then
! don't ask for a file name 
   open(unit=self%psunit,file=trim(EMsoft%toNativePath(self%psname)),status='unknown',action='write',form='formatted')
   call Message%printMessage('Opening temporary file for PostScript output', frm = "(A)")
 else
! do ask for a file name 
   call Message%ReadValue(' Enter Postscript file name : ', gname,"(A)")
   self%psname = trim(gname)
   open(unit=self%psunit,file=trim(EMsoft%toNativePath(self%psname)),status='unknown',form='formatted')
 end if

! write the preamble
 write (self%psunit,"(A)") PSpreamble(1)
 write (self%psunit,"(A,' ',A)") trim(PSpreamble(2)), trim(EMsoft%getConfigParameter('UserName'))
 write (self%psunit,"(A,' ',A)") trim(PSpreamble(3)), trim(progdesc)
 do i=4,23
  write (self%psunit,"(A)") PSpreamble(i)
 end do 

! determine lower left corner and translate to that point
 fw=0.5*(8.50-self%psscale*self%psfigwidth)
 fh=0.5*(11.0-self%psscale*self%psfigheight)
 write (self%psunit,"(F12.7,' ',F12.7,' T')") fw,fh
 write (self%psunit,"(F12.7,' setlinewidth')") 0.01
 write (self%psunit,"(F12.7,' ',F12.7,' scale')") self%psscale,self%psscale

! set page number counter to zero
 self%pspage = 0

end subroutine openfile_

!--------------------------------------------------------------------------
recursive subroutine closefile_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! close and save postscript file

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)    :: self 

! write the trailer to the file 
 write (self%psunit,*) 'showpage'
 write (self%psunit,"(' %%Pages: ',i3)") self%pspage
 write (self%psunit,"(' %%EOF')")
 
! and close it
  close(unit=self%psunit,status='keep')

end subroutine closefile_

!--------------------------------------------------------------------------
recursive subroutine newpage_(self, frm, btxt)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! start a new page in the PS file

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self
logical,INTENT(IN)        	        :: frm		
 !! logical draw frame or not
character(*),INTENT(IN)   	        :: btxt		
 !! character string for header balloon

 if (self%pspage.ne.0) then
  write (self%psunit,*) 'showpage saveobj restore'
 end if

! update the page counter
 self%pspage = self%pspage + 1
 write (self%psunit,"(' %%Page: ',i3,i3)") self%pspage-1,self%pspage
 write (self%psunit,*) '/saveobj save def'
 
! prepare to draw a header balloon 
 call self%setfont(PSfonts(3),0.18)
 write (self%psunit,"(1x,F12.7,' ',F12.7,' M (',I8,') show')") 6.75,self%psfigheight-0.2,self%pspage
 if (frm.eqv..TRUE.) then  ! we need a frame
  call self%drawframe(6.75,self%psfigheight)
 endif
 
! output the text balloon
 call self%setlinewidth(0.012)
 call self%textballoon(2.0,9.2,btxt,PSfonts(2),0.25)
 call self%setfont(PSfonts(5),0.07)
 call self%text(0.1,-0.1,PSlbl)
 
end subroutine newpage_

!--------------------------------------------------------------------------
recursive subroutine cellinfo_(self, cell, xo, yo)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! write unit cell information (for drawing programs)

use mod_crystallography 

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)        :: self
type(Cell_T),INTENT(INOUT)               :: cell
real(kind=sgl),INTENT(IN)                :: xo,yo		
 !! starting location for output 

 call self%setfont(PSfonts(2),0.12/self%psscale)
 call self%text(xo,yo,'Filename: '//trim(cell%getFileName()))
 call self%text(xo,yo-0.16,'a [nm]          ')
 call self%text(xo,yo-0.30,'b [nm]          ')
 call self%text(xo,yo-0.44,'c [nm]          ')
 call self%text(xo,yo-0.58,'alpha [deg]     ')
 call self%text(xo,yo-0.72,'beta  [deg]     ')
 call self%text(xo,yo-0.86,'gamma [deg]     ')
 call self%textvardbl(xo+0.75,yo-0.16,': ',cell%getLatParm('a'))
 call self%textvardbl(xo+0.75,yo-0.30,': ',cell%getLatParm('b'))
 call self%textvardbl(xo+0.75,yo-0.44,': ',cell%getLatParm('c'))
 call self%textvardbl(xo+0.75,yo-0.58,': ',cell%getLatParm('alpha'))
 call self%textvardbl(xo+0.75,yo-0.72,': ',cell%getLatParm('beta'))
 call self%textvardbl(xo+0.75,yo-0.86,': ',cell%getLatParm('gamma'))

end subroutine cellinfo_


!--------------------------------------------------------------------------
recursive function getpsscale_(self) result(psscale)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! make the last path the clippath

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)  :: self
real(kind=sgl)                     :: psscale

psscale = self%psscale

end function getpsscale_

!--------------------------------------------------------------------------
recursive subroutine clippath_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! make the last path the clippath

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)        :: self

 write (self%psunit,"('Cl clip')")

end subroutine clippath_

!--------------------------------------------------------------------------
recursive subroutine translate_(self,x,y)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! redefine the origin of the current coordinate frame 

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)     :: self
real(kind=sgl),INTENT(IN) 	          :: x,y	
 !! coordinates of new origin

 write (self%psunit,"(F18.7,' ',F18.7,' T')") x,y

end subroutine translate_

!--------------------------------------------------------------------------
recursive subroutine move_(self,x,y)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! move to a given location

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self
real(kind=sgl),INTENT(IN)           :: x,y	
 !! move to this location

 write (self%psunit,"(F18.7,' ',F18.7,' M')") x,y

end subroutine move_

!--------------------------------------------------------------------------
recursive subroutine draw_(self,x,y)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a line from the current point to the new point

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)         :: x,y	
 !! end coordinates of draw

  write (self%psunit,"(F18.7,' ',F18.7,' L')") x,y

end subroutine draw_

!--------------------------------------------------------------------------
recursive subroutine line_gray_(self,x1,y1,x2,y2,gray)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a line with a given gray level from the current point to the new point        

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x1,y1		
 !! starting point
real(kind=sgl),INTENT(IN)  	      :: x2,y2		
 !! end point
real(kind=sgl),INTENT(IN)  	      :: gray		  
 !! gray level

  write (self%psunit,"(F18.7,' setgray ')") gray  
  call self%move(x1,y1)
  call self%draw(x2,y2)

! and reset the gray level to black
  write (self%psunit,"('S  0.0 setgray ')")      

end subroutine line_gray_

!--------------------------------------------------------------------------
recursive subroutine setlinewidth_(self,x)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !!  set the line width

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x		
 !! line width parameter

 write (self%psunit,"(F12.7,' setlinewidth')") x

end subroutine setlinewidth_

!--------------------------------------------------------------------------
recursive subroutine square_(self,x,y,edge)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a square 

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)       	:: x,y	
 !! center coordinates
real(kind=sgl),INTENT(IN)       	:: edge	
 !! edge length

real(kind=sgl)  		              :: ed	

 ed=0.5*edge
 write (self%psunit,"('0.0 setgray')")
 write (self%psunit,"('newpath')")
 write (self%psunit,"(2(F12.7,' '),'moveto')") x-ed,y-ed
 write (self%psunit,"(2(F12.7,' '),'lineto')") x-ed,y+ed
 write (self%psunit,"(2(F12.7,' '),'lineto')") x+ed,y+ed
 write (self%psunit,"(2(F12.7,' '),'lineto closepath S')") x+ed,y-ed

end subroutine square_

!--------------------------------------------------------------------------
recursive subroutine filledsquare_(self,x,y,edge,graylevel)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a filled square       

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x,y		    
 !! center coordinates
real(kind=sgl),INTENT(IN)  	      :: edge		    
 !! edge length
real(kind=sgl),INTENT(IN)  	      :: graylevel	
 !! gray level for filling

real(kind=sgl)  		:: ed

 ed=0.5*edge
 write (self%psunit,"(F12.7,' setgray')") graylevel
 call self%newpath()
 call self%move(x-ed,y-ed)
 call self%draw(x-ed,y+ed)
 call self%draw(x+ed,y+ed)
 write (self%psunit,"(2(F12.7,' '),'lineto closepath fill S')") x+ed,y-ed

end subroutine filledsquare_

!--------------------------------------------------------------------------
recursive subroutine cross_(self,x,y,edge,lw)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a small cross       

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x,y	
 !! center coordinates
real(kind=sgl),INTENT(IN)  	      :: edge	
 !! edge length
real(kind=sgl),INTENT(IN)  	      :: lw		
 !! line width

real(kind=sgl)  		              :: ed	

 ed=0.5*edge
 call self%setlinewidth(lw)
 call self%newpath()
 call self%move(x-ed,y-ed)
 call self%draw(x+ed,y+ed)
 call self%move(x-ed,y+ed)
 call self%draw(x+ed,y-ed)
 call self%stroke()

end subroutine cross_

!--------------------------------------------------------------------------
recursive subroutine sphere_(self,x,y,r,clr)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a colored sphere
  !!
  !! method modified from Earl J. Kirkland''s book, page 226, adapted for 
  !! color PostScript

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)         :: x,y    
 !! center coordinates
real(kind=sgl),INTENT(IN)         :: r      
 !! radius
integer(kind=irg),INTENT(IN)      :: clr    
 !! atomic number

write (self%psunit,"(1x,7(f12.5,1x),'sp')") ATOM_colors(1:3,clr),r,r,x,y

end subroutine sphere_

!--------------------------------------------------------------------------
recursive subroutine arc_(self,x0,y0,x,y,radius,ang1,ang2)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw an arc of a circle (see PostScript 'arc' command for details)

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x0,y0		   
 !! new origin coordinates
real(kind=sgl),INTENT(IN)  	      :: x,y			   
 !! center coordinates
real(kind=sgl),INTENT(IN)  	      :: radius		   
 !! radius
real(kind=sgl),INTENT(IN)  	      :: ang1,ang2	 
 !! start and end angles


 write (self%psunit,"('N ',2(F16.10,' '),' moveto ',5(E16.8,' '),' arc S')") x0,y0,x,y,radius,ang1,ang2

end subroutine arc_

!--------------------------------------------------------------------------
recursive subroutine circle_(self,x,y,radius)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a circle

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x,y			
  !! center coordinates
real(kind=sgl),INTENT(IN)  	      :: radius		
  !! radius

 write (self%psunit,"('N ',3(F16.10,' '),'0 360 arc Cl S')") x,y,radius

end subroutine circle_

!--------------------------------------------------------------------------
recursive subroutine filledcircle_(self,x,y,radius,graylevel)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a filled circle

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x,y		    
 !! center coordinates
real(kind=sgl),INTENT(IN)  	      :: radius	    
 !! radius
real(kind=sgl),INTENT(IN)  	      :: graylevel	
 !! gray level

 write (self%psunit,"(F12.7,' setgray')") graylevel
 write (self%psunit,"('N ',3(F12.7,' '),'0 360 arc Cl F')") x,y,radius

end subroutine filledcircle_

!--------------------------------------------------------------------------
recursive subroutine drawframe_(self,x,y)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw the main frame

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x		
 !! frame width
real(kind=sgl),INTENT(IN)  	      :: y 		
 !! frame height

call self%drawrect(0.0,0.0,x,y)

end subroutine drawframe_

!--------------------------------------------------------------------------
recursive subroutine drawrect_(self,x1,y1,x2,y2)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a rectangle

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x1, y1		
 !! lower left
real(kind=sgl),INTENT(IN)  	      :: x2, y2 	
 !! upper right

 write (self%psunit,"('N')") 
 call self%move(x1,y1)
 call self%draw(x1,y2)
 call self%draw(x2,y2)
 call self%draw(x2,y1)
 call self%draw(x1,y1)
 call self%closepathS()
 write (self%psunit,"('[0.15 0.03 0.02 0.03] 0 setdash')") 
 write (self%psunit,"('[] 0 setdash')") 

end subroutine drawrect_

!--------------------------------------------------------------------------
recursive subroutine line_(self,x1,y1,x2,y2)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a line between two points

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)  	      :: x1, y1		
 !! first point
real(kind=sgl),INTENT(IN)  	      :: x2, y2 	
 !! second point

  call self%move(x1,y1)
  call self%draw(x2,y2)
  write (self%psunit,"('S')")      

end subroutine line_

!--------------------------------------------------------------------------
recursive subroutine setdash_(self, num)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! define a dash pattern
  !!
  !! @details Note that the dash pattern must be defined in the calling program.

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)  	  :: num	
 !! dash pattern number of components/segments

integer(kind=irg)  		            :: i	! loop counter

 write (self%psunit,"('[')")
 do i=1,num
   write (self%psunit,"(F12.7,' ')") self%psdash(i)
 end do
 write (self%psunit,"('] ',I4,' setdash')") int(self%psdash(num+1))

end subroutine setdash_

!--------------------------------------------------------------------------
recursive subroutine closepathS_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! close current path and Stroke

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self

write (self%psunit,"('Cl S')")

end subroutine closepathS_

!--------------------------------------------------------------------------
recursive subroutine stroke_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! stroke the current path

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self

write (self%psunit,"('S ')")

end subroutine stroke_

!--------------------------------------------------------------------------
recursive subroutine gsave_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! save the current graphics settings

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self

write (self%psunit,"('gsave ')")

end subroutine gsave_

!--------------------------------------------------------------------------
recursive subroutine grestore_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! restore the previous graphics settings

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self

write (self%psunit,"('grestore ')")

end subroutine grestore_

!--------------------------------------------------------------------------
recursive subroutine closepath_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! close the current path

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self

write (self%psunit,"('Cl ')")

end subroutine closepath_

!--------------------------------------------------------------------------
recursive subroutine newpath_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! start a new path

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self

write (self%psunit,"('newpath ')")

end subroutine newpath_

!--------------------------------------------------------------------------
recursive subroutine text_(self,x,y,line)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw text at a given location

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y	
 !! text start coordinates
character(*),INTENT(IN)	          :: line	
 !! output string

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(') show')") 

end subroutine text_

!--------------------------------------------------------------------------
recursive subroutine textv_(self,x,y,line)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw text rotated counterclockwise by 90 degrees

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y	
 !! text start coordinates
character(*),INTENT(IN)	          :: line	
 !! output string

 write (self%psunit,"('gsave ')") 
 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('90.0 rotate')") 
 write (self%psunit,"('( ',A,' ) show')") line
 write (self%psunit,"('-90.0 rotate grestore')") 

end subroutine textv_

!--------------------------------------------------------------------------
recursive subroutine texttitle_(self,x,y,line,q)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw the title

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y	
 !! text start coordinates
character(*),INTENT(IN)	          :: line	
 !! output string
real(kind=sgl),INTENT(IN)	        :: q	  
 !! 

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"('  [x',1PE8.0,'] ) show')") q

end subroutine texttitle_

!--------------------------------------------------------------------------
recursive subroutine textvtitle_(self,x,y,line,q)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw a vertical title

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)         :: x,y  
 !! text start coordinates
character(*),INTENT(IN)           :: line 
 !! output string
real(kind=sgl),INTENT(IN)         :: q    
 !! 

 write (self%psunit,"('gsave ')") 
 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('90.0 rotate')") 
 write (self%psunit,"('(')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"('  [x',1PE8.0,'] ) show')") q
 write (self%psunit,"('-90.0 rotate grestore')") 

end subroutine textvtitle_

!--------------------------------------------------------------------------
recursive subroutine textint_(self,x,y,line,vl)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! text followed by an integer number

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y		  
 !! text start coordinates
character(*),INTENT(IN)	          :: line	      
 !! output string
integer(kind=irg),INTENT(IN)	    :: vl		
 !! integer output value

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(I4)",advance="no") vl
 write (self%psunit,"(') show')") 

end subroutine textint_

!--------------------------------------------------------------------------
recursive subroutine textvar_(self,x,y,line,vl)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! text followed by a real number

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y	
 !! text start coordinates
character(*),INTENT(IN)	          :: line	
 !! output string
real(kind=sgl),INTENT(IN)	        :: vl		
 !! real output value

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(F14.4)",advance="no") vl
 write (self%psunit,"(') show')") 

end subroutine textvar_

!--------------------------------------------------------------------------
recursive subroutine textvardbl_(self,x,y,line,vl)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! text followed by a double precision real number

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y	
 !! text start coordinates
character(*),INTENT(IN)	          :: line	
 !! output string
real(kind=dbl),INTENT(IN)	        :: vl		
 !! double output value

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(F12.6)",advance="no") vl
 write (self%psunit,"(') show')") 

end subroutine textvardbl_

!--------------------------------------------------------------------------
recursive subroutine textballoon_(self,x,y,line,font,sc)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! text inside a rounded balloon

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y	
 !! text start coordinates
character(*),INTENT(IN)	          :: line	
 !! output string
character(*),INTENT(IN)	          :: font	
 !! font string
real(kind=sgl),INTENT(IN)	        :: sc	  
 !! scale factor

 call self%setfont(font,sc)
 write (self%psunit,"('/length (')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(') stringwidth pop def')") 
 write (self%psunit,"('/height ',F6.4,' def /border ',F6.4,' def')") 0.11*sc/0.2,0.06*sc/0.2
 write (self%psunit,"('/bottom ',F12.7,' border sub def')") y
 write (self%psunit,"('/top ',F12.7,' height add border add def')") y
 write (self%psunit,"('/left ',F12.7,' border sub def')") x
 write (self%psunit,"('/right ',F12.7,' length add border add def')") x
 write (self%psunit,"('/rad 0.04 def frame')") 
 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no") 
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(') show')") 
 
end subroutine textballoon_

!--------------------------------------------------------------------------
recursive subroutine balloon_(self,x,y,le,he,w)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw an empty balloon

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: x,y		
 !! start coordinates
real(kind=sgl),INTENT(IN)	        :: le, he	
 !! length and height
real(kind=sgl),INTENT(IN)	        :: w		  
 !! width parameter

 write (self%psunit,"('/he ',F6.4,' def /bo ',F6.4,' def /wi ',F6.4,' def')") he,0.5*w,le
 write (self%psunit,"('/bottom ',F12.7,' bo add def')") y
 write (self%psunit,"('/top bottom he add bo sub bo sub def')")  
 write (self%psunit,"('/left ',F12.7,' bo add def')") x
 write (self%psunit,"('/right left wi add bo sub def')")  
 write (self%psunit,"('/rad 0.04 def frame')") 

end subroutine balloon_

!--------------------------------------------------------------------------
recursive subroutine setfont_(self,line,sc)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! select a font and make it active

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)	        :: sc	  
 !! font scale factor  
character(*),INTENT(IN)	          :: line	
 !! font string

 write (self%psunit,"()",advance="no") 
 write (self%psunit,"('/',A)",advance="no") line
 write (self%psunit,"(' findfont')") 
 write (self%psunit,"(F6.4,' scalefont ')") sc
 write (self%psunit,"('setfont')")

end subroutine setfont_

!--------------------------------------------------------------------------
recursive subroutine Printhkl_(self,x,y,h,k,l)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! print hkl indices in PostScript format

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)      :: h,k,l		!< Miller index triplet
real(kind=sgl),INTENT(IN)		      :: x,y			!< starting position of indices
character(1),parameter 		        :: numbers(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
real(kind=sgl)         		        :: xo,yo,dx,dy,x1,y1
character(1)           		        :: line

 call self%setfont(PSfonts(5),0.065)
 call self%setlinewidth(0.004)
 xo = 0.050
 yo = 0.050
 dx = 0.050
 dy = 0.065
! THIS ONLY WORKS FOR INDICES -9 <= ... <= 9  !!!
! first index
 x1=x+xo
 y1=y+yo
 line=numbers(abs(h))
 call self%text(x1,y1,line)
 if (h.lt.0) then
  call self%line(x1,y1+dy,x1+0.5*dx,y1+dy)
 end if
! second index
 x1=x1+dx
 line=numbers(abs(k))
 call self%text(x1,y1,line)
 if (k.lt.0) then
  call self%line(x1,y1+dy,x1+0.5*dx,y1+dy)
 end if
! third index
 x1=x1+dx
 line=numbers(abs(l))
 call self%text(x1,y1,line)
 if (l.lt.0) then
  call self%line(x1,y1+dy,x1+0.5*dx,y1+dy)
 end if

end subroutine Printhkl_

!--------------------------------------------------------------------------
recursive subroutine DumpIndices_(self,hexset,S,h,k,l,c,x,y,n)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! version of Printhkl used by stereographic projection program

use mod_misc

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)       :: self
!f2py intent(in,out) ::  PS
logical,INTENT(IN)                      :: hexset
character(1),INTENT(IN)           	    :: S			
 !! space character 'd' or 'r'
integer(kind=irg),INTENT(IN)      	    :: h,k,l	
 !! Miller index triplet
integer(kind=irg),INTENT(IN)		        :: c			
 !! positioning parameter
real(kind=sgl),INTENT(IN)		            :: x,y		
 !! starting position of indices
logical,INTENT(IN)		                  :: n			
 !! logical

character(1),parameter :: numbers(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
character(1)           :: line
integer(kind=irg)      :: uvw(3),uvtw(4)
real(kind=sgl)         :: xo,yo,dx,dy,x1,y1

 call self%setfont(PSfonts(5),0.08/self%psscale)
 xo = 0.050/self%psscale
 yo = 0.050/self%psscale
 dx = 0.050/self%psscale
 dy = 0.075/self%psscale

 if (n.eqv..FALSE.) then
  xo = -0.30/self%psscale
  yo = -0.10/self%psscale
 endif

 if (c.eq.2) then
  xo = -0.30/self%psscale
  yo = 0.05/self%psscale
 end if

 uvw =(/ h,k,l /)

 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  call MilBrav(uvw,uvtw,'34')
  call IndexReduceMB(uvtw)
 end if

! opening bracket
 if (S.eq.'d') then
  call self%text(x+xo,y+yo,'[')
 else
  call self%text(x+xo,y+yo,'\(')
 end if

!first index
 x1=x+xo+dx
 y1=y+yo
 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  line=numbers(abs(uvtw(1)))
 else
  line=numbers(abs(uvw(1)))
 end if
 call self%text(x1,y1,line)
 if (h.lt.0) then
  call self%line(x1,y1+dy,x1+0.8*dx,y1+dy)
 end if

! second index
 x1=x1+dx
 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  line=numbers(abs(uvtw(2)))
 else
  line=numbers(abs(uvw(2)))
 end if
 call self%text(x1,y1,line)
 if (k.lt.0) then
  call self%line(x1,y1+dy,x1+0.8*dx,y1+dy)
 end if

! third index (if hexset = .TRUE.) -> put a period
 if (hexset.eqv..TRUE.) then
  x1=x1+dx
  line='.'
  call self%text(x1,y1,line)
 end if

! last index
 x1=x1+dx
 if ((S.eq.'d').AND.(hexset.eqv..TRUE.)) then
  line=numbers(abs(uvtw(4)))
 else
  line=numbers(abs(uvw(3)))
 end if
 call self%text(x1,y1,line)
 if (l.lt.0) then
  call self%line(x1,y1+dy,x1+0.8*dx,y1+dy)
 end if

! closing bracket
 x1=x1+dx
 if (S.eq.'d') then
  call self%text(x1,y+yo,']')
 else
  call self%text(x1,y+yo,'\)')
 end if

 dx=dx*0.3
 if (c.eq.2) then
  call self%line(x+xo,y1-0.02/self%psscale,x1+dx,y1-0.02/self%psscale)
 end if

end subroutine DumpIndices_

!--------------------------------------------------------------------------
recursive subroutine PrintIndices_(self,S,hexset,h,k,l,x,y)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! version of Printhkl used by stereographic projection program

use mod_misc

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self
character(1),INTENT(IN)           	:: S			  
 !! space character 'd' or 'r'
logical,INTENT(IN)                  :: hexset
 !! hexagonal setting logical
integer(kind=irg),INTENT(IN)      	:: h,k,l		
 !! Miller index triplet
real(kind=sgl),INTENT(IN)		        :: x,y			
 !! starting position of indices

character(12)    			              :: line
integer(kind=irg)			              :: hkl(3)

 hkl = (/ h,k,l /)
 call IndexString(hexset,line,hkl,S)
 call self%text(x,y,line)

end subroutine PrintIndices_

!--------------------------------------------------------------------------
recursive subroutine DumpImage_(self,imaint,imanum,x0,y0,npx,npy,scl)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw integer image (512x512 maximum size) at given location with given scale

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN) 	        :: imaint(npx,npy)
 !! image array
integer(kind=irg),INTENT(INOUT)       :: imanum
 !! image number 
real(kind=sgl),INTENT(IN)             :: x0,y0		
 !! image position
integer(kind=irg),INTENT(IN)         	:: npx,npy	
 !! image size
real(kind=sgl),INTENT(IN)             :: scl			
 !! image scale factor

 call self%DumpImageDistort(imaint,imanum,x0,y0,npx,npy,scl,scl)
 
end subroutine DumpImage_

!--------------------------------------------------------------------------
recursive subroutine DumpImageDistort_(self,imaint,imanum,x0,y0,npx,npy,sclx,scly)   
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! draw integer image (512x512 maximum size) at given location with given scale which may be different along x and y

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN) 	        :: imaint(npx,npy)
integer(kind=irg),INTENT(INOUT)       :: imanum
real(kind=sgl),INTENT(IN)             :: x0,y0		  
 !! image position
integer(kind=irg),INTENT(IN)         	:: npx,npy		
 !! image size
real(kind=sgl),INTENT(IN)             :: sclx,scly	
 !! image scale factors

integer(kind=irg)                 	  :: iq,i,j,ir,iq1,iq2,k
integer(kind=irg),parameter       	  :: bpp=8
character(2*npx)                  	  :: bigone
character(3),parameter            	  :: imnm(20) = (/'i01','i02','i03','i04','i05','i06', &
                                                      'i07','i08','i09','i10','i11','i12','i13','i14', &
                                                      'i15','i16','i17','i18','i19','i20'/)
character(1),parameter            	  :: hd(0:15) = (/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'/)

 imanum = imanum + 1

! integer image
 call self%translate(x0,y0)
 write (self%psunit,"('/picstr ',i3,' string def')") npx
 write (self%psunit,"(' ')")
 write (self%psunit,"('/',3A)") imnm(imanum)
 write (self%psunit,"('{ ',i3,' ',i3,' ',i1,' [',i3,' 0 0 ',i3,' 0 0]')") npx,npy,bpp,npx,npy
 write (self%psunit,"(' { currentfile picstr readhexstring pop } image } def')")
 write (self%psunit,"(' gsave ',f7.4,' ',f7.4,' scale ')") sclx,scly*float(npy)/float(npx)
 write (self%psunit,"(3A)") imnm(imanum)
 
 do j=1,npy
  do i=1,npx
   ir=2*i-1
   iq=imaint(i,j)
   iq1=iq/16
   iq2=mod(iq,16)
   bigone(ir:ir)=hd(iq1)
   ir=ir+1
   bigone(ir:ir)=hd(iq2)
  end do
  k=2*npx
  write (self%psunit,"(A)") bigone
 end do

 write (self%psunit,"('grestore')")
 call self%translate(-x0,-y0)

end subroutine DumpImageDistort_



!--------------------------------------------------------------------------
recursive subroutine DrawSPFrame_(self,cell,SG,CX,CY,CRad,iview,sp)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !!  draw a stereographic projection layout 

use mod_crystallography 
use mod_symmetry
use mod_misc, only: IndexString

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)        :: cell
type(SpaceGroup_T),INTENT(INOUT)  :: SG
real(kind=sgl),INTENT(IN)    			:: CX, CY		
 !! circle center coordinates 
real(kind=sgl),INTENT(IN)    			:: CRad 		
 !! circle radius
integer(kind=irg),INTENT(INOUT)		:: iview(3)	
 !! zone axis indices
character(1),INTENT(IN)	      		:: sp			  
 !! drawing space

character(12)     :: instr
character(17)     :: str

 call self%newpage(.FALSE.,'Stereographic Projection')

 call self%setlinewidth(0.012)
 call self%circle(CX,CY,CRad)
 call self%setlinewidth(0.008)
 call self%line(CX-CRad,CY,CX+CRad,CY)
 call self%line(CX,CY-CRad,CX,CY+CRad)
 call self%text(CX-CRad-0.07,CY-0.025,'A')
 call self%text(CX+CRad+0.03,CY-0.025,'B')
 call self%text(CX-0.03,CY-CRad-0.08,'M''')
 call self%text(CX-0.03,CY+CRad+0.05,'M"')
 call self%cellinfo(cell,0.00,8.30)
 call IndexString(SG%getSpaceGrouphexset(),instr,iview,'d')
 call self%text(CX,8.14,'Viewing Direction '//instr)

 if (sp.eq.'d') then 
  str='direct space'
 else
  str='reciprocal space'
 end if

 call self%text(CX,8.00,'Projection of '//str)

end subroutine DrawSPFrame_

!--------------------------------------------------------------------------
subroutine DrawcellFrame_(self, cell, iview, sp, CX, CY, hexset)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! format the output page for the EMdrawcell program

use mod_io 
use mod_crystallography 
use mod_misc 

IMPLICIT NONE

class(PostScript_T), INTENT(INOUT)  :: self
type(Cell_T), INTENT(INOUT)         :: cell 
integer(kind=irg),INTENT(INOUT)     :: iview(3)             
 !! viewing direction
character(1),INTENT(IN)             :: sp                   
 !! drawing space character
real(kind=sgl), INTENT(IN)          :: CX, CY               
 !! center of page   
logical, INTENT(IN)                 :: hexset

character(12)                       :: instr
character(17)                       :: str

 if (sp.eq.'d') then 
  call self%newpage(.FALSE.,'Crystal Structure Drawing')
 else
  call self%newpage(.FALSE.,'Reciprocal Lattice Drawing')
 endif
 call self%setlinewidth(0.012)
 call self%drawrect(0.0,0.0,CX,CY)
 call self%setlinewidth(0.008)
 call self%setfont(PSfonts(2),0.12/self % psscale)
 call self%cellinfo(cell,0.0,8.3)
 call IndexString(hexset, instr, iview, 'd')
 call self%text(CX*0.5,8.14,'Viewing Direction '//instr)
 if (sp.eq.'d') then 
  str='direct space'
 else
  str='reciprocal space'
 endif
 call self%text(CX*0.5,8.00,'Drawing of '//str)

end subroutine DrawcellFrame_

!--------------------------------------------------------------------------
recursive subroutine StereoProj_(self, cell, SG, sp, iview, hm, km, lm, topbot)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! draw a stereographic projection

use mod_crystallography
use mod_symmetry 
use mod_math 
use mod_io 
use mod_misc 

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self 
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
character(1),INTENT(IN)             :: sp         
 !! space character 'd' or 'r'
integer(kind=irg),INTENT(INOUT)     :: iview(3)   
 !! viewing direction
integer(kind=irg),INTENT(IN)        :: hm, km, lm 
 !! maximum h,k,l indices to be included in drawing
logical,INTENT(IN)                  :: topbot     
 !! logical 

logical                             :: nn
logical,allocatable                 :: z(:,:,:)
real(kind=sgl)                      :: rr(3),g(3),r(3),M(3,3), CX, CY, CRad,xst,yst 
real(kind=sgl),parameter            :: negthresh = -0.000001
integer(kind=irg)                   :: i,h,k,l,hkl(3),hkil(4),cr,hh,kk,ll,num,hkm 
integer(kind=irg),allocatable       :: itmp(:,:) 

! 20cm radius projection circle [inches]
 CRad = 3.937
 CX = 3.25
 CY = 3.5
 
! create transformation matrix
 call ProjectionMatrix(cell, iview, M)
 
! write text and draw projection circle
 call self%DrawSPFrame(cell, SG, CX, CY, CRad, iview, sp)

! loop over families
! make sure that the arrays are big enough for the hexagonal case...
if (SG%getSpaceGrouphexset().eqv..TRUE.) then 
  hkm = abs(hm)+abs(km)
  allocate(z(-hkm:hkm,-hkm:hkm,-lm:lm))
else
 allocate(z(-hm:hm,-km:km,-lm:lm))
end if

z = .FALSE.
 do hh=-hm,hm
  do kk=-km,km
   do ll=-lm,lm
    if (z(hh,kk,ll).eqv..TRUE.) cycle
! determine the family members
    if ((SG%getSpaceGrouphexset().eqv..TRUE.).AND.(sp.eq.'d')) then
     hkil= (/ hh,kk,-(hh+kk),ll /)
     call MilBrav(hkl,hkil,'43')
    else
     hkl= (/ hh,kk,ll /)
    end if
    if ((hh**2+kk**2+ll**2).ne.0) then
     call IndexReduce(hkl)
     call SG%CalcFamily(hkl, num, sp, itmp)
! loop over all points and draw projection+label
      do i=0,num-1
       h=itmp(i,1)
       k=itmp(i,2)
       l=itmp(i,3)
       hkl(1:3)=itmp(i,1:3)
       if ((h.le.hm).and.(k.le.km).and.(l.le.lm)) then
! reduce to smallest integers to avoid overlap
! of indices, such as (111) and (222)
         call IndexReduce(hkl)
         h=hkl(1)
         k=hkl(2)
         l=hkl(3)
         g=float( (/h,k,l/) )
         call cell%TransSpace(g,r,sp,'c')
         call cell%NormVec(r,'c')
! apply viewing tansformation
         rr = matmul(M,r)
! compute stereographic projection coordinates
         xst=CX+CRad*rr(1)/(1.0+abs(rr(3)))
         yst=CY+CRad*rr(2)/(1.0+abs(rr(3)))
         cr=1
         if (z(h,k,l).eqv..FALSE.) then
          if (rr(3).gt.negthresh) then
           call self%filledcircle(xst,yst,0.015/self%getpsscale(),0.0)
           nn = .TRUE.
           call self%DumpIndices(SG%getSpaceGrouphexset(),sp,h,k,l,cr,xst,yst,nn)
          else if (topbot) then
           call self%circle(xst,yst,0.035/self%getpsscale())
           nn = .FALSE.
           call self%DumpIndices(SG%getSpaceGrouphexset(),sp,h,k,l,cr,xst,yst,nn)
          end if
         end if
         z(h,k,l) = .TRUE.
       end if
      end do
    end if
   end do 
  end do 
 end do 

end subroutine StereoProj_

!--------------------------------------------------------------------------
recursive subroutine DumpZAP_(self,cell,SG,xo,yo,u,v,w,p,np,first,indi,laL,icnt,dbdiff,Vg,Vgsave,rg,rfamily,rnumfam,hhcc)
  !! author: MDG
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! draw a single zone axis diffraction pattern

use mod_io
use mod_crystallography
use mod_symmetry

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self 
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
real(kind=sgl),INTENT(IN)           :: xo, yo               
 !! lower left position
integer(kind=irg),INTENT(IN)        :: u, v, w              
 !! zone axis components
integer(kind=irg),INTENT(IN)        :: p                    
 !! ??
logical,INTENT(IN)                  :: np                   
 !! logical for new page
logical,INTENT(IN)                  :: first                
 !! logical
integer(kind=irg),INTENT(IN)        :: indi                 
 !! ??
real(kind=sgl),INTENT(IN)           :: laL                  
 !! camera length
integer(kind=irg),INTENT(IN)        :: icnt                 
 !! counter
logical,INTENT(IN)                  :: dbdiff(icnt)         
 !! array to deal with double diffraction spots
real(kind=sgl),INTENT(IN)           :: Vg(*),rg(*),Vgsave(*)
integer(kind=irg),INTENT(IN)        :: hhcc
integer(kind=irg),INTENT(IN)        :: rfamily(hhcc,48,3),rnumfam(*)

type(IO_T)                          :: Message
! nref is the anticipated maximum number of reflections per pattern
integer(kind=irg),parameter         :: nref = 2000
integer(kind=irg)                   :: dp,i,j,jcnt,ui,vi,wi,pp,locg(nref,3),ier, io_int(1)
integer(kind=irg),allocatable       :: idx(:)
real(kind=sgl)                      :: sc,gmax,leng(nref),PX,PY,qx,qy,locv(nref),locvsave(nref),t(3),c(3),gg(3),gx(3),gy(3)
real(kind=sgl),allocatable          :: lng(:)
real(kind=sgl),parameter            :: le=3.25,he=2.9375,thr=1.0E-4
logical                             :: dbd(nref)

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
 if (np) then
  call self%newpage(.FALSE.,'Kinematical Zone Axis Patterns')
  call self%text(5.25,-0.05,'scale bar in reciprocal nm')
  gmax = laL
  call self%textvar(5.25,self % psfigheight+0.02,'Camera Constant [nm mm]',gmax)
  call self%setfont(PSfonts(2),0.15)
  call self%text(-0.25,self % psfigheight+0.02,'Structure File : '//cell %getFileName())
 end if

! draw frame and related stuff
 call self%setlinewidth(0.012)
 call self%balloon(xo,yo,le,he,0.0312)
 call self%setlinewidth(0.001)
 PX = xo+1.8125
 PY = yo+1.0+15.0/32.0
 call self%circle(PX,PY,1.375)

! zone axis
 call self%setfont(PSfonts(2),0.12)
 call self%text(xo+0.05,yo+he-0.15,'Zone axis ')
 ui=u
 vi=v
 wi=w
 call self%PrintIndices('d',SG%getSpaceGrouphexset(),ui,vi,wi,xo+0.6,yo+he-0.15)

! multiplicity
 call self%setfont(PSfonts(2),0.12)
 pp=p
 call self%textint(xo+0.05,yo+he-0.30,'Multiplicity ',pp)

! scale bar (sc is the conversion factor from nm-1 to inches)
 sc = laL/25.4
 call self%setlinewidth(0.020)
 call self%line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
 call self%setfont(PSfonts(2),0.15)
 call self%text(xo+0.05+2.5*sc,yo+0.10,'5 ')

! select all reflections belonging to this zone axis pattern
 leng(1:nref)=0.0
 jcnt=0
 do i=1,icnt
  do j=1,rnumfam(i)
   dp=u*rfamily(i,j,1)+v*rfamily(i,j,2)+w*rfamily(i,j,3)
   if (dp.eq.0) then
    jcnt=jcnt+1
    if (jcnt.gt.nref) call Message%printError('DumpZAP ',' too many reflections (<2000)' )
    locg(jcnt,1:3)=rfamily(i,j,1:3)
    leng(jcnt)=rg(i)
    locv(jcnt)=Vg(i)
    locvsave(jcnt)=Vgsave(i)
    dbd(jcnt)=.FALSE.
! take care of potential double diffraction reflections
    if ((SG%getSpaceGroupSymmorphic()).and.(dbdiff(i))) dbd(jcnt) = .TRUE.
   end if
  end do
 end do

! rank them by length (use SLATEC routine)
 allocate(idx(jcnt))
 allocate(lng(jcnt))
 lng(1:jcnt) = leng(1:jcnt)
 call SPSORT(lng,jcnt,idx,1,ier)
 io_int(1) = jcnt
 call Message%WriteValue(' Number of reflections : ', io_int, 1)

! normalize the zone axis in cartesian components; this is the z-axis
 t(1)=-float(u)
 t(2)=-float(v)
 t(3)=-float(w)
 call cell%TransSpace(t,c,'d','c')
 call cell%NormVec(c,'c')

! take the first reflection in the list and make that the x-axis
! skip the zero reflection !!
 j=idx(2)

! normalize the first reciprocal vector in cartesian components
! this will be the x-axis of the diffraction pattern
 gg(1:3)=float(locg(j,1:3))
 call cell%TransSpace(gg,gx,'r','c')
 call cell%NormVec(gx,'c')

! then get the cross product between t and g; this is the y-axis
 call cell%CalcCross(c,gx,gy,'c','c',0)

! plot origin of reciprocal space 
 call self%filledcircle(PX,PY,0.05,0.0)

! then plot the remaining reflections
 do i=1,jcnt
  j=idx(i)
  gg(1:3)=float(locg(j,1:3))
  call cell%TransSpace(gg,c,'r','c')
  qx=PX-cell%CalcDot(c,gx,'c')*sc
  qy=PY+cell%CalcDot(c,gy,'c')*sc

! first check for systematic absence due to lattice centering
  if ((locvsave(j).eq.-100.0).and.(indi.ge.2)) then
    call self%cross(qx,qy,0.03,0.001)
  end if

! could it be a double diffraction spot ?
  if ((SG%getSpaceGroupSymmorphic()).and.(dbd(j))) call self%square(qx,qy,0.04)

! is it a regular reflection ?
  if (locv(j).ge.thr) then
   call self%filledcircle(qx,qy,locv(j),0.0)
   if ((indi.eq.1).or.(indi.eq.3)) then 
    call self%Printhkl(qx,qy,locg(j,1),locg(j,2),locg(j,3))
   end if
  end if

 end do

 deallocate(idx, lng)
 
end subroutine DumpZAP_

!--------------------------------------------------------------------------
recursive subroutine DumpPP_(self,cell,xo,yo,np,laL,icnt,Vgsave,rg,rnumfam)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! draw a kinematical powder pattern

use mod_crystallography 

IMPLICIT NONE 

class(PostScript_T),INTENT(INOUT)   :: self 
type(Cell_T),INTENT(INOUT)          :: cell
real(kind=sgl),INTENT(IN)           :: xo, yo               
 !! lower left position   
logical,INTENT(IN)                  :: np                   
 !! logical for new page    
real(kind=sgl),INTENT(IN)           :: laL                  
 !! camera length   
integer(kind=irg),INTENT(IN)        :: icnt                 
 !! counter   
real(kind=sgl),INTENT(IN)           :: rg(*),Vgsave(*)
integer(kind=irg),INTENT(IN)        :: rnumfam(*)
    
! nref = max number of rings    
integer(kind=irg),parameter         :: nref = 500
integer(kind=irg)                   :: i,j
real(kind=sgl)                      :: sc,gmax,leng(nref),PX,PY,locv(nref),grad,w,Vmax
real(kind=sgl),parameter            :: le=3.25,he=2.9375,thr=1.0E-4

! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
 if (np) then
  call self%newpage(.FALSE.,'Kinematical Zone Axis Patterns')
  call self%text(5.25,-0.05,'scale bar in reciprocal nm')
  gmax = laL
  call self%textvar(5.25,self%psfigheight+0.02,'Camera Constant [nm mm]',gmax)
  call self%setfont(PSfonts(2),0.15)
  call self%text(-0.25,self%psfigheight+0.02,'Structure File : '//cell%getFileName())
 end if

! draw frame and related stuff
 call self%setlinewidth(0.012)
 call self%balloon(xo,yo,le,he,0.0312)
 call self%setlinewidth(0.001)
 PX = xo+1.8125
 PY = yo+1.0+15.0/32.0
 call self%circle(PX,PY,1.375)

! frame title
 call self%setfont(PSfonts(2),0.12)
 call self%text(xo+0.05,yo+he-0.15,'Powder Pattern')

! scale bar (sc is the conversion factor from nm-1 to inches)
 sc = laL/25.4
 call self%setlinewidth(0.020)
 call self%line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
 call self%setfont(PSfonts(2),0.15)
 call self%text(xo+0.05+2.5*sc,yo+0.10,'5 ')

! scale all reflection intensities by the multiplicity
 leng(1:nref)=0.0
 Vmax = 0.0
 do i=1,icnt-1
  leng(i)=rg(i)
  locv(i)=Vgsave(i)*rnumfam(i)
  if (locv(i).gt.Vmax) Vmax=locv(i)
 end do

! plot origin of reciprocal space 
 call self%filledcircle(PX,PY,0.05,0.0)

! then plot the diffraction circles 
 do i=1,icnt
  j=icnt+1-i
! get the circle radius and intensity
  grad = leng(j)*sc
  w = locv(j)*0.03/Vmax
! draw circle if radius fits in drawing frame and intensity large enough
  if ((w.gt.0.0001).AND.(grad.le.1.375)) then 
   call self%setlinewidth(w)
   call self%circle(PX,PY,grad)
  end if
 end do

end subroutine DumpPP_


end module mod_postscript
