! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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
      procedure, pass(self) :: setpspage_
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
      procedure, pass(self) :: getpsunit_
      procedure, pass(self) :: getpsscale_
      procedure, pass(self) :: getpsfigwidth_
      procedure, pass(self) :: getpsfigheight_
      procedure, pass(self) :: StereoProj_
      procedure, pass(self) :: DumpZAP_
      procedure, pass(self) :: DumpPP_
      procedure, pass(self) :: DumpLine_
      procedure, pass(self) :: DumpMatrix_
      procedure, pass(self) :: DumpAtom_
      procedure, pass(self) :: InfoPage_
      procedure, pass(self) :: StrucFacPage_
      procedure, pass(self) :: StereoPage_
      procedure, pass(self) :: DiffPage_
      final :: PS_destructor

      generic, public :: openFile => openFile_
      generic, public :: closeFile => closeFile_
      generic, public :: newpage => newpage_
      generic, public :: setpspage => setpspage_
      generic, public :: cellinfo => cellinfo_
      generic, public :: clippath => clippath_
      generic, public :: translate => translate_
      generic, public :: move => move_
      generic, public :: draw => draw_
      generic, public :: line_gray => line_gray_
      generic, public :: setlinewidth => setlinewidth_
      generic, public :: square => square_
      generic, public :: filledsquare => filledsquare_
      generic, public :: cross => cross_
      generic, public :: sphere => sphere_
      generic, public :: arc => arc_
      generic, public :: circle => circle_
      generic, public :: filledcircle => filledcircle_
      generic, public :: drawframe => drawframe_
      generic, public :: drawrect => drawrect_
      generic, public :: line => line_
      generic, public :: setdash => setdash_
      generic, public :: closepathS => closepathS_
      generic, public :: stroke => stroke_
      generic, public :: gsave => gsave_
      generic, public :: grestore => grestore_
      generic, public :: closepath => closepath_
      generic, public :: newpath => newpath_
      generic, public :: text => text_
      generic, public :: textv => textv_
      generic, public :: texttitle => texttitle_
      generic, public :: textvtitle => textvtitle_
      generic, public :: textint => textint_
      generic, public :: textvar => textvar_
      generic, public :: textvardbl => textvardbl_
      generic, public :: textballoon => textballoon_
      generic, public :: balloon => balloon_
      generic, public :: setfont => setfont_
      generic, public :: Printhkl => Printhkl_

      generic, public :: DumpIndices => DumpIndices_
      generic, public :: PrintIndices => PrintIndices_
      generic, public :: DumpImage => DumpImage_
      generic, public :: DumpImageDistort => DumpImageDistort_
      generic, public :: DrawSPFrame => DrawSPFrame_
      generic, public :: DrawcellFrame => DrawcellFrame_
      generic, public :: getpsunit => getpsunit_
      generic, public :: getpsscale => getpsscale_
      generic, public :: getpsfigwidth => getpsfigwidth_
      generic, public :: getpsfigheight => getpsfigheight_
      generic, public :: StereoProj => StereoProj_
      generic, public :: DumpZAP => DumpZAP_
      generic, public :: DumpPP => DumpPP_
      generic, public :: dumpline => dumpline_
      generic, public :: DumpMatrix => DumpMatrix_
      generic, public :: DumpAtom => DumpAtom_
      generic, public :: InfoPage => InfoPage_
      generic, public :: StrucFacPage => StrucFacPage_
      generic, public :: StereoPage => StereoPage_
      generic, public :: DiffPage => DiffPage_


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
!DEC$ ATTRIBUTES DLLEXPORT :: PS_constructor
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
subroutine PS_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: PS_destructor
!! author: MDG
!! version: 1.0
!! date: 02/02/20
!!
!! destructor for the PostScript_T Class

use mod_io

IMPLICIT NONE

type(PostScript_T), INTENT(INOUT)     :: self

type(IO_T)                            :: Message
logical                               :: itsopen

call reportDestructor('PostScript_T')

! if the MRC uit is still open, close it here
inquire(unit=self%psunit, opened=itsopen)

if (itsopen.eqv..TRUE.) then
  close(unit=self%psunit, status='keep')
  call Message%printMessage(' Closed PostScript file '//trim(self%psname))
end if

end subroutine PS_destructor

!--------------------------------------------------------------------------
recursive subroutine openfile_(self, progdesc, EMsoft, dontask)
!DEC$ ATTRIBUTES DLLEXPORT :: openFile_
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

real(kind=sgl)                            :: fw, fh        !< page format parameters
integer(kind=irg)                         :: i            !< loop counter
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
recursive subroutine closeFile_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: closeFile_
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

end subroutine closeFile_

!--------------------------------------------------------------------------
recursive subroutine newpage_(self, frm, btxt)
!DEC$ ATTRIBUTES DLLEXPORT :: newpage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! start a new page in the PS file

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self
logical,INTENT(IN)                    :: frm
 !! logical draw frame or not
character(*),INTENT(IN)               :: btxt
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
recursive subroutine setpspage_(self, pgnum)
!DEC$ ATTRIBUTES DLLEXPORT :: setpspage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! set a page number

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)       :: pgnum

self%pspage = pgnum

end subroutine setpspage_

!--------------------------------------------------------------------------
recursive subroutine cellinfo_(self, cell, xo, yo)
!DEC$ ATTRIBUTES DLLEXPORT :: cellinfo_
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
recursive function getpsunit_(self) result(psu)
!DEC$ ATTRIBUTES DLLEXPORT :: getpsunit_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! get the postscript output device id

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)  :: self
integer(kind=irg)                  :: psu

psu = self%psunit

end function getpsunit_

!--------------------------------------------------------------------------
recursive function getpsscale_(self) result(psscale)
!DEC$ ATTRIBUTES DLLEXPORT :: getpsscale_
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
recursive function getpsfigwidth_(self) result(psfw)
!DEC$ ATTRIBUTES DLLEXPORT :: getpsfigwidth_
  !! author: MDG
  !! version: 1.0
  !! date: 01/28/20
  !!
  !! return ps figure width

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)  :: self
real(kind=sgl)                     :: psfw

psfw = self%psfigwidth

end function getpsfigwidth_

!--------------------------------------------------------------------------
recursive function getpsfigheight_(self) result(psfh)
!DEC$ ATTRIBUTES DLLEXPORT :: getpsfigheight_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! make the last path the clippath

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)  :: self
real(kind=sgl)                     :: psfh

psfh = self%psfigheight

end function getpsfigheight_


!--------------------------------------------------------------------------
recursive subroutine clippath_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: clippath_
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
!DEC$ ATTRIBUTES DLLEXPORT :: translate_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! redefine the origin of the current coordinate frame

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)     :: self
real(kind=sgl),INTENT(IN)               :: x,y
 !! coordinates of new origin

 write (self%psunit,"(F18.7,' ',F18.7,' T')") x,y

end subroutine translate_

!--------------------------------------------------------------------------
recursive subroutine move_(self,x,y)
!DEC$ ATTRIBUTES DLLEXPORT :: move_
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
!DEC$ ATTRIBUTES DLLEXPORT :: draw_
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
!DEC$ ATTRIBUTES DLLEXPORT :: line_gray_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a line with a given gray level from the current point to the new point

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x1,y1
 !! starting point
real(kind=sgl),INTENT(IN)            :: x2,y2
 !! end point
real(kind=sgl),INTENT(IN)            :: gray
 !! gray level

  write (self%psunit,"(F18.7,' setgray ')") gray
  call self%move(x1,y1)
  call self%draw(x2,y2)

! and reset the gray level to black
  write (self%psunit,"('S  0.0 setgray ')")

end subroutine line_gray_

!--------------------------------------------------------------------------
recursive subroutine setlinewidth_(self,x)
!DEC$ ATTRIBUTES DLLEXPORT :: setlinewidth_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !!  set the line width

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x
 !! line width parameter

 write (self%psunit,"(F12.7,' setlinewidth')") x

end subroutine setlinewidth_

!--------------------------------------------------------------------------
recursive subroutine square_(self,x,y,edge)
!DEC$ ATTRIBUTES DLLEXPORT :: square_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a square

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)           :: x,y
 !! center coordinates
real(kind=sgl),INTENT(IN)           :: edge
 !! edge length

real(kind=sgl)                        :: ed

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
!DEC$ ATTRIBUTES DLLEXPORT :: filledsquare_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a filled square

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! center coordinates
real(kind=sgl),INTENT(IN)            :: edge
 !! edge length
real(kind=sgl),INTENT(IN)            :: graylevel
 !! gray level for filling

real(kind=sgl)          :: ed

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
!DEC$ ATTRIBUTES DLLEXPORT :: cross_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a small cross

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! center coordinates
real(kind=sgl),INTENT(IN)            :: edge
 !! edge length
real(kind=sgl),INTENT(IN)            :: lw
 !! line width

real(kind=sgl)                        :: ed

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
!DEC$ ATTRIBUTES DLLEXPORT :: sphere_
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
!DEC$ ATTRIBUTES DLLEXPORT :: arc_

  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw an arc of a circle (see PostScript 'arc' command for details)

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x0,y0
 !! new origin coordinates
real(kind=sgl),INTENT(IN)            :: x,y
 !! center coordinates
real(kind=sgl),INTENT(IN)            :: radius
 !! radius
real(kind=sgl),INTENT(IN)            :: ang1,ang2
 !! start and end angles


 write (self%psunit,"('N ',2(F16.10,' '),' moveto ',5(E16.8,' '),' arc S')") x0,y0,x,y,radius,ang1,ang2

end subroutine arc_

!--------------------------------------------------------------------------
recursive subroutine circle_(self,x,y,radius)
!DEC$ ATTRIBUTES DLLEXPORT :: circle_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a circle

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
  !! center coordinates
real(kind=sgl),INTENT(IN)            :: radius
  !! radius

 write (self%psunit,"('N ',3(F16.10,' '),'0 360 arc Cl S')") x,y,radius

end subroutine circle_

!--------------------------------------------------------------------------
recursive subroutine filledcircle_(self,x,y,radius,graylevel)
!DEC$ ATTRIBUTES DLLEXPORT :: filledcircle_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a filled circle

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! center coordinates
real(kind=sgl),INTENT(IN)            :: radius
 !! radius
real(kind=sgl),INTENT(IN)            :: graylevel
 !! gray level

 write (self%psunit,"(F12.7,' setgray')") graylevel
 write (self%psunit,"('N ',3(F12.7,' '),'0 360 arc Cl F')") x,y,radius

end subroutine filledcircle_

!--------------------------------------------------------------------------
recursive subroutine drawframe_(self,x,y)
!DEC$ ATTRIBUTES DLLEXPORT :: drawframe_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw the main frame

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x
 !! frame width
real(kind=sgl),INTENT(IN)            :: y
 !! frame height

call self%drawrect(0.0,0.0,x,y)

end subroutine drawframe_

!--------------------------------------------------------------------------
recursive subroutine drawrect_(self,x1,y1,x2,y2)
!DEC$ ATTRIBUTES DLLEXPORT :: drawrect_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a rectangle

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x1, y1
 !! lower left
real(kind=sgl),INTENT(IN)            :: x2, y2
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
!DEC$ ATTRIBUTES DLLEXPORT :: line_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw a line between two points

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x1, y1
 !! first point
real(kind=sgl),INTENT(IN)            :: x2, y2
 !! second point

  call self%move(x1,y1)
  call self%draw(x2,y2)
  write (self%psunit,"('S')")

end subroutine line_

!--------------------------------------------------------------------------
recursive subroutine setdash_(self, num)
!DEC$ ATTRIBUTES DLLEXPORT :: setdash_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! define a dash pattern
  !!
  !! @details Note that the dash pattern must be defined in the calling program.

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)        :: num
 !! dash pattern number of components/segments

integer(kind=irg)                      :: i    ! loop counter

 write (self%psunit,"('[')")
 do i=1,num
   write (self%psunit,"(F12.7,' ')") self%psdash(i)
 end do
 write (self%psunit,"('] ',I4,' setdash')") int(self%psdash(num+1))

end subroutine setdash_

!--------------------------------------------------------------------------
recursive subroutine closepathS_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: closepathS_
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
!DEC$ ATTRIBUTES DLLEXPORT :: stroke_

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
!DEC$ ATTRIBUTES DLLEXPORT :: gsave_
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
!DEC$ ATTRIBUTES DLLEXPORT :: grestore_
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
!DEC$ ATTRIBUTES DLLEXPORT :: closepath_
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
!DEC$ ATTRIBUTES DLLEXPORT :: newpath_

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
!DEC$ ATTRIBUTES DLLEXPORT :: text_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw text at a given location

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! text start coordinates
character(*),INTENT(IN)              :: line
 !! output string

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no")
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(') show')")

end subroutine text_

!--------------------------------------------------------------------------
recursive subroutine textv_(self,x,y,line)
!DEC$ ATTRIBUTES DLLEXPORT :: textv_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw text rotated counterclockwise by 90 degrees

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! text start coordinates
character(*),INTENT(IN)              :: line
 !! output string

 write (self%psunit,"('gsave ')")
 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('90.0 rotate')")
 write (self%psunit,"('( ',A,' ) show')") line
 write (self%psunit,"('-90.0 rotate grestore')")

end subroutine textv_

!--------------------------------------------------------------------------
recursive subroutine texttitle_(self,x,y,line,q)
!DEC$ ATTRIBUTES DLLEXPORT :: texttitle_

  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw the title

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! text start coordinates
character(*),INTENT(IN)              :: line
 !! output string
real(kind=sgl),INTENT(IN)            :: q
 !!

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no")
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"('  [x',1PE8.0,'] ) show')") q

end subroutine texttitle_

!--------------------------------------------------------------------------
recursive subroutine textvtitle_(self,x,y,line,q)
!DEC$ ATTRIBUTES DLLEXPORT :: textvtitle_
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
!DEC$ ATTRIBUTES DLLEXPORT :: textint_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! text followed by an integer number

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! text start coordinates
character(*),INTENT(IN)              :: line
 !! output string
integer(kind=irg),INTENT(IN)        :: vl
 !! integer output value

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no")
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(I4)",advance="no") vl
 write (self%psunit,"(') show')")

end subroutine textint_

!--------------------------------------------------------------------------
recursive subroutine textvar_(self,x,y,line,vl)
!DEC$ ATTRIBUTES DLLEXPORT :: textvar_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! text followed by a real number

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! text start coordinates
character(*),INTENT(IN)              :: line
 !! output string
real(kind=sgl),INTENT(IN)            :: vl
 !! real output value

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no")
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(F14.4)",advance="no") vl
 write (self%psunit,"(') show')")

end subroutine textvar_

!--------------------------------------------------------------------------
recursive subroutine textvardbl_(self,x,y,line,vl)
!DEC$ ATTRIBUTES DLLEXPORT :: textvardbl_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! text followed by a double precision real number

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! text start coordinates
character(*),INTENT(IN)              :: line
 !! output string
real(kind=dbl),INTENT(IN)            :: vl
 !! double output value

 write (self%psunit,"(F12.7,' ',F12.7,' M')") x,y
 write (self%psunit,"('(')",advance="no")
 write (self%psunit,"(A)",advance="no") line
 write (self%psunit,"(F12.6)",advance="no") vl
 write (self%psunit,"(') show')")

end subroutine textvardbl_

!--------------------------------------------------------------------------
recursive subroutine textballoon_(self,x,y,line,font,sc)
!DEC$ ATTRIBUTES DLLEXPORT :: textballoon_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! text inside a rounded balloon

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! text start coordinates
character(*),INTENT(IN)              :: line
 !! output string
character(*),INTENT(IN)              :: font
 !! font string
real(kind=sgl),INTENT(IN)            :: sc
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
!DEC$ ATTRIBUTES DLLEXPORT :: balloon_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw an empty balloon

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: x,y
 !! start coordinates
real(kind=sgl),INTENT(IN)            :: le, he
 !! length and height
real(kind=sgl),INTENT(IN)            :: w
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
!DEC$ ATTRIBUTES DLLEXPORT :: setfont_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! select a font and make it active

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
real(kind=sgl),INTENT(IN)            :: sc
 !! font scale factor
character(*),INTENT(IN)              :: line
 !! font string

 write (self%psunit,"()",advance="no")
 write (self%psunit,"('/',A)",advance="no") line
 write (self%psunit,"(' findfont')")
 write (self%psunit,"(F6.4,' scalefont ')") sc
 write (self%psunit,"('setfont')")

end subroutine setfont_

!--------------------------------------------------------------------------
recursive subroutine Printhkl_(self,x,y,h,k,l)
!DEC$ ATTRIBUTES DLLEXPORT :: Printhkl_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! print hkl indices in PostScript format

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)      :: h,k,l        !< Miller index triplet
real(kind=sgl),INTENT(IN)              :: x,y            !< starting position of indices
character(1),parameter                 :: numbers(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
real(kind=sgl)                         :: xo,yo,dx,dy,x1,y1
character(1)                           :: line

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
!DEC$ ATTRIBUTES DLLEXPORT :: DumpIndices_
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
character(1),INTENT(IN)                   :: S
 !! space character 'd' or 'r'
integer(kind=irg),INTENT(IN)              :: h,k,l
 !! Miller index triplet
integer(kind=irg),INTENT(IN)                :: c
 !! positioning parameter
real(kind=sgl),INTENT(IN)                    :: x,y
 !! starting position of indices
logical,INTENT(IN)                          :: n
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
!DEC$ ATTRIBUTES DLLEXPORT :: PrintIndices_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! version of Printhkl used by stereographic projection program

use mod_misc

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self
character(1),INTENT(IN)               :: S
 !! space character 'd' or 'r'
logical,INTENT(IN)                  :: hexset
 !! hexagonal setting logical
integer(kind=irg),INTENT(IN)          :: h,k,l
 !! Miller index triplet
real(kind=sgl),INTENT(IN)                :: x,y
 !! starting position of indices

character(12)                              :: line
integer(kind=irg)                          :: hkl(3)

 hkl = (/ h,k,l /)
 call IndexString(hexset,line,hkl,S)
 call self%text(x,y,line)

end subroutine PrintIndices_

!--------------------------------------------------------------------------
recursive subroutine DumpImage_(self,imaint,imanum,x0,y0,npx,npy,scl)
!DEC$ ATTRIBUTES DLLEXPORT :: DumpImage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw integer image (512x512 maximum size) at given location with given scale

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)             :: imaint(npx,npy)
 !! image array
integer(kind=irg),INTENT(INOUT)       :: imanum
 !! image number
real(kind=sgl),INTENT(IN)             :: x0,y0
 !! image position
integer(kind=irg),INTENT(IN)             :: npx,npy
 !! image size
real(kind=sgl),INTENT(IN)             :: scl
 !! image scale factor

 call self%DumpImageDistort(imaint,imanum,x0,y0,npx,npy,scl,scl)

end subroutine DumpImage_

!--------------------------------------------------------------------------
recursive subroutine DumpImageDistort_(self,imaint,imanum,x0,y0,npx,npy,sclx,scly)
!DEC$ ATTRIBUTES DLLEXPORT :: DumpImageDistort_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! draw integer image (512x512 maximum size) at given location with given scale which may be different along x and y

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)     :: self
integer(kind=irg),INTENT(IN)             :: imaint(npx,npy)
integer(kind=irg),INTENT(INOUT)       :: imanum
real(kind=sgl),INTENT(IN)             :: x0,y0
 !! image position
integer(kind=irg),INTENT(IN)             :: npx,npy
 !! image size
real(kind=sgl),INTENT(IN)             :: sclx,scly
 !! image scale factors

integer(kind=irg)                       :: iq,i,j,ir,iq1,iq2,k
integer(kind=irg),parameter             :: bpp=8
character(2*npx)                        :: bigone
character(3),parameter                  :: imnm(20) = (/'i01','i02','i03','i04','i05','i06', &
                                                      'i07','i08','i09','i10','i11','i12','i13','i14', &
                                                      'i15','i16','i17','i18','i19','i20'/)
character(1),parameter                  :: hd(0:15) = (/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'/)

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
!DEC$ ATTRIBUTES DLLEXPORT :: DrawSPFrame_
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
real(kind=sgl),INTENT(IN)                :: CX, CY
 !! circle center coordinates
real(kind=sgl),INTENT(IN)                :: CRad
 !! circle radius
integer(kind=irg),INTENT(INOUT)        :: iview(3)
 !! zone axis indices
character(1),INTENT(IN)                  :: sp
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
!DEC$ ATTRIBUTES DLLEXPORT :: DrawcellFrame_
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
!DEC$ ATTRIBUTES DLLEXPORT :: StereoProj_
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
!DEC$ ATTRIBUTES DLLEXPORT :: DumpZAP_
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
!DEC$ ATTRIBUTES DLLEXPORT :: DumpPP_
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

!--------------------------------------------------------------------------
subroutine InfoPage_(self, cell, SG)
!dec$ attributes dllexport :: InfoPage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! creates a page with crystallographic unit cell info

use mod_io
use mod_crystallography
use mod_symmetry

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT) :: self
type(Cell_T),INTENT(INOUT)        :: cell
type(SpaceGroup_T),INTENT(INOUT)  :: SG

real(kind=sgl)                    :: ast,bst,cst,alphast,betast,gammast,coor(5),dx,dy,db,topx,topy,x,y
integer(kind=irg)                 :: iap,i,j
character(1),parameter            :: xsys(7) = (/'c','t','o','h','R','m','a'/)
character(2)                      :: brvs
character(11)                     :: sgname
real(kind=dbl)                    :: dmt(3,3), rmt(3,3), dsm(3,3), rsm(3,3)
real(kind=dbl),allocatable        :: apos(:,:)
integer(kind=irg),allocatable     :: atp(:)

! start first page
 call self%newpage(.TRUE.,trim(cell%getFileName()))

! write text
 call self%text(5.25,0.05,'Distances [nm], angles [degrees]')
 call self%setlinewidth(0.012)

! direct space information
 call self%textballoon(0.75,8.45,'Direct Space',PSfonts(2),0.20)
 call self%setfont(PSfonts(4),0.14)
 call self%text(1.0,8.00,'a :')
 call self%text(1.0,7.84,'b :')
 call self%text(1.0,7.68,'c :')
 write (self%psunit,900) 1.25,8.00,cell%getLatParm('a')
 write (self%psunit,900) 1.25,7.84,cell%getLatParm('b')
 write (self%psunit,900) 1.25,7.68,cell%getLatParm('c')
 call self%setfont(PSfonts(1),0.14)
 call self%text(1.0,7.52,'a :')
 call self%text(1.0,7.36,'b :')
 call self%text(1.0,7.20,'g :')
 call self%setfont(PSfonts(4),0.14)
 write (self%psunit,900) 1.25,7.52,cell%getLatParm('alpha')
 write (self%psunit,900) 1.25,7.36,cell%getLatParm('beta')
 write (self%psunit,900) 1.25,7.20,cell%getLatParm('gamma')

! space group information and unit cell volume
 call self%text(3.0,8.00,'Space Group :')
 call self%setfont(PSfonts(2),0.14)
 write (self%psunit,905) 4.2,8.00,SG%getSpaceGroupName(), SG%getSpaceGroupNumber()
 call self%setfont(PSfonts(4),0.14)
 call self%text(3.0,7.84,'Bravais Lattice :')
 call self%setfont(PSfonts(2),0.14)
 brvs(1:1)=xsys(SG%getSpaceGroupXtalSystem())
 sgname = SG%getSpaceGroupName()
 if (SG%getSpaceGroupXtalSystem().ne.5) then
  brvs(2:2)=sgname(2:2)
 else
  brvs(2:2)=' '
 endif
 write (self%psunit,902) 4.2,7.84,brvs
 call self%setfont(PSfonts(4),0.14)
 call self%text(3.0,7.52,'Volume V [nm^3] :')
 write (self%psunit,900) 4.2,7.52,cell%getVolume()

 dsm = cell%getdsm()
 rsm = cell%getrsm()
 dmt = cell%getdmt()
 rmt = cell%getrmt()

! reciprocal space information
 call self%textballoon(0.75,6.75,'Reciprocal Space',PSfonts(2),0.20)
 ast=sqrt(rmt(1,1))
 bst=sqrt(rmt(2,2))
 cst=sqrt(rmt(3,3))
 alphast=acos(rmt(2,3)/bst/cst)/dtor
 betast =acos(rmt(1,3)/ast/cst)/dtor
 gammast=acos(rmt(1,2)/ast/bst)/dtor
 call self%setfont(PSfonts(4),0.14)
 call self%text(1.0,6.30,'a* :')
 call self%text(1.0,6.14,'b* :')
 call self%text(1.0,5.98,'c* :')
 write (self%psunit,900) 1.25,6.30,ast
 write (self%psunit,900) 1.25,6.14,bst
 write (self%psunit,900) 1.25,5.98,cst
 call self%setfont(PSfonts(1),0.14)
 call self%text(1.0,5.82,'a* :')
 call self%text(1.0,5.66,'b* :')
 call self%text(1.0,5.50,'g* :')
 call self%setfont(PSfonts(4),0.14)
 write (self%psunit,900) 1.25,5.82,alphast
 write (self%psunit,900) 1.25,5.66,betast
 write (self%psunit,900) 1.25,5.50,gammast

! unit cell volume and reciprocal basis vectors
 call self%text(3.0,6.30,'Volume V* [nm^-3] :')
 write (self%psunit,900) 4.2,6.30,1.0/cell%getVolume()

! metric tensors and structure matrices
 call self%textballoon(0.75,5.05,'Important Matrices',PSfonts(2),0.20)
 call self%DumpMatrix(1.0,4.2,sngl(dmt),'g =')
 call self%DumpMatrix(3.5,4.2,sngl(rmt),'g*=')
 call self%DumpMatrix(1.0,3.2,sngl(dsm),'a =')
 call self%DumpMatrix(3.5,3.2,sngl(rsm),'b =')

apos = cell%getAsymPosData()
atp = cell%getAtomtype()

! asymmetric unit
 call self%textballoon(0.75,2.80,'Asymmetric Unit',PSfonts(2),0.20)
 dx= 2.75
 dy=-0.16
 topx =-1.75
 topy = 2.3
 x=topx
 y=topy
 db = 0.2
 do i=1,cell%getNatomtype()
  if (mod(i-1,10).eq.0) then
   x=x+dx
   call self%setfont(PSfonts(3),0.14)
   call self%text(x,y+0.2,'atom ')
   call self%text(x+0.55,y+0.2,'x ')
   call self%text(x+0.90,y+0.2,'y ')
   call self%text(x+1.25,y+0.2,'z ')
   call self%text(x+1.50,y+0.2,'occ. ')
   call self%text(x+1.95,y+0.2,'B ')
   call self%setfont(PSfonts(4),0.12)
  endif
  y=topy + mod(i-1,10)*dy
  do j=1,3
   coor(j) = apos(i,j)
  end do
  coor(4) = apos(i,4)
  coor(5) = apos(i,5)
  iap = atp(i)
  call self%DumpAtom(x,y,ATOM_sym(iap),iap,coor)
 end do

! real number output
 900    format (1x,F12.7,' ',F12.7,' M (',F12.4,') show')
! integer number output
 901    format (1x,F12.7,' ',F12.7,' M (',I8,') show')
! character output
 902    format (1x,F12.7,' ',F12.7,' M (',A,') show')
 905    format (1x,F12.7,' ',F12.7,' M (',A11,'  [#',I3,']) show')

end subroutine InfoPage_

!--------------------------------------------------------------------------
subroutine DumpMatrix_(self, x, y, g, tt)
!dec$ attributes dllexport :: DumpMatrix_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! format a 3x3 matrix in PostScript

class(PostScript_T),INTENT(INOUT)      :: self
real(kind=sgl),INTENT(IN)              :: g(3,3)
real(kind=sgl),INTENT(IN)              :: x
real(kind=sgl),INTENT(IN)              :: y
character(3),INTENT(IN)                :: tt

real(kind=sgl)                         :: dy,dx

! set the matrix label
 call self%setlinewidth(0.012)
 call self%setfont(PSfonts(4),0.12)
 call self%text(x-0.1,y+0.3,tt)

! draw the left bracket
 dy=0.26
 dx=x+0.1
 call self%setlinewidth(0.004)
 call self%line(dx+0.1,y-0.1,dx,y-0.1)
 call self%line(dx,y-0.1,dx,y+2.5*dy)
 call self%line(dx,y+2.5*dy,dx+0.1,y+2.5*dy)
 call self%setlinewidth(0.012)

! write all the entries
 write (self%psunit,900) dx,y+2.0*dy,g(1,1)
 write (self%psunit,900) dx,y+1.0*dy,g(2,1)
 write (self%psunit,900) dx,y+0.0*dy,g(3,1)
 dx=dx+0.5
 write (self%psunit,900) dx,y+2.0*dy,g(1,2)
 write (self%psunit,900) dx,y+1.0*dy,g(2,2)
 write (self%psunit,900) dx,y+0.0*dy,g(3,2)
 dx=dx+0.5
 write (self%psunit,900) dx,y+2.0*dy,g(1,3)
 write (self%psunit,900) dx,y+1.0*dy,g(2,3)
 write (self%psunit,900) dx,y+0.0*dy,g(3,3)
 dx=dx+0.65

! and conclude with the right bracket
 call self%setlinewidth(0.004)
 call self%line(dx-0.1,y-0.1,dx,y-0.1)
 call self%line(dx,y-0.1,dx,y+2.5*dy)
 call self%line(dx,y+2.5*dy,dx-0.1,y+2.5*dy)
 call self%setlinewidth(0.012)

 900    format (1x,F12.7,' ',F12.7,' M (',F12.5,') show')

end subroutine DumpMatrix_

!--------------------------------------------------------------------------
subroutine DumpAtom_(self, x, y, A, AT, coor)
!dec$ attributes dllexport :: DumpAtom_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! format atom info for asymmetric unit

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)       :: self

real(kind=sgl),INTENT(IN)               :: x,y,coor(5)
character(2),INTENT(IN)                 :: A
integer(kind=irg),INTENT(IN)            :: AT
real(kind=sgl)                          :: dx

 write (self%psunit,900) x,y,A,AT
 dx=x+0.40
 write (self%psunit,910) dx,y,coor(1)
 dx=dx+0.36
 write (self%psunit,910) dx,y,coor(2)
 dx=dx+0.36
 write (self%psunit,910) dx,y,coor(3)
 dx=dx+0.36
 write (self%psunit,910) dx,y,coor(4)
 dx=dx+0.36
 write (self%psunit,910) dx,y,coor(5)

! formats
 900    format (1x,F12.7,' ',F12.7,' M (',A2,' [',I2,'] ) show')
 910    format (1x,F12.7,' ',F12.7,' M (',f6.4,') show')

end subroutine DumpAtom_

!--------------------------------------------------------------------------
subroutine StrucFacPage_(self, cell, SG, Diff)
!dec$ attributes dllexport :: StrucFacPage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! Structure Factor Information Page

use mod_io
use mod_symmetry
use mod_crystallography
use mod_diffraction

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)      :: self
type(Cell_T),INTENT(INOUT)             :: cell
type(SpaceGroup_T),INTENT(INOUT)       :: SG
type(Diffraction_T),INTENT(INOUT)      :: Diff

type(IO_T)                             :: Message
integer(kind=irg),parameter            :: inm = 4
integer(kind=irg),allocatable          :: family(:,:,:),numfam(:),idx(:), itmp(:,:)
integer(kind=irg)                      :: i,j,ii,ier,h,k,l,totfam,ind(3),icnt, oi_int(1),num
logical                                :: first
logical,allocatable                    :: z(:,:,:)
real(kind=sgl)                         :: g(3),twopi,xs,rag, thr
real(kind=sgl),allocatable             :: Vgg(:,:),ddg(:),gg(:),xi(:),xip(:),twth(:)
character(1)                           :: space
type(gnode)                            :: rlp

! start page
 call self%newpage(.TRUE.,'Structure Factor Information')

! write text
 call self%text(4.0,0.05,'Distances [nm], angles [mrad], potential [V,rad]')

! label columns
 xs = 0.0
 call self%setfont(PSfonts(3),0.14)
 call self%text(xs+1.03,8.15,'h ')
 call self%text(xs+1.18,8.15,'k ')
 call self%text(xs+1.33,8.15,'l ')
 call self%text(xs+1.53,8.15,'p ')
 call self%text(xs+1.93,8.15,'d ')
 call self%text(xs+2.38,8.15,'g ')
 call self%setfont(PSfonts(1),0.14)
 call self%text(xs+2.83,8.15,'2q ')
 call self%setfont(PSfonts(3),0.14)
 call self%text(xs+3.33,8.15,'|V| ')
 call self%text(xs+3.68,8.15,'Vphase ')
 call self%text(xs+4.18,8.15,'|V''| ')
 call self%text(xs+4.58,8.15,'V''phase ')
 call self%setfont(PSfonts(1),0.14)
 call self%text(xs+5.23,8.15,'x ')
 call self%setfont(PSfonts(3),0.14)
 call self%text(xs+5.31,8.13,'g ')
 call self%setfont(PSfonts(1),0.14)
 call self%text(xs+5.58,8.15,'x''')
 call self%setfont(PSfonts(3),0.14)
 call self%text(xs+5.66,8.13,'g ')

! initialize parameters
 thr = 1.E-5
 twopi = 2.0*cPi
 space = 'r'

! allocate all arrays
 allocate(z(-inm:inm,-inm:inm,-inm:inm))
 ii = (2*inm+1)**3
 allocate(family(ii,48,3))
 allocate(numfam(ii))
 allocate(Vgg(ii,6))
 allocate(ddg(ii))
 allocate(gg(ii))
 allocate(xi(ii))
 allocate(xip(ii))
 allocate(twth(ii))
! determine the families of reflections with (hkl)<=(inm,inm,inm)
! first initialize the boolean array z
 z(-inm:inm,-inm:inm,-inm:inm) = .FALSE.
! then loop through all (hkl) values
 first = .TRUE.
 icnt = 1
 totfam=0
 do h=-inm,inm
  ind(1)=-h
  do k=-inm,inm
   ind(2)=-k
   do l=-inm,inm
    ind(3)=-l

! make sure we have not already done this one in another family
    if (.not.z(-h,-k,-l)) then

! if it is a new one, then determine the entire family
     call Diff%CalcUcg(cell,ind)
     rlp = Diff%getrlp()

! but ignore the reciprocal lattice point if Vgg is small
     if (rlp%Vmod.ge.thr) then

! copy family in array and label all its members in z-array
      call SG%CalcFamily(ind,num,space,itmp)
      do i=1,num
       do j=1,3
        family(icnt,i,j)=itmp(i,j)
       end do
       z(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
      end do

! store the Fourier coefficient of the lattice potential
      Vgg(icnt,1)=rlp%Vmod
      Vgg(icnt,2)=rlp%Vphase
      Vgg(icnt,3)=rlp%Vpmod
      Vgg(icnt,4)=rlp%Vpphase
      Vgg(icnt,5)=rlp%xgp
      Vgg(icnt,6)=rlp%xg

! increment family counter
      numfam(icnt)=num
      totfam=totfam+num-1
      icnt=icnt+1
     end if
    end if
   end do
  end do
 end do

 icnt=icnt-1
 oi_int(1)=icnt
 call Message%WriteValue(' Total number of families        = ', oi_int, 1, "(I6)")
 oi_int(1)=totfam
 call Message%WriteValue(' Total number of family members  = ', oi_int, 1, "(I6)")

! compute d-spacings, g-spacings, two-theta, and extinction distance
 do k=1,icnt
  g(1:3)=float(family(k,1,1:3))
  gg(k)=cell%CalcLength(g,'r')
  ddg(k)=1.0/gg(k)
  twth(k)=2.0*asin(0.5*Diff%getWaveLength()*gg(k))
  xip(k) = Vgg(k,5)
  xi(k)  = Vgg(k,6)
 end do

! take care of (0,0,0) reflection
 ind(1)=0
 ind(2)=0
 ind(3)=0
 call Diff%CalcUcg(cell,ind)
 rlp = Diff%getrlp()
 Vgg(icnt,1)=rlp%Vmod
 Vgg(icnt,2)=rlp%Vphase
 Vgg(icnt,3)=rlp%Vpmod
 Vgg(icnt,4)=rlp%Vpphase
 gg(icnt)=0.0
 twth(icnt)=0.0
 xi(icnt)=0.0
 xip(icnt) = rlp%xgp

! use the spsort.f routine from SLATEC
! to rank by increasing value of gg
 allocate(idx(ii))
 call SPSORT(gg,icnt,idx,1,ier)

! and create output table
 i=1
 j=icnt
 rag = 0.0
 call self%DumpLine(i,family(j,1,1),family(j,1,2),family(j,1,3),numfam(j),ddg(j),rag,twth(j), &
               Vgg(j,1),Vgg(j,2),Vgg(j,3),Vgg(j,4),xi(j),xip(j))
 do i=2,icnt
  j=idx(i)
  if (ddg(j).ne.0.0) then
   rag = 1.0/ddg(j)
  else
   rag = 0.0
  endif
 call self%DumpLine(i,family(j,1,1),family(j,1,2),family(j,1,3),numfam(j),ddg(j),rag,twth(j), &
               Vgg(j,1),Vgg(j,2),Vgg(j,3),Vgg(j,4),xi(j),xip(j))
 end do


deallocate(idx)
deallocate(z)
deallocate(family)
deallocate(numfam)
deallocate(Vgg)
deallocate(ddg)
deallocate(gg)
deallocate(xi)
deallocate(xip)
deallocate(twth)

! integer number output
 901    format (1x,F12.7,' ',F12.7,' M (',I8,') show')

end subroutine StrucFacPage_

!--------------------------------------------------------------------------
subroutine DumpLine_(self,i,h,k,l,p,d,g,th,vm,vp,vpm,vpp,xi,xip)
!dec$ attributes dllexport :: DumpLine_

  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! formats a line of structure factor information

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)      :: self
integer                                :: h,k,l,p,i
real                                   :: d,g,th,xi,vm,vp,vpm,vpp,xip, x, y, xs

! start a new page (if needed) and label columns
 if ((mod(i-1,50).eq.0).AND.(i.ne.1)) then
  call self%newpage(.TRUE.,'Structure Factor Information')

! write text
  call self%text(4.0,0.05,'Distances [nm], angles [mrad], potential [V,rad]')

! label columns
  xs = 0.0
  call self%setfont(PSfonts(3),0.14)
  call self%text(xs+1.03,8.15,'h ')
  call self%text(xs+1.18,8.15,'k ')
  call self%text(xs+1.33,8.15,'l ')
  call self%text(xs+1.53,8.15,'p ')
  call self%text(xs+1.93,8.15,'d ')
  call self%text(xs+2.38,8.15,'g ')
  call self%setfont(PSfonts(1),0.14)
  call self%text(xs+2.83,8.15,'2q ')
  call self%setfont(PSfonts(3),0.14)
  call self%text(xs+3.33,8.15,'|V| ')
  call self%text(xs+3.68,8.15,'Vphase ')
  call self%text(xs+4.18,8.15,'|V''| ')
  call self%text(xs+4.58,8.15,'V''phase ')
  call self%setfont(PSfonts(1),0.14)
  call self%text(xs+5.23,8.15,'x ')
  call self%setfont(PSfonts(3),0.14)
  call self%text(xs+5.31,8.13,'g ')
  call self%setfont(PSfonts(1),0.14)
  call self%text(xs+5.58,8.15,'x''')
  call self%setfont(PSfonts(3),0.14)
  call self%text(xs+5.66,8.13,'g ')
 end if

! put all entries in the correct positions
 x = xs+1.0
 y = 8.0 - (mod(i-1,50))*0.15
 call self%setfont(PSfonts(4),0.12)
 write (self%psunit,900) x,y,h
 x=x+0.15
 write (self%psunit,900) x,y,k
 x=x+0.15
 write (self%psunit,900) x,y,l
 x=x+0.20
 write (self%psunit,900) x,y,p
 x=x+0.25
 write (self%psunit,901) x,y,d
 x=x+0.45
 write (self%psunit,901) x,y,g
 x=x+0.45
 write (self%psunit,902) x,y,th*1000.0
 x=x+0.55
 write (self%psunit,901) x,y,vm
 x=x+0.45
 write (self%psunit,901) x,y,vp
 x=x+0.45
 write (self%psunit,901) x,y,vpm
 x=x+0.45
 write (self%psunit,901) x,y,vpp
 x=x+0.45
 if (xi.gt.100000.0) then
  write (self%psunit,905) x,y
 else
  write (self%psunit,906) x,y,xi
 endif
 x=x+0.45
 if (xip.gt.100000.0) then
  write (self%psunit,905) x,y
 else
  write (self%psunit,906) x,y,xip
 endif
!
 900    format (1x,F12.7,' ',F12.7,' M (',I2,') show')
 901    format (1x,F12.7,' ',F12.7,' M (',F8.5,') show')
 902    format (1x,F12.7,' ',F12.7,' M (',F9.4,') show')
 903    format (1x,F12.7,' ',F12.7,' M (',F12.4,') show')
 904    format (1x,F12.7,' ',F12.7,' M (',I8,') show')
 905    format (1x,F12.7,' ',F12.7,' M ( >100,000 ) show')
 906    format (1x,F12.7,' ',F12.7,' M (',F8.2,') show')

end subroutine DumpLine_

!--------------------------------------------------------------------------
subroutine StereoPage_(self, cell, SG)
!dec$ attributes dllexport :: StereoPage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! creates pages with stereographic projections

use mod_io
use mod_crystallography
use mod_symmetry

class(PostScript_T),INTENT(INOUT)     :: self
type(Cell_T),INTENT(INOUT)            :: cell
type(SpaceGroup_T),INTENT(INOUT)      :: SG

type(IO_T)                            :: Message
character(1)                          :: sp
logical                               :: topbot
integer(kind=irg)                     :: hm,km,lm,iview(3), io_int(3)

 topbot=.TRUE.
 call Message%printMessage( (/ 'Enter the maximum index for h,k and l, or for ',&
                               'u,v, and w. For a hexagonal system, please use',&
                               '4-index notation [uv.w] or (hk.l) to determine',&
                               'the largest index.                            ' /), frm = "(A)")
 call Message%ReadValue(' Enter maximum indices : ', io_int, 3)
 hm = io_int(1)
 km = io_int(2)
 lm = io_int(3)
! first [001] projection of real space
 sp = 'd'
 iview(1)=0
 iview(2)=0
 iview(3)=1
! call the drawing routine
 call self%StereoProj(cell,SG,sp,iview,hm,km,lm,topbot)
! then [001] projection of reciprocal space
 sp = 'r'
! call the drawing routine
 call self%StereoProj(cell,SG,sp,iview,hm,km,lm,topbot)

end subroutine StereoPage_

!--------------------------------------------------------------------------
recursive subroutine DiffPage_(self, cell, SG, Diff, camlen)
!dec$ attributes dllexport :: DiffPage_
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !!  draw kinematical zone axis electron diffraction patterns

use mod_crystallography
use mod_symmetry
use mod_io
use mod_diffraction
use mod_misc

IMPLICIT NONE

class(PostScript_T),INTENT(INOUT)   :: self
type(Cell_T),INTENT(INOUT)          :: cell
type(SpaceGroup_T),INTENT(INOUT)    :: SG
class(Diffraction_T),INTENT(INOUT)  :: Diff
real(kind=sgl),INTENT(IN)           :: camlen

type(IO_T)                          :: Message

integer(kind=irg),parameter         :: inm = 5
character(1)                        :: list(256)
logical                             :: first,np,ppat,a
logical,allocatable                 :: z(:,:,:),zr(:,:,:)
integer(kind=irg)                   :: i,j,h,k,l,m,totfam,hh,ll,fmax,inmhkl(3),ricnt,icnt,ind(3),uu,vv,ww,slect(256), &
                                       js,ii,num,hc,hhcc,iinm,dpcnt,imo,ih,ik,il,ier,iref, io_int(4)
integer(kind=irg),allocatable       :: idx(:)
integer(kind=irg),allocatable       :: family(:,:),numfam(:)
real(kind=sgl)                      :: twopi,ggl,g(3),Vmax,laL,gmax,RR,thr, oi_real(1)
real(kind=sgl),allocatable          :: gg(:)
real(kind=sgl)                      :: rmt(3,3)
real(kind=sgl),parameter            :: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/), &
                                       eps = 1.0E-3
logical,allocatable                 :: dbdiff(:)
integer(kind=irg),allocatable       :: itmp(:,:)   !< array used for family computations etc

real(kind=sgl),allocatable          :: Vg(:),rg(:),Vgsave(:)
integer(kind=irg),allocatable       :: rfamily(:,:,:),rnumfam(:)
type(gnode)                         :: rlp

! set some parameters
 thr = 1.E-4
 twopi = 2.0*cPi
 Vmax = 0.0

! gmax is the radius of the sphere whose intersection with the
! back focal plane is the circle on the output zone axis patterns
 laL = sngl(Diff%getWaveLength()) * camlen
 RR = 1.375 * 25.4
 gmax = RR/laL

 oi_real(1) = sngl(Diff%getWaveLength())
 call Message%WriteValue('wavelength [nm] = ', oi_real, 1)
 oi_real(1) = camlen
 call Message%WriteValue(' L         [mm] = ', oi_real, 1)
 oi_real(1) = laL
 call Message%WriteValue('camera length lambda*L [mm nm] = ', oi_real, 1)

! determine the families of reciprocal lattice points
! in a region of reciprocal space.
! set the index boundaries

rmt = cell%getrmt()

 do i=1,3
  inmhkl(i) = 2*int(gmax/sqrt(rmt(i,i)))
 end do
 hc = maxval(inmhkl)

! allocate all the necessary arrays
 allocate(zr(-hc:hc,-hc:hc,-hc:hc))
 allocate(z(-inm:inm,-inm:inm,-inm:inm))
 hhcc = (2*hc+1)**3
 allocate(Vg(hhcc))
 allocate(Vgsave(hhcc))
 allocate(rfamily(hhcc,48,3))
 allocate(rnumfam(hhcc))
 allocate(rg(hhcc))
 iinm = (2*inm+1)**3
 allocate(family(iinm,3))
 allocate(numfam(iinm))

! if this is a non-symmorphic space group, then also
! allocate the dbdiff array to tag potential double
! diffraction reflections
 if (.not.SG%getSpaceGroupSymmorphic()) then
   allocate(dbdiff(hhcc))
   dbdiff(1:hhcc) = .FALSE.
 endif

! and initialize the ones that need to be initialized
 zr(-hc:hc,-hc:hc,-hc:hc) = .FALSE.
 z(-inm:inm,-inm:inm,-inm:inm) = .FALSE.

! here we go:
 first = .TRUE.
 icnt = 1
 totfam=0
 do h=-inmhkl(1),inmhkl(1)
  ind(1)=h
  do k=-inmhkl(2),inmhkl(2)
   ind(2)=k
   do l=-inmhkl(3),inmhkl(3)
    ind(3)=l
! make sure we have not already done this one in another family
    if (.not.zr(h,k,l)) then
! check the length, to make sure it lies within the sphere gmax
     ggl = cell%CalcLength(float(ind),'r')
! if it is larger than gmax, then compute the entire family
     if (ggl.ge.gmax) then
      call SG%CalcFamily(ind,num,'r',itmp)
! and label the family members in the zr array so that we
! do not include those points later on in the loop
      do i=1,num
       zr(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
      end do
     else
! if the length is smaller than gmax, then compute the
! Fourier coefficient Vg and determine the entire family
! [recall that all members in a family have the same Vg]
! Do this only for those reflections that are allowed by
! the lattice centering !
      a = SG%IsGAllowed((/h,k,l/))
      if (a) then
       call Diff%CalcUcg(cell,ind)
       rlp = Diff%getrlp()

! check for nonsymmorphic systematic absences
       if ((.not.SG%getSpaceGroupSymmorphic()).and.(rlp%Vmod.lt.eps)) then
        io_int = (/ h, k, l, 0 /)
        call Message%WriteValue(' potential double diffraction family :', io_int,3,"('{',I3,I3,I3,'}')")
        dbdiff(icnt) = .TRUE.
        rlp%Vmod = 0.0
       endif
! compute the entire family
       deallocate(itmp)
       call SG%CalcFamily(ind,num,'r',itmp)
       rg(icnt)=ggl
! copy family in array
       do i=1,num
        rfamily(icnt,i,1:3)=itmp(i,1:3)
        zr(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
       end do
! and take the modulus squared for the intensity
       Vg(icnt)=rlp%Vmod**2
       Vgsave(icnt)=Vg(icnt)
! update the maximum value
       Vmax = max(Vg(icnt),Vmax)
      else
! remove the equivalent systematic absences
       deallocate(itmp)
       call SG%CalcFamily(ind,num,'r',itmp)
       rg(icnt)=ggl
       do i=1,num
        rfamily(icnt,i,1:3)=itmp(i,1:3)
        zr(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
       end do
! and put the intensity to a negative value
! that way we will know whether or not to draw them
       Vg(icnt)=-100.0
       Vgsave(icnt)=Vg(icnt)
      end if
! and increment the family counter
      rnumfam(icnt)=num
      totfam=totfam+num-1
      icnt=icnt+1
     end if
    end if
   end do
  end do
 end do

 icnt=icnt-1
! normalize potential coefficients to largest one
! and scale in a non-linear way to mimic density on
! an electron micrograph [Gonzalez & Windtz]
! Use the where operator to avoid the negative intensities
 call Message%ReadValue('logarithmic[0] or exponential[1] intensity scale ', io_int, 1)
 ll = io_int(1)
 if (ll.eq.0) then
  where(Vg.gt.0.0) Vg=0.01*alog(1.0+0.1*Vg)
 else
  where(Vg.gt.0.0) Vg=0.05*(Vg/Vmax)**0.2
 end if
 ricnt=icnt

! determine families of directions
 first = .TRUE.
 icnt = 1
 totfam=0
 do uu=-inm,inm
  do vv=-inm,inm
   do ww=-inm,inm
    if ((uu**2+vv**2+ww**2).ne.0) then
! make sure we have not already done this one in another family
     ind= (/ -uu, -vv, -ww /)
     call IndexReduce(ind)
     if (.not.z(ind(1),ind(2),ind(3))) then
! determine the family <uvw>
      deallocate(itmp)
      call SG%CalcFamily(ind,num,'d',itmp)
! and keep only one family member, namely the one with the
! largest sum of the three integers, i.e. u+v+w
! [this is a simple way to get mostly positive indices as
! the zone axis indices]
      js = -100
      ii = 0
      do i=1,num
       hh = itmp(i,1)+itmp(i,2)+itmp(i,3)
       if (hh.gt.js) then
        ii = i
        js = hh
       end if
! then remove the multiples of those direction indices from list
       do m=-inm,inm
        ih=itmp(i,1)*m
        ik=itmp(i,2)*m
        il=itmp(i,3)*m
        if (((abs(ih).le.inm).and.(abs(ik).le.inm)).and.(abs(il).le.inm)) then
         z(ih,ik,il)=.TRUE.
        end if
       end do
      end do
      family(icnt,1:3)=itmp(ii,1:3)
! increment family counter
      numfam(icnt)=num
      totfam=totfam+num-1
      icnt=icnt+1
     end if
    end if
   end do
  end do
 end do
 icnt=icnt-1
 io_int(1) = icnt
 call Message%WriteValue('->Total number of direction families = ', io_int, 1)

! compute length of direction vectors and rank by increasing length
 allocate(idx(icnt))
 allocate(gg(icnt))
 gg(1:icnt) = 0.0
 do k=1,icnt
  g(1:3)=float(family(k,1:3))
  gg(k)=cell%CalcLength(g,'d')
 end do

! rank by increasing value of gg (use SLATEC routine)
 call SPSORT(gg,icnt,idx,1,ier)

! ask for number to be included in output
 call Message%printMessage('List of available zone axis patterns', frm = "(A)")
 do i=1,icnt
  j=idx(i)
  io_int(1)=i
  do k=1,3
   io_int(k+1) = family(j,k)
  end do
  if (mod(i,4).eq.0) then
   call Message%WriteValue('', io_int,4,"(I3,' [',3I3,'];')")
  else
   call Message%WriteValue('', io_int,4,"(I3,' [',3I3,'];')",advance="no")
  endif
 end do
 call Message%printMessage('Enter selection (e.g. 4,10-20, ... ) ', frm = "(//,A)")
 call Message%printMessage('[Include 0 to also draw a powder pattern] ', frm = "(A)")
 list = (/ (' ',j=1,256) /)
 call Message%printMessage(' -> ', frm = "(A,' ')",advance="no")
 read (5,"(256A)") list
 call studylist(list,slect,fmax,ppat)

 if (.not.SG%getSpaceGroupSymmorphic()) then
  call Message%printMessage('Potential double diffraction reflections will be indicated by open squares.', frm = "(A,/)")
 end if
 call Message%ReadValue('No indices (0), labels (1), extinctions (2), labels + extinctions (3): ', io_int, 1)
 iref = io_int(1)

! and create output in 2 columns, 3 rows
 do i=1,fmax
  dpcnt=dpcnt+1
  j=idx(slect(i))
  imo = mod(i-1,6)
  if (imo.eq.0) then
   np=.TRUE.
  else
   np=.FALSE.
  endif
  if (i.eq.1) then
   first=.TRUE.
  else
   first=.FALSE.
  endif
  if (slect(i).eq.0) then
   call Message%printMessage('Creating Powder Pattern ', frm = "(A)")
   call self%DumpPP(cell,xoff(imo),yoff(imo),np,laL,ricnt,Vgsave,rg,rnumfam)
   ppat=.FALSE.
  else
   io_int(1:3) = family(j,1:3)
   call Message%WriteValue('Creating ZAP ', io_int,3, "('[',3i3,'] : ')",advance="no")
   call self%DumpZAP(cell,SG,xoff(imo),yoff(imo),family(j,1),family(j,2),family(j,3),numfam(j),np,first,iref,laL,ricnt,dbdiff, &
        Vg, Vgsave, rg, rfamily, rnumfam, hhcc)
  endif
 end do

! and clean up all variables
 deallocate(zr,z,Vg, Vgsave, rfamily, rnumfam, rg, family, numfam, idx, gg)
 if (.not.SG%getSpaceGroupSymmorphic()) deallocate(dbdiff)

end subroutine DiffPage_

!--------------------------------------------------------------------------
recursive subroutine studylist(list,slect,np,ppat)
!DEC$ ATTRIBUTES DLLEXPORT :: studylist
  !! author: MDG
  !! version: 1.0
  !! date: 01/27/20
  !!
  !! determine which zone axis pattern to draw from user input
  !! the parameter ppat is returned as true if the user also requested a powder pattern

IMPLICIT NONE

character(1),INTENT(IN)                 :: list(256)            !< input string
integer(kind=irg),INTENT(OUT)           :: slect(256)           !< list of patterns to be drawn
integer(kind=irg),INTENT(OUT)           :: np                           !< number of patterns
logical,INTENT(INOUT)                   :: ppat                 !< powder pattern included ?
!f2py intent(in,out) ::  ppat                 !< powder pattern included ?

integer(kind=irg)                       :: comma(100),hyphen(100),ccnt,hcnt,i,j,k,ip,icnt,nd,n,istart,istop
integer(kind=irg),parameter             :: nmb(48:57)=(/0,1,2,3,4,5,6,7,8,9/)

! initialize a few parameters
 ccnt = 0
 hcnt = 0
 ppat = .FALSE.
 slect = 0
 comma = 0
 hyphen= 0
 j = 0

! count characters and search for , and -
 do i=1,256
  if (list(i)(1:1).ne.' ') j=j+1
  if (list(i)(1:1).eq.',') then
   ccnt = ccnt+1
   comma(ccnt)=i
  end if
  if (list(i)(1:1).eq.'-') then
   hcnt = hcnt+1
   hyphen(hcnt)=i
  end if
 end do
 ccnt = ccnt+1
 comma(ccnt) = j+1

! interpret the string
 j = 1
 ip = 1
 icnt = 0
 do i=1,ccnt
! is it a range ?
  if (((hyphen(j).lt.comma(i)).and.(hcnt.gt.0)).and.(j.le.hcnt)) then

! yes, it is;  get the first number
   nd = hyphen(j)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   istart = n
   ip = hyphen(j)+1

! and then the second number
   nd = comma(i)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   istop = n

! and fill in the entire range
   do k=istart,istop
    icnt=icnt+1
    slect(icnt)=k
    if (k.eq.0) then
     ppat = .TRUE.
    end if
   end do
   ip = comma(i)+1
   j=j+1
  else

! no, it is not a range; determine number of digits
   nd = comma(i)-ip
   n = 0
   do k=0,nd-1
    n = 10*n+nmb(ichar(list(ip+k)(1:1)))
   end do
   icnt=icnt+1
   slect(icnt)=n
   if (n.eq.0) then
    ppat = .TRUE.
   end if
   ip = comma(i)+1
  end if
 end do
 np = icnt

end subroutine studylist


end module mod_postscript
