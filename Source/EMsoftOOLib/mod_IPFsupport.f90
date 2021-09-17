! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_IPFsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 09/03/21
  !!
  !! module for inverse pole figure maps

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: IPFmap_T
private 
  integer(kind=irg)   :: ipf_wd
  integer(kind=irg)   :: ipf_ht
  integer(kind=irg)   :: ipf_LaueClass
  integer(kind=irg)   :: ipf_nthreads
  character(fnlen)    :: ipf_mode
  character(fnlen)    :: ipf_filename

contains
private 
  procedure, pass(self) :: get_ipf_wd_
  procedure, pass(self) :: get_ipf_ht_
  procedure, pass(self) :: get_ipf_mode_
  procedure, pass(self) :: get_ipf_filename_
  procedure, pass(self) :: get_ipf_LaueClass_
  procedure, pass(self) :: get_ipf_nthreads_
  procedure, pass(self) :: set_ipf_wd_
  procedure, pass(self) :: set_ipf_ht_
  procedure, pass(self) :: set_ipf_mode_
  procedure, pass(self) :: set_ipf_filename_
  procedure, pass(self) :: set_ipf_LaueClass_
  procedure, pass(self) :: set_ipf_nthreads_
  procedure, pass(self) :: get_ipf_RGB_
  procedure, pass(self) :: get_IPFMap_

  generic, public :: get_ipf_wd => get_ipf_wd_
  generic, public :: get_ipf_ht => get_ipf_ht_
  generic, public :: get_ipf_mode => get_ipf_mode_
  generic, public :: get_ipf_filename => get_ipf_filename_
  generic, public :: get_ipf_LaueClass => get_ipf_LaueClass_
  generic, public :: get_ipf_nthreads => get_ipf_nthreads_
  generic, public :: set_ipf_wd => set_ipf_wd_
  generic, public :: set_ipf_ht => set_ipf_ht_
  generic, public :: set_ipf_mode => set_ipf_mode_
  generic, public :: set_ipf_filename => set_ipf_filename_
  generic, public :: set_ipf_LaueClass => set_ipf_LaueClass_
  generic, public :: set_ipf_nthreads => set_ipf_nthreads_
  generic, public :: get_ipf_RGB => get_ipf_RGB_
  generic, public :: get_IPFMap => get_IPFMap_

end type IPFmap_T

! the constructor routine for this class 
interface IPFmap_T
  module procedure IPF_constructor
end interface IPFmap_T

contains

!--------------------------------------------------------------------------
type(IPFmap_T) function IPF_constructor( ) result(IPF)
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! constructor for the IPFmap_T Class; reads the name list 
 
IMPLICIT NONE


end function IPF_constructor

!--------------------------------------------------------------------------
subroutine IPF_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! destructor for the IPFmap_T Class
 
IMPLICIT NONE

type(IPFmap_T), INTENT(INOUT)  :: self 

call reportDestructor('IPFmap_T')

end subroutine IPF_destructor


!--------------------------------------------------------------------------
function get_ipf_wd_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_wd_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_wd from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%ipf_wd

end function get_ipf_wd_

!--------------------------------------------------------------------------
subroutine set_ipf_wd_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_wd_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_wd in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_wd = inp

end subroutine set_ipf_wd_

!--------------------------------------------------------------------------
function get_ipf_ht_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_ht_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_ht from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%ipf_ht

end function get_ipf_ht_

!--------------------------------------------------------------------------
subroutine set_ipf_ht_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_ht_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_ht in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_ht = inp

end subroutine set_ipf_ht_

!--------------------------------------------------------------------------
function get_ipf_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_nthreads from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg)               :: out

out = self%ipf_nthreads

end function get_ipf_nthreads_

!--------------------------------------------------------------------------
subroutine set_ipf_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_nthreads in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_nthreads = inp

end subroutine set_ipf_nthreads_

!--------------------------------------------------------------------------
function get_ipf_filename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_filename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_filename from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen)                :: out

out = self%ipf_filename

end function get_ipf_filename_

!--------------------------------------------------------------------------
subroutine set_ipf_filename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_filename_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_filename in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen), INTENT(IN)    :: inp

self%ipf_filename = inp

end subroutine set_ipf_filename_

!--------------------------------------------------------------------------
function get_ipf_mode_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_mode_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_mode from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen)                :: out

out = self%ipf_mode

end function get_ipf_mode_

!--------------------------------------------------------------------------
subroutine set_ipf_mode_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_mode_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_mode in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
character(fnlen), INTENT(IN)    :: inp

self%ipf_mode = inp

end subroutine set_ipf_mode_

!--------------------------------------------------------------------------
function get_ipf_LaueClass_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_LaueClass_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! get ipf_LaueClass from the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
integer(kind=irg)               :: out

out = self%ipf_LaueClass

end function get_ipf_LaueClass_

!--------------------------------------------------------------------------
subroutine set_ipf_LaueClass_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ipf_LaueClass_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! set ipf_LaueClass in the IPFmap_T class

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)  :: self
integer(kind=irg), INTENT(IN)   :: inp

self%ipf_LaueClass = inp

end subroutine set_ipf_LaueClass_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive function in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax ) result(inside)
!! author: MDG 
!! version: 1.0 
!! date: 09/06/21
!!
!! generic unit triangle test 

IMPLICIT NONE 

real(kind=dbl),INTENT(IN)       :: etaDeg
real(kind=dbl),INTENT(IN)       :: chiDeg
real(kind=dbl),INTENT(IN)       :: etaMin 
real(kind=dbl),INTENT(IN)       :: etaMax 
real(kind=dbl),INTENT(IN)       :: chiMin 
real(kind=dbl),INTENT(IN)       :: chiMax 
logical                         :: inside 

real(kind=dbl),parameter        :: eps = 1.0D-12

inside = .FALSE.

if ( (etaDeg.ge.etaMin) .and. (etaDeg.le.(etaMax+eps)) ) then 
  if ( (chiDeg.ge.chiMin) .and. (chiDeg.le.(chiMax+eps)) )  then 
    inside = .TRUE.
  end if 
end if 

end function in_generic_unit_triangle

!--------------------------------------------------------------------------
recursive function in_cubic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMax ) result(inside)
!! author: MDG 
!! version: 1.0 
!! date: 09/06/21
!!
!! cubic unit triangle test 

IMPLICIT NONE 

real(kind=dbl),INTENT(IN)       :: etaDeg
real(kind=dbl),INTENT(IN)       :: chiDeg
real(kind=dbl),INTENT(IN)       :: etaMin
real(kind=dbl),INTENT(IN)       :: etaMax 
real(kind=dbl),INTENT(OUT)      :: chiMax 
logical                         :: inside 

real(kind=dbl),parameter        :: eps = 1.0D-12

inside = .FALSE.

! does this lie in the cubic high unit triangle ? 
if (etaDeg.ge.45D0) then 
  chiMax = dsqrt(1.D0 / (2.D0 + tan((90.D0 - etaDeg)*dtor)**2 ))
else 
  chiMax = dsqrt(1.D0 / (2.D0 + tan(etaDeg*dtor)**2 ))
end if 
if (chiMax.gt.1.D0) chiMax = 1.D0
chiMax = dacos(chiMax)
if ( (etaDeg.ge.0.D0) .and. (etaDeg.le.(etaMax + eps)) ) then 
  if ( (chiDeg.ge.0.D0) .and. (chiDeg.le.((chiMax + eps ) * rtod) ) )  then 
    inside = .TRUE.
  end if 
end if 

end function in_cubic_unit_triangle

!--------------------------------------------------------------------------
recursive function get_ipf_RGB_(self, sampleDir, qu, sym, Pm, clr ) result(RGB)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ipf_RGB_
!! author: MDG 
!! version: 1.0 
!! date: 09/03/21
!!
!! returns an RGB triplet for the specific Laue Class, sample direction, and grain orientation
!! This code is based on the LaueOps routines from DREAM.3D, but combines all Laue classes into
!! a single routine.
!!
!! At present, only Euler coloring and standard Laue Class RGB coloring are available.
!! The more extensive color schemes proposed in [Nolze, G., & Hielscher, R. (2016). Orientations–perfectly colored. 
!! Journal of Applied Crystallography, 49(5), 1786-1802.], which have minimal color discontinuities
!! will be implemented in a future version of this module.

use mod_quaternions
use mod_rotations 
use mod_io
use mod_colorspace

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)        :: self
real(kind=dbl), INTENT(IN)            :: sampleDir(3)
type(Quaternion_T), INTENT(IN)        :: qu
type(QuaternionArray_T), INTENT(IN)   :: sym 
integer(kind=irg), INTENT(IN)         :: Pm
integer(kind=ish)                     :: RGB(3)
type(colorspace_T),INTENT(IN),OPTIONAL:: clr

type(e_T)                             :: eu
type(q_T)                             :: qq
type(IO_T)                            :: Message
type(Quaternion_T)                    :: qs, qc, qtest 

real(kind=dbl)                        :: e(3), tp, hsl(0:2) 
real(kind=dbl)                        :: refDir(3), chi, eta, etaDeg, chiDeg, RGBd(3), ma 
real(kind=dbl)                        :: chiMin, chiMax, etaMin, etaMax
integer(kind=irg)                     :: i 
logical                               :: inside

! do a simple test 
! qtest = Quaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0 /) )
! RGB = self%get_RGB_cubic_high_( sampleDir, qtest, sym, Pm )
! write (*,*) ' RGB [001] = ', RGB 
! qtest = Quaternion_T( qd = (/ cos(0.5D0*45.D0*dtor), 0.D0, sin(0.5D0*45.D0*dtor), 0.D0 /) )
! RGB = self%get_RGB_cubic_high_( sampleDir, qtest, sym, Pm )
! write (*,*) ' RGB [101] = ', RGB 
! stop

! use a simple Euler angle color scheme, mapping all of Euler space linearly onto RGB color space
if (trim(self%ipf_mode).eq.'Euler') then 
  tp = 2.D0 * cPi 
! scale the full Euler space and map onto the RGB space
  qq = q_T( qdinp = qu%get_quatd() )
  eu = qq%qe()
  e = eu%e_copyd()
  RGB = (/ int(e(1) * 255/ tp), int(e(2) * 255 / cPi), int(e(3) * 255/ tp) /) 
end if

! use the standard IPF colors based on the TSL Euler angles
if ( (trim(self%ipf_mode).eq.'TSL') .or. (trim(self%ipf_mode).eq.'OPC') ) then 
  RGB = (/ 0, 0, 0 /)

  findloop: do i=1,Pm 
! apply the symmetry operation
    qs = sym%getQuatfromArray(i)
    qc = qs * qu
    call qc%quat_pos()
! convert the sampleDir to refDir
    refDir = qc%quat_Lp( sampleDir )
    if (refDir(3).lt.0.D0) refDir = -refDir
! convert to stereographic projection angles 
    chi = acos(refDir(3))
    eta = atan2(refDir(2),refDir(1))
    etaDeg = eta * rtod
    chiDeg = chi * rtod
! does this lie in the orthorhombic unit triangle ? 
    select case (self%ipf_LaueClass)
      case(1,2)     ! triclinic, monoclinic
          etaMin=0.D0
          etaMax=180.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(3,4)     ! orthorhombic, tetragonal-low
          etaMin=0.D0
          etaMax=90.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(5)       ! tetragonal-high
          etaMin=0.D0
          etaMax=45.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(6)       ! trigonal-low
          etaMin=-120.D0
          etaMax=0.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(7)       ! trigonal-high
          etaMin=-90.D0
          etaMax=-30.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(8)       ! hexagonal-low
          etaMin=0.D0
          etaMax=60.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(9)       ! hexagonal-high
          etaMin=0.D0
          etaMax=30.D0
          chiMin=0.D0
          chiMax=90.D0
          inside = in_generic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMin, chiMax )
      case(10)      ! cubic-low
          etaMin=0.D0
          etaMax=90.D0
          inside = in_cubic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMax )
      case(11)      ! cubic-high
          etaMin=0.D0
          etaMax=45.D0
          inside = in_cubic_unit_triangle( etaDeg, chiDeg, etaMin, etaMax, chiMax )
      case default 
        call Message%printError('get_ipf_RGB_',' Unknown Laue class')
    end select 

    if (inside.eqv..TRUE.) then ! if .TRUE., then get the correct color and exit the loop
      chiDeg = chiDeg / chiMax
      RGBd(1) = 1.D0 - chiDeg 
      RGBd(3) = abs(etaDeg - etaMin) / (etaMax-etaMin)
      RGBd(2) = (1.D0 - RGBd(3)) * chiDeg
      RGBd(3) = RGBd(3) * chiDeg
      RGBd = sqrt(RGBd) 
      ma = maxval(RGBd)
      RGBd = RGBd / ma 
      if (trim(self%ipf_mode).eq.'OPC') then 
        hsl = clr%rgb2hsl(RGBd)
        RGBd = clr%sph2rgb(hsl)
      end if 
      if (trim(self%ipf_mode).eq.'PUC') then 
        hsl = clr%rgb2hsl(RGBd)
        RGBd = clr%sph2rgb(hsl)
      end if 
      RGB = int(RGBd * 255)
      exit findloop
    end if 
  end do findloop
end if 

end function get_ipf_RGB_

!--------------------------------------------------------------------------
subroutine get_IPFMap_( self, EMsoft, sampleDir, Orientations, sym ) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_IPFMap_
!! author: MDG
!! version: 1.0
!! date: 09/03/21
!!
!! This takes an orientation list in quaternion representation and turns it into 
!! an inverse pole figure map (using one of a few color schemes) that is then stored 
!! in a .tiff file.
!!
!! This routine assumes that the IPF map size and filename as well as the Laue Class 
!! have been set by the calling program using the appropriate methods; the number of orientations
!! in the Orientations array must equal the product self%ipf_wd * self%ipf_ht.

use mod_EMsoft
use mod_quaternions
use mod_OMPsupport
use omp_lib
use mod_image
use mod_io
use mod_colorspace
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(IPFmap_T), INTENT(INOUT)              :: self
type(EMsoft_T), INTENT(INOUT)               :: EMsoft
integer(kind=irg), INTENT(IN)               :: sampleDir(3)
type(QuaternionArray_T), INTENT(IN)         :: Orientations
type(QuaternionArray_T), INTENT(INOUT)      :: sym

integer(kind=irg),allocatable               :: IPFmap(:,:,:)
integer(kind=irg)                           :: ix, iy, iq, TID 
type(Quaternion_T)                          :: qu
type(IO_T)                                  :: Message
type(colorspace_T)                          :: clr

character(fnlen)                            :: fname, TIFF_filename
real(kind=dbl)                              :: sDir(3)

! declare variables for use in object oriented image module
integer                                     :: iostat, io_int(2)
character(len=128)                          :: iomsg
logical                                     :: isInteger, OPC
type(image_t)                               :: im
integer(int8), allocatable                  :: TIFF_image(:,:)
integer                                     :: dim2(2), Pm
integer(c_int32_t)                          :: result

allocate(IPFmap(3, self%ipf_wd, self%ipf_ht))

! make sure sampleDir is normalized 
sDir = dble(sampleDir)/NORM2(dble(sampleDir))

io_int = (/ self%ipf_wd, self%ipf_ht /)
call Message%printMessage(' ')
call Message%WriteValue(' Generating IPF map of dimensions : ', io_int,2,frm="(I5,' x ',I5)")

Pm = sym%getQnumber()

OPC = .FALSE.
if (trim(self%ipf_mode).eq.'OPC') OPC = .TRUE.

io_int(1) = Pm
call Message%WriteValue(' # symops = ', io_int, 1) 

io_int(1) = self%get_ipf_nthreads_()
call Message%WriteValue(' # threads = ', io_int, 1)

call OMP_SET_NUM_THREADS( self%get_ipf_nthreads_() )

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, iy, iq, qu, clr)

TID = OMP_GET_THREAD_NUM()

if (OPC.eqv..TRUE.) clr = colorspace_T()

!$OMP DO SCHEDULE(DYNAMIC)
  do iy = 1, self%ipf_ht 
    do ix = 1, self%ipf_wd 
      iq = (iy-1) * self%ipf_wd + ix 
      qu = Orientations%getQuatfromArray(iq)
      call qu%quat_pos()
      if (OPC.eqv..TRUE.) then 
        IPFmap(1:3, ix, iy) = self%get_ipf_RGB_( sDir, qu, sym, Pm, clr )
      else 
        IPFmap(1:3, ix, iy) = self%get_ipf_RGB_( sDir, qu, sym, Pm )
      end if 
    end do 
  end do
!$OMP END DO

!$OMP END PARALLEL

! finally, store the IPFmap in a color TIFF file.
fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(self%get_ipf_filename_())
TIFF_filename = trim(fname)

! allocate memory for a color image; each pixel has 3 bytes [RGB]
allocate(TIFF_image(3*self%ipf_wd,self%ipf_ht))
TIFF_image = reshape( IPFmap, (/ 3*self%ipf_wd, self%ipf_ht /) )

! set up the image_t structure
im = image_t(TIFF_image)
im%dims = (/ self%ipf_wd, self%ipf_ht /)
im%samplesPerPixel = 3
if(im%empty()) call Message%printMessage("EMdpmerge: failed to convert array to rgb image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" Failed to write image to file : "//iomsg)
else
  call Message%printMessage(' IPF map written to '//trim(TIFF_filename),"(/A)")
end if

deallocate(TIFF_image, IPFmap)

end subroutine get_IPFMap_

end module mod_IPFsupport