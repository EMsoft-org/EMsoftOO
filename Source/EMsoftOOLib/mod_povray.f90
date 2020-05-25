! ###################################################################
! Copyright (c) 2016-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_povray
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/21/20
  !!
  !! All povray-related routines
  !!
  !! Routines to generate output for PoVray file creation; this is mostly
  !! used for the visualization of orientation data sets in one of the many
  !! representations and fundamental zones...
  !!
  !! The idea is that the program using this module will generate a skeleton 
  !! PoVRay file that can then be customized by the user.

use mod_kinds
use mod_global

IMPLICIT NONE 

private

type, public :: PoVRay_T 
  private

    character(fnlen)    :: filename = 'undefined'
    character(fnlen)    :: nmlfile = 'undefined' 
    character(fnlen)    :: locationline    ! default < 1.0, 0.0, 0.0 >
    character(fnlen)    :: skyline         ! default < 0.0, 0.0, 1.0>
    character(fnlen)    :: lightline       ! default <1, 2, -2>*50
    real(kind=sgl)      :: eyepos(3)
    logical             :: background
    integer(kind=irg)   :: dunit = 0       ! default value = 90
    integer(kind=irg)   :: nmlunit = 88    ! default value 
    integer(kind=irg)   :: df3unit = 86    ! default value 

  contains
  private

    procedure, pass(self) :: openFile_
    procedure, pass(self) :: setCamera_
    procedure, pass(self) :: setLightSource_
    procedure, pass(self) :: addEulerBox_
    procedure, pass(self) :: declare_DF3file_
    procedure, pass(self) :: write_DF3file_
    procedure, pass(self) :: addWireFrameSphere_
    procedure, pass(self) :: addReferenceFrame_
    procedure, pass(self) :: addSphere_
    procedure, pass(self) :: addCylinder_
    procedure, pass(self) :: addCubochoricCube_
    procedure, pass(self) :: getpos_FZ432_
    procedure, pass(self) :: getpos_FZ23_
    procedure, pass(self) :: getpos_FZ622_
    procedure, pass(self) :: getpos_FZ422_
    procedure, pass(self) :: getpos_FZ32_
    procedure, pass(self) :: getpos_FZ222_
    procedure, pass(self) :: drawFZ_
    procedure, pass(self) :: initFZCyclic_
    procedure, pass(self) :: fliprotationmatrix_
    procedure, pass(self) :: closeFile_
    final :: PoVRay_destructor

! general utility procedures
    generic, public :: openFile => openFile_
    generic, public :: setCamera => setCamera_
    generic, public :: setLightSource => setLightSource_
! routines for handling DF3 volume files 
    generic, public :: declare_DF3file => declare_DF3file_
    generic, public :: write_DF3file => write_DF3file_
! routines for drawing selected primitives and other scene objects 
    generic, public :: addWireFrameSphere => addWireFrameSphere_
    generic, public :: addReferenceFrame => addReferenceFrame_
    generic, public :: addSphere => addSphere_
    generic, public :: addCylinder => addCylinder_
! routines for drawing fundamental zones
    generic, public :: addEulerBox => addEulerBox_
    generic, public :: addCubochoricCube => addCubochoricCube_
    generic, public :: getpos_FZ432 => getpos_FZ432_
    generic, public :: getpos_FZ23 => getpos_FZ23_
    generic, public :: getpos_FZ622 => getpos_FZ622_
    generic, public :: getpos_FZ422 => getpos_FZ422_
    generic, public :: getpos_FZ32 => getpos_FZ32_
    generic, public :: getpos_FZ222 => getpos_FZ222_
    generic, public :: drawFZ => drawFZ_
    generic, public :: initFZCyclic => initFZCyclic_
    generic, public :: closeFile => closeFile_
! PoVRay uses a left-handed reference frame, so we provide a conversion routine
! for the EMsoft right-handed convention
    generic, public :: fliprotationmatrix => fliprotationmatrix_

end type PoVRay_T

! the constructor routine for this class 
interface PoVRay_T
  module procedure PoVRay_constructor
end interface PoVRay_T

!DEC$ ATTRIBUTES DLLEXPORT :: openFile
!DEC$ ATTRIBUTES DLLEXPORT :: setCamera
!DEC$ ATTRIBUTES DLLEXPORT :: setLightSource
!DEC$ ATTRIBUTES DLLEXPORT :: addEulerBox
!DEC$ ATTRIBUTES DLLEXPORT :: declare_DF3file
!DEC$ ATTRIBUTES DLLEXPORT :: write_DF3file
!DEC$ ATTRIBUTES DLLEXPORT :: addWireFrameSphere
!DEC$ ATTRIBUTES DLLEXPORT :: addReferenceFrame
!DEC$ ATTRIBUTES DLLEXPORT :: addSphere
!DEC$ ATTRIBUTES DLLEXPORT :: addCylinder
!DEC$ ATTRIBUTES DLLEXPORT :: addCubochoricCube
!DEC$ ATTRIBUTES DLLEXPORT :: getpos_FZ432
!DEC$ ATTRIBUTES DLLEXPORT :: getpos_FZ23
!DEC$ ATTRIBUTES DLLEXPORT :: getpos_FZ622
!DEC$ ATTRIBUTES DLLEXPORT :: getpos_FZ422
!DEC$ ATTRIBUTES DLLEXPORT :: getpos_FZ32
!DEC$ ATTRIBUTES DLLEXPORT :: getpos_FZ222
!DEC$ ATTRIBUTES DLLEXPORT :: drawFZ
!DEC$ ATTRIBUTES DLLEXPORT :: initFZCyclic
!DEC$ ATTRIBUTES DLLEXPORT :: closeFile
!DEC$ ATTRIBUTES DLLEXPORT :: PoVRay_fliprotationmatrix

contains 

!--------------------------------------------------------------------------
type(PoVRay_T) function PoVRay_constructor( EMsoft, fname, dunit, nmlfile, locationline, lightline, skyline) result(PV)
!! author: MDG 
!! version: 1.0 
!! date: 01/21/20
!!
!! constructor for the PoVRay_T Class 
 
use mod_EMsoft 

IMPLICIT NONE

type(EMsoft_T), INTENT(INOUT)          :: EMsoft
 !! used to set file path etc...
character(fnlen), INTENT(IN)           :: fname   
 !! PoVRay file name (full path)
integer(kind=irg),INTENT(IN), OPTIONAL :: dunit 
 !! optional output unit number (default = 90)
character(fnlen), INTENT(IN), OPTIONAL :: nmlfile
 !! full path to an nml file to be included in the top comments
character(fnlen), INTENT(IN), OPTIONAL :: locationline
character(fnlen), INTENT(IN), OPTIONAL :: lightline
 !! position of first light source (default <1, 2, -2>*50)
character(fnlen), INTENT(IN), OPTIONAL :: skyline
 !! position of first light source (default <1, 2, -2>*50)

if (present(dunit)) then 
  PV%dunit = dunit 
else 
  PV%dunit = 90 
end if 

if (present(skyline)) then 
  PV%skyline = trim(skyline)
else 
  PV%skyline = 'sky < 0.0, 0.0, 1.0>'
end if

if (present(locationline)) then 
  PV%locationline = trim(locationline)
else 
  PV%locationline = "location < 1.0, 0.0, 0.0 >"
end if 

if (present(lightline)) then 
  PV%lightline = trim(lightline)
else 
  PV%lightline = "<1, 2, -2>*50"
end if 

PV%filename = trim(fname) 

if (present(nmlfile)) then
  PV%nmlfile = trim(nmlfile)
else 
  PV%nmlfile = 'undefined'
end if 

call PV%openFile(EMsoft)
call PV%setCamera()
call PV%setLightSource()

! from here on, the output file is ready to receive user scene commands

end function PoVRay_constructor

!--------------------------------------------------------------------------
subroutine PoVRay_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/02/20
!!
!! destructor for the PoVRay_T Class

use mod_io 

IMPLICIT NONE

type(PoVRay_T), INTENT(INOUT)   :: self 

type(IO_T)                      :: Message
logical                         :: itsopen 

call reportDestructor('PoVRay_T')

! if the file unit is still open, close it here 
if (self%dunit.ne.0) then 
  inquire(unit=self%dunit, opened=itsopen)

  if (itsopen.eqv..TRUE.) then 
    close(unit=self%dunit, status='keep')
    call Message%printMessage(' Closed PoVray file '//trim(self%filename))
  end if 
end if 

end subroutine PoVRay_destructor

!--------------------------------------------------------------------------
recursive subroutine closeFile_(self)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! clean up routine 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 

close(unit=self%dunit, status = 'keep')

end subroutine closeFile_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! general utilities to open/close file and add code blocks for standard objects
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine openFile_(self, EMsoft)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! open a PoVRay file and add some information as comments

use mod_EMsoft 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 
type(EMsoft_T),INTENT(INOUT)          :: EMsoft 

character(fnlen)                      :: fname, line, cwd
integer(kind=irg)                     :: io

open(unit=self%dunit,file=trim(self%filename),form='formatted',status='unknown')

write (self%dunit,"(A)") "//Persistence of Vision Ray Tracer Scene Description File"
write (self%dunit,"(A)") "//Created by EMsoft package"

! should we copy the complete namelist file as a comment at this point ?
if (trim(self%nmlfile).ne.'undefined') then
  call getcwd(cwd)
  fname = trim(cwd)//'/'//trim(self%nmlfile)
  fname = EMsoft%toNativePath(fname)
  write(self%dunit,"(A)") "// "
  write(self%dunit,"(A)") "// Contents of the namelist used to create this file "//trim(fname)
  open(unit=self%nmlunit,file=trim(fname),status='unknown',form='formatted')
  do
    read(self%nmlunit,"(A)",iostat=io) line
    if (io.eq.0) then
      write(self%dunit,"(A)") "// "//trim(line)
    else 
      EXIT 
    end if
  end do
  close(unit=self%nmlunit,status='keep')
  write(self%dunit,"(A)") "// "
end if 

! and write the include statements to the file
write (self%dunit,"(A)") "#include ""colors.inc"""
write (self%dunit,"(A)") "#include ""textures.inc"""
write (self%dunit,"(A)") "#include ""glass.inc"""
write (self%dunit,"(A)") "// "
write (self%dunit,"(A)") "global_settings"
write (self%dunit,"(A)") "{  ambient_light <1,1,1>"
write (self%dunit,"(A)") "   assumed_gamma 1"
write (self%dunit,"(A)") "}"

end subroutine openFile_

!--------------------------------------------------------------------------
recursive subroutine setCamera_(self)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add the camera command to the current PoVRay file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 

write (self%dunit,"(A)") " "
write (self%dunit,"(A)") "camera {"
write (self%dunit,"(A)") "perspective "
write (self%dunit,"(A)") trim(self%locationline)
write (self%dunit,"(A)") trim(self%skyline)
write (self%dunit,"(A)") "right y * 1"
write (self%dunit,"(A)") "up z"
write (self%dunit,"(A)") "angle 50"
write (self%dunit,"(A)") "look_at < 0.0, 0.0, 0.0>"
write (self%dunit,"(A)") "}"


end subroutine setCamera_

!--------------------------------------------------------------------------
recursive subroutine setLightSource_(self)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add a light command with default properties to the current PoVRay file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 

write (self%dunit,"(A)") " "
write (self%dunit,"(A)") "light_source {"
write (self%dunit,"(A)") trim(self%lightline)
write (self%dunit,"(A)") "color White"
write (self%dunit,"(A)") "media_interaction on"
write (self%dunit,"(A)") "media_attenuation on"
write (self%dunit,"(A)") "shadowless"
write (self%dunit,"(A)") "}"
if (self%background) then
  write (self%dunit,"(A)") "background { color White }"
  write (self%dunit,"(A)") " "
end if

end subroutine setLightSource_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! routines to create and handle DF3 volume files 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine declare_DF3file_(self, df3name, levelset)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! set up the PoVRay code to include a 3D density field (df3) file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 

character(fnlen),INTENT(IN)     :: df3name
 !! DF3 output file name  (full path)
logical,OPTIONAL,INTENT(IN)     :: levelset
 !! if TRUE, then different colors are used for level contours

logical                         :: levels

write (self%dunit,"(A)") ''
write (self%dunit,"(A)") '#declare DENS = interior'
write (self%dunit,"(A)") '  {  media'
write (self%dunit,"(A)") '    {  intervals 100                   // number of ray-intervals'
write (self%dunit,"(A)") '      ratio 0.5'
write (self%dunit,"(A)") '      samples 4,4                     // maximum,minimum of samples per voxel'
write (self%dunit,"(A)") '      method 2                        // 1, 2 or 3 (3 requires samples 3,3)'
write (self%dunit,"(A)") '      emission 20*<1,1,1>'
write (self%dunit,"(A)") '      absorption <1,1,1>'
write (self%dunit,"(A)") '      scattering { 1, <0,0,0> }'
write (self%dunit,"(A)") '      confidence 0.9                // default: 0.9'
write (self%dunit,"(A)") '      variance 1/100                // default: 1/128'
write (self%dunit,"(A)") '      density'
     write (self%dunit,"('      {  density_file df3 ""',A,'""')") trim(df3name)
write (self%dunit,"(A)") '        interpolate 1'
write (self%dunit,"(A)") '        color_map                    // colour map with (smooth) linear transition(s)'

levels = .FALSE.
if (PRESENT(levelset)) then
  if (levelset.eqv..TRUE.) then
    levels = .TRUE.
  end if
end if

! these are experimental level values; the user may want to modifiy those in the 
! output PoVRay file ...

if (levels.eqv..TRUE.) then 
  write (self%dunit,"(A)") '        {  [0.0 rgb <0.0,0.0,0.0>] '
  write (self%dunit,"(A)") '           [0.1 rgb <0.6,0.0,0.0>]'
  write (self%dunit,"(A)") '           [0.2 rgb <0.0,0.6,0.0>]'
  write (self%dunit,"(A)") '           [0.4 rgb <0.2,0.3,0.8>]'
  write (self%dunit,"(A)") '           [0.6 rgb <0.6,0.5,0.8>]'
  write (self%dunit,"(A)") '           [1.0 rgb <1.0,1.0,1.0>]'
else
  write (self%dunit,"(A)") '        {  [0.1 rgb <0.0,0.0,0.0>] '
  write (self%dunit,"(A)") '           [0.2 rgb <0.5,0.0,0.0>]'
  write (self%dunit,"(A)") '           [0.3 rgb <1.0,0.0,0.0>]'
  write (self%dunit,"(A)") '           [0.4 rgb <1.0,0.3,0.0>]'
  write (self%dunit,"(A)") '           [0.5 rgb <1.0,0.6,0.0>]'
  write (self%dunit,"(A)") '           [0.6 rgb <1.0,1.0,0.0>]'
  write (self%dunit,"(A)") '           [0.7 rgb <1.0,1.0,0.3>]'
  write (self%dunit,"(A)") '           [0.8 rgb <1.0,1.0,0.6>]'
  write (self%dunit,"(A)") '           [0.9 rgb <1.0,1.0,1.0>]'
end if

write (self%dunit,"(A)") '        }'
write (self%dunit,"(A)") '      }'
write (self%dunit,"(A)") '    }'
write (self%dunit,"(A)") '  }'
write (self%dunit,"(A)") ''
write (self%dunit,"(A)") '#declare renderbox = box'
write (self%dunit,"(A)") '  {  <0,0,0>, <1,1,1>'
write (self%dunit,"(A)") '     pigment { rgbt <0,0,0,1> }'
write (self%dunit,"(A)") '     hollow'
write (self%dunit,"(A)") '     interior { DENS }'
write (self%dunit,"(A)") '  }'

end subroutine declare_DF3file_

!--------------------------------------------------------------------------
recursive subroutine write_DF3file_(self, df3name, volume, ndims, scalingmode)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! set up the PoVRay code to include a 3D density field (df3) file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 
character(fnlen),INTENT(IN)     :: df3name
 !! output file name 
integer(kind=irg),INTENT(IN)    :: ndims(3)
 !! dimensions of the volume array 
real(kind=sgl),INTENT(INOUT)    :: volume(-ndims(1):ndims(1),-ndims(2):ndims(2),-ndims(3):ndims(3))
 !! volume array to be written to DF3 file 
character(3),INTENT(IN)         :: scalingmode
 !! scaling type: 'lin' or 'log' or 'lev'

integer(kind=ish)                        :: ivol(-ndims(1):ndims(1),-ndims(2):ndims(2),-ndims(3):ndims(3))
integer(kind=ish)                          :: idims(3), mval
real(kind=sgl)                              :: mi, ma, levels(6)
integer(kind=irg)                          :: recno, i, j, k

! This format is described on the following web pages:
! http://wwwmpa.mpa-garching.mpg.de/~mselig/povray/povray.html
! http://www.povray.org/documentation/view/3.6.1/374/
!
! Essentially, this is a binary file with first 3 2-byte integers
! with the array dimensions, and then the array itself as 4-byte integers.
! Note that this file *must* be written in big-endian format (most significant
! byte first).  This can be achieved by using the CONVERT='BIG_ENDIAN' qualifier
! when the file is opened.

! we will write this as a file of 16 bit integers (kind=ish), so we need to rescale them all
! to the correct range
mval = int(2**15-1,kind=ish)
if ((scalingmode.eq.'log').or.(scalingmode.eq.'lev')) volume = alog10(volume+1.0)
mi = minval(volume)
ma = maxval(volume)
volume = (volume - mi) / (ma - mi)

if (scalingmode.eq.'lev') then
! we discretize the logarithmic values into bins (levels) which are then displayed
! as nested contour surfaces (ideally) by PoVRay's Density_field routine...
! Since we are testing this capability, we'll only include a few different contour levels
! for now; this can be modified later on...
  levels = (/ 0.0, 0.4, 0.6, 0.8, 0.9, 1.0 /)
  do i=2,6
    where ((volume.ge.levels(i-1)).and.(volume.lt.levels(i))) volume = levels(i-1)
  end do
endif

volume = volume * mval
ivol = int(volume, kind=ish)

! open the file and write the rescaled integer data
open(unit=self%df3unit,file=trim(df3name),status='unknown',access="DIRECT",action="WRITE", &
     recl=2,form='unformatted',convert='big_endian')
idims = 2*ndims+1
recno = 1
do i=1,3
  write (self%df3unit,rec=recno) idims(i)
  recno = recno+1
end do
do k=-ndims(3),ndims(3)
  do j=-ndims(2),ndims(2)
    do i=-ndims(1),ndims(1)
      write (self%df3unit,rec=recno) ivol(i,j,k)
      recno = recno+1
    end do
  end do
end do
close(self%df3unit,status='keep')

end subroutine write_DF3file_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! routines to create a few standard scene objects (there are many more in PoVRay...)
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine addWireFrameSphere_(self, sphereRadius)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add a wireframe sphere to the current PoVRay file

use mod_io 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 
real(kind=dbl),INTENT(IN)       :: sphereRadius

type(IO_T)                      :: Message 
character(12)                   :: px

call Message%printMessage( (/ &
     "#macro WireFrameSphere(NrLongitudes, NrLatitudes, Rmaj, Rmin)", &
     "                                                             ", & 
     "  #local dLongitude = 360/NrLongitudes;                      ", & 
     "  #local dLatitude = 180/NrLatitudes;                        ", & 
     "  #local Cnt = 0;                                            ", & 
     "  #while (Cnt < NrLongitudes)                                ", & 
     "    #local Longitude = Cnt*dLongitude;                       ", & 
     "    difference {                                             ", & 
     "      torus { Rmaj, Rmin }                                   ", & 
     "      plane { -z, 0 }                                        ", & 
     "      rotate -90*z                                           ", & 
     "      rotate Longitude*y                                     ", & 
     "    }                                                        ", & 
     "    #local Cnt = Cnt + 1;                                    ", & 
     "  #end // while                                              ", & 
     "                                                             ", & 
     "  #local Cnt = 1;                                            ", & 
     "  #while (Cnt < NrLatitudes)                                 ", & 
     "    #local Latitude = radians(Cnt*dLatitude - 90);           ", & 
     "    torus {                                                  ", & 
     "      Rmaj*cos(Latitude), Rmin                               ", & 
     "      translate Rmaj*sin(Latitude)*y                         ", & 
     "    }                                                        ", & 
     "    #local Cnt = Cnt + 1;                                    ", & 
     "  #end // while                                              ", & 
     "                                                             ", & 
     "#end // macro WireFrameSphere                                ", & 
     "                                                             " /), redirect = self%dunit)

write (px,"(F12.6)") sphereRadius
write (self%dunit,"(A)") "#declare Rglobe = "//px//";"

call Message%printMessage( (/ &
     "#declare Rwireframe = 0.0033;                               ", &
     "#declare Rspheres = Rwireframe*2;                           ", &
     "                                                            ", &
     "// Number of longitude intervals                            ", &
     "#declare Longitudes = 12;                                   ", &
     "                                                            ", &
     "// Number of latitude intervals                             ", &
     "#declare Latitudes = 6;                                     ", &
     "                                                            ", &
     "                                                            ", &
     "union {                                                     ", &
     "  WireFrameSphere(Longitudes, Latitudes, Rglobe, Rwireframe)", &
     "  rotate 90.0*x                                             ", &
     "  pigment { color Black*0.7 }                               ", &
     "}                                                           " /), redirect = self%dunit)

end subroutine addWireFrameSphere_

!--------------------------------------------------------------------------
recursive subroutine addReferenceFrame_(self, ac, cylr)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add a reference frame to the current PoVRay file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 
real(kind=dbl),INTENT(IN)             :: ac
 !! maximum semi-length of axis
real(kind=dbl),OPTIONAL,INTENT(IN)    :: cylr
 !! optional cylinder radius 

real(kind=dbl)                        :: cc

if (present(cylr).eqv..TRUE.) then
  cc = cylr
else
  cc = 0.005D0
end if

write (self%dunit,"('cylinder { <',2(F9.6,','),F9.6,'>,<',2(F9.6,','),F9.6,'>,',F9.6,' pigment { color Red*0.7 } }')") &
      -ac, 0.0, 0.0,  ac, 0.0, 0.0, cc
write (self%dunit,"('cylinder { <',2(F9.6,','),F9.6,'>,<',2(F9.6,','),F9.6,'>,',F9.6,' pigment { color Green*0.7 } }')") &
      0.0, -ac, 0.0, 0.0,  ac, 0.0, cc
write (self%dunit,"('cylinder { <',2(F9.6,','),F9.6,'>,<',2(F9.6,','),F9.6,'>,',F9.6,' pigment { color Blue*0.7 } }')") &
      0.0, 0.0, -ac, 0.0, 0.0,  ac, cc

end subroutine addReferenceFrame_

!--------------------------------------------------------------------------
recursive subroutine addSphere_(self, ctr, radius, rgb)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add a sphere to the current PoVRay file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 
real(kind=dbl),INTENT(IN)             :: ctr(3)
 !! center coordinates of the sphere
real(kind=dbl),INTENT(IN)             :: radius 
 !! sphere radius
real(kind=sgl),INTENT(IN)             :: rgb(3)
 !! rgb color triplet

write (self%dunit,"('sphere { <',2(F9.6,','),F9.6,'>,',F9.6,' material { texture { pigment { rgb <', &
                   &2(F9.6,','),F9.6,'>}}}}')") ctr(1:3), radius, rgb(1:3)

end subroutine addSphere_

!--------------------------------------------------------------------------
recursive subroutine addCylinder_(self, p1, p2, radius, rgb)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add a cylinder to the current PoVRay file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 
real(kind=dbl),INTENT(IN)             :: p1(3)
 !! starting point of cylinder (on axis)
real(kind=dbl),INTENT(IN)             :: p2(3)
 !! end point of cylinder (on axis)
real(kind=dbl),INTENT(IN)             :: radius 
 !! cylinder radius
real(kind=sgl),INTENT(IN)             :: rgb(3)
 !! color triplet (RGB)

write (self%dunit,"('cylinder { <',2(F9.6,','),F9.6,'>,<',2(F9.6,','),F9.6,'>,', F9.6,' pigment { ', &
                  &'rgb <',2(F9.6,','),F9.6,'>}}')") p1(1:3), p2(1:3), radius, rgb(1:3)

end subroutine addCylinder_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! code to draw a variety of fundamental zones 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine addEulerBox_(self)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add a wireframe Euler Box to the current PoVRay file

use mod_io 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 

type(IO_T)                      :: Message 

call Message%printMessage( (/ &
  "cylinder {<-3.141593,-1.570796,-3.141593>,<-3.141593, 1.570796,-3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {<-3.141593,-1.570796, 3.141593>,<-3.141593, 1.570796, 3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {< 3.141593,-1.570796,-3.141593>,< 3.141593, 1.570796,-3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {< 3.141593,-1.570796, 3.141593>,< 3.141593, 1.570796, 3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {<-3.141593,-1.570796,-3.141593>,<-3.141593,-1.570796, 3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {<-3.141593,-1.570796, 3.141593>,< 3.141593,-1.570796, 3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {< 3.141593,-1.570796, 3.141593>,< 3.141593,-1.570796,-3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {< 3.141593,-1.570796,-3.141593>,<-3.141593,-1.570796,-3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {<-3.141593, 1.570796,-3.141593>,<-3.141593, 1.570796, 3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {<-3.141593, 1.570796, 3.141593>,< 3.141593, 1.570796, 3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {< 3.141593, 1.570796, 3.141593>,< 3.141593, 1.570796,-3.141593>, 0.005 pigment {color Green*0.7}}", &
  "cylinder {< 3.141593, 1.570796,-3.141593>,<-3.141593, 1.570796,-3.141593>, 0.005 pigment {color Green*0.7}}"/), &
  redirect = self%dunit)

end subroutine addEulerBox_

!--------------------------------------------------------------------------
recursive subroutine addCubochoricCube_(self)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! add a wireframe cubochoric Box to the current PoVRay file

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 

real(kind=dbl)                  :: ac

! create the cubochoric cube
ac = 0.5 * LPs%ap
call self%addCylinder((/ -ac,  ac,  ac /), (/ ac,  ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/ -ac,  ac, -ac /), (/ ac,  ac, -ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/ -ac, -ac,  ac /), (/ ac, -ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/ -ac, -ac, -ac /), (/ ac, -ac, -ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))

call self%addCylinder((/  ac, -ac,  ac /), (/ ac,  ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/ -ac, -ac,  ac /), (/-ac,  ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/  ac, -ac, -ac /), (/ ac,  ac, -ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/ -ac, -ac, -ac /), (/-ac,  ac, -ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))

call self%addCylinder((/  ac,  ac, -ac /), (/ ac,  ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/ -ac,  ac, -ac /), (/-ac,  ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/  ac, -ac, -ac /), (/ ac, -ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))
call self%addCylinder((/ -ac, -ac, -ac /), (/-ac, -ac,  ac /), 0.005D0, (/ 0.7, 0.7, 0.7 /))

end subroutine addCubochoricCube_

!--------------------------------------------------------------------------
recursive function fliprotationmatrix_(self, M) result(O)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! change a rotation matrix to the POVray axes convention (left handed)

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)   :: self 

real(kind=dbl),INTENT(IN)       :: M(3,3)
real(kind=dbl)                  :: O(3,3)

O = reshape( (/ M(1,1), M(1,3), M(1,2), M(3,1), M(3,3), M(3,2), M(2,1), M(2,3), M(2,2) /), (/ 3,3 /) )

end function fliprotationmatrix_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! below this line are a number of routines for different rotation representations
! and different crystallographic fundamental zones... below this line, there
! should be no direct writing to the dunit file; all writes should pass 
! through the routines above.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine getpos_FZ432_(self, dims, cpos, s_edge, t_edge, ns, d, nt, MFZ)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! initialize the PoVRay output for rotational group 432 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 

integer(kind=irg),INTENT(IN)          :: dims(3)
 !! array dimensions
real(kind=dbl),INTENT(INOUT)          :: cpos(3,dims(1))
 !! vertex coordinates
integer(kind=irg),INTENT(INOUT)       :: s_edge(2,dims(2))
 !! first set of edge connectivities
integer(kind=irg),INTENT(INOUT)       :: t_edge(2,dims(3))
 !! second set of edge connectivities
integer(kind=irg),INTENT(OUT)         :: ns
 !! aux parameter
real(kind=dbl),INTENT(OUT)            :: d
 !! aux parameter
integer(kind=irg),INTENT(OUT)         :: nt
 !! aux parameter
logical,OPTIONAL,INTENT(IN)           :: MFZ
 !! (optional) return coordinates for Mackenzie FZ instead of regular FZ

real(kind=dbl)  :: a = 0.41421356237D0, c = 0.17157287525D0, dt = 0.34314575050D0, ds = 0.6340506711D0, &
                   dd, e = 0.29289323D0, f = 0.333333333D0

d = 0.610395774912

if (present(MFZ)) then
  if (MFZ) then 
    d = 1.0
! define the coordinates of the cubic Mackenzie FZ in Rodrigues Space
    cpos(1:3, 1) = (/  0.D0,  0.D0,  0.D0 /)
    cpos(1:3, 2) = (/  a,  0.D0,  0.D0 /)
    cpos(1:3, 3) = (/  a,  a,  0.D0 /)

    cpos(1:3, 4) = (/  a,  a,  c /)
    cpos(1:3, 5) = (/  f,  f,  f /)
    cpos(1:3, 6) = (/  a,  e,  e /)

    ! cpos = cpos/d

    ns = 200
    nt = nint( ns * dt/ds )

! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1, 2 /)
    s_edge(1:2, 2) = (/  1, 3 /)
    s_edge(1:2, 3) = (/  1, 5 /)
    s_edge(1:2, 4) = (/  2, 3 /)
    s_edge(1:2, 5) = (/  2, 6 /)
    s_edge(1:2, 6) = (/  3, 4 /)
    s_edge(1:2, 7) = (/  4, 5 /)
    s_edge(1:2, 8) = (/  4, 6 /)
    s_edge(1:2, 9) = (/  5, 6 /)
  end if 
else ! define the coordinates of the cubic FZ in Rodrigues Space
    cpos(1:3, 1) = (/  a,  a,  c /)
    cpos(1:3, 2) = (/  c,  a,  a /)
    cpos(1:3, 3) = (/  a,  c,  a /)

    cpos(1:3, 4) = (/ -a,  a,  c /)
    cpos(1:3, 5) = (/ -c,  a,  a /)
    cpos(1:3, 6) = (/ -a,  c,  a /)

    cpos(1:3, 7) = (/ -a, -a,  c /)
    cpos(1:3, 8) = (/ -c, -a,  a /)
    cpos(1:3, 9) = (/ -a, -c,  a /)

    cpos(1:3,10) = (/  a, -a,  c /)
    cpos(1:3,11) = (/  c, -a,  a /)
    cpos(1:3,12) = (/  a, -c,  a /)

    cpos(1:3,13) = (/  a,  a, -c /)
    cpos(1:3,14) = (/  a,  c, -a /)
    cpos(1:3,15) = (/  c,  a, -a /)

    cpos(1:3,16) = (/ -a,  a, -c /)
    cpos(1:3,17) = (/ -a,  c, -a /)
    cpos(1:3,18) = (/ -c,  a, -a /)

    cpos(1:3,19) = (/ -a, -a, -c /)
    cpos(1:3,20) = (/ -a, -c, -a /)
    cpos(1:3,21) = (/ -c, -a, -a /)

    cpos(1:3,22) = (/  a, -a, -c /)
    cpos(1:3,23) = (/  a, -c, -a /)
    cpos(1:3,24) = (/  c, -a, -a /)

    cpos = cpos/d

    ns = 200
    nt = nint( ns * dt/ds )

    ! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  3, 12 /)
    s_edge(1:2, 2) = (/ 10, 22 /)
    s_edge(1:2, 3) = (/ 14, 23 /)
    s_edge(1:2, 4) = (/  1, 13 /)
    s_edge(1:2, 5) = (/  2,  5 /)
    s_edge(1:2, 6) = (/  8, 11 /)
    s_edge(1:2, 7) = (/ 21, 24 /)
    s_edge(1:2, 8) = (/ 15, 18 /)
    s_edge(1:2, 9) = (/  6,  9 /)
    s_edge(1:2,10) = (/  7, 19 /)
    s_edge(1:2,11) = (/ 17, 20 /)
    s_edge(1:2,12) = (/  4, 16 /)

    t_edge(1:2, 1) = (/  1,  2 /)
    t_edge(1:2, 2) = (/  2,  3 /)
    t_edge(1:2, 3) = (/  3,  1 /)

    t_edge(1:2, 4) = (/  4,  5 /)
    t_edge(1:2, 5) = (/  5,  6 /)
    t_edge(1:2, 6) = (/  6,  4 /)

    t_edge(1:2, 7) = (/  7,  8 /)
    t_edge(1:2, 8) = (/  8,  9 /)
    t_edge(1:2, 9) = (/  9,  7 /)

    t_edge(1:2,10) = (/ 10, 11 /)
    t_edge(1:2,11) = (/ 11, 12 /)
    t_edge(1:2,12) = (/ 12, 10 /)

    t_edge(1:2,13) = (/ 13, 14 /)
    t_edge(1:2,14) = (/ 14, 15 /)
    t_edge(1:2,15) = (/ 15, 13 /)

    t_edge(1:2,16) = (/ 16, 17 /)
    t_edge(1:2,17) = (/ 17, 18 /)
    t_edge(1:2,18) = (/ 18, 16 /)

    t_edge(1:2,19) = (/ 19, 20 /)
    t_edge(1:2,20) = (/ 20, 21 /)
    t_edge(1:2,21) = (/ 21, 19 /)

    t_edge(1:2,22) = (/ 22, 23 /)
    t_edge(1:2,23) = (/ 23, 24 /)
    t_edge(1:2,24) = (/ 24, 22 /)
end if

end subroutine getpos_FZ432_

!--------------------------------------------------------------------------
recursive subroutine getpos_FZ23_(self, dims, cpos, s_edge, t_edge, ns, d, nt, MFZ)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! initialize the PoVRay output for rotational group 23 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 

integer(kind=irg),INTENT(IN)          :: dims(3)
 !! array dimensions
real(kind=dbl),INTENT(INOUT)          :: cpos(3,dims(1))
 !! vertex coordinates
integer(kind=irg),INTENT(INOUT)       :: s_edge(2,dims(2))
 !! first set of edge connectivities
integer(kind=irg),INTENT(INOUT)       :: t_edge(2,dims(3))
 !! second set of edge connectivities
integer(kind=irg),INTENT(OUT)         :: ns
 !! aux parameter
real(kind=dbl),INTENT(OUT)            :: d
 !! aux parameter
integer(kind=irg),INTENT(OUT)         :: nt
 !! aux parameter
logical,OPTIONAL,INTENT(IN)           :: MFZ
 !! (optional) return coordinates for Mackenzie FZ instead of regular FZ

real(kind=dbl)  :: a = 1.D0, b = 0.0D0, c = 0.5773502692D0, e = 0.333333333D0, &
                   ds = 0.6340506711D0, dt = 1.4142135623730D0, dd, zz = 0.D0, oo = 1.D0

d = 1.0

if (present(MFZ)) then
  if (MFZ) then 
    d = 1.0

! define the coordinates of the cubic Mackenzie FZ in Rodrigues Space
    cpos(1:3, 1) = (/  b,  b,  b /)
    cpos(1:3, 2) = (/  a,  b,  b /)
    cpos(1:3, 3) = (/  b,  a,  b /)
    cpos(1:3, 4) = (/  e,  e,  e /)

    cpos = cpos/d

    ns = 200
    nt = nint( ns * dt/ds )

! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1, 2 /)
    s_edge(1:2, 2) = (/  1, 3 /)
    s_edge(1:2, 3) = (/  1, 4 /)
    s_edge(1:2, 4) = (/  2, 3 /)
    s_edge(1:2, 5) = (/  2, 4 /)
    s_edge(1:2, 6) = (/  3, 4 /)
  end if 
else ! define the coordinates of the cubic FZ in Rodrigues Space
    cpos(1:3, 1) = (/  a,  b,  b /)
    cpos(1:3, 2) = (/  b,  a,  b /)
    cpos(1:3, 3) = (/ -a,  b,  b /)
    cpos(1:3, 4) = (/  b, -a,  b /)
    cpos(1:3, 5) = (/  b,  b,  a /)
    cpos(1:3, 6) = (/  b,  b, -a /)

    cpos = cpos / d

    ns = 200
    nt = nint( ns * dt/ds )

    ! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1,  2 /)
    s_edge(1:2, 2) = (/  2,  3 /)
    s_edge(1:2, 3) = (/  3,  4 /)
    s_edge(1:2, 4) = (/  4,  1 /)
    s_edge(1:2, 5) = (/  1,  5 /)
    s_edge(1:2, 6) = (/  2,  5 /)
    s_edge(1:2, 7) = (/  3,  5 /)
    s_edge(1:2, 8) = (/  4,  5 /)
    s_edge(1:2, 9) = (/  1,  6 /)
    s_edge(1:2,10) = (/  2,  6 /)
    s_edge(1:2,11) = (/  3,  6 /)
    s_edge(1:2,12) = (/  4,  6 /)
 
end if

end subroutine getpos_FZ23_

!--------------------------------------------------------------------------
recursive subroutine getpos_FZ622_(self, dims, cpos, s_edge, t_edge, ns, d, nt, MFZ)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! initialize the PoVRay output for rotational group 622 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 

integer(kind=irg),INTENT(IN)          :: dims(3)
 !! array dimensions
real(kind=dbl),INTENT(INOUT)          :: cpos(3,dims(1))
 !! vertex coordinates
integer(kind=irg),INTENT(INOUT)       :: s_edge(2,dims(2))
 !! first set of edge connectivities
integer(kind=irg),INTENT(INOUT)       :: t_edge(2,dims(3))
 !! second set of edge connectivities
integer(kind=irg),INTENT(OUT)         :: ns
 !! aux parameter
real(kind=dbl),INTENT(OUT)            :: d
 !! aux parameter
integer(kind=irg),INTENT(OUT)         :: nt
 !! aux parameter
logical,OPTIONAL,INTENT(IN)           :: MFZ
 !! (optional) return coordinates for Mackenzie FZ instead of regular FZ

real(kind=dbl)  :: a = 1.0D0, b = 0.267949192431D0, c = 0.732050807569D0, &
                   dt = 0.5358983848622454D0, ds = 0.5358983848622454D0, di =1.069389330154823D0, dd, &
                   z = 0.D0, o = 0.86602540378443D0, p = 0.5D0


d = 1.0693893290743279D0
if (present(MFZ)) then
  if (MFZ) then ! define the coordinates of the hexagonal Mackenzie FZ in Rodrigues Space
    d = 1.0
    cpos(1:3, 1) = (/  z,  z,  z /)
    cpos(1:3, 2) = (/  z,  z,  b /)
    cpos(1:3, 3) = (/  a,  z,  z /)
    cpos(1:3, 4) = (/  a,  z,  b /)
    cpos(1:3, 5) = (/  a,  b,  z /)
    cpos(1:3, 6) = (/  a,  b,  b /)
    cpos(1:3, 7) = (/  o,  p,  z /)
    cpos(1:3, 8) = (/  o,  p,  b /)

    cpos = cpos/d

    ns = 200
    nt = nint( ns * dt/ds )

! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1, 2 /)
    s_edge(1:2, 2) = (/  1, 3 /)
    s_edge(1:2, 3) = (/  1, 7 /)
    s_edge(1:2, 4) = (/  2, 4 /)
    s_edge(1:2, 5) = (/  2, 8 /)
    s_edge(1:2, 6) = (/  3, 4 /)
    s_edge(1:2, 7) = (/  3, 5 /)
    s_edge(1:2, 8) = (/  4, 6 /)
    s_edge(1:2, 9) = (/  5, 6 /)
    s_edge(1:2, 10) = (/  5, 7 /)
    s_edge(1:2, 11) = (/  6, 8 /)
    s_edge(1:2, 12) = (/  7, 8 /)
  end if 

else ! define the coordinates of the hexagonal FZ in Rodrigues Space
    cpos(1:3, 1) = (/  a,  b,  b /)
    cpos(1:3, 2) = (/  c,  c,  b /)
    cpos(1:3, 3) = (/  b,  a,  b /)

    cpos(1:3, 4) = (/ -b,  a,  b /)
    cpos(1:3, 5) = (/ -c,  c,  b /)
    cpos(1:3, 6) = (/ -a,  b,  b /)

    cpos(1:3, 7) = (/ -a, -b,  b /)
    cpos(1:3, 8) = (/ -c, -c,  b /)
    cpos(1:3, 9) = (/ -b, -a,  b /)

    cpos(1:3,10) = (/  b, -a,  b /)
    cpos(1:3,11) = (/  c, -c,  b /)
    cpos(1:3,12) = (/  a, -b,  b /)

    cpos(1:3,13) = (/  a,  b, -b /)
    cpos(1:3,14) = (/  c,  c, -b /)
    cpos(1:3,15) = (/  b,  a, -b /)

    cpos(1:3,16) = (/ -b,  a, -b /)
    cpos(1:3,17) = (/ -c,  c, -b /)
    cpos(1:3,18) = (/ -a,  b, -b /)

    cpos(1:3,19) = (/ -a, -b, -b /)
    cpos(1:3,20) = (/ -c, -c, -b /)
    cpos(1:3,21) = (/ -b, -a, -b /)

    cpos(1:3,22) = (/  b, -a, -b /)
    cpos(1:3,23) = (/  c, -c, -b /)
    cpos(1:3,24) = (/  a, -b, -b /)

    cpos = cpos / d

    ns = 200
    nt = nint( ns * dt/ds )

    ! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1,  2 /)
    s_edge(1:2, 2) = (/  2,  3 /)
    s_edge(1:2, 3) = (/  3,  4 /)
    s_edge(1:2, 4) = (/  4,  5 /)
    s_edge(1:2, 5) = (/  5,  6 /)
    s_edge(1:2, 6) = (/  6,  7 /)
    s_edge(1:2, 7) = (/  7,  8 /)
    s_edge(1:2, 8) = (/  8,  9 /)
    s_edge(1:2, 9) = (/  9, 10 /)
    s_edge(1:2,10) = (/ 10, 11 /)
    s_edge(1:2,11) = (/ 11, 12 /)
    s_edge(1:2,12) = (/ 12,  1 /)
    s_edge(1:2,13) = (/ 13, 14 /)
    s_edge(1:2,14) = (/ 14, 15 /)
    s_edge(1:2,15) = (/ 15, 16 /)
    s_edge(1:2,16) = (/ 16, 17 /)
    s_edge(1:2,17) = (/ 17, 18 /)
    s_edge(1:2,18) = (/ 18, 19 /)
    s_edge(1:2,19) = (/ 19, 20 /)
    s_edge(1:2,20) = (/ 20, 21 /)
    s_edge(1:2,21) = (/ 21, 22 /)
    s_edge(1:2,22) = (/ 22, 23 /)
    s_edge(1:2,23) = (/ 23, 24 /)
    s_edge(1:2,24) = (/ 24, 13 /)

    t_edge(1:2, 1) = (/  1, 13 /)
    t_edge(1:2, 2) = (/  2, 14 /)
    t_edge(1:2, 3) = (/  3, 15 /)
    t_edge(1:2, 4) = (/  4, 16 /)
    t_edge(1:2, 5) = (/  5, 17 /)
    t_edge(1:2, 6) = (/  6, 18 /)
    t_edge(1:2, 7) = (/  7, 19 /)
    t_edge(1:2, 8) = (/  8, 20 /)
    t_edge(1:2, 9) = (/  9, 21 /)
    t_edge(1:2,10) = (/ 10, 22 /)
    t_edge(1:2,11) = (/ 11, 23 /)
    t_edge(1:2,12) = (/ 12, 24 /)
end if

end subroutine getpos_FZ622_

!--------------------------------------------------------------------------
recursive subroutine getpos_FZ422_(self, dims, cpos, s_edge, t_edge, ns, d, nt, MFZ)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! initialize the PoVRay output for rotational group 622 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 

integer(kind=irg),INTENT(IN)          :: dims(3)
 !! array dimensions
real(kind=dbl),INTENT(INOUT)          :: cpos(3,dims(1))
 !! vertex coordinates
integer(kind=irg),INTENT(INOUT)       :: s_edge(2,dims(2))
 !! first set of edge connectivities
integer(kind=irg),INTENT(INOUT)       :: t_edge(2,dims(3))
 !! second set of edge connectivities
integer(kind=irg),INTENT(OUT)         :: ns
 !! aux parameter
real(kind=dbl),INTENT(OUT)            :: d
 !! aux parameter
integer(kind=irg),INTENT(OUT)         :: nt
 !! aux parameter
logical,OPTIONAL,INTENT(IN)           :: MFZ
 !! (optional) return coordinates for Mackenzie FZ instead of regular FZ

real(kind=dbl)    :: a = 1.0D0, b = 0.41421354D0, c = 0.41421354D0, dt = 0.8284270763397216D0, &
                     ds = 0.8284270763397216D0, dd, z = 0.D0, o = 0.70710678118654746D0

d = 1.158941651036677D0
if (present(MFZ)) then
  if (MFZ) then ! define the coordinates of the tetragonal Mackenzie FZ in Rodrigues Space
    d = 1.0
    cpos(1:3, 1) = (/  z,  z,  z /)
    cpos(1:3, 2) = (/  z,  z,  c /)
    cpos(1:3, 3) = (/  a,  z,  z /)
    cpos(1:3, 4) = (/  a,  z,  c /)
    cpos(1:3, 5) = (/  a,  b,  z /)
    cpos(1:3, 6) = (/  a,  b,  c /)
    cpos(1:3, 7) = (/  o,  o,  z /)
    cpos(1:3, 8) = (/  o,  o,  c /)

    cpos = cpos/d

    ns = 200
    nt = nint( ns * dt/ds )

! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1, 2 /)
    s_edge(1:2, 2) = (/  1, 3 /)
    s_edge(1:2, 3) = (/  1, 7 /)
    s_edge(1:2, 4) = (/  2, 4 /)
    s_edge(1:2, 5) = (/  2, 8 /)
    s_edge(1:2, 6) = (/  3, 4 /)
    s_edge(1:2, 7) = (/  3, 5 /)
    s_edge(1:2, 8) = (/  4, 6 /)
    s_edge(1:2, 9) = (/  5, 6 /)
    s_edge(1:2, 10) = (/  5, 7 /)
    s_edge(1:2, 11) = (/  6, 8 /)
    s_edge(1:2, 12) = (/  7, 8 /)
  end if 
else ! define the coordinates of the tetragonal 422 FZ in Rodrigues Space
    cpos(1:3, 1) = (/  a,  b,  c /)
    cpos(1:3, 2) = (/  b,  a,  c /)
    cpos(1:3, 3) = (/ -b,  a,  c /)
    cpos(1:3, 4) = (/ -a,  b,  c /)
    cpos(1:3, 5) = (/ -a, -b,  c /)
    cpos(1:3, 6) = (/ -b, -a,  c /)
    cpos(1:3, 7) = (/  b, -a,  c /)
    cpos(1:3, 8) = (/  a, -b,  c /)
    cpos(1:3, 9) = (/  a,  b, -c /)
    cpos(1:3,10) = (/  b,  a, -c /)
    cpos(1:3,11) = (/ -b,  a, -c /)
    cpos(1:3,12) = (/ -a,  b, -c /)
    cpos(1:3,13) = (/ -a, -b, -c /)
    cpos(1:3,14) = (/ -b, -a, -c /)
    cpos(1:3,15) = (/  b, -a, -c /)
    cpos(1:3,16) = (/  a, -b, -c /)

    cpos = cpos / d

    ns = 200
    nt = nint( ns * dt/ds )

! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1,  2 /)
    s_edge(1:2, 2) = (/  2,  3 /)
    s_edge(1:2, 3) = (/  3,  4 /)
    s_edge(1:2, 4) = (/  4,  5 /)
    s_edge(1:2, 5) = (/  5,  6 /)
    s_edge(1:2, 6) = (/  6,  7 /)
    s_edge(1:2, 7) = (/  7,  8 /)
    s_edge(1:2, 8) = (/  8,  1 /)
    s_edge(1:2, 9) = (/  9, 10 /)
    s_edge(1:2,10) = (/ 10, 11 /)
    s_edge(1:2,11) = (/ 11, 12 /)
    s_edge(1:2,12) = (/ 12, 13 /)
    s_edge(1:2,13) = (/ 13, 14 /)
    s_edge(1:2,14) = (/ 14, 15 /)
    s_edge(1:2,15) = (/ 15, 16 /)
    s_edge(1:2,16) = (/ 16,  9 /)

    t_edge(1:2, 1) = (/  1,  9 /)
    t_edge(1:2, 2) = (/  2, 10 /)
    t_edge(1:2, 3) = (/  3, 11 /)
    t_edge(1:2, 4) = (/  4, 12 /)
    t_edge(1:2, 5) = (/  5, 13 /)
    t_edge(1:2, 6) = (/  6, 14 /)
    t_edge(1:2, 7) = (/  7, 15 /)
    t_edge(1:2, 8) = (/  8, 16 /)

end if

end subroutine getpos_FZ422_

!--------------------------------------------------------------------------
recursive subroutine getpos_FZ32_(self, dims, cpos, s_edge, t_edge, ns, d, nt, MFZ)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! initialize the PoVRay output for rotational group 32 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 

integer(kind=irg),INTENT(IN)          :: dims(3)
 !! array dimensions
real(kind=dbl),INTENT(INOUT)          :: cpos(3,dims(1))
 !! vertex coordinates
integer(kind=irg),INTENT(INOUT)       :: s_edge(2,dims(2))
 !! first set of edge connectivities
integer(kind=irg),INTENT(INOUT)       :: t_edge(2,dims(3))
 !! second set of edge connectivities
integer(kind=irg),INTENT(OUT)         :: ns
 !! aux parameter
real(kind=dbl),INTENT(OUT)            :: d
 !! aux parameter
integer(kind=irg),INTENT(OUT)         :: nt
 !! aux parameter
logical,OPTIONAL,INTENT(IN)           :: MFZ
 !! (optional) return coordinates for Mackenzie FZ instead of regular FZ

real(kind=dbl)    :: a = 0.8660254038D0, b = 0.5D0, c = 0.5773502692D0, dt = 0.34314575050D0, &
                     ds = 0.6340506711D0, dd, z = 0.D0, oo = 1.D0, o = 0.86602540378443D0, p = 0.5D0

d = 1.1547005384D0
if (present(MFZ)) then
  if (MFZ) then ! define the coordinates of the tetragonal Mackenzie FZ in Rodrigues Space
    d = 1.0
    cpos(1:3, 1) = (/  z,  z,  z /)
    cpos(1:3, 2) = (/  z,  z,  c /)
    cpos(1:3, 3) = (/  a, -p,  z /)
    cpos(1:3, 4) = (/  a, -p,  c /)
    cpos(1:3, 5) = (/  a,  p,  z /)
    cpos(1:3, 6) = (/  a,  p,  c /)

    cpos = cpos/d

    ns = 200
    nt = nint( ns * dt/ds )

! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1, 2 /)
    s_edge(1:2, 2) = (/  1, 3 /)
    s_edge(1:2, 3) = (/  1, 5 /)
    s_edge(1:2, 4) = (/  2, 4 /)
    s_edge(1:2, 5) = (/  2, 6 /)
    s_edge(1:2, 6) = (/  3, 4 /)
    s_edge(1:2, 7) = (/  3, 5 /)
    s_edge(1:2, 8) = (/  4, 6 /)
    s_edge(1:2, 9) = (/  5, 6 /)
  end if 
else ! define the coordinates of the cubic FZ in Rodrigues Space
    cpos(1:3, 1) = (/  a,  b,  c /)
    cpos(1:3, 2) = (/  z, oo,  c /)
    cpos(1:3, 3) = (/ -a,  b,  c /)

    cpos(1:3, 4) = (/ -a, -b,  c /)
    cpos(1:3, 5) = (/  z,-oo,  c /)
    cpos(1:3, 6) = (/  a, -b,  c /)

    cpos(1:3, 7) = (/  a,  b, -c /)
    cpos(1:3, 8) = (/  z, oo, -c /)
    cpos(1:3, 9) = (/ -a,  b, -c /)

    cpos(1:3,10) = (/ -a, -b, -c /)
    cpos(1:3,11) = (/  z,-oo, -c /)
    cpos(1:3,12) = (/  a, -b, -c /)

    cpos = cpos / d

    ns = 200
    nt = nint( ns * dt/ds )

    ! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1,  2 /)
    s_edge(1:2, 2) = (/  2,  3 /)
    s_edge(1:2, 3) = (/  3,  4 /)
    s_edge(1:2, 4) = (/  4,  5 /)
    s_edge(1:2, 5) = (/  5,  6 /)
    s_edge(1:2, 6) = (/  6,  1 /)
    s_edge(1:2, 7) = (/  7,  8 /)
    s_edge(1:2, 8) = (/  8,  9 /)
    s_edge(1:2, 9) = (/  9, 10 /)
    s_edge(1:2,10) = (/ 10, 11 /)
    s_edge(1:2,11) = (/ 11, 12 /)
    s_edge(1:2,12) = (/ 12,  7 /)

    t_edge(1:2, 1) = (/  1,  7 /)
    t_edge(1:2, 2) = (/  2,  8 /)
    t_edge(1:2, 3) = (/  3,  9 /)
    t_edge(1:2, 4) = (/  4, 10 /)
    t_edge(1:2, 5) = (/  5, 11 /)
    t_edge(1:2, 6) = (/  6, 12 /)

end if

end subroutine getpos_FZ32_

!--------------------------------------------------------------------------
recursive subroutine getpos_FZ222_(self, dims, cpos, s_edge, t_edge, ns, d, nt, MFZ)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! initialize the PoVRay output for rotational group 222 

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 

integer(kind=irg),INTENT(IN)          :: dims(3)
 !! array dimensions
real(kind=dbl),INTENT(INOUT)          :: cpos(3,dims(1))
 !! vertex coordinates
integer(kind=irg),INTENT(INOUT)       :: s_edge(2,dims(2))
 !! first set of edge connectivities
integer(kind=irg),INTENT(INOUT)       :: t_edge(2,dims(3))
 !! second set of edge connectivities
integer(kind=irg),INTENT(OUT)         :: ns
 !! aux parameter
real(kind=dbl),INTENT(OUT)            :: d
 !! aux parameter
integer(kind=irg),INTENT(OUT)         :: nt
 !! aux parameter
logical,OPTIONAL,INTENT(IN)           :: MFZ
 !! (optional) return coordinates for Mackenzie FZ instead of regular FZ

real(kind=dbl)    :: a = 1.0D0, b = 1.0D0, c = 1D0, dt = 2.0D0, &
                     ds = 2.0D0, dd, z = 0.D0, oo = 1.D0

d = 1.7320508075688772D0
if (present(MFZ)) then
  if (MFZ) then ! define the coordinates of the tetragonal Mackenzie FZ in Rodrigues Space
    d = 1.0
    cpos(1:3, 1) = (/  z, -a,  z /)
    cpos(1:3, 2) = (/  z, -a,  c /)
    cpos(1:3, 3) = (/  a, -a,  z /)
    cpos(1:3, 4) = (/  a, -a,  c /)
    cpos(1:3, 5) = (/  a,  a,  z /)
    cpos(1:3, 6) = (/  a,  a,  c /)
    cpos(1:3, 7) = (/  z,  a,  z /)
    cpos(1:3, 8) = (/  z,  a,  c /)

    cpos = cpos/d

    ns = 200
    nt = nint( ns * dt/ds )

! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1, 2 /)
    s_edge(1:2, 2) = (/  1, 3 /)
    s_edge(1:2, 3) = (/  1, 7 /)
    s_edge(1:2, 4) = (/  2, 4 /)
    s_edge(1:2, 5) = (/  2, 8 /)
    s_edge(1:2, 6) = (/  3, 4 /)
    s_edge(1:2, 7) = (/  3, 5 /)
    s_edge(1:2, 8) = (/  4, 6 /)
    s_edge(1:2, 9) = (/  5, 6 /)
    s_edge(1:2, 10) = (/  5, 7 /)
    s_edge(1:2, 11) = (/  6, 8 /)
    s_edge(1:2, 12) = (/  7, 8 /)
  end if 
else ! define the coordinates of the FZ in Rodrigues Space
    cpos(1:3, 1) = (/  a,  b,  c /)
    cpos(1:3, 2) = (/ -a,  b,  c /)
    cpos(1:3, 3) = (/ -a, -b,  c /)
    cpos(1:3, 4) = (/  a, -b,  c /)
    cpos(1:3, 5) = (/  a,  b, -c /)
    cpos(1:3, 6) = (/ -a,  b, -c /)
    cpos(1:3, 7) = (/ -a, -b, -c /)
    cpos(1:3, 8) = (/  a, -b, -c /)

    cpos = cpos / d

    ns = 200
    nt = nint( ns * dt/ds )

    ! define the connectivity of all the edges
    s_edge(1:2, 1) = (/  1,  2 /)
    s_edge(1:2, 2) = (/  2,  3 /)
    s_edge(1:2, 3) = (/  3,  4 /)
    s_edge(1:2, 4) = (/  4,  1 /)
    s_edge(1:2, 5) = (/  5,  6 /)
    s_edge(1:2, 6) = (/  6,  7 /)
    s_edge(1:2, 7) = (/  7,  8 /)
    s_edge(1:2, 8) = (/  8,  5 /)

    t_edge(1:2, 1) = (/  1,  5 /)
    t_edge(1:2, 2) = (/  2,  6 /)
    t_edge(1:2, 3) = (/  3,  7 /)
    t_edge(1:2, 4) = (/  4,  8 /)
end if

end subroutine getpos_FZ222_

!--------------------------------------------------------------------------
recursive subroutine drawFZ_(self, SO, rmode, cylr)
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/21/20
 !!
 !! initialize the PoVRay output for any rotation group
 !!
 !! This routine draws the outline of either the Rodrigues Fundamental
 !! zone, or the Mackenzie Fundamental Zone (if SO%getMK() is true).

use mod_rotations 
use mod_so3

IMPLICIT NONE

class(PoVRay_T),INTENT(INOUT)         :: self 
type(so3_T),INTENT(INOUT)             :: SO
integer(kind=irg),INTENT(IN)          :: rmode
 !! 1(cubochoric)|2(homochoric)|3(stereographic)|4(Rodrigues)|5(Euler)
real(kind=dbl),INTENT(IN)             :: cylr
 !! cylinder radius

type(e_T)                             :: eul, eu, euld, eulast 
type(r_T)                             :: ro1, ro2, ro, rolast, ron
type(q_T)                             :: qu
type(s_T)                             :: sp, splast 
type(h_T)                             :: h, ho, holast
type(o_T)                             :: om 
type(c_T)                             :: cu, culast
type(a_T)                             :: axang
type(orientation_T)                   :: ot

real(kind=dbl)                        :: rmax, dx, r, xmax, x, y, z, zsmall, ac, sh(3), xx, d, dd, &
                                         tpi, hpi, aux3(3), aux(4), aux4a(4), aux4b(4)

integer(kind=irg),allocatable         :: s_edge(:,:), t_edge(:,:)
real(kind=dbl),allocatable            :: cpos(:,:)

logical                               :: doMFZ, twostep
integer(kind=irg)                     :: i,j,k, icnt, imax, nt, ns, dims(3), FZtype, FZorder 

call setRotationPrecision('Double')

tpi = 2.D0 * cPi
hpi = 0.5D0 * cPi
doMFZ = SO%getMK() 
call SO%getFZtypeandorder(FZtype, FZorder)

if (FZtype.eq.2) then
    if (FZorder.eq.6) then
      if (doMFZ) then
        twostep = .FALSE.
        dims = (/ 8, 12, 1 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ622(dims, cpos, s_edge, t_edge, ns, d, nt, doMFZ)
      else
        twostep = .TRUE.
        dims = (/ 24, 24, 12 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ622(dims, cpos, s_edge, t_edge, ns, d, nt)
      end if
    end if
    if (FZorder.eq.4) then
      if (doMFZ) then
        twostep = .FALSE.
        dims = (/ 8, 12, 1 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ422(dims, cpos, s_edge, t_edge, ns, d, nt, doMFZ)
      else
        twostep = .TRUE.
        dims = (/ 16, 16, 8 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ422(dims, cpos, s_edge, t_edge, ns, d, nt)
      end if
    end if
    if (FZorder.eq.3) then
      if (doMFZ) then
        twostep = .FALSE.
        dims = (/ 6, 9, 1 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ32(dims, cpos, s_edge, t_edge, ns, d, nt, doMFZ)
      else
        twostep = .TRUE.
        dims = (/ 12, 12, 6 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ32(dims, cpos, s_edge, t_edge, ns, d, nt)
      end if
    end if
    if (FZorder.eq.2) then
      if (doMFZ) then
        twostep = .FALSE.
        dims = (/ 8, 12, 1 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ222(dims, cpos, s_edge, t_edge, ns, d, nt, doMFZ)
      else
        twostep = .TRUE.
        dims = (/ 8, 8, 4 /)
        allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
        call self%getpos_FZ222(dims, cpos, s_edge, t_edge, ns, d, nt)
      end if
    end if
end if

if (FZtype.eq.3) then
! rotational group 23
    if (doMFZ) then
      twostep = .FALSE.
      dims = (/ 4, 6, 1 /)
      allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
      call self%getpos_FZ23(dims, cpos, s_edge, t_edge, ns, d, nt, doMFZ)
    else
      twostep = .FALSE.
      dims = (/ 6, 12, 1 /)
      allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
      call self%getpos_FZ23(dims, cpos, s_edge, t_edge, ns, d, nt)
    end if
end if

if (FZtype.eq.4) then
! rotational group 432
    if (doMFZ) then
      twostep = .FALSE.
      dims = (/ 6, 9, 1 /)
      allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
      call self%getpos_FZ432(dims, cpos, s_edge, t_edge, ns, d, nt, doMFZ)
    else
      twostep = .TRUE.
      dims = (/ 24, 12, 24 /)
      allocate(cpos(3,dims(1)), s_edge(2,dims(2)), t_edge(2,dims(3)))
      call self%getpos_FZ432(dims, cpos, s_edge, t_edge, ns, d, nt)
    end if
end if

! add the reference frame and any necessary wireframes 
if (rmode.eq.1) then
  ac = 0.5D0 * LPs%ap
  call self%addReferenceFrame(ac, cylr)
  call self%addCubochoricCube()
end if
if (rmode.eq.2) then
  ac = 1.33067D0
  call self%addReferenceFrame(ac, cylr)
  call self%addWireFrameSphere(ac)
end if
if (rmode.eq.3) then 
  ac = 1.0D0
  call self%addReferenceFrame(ac, cylr)
  call self%addWireFrameSphere(ac)
end if
if (rmode.eq.4) then 
  ac = 1.0D0
  call self%addReferenceFrame(ac, cylr)
end if
if (rmode.eq.5) then 
  call self%addEulerBox()
end if

! and next, draw the outline of the FZ or MFZ
if ((rmode.eq.1).or.(rmode.eq.2)) then
! create the square edges first
 dx = 1.D0/dble(ns)
 do i=1,dims(2)
  ro1 = r_T( rdinp = (/ cpos(1:3,s_edge(1,i)), d/) )
  ro2 = r_T( rdinp = (/ cpos(1:3,s_edge(2,i)), d/) )
  culast = ro1%rc()
  holast = ro1%rh()
  do j=1,ns+1
    aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
    xx = dsqrt( sum (aux(1:3)**2) )
    ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
    cu = ro%rc()
    ho = ro%rh()
! and create a cylinder with these points
    if (rmode.eq.1) then
      call self%addCylinder(culast%c_copyd(),cu%c_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
    else
      call self%addCylinder(holast%h_copyd(),ho%h_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
    end if
    culast = cu
    holast = ho
  end do 
 end do

 if (twostep) then
   dx = 1.D0/dble(nt)
   do i=1,dims(3)
    ro1 = r_T( rdinp = (/ cpos(1:3,t_edge(1,i)), d/) )
    ro2 = r_T( rdinp = (/ cpos(1:3,t_edge(2,i)), d/) )
    culast = ro1%rc()
    holast = ro1%rh()
    do j=1,nt+1
      aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
      xx = dsqrt( sum (aux(1:3)**2) )
      ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
      cu = ro%rc()
      ho = ro%rh()
  ! and create a cylinder with these points
      if (rmode.eq.1) then
        call self%addCylinder(culast%c_copyd(),cu%c_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      else
        call self%addCylinder(holast%h_copyd(),ho%h_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      end if
      culast = cu
      holast = ho
    end do 
   end do
  end if 
end if

if ((rmode.eq.3).or.(rmode.eq.4)) then
 dx = 1.D0/dble(ns)
 do i=1,dims(2)
  ro1 = r_T( rdinp = (/ cpos(1:3,s_edge(1,i)), d /) )
  ro2 = r_T( rdinp = (/ cpos(1:3,s_edge(2,i)), d /) )
  culast = ro1%rc()
  holast = ro1%rh()
  rolast = ro1
  qu = ro1%rq()
  splast = qu%qs() 
  do j=1,ns+1
    aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
    xx = dsqrt( sum (aux(1:3)**2) )
    ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
    qu = ro%rq()
    sp = qu%qs() 
! and create a cylinder with these points
    if (rmode.eq.3) then
      call self%addCylinder(splast%s_copyd(),sp%s_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
    else
      aux4a = rolast%r_copyd()
      aux4b = ro%r_copyd()
      call self%addCylinder(aux4a(1:3)*aux4a(4),aux4b(1:3)*aux4b(4),cylr,(/ 0.0, 0.0, 1.0 /))
    end if
    rolast = ro
    splast = sp
  end do 
 end do

 if (twostep) then
   dx = 1.D0/dble(nt)
   do i=1,dims(3)
    ro1 = r_T( rdinp = (/ cpos(1:3,t_edge(1,i)), d /) )
    ro2 = r_T( rdinp = (/ cpos(1:3,t_edge(2,i)), d /) )
    rolast = ro1
    qu = ro1%rq()
    splast = qu%qs() 
    do j=1,nt+1
      aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
      xx = dsqrt( sum (aux(1:3)**2) )
      ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
      qu = ro%rq()
      sp = qu%qs() 
  ! and create a cylinder with these points
      if (rmode.eq.3) then
        call self%addCylinder(splast%s_copyd(),sp%s_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      else
        aux4a = rolast%r_copyd()
        aux4b = ro%r_copyd()
        call self%addCylinder(aux4a(1:3)*aux4a(4),aux4b(1:3)*aux4b(4),cylr,(/ 1.0, 0.0, 0.0 /))
      end if
      rolast = ro
      splast = sp
    end do 
   end do
 end if
end if
  
! and next, draw the outline of the FZ or MFZ
if (rmode.eq.5) then
  sh = (/ cPi, cPi/2.D0, cPi /)
! create the square edges first
 dx = 1.D0/dble(ns)
 do i=1,dims(2)
  ro1 = r_T( rdinp = (/ cpos(1:3,s_edge(1,i)), d /) )
  ro2 = r_T( rdinp = (/ cpos(1:3,s_edge(2,i)), d /) )
  eulast = ro1%re()
  aux3 = eulast%e_copyd()
  aux3(1) = mod(aux3(1)+10.D0*cPi,2.D0*cPi)
  aux3(2) = mod(aux3(2)+10.D0*cPi,cPi)
  aux3(3) = mod(aux3(3)+10.D0*cPi,2.D0*cPi)
  eulast = e_T( edinp = aux3 - sh )
  do j=1,ns+1
    aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
    xx = dsqrt( sum (aux(1:3)**2) )
    ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
    eu = ro%re()
    aux3 = eu%e_copyd()
    aux3(1) = mod(aux3(1)+10.D0*cPi,2.D0*cPi)
    aux3(2) = mod(aux3(2)+10.D0*cPi,cPi)
    aux3(3) = mod(aux3(3)+10.D0*cPi,2.D0*cPi)
    eu = e_T( edinp = aux3 - sh )
! and create a cylinder with these points
    aux3 = eulast%e_copyd() - eu%e_copyd()
    if (maxval(abs(aux3)).lt.cPi) then 
      call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
    end if
    eulast = eu
  end do 
 end do

 if (twostep) then
   dx = 1.D0/dble(nt)
   do i=1,dims(3)
    ro1 = r_T( rdinp = (/ cpos(1:3,t_edge(1,i)), d /) )
    ro2 = r_T( rdinp = (/ cpos(1:3,t_edge(2,i)), d /) )
    eulast = ro1%re()
    aux3 = eulast%e_copyd()
    aux3(1) = mod(aux3(1)+10.D0*cPi,2.D0*cPi)
    aux3(2) = mod(aux3(2)+10.D0*cPi,cPi)
    aux3(3) = mod(aux3(3)+10.D0*cPi,2.D0*cPi)
    eulast = e_T( edinp = aux3 - sh )
    do j=1,nt+1
      aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
      xx = dsqrt( sum (aux(1:3)**2) )
      ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
      eu = ro%re()
      aux3 = eu%e_copyd()
      aux3(1) = mod(aux3(1)+10.D0*cPi,2.D0*cPi)
      aux3(2) = mod(aux3(2)+10.D0*cPi,cPi)
      aux3(3) = mod(aux3(3)+10.D0*cPi,2.D0*cPi)
      eu = e_T( edinp = aux3 - sh )
  ! and create a cylinder with these points
      aux3 = eulast%e_copyd() - eu%e_copyd()
      if (maxval(abs(aux3)).lt.cPi) then 
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
      end if
      eulast = eu
    end do 
   end do
  end if 
! for the Euler representation of the Rodrigues fundamental zones we also need to draw
! a few additional lines to complete the volume appearance of the FZ; this depends on
! the order of the FZ, naturally, so we have a couple of possible cases...

! the diagonal lines in the Phi=0 plane are the most important lines; they should be drawn
! for all the rotation groups except for the identity
  if ((FZtype.ge.1).and.(FZtype.le.2)) then
    if (FZorder.ne.0) then  ! we have cyclic or dihedral
      xx = cPi/dble(FZorder)
! draw four diagonal lines
      eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
      eulast = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
      call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
      eulast = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
      call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
      eulast = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
      call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
      eulast = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
      call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
! for cyclic groups we also need to draw the diagonals in the top surface
! and the vertical lines connecting bottom and top planes
      if (FZtype.eq.1) then
! top plane
        eu = e_T( edinp = (/ xx, cPi, 0.D0 /) - sh )
        eulast = e_T( edinp = (/ 0.D0, cPi, xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi-xx, cPi, 0.D0 /) - sh )
        eulast = e_T( edinp = (/ 0.D0, cPi, tpi-xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi, cPi, xx /) - sh )
        eulast = e_T( edinp = (/ xx, cPi, tpi /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi, cPi, tpi-xx /) - sh )
        eulast = e_T( edinp = (/ tpi-xx, cPi, tpi /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
! verticals
        eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
        eulast = e_T( edinp = (/ xx, cPi, 0.D0 /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
        eulast = e_T( edinp = (/ 0.D0, cPi, xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
        eulast = e_T( edinp = (/ tpi-xx, cPi, 0.D0 /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
        eulast = e_T( edinp = (/ 0.D0, cPi, tpi-xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
        eulast = e_T( edinp = (/ tpi, cPi, xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
        eulast = e_T( edinp = (/ xx, cPi, tpi /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
        eulast = e_T( edinp = (/ tpi, cPi, tpi-xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
        eulast = e_T( edinp = (/ tpi-xx, cPi, tpi /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      end if
      if (FZtype.eq.2) then
! in this case, the verticals need to be drawn but only up to the level of the FZ surface
        eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
        eulast = e_T( edinp = (/ xx, hPi, 0.D0 /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
        eulast = e_T( edinp = (/ 0.D0, hPi, xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
        eulast = e_T( edinp = (/ tpi-xx, hPi, 0.D0 /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
        eulast = e_T( edinp = (/ 0.D0, hPi, tpi-xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
        eulast = e_T( edinp = (/ tpi, hPi, xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
        eulast = e_T( edinp = (/ xx, hPi, tpi /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
        eulast = e_T( edinp = (/ tpi, hPi, tpi-xx /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
        eu = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
        eulast = e_T( edinp = (/ tpi-xx, hPi, tpi /) - sh )
        call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
      end if
    end if
  end if
end if

if (FZtype.eq.3) then
    xx = cPi/dble(2)
! draw four diagonal lines
    eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
    eulast = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
    eulast = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 

! and finally the corner posts
    eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ hPi, 0.D0, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 

    eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi - hPi, 0.D0, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi, 0.D0, hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 

    eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi - hPi, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi, 0.D0, tpi-hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 

    eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ 0.D0 + hPi, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, tpi-hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
end if

if ((FZtype.eq.4).and.(rmode.eq.5)) then
    xx = cPi/dble(4)
! draw four diagonal lines
    eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
    eulast = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
    eulast = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
! and verticals
    hPi = hPi * 0.5D0
    eu = e_T( edinp = (/ xx, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ xx, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ 0.D0, 0.D0, xx /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi-xx, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi-xx, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ 0.D0, 0.D0, tpi-xx /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, tpi-xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, 0.D0, xx /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ xx, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ xx, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, 0.D0, tpi-xx /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, tpi-xx /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi-xx, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi-xx, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
! and the closing segments
    eu = e_T( edinp = (/ xx, hPi, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ 0.D0, hPi, xx /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi-xx, hPi, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ 0.D0, hPi, tpi-xx /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, hPi, xx /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ xx, hPi, tpi /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi, hPi, tpi-xx /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
    eu = e_T( edinp = (/ tpi-xx, hPi, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /))
! and finally the corner posts
    eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ hPi, 0.D0, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 

    eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi - hPi, 0.D0, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, 0.D0 /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, 0.D0 /) - sh )
    eulast = e_T( edinp = (/ tpi, 0.D0, hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 

    eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi - hPi, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ tpi, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ tpi, 0.D0, tpi-hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 

    eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ 0.D0 + hPi, 0.D0, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ 0.D0, hPi, tpi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
    eu = e_T( edinp = (/ 0.D0, 0.D0, tpi /) - sh )
    eulast = e_T( edinp = (/ 0.D0, 0.D0, tpi-hPi /) - sh )
    call self%addCylinder(eulast%e_copyd(),eu%e_copyd(),cylr,(/ 1.0, 0.0, 0.0 /)) 
end if

end subroutine drawFZ_

!--------------------------------------------------------------------------
recursive subroutine initFZCyclic_(self, FZorder, cylr, rmode)
!! author: MDG 
!! version: 1.0 
!! date: 01/21/20
!!
!! generate the PoVRay output for the cyclic rotational groups

use mod_rotations

IMPLICIT NONE

class(PoVRay_T), INTENT(INOUT)        :: self 
integer(kind=irg),INTENT(IN)          :: FZorder
 !! 2, 3, 4, or 6
real(kind=dbl),INTENT(IN)             :: cylr
 !! cylinder radius
integer(kind=irg),INTENT(IN)          :: rmode
 !! 1(cubochoric)|2(homochoric)|3(stereographic)|4(Rodrigues)

type(e_T)                             :: eul, eu, euld, eulast 
type(r_T)                             :: ro1, ro2, ro, rolast, ron
type(q_T)                             :: qu
type(s_T)                             :: sp, splast 
type(h_T)                             :: h, ho, holast
type(o_T)                             :: om 
type(c_T)                             :: cu, culast
type(a_T)                             :: axang
type(orientation_T)                   :: ot

real(kind=dbl)                        :: rmax, dx, r, xmax, x, y, z, zsmall, ac, sh(3), xx, &
                                         tpi, hpi, aux(3), aux4a(4), aux4b(4)

integer(kind=irg)                     :: i,j,k, icnt, imax, nt, ns, idpos, icpos, ihedge 
integer(kind=irg),allocatable         :: h_edge(:,:)
real(kind=dbl),allocatable            :: cpos(:,:), dpos(:)
! parameters that depend on the cyclic group
real(kind=dbl)                        :: a, b, c, dt, ds, d, dd, zz, oo, c2, tmp

select case(FZorder)
  case(2) ! define the coordinates of the monoclinic C2 (2) FZ in Rodrigues Space
    a = 570.289922125538D0 
    b = 1.0D0 
    c = 1.0D0 
    dt = 114.57984425107713D0
    ds = 2.0D0 
    d = 1.7320508075688772D0 
    zz = 0.D0 
    oo = 1.D0
    idpos = 200
    icpos = 200
    ihedge = 100
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) ) 

    do i=-12,12
      if (abs(i).ne.12) then 
        cpos(1:3, 13+i) = (/ -a, dtan(dble(i)*15.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 13+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,13+i) = -cpos(2,13+i)
      end if
    end do

    do i=1,25
      cpos(1:3,25+i) = cpos(1:3,i)
      cpos(1,25+i) = -cpos(1,25+i)
    end do

    do i=1,50
      cpos(1:3,50+i) = cpos(1:3,i)
      tmp = cpos(1,50+i)
      cpos(1,50+i) = cpos(2,50+i)
      cpos(2,50+i) = tmp
    end do

    do i=1,100
      cpos(1:3,100+i) = cpos(1:3,i)
      cpos(3,100+i) = -cpos(3,100+i)
    end do

    ! and normalize
    do i=1,200
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

    ns = 2000
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,25
      h_edge(1:2,i)    = (/    i, 25+i /)
      h_edge(1:2,25+i) = (/ 50+i, 75+i /)
      h_edge(1:2,50+i) = (/100+i,125+i /)
      h_edge(1:2,75+i) = (/150+i,175+i /)
    end do
  case(3) ! define the coordinates of the trigonal C3 (3) FZ in Rodrigues Space
    a = 57.289922125538D0
    b = 1.0D0
    c = 0.577350269120D0
    dt = 114.57984425107713D0
    ds = 2.0D0
    d = 1.7320508075688772D0
    zz = 0.D0
    oo = 1.D0
    c2 = 1.7320508075688767D0
    icpos = 208
    idpos = 208
    ihedge = 104
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) ) 

    do i=-6,6 
      if (abs(i).ne.6) then 
        cpos(1:3, 7+i) = (/ -a, dtan(dble(i)*30.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 7+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,7+i) = -cpos(2,7+i)
      end if
    end do
    do i=1,13
      cpos(1:3,13+i) = cpos(1:3,i)
      cpos(1,13+i) = -cpos(1,13+i)
    end do

    do i=1,26
      cpos(1:3,26+i) = cpos(1:3,i)
      tmp = cpos(1,26+i)
      cpos(1,26+i) = cpos(2,26+i)
      cpos(2,26+i) = tmp
    end do

    do i=1,52
      cpos(1:3,52+i) = cpos(1:3,i)
      cpos(3,52+i) = -cpos(3,52+i)
    end do

    do i=1,104
      cpos(1:3,104+i) = cpos(1:3,i)
      if (cpos(3,104+i).lt.0.D0) then
        cpos(3,104+i) = -c2
      else 
        cpos(3,104+i) = c2
      end if
    end do

    ! and normalize
    do i=1,208
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

    ns = 200
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,13
      h_edge(1:2,i)    = (/    i, 13+i /)
      h_edge(1:2,13+i) = (/ 26+i, 39+i /)
      h_edge(1:2,26+i) = (/ 52+i, 65+i /)
      h_edge(1:2,39+i) = (/ 78+i, 91+i /)

      h_edge(1:2,52+i) = (/104+i,117+i /)
      h_edge(1:2,65+i) = (/130+i,143+i /)
      h_edge(1:2,78+i) = (/156+i,169+i /)
      h_edge(1:2,91+i) = (/182+i,195+i /)
    end do
  case(4) ! define the coordinates of the tetragonal C4 (4) FZ in Rodrigues Space

    ! NEEDS TO BE FIXED !!!
    a = 57.289922125538D0
    b = 1.0D0
    c = 1.0D0
    dt = 114.57984425107713D0
    ds = 2.0D0
    d = 1.7320508075688772D0
    zz = 0.D0
    oo = 1.D0
    icpos = 104
    idpos = 104
    ihedge = 52
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) ) 

    do i=-6,6 
      if (abs(i).ne.6) then 
        cpos(1:3, 7+i) = (/ -a, dtan(dble(i)*30.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 7+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,7+i) = -cpos(2,7+i)
      end if
    end do
    do i=1,13
      cpos(1:3,13+i) = cpos(1:3,i)
      cpos(1,13+i) = -cpos(1,13+i)
    end do

    do i=1,26
      cpos(1:3,26+i) = cpos(1:3,i)
      tmp = cpos(1,26+i)
      cpos(1,26+i) = cpos(2,26+i)
      cpos(2,26+i) = tmp
    end do

    do i=1,52
      cpos(1:3,52+i) = cpos(1:3,i)
      cpos(3,52+i) = -cpos(3,52+i)
    end do

    ! and normalize
    do i=1,104
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

    ns = 200
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,13
      h_edge(1:2,i)    = (/    i, 13+i /)
      h_edge(1:2,13+i) = (/ 26+i, 39+i /)
      h_edge(1:2,26+i) = (/ 52+i, 65+i /)
      h_edge(1:2,39+i) = (/ 78+i, 91+i /)
    end do
  case(6) ! define the coordinates of the hexagonal C6 (6) FZ in Rodrigues Space
    a = 57.289922125538D0
    b = 1.0D0
    c = 1.0D0
    dt = 114.57984425107713D0
    ds = 2.0D0
    d = 1.7320508075688772D0
    zz = 0.D0
    oo = 1.D0
    icpos = 104
    idpos = 104
    ihedge = 52
    allocate(cpos(3,icpos), h_edge(2,ihedge), dpos(idpos) ) 

    do i=-6,6 
      if (abs(i).ne.6) then 
        cpos(1:3, 7+i) = (/ -a, dtan(dble(i)*30.D0*dtor*0.5D0) ,  c /)
      else
        cpos(1:3, 7+i) = (/ -a, a ,  c /)
        if (i.lt.0) cpos(2,7+i) = -cpos(2,7+i)
      end if
    end do
    do i=1,13
      cpos(1:3,13+i) = cpos(1:3,i)
      cpos(1,13+i) = -cpos(1,13+i)
    end do

    do i=1,26
      cpos(1:3,26+i) = cpos(1:3,i)
      tmp = cpos(1,26+i)
      cpos(1,26+i) = cpos(2,26+i)
      cpos(2,26+i) = tmp
    end do

    do i=1,52
      cpos(1:3,52+i) = cpos(1:3,i)
      cpos(3,52+i) = -cpos(3,52+i)
    end do

    ! and normalize
    do i=1,104
      dpos(i) = dsqrt(sum(cpos(1:3,i)*cpos(1:3,i)))
    end do

    ns = 200
    dx = dt/float(ns-1)

    ! define the connectivity of all the edges
    do i=1,13
      h_edge(1:2,i)    = (/    i, 13+i /)
      h_edge(1:2,13+i) = (/ 26+i, 39+i /)
      h_edge(1:2,26+i) = (/ 52+i, 65+i /)
      h_edge(1:2,39+i) = (/ 78+i, 91+i /)
    end do
  case default
end select 

! add the reference frame and any necessary wireframes 
if (rmode.eq.1) then
  ac = 0.5D0 * LPs%ap
  call self%addReferenceFrame(ac, cylr)
  call self%addCubochoricCube()
end if
if (rmode.eq.2) then
  ac = 1.33067D0
  call self%addReferenceFrame(ac, cylr)
  call self%addWireFrameSphere(ac)
end if
if (rmode.eq.3) then 
  ac = 1.0D0
  call self%addReferenceFrame(ac, cylr)
  call self%addWireFrameSphere(ac)
end if
if (rmode.eq.4) then 
  ac = 1.0D0
  call self%addReferenceFrame(ac, cylr)
end if
if (rmode.eq.5) then 
  call self%addEulerBox()
end if



if ((rmode.eq.1).or.(rmode.eq.2)) then
  ! create the square edges first
 dx = 1.D0/dble(ns)
 do i=1, ihedge
  ro1 = r_T( rdinp = (/ cpos(1:3,h_edge(1,i))/dpos(i), dpos(i) /) )
  ro2 = r_T( rdinp = (/ cpos(1:3,h_edge(2,i))/dpos(i), dpos(i) /) )
  culast = ro1%rc()
  holast = ro1%rh()
  do j=1,ns+1
    aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
    xx = dsqrt( sum (aux(1:3)**2) )
    ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
    cu = ro%rc()
    ho = ro%rh()
! and create a cylinder with these points
    if (rmode.eq.1) then
      call self%addCylinder(culast%c_copyd(),cu%c_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
    else
      call self%addCylinder(holast%h_copyd(),ho%h_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
    end if
    culast = cu
    holast = ho
  end do 
 end do
end if
  
if ((rmode.eq.3).or.(rmode.eq.4)) then

 dx = 1.D0/dble(ns)
 do i=1, ihedge
  ro1 = r_T( rdinp = (/ cpos(1:3,h_edge(1,i))/dpos(i), dpos(i) /) )
  ro2 = r_T( rdinp = (/ cpos(1:3,h_edge(2,i))/dpos(i), dpos(i) /) )
  rolast = ro1
  qu = ro1%rq()
  splast = qu%qs() 
  do j=1,ns+1
    aux = d*ro1%r_copyd() + d*(ro2%r_copyd() - ro1%r_copyd()) * j * dx
    xx = dsqrt( sum (aux(1:3)**2) )
    ro = r_T( rdinp = (/ aux(1:3)/xx, xx /) )
    qu = ro%rq()
    sp = qu%qs() 
! and create a cylinder with these points
    if (rmode.eq.3) then
      call self%addCylinder(splast%s_copyd(),sp%s_copyd(),cylr,(/ 0.0, 0.0, 1.0 /))
    else
      aux4a = rolast%r_copyd()
      aux4b = ro%r_copyd()
      call self%addCylinder(aux4a(1:3)*aux4a(4),aux4b(1:3)*aux4b(4),cylr,(/ 0.0, 0.0, 1.0 /))
    end if
    rolast = ro
    splast = sp
  end do 
 end do

end if
  
end subroutine initFZCyclic_

end module mod_povray
