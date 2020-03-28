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

module mod_OrientationViz
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/27/20
  !!
  !! class definition for the EMOrientationViz program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMOrientationViz program
type, public :: OrientationVizNameListType
  integer(kind=irg) :: cubochoric
  integer(kind=irg) :: homochoric
  integer(kind=irg) :: rodrigues
  integer(kind=irg) :: stereographic
  integer(kind=irg) :: eulerspace
  integer(kind=irg) :: reducetoRFZ
  integer(kind=irg) :: nx
  integer(kind=irg) :: ny
  integer(kind=irg) :: nz
  integer(kind=irg) :: overridepgnum
  integer(kind=irg) :: MacKenzieCell
  real(kind=sgl)    :: rgb(3)
  real(kind=sgl)    :: sphrad
  real(kind=sgl)    :: distance
  character(3)      :: scalingmode
  character(3)      :: mrcmode
  character(fnlen)  :: df3file
  character(fnlen)  :: mrcfile
  character(fnlen)  :: framemrcfile
  character(fnlen)  :: xtalname
  character(fnlen)  :: povrayfile
  character(fnlen)  :: anglefile
end type OrientationVizNameListType

! class definition
type, public :: OrientationViz_T
private 
  character(fnlen)                  :: nmldeffile = 'EMOrientationViz.nml'
  type(OrientationVizNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: OrientationViz_
  procedure, pass(self) :: get_cubochoric_
  procedure, pass(self) :: get_homochoric_
  procedure, pass(self) :: get_rodrigues_
  procedure, pass(self) :: get_stereographic_
  procedure, pass(self) :: get_eulerspace_
  procedure, pass(self) :: get_reducetoRFZ_
  procedure, pass(self) :: get_nx_
  procedure, pass(self) :: get_ny_
  procedure, pass(self) :: get_nz_
  procedure, pass(self) :: get_overridepgnum_
  procedure, pass(self) :: get_MacKenzieCell_
  procedure, pass(self) :: get_rgb_
  procedure, pass(self) :: get_sphrad_
  procedure, pass(self) :: get_distance_
  procedure, pass(self) :: get_scalingmode_
  procedure, pass(self) :: get_mrcmode_
  procedure, pass(self) :: get_df3file_
  procedure, pass(self) :: get_mrcfile_
  procedure, pass(self) :: get_framemrcfile_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: get_povrayfile_
  procedure, pass(self) :: get_anglefile_
  procedure, pass(self) :: set_cubochoric_
  procedure, pass(self) :: set_homochoric_
  procedure, pass(self) :: set_rodrigues_
  procedure, pass(self) :: set_stereographic_
  procedure, pass(self) :: set_eulerspace_
  procedure, pass(self) :: set_reducetoRFZ_
  procedure, pass(self) :: set_nx_
  procedure, pass(self) :: set_ny_
  procedure, pass(self) :: set_nz_
  procedure, pass(self) :: set_overridepgnum_
  procedure, pass(self) :: set_MacKenzieCell_
  procedure, pass(self) :: set_rgb_
  procedure, pass(self) :: set_sphrad_
  procedure, pass(self) :: set_distance_
  procedure, pass(self) :: set_scalingmode_
  procedure, pass(self) :: set_mrcmode_
  procedure, pass(self) :: set_df3file_
  procedure, pass(self) :: set_mrcfile_
  procedure, pass(self) :: set_framemrcfile_
  procedure, pass(self) :: set_xtalname_
  procedure, pass(self) :: set_povrayfile_
  procedure, pass(self) :: set_anglefile_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: OrientationViz => OrientationViz_
  generic, public :: get_cubochoric => get_cubochoric_
  generic, public :: get_homochoric => get_homochoric_
  generic, public :: get_rodrigues => get_rodrigues_
  generic, public :: get_stereographic => get_stereographic_
  generic, public :: get_eulerspace => get_eulerspace_
  generic, public :: get_reducetoRFZ => get_reducetoRFZ_
  generic, public :: get_nx => get_nx_
  generic, public :: get_ny => get_ny_
  generic, public :: get_nz => get_nz_
  generic, public :: get_overridepgnum => get_overridepgnum_
  generic, public :: get_MacKenzieCell => get_MacKenzieCell_
  generic, public :: get_rgb => get_rgb_
  generic, public :: get_sphrad => get_sphrad_
  generic, public :: get_distance => get_distance_
  generic, public :: get_scalingmode => get_scalingmode_
  generic, public :: get_mrcmode => get_mrcmode_
  generic, public :: get_df3file => get_df3file_
  generic, public :: get_mrcfile => get_mrcfile_
  generic, public :: get_framemrcfile => get_framemrcfile_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: get_povrayfile => get_povrayfile_
  generic, public :: get_anglefile => get_anglefile_
  generic, public :: set_cubochoric => set_cubochoric_
  generic, public :: set_homochoric => set_homochoric_
  generic, public :: set_rodrigues => set_rodrigues_
  generic, public :: set_stereographic => set_stereographic_
  generic, public :: set_eulerspace => set_eulerspace_
  generic, public :: set_reducetoRFZ => set_reducetoRFZ_
  generic, public :: set_nx => set_nx_
  generic, public :: set_ny => set_ny_
  generic, public :: set_nz => set_nz_
  generic, public :: set_overridepgnum => set_overridepgnum_
  generic, public :: set_MacKenzieCell => set_MacKenzieCell_
  generic, public :: set_rgb => set_rgb_
  generic, public :: set_sphrad => set_sphrad_
  generic, public :: set_distance => set_distance_
  generic, public :: set_scalingmode => set_scalingmode_
  generic, public :: set_mrcmode => set_mrcmode_
  generic, public :: set_df3file => set_df3file_
  generic, public :: set_mrcfile => set_mrcfile_
  generic, public :: set_framemrcfile => set_framemrcfile_
  generic, public :: set_xtalname => set_xtalname_
  generic, public :: set_povrayfile => set_povrayfile_
  generic, public :: set_anglefile => set_anglefile_
end type OrientationViz_T

!DEC$ ATTRIBUTES DLLEXPORT :: getNameList
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList
!DEC$ ATTRIBUTES DLLEXPORT :: OrientationViz
!DEC$ ATTRIBUTES DLLEXPORT :: get_cubochoric
!DEC$ ATTRIBUTES DLLEXPORT :: set_cubochoric
!DEC$ ATTRIBUTES DLLEXPORT :: get_homochoric
!DEC$ ATTRIBUTES DLLEXPORT :: set_homochoric
!DEC$ ATTRIBUTES DLLEXPORT :: get_rodrigues
!DEC$ ATTRIBUTES DLLEXPORT :: set_rodrigues
!DEC$ ATTRIBUTES DLLEXPORT :: get_stereographic
!DEC$ ATTRIBUTES DLLEXPORT :: set_stereographic
!DEC$ ATTRIBUTES DLLEXPORT :: get_eulerspace
!DEC$ ATTRIBUTES DLLEXPORT :: set_eulerspace
!DEC$ ATTRIBUTES DLLEXPORT :: get_reducetoRFZ
!DEC$ ATTRIBUTES DLLEXPORT :: set_reducetoRFZ
!DEC$ ATTRIBUTES DLLEXPORT :: get_nx
!DEC$ ATTRIBUTES DLLEXPORT :: set_nx
!DEC$ ATTRIBUTES DLLEXPORT :: get_ny
!DEC$ ATTRIBUTES DLLEXPORT :: set_ny
!DEC$ ATTRIBUTES DLLEXPORT :: get_nz
!DEC$ ATTRIBUTES DLLEXPORT :: set_nz
!DEC$ ATTRIBUTES DLLEXPORT :: get_overridepgnum
!DEC$ ATTRIBUTES DLLEXPORT :: set_overridepgnum
!DEC$ ATTRIBUTES DLLEXPORT :: get_MacKenzieCell
!DEC$ ATTRIBUTES DLLEXPORT :: set_MacKenzieCell
!DEC$ ATTRIBUTES DLLEXPORT :: get_rgb
!DEC$ ATTRIBUTES DLLEXPORT :: set_rgb
!DEC$ ATTRIBUTES DLLEXPORT :: get_sphrad
!DEC$ ATTRIBUTES DLLEXPORT :: set_sphrad
!DEC$ ATTRIBUTES DLLEXPORT :: get_distance
!DEC$ ATTRIBUTES DLLEXPORT :: set_distance
!DEC$ ATTRIBUTES DLLEXPORT :: get_scalingmode
!DEC$ ATTRIBUTES DLLEXPORT :: set_scalingmode
!DEC$ ATTRIBUTES DLLEXPORT :: get_mrcmode
!DEC$ ATTRIBUTES DLLEXPORT :: set_mrcmode
!DEC$ ATTRIBUTES DLLEXPORT :: get_df3file
!DEC$ ATTRIBUTES DLLEXPORT :: set_df3file
!DEC$ ATTRIBUTES DLLEXPORT :: get_mrcfile
!DEC$ ATTRIBUTES DLLEXPORT :: set_mrcfile
!DEC$ ATTRIBUTES DLLEXPORT :: get_framemrcfile
!DEC$ ATTRIBUTES DLLEXPORT :: set_framemrcfile
!DEC$ ATTRIBUTES DLLEXPORT :: get_xtalname
!DEC$ ATTRIBUTES DLLEXPORT :: set_xtalname
!DEC$ ATTRIBUTES DLLEXPORT :: get_povrayfile
!DEC$ ATTRIBUTES DLLEXPORT :: set_povrayfile
!DEC$ ATTRIBUTES DLLEXPORT :: get_anglefile
!DEC$ ATTRIBUTES DLLEXPORT :: set_anglefile

! the constructor routine for this class 
interface OrientationViz_T
  module procedure OrientationViz_constructor
end interface OrientationViz_T

contains

!--------------------------------------------------------------------------
type(OrientationViz_T) function OrientationViz_constructor( nmlfile ) result(OrientationViz)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! constructor for the OrientationViz_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call OrientationViz%readNameList(nmlfile)

end function OrientationViz_constructor

!--------------------------------------------------------------------------
subroutine OrientationViz_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! destructor for the OrientationViz_T Class
 
IMPLICIT NONE

type(OrientationViz_T), INTENT(INOUT)  :: self 

call reportDestructor('OrientationViz_T')

end subroutine OrientationViz_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! read the namelist from an nml file for the OrientationViz_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)  :: self
character(fnlen),INTENT(IN)             :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)             :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                          :: EMsoft 
type(IO_T)                              :: Message       
logical                                 :: skipread = .FALSE.

integer(kind=irg) :: cubochoric
integer(kind=irg) :: homochoric
integer(kind=irg) :: rodrigues
integer(kind=irg) :: stereographic
integer(kind=irg) :: eulerspace
integer(kind=irg) :: reducetoRFZ
integer(kind=irg) :: nx
integer(kind=irg) :: ny
integer(kind=irg) :: nz
integer(kind=irg) :: overridepgnum
integer(kind=irg) :: MacKenzieCell
real(kind=sgl)    :: rgb(3)
real(kind=sgl)    :: sphrad
real(kind=sgl)    :: distance
character(3)      :: scalingmode
character(3)      :: mrcmode
character(fnlen)  :: df3file
character(fnlen)  :: mrcfile
character(fnlen)  :: framemrcfile
character(fnlen)  :: xtalname
character(fnlen)  :: povrayfile
character(fnlen)  :: anglefile

! define the IO namelist to facilitate passing variables to the program.
namelist  / EMOrientationViz / cubochoric, homochoric, rodrigues, stereographic, eulerspace, &
                               xtalname, povrayfile, anglefile, reducetoRFZ, rgb, sphrad, df3file, &
                               mrcfile, framemrcfile, mrcmode, &
                               nx, ny, nz, distance, scalingmode, overridepgnum, MacKenzieCell

! initialize
cubochoric = 0
homochoric = 0
rodrigues = 0
stereographic = 0
eulerspace = 0
reducetoRFZ = 1
overridepgnum = 0
MacKenzieCell = 0
rgb = (/ 0.0, 0.0, 1.0 /)
sphrad = 0.015
distance = 4.0
nx = 64
ny = 64
nz = 64
scalingmode = 'lev'   ! or 'log' or 'lev' (for equi-level contours)
mrcmode = 'off'       ! 'off', 'reg', 'frm'
df3file = 'undefined'
mrcfile = 'undefined'
framemrcfile = 'undefined'
xtalname = 'undefined'
povrayfile = 'undefined'
anglefile = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EMOrientationViz)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' structure file name is undefined in '//nmlfile)
 end if
 if (mrcmode.eq.'off') then
   if (trim(povrayfile).eq.'undefined') then
    call Message%printError('readNameList:',' povray file name is undefined in '//nmlfile)
   end if
 else
  if (mrcmode.eq.'reg') then
   if (trim(mrcfile).eq.'undefined') then
    call Message%printError('readNameList:',' mrc file name is undefined in '//nmlfile)
   end if
  else
   if (trim(framemrcfile).eq.'undefined') then
    call Message%printError('readNameList:',' frame mrc file name is undefined in '//nmlfile)
   end if
  end if
 end if
 if (trim(anglefile).eq.'undefined') then
  call Message%printError('readNameList:',' angle file name is undefined in '//nmlfile)
 end if
end if

self%nml%cubochoric = cubochoric
self%nml%homochoric = homochoric
self%nml%rodrigues = rodrigues
self%nml%stereographic = stereographic
self%nml%eulerspace = eulerspace
self%nml%reducetoRFZ = reducetoRFZ
self%nml%nx = nx
self%nml%ny = ny
self%nml%nz = nz
self%nml%overridepgnum = overridepgnum
self%nml%MacKenzieCell = MacKenzieCell
self%nml%rgb = rgb
self%nml%sphrad = sphrad
self%nml%distance = distance
self%nml%scalingmode = scalingmode
self%nml%mrcmode = mrcmode
self%nml%df3file = df3file
self%nml%mrcfile = mrcfile
self%nml%framemrcfile = framemrcfile
self%nml%xtalname = trim(xtalname)
self%nml%povrayfile = trim(povrayfile)
self%nml%anglefile = trim(anglefile)

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! pass the namelist for the OrientationViz_T Class to the calling program

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)          :: self
type(OrientationVizNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
function get_cubochoric_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get cubochoric from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%cubochoric

end function get_cubochoric_

!--------------------------------------------------------------------------
subroutine set_cubochoric_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set cubochoric in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%cubochoric = inp

end subroutine set_cubochoric_

!--------------------------------------------------------------------------
function get_homochoric_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get homochoric from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%homochoric

end function get_homochoric_

!--------------------------------------------------------------------------
subroutine set_homochoric_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set homochoric in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%homochoric = inp

end subroutine set_homochoric_

!--------------------------------------------------------------------------
function get_rodrigues_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get rodrigues from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%rodrigues

end function get_rodrigues_

!--------------------------------------------------------------------------
subroutine set_rodrigues_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set rodrigues in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%rodrigues = inp

end subroutine set_rodrigues_

!--------------------------------------------------------------------------
function get_stereographic_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get stereographic from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%stereographic

end function get_stereographic_

!--------------------------------------------------------------------------
subroutine set_stereographic_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set stereographic in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%stereographic = inp

end subroutine set_stereographic_

!--------------------------------------------------------------------------
function get_eulerspace_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get eulerspace from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%eulerspace

end function get_eulerspace_

!--------------------------------------------------------------------------
subroutine set_eulerspace_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set eulerspace in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%eulerspace = inp

end subroutine set_eulerspace_

!--------------------------------------------------------------------------
function get_reducetoRFZ_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get reducetoRFZ from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%reducetoRFZ

end function get_reducetoRFZ_

!--------------------------------------------------------------------------
subroutine set_reducetoRFZ_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set reducetoRFZ in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%reducetoRFZ = inp

end subroutine set_reducetoRFZ_

!--------------------------------------------------------------------------
function get_nx_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get nx from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%nx

end function get_nx_

!--------------------------------------------------------------------------
subroutine set_nx_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set nx in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%nx = inp

end subroutine set_nx_

!--------------------------------------------------------------------------
function get_ny_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get ny from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%ny

end function get_ny_

!--------------------------------------------------------------------------
subroutine set_ny_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set ny in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%ny = inp

end subroutine set_ny_

!--------------------------------------------------------------------------
function get_nz_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get nz from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%nz

end function get_nz_

!--------------------------------------------------------------------------
subroutine set_nz_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set nz in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%nz = inp

end subroutine set_nz_

!--------------------------------------------------------------------------
function get_overridepgnum_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get overridepgnum from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%overridepgnum

end function get_overridepgnum_

!--------------------------------------------------------------------------
subroutine set_overridepgnum_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set overridepgnum in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%overridepgnum = inp

end subroutine set_overridepgnum_

!--------------------------------------------------------------------------
function get_MacKenzieCell_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get MacKenzieCell from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg)                          :: out

out = self%nml%MacKenzieCell

end function get_MacKenzieCell_

!--------------------------------------------------------------------------
subroutine set_MacKenzieCell_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set MacKenzieCell in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)              :: inp

self%nml%MacKenzieCell = inp

end subroutine set_MacKenzieCell_

!--------------------------------------------------------------------------
function get_rgb_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get rgb from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
real(kind=sgl)                             :: out(3)

out = self%nml%rgb

end function get_rgb_

!--------------------------------------------------------------------------
subroutine set_rgb_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set rgb(3) in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)                 :: inp(3)

self%nml%rgb = inp

end subroutine set_rgb_

!--------------------------------------------------------------------------
function get_sphrad_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get sphrad from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
real(kind=sgl)                             :: out

out = self%nml%sphrad

end function get_sphrad_

!--------------------------------------------------------------------------
subroutine set_sphrad_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set sphrad in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)                 :: inp

self%nml%sphrad = inp

end subroutine set_sphrad_

!--------------------------------------------------------------------------
function get_distance_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get distance from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
real(kind=sgl)                             :: out

out = self%nml%distance

end function get_distance_

!--------------------------------------------------------------------------
subroutine set_distance_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set distance in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)                 :: inp

self%nml%distance = inp

end subroutine set_distance_

!--------------------------------------------------------------------------
function get_scalingmode_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get scalingmode from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(3)                               :: out

out = self%nml%scalingmode

end function get_scalingmode_

!--------------------------------------------------------------------------
subroutine set_scalingmode_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set scalingmode in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)                   :: inp

self%nml%scalingmode = inp

end subroutine set_scalingmode_

!--------------------------------------------------------------------------
function get_mrcmode_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get mrcmode from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(3)                               :: out

out = self%nml%mrcmode

end function get_mrcmode_

!--------------------------------------------------------------------------
subroutine set_mrcmode_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set mrcmode in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(3), INTENT(IN)                   :: inp

self%nml%mrcmode = inp

end subroutine set_mrcmode_

!--------------------------------------------------------------------------
function get_df3file_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get df3file from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%df3file

end function get_df3file_

!--------------------------------------------------------------------------
subroutine set_df3file_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set df3file in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%df3file = inp

end subroutine set_df3file_

!--------------------------------------------------------------------------
function get_mrcfile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get mrcfile from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%mrcfile

end function get_mrcfile_

!--------------------------------------------------------------------------
subroutine set_mrcfile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set mrcfile in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%mrcfile = inp

end subroutine set_mrcfile_

!--------------------------------------------------------------------------
function get_framemrcfile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get framemrcfile from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%framemrcfile

end function get_framemrcfile_

!--------------------------------------------------------------------------
subroutine set_framemrcfile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set framemrcfile in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%framemrcfile = inp

end subroutine set_framemrcfile_

!--------------------------------------------------------------------------
function get_xtalname_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get xtalname from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%xtalname

end function get_xtalname_

!--------------------------------------------------------------------------
subroutine set_xtalname_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set xtalname in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
function get_povrayfile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get povrayfile from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%povrayfile

end function get_povrayfile_

!--------------------------------------------------------------------------
subroutine set_povrayfile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set povrayfile in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%povrayfile = inp

end subroutine set_povrayfile_

!--------------------------------------------------------------------------
function get_anglefile_(self) result(out)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! get anglefile from the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen)                           :: out

out = self%nml%anglefile

end function get_anglefile_

!--------------------------------------------------------------------------
subroutine set_anglefile_(self,inp)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! set anglefile in the OrientationViz_T class

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)               :: inp

self%nml%anglefile = inp

end subroutine set_anglefile_

!--------------------------------------------------------------------------
subroutine OrientationViz_(self, EMsoft, progname, progdesc)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! perform the computations

use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry
use mod_rotations
use mod_Lambert
use mod_quaternions
use mod_povray
use mod_so3
use mod_dirstats
use mod_MRC
use mod_math
use mod_so3
use mod_HDFsupport, only: openFortranHDFInterface, closeFortranHDFInterface

IMPLICIT NONE 

class(OrientationViz_T), INTENT(INOUT)  :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 
character(fnlen), INTENT(INOUT)         :: progdesc 

type(IO_T)              :: Message 
type(PoVRay_T)          :: PoVcu, PoVho, PoVro, PoVst, PoVeu 
type(MRC_T)             :: MRC
type(DirStat_T)         :: dict
type(so3_T)             :: SO
type(SpaceGroup_T)      :: SG
type(cell_T)            :: cell
type(QuaternionArray_T) :: Pm, dummy 
type(r_T)               :: ro 
type(h_T)               :: ho 
type(s_T)               :: st  
type(e_T)               :: eu 
type(c_T)               :: cu 

real(kind=dbl)          :: rod(4), sh(3), xyz(3), xyz4(4), XY(2), euFZ(3), rstep, ac, dd, qur(4)
integer(kind=irg)       :: i,j,k, icnt, imax, nt, npx, ngroups, groups(10), dataunit4=25, dataunit5=40, &
                           ierr, ig, ix, iy, iz, num, ixyz(3), pgnum
real(kind=dbl)          :: delta, eps = 1.0D-2
character(fnlen)        :: locationline, fname, dataname, outname, lightline, skyline, rgbstring, colorstring, df3name, mrcname
character(11)           :: p0
character(21)           :: p1, p2
character(5)            :: px, py, pz, pd
character(7)            :: pd2
character(3)            :: gid(10)
character(2)            :: angleformat
integer(kind=irg)       :: FZtype, FZorder, numpoints, seed, istat, FZtype_override, FZorder_override
logical                 :: levelset, skipthis, drawMK
real(kind=dbl),allocatable :: samples(:,:)
real(kind=dbl)          :: muhat(4), kappahat
! rendering volumes
real(kind=sgl),allocatable :: rovol(:,:,:), spvol(:,:,:), euvol(:,:,:), cuvol(:,:,:), hovol(:,:,:)
type(MRCstruct)         :: MRCheader
type(FEIstruct)         :: FEIheaders(1024)
real(kind=dbl),allocatable :: psum(:)
real(kind=dbl),allocatable :: volume(:,:,:)  ! you'll need to fill this array with values ... 
integer(kind=irg)       :: numx, numy, numz       ! set these to the size of the volume array
real(kind=sgl)          :: maxRFZdis(5), rodx, rody, rodz, eudx, eudy, eudz, spdx, spdy, spdz, cudx, cudy, cudz, &
                           hodx, hody, hodz, scalefactors(3,5), acubo, ahomo, grid3(3,3,3), eyepos(3)
type(FZpointd),pointer  :: FZtmp 

associate( enl=>self%nml )

call setRotationPrecision('double')

EMsoft = EMsoft_T( progname, progdesc, silent = .TRUE.)

! init a few parameters and arrays
grid3 = trilinear_splat( (/ 0.0, 0.0, 0.0/), (/ 0.0, 0.0, 0.0/), init=.TRUE.)

! define some parameters 
sh = (/ sngl(cPi), sngl(cPi/2.D0), sngl(cPi) /)   ! offset parameter for primary Euler cell

! get the space group symmetry 
call openFortranHDFInterface()
call cell%getCrystalData(enl%xtalname, SG, EMsoft)
call closeFortranHDFInterface()
pgnum = SG%getPGnumber()

! set the FZtype and FZorder parameters in the SO class 
if (enl%overridepgnum.eq.0) then 
  pgnum = enl%overridepgnum
end if 
SO = so3_T( pgnum, zerolist='FZ')
call SO%getFZtypeandorder(FZtype, FZorder)

! Regular or MacKenzie FZ ?
if (enl%MacKenzieCell.eq.1) then 
  call SO%setMK(.TRUE.)
  call SO%setMFZtypeandorder(pgnum)
  call SO%getMFZtypeandorder(FZtype, FZorder)
end if

! get the symmetry operator quaternions for the point group 
call dummy%QSym_Init(pgnum, Pm)

! make sure that the outname does not have an extension (no . in the string)
outname = trim(enl%povrayfile)
if ((index(trim(outname),'.').ne.0).and.(index(trim(outname),'.').gt.20)) then
  call Message%printError('EMOrientationViz','the output filename may not contain any periods')
end if 
outname = EMsoft%generateFilePath('EMdatapathname', enl%povrayfile)

! PoVRay/DF3 initializations  (this opens the file, sets the camera and light source
! and allocates the rendering volume)
!===========
! cubochoric 
!===========
if (enl%cubochoric.ne.0) then
  if (enl%mrcmode.eq.'off') then
    call initFiles(EMsoft, enl, PoVcu, SO, 'cu', outname, dataunit )
  end if
! create the rendering volume
  allocate(cuvol(-enl%nx:enl%nx,-enl%ny:enl%ny,-enl%nz:enl%nz),stat=istat)
  cuvol = 0.0
! generate the x/y/z scaling factors
! these should be set so that the render box has the outer dimensions of the cubochoric cube,
  acubo = sngl(cPi)**0.6666666 * 0.5
  cudx = float(enl%nx) / acubo
  cudy = float(enl%nx) / acubo
  cudz = float(enl%nz) / acubo
end if

!===========
! homochoric 
!===========
if (enl%homochoric.ne.0) then
  if (enl%mrcmode.eq.'off') then
    call initFiles(EMsoft, enl, PoVcu, SO, 'ho', outname, dataunit2 )
  end if
! create the rendering volume
  allocate(hovol(-enl%nx:enl%nx,-enl%ny:enl%ny,-enl%nz:enl%nz),stat=istat)
  hovol = 0.0
! generate the x/y/z scaling factors
! these should be set so that the render box has the outer dimensions of the cubochoric cube,
  ahomo = (3.0*sngl(cPi)/4.0)**0.3333333
  hodx = float(enl%nx) / ahomo
  hody = float(enl%nx) / ahomo
  hodz = float(enl%nz) / ahomo
end if

!===========
! Rodrigues 
!===========
if (enl%rodrigues.ne.0) then
  if (enl%mrcmode.eq.'off') then
    call initFiles(EMsoft, enl, PoVro, SO, 'ro', outname, dataunit3 )
  end if
! create the rendering volume
  allocate(rovol(-enl%nx:enl%nx,-enl%ny:enl%ny,-enl%nz:enl%nz),stat=istat)
  rovol = 0.0
! generate the x/y/z scaling factors
! these should be set so that the render box has the outer dimensions of the RFZ,
! with the exception of the cyclic groups for which we pick some value that is 
! reasonable.  For the dihedral groups the box should have an edge length of 2,
! for tetrahedral 1/3, and for octahedral sqrt(2)-1.
  rodx = float(enl%nx) / maxRFZdis(FZtype+1)
  rody = float(enl%nx) / maxRFZdis(FZtype+1)
  if ((FZtype.eq.1).or.(FZtype.eq.2)) then
    rodz = float(enl%nz) / tan(cPi*0.5/float(FZorder))
  else
    rodz = float(enl%nz) / maxRFZdis(FZtype+1) 
  end if
end if

!==============
! stereographic
!==============
if (enl%stereographic.ne.0) then
  if (enl%mrcmode.eq.'off') then
    call initFiles(EMsoft, enl, PoVst, SO, 'st', outname, dataunit4 )
  end if
! create the rendering volume
  allocate(spvol(-enl%nx:enl%nx,-enl%ny:enl%ny,-enl%nz:enl%nz),stat=istat)
  spvol = 0.0
! generate the x/y/z scaling factors
! these should be set so that the render box has the outer dimensions of the stereographic sphere,
! which has unit radius.
  spdx = float(enl%nx)
  spdy = float(enl%nx)
  spdz = float(enl%nz)
end if

!==============
! stereographic
!==============
if (enl%eulerspace.ne.0) then
  if (enl%mrcmode.eq.'off') then
! Euler space has a slightly different eye position so we need to redefine the locationline and skyline strings
    fname = trim(outname)//'-eu.pov'
    eyepos(1:3) = (/ 1.50, 1.60, 1.60 /)
    PoVeu = PoVRay_T( EMsoft, fname, dunit=dataunit5, nmlfile=EMsoft%nmldeffile, eyepos=eyepos,  distance = dble(enl%distance) )
    call initFiles(EMsoft, enl, PoVeu, SO, 'eu', outname, dataunit5 )
  end if
! create the rendering volume
  allocate(euvol(1:2*enl%nx+1,1:2*enl%ny+1,1:2*enl%nz+1),stat=istat)
  euvol = 0.0
! generate the x/y/z scaling factors
! these should be set so that the render box has the outer dimensions of the primary Euler cell
! which has unit radius.
  eudx = float(2*enl%nx+1) / (2.0*cPi)
  eudy = float(2*enl%ny+1) / cPi
  eudz = float(2*enl%nz+1) / (2.0*cPi)
end if

! data file of orientations, convert the orientations to the RFZ
fname = EMsoft%generateFilePath('EMdatapathname', enl%anglefile)
call SO%getOrientationsfromFile(fname)
! next, we reduce these orientations to the RFZ or MacKenzie cell
if (enl%MacKenzieCell.eq.1) then 
  call SO%ReducelisttoRFZ(Pm)
else
  call SO%ReducelisttoMFZ(SG)  ! this requires the full crystal point group
end if

FZtmp => SO%getListHead('FZ')          ! point to the top of the list
numpoints = SO%getListCount('FZ')

pointloop: do ix = 1,numpoints 
  ro = FZtmp%rod 
! are we drawing a point or splatting it into a volume ? 
  if ((trim(enl%df3file).eq.'undefined').and.(enl%mrcmode.eq.'off')) then 
! we have the point, so now transform it to all alternative representations and draw the 
! corresponding spheres 
! cubochoric
    if (enl%cubochoric.ne.0) then
      cu = ro%rc()
      xyz = cu%c_copy()
      write (dataunit,"('sphere { <',2(F14.6,','),F14.6,'>,',F6.4,'  }')") xyz(1:3), enl%sphrad
    end if
! homochoric
    if (enl%homochoric.ne.0) then
      ho = ro%rh()
      xyz = ho%h_copy()
      write (dataunit2,"('sphere { <',2(F14.6,','),F14.6,'>,',F6.4,' }')") xyz(1:3), enl%sphrad
    end if
! rodrigues
    if (enl%rodrigues.ne.0) then
      xyz4 = ro%r_copy()
      if (xyz4(4).le.1000.0) then 
        write (dataunit3,"('sphere { <',2(F20.6,','),F20.6,'>,',F6.4,' }')") xyz4(1:3)*xyz4(4), enl%sphrad
      end if
    end if
! stereographic
    if (enl%stereographic.ne.0) then
      st = ro%rs()
      xyz = st%s_copy()
      write (dataunit4,"('sphere { <',2(F14.6,','),F14.6,'>,',F6.4,' }')") xyz(1:3), enl%sphrad
    end if
! Euler box
    if (enl%eulerspace.ne.0) then
      eu = ro%re()
      xyz = eu%e_copy()
      write (dataunit5,"('sphere { <',2(F14.6,','),F14.6,'>,',F6.4,' }')") xyz(1:3) - sh(1:3), enl%sphrad
    end if

  else  ! we're filling up the rendering volume(s)

! cubochoric rendering volume
    if (enl%cubochoric.ne.0) then
      cu = ro%rc()
      xyz = cu%c_copy()
      grid3 = trilinear_splat( sngl(xyz(1:3)), (/ cudx, cudy, cudz /), init=.FALSE.)
      ixyz(1:3) = nint(xyz(1:3) * (/ cudx, cudy, cudz /))
      if ((abs(ixyz(1)).lt.enl%nx).and.(abs(ixyz(2)).lt.enl%ny).and.(abs(ixyz(3)).lt.enl%nz) ) then
        cuvol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) = &
          cuvol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) + grid3
      end if
    end if

! homochoric rendering volume
    if (enl%homochoric.ne.0) then
      ho = ro%rh()
      xyz = ho%h_copy()
      grid3 = trilinear_splat( sngl(xyz(1:3)), (/ hodx, hody, hodz /), init=.FALSE.)
      ixyz(1:3) = nint(xyz(1:3) * (/ hodx, hody, hodz /))
      if ((abs(ixyz(1)).lt.enl%nx).and.(abs(ixyz(2)).lt.enl%ny).and.(abs(ixyz(3)).lt.enl%nz) ) then
        hovol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) = &
          hovol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) + grid3
      end if
    end if

! rodrigues rendering volume
    if (enl%rodrigues.ne.0) then
      xyz4 = ro%r_copy()
      grid3 = trilinear_splat( sngl(xyz4(1:3)*xyz4(4)), (/ rodx, rody, rodz /), init=.FALSE.)
      ixyz(1:3) = nint((xyz4(1:3)*xyz4(4)) * (/ rodx, rody, rodz /))
      if ((abs(ixyz(1)).lt.enl%nx).and.(abs(ixyz(2)).lt.enl%ny).and.(abs(ixyz(3)).lt.enl%nz) ) then
        rovol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) = &
          rovol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) + grid3
      end if
    end if

! stereographic rendering volume
    if (enl%stereographic.ne.0) then
      st = ro%rs()
      xyz = st%s_copy()
      grid3 = trilinear_splat( sngl(xyz), (/ spdx, spdy, spdz /), init=.FALSE.)
      ixyz(1:3) = nint(xyz(1:3) * (/ spdx, spdy, spdz /))
      ! in the following test, we are potentially eliminating a few points that lie on
      ! the surface of the volume array...
      if ((abs(ixyz(1)).lt.enl%nx).and.(abs(ixyz(2)).lt.enl%ny).and.(abs(ixyz(3)).lt.enl%nz) ) then
        spvol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) = &
          spvol(ixyz(1)-1:ixyz(1)+1,ixyz(2)-1:ixyz(2)+1,ixyz(3)-1:ixyz(3)+1) + grid3
      end if
    end if

! Euler primary cell rendering volume  [needs to be updated with trilinear splatting]
    if (enl%eulerspace.ne.0) then
      eu = ro%re()
      xyz = eu%e_copy()
      xyz(1) = mod(xyz(1)+10.0*cPi,2.0*cPi)
      xyz(2) = mod(xyz(2)+5.0*cPi,cPi)
      xyz(3) = mod(xyz(3)+10.0*cPi,2.0*cPi)
      ixyz(1:3) = nint(xyz(1:3) * (/ eudx, eudy, eudz /))
      if ((ixyz(1).le.2*enl%nx).and.(ixyz(2).le.2*enl%ny).and.(ixyz(3).le.2*enl%nz) ) then
        euvol(ixyz(1)+1,ixyz(2)+1,ixyz(3)+1) = euvol(ixyz(1)+1,ixyz(2)+1,ixyz(3)+1) + 1.0
      end if
    end if
  end if 
  FZtmp => FZtmp%next
end do pointloop

! final stage of the program so close the PoVray files  
if (enl%mrcmode.eq.'off') then 
 if (trim(enl%df3file).eq.'undefined') then ! we're drawing spheres, so close the union and add color information
  write(rgbstring,"(F8.6,',',F8.6,',',F8.6)") enl%rgb(1:3)
  colorstring = 'material { texture { pigment { rgb <'//trim(rgbstring)//'> filter 0.95 }'
  if (enl%cubochoric.ne.0) then
    write (dataunit,"(A)") trim(colorstring)
    write (dataunit,"(' finish { diffuse 0.6, 0.6 brilliance 1.0 } } } } ')")
    write (dataunit,"(A)") 'background { color rgb <0.9, 0.9, 0.9> }'
    call PoVcu%closeFile() 
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-cu.pov')
  end if

  if (enl%homochoric.ne.0) then
    write (dataunit2,"(A)") trim(colorstring)
    write (dataunit2,"(' finish { diffuse 0.6, 0.6 brilliance 1.0 }  } } }')")
    write (dataunit2,"(A)") 'background { color rgb <0.9, 0.9, 0.9> }'
    call PoVho%closeFile() 
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-ho.pov')
  end if

  if (enl%rodrigues.ne.0) then
    write (dataunit3,"(A)") trim(colorstring)
    write (dataunit3,"('finish { diffuse 0.6, 0.6 brilliance 1.0 }  } } }')")
    write (dataunit3,"(A)") 'background { color rgb <0.9, 0.9, 0.9> }'
    call PoVro%closeFile() 
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-ro.pov')
  end if

  if (enl%stereographic.ne.0) then
    write (dataunit4,"(A)") trim(colorstring)
    write (dataunit4,"(' finish { diffuse 0.6, 0.6 brilliance 1.0 }  } } }')")
    write (dataunit4,"(A)") 'background { color rgb <0.9, 0.9, 0.9> }'
    call PoVst%closeFile() 
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-sp.pov')
  end if

  if (enl%eulerspace.ne.0) then
    write (dataunit5,"(A)") trim(colorstring)
    write (dataunit5,"(' finish { diffuse 0.6, 0.6 brilliance 1.0 }  } } }')")
    write (dataunit5,"(A)") 'background { color rgb <0.9, 0.9, 0.9> }'
    call PoVeu%closeFile() 
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-eu.pov')
  end if

 else ! we're writing the rendering volume to a DF3 file and closing the povray file

  if (enl%cubochoric.ne.0) then
! output the final rendering commands
    write (dataunit,"(A)") 'background { color rgb <0.2, 0.2, 0.2> }'
    write (dataunit,"('object { renderbox translate <-0.5, -0.5, -0.5>')")
    write (dataunit,"(' scale <',F10.6,',',F10.6,',',F10.6,' > }')") 2.0 * (/ acubo, acubo, acubo /)
! and close the file
    close(UNIT=dataunit,STATUS='keep')
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-cu.pov')
    df3name = EMsoft%generateFilePath('EMdatapathname', enl%df3file//'-cu.df3')
    call PoVcu%write_DF3file(df3name, cuvol, (/ enl%nx, enl%ny, enl%nz /), enl%scalingmode)
  end if

  if (enl%homochoric.ne.0) then
! output the final rendering commands
    write (dataunit2,"(A)") 'background { color rgb <0.2, 0.2, 0.2> }'
    write (dataunit2,"('object { renderbox translate <-0.5, -0.5, -0.5>')")
    write (dataunit2,"(' scale <',F10.6,',',F10.6,',',F10.6,' > }')") 2.0 * (/ ahomo, ahomo, ahomo /)
! and close the file
    close(UNIT=dataunit2,STATUS='keep')
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-ho.pov')
    df3name = EMsoft%generateFilePath('EMdatapathname', enl%df3file//'-ho.df3')
    call PoVho%write_DF3file(df3name, hovol, (/ enl%nx, enl%ny, enl%nz /), enl%scalingmode)
  end if

  if (enl%rodrigues.ne.0) then
! output the final rendering commands
    write (dataunit3,"(A)") 'background { color rgb <0.2, 0.2, 0.2> }'
    write (dataunit3,"('object { renderbox translate <-0.5, -0.5, -0.5>')")
    if (enl%overridepgnum.eq.0) then
      write (dataunit3,"(' scale <',F10.6,',',F10.6,',',F10.6,' > }')") scalefactors(1:3,FZtype+1)
    else
      write (dataunit3,"(' scale <',F10.6,',',F10.6,',',F10.6,' > }')") scalefactors(1:3,FZtype_override+1)
    end if
! and close the file
    close(UNIT=dataunit3,STATUS='keep')
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-ro.pov')
    df3name = EMsoft%generateFilePath('EMdatapathname', enl%df3file//'-ro.df3')
    call PoVro%write_DF3file(df3name, rovol, (/ enl%nx, enl%ny, enl%nz /), enl%scalingmode)
  end if

  if (enl%stereographic.ne.0) then
! output the final rendering commands
    write (dataunit4,"(A)") 'background { color rgb <0.2, 0.2, 0.2> }'
    write (dataunit4,"('object { renderbox translate <-0.5, -0.5, -0.5>')")
    write (dataunit4,"(' scale <',F10.6,',',F10.6,',',F10.6,' > }')") (/ 2.0, 2.0, 2.0 /)
! and close the file
    close(UNIT=dataunit4,STATUS='keep')
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-sp.pov')
    df3name = EMsoft%generateFilePath('EMdatapathname', enl%df3file//'-sp.df3')
    call PoVst%write_DF3file(df3name, spvol, (/ enl%nx, enl%ny, enl%nz /), enl%scalingmode)
  end if

  if (enl%eulerspace.ne.0) then
! output the final rendering commands
    write (dataunit5,"(A)") 'background { color rgb <0.2, 0.2, 0.2> }'
    write (dataunit5,"('object { renderbox scale < 6.2831855, 3.1415927, 6.2831855 > ')")
    write (dataunit5,"(' translate <  -3.1415927, -1.5707964, -3.1415927 > } ')")
! and close the file
    close(UNIT=dataunit5,STATUS='keep')
    call Message%printMessage('PoVray rendering script stored in '//trim(outname)//'-eu.pov')
    df3name = EMsoft%generateFilePath('EMdatapathname', enl%df3file//'-eu.df3')
    call PoVeu%write_DF3file(df3name, euvol, (/ enl%nx, enl%ny, enl%nz /), enl%scalingmode)
  end if
 end if
else  ! we're creating an .mrc file, so we do not need any of the povray commands...
 !write(*,*) 'starting creation of .mrc file'
! parameters that are generic to all the volumes
  numx = 2*enl%nx+1
  numy = 2*enl%ny+1
  numz = 2*enl%nz+1
  allocate(volume(numx,numy,numz))
  MRCheader%nx = numx
  MRCheader%ny = numy
  MRCheader%nz = numz
  MRCheader%mode = 2    ! for floating point output
  MRCheader%mx = numx
  MRCheader%my = numy
  MRCheader%mz = numz
  MRCheader%xlen = numx
  MRCheader%ylen = numy
  MRCheader%zlen = numz
  do i=1,numz
    FEIheaders(i)%b_tilt = 0.0
    FEIheaders(i)%defocus = 0.0
    FEIheaders(i)%pixelsize = 1.0e-9
    FEIheaders(i)%magnification = 1000.0
    FEIheaders(i)%voltage = 0.0
  end do
  allocate(psum(numz))
! write (*,*) 'dimensions :',numx, numy, numz

  if (enl%cubochoric.ne.0) then
! copy the cubochoric array into the volume array
    do ix = -enl%nx,enl%nx
     do iy = -enl%ny,enl%ny
      do iz = -enl%nz,enl%nz
        volume(ix+enl%nx+1, iy+enl%ny+1, iz+enl%nz+1) = dble(cuvol(ix,iy,iz))
      end do
     end do 
    end do    
! parameters specific to this volume
    psum = sum(sum(volume,1),1)
    do iz=1,numz
      FEIheaders(iz)%mean_int = psum(iz)/float(numx)/float(numy)
    end do
    MRCheader%amin = minval(volume)
    MRCheader%amax = maxval(volume)
    MRCheader%amean = sum(volume)/float(numx)/float(numy)/float(numz)
! set the filename
    mrcname = EMsoft%generateFilePath('EMdatapathname', enl%mrcfile)//'-cu.mrc'
! and write the volume to file
    call MRC%setMRCFileName(mrcname)
    call MRC%write_3Dvolume(volume,verbose=.TRUE.) 
  end if


  if (enl%homochoric.ne.0) then
! copy the homochoric array into the volume array
    do ix = -enl%nx,enl%nx
     do iy = -enl%ny,enl%ny
      do iz = -enl%nz,enl%nz
        volume(ix+enl%nx+1, iy+enl%ny+1, iz+enl%nz+1) = dble(hovol(ix,iy,iz))
      end do
     end do 
    end do    
! parameters specific to this volume
    psum = sum(sum(volume,1),1)
    do iz=0,numz-1
      FEIheaders(iz)%mean_int = psum(iz)/float(numx)/float(numy)
    end do
    MRCheader%amin = minval(volume)
    MRCheader%amax = maxval(volume)
    MRCheader%amean = sum(volume)/float(numx)/float(numy)/float(numz)
! set the filename
    mrcname = EMsoft%generateFilePath('EMdatapathname', enl%mrcfile)//'-ho.mrc'
! and write the volume to file
    call MRC%setMRCFileName(mrcname)
    call MRC%write_3Dvolume(volume,verbose=.TRUE.) 
  end if


  if (enl%rodrigues.ne.0) then
! copy the rodrigues array into the volume array
    do ix = -enl%nx,enl%nx
     do iy = -enl%ny,enl%ny
      do iz = -enl%nz,enl%nz
        volume(ix+enl%nx+1, iy+enl%ny+1, iz+enl%nz+1) = dble(rovol(ix,iy,iz))
      end do
     end do 
    end do    
! parameters specific to this volume
    psum = sum(sum(volume,1),1)
    do iz=1,numz
      FEIheaders(iz)%mean_int = psum(iz)/float(numx)/float(numy)
    end do
    MRCheader%amin = minval(volume)
    MRCheader%amax = maxval(volume)
    MRCheader%amean = sum(volume)/float(numx)/float(numy)/float(numz)
! set the filename
    mrcname = EMsoft%generateFilePath('EMdatapathname', enl%mrcfile)//'-ro.mrc'
! and write the volume to file
    call MRC%setMRCFileName(mrcname)
    call MRC%write_3Dvolume(volume,verbose=.TRUE.) 
  end if


  if (enl%stereographic.ne.0) then
  ! write (*,*) 'starting volume array copy'
! copy the stereographic array into the volume array
    do ix = -enl%nx,enl%nx
     do iy = -enl%ny,enl%ny
      do iz = -enl%nz,enl%nz
        volume(ix+enl%nx+1, iy+enl%ny+1, iz+enl%nz+1) = dble(spvol(ix,iy,iz))
      end do
     end do 
    end do    
  ! write (*,*) '  --> done'
! parameters specific to this volume
    psum = sum(sum(volume,1),1)
    do iz=1,numz
      FEIheaders(iz)%mean_int = psum(iz)/float(numx)/float(numy)
    end do
    MRCheader%amin = minval(volume)
    MRCheader%amax = maxval(volume)
    MRCheader%amean = sum(volume)/float(numx)/float(numy)/float(numz)
  ! write(*,*) MRCheader%amin, MRCheader%amax, MRCheader%amean
! set the filename
    mrcname = EMsoft%generateFilePath('EMdatapathname', enl%mrcfile)//'-sp.mrc'
! and write the volume to file
    call MRC%setMRCFileName(mrcname)
    call MRC%write_3Dvolume(volume,verbose=.TRUE.) 
  end if


  if (enl%eulerspace.ne.0) then
! copy the eulerspace array into the volume array
    volume = dble(euvol)
! parameters specific to this volume
    psum = sum(sum(volume,1),1)
    do iz=1,numz
      FEIheaders(iz)%mean_int = psum(iz)/float(numx)/float(numy)
    end do
    MRCheader%amin = minval(volume)
    MRCheader%amax = maxval(volume)
    MRCheader%amean = sum(volume)/float(numx)/float(numy)/float(numz)
! set the filename
    mrcname = EMsoft%generateFilePath('EMdatapathname', enl%mrcfile)//'-eu.mrc'
! and write the volume to file
    call MRC%setMRCFileName(mrcname)
    call MRC%write_3Dvolume(volume,verbose=.TRUE.) 
  end if

end if

end associate

end subroutine OrientationViz_


!--------------------------------------------------------------------------
subroutine initFiles(EMsoft, enl, PoV, SO, rep, outname, dunit)
!! author: MDG 
!! version: 1.0 
!! date: 03/27/20
!!
!! perform the file initializations

use mod_EMsoft
use mod_io
use mod_povray
use mod_so3

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT)                      :: EMsoft 
type(OrientationVizNameListType),INTENT(INOUT)    :: enl 
type(PoVRay_T),INTENT(INOUT)                      :: PoV
character(2),INTENT(IN)                           :: rep 
type(so3_T),INTENT(INOUT)                         :: SO
character(fnlen),INTENT(IN)                       :: outname
integer(kind=irg),INTENT(IN)                      :: dunit

type(IO_T)                                        :: Message
character(fnlen)                                  :: fname, DF3name 
logical                                           :: drawMK = .FALSE.
integer(kind=irg)                                 :: dFZ 

drawMK = SO%getMK()

select case(rep)
  case('cu')
    dFZ = 1 
  case('ho')
    dFZ = 2 
  case('st')
    dFZ = 3
  case('ro')
    dFZ = 4 
  case('eu')
    dFZ = 5
end select

fname = trim(outname)//'-'//rep//'.pov'
call Message%printMessage('opening '//trim(fname))
if (rep.ne.'eu') then
  PoV = PoVRay_T( EMsoft, fname, dunit=dunit, nmlfile=EMsoft%nmldeffile, distance = dble(enl%distance) )
end if 
if (trim(enl%df3file).eq.'undefined') then 
! we're just going to draw a bunch of spheres, so put them together in a PoVRay union
  call PoV%drawFZ(SO, dFZ, 0.005D0)
! open the union here
  write (dunit,"('union { ')")
else
! insert code to read in a 3D Density File (df3) containing the object to be rendered
  df3name = EMsoft%generateFilePath('EMdatapathname',enl%df3file)//'-'//rep//'.df3'
  if (enl%scalingmode.eq.'lev') then
    call PoV%declare_DF3file(df3name,levelset=.TRUE.)
  else
    call PoV%declare_DF3file(df3name)
  end if
  call PoV%drawFZ(SO, dFZ, 0.005D0)
end if

end subroutine initFiles

end module mod_OrientationViz