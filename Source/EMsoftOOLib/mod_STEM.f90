! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_STEM
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/13/24
  !!
  !! STEM class definition

use mod_kinds
use mod_global

IMPLICIT NONE 

type STEMGeometryNameListType
  integer(kind=irg)               :: numberofsvalues
  integer(kind=irg)               :: numCL
  real(kind=sgl)                  :: BFradius
  real(kind=sgl)                  :: ADFinnerradius
  real(kind=sgl)                  :: ADFouterradius
  real(kind=sgl)                  :: kt
  real(kind=sgl)                  :: beamconvergence
  real(kind=sgl)                  :: diffaprad
  real(kind=sgl)                  :: diffapcenter
  real(kind=sgl)                  :: CLarray(20)
  character(2)                    :: geometry
end type STEMGeometryNameListType

! class definition
type, public :: STEM_T
private 
  type(STEMGeometryNameListType)  :: nml 
  integer(kind=irg)               :: numk
  real(kind=sgl)                  :: BFmrad
  real(kind=sgl)                  :: ADFimrad
  real(kind=sgl)                  :: ADFomrad
  real(kind=sgl)                  :: diffapmrad
  real(kind=sgl)                  :: diffapmcenter
  logical,allocatable             :: ZABFweightsarray(:,:,:)
  logical,allocatable             :: ZAADFweightsarray(:,:,:) ! only used for the zone axis case
  real(kind=sgl),allocatable      :: sgarray(:,:)
  real(kind=sgl),allocatable      :: BFweightsarray(:,:,:)
  real(kind=sgl),allocatable      :: ADFweightsarray(:,:,:)   ! only used for the systematic row case

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: init_STEM_
  procedure, pass(self) :: init_STEM_ZA_
  procedure, pass(self) :: setnumberofsvalues_
  procedure, pass(self) :: getnumberofsvalues_
  procedure, pass(self) :: setnumCL_
  procedure, pass(self) :: getnumCL_
  procedure, pass(self) :: setnumk_
  procedure, pass(self) :: getnumk_
  procedure, pass(self) :: setBFradius_
  procedure, pass(self) :: getBFradius_
  procedure, pass(self) :: setADFinnerradius_
  procedure, pass(self) :: getADFinnerradius_
  procedure, pass(self) :: setADFouterradius_
  procedure, pass(self) :: getADFouterradius_
  procedure, pass(self) :: setkt_
  procedure, pass(self) :: getkt_
  procedure, pass(self) :: setbeamconvergence_
  procedure, pass(self) :: getbeamconvergence_
  procedure, pass(self) :: setdiffaprad_
  procedure, pass(self) :: getdiffaprad_
  procedure, pass(self) :: setdiffapcenter_
  procedure, pass(self) :: getdiffapcenter_
  procedure, pass(self) :: setCLarray_
  procedure, pass(self) :: getCLarray_
  procedure, pass(self) :: setgeometry_
  procedure, pass(self) :: getgeometry_            
  procedure, pass(self) :: getsgarray_
  procedure, pass(self) :: getBFweightsarray_
  procedure, pass(self) :: getADFweightsarray_
  procedure, pass(self) :: getZABFweightsarray_
  procedure, pass(self) :: getZAADFweightsarray_

  generic, public :: readNameList => readNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: getNameList => getNameList_
  generic, public :: init_STEM => init_STEM_
  generic, public :: init_STEM_ZA => init_STEM_ZA_
  generic, public :: setnumberofsvalues => setnumberofsvalues_
  generic, public :: getnumberofsvalues => getnumberofsvalues_
  generic, public :: setnumCL => setnumCL_
  generic, public :: getnumCL => getnumCL_
  generic, public :: setnumk => setnumk_
  generic, public :: getnumk => getnumk_
  generic, public :: setBFradius => setBFradius_
  generic, public :: getBFradius => getBFradius_
  generic, public :: setADFinnerradius => setADFinnerradius_
  generic, public :: getADFinnerradius => getADFinnerradius_
  generic, public :: setADFouterradius => setADFouterradius_
  generic, public :: getADFouterradius => getADFouterradius_
  generic, public :: setkt => setkt_
  generic, public :: getkt => getkt_
  generic, public :: setbeamconvergence => setbeamconvergence_
  generic, public :: getbeamconvergence => getbeamconvergence_
  generic, public :: setdiffaprad => setdiffaprad_
  generic, public :: getdiffaprad => getdiffaprad_
  generic, public :: setdiffapcenter => setdiffapcenter_
  generic, public :: getdiffapcenter => getdiffapcenter_
  generic, public :: setCLarray => setCLarray_
  generic, public :: getCLarray => getCLarray_
  generic, public :: setgeometry => setgeometry_
  generic, public :: getgeometry => getgeometry_
  generic, public :: getsgarray => getsgarray_
  generic, public :: getBFweightsarray => getBFweightsarray_
  generic, public :: getADFweightsarray => getADFweightsarray_
  generic, public :: getZABFweightsarray => getZABFweightsarray_
  generic, public :: getZAADFweightsarray => getZAADFweightsarray_
end type STEM_T

! the constructor routine for this class 
interface STEM_T
  module procedure STEM_constructor
end interface STEM_T

contains

!--------------------------------------------------------------------------
type(STEM_T) function STEM_constructor( nmlfile ) result(STEM)
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! constructor for the STEM_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call STEM%readNameList_(nmlfile)

end function STEM_constructor

!--------------------------------------------------------------------------
subroutine STEM_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! destructor for the STEM_T Class
 
IMPLICIT NONE

type(STEM_T), INTENT(INOUT)  :: self 

call reportDestructor('STEM_T')

end subroutine STEM_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! read the namelist from an nml file for the SRdefect_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(STEM_T), INTENT(INOUT)         :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: numberofsvalues 
integer(kind=irg)                    :: numCL
real(kind=sgl)                       :: BFradius
real(kind=sgl)                       :: ADFinnerradius
real(kind=sgl)                       :: ADFouterradius
real(kind=sgl)                       :: kt
real(kind=sgl)                       :: beamconvergence
real(kind=sgl)                       :: diffaprad
real(kind=sgl)                       :: diffapcenter
real(kind=sgl)                       :: CLarray(20)
character(2)                         :: geometry

! define the IO namelist to facilitate passing variables to the program.
namelist /STEMGeometrydata/ numberofsvalues, numCL, BFradius, ADFinnerradius, ADFouterradius, kt, &
                            beamconvergence, CLarray, geometry, diffaprad, diffapcenter

! SET INPUT PARAMETERS TO DEFAULT VALUES (EXCEPT XTALNAME, WHICH MUST BE PRESENT)
numberofsvalues = 21
numCL = 2
BFradius = 3.5
ADFinnerradius = 3.5
ADFouterradius = 10.0
kt = 1.0
beamconvergence = 8.2
diffaprad = 0.0
diffapcenter = 0.0
CLarray = 0.0
CLarray(1) = 150.0
CLarray(2) = 250.0
geometry = 'ZA'

if (.not.skipread) then
! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=STEMGeometrydata)
  close(UNIT=dataunit,STATUS='keep')
end if

if ( (trim(geometry).eq.'SR').and.(mod(numberofsvalues,2).eq.0) ) numberofsvalues=numberofsvalues+1
self%nml%numberofsvalues = numberofsvalues 
self%nml%numCL = numCL
self%nml%BFradius = BFradius 
self%nml%ADFinnerradius = ADFinnerradius
self%nml%ADFouterradius = ADFouterradius
self%nml%kt = kt 
self%nml%beamconvergence = beamconvergence 
self%nml%diffaprad = diffaprad
self%nml%diffapcenter = diffapcenter
self%nml%CLarray = CLarray 
self%nml%geometry = geometry 

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! pass the namelist for the SRdefect_T Class to the calling program

IMPLICIT NONE 

class(STEM_T), INTENT(INOUT)              :: self
type(STEMGeometryNameListType)            :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 02/13/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)            :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 2, n_real = 7
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%numberofsvalues, enl%numCL /)
intlist(1) = 'numberofsvalues'
intlist(2) = 'numCL'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%BFradius, enl%ADFinnerradius, enl%ADFouterradius, enl%kt, enl%beamconvergence, enl%diffaprad, &
             enl%diffapcenter /)
reallist(1) = 'BFradius'
reallist(2) = 'ADFinnerradius'
reallist(3) = 'ADFouterradius'
reallist(4) = 'kt'
reallist(5) = 'beamconvergence'
reallist(6) = 'diffaprad'
reallist(7) = 'diffapcenter'
call HDF%writeNMLreals(io_real, reallist, n_real)

! a 20-vector
dataset = 'CLarray'
hdferr = HDF%writeDatasetFloatArray(dataset, enl%CLarray, 20)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create CLarray dataset', hdferr)

! write all the strings
dataset = 'geometry'
line2(1) = trim(enl%geometry)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create geometry dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine setnumberofsvalues_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumberofsvalues_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set numberofsvalues in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numberofsvalues = inp

end subroutine setnumberofsvalues_

!--------------------------------------------------------------------------
function getnumberofsvalues_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumberofsvalues_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get numberofsvalues from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numberofsvalues

end function getnumberofsvalues_

!--------------------------------------------------------------------------
subroutine setnumCL_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumCL_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set numCL in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%nml%numCL = inp

end subroutine setnumCL_

!--------------------------------------------------------------------------
function getnumCL_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumCL_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get numCL from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%nml%numCL

end function getnumCL_

!--------------------------------------------------------------------------
subroutine setnumk_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setnumk_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set numk in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%numk = inp

end subroutine setnumk_

!--------------------------------------------------------------------------
function getnumk_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getnumk_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get numk from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%numk

end function getnumk_

!--------------------------------------------------------------------------
subroutine setBFradius_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setBFradius_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set BFradius in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%BFradius = inp

end subroutine setBFradius_

!--------------------------------------------------------------------------
function getBFradius_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getBFradius_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get BFradius from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%BFradius

end function getBFradius_

!--------------------------------------------------------------------------
subroutine setADFinnerradius_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setADFinnerradius_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set ADFinnerradius in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%ADFinnerradius = inp

end subroutine setADFinnerradius_

!--------------------------------------------------------------------------
function getADFinnerradius_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getADFinnerradius_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get ADFinnerradius from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%ADFinnerradius

end function getADFinnerradius_

!--------------------------------------------------------------------------
subroutine setADFouterradius_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setADFouterradius_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set ADFouterradius in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%ADFouterradius = inp

end subroutine setADFouterradius_

!--------------------------------------------------------------------------
function getADFouterradius_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getADFouterradius_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get ADFouterradius from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%ADFouterradius

end function getADFouterradius_

!--------------------------------------------------------------------------
subroutine setkt_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setkt_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set kt in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%kt = inp

end subroutine setkt_

!--------------------------------------------------------------------------
function getkt_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getkt_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get kt from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%kt

end function getkt_

!--------------------------------------------------------------------------
subroutine setbeamconvergence_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setbeamconvergence_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set beamconvergence in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%beamconvergence = inp

end subroutine setbeamconvergence_

!--------------------------------------------------------------------------
function getbeamconvergence_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getbeamconvergence_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get beamconvergence from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%beamconvergence

end function getbeamconvergence_

!--------------------------------------------------------------------------
subroutine setdiffaprad_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdiffaprad_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set diffaprad in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%diffaprad = inp

end subroutine setdiffaprad_

!--------------------------------------------------------------------------
function getdiffaprad_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdiffaprad_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get diffaprad from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%diffaprad

end function getdiffaprad_

!--------------------------------------------------------------------------
subroutine setdiffapcenter_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setdiffapcenter_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set diffapcenter in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp

self%nml%diffapcenter = inp

end subroutine setdiffapcenter_

!--------------------------------------------------------------------------
function getdiffapcenter_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getdiffapcenter_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get diffapcenter from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out

out = self%nml%diffapcenter

end function getdiffapcenter_

!--------------------------------------------------------------------------
subroutine setCLarray_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setCLarray_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set CLarray in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)       :: inp(20)

self%nml%CLarray = inp

end subroutine setCLarray_

!--------------------------------------------------------------------------
function getCLarray_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getCLarray_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get CLarray from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl)                   :: out(20)

out = self%nml%CLarray

end function getCLarray_

!--------------------------------------------------------------------------
subroutine setgeometry_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setgeometry_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! set geometry in the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
character(2), INTENT(IN)       :: inp

self%nml%geometry = trim(inp)

end subroutine setgeometry_

!--------------------------------------------------------------------------
function getgeometry_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getgeometry_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get geometry from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
character(2)                   :: out

out = trim(self%nml%geometry)

end function getgeometry_

!--------------------------------------------------------------------------
function getsgarray_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsgarray_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get sgarray from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl),allocatable       :: out(:,:)

integer(kind=irg)                :: sz(2) 

sz = shape(self%sgarray)
allocate(out(sz(1), sz(2)))
out = self%sgarray

end function getsgarray_

!--------------------------------------------------------------------------
function getBFweightsarray_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getBFweightsarray_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get BFweightsarray from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl),allocatable       :: out(:,:,:)

integer(kind=irg)                :: sz(3) 

sz = shape(self%BFweightsarray)
allocate(out(sz(1), sz(2), sz(3)))
out = self%BFweightsarray

end function getBFweightsarray_

!--------------------------------------------------------------------------
function getADFweightsarray_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getADFweightsarray_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get ADFweightsarray from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
real(kind=sgl),allocatable       :: out(:,:,:)

integer(kind=irg)                :: sz(3) 

sz = shape(self%ADFweightsarray)
allocate(out(sz(1), sz(2), sz(3)))
out = self%ADFweightsarray

end function getADFweightsarray_

!--------------------------------------------------------------------------
function getZABFweightsarray_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getZABFweightsarray_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get ZABFweightsarray from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg),allocatable    :: out(:,:,:)

integer(kind=irg)                :: sz(3) 

sz = shape(self%ZABFweightsarray)
allocate(out(sz(1), sz(2), sz(3)))
out = 0
where(self%ZABFweightsarray)
  out = 1
end where

end function getZABFweightsarray_

!--------------------------------------------------------------------------
function getZAADFweightsarray_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getZAADFweightsarray_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! get ZAADFweightsarray from the STEM_T class

IMPLICIT NONE

class(STEM_T), INTENT(INOUT)     :: self
integer(kind=irg),allocatable    :: out(:,:,:)

integer(kind=irg)                :: sz(3) 

sz = shape(self%ZAADFweightsarray)
allocate(out(sz(1), sz(2), sz(3)))
out = 0
where(self%ZAADFweightsarray)
  out = 1
end where

end function getZAADFweightsarray_

!--------------------------------------------------------------------------
recursive subroutine init_STEM_(self,cell,Diff,nn,g)
!DEC$ ATTRIBUTES DLLEXPORT :: init_STEM_
!! author: MDG
!! version: 1.0
!! date: 02/13/24
!!
!! initialize a bunch of things ... 

use mod_diffraction
use mod_io 
use mod_crystallography

IMPLICIT NONE

class(STEM_T),INTENT(INOUT)             :: self
type(Cell_T),INTENT(INOUT)              :: cell
type(Diffraction_T),INTENT(INOUT)       :: Diff
integer(kind=irg),INTENT(IN)            :: nn
integer(kind=irg),INTENT(IN)            :: g(3)

type(IO_T)                              :: Message

integer(kind=irg)                       :: i,j,n,ira,jj,k,kk, iCL
real(kind=sgl)                          :: glen, thb, alp, omega_c, omega_min, omega_max,omega,a,b,c,th,dom,p,q,dr,dx,lambda

! these are only used to debug this routine
real(kind=sgl),allocatable              :: thetar(:),outar(:,:,:)
logical                                 :: debug = .FALSE., diffappresent = .FALSE., apinBF=.FALSE. , apinADF = .FALSE.

! this routine initializes the excitation error arrays and the weight-factor arrays for systematic row STEM signals
! we'll assume that the central beam is centered on the BF detector; then we can 
! compute the complete geometry by working in mrad units throughout.

! allocate the excitation error array areal(1..nn,1..self%nml%numberofsvalues)
allocate(self%sgarray(nn,self%nml%numberofsvalues))

! determine the lower and upper bounds of the excitation error for the fundamental reflection G
thb = Diff%CalcDiffAngle( cell, g ) * 0.5  ! Bragg angle in radians
lambda = sngl(Diff%getWaveLength())

! convert k_t to the alp and omega angles (in radians)
glen = cell%CalcLength( float(g),'r' )
alp = -2.0 * self%nml%kt * thb
omega_c = cPi*0.5+alp
omega_min = omega_c - self%nml%beamconvergence/1000.0
omega_max = omega_c + self%nml%beamconvergence/1000.0

! step size
dom = (omega_max - omega_min)/float(self%nml%numberofsvalues-1)

! and for each value in between, compute each reflection's excitation error
ira = (nn+1)/2

do j=1,self%nml%numberofsvalues
! set omega angle
   omega = omega_min+float(j-1)*dom
   do i=1,nn
    n = -ira+i
! excitation error
    self%sgarray(nn+1-i,j) = -n*glen*cos(omega)-(1.0-sqrt(1.0-(n*lambda*glen*sin(omega))**2))/lambda
   end do
end do

! if (debug) then
!   allocate(thetar(25))
!   thetar = 0.0
!   thetar(1:7) = (/ a,b,c,th,dom,dr,dx /)
! end if 

! next, we compute the weightfactors, i.e., how much does each excitation error value contribute
! to the BF or ADF signal?  The weight factor is basically the length of the chord across the overlap
! area of the diffraction disk and the detector, which requires a little bit of math to figure out;
! the math employs the concept of the radical line (see mathworld.com section on circle-circle intersections)

! this computation is carried out in mrad units !
allocate(self%BFweightsarray(nn,self%nml%numberofsvalues,self%nml%numCL), &
         self%ADFweightsarray(nn,self%nml%numberofsvalues,self%nml%numCL))

self%BFweightsarray = 0.0
self%ADFweightsarray = 0.0


outerCLloop: do iCL=1,self%nml%numCL    ! this is the outer loop over the microscope camera lengths (very long loop !!!)

! fist, convert the detector parameters to mrad units
  self%BFmrad = atan(self%nml%BFradius/self%nml%CLarray(iCL))*1000.0
  self%ADFimrad = atan(self%nml%ADFinnerradius/self%nml%CLarray(iCL))*1000.0
  self%ADFomrad = atan(self%nml%ADFouterradius/self%nml%CLarray(iCL))*1000.0
  if (self%nml%diffaprad.ne.0.0) diffappresent = .TRUE.

! then, for each point inside each diffraction disk, determine where it falls with respect to
! the BD and ADF detectors ... Also, look for disk overlaps as they might require amplitudes
! to be added in instead of intensities (for starters, we could just not allow that to happen...)

! rename some variables to short symbols
  a = self%ADFimrad
  b = self%ADFomrad
  c = self%BFmrad
  th = self%nml%beamconvergence
  n = self%nml%numberofsvalues
  if (diffappresent) then
    dr = self%nml%diffaprad
    dx = self%nml%diffapcenter
  end if
  omega_min = -th
  omega_max = th
  dom = 2.0*th/float(n-1)

  if (.not.diffappresent) then   ! there is no diffraction aperture, so compute the regular weight factors
    ! first, do the math for the g=0 disk  (we're dropping common factors of 2 for the weights) 
    i = ira
    do j=(n+1)/2,n
      omega = omega_min+float(j-1)*dom
      if (th.gt.c) then       ! the zero disk is larger than the BF detector, so it (potentially) gives two signals
        if (omega.le.c) self%BFweightsarray(i,j,iCL) = sqrt(c**2-omega**2)
        if (th.ge.a) then    ! there's overlap with the ADF detector
          if (omega.le.a) then  ! the first part needs to have a bit subtracted
            self%ADFweightsarray(i,j,iCL) = sqrt((th**2-omega**2)) - sqrt((a**2-omega**2))
          else   ! the second part does not
            self%ADFweightsarray(i,j,iCL) = sqrt((th**2-omega**2)) 
          end if
        end if
      else  ! the zero disk is smaller than the BF detector, so only a BF signal
       self%BFweightsarray(i,j,iCL) = sqrt((th**2-omega**2))
      end if
    ! then apply symmetry for the other half of the g=0 disk
      if (j.ne.(n+1)/2) then
        jj = n+1 - j
        self%BFweightsarray(i,jj,iCL) = self%BFweightsarray(i,j,iCL)
        self%ADFweightsarray(i,jj,iCL) = self%ADFweightsarray(i,j,iCL)
      end if  
    end do  ! that completes the central disk weight factors

! the other disks are quite a bit more difficult to deal with ... there are a lot of possible cases to consider ...
    do i=ira+1,nn      ! loop over the positive reflections of the systematic row (the rest follows by symmetry)
    ! redefine a couple of parameters
      j = i - ira
      thb = Diff%CalcDiffAngle(cell,j*g)*1000.0  ! diffraction angle in mrad
      omega_min = thb - th
      omega_max = thb + th
    ! only used for debugging
     ! if (debug)  thetar(7+j) = thb
    ! first check if a part of this disk lies inside the BF detector
      if (omega_min.lt.c) then     ! yes, it does, so determine the BF weight factors for this disk
        if (omega_max.le.c) then  ! does it lie completely inside the BF detector?
          do j=1,n   ! yes it does, so compute all weight factors
            omega = omega_min + float(j-1)*dom
            self%BFweightsarray(i,j,iCL) = sqrt(th**2 - (thb-omega)**2)
            self%BFweightsarray(2*ira-i,n+1-j,iCL) = self%BFweightsarray(i,j,iCL)
          end do
        else  ! no, there's some potential overlap with the ADF detector 
          do j=1,n   ! once again, there are a few cases
            omega = omega_min + float(j-1)*dom    ! this is the position
            p = (thb**2-th**2+a**2)*0.5/thb             ! this is the location of the radical line for the ADF detector
            q = (thb**2-th**2+c**2)*0.5/thb             ! this is the location of the radical line for the BF detector
            if (omega.le.q) then   ! this point contributes to the BF detector
              self%BFweightsarray(i,j,iCL) = sqrt(th**2 - (thb-omega)**2)
            end if
            if ((omega.gt.q).and.(omega.le.c)) then   ! this point contributes to the BF detector
              self%BFweightsarray(i,j,iCL) = sqrt(c**2 - omega**2)
            end if
            if ((omega_max.ge.a).and.(omega.ge.p).and.(omega.le.a)) then ! this point contributes to the ADF detector (using radical line position)
              self%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)  - sqrt(a**2 - omega**2)
            end if
             if ((omega_max.ge.a).and.(omega.gt.a)) then ! this point lies on the ADF detector 
              self%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)
            end if
           self%BFweightsarray(2*ira-i,n+1-j,iCL) = self%BFweightsarray(i,j,iCL)       
           self%ADFweightsarray(2*ira-i,n+1-j,iCL) = self%ADFweightsarray(i,j,iCL)
          end do
        end if 
      else    ! no, it does not intersect the BF detector, so this disk can only contribute to the ADF weight factors
    ! once more there are several cases which we'll treat in increasing value of the position...
          do j=1,n 
            omega = omega_min + float(j-1)*dom    ! this is the position
            p = (thb**2-th**2+a**2)*0.5/thb             ! this is the location of the radical line for the inner ADF detector edge
            q = (thb**2-th**2+b**2)*0.5/thb             ! this is the location of the radical line for the outer ADF detector edge
            if ((omega.lt.a).and.(omega.ge.p))  then    ! inside the inner ADF edge, but close enough to contribute
              self%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)  - sqrt(a**2 - omega**2)
            end if
             if ((omega.ge.a).and.(omega_max.le.b)) then ! this point lies on the ADF detector 
              self%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)
            end if
            if ((omega_max.gt.b).and.(omega.le.q)) then   ! this point lies on the ADF detector
              self%ADFweightsarray(i,j,iCL) = sqrt(th**2 - (thb-omega)**2)
            end if
            if ((omega_max.gt.b).and.(omega.gt.q).and.(omega.le.b))  then   ! this point contributes to the ADF detector
              self%ADFweightsarray(i,j,iCL) = sqrt(b**2 - omega**2)
            end if
            self%ADFweightsarray(2*ira-i,n+1-j,iCL) = self%ADFweightsarray(i,j,iCL)
          end do
      end if

    end do
  end if ! end of regular weight factors without a diffraction aperture



  if (diffappresent) then   ! there is a diffraction aperture, so revisit the weight factors.
  ! once again, there are many different cases that need to be addressed...
    
  ! we do not allow for a diffraction aperture that overlaps the boundary between BF and ADF detectors,
  ! nor an aperture that lies entirely beyond the ADF detector
  ! first the BF test
    if ( ((dx-dr).gt.-c).and.((dx+dr).lt.c))  apinBF = .TRUE.

  ! then the ADF detector
    if ( (((dx-dr).gt.-b).and.((dx+dr).lt.-c)) .or.(((dx-dr).gt.a).and.((dx+dr).lt.b)) )  apinADF = .TRUE. 

  ! if the aperture is outside the ADF detector, or it overlaps the space between the detectors, then abort
    if ( .not.apinBF .and. .not.apinADF ) then
      call Message%printError('init_STEM_','Please fix input: Diffraction aperture outside BF detector disk or ADF ring !')
    end if

    if (apinBF) then
    ! figure out which diffraction disk(s) contribute to the BF detector
     do i=1,nn      ! loop over all reflections of the systematic row
    ! redefine a couple of parameters
      j = -(nn-1)/2-1+i
      if (j.ne.0) then 
        thb = (j/abs(j)) * Diff%CalcDiffAngle(cell,j*g)*1000.0  ! diffraction angle in mrad
      else
        thb = 0.0
      end if  
     ! only used for debugging
     ! if (debug)  thetar(7+i) = thb
      omega_min = thb - th
      omega_max = thb + th
    ! check whether or not there is any overlap between this disk and the diffraction aperture opening
      if ((omega_max.lt.(dx-dr)).or.(omega_min.gt.(dx+dr))) then  ! this disks does not fall inside the diffraction aperture 
        self%BFweightsarray(i,1:n,iCL) = 0.0
      else
    ! case 1: dx-dr < omega_min < omega_max < dx+dr
        if ((omega_max.lt.(dx+dr)).and.(omega_min.gt.(dx-dr))) then
          do k=1,n
            omega = omega_min+float(k-1)*dom
            kk = k
            if (j.lt.0) kk = n+1-k 
            self%BFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
         end do 
        end if
    ! case 2: omega_min < dx-dr  < dx+dr < omega_max 
        if ((omega_max.gt.(dx+dr)).and.(omega_min.lt.(dx-dr))) then
          do k=1,n
            omega = omega_min+float(k-1)*dom
            kk = k
            if (j.lt.0) kk = n+1-k 
            if ((omega.gt.dx-dr).and.(omega.lt.dx+dr)) then
               self%BFweightsarray(i,kk,iCL) =  sqrt((dr**2-(omega-dx)**2))
            end if
          end do 
        end if
    ! case 3: omega_min < dx-dr   < omega_max < dx+dr
        if ((omega_min.lt.(dx-dr)).and.(omega_max.lt.(dx+dr))) then
           p = ((dx-thb)**2-dr**2+th**2)/2.0/(dx-thb)+thb
           do k=1,n
             omega = omega_min+float(k-1)*dom
             kk = k
    !         if (j.lt.0) kk = n+1-k 
            if ((omega.gt.dx-dr).and.(omega.le.p)) then
              self%BFweightsarray(i,kk,iCL) =  sqrt(dr**2-(omega-dx)**2)
            end if
            if ((omega.gt.p).and.(omega.le.omega_max)) then
              self%BFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))       
            end if
          end do 
        end if
     ! case 4:  dx-dr   < omega_min < dx+dr < omega_max
        if ((omega_min.gt.(dx-dr)).and.(omega_max.gt.(dx+dr))) then
           p = ((dx-thb)**2-th**2+dr**2)/2.0/(thb-dx) + dx
            do k=1,n
            omega = omega_min+float(k-1)*dom
            kk = k
    !        if (j.lt.0) kk = n+1-k 
            if ((omega.gt.p).and.(omega.le.dx+dr)) then
              self%BFweightsarray(i,kk,iCL) = sqrt((dr**2-(omega-dx)**2)) 
            end if
            if ((omega.gt.omega_min).and.(omega.le.p)) then
               self%BFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
            end if
          end do 
        end if 
      end if
   end do  ! this completes the BF weight factors when a diffraction aperture is present and apinBF=.TRUE.
  end if 


  ! next determine the ADF weight factors in the presence of an aperture
  if (apinADF) then
  ! figure out which diffraction disk(s) contribute to the ADF detector
   do i=1,nn      ! loop over all reflections of the systematic row
  ! redefine a couple of parameters
      j = -(nn-1)/2-1+i
      if (j.ne.0) then 
        thb = (j/abs(j)) * Diff%CalcDiffAngle(cell,j*g)*1000.0  ! diffraction angle in mrad
      else
        thb = 0.0
     end if  
   ! only used for debugging
   ! if (debug)  thetar(7+i) = thb
    omega_min = thb - th
    omega_max = thb + th
  ! check whether or not there is any overlap between this disk and the diffraction aperture opening
    if ((omega_max.lt.(dx-dr)).or.(omega_min.gt.(dx+dr))) then  ! this disks does not fall inside the diffraction aperture 
      self%ADFweightsarray(i,1:n,iCL) = 0.0
    else
  ! case 1: dx-dr < omega_min < omega_max < dx+dr
      if ((omega_max.lt.(dx+dr)).and.(omega_min.gt.(dx-dr))) then
        do k=1,n
          omega = omega_min+float(k-1)*dom
          kk = k
          if (j.lt.0) kk = n+1-k 
          self%ADFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
       end do 
      end if
  ! case 2: omega_min < dx-dr  < dx+dr < omega_max 
      if ((omega_max.gt.(dx+dr)).and.(omega_min.lt.(dx-dr))) then
        do k=1,n
          omega = omega_min+float(k-1)*dom
          kk = k
          if (j.lt.0) kk = n+1-k 
          if ((omega.gt.dx-dr).and.(omega.lt.dx+dr)) then
             self%ADFweightsarray(i,kk,iCL) =  sqrt((dr**2-(omega-dx)**2))
          end if
        end do 
      end if
  ! case 3: omega_min < dx-dr   < omega_max < dx+dr
      if ((omega_min.lt.(dx-dr)).and.(omega_max.lt.(dx+dr))) then
         p = ((dx-thb)**2-dr**2+th**2)/2.0/(dx-thb)+thb
         do k=1,n
           omega = omega_min+float(k-1)*dom
           kk = k
  !         if (j.lt.0) kk = n+1-k 
          if ((omega.gt.dx-dr).and.(omega.le.p)) then
            self%ADFweightsarray(i,kk,iCL) =  sqrt(dr**2-(omega-dx)**2)
          end if
          if ((omega.gt.p).and.(omega.le.omega_max)) then
            self%ADFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))       
          end if
        end do 
      end if
   ! case 4:  dx-dr   < omega_min < dx+dr < omega_max
      if ((omega_min.gt.(dx-dr)).and.(omega_max.gt.(dx+dr))) then
         p = ((dx-thb)**2-th**2+dr**2)/2.0/(thb-dx) + dx
          do k=1,n
          omega = omega_min+float(k-1)*dom
          kk = k
  !        if (j.lt.0) kk = n+1-k 
          if ((omega.gt.p).and.(omega.le.dx+dr)) then
            self%ADFweightsarray(i,kk,iCL) = sqrt((dr**2-(omega-dx)**2)) 
          end if
          if ((omega.gt.omega_min).and.(omega.le.p)) then
             self%ADFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
          end if
        end do 
      end if 
    end if
   end do  ! this completes the ADF weight factors when a diffraction aperture is present and apinADF=.TRUE.
   
  end if

end if ! if aperture is present

end do outerCLloop   ! see line 814


! and the rest is also only used for debugging purposes
! if (debug) then 
!   allocate(outar(2*nn,self%nml%numberofsvalues,self%nml%numCL))
!   outar(1:nn,1:self%nml%numberofsvalues,1:self%nml%numCL) = self%BFweightsarray(1:nn,1:self%nml%numberofsvalues,1:self%nml%numCL)
!   outar(nn+1:2*nn,1:self%nml%numberofsvalues,1:self%nml%numCL) = self%ADFweightsarray(1:nn,1:self%nml%numberofsvalues,1:self%nml%numCL)
! ! to make sure that everything is correct, let's export this array so that we can display it in IDL
!   open(unit=dataunit,file='STEMprofiles.data',status='unknown',form='unformatted')
!   write(unit=dataunit) nn,self%nml%numberofsvalues,self%nml%numCL
!   write(unit=dataunit) thetar
!   write(unit=dataunit) outar
!   close(unit=dataunit,status='keep')
! end if

end subroutine init_STEM_


!--------------------------------------------------------------------------
!
! SUBROUTINE: init_STEM_ZA
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize weight factors for zone-axis STEM case
! 
!> @note This will need to be reconsidered when we implement sectored detectors ... 
!
!> @param STEM STEM structure
!> @param stemnl STEM namelist
!> @param cell unit cell pointer
!> @param F foil normal
!> @param khead top of kvector list
!> @param reflist top of reflection list
!> @param nn number of reflections
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/12/13 MDG 2.0 rewrite 
!> @date   06/09/14 MDG 3.0 added STEM and cell structures and khead+reflist linked lists
!> @date   06/10/14 MDG 3.1 added F, Dyn argument
!> @date   07/02/17 MDG 3.2 split STEM into STEM and stemnl
!--------------------------------------------------------------------------
recursive subroutine init_STEM_ZA_(self, stemnl, cell, F, Diff, khead, reflist, nn)
!DEC$ ATTRIBUTES DLLEXPORT :: init_STEM_ZA_

use mod_crystallography
use mod_diffraction
use mod_kvectors
use mod_gvectors

IMPLICIT NONE

class(STEM_T),INTENT(INOUT)         :: self
type(STEMGeometryNameListType),INTENT(INOUT)    :: stemnl
type(Cell_T),INTENT(INOUT)          :: cell
real(kind=dbl),INTENT(INOUT)        :: F(3)
type(Diffraction_T),INTENT(INOUT)   :: Diff
type(kvectorlist),pointer           :: khead
type(reflisttype),pointer           :: reflist
integer(kind=irg),INTENT(IN)        :: nn

integer(kind=irg)                   :: ik,ig, iCL
real(kind=sgl)                      :: ll(3), lpg(3), gg(3), glen, gplen, kpg, lambda
type(kvectorlist),pointer           :: ktmp
type(reflisttype),pointer           :: rltmpa

! this routine initializes the excitation error arrays and the weight-factor arrays for zone axis STEM signals
! the weightfactors are quite a bit different from the ones for the systematic row case;
! they are simpler in principle, since each point in the diffracted disk can only lie in one
! place, and hence only contributes to one detector.  However, not all points in a disk
! contribute to the same detector...  The length of the vector k_t+g, expressed in mrad,
! is what needs to be compared to the radii of the BF and ADF detectors.  For each incident 
! beam direction, we take the tangential component of the wave vector and loop over all
! reflections to compute the relevant angle; this then allows us to assign the weight factors
! which are now either 1 or 0 (so they can be stored as logicals).

  lambda = sngl(Diff%getWaveLength())

! allocate the excitation error array areal(1..nn,1..STEM%numk)
  allocate(self%sgarray(nn,self%numk))
  
! transform the foil normal to real space and normalize
  call cell%TransSpace(sngl(F),Diff%Dyn%FN,'d','r')
  call cell%NormVec(Diff%Dyn%FN,'r')

! allocate the weight factor arrays, one entry for each beam direction, reflection, and camera length
  allocate(self%ZABFweightsarray(nn,self%numk,stemnl%numCL))
  allocate(self%ZAADFweightsarray(nn,self%numk,stemnl%numCL))
  self%ZABFweightsarray = .FALSE.
  self%ZAADFweightsarray = .FALSE.

! loop over the wave vector linked list
  ktmp => khead
  beamloopCL: do ik=1,self%numk
    ll = ktmp%kt        ! this is the tangential component of the wave vector
! and loop over all reflections
    rltmpa => reflist%next
    reflectionloopCL: do ig=1,nn
      gg = float(rltmpa%hkl)
      glen = cell%CalcLength(gg,'r')
      lpg = ll + gg                ! Laue + g
      gplen = cell%CalcLength(lpg,'r')
      kpg = 2000.0*asin(0.50*lambda*gplen)    ! 2theta in mrad
      do iCL=1,stemnl%numCL
        self%BFmrad = atan(stemnl%BFradius/stemnl%CLarray(iCL))*1000.0
        self%ADFimrad = atan(stemnl%ADFinnerradius/stemnl%CLarray(iCL))*1000.0
        self%ADFomrad = atan(stemnl%ADFouterradius/stemnl%CLarray(iCL))*1000.0
        if (kpg.le.self%BFmrad) self%ZABFweightsarray(ig,ik,iCL) = .TRUE.
        if ((kpg.ge.self%ADFimrad).AND.(kpg.le.self%ADFomrad)) self%ZAADFweightsarray(ig,ik,iCL) = .TRUE.
      end do  ! loop over camera lengths
      self%sgarray(ig,ik) = Diff%Calcsg(cell,gg,sngl(ktmp%k),Diff%Dyn%FN)
 ! and we move to the next reflection in the list
      rltmpa => rltmpa%next
    end do reflectionloopCL  
    ktmp => ktmp%next
  end do beamloopCL

! that's it folks!
end subroutine init_STEM_ZA_

end module mod_STEM