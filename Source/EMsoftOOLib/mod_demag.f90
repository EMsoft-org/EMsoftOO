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

module mod_demag
  !! author: MDG 
  !! version: 1.0 
  !! date: 08/16/21
  !!
  !! class definition for the EMdemag program
  !!
  !! this class allows for the computation of the following quantities:
  !! - demag tensor field                                   -> demag%getDemagTensorField
  !! - volume-averaged demag tensor                         -> demag%getDemagFactors
  !! - demag energy field                                   -> demag%getEnergyField
  !! - magnetic induction configuration                     -> demag%getMagneticInduction
  !! - magnetic phase shift for a given shape amplitude     -> demag%getMagneticPhaseShift
  !! - generate an STL file to represent the demag tensor field with ellipsoids -> demag%dtf2stl
  !!

use mod_kinds
use mod_global
! use mod_ShapeAmplitude

IMPLICIT NONE 

! namelist for the EMdemag program
type, public :: demagNameListType
  integer(kind=irg)         :: dims         ! all 3 dimensions are the same to keep things simple
  logical                   :: dtf
  character(fnlen)          :: shampFilename
  character(fnlen)          :: outname
end type demagNameListType

! class definition
type, public :: demag_T
private 
  character(fnlen)                  :: nmldeffile = 'EMdemag.nml'
  type(demagNameListType)           :: nml 
  real(kind=dbl),allocatable        :: kx(:,:,:), ky(:,:,:), kz(:,:,:), kk(:,:,:), gfilter(:,:,:)  
  real(kind=dbl),allocatable        :: Nij(:,:,:,:)   ! use Voigt notation for 6 tensor components
  real(kind=dbl)                    :: Nav(6)
  complex(kind=dbl),allocatable     :: shamp(:,:,:)
  real(kind=dbl)                    :: dk, sc 
  ! type(ShapeAmplitudeNameListType)  :: shampnml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: Demag_
  procedure, pass(self) :: getDemagTensorField_ 
  procedure, pass(self) :: getDemagFactors_
  procedure, pass(self) :: get_dims_
  procedure, pass(self) :: get_dtf_
  procedure, pass(self) :: get_shampFilename_
  procedure, pass(self) :: get_outname_
  procedure, pass(self) :: set_dims_
  procedure, pass(self) :: set_dtf_
  procedure, pass(self) :: set_shampFilename_
  procedure, pass(self) :: set_outname_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: Demag => Demag_
  generic, public :: getDemagFactors => getDemagFactors_
  generic, public :: getDemagTensorField => getDemagTensorField_
  generic, public :: get_dims => get_dims_
  generic, public :: get_dtf => get_dtf_
  generic, public :: get_shampFilename => get_shampFilename_
  generic, public :: get_outname => get_outname_
  generic, public :: set_dims => set_dims_
  generic, public :: set_dtf => set_dtf_
  generic, public :: set_shampFilename => set_shampFilename_
  generic, public :: set_outname => set_outname_

end type demag_T

! the constructor routine for this class 
interface demag_T
  module procedure demag_constructor
end interface demag_T

contains

!--------------------------------------------------------------------------
type(demag_T) function demag_constructor( nmlfile ) result(demag)
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! constructor for the demag_T Class; reads the name list 
 
use mod_io 
use mod_ShapeAmplitude

IMPLICIT NONE

character(fnlen), OPTIONAL    :: nmlfile 

type(IO_T)                    :: Message 
type(ShapeAmplitude_T)        :: shamp 

integer(kind=irg)             :: d, i, j, shd(3), io_int(2) 
real(kind=dbl),allocatable    :: line(:) 

call demag%readNameList(nmlfile)

call Message%printMessage(' Generating k-vector coordinate arrays... ',frm="(A$)")

! create the reciprocal space k-vector arrays;
! to facilitate fft computations, all arrays will have the k-space origin in position (1,1,1)
d = demag%nml%dims
demag%dk = 1.D0/dble(d)
demag%sc = 2.D0*cPi*demag%dk

allocate(demag%kx(d,d,d), demag%ky(d,d,d), demag%kz(d,d,d), demag%kk(d,d,d), demag%gfilter(d,d,d))
allocate(line(d))

! generate a single line
line = cshift( demag%sc * ( (/ (dble(i),i=0,d-1) /) - dble(d/2) ), d/2, 1) 

! and populate the coordinate arrays
do i=1,d 
  do j=1,d 
    demag%kx(:,i,j) = line(:)
    demag%ky(i,:,j) = line(:)
    demag%kz(i,j,:) = line(:)
  end do 
end do 

! compute the distance array kk 
demag%kk = sqrt( demag%kx*demag%kx+demag%ky*demag%ky+demag%kz*demag%kz )
! and set the origin to 1 to prevent division by zero 
demag%kk(1,1,1) = 1.D0
! scale the component arrays
demag%kx = demag%kx/demag%kk
demag%ky = demag%ky/demag%kk
demag%kz = demag%kz/demag%kk

demag%gfilter = (sin(demag%kk)/demag%kk)**2 - sin(2.D0*demag%kk)/(2.0*demag%kk)
demag%gfilter = 3.D0 * demag%gfilter/demag%kk**2/dble(d)**3
demag%gfilter(1,1,1) = 0.D0

call Message%printMessage(' done.')

call Message%printMessage(' Reading shape amplitude from file... ', advance='no')

shamp = ShapeAmplitude_T()
call shamp%readShapeAmplitude(demag%nml%shampFilename, demag%shamp) ! demag%shampnml

! check to make sure the dimensions of the shamp array are compatible with the coordinate arrays
shd = shape(demag%shamp) 
if (shd(1).ne.d) then 
  io_int = (/ shd(1), d /)
  call Message%printMessage('')
  call Message%WriteValue(' Array dimension mismatch: ', io_int, 2, frm="(I5,'<---->',I5)")
  call Message%printError(' EMdemag class constructor', 'Incompatible array dimension in shape amplitude file')
end if 

call Message%printMessage(' done.')

end function demag_constructor

!--------------------------------------------------------------------------
subroutine demag_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! destructor for the demag_T Class
 
IMPLICIT NONE

type(demag_T), INTENT(INOUT)  :: self 

call reportDestructor('demag_T')

end subroutine demag_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! read the namelist from an nml file for the demag_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)        :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)         :: dims
logical                   :: dtf
character(fnlen)          :: shampFilename
character(fnlen)          :: outname

namelist / Demag / dims, dtf, shampFilename, outname 

dims = 128 
dtf = .TRUE. 
shampFilename = 'undefined'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=Demag)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(shampFilename).eq.'undefined') then
        call Message%printError('readNameList:',' shape amplitude file name is undefined in '//nmlfile)
    end if

    if (trim(outname).eq.'undefined') then
        call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
    end if
end if 

self%nml%dims = dims 
self%nml%dtf = dtf 
self%nml%shampFilename = shampFilename
self%nml%outname = outname

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! pass the namelist for the demag_T Class to the calling program

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)          :: self
type(demagNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(demag_T), INTENT(INOUT)           :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 2, n_real = 9
integer(kind=irg)                       :: hdferr,  io_int(n_int), d
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( emnl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

d = 0
if (emnl%dtf.eqv..TRUE.) d = 1
io_int = (/ emnl%dims, d /)
intlist(1) = 'dims'
intlist(2) = 'dtf'
call HDF%writeNMLintegers(io_int, intlist, n_int)

dataset = 'shampFilename'
line2(1) = emnl%shampFilename
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create shampFilename dataset', hdferr)

dataset = 'outname'
line2(1) = emnl%outname
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create outname dataset', hdferr)

! pop this group off the stack
call HDF%pop()
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
function get_dims_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dims_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! get dims from the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
integer(kind=irg)                 :: out

out = self%nml%dims

end function get_dims_

!--------------------------------------------------------------------------
subroutine set_dims_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dims_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! set dims in the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)     :: inp

self%nml%dims = inp

end subroutine set_dims_

!--------------------------------------------------------------------------
function get_dtf_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dtf_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! get dtf from the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
logical                           :: out

out = self%nml%dtf

end function get_dtf_

!--------------------------------------------------------------------------
subroutine set_dtf_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_dtf_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! set dtf in the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
logical, INTENT(IN)               :: inp

self%nml%dtf = inp

end subroutine set_dtf_

!--------------------------------------------------------------------------
function get_shampFilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_shampFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! get shampFilename from the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
character(fnlen)                  :: out

out = self%nml%shampFilename

end function get_shampFilename_

!--------------------------------------------------------------------------
subroutine set_shampFilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_shampFilename_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! set shampFilename in the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)      :: inp

self%nml%shampFilename = inp

end subroutine set_shampFilename_

!--------------------------------------------------------------------------
function get_outname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_outname_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! get outname from the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
character(fnlen)                  :: out

out = self%nml%outname

end function get_outname_

!--------------------------------------------------------------------------
subroutine set_outname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_outname_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/2021
!!
!! set outname in the demag_T class

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)      :: inp

self%nml%outname = inp

end subroutine set_outname_

!--------------------------------------------------------------------------
subroutine getDemagTensorField_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getDemagTensorField_
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! compute the 3D demag tensor field; uses Voigt notation due to tensor symmetry 

use mod_fftw3
use mod_EMsoft

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)             :: self 

type(EMsoft_T)                            :: EMsoft
integer(kind=irg)                         :: d
type(C_PTR)                               :: plan
character(fnlen)                          :: p
complex(C_DOUBLE_COMPLEX),pointer         :: inp(:,:,:)
complex(C_DOUBLE_COMPLEX),pointer         :: outp(:,:,:)
real(kind=dbl),allocatable                :: dtfcomp(:,:,:)
type(C_PTR)                               :: i, o

p = ''
EMsoft = EMsoft_T(p, p, silent=.TRUE.)
d = self%nml%dims 

! allocate fftw arrays
i = fftw_alloc_complex(int(d*d*d,C_SIZE_T))
call c_f_pointer(i, inp, [d,d,d])

o = fftw_alloc_complex(int(d*d*d,C_SIZE_T))
call c_f_pointer(o, outp, [d,d,d])

! allocate the output array 
allocate( self%Nij(d,d,d,6) )

inp = cmplx(0.D0,0D0)
outp = cmplx(0.D0,0.D0)

! set up the fftw plan for the forward transform
plan = fftw_plan_dft_3d(d,d,d,inp,outp,FFTW_BACKWARD,FFTW_ESTIMATE)

! compute the inverse ffts of the shape amplitude multiplied with a combination
! of k-vector coordinate arrays and then shift each to put the origin at the center
! of the array
allocate(dtfcomp(d,d,d))

! Nxx
inp = self%shamp * cmplx(self%kx**2 * self%gfilter) 
call fftw_execute_dft(plan, inp, outp)
dtfcomp = cshift(real(outp),d/2,1)
dtfcomp = cshift(dtfcomp,d/2,2)
dtfcomp = cshift(dtfcomp,d/2,3)
self%Nij(:,:,:,1) = dtfcomp  

! Nyy
inp = self%shamp * cmplx(self%ky**2 * self%gfilter) 
call fftw_execute_dft(plan, inp, outp)
dtfcomp = cshift(real(outp),d/2,1)
dtfcomp = cshift(dtfcomp,d/2,2)
dtfcomp = cshift(dtfcomp,d/2,3)
self%Nij(:,:,:,2) = dtfcomp  

! Nzz
inp = self%shamp * cmplx(self%kz**2 * self%gfilter) 
call fftw_execute_dft(plan, inp, outp)
dtfcomp = cshift(real(outp),d/2,1)
dtfcomp = cshift(dtfcomp,d/2,2)
dtfcomp = cshift(dtfcomp,d/2,3)
self%Nij(:,:,:,3) = dtfcomp  

! Nyz, Nzy
inp = self%shamp * cmplx(self%kz*self%ky * self%gfilter) 
call fftw_execute_dft(plan, inp, outp)
dtfcomp = cshift(real(outp),d/2,1)
dtfcomp = cshift(dtfcomp,d/2,2)
dtfcomp = cshift(dtfcomp,d/2,3)
self%Nij(:,:,:,4) = dtfcomp  

! Nxz, Nzx
inp = self%shamp * cmplx(self%kx*self%kz * self%gfilter) 
call fftw_execute_dft(plan, inp, outp)
dtfcomp = cshift(real(outp),d/2,1)
dtfcomp = cshift(dtfcomp,d/2,2)
dtfcomp = cshift(dtfcomp,d/2,3)
self%Nij(:,:,:,5) = dtfcomp  

! Nxy, Nyx
inp = self%shamp * cmplx(self%kx*self%ky * self%gfilter) 
call fftw_execute_dft(plan, inp, outp)
dtfcomp = cshift(real(outp),d/2,1)
dtfcomp = cshift(dtfcomp,d/2,2)
dtfcomp = cshift(dtfcomp,d/2,3)
self%Nij(:,:,:,6) = dtfcomp  

call fftw_free(i)
call fftw_free(o)
deallocate(dtfcomp)

end subroutine getDemagTensorField_

!--------------------------------------------------------------------------
subroutine getDemagFactors_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getDemagFactors_
!! author: MDG 
!! version: 1.0 
!! date: 08/17/21
!!
!! compute the volume averaged demag factors; uses Voigt notation due to tensor symmetry 

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)             :: self 

integer(kind=irg)                         :: d
real(kind=dbl),allocatable                :: shint(:,:,:)
real(kind=dbl)                            :: s

d = self%nml%dims 

! shape intensity
allocate(shint(d,d,d))
shint = abs(self%shamp)**2

! Nxx
self%Nav(1) = sum(shint * self%kx**2) 

! Nyy
self%Nav(2) = sum(shint * self%ky**2)

! Nzz
self%Nav(3) = sum(shint * self%kz**2)

! Nyz, Nzy
self%Nav(4) = sum(shint * self%kz*self%ky)

! Nxz, Nzx
self%Nav(5) = sum(shint * self%kz*self%kx)

! Nxy, Nyx
self%Nav(6) = sum(shint * self%kx*self%ky)

! divide by volume and d^3
self%Nav = self%Nav/self%shamp(1,1,1)/dble(d)**3

! Normalize to force the trace equal to one ()
s = sum(self%Nav(1:3))
self%Nav = self%Nav/s 

deallocate(shint)

end subroutine getDemagFactors_

!--------------------------------------------------------------------------
subroutine Demag_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: Demag_
!! author: MDG 
!! version: 1.0 
!! date: 08/16/21
!!
!! perform the computations and store the output in an HDF5 file 

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_timing
use stringconstants

IMPLICIT NONE 

class(demag_T), INTENT(INOUT)       :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 
type(HDFnames_T), INTENT(INOUT)     :: HDFnames

type(Timing_T)                      :: timer 
type(HDF_T)                         :: HDF 
character(fnlen,kind=c_char)        :: line2(1)
character(fnlen)                    :: datafile, dataset, datagroupname, attributename, &
                                       HDF_FileVersion 
character(80)                       :: header
character(11)                       :: dstr
character(15)                       :: tstrb
character(15)                       :: tstre
integer(kind=irg)                   :: hdferr, d 

associate(emnl => self%nml)

d = emnl%dims 

! do we need the demagnetization tensor field ?
if (emnl%dtf.eqv..TRUE.) then 
  call self%getDemagTensorField_()
  write (*,*) ' central trace of Nij = ', sum(self%Nij(d/2,d/2,d/2,1:3))
end if 

call self%getDemagFactors_()
write (*,*) ' volume-averaged demag factors : ', self%Nav

!====================================
!====================================
!====================================
! generate an HDF5 file with all the necessary arrays ... 
call openFortranHDFInterface()
HDF = HDF_T()

timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()

! Create a new file using the default properties.
datafile = EMsoft%generateFilePath('EMdatapathname', emnl%outname)

hdferr =  HDF%createFile(datafile)
if (hdferr.ne.0) call HDF%error_check('HDF_createFile ', hdferr)

dataset = SC_Manufacturer
line2(1) = 'EMsoftOO'
line2(1) = cstringify(line2(1))
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)

! write the EMheader to the file
datagroupname = trim(HDFnames%get_ProgramData()) 
call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(HDFnames%get_NMLfiles())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLfiles', hdferr)

call HDF%pop()

! create a NMLparameters group to write all the namelist entries into
hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup NMLparameters', hdferr)

! to be written ... 

! and leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(HDFnames%get_EMData())
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup EMData', hdferr)

! create the sub group and add a HDF_FileVersion attribute to it
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('HDF_createGroup ShapeAmplitude', hdferr)
HDF_FileVersion = '4.1'
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

dataset = 'demagtensorfield'
hdferr = HDF%writeDatasetDoubleArray(dataset, self%Nij, d, d, d, 6)

dataset = 'kx'
hdferr = HDF%writeDatasetDoubleArray(dataset, self%kx, d, d, d)

call HDF%pop(.TRUE.)
call closeFortranHDFInterface()

end associate 

end subroutine Demag_



end module mod_demag