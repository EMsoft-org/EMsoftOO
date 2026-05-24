! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_EBSDFull
  !! author: MDG/SS + EMsoftOO port
  !! version: 1.0
  !! date: 04/13/26
  !!
  !! Direct/full EBSD pattern simulator.
  !! This is a repo-native port of the corrected legacy EMsoft EMEBSDFull implementation.
  !! It runs a full Monte Carlo calculation for the detector geometry and then
  !! computes the dynamical EBSD patterns directly for each orientation and
  !! energy bin, without going through EMEBSDmaster.

use mod_kinds
use mod_global
use mod_diffraction, only: Diffraction_T

implicit none
private

real(kind=dbl),parameter :: nAmpere = 6.241D18
integer,parameter        :: source_length = 50000

type(Diffraction_T),save :: Diff
!$OMP THREADPRIVATE(Diff)

type :: EBSDFullPixel
  real(kind=sgl),allocatable  :: lambdaEZ(:,:)
  real(kind=dbl)              :: dc(3)
  real(kind=dbl)              :: cfactor
end type EBSDFullPixel

type :: EBSDFullDetectorType
  type(EBSDFullPixel),allocatable :: detector(:,:)
end type EBSDFullDetectorType

type, public :: EBSDFullNameListType
  character(fnlen)        :: xtalname
  real(kind=dbl)          :: dmin
  integer(kind=ill)       :: totnum_el
  integer(kind=irg)       :: multiplier
  real(kind=dbl)          :: EkeV
  real(kind=dbl)          :: Ehistmin
  real(kind=dbl)          :: Ebinsize
  real(kind=dbl)          :: depthmax
  real(kind=dbl)          :: depthstep
  real(kind=dbl)          :: beamcurrent
  real(kind=dbl)          :: dwelltime
  real(kind=dbl)          :: sig
  real(kind=dbl)          :: omega
  real(kind=sgl)          :: L
  real(kind=sgl)          :: thetac
  real(kind=sgl)          :: delta
  integer(kind=irg)       :: numsx
  integer(kind=irg)       :: numsy
  real(kind=sgl)          :: xpc
  real(kind=sgl)          :: ypc
  integer(kind=irg)       :: binning
  character(3)            :: scalingmode
  real(kind=sgl)          :: gammavalue
  character(1)            :: maskpattern
  integer(kind=irg)       :: nthreads
  integer(kind=irg)       :: platid
  integer(kind=irg)       :: devid
  integer(kind=irg)       :: globalworkgrpsz
  integer(kind=irg)       :: num_el
  character(3)            :: eulerconvention
  character(fnlen)        :: anglefile
  character(fnlen)        :: datafile
end type EBSDFullNameListType

type, public :: EBSDFull_T
private
  character(fnlen)           :: nmldeffile = 'EMEBSDFull.nml'
  type(EBSDFullNameListType) :: nml
  type(EBSDFullDetectorType) :: det

contains
private
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: readAngles_
  procedure, pass(self) :: GenerateDetector_
  procedure, pass(self) :: ComputeFullEBSDPatterns_
  procedure, pass(self) :: EBSDFull_

  generic, public :: readNameList => readNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readAngles => readAngles_
  generic, public :: GenerateDetector => GenerateDetector_
  generic, public :: ComputeFullEBSDPatterns => ComputeFullEBSDPatterns_
  generic, public :: EBSDFull => EBSDFull_
end type EBSDFull_T

interface EBSDFull_T
  module procedure EBSDFull_constructor
end interface EBSDFull_T

contains

type(EBSDFull_T) function EBSDFull_constructor(nmlfile) result(EBSDFull)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDFull_constructor

implicit none

character(fnlen),optional :: nmlfile

if (present(nmlfile)) EBSDFull%nmldeffile = trim(nmlfile)
call EBSDFull%readNameList(nmlfile)

end function EBSDFull_constructor

subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_

use mod_io

implicit none

class(EBSDFull_T),intent(inout)      :: self
character(fnlen),intent(in),optional :: nmlfile
logical,optional,intent(in)          :: initonly

type(IO_T)                           :: Message
logical                              :: skipread = .FALSE.
character(fnlen)                     :: localnmlfile

character(fnlen)                     :: xtalname
real(kind=dbl)                       :: dmin
integer(kind=ill)                    :: totnum_el
integer(kind=irg)                    :: multiplier
real(kind=dbl)                       :: EkeV
real(kind=dbl)                       :: Ehistmin
real(kind=dbl)                       :: Ebinsize
real(kind=dbl)                       :: depthmax
real(kind=dbl)                       :: depthstep
real(kind=dbl)                       :: beamcurrent
real(kind=dbl)                       :: dwelltime
real(kind=dbl)                       :: sig
real(kind=dbl)                       :: omega
real(kind=sgl)                       :: L
real(kind=sgl)                       :: thetac
real(kind=sgl)                       :: delta
integer(kind=irg)                    :: numsx
integer(kind=irg)                    :: numsy
real(kind=sgl)                       :: xpc
real(kind=sgl)                       :: ypc
integer(kind=irg)                    :: binning
character(3)                         :: scalingmode
real(kind=sgl)                       :: gammavalue
character(1)                         :: maskpattern
integer(kind=irg)                    :: nthreads
integer(kind=irg)                    :: platid
integer(kind=irg)                    :: devid
integer(kind=irg)                    :: globalworkgrpsz
integer(kind=irg)                    :: num_el
character(3)                         :: eulerconvention
character(fnlen)                     :: anglefile
character(fnlen)                     :: datafile

namelist /EBSDFulldata/ xtalname, dmin, totnum_el, multiplier, EkeV, Ehistmin, Ebinsize, depthmax, depthstep, &
                        beamcurrent, dwelltime, sig, omega, L, thetac, delta, numsx, numsy, xpc, ypc, binning, &
                        scalingmode, gammavalue, maskpattern, nthreads, platid, devid, globalworkgrpsz, num_el, &
                        eulerconvention, anglefile, datafile

xtalname = 'undefined'
dmin = 0.04D0
totnum_el = 2000000000_ill
multiplier = 1
EkeV = 30.D0
Ehistmin = 15.D0
Ebinsize = 1.D0
depthmax = 100.D0
depthstep = 1.D0
beamcurrent = 150.D0
dwelltime = 100.D0
sig = 70.D0
omega = 0.D0
L = 15000.0
thetac = 10.0
delta = 50.0
numsx = 0
numsy = 0
xpc = 0.0
ypc = 0.0
binning = 1
scalingmode = 'lin'
gammavalue = 0.34
maskpattern = 'n'
nthreads = 1
platid = 2
devid = 1
globalworkgrpsz = 512
num_el = 10
eulerconvention = 'tsl'
anglefile = 'undefined'
datafile = 'EBSDout.h5'

if (present(initonly)) then
  if (initonly.eqv..TRUE.) skipread = .TRUE.
end if

localnmlfile = self%nmldeffile
if (present(nmlfile)) localnmlfile = trim(nmlfile)
self%nmldeffile = trim(localnmlfile)

if (.not.skipread) then
  open(unit=dataunit,file=trim(localnmlfile),status='old',delim='apostrophe')
  read(unit=dataunit,nml=EBSDFulldata)
  close(unit=dataunit,status='keep')

  if (trim(xtalname).eq.'undefined') then
    call Message%printError('readNameList',' xtalname is undefined in '//trim(localnmlfile))
  end if
  if (trim(anglefile).eq.'undefined') then
    call Message%printError('readNameList',' anglefile is undefined in '//trim(localnmlfile))
  end if
  if (numsx.le.0 .or. numsy.le.0) then
    call Message%printError('readNameList',' numsx and numsy must be positive in '//trim(localnmlfile))
  end if
end if

self%nml%xtalname = xtalname
self%nml%dmin = dmin
self%nml%totnum_el = totnum_el
self%nml%multiplier = multiplier
self%nml%EkeV = EkeV
self%nml%Ehistmin = Ehistmin
self%nml%Ebinsize = Ebinsize
self%nml%depthmax = depthmax
self%nml%depthstep = depthstep
self%nml%beamcurrent = beamcurrent
self%nml%dwelltime = dwelltime
self%nml%sig = sig
self%nml%omega = omega
self%nml%L = L
self%nml%thetac = thetac
self%nml%delta = delta
self%nml%numsx = numsx
self%nml%numsy = numsy
self%nml%xpc = xpc
self%nml%ypc = ypc
self%nml%binning = binning
self%nml%scalingmode = scalingmode
self%nml%gammavalue = gammavalue
self%nml%maskpattern = maskpattern
self%nml%nthreads = nthreads
self%nml%platid = platid
self%nml%devid = devid
self%nml%globalworkgrpsz = globalworkgrpsz
self%nml%num_el = num_el
self%nml%eulerconvention = eulerconvention
self%nml%anglefile = anglefile
self%nml%datafile = datafile

end subroutine readNameList_

subroutine writeHDFNameList_(self, HDF)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_

use mod_HDFsupport
use ISO_C_BINDING

implicit none

class(EBSDFull_T),intent(inout) :: self
type(HDF_T),intent(inout)       :: HDF

integer(kind=irg),parameter     :: n_int = 7, n_real = 6
integer(kind=irg)               :: hdferr, io_int(n_int)
real(kind=sgl)                  :: io_real(n_real)
character(20)                   :: intlist(n_int), reallist(n_real)
character(fnlen,kind=c_char)    :: line2(1)
character(fnlen)                :: dataset, groupname

associate(enl => self%nml)

groupname = 'EBSDFulldata'
hdferr = HDF%createGroup(groupname)

io_int = (/ enl%multiplier, enl%numsx, enl%numsy, enl%binning, enl%nthreads, enl%platid, enl%devid /)
intlist(1) = 'multiplier'
intlist(2) = 'numsx'
intlist(3) = 'numsy'
intlist(4) = 'binning'
intlist(5) = 'nthreads'
intlist(6) = 'platid'
intlist(7) = 'devid'
call HDF%writeNMLintegers(io_int, intlist, n_int)

io_real = (/ enl%L, enl%thetac, enl%delta, enl%xpc, enl%ypc, enl%gammavalue /)
reallist(1) = 'L'
reallist(2) = 'thetac'
reallist(3) = 'delta'
reallist(4) = 'xpc'
reallist(5) = 'ypc'
reallist(6) = 'gammavalue'
call HDF%writeNMLreals(io_real, reallist, n_real)

dataset = 'dmin'
hdferr = HDF%writeDatasetDouble(dataset, enl%dmin)
dataset = 'totnum_el'
hdferr = HDF%writeDatasetDouble(dataset, dble(enl%totnum_el))
dataset = 'EkeV'
hdferr = HDF%writeDatasetDouble(dataset, enl%EkeV)
dataset = 'Ehistmin'
hdferr = HDF%writeDatasetDouble(dataset, enl%Ehistmin)
dataset = 'Ebinsize'
hdferr = HDF%writeDatasetDouble(dataset, enl%Ebinsize)
dataset = 'depthmax'
hdferr = HDF%writeDatasetDouble(dataset, enl%depthmax)
dataset = 'depthstep'
hdferr = HDF%writeDatasetDouble(dataset, enl%depthstep)
dataset = 'beamcurrent'
hdferr = HDF%writeDatasetDouble(dataset, enl%beamcurrent)
dataset = 'dwelltime'
hdferr = HDF%writeDatasetDouble(dataset, enl%dwelltime)
dataset = 'sig'
hdferr = HDF%writeDatasetDouble(dataset, enl%sig)
dataset = 'omega'
hdferr = HDF%writeDatasetDouble(dataset, enl%omega)
dataset = 'globalworkgrpsz'
hdferr = HDF%writeDatasetInteger(dataset, enl%globalworkgrpsz)
dataset = 'num_el'
hdferr = HDF%writeDatasetInteger(dataset, enl%num_el)

dataset = 'xtalname'
line2(1) = trim(enl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
dataset = 'scalingmode'
line2(1) = trim(enl%scalingmode)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
dataset = 'maskpattern'
line2(1) = trim(enl%maskpattern)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
dataset = 'eulerconvention'
line2(1) = trim(enl%eulerconvention)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
dataset = 'anglefile'
line2(1) = trim(enl%anglefile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
dataset = 'datafile'
line2(1) = trim(enl%datafile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)

call HDF%pop()

end associate

end subroutine writeHDFNameList_

subroutine readAngles_(self, EMsoft, numangles, angles, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: readAngles_

use mod_EMsoft
use mod_io
use mod_rotations
use mod_quaternions

implicit none

class(EBSDFull_T),intent(inout)           :: self
type(EMsoft_T),intent(inout)              :: EMsoft
integer(kind=irg),intent(out)             :: numangles
type(QuaternionArray_T),intent(inout)     :: angles
logical,optional,intent(in)               :: verbose

type(IO_T)                                :: Message
type(e_T)                                 :: e
type(q_T)                                 :: q
type(Quaternion_T)                        :: qq

integer(kind=irg)                         :: io_int(1), i
character(2)                              :: atype
real(kind=dbl)                            :: eulang(3), quatd(4)
character(fnlen)                          :: fname

fname = EMsoft%generateFilePath('EMdatapathname', trim(self%nml%anglefile))
open(unit=dataunit,file=trim(fname),status='old',action='read')

read(dataunit,*) atype
read(dataunit,*) numangles

if (present(verbose)) then
  if (verbose.eqv..TRUE.) then
    io_int(1) = numangles
    call Message%WriteValue(' Number of angle entries = ', io_int, 1)
  end if
end if

call setRotationPrecision('d')
angles = QuaternionArray_T(n=numangles, s='d')

select case (trim(atype))
  case('eu')
    if (present(verbose)) then
      if (verbose.eqv..TRUE.) call Message%printMessage('  -> converting Euler angles to quaternions', frm='(A/)')
    end if
    do i=1,numangles
      read(dataunit,*) eulang
      if (self%nml%eulerconvention.eq.'hkl') eulang(1) = eulang(1) + 90.D0
      call e%e_setd(eulang*dtoR)
      q = e%eq()
      qq = Quaternion_T(qd=q%q_copyd())
      call angles%insertQuatinArray(i, qq)
    end do
  case('qu')
    do i=1,numangles
      read(dataunit,*) quatd
      qq = Quaternion_T(qd=quatd)
      call angles%insertQuatinArray(i, qq)
    end do
  case default
    call Message%printError('readAngles',' unsupported angle file type '//trim(atype))
end select

close(unit=dataunit,status='keep')
call Message%printMessage(' Completed reading orientations')

end subroutine readAngles_

subroutine GenerateDetector_(self, numEbins, numzbins, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: GenerateDetector_

use mod_io
use mod_math

implicit none

class(EBSDFull_T),intent(inout) :: self
integer(kind=irg),intent(in)    :: numEbins, numzbins
logical,optional,intent(in)     :: verbose

  type(IO_T)                      :: Message
  real(kind=sgl),allocatable      :: scin_x(:), scin_y(:)
  real(kind=sgl)                  :: alp, ca, sa, cw, sw
  real(kind=sgl)                  :: L2, Ls, Lc, calpha
  integer(kind=irg)               :: i, j, ipx, ipy, istat
  real(kind=sgl)                  :: dc(3), pcvec(3), alpha, theta, dp

associate(enl => self%nml, det => self%det)

if (.not.allocated(det%detector)) then
  allocate(det%detector(enl%numsx, enl%numsy), stat=istat)
  if (istat.ne.0) call Message%printError('GenerateDetector',' unable to allocate detector pixel array')
end if

allocate(scin_x(enl%numsx), scin_y(enl%numsy), stat=istat)
if (istat.ne.0) call Message%printError('GenerateDetector',' unable to allocate scintillator coordinate arrays')

scin_x = - ( -enl%xpc - (1.0 - enl%numsx) * 0.5 - (/ (real(i-1,kind=sgl), i=1,enl%numsx) /) ) * enl%delta
scin_y =   (  enl%ypc - (1.0 - enl%numsy) * 0.5 - (/ (real(i-1,kind=sgl), i=1,enl%numsy) /) ) * enl%delta

alp = 0.5 * sngl(cPi) - sngl(enl%sig - enl%thetac) * sngl(dtoR)
ca = cos(alp)
sa = sin(alp)
cw = cos(sngl(enl%omega) * sngl(dtoR))
sw = sin(sngl(enl%omega) * sngl(dtoR))

L2 = enl%L * enl%L
do j=1,enl%numsx
  Ls = -sw * scin_x(j) + enl%L * cw
  Lc =  cw * scin_x(j) + enl%L * sw
  do i=1,enl%numsy
    if (.not.allocated(det%detector(j,i)%lambdaEZ)) then
      allocate(det%detector(j,i)%lambdaEZ(1:numEbins,1:numzbins), stat=istat)
      if (istat.ne.0) call Message%printError('GenerateDetector',' unable to allocate lambdaEZ array')
    end if
    det%detector(j,i)%lambdaEZ = 0.0
    dc = (/ (scin_y(i) * ca + sa * Ls), Lc, (-sa * scin_y(i) + ca * Ls) /)
    dc = dc / sngl(vecnorm(dble(dc)))
    det%detector(j,i)%dc = dble(dc)
  end do
end do
deallocate(scin_x, scin_y)

alpha = atan(enl%delta / enl%L / sqrt(sngl(cPi)))
ipx = nint(real(enl%numsx,kind=sgl) * 0.5 + enl%xpc)
ipy = nint(real(enl%numsy,kind=sgl) * 0.5 + enl%ypc)
if ((ipx.lt.1).or.(ipx.gt.enl%numsx).or.(ipy.lt.1).or.(ipy.gt.enl%numsy)) then
  call Message%printError('GenerateDetector',' pattern center lies outside the detector area')
end if

pcvec = real(det%detector(ipx,ipy)%dc, kind=sgl)
calpha = cos(alpha)
do i=1,enl%numsx
  do j=1,enl%numsy
    dc = real(det%detector(i,j)%dc, kind=sgl)
    dp = dot_product(pcvec, dc)
    theta = acos(dp)
    if ((i.eq.ipx).and.(j.eq.ipy)) then
      det%detector(i,j)%cfactor = 0.25D0
    else
      det%detector(i,j)%cfactor = dble(((calpha*calpha + dp*dp - 1.0)**1.5)/(calpha**3))
    end if
  end do
end do

if (present(verbose)) then
  if (verbose.eqv..TRUE.) call Message%printMessage(' --> completed detector generation', frm='(A)')
end if

end associate

end subroutine GenerateDetector_

subroutine EBSDFull_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: EBSDFull_

use mod_EMsoft
use mod_io
use mod_quaternions

implicit none

class(EBSDFull_T),intent(inout) :: self
type(EMsoft_T),intent(inout)    :: EMsoft
character(fnlen),intent(inout)  :: progname

type(IO_T)                      :: Message
type(QuaternionArray_T)         :: angles
integer(kind=irg)               :: numangles, numEbins, numzbins

call self%readAngles(EMsoft, numangles, angles, verbose=.TRUE.)

numEbins = int((self%nml%EkeV - self%nml%Ehistmin)/self%nml%Ebinsize) + 1
numzbins = int(self%nml%depthmax/self%nml%depthstep) + 1
call self%GenerateDetector(numEbins, numzbins, verbose=.TRUE.)
call self%ComputeFullEBSDPatterns(EMsoft, numangles, angles, progname)

call Message%printMessage(' EMEBSDFull completed')

end subroutine EBSDFull_

subroutine ComputeFullEBSDPatterns_(self, EMsoft, numangles, angles, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: ComputeFullEBSDPatterns_

use mod_EMsoft
use mod_symmetry
use mod_crystallography
use mod_io
use mod_initializers
use mod_diffraction, only: DynType, gnode
use mod_gvectors
use mod_quaternions
use mod_rotations
use mod_Lambert
use mod_math
use HDF5
use mod_HDFsupport
use mod_CLsupport
use mod_timing
use clfortran
use omp_lib
use ISO_C_BINDING

implicit none

class(EBSDFull_T),intent(inout)           :: self
type(EMsoft_T),intent(inout)              :: EMsoft
integer(kind=irg),intent(in)              :: numangles
type(QuaternionArray_T),intent(inout)     :: angles
character(fnlen),intent(in)               :: progname

type(SpaceGroup_T)                        :: SG
type(IO_T)                                :: Message
type(Timing_T)                            :: timer
type(Cell_T)                              :: cell
type(DynType)                             :: Dyn
type(HDF_T)                               :: HDF
type(OpenCL_T)                            :: CL
type(gvectors_T)                          :: reflist
type(reflisttype),pointer                 :: firstw
type(Quaternion_T)                        :: quat, quinv
type(q_T)                                 :: qq
type(e_T)                                 :: eu
type(gnode)                               :: rlp
type(Lambert_T)                           :: L

logical                                   :: verbose
integer(kind=irg)                         :: i, j, k, iE, iz, iang, numEbins, numzbins, nseeds, ierrF, hdferr, &
                                             nns, nnw, numset, izz, istat, ipos
integer(kind=irg)                         :: idxy(2), io_int(3)
integer(kind=ill)                         :: num_max, totnum_el, val
integer(kind=8)                           :: size_in_bytes, size_in_bytes_seeds
integer(kind=8),target                    :: globalsize(2)
real(kind=sgl)                            :: io_real(3), xyz(3), r1, r2, r3, rho, edis, alpha, tana, cota, sa, ca
real(kind=sgl),target                     :: EkeV, Ze, density, at_wt, sig, omega
integer(kind=4),target                    :: globalworkgrpsz, num_el, steps
real(kind=sgl),allocatable,target         :: Lamresx(:), Lamresy(:), depthres(:), energyres(:)
integer(kind=4),allocatable,target        :: init_seeds(:)
integer(kind=4),allocatable               :: rnseeds(:)
real(kind=sgl),allocatable                :: eulerangles(:,:)
character(fnlen)                          :: datafile, dataset, dstr, tstrb, tstre, bethefile, fname, sourcefile, &
                                             groupname, datagroupname
character(fnlen,kind=c_char)              :: line2(1)
character(len=source_length),target       :: source
character(len=source_length,kind=c_char),target :: csource
type(c_ptr),target                        :: psource
integer(c_intptr_t),allocatable,target    :: platform(:), device(:)
integer(c_intptr_t),target                :: context, command_queue, prog, kernel, LamX, LamY, depth, energy, seeds
integer(c_int32_t)                        :: ierr, ierr2, pcnt, nump, numd
integer(c_size_t),target                  :: slength
integer(c_size_t)                         :: cnum
character(3),target                       :: kernelname
character(fnlen),target                   :: info
integer                                   :: loopcount, batch

associate(enl => self%nml, det => self%det)

verbose = .TRUE.
call setRotationPrecision('d')

call openFortranHDFInterface()
HDF = HDF_T()

call cell%setFileName(enl%xtalname)
call Diff%setV(enl%EkeV)
call Diff%setrlpmethod('WK')
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, sngl(enl%dmin), verbose=verbose, useHDF=HDF)

call cell%calcDensity()
io_real(1:3) = real(cell%getDensity(), kind=sgl)
density = io_real(1)
at_wt = io_real(2)
Ze = io_real(3)
call Message%WriteValue(' Density, avA, avZ = ', io_real, 3, "(/2f10.5,',',f10.5)")

timer = Timing_T()
tstrb = timer%getTimeString()
dstr = timer%getDateString()
tstre = tstrb

datafile = EMsoft%generateFilePath('EMdatapathname', enl%datafile)
hdferr = HDF%createFile(datafile)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_createFile', hdferr)

datagroupname = 'EBSD'
call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)
call cell%addXtalDataGroup(SG, EMsoft, HDF)

groupname = 'NMLfiles'
hdferr = HDF%createGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_createGroup NMLfiles', hdferr)
dataset = 'EMEBSDFullNML'
hdferr = HDF%writeDatasetTextFile(dataset, self%nmldeffile)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetTextFile', hdferr)
call HDF%pop()

groupname = 'NMLparameters'
hdferr = HDF%createGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_createGroup NMLparameters', hdferr)
call self%writeHDFNameList(HDF)

bethefile = 'BetheParameters.nml'
ipos = scan(trim(datafile), '/\', back=.TRUE.)
if (ipos.gt.0) bethefile = trim(datafile(1:ipos))//'BetheParameters.nml'
call Diff%SetBetheParameters(EMsoft, .FALSE., bethefile)
call Diff%writeBetheparameterNameList(HDF)
call HDF%pop()

groupname = 'EMData'
hdferr = HDF%createGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_createGroup EMData', hdferr)
hdferr = HDF%createGroup(datagroupname)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_createGroup EBSD', hdferr)

dataset = 'xtalname'
line2(1) = trim(enl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetStringArray xtalname', hdferr)

dataset = 'numangles'
hdferr = HDF%writeDatasetInteger(dataset, numangles)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetInteger numangles', hdferr)

allocate(eulerangles(3,numangles))
do i=1,numangles
  quat = angles%getQuatfromArray(i)
  call qq%q_setd(quat%get_quatd())
  eu = qq%qe()
  eulerangles(1:3,i) = real(eu%e_copyd()/dtoR, kind=sgl)
end do
dataset = 'Eulerangles'
hdferr = HDF%writeDatasetFloatArray(dataset, eulerangles, 3, numangles)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetFloatArray Eulerangles', hdferr)
deallocate(eulerangles)

EkeV = sngl(enl%EkeV)
sig = sngl(enl%sig*dtoR)
omega = sngl(enl%omega*dtoR)
globalworkgrpsz = enl%globalworkgrpsz
num_el = enl%num_el
num_max = int(globalworkgrpsz,kind=ill) * int(globalworkgrpsz,kind=ill) * int(num_el,kind=ill)
totnum_el = enl%totnum_el * int(enl%multiplier,kind=ill)
steps = 300
globalsize = (/ int(globalworkgrpsz,kind=8), int(globalworkgrpsz,kind=8) /)
numEbins = int((enl%EkeV - enl%Ehistmin)/enl%Ebinsize) + 1
numzbins = int(enl%depthmax/enl%depthstep) + 1

alpha = 0.5 * sngl(cPi) - sngl(enl%sig - enl%thetac) * sngl(dtoR)
tana = tan(alpha)
cota = 1.0 / tana
sa = sin(alpha)
ca = cos(alpha)

allocate(Lamresx(num_max), Lamresy(num_max), depthres(num_max), energyres(num_max), stat=istat)
if (istat.ne.0) call Message%printError('ComputeFullEBSDPatterns',' unable to allocate Monte Carlo result arrays')
Lamresx = 0.0
Lamresy = 0.0
depthres = 0.0
energyres = 0.0

size_in_bytes = int(num_max,kind=8) * sizeof(EkeV)
size_in_bytes_seeds = int(4*globalworkgrpsz*globalworkgrpsz,kind=8) * sizeof(EkeV)

CL = OpenCL_T()
call CL%init_PDCCQ(platform, nump, enl%platid, device, numd, enl%devid, info, context, command_queue)

sourcefile = 'EMMC.cl'
call Message%printMessage(' OpenCL source file set to : '//trim(sourcefile))
call CL%read_source_file(EMsoft, sourcefile, csource, slength)

pcnt = 1
psource = c_loc(csource)
prog = clCreateProgramWithSource(context, pcnt, c_loc(psource), c_loc(slength), ierr)
call CL%error_check('ComputeFullEBSDPatterns:clCreateProgramWithSource', ierr)

ierr = clBuildProgram(prog, numd, c_loc(device), c_null_ptr, c_null_funptr, c_null_ptr)
ierr2 = clGetProgramBuildInfo(prog, device(enl%devid), CL_PROGRAM_BUILD_LOG, sizeof(source), c_loc(source), cnum)
if (len(trim(source)).gt.0) call Message%printMessage(trim(source(1:cnum)), frm='(A)')
call CL%error_check('ComputeFullEBSDPatterns:clBuildProgram', ierr)
call CL%error_check('ComputeFullEBSDPatterns:clGetProgramBuildInfo', ierr2)

call Message%printMessage(' Program Build Successful... Creating kernel')
kernelname = 'MC'//char(0)
kernel = clCreateKernel(prog, c_loc(kernelname), ierr)
call CL%error_check('ComputeFullEBSDPatterns:clCreateKernel:MC', ierr)
ierr = clReleaseProgram(prog)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseProgram', ierr)

fname = EMsoft%generateFilePath('Randomseedfilename')
open(unit=10,file=trim(fname),form='unformatted',status='old')
read(10) nseeds
allocate(rnseeds(nseeds))
read(10) rnseeds
close(unit=10,status='keep')

if (4*globalworkgrpsz*globalworkgrpsz.gt.nseeds) then
  call Message%printMessage('------------------------------')
  io_int(1) = nseeds
  call Message%WriteValue('Total number of prime number seeds available = ', io_int, 1)
  io_int(1) = 4*globalworkgrpsz*globalworkgrpsz
  call Message%WriteValue('Total number of prime number seeds needed    = ', io_int, 1)
  call Message%printMessage('Please reduce the globalworkgrpsz parameter or increase the number of seeds')
  call Message%printMessage('in the '//trim(fname)//' file.')
  call Message%printError('ComputeFullEBSDPatterns',' insufficient prime number seeds')
end if

allocate(init_seeds(4*globalworkgrpsz*globalworkgrpsz), stat=istat)
if (istat.ne.0) call Message%printError('ComputeFullEBSDPatterns',' unable to allocate OpenCL seed array')
init_seeds = 0
do i=1,globalworkgrpsz
  do j=1,globalworkgrpsz
    do k=1,4
      init_seeds(4*((i-1)*globalworkgrpsz + (j-1)) + k) = rnseeds(4*((i-1)*globalworkgrpsz + j) + k)
    end do
  end do
end do

LamX = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, c_null_ptr, ierr)
call CL%error_check('ComputeFullEBSDPatterns:clCreateBuffer:LamX', ierr)
LamY = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, c_null_ptr, ierr)
call CL%error_check('ComputeFullEBSDPatterns:clCreateBuffer:LamY', ierr)
depth = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, c_null_ptr, ierr)
call CL%error_check('ComputeFullEBSDPatterns:clCreateBuffer:depth', ierr)
energy = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, c_null_ptr, ierr)
call CL%error_check('ComputeFullEBSDPatterns:clCreateBuffer:energy', ierr)
seeds = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, c_null_ptr, ierr)
call CL%error_check('ComputeFullEBSDPatterns:clCreateBuffer:seeds', ierr)

ierr = clEnqueueWriteBuffer(command_queue, seeds, CL_TRUE, 0_8, size_in_bytes_seeds, c_loc(init_seeds(1)), 0, c_null_ptr, c_null_ptr)
call CL%error_check('ComputeFullEBSDPatterns:clEnqueueWriteBuffer', ierr)
call Message%printMessage(' Monte Carlo mode set to full. Performing full calculation...', frm='(A/)')

ierr = clSetKernelArg(kernel, 0, sizeof(LamX), c_loc(LamX))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:LamX', ierr)
ierr = clSetKernelArg(kernel, 1, sizeof(LamY), c_loc(LamY))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:LamY', ierr)
ierr = clSetKernelArg(kernel, 2, sizeof(EkeV), c_loc(EkeV))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:EkeV', ierr)
ierr = clSetKernelArg(kernel, 3, sizeof(globalworkgrpsz), c_loc(globalworkgrpsz))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:globalworkgrpsz', ierr)
ierr = clSetKernelArg(kernel, 4, sizeof(Ze), c_loc(Ze))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:Ze', ierr)
ierr = clSetKernelArg(kernel, 5, sizeof(density), c_loc(density))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:density', ierr)
ierr = clSetKernelArg(kernel, 6, sizeof(at_wt), c_loc(at_wt))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:at_wt', ierr)
ierr = clSetKernelArg(kernel, 7, sizeof(num_el), c_loc(num_el))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:num_el', ierr)
ierr = clSetKernelArg(kernel, 8, sizeof(seeds), c_loc(seeds))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:seeds', ierr)
ierr = clSetKernelArg(kernel, 9, sizeof(sig), c_loc(sig))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:sig', ierr)
ierr = clSetKernelArg(kernel, 10, sizeof(omega), c_loc(omega))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:omega', ierr)
ierr = clSetKernelArg(kernel, 11, sizeof(depth), c_loc(depth))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:depth', ierr)
ierr = clSetKernelArg(kernel, 12, sizeof(energy), c_loc(energy))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:energy', ierr)
ierr = clSetKernelArg(kernel, 13, sizeof(steps), c_loc(steps))
call CL%error_check('ComputeFullEBSDPatterns:clSetKernelArg:steps', ierr)

call timer%Time_tick()
val = 0_ill
loopcount = int(totnum_el/num_max + 1_ill, kind=irg)

do batch=1,loopcount
  ierr = clEnqueueNDRangeKernel(command_queue, kernel, 2, c_null_ptr, c_loc(globalsize), c_null_ptr, 0, c_null_ptr, c_null_ptr)
  call CL%error_check('ComputeFullEBSDPatterns:clEnqueueNDRangeKernel', ierr)
  ierr = clFinish(command_queue)
  call CL%error_check('ComputeFullEBSDPatterns:clFinish', ierr)

  ierr = clEnqueueReadBuffer(command_queue, LamX, CL_TRUE, 0_8, size_in_bytes, c_loc(Lamresx(1)), 0, c_null_ptr, c_null_ptr)
  call CL%error_check('ComputeFullEBSDPatterns:clEnqueueReadBuffer:Lamresx', ierr)
  ierr = clEnqueueReadBuffer(command_queue, LamY, CL_TRUE, 0_8, size_in_bytes, c_loc(Lamresy(1)), 0, c_null_ptr, c_null_ptr)
  call CL%error_check('ComputeFullEBSDPatterns:clEnqueueReadBuffer:Lamresy', ierr)
  ierr = clEnqueueReadBuffer(command_queue, depth, CL_TRUE, 0_8, size_in_bytes, c_loc(depthres(1)), 0, c_null_ptr, c_null_ptr)
  call CL%error_check('ComputeFullEBSDPatterns:clEnqueueReadBuffer:depthres', ierr)
  ierr = clEnqueueReadBuffer(command_queue, energy, CL_TRUE, 0_8, size_in_bytes, c_loc(energyres(1)), 0, c_null_ptr, c_null_ptr)
  call CL%error_check('ComputeFullEBSDPatterns:clEnqueueReadBuffer:energyres', ierr)

  do j=1,int(num_max,kind=irg)
    if ((Lamresx(j).ne.-10.0).and.(Lamresy(j).ne.-10.0).and.(depthres(j).ne.10.0).and.(energyres(j).ne.0.0)) then
      L = Lambert_T(xy=(/ Lamresx(j), Lamresy(j) /))
      ierrF = L%LambertSquareToSphere(xyz)
      if (ierrF.eq.0) then
        xyz = xyz / sngl(vecnorm(dble(xyz)))
        r1 = xyz(1)
        r2 = xyz(2)
        r3 = xyz(3)
        rho = (enl%L/enl%delta) * (tana + cota) / (r3/sa + r1/ca)
        r1 = r1 * rho
        r2 = r2 * rho
        idxy(1) = nint(enl%xpc - r2 + real(enl%numsx,kind=sgl)*0.5)
        idxy(2) = nint(enl%ypc - (r1 - (enl%L/enl%delta)*sa)/ca + real(enl%numsy,kind=sgl)*0.5)
        if ((idxy(1).ge.1).and.(idxy(1).le.enl%numsx).and.(idxy(2).ge.1).and.(idxy(2).le.enl%numsy)) then
          if (energyres(j).gt.enl%Ehistmin) then
            iE = nint((energyres(j) - sngl(enl%Ehistmin))/sngl(enl%Ebinsize)) + 1
            edis = abs(depthres(j))
            iz = nint(edis/sngl(enl%depthstep)) + 1
            if ((iE.ge.1).and.(iE.le.numEbins).and.(iz.ge.1).and.(iz.le.numzbins)) then
              val = val + 1_ill
              det%detector(idxy(1),idxy(2))%lambdaEZ(iE,iz) = det%detector(idxy(1),idxy(2))%lambdaEZ(iE,iz) + 1.0
            end if
          end if
        end if
      end if
    end if
  end do

  if (mod(batch,50).eq.0) then
    io_int(1) = batch * int(num_max,kind=irg)
    call Message%WriteValue('Total number of incident electrons = ', io_int, 1, '(I15)')
    io_int(1) = int(val,kind=irg)
    call Message%WriteValue('Number of BSE electrons intercepted by detector = ', io_int, 1, '(I15)')
  end if
end do

io_int(1) = int(totnum_el,kind=irg)
call Message%WriteValue('Total number of incident electrons = ', io_int, 1, '(I15)')
io_int(1) = int(val,kind=irg)
call Message%WriteValue('Total number of BSE electrons intercepted by detector = ', io_int, 1, '(I15)')
io_real(1) = real(val,kind=sgl) / real(totnum_el,kind=sgl)
call Message%WriteValue('Backscatter yield on detector = ', io_real, 1, '(F15.6)')

ierr = clReleaseKernel(kernel)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseKernel', ierr)
ierr = clReleaseCommandQueue(command_queue)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseCommandQueue', ierr)
ierr = clReleaseContext(context)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseContext', ierr)
ierr = clReleaseMemObject(LamX)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseMemObject:LamX', ierr)
ierr = clReleaseMemObject(LamY)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseMemObject:LamY', ierr)
ierr = clReleaseMemObject(depth)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseMemObject:depth', ierr)
ierr = clReleaseMemObject(energy)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseMemObject:energy', ierr)
ierr = clReleaseMemObject(seeds)
call CL%error_check('ComputeFullEBSDPatterns:clReleaseMemObject:seeds', ierr)

call ComputeFullDynamicalPatterns(self, EMsoft, HDF, cell, SG, timer, datafile, numangles, angles, &
                                  numEbins, numzbins, totnum_el, num_max)

call closeFortranHDFInterface()

end associate

end subroutine ComputeFullEBSDPatterns_

subroutine ComputeFullDynamicalPatterns(self, EMsoft, HDF, cell, SG, timer, datafile, numangles, angles, &
                                        numEbins, numzbins, totnum_el, num_max)

use mod_EMsoft
use mod_symmetry
use mod_crystallography
use mod_io
use mod_gvectors
use mod_diffraction, only: gnode
use mod_quaternions
use mod_rotations
use mod_HDFsupport
use mod_timing
use ISO_C_BINDING
use omp_lib

implicit none

class(EBSDFull_T),intent(inout)       :: self
type(EMsoft_T),intent(inout)          :: EMsoft
type(HDF_T),intent(inout)             :: HDF
type(Cell_T),intent(inout)            :: cell
type(SpaceGroup_T),intent(inout)      :: SG
type(Timing_T),intent(inout)          :: timer
character(fnlen),intent(in)           :: datafile
integer(kind=irg),intent(in)          :: numangles
type(QuaternionArray_T),intent(inout) :: angles
integer(kind=irg),intent(in)          :: numEbins, numzbins
integer(kind=ill),intent(in)          :: totnum_el, num_max

type(IO_T)                            :: Message
type(gvectors_T)                      :: reflist
type(reflisttype),pointer             :: firstw
type(Quaternion_T)                    :: quat, quinv
type(gnode)                           :: rlp

integer(kind=irg)                     :: i, j, k, ix, iE, iang, numset, izz, nns, nnw, istat, hdferr
integer(kind=irg)                     :: nat(maxpasym), io_int(2), badrow, badcol
integer(kind=irg),allocatable         :: numat(:)
real(kind=sgl)                        :: prefactor, fnat, nabsl
real(kind=sgl)                        :: kk(3), FN(3), kkk(3)
real(kind=sgl),allocatable            :: EkeVs(:), nabsfact(:), lambdaZ(:), svals(:)
real(kind=sgl),allocatable            :: EBSDPatterns(:,:,:,:), lambdas(:,:,:,:)
real(kind=dbl),allocatable            :: thick(:,:,:)
complex(kind=dbl),allocatable         :: DynMat(:,:), Sgh(:,:,:), Lghtmp(:,:,:)
character(fnlen)                      :: tstre, dataset, groupname
character(fnlen,kind=c_char)          :: line2(1)
integer,save                          :: nonfinite_trace_count = 0

associate(enl => self%nml, det => self%det)

prefactor = sngl(0.25D0 * nAmpere * enl%beamcurrent * enl%dwelltime * 1.0D-15 / dble(totnum_el))

do i=1,enl%numsx
  do j=1,enl%numsy
    det%detector(i,j)%lambdaEZ = det%detector(i,j)%lambdaEZ / dble(totnum_el + num_max)
  end do
end do

allocate(EkeVs(numEbins), thick(numEbins,enl%numsx,enl%numsy), lambdas(enl%numsx,enl%numsy,numEbins,numzbins), &
         nabsfact(numzbins), stat=istat)
if (istat.ne.0) call Message%printError('ComputeFullDynamicalPatterns',' unable to allocate lambda storage arrays')
thick = 0.D0
lambdas = 0.0
EkeVs = 0.0
nabsfact = 0.0

call Diff%CalcUcg(cell, (/ 0,0,0 /))
rlp = Diff%getrlp()
nabsl = rlp%xgp
do i=1,numzbins
  nabsfact(i) = exp(2.0*sngl(cPi) * real(i-1,kind=sgl) * sngl(enl%depthstep) / nabsl)
end do

do i=1,numEbins
  EkeVs(i) = sngl(enl%Ehistmin + dble(i-1) * enl%Ebinsize)
end do

do i=1,enl%numsx
  do j=1,enl%numsy
    do iE=1,numEbins
      lambdas(i,j,iE,1:numzbins) = det%detector(i,j)%lambdaEZ(iE,1:numzbins) * nabsfact
      izz = 1
      do while ((izz.le.numzbins) .and. (sum(det%detector(i,j)%lambdaEZ(iE,izz:numzbins)).gt.0.0))
        izz = izz + 1
      end do
      thick(iE,i,j) = dble(izz)
    end do
  end do
end do

allocate(EBSDPatterns(enl%numsx,enl%numsy,numEbins,numangles), stat=istat)
if (istat.ne.0) call Message%printError('ComputeFullDynamicalPatterns',' unable to allocate EBSD pattern array')
EBSDPatterns = 0.0

numset = cell%getNatomtype()
numat = cell%getnumat()
nat = 0
do i=1,numset
  nat(i) = numat(i)
end do
if (sum(nat(1:numset)).le.0) call Message%printError('ComputeFullDynamicalPatterns',' invalid nat normalization factor')
fnat = 1.0 / real(sum(nat(1:numset)),kind=sgl)

call Diff%Initialize_SghLUT(cell, SG, sngl(enl%dmin), numset, nat, .FALSE.)

call Message%printMessage(' Starting direct dynamical pattern calculation')
call omp_set_num_threads(enl%nthreads)
io_int(1) = enl%nthreads
call Message%WriteValue(' Attempting to set number of threads to ', io_int, 1, frm='(I4)')

do iang=1,numangles
  quat = angles%getQuatfromArray(iang)
  quinv = conjg(quat)

  do iE=numEbins,1,-1
    call Diff%setV(dble(EkeVs(iE)))
    call Diff%CalcWaveLength(cell)

!$OMP PARALLEL default(SHARED) COPYIN(Diff) &
!$OMP& PRIVATE(i,j,k,ix,kk,kkk,FN,reflist,firstw,nns,nnw,DynMat,Sgh,Lghtmp,lambdaZ,svals,istat,rlp,badrow,badcol) &
!$OMP& PRIVATE(io_int)

    allocate(lambdaZ(numzbins), svals(numset), stat=istat)
    if (istat.ne.0) call Message%printError('ComputeFullDynamicalPatterns',' unable to allocate per-thread arrays')

!$OMP DO SCHEDULE(DYNAMIC)
    do k=1,enl%numsx*enl%numsy
      j = mod(k-1,enl%numsy) + 1
      i = (k-1)/enl%numsy + 1

      lambdaZ = det%detector(i,j)%lambdaEZ(iE,1:numzbins) * nabsfact
      kk = real(quinv%quat_Lp(det%detector(i,j)%dc), kind=sgl)
      kk = kk / sngl(Diff%getWaveLength())

      call cell%TransSpace(kk, kkk, 'c', 'r')
      FN = kkk

      reflist = gvectors_T()
      call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, kkk, sngl(enl%dmin), .FALSE.)
      nullify(firstw)
      nns = 0
      nnw = 0
      call reflist%Apply_BethePotentials(Diff, firstw, nns, nnw)

      allocate(DynMat(nns,nns), stat=istat)
      if (istat.ne.0) call Message%printError('ComputeFullDynamicalPatterns',' unable to allocate dynamical matrix')
      DynMat = cmplx(0.D0,0.D0,dbl)
      call reflist%GetDynMat(cell, Diff, firstw, DynMat, nns, nnw)

      allocate(Sgh(nns,nns,numset), Lghtmp(nns,nns,numzbins), stat=istat)
      if (istat.ne.0) call Message%printError('ComputeFullDynamicalPatterns',' unable to allocate scattering matrices')
      Sgh = cmplx(0.D0,0.D0,dbl)
      Lghtmp = cmplx(0.D0,0.D0,dbl)

      call reflist%getSghfromLUT(Diff, nns, numset, Sgh)
      call CalcLghSM(DynMat, nns, lambdaZ, numzbins, sngl(enl%depthstep), Diff%getWaveLength(), Lghtmp, &
                     badrow, badcol, nonfinite_trace_count)

      svals = 0.0
      do ix=1,numset
        svals(ix) = real(sum(Lghtmp(1:nns,1:nns,numzbins) * Sgh(1:nns,1:nns,ix)))
      end do
      svals = svals * fnat
      EBSDPatterns(i,j,iE,iang) = sum(svals) * prefactor

      call reflist%Delete_gvectorlist()
      deallocate(Sgh, Lghtmp, DynMat)

      if (mod(k,10000).eq.0) then
        io_int(1) = k
        io_int(2) = enl%numsx*enl%numsy
        call Message%WriteValue(' completed ', io_int, 2, "(I8,' of ',I8)")
      end if
    end do
!$OMP END DO

    deallocate(lambdaZ, svals)
!$OMP END PARALLEL

    io_int(1) = iE
    io_int(2) = iang
    call Message%WriteValue('  completed energy bin for angle ', io_int, 2, "(2I8)")
  end do
end do

call timer%Time_tock()
dataset = 'EBSDPatterns'
hdferr = HDF%writeDatasetFloatArray(dataset, EBSDPatterns, enl%numsx, enl%numsy, numEbins, numangles)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetFloatArray4D EBSDPatterns', hdferr)
dataset = 'Lambdas'
hdferr = HDF%writeDatasetFloatArray(dataset, lambdas, enl%numsx, enl%numsy, numEbins, numzbins)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetFloatArray4D Lambdas', hdferr)
dataset = 'EkeVs'
hdferr = HDF%writeDatasetFloatArray(dataset, EkeVs, numEbins)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetFloatArray EkeVs', hdferr)

call HDF%pop()
call HDF%pop()

call timer%makeTimeStamp()
tstre = timer%getTimeString()
groupname = 'EMheader'
hdferr = HDF%openGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_openGroup EMheader', hdferr)
groupname = 'EBSD'
hdferr = HDF%openGroup(groupname)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_openGroup EBSD header', hdferr)
dataset = 'StopTime'
line2(1) = trim(timer%getDateString())//', '//trim(tstre)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite=.TRUE.)
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetStringArray StopTime', hdferr)
dataset = 'Duration'
hdferr = HDF%writeDatasetFloat(dataset, timer%getInterval())
if (hdferr.ne.0) call HDF%error_check('EMEBSDFull:HDF_writeDatasetFloat Duration', hdferr)
call HDF%popall()

end associate

end subroutine ComputeFullDynamicalPatterns

recursive subroutine CalcLghSM(DynMat, nn, lambdaZ, nt, dthick, wavelength, Lgh, badrow, badcol, tracecount)

use mod_math
use, intrinsic :: ieee_arithmetic

implicit none

complex(kind=dbl),intent(in)           :: DynMat(nn,nn)
integer(kind=irg),intent(in)           :: nn
real(kind=sgl),intent(in)              :: lambdaZ(1:nt)
integer(kind=irg),intent(in)           :: nt
real(kind=sgl),intent(in)              :: dthick
real(kind=dbl),intent(in)              :: wavelength
complex(kind=dbl),intent(out)          :: Lgh(nn,nn,nt)
integer(kind=irg),intent(out)          :: badrow, badcol
integer,intent(inout)                  :: tracecount

integer                                :: i, ii, jj
complex(kind=dbl),allocatable          :: Minp(:,:), Azz(:,:), ampl(:), ampl2(:)
logical                                :: badMinp, badAzz

allocate(Minp(nn,nn), Azz(nn,nn), ampl(nn), ampl2(nn))

Minp = DynMat * cmplx(0.D0, cPi * wavelength, dbl)
badMinp = .FALSE.
badrow = 0
badcol = 0
do jj=1,nn
  do ii=1,nn
    if ((.not.ieee_is_finite(real(Minp(ii,jj)))).or.(.not.ieee_is_finite(aimag(Minp(ii,jj))))) then
      badMinp = .TRUE.
      badrow = ii
      badcol = jj
      exit
    end if
  end do
  if (badMinp.eqv..TRUE.) exit
end do
if ((badMinp.eqv..TRUE.).and.(tracecount.lt.10)) then
  tracecount = tracecount + 1
  write (*,'(A)') 'EMEBSDFull trace: non-finite Minp before MatrixExponential'
  write (*,'(A,I8)') '  nn        = ', nn
  write (*,'(A,I8)') '  row       = ', badrow
  write (*,'(A,I8)') '  col       = ', badcol
end if

call MatrixExponential(Minp, Azz, dble(dthick), 'Pade', nn)

badAzz = .FALSE.
badrow = 0
badcol = 0
do jj=1,nn
  do ii=1,nn
    if ((.not.ieee_is_finite(real(Azz(ii,jj)))).or.(.not.ieee_is_finite(aimag(Azz(ii,jj))))) then
      badAzz = .TRUE.
      badrow = ii
      badcol = jj
      exit
    end if
  end do
  if (badAzz.eqv..TRUE.) exit
end do
if ((badAzz.eqv..TRUE.).and.(tracecount.lt.10)) then
  tracecount = tracecount + 1
  write (*,'(A)') 'EMEBSDFull trace: non-finite Azz after MatrixExponential'
  write (*,'(A,I8)') '  nn        = ', nn
  write (*,'(A,I8)') '  row       = ', badrow
  write (*,'(A,I8)') '  col       = ', badcol
end if

ampl = cmplx(0.D0,0.D0,dbl)
ampl(1) = cmplx(1.D0,0.D0,dbl)
Lgh = cmplx(0.D0,0.D0,dbl)

do i=1,nt
  ampl2 = matmul(Azz, ampl)
  if (i.eq.1) then
    Lgh(1:nn,1:nn,i) = lambdaZ(i) * spread(ampl2(1:nn),dim=2,ncopies=nn) * spread(conjg(ampl2(1:nn)),dim=1,ncopies=nn)
  else
    Lgh(1:nn,1:nn,i) = Lgh(1:nn,1:nn,i-1) + lambdaZ(i) * spread(ampl2(1:nn),dim=2,ncopies=nn) * &
                       spread(conjg(ampl2(1:nn)),dim=1,ncopies=nn)
  end if
  ampl = ampl2
end do

deallocate(Minp, Azz, ampl, ampl2)

end subroutine CalcLghSM

end module mod_EBSDFull
