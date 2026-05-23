! ###################################################################
! Copyright (c) 2013-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_Lorentz
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/01/24
  !!
  !! class definition for the EMLorentz program
  !!
  !! Many of these routines were developed with DOE Basic Energy Sciences 
  !! support during the years 2003-2012; most of them were originally written in 
  !! IDL (Interactive Data Language) and have been translated into f90.

use mod_kinds
use mod_global
use mod_fftw3

IMPLICIT NONE 

! namelist for the EMLorentz program
type, public :: LorentzNameListType
end type LorentzNameListType

! class definition
type, public :: Lorentz_T
private 
  character(fnlen)                            :: nmldeffile = 'EMLorentz.nml'
  type(LorentzNameListType)                   :: nml 
! fftw variables
  type(C_PTR)                                 :: planf, planb
  complex(C_DOUBLE_COMPLEX),allocatable       :: inp(:,:), outp(:,:)

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: getPhaseMapMansuripur_
  procedure, pass(self) :: getFFTWplans_1D_
  procedure, pass(self) :: getFFTWplans_2D_
  procedure, pass(self) :: getFFTWplans_3D_
  procedure, pass(self) :: Lorentz_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: getPhaseMapMansuripur => getPhaseMapMansuripur_
  generic, public :: getFFTWplans => getFFTWplans_1D_, getFFTWplans_2D_, getFFTWplans_3D_
  generic, public :: Lorentz => Lorentz_

end type Lorentz_T

! the constructor routine for this class 
interface Lorentz_T
  module procedure Lorentz_constructor
end interface Lorentz_T

contains

!--------------------------------------------------------------------------
type(Lorentz_T) function Lorentz_constructor( nmlfile ) result(Lorentz)
!! author: MDG 
!! version: 1.0 
!! date: 03/01/24
!!
!! constructor for the Lorentz_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call Lorentz%readNameList(nmlfile)

end function Lorentz_constructor

!--------------------------------------------------------------------------
subroutine Lorentz_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 03/01/24
!!
!! destructor for the Lorentz_T Class
 
IMPLICIT NONE

type(Lorentz_T), INTENT(INOUT)  :: self 

call reportDestructor('Lorentz_T')

end subroutine Lorentz_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 03/01/24
!!
!! read the namelist from an nml file for the Lorentz_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(Lorentz_T), INTENT(INOUT) :: self
character(fnlen),INTENT(IN)     :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)     :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                  :: EMsoft 
type(IO_T)                      :: Message       
logical                         :: skipread = .FALSE.



end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 03/01/24
!!
!! pass the namelist for the Lorentz_T Class to the calling program

IMPLICIT NONE 

class(Lorentz_T), INTENT(INOUT)          :: self
type(LorentzNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 03/01/24
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(Lorentz_T), INTENT(INOUT) :: self 
type(HDF_T), INTENT(INOUT)      :: HDF
type(HDFnames_T), INTENT(INOUT) :: HDFnames

integer(kind=irg),parameter     :: n_int = 11, n_real = 9
integer(kind=irg)               :: hdferr,  io_int(n_int)
real(kind=sgl)                  :: io_real(n_real)
character(20)                   :: intlist(n_int), reallist(n_real)
character(fnlen)                :: dataset, sval(1),groupname
character(fnlen,kind=c_char)    :: line2(1)

associate( mcnl => self%nml )

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
recursive subroutine getPhaseMapMansuripur_(self, EMsoft, B0t, beam, nx, ny, Mag, phasemap, init, destroy)
!DEC$ ATTRIBUTES DLLEXPORT :: getPhaseMapMansuripur_

use mod_io
use mod_EMsoft

use, intrinsic                        :: iso_c_binding

class(Lorentz_T),INTENT(INOUT)        :: self
type(EMsoft_T),INTENT(INOUT)          :: EMsoft
real(kind=dbl),INTENT(IN)             :: B0t
real(kind=dbl),INTENT(IN)             :: beam(3)
integer(kind=irg),INTENT(IN)          :: nx
integer(kind=irg),INTENT(IN)          :: ny
complex(kind=dbl),INTENT(IN)          :: Mag(3,nx,ny)
real(kind=dbl),INTENT(OUT)            :: phasemap(nx,ny)
logical,INTENT(IN),OPTIONAL           :: init
logical,INTENT(IN),OPTIONAL           :: destroy

integer(kind=irg)                     :: ix, iy, iloc
real(kind=dbl)                        :: prefac, b(3), sx, sy, fnx, fny, s, sigx, sigy, gp, pre, psig
complex(C_DOUBLE_COMPLEX),allocatable :: bx(:,:), by(:,:), bz(:,:)
complex(kind=dbl)                     :: cone = cmplx(1.D0,0.D0), czero = cmplx(0.D0,0.D0), d, prex, prey
logical                               :: igp


! are we just destroying the fftw plans ?
if (present(destroy)) then
  if (destroy.eqv..TRUE.) then
    deallocate(self%inp, self%outp)
    call fftw_destroy_plan(self%planf)
    call fftw_destroy_plan(self%planb)
    return
  end if
end if

! if init=.TRUE. then initialize the fftw plans
if (present(init)) then
  if (init) then
    allocate(self%inp(nx, ny), self%outp(nx, ny)) 
! initialize arrays
  self%inp = czero
  self%outp = czero

! then we set up the fftw plans for forward and reverse transforms
  call self%getFFTWplans_2D_( (/ nx, ny /) )
   return
  end if
end if

! prefactor e/hbar (factor of pi cancels out in the end)
prefac = 2.D0*cCharge*1.0D-18/cPlanck
fnx = 1.D0/dble(nx)
fny = 1.D0/dble(ny)

! normalize the incident beam direction
b = beam / sqrt(sum(beam**2))

! if the beam is along the z-axis, we do not need to compute the
! function G_p(ts) (controlled by variable igp)
igp=.TRUE.
if ((beam(1).eq.0.D0).AND.(beam(2).eq.0.D0)) then
  igp=.FALSE.
endif

! if neither init nor destroy are present then we do an actual computation of the phase shift
allocate(bx(ny,ny),by(nx,ny),bz(nx,ny)) 
bx = czero
by = czero
bz = czero

! transform the magnetization to Fourier space (componentwise)
inp(1:nx,1:ny) = Mag(1,1:nx,1:ny)
call fftw_execute_dft(self%planf, self%inp, self%outp)
bx = outp

inp(1:nx,1:ny) = Mag(2,1:nx,1:ny)
call fftw_execute_dft(self%planf, self%inp, self%outp)
by = outp

inp(1:nx,1:ny) = Mag(3,1:nx,1:ny)
call fftw_execute_dft(self%planf, self%inp, self%outp)
bz = outp

! Compute eqn. 13a-b in the Mansuripur paper.
d = czero
! loop over the y-axis
do iy=1,ny
  if (iy.lt.ny/2) then
    sy = dble(iy-1)*fny
  else
    sy = dble(iy-1)*fny - 1.D0
  end if
! loop over the x-axis
  do ix=1,nx
    if (ix.lt.nx/2) then
      sx = dble(ix-1)*fnx
    else
      sx = dble(ix-1)*fnx - 1.D0
    end if
! compute normalized frequency components
    s=sqrt(sx**2+sy**2)
    if (s.eq.0.D0) then
     sigx=0.D0
     sigy=0.D0
     iloc = 1
    else 
     sigx=sx/s
     sigy=sy/s
     iloc = 0
    endif
! compute the products of various vectors
    prex = by(ix,iy)*cmplx(beam(1)**2+beam(3)**2,0.D0) &
          -bx(ix,iy)*cmplx(beam(1)*beam(2),0.D0) &
          -bz(ix,iy)*cmplx(beam(2)*beam(3),0.D0)
    prey = -bx(ix,iy)*cmplx(beam(2)**2+beam(3)**2,0.D0) &
          +by(ix,iy)*cmplx(beam(1)*beam(2),0.D0) &
          +bz(ix,iy)*cmplx(beam(1)*beam(3),0.D0)
    d = prex * cmplx(sigx,0.D0) + prey * cmplx(sigy,0.D0)
! compute G_p(ts) (or not)
    gp=1.D0
    if (igp) then 
      psig=(beam(1)*sigx+beam(2)*sigy)/beam(3)
      pre=1.D0/(psig**2+1.D0)/beam(3)**2
      if ((psig.ne.0.D0).and.(iloc.ne.1)) then
        psig=psig*cPi*s*B0t
        gp=pre*sin(psig)/psig
      else
        gp=pre
      end if
    end if
! multiply the whole thing and store in bx
    if (iloc.ne.1) then 
      bx(ix,iy) = cmplx(0.D0,gp*B0t/s) * d
    else
      bx(ix,iy) = czero
    endif

  end do
end do

! perform the inverse transform
inp = bx
call fftw_execute_dft(self%planb, self%inp, self%outp)
phasemap = real(self%outp) * prefac / float(nx) / float(ny)


end subroutine getPhaseMapMansuripur_

!--------------------------------------------------------------------------
recursive subroutine getFFTWplans_1D_(self, dims)
!DEC$ ATTRIBUTES DLLEXPORT :: getFFTWplans_1D_

use, intrinsic               :: iso_c_binding
use mod_FFTW3

integer(kind=irg),INTENT(IN) :: dims(1)

complex(C_DOUBLE_COMPLEX)    :: inp(dims(1)), outp(dims(1))

inp = cmplx(0.D0,0.D0)
outp = cmplx(0.D0,0.D0)

! then we set up the fftw plans for forward and reverse transforms
self%planf = fftw_plan_dft_1d(dims(1), inp, outp, FFTW_FORWARD, FFTW_ESTIMATE)
self%planb = fftw_plan_dft_1d(dims(1), inp, outp, FFTW_BACKWARD, FFTW_ESTIMATE)

end subroutine getFFTWplans_1D_

!--------------------------------------------------------------------------
recursive subroutine getFFTWplans_2D_(self, dims)
!DEC$ ATTRIBUTES DLLEXPORT :: getFFTWplans_2D_

use, intrinsic               :: iso_c_binding
use mod_FFTW3

integer(kind=irg),INTENT(IN) :: dims(2)

complex(C_DOUBLE_COMPLEX)    :: inp(dims(1),dims(2)), outp(dims(1),dims(2))

inp = cmplx(0.D0,0.D0)
outp = cmplx(0.D0,0.D0)

! then we set up the fftw plans for forward and reverse transforms
self%planf = fftw_plan_dft_2d(dims(2), dims(1), inp, outp, FFTW_FORWARD, FFTW_ESTIMATE)
self%planb = fftw_plan_dft_2d(dims(2), dims(1), inp, outp, FFTW_BACKWARD, FFTW_ESTIMATE)

end subroutine getFFTWplans_2D_

!--------------------------------------------------------------------------
recursive subroutine getFFTWplans_3D_(self, dims)
!DEC$ ATTRIBUTES DLLEXPORT :: getFFTWplans_3D_

use, intrinsic               :: iso_c_binding
use mod_FFTW3

integer(kind=irg),INTENT(IN) :: dims(3)

complex(C_DOUBLE_COMPLEX)    :: inp(dims(1),dims(2),dims(3)), outp(dims(1),dims(2),dims(3))

inp = cmplx(0.D0,0.D0)
outp = cmplx(0.D0,0.D0)

! then we set up the fftw plans for forward and reverse transforms
self%planf = fftw_plan_dft_3d(dims(3), dims(2), dims(1), inp, outp, FFTW_FORWARD, FFTW_ESTIMATE)
self%planb = fftw_plan_dft_3d(dims(3), dims(2), dims(1), inp, outp, FFTW_BACKWARD, FFTW_ESTIMATE)

end subroutine getFFTWplans_3D_




!--------------------------------------------------------------------------
subroutine Lorentz_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: Lorentz_
!! author: MDG 
!! version: 1.0 
!! date: 03/01/24
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames

IMPLICIT NONE 

class(Lorentz_T), INTENT(INOUT) :: self
type(EMsoft_T), INTENT(INOUT)   :: EMsoft
character(fnlen), INTENT(INOUT) :: progname 
type(HDFnames_T), INTENT(INOUT) :: HDFnames

end subroutine Lorentz_



end module mod_Lorentz