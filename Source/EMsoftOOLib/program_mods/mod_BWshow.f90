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

module mod_BWshow
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/23/24
  !!
  !! class definition for the EMBWshow program

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: BWshow_T
  character(fnlen)                :: outname 
  character(fnlen)                :: plotprefix
  character(2)                    :: TBSR 
  integer(kind=irg)               :: nn, ns, numthick 
  integer(kind=irg),allocatable   :: g(:), k(:), fn(:), PIVOT(:,:)
  real(kind=sgl)                  :: kz, startthick, thickinc
  real(kind=sgl),allocatable      :: kttb(:), kn(:)
  real(kind=dbl),allocatable      :: PVarrays(:,:,:)
  complex(kind=dbl),allocatable   :: W(:,:), CG(:,:,:), alpha(:,:)

contains
private 
  procedure, pass(self) :: BWshow_
  procedure, pass(self) :: readTBSRBWfile_
  procedure, pass(self) :: pivot2_
  procedure, pass(self) :: pivot3_
  procedure, pass(self) :: makePIVOTs_
  procedure, pass(self) :: BWtoI_ 

  generic, public :: BWshow => BWshow_

end type BWshow_T

! the constructor routine for this class 
interface BWshow_T
  module procedure BWshow_constructor
end interface BWshow_T

contains

!--------------------------------------------------------------------------
type(BWshow_T) function BWshow_constructor( ) result(BWshow)
!! author: MDG 
!! version: 1.0 
!! date: 02/23/24
!!
!! constructor for the BWshow_T Class; reads the name list 
 
IMPLICIT NONE


end function BWshow_constructor

!--------------------------------------------------------------------------
subroutine BWshow_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 02/23/24
!!
!! destructor for the BWshow_T Class
 
IMPLICIT NONE

type(BWshow_T), INTENT(INOUT)  :: self 

call reportDestructor('BWshow_T')

end subroutine BWshow_destructor

!--------------------------------------------------------------------------
subroutine makePIVOTs_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: makePIVOTs_
!! author: MDG 
!! version: 1.0 
!! date: 02/26/24
!!
!! make a LAPACK pivot table to an array of arrays

IMPLICIT NONE

class(BWshow_T),INTENT(INOUT)     :: self

real(kind=dbl)                    :: M(self%nn,self%nn), ident(self%nn,self%nn), row(self%nn)
integer(kind=irg)                 :: i, j 

allocate(self%PVarrays(self%nn,self%nn,self%ns))

ident = 0.D0
do i=1,self%nn
  ident(i,i) = 1.D0
end do 

do i=1,self%ns 
  M = ident 
  do j=1,self%nn 
    if (j.ne.self%PIVOT(j,i)) then ! we need to swap two rows
      row(:) = M(self%PIVOT(j,i),:)
      M(self%PIVOT(j,i),:) = M(j,:)
      M(j,:) = row(:)
    end if 
  end do
  self%PVarrays(:,:,i) = transpose(M)
end do 

end subroutine makePIVOTs_

!--------------------------------------------------------------------------
subroutine pivot2_(self, A )
!DEC$ ATTRIBUTES DLLEXPORT :: pivot2_
!! author: MDG 
!! version: 1.0 
!! date: 02/26/24
!!
!! apply a LAPACK pivot table to an array of column vectors

IMPLICIT NONE

class(BWshow_T),INTENT(INOUT)     :: self
real(kind=dbl),INTENT(INOUT)      :: A(self%nn,self%ns)

integer(kind=irg)                 :: i,j
real(kind=dbl)                    :: tmp(self%nn), row, M(self%nn,self%nn)

do i=1,self%ns 
  tmp(:) = A(:,i)
  M = self%PVarrays(:,:,i)
  tmp = matmul(M, tmp)
  ! do j=1,self%nn
  !   if (j.ne.self%PIVOT(j,i)) then ! we need to swap two rows
  !     row = tmp(self%PIVOT(j,i))
  !     tmp(self%PIVOT(j,i)) = tmp(j)
  !     tmp(j) = row
  !   end if 
  ! end do
  A(:,i) = tmp(:)
end do

end subroutine pivot2_

!--------------------------------------------------------------------------
subroutine pivot3_(self, A )
!DEC$ ATTRIBUTES DLLEXPORT :: pivot3_
!! author: MDG 
!! version: 1.0 
!! date: 02/26/24
!!
!! apply a LAPACK pivot table to an array of arrays

IMPLICIT NONE

class(BWshow_T),INTENT(INOUT)     :: self
real(kind=dbl),INTENT(INOUT)      :: A(self%nn,self%nn,self%ns)

integer(kind=irg)                 :: i,j
real(kind=dbl)                    :: tmp(self%nn,self%nn), row(self%nn), M(self%nn,self%nn)

do i=1,self%ns 
  tmp(:,:) = A(:,:,i)
  M = self%PVarrays(:,:,i)
  tmp = matmul(M, tmp)
  ! do j=1,self%nn
  !   if (j.ne.self%PIVOT(j,i)) then ! we need to swap two rows
  !     row(:) = tmp(self%PIVOT(j,i),:)
  !     tmp(self%PIVOT(j,i),:) = tmp(j,:)
  !     tmp(j,:) = row(:)
  !   end if 
  ! end do
  A(:,:,i) = tmp(:,:)
end do

end subroutine pivot3_

!--------------------------------------------------------------------------
subroutine readTBSRBWfile_(self, EMsoft, HDFnames, noPIVOT)
!DEC$ ATTRIBUTES DLLEXPORT :: readTBSRBWfile_
!! author: MDG 
!! version: 1.0 
!! date: 02/23/24
!!
!! read data from HDF5 Bloch wave file 

use mod_EMsoft
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_io
use stringconstants

use ISO_C_BINDING

IMPLICIT NONE

class(BWshow_T),INTENT(INOUT)     :: self
type(EMsoft_T),INTENT(INOUT)      :: EMsoft
type(HDFnames_T),INTENT(INOUT)    :: HDFnames
logical,INTENT(IN),OPTIONAL       :: noPIVOT

type(HDF_T)                       :: HDF
type(IO_T)                        :: Message

character(fnlen)                  :: fname, groupname, dataset 
integer(kind=irg)                 :: hdferr 
integer(HSIZE_T)                  :: dims1(1), dims2(2), dims3(3)
real(kind=dbl),allocatable        :: W_R(:,:), CG_R(:,:,:), alpha_R(:,:), W_I(:,:), CG_I(:,:,:), alpha_I(:,:)
logical                           :: readonly=.TRUE., applyPIVOT = .TRUE. 

if (present(noPIVOT)) then 
  if (noPIVOT.eqv..TRUE.) applyPIVOT=.FALSE. 
end if 

call openFortranHDFInterface()
HDF = HDF_T() 

fname = EMsoft%generateFilePath('EMdatapathname', self%outname)
hdferr =  HDF%openFile(fname, readonly)

groupname = SC_EMData
  hdferr = HDF%openGroup(groupname)
groupname = trim(HDFnames%get_ProgramData())
  hdferr = HDF%openGroup(groupname)

dataset = 'nn'
  call HDF%readDatasetInteger(dataset, hdferr, self%nn)

if (self%nn.eq.2) then 
  self%TBSR = 'TB'
else
  self%TBSR = 'SR'
end if

dataset = 'ns'
  call HDF%readDatasetInteger(dataset, hdferr, self%ns)

dataset = 'g'
  call HDF%readDatasetIntegerArray(dataset, dims1, hdferr, self%g)

dataset = 'k'
  call HDF%readDatasetIntegerArray(dataset, dims1, hdferr, self%k)

dataset = 'f'
  call HDF%readDatasetIntegerArray(dataset, dims1, hdferr, self%fn)

dataset = 'kttb'
  call HDF%readDatasetFloatArray(dataset, dims2, hdferr, self%kttb)

dataset = 'kn'
  call HDF%readDatasetFloatArray(dataset, dims2, hdferr, self%kn)

dataset = 'kz'
  call HDF%readDatasetFloat(dataset, hdferr, self%kz)

! read the pivot array; the LAPACK routines use row pivoting to compute the 
! matrix inverse, so we need to undo the pivots by multiplying by the 
! transpose of the permutation matrix.
if (applyPIVOT.eqv..TRUE.) then 
dataset = 'PIVOT'
  if (allocated(self%PIVOT)) deallocate(self%PIVOT)
  call HDF%readDatasetIntegerArray(dataset, dims2, hdferr, self%PIVOT)
  call Message%printMessage(' applying pivot table to parameter arrays...')
  call self%makePIVOTs_()
end if 

dataset = 'W_R'
  if (allocated(self%W)) deallocate(self%W)
  call HDF%readDatasetDoubleArray(dataset, dims2, hdferr, W_R)
  if (applyPIVOT.eqv..TRUE.) call self%pivot2_(W_R)

dataset = 'W_I'
  call HDF%readDatasetDoubleArray(dataset, dims2, hdferr, W_I)
  if (applyPIVOT.eqv..TRUE.) call self%pivot2_(W_I)

allocate(self%W(dims2(1),dims2(2)))
self%W = cmplx(W_R,W_I,dbl)
deallocate(W_R,W_I)

dataset = 'CG_R'
  if (allocated(self%CG)) deallocate(self%CG)
  call HDF%readDatasetDoubleArray(dataset, dims3, hdferr, CG_R)
  if (applyPIVOT.eqv..TRUE.) call self%pivot3_(CG_R)

dataset = 'CG_I'
  call HDF%readDatasetDoubleArray(dataset, dims3, hdferr, CG_I)
  if (applyPIVOT.eqv..TRUE.) call self%pivot3_(CG_I)

allocate(self%CG(dims3(1),dims3(2),dims3(3)))
self%CG = cmplx(CG_R,CG_I,dbl)
deallocate(CG_R,CG_I)

dataset = 'alpha_R'
  if (allocated(self%alpha)) deallocate(self%alpha)
  call HDF%readDatasetDoubleArray(dataset, dims2, hdferr, alpha_R)
  if (applyPIVOT.eqv..TRUE.) call self%pivot2_(alpha_R)

dataset = 'alpha_I'
  call HDF%readDatasetDoubleArray(dataset, dims2, hdferr, alpha_I)
  if (applyPIVOT.eqv..TRUE.) call self%pivot2_(alpha_I)

allocate(self%alpha(dims2(1),dims2(2)))
self%alpha = cmplx(alpha_R,alpha_I,dbl)
deallocate(alpha_R,alpha_I)

call HDF%pop() 
call HDF%pop() 

! get the thickness parameters
  hdferr = HDF%openGroup(HDFnames%get_NMLparameters())
  hdferr = HDF%openGroup(HDFnames%get_NMLlist())

dataset = SC_numthick
  call HDF%readDatasetInteger(dataset, hdferr, self%numthick)

dataset = SC_startthick
  call HDF%readDatasetFloat(dataset, hdferr, self%startthick)

dataset = SC_thickinc
  call HDF%readDatasetFloat(dataset, hdferr, self%thickinc)

call HDF%popall()

call closeFortranHDFInterface()

 if (applyPIVOT.eqv..TRUE.) deallocate(self%PVarrays)

end subroutine readTBSRBWfile_

!--------------------------------------------------------------------------
subroutine BWtoI_(self, EMsoft, minthick, maxthick)
!DEC$ ATTRIBUTES DLLEXPORT :: BWtoI_
!! author: MDG 
!! version: 1.0 
!! date: 02/23/24
!!
!! read data from HDF5 Bloch wave file 

use mod_EMsoft
use mod_image
use mod_io

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

class(BWshow_T),INTENT(INOUT)     :: self
type(EMsoft_T),INTENT(INOUT)      :: EMsoft
real(kind=sgl),INTENT(IN)         :: minthick
real(kind=sgl),INTENT(IN)         :: maxthick

type(IO_T)                        :: Message 

integer(kind=irg)                 :: i, j, k     
real(kind=sgl),allocatable        :: images(:,:,:)
real(kind=dbl)                    :: dz, z, arg, Wr(self%nn), Wi(self%nn)
real(kind=sgl)                    :: mi, ma, ima(self%ns, self%ns) 
complex(kind=dbl)                 :: amp, q(self%nn), diag(self%nn) 
character(fnlen)                  :: TIFF_filename

! declare variables for use in object oriented image module
integer                           :: iostat
character(len=128)                :: iomsg
logical                           :: isInteger
type(image_t)                     :: im
integer(int8)                     :: i8 (3,4)
integer(int8), allocatable        :: TIFF_image(:,:)

! all images are square
allocate(images(self%ns,self%ns,self%nn))

dz = ( maxthick - minthick )/float(self%ns)

do k=1,self%ns        ! horizontal image dimension
  do i=1,self%ns      ! vertical image dimension
    z = minthick + dz*float(i)
    arg = 2.0*sngl(cPi)*z
    Wr(:) = arg * self%W(:,k)%re
    Wi(:) = arg * self%W(:,k)%im
    q(:) = cmplx(cos(Wr(:)),sin(Wr(:)))
    diag(:) = exp(Wi(:)) * q(:) * self%alpha(:,k)
    do j=1,self%nn    ! individual diffracted beams
     amp = sum(self%CG(j,:,k)*diag(:))
     images(k,i,j) = sngl(abs(amp)**2)
    end do 
  end do
end do

allocate(TIFF_image(self%nn*self%ns,self%ns))

! scale each image of the systematic row separately
do k=1,self%nn 
  ima = images(:,:,k)
  mi = minval(ima)
  ma = maxval(ima)
  TIFF_image((k-1)*self%ns:k*self%ns-1, self%ns:1:-1) = nint(255*(ima(1:self%ns,1:self%ns)-mi)/(ma-mi))
end do 

! delineate the edges between the images by means of gray 2-pixel wide lines
do k=1,self%nn-1
  TIFF_image(k*self%ns:k*self%ns+1, 1:self%ns) = nint(126.5)
end do

TIFF_filename = EMsoft%generateFilePath('EMdatapathname',trim(self%plotprefix)//'.tiff') 

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMgetADP","failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else
  call Message%printMessage(' Systematic row images written to '//trim(TIFF_filename))
end if

end subroutine BWtoI_

!--------------------------------------------------------------------------
subroutine BWshow_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: BWshow_
!! author: MDG 
!! version: 1.0 
!! date: 02/23/24
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_io
use mod_axis
use mod_postscript

IMPLICIT NONE 

class(BWshow_T), INTENT(INOUT)      :: self
type(EMsoft_T), INTENT(INOUT)       :: EMsoft
character(fnlen), INTENT(INOUT)     :: progname 
type(HDFnames_T), INTENT(INOUT)     :: HDFnames

type(IO_T)                          :: Message
type(axis_T)                        :: axis
type(PostScript_T)                  :: PS 
type(axis_T)                        :: AX

integer(kind=irg)                   :: io_int(2) 
real(kind=sgl)                      :: io_real(2)
character(12),parameter             :: yt(4) = ['gamma^(j)   ', &
                                                'k^(j)_z-k_0 ', &
                                                'alpha^(j)   ', &
                                                'q^(j)       ']
character(40),parameter             :: gt(4) = ['Bloch wave eigenvalues                  ', &
                                                'Bloch wave eigenvalues                  ', &
                                                'Bloch wave excitation amplitudes        ', &
                                                'Bloch wave absorption parameters        ']
real,parameter                      :: xo(4)=[0.5,4.25,0.5,4.25],yo(4)=[5.5,5.5,1.0,1.0]
real,parameter                      :: xx(11) = [0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35,4.45,0.25,2.35], &
                                       yy(11) = [6.6,6.6,6.6,4.5,4.5,4.5,2.4,2.4,2.4,0.3,0.3]
integer(kind=irg)                   :: n(1),nn, nf, i, j 
real(kind=sgl),allocatable          :: images(:,:,:),y(:),worky(:,:)
real(kind=sgl)                      :: xmin,xmax,ymin,ymax,imax,imin,tmin,tmax
integer(kind=irg)                   :: ntx,nty,tob
logical                             :: npg, readonly=.TRUE.
character(fnlen)                    :: oname

! get the input filename
call Message%printMessage(' BWshow: generating plots ')

! read all the Bloch wave data from the HDF5 file 
call self%readTBSRBWfile_(EMsoft, HDFnames)

if (self%TBSR.eq.'TB') then 
  call Message%printMessage(' Read two-beam data from HDF5 file')
  io_int(1) = self%ns
  call Message%WriteValue(' Number of beam directions  : ', io_int, 1)
else 
  call Message%printMessage(' Read systematic row data from HDF5 file')
  io_int(1) = self%nn
  call Message%WriteValue(' Number of diffracted beams : ', io_int, 1)
  io_int(1) = self%ns
  call Message%WriteValue(' Number of beam directions  : ', io_int, 1)
end if 

io_real = (/ self%startthick, self%startthick + (self%numthick-1)*self%thickinc /)
call Message%WriteValue(' Minimum and maximum foil thickness : ', io_real, 2)
io_int(1) = self%numthick
call Message%WriteValue(' Number of thickness values         : ', io_int, 1)

! we'll generate PostScript output for the four plots, and then tiff output for
! the images... All four plots will fit on one page

! original workx is equal to self%kttb array
allocate(worky(self%ns,self%nn), y(self%ns))
xmin = minval(self%kttb)
xmax = maxval(self%kttb)

oname = trim(self%plotprefix)//'.eps'
write (*,*) ' oname = ', trim(oname)
oname = EMsoft%generateFilePath('EMdatapathname', oname)
write (*,*) ' oname = ', trim(oname)

! the first time we call the axis class explicitly
! AX = axis_T( PS, 4.2, xo(1), yo(1) )
AX = axis_T( 4.2, xo(1), yo(1) )
AX%PS = PostScript_T(progname, EMsoft, 0, dontask=.TRUE., psname=oname)

! call AX%PS%newpage(.FALSE.,'Bloch Wave Results')

do nf = 1, 4 

  select case (nf)
    case(1) ! plot 1 Bloch wave eigenvalues (gamma^(j))
      worky = -self%W%re

    case(2) ! plot 2 Bloch wave eigenvalues (k^(j)_z-k_0)
      do i=1,self%ns
        worky(:,i) = -(self%kn(i)+self%W(:,i)%re+self%kz)
      end do 

    case(3) ! plot 3 Bloch wave excitation amplitudes (alpha^(j))
      do i=1,self%ns
        worky(:,i) = abs(self%alpha(:,i))**2
      end do

    case(4) ! plot 4 Bloch wave absorption parameters (q^(j))
      worky = self%W%im

    case default
  end select

  ymin = minval(worky)
  ymax = maxval(worky)

! create the axis class and pass the PostScript class into it
  if (nf.ne.1) then 
    AX%axw = 4.2
    AX%xll = xo(nf)
    AX%yll = yo(nf)
  end if 
  y(:) = worky(1,:)
  call AX%axis(self%ns,self%kttb,y,xmin,xmax,ymin,ymax,.FALSE.,.FALSE.,'lin','lin','CON',1, &
               'BOT','LEF',.FALSE.,.TRUE.,gt(nf),'kt/g',yt(nf))
! superimpose more curves
  do j=2,self%nn
    y(1:self%ns) = worky(j,1:self%ns)
    call AX%axis(self%ns,self%kttb,y,xmin,xmax,ymin,ymax,.FALSE.,.FALSE.,'lin','lin','CON',1, &
                'BOT','LEF',.TRUE.,.FALSE.,' ',' ',' ')
  end do

end do 

call Ax%PS%closefile()

deallocate(worky)

! compute the output image for the systematic row
! to do so, we need to reread the Bloch wave arrays from the h5 file, but this time
! without application of the PIVOT array...

! read all the Bloch wave data from the HDF5 file 
call self%readTBSRBWfile_(EMsoft, HDFnames, noPIVOT=.TRUE.)

call self%BWtoI_(EMsoft, self%startthick, self%startthick+self%numthick*self%thickinc )

! and, finally, delete the HDF5 file since we really no longer need it
! and it is very fast to redo the computations if necessary


end subroutine BWshow_



end module mod_BWshow