! ###################################################################
! Copyright (c) 2016-2024, Marc De Graef Research Group/Carnegie Mellon University
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

program EMTBBFDF
  !! author: MDG
  !! version: 1.0 
  !! date: 02/21/24
  !!
  !! Two beam bright field dark field image pairs (from CTEM book)
  !!
  !! This is a legacy program to reproduce some of the figures in chapter 6 of the CTEM book
  !! The original code has been updated using the modern object oriented fortran 2018 approach. 

use mod_kinds
use mod_global
use mod_EMsoft
use mod_crystallography
use mod_symmetry
use mod_diffraction
use mod_BFDF
use mod_image
use ISO_C_BINDING
use mod_memory
use mod_HDFsupport
use mod_TB
use mod_misc

use, intrinsic :: iso_fortran_env

IMPLICIT NONE

character(fnlen)                    :: progname = 'EMTBBFDF.f90'
character(fnlen)                    :: progdesc = 'Two-beam bright field-dark field images [from CTEM book]'

type(EMsoft_T)                      :: EMsoft
type(memory_T)                      :: mem 
type(IO_T)                          :: Message
type(Cell_T)                        :: cell 
type(SpaceGroup_T)                  :: SG
type(Diffraction_T)                 :: Diff
type(gnode)                         :: rlp

real(kind=sgl),parameter            :: thr = 1.0E-6
real(kind=sgl),allocatable          :: BF(:,:), DF(:,:), exe(:,:), z(:,:)
real(kind=sgl)                      :: dz, zmin, zmax
integer(kind=irg)                   :: jdim, zero(3), ind(3), io_int(3), iflag, j,i, ihole, jhole, irad, jrad, nw, &
                                       npx, npy, k, TIFF_nx, TIFF_ny
real(kind=sgl)                      :: r(200),p(200),xig,xigp,betag,xizero, io_real(3), ss, sm, f1, dd, T0, T1, &
                                       av, zdev, ff, fi, fj, ztot, q, zmi, zma, tps, x0, y0, scl
integer(kind=irg)                   :: values(1:8), kk
integer, dimension(:), allocatable  :: seed
real(kind=dbl)                      :: zz
character(fnlen)                    :: TIFF_filename, xtalname, fname
logical                             :: isallowed
real(kind=dbl)                      :: eps = 1.0D-6

! declare variables for use in object oriented image module
integer                             :: iostat
character(len=128)                  :: iomsg
logical                             :: isInteger
type(image_t)                       :: im
integer(int8)                       :: i8 (3,4)
integer(int8), allocatable          :: TIFF_image(:,:)

! print the EMsoft header, no command line argument handling 
EMsoft = EMsoft_T( progname, progdesc, noCLA=.TRUE. )

! memory class 
mem = memory_T()

call openFortranHDFInterface()

! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft)
call cell%calcPositions(SG, 'v')

call Diff%GetVoltage(cell)
call Diff%setrlpmethod('WK')

call Message%printMessage('Enter Miller indices :', frm = "(/A)")
call GetIndex(SG%getSpaceGrouphexset(), ind, 'r')
call Diff%CalcUcg(cell, ind)
isallowed = SG%IsGAllowed(ind)
rlp = Diff%getrlp() 
if ((rlp%Umod.lt.eps).and.(isallowed.eqv..FALSE.)) then
  call Message%printError('EMTBBFDF',' This reflection is absent due to lattice centering; please try another one.')
end if
if ((rlp%Umod.lt.eps).and.(isallowed.eqv..TRUE.)) then
  call Message%printError('EMTBBFDF','This reflection is absent due to glide/screw symmetry; please try another one.')
end if
xigp = rlp%xgp
xig  = rlp%xg
betag = rlp%Vpphase-rlp%Vphase
io_real(1) = xig
call Message%WriteValue(' Extinction distance [nm] = ', io_real, 1, "(f10.4)")
io_real(1) = xigp
call Message%WriteValue(' Absorption length   [nm] = ', io_real, 1, "(f10.4)")

! initialize all arrays
call Message%ReadValue(' Enter dimension of output image (square)', io_int, 1)
jdim = io_int(1)
call mem%alloc(BF, (/ jdim, jdim /), 'DF', initval = 0.0)
call mem%alloc(DF, (/ jdim, jdim /), 'BF', initval = 0.0)
call mem%alloc(exe, (/ jdim, jdim /), 'exe', initval = 0.0)
call mem%alloc(z, (/ jdim, jdim /), 'z', initval = 0.0)

! normal absorption length
call Diff%CalcUcg(cell, (/ 0, 0, 0 /) )
rlp = Diff%getrlp()
xizero = rlp%xgp
io_real(1) = xizero
call Message%WriteValue(' Normal Absorption   [nm] = ', io_real, 1, "(f10.4)")

! the following used to be in a separate file (BFDF.routines); in the current version,
! the user is prompted for a selection.
call Message%printMessage( (/ ' Which foil geometry do you wish to use ? ', &
                              ' [1] Straight wedge                       ', &
                              ' [2] Bent wedge                           ', &
                              ' [3] Bent wedge with a hole               ', &
                              ' [4] Random superposition of cosine waves ' /) )
call Message%ReadValue('   -> your choice : ', io_int, 1)

select case (io_int(1))
  case (1) ! Case 1 straight wedge
    call BFDF_straight_wedge(jdim, z, exe, xig)

  case (2) ! Case 2 bent wedge
    call BFDF_bent_wedge(jdim, z, exe, xig)

  case(3) ! Case 3 bent foil with circular hole
    call BFDF_bent_hole(jdim, z, exe, xig)

  case(4) ! Case 4 Random superposition of cosine waves
    call BFDF_random_thickness(jdim, z, exe, xig)

  case default
    call Message%printError('EMTBBFDF',' Option not implemented')
end select  

! compute intensities
 T0= SECNDS(0.0)
 call Message%printMessage('computation start ')
 iflag=0
 do i=1,jdim
  do j=1,jdim
   call TBCalcInten(BF(i,j),DF(i,j),exe(i,j),z(i,j),xig,xigp,xizero,betag)
   if (iflag.eq.0) then 
    if ((BF(i,j).gt.1.D0).or.(DF(i,j).gt.1.D0)) iflag=1
   endif
  end do
 end do
 T1= SECNDS(T0)

 io_real(1)=T1
 call Message%WriteValue(' computation end   -> ', io_real, 1, "(F,' seconds')")
 tps = T1/float(jdim)**2
 io_real(1) = tps 
 call Message%WriteValue(' time per pixel    -> ', io_real, 1, "(F,' seconds')")
 if (iflag.eq.1) then 
  call Message%printWarning('EMTBBFDF',(/' WARNING: one or more intensities larger than 1 !'/))
  call Message%printWarning('EMTBBFDF',(/' Normal absorption length may be too large...'/))
 endif

! finally, produce output images (BF and DF next to each other) 
! output image is in the current folder and called BFDF.tiffl it is up to the 
! user to change the filename if multiple runs are needed.

fname = trim(EMsoft%generateFilePath('EMdatapathname'))//'BFDF.tiff'
TIFF_filename = trim(fname)

! allocate memory for image
TIFF_nx = 2*jdim+2
TIFF_ny = jdim
allocate(TIFF_image(TIFF_nx,TIFF_ny))

! fill the image with whatever data you have (between 0 and 255)
do j=1,jdim
 do i=1,jdim
  TIFF_image(i,jdim+1-j) = int(255 * BF(i,j))
  TIFF_image(jdim+2+i,jdim+1-j) = int(255 * DF(i,j))
 end do
end do

! set up the image_t structure
im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMTBBFDF"," failed to convert array to image")

! create the file
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write image to file : "//iomsg)
else
  call Message%printMessage(' BFDF images written to '//trim(TIFF_filename))
end if

call mem%dealloc(BF, 'BF')
call mem%dealloc(DF, 'DF')
call mem%dealloc(z, 'z')
call mem%dealloc(exe, 'exe')
deallocate(TIFF_image)

end program EMTBBFDF
