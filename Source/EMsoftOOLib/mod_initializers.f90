! ###################################################################
! Copyright (c) 2014-2021, Marc De Graef Research Group/Carnegie Mellon University
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
module mod_initializers
  !! author: MDG
  !! version: 1.0
  !! date: 02/05/20
  !!
  !! various initializer routines.
  !!
  !! Due to some complicated module interdependencies these routines are in
  !! this module rather than elsewhere...

use mod_kinds
use mod_global

IMPLICIT NONE

contains

!--------------------------------------------------------------------------
recursive subroutine Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, &
                                      useHDF, initLUT, noLUT, interpolate)
!DEC$ ATTRIBUTES DLLEXPORT :: Initialize_Cell
 !! author: MDG
 !! version: 1.0
 !! date: 02/05/20
 !!
 !! perform all steps to initialize a unit Diff type variable

use mod_EMsoft
use mod_symmetry
use mod_io
use mod_gvectors
use mod_diffraction
use mod_crystallography
use mod_HDFsupport

IMPLICIT NONE

type(Cell_T),INTENT(INOUT)                 :: cell
type(Diffraction_T),INTENT(INOUT)          :: Diff
type(SpaceGroup_T),INTENT(INOUT)           :: SG
type(DynType),INTENT(INOUT)                :: Dyn
type(EMsoft_T),INTENT(INOUT)               :: EMsoft
real(kind=sgl),INTENT(IN)                  :: dmin
logical,INTENT(IN),OPTIONAL                :: verbose
type(HDF_T),OPTIONAL,INTENT(INOUT)         :: useHDF
logical,INTENT(IN),OPTIONAL                :: initLUT
logical,INTENT(IN),OPTIONAL                :: noLUT
logical,INTENT(IN),OPTIONAL                :: interpolate

type(IO_T)                                 :: Message
type(gnode)                                :: rlp
integer(kind=irg)                          :: istat, io_int(3), skip
integer(kind=irg)                          :: imh, imk, iml, gg(3), ix, iy, iz
real(kind=sgl)                             :: dhkl, io_real(3), ddt
logical                                    :: loadingfile, justinit, interp, compute
real(kind=sgl),parameter                   :: gstepsize = 0.001  ! [nm^-1] interpolation stepsize
character(fnlen)                           :: xtalname

interp = .FALSE.
if (present(interpolate)) then
  if (interpolate) then
    interp = .TRUE.
    call Diff%setrlpmethod('IP')
  end if
end if

justinit = .FALSE.
if(present(initLUT)) then
  if(initLUT) justinit = .TRUE.
end if

compute = .TRUE.
if(present(noLUT)) then
  if(noLUT) compute= .FALSE.
end if

if(.not. justinit) then
! load the crystal structure file, which also computes all the important
! matrices as well as all the symmetry arrays
 xtalname = trim(cell%getFileName())
 if (present(useHDF)) then
   call cell%getCrystalData(xtalname, SG, EMsoft, verbose, useHDF)
 else
   call cell%getCrystalData(xtalname, SG, EMsoft, verbose)
 end if
end if

! we assume that the wavelength has already been set in the Diff class...
call Diff%CalcWaveLength(cell)

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
if (compute) then
 imh = 1
 do
   dhkl = 1.0/cell%CalcLength( (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imh = imh + 1
 end do
 imk = 1
 do
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,float(imk),0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imk = imk + 1
 end do
 iml = 1
 do
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
   iml = iml + 1
 end do

 if (present(verbose)) then
  if (verbose) then
    io_int = (/ imh, imk, iml /)
    call Message%WriteValue(' Range of reflections along a*, b* and c* = ',io_int,3)
  end if
 end if

! do we need to pre-compute the scattering factors and store them in Diff%scatfac ?
 if (interp.eqv..TRUE.) then
   call Diff%PreCalcFSCATT(cell, dmin, gstepsize)
 end if

if(.not. justinit) then
  call Diff%allocateLUT( (/ imh, imk, iml /) )
end if

ddt = 1.0e-5

! next, we compute the overall lookup table Diff%LUT; we do not, at this point, create a
! list of linked reflections; in the old code, this was done at the same time, but it appears
! it is better to decouple these two computations. In this new approach, we'll compute a much
! shorter linked list based on the incident wave vector direction.


! first, we deal with the transmitted beam
 gg = (/ 0,0,0 /)
 if (interp.eqv..TRUE.) then
   call Diff%CalcUcg(cell,gg,applyqgshift=.TRUE.,interpolate=.TRUE.)
 else
   call Diff%CalcUcg(cell,gg,applyqgshift=.TRUE.)
 end if
 rlp = Diff%getrlp()

 Dyn%Upz = rlp%Vpmod         ! U'0 normal absorption parameter
 if (present(verbose)) then
  if (verbose) then
   io_real(1) = rlp%xgp
   call Message%WriteValue(' Normal absorption length [nm] = ', io_real, 1)
  end if
 end if

! and add this reflection to the look-up table
 call Diff%setLUT( gg, rlp%Ucg )
 call Diff%setLUTqg( gg, rlp%qg )

 if (present(verbose)) then
  if (verbose) then
   call Message%printMessage(' Generating Fourier coefficient lookup table ... ', frm = "(/A)",advance="no")
  end if
 end if

! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!
ixl: do ix=-2*imh,2*imh
iyl:  do iy=-2*imk,2*imk
izl:   do iz=-2*iml,2*iml
        gg = (/ ix, iy, iz /)
        if (SG%IsGAllowed(gg)) then  ! is this reflection allowed by lattice centering ?
! add the reflection to the look up table
           if (interp.eqv..TRUE.) then
             call Diff%CalcUcg(cell,gg,applyqgshift=.TRUE.,interpolate=.TRUE.)
           else
             call Diff%CalcUcg(cell,gg,applyqgshift=.TRUE.)
           end if
           rlp = Diff%getrlp()
           call Diff%setLUT( gg, rlp%Ucg )
           call Diff%setLUTqg( gg, rlp%qg )
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
           if (cabs(rlp%Ucg).le.ddt) then
             call Diff%setdbdiff( gg, .TRUE. )
           end if
        end if ! IsGAllowed
       end do izl
      end do iyl
    end do ixl

  if (present(verbose)) then
   if (verbose) then
    call Message%printMessage('Done', frm = "(A/)")
   end if
  end if
 end if

! that's it
end subroutine Initialize_Cell



!--------------------------------------------------------------------------
recursive subroutine Initialize_Cell_noHDF(cell, Diff, SG, Dyn, dmin)
!DEC$ ATTRIBUTES DLLEXPORT :: Initialize_Cell_NoHDF
 !! author: MDG
 !! version: 1.0
 !! date: 04/15/21
 !!
 !! perform all steps to initialize a unit cell without needing to read an .xtal file
 !! we assume that all the atom coordinates and symmetry operators have already been created

use mod_symmetry
use mod_io
use mod_gvectors
use mod_diffraction
use mod_crystallography

IMPLICIT NONE

type(Cell_T),INTENT(INOUT)                 :: cell
type(Diffraction_T),INTENT(INOUT)          :: Diff
type(SpaceGroup_T),INTENT(INOUT)           :: SG
type(DynType),INTENT(INOUT)                :: Dyn
real(kind=sgl),INTENT(IN)                  :: dmin

type(IO_T)                                 :: Message
type(gnode)                                :: rlp
integer(kind=irg)                          :: istat, io_int(3), skip
integer(kind=irg)                          :: imh, imk, iml, gg(3), ix, iy, iz
real(kind=sgl)                             :: dhkl, io_real(3), ddt
real(kind=sgl),parameter                   :: gstepsize = 0.001  ! [nm^-1] interpolation stepsize

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
 imh = 1
 do
   dhkl = 1.0/cell%CalcLength( (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imh = imh + 1
 end do
 imk = 1
 do
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,float(imk),0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
   imk = imk + 1
 end do
 iml = 1
 do
   dhkl = 1.0/cell%CalcLength((/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
   iml = iml + 1
 end do

! we assume that the wavelength has already been set in the Diff class...
 call Diff%CalcWaveLength(cell)

 call Diff%allocateLUT( (/ imh, imk, iml /) )

ddt = 1.0e-5

! next, we compute the overall lookup table Diff%LUT; we do not, at this point, create a
! list of linked reflections; in the old code, this was done at the same time, but it appears
! it is better to decouple these two computations. In this new approach, we'll compute a much
! shorter linked list based on the incident wave vector direction.

! first, we deal with the transmitted beam
 gg = (/ 0,0,0 /)
 call Diff%CalcUcg(cell,gg,applyqgshift=.TRUE.)
 rlp = Diff%getrlp()
 Dyn%Upz = rlp%Vpmod         ! U'0 normal absorption parameter

! and add this reflection to the look-up table
 call Diff%setLUT( gg, rlp%Ucg )
 call Diff%setLUTqg( gg, rlp%qg )

! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!
ixl: do ix=-2*imh,2*imh
iyl:  do iy=-2*imk,2*imk
izl:   do iz=-2*iml,2*iml
        gg = (/ ix, iy, iz /)
        if (SG%IsGAllowed(gg)) then  ! is this reflection allowed by lattice centering ?
! add the reflection to the look up table
           call Diff%CalcUcg(cell,gg,applyqgshift=.TRUE.)
           rlp = Diff%getrlp()
           call Diff%setLUT( gg, rlp%Ucg )
           call Diff%setLUTqg( gg, rlp%qg )
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
           if (cabs(rlp%Ucg).le.ddt) then
             call Diff%setdbdiff( gg, .TRUE. )
           end if
        end if ! IsGAllowed
       end do izl
      end do iyl
    end do ixl

! that's it
end subroutine Initialize_Cell_NoHDF



end module mod_initializers
