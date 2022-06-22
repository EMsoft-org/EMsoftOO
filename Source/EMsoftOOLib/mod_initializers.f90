! ###################################################################
! Copyright (c) 2014-2022, Marc De Graef Research Group/Carnegie Mellon University
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

ddt = 1.0e-10

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

!--------------------------------------------------------------------------
recursive subroutine Initialize_QCCell(QCcell, QCDiff, QCSG, EMsoft, xtalname, dmin_qc, dmin_p, &
                                       voltage, nthreads, verbose, initLUT)
!DEC$ ATTRIBUTES DLLEXPORT :: Initialize_QCCell

use mod_io
use mod_EMsoft
use mod_QCcrystallography
use mod_QCsymmetry
use mod_QCdiffraction
use omp_lib

IMPLICIT NONE

class(QCcell_T),INTENT(INOUT)              :: QCcell
type(QCDiffraction_T),INTENT(INOUT)        :: QCDiff
type(QCspacegroup_T),INTENT(INOUT)         :: QCSG
type(EMsoft_T),INTENT(INOUT)               :: EMsoft
character(fnlen),INTENT(IN)                :: xtalname
real(kind=sgl),INTENT(IN)                  :: dmin_qc
real(kind=sgl),INTENT(IN)                  :: dmin_p
real(kind=dbl),INTENT(IN)                  :: voltage
integer(kind=irg),INTENT(IN)               :: nthreads
logical,INTENT(IN),OPTIONAL                :: verbose
logical,INTENT(IN),OPTIONAL                :: initLUT

type(IO_T)                                 :: Message 

integer(kind=irg)                          :: m, istat, skip, QCindex, nLUT
integer(kind=irg)                          :: imh, imk, ia1, ia2, ia3, ia4, ia5, ia6, id, TID
integer(kind=irg),allocatable              :: io_int(:), gg(:)
real(kind=sgl)                             :: dhkl, ddt
real(kind=sgl),allocatable                 :: io_real(:)
logical                                    :: loadingfile, justinit
complex(kind=dbl)                          :: Ucg
complex(kind=dbl)                          :: qg
real(kind=dbl)                             :: Vmod
real(kind=dbl)                             :: Vpmod, Upz
real(kind=dbl)                             :: xig
real(kind=dbl)                             :: xgp

justinit = .FALSE.
if(present(initLUT)) then
  if(initLUT) justinit = .TRUE.
end if

if(.not. justinit) then
! load the crystal structure file, which also computes all the important 
! matrices as well as all the symmetry arrays
 call QCcell%setfname(xtalname)
 call QCcell%ReadQCDataHDF(QCSG, EMsoft)
end if

! initialize the accelerating voltage
QCDiff = QCDiffraction_T( dble(voltage), QCcell, QCSG, verbose=.TRUE. )

! always use Weickenmeier&Kohl scattering coefficients, including absorptive form factors
call QCDiff%QCCalcWaveLength(QCcell, QCSG, verbose)

! compute the range of reflections for the lookup table and allocate the table
! The master list is easily created by brute force
select type (QCcell)
  class is (QCcell_axial_T)
    m = 5
    call QCcell%setdmin( dmin_qc, dmin_p )
    call QCcell%get_imax( imh, imk )
    imh = imh/2
    imk = imk/2

    if(.not. justinit) then
       imh = 1
       do 
          dhkl = 1.D0/QCcell%getvectorLength((/imh,0,0,0,0/), 'P', 'r')
          if (dhkl.lt.dmin_qc) EXIT
          imh = imh + 1
        end do

       imk = 1
       do 
          dhkl = 1.D0/QCcell%getvectorLength((/0,0,0,0,imk/), 'P', 'r')
          if (dhkl.lt.dmin_p) EXIT
          imk = imk + 1
       end do
       call QCcell%set_imax( 2*imh, 2*imk )
    end if

    if (present(verbose)) then
      if (verbose) then
        io_int(1:2) = (/imh,imk/)
        call Message%WriteValue(' Number of reflections along a*_i | i = {1,2,3,4} and a*_5 = ',io_int,2)
      end if
    end if

    nLUT = QCcell%getnDIndex( (/2*imh, 2*imh, 2*imh, 2*imh, 2*imk/) ) 

  class is (QCcell_icosahedral_T)
    m = 6
    call QCcell%setdmin( dmin_qc )
    imh = QCcell%get_imax()/2

    if(.not. justinit) then
       imh = 1
       do 
         dhkl = 1.D0/QCcell%getvectorLength( (/imh,0,0,0,0,0/), 'P', 'r')
         if (dhkl.lt.dmin_qc) EXIT
         imh = imh + 1
       end do
       call QCcell%set_imax( 2*imh )
    end if     

    nLUT = QCcell%getnDIndex( (/2*imh, 2*imh, 2*imh, 2*imh, 2*imh, 2*imh /) ) 

end select

allocate( gg(m), io_int(m), io_real(m) )

if(.not. justinit) then
  call QCDiff%allocateLUT(nLUT)
  select type (QCcell)
    class is (QCcell_axial_T)
     call QCcell%allocateinverseIndex(nLUT)
    class is (QCcell_icosahedral_T)
     call QCcell%allocateinverseIndex(nLUT)
   end select
end if 

ddt = 1.0e-5 

if (present(verbose)) then
  if (verbose) then
   call Message%printMessage('Generating Fourier coefficient lookup table ... ', frm = "(/A)",advance="no")
  end if
end if

! first, we deal with the transmitted beam
select type (QCcell)
  class is (QCcell_axial_T)
    gg  = (/ 0,0,0,0,0 /)
    Ucg = QCDiff%CalcQCUcg(QCcell, QCSG, gg, qg, Vmod, Vpmod, xig, xgp) 
    Upz = Vpmod         ! U'0 normal absorption parameter 
    id  = QCcell%getnDIndex(gg)
! and add this reflection to the look-up table
    call QCDiff%setLUT(id, Ucg)
    call QCDiff%setLUTqg(id, qg)
 
    if(.not. justinit) then
! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!  
     ia1l: do ia1 = -2*imh,2*imh
       ia2l: do ia2 = -2*imh,2*imh
         ia3l: do ia3 = -2*imh,2*imh
           ia4l: do ia4 = -2*imh, 2*imh
             ia5l: do ia5 = -2*imk, 2*imk
                gg  = (/ ia1, ia2, ia3, ia4, ia5 /)
                id  = QCcell%getnDIndex(gg)
                call QCcell%setinverseIndex(id, gg(1:5))
              end do ia5l
            end do ia4l
         end do ia3l
        end do ia2l
      end do ia1l
    end if

  class is (QCcell_icosahedral_T)
    gg  = (/ 0,0,0,0,0,0 /)
    Ucg = QCDiff%CalcQCUcg(QCcell, QCSG, gg, qg, Vmod, Vpmod, xig, xgp) 
    Upz = Vpmod         ! U'0 normal absorption parameter 
    id  = QCcell%getnDIndex(gg)
! and add this reflection to the look-up table
    call QCDiff%setLUT(id, Ucg)
    call QCDiff%setLUTqg(id, qg)

    if(.not. justinit) then
! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!  
     ia1r: do ia1 = -2*imh,2*imh
       ia2r: do ia2 = -2*imh,2*imh
         ia3r: do ia3 = -2*imh,2*imh
           ia4r: do ia4 = -2*imh, 2*imh
             ia5r: do ia5 = -2*imh, 2*imh
               ia6r: do ia6 = -2*imh, 2*imh
                  gg = (/ ia1, ia2, ia3, ia4, ia5, ia6 /)
                  id = QCcell%getnDIndex(gg)
                  call QCcell%setinverseIndex(id, gg(1:6))
                end do ia6r
              end do ia5r
            end do ia4r
         end do ia3r
        end do ia2r
      end do ia1r
    end if

end select

! set the number of OpenMP threads 
call OMP_SET_NUM_THREADS(nthreads)
if (present(verbose)) then
     if (verbose) then
        io_int(1) = nthreads
        call Message%WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")
     end if
end if

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(TID, gg, Ucg, io_real, qg, Vmod, VPmod, xig, xgp)

TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC)

do id = 1,nLUT
  select type (QCcell)
    class is (QCcell_axial_T)
     gg = QCCell%getinverseIndex(id)
    class is (QCcell_icosahedral_T)
     gg = QCCell%getinverseIndex(id)
  end select

  Ucg = QCDiff%CalcQCUcg(QCcell, QCSG, gg, qg, Vmod, Vpmod, xig, xgp)
  call QCDiff%setLUT(id, Ucg)
  call QCDiff%setLUTqg(id, qg)

! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
   if (abs(Ucg).le.ddt) then 
     call QCDiff%setdbdiff(id, .TRUE.)
   end if

   if(present(verbose)) then
     if(verbose) then
       if(mod(id,250000) .eq. 0) then
         io_real(1) = 100.D0 * dble(id)/dble(nLUT)
         call Message%WriteValue(' Finished computing ',io_real, 1, '(F8.2, " % of the total coefficients ")')
       end if
     end if
   end if
end do
! end of OpenMP portion
!$OMP END DO
!$OMP END PARALLEL

if (present(verbose)) then
  if (verbose) then
    call Message%printMessage('Done', frm = "(A/)")
  end if
end if

! that's it
end subroutine Initialize_QCCell

end module mod_initializers
