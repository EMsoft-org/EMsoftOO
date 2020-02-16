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

module MODkvectors
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/13/20
  !!
  !! perform a series of unit tests on the kvectors module 

use stringconstants
use mod_global
use mod_kinds

contains 

subroutine MODkvectorsExecuteTest(res) &
           bind(c, name='MODkvectorsExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODkvectorsExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_kinds
use mod_kvectors
use mod_diffraction 
use mod_EMsoft
use mod_crystallography
use mod_symmetry
use HDF5
use mod_HDFsupport

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

type(kvectors_T)                :: kvec 
type(kvectorlist), pointer      :: klist
type(Diffraction_T)             :: Diff
type(EMsoft_T)                  :: EMsoft
type(cell_T)                    :: cell 
type(SpaceGroup_T)              :: SG

character(fnlen)                :: progname = 'MODkvectorsTest.f90'
character(fnlen)                :: progdesc = 'test program for the mod_kvectors module'
integer(kind=irg)               :: numk, i
character(fnlen)                :: mapmode, xtalname 
logical                         :: ok
real(kind=dbl)                  :: kV 

!===================================================
! set the reference values (verified with Mathematica scripts)




! initialize the error identifier to zero (should remain zero upon successful exit)
res = 0

!===================================================
!===================Start Tests=====================
!===================================================

! execute the constructor
kvec = kvectors_T() 
klist => kvec%get_ListHead()
numk = kvec%get_numk() 

! make sure that klist%next is nullified and empty
if (associated(klist%next)) then 
  res = 1
  write (*,"(' klist%next is associated; should be nullified upon construction ')")
  return
end if

if (numk.ne.0) then 
  res = 2
  write (*,"(' initial list should be empty ')")
  return
end if

! a few mapmode tests 
call kvec%set_mapmode('RoscaLambert')
mapmode = kvec%get_mapmode()
if (trim(mapmode).ne.'RoscaLambert') then 
  res = 3
  write (*,"(' mapmode incorrectly returned ')")
  return
end if

mapmode = 'incorrect'
ok = kvec%check_mapmode(mapmode)
if (ok) then 
  res = 4
  write (*,"(' nonexisting mapmode accepted as correct ')")
  return
end if

mapmode = 'Conical'
ok = kvec%check_mapmode(mapmode)
if (.not.ok) then 
  res = 5
  write (*,"(' existing mapmode returned as invalid ')")
  return
end if

!=======================================================
!=======================================================
!=======================================================
! next we do some tests that require an actual crystal structure to be present...
! open the HDF interface
call openFortranHDFInterface()

! set the XtalFolder to the resource folder location 
EMsoft = EMsoft_T( progname, progdesc, silent = .TRUE.)
call EMsoft%setConfigParameter('EMXtalFolderpathname', trim(EMsoft%getConfigParameter('Resourcepathname')) ) 

! read the test crystal structure 
xtalname = 'test_triclinic.xtal'
call cell%getCrystalData(xtalname, SG, EMsoft, verbose = .TRUE.)

! set the accelerating voltage
Diff = Diffraction_T( 20.D0, cell )

! and close the HDF interface 
call closeFortranHDFInterface() 
!=======================================================
!=======================================================
!=======================================================

! now we are ready to some more testing

! add a single k-vector to the list; we'll use the triclinic test structure in the resources folder 
call kvec%set_mapmode('RoscaLambert')
! call kvec%Calckvectors(cell, SG, Diff, )






end subroutine MODkvectorsExecuteTest

end module MODkvectors