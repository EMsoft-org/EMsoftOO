! ###################################################################
! Copyright (c) 2016-2023, Marc De Graef Research Group/Carnegie Mellon University
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

!--------------------------------------------------------------------------
! EMsoft:HDFcharTest.f90
!--------------------------------------------------------------------------
!
! MODULE: HDFcharTest
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief test module for writing and reading of real(sgl) to/from HDF5 files
!
!> @date 10/29/16   MDG 1.0 original
!--------------------------------------------------------------------------

module HDFcharTest

use stringconstants
use mod_global
use mod_kinds
contains 

subroutine HDFcharExecuteTest(res) &
           bind(c, name='HDFcharExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: HDFcharExecuteTest

use,INTRINSIC :: ISO_C_BINDING
use mod_EMsoft
use HDF5
use mod_HDFsupport
IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res

character(fnlen)                :: HDFfilename, tmppath, groupname, dataset, textfile, progname, progdesc
character(1)                    :: EMsoftnativedelimiter

integer(kind=irg)               :: i1, i2, i3, i4, dim1, dim2, dim3, dim4, hdferr, isum 
real(kind=sgl)                  :: fval, fval_save
integer(HSIZE_T)                :: dims1(1), dims2(2), dims3(3), dims4(4)
character(len=1),allocatable    :: carr1(:), carr2(:,:), carr3(:,:,:), carr4(:,:,:,:)
character(len=1),allocatable    :: carr1_save(:), carr2_save(:,:), carr3_save(:,:,:), carr4_save(:,:,:,:)

type(HDF_T)                     :: HDF
type(EMsoft_T)                  :: EMsoft

progname = 'HDFcharTest'
progdesc = 'Test program for character array reading and writing using the HDF class'
EMsoft = EMsoft_T(progname, progdesc)

!====================================
! generate the real arrays
dim1 = 5
dim2 = 10
dim3 = 15
dim4 = 20

dims1 = (/ dim1 /)
dims2 = (/ dim1, dim2 /)
dims3 = (/ dim1, dim2, dim3 /)
dims4 = (/ dim1, dim2, dim3, dim4 /)


ALLOCATE (carr1(dim1))
ALLOCATE (carr2(dim1,dim2))
ALLOCATE (carr3(dim1,dim2,dim3))
ALLOCATE (carr4(dim1,dim2,dim3,dim4))

do i1=1,dim1
  carr1(i1) = char(i1)
  do i2=1,dim2
    carr2(i1,i2) = char(mod(i1 * i2,128))
    do i3=1,dim3
      carr3(i1,i2,i3) = char(mod(i1 * i2 * i3,128))
      do i4=1,dim4
        carr4(i1,i2,i3,i4) = char(mod(i1 * i2 * i3 * i4,128))
      end do
    end do
  end do
end do

ALLOCATE (carr1_save(dim1))
ALLOCATE (carr2_save(dim1,dim2))
ALLOCATE (carr3_save(dim1,dim2,dim3))
ALLOCATE (carr4_save(dim1,dim2,dim3,dim4))

carr1_save = carr1
carr2_save = carr2
carr3_save = carr3
carr4_save = carr4

!====================================
! initialize the HDF class
HDF = HDF_T()

! determine the pathname delimiter character
EMsoftnativedelimiter = EMsoft%getConfigParameter('EMsoftnativedelimiter')

! get the location of the Temporary folder inside the Build folder (it always exists)
tmppath = EMsoft%getConfigParameter('EMsofttestpath')

! create and open the test file
HDFfilename = trim(tmppath)//EMsoftnativedelimiter//'HDFtest_char.h5'

write(*,*) 'writing filename = <'//trim(HDFfilename)//'>, has length ',len(trim(HDFfilename))
write (*,*) 'cstring : <'//cstringify(HDFfilename)//'> has length ',len(cstringify(HDFfilename))

hdferr =  HDF%createFile(HDFfilename)
if (hdferr.ne.0) then
  res = 1
  return
end if

! write the char arrays to the file
dataset = SC_char1D
hdferr = HDF%writeDatasetCharArray(dataset, carr1, dims1)
if (hdferr.ne.0) then
  res = 3
  return
end if

dataset = SC_char2D
hdferr = HDF%writeDatasetCharArray(dataset, carr2, dims2)
if (hdferr.ne.0) then
  res = 4
  return
end if

dataset = SC_char3D
hdferr = HDF%writeDatasetCharArray(dataset, carr3, dims3)
if (hdferr.ne.0) then
  res = 5
  return
end if

dataset = SC_char4D
hdferr = HDF%writeDatasetCharArray(dataset, carr4, dims4)
if (hdferr.ne.0) then
  res = 6
  return
end if

call HDF%pop(.TRUE.)

!====================================


!====================================
! deallocate the integer arrays (they will be recreated upon reading)
deallocate( carr1, carr2, carr3, carr4)
!====================================

!====================================
! next, we read the data sets from the HDF5 file

! open the test file
HDFfilename = trim(tmppath)//EMsoftnativedelimiter//'HDFtest_char.h5'

write(*,*) 'reading filename = <'//trim(HDFfilename)//'>, has length ',len(trim(HDFfilename))
write (*,*) 'cstring : <'//cstringify(HDFfilename)//'> has length ',len(cstringify(HDFfilename))

hdferr =  HDF%openFile(HDFfilename)
if (hdferr.ne.0) then
  res = 7
  return
end if

! read the char arrays
dataset = SC_char1D
call HDF%readDatasetCharArray(dataset, dims1, hdferr, carr1)
if (hdferr.ne.0) then
  res = 9
  return
end if

dataset = SC_char2D
call HDF%readDatasetCharArray(dataset, dims2, hdferr, carr2)
if (hdferr.ne.0) then
  res = 10 
  return
end if

dataset = SC_char3D
call HDF%readDatasetCharArray(dataset, dims3, hdferr, carr3)
if (hdferr.ne.0) then
  res = 11
  return
end if

dataset = SC_char4D
call HDF%readDatasetCharArray(dataset, dims4, hdferr, carr4)
if (hdferr.ne.0) then
  res = 12
  return
end if

call HDF%pop(.TRUE.)

!====================================


!====================================
! compare the entries with the stored values
isum = 20
if (all(carr1.eq.carr1_save).eqv..FALSE.) isum = isum + 2
if (all(carr2.eq.carr2_save).eqv..FALSE.) isum = isum + 4
if (all(carr3.eq.carr3_save).eqv..FALSE.) isum = isum + 8
if (all(carr4.eq.carr4_save).eqv..FALSE.) isum = isum + 16

if (isum.eq.20) then
  res = 0
else
  res = isum
end if

!====================================

end subroutine HDFcharExecuteTest


end module HDFcharTest
