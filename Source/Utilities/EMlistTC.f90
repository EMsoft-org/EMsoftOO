! ###################################################################
! Copyright (c) 2015-2020, Marc De Graef Research Group/Carnegie Mellon University
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

program EMlistTC
  !! author: MDG
  !! version: 1.0 
  !! date: 01/23/20
  !!
  !! List information for texture components in an arbitrary crystal system

use mod_kinds
use mod_global
use mod_EMsoft
use mod_symmetry
use mod_io
use mod_crystallography
use mod_rotations
use mod_quaternions
use mod_so3
use mod_misc 

IMPLICIT NONE

character(fnlen)               :: progname = 'EMlistTC.f90'
character(fnlen)               :: progdesc = 'List information about texture components'

type(EMsoft_T)                 :: EMsoft
type(IO_T)                     :: Message
type(Cell_T)                   :: cell 
type(SpaceGroup_T)             :: SG 
type(so3_T)                    :: SO
type(QuaternionArray_T)        :: Pm, dummy
type(Orientation_T)            :: ot
type(o_T)                      :: o
type(r_T)                      :: r 

integer(kind=irg)              :: io_int3(3), io_int4(4), io_int7(7), nvec(3), bvec4(4), bvec(3)
real(kind=dbl)                 :: cnvec(3), cbvec(3), t(3), om(3,3)
character(fnlen)               :: xtalname
integer(kind=irg), allocatable :: ntmp(:,:), btmp(:,:)
integer(kind=irg)              :: i, j, k, nnum, bnum, ntc, pgnum
integer(kind=irg), allocatable :: ortho(:,:)


! print some information
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 916 /) )

! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft)

! allocate the symmetry operations for Fundamental Zone handling.
pgnum = SG%getPGnumber()
SO = so3_T( pgnum )
call dummy%QSym_Init( pgnum, Pm )

! ask for the texture component (distinguish between hexagonal indices and regular)
call Message%printMessage('  ')
if (SG%getSpaceGrouphexset().eqv..TRUE.) then 
  call Message%printMessage(' Enter the texture component symbol {hk.l}<uvtw>')
  call Message%printMessage('  (This structure uses the hexagonal indexing system!) ')
  call Message%ReadValue('   Enter the three-index planar indices {hkl} ',io_int3,3)
  nvec = io_int3(1:3)
  call Message%ReadValue('   Enter the four-index direction indices <uvtw> ',io_int4,4)
  bvec4 = io_int4(1:4)
! convert four-index to three-index 
  call MilBrav(bvec,bvec4,'43')
  io_int3(1:3) = bvec(1:3)
  call Message%WriteValue('   ---> in three index notation: ',io_int3,3,"('[',I4,I4,I4,']')")
else
  call Message%printMessage(' Enter the texture component symbol {hkl}<uvw>')
  call Message%ReadValue('   Enter the planar indices {hkl} ',io_int3,3)
  nvec = io_int3(1:3)
  call Message%ReadValue('   Enter the direction indices <uvw> ',io_int3,3)
  bvec = io_int3(1:3)
end if 

! In determining the list of equivalent texture components, we need to make sure that we only keep
! those for which the plane vector is perpendicular to the direction vector !
call Message%printMessage('  ')
call Message%printMessage(' Determining unique texture components ')
call SG%CalcFamily(nvec,nnum,'r',ntmp)
call SG%CalcFamily(bvec,bnum,'d',btmp)
io_int3(1) = nnum
call Message%WriteValue('    Multiplicity of {hkl} : ', io_int3,1,"(I2)")
io_int3(1) = bnum
call Message%WriteValue('    Multiplicity of <uvw> : ', io_int3,1,"(I2)")

! next, take every single combination of {hkl} and <uvw> and keep only those that are orthogonal
allocate(ortho(nnum,bnum))
ortho = 0
do i=1,nnum
  do j=1,bnum
    if (sum(ntmp(i,1:3)*btmp(j,1:3)).eq.0) ortho(i,j) = 1
  end do
end do
! how many are there ?
ntc = sum(ortho)
io_int3(1) = ntc
call Message%WriteValue(' There are ',io_int3,1,"(I4, ' equivalent texture components.')")
call Message%printMessage('  ')

! next, get the rotation representation for all of them in all the rotation representations
k = 0
do i=1,nnum
  do j=1,bnum
    if (ortho(i,j).eq.1) then 
      k = k+1
      io_int3(1) = k
      call Message%printMessage('-------')
      call Message%WriteValue(' TC ',io_int3,1,"(I3,'  -->  ')",advance="no")
! transform the two vectors to the Cartesian crystallographic reference frame
      call cell%TransSpace(dble(ntmp(i,1:3)), cnvec, 'r', 'c') 
      call cell%TransSpace(dble(btmp(j,1:3)), cbvec, 'd', 'c') 
! normalize them
      call cell%NormVec(cnvec, 'c')
      call cell%NormVec(cbvec, 'c')
! get the cross product n x b and normalize
      call cell%CalcCross(cnvec,cbvec,t,'c','c',0)
      call cell%NormVec(t, 'c')
! get these into a rotation matrix with b,t,n in the columns (in that order)
      om(1:3,1) = cbvec(1:3)
      om(1:3,2) = t(1:3)
      om(1:3,3) = cnvec(1:3)
      o = o_T( odinp = om )
      ot = Orientation_T(o)
! and print the information to the command line
      if (SG%getSpaceGrouphexset().eqv..FALSE.) then
        io_int7(1:3) = ntmp(i,1:3)
        io_int7(4:6) = btmp(j,1:3)
        call Message%WriteValue('',io_int7,6,"('(',I3,I3,I3,') [',I3, I3, I3,']')")
      else
        io_int7(1:3) = ntmp(i,1:3)
        call MilBrav(btmp(j,1:3),io_int7(4:7),'34')
        call Message%WriteValue('',io_int7,7,"('(',I3,I3,' . ',I3,') [',I3, I3, I3, I3,']')")
      end if
      call ot%print_orientation('d')
! is this one inside the fundamental zone ?
      r = ot%getClass_r()
      if (SO%IsinsideFZ(r).eqv..TRUE.) then 
        call Message%printMessage('  ---> THIS ROTATION LIES INSIDE THE FUNDAMENTAL ZONE.',"(A/)")
      end if 
    end if
  end do
end do

end program EMlistTC
