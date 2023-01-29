! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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

program EMstar
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! Prints the star of a reciprocal lattice position for a given space group;
  !!
  !! also lists the Fourier coefficients of the lattice potential

use mod_EMsoft
use mod_global
use mod_io
use mod_symmetry
use mod_crystallography
use mod_diffraction
use mod_HDFsupport

IMPLICIT NONE

type(EMsoft_T)              :: EMsoft 
type(IO_T)                  :: Message
type(Cell_T)                :: cell
type(SpaceGroup_T)          :: SG
type(Diffraction_T)         :: Diff 

character(fnlen)            :: progname = 'EMstar.f90'
character(fnlen)            :: progdesc = 'Computes the star of a reciprocal lattice vector'

integer(kind=irg)           :: g(3),gg(3),ans,n,i, io_int(3)
real(kind=dbl)              :: kk(3)
real(kind=dbl),allocatable  :: stmp(:,:)
logical                     :: first
character(1)                :: space
character(fnlen)            :: xtalname

EMsoft = EMsoft_T(progname, progdesc, tpl = (/ 924 /) )

call openFortranHDFInterface()

! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft)

call Diff%setrlpmethod('WK')
call Diff%getVoltage(cell, verbose=.TRUE.) 

ans = 1

do while (ans.eq.1)
  call Message%printMessage(' ',"(/A)")
  call Message%ReadValue(' Enter reciprocal lattice coordinates [I] : ', io_int, 3)
  g(1:3) = io_int(1:3)
  kk = dble(g)
  call SG%CalcStar(kk, n, stmp, space)
  io_int(1) = n
  call Message%WriteValue(' Number of equivalent planes in star = ', io_int, 1, "(I3)")
! compute and display the structure factor
  first = .TRUE.
  do i=1,n
   gg(1:3) = nint(stmp(i,1:3))
   call Diff%CalcUcg(cell, gg)
   call Diff%Printrlp(first)
  end do
  call Message%ReadValue(' Another star ? (1/0) ', io_int, 1)
  ans = io_int(1)
end do

end program EMstar
