! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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
program EMqg
  !! author: MDG
  !! version: 1.0 
  !! date: 01/27/20
  !!
  !! compute the complex extinction distance q_g

use mod_kinds
use mod_global
use mod_EMsoft
use mod_io
use mod_symmetry
use mod_crystallography
use mod_diffraction 
use mod_misc 

IMPLICIT NONE

character(fnlen)               :: progname = 'EMqg.f90'
character(fnlen)               :: progdesc = 'Display potential coefficient values'

type(EMsoft_T)                 :: EMsoft 
type(IO_T)                     :: Message
type(Cell_T)                   :: cell 
type(SpaceGroup_T)             :: SG
type(Diffraction_T)            :: Diff

integer(kind=irg)              :: ind(3),ans, oi_int(3)
real(kind=sgl)                 :: oi_real(7)
complex(kind=sgl)              :: oi_cmplx(1)
real(kind=sgl)                 :: preg, dmin, gstepsize
real(kind=dbl)                         :: eps = 1.0D-6
character(fnlen)               :: xtalname
character(200)                 :: parta
logical                        :: isallowed
type(gnode)                    :: rlp

EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 920 /) )

! ask for the crystal structure file
call Message%ReadValue(' Enter xtal file name : ', xtalname,"(A)")
call cell%getCrystalData(xtalname, SG, EMsoft)
call cell%calcPositions(SG, 'v')

call Diff%GetVoltage(cell)
preg = 2.0 * sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
call Diff%setrlpmethod('WK')

 ans = 1
 do while (ans.eq.1)
  call Message%printMessage('Enter Miller indices :', frm = "(/A)")
  call GetIndex(SG%getSpaceGrouphexset(), ind, 'r')
  call Diff%CalcUcg(cell, ind)
  isallowed = SG%IsGAllowed(ind)
  rlp = Diff%getrlp() 

! check whether this is a lattice extinction, a symmetry extinction, or an allowed reflection
  if ((rlp%Umod.lt.eps).and.(isallowed.eqv..FALSE.)) then
    call Message%printMessage('This reflection is absent due to lattice centering.')
    CYCLE
  end if
  if ((rlp%Umod.lt.eps).and.(isallowed.eqv..TRUE.)) then
    call Message%printMessage('This reflection is absent due to glide or screw symmetry elements.')
    CYCLE
  end if
  if (rlp%Umod.ne.0.D0) then
    parta = '   h  k  l    |g|    Ucg_r     Ucg_i      |Ug|      phase'// &
            '     |Ugp|     phase     xi_g   xi_gp    ratio    Re-1/q_g-Im'  
    call Message%printMessage(parta, frm = "(200A)") 

    oi_int(1:3) = rlp%hkl(1:3)
    call Message%WriteValue('',oi_int, 3, "(1x,3I3,1x)",advance="no")
    oi_real(1) = rlp%g
    call Message%WriteValue('',oi_real, 1, "(F9.4)",advance="no")
    oi_cmplx(1) = rlp%Ucg
    call Message%WriteValue('',oi_cmplx, 1, "(2F10.6,1x)",advance="no")
    oi_real(1:7)  = (/ rlp%Umod,rlp%Vphase*180.0/sngl(cPi),rlp%Upmod,rlp%Vpphase*180.0/sngl(cPi),rlp%xg,rlp%xgp,rlp%ar /)
    call Message%WriteValue('',oi_real, 7, "(4F10.5,3F8.1)",advance="no")
    oi_cmplx(1) = rlp%qg
    call Message%WriteValue('',oi_cmplx, 1, "(2F8.5)")
  end if

  call Message%ReadValue(' Another one ? (1/0) : ', oi_int, 1)
  ans = oi_int(1)
 end do 

end program EMqg
       



