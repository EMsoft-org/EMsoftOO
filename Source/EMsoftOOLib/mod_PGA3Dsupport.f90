! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_PGA3Dsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/21/21
  !!
  !! support routines for Projective Geometric Algebra in 3D 
  !!

use mod_kinds
use mod_global
use mod_PGA3D 

IMPLICIT NONE 

! define the basis components of the algebra
type(PGA3D_T)             :: ONE, E0, E1, E2, E3, E01, E02, E03, E12, E23, E31, E123, E032, E021, E013, E0123, &
                             E023, E012, E031 

contains 

!--------------------------------------------------------------------------
subroutine PGA3D_initialize(verbose) 
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! define the basis elements of the projective geometric algebra in 3D 

IMPLICIT NONE

logical, INTENT(IN), OPTIONAL     :: verbose 


ONE = PGA3D_T(val=1.D0, ind=0) 

E0 = PGA3D_T(val=1.D0, ind=1)           ! ideal plane
E1 = PGA3D_T(val=1.D0, ind=2)           ! x=0 plane
E2 = PGA3D_T(val=1.D0, ind=3)           ! y=0 plane
E3 = PGA3D_T(val=1.D0, ind=4)           ! z=0 plane

E01 = E0 .wedge. E1 
E02 = E0 .wedge. E2 
E03 = E0 .wedge. E3 
E12 = E1 .wedge. E2 
E31 = E3 .wedge. E1 
E23 = E2 .wedge. E3 

E123 = E1 .wedge. E2 .wedge. E3
E032 = E0 .wedge. E3 .wedge. E2
E013 = E0 .wedge. E1 .wedge. E3
E021 = E0 .wedge. E2 .wedge. E1
E023 = E032 .muls. (-1.D0)
E031 = E013 .muls. (-1.D0)
E012 = E021 .muls. (-1.D0)

E0123 = E0 .wedge. E1 .wedge. E2 .wedge. E3 

if (present(verbose)) then 
  if (verbose.eqv..TRUE.) then 
    call ONE%log( pre = 'ONE')
    call E0%log( pre = 'E0')
    call E1%log( pre = 'E1')
    call E2%log( pre = 'E2')
    call E3%log( pre = 'E3')
    call E01%log( pre = 'E01')
    call E02%log( pre = 'E02')
    call E03%log( pre = 'E03')
    call E12%log( pre = 'E12')
    call E23%log( pre = 'E23')
    call E31%log( pre = 'E31')
    call E123%log( pre = 'E123')
    call E032%log( pre = 'E032')
    call E013%log( pre = 'E013')
    call E021%log( pre = 'E021')
    call E0123%log( pre = 'E0123')
  end if
end if 

end subroutine PGA3D_initialize

!--------------------------------------------------------------------------
recursive function rotor(angle, line) result(rtr)
!DEC$ ATTRIBUTES DLLEXPORT :: rotor
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! create a rotor multivector

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: angle
type(PGA3D_T), INTENT(IN)           :: line
type(PGA3D_T)                       :: rtr

type(PGA3D_T)                       :: b

b = line%normalized()
b = b.muls.sin(angle*0.5D0)
rtr = b.adds.cos(angle*0.5D0)

end function rotor

!--------------------------------------------------------------------------
recursive function translator(dist, line) result(trs)
!DEC$ ATTRIBUTES DLLEXPORT :: translator
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! create a translator multivector

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: dist
type(PGA3D_T), INTENT(IN)           :: line
type(PGA3D_T)                       :: trs

type(PGA3D_T)                       :: b

b = line
b = b.muls.(dist*0.5D0)
trs = b.adds.(1.D0)

end function translator

!--------------------------------------------------------------------------
recursive function plane(a, b, c, d) result(pl)
!DEC$ ATTRIBUTES DLLEXPORT :: plane
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! create a plane

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: a, b, c, d 
type(PGA3D_T)                       :: pl

pl = (E1.muls.a) + (E2.muls.b) + (E3.muls.c) + (E0.muls.d)

end function plane

!--------------------------------------------------------------------------
recursive function point(x, y, z) result(pt)
!DEC$ ATTRIBUTES DLLEXPORT :: point
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! create a point

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: x, y, z
type(PGA3D_T)                       :: pt

pt = (E032.muls.x) + (E013.muls.y) + (E021.muls.z) + E123

end function point



end module mod_PGA3Dsupport