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

module mod_PGA3Dsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/21/21
  !!
  !! support routines for Projective Geometric Algebra in 3D; this module is based on 
  !! https://bivector.net/3DPGA.pdf 
  !!
  !! as of 7/23/2021 the following operations are available
  !!
  !! returning a multivector:
  !! function rotor(angle, line) result(rtr)
  !! function translator(dist, line) result(trs)
  !! function plane(a, b, c, d) result(pl)
  !! function point(x, y, z, ideal) result(pt)
  !! function line(a, b, c, ideal) result(l)
  !! function circle(t, radius, line) result(cr)
  !! function torus(s, t, r1, l1, r2, l2) result(to)
  !! function pointoncircle(t, radius, line) result(poc)
  !! function join(mv1, mv2, mv3, mv4) result(j)
  !! function meet(mv1, mv2, mv3, mv4) result(j)
  !! function normalthru(mv1, mv2) result(p)
  !! function normalto(mv1) result(n)
  !! function projectonto(mv1, mv2) result(p)
  !!
  !! returnig a scalar result:
  !! function distpoints(mv1, mv2) result(d)
  !! function distplanes(mv1, mv2) result(d)
  !! function distplaneline(mv1, mv2) result(d)
  !! function ordisttoplane(mv1, mv2) result(d)
  !! function ordisttoline(mv1, mv2) result(d)
  !! function angleplanes(mv1, mv2) result(d)
  !! function angleplaneline(mv1, mv2) result(d)
  !!
  !! miscellaneous:
  !! subroutine getplane(pl, a, b, c, d) 
  !! subroutine getpoint(pt, x, y, z)


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
!DEC$ ATTRIBUTES DLLEXPORT :: PGA3D_initialize
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
recursive subroutine getplane(pl, a, b, c, d) 
!DEC$ ATTRIBUTES DLLEXPORT :: getplane
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! get the components of a plane

IMPLICIT NONE 

type(PGA3D_T), INTENT(INOUT)         :: pl
real(kind=dbl), INTENT(OUT)          :: a, b, c, d 

a = pl%getcomp(2)
b = pl%getcomp(3)
c = pl%getcomp(4)
d = pl%getcomp(1)

end subroutine getplane

!--------------------------------------------------------------------------
recursive function point(x, y, z, ideal) result(pt)
!DEC$ ATTRIBUTES DLLEXPORT :: point
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! create a regular or an ideal point

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: x, y, z
logical, INTENT(IN), OPTIONAL       :: ideal 
type(PGA3D_T)                       :: pt

pt = (E032.muls.x) + (E013.muls.y) + (E021.muls.z) 
if (.not.present(ideal)) pt = pt + E123

end function point

!--------------------------------------------------------------------------
recursive subroutine getpoint(pt, x, y, z)
!DEC$ ATTRIBUTES DLLEXPORT :: getpoint
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! get the coordinates of a point

IMPLICIT NONE 

type(PGA3D_T), INTENT(INOUT)         :: pt
real(kind=dbl), INTENT(OUT)          :: x, y, z

x = pt%getcomp(13)
y = pt%getcomp(12)
z = pt%getcomp(11)

end subroutine getpoint

!--------------------------------------------------------------------------
recursive function line(a, b, c, ideal) result(l)
!DEC$ ATTRIBUTES DLLEXPORT :: line
!! author: MDG 
!! version: 1.0 
!! date: 07/21/21
!!
!! create a regular or an ideal line

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: a, b, c
logical, INTENT(IN), OPTIONAL       :: ideal 
type(PGA3D_T)                       :: l

if (present(ideal)) then 
  l = (E01.muls.a) + (E02.muls.b) + (E03.muls.c) 
else
  l = (E12.muls.c) + (E31.muls.b) + (E23.muls.a) 
end if 

end function line

!--------------------------------------------------------------------------
recursive function circle(t, radius, line) result(cr)
!DEC$ ATTRIBUTES DLLEXPORT :: circle
!! author: MDG 
!! version: 1.0 
!! date: 07/22/21
!!
!! create a circle

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: t
real(kind=dbl), INTENT(IN)          :: radius 
type(PGA3D_T), INTENT(IN)           :: line
type(PGA3D_T)                       :: cr

cr = rotor(t * cPi * 2.D0, line) * translator( radius, E1*E0 )

end function circle

!--------------------------------------------------------------------------
recursive function torus(s, t, r1, l1, r2, l2) result(to)
!DEC$ ATTRIBUTES DLLEXPORT :: torus
!! author: MDG 
!! version: 1.0 
!! date: 07/22/21
!!
!! create a circle

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: s, t, r1, r2 
type(PGA3D_T), INTENT(IN)           :: l1, l2
type(PGA3D_T)                       :: to

to = circle(s, r2, l2) * circle(t, r1, l1)

end function torus

!--------------------------------------------------------------------------
recursive function pointoncircle(t, radius, line) result(poc)
!DEC$ ATTRIBUTES DLLEXPORT :: pointoncircle
!! author: MDG 
!! version: 1.0 
!! date: 07/22/21
!!
!! temporary routine

IMPLICIT NONE 

real(kind=dbl), INTENT(IN)          :: t
real(kind=dbl), INTENT(IN)          :: radius 
type(PGA3D_T), INTENT(IN)           :: line
type(PGA3D_T)                       :: poc

type(PGA3D_T)                       :: cr

cr = circle(t, radius, line)
poc = cr .sandwich. E123

end function pointoncircle

!--------------------------------------------------------------------------
recursive function join(mv1, mv2, mv3, mv4) result(j)
!DEC$ ATTRIBUTES DLLEXPORT :: join
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! join two, three, or four multivectors

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2 
type(PGA3D_T), INTENT(IN),OPTIONAL  :: mv3, mv4 
type(PGA3D_T)                       :: j

j = mv1.vee.mv2 

if (present(mv3)) j = j.vee.mv3
if (present(mv4)) j = j.vee.mv4

end function join

!--------------------------------------------------------------------------
recursive function meet(mv1, mv2, mv3, mv4) result(j)
!DEC$ ATTRIBUTES DLLEXPORT :: meet
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! meet two, three, or four multivectors 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2 
type(PGA3D_T), INTENT(IN),OPTIONAL  :: mv3, mv4 
type(PGA3D_T)                       :: j

j = mv1.wedge.mv2 

if (present(mv3)) j = j.wedge.mv3
if (present(mv4)) j = j.wedge.mv4

end function meet

!--------------------------------------------------------------------------
recursive function normalthru(mv1, mv2) result(p)
!DEC$ ATTRIBUTES DLLEXPORT :: normalthru
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! find multivector normal to mv1 through mv2

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2 
type(PGA3D_T)                       :: p

p = mv1.inner.mv2 

end function normalthru

!--------------------------------------------------------------------------
recursive function normalto(mv1) result(n)
!DEC$ ATTRIBUTES DLLEXPORT :: normalto
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! find multivector normal to mv1 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1 
type(PGA3D_T)                       :: n

n = mv1 * E0123

end function normalto

!--------------------------------------------------------------------------
recursive function projectonto(mv1, mv2) result(p)
!DEC$ ATTRIBUTES DLLEXPORT :: projectonto
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! project mv1 onto mv2 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
type(PGA3D_T)                       :: p

p = (mv2.inner.mv1) * mv2

end function projectonto

!--------------------------------------------------------------------------
recursive function distpoints(mv1, mv2) result(d)
!DEC$ ATTRIBUTES DLLEXPORT :: distpoints
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! distance bwtween to points mv1 and mv2 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
real(kind=dbl)                      :: d

type(PGA3D_T)                       :: p

p = (mv1%normalized()).vee.(mv2%normalized())
d = p%norm()

end function distpoints

!--------------------------------------------------------------------------
recursive function distplanes(mv1, mv2) result(d)
!DEC$ ATTRIBUTES DLLEXPORT :: distplanes
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! distance between parallel planes mv1 and mv2 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
real(kind=dbl)                      :: d

type(PGA3D_T)                       :: p

p = (mv1%normalized()).wedge.(mv2%normalized())
d = p%inorm()

end function distplanes

!--------------------------------------------------------------------------
recursive function distplaneline(mv1, mv2) result(d)
!DEC$ ATTRIBUTES DLLEXPORT :: distplaneline
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! distance between parallel plane mv1 and line mv2 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
real(kind=dbl)                      :: d

type(PGA3D_T)                       :: p
real(kind=dbl), allocatable         :: v(:)
integer(kind=irg)                   :: i

p = mv1%normalized() * mv2%normalized() 
v = p%getgrade(3)

d = norm2(v(1:3))

end function distplaneline

!--------------------------------------------------------------------------
recursive function ordisttoplane(mv1, mv2) result(d)
!DEC$ ATTRIBUTES DLLEXPORT :: ordisttoplane
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! oriented distance between point mv1 and plane mv2 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
real(kind=dbl)                      :: d

type(PGA3D_T)                       :: p

p = (mv1%normalized()) .wedge. (mv2%normalized()) 

d = p%getcomp(15)

end function ordisttoplane

!--------------------------------------------------------------------------
recursive function ordisttoline(mv1, mv2) result(d)
!DEC$ ATTRIBUTES DLLEXPORT :: ordisttoline
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! oriented distance between point mv1 and line mv2 

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
real(kind=dbl)                      :: d

type(PGA3D_T)                       :: p

p = (mv1%normalized()) .vee. (mv2%normalized()) 

d = p%norm()

end function ordisttoline

!--------------------------------------------------------------------------
recursive function angleplanes(mv1, mv2) result(d)
!DEC$ ATTRIBUTES DLLEXPORT :: angleplanes
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! angle between intersecting planes

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
real(kind=dbl)                      :: d

type(PGA3D_T)                       :: p

p = (mv1%normalized()).inner.(mv2%normalized())
d = acos(p%getcomp(0))

end function angleplanes

!--------------------------------------------------------------------------
recursive function angleplaneline(mv1, mv2) result(d)
!DEC$ ATTRIBUTES DLLEXPORT :: angleplaneline
!! author: MDG 
!! version: 1.0 
!! date: 07/23/21
!!
!! angle between plane and non-parallel line

IMPLICIT NONE 

type(PGA3D_T), INTENT(IN)           :: mv1, mv2
real(kind=dbl)                      :: d

type(PGA3D_T)                       :: p
real(kind=dbl), allocatable         :: v(:)

p = mv1%normalized() * mv2%normalized() 
v = p%getgrade(3)

d = asin( norm2(v(2:4)) )

end function angleplaneline




end module mod_PGA3Dsupport
