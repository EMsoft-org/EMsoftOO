! ###################################################################
! Copyright (c) 2013-2026, Marc De Graef Research Group/Carnegie Mellon University
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


module mod_productgrid
  !! author: ZTV
  !! version: 1.0
  !! date: 11/02/26
  !!
  !! product grids on the Hopf fibration for every proper point group
  !!
  !! Z.T. Varley and M. De Graef, "Uniform crystal orientation grids via the Hopf fibration"
  !! (submitted to J. Appl. Cryst.)
  !!
  !! The orientation matrix of an EMsoft quaternion q takes sample to crystal coordinates, so its third
  !! column pi(q) = q z q* is the sample normal in the crystal frame (the inverse pole figure direction),
  !! and the crystal symmetry acts on the left, pi(g q) = g pi(q). The projection pi is a Hopf fibration
  !! that is compatible with every point group K; its fibres are the circles { q R_z(t) } of rotations
  !! about the sample normal. A node v of the sphere with a stabilizer of order k carries a fibre on which
  !! the stabilizer identifies angles that differ by 2 pi / k.
  !!
  !! Three root grids (O, I, D6) are built directly from a K-invariant mesh of the sphere:
  !!
  !! - O and I: the faces of a regular octahedron or icosahedron subdivided at frequency f and mapped to
  !!   the sphere by an area-preserving Knothe-Rosenblatt map in the symmetry chambers of the polyhedron;
  !!   D6: equal-area latitude rings;
  !! - one representative v of every K-orbit of nodes, taken in the sector at z so that the sections of
  !!   neighbouring representatives are nearly aligned, carries n / k evenly spaced rotations;
  !! - alternate layers of nodes are staggered by half a fibre step.
  !!
  !! There are then n rotations on the fibre over every node and no orientation is listed twice:
  !! N_FZ = n |V| / |K|. Every other group H inherits the grid of a root group P that contains it: the
  !! P-extended set is H-invariant in a suitable setting, and one orientation per H-orbit is kept, so that
  !! N_FZ grows by |P| / |H|.
  !!
  !! Internally the routines use the active quaternion (the conjugate of the EMsoft quaternion), for which
  !! the symmetry acts on the right; productGrid converts on output. All quaternion products use the
  !! standard (Hamilton) definition, which is the one mod_quaternions uses for epsijk = +1.

use mod_kinds
use mod_global

IMPLICIT NONE
private

! number of azimuths in the tables of the area-preserving map
integer(kind=irg), parameter :: ntab = 65537

type :: areaMap_T
  !! the symmetry chamber of a polyhedron and the azimuthal tables of its area-preserving map
  real(kind=dbl), allocatable :: ops(:,:,:)      ! symmetry operations of the polyhedron
  real(kind=dbl), allocatable :: rotp(:,:)       ! the chamber's base point, rotated by each operation
  real(kind=dbl), allocatable :: planes(:,:)     ! chamber planes that do not contain the pole
  real(kind=dbl)              :: normal(3)       ! normal of the face that holds the chamber
  real(kind=dbl)              :: lo, hi          ! azimuthal range of the chamber
  real(kind=dbl), allocatable :: phi(:)          ! azimuths
  real(kind=dbl), allocatable :: scdf(:)         ! cumulative area of the flat face
  real(kind=dbl), allocatable :: sder(:)         ! its derivative at the azimuths
  real(kind=dbl), allocatable :: tcdf(:)         ! cumulative solid angle of the sphere
  real(kind=dbl), allocatable :: tder(:)         ! derivative of the azimuth with respect to tcdf
end type areaMap_T

public :: productGridCount, productGrid

contains

!--------------------------------------------------------------------------
pure function qmul(a, b) result(c)
  !! Hamilton product of two quaternions

IMPLICIT NONE

real(kind=dbl), INTENT(IN)  :: a(4), b(4)
real(kind=dbl)              :: c(4)

c(1) = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4)
c(2) = a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3)
c(3) = a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2)
c(4) = a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1)

end function qmul

!--------------------------------------------------------------------------
pure function qconj(a) result(c)
  !! quaternion conjugate

IMPLICIT NONE

real(kind=dbl), INTENT(IN)  :: a(4)
real(kind=dbl)              :: c(4)

c = (/ a(1), -a(2), -a(3), -a(4) /)

end function qconj

!--------------------------------------------------------------------------
pure function qrot(axis, turns) result(q)
  !! quaternion of a rotation by 2 pi / turns around axis

IMPLICIT NONE

real(kind=dbl), INTENT(IN)  :: axis(3), turns
real(kind=dbl)              :: q(4)

q(1) = cos(cPi/turns)
q(2:4) = axis * sin(cPi/turns) / sqrt(sum(axis*axis))

end function qrot

!--------------------------------------------------------------------------
pure function qmat(q) result(m)
  !! rotation matrix of an active quaternion

IMPLICIT NONE

real(kind=dbl), INTENT(IN)  :: q(4)
real(kind=dbl)              :: m(3,3)

m(1,1) = 1.D0 - 2.D0*(q(3)*q(3) + q(4)*q(4))
m(1,2) = 2.D0*(q(2)*q(3) - q(1)*q(4))
m(1,3) = 2.D0*(q(2)*q(4) + q(1)*q(3))
m(2,1) = 2.D0*(q(2)*q(3) + q(1)*q(4))
m(2,2) = 1.D0 - 2.D0*(q(2)*q(2) + q(4)*q(4))
m(2,3) = 2.D0*(q(3)*q(4) - q(1)*q(2))
m(3,1) = 2.D0*(q(2)*q(4) - q(1)*q(3))
m(3,2) = 2.D0*(q(3)*q(4) + q(1)*q(2))
m(3,3) = 1.D0 - 2.D0*(q(2)*q(2) + q(3)*q(3))

end function qmat

!--------------------------------------------------------------------------
subroutine closeGroup(gens, el)
  !! close a set of quaternion generators under multiplication; identity first, scalar part >= 0

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: gens(:,:)
real(kind=dbl), allocatable, INTENT(OUT)  :: el(:,:)

real(kind=dbl)                            :: found(4,128), p(4)
integer(kind=irg)                         :: i, j, cnt

cnt = 1
found(1:4,1) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
i = 1
do while (i.le.cnt)
  do j=1,size(gens,2)
    p = qmul(found(1:4,i), gens(1:4,j))
    if (maxval(abs(matmul(p, found(1:4,1:cnt)))).lt.1.D0-1.D-9) then
      cnt = cnt + 1
      found(1:4,cnt) = p
    end if
  end do
  i = i + 1
end do

allocate(el(4,cnt))
do i=1,cnt
  el(1:4,i) = found(1:4,i) * merge(1.D0, -1.D0, found(1,i).ge.0.D0)
end do

end subroutine closeGroup

!--------------------------------------------------------------------------
subroutine groupElements(group, el)
  !! the elements of a rotation group: C_k with the k-fold axis along z, D_k in addition with a
  !! two-fold axis along x, and I with five-fold axes along z and (2,0,1) (the IcoVertices setting)

IMPLICIT NONE

character(2), INTENT(IN)                  :: group
real(kind=dbl), allocatable, INTENT(OUT)  :: el(:,:)

real(kind=dbl)                            :: g(4,3), z(3), x(3), k

z = (/ 0.D0, 0.D0, 1.D0 /)
x = (/ 1.D0, 0.D0, 0.D0 /)
select case (group(1:1))
  case ('C')
    read (group(2:2),"(F1.0)") k
    g(1:4,1) = qrot(z, k)
    call closeGroup(g(1:4,1:1), el)
  case ('D')
    read (group(2:2),"(F1.0)") k
    g(1:4,1) = qrot(z, k)
    g(1:4,2) = qrot(x, 2.D0)
    call closeGroup(g(1:4,1:2), el)
  case ('T')
    g(1:4,1) = qrot(z, 2.D0)
    g(1:4,2) = qrot(x, 2.D0)
    g(1:4,3) = qrot((/ 1.D0, 1.D0, 1.D0 /), 3.D0)
    call closeGroup(g, el)
  case ('O')
    g(1:4,1) = qrot(z, 4.D0)
    g(1:4,2) = qrot((/ 1.D0, 1.D0, 1.D0 /), 3.D0)
    call closeGroup(g(1:4,1:2), el)
  case ('I')
    g(1:4,1) = qrot(z, 5.D0)
    g(1:4,2) = qrot((/ 2.D0, 0.D0, 1.D0 /), 5.D0)
    call closeGroup(g(1:4,1:2), el)
end select

end subroutine groupElements

!--------------------------------------------------------------------------
subroutine groupSetting(FZtype, FZorder, FZordersign, group, root, embed, setting, info)
  !! the rotation group of a fundamental zone, the root grid it inherits, and two rotations:
  !! embed places the group inside its root (the group's axes on rotation axes of the root; the
  !! icosahedral grid is used wherever the group embeds in I), and setting takes the standard
  !! setting of the group to the one of the EMsoft fundamental zone

IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: FZtype, FZorder, FZordersign
character(2), INTENT(OUT)       :: group, root
real(kind=dbl), INTENT(OUT)     :: embed(4), setting(4)
integer(kind=irg), INTENT(OUT)  :: info

real(kind=dbl), parameter       :: e2(4) = (/ 0.96193835778391745D0, 0.D0, 0.27326652891267172D0, 0.D0 /), &
                                   e3(4) = (/ 0.94597537489177308D0, 0.17451645846277514D0, &
                                              0.26873177999583048D0, -0.049576468654520228D0 /), &
                                   id(4) = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
real(kind=dbl)                  :: z(3)

info = 0
z = (/ 0.D0, 0.D0, 1.D0 /)
embed = id
setting = id
select case (FZtype)
  case (0)
    group = 'C1'
    root = 'I '
  case (1)
    select case (abs(FZorder))
      case (2)
        group = 'C2'
        root = 'I '
        embed = e2
! monoclinic groups have the two-fold axis along y, and mm2 along z
        if (FZordersign.eq.1) setting = qrot((/ 1.D0, 0.D0, 0.D0 /), -4.D0)
      case (3)
        group = 'C3'
        root = 'I '
        embed = e3
      case (4)
        group = 'C4'
        root = 'O '
      case (6)
        group = 'C6'
        root = 'D6'
      case default
        info = 1
    end select
  case (2)
    select case (abs(FZorder))
      case (2)
        group = 'D2'
        root = 'I '
        embed = e2
! -4m2 is 222 rotated by 45 degrees around z
        if (FZordersign.eq.-1) setting = qrot(z, 8.D0)
      case (3)
        group = 'D3'
        root = 'I '
        embed = e3
! 312, -31m and -62m are 32 rotated by 30 degrees around z
        if (FZordersign.eq.-1) setting = qrot(z, 12.D0)
      case (4)
        group = 'D4'
        root = 'O '
      case (6)
        group = 'D6'
        root = 'D6'
      case default
        info = 1
    end select
  case (3)
    group = 'T '
    root = 'I '
    embed = e2
  case (4)
    group = 'O '
    root = 'O '
  case (5)
    group = 'I '
    root = 'I '
  case default
    info = 1
end select

end subroutine groupSetting

!--------------------------------------------------------------------------
pure function fibreCount(root, f) result(n)
  !! rotations per full fibre of a root grid at frequency f

IMPLICIT NONE

character(2), INTENT(IN)        :: root
integer(kind=irg), INTENT(IN)   :: f
integer(kind=irg)               :: n

select case (root)
  case ('O ')
    n = 4*f
    if (f.ge.6) n = 4*f - 12*((2*f+9)/18)
  case ('I ')
    n = 5*f
  case default
    n = 6*((2*f+2)/3)
end select

end function fibreCount

!--------------------------------------------------------------------------
pure function ringSize(f, j) result(m)
  !! population of ring j of the hexagonal sphere grid (j = f is the equator)

IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: f, j
integer(kind=irg)               :: m, mult

mult = merge(12, 6, j.eq.f)
m = max(mult, mult*floor(4.D0*dble(f)*sin(cPi*dble(j)/dble(2*f))/dble(mult) + 0.5D0))

end function ringSize

!--------------------------------------------------------------------------
pure function nodeCount(root, f) result(nv)
  !! number of nodes of a root's sphere grid

IMPLICIT NONE

character(2), INTENT(IN)        :: root
integer(kind=irg), INTENT(IN)   :: f
integer(kind=ill)               :: nv

integer(kind=irg)               :: j

select case (root)
  case ('O ')
    nv = 4_ill*f*f + 2_ill
  case ('I ')
    nv = 10_ill*f*f + 2_ill
  case default
    nv = 2_ill + int(ringSize(f, f), ill)
    do j=1,f-1
      nv = nv + 2_ill*ringSize(f, j)
    end do
end select

end function nodeCount

!--------------------------------------------------------------------------
recursive function productGridCount(FZtype, FZorder, FZordersign, f) result(cnt)
!DEC$ ATTRIBUTES DLLEXPORT :: productGridCount
  !! author: ZTV
  !! version: 1.0
  !! date: 11/02/26
  !!
  !! number of zone points N_FZ = n |V| / |K| of the product grid at frequency f; 0 if there is no such grid

IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: FZtype, FZorder, FZordersign, f
integer(kind=ill)               :: cnt

character(2)                    :: group, root
real(kind=dbl)                  :: embed(4), setting(4)
real(kind=dbl), allocatable     :: el(:,:)
integer(kind=irg)               :: info

cnt = 0_ill
call groupSetting(FZtype, FZorder, FZordersign, group, root, embed, setting, info)
if ((info.ne.0).or.(f.lt.1)) return
call groupElements(group, el)
cnt = int(fibreCount(root, f), ill) * nodeCount(root, f) / int(size(el,2), ill)

end function productGridCount

!--------------------------------------------------------------------------
subroutine pchipDerivatives(x, y, d)
  !! derivatives of the monotone piecewise cubic Hermite interpolant (Fritsch-Butland)

IMPLICIT NONE

real(kind=dbl), INTENT(IN)      :: x(:), y(:)
real(kind=dbl), INTENT(OUT)     :: d(:)

real(kind=dbl), allocatable     :: h(:), m(:)
real(kind=dbl)                  :: w1, w2
integer(kind=irg)               :: i, n

n = size(x)
allocate(h(n-1), m(n-1))
h = x(2:n) - x(1:n-1)
m = (y(2:n) - y(1:n-1)) / h

d(1) = pchipEdge(h(1), h(2), m(1), m(2))
d(n) = pchipEdge(h(n-1), h(n-2), m(n-1), m(n-2))
do i=2,n-1
  if ((m(i-1)*m(i)).le.0.D0) then
    d(i) = 0.D0
  else
    w1 = 2.D0*h(i) + h(i-1)
    w2 = h(i) + 2.D0*h(i-1)
    d(i) = (w1 + w2) / (w1/m(i-1) + w2/m(i))
  end if
end do

contains

  pure function pchipEdge(h0, h1, m0, m1) result(e)
  real(kind=dbl), INTENT(IN)    :: h0, h1, m0, m1
  real(kind=dbl)                :: e
  e = ((2.D0*h0 + h1)*m0 - h0*m1) / (h0 + h1)
  if ((e*m0).le.0.D0) then
    e = 0.D0
  else if (((m0*m1).le.0.D0).and.(abs(e).gt.3.D0*abs(m0))) then
    e = 3.D0*m0
  end if
  end function pchipEdge

end subroutine pchipDerivatives

!--------------------------------------------------------------------------
pure function pchipValue(x, y, d, xi) result(v)
  !! value of the interpolant at xi, which must lie in the range of x

IMPLICIT NONE

real(kind=dbl), INTENT(IN)      :: x(:), y(:), d(:), xi
real(kind=dbl)                  :: v

real(kind=dbl)                  :: h, t
integer(kind=irg)               :: lo, hi, mid

lo = 1
hi = size(x)
do while (hi-lo.gt.1)
  mid = (lo+hi)/2
  if (xi.ge.x(mid)) then
    lo = mid
  else
    hi = mid
  end if
end do

h = x(hi) - x(lo)
t = (xi - x(lo)) / h
v = (1.D0+2.D0*t)*(1.D0-t)**2 * y(lo) + t*(1.D0-t)**2 * h * d(lo) + &
    t*t*(3.D0-2.D0*t) * y(hi) + t*t*(t-1.D0) * h * d(hi)

end function pchipValue

!--------------------------------------------------------------------------
pure subroutine cumulativeSimpson(y, h, c)
  !! cumulative integral of samples at spacing h, for an even number of intervals, with
  !! Simpson's rule on every pair of intervals, split between them by the same parabola

IMPLICIT NONE

real(kind=dbl), INTENT(IN)      :: y(:), h
real(kind=dbl), INTENT(OUT)     :: c(:)

integer(kind=irg)               :: i

c(1) = 0.D0
do i=1,size(y)-2,2
  c(i+1) = c(i) + h/12.D0 * (5.D0*y(i) + 8.D0*y(i+1) - y(i+2))
  c(i+2) = c(i+1) + h/12.D0 * (-y(i) + 8.D0*y(i+1) + 5.D0*y(i+2))
end do

end subroutine cumulativeSimpson

!--------------------------------------------------------------------------
pure function thetaLimit(M, phi) result(th)
  !! largest polar angle of the chamber at azimuth phi

IMPLICIT NONE

type(areaMap_T), INTENT(IN)     :: M
real(kind=dbl), INTENT(IN)      :: phi
real(kind=dbl)                  :: th

real(kind=dbl)                  :: cp, sp, psi
integer(kind=irg)               :: i

cp = cos(phi)
sp = sin(phi)
psi = 0.D0
do i=1,size(M%planes,2)
  psi = max(psi, atan2(M%planes(3,i), M%planes(1,i)*cp + M%planes(2,i)*sp))
end do
th = cPi - psi

end function thetaLimit

!--------------------------------------------------------------------------
subroutine polyhedron(root, vv, faces)
  !! vertices of the regular polyhedron of a root grid and the vertex triplets of its faces

IMPLICIT NONE

character(2), INTENT(IN)                     :: root
real(kind=dbl), allocatable, INTENT(OUT)     :: vv(:,:)
integer(kind=irg), allocatable, INTENT(OUT)  :: faces(:,:)

real(kind=dbl)                               :: adjacent
integer(kind=irg)                            :: i, j, k, cnt, nf

if (root.eq.'O ') then
  allocate(vv(3,6))
  vv = 0.D0
  do i=1,3
    vv(i,i) = 1.D0
    vv(i,i+3) = -1.D0
  end do
  nf = 8
  adjacent = -0.1D0
else
  allocate(vv(3,12))
  vv = IcoVertices
  nf = 20
  adjacent = 0.4D0
end if

allocate(faces(3,nf))
cnt = 0
do i=1,size(vv,2)
  do j=i+1,size(vv,2)
    if (dot_product(vv(1:3,i),vv(1:3,j)).lt.adjacent) cycle
    do k=j+1,size(vv,2)
      if ((dot_product(vv(1:3,i),vv(1:3,k)).lt.adjacent).or.(dot_product(vv(1:3,j),vv(1:3,k)).lt.adjacent)) cycle
      cnt = cnt + 1
      faces(1:3,cnt) = (/ i, j, k /)
    end do
  end do
end do

end subroutine polyhedron

!--------------------------------------------------------------------------
subroutine areaMapInit(root, vv, faces, M)
  !! the chamber of the polyhedron and the azimuthal cumulatives of the flat face and the sphere

IMPLICIT NONE

character(2), INTENT(IN)        :: root
real(kind=dbl), INTENT(IN)      :: vv(:,:)
integer(kind=irg), INTENT(IN)   :: faces(:,:)
type(areaMap_T), INTENT(OUT)    :: M

real(kind=dbl), allocatable     :: el(:,:), tmp(:,:)
real(kind=dbl), allocatable     :: src(:), tgt(:)
real(kind=dbl)                  :: p(3), h(3), inside(3), n(3), a, t, th, dphi, tmid, mid, best, val
integer(kind=irg)               :: i, j, cnt, ng

call groupElements(root, el)
ng = size(el,2)
allocate(M%ops(3,3,2*ng), M%rotp(3,2*ng), tmp(3,2*ng))

! the base point of the chamber, whose spherical Voronoi cell under the symmetry operations it is
if (root.eq.'O ') then
  p = (/ 1.D0, 2.D0, 3.D0 /) / sqrt(14.D0)
  M%lo = 0.25D0*cPi
  M%hi = 0.5D0*cPi
else
  p = (/ sin(0.2D0)*cos(0.1D0*cPi), sin(0.2D0)*sin(0.1D0*cPi), cos(0.2D0) /)
  M%lo = 0.D0
  M%hi = 0.2D0*cPi
end if

! the rotations and their negatives, with near-integer entries snapped to the integers
do i=1,ng
  M%ops(1:3,1:3,i) = qmat(el(1:4,i))
  where (abs(M%ops(1:3,1:3,i)-anint(M%ops(1:3,1:3,i))).lt.1.D-10) M%ops(1:3,1:3,i) = anint(M%ops(1:3,1:3,i))
  M%ops(1:3,1:3,ng+i) = -M%ops(1:3,1:3,i)
end do

! each operation contributes the inward normal p - M^T p; only those that do not contain the pole limit theta
cnt = 0
do i=1,2*ng
  M%rotp(1:3,i) = matmul(transpose(M%ops(1:3,1:3,i)), p)
  h = p - M%rotp(1:3,i)
  a = sqrt(sum(h*h))
  if (a.le.1.D-9) cycle
  h = h/a
  where (abs(h).lt.1.D-12) h = 0.D0
  if (h(3).eq.0.D0) cycle
  cnt = cnt + 1
  tmp(1:3,cnt) = h
end do
allocate(M%planes(3,cnt))
M%planes = tmp(1:3,1:cnt)

! the face that holds the chamber
mid = 0.5D0*(M%lo+M%hi)
tmid = 0.5D0*thetaLimit(M, mid)
inside = (/ sin(tmid)*cos(mid), sin(tmid)*sin(mid), cos(tmid) /)
best = -2.D0
do i=1,size(faces,2)
  n = 0.D0
  do j=1,3
    n = n + vv(1:3,faces(j,i))
  end do
  n = n / sqrt(sum(n*n))
  val = dot_product(n, inside)
  if (val.gt.best) then
    best = val
    M%normal = n
  end if
end do

! flat-face area and solid angle per unit azimuth
allocate(src(ntab), tgt(ntab), M%phi(ntab), M%scdf(ntab), M%sder(ntab), M%tcdf(ntab), M%tder(ntab))
dphi = (M%hi-M%lo)/dble(ntab-1)
do i=1,ntab
  M%phi(i) = M%lo + dble(i-1)*dphi
end do
M%phi(ntab) = M%hi
do i=1,ntab
  th = thetaLimit(M, M%phi(i))
  a = M%normal(1)*cos(M%phi(i)) + M%normal(2)*sin(M%phi(i))
  t = tan(th)
  src(i) = (t/(M%normal(3) + a*t))**2
  tgt(i) = 2.D0*sin(0.5D0*th)**2
end do
call cumulativeSimpson(src, dphi, M%scdf)
call cumulativeSimpson(tgt, dphi, M%tcdf)
M%scdf = M%scdf/M%scdf(ntab)
M%tcdf = M%tcdf/M%tcdf(ntab)
call pchipDerivatives(M%phi, M%scdf, M%sder)
call pchipDerivatives(M%tcdf, M%phi, M%tder)

end subroutine areaMapInit

!--------------------------------------------------------------------------
pure function areaPreserving(M, v) result(out)
  !! Knothe-Rosenblatt map from a radially projected point of a flat face to the sphere

IMPLICIT NONE

type(areaMap_T), INTENT(IN)     :: M
real(kind=dbl), INTENT(IN)      :: v(3)
real(kind=dbl)                  :: out(3)

real(kind=dbl)                  :: w(3), best, score, phi, theta, top, a, t, tm, fraction, newphi, newtop, newtheta
integer(kind=irg)               :: g, ig

! carry the point into the chamber with the operation that brings it closest to the base point
ig = 1
best = dot_product(v, M%rotp(1:3,1))
do g=2,size(M%rotp,2)
  score = dot_product(v, M%rotp(1:3,g))
  if (score.gt.best) then
    best = score
    ig = g
  end if
end do
w = matmul(M%ops(1:3,1:3,ig), v)

phi = min(max(atan2(w(2), w(1)), M%lo), M%hi)
theta = atan2(sqrt(w(1)*w(1) + w(2)*w(2)), w(3))
top = thetaLimit(M, phi)
a = M%normal(1)*cos(phi) + M%normal(2)*sin(phi)
t = tan(theta)
tm = tan(top)

! the area fraction of the column is kept, and the column goes to the azimuth with the same cumulative
fraction = (t/(M%normal(3) + a*t) / (tm/(M%normal(3) + a*tm)))**2
newphi = pchipValue(M%tcdf, M%phi, M%tder, min(max(pchipValue(M%phi, M%scdf, M%sder, phi), 0.D0), 1.D0))
newtop = thetaLimit(M, newphi)
newtheta = 2.D0*asin(sqrt(min(max(fraction, 0.D0), 1.D0))*sin(0.5D0*newtop))

out = (/ sin(newtheta)*cos(newphi), sin(newtheta)*sin(newphi), cos(newtheta) /)
out = matmul(transpose(M%ops(1:3,1:3,ig)), out)
out = out / sqrt(sum(out*out))

end function areaPreserving

!--------------------------------------------------------------------------
subroutine polyhedronNodes(root, f, v, layer)
  !! nodes of the subdivided polyhedron on the sphere; the layer of a node with barycentric
  !! indices (i,j,k) counts rows of the flat face lattice from the nearest vertex, f - max(i,j,k)

IMPLICIT NONE

character(2), INTENT(IN)        :: root
integer(kind=irg), INTENT(IN)   :: f
real(kind=dbl), INTENT(OUT)     :: v(:,:)
integer(kind=irg), INTENT(OUT)  :: layer(:)

real(kind=dbl), allocatable     :: vv(:,:)
integer(kind=irg), allocatable  :: faces(:,:)
type(areaMap_T)                 :: M
real(kind=dbl)                  :: u(3)
integer(kind=irg)               :: a, b, c, i, j, n

call polyhedron(root, vv, faces)
call areaMapInit(root, vv, faces, M)

n = 0

! the polyhedron vertices
do a=1,size(vv,2)
  n = n + 1
  v(1:3,n) = areaPreserving(M, vv(1:3,a))
  layer(n) = 0
end do

! the interior nodes of the edges
do a=1,size(vv,2)
  do b=a+1,size(vv,2)
    if (.not.any((faces(1,:).eq.a.or.faces(2,:).eq.a.or.faces(3,:).eq.a).and. &
                 (faces(1,:).eq.b.or.faces(2,:).eq.b.or.faces(3,:).eq.b))) cycle
    do i=1,f-1
      u = dble(i)*vv(1:3,a) + dble(f-i)*vv(1:3,b)
      n = n + 1
      v(1:3,n) = areaPreserving(M, u/sqrt(sum(u*u)))
      layer(n) = f - max(i, f-i)
    end do
  end do
end do

! the interior nodes of the faces
do c=1,size(faces,2)
  do i=1,f-2
    do j=1,f-1-i
      u = dble(i)*vv(1:3,faces(1,c)) + dble(j)*vv(1:3,faces(2,c)) + dble(f-i-j)*vv(1:3,faces(3,c))
      n = n + 1
      v(1:3,n) = areaPreserving(M, u/sqrt(sum(u*u)))
      layer(n) = f - max(i, j, f-i-j)
    end do
  end do
end do

end subroutine polyhedronNodes

!--------------------------------------------------------------------------
subroutine ringNodes(f, v, layer)
  !! D6-invariant equal-area rings: both poles, 2f-2 rings and the equator. Every ring sits at the
  !! middle height of its band of area 4 pi m / |V|, and alternate rings are staggered by half a step;
  !! the layer of a node counts rings from the nearest pole

IMPLICIT NONE

integer(kind=irg), INTENT(IN)   :: f
real(kind=dbl), INTENT(OUT)     :: v(:,:)
integer(kind=irg), INTENT(OUT)  :: layer(:)

real(kind=dbl)                  :: top, width, z, r, phi
integer(kind=irg)               :: i, j, m, n, cnt

cnt = size(v,2)
v(1:3,1) = (/ 0.D0, 0.D0, 1.D0 /)
v(1:3,2) = (/ 0.D0, 0.D0, -1.D0 /)
layer(1:2) = 0
n = 2
top = 1.D0 - 2.D0/dble(cnt)
do j=1,f-1
  m = ringSize(f, j)
  width = 2.D0*dble(m)/dble(cnt)
  z = top - 0.5D0*width
  r = sqrt(1.D0 - z*z)
  do i=0,m-1
    phi = 2.D0*cPi*(dble(i) + 0.5D0*dble(mod(j,2)))/dble(m)
    v(1:3,n+1) = (/ r*cos(phi), r*sin(phi), z /)
    v(1:3,n+2) = (/ r*cos(phi), r*sin(phi), -z /)
    layer(n+1:n+2) = j
    n = n + 2
  end do
  top = top - width
end do
m = ringSize(f, f)
do i=0,m-1
  phi = 2.D0*cPi*dble(i)/dble(m)
  n = n + 1
  v(1:3,n) = (/ cos(phi), sin(phi), 0.D0 /)
  layer(n) = f
end do

end subroutine ringNodes

!--------------------------------------------------------------------------
subroutine rootGrid(root, f, q, info)
  !! the active quaternions of a root grid at frequency f, before reduction to the fundamental zone

IMPLICIT NONE

character(2), INTENT(IN)                  :: root
integer(kind=irg), INTENT(IN)             :: f
real(kind=dbl), allocatable, INTENT(OUT)  :: q(:,:)
integer(kind=irg), INTENT(OUT)            :: info

real(kind=dbl), allocatable               :: v(:,:), el(:,:), mats(:,:,:), poles(:,:), images(:,:)
integer(kind=irg), allocatable            :: layer(:), faces(:,:)
real(kind=dbl)                            :: w(3), s(4), rz(4), t, phase, az, zmax, sn, kx, ky, bx, by, sectorAngle
integer(kind=irg)                         :: i, g, j, k, n, nv, ng, zorder, cnt, ip
integer(kind=ill)                         :: total
logical                                   :: found

info = 0
n = fibreCount(root, f)
total = int(n, ill)*nodeCount(root, f)
nv = int(nodeCount(root, f))
allocate(v(3,nv), layer(nv))

if (root.eq.'D6') then
  call ringNodes(f, v, layer)
  allocate(poles(3,2))
  poles(1:3,1) = (/ 0.D0, 0.D0, 1.D0 /)
  poles(1:3,2) = (/ 0.D0, 0.D0, -1.D0 /)
  zorder = 6
else
  call polyhedronNodes(root, f, v, layer)
  call polyhedron(root, poles, faces)
  zorder = merge(4, 5, root.eq.'O ')
end if
sectorAngle = 2.D0*cPi/dble(zorder)

call groupElements(root, el)
ng = size(el,2)
total = total/ng
allocate(mats(3,3,ng), images(3,ng), q(4,total))
do g=1,ng
  mats(1:3,1:3,g) = transpose(qmat(el(1:4,g)))
end do

cnt = 0
do i=1,nv

! the representative of an orbit is its node in the sector at z, the cell of the vertex z among the
! polyhedron's vertices (for D6 the northern hemisphere) at azimuths 0 .. 2 pi / zorder; on the edges of
! the sector it is the one with the largest x, then y
  found = .FALSE.
  kx = 0.D0
  ky = 0.D0
  do g=1,ng
    images(1:3,g) = matmul(mats(1:3,1:3,g), v(1:3,i))
    zmax = -2.D0
    do ip=1,size(poles,2)
      zmax = max(zmax, dot_product(images(1:3,g), poles(1:3,ip)))
    end do
    if (images(3,g).lt.zmax-1.D-11) cycle
    az = atan2(images(2,g), images(1,g))
    if (.not.(((az.ge.-1.D-11).and.(az.le.sectorAngle+1.D-11)).or. &
              (sqrt(images(1,g)**2+images(2,g)**2).lt.1.D-11))) cycle
    bx = anint(images(1,g)*1.D11)
    by = anint(images(2,g)*1.D11)
    if ((.not.found).or.(bx.gt.kx).or.((bx.eq.kx).and.(by.gt.ky))) then
      found = .TRUE.
      kx = bx
      ky = by
      w = images(1:3,g)
    end if
  end do
  if (.not.found) then
    info = 2
    return
  end if
  if (sum(abs(w-v(1:3,i))).gt.1.D-9) cycle

! the stabilizer order of the node
  k = 0
  do g=1,ng
    if (sum(abs(images(1:3,g)-w)).lt.1.D-9) k = k + 1
  end do
  if (mod(n,k).ne.0) then
    info = 3
    return
  end if

! a rotation s taking v to z; a fixed 180 degree turn about x at -z
  s = (/ 1.D0 + w(3), w(2), -w(1), 0.D0 /)
  sn = sqrt(sum(s*s))
  if (sn.lt.1.D-12) then
    s = (/ 0.D0, 1.D0, 0.D0, 0.D0 /)
  else
    s = s/sn
  end if

! the rotations on the fibre, with odd layers staggered by half a step
  phase = cPi*dble(mod(layer(i),2))
  do j=0,n/k-1
    t = (2.D0*cPi*dble(j) + phase)/dble(n)
    rz = (/ cos(0.5D0*t), 0.D0, 0.D0, sin(0.5D0*t) /)
    cnt = cnt + 1
    if (cnt.gt.total) then
      info = 2
      return
    end if
    q(1:4,cnt) = qmul(rz, s)
    q(1:4,cnt) = q(1:4,cnt)/sqrt(sum(q(1:4,cnt)**2))
  end do
end do

if (int(cnt,ill).ne.total) info = 2

end subroutine rootGrid

!--------------------------------------------------------------------------
subroutine productGrid(FZtype, FZorder, FZordersign, f, qFZ, info)
!DEC$ ATTRIBUTES DLLEXPORT :: productGrid
  !! author: ZTV
  !! version: 1.0
  !! date: 11/02/26
  !!
  !! the product grid at frequency f for a fundamental zone given by its so3_T type and order
  !!
  !! returns EMsoft unit quaternions (4,N_FZ), reduced to the fundamental zone with scalar part >= 0.
  !! info = 0 on success, 1 for a zone without a product grid (822, 1022, 1222) or f < 1, 2 for an
  !! internal inconsistency, and 3 if the fibre count is not divisible by a stabilizer order.

IMPLICIT NONE

integer(kind=irg), INTENT(IN)             :: FZtype, FZorder, FZordersign, f
real(kind=dbl), allocatable, INTENT(OUT)  :: qFZ(:,:)
integer(kind=irg), INTENT(OUT)            :: info

character(2)                              :: group, root
real(kind=dbl)                            :: embed(4), setting(4), p(4), best(4), c(4), s(4)
real(kind=dbl), allocatable               :: q(:,:), parent(:,:), sub(:,:), reps(:,:), K(:,:)
logical, allocatable                      :: seen(:)
integer(kind=irg)                         :: i, j, g, a, nc
integer(kind=ill)                         :: total, cnt

call groupSetting(FZtype, FZorder, FZordersign, group, root, embed, setting, info)
if ((info.ne.0).or.(f.lt.1)) then
  info = 1
  return
end if

call rootGrid(root, f, q, info)
if (info.ne.0) return

! the subgroup R H R* inside the root group, and representatives of its left cosets
call groupElements(root, parent)
call groupElements(group, sub)
do i=1,size(sub,2)
  sub(1:4,i) = qmul(qmul(embed, sub(1:4,i)), qconj(embed))
end do
allocate(seen(size(parent,2)), reps(4,size(parent,2)))
seen = .FALSE.
nc = 0
do i=1,size(parent,2)
  if (seen(i)) cycle
  nc = nc + 1
  reps(1:4,nc) = parent(1:4,i)
  do j=1,size(sub,2)
    c = qmul(parent(1:4,i), sub(1:4,j))
    do a=1,size(parent,2)
      if (abs(dot_product(parent(1:4,a), c)).gt.1.D0-1.D-9) seen(a) = .TRUE.
    end do
  end do
end do
if (nc*size(sub,2).ne.size(parent,2)) then
  info = 2
  return
end if

! the symmetry operators of the EMsoft fundamental zone
call groupElements(group, K)
do i=1,size(K,2)
  K(1:4,i) = qmul(qmul(setting, K(1:4,i)), qconj(setting))
end do

total = int(size(q,2),ill)*nc
allocate(qFZ(4,total))
cnt = 0_ill
do i=1,size(q,2)
  do j=1,nc

! the active orientation q c R goes over to the EMsoft quaternion in the setting of the zone
    p = qmul(setting, qconj(qmul(qmul(q(1:4,i), reps(1:4,j)), embed)))

! its symmetry equivalent g p closest to the identity is in the zone
    best = p*merge(1.D0, -1.D0, p(1).ge.0.D0)
    do g=1,size(K,2)
      s = qmul(K(1:4,g), p)
      if (abs(s(1)).gt.best(1)) best = s*merge(1.D0, -1.D0, s(1).ge.0.D0)
    end do

! rounding must not leave a residue in a component that vanishes, as for a 180 degree rotation in a cyclic zone
    where (abs(best).lt.1.D-12) best = 0.D0
    cnt = cnt + 1_ill
    qFZ(1:4,cnt) = best/sqrt(sum(best*best))
  end do
end do

end subroutine productGrid

end module mod_productgrid
