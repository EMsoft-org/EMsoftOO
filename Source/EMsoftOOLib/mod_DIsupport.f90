! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_DIsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/31/20
  !!
  !! EMDIsupport contains a series of routines that are used by the dictionary indexing programs
  !! Some of these routines used to be in the commonmod module.

use mod_kinds
use mod_global
use mod_HDFsupport 
use HDF5
use stringconstants

IMPLICIT NONE 

public :: DIgetAverageOrientations, getOrientationSimilarityMap, getIndexingSuccessMap, &
          getKAMMap, getEMsoftPCcoordinates

!DEC$ ATTRIBUTES DLLEXPORT :: DIgetAverageOrientations
!DEC$ ATTRIBUTES DLLEXPORT :: getOrientationSimilarityMap
!DEC$ ATTRIBUTES DLLEXPORT :: getIndexingSuccessMap
!DEC$ ATTRIBUTES DLLEXPORT :: getKAMMap
!DEC$ ATTRIBUTES DLLEXPORT :: getEMsoftPCcoordinates
  
contains

!--------------------------------------------------------------------------
recursive subroutine DIgetAverageOrientations(ipar, Eulers, tmi, dplist, avEuler, disorient)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! Use the top near-matches list to compute the averaged orientations

use mod_rotations
use mod_quaternions
use mod_so3

IMPLICIT NONE

! ipar(1) = pgnum
! ipar(2) = FZcnt
! ipar(3) = Nexp
! ipar(4) = nnk
! ipar(5) = indi (=Ne*ceiling(float(totnumexpt)/float(Ne)))
! ipar(6) = nmuse

integer(kind=irg),INTENT(IN)            :: ipar(6)
real(kind=sgl),INTENT(IN)               :: Eulers(3,ipar(2))
integer(kind=irg),INTENT(IN)            :: tmi(ipar(4),ipar(5))
real(kind=sgl),INTENT(IN)               :: dplist(ipar(4),ipar(5))
real(kind=sgl),INTENT(OUT)              :: avEuler(3,ipar(3))
real(kind=sgl),INTENT(OUT),OPTIONAL     :: disorient(ipar(3),ipar(6))

type(QuaternionArray_T)                 :: Pm, dummy
type(q_T)                               :: q1, q2
type(a_T)                               :: ax
type(e_T)                               :: eu
type(so3_T)                             :: SO

integer(kind=irg)                       :: Nexp, nnm, nnk, Pmdims, i, j, k, ipat, pgnum, FZcnt, indi, nmuse
real(kind=dbl)                          :: aa(4), accum(3), qus(4), a, p(4), qsmall(4), theta, vec(3)

real(kind=dbl),allocatable              :: logq(:,:), w(:)
integer(kind=irg),allocatable           :: EulerIDs(:)
logical                                 :: store

pgnum = ipar(1)
FZcnt = ipar(2)
Nexp = ipar(3)
nnk = ipar(4)
indi = ipar(5)
nmuse = ipar(6)

store = .FALSE.
if (PRESENT(disorient)) store = .TRUE.

! get the symmetry operator quaternions for the point group 
call dummy%QSym_Init(pgnum, Pm)
Pmdims = Pm%getQnumber()

!===================================
! ok, so now we have all the necessary data
! next, we convert all the dictionary Euler angles into axis-angle triplets 
! but with half the angle so that they become the logarithm of the corresponding 
! quaternions (we omit the scalar part which is always zero)
allocate(logq(3,FZcnt))
do i=1,FZcnt
  eu = e_T( edinp = dble(Eulers(1:3,i)) )
  ax = eu%ea()
  aa = ax%a_copyd()
  logq(1:3,i) = aa(1:3) * aa(4) * 0.5
end do

! next, loop over all the experimental points, determine the weight factors
! and perform the averaging over those symmetrically equivalent orientations
! that have the smallest possible misorientation with respect to the best match.
allocate(EulerIDs(nmuse),w(nmuse))

do i=1,Nexp
! get the Euler angle IDs in the short list
  EulerIDs(1:nmuse) = tmi(1:nmuse,i)

! get the weight factors and normalize them to sum to 1
  w(1:nmuse) = dplist(1:nmuse,i)
  w = w - w(nmuse)
  w = w/sum(w)

! for each orientation in the list, determine the symmetrycally equivalent one
! that has the smallest misorientation with respect to the first entry on the list
! and add its weighted logarithm to the accum sum
  eu = e_T( edinp = dble(Eulers(1:3,EulerIDs(1))) )
  q1 = eu%eq()
  accum = logq(1:3,EulerIDs(1)) * w(1)
  do j=2,nmuse!  -1   ! -1 because the last one has weight factor 0
    eu = e_T( edinp = dble(Eulers(1:3,EulerIDs(j))) )
    q2 = eu%eq()
! determine the orientation with the smallest misorientation w.r.t. q1 and store it in qsmall
    call SO%getDisorientation(Pm, q1, q2, ax, fix1=.TRUE.)
    p = ax%a_copyd()
    if (store.eqv..TRUE.) disorient(i,j) = p(4)
! add p with the appropriate weight factor to accum
    accum(1:3) = accum(1:3) + p(1:3) * p(4) * 0.5 * w(j)
  end do
  
! accum is now the logarithm of the desired orientation quaternion, so we need to convert
! this back to an Euler angle triplet
  theta = sqrt(sum(accum**2))
  theta = mod(theta, 2.0*sngl(cPi))
  if (theta.ne.0.0) then
    vec = accum/theta
  else 
    vec = (/ 0.0, 0.0, 1.0 /)
  end if
  p(1:3) = vec(1:3)
  p(4) = theta  * 2.D0
  ax = a_T( adinp = p )
  eu = ax%ae()
  avEuler(1:3,i) = eu%e_copyd()
end do

! and put this array back in degrees
avEuler = avEuler * rtod

end subroutine DIgetAverageOrientations

!--------------------------------------------------------------------------
recursive subroutine getOrientationSimilarityMap(idims, tmi, nm, ipf_wd, ipf_ht, osm)
!! author: MDG 
!! version: 1.0 
!! date: 04/01/20
!!
!! compute the OSM (Orientation Similarity Map) given a set of near matches

use mod_math
use mod_io

IMPLICIT NONE

integer(kind=irg),INTENT(IN)     :: idims(2)
integer(kind=irg),INTENT(IN)     :: tmi(idims(1),idims(2))
integer(kind=irg),INTENT(IN)     :: nm
integer(kind=irg),INTENT(IN)     :: ipf_wd
integer(kind=irg),INTENT(IN)     :: ipf_ht
real(kind=sgl),INTENT(OUT)       :: osm(ipf_wd,ipf_ht)

type(IO_T)                       :: Message
real(kind=sgl),allocatable       :: localosm(:)
integer(kind=irg),allocatable    :: lstore(:,:), pstore(:,:), cp(:), lp(:)
integer(kind=irg)                :: ii, jj, iii , dp, lnm, io_int(2), totnumexpt

osm = 0.0
totnumexpt = ipf_wd * ipf_ht

! make sure that the requested number of near-matches is smaller than/equal to the available number
if (nm.gt.idims(1)) then
  io_int(1) = nm
  io_int(2) = idims(1)
  call Message%WriteValue('Requested number of near matches is too large: ',io_int,2,"(I4,' > ',I4)")
  call Message%printMessage(' --> Resetting requested number to maximum available')
  lnm = idims(1)
else
  lnm = nm
end if

allocate(lstore(lnm, ipf_wd), pstore(lnm, ipf_wd), cp(lnm), lp(lnm), localosm(ipf_wd*ipf_ht))

localosm = 0.0
lstore = 0
pstore = 0
cp = 0
lp = 0

! we'll do this computation on the 1D array, in the same way 
! as the ADP (Average Dot Product) map in the EBSDDI.f90 program
do iii = 1,totnumexpt
    ii = mod(iii,ipf_wd)
    if (ii.eq.0) ii = ipf_wd
    jj = iii/ipf_wd+1
! do we need to copy pstore into lstore ?
    if ((ii.eq.1).and.(jj.gt.1)) lstore = pstore
! determine to which osm entries we need to add the similarity count
    if (ii.eq.1) then
      cp(1:lnm) = tmi(1:lnm, iii)
      pstore(1:lnm,ii) = cp(1:lnm)
    else
      lp = cp
      cp(1:lnm) = tmi(1:lnm, iii)
      pstore(1:lnm,ii) = cp(1:lnm)
      dp = vectormatch(lnm, cp, lp)
      localosm(iii-1) = localosm(iii-1) + dp
      localosm(iii) = localosm(iii) + dp
    end if
    if (jj.gt.1) then
      dp = vectormatch(lnm,lstore(1:lnm,ii),cp)
      localosm(iii-ipf_wd+1) = localosm(iii-ipf_wd+1) + dp
      localosm(iii) = localosm(iii) + dp
    end if
end do

! correct the osm values depending on inside, edge, or corner pixels
! divide by 4
localosm = localosm*0.25

! correct the straight segments 
localosm(2:ipf_wd-1) = localosm(2:ipf_wd-1) * 4.0/3.0
localosm(totnumexpt-ipf_wd+2:totnumexpt-1) = localosm(totnumexpt-ipf_wd+2:totnumexpt-1) * 4.0/3.0
do jj=1,ipf_ht-2
  localosm(ipf_wd*jj+1) = localosm(ipf_wd*jj+1) * 4.0/3.0
end do
do jj=2,ipf_ht-1
  localosm(ipf_wd*jj) = localosm(ipf_wd*jj) * 4.0/3.0
end do

! and the corners
localosm(1) = localosm(1) * 4.0
localosm(ipf_wd) = localosm(ipf_wd) * 2.0
localosm(totnumexpt) = localosm(totnumexpt) * 2.0
localosm(totnumexpt-ipf_wd+1) = localosm(totnumexpt-ipf_wd+1) * 4.0/3.0

! and we deallocate the auxiliary variables 
deallocate(lstore,pstore,lp,cp)

do ii=1,ipf_wd
  do jj=1,ipf_ht
    osm(ii,jj) = localosm(ipf_wd*(jj-1)+ii)
  end do
end do

end subroutine getOrientationSimilarityMap

!--------------------------------------------------------------------------
recursive subroutine getIndexingSuccessMap(ipar, tmi, ea, nism, nnk, nt, ism)
!! author: MDG 
!! version: 1.0 
!! date: 04/01/20
!!
!! compute the ISM (Indexing Success Map) given a set of near matches

use omp_lib
use mod_io
use mod_quaternions 
use mod_rotations
use mod_so3

IMPLICIT NONE

integer(kind=irg),INTENT(IN)                      :: ipar(10)
integer(kind=irg),INTENT(IN)                      :: tmi(ipar(1),ipar(2))
real(kind=sgl),INTENT(INOUT)                      :: ea(3,ipar(4))
!f2py intent(in,out) ::  ea
integer(kind=irg),INTENT(IN)                      :: nism 
integer(kind=irg),INTENT(IN)                      :: nnk
integer(kind=irg),INTENT(IN)                      :: nt
real(kind=sgl),INTENT(OUT)                        :: ism(ipar(7)*ipar(8))
       
type(IO_T)                                        :: Message
type(a_T)                                         :: disax
type(e_T)                                         :: ea1, ea2
type(QuaternionArray_T)                           :: Pm, dummy
type(so3_T)                                       :: SO 

integer(kind=irg)                                 :: io_int(2), lnism, i, j
real(kind=dbl),allocatable                        :: angles(:)
real(kind=dbl)                                    :: angle, a(4)

ism = 0.0

! make sure that the requested number of near-matches is smaller than/equal to the available number
if (nism.gt.nnk-1) then
  io_int(1) = nism
  io_int(2) = nnk-1
  call Message%WriteValue('Requested number of near matches is too large: ',io_int,2,"(I4,' > ',I4)")
  call Message%printMessage(' --> Resetting requested number to maximum available')
  lnism = nnk-1
else
  lnism = nism
end if

! set up the correct symmetry variables 
call dummy%QSym_Init(ipar(6), Pm)

! next we go through the entire list of points in tmi and compute the misorientation angle
! for the best match with respect to the next nism matches

! this should be done in parallel ... 
call OMP_SET_NUM_THREADS(nt)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,angles,angle)
allocate(angles(lnism))
!$OMP DO SCHEDULE(DYNAMIC)
do i=1,ipar(3)
  ea1 = e_T( edinp = dble(ea(1:3,tmi(1,i))) )
  do j=2,lnism+1
    ea2 = e_T( edinp = dble(ea(1:3,tmi(j,i))) )
    call SO%getDisorientation(Pm, ea1, ea2, disax)
    a = disax%a_copyd()
    angles(j-1) = a(4)
  end do
  ism(i) = minval(angles)
end do
!$OMP END DO
deallocate(angles)
!$OMP END PARALLEL

ism = ism * rtod

! that's it.

end subroutine getIndexingSuccessMap

!--------------------------------------------------------------------------
recursive subroutine getKAMMap(numeu, eulers, ipf_wd, ipf_ht, pgnum, kam)
!! author: MDG 
!! version: 1.0 
!! date: 04/01/20
!!
!! compute the KAM (Kernel Average Misorientation Map) given a set of orientations

use mod_math
use mod_quaternions
use mod_rotations
use mod_so3

IMPLICIT NONE

integer(kind=irg),INTENT(IN)      :: numeu
real(kind=sgl),INTENT(IN)         :: eulers(3,numeu)
integer(kind=irg),INTENT(IN)      :: ipf_wd
integer(kind=irg),INTENT(IN)      :: ipf_ht
integer(kind=irg),INTENT(IN)      :: pgnum
real(kind=sgl),INTENT(OUT)        :: kam(ipf_wd,ipf_ht)

type(a_T)                         :: disax
type(e_T)                         :: ea1, ea2
type(QuaternionArray_T)           :: Pm, dummy
type(so3_T)                       :: SO 

real(kind=dbl),allocatable        :: localkam(:)
real(kind=dbl),allocatable        :: lstore(:,:), pstore(:,:)
real(kind=dbl)                    :: cp(3), lp(3)
integer(kind=irg)                 :: ii, jj, iii 
real(kind=dbl)                    :: dp, a(4)

kam = 0.0

allocate(lstore(3,ipf_wd), pstore(3,ipf_wd), localkam(numeu))

localkam = 0.0
lstore = 0
pstore = 0
cp = 0
lp = 0

! set up the correct symmetry variables 
call dummy%QSym_Init(pgnum, Pm)

! we'll do this computation on the 1D array, in the same way 
! as the ADP (Average Dot Product) routine.
do iii = 1,numeu
    ii = mod(iii,ipf_wd)
    if (ii.eq.0) ii = ipf_wd
    jj = iii/ipf_wd+1
! do we need to copy pstore into lstore ?
    if ((ii.eq.1).and.(jj.gt.1)) lstore = pstore
! determine to which kam entries we need to add the next disorientation value
    if (ii.eq.1) then
      cp = eulers(1:3,iii)
      pstore(1:3,ii) = cp
    else
      lp = cp
      cp = eulers(1:3,iii)
      pstore(1:3,ii) = cp
      ea1 = e_T( edinp = dble(lp) )
      ea2 = e_T( edinp = dble(cp) )
      call SO%getDisorientation(Pm, ea1, ea2, disax)
      a = disax%a_copyd()
      localkam(iii-1) = localkam(iii-1) + a(4)
      localkam(iii) = localkam(iii) + a(4)
    end if
    if (jj.gt.1) then
      ea1 = e_T( edinp = dble(lstore(1:3,ii)) )
      ea2 = e_T( edinp = dble(cp) )
      call SO%getDisorientation(Pm, ea1, ea2, disax)
      a = disax%a_copyd()
      localkam(iii-ipf_wd+1) = localkam(iii-ipf_wd+1) + a(4)
      localkam(iii) = localkam(iii) + a(4)
    end if
end do

! correct the kam values depending on inside, edge, or corner pixels
! divide by 4
localkam = localkam*0.25

! correct the straight segments  
localkam(2:ipf_wd-1) = localkam(2:ipf_wd-1) * 4.0/3.0
localkam(numeu-ipf_wd+2:numeu-1) = localkam(numeu-ipf_wd+2:numeu-1) * 4.0/3.0
do jj=1,ipf_ht-2
  localkam(ipf_wd*jj+1) = localkam(ipf_wd*jj+1) * 4.0/3.0
end do
do jj=2,ipf_ht-1
  localkam(ipf_wd*jj) = localkam(ipf_wd*jj) * 4.0/3.0
end do

! and the corners
localkam(1) = localkam(1) * 4.0
localkam(ipf_wd) = localkam(ipf_wd) * 2.0
localkam(numeu) = localkam(numeu) * 2.0
localkam(numeu-ipf_wd+1) = localkam(numeu-ipf_wd+1) * 4.0/3.0

! and we deallocate the auxiliary variables 
deallocate(lstore,pstore)

do ii=1,ipf_wd
  do jj=1,ipf_ht
    kam(ii,jj) = sngl(localkam(ipf_wd*(jj-1)+ii))
  end do
end do

end subroutine getKAMMap


!--------------------------------------------------------------------------
!
! subroutine: getEMsoftPCcoordinates
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert pattern center coordinates to EMsoft units for each vendor
!
!> @param pctr array of 3 PC coordinates
!> @param vendor vendor string
!> @param delta pixel size [microns]
!> @param Nx number of detector pixels along x
!> @param Ny number of detector pixels along y
!
!> @date 08/20/19 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function getEMsoftPCcoordinates(pctr, vendor, delta, Nx, Ny) result(EMsoftpc)
!! author: MDG 
!! version: 1.0 
!! date: 04/01/20
!!
!! convert pattern center coordinates to EMsoft units for each vendor

use mod_io

IMPLICIT NONE 

real(kind=sgl),INTENT(IN)           :: pctr(3) 
character(fnlen),INTENT(IN)         :: vendor 
real(kind=sgl),INTENT(IN)           :: delta
integer(kind=irg),INTENT(IN)        :: Nx
integer(kind=irg),INTENT(IN)        :: Ny
real(kind=sgl)                      :: EMsoftpc(3)

type(IO_T)                          :: Message
real(kind=sgl)                      :: io_real(3)

if (trim(vendor).eq.'EMsoft') then 
  EMsoftpc = pctr 
end if 

if (trim(vendor).eq.'EDAX/TSL') then 
  EMsoftpc = (/ Nx * (pctr(1) - 0.5), Nx * pctr(2) - Ny*0.5, Nx * delta * pctr(3) /)
end if 

if (trim(vendor).eq.'Oxford') then 
  EMsoftpc = (/ Nx * (pctr(1) - 0.5), Ny * (pctr(2) - 0.5), Nx * delta * pctr(3) /)
end if 

if (trim(vendor).eq.'Bruker') then 
  EMsoftpc = (/ Nx * (pctr(1) - 0.5), Ny * (0.5 - pctr(2)), Nx * delta * pctr(3) /)
end if 

if (trim(vendor).ne.'EMsoft') then 
  io_real = pctr 
  call Message%WriteValue('Input pattern center coordinates in '//trim(vendor)//' convention : ',io_real,3)
  io_real = EMsoftpc
  call Message%WriteValue('  --> pattern center coordinates in EMsoft convention : ',io_real,3)
end if

end function getEMsoftPCcoordinates



end module mod_DIsupport