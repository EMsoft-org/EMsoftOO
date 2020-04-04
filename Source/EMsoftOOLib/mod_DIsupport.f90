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

public :: h5_writeFile, DIgetAverageOrientations, getOrientationSimilarityMap, getIndexingSuccessMap, &
          getKAMMap, getEMsoftPCcoordinates

private :: h5_writeInfo, h5_write2DImageFromVector, h5_writeCoordinateSystemGroup, h5_writePatternCenterGroup, &
           h5_writePhaseGroup

!DEC$ ATTRIBUTES DLLEXPORT :: h5_writeFile
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
recursive subroutine getIndexingSuccessMap(ipar, tmi, ea, nml, ism)
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
use mod_DIfiles

IMPLICIT NONE

integer(kind=irg),INTENT(IN)                      :: ipar(10)
integer(kind=irg),INTENT(IN)                      :: tmi(ipar(1),ipar(2))
real(kind=sgl),INTENT(INOUT)                      :: ea(3,ipar(4))
!f2py intent(in,out) ::  ea
class(DictionaryIndexingNameListType),INTENT(IN)  :: nml
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
if (nml%nism.gt.nml%nnk-1) then
  io_int(1) = nml%nism
  io_int(2) = nml%nnk-1
  call Message%WriteValue('Requested number of near matches is too large: ',io_int,2,"(I4,' > ',I4)")
  call Message%printMessage(' --> Resetting requested number to maximum available')
  lnism = nml%nnk-1
else
  lnism = nml%nism
end if

! set up the correct symmetry variables 
call dummy%QSym_Init(ipar(6), Pm)

! next we go through the entire list of points in tmi and compute the misorientation angle
! for the best match with respect to the next nism matches

! this should be done in parallel ... 
call OMP_SET_NUM_THREADS(nml%nthreads)
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

ism = ism * 180.0/sngl(cPi)

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


!--------------------------------------------------------------------------
subroutine h5_writeFile(EMsoft, HDF, HDFnames, vendor, ebsdnl, mcnl, xtalname, dstr, tstrb, ipar, resultmain, exptIQ, &
                         indexmain, dicteulerarray, dpmap, progname, nmldeffile, OSMmap)
!! author: MDG 
!! version: 1.0 
!! date: 04/03/20
!!
!! perform the computations

use mod_EMsoft
use mod_HDFnames
use mod_io
use mod_timing
use mod_MCfiles
use mod_DIfiles 

IMPLICIT NONE 

type(EMsoft_T), INTENT(INOUT)                       :: EMsoft
type(HDF_T), INTENT(INOUT)                          :: HDF
type(HDFnames_T), INTENT(INOUT)                     :: HDFnames
character(3),INTENT(IN)                             :: vendor   ! 'TSL' 'HKL' 'BRU'
type(DictionaryIndexingNameListType),INTENT(INOUT)  :: ebsdnl
!f2py intent(in,out) ::  ebsdnl
type(MCOpenCLNameListType),INTENT(INOUT)            :: mcnl
character(fnlen),INTENT(IN)                         :: xtalname
character(11),INTENT(INOUT)                         :: dstr
!f2py intent(in,out) ::  dstr
character(15),INTENT(IN)                            :: tstrb
integer(kind=irg),INTENT(INOUT)                     :: ipar(10)
!f2py intent(in,out) ::  ipar
real(kind=sgl),INTENT(IN)                           :: resultmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(IN)                           :: exptIQ(ipar(3))
integer(kind=irg),INTENT(IN)                        :: indexmain(ipar(1),ipar(2))
real(kind=sgl),INTENT(INOUT)                        :: dicteulerarray(3,ipar(4))
!f2py intent(in,out) ::  dicteulerarray
real(kind=sgl),INTENT(IN)                           :: dpmap(ipar(3))
character(fnlen),INTENT(IN)                         :: progname
character(fnlen),INTENT(IN)                         :: nmldeffile
real(kind=sgl),INTENT(OUT)                          :: OSMmap(ipar(7),ipar(8))

type(IO_T)                                          :: Message 
type(timing_T)                                      :: timer

character(15)                                       :: tstre
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
integer(kind=irg)                                   :: hdferr, filetype, i, j, ii, jj,indx, istat, ipar2(6), L
character(fnlen)                                    :: groupname, dataset, h5ebsdfile, savefile
logical                                             :: noindex, g_exists, overwrite=.TRUE.
real(kind=sgl)                                      :: eulerarray(3,ipar(4)), WD

real(kind=sgl),allocatable                          :: kam(:,:), ISMap(:)

real(kind=sgl),allocatable                          :: exptCI(:), eangle(:), eangles(:,:), results(:), avEuler(:,:), &
                                                       lresultmain(:,:), eulers(:,:) 
integer(kind=1),allocatable                         :: iPhase(:), valid(:)
integer(kind=irg),allocatable                       :: SEMsignal(:), lindexmain(:,:)
real(kind=sgl)                                      :: isratio, io_real(1)

! copy the dictionary euler angle array 
eulerarray = dicteulerarray

!=====================================================
! write the output in the format of an h5ebsd file
!!!! THIS PART IS STILL UNDER DEVELOPMENT !!!!
! we use the TSL h5ebsd file as a template for now; this 
! can be extended later other vendor formats
!=====================================================

if (vendor.ne.'TSL') then
  call Message%printMessage(' Only TSL h5ebsd file format is implemented in this version.')
  call Message%printMessage(' Program results will be saved in this format.')
end if

allocate(stringarray(1))

call timer%makeTimeStamp()
tstre = timer%getTimeString()

! Create a new file using the default properties.
  hdferr =  HDF%createFile(ebsdnl%datafile)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening file', hdferr)
  filetype = 1

  call h5_writeInfo(EMsoft, HDF, HDFnames, filetype, dstr, tstrb, tstre, progname, ebsdnl, nmldeffile)

! here we start with the h5ebsd-specific stuff
  groupname = 'Scan 1'
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening group Scan 1', hdferr)

groupname = SC_EBSD
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening group EBSD', hdferr)

!=====================================================
!=====================================================
groupname = SC_Data
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error opening group Data', hdferr)

! there are 15 datasets in this Data group: CI, Fit, IQ, PRIAS Bottom Strip,
! PRIAS Center Square, PRIAS Top Strip, Pattern, Phase, Phi, Phi1, Phi2,
! SEM Signal, Valid, X Position, Y Position

!=====================================================
! CI Confidence Index: real(kind=sgl), one for each pattern... we take this
! to be the largest dot product
dataset = SC_CI
  allocate(exptCI(ipar(3)))
  exptCI = resultmain(1,1:ipar(3))
  hdferr = HDF%writeDatasetFloatArray(dataset, exptCI, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset CI', hdferr)

! we also insert a visual map of the Confidence Index, resampled on a rectangular array
dataset = SC_CIMap
  call h5_write2DImageFromVector(HDF, dataset, exptCI, ipar(3), ebsdnl)
  deallocate(exptCI)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset CImap', hdferr)

!=====================================================
! Fit 
dataset = SC_Fit
  allocate(eangle(ipar(3)),stat=istat)
  eangle = 1.0
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Fit', hdferr)
  deallocate(eangle)

!=====================================================
! this option is disabled starting in version 4.3 [05/19/19, MDG]; code block can be deleted
! Averaged Orientation Map (using near-match list and quaternion logarithm averaging)
! define the ipar2 entries
!   allocate(avEuler(3,ipar(3)))
  ipar2(1) = ipar(6)
  ipar2(2) = ipar(5)
  ipar2(3) = ipar(3)
  ipar2(4) = ipar(1)
  ipar2(5) = ipar(2)
  ipar2(6) = ebsdnl%nnav
! ! get the avEuler array
!   eulerarray = eulerarray * sngl(cPi)/180.0
!   call EBSDgetAverageOrientations(ipar2, eulerarray, indexmain(1:ipar2(4),1:ipar2(5)), resultmain(1:ipar2(4),1:ipar2(5)), &
!                                   avEuler)
!   eulerarray = eulerarray * 180.0/sngl(cPi)

! ! and write it to the HDF file
! dataset = SC_AverageOrientations
!   hdferr = HDF_writeDatasetFloatArray2D(dataset, avEuler, 3, ipar(3))
!   if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset AverageOrientations')

! get the nearest neighbor Kernel Average Misorientation Map (KAM)
  allocate(eulers(3,ipar(3)))
  do i=1,ipar(3)
    eulers(1:3,i) = eulerarray(1:3,indexmain(1,i))
  end do
  eulers = eulers*dtor
dataset = SC_KAM
  if (sum(ebsdnl%ROI).ne.0) then
    allocate(kam(ebsdnl%ROI(3),ebsdnl%ROI(4)))
    call getKAMMap(ipar(3), eulers, ebsdnl%ROI(3), ebsdnl%ROI(4), ipar(6), kam)
    kam = kam*rtod
    hdferr = HDF%writeDatasetFloatArray(dataset, kam, ebsdnl%ROI(3), ebsdnl%ROI(4))
  else
    allocate(kam(ebsdnl%ipf_wd,ebsdnl%ipf_ht))
    call getKAMMap(ipar(3), eulers, ebsdnl%ipf_wd, ebsdnl%ipf_ht, ipar(6), kam)
    kam = kam*rtod
    hdferr = HDF%writeDatasetFloatArray(dataset, kam, ebsdnl%ipf_wd, ebsdnl%ipf_ht)
  end if
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset KAM', hdferr)
  deallocate(kam, eulers)

! get the Orientation Similarity Map (OSM); map is now returned to calling routine [MDG, 3/5/18]
dataset = SC_OSM
  if (sum(ebsdnl%ROI).ne.0) then
!   allocate(osm(ebsdnl%ROI(3),ebsdnl%ROI(4)))
    call getOrientationSimilarityMap( (/ipar(1), ipar(2)/), indexmain, ebsdnl%nosm, ebsdnl%ROI(3), ebsdnl%ROI(4), OSMmap)
    hdferr = HDF%writeDatasetFloatArray(dataset, OSMmap, ebsdnl%ROI(3), ebsdnl%ROI(4))
  else
!   allocate(osm(ebsdnl%ipf_wd,ebsdnl%ipf_ht))
    call getOrientationSimilarityMap( (/ipar(1), ipar(2)/), indexmain, ebsdnl%nosm, ebsdnl%ipf_wd, ebsdnl%ipf_ht, OSMmap)
    hdferr = HDF%writeDatasetFloatArray(dataset, OSMmap, ebsdnl%ipf_wd, ebsdnl%ipf_ht)
  end if
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset OSM', hdferr)

! we also insert a visual map of the Confidence Index, resampled on a rectangular array
! dataset = 'FitMap'
! call h5ebsd_write2DImageFromVector(dataset, totnumexpt, exptCI, ebsdnl)

!=====================================================
! IQ Image Quality; computed using the second moment of the pattern power spectrum
dataset = SC_IQ
  hdferr = HDF%writeDatasetFloatArray(dataset, exptIQ, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset IQ', hdferr)
  
! we also insert a visual map of the Image Quality, resampled on a rectangular array
dataset = SC_IQMap
  call h5_write2DImageFromVector(HDF, dataset, exptIQ, ipar(3), ebsdnl)

!=====================================================
! generate the indexing success map (ISM) 
  eulerarray = eulerarray * dtor
  allocate(ISMap(ipar(7)*ipar(8)))
  call getIndexingSuccessMap(ipar, indexmain, eulerarray, ebsdnl, ISMap)
  eulerarray = eulerarray * rtod

dataset = SC_ISM
  hdferr = HDF%writeDatasetFloatArray(dataset, ISMap, ipar(7)*ipar(8))

dataset = SC_ISMap
  call h5_write2DImageFromVector(HDF, dataset, ISMap, ipar(7)*ipar(8), ebsdnl, binary=ebsdnl%isangle)
  j = 0
  do i=1,ipar(7)*ipar(8)
    if (ISMap(i).le.ebsdnl%isangle) j = j+1
  end do 
  isratio = 100.0 * real(j) / real(ipar(7)*ipar(8))
  io_real(1) = isratio 
  call Message%WriteValue('Indexing Success Rate (%) : ',io_real,1)
  deallocate(ISMap)

dataset = SC_ISR 
  hdferr = HDF%writeDatasetFloat(dataset, isratio)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset ISR', hdferr)
 

!=====================================================
! PRIAS Bottom Strip: to be implemented
!   call Message('h5ebsd_writeFile: writing of ->PRIAS Bottom Strip<- data not yet implemented.')

!=====================================================
! PRIAS Center Square: to be implemented
!   call Message('h5ebsd_writeFile: writing of ->PRIAS Center Strip<- data not yet implemented.')

!=====================================================
! PRIAS Top Strip: to be implemented
!   call Message('h5ebsd_writeFile: writing of ->PRIAS Top Strip<- data not yet implemented.')

!=====================================================
! Pattern: in principle, this is where the fitted patterns could be stored
! This will require re-computing them for the best match orientations; we 
! could leave this as an option for the user, to be implemented.
!   call Message('h5ebsd_writeFile: writing of ->Pattern<- data not yet implemented.')

!=====================================================
! Phase: Phase identifier (all zero for now)
dataset = SC_Phase
  allocate(iPhase(ipar(3)),stat=istat)
  iPhase = 0
  hdferr = HDF%writeDatasetInteger1byteArray1D(dataset, iPhase, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phase', hdferr)
  deallocate(iPhase)

!=====================================================
! SEM Signal: all 0 for now
dataset = SC_SEMSignal
  allocate(SEMsignal(ipar(3)),stat=istat)
  SEMsignal = 10000
  hdferr = HDF%writeDatasetIntegerArray(dataset, SEMsignal, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset SEM Signal', hdferr)
  deallocate(SEMsignal)

!=====================================================
! Valid : all 0 for now
dataset = SC_Valid
  allocate(valid(ipar(3)),stat=istat)
  valid = 0
  hdferr = HDF%writeDatasetInteger1byteArray1D(dataset, valid, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Valid', hdferr)
  deallocate(valid)

dataset = SC_EulerAngles
  call H5Lexists_f(HDF%getObjectID(),trim(dataset),g_exists, hdferr)
  allocate(eangles(3,ipar(3)),stat=istat)
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangles(1:3,ii) = eulerarray(1:3,indx)
  end do
  eangles = eangles * dtor
  if (g_exists) then 
     hdferr = HDF%writeDatasetFloatArray(dataset, eangles, 3, ipar(3), overwrite)
  else
     hdferr = HDF%writeDatasetFloatArray(dataset, eangles, 3, ipar(3))
  end if
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset EulerAngles', hdferr)

!=====================================================
! Euler angles: Phi 
dataset = SC_Phi
  allocate(eangle(ipar(3)),stat=istat)
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangle(ii) = eulerarray(2,indx)
  end do
  eangle = eangle * dtor
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phi', hdferr)

!=====================================================
! Euler angles: Phi1
dataset = SC_Phi1
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangle(ii) = eulerarray(1,indx)
  end do
  eangle = eangle * dtor
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phi1', hdferr)

!=====================================================
! Euler angles: Phi2 
dataset = SC_Phi2
  do ii = 1,ipar(3)
    indx = indexmain(1,ii)
    eangle(ii) = eulerarray(3,indx)
  end do
  eangle = eangle * dtor
  hdferr = HDF%writeDatasetFloatArray(dataset, eangle, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Phi2', hdferr)
  deallocate(eangle)
 
!=====================================================
! X Position: list of x positions for sampling points; requires knowledge of step size
! from Header
dataset = SC_XPos
  allocate(results(ipar(3)),stat=istat)
  if (sum(ebsdnl%ROI).eq.0) then
    do jj=1,ebsdnl%ipf_ht
      do ii=1,ebsdnl%ipf_wd
        results(ebsdnl%ipf_wd*(jj-1)+ii) = (ii-1)*ebsdnl%StepX 
      end do
    end do
  else
    do jj=1,ebsdnl%ROI(4)
      do ii=1,ebsdnl%ROI(3)
        results(ebsdnl%ROI(3)*(jj-1)+ii) = (ii-1)*ebsdnl%StepX 
      end do
    end do
  end if
  hdferr = HDF%writeDatasetFloatArray(dataset, results, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset X Position', hdferr)

 
!=====================================================
! Y Position: list of y positions for sampling points; requires knowledge of step size
! from Header
dataset = SC_YPos
  if (sum(ebsdnl%ROI).eq.0) then
    do jj=1,ebsdnl%ipf_ht
      do ii=1,ebsdnl%ipf_wd
        results(ebsdnl%ipf_wd*(jj-1)+ii) = (ii-1)*ebsdnl%StepY 
      end do
    end do
  else
    do jj=1,ebsdnl%ROI(4)
      do ii=1,ebsdnl%ROI(3)
        results(ebsdnl%ROI(3)*(jj-1)+ii) = (ii-1)*ebsdnl%StepY 
      end do
    end do
  end if
  hdferr = HDF%writeDatasetFloatArray(dataset, results, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Y position', hdferr)
  deallocate(results)
 

!=====================================================
!=====================================================
! this concludes the standard data sets in the Data group
! here, we have additional data sets based on results from the 
! dictionary indexing program; these are not part of the standard
! TSL h5ebsd file format.
!=====================================================
! EBSD average dot product map 
dataset = SC_AvDotProductMap
  call h5_write2DImageFromVector(HDF, dataset, dpmap, ipar(3), ebsdnl)

! number of samples in dictionary
dataset = SC_FZcnt
  hdferr = HDF%writeDatasetInteger(dataset, ipar(5))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset FZcnt', hdferr)

! point group number
dataset = SC_PointGroupNumber
  hdferr = HDF%writeDatasetInteger(dataset, ipar(6))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset PointGroupNumber', hdferr)

! Ncubochoric
dataset = SC_Ncubochoric
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%ncubochoric)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Ncubochoric', hdferr)

! write the list of sampled Euler angles
dataset = SC_DictionaryEulerAngles
  hdferr = HDF%writeDatasetFloatArray(dataset, dicteulerarray, 3, ipar(4))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset EulerAngles', hdferr)

! number of experimental patterns 
dataset = SC_NumExptPatterns
  hdferr = HDF%writeDatasetInteger(dataset, ipar(3))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset NumExptPatterns', hdferr)

! list of top nnk dot product values
dataset = SC_TopDotProductList
  hdferr = HDF%writeDatasetFloatArray(dataset, resultmain, ipar(1), ipar(2))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset TopDotProductList', hdferr)

! indices of top matches into Euler angle list
dataset = SC_TopMatchIndices
  hdferr = HDF%writeDatasetIntegerArray(dataset, indexmain, ipar(1), ipar(2))
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset TopMatchIndices', hdferr)

! leave this group
  call HDF%pop()
!=====================================================
!=====================================================

!=====================================================
!=====================================================
! create the Header group
groupname = SC_Header
  hdferr = HDF%createGroup(groupname)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error creating group Header', hdferr)

! there are 15 datasets in this group: Camera Azimuth Angle, Camera Elevation Angle,
! Grid Type, Notes, Operator, Pattern Height, Pattern Width, Sample ID, Sample Tilt,
! Scan ID, Step X, Step Y, Working Distance, nColumns, nRows
! there are also 3 groups: Coordinate System, Pattern Center Calibration, and Phase

!=====================================================
! Camera Azimuthal Angle
dataset = SC_CameraAzimuthalAngle
  hdferr = HDF%writeDatasetFloat(dataset, 0.0)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Camera Azimuthal Angle', hdferr)
 
!=====================================================
! Camera Elevation Angle
dataset = SC_CameraElevationAngle
  hdferr = HDF%writeDatasetFloat(dataset, ebsdnl%thetac)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset TopMatchIndices', hdferr)
 
!=====================================================
! Coordinate System group
  call h5_writeCoordinateSystemGroup(EMsoft, HDF)

!=====================================================
! Grid Type
dataset = SC_GridType
  stringarray(1) = 'SqrGrid'
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Grid Type', hdferr)

!=====================================================
! Notes
dataset = SC_Notes
  stringarray(1) = ''
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Notes', hdferr)

!=====================================================
! Operator
dataset = SC_Operator
  stringarray(1) = trim(EMsoft%getConfigParameter('Username'))//' ['// &
                   trim(EMsoft%getConfigParameter('Useremail'))//']'
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Operator', hdferr)

!=====================================================
! Pattern Center Calibration group
  call h5_writePatternCenterGroup(HDF, ebsdnl%xpc, ebsdnl%ypc, ebsdnl%L, ebsdnl%delta, (/ebsdnl%numsx, ebsdnl%numsy/))

!=====================================================
! Pattern height
dataset = SC_PatternHeight
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%numsx/ebsdnl%binning)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Pattern Height', hdferr)

!=====================================================
! Pattern width
dataset = SC_PatternWidth
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%numsy/ebsdnl%binning)
  if (hdferr.ne.0) call HDF%error_check('h5_writeFile:Error writing dataset Pattern Width', hdferr)

!=====================================================
! Phase group
groupname = SC_Phase
  hdferr = HDF%createGroup(groupname)
groupname = "1"
  call h5_writePhaseGroup(EMsoft, HDF, groupname, xtalname)

! close the Phase group
  call HDF%pop()

!=====================================================
! Sample ID
dataset = SC_SampleID
  stringarray(1) = ''
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

!=====================================================
! Sample Tilt
dataset = SC_SampleTilt
  hdferr = HDF%writeDatasetFloat(dataset, sngl(mcnl%sig))
 
!=====================================================
! Scan ID
dataset = SC_ScanID
  stringarray(1) = ''
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

!=====================================================
! Step X
dataset = SC_StepX
  hdferr = HDF%writeDatasetFloat(dataset, ebsdnl%StepX)
 
!=====================================================
! Step Y
dataset = SC_StepY
  hdferr = HDF%writeDatasetFloat(dataset, ebsdnl%StepY)

!=====================================================
! Working Distance
dataset = SC_WorkingDistance
  WD = 0.0
  hdferr = HDF%writeDatasetFloat(dataset, WD)

!=====================================================
! nColumns
dataset = SC_nColumns
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%ipf_wd)

!=====================================================
! nRows
dataset = SC_nRows
  hdferr = HDF%writeDatasetInteger(dataset, ebsdnl%ipf_ht)

!=====================================================
!=====================================================

! once all these have been written, we simply pop all the way to the top and close the file
  call HDF%pop(.TRUE.)

end subroutine h5_writeFile

!--------------------------------------------------------------------------
subroutine h5_writeInfo(EMsoft, HDF, HDFnames, filetype, dstr, tstrb, tstre, progname, ebsdnl, nmldeffile)
!! author: MDG 
!! version: 1.0 
!! date: 04/03/20
!!
!! write general information fields to the h5ebsd file, including EMsoft specific fields

use mod_EMsoft 
use mod_HDFnames
use mod_DIfiles

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT)                        :: EMsoft
type(HDF_T),INTENT(INOUT)                           :: HDF
type(HDFnames_T),INTENT(INOUT)                      :: HDFnames
integer(kind=irg),INTENT(IN)                        :: filetype
character(11),INTENT(INOUT)                         :: dstr
!f2py intent(in,out) ::  dstr
character(15),INTENT(IN)                            :: tstrb
character(15),INTENT(IN)                            :: tstre
character(fnlen),INTENT(IN)                         :: progname
type(DictionaryIndexingNameListType),INTENT(INOUT)  :: ebsdnl
!f2py intent(in,out) ::  ebsdnl
character(fnlen),INTENT(IN)                         :: nmldeffile

type(DIfile_T)                                      :: DIFT
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                                    :: groupname, dataset, nmlname, manufacturer
integer(kind=irg)                                   :: hdferr

! to be replaced by HDFnames code
if (filetype.eq.1) then ! EBSDDictionarIndexing file
  manufacturer = 'EMDI.f90'
  nmlname = 'EBSDDictionaryIndexingNML'
else
  manufacturer = ''
  nmlname = ''
end if

allocate(stringarray(1))

! set the Manufacturer and Version data sets
dataset = SC_Manufacturer
  stringarray(1)= trim(manufacturer)
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

dataset = SC_Version
  stringarray(1)= 'EMsoft '//EMsoft%getConfigParameter('EMsoftversion')
  hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! add the EMsoft header group
! write the EMheader to the file
groupname = SC_h5EBSD
  call HDF%writeEMheader(dstr, tstrb, tstre, progname)

! create a namelist group to write all the namelist files into
groupname = SC_NMLfiles
  hdferr = HDF%createGroup(groupname)

! and write the nml file for this program to the HDF5 file
! read the text file and write the array to the file
  dataset = trim(nmlname)
  hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)

! leave this group
  call HDF%pop()
  
! create a namelist group to write all the namelist files into
groupname = SC_NMLparameters
  hdferr = HDF%createGroup(groupname)
  if (filetype.eq.1) then 
    call DIFT%writeHDFNameList(HDF, HDFnames, ebsdnl)
  end if

! leave this group
  call HDF%pop()

end subroutine h5_writeInfo

!--------------------------------------------------------------------------
subroutine h5_write2DImageFromVector(HDF, dataset, inpvec, nump, ebsdnl, binary)
!! author: MDG 
!! version: 1.0 
!! date: 04/03/20
!!
!! write a gray scale image to the HDF5 file starting from a 1D vector

use mod_io
use mod_DIfiles

IMPLICIT NONE

type(HDF_T),INTENT(INOUT)               :: HDF
character(fnlen),INTENT(IN)             :: dataset
integer(kind=irg),INTENT(IN)            :: nump
real(kind=sgl),INTENT(IN)               :: inpvec(nump)
type(DictionaryIndexingNameListType),INTENT(IN)  :: ebsdnl
real(kind=sgl),OPTIONAL,INTENT(IN)      :: binary

type(IO_T)                              :: Message
real(kind=sgl)                          :: mi, ma
integer(kind=irg)                       :: istat, ii, jj, hdferr
real(kind=sgl),allocatable              :: newvec(:)
integer(kind=irg),allocatable           :: image(:,:)
integer(HSIZE_T)                        :: width, height
logical                                 :: isbinary

isbinary = .FALSE.
if (present(binary)) isbinary=.TRUE.

allocate(newvec(nump),stat=istat)
if (istat.ne.0) call Message%printError('h5_write2DImageFromVector','Could not allocate array for copy of input image')

newvec = inpvec

if (sum(ebsdnl%ROI).ne.0) then
  width = ebsdnl%ROI(3)
  height = ebsdnl%ROI(4)
else
  width = ebsdnl%ipf_wd
  height = ebsdnl%ipf_ht
end if
allocate(image(width,height),stat=istat)
if (istat.ne.0) call Message%printError('h5_write2DImageFromVector','Could not allocate array for output image')

if (isbinary.eqv..TRUE.) then 
  do jj = 1,height
    do ii = 1, width
      if (newvec((jj-1)*width+ii).gt.ebsdnl%isangle) then
        image(ii,jj) = 0
      else
        image(ii,jj) = 255
      end if 
    end do 
  end do
else
  mi = minval(newvec)
  newvec = newvec - mi
  ma = maxval(newvec)

  do jj = 1,height
    image(1:width,jj) = int(255.0*newvec((jj-1)*width+1:jj*width)/ma)
  end do
end if 

call h5immake_image_8bit_f(HDF%getObjectID(),dataset,width,height,image,hdferr)
deallocate(image, newvec)

end subroutine h5_write2DImageFromVector

!--------------------------------------------------------------------------
subroutine h5_writeCoordinateSystemGroup(EMsoft, HDF)
!! author: MDG 
!! version: 1.0 
!! date: 02/17/20
!!
!! write information about the sample/detector coordinate frames 

use mod_EMsoft
use mod_io
use h5im

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT)       :: EMsoft
type(HDF_T),INTENT(INOUT)          :: HDF

type(IO_T)                         :: Message
character(fnlen)                   :: groupname, dataset, fname, resourcepathname
integer(kind=irg)                  :: hdferr
integer(HSIZE_T),allocatable       :: EBSDview(:,:,:), schematic(:,:,:)
character(1),allocatable           :: chararr(:,:,:)
integer(kind=irg)                  :: dims(3), istat
integer(HSIZE_T)                   :: width, height

! create the Coordinate System group
groupname = 'Coordinate System'
hdferr = HDF%createGroup(groupname)

!=====================================================
! EBSD View Reference Frame
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'EBSDview.data'
open(unit=50,file=trim(fname),status='old',form='unformatted')
read(50) dims
allocate(EBSDview(dims(1),dims(2),dims(3)),chararr(dims(1),dims(2),dims(3)),stat=istat)
if (istat.ne.0) call Message%printError('h5_writeCoordinateSystemGroup','Could not allocate array for EBSD view output image')
read(50) chararr
close(unit=50,status='keep')
EBSDview = ichar(chararr)

dataset = 'EBSD View Reference Frame'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(EBSDview),hdferr)
deallocate(EBSDview,chararr)

!=====================================================
! Schematic 1
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic1.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
allocate(schematic(dims(1),dims(2),dims(3)),chararr(dims(1),dims(2),dims(3)),stat=istat)
if (istat.ne.0) call Message%printError('h5_writeCoordinateSystemGroup','Could not allocate array for Schematic output image')
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 1'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

!=====================================================
! Schematic 2
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic2.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 2'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

!=====================================================
! Schematic 3
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic3.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 3'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

!=====================================================
! Schematic 4
fname = trim(EMsoft%generateFilePath('Resourcepathname'))//'Schematic4.data'
open(unit=50,file=fname,status='old',form='unformatted')
read(50) dims
read(50) chararr
close(unit=50,status='keep')
schematic = ichar(chararr)

dataset = 'Schematic 4'
width = dims(2)
height = dims(3)
call h5immake_image_24bit_f(HDF%getObjectID(),dataset,width,height,'INTERLACE_PIXEL',int(schematic),hdferr)

deallocate(schematic,chararr)
!=====================================================
! and finally the selected type
dataset = SC_ID
hdferr = HDF%writeDatasetInteger(dataset, 2)
 
call HDF%pop()

end subroutine h5_writeCoordinateSystemGroup

!--------------------------------------------------------------------------
subroutine h5_writePatternCenterGroup(HDF, xpc, ypc, L, delta, scdim)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! write the pattern center group

use mod_io

IMPLICIT NONE

type(HDF_T),INTENT(INOUT)            :: HDF
real(kind=sgl),INTENT(IN)            :: xpc      ! pattern center x [pixels]
real(kind=sgl),INTENT(IN)            :: ypc      ! pattern center y [pixels]
real(kind=sgl),INTENT(IN)            :: L        ! sample-scintillator distance [micron]
real(kind=sgl),INTENT(IN)            :: delta    ! scintillator pixel size [micron]
integer(kind=irg),INTENT(IN)         :: scdim(2) ! scintillator dimensions [pixels]

character(fnlen)                     :: groupname, dataset
integer(kind=irg)                    :: hdferr
real(kind=sgl)                       :: xstar, ystar, zstar

! create the Coordinate System group
groupname = 'Pattern Center Calibration'
hdferr = HDF%createGroup(groupname)

! we assume that we are writing a TSL file

! in EMsoft, the pattern center is measured in units of pixels from the 
! center of the scintillator.  For TSL, the pattern center is measured 
! from the bottom left of the scintillator (when looking towards it from the 
! sample) and in units of the width of the scintillator.

xstar = ( float(scdim(1))*0.5 + xpc ) / float(scdim(1)) 
ystar = ( float(scdim(2))*0.5 + ypc ) / float(scdim(2)) 
zstar = L / ( delta * float(scdim(1)) )

dataset = SC_xstar
hdferr = HDF%writeDatasetFloat(dataset, xstar)

dataset = SC_ystar
hdferr = HDF%writeDatasetFloat(dataset, ystar)

dataset = SC_zstar
hdferr = HDF%writeDatasetFloat(dataset, zstar)

call HDF%pop()

end subroutine h5_writePatternCenterGroup

!--------------------------------------------------------------------------
subroutine h5_writePhaseGroup(EMsoft, HDF, groupname, xtalname)
!DEC$ ATTRIBUTES DLLEXPORT :: h5ebsd_writePhaseGroup
 !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! write the phase group, describing the crystal structure

use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT)                      :: EMsoft
type(HDF_T),INTENT(INOUT)                         :: HDF
character(fnlen),intent(IN)                       :: groupname
character(fnlen),intent(IN)                       :: xtalname
                
type(HDF_T)                                       :: localHDF
type(cell_T)                                      :: cell 
type(SpaceGroup_T)                                :: SG
                
character(fnlen)                                  :: dataset, grname, filename
integer(kind=irg)                                 :: istat, SGnum, hdferr
real(kind=dbl),allocatable                        :: cellparams(:)
integer(HSIZE_T)                                  :: dims(1)
logical                                           :: readonly, stat
                
                
integer(kind=irg)                                 :: i, pgnum
character(fnlen, KIND=c_char),allocatable,TARGET  :: stringarray(:)

! TSL point group labels [courtesy of S. Wright]
character(26),parameter       :: TSLpgname(32) = (/ "Triclinic (C1) [1]        ", "Triclinic (S2, Ci) [-1]   ",&
                        "Monoclinic b (C2)[2]      ", "Monoclinic b (C1h, Cs) [m]", "Monoclinic b (C2h) [2/m]  ",&
                        "Orthorhombic (D2) [222]   ", "Orthorhombic (C2v) [mm2]  ", "Orthorhombic (D2h) [mmm]  ",&
                        "Tetragonal (C4) [4]       ", "Tetragonal (S4) [-4]      ", "Tetragonal (C4h) [4/m]    ",&
                        "Tetragonal (D4) [422]     ", "Tetragonal (C4v) [4mm]    ", "Tetragonal (D2d) [-42m]   ",&
                        "Tetragonal (D4h) [4/mmm]  ", "Trigonal (C3) [3]         ", "Trigonal (S6, C3i) [-3]   ",&
                        "Trigonal (D3) [32]        ", "Trigonal (C3v) [3m]       ", "Trigonal (D3d) [-3m]      ",&
                        "Hexagonal (C6) [6]        ", "Hexagonal (C3h) [-6]      ", "Hexagonal (C6h) [6/m]     ",&
                        "Hexagonal (D6) [622]      ", "Hexagonal (C6v) [6mm]     ", "Hexagonal (D3h) [-6m2]    ",&
                        "Hexagonal (D6h) [6/mmm]   ", "Cubic (T) [23]            ", "Cubic (Th) [m3]           ",&
                        "Cubic (O) [432]           ", "Cubic (Td) [-43m]         ", "Cubic (Oh) [m3m]          " /)

! TSL old symmetry identifiers [courtesy of S. Wright]
integer(kind=irg),parameter    :: TSLoldID(32) = (/ 1,1,2,2,2,22,22,22,4,4,4,42,42,42,42,3,3, &
                                                  32,32,32,6,6,6,62,62,62,62,23,23,43,43,43 /)


! this routine first extracts information from the xtal file and then
! puts it in the right format for the h5ebsd file format.
! This is organized by phase, so each phase is a separate numbered
! subgroup; the subgroupname is passed in as groupname

! test to make sure the input file exists and is HDF5 format
cell = cell_T( )
call cell%setFileName(xtalname)
call cell%readDataHDF(SG, EMsoft, useHDF=localHDF)

! create the subgroup [now we are back in the original HDF5 file]
hdferr = HDF%createGroup(groupname)

! the following data sets need to be created: Formula, Info, Lattice Constant a,
! b, c, alpha, beta, gamma, Laue Group, MaterialName, NumberFamilies, Point Group,
! Symmetry, hkl Families.  These last ones are typically used by the EDAX/TSL 
! software, so we do not necessarily have to fill them in here.

cellparams = cell%getLatParm()
sgnum = SG%getSpaceGroupNumber() 

! lattice parameters [in Angstrom]
dataset = 'Lattice Constant a'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(1))*10.0)
dataset = 'Lattice Constant b'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(2))*10.0)
dataset = 'Lattice Constant c'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(3))*10.0)
dataset = 'Lattice Constant alpha'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(4)))
dataset = 'Lattice Constant beta'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(5)))
dataset = 'Lattice Constant gamma'
hdferr = HDF%writeDatasetFloat(dataset, sngl(cellparams(6)))

allocate(stringarray(1))

! point group
pgnum = 0
do i=1,32
  if (SGPG(i).le.sgnum) pgnum = i
end do
dataset = 'Point Group'
stringarray(1)= trim(TSLpgname(pgnum))
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Laue group
dataset = 'Laue Group'
stringarray(1)= trim(TSLpgname(PGrot(pgnum)))
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Symmetry
dataset = SC_Symmetry
hdferr = HDF%writeDatasetInteger(dataset, TSLoldID(pgnum))

! various other strings

! Formula [extract this from the first part of xtalname]
dataset = SC_Formula
i = scan(trim(xtalname),'.')
stringarray(1) = xtalname(1:i-1)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Material name [same as Formula; this would require adding a field to the .xtal files]
dataset = SC_MaterialName
stringarray(1) = xtalname(1:i-1)
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! Info [empty string most of the time]
dataset = SC_Info
stringarray(1) = ''
hdferr = HDF%writeDatasetStringArray(dataset, stringarray, 1)

! hkl Families [this will require a bit of work !!!!!]
! this item uses the Compound data type; we will need to generate the 
! families of unique planes, and compute structure factors ... 

! in this version of the software [EMsoft 3.1], we leave these datasets empty
dataset = SC_NumberFamilies
i = 0
hdferr = HDF%writeDatasetInteger(dataset, i)
! call Message('h5ebsd_writePhaseGroup: writing of ->NumberFamilies<- data not yet implemented.')

dataset = 'hkl Families'
hdferr = HDF%writeDatasetInteger(dataset, i)
! call Message('h5ebsd_writePhaseGroup: writing of ->hkl Families<- data not yet implemented.')


! and leave this group
call HDF%pop()

end subroutine h5_writePhaseGroup




end module mod_DIsupport