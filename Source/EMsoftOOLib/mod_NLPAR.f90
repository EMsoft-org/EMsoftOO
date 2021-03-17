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
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIGBILITY OF SUCH DAMAGE.
! ###################################################################

module mod_NLPAR
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/15/21
  !!
  !! class definition for the NLPAR module
  !!
  !! This module is based on the NLPAR paper:
  !! 
  !! "NLPAR: Non-local smoothing for enhanced EBSD pattern indexing"
  !! P.T. Brewick, S.I. Wright, D.J. Rowenhorst, Ultramicroscopy 200 (2019) pp. 50-61

use mod_kinds
use mod_global

IMPLICIT NONE 

! class definition
type, public :: NLPAR_T
private 
  integer(kind=irg)     :: searchWindow
  real(kind=sgl)        :: lambda


contains
private 
  procedure, pass(self) :: setSearchWindow_
  procedure, pass(self) :: getSearchWindow_
  procedure, pass(self) :: setLambda_
  procedure, pass(self) :: getLambda_
  procedure, pass(self) :: estimateSigma_
  procedure, pass(self) :: getWeightFactors_
  procedure, pass(self) :: averagePatterns_
  procedure, pass(self) :: doNLPAR_

  generic, public :: setSearchWindow => setSearchWindow_
  generic, public :: getSearchWindow => getSearchWindow_
  generic, public :: setLambda => setLambda_
  generic, public :: getLambda => getLambda_
  generic, public :: estimateSigma => estimateSigma_
  generic, public :: getWeightFactors => getWeightFactors_
  generic, public :: averagePatterns => averagePatterns_
  generic, public :: doNLPAR => doNLPAR_

end type NLPAR_T

! the constructor routine for this class 
interface NLPAR_T
  module procedure NLPAR_constructor
end interface NLPAR_T

contains

!--------------------------------------------------------------------------
type(NLPAR_T) function NLPAR_constructor( ) result(NLPAR)
!DEC$ ATTRIBUTES DLLEXPORT :: NLPAR_constructor
!! author: MDG 
!! version: 1.0 
!! date: 03/15/21
!!
!! constructor for the NLPAR_T Class; reads the name list 
 
IMPLICIT NONE

end function NLPAR_constructor

!--------------------------------------------------------------------------
subroutine NLPAR_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: NLPAR_destructor
!! author: MDG 
!! version: 1.0 
!! date: 03/15/21
!!
!! destructor for the NLPAR_T Class
 
IMPLICIT NONE

type(NLPAR_T), INTENT(INOUT)  :: self 

call reportDestructor('NLPAR_T')

end subroutine NLPAR_destructor

!--------------------------------------------------------------------------
subroutine setSearchWindow_(self, inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setSearchWindow_
!! author: MDG
!! version: 1.0
!! date: 03/15/21
!!
!! set searchWindow parameter

IMPLICIT NONE

class(NLPAR_T), INTENT(INOUT)   :: self
integer(kind=irg), INTENT(IN)   :: inp

self%searchWindow = inp

end subroutine setSearchWindow_

!--------------------------------------------------------------------------
function getSearchWindow_(self) result(outp)
!DEC$ ATTRIBUTES DLLEXPORT :: getSearchWindow_
!! author: MDG
!! version: 1.0
!! date: 03/15/21
!!
!! get searchWindow parameter

IMPLICIT NONE

class(NLPAR_T), INTENT(INOUT)   :: self
integer(kind=irg)               :: outp

outp = self%searchWindow 

end function getSearchWindow_

!--------------------------------------------------------------------------
subroutine setLambda_(self, inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setLambda_
!! author: MDG
!! version: 1.0
!! date: 03/15/21
!!
!! set lambda parameter

IMPLICIT NONE

class(NLPAR_T), INTENT(INOUT)   :: self
real(kind=sgl), INTENT(IN)      :: inp

self%lambda= inp

end subroutine setLambda_

!--------------------------------------------------------------------------
function getLambda_(self) result(outp)
!DEC$ ATTRIBUTES DLLEXPORT :: getLambda_
!! author: MDG
!! version: 1.0
!! date: 03/17/21
!!
!! get lambda parameter

IMPLICIT NONE

class(NLPAR_T), INTENT(INOUT)   :: self
real(kind=sgl)                  :: outp

outp = self%lambda

end function getLambda_

!--------------------------------------------------------------------------
function estimateSigma_(self, pb, ps, wd, sw, jrow) result(sigEst)
!DEC$ ATTRIBUTES DLLEXPORT :: estimateSigma_
!! author: MDG 
!! version: 1.0 
!! date: 03/15/21
!!
!! estimate the sigma parameter for a given pattern (needs 3x3 surrounding patterns)
!!
!! we do this for a single row and duplicate the first and last values; this is done 
!! using multiple threads 

use omp_lib

IMPLICIT NONE 

class(NLPAR_T), INTENT(INOUT)          :: self
integer(kind=irg), INTENT(IN)          :: ps 
integer(kind=irg), INTENT(IN)          :: wd 
integer(kind=irg), INTENT(IN)          :: sw 
real(kind=sgl), INTENT(IN)             :: pb( ps * wd * (2*sw+2) )
integer(kind=irg), INTENT(IN)          :: jrow
real(kind=sgl)                         :: sigEst( wd ) 

real(kind=sgl)                         :: dpvals(8), cp(ps), nv
integer(kind=irg)                      :: ip, i, j, pstart(8), ik, pc

sigEst = 0.0

nv = 0.5/float(wd)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, ik, pstart, dpvals, pc, cp)
!$OMP DO SCHEDULE(DYNAMIC,5)
  do ip=2,wd-1
    dpvals = 0.0 
    pc = ( (jrow-1)*wd+ip-1) * ps + 1
    cp = pb(pc:pc+ps-1)
    ik = 0
    do i=-1,1
      do j=-1,1
        if ((abs(i)+abs(j)).ne.0)  then 
          ik = ik+1 
          pstart(ik) = pc + (i + j * wd) * ps 
        end if
      end do 
    end do 
    do ik=1,8 
      dpvals(ik) = sum( ( cp - pb(pstart(ik):pstart(ik)+ps-1) )**2 )
    end do
    sigEst(ip) = minval(dpvals)*nv  ! we only use this in squared form so no need for sqrt
  end do
!$OMP END DO 
!$OMP END PARALLEL

! copy the edge points; should be a good estimate
sigEst(1) = sigEst(2)
sigEst(wd) = sigEst(wd-1)

end function estimateSigma_

!--------------------------------------------------------------------------
function getWeightFactors_(self, pb, se, ps, wd, sw, lambda, row) result(wt)
!DEC$ ATTRIBUTES DLLEXPORT :: getWeightFactors_
!! author: MDG 
!! version: 1.0 
!! date: 03/15/21
!!
!! compute the (2*SW+1)x(2*SW+1) pairwise distance arrays

use omp_lib

IMPLICIT NONE 

class(NLPAR_T), INTENT(INOUT)          :: self
integer(kind=irg), INTENT(IN)          :: ps 
integer(kind=irg), INTENT(IN)          :: wd 
integer(kind=irg), INTENT(IN)          :: sw 
real(kind=sgl), INTENT(IN)             :: pb( ps * wd * (2*sw+1) )
real(kind=sgl), INTENT(IN)             :: se( wd * (2*sw+1) )
real(kind=sgl), INTENT(IN)             :: lambda
integer(kind=sgl), INTENT(IN), OPTIONAL:: row 
real(kind=sgl)                         :: wt( 2*sw+1, 2*sw+1, wd ) 

real(kind=sgl)                         :: cp(ps), diff(ps)
integer(kind=irg)                      :: ip, i, j, pstart, ik, jrow, icol, pc, fwd, sfwd, sigi, sigj

fwd = float(wd)
sfwd = sqrt(2.0*fwd)

! we need to make sure that we handle the bottom and top SW patterns differently!
! The point with indices (icol, jrow) is the reference pattern for the NLPAR algorithm
jrow = sw+1
if (present(row)) jrow=row

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, ik, ip, icol, pc, cp, pstart, sigi, sigj, diff)
!$OMP DO SCHEDULE(DYNAMIC,10)
  do ik=1,wd  ! loop over all the reference patterns 
! icol is the starting column for the double loop
    icol = ik-sw
    if (ik.le.sw) icol = 1
    if (ik.gt.wd-sw) icol=wd-2*sw
    pc = (jrow-1)*wd*ps + (icol-1)*ps + 1
    cp = pb(pc:pc+ps-1)  ! this is the reference pattern vector 
    sigi = se( (jrow-1)*wd + icol )  ! sigma_i^2 for the reference pattern 
! loop over all the patterns in the box
    do i=icol,icol+2*sw
      ip = i-icol+1
      do j=1,2*sw+1
        if ( (ip.eq.icol) .and. (j.eq.jrow)) then
          wt(ip,j,ik) = 1.0
        else
          pstart = (j-1) * wd * ps + (ip-1) * ps + 1
          sigj = sigi + se( (j-1) * wd  + ip )
          diff = cp - pb(pstart:pstart+ps-1)
          wt(ip,j,ik) = exp( - lambda*maxval( (/0.0, (sum(diff*diff)-fwd*sigj)/(sfwd*sigj) /)) ) 
        end if
      end do
    end do     
    wt(:,:,ik) = wt(:,:,ik) / sum(wt(:,:,ik))
  end do 
!$OMP END DO 
!$OMP END PARALLEL

end function getWeightFactors_


!--------------------------------------------------------------------------
function averagePatterns_(self, pb, wt, ps, wd, sw) result(pat)
!DEC$ ATTRIBUTES DLLEXPORT :: averagePatterns_
!! author: MDG 
!! version: 1.0 
!! date: 03/15/21
!!
!! compute the (2*SW+1)x(2*SW+1) pairwise distance arrays

use omp_lib

IMPLICIT NONE 

class(NLPAR_T), INTENT(INOUT)          :: self
integer(kind=irg), INTENT(IN)          :: ps 
integer(kind=irg), INTENT(IN)          :: wd 
integer(kind=irg), INTENT(IN)          :: sw 
real(kind=sgl), INTENT(IN)             :: pb( ps * wd * (2*sw+1) )
real(kind=sgl), INTENT(IN)             :: wt( 2*sw+1, 2*sw+1, wd )
real(kind=sgl)                         :: pat( ps * wd ) 

real(kind=sgl)                         :: psum(ps)
integer(kind=irg)                      :: ip, i, j, pstart, ik, jrow, icol

pat = 0.0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, ik, ip, icol, psum, pstart)
!$OMP DO SCHEDULE(DYNAMIC,10)
  do ik=1,wd  ! loop over all the reference patterns 
! icol is the starting column for the double loop
    icol = ik-sw
    if (ik.le.sw) icol = 1
    if (ik.gt.wd-sw) icol=wd-2*sw
! next we sum all the weighted patterns in each block 
    psum = 0.0
    do i=icol,icol+2*sw
      ip = i-icol+1
      do j=1,2*sw+1
        pstart = (j-1) * wd * ps + (icol-1) * ps + 1
        psum(1:ps) = psum(1:ps) + wt(ip,j,ik) * pb(pstart:pstart+ps-1)
      end do 
    end do 
    pat((ik-1)*ps+1:ik*ps) = psum(1:ps)
  end do 
!$OMP END DO 
!$OMP END PARALLEL


end function averagePatterns_

!--------------------------------------------------------------------------
recursive subroutine doNLPAR_(self, EMsoft, HDF, inRAM, nml, binx, biny, masklin, correctsize, totnumexpt, &
                              epatterns, exptIQ)
!DEC$ ATTRIBUTES DLLEXPORT :: doNLPAR_
!! author: MDG
!! version: 1.0
!! date: 03/15/21
!!
!! apply the NLPAR averaging method to the pattern file and generate a preprocessed pattern file 
!! this includes the normal pattern preprocessing step performed for DI etc... programs, in other 
!! words, this is a more complicated version of the PreProcessPatterns routine in mod_patterns, with 
!! some duplication of code; this is unavoidable because NLPAR needs to read multiple lines of patterns
!! instead of just one, and it would needlessly complicate the PreProcessPatterns routine to do it the 
!! other way around. 

use mod_EMsoft
use mod_io
use omp_lib
use mod_filters
use mod_timing
use mod_vendors
use mod_DIfiles
use mod_patterns
use HDF5
use mod_OMPsupport
use ISO_C_BINDING
use mod_memory

class(NLPAR_T), INTENT(INOUT)                     :: self
type(EMsoft_T),INTENT(INOUT)                      :: EMsoft
type(HDF_T),INTENT(INOUT)                         :: HDF
logical,INTENT(IN)                                :: inRAM
class(DictionaryIndexingNameListType),INTENT(IN)  :: nml
integer(kind=irg),INTENT(IN)                      :: binx
integer(kind=irg),INTENT(IN)                      :: biny
real(kind=sgl),INTENT(IN)                         :: masklin(binx*biny)
integer(kind=irg),INTENT(IN)                      :: correctsize
integer(kind=irg),INTENT(IN)                      :: totnumexpt
real(kind=sgl),INTENT(INOUT),OPTIONAL             :: epatterns(correctsize, totnumexpt)
!f2py intent(in,out) ::  epatterns
real(kind=sgl),INTENT(INOUT),OPTIONAL             :: exptIQ(totnumexpt)
!f2py intent(in,out) ::  exptIQ

type(IO_T)                                        :: Message
type(Vendor_T)                                    :: VT
type(timing_T)                                    :: timer
type(memory_T)                                    :: mem 

integer(kind=irg)                                 :: itype, istat, L, recordsize, patsz, iiistart, iiiend, jjend, ierr, &
                                                     SW, i, iii, ii, j, jj, jrow, kk 
integer(kind=irg)                                 :: tickstart, tstop, io_int(2)
integer(HSIZE_T)                                  :: dims3(3), offset3(3)
integer(kind=irg),parameter                       :: iunitexpt = 41, itmpexpt = 42
real(kind=dbl)                                    :: w, Jres
real(kind=sgl)                                    :: vlen, tmp, mi, ma, lambda, fwd, sfwd 
real(kind=sgl),allocatable                        :: Pattern(:,:), imageexpt(:), exppatarray(:), Pat(:,:), exptblock(:), &
                                                     sigvals(:,:), wtfactors(:,:,:), sigEst(:), dpp(:,:,:), tmpimageexpt(:)
real(kind=dbl),allocatable                        :: rrdata(:,:), ffdata(:,:), ksqarray(:,:)
integer(kind=irg),allocatable                     :: pint(:,:)
complex(kind=dbl),allocatable                     :: hpmask(:,:)
complex(C_DOUBLE_COMPLEX),allocatable             :: inp(:,:), outp(:,:)
type(C_PTR)                                       :: planf, HPplanf, HPplanb
logical                                           :: isEBSD = .FALSE., isTKD = .FALSE., ROIselected, f_exists, verbose
character(fnlen)                                  :: fname 

SW = self%getSearchWindow_()
lambda = 1.0/self%getLambda_()**2
mem = memory_T()

verbose = .FALSE.

call Message%printMessage(' Performing NLPAR averaging on experimental patterns')

if (nml%DIModality.eq.'EBSD') then
  isEBSD = .TRUE.
  else if (nml%DIModality.eq.'TKD') then
    isTKD = .TRUE.
  else
    call Message%printError('PreProcessPatterns', 'unknown name list type requested')
  end if

!===================================================================================
! define a bunch of mostly integer parameters
recordsize = correctsize*4
L = binx*biny
patsz = correctsize
w = nml%hipassw
fwd = float(L)
sfwd = sqrt(2.0*float(L))

if (sum(nml%ROI).ne.0) then
  ROIselected = .TRUE.
  iiistart = nml%ROI(2)
  iiiend = nml%ROI(2)+nml%ROI(4)-1
  jjend = nml%ROI(3)
else
  ROIselected = .FALSE.
  iiistart = 1
  iiiend = nml%ipf_ht
  jjend = nml%ipf_wd
end if

fname = trim(EMsoft%generateFilePath('EMtmppathname'))//trim(nml%tmpfile)

call Message%WriteValue('Creating temporary file :',trim(fname))

if (f_exists) then  ! delete the file if it already exists
  open(unit=itmpexpt,file=trim(fname),&
       status='unknown',form='unformatted',access='direct',recl=recordsize,iostat=ierr)
  close(unit=itmpexpt,status='delete')
end if
open(unit=itmpexpt,file=trim(fname),&
    status='unknown',form='unformatted',access='direct',recl=recordsize,iostat=ierr)


!===================================================================================
! set the vendor inputtype for the pattern file and, if necessary, declare the HDF class
VT = Vendor_T( nml%inputtype )
itype = VT%get_itype()
call VT%set_filename(nml%exptfile)

!===================================================================================
! open the file with experimental patterns; depending on the inputtype parameter, this
! can be a regular binary file, as produced by a MatLab or IDL script (default); a
! pattern file produced by EMEBSD.f90 etc.; or a vendor binary or HDF5 file... in each case we need to
! open the file and leave it open, then use the getExpPatternRow() routine to read a row
! of patterns into the exppatarray variable ...  at the end, we use closeExpPatternFile() to
! properly close the experimental pattern file
if ( (itype.eq.4) .or. (itype.eq.7) .or. (itype.eq.8) ) then
  istat = VT%openExpPatternFile(EMsoft, nml%ipf_wd, L, recordsize, nml%HDFstrings, HDF)
else
  istat = VT%openExpPatternFile(EMsoft, nml%ipf_wd, L, recordsize)
end if

if (istat.ne.0) then
    call Message%printError("PreProcessPatterns:", "Fatal error handling experimental pattern file")
end if

! this next part is done with OpenMP, with only thread 0 doing the reading;
! Thread 0 reads one or more lines worth of patterns from the input file, then all threads do
! the work, and thread 0 adds them to the epatterns array in RAM or the temporary output file; 
! repeat until all patterns have been processed.
call OMP_setNThreads(nml%nthreads)

! allocate the arrays that holds the experimental patterns from (2*SW+1) rows of the region of interest
call mem%alloc1(exppatarray, (/patsz * nml%ipf_wd/), 'exppatarray')
call mem%alloc1(exptblock, (/patsz * nml%ipf_wd * (2*SW+1)/), 'exptblock')
call mem%alloc1(sigEst, (/nml%ipf_wd * (2*SW+1)/), 'sigEst')

if (present(exptIQ)) then
! prepare the fftw plan for this pattern size to compute pattern quality (pattern sharpness Q)
  call mem%alloc2(Pat, (/binx,biny/),'Pat')
  call mem%alloc2(ksqarray, (/binx,biny/),'ksqarray')
  Jres = 0.0
  call init_getEBSDIQ(binx, biny, Pat, ksqarray, Jres, planf)
  call mem%dealloc2(Pat,'Pat')
end if

! initialize the HiPassFilter routine (has its own FFTW plans)
call mem%alloc2(hpmask, (/binx,biny/),'hpmask')
call mem%alloc2(inp, (/binx,biny/),'inp')
call mem%alloc2(outp, (/binx,biny/),'outp')
call init_HiPassFilter(w, (/ binx, biny /), hpmask, inp, outp, HPplanf, HPplanb)
call mem%dealloc2(inp,'inp')
call mem%dealloc2(outp,'outp')
! keep the hpmask array 

call mem%allocated_memory_use()

call Message%printMessage('Starting processing of experimental patterns')
timer = Timing_T()
call timer%Time_tick(1)

!==================================================
! the bottom SW rows require some special treatment 
!==================================================

! Loop over all the rows and make sure that we always have 2*SW+1 rows in the exptblock array; so,
! we read the first 2*SW+1 rows and then start the row loop.  To keep things simple, we read complete rows 
! and apply the ROI afterwards if it is defined.
dims3 = (/ binx, biny, nml%ipf_wd /)
do jrow = 1, 2*SW+1+1  ! the second +1 is to have enough rows for the sigEst routine
  offset3 = (/ 0, 0, (jrow-1)*nml%ipf_wd /)
  if ( (itype.eq.4) .or. (itype.eq.7) .or. (itype.eq.8) ) then
    call VT%getExpPatternRow(jrow, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray, &
                             HDFstrings=nml%HDFstrings, HDF=HDF)
  else
    call VT%getExpPatternRow(jrow, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray)
  end if
  exptblock( (jrow-1) * patsz * nml%ipf_wd+1: jrow * patsz * nml%ipf_wd ) = exppatarray(1:patsz * nml%ipf_wd)
  if (verbose.eqv..TRUE.) write (*,*) ' read line ',jrow,' of ',nml%ipf_ht
end do 

! we need to compute the sigma values for the bottom 2 rows 
! and we'll take the bottom row equal to the first row to avoid 
! having to deal with an incomplete nearest neighbor set.
! We send the first three concecutive rows to the estimateSigma_ routine.
do jrow=2,2*SW+1
  sigEst((jrow-1)*nml%ipf_wd+1:jrow*nml%ipf_wd) = self%estimateSigma_(exptblock, patsz, nml%ipf_wd, SW, jrow)
  if (verbose.eqv..TRUE.) write (*,*) ' computed sigma^2 for row ', jrow 
end do 
! make sure to copy the second row into the first row 
sigEst(1:nml%ipf_wd) = sigEst(nml%ipf_wd+1:2*nml%ipf_wd)
if (verbose.eqv..TRUE.) write (*,*) ' copied sigma^2 for row 1'

!==================================================
! The following is a very complicated main loop!  
!==================================================
call mem%alloc3(wtfactors, (/2*SW+1,2*SW+1,nml%ipf_wd/), 'wtfactors')
do jrow=1,nml%ipf_ht
! loop over all the experimental rows.  First we read the next line but only
! if we are in the central region.
  if ( (jrow.gt.SW+1) .and. (jrow.lt.(nml%ipf_ht-SW+1)) ) then   ! we need to shift arrays and read the following row of patterns 
    exptblock = cshift(exptblock, nml%ipf_wd * patsz)  ! experimental patterns 
    sigEst = cshift(sigEst, nml%ipf_wd)                ! estimated sigma^2 values
    offset3 = (/ 0, 0, (jrow+SW)*nml%ipf_wd /)
    if ( (itype.eq.4) .or. (itype.eq.7) .or. (itype.eq.8) ) then
      call VT%getExpPatternRow(jrow+SW+1, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray, &
                               HDFstrings=nml%HDFstrings, HDF=HDF)
    else
      call VT%getExpPatternRow(jrow+SW+1, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray)
    end if
    ! we add this as the last row of the exptblock array
    exptblock( (2 * SW +1) * patsz * nml%ipf_wd+1: (2 * SW + 2) * patsz * nml%ipf_wd ) = exppatarray(1:patsz * nml%ipf_wd)
    if (verbose.eqv..TRUE.) write (*,*) ' read line ',jrow+SW+1,' of ',nml%ipf_ht
! at this point we might as well compute the estimated sigma^2 values for the penultimate 
! row in the exptblock array ... except for when we are at the end for which we need to 
! copy a row ... if jrow = nml%ipf_ht-SW
    if (jrow.eq.(nml%ipf_ht-SW)) then  ! copy the last line
      sigEst((2*SW)*nml%ipf_wd+1:(2*SW+1)*nml%ipf_wd) = sigEst((2*SW-1)*nml%ipf_wd+1:(2*SW)*nml%ipf_wd)
      if (verbose.eqv..TRUE.) write (*,*) ' copied sigma^2 for row ', jrow 
    else 
      sigEst((2*SW)*nml%ipf_wd+1:(2*SW+1)*nml%ipf_wd) = self%estimateSigma_(exptblock, patsz, nml%ipf_wd, SW, 2*SW+1)
      if (verbose.eqv..TRUE.) write (*,*) ' computed sigma^2 for row ', jrow 
    end if 
  end if 

! get the weightfactor array (this includes the normalized distance computation)
  if (jrow.le.SW+1) then ! we are in the bottom block of SW+1 patterns without loading
    wtfactors = self%getWeightFactors_(exptblock, sigEst, patsz, nml%ipf_wd, SW, lambda, row=jrow)
    if (verbose.eqv..TRUE.) write (*,*) ' getWeightFactors_ fixed',jrow 
  else 
    if (jrow.ge.nml%ipf_ht-SW) then ! same for the top block of SW+1 patterns
      wtfactors = self%getWeightFactors_(exptblock, sigEst, patsz, nml%ipf_wd, SW, lambda, row=2*SW+1-(nml%ipf_ht-jrow))
      if (verbose.eqv..TRUE.) write (*,*) ' getWeightFactors_ fixed',2*SW+1-(nml%ipf_ht-jrow) 
    else  ! we just loaded a new pattern row so we recompute distances
      wtfactors = self%getWeightFactors_(exptblock, sigEst, patsz, nml%ipf_wd, SW, lambda)
      if (verbose.eqv..TRUE.) write (*,*) ' getWeightFactors_ ',jrow 
    end if
  end if 

! and compute the weighted sum of patterns 
  exppatarray = self%averagePatterns(exptblock, wtfactors, patsz, nml%ipf_wd, SW)
  if (verbose.eqv..TRUE.) write (*,*) ' averagePatterns ',jrow

! next we perform the normal pre-processing steps on these averaged pattern vectors 
! these lines are taken from the PreProcessPatterns routine; maybe in a later version
! we can avoid duplicating this code segment (it also appears twice further down)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jj, kk, mi, ma, istat) &
!$OMP& PRIVATE(imageexpt, tmpimageexpt, Pat, rrdata, ffdata, pint, vlen, tmp, inp, outp)

  allocate(Pat(binx,biny),rrdata(binx,biny),ffdata(binx,biny),tmpimageexpt(binx*biny),stat=istat)
  if (istat .ne. 0) stop 'could not allocate arrays for Hi-Pass filter'

  allocate(pint(binx,biny),stat=istat)
  if (istat .ne. 0) stop 'could not allocate pint array'

  allocate(inp(binx,biny),outp(binx,biny),stat=istat)
  if (istat .ne. 0) stop 'could not allocate inp, outp arrays'

  tmpimageexpt = 0.0
  rrdata = 0.D0
  ffdata = 0.D0

!$OMP DO SCHEDULE(DYNAMIC)
  do jj=1,nml%ipf_wd
! convert imageexpt to 2D EBSD Pattern array
      do kk=1,biny
        Pat(1:binx,kk) = exppatarray((jj-1)*patsz+(kk-1)*binx+1:(jj-1)*patsz+kk*binx)
      end do

      if (present(exptIQ)) then
! compute the pattern Image Quality
        exptIQ((iii-iiistart)*jjend + jj) = sngl(computeEBSDIQ(binx, biny, Pat, ksqarray, Jres, planf))
      end if

! Hi-Pass filter
      rrdata = dble(Pat)
      ffdata = applyHiPassFilter(rrdata, (/ binx, biny /), w, hpmask, inp, outp, HPplanf, HPplanb)
      Pat = sngl(ffdata)

! adaptive histogram equalization
      ma = maxval(Pat)
      mi = minval(Pat)

      pint = nint(((Pat - mi) / (ma-mi))*255.0)
      Pat = float(adhisteq(nml%nregions,binx,biny,pint))

! convert back to 1D vector
      do kk=1,biny
        exppatarray((jj-1)*patsz+(kk-1)*binx+1:(jj-1)*patsz+kk*binx) = Pat(1:binx,kk)
      end do

! apply circular mask and normalize for the dot product computation
      exppatarray((jj-1)*patsz+1:(jj-1)*patsz+L) = exppatarray((jj-1)*patsz+1:(jj-1)*patsz+L) * masklin(1:L)
      vlen = vecnorm(exppatarray((jj-1)*patsz+1:(jj-1)*patsz+L))
      if (vlen.ne.0.0) then
        exppatarray((jj-1)*patsz+1:(jj-1)*patsz+L) = exppatarray((jj-1)*patsz+1:(jj-1)*patsz+L)/vlen
      else
        exppatarray((jj-1)*patsz+1:(jj-1)*patsz+L) = 0.0
      end if
  end do
!$OMP END DO
deallocate(tmpimageexpt, Pat, rrdata, ffdata, pint, inp, outp)
!$OMP BARRIER
!$OMP END PARALLEL

! and we either write the resulting patterns to the temp file or we keep them in RAM 
      if (inRAM.eqv..TRUE.) then
        do jj=1,jjend
          epatterns(1:patsz,(jrow-iiistart)*jjend + jj) = exppatarray((jj-1)*patsz+1:jj*patsz)
        end do
      else
        do jj=1,jjend
          write(itmpexpt,rec=(jrow-iiistart)*jjend + jj) exppatarray((jj-1)*patsz+1:jj*patsz)
        end do
      end if

! print an update of progress
    if (mod(jrow-iiistart+1,5).eq.0) then
      ! if (ROIselected.eqv..TRUE.) then
      !   io_int(1:2) = (/ jrow-iiistart+1, nml%ROI(4) /)
      !   call Message%WriteValue('Completed row ',io_int,2,"(I4,' of ',I4,' rows')")
      ! else
        io_int(1:2) = (/ jrow-iiistart+1, nml%ipf_ht /)
        call Message%WriteValue('Completed row ',io_int,2,"(I4,' of ',I4,' rows')")
      ! end if
    end if

end do 

close(unit=itmpexpt,status = 'keep')

!====================================
! that should be it... some clean up and we return to the calling program
!====================================
call mem%dealloc1(exppatarray, 'exppatarray')
call mem%dealloc1(exptblock, 'exptblock')
call mem%dealloc1(sigEst, 'sigEst')
call mem%dealloc2(hpmask, 'hpmask')
call mem%dealloc3(wtfactors, 'wtfactors')

if (present(exptIQ)) then
  call mem%dealloc2(ksqarray, 'ksqarray')
end if

call mem%allocated_memory_use()

end subroutine doNLPAR_


end module mod_NLPAR