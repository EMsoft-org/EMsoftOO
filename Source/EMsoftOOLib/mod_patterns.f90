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

module mod_patterns
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/31/20
  !!
  !! a collection of routines for operations on diffraction patterns; these are 
  !! mostly standalone routines that are called by a variety of programs. They 
  !! belong together in a loose logical way... In the older EMsoft version, they
  !! were scattered over multiple files: patternmod, commonmod, etc...

use mod_kinds
use mod_global

IMPLICIT NONE 

contains

!DEC$ ATTRIBUTES DLLEXPORT :: PreProcessPatterns
!DEC$ ATTRIBUTES DLLEXPORT :: init_getEBSDIQ
!DEC$ ATTRIBUTES DLLEXPORT :: computeEBSDIQ

!--------------------------------------------------------------------------
recursive subroutine init_getEBSDIQ(dimx, dimy, pattern, ksqarray, Jres, planf) 
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! initialize variables for the EBSD Image Quality using the second moment of the power spectrum
!!
!! This is based on Krieger Lassen's pattern sharpness definition: Q = 1 - J / Jres wtot
!! more details page 93 of thesis of Farangis Ram.

use mod_FFTW3

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: dimx
integer(kind=irg),INTENT(IN)            :: dimy
real(kind=sgl),INTENT(IN)               :: pattern(dimx,dimy)
real(kind=dbl),INTENT(OUT)              :: ksqarray(dimx,dimy)
real(kind=dbl),INTENT(OUT)              :: Jres
type(C_PTR),INTENT(OUT)                 :: planf

complex(C_DOUBLE_COMPLEX),pointer       :: inp(:,:)
complex(C_DOUBLE_COMPLEX),pointer       :: outp(:,:)
type(C_PTR)                             :: p, o
real(kind=dbl)                          :: linex(dimx), liney(dimy)
integer(kind=irg)                       :: i

p = fftw_alloc_complex(int(dimx*dimy,C_SIZE_T))
call c_f_pointer(p, inp, [dimx,dimy])

o = fftw_alloc_complex(int(dimx*dimy,C_SIZE_T))
call c_f_pointer(o, outp, [dimx,dimy])

inp = cmplx(0.D0,0D0)
outp = cmplx(0.D0,0.D0)

! set up the fftw plan for the forward transform
planf = fftw_plan_dft_2d(dimy,dimx,inp,outp, FFTW_FORWARD,FFTW_ESTIMATE)

! generate the parameter/array needed by the getEBSDIQ function
ksqarray = 0.D0
Jres = 0.D0

linex = (/ (dble(i),i=0,dimx-1) /) 
linex(dimx/2+1:dimx) = linex(dimx/2+1:dimx) - dble(dimx)
linex = linex**2
liney = (/ (dble(i),i=0,dimy-1) /) 
liney(dimy/2+1:dimy) = liney(dimy/2+1:dimy) - dble(dimy)
liney = liney**2

do i=1,dimx
    ksqarray(i,1:dimy) = linex(i) + liney(1:dimy)
end do
Jres = sum(ksqarray) / dble(dimx) / dble(dimy)

call fftw_free(p)
call fftw_free(o)
call fftw_cleanup()

end subroutine init_getEBSDIQ

!--------------------------------------------------------------------------
recursive function computeEBSDIQ(dimx, dimy, pattern, ksqarray, Jres, planf) result(Q)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! compute the EBSD Image Quality using the second moment of the power spectrum
!!
!! This is based on Krieger Lassen's pattern sharpness definition: Q = 1 - J / Jres wtot
!! more details page 93 of thesis of Farangis Ram.

use mod_FFTW3

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: dimx
integer(kind=irg),INTENT(IN)            :: dimy
real(kind=sgl),INTENT(IN)               :: pattern(dimx,dimy)
real(kind=dbl),INTENT(IN)               :: ksqarray(dimx,dimy)
real(kind=dbl),INTENT(IN)               :: Jres
type(C_PTR),INTENT(IN)                  :: planf
real(kind=dbl)                          :: Q

real(kind=dbl)                          :: J, wtot
complex(C_DOUBLE_COMPLEX),pointer       :: inp(:,:)
complex(C_DOUBLE_COMPLEX),pointer       :: outp(:,:)
type(C_PTR)                             :: p, o
real(kind=dbl)                          :: w(dimx,dimy), linex(dimx), liney(dimy)
integer(kind=irg)                       :: i

p = fftw_alloc_complex(int(dimx*dimy,C_SIZE_T))
call c_f_pointer(p, inp, [dimx,dimy])

o = fftw_alloc_complex(int(dimx*dimy,C_SIZE_T))
call c_f_pointer(o, outp, [dimx,dimy])

inp = pattern
outp = cmplx(0.D0,0.D0)

! compute the Fourier transform
call fftw_execute_dft(planf, inp, outp)

w = sqrt(real(outp)**2 + aimag(outp)**2)

! sum over the arrays
J = sum(w*ksqarray)
wtot = sum(w)

! and return the quality parameter
Q = 1.0 - J/Jres/wtot

call fftw_free(p)
call fftw_free(o)
call fftw_cleanup()

end function computeEBSDIQ

!--------------------------------------------------------------------------
recursive subroutine PreProcessPatterns(EMsoft, HDF, inRAM, nml, binx, biny, masklin, correctsize, totnumexpt, &
                                        epatterns, exptIQ)
!! author: MDG 
!! version: 1.0 
!! date: 03/31/20
!!
!! standard preprocessing of EBSD/TKD patterns (hi-pass filter+adaptive histogram equalization) 
!!
!! This is one of the core routines used to pre-process patterns for dictionary indexing.
!! The routine reads the experimental patterns row by row, and performs a hi-pass filter and adaptive
!! histogram equalization.  Then the preprocessed patterns are either stored in a binary direct access file,
!! or they are kept in RAM, depending on the user parameter setting. 

use mod_EMsoft
use mod_io
use omp_lib
use mod_filters
use mod_timing
use mod_DIfiles
use mod_vendors
use HDF5
use mod_OMPsupport
use ISO_C_BINDING

IMPLICIT NONE

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

logical                                           :: ROIselected, f_exists
character(fnlen)                                  :: fname
integer(kind=irg)                                 :: istat, L, recordsize, io_int(2), patsz, iii, &
                                                     iiistart, iiiend, jjend, TID, jj, kk, ierr, itype
integer(HSIZE_T)                                  :: dims3(3), offset3(3)
integer(kind=irg),parameter                       :: iunitexpt = 41, itmpexpt = 42
integer(kind=irg)                                 :: tickstart, tstop 
real(kind=sgl)                                    :: vlen, tmp, ma, mi, io_real(1)
real(kind=dbl)                                    :: w, Jres
integer(kind=irg),allocatable                     :: pint(:,:)
real(kind=sgl),allocatable                        :: tmpimageexpt(:), Pattern(:,:), imageexpt(:), exppatarray(:), Pat(:,:)
real(kind=dbl),allocatable                        :: rrdata(:,:), ffdata(:,:), ksqarray(:,:)
complex(kind=dbl),allocatable                     :: hpmask(:,:)
complex(C_DOUBLE_COMPLEX),allocatable             :: inp(:,:), outp(:,:)
type(C_PTR)                                       :: planf, HPplanf, HPplanb
logical                                           :: isEBSD = .FALSE., isTKD = .FALSE.

call Message%printMessage(' Preprocessing experimental patterns')

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

!===================================================================================
! open the output file if the patterns are not to be kept in memory
if (inRAM.eqv..FALSE.) then
! first, make sure that this file does not already exist
   f_exists = .FALSE.
   if (nml%tmpfile(1:1).ne.trim(EMsoft%getConfigParameter('EMsoftnativedelimiter'))) then 
     fname = trim(EMsoft%generateFilePath('EMtmppathname'))//trim(nml%tmpfile)
   else 
     fname = trim(nml%tmpfile)
   end if
   inquire(file=trim(fname), exist=f_exists)

   call Message%WriteValue('Creating temporary file :',trim(fname))

   if (f_exists) then  ! delete the file if it already exists 
      open(unit=itmpexpt,file=trim(fname),&
           status='unknown',form='unformatted',access='direct',recl=recordsize,iostat=ierr)
      close(unit=itmpexpt,status='delete')
   end if
   open(unit=itmpexpt,file=trim(fname),&
        status='unknown',form='unformatted',access='direct',recl=recordsize,iostat=ierr)
end if
  
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
! Thread 0 reads one line worth of patterns from the input file, then all threads do 
! the work, and thread 0 adds them to the epatterns array in RAM; repeat until all patterns have been processed.
call OMP_setNThreads(nml%nthreads)

! allocate the arrays that holds the experimental patterns from a single row of the region of interest
allocate(exppatarray(patsz * nml%ipf_wd),stat=istat)
if (istat .ne. 0) stop 'could not allocate exppatarray'

if (present(exptIQ)) then
! prepare the fftw plan for this pattern size to compute pattern quality (pattern sharpness Q)
  allocate(Pat(binx,biny),stat=istat)
  if (istat .ne. 0) stop 'could not allocate arrays for Pat filter'
  Pat = 0.0
  allocate(ksqarray(binx,biny),stat=istat)
  if (istat .ne. 0) stop 'could not allocate ksqarray array'
  Jres = 0.0
  call init_getEBSDIQ(binx, biny, Pat, ksqarray, Jres, planf)
  deallocate(Pat)
end if

! initialize the HiPassFilter routine (has its own FFTW plans)
allocate(hpmask(binx,biny),inp(binx,biny),outp(binx,biny),stat=istat)
if (istat .ne. 0) stop 'could not allocate hpmask array'
call init_HiPassFilter(w, (/ binx, biny /), hpmask, inp, outp, HPplanf, HPplanb) 
deallocate(inp, outp)

call Message%printMessage('Starting processing of experimental patterns')
timer = Timing_T()
call timer%Time_tick(1)

dims3 = (/ binx, biny, nml%ipf_wd /)

!===================================================================================
! we do one row at a time
prepexperimentalloop: do iii = iiistart,iiiend

! start the OpenMP portion
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, jj, kk, mi, ma, istat) &
!$OMP& PRIVATE(imageexpt, tmpimageexpt, Pat, rrdata, ffdata, pint, vlen, tmp, inp, outp)

! set the thread ID
    TID = OMP_GET_THREAD_NUM()

! initialize thread private variables
    allocate(Pat(binx,biny),rrdata(binx,biny),ffdata(binx,biny),tmpimageexpt(binx*biny),stat=istat)
    if (istat .ne. 0) stop 'could not allocate arrays for Hi-Pass filter'

    allocate(pint(binx,biny),stat=istat)
    if (istat .ne. 0) stop 'could not allocate pint array'

    allocate(inp(binx,biny),outp(binx,biny),stat=istat)
    if (istat .ne. 0) stop 'could not allocate inp, outp arrays'

    tmpimageexpt = 0.0
    rrdata = 0.D0
    ffdata = 0.D0

! thread 0 reads the next row of patterns from the input file
! we have to allow for all the different types of input files here...
    if (TID.eq.0) then
        offset3 = (/ 0, 0, (iii-1)*nml%ipf_wd /)
        if (ROIselected.eqv..TRUE.) then
          if ( (itype.eq.4) .or. (itype.eq.7) .or. (itype.eq.8) ) then
            call VT%getExpPatternRow(iii, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray, nml%ROI, &
                                     HDFstrings=nml%HDFstrings, HDF=HDF)
          else 
            call VT%getExpPatternRow(iii, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray, nml%ROI)
          end if 
        else
         if ( (itype.eq.4) .or. (itype.eq.7) .or. (itype.eq.8) ) then
            call VT%getExpPatternRow(iii, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray, &
                                     HDFstrings=nml%HDFstrings, HDF=HDF)
          else 
            call VT%getExpPatternRow(iii, nml%ipf_wd, patsz, L, dims3, offset3, exppatarray)
          end if 
        end if
    end if

! other threads must wait until T0 is ready
!$OMP BARRIER
    jj=0

! then loop in parallel over all patterns to perform the preprocessing steps
!$OMP DO SCHEDULE(DYNAMIC)
    do jj=1,jjend
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

! thread 0 writes the row of patterns to the epatterns array or to the file, depending on the enl%inRAM parameter
    if (TID.eq.0) then
      if (inRAM.eqv..TRUE.) then
        do jj=1,jjend
          epatterns(1:patsz,(iii-iiistart)*jjend + jj) = exppatarray((jj-1)*patsz+1:jj*patsz)
        end do
      else
        do jj=1,jjend
          write(itmpexpt,rec=(iii-iiistart)*jjend + jj) exppatarray((jj-1)*patsz+1:jj*patsz)
        end do
      end if
    end if

deallocate(tmpimageexpt, Pat, rrdata, ffdata, pint, inp, outp)
!$OMP BARRIER
!$OMP END PARALLEL

! print an update of progress
    if (mod(iii-iiistart+1,5).eq.0) then
      if (ROIselected.eqv..TRUE.) then
        io_int(1:2) = (/ iii-iiistart+1, nml%ROI(4) /)
        call Message%WriteValue('Completed row ',io_int,2,"(I4,' of ',I4,' rows')")
      else
        io_int(1:2) = (/ iii-iiistart+1, nml%ipf_ht /)
        call Message%WriteValue('Completed row ',io_int,2,"(I4,' of ',I4,' rows')")
      end if
    end if
end do prepexperimentalloop

call Message%printMessage(' -> experimental patterns preprocessed')

call VT%closeExpPatternFile()

if (inRAM.eqv..FALSE.) then
  close(unit=itmpexpt,status='keep')
end if

! print some timing information
call timer%Time_tock(1) 
tstop = timer%getInterval(1)
if (tstop.eq.0.0) then 
  call Message%printMessage(' # experimental patterns processed per second : ? [time shorter than system time resolution] ')
else
  io_real(1) = float(nml%nthreads) * float(totnumexpt)/tstop
  call Message%WriteValue(' # experimental patterns processed per second : ',io_real,1,"(F10.1,/)")
end if

end subroutine PreProcessPatterns

end module mod_patterns