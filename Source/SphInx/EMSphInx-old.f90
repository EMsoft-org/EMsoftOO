!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!*                                                                     *
!* Copyright (c) 2019-2024, De Graef Group, Carnegie Mellon University *
!* All rights reserved.                                                *
!*                                                                     *
!* Author: William C. Lenthe                                           *
!*                                                                     *
!* EMSphInx is available for academic or non-profit non-commercial     *
!* research use. Please, see the license.txt file in this distribution *
!* for further details.                                                *
!*                                                                     *
!* Interested in a commercial license? Contact:                        *
!*                                                                     *
!* Center for Technology Transfer and Enterprise Creation              *
!* 4615 Forbes Avenue, Suite 302                                       *
!* Pittsburgh, PA 15213                                                *
!*                                                                     *
!* phone. : 412.268.7393                                               *
!* email  : innovation@cmu.edu                                         *
!* website: https://www.cmu.edu/cttec/                                 *
!*                                                                     *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!--------------------------------------------------------------------------
! EMsoft:EMSphInx.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMSphInx
!
!> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
!
!> @brief Indexing of EBSD patterns using the spherical indexing approach. 
!
!> @date 06/11/19 MDG 1.0 original, starting from EMEBSDDI.f90
!> @date 07/30/19 WCL 1.1 cleaning up compilation issues
!> @date 08/04/19 MDG 1.2 added OpenMP handling
!> @date 08/05/19 MDG 1.3 added HDF5 output
!--------------------------------------------------------------------------
program EMSphInx

  use mod_kinds
  use files
  use mod_io
  use mod_EBSDmod
  use HDF5
  use mod_HDFsupport

IMPLICIT NONE

  character(fnlen                 )             :: nmldeffile, progname, progdesc
  type     (SphInxNameListType    )             :: sinl
  type     (MCCLNameListType      )             :: mcnl
  type     (EBSDMasterNameListType)             :: mpnl
  type     (EBSDMCdataType        )             :: EBSDMCdata
  type     (EBSDMPdataType        )             :: EBSDMPdata
  logical                                       :: verbose
  real     (kind=dbl              ),allocatable :: mLPNH(:,:), mLPSH(:,:), weights(:)
  integer                                       :: d, hdferr, res, i, n

  nmldeffile = 'EMSphInx.nml'
  progname   = 'EMSphInx.f90'
  progdesc   = 'Program to index EBSD patterns using Spherical Indexing'
  verbose    = .TRUE.

  ! print some information
  call EMsoft(progname, progdesc)

  ! deal with the command line arguments, if any
  call Interpret_Program_Arguments(nmldeffile,1,(/ 120 /), progname)

  ! deal with the namelist stuff
  res = index(nmldeffile,'.nml',kind=irg)
  if (res.eq.0) then
    call FatalError('EMSphInx','JSON input not yet implemented')
  !  call JSONreadEBSDIndexingNameList(sinl, nmldeffile, error_cnt)
  else
    call GetSphInxNameList(nmldeffile,sinl)
  end if

  ! here we read the master pattern data and everything else that is needed
  ! for the spherical indexing approach; we also pre-compute the SHT of the 
  ! master pattern, so that we no longer need to keep any of these data arrays.

  ! 1. read the Monte Carlo data file
  call h5open_EMsoft(hdferr)
  call readEBSDMonteCarloFile(sinl%masterfile, mcnl, hdferr, EBSDMCdata, getAccume=.TRUE.)
  ! copy the sample tilt angle into the correct variable for writing to the dot product file
  sinl%sig = mcnl%sig

  ! 2. read EBSD master pattern file
  call readEBSDMasterPatternFile(sinl%masterfile, mpnl, hdferr, EBSDMPdata, getmLPNH=.TRUE., getmLPSH=.TRUE.)
  call h5close_EMsoft(hdferr)

  ! 3. build the energy averaged master pattern
  ! ! first make sure that we do some appropriate energy weighting to make the master patterns 2D instead of 3D
  ! ! sum MC counts over rectangle [-1/2,1/2] and [-1,-1/3] in square Lambert, then normalize and use as 
  ! ! weight factors, then store in 2D MPs
  n = size(EBSDMCdata%accum_e, 1) ! get number of energy bins
  allocate(weights(n)) ! allocate space for energy histogram
  do i = 1, n
    weights(i) = sum(EBSDMCdata%accum_e(i,:,:)) ! this could be modified to sum over partial rectangle
  enddo
  weights = weights / sum(weights) ! this is currently wieghted over the full square Lambert

  ! build energy weighted master pattern
  d = mpnl%npx
  allocate(mLPNH(-d:d,-d:d))
  allocate(mLPSH(-d:d,-d:d))
  mLPNH = 0.D0
  mLPSH = 0.D0
  do i = 1, n
    mLPNH = mLPNH + EBSDMPdata%mLPNH(:,:,i) * weights(i)
    mLPSH = mLPSH + EBSDMPdata%mLPSH(:,:,i) * weights(i)
  enddo

! we moved the indexer%init part into the indexing subroutine to facilitate OpenMP implementation

  ! perform the spherical indexing computations
  call EBSDSISubroutine(sinl, mcnl, mpnl, mLPNH, mLPSH, progname, nmldeffile)

end program EMSphInx

!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDSISubroutine
!
!> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
!
!> @brief Master subroutine to control Spherical Indexing
!
!> @param sinl SphInx namelist pointer
!> @param mcnl Monte Carlo data
!> @param mpnl Master Pattern data 
!> @param mLPNH Northern hemisphere master pattern
!> @param mLPSH Southern hemisphere master pattern
!> @param progname name of the program
!> @param nmldeffile namelist filename
!
!> @date 06/11/19 MDG 1.0 original, copied from EMEBSDDI.f90
!> @date 08/03/19 MDG 1.1 read row of patterns instead of single pattern
!> @date 08/04/19 MDG 1.2 reorganizaton of code to enable OpenMP implementation
!> @date 08/05/19 MDG 1.3 added HDF5 output; added optional .ctf and/or .ang output; added IQ map computation
!> @date 08/20/19 MDG 1.4 updated for modified sinl name list
!--------------------------------------------------------------------------
subroutine EBSDSISubroutine(sinl, mcnl, mpnl, mLPNH, mLPSH, progname, nmldeffile)

  use local
  use typedefs
  use NameListTypedefs
  use NameListHandlers
  use sphnamelisthandler
  use io
  use HDFsupport
  use fft_wrap
  use ECPmod, only: GetPointGroup
  use notifications
  ! use timing
  use sphtypedefs
  use sphh5ebsd
  use sphEBSDiomod
  use SphereIndexerMod
  use patternmod
  use filters
  use omp_lib


IMPLICIT NONE

  type     (SphInxNameListType      ),INTENT(INOUT) :: sinl
  type     (MCCLNameListType        ),INTENT(INOUT) :: mcnl
  type     (EBSDMasterNameListType  ),INTENT(INOUT) :: mpnl
  real     (kind=dbl                ),INTENT(INOUT) :: mLPNH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx)
  real     (kind=dbl                ),INTENT(INOUT) :: mLPSH(-mpnl%npx:mpnl%npx,-mpnl%npx:mpnl%npx)
  character(fnlen                   ),INTENT(IN   ) :: progname
  character(fnlen                   ),INTENT(IN   ) :: nmldeffile

  type     (SphereIndexer           )               :: threadindexer 
  integer  (kind=irg                )               :: L, totnumexpt, hdferr, TID, NUMTHREADS, io_int(3), rcnt, ipos
  integer  (kind=irg                )               :: dims(2), dim2(2), dim3(3)
  character(11                      )               :: dstr 
  character(15                      )               :: tstrb, tstre
  character(3                       )               :: vendor
  integer  (kind=irg                )               :: i,j, std=6, kk
  integer  (kind=irg                )               :: pgnum
  ! integer  (kind=irg                ),allocatable   :: indexmain(:,:)
  real     (kind=sgl                )               :: tstart, tstop, io_real(2), ma, mi, pcvals(3)
  character(fnlen                   )               :: xtalname
  integer  (kind=irg                )               :: istart, iend, jstart, jend, ipf_wd, ipf_ht
  logical                                           :: ROIselected 
  integer  (kind=irg                )               :: ipar(4)
  character(fnlen                   ),allocatable   :: MessageLines(:)
  integer  (kind=irg)                               :: NumLines, numsx, numsy
  character(fnlen                   )               :: TitleMessage, exectime, wisdomFile
  character(100                     )               :: c
  real(kind=dbl)                                    :: w, Jres


  real     (kind=dbl                ),allocatable   :: quats(:,:), xcorr(:), pat(:,:), ksqarray(:,:)
  real     (kind=sgl                ),allocatable   :: pat32(:), patrow(:), exptIQ(:), EBSDpat(:,:)
  type     (IdxRes                  )               :: ires
  integer  (kind=irg                )               :: istat, iunitexpt = 41
  integer  (HSIZE_T                 )               :: dims3(3), offset3(3)
  type(C_PTR)                                       :: planf

  call timestamp(datestring=dstr, timestring=tstrb)

  ! convert some parameters from the name list to the regular EMsoft variables
  numsx = sinl%patdims(1)
  numsy = sinl%patdims(2)
  ipf_wd = int(sinl%scandims(1)) 
  ipf_ht = int(sinl%scandims(2))

  pcvals = getEMsoftPCcoordinates(sinl%pctr, sinl%vendor, sinl%delta, numsx, numsy)
  sinl%xpc = pcvals(1)
  sinl%ypc = pcvals(2)
  sinl%L   = pcvals(3)

  ! do we need a region of interest (ROI) ?
  if (sum(sinl%ROImask).ne.0) then
    ROIselected = .TRUE.
    istart = sinl%ROImask(1)
    jstart = sinl%ROImask(2)
    iend   = istart+sinl%ROImask(3)-1
    jend   = jstart+sinl%ROImask(4)-1
  else
    ROIselected = .FALSE.
    istart = 1
    jstart = 1
    iend   = ipf_wd
    jend   = ipf_ht
  end if

  L = numsx*numsy/sinl%binning**2
  if (ROIselected.eqv..TRUE.) then 
      totnumexpt = sinl%ROImask(3)*sinl%ROImask(4)
  else
      totnumexpt = ipf_wd*ipf_ht
  end if

  call Message(' reading from xtalfile '//trim(mcnl%xtalname))
  pgnum = GetPointGroup(mcnl%xtalname)

  !=====================================================
  ! open the pattern file for reading
  !=====================================================
  call h5open_EMsoft(hdferr)

  istat = openExpPatternFile(sinl%patfile, ipf_wd, L, sinl%inputtype, L*4, iunitexpt, sinl%HDFstrings) 
  if (istat.ne.0) then
    call patternmod_errormessage(istat)
    call FatalError("EMSphInx:", "Fatal error handling experimental pattern file")
  end if
  
  if (sinl%flipy.eqv..TRUE.) then 
    call Message(' EBSD patterns will be flipped vertically ')
  else
    call Message(' EBSD patterns will not be flipped vertically ')
  end if
  
  dims3 = (/ numsx, numsy, ipf_wd /)
  offset3 = (/ 0, 0, 0 /)

  !=====================================================
  !=====================================================
  !=====================================================
  ! indexing loop goes here.
  allocate(quats (4,totnumexpt))
  allocate(xcorr (  totnumexpt))

  if (sinl%nthread.eq.0) then 
    call OMP_SET_NUM_THREADS(OMP_GET_MAX_THREADS())
  else
    call OMP_SET_NUM_THREADS(sinl%nthread)
  end if

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(TID, threadindexer, i, j, offset3, patrow, pat32, mi, ma, pat, ires, rcnt, ipos)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()
!$OMP CRITICAL  
  if (TID.eq.0) then 
    io_int(1) = NUMTHREADS
    call WriteValue('  --> number of threads allocated : ',io_int,1)
    rcnt = 1
  end if
!$OMP END CRITICAL  

  allocate(pat   ( numsx ,numsy ))
  allocate(patrow(    L * ipf_wd))
  allocate(pat32 (             L))

!$OMP BARRIER
! this should be possible in a parallel fashion if we remove accessing the EMsoftConfig.json and fftw wisdom files 
! from the DiscreteSHTConstants_Init routine...
! 
  if (TID.eq.0) then 
    io_int(1) = TID
    call WriteValue(' loading wisdom for thread ',io_int,1)
    call FFTWisdom%load()
    call Message(' initializing all other threads ... ',"(A/)")
  end if
!$OMP BARRIER
!$OMP CRITICAL  
    call threadindexer%init(sinl%bw, dble(sinl%sig), dble(sinl%L), dble(sinl%thetac), dble(sinl%delta), numsx, numsy, &
                      mpnl%npx, mLPNH, mLPSH, sinl%circmask)
!$OMP END CRITICAL  

!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC)
  do j = jstart, jend ! loop over rows
    if (TID.eq.0) then 
      io_int(1:3) = (/ rcnt, minval( (/ rcnt+NUMTHREADS-1, jend /) ), jend /)
      call WriteValue(" reading/indexing rows ", io_int, 3,"(I4,' through ',I4,' of',I4)")
    end if

! get a pattern row (slightly faster than reading one pattern at a time); make sure threads wait in line to acces the file
    offset3 = (/ 0, 0, (j-1) * ipf_wd /) ! get hyperslab offset for this pattern row
!$OMP CRITICAL
    call getExpPatternRow(j, ipf_wd, L, L, dims3, offset3, iunitexpt, sinl%inputtype, sinl%HDFstrings, patrow, flipy=sinl%flipy) 
!$OMP END CRITICAL

    do i = istart, iend ! loop over columns
      ! get a pattern from the row
      pat32 = patrow(L*(i-1)+1:L*i)
      
      ! do adaptive histogram equalization if needed
      if(sinl%nregions.gt.1) then
        ma = maxval(pat32)
        mi = minval(pat32)
        pat32 = ( (pat32 - mi) / (ma - mi) ) * 255.0 ! rescale from 0-255
        pat = adhisteq(sinl%nregions, numsx, numsy, nint(pat32)) ! do AHE
      else
        pat = reshape(pat32, (/ numsx, numsy /) ) ! convert from vectorized to 2d (and float to double)
      endif

      ! index the pattern
      ires = threadindexer%index(pat, dble(sinl%xpc), dble(sinl%ypc), sinl%refine) ! for now just use a single pattern center
      if (ROIselected.eqv..TRUE.) then 
        ipos = i + (j-1) * sinl%ROImask(3)
      else
        ipos = i + (j-1) * ipf_wd
      end if 
      quats(:,ipos) = ires%qu
      xcorr(  ipos) = ires%xc
    enddo
    if (TID.eq.0) rcnt = rcnt + NUMTHREADS
  enddo
!$OMP END DO
  deallocate(pat32 )
  deallocate(pat   )
  deallocate(patrow)
  !$OMP BARRIER
  if (TID.eq.0) call FFTWisdom%save()
!$OMP END PARALLEL

  !=====================================================
  !=====================================================
  !=====================================================

  ! perform some timing stuff
  call CPU_TIME(tstop)
  tstop = tstop - tstart
  io_real(1) = float(sinl%nthread * totnumexpt) / tstop
  call WriteValue(' Number of experimental patterns indexed per second : ',io_real,1,"(/,F10.2,/)")

  !=====================================================
  !=====================================================
  !=====================================================

  ! while the pattern file is still open, compute the pattern quality array exptIQ 
  call Message(' Computing pattern quality array ... ',"(A)",advance="no")
  allocate(exptIQ(totnumexpt))

  ! prepare the fftw plan for this pattern size to compute pattern quality (pattern sharpness Q)
  allocate(EBSDPat(numsx,numsy),stat=istat)
  if (istat .ne. 0) stop '  Could not allocate arrays for EBSDPat filter'
  EBSDPat = 0.0
  allocate(ksqarray(numsx,numsy),stat=istat)
  if (istat .ne. 0) stop '  Could not allocate ksqarray array'
  Jres = 0.0
  call init_getEBSDIQ(numsx, numsy, EBSDPat, ksqarray, Jres, planf)
  deallocate(EBSDPat)


!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(TID, i, j, kk, offset3, EBSDpat, patrow, pat, ipos, rcnt)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()
  allocate(EBSDPat(numsx, numsy))
  allocate(patrow(   L * ipf_wd))

!$OMP DO SCHEDULE(DYNAMIC)
  do j = jstart, jend ! loop over rows

! get a pattern row (slightly faster than reading one pattern at a time); make sure threads wait in line to acces the file
    offset3 = (/ 0, 0, (j-1) * ipf_wd /) ! get hyperslab offset for this pattern row
!$OMP CRITICAL
    call getExpPatternRow(j, ipf_wd, L, L, dims3, offset3, iunitexpt, sinl%inputtype, sinl%HDFstrings, patrow) 
!$OMP END CRITICAL

    do i = istart, iend ! loop over columns
      ! get a pattern from the row
      do kk=1,numsy
          EBSDPat(1:numsx,kk) = patrow((i-1)*L+(kk-1)*numsx+1:(i-1)*L+kk*numsx)
      end do
      ipos = i+(j-1)*ipf_wd
      exptIQ(ipos) = sngl(computeEBSDIQ(numsx, numsy, EBSDPat, ksqarray, Jres, planf))
    end do
  enddo
!$OMP END DO
  deallocate(EBSDpat, patrow)
!$OMP END PARALLEL
  call Message('   done.',"(A/)")

  call closeExpPatternFile(sinl%inputtype, iunitexpt)

 ! close the fortran HDF5 interface
  call h5close_EMsoft(hdferr)

  !=====================================================
  !=====================================================
  !=====================================================

  ! ===================
  ! MAIN OUTPUT SECTION
  ! ===================

  ! fill the ipar array with integer parameters that are needed to write the h5ebsd file
  ! (anything other than what is already in the sinl structure)
  ipar = 0
  ipar(1) = totnumexpt
  ipar(2) = pgnum
  if (ROIselected.eqv..TRUE.) then
    ipar(3) = sinl%ROImask(3)
    ipar(4) = sinl%ROImask(4)
  else
    ipar(3) = ipf_wd
    ipar(4) = ipf_ht
  end if 

  ! Initialize FORTRAN interface.
  call h5open_EMsoft(hdferr)

  if (sinl%datafile.ne.'undefined') then 
    vendor = 'TSL'
    call sphh5ebsd_writeFile(vendor, sinl, mcnl%xtalname, dstr, tstrb, ipar, sngl(mcnl%sig), xcorr, quats, &
                             exptIQ, progname, nmldeffile)
    call Message(' Data stored in h5ebsd file : '//trim(sinl%datafile),"(A)")
  end if

! we will need to also write .ang and .ctf files, if requested 
  if (trim(sinl%ctffile).ne.'undefined') then 
    call sphctfebsd_writeFile(sinl, mcnl%xtalname, ipar, sngl(mcnl%EkeV), sngl(mcnl%sig), xcorr, quats, exptIQ)
    call Message(' Data stored in ctf file : '//trim(sinl%ctffile),"(/A)")
  end if
  
  if (trim(sinl%angfile).ne.'undefined') then 
      call sphangebsd_writeFile(sinl,mcnl%xtalname,ipar,xcorr,quats, exptIQ)
      call Message(' Data stored in ang file : '//trim(sinl%angfile),"(/A)")
  end if
  

  ! close the fortran HDF5 interface
  call h5close_EMsoft(hdferr)

  ! if requested, we notify the user that this program has completed its run
  ! if (trim(EMsoft_getNotify()).ne.'Off') then
  !   ! if (trim(sinl%Notify).eq.'On') then ! this don't exist in the spherical namelist !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   if (.false.) then 
  !     NumLines = 3
  !     allocate(MessageLines(NumLines))

  !     call hostnm(c)
   
  !     MessageLines(1) = 'EMEBSDDI program has ended successfully'
  !     MessageLines(2) = 'Indexed data stored in '//trim(sinl%datafile)
  !     write (exectime,"(F15.0)") tstop  
  !     MessageLines(3) = 'Total execution time [s]: '//trim(exectime)
  !     TitleMessage = 'EMsoft on '//trim(c)
  !     i = PostMessage(MessageLines, NumLines, TitleMessage)
  !   end if
  ! end if

end subroutine EBSDSISubroutine
