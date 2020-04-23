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

module mod_DI
  !! author: MDG 
  !! version: 1.0 
  !! date: 03/31/20
  !!
  !! routines for the EMDI program

use mod_kinds
use mod_global
use mod_DIfiles

IMPLICIT NONE 

public :: DIdriver 

!DEC$ ATTRIBUTES DLLEXPORT :: DIdriver

contains

!--------------------------------------------------------------------------
subroutine DIdriver(Cnmldeffile, Cprogname, cproc, ctimeproc, cerrorproc, objAddress, cancel) &
           bind(c, name='DIdriver') 
!! author: MDG 
!! version: 1.0 
!! date: 04/02/20
!!
!! perform the DI computations.
!!
!! this routine must be callable from C++ as well, so the parameter list 
!! is a bit different from that of the most other programs.  Furthermore, this 
!! driver routine must be able to handle EBSD, ECP, and TKD patterns (at least).

use mod_EMsoft
use mod_io
use mod_initializers
use HDF5
use mod_HDFsupport
use mod_patterns
use mod_Lambert
use mod_others
use mod_crystallography
use mod_gvectors
use mod_filters
use mod_diffraction
use mod_symmetry
use mod_quaternions
use mod_rotations
use mod_so3
use mod_math
use clfortran
use mod_CLsupport
use omp_lib
use mod_OMPsupport
use h5im
use h5lt
use ISO_C_BINDING
use mod_notifications
use mod_timing
use mod_MCfiles
use mod_MPfiles
use mod_DIfiles
use mod_DIsupport
use mod_HDFnames
use mod_EBSD
use mod_ECP
use mod_so3
use mod_vendors

IMPLICIT NONE 

! interface for the callback routines
ABSTRACT INTERFACE
   SUBROUTINE ProgCallBackTypeTimingdriver(objAddress, loopCompleted, totalLoops, timeRemaining) bind(C)
    USE, INTRINSIC :: ISO_C_BINDING
    INTEGER(c_size_t),INTENT(IN), VALUE             :: objAddress
    INTEGER(KIND=4), INTENT(IN), VALUE              :: loopCompleted
    INTEGER(KIND=4), INTENT(IN), VALUE              :: totalLoops
    REAL(KIND=4),INTENT(IN), VALUE                  :: timeRemaining
   END SUBROUTINE ProgCallBackTypeTimingdriver

   SUBROUTINE ProgCallBackTypeDIdriver(objAddress, Ndict, euarr_cptr, dparr_cptr, indarr_cptr) bind(C)
    USE, INTRINSIC :: ISO_C_BINDING
    INTEGER(c_size_t),INTENT(IN), VALUE             :: objAddress
    INTEGER(KIND=4), INTENT(IN), VALUE              :: Ndict 
    type(c_ptr), INTENT(OUT)                        :: euarr_cptr
    type(c_ptr), INTENT(OUT)                        :: dparr_cptr
    type(c_ptr), INTENT(OUT)                        :: indarr_cptr
   END SUBROUTINE ProgCallBackTypeDIdriver

   SUBROUTINE ProgCallBackTypeErrorDIdriver(objAddress, errorCode) bind(C)
    USE, INTRINSIC :: ISO_C_BINDING
    INTEGER(c_size_t),INTENT(IN), VALUE             :: objAddress
    INTEGER(KIND=4), INTENT(IN), VALUE              :: errorCode
   END SUBROUTINE ProgCallBackTypeErrorDIdriver

END INTERFACE

character(kind=c_char), INTENT(IN)                  :: Cnmldeffile(fnlen)
character(kind=c_char), INTENT(IN)                  :: Cprogname(fnlen)
TYPE(C_FUNPTR), INTENT(IN), VALUE                   :: cproc
TYPE(C_FUNPTR), INTENT(IN), VALUE                   :: ctimeproc
TYPE(C_FUNPTR), INTENT(IN), VALUE                   :: cerrorproc
integer(c_size_t),INTENT(IN), VALUE                 :: objAddress
character(len=1),INTENT(IN), OPTIONAL               :: cancel

! callback procedure pointer definitions
PROCEDURE(ProgCallBackTypeDIdriver), POINTER        :: proc
PROCEDURE(ProgCallBackTypeTimingdriver), POINTER    :: timeproc
PROCEDURE(ProgCallBackTypeErrorDIdriver), POINTER   :: errorproc
type(c_ptr)                                         :: dparr_cptr=c_null_ptr, euarr_cptr=c_null_ptr, indarr_cptr=c_null_ptr

type(MCfile_T)                                      :: MCFT 
type(MPfile_T)                                      :: MPFT 
type(DIfile_T)                                      :: DIFT 
type(EMsoft_T)                                      :: EMsoft 
type(cell_T)                                        :: cell 
type(HDF_T)                                         :: HDF
type(HDFnames_T)                                    :: HDFnames
type(EBSD_T)                                        :: EBSD
type(ECP_T)                                         :: ECP
type(Timing_T)                                      :: timer 
type(IO_T)                                          :: Message
type(OpenCL_T)                                      :: CL
type(SpaceGroup_T)                                  :: SG
type(so3_T)                                         :: SO 
type(q_T)                                           :: quat
type(e_T)                                           :: eu 
type(r_T)                                           :: ro
type(Vendor_T)                                      :: VT
type(Quaternion_T)                                  :: qu 
type(IncidentListECP),pointer                       :: ktmp

type(MCOpenCLNameListType)                          :: mcnl
type(SEMmasterNameListType)                         :: mpnl
        
logical                                             :: verbose

type(DynType)                                       :: Dyn
type(gnode)                                         :: rlp


integer(c_intptr_t),allocatable, target             :: platform(:)
integer(c_intptr_t),allocatable, target             :: device(:)
integer(c_intptr_t),target                          :: context
integer(c_intptr_t),target                          :: command_queue
integer(c_intptr_t),target                          :: cl_expt,cl_dict
character(len = 50000), target                      :: source
integer(kind=irg), parameter                        :: source_length = 50000
integer(kind=irg), target                           :: source_l
character(len=source_length, KIND=c_char),TARGET    :: csource
type(c_ptr), target                                 :: psource
integer(c_int32_t)                                  :: ierr2, pcnt
integer(c_intptr_t),target                          :: prog
integer(c_intptr_t),target                          :: kernel
integer(c_size_t)                                   :: cnum
character(9),target                                 :: kernelname
character(10, KIND=c_char),target                   :: ckernelname

integer(kind=irg)                                   :: num,ierr,irec,istat, jpar(7), SGnum, nlines
integer(kind=irg),parameter                         :: iunit = 40
integer(kind=irg),parameter                         :: iunitexpt = 41
integer(kind=irg),parameter                         :: iunitdict = 42
character(fnlen)                                    :: info ! info about the GPU
real(kind=dbl),parameter                            :: nAmpere = 6.241D+18   ! Coulomb per second


integer(kind=irg)                                   :: Ne,Nd,L,totnumexpt,numdictsingle,numexptsingle,imght,imgwd,nnk,numE,&
                                                       recordsize, fratio, cratio, fratioE, cratioE, iii, itmpexpt, hdferr, &
                                                       nsig, numk, recordsize_correct, patsz, tickstart, tickstart2, tock, &
                                                       npy, sz(3), jjj
integer(kind=8)                                     :: size_in_bytes_dict,size_in_bytes_expt, Nres
real(kind=sgl),pointer                              :: dict(:), T0dict(:)
real(kind=sgl),allocatable,TARGET                   :: dict1(:), dict2(:), eudictarray(:)
real(kind=sgl),allocatable                          :: imageexpt(:),imagedict(:), mask(:,:),masklin(:), exptIQ(:), &
                                                       exptCI(:), exptFit(:), exppatarray(:), tmpexppatarray(:)
real(kind=sgl),allocatable                          :: imageexptflt(:),binned(:,:),imagedictflt(:),imagedictfltflip(:), &
                                                       tmpimageexpt(:), OSMmap(:,:)
real(kind=sgl),allocatable, target                  :: results(:),expt(:),dicttranspose(:),resultarray(:), dparray(:), &
                                                       eulerarray(:,:),eulerarray2(:,:),resultmain(:,:),resulttmp(:,:)
integer(kind=irg),allocatable                       :: acc_array(:,:), ppend(:), ppendE(:)
integer(kind=irg),allocatable,target                :: indarray(:) 
integer*4,allocatable                               :: iexptCI(:,:), iexptIQ(:,:)
real(kind=sgl),allocatable                          :: meandict(:),meanexpt(:),wf(:),mLPNH(:,:,:),mLPSH(:,:,:),accum_e_MC(:,:,:)
real(kind=sgl),allocatable                          :: mLPNH_simple(:,:), mLPSH_simple(:,:), eangle(:), mLPNH2D(:,:), mLPSH2D(:,:)
real(kind=sgl),allocatable                          :: pattern(:,:), FZarray(:,:), dpmap(:), lstore(:,:), pstore(:,:)
real(kind=sgl),allocatable                          :: patternintd(:,:), lp(:), cp(:), EBSDpat(:,:)
integer(kind=irg),allocatable                       :: patterninteger(:,:), patternad(:,:), EBSDpint(:,:), kij(:,:)
character(kind=c_char),allocatable                  :: EBSDdictpat(:,:,:)
real(kind=sgl),allocatable                          :: dictpatflt(:,:), anglewf(:)
real(kind=dbl),allocatable                          :: rdata(:,:), fdata(:,:), rrdata(:,:), ffdata(:,:), ksqarray(:,:), klist(:,:)
complex(kind=dbl),allocatable                       :: hpmask(:,:)
complex(C_DOUBLE_COMPLEX),allocatable               :: inp(:,:), outp(:,:)
real(kind=dbl)                                      :: w, Jres
integer(kind=irg)                                   :: dims(2)
character(11)                                       :: dstr
character(15)                                       :: tstrb
character(15)                                       :: tstre
character(3)                                        :: vendor
character(fnlen, KIND=c_char),allocatable,TARGET    :: stringarray(:)
character(fnlen)                                    :: groupname, dataset, fname, clname, ename, sourcefile, &
                                                       datagroupname, dictfile, attname
integer(hsize_t)                                    :: expwidth, expheight
integer(hsize_t),allocatable                        :: iPhase(:), iValid(:)
integer(c_size_t),target                            :: slength
integer(c_int)                                      :: numd, nump
type(C_PTR)                                         :: planf, HPplanf, HPplanb
integer(HSIZE_T)                                    :: dims2(2), offset2(2), dims3(3), offset3(3)

integer(kind=irg)                                   :: i,j,ii,jj,kk,ll,mm,pp,qq, cn, dn, totn
integer(kind=irg)                                   :: FZcnt, pgnum, io_int(4), ncubochoric, pc, ecpipar(4)
type(FZpointd),pointer                              :: FZlist, FZtmp
integer(kind=irg),allocatable                       :: indexlist(:),indexarray(:),indexmain(:,:),indextmp(:,:)
real(kind=sgl)                                      :: dmin,voltage,scl,ratio, mi, ma, ratioE, io_real(2), tstart, tmp, &
                                                       totnum_el, vlen, tstop, ttime
real(kind=dbl)                                      :: prefactor
character(fnlen)                                    :: xtalname
integer(kind=irg)                                   :: binx,biny,TID,nthreads,Emin,Emax, iiistart, iiiend, jjend
real(kind=sgl)                                      :: sx,dx,dxm,dy,dym,rhos,x,projweight, dp, mvres, nel, emult
real(kind=sgl)                                      :: dc(3),ixy(2),bindx, MCsig, WD, fpar1(1), fpar2(2)
integer(kind=irg)                                   :: nix,niy,nixp,niyp
real(kind=sgl)                                      :: euler(3)
integer(kind=irg)                                   :: indx
integer(kind=irg)                                   :: correctsize
logical                                             :: f_exists, init, ROIselected, Clinked, cancelled, isTKD = .FALSE., &
                                                       isEBSD = .FALSE., isECP = .FALSE., switchwfoff

integer(kind=irg)                                   :: ipar(10)

character(fnlen),ALLOCATABLE                        :: MessageLines(:)
integer(kind=irg)                                   :: NumLines
character(fnlen)                                    :: TitleMessage, exectime, progname, nmldeffile
character(100)                                      :: c
character(1000)                                     :: charline
character(3)                                        :: stratt
character(fnlen)                                    :: progdesc 

! convert the input strings from C to fortran format
nmldeffile = trim(fstringify(Cnmldeffile))
progname = trim(fstringify(Cprogname))
progdesc = 'Indexing of EBSD/ECP/TKD patterns using a dynamically calculated pattern dictionary'

! open the HDF interface
call openFortranHDFInterface()
HDF = HDF_T() 

! we've already shown the standard splash screen, so we do this one silently
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 80 /), silent=.TRUE.)

! link the proc procedure to the cproc argument
Clinked = .FALSE.
if (present(cancel)) then 
  Clinked = .TRUE.
  nullify(proc, errorproc, timeproc)
  CALL C_F_PROCPOINTER (cproc, proc)
  CALL C_F_PROCPOINTER (ctimeproc, timeproc)
  CALL C_F_PROCPOINTER (cerrorproc, errorproc)
end if

! deal with the namelist stuff
DIFT = DIfile_T(nmldeffile)

! set the HDF group names for this program
HDFnames = HDFnames_T() 

call setRotationPrecision('d')

associate( dinl=>DIFT%nml, MPDT=>MPFT%MPDT, MCDT=>MCFT%MCDT, det=>EBSD%det, enl=>EBSD%nml, ecpnl=>ECP%nml )

! make sure that nthreads is at least 2 
if (dinl%nthreads.lt.2) then 
  call Message%printError('DIdriver:', 'Dictionary Indexing requires at least 2 compute threads')
end if 

! determine the modality from the master pattern file, and also set it in the dinl name list
fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
call MPFT%determineModality(HDF, fname)
call Message%printMessage(' Master Pattern modality : '//trim(MPFT%getModality()))
call DIFT%setModality(MPFT%getModality())

if (trim(MPFT%getModality()).eq.'EBSD') then 
  isEBSD = .TRUE.
else if (trim(MPFT%getModality()).eq.'TKD') then 
  isTKD = .TRUE.
else if (trim(MPFT%getModality()).eq.'ECP') then
  isECP = .TRUE.
  end if  

! is this a dynamic calculation (i.e., do we actually compute the diffraction patterns)?
if (trim(dinl%indexingmode).eq.'dynamic') then 

    ! 1. read the Monte Carlo data file
    call HDFnames%set_ProgramData(SC_MCOpenCL) 
    call HDFnames%set_NMLlist(SC_MCCLNameList) 
    call HDFnames%set_NMLfilename(SC_MCOpenCLNML)
    fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
    call MCFT%setFileName(fname)
    call MCFT%readMCfile(HDF, HDFnames, getAccume=.TRUE.)
    mcnl = MCFT%getnml()
    xtalname = trim(mcnl%xtalname)

    ! 2. read the master pattern file
    if (isTKD.eqv..TRUE.) then
      call HDFnames%set_ProgramData(SC_TKDmaster) 
      call HDFnames%set_NMLlist(SC_TKDmasterNameList) 
      call HDFnames%set_NMLfilename(SC_TKDmasterNML) 
    end if 
    if (isEBSD.eqv..TRUE.) then
      call HDFnames%set_ProgramData(SC_EBSDmaster) 
      call HDFnames%set_NMLlist(SC_EBSDmasterNameList) 
      call HDFnames%set_NMLfilename(SC_EBSDmasterNML) 
    end if 
    if (isECP.eqv..TRUE.) then
      call HDFnames%set_ProgramData(SC_ECPmaster) 
      call HDFnames%set_NMLlist(SC_ECPmasterNameList) 
      call HDFnames%set_NMLfilename(SC_ECPmasterNML) 
    end if 
    call HDFnames%set_Variable(SC_MCOpenCL) 

    fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
    call MPFT%setFileName(fname)
    call MPFT%readMPfile(HDF, HDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

! set the HDFnames for the current program (same for all modalities)
    call HDFnames%set_ProgramData(SC_EMDI) 
    call HDFnames%set_NMLlist(SC_EMDINameList) 
    call HDFnames%set_NMLfilename(SC_EMDI) 

! we know that the master pattern file exists, and it also has all the 
! crystallographic data in it, so we read that here instead of assuming 
! that the actual .xtal file exists on this system ...
    call cell%setFileName(xtalname)
    call cell%readDataHDF(SG, EMsoft, useXtalName=fname)
! extract the point group number 
    pgnum = SG%getPGnumber()
    io_int = pgnum 
    call Message%WriteValue(' Setting point group number to ',io_int,1)

    ! 3. for EBSD/TKD copy a few parameters from dinl to enl
    ! and then generate the detector arrays
    if ( (isEBSD.eqv..TRUE.) .or. (isTKD.eqv..TRUE.)) then 
      allocate(det%rgx(dinl%numsx,dinl%numsy), &
               det%rgy(dinl%numsx,dinl%numsy), &
               det%rgz(dinl%numsx,dinl%numsy), &
               det%accum_e_detector(MCDT%numEbins,dinl%numsx,dinl%numsy), stat=istat)
      enl%numsx = dinl%numsx
      enl%numsy = dinl%numsy
      enl%xpc = dinl%xpc
      enl%ypc = dinl%ypc
      enl%delta = dinl%delta
      enl%thetac = dinl%thetac
      enl%L = dinl%L
      enl%energymin = dinl%energymin
      enl%energymax = dinl%energymax

      if (isTKD.eqv..TRUE.) then 
        call EBSD%GenerateDetector(MCFT, verbose, isTKD)
      end if 
      if (isEBSD.eqv..TRUE.) then 
        call EBSD%GenerateDetector(MCFT, verbose)
      end if 
    else  ! this must be an ECP indexing run so we initialize the appropriate detector arrays
      if (isECP.eqv..TRUE.) then 
        ECP = ECP_T()
      ! copy a few parameters
        ecpnl%conesemiangle = dinl%conesemiangle
        ecpnl%sampletilt = dinl%sampletilt
        ecpnl%npix = dinl%npix 
        ecpnl%workingdistance = dinl%workingdistance
        ecpnl%Rin = dinl%Rin 
        ecpnl%Rout = dinl%Rout

        call ECP%ECPGenerateDetector(verbose=.TRUE.)
        nsig = nint((ecpnl%conesemiangle) + abs(ecpnl%sampletilt)) + 1
        allocate(anglewf(1:nsig),stat=istat)

        call Message%printMessage(' --> Calculating weight factors', frm = "(A)" )
        call ECP%ECPGetWeightFactors(mcnl, MCFT, anglewf, nsig, verbose=.TRUE.)

        !=================================================================
        ! check if there are enough angles in MC for detector geometry
        !=================================================================
        if (mcnl%sigend .lt. (abs(ecpnl%sampletilt) + ecpnl%conesemiangle)) then
          call Message%printMessage(' Not enough angles in Monte carlo file...interpolation will be done without &
          appropriate weight factors',frm = "(A)")
          switchwfoff = .TRUE.
        end if

        if ((-mcnl%sigend .gt. (ecpnl%conesemiangle - abs(ecpnl%sampletilt))) .and. (switchwfoff .eqv. .FALSE.)) then
          call Message%printMessage(' Not enough angles in Monte carlo file...interpolation will be done without &
          appropriate weight factors',frm = "(A)")
          switchwfoff = .TRUE.
        end if

        !=================================================================
        ! generate list of incident vectors
        !=================================================================
        numk = 0
        call ECP%GetVectorsCone()
        numk = ECP%get_numk()
        allocate(kij(2,numk),klist(3,numk),stat=istat)

        io_int(1) = numk
        call Message%WriteValue(' Number of beams for which interpolation will be done = ',io_int,1) 

        ktmp => ECP%get_ListHead()
        ! converting to array for OpenMP parallelization
        do i = 1,numk
           klist(1:3,i) = ktmp%k(1:3)
           kij(1:2,i) = (/ktmp%i,ktmp%j/)
           ktmp => ktmp%next
        end do
        ecpipar(1) = nsig 
        ecpipar(2) = numk 
        ecpipar(3) = ecpnl%npix
        ecpipar(4) = mpnl%npx
      end if 
    end if 

    ! also copy the sample tilt angle into the correct variable for writing to the dot product file
    MCsig = mcnl%sig

    call Message%printMessage(' Completed reading all MC/MP input data; generated detector ')
end if

! set the timer 
timer = Timing_T()
dstr = timer%getDateString()
tstrb = timer%getTimeString()
tstre = ''

if (trim(dinl%indexingmode).eq.'static') then

    ! get the full filename
    if (dinl%dictfile(1:1).ne.EMsoft%getConfigParameter('EMsoftnativedelimiter')) then 
      dictfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(dinl%dictfile)
    else 
      dictfile = trim(dinl%dictfile)
    end if 

    call Message%printMessage('-->  '//'Opening HDF5 dictionary file '//trim(dinl%dictfile))

    hdferr =  HDF%openFile(dictfile)
    if (hdferr.ne.0) call HDF%error_check('HDF_openFile ', hdferr)

    ! we need the point group number (derived from the space group number)
    ! if MPDT%newSGnumber is set to 2, then pgnum must be set to 1 for 
    ! overlap master patterns  [ added by MDG, 06/19/19 ]
    if (MPDT%AveragedMP.eqv..TRUE.) then 
        pgnum = MPDT%newPGnumber
        io_int = pgnum 
        call Message%WriteValue(' Setting point group number to ',io_int,1)
    else
        groupname = SC_CrystalData
        hdferr = HDF%openGroup(groupname)
        if (hdferr.ne.0) call HDF%error_check('HDF_openGroup:CrystalData', hdferr)

        dataset = SC_SpaceGroupNumber
        call HDF%readDatasetInteger(dataset, hdferr, SGnum)
        if (hdferr.ne.0) call HDF%error_check('HDF_readDatasetInteger:SpaceGroupNumber', hdferr)
        call HDF%pop()
    ! get the point group number    
        if (SGnum.ge.221) then
          pgnum = 32
        else
          i=0
          do while (SGPG(i+1).le.SGnum) 
            i = i+1
          end do
          pgnum = i
        end if
        io_int = pgnum 
        call Message%WriteValue(' Setting point group number to ',io_int,1)
    end if

    ! then read some more data from the EMData group
    hdferr = HDF%openGroup(HDFnames%get_EMData())
    if (hdferr.ne.0) call HDF%error_check('HDF_openGroup:EMData', hdferr)

    datagroupname = 'EBSD'
    hdferr = HDF%openGroup(datagroupname)
    if (hdferr.ne.0) call HDF%error_check('HDF_openGroup:EBSD', hdferr)

! test the HDF_FileVersion to make sure that the dictionary file is recent enough
    attname = 'HDF_FileVersion'
    hdferr = HDF%getStringAttributeFromGroup(attname, stratt, 3_SIZE_T)

    if (stratt.eq.'4.0') then
        call Message%printMessage('The dictionary file was created with an older version of the EMEBSD program.')
        call Message%printMessage('This file can not be used by the present program; must be version 4.1 or higher.')
        call Message%printMessage('')
        call Message%printError('DIdriver','Incompatible dictionary file; please rerun the EMEBSD program.')
    end if

    ! we already have the xtalname string from the Monte Carlo name list
    xtalname = trim(mcnl%xtalname)

    ! number of Eulerangles numangles
    dataset = SC_numangles
    call HDF%readDatasetInteger(dataset, hdferr, FZcnt)
    if (hdferr.ne.0) call HDF%error_check('HDF_readDatasetInteger:numangles', hdferr)

    ! euler angle list Eulerangles
    dataset = SC_Eulerangles
    call HDF%readDatasetFloatArray(dataset, dims2, hdferr, eulerarray2)
    eulerarray2 = eulerarray2 * rtod
    if (hdferr.ne.0) call HDF%error_check('HDF_readDatasetFloatArray2D:Eulerangles', hdferr)

    ! we leave this file open since we still need to read all the patterns...
    !=====================================================
    call Message%printMessage('-->  completed initial reading of dictionary file ')
end if

if (sum(dinl%ROI).ne.0) then
  ROIselected = .TRUE.
  iiistart = dinl%ROI(2)
  iiiend = dinl%ROI(2)+dinl%ROI(4)-1
  jjend = dinl%ROI(3)
else
  ROIselected = .FALSE.
  iiistart = 1
  iiiend = dinl%ipf_ht
  jjend = dinl%ipf_wd
end if

verbose = .FALSE.
init = .TRUE.
Ne = dinl%numexptsingle
Nd = dinl%numdictsingle
L = dinl%numsx*dinl%numsy/dinl%binning**2
if (ROIselected.eqv..TRUE.) then 
    totnumexpt = dinl%ROI(3)*dinl%ROI(4)
else
    totnumexpt = dinl%ipf_wd*dinl%ipf_ht
end if
imght = dinl%numsx/dinl%binning
imgwd = dinl%numsy/dinl%binning
dims = (/imght, imgwd/)
nnk = dinl%nnk
ncubochoric = dinl%ncubochoric
recordsize = L*4
itmpexpt = 43
w = dinl%hipassw
source_l = source_length

! these will eventually need to be read from an experimental data file but we'll set default values here.
WD = 10.0

! nullify the dict and T0dict pointers
nullify(dict,T0dict)

! make sure that correctsize is a multiple of 16; if not, make it so
if (mod(L,16) .ne. 0) then
    correctsize = 16*ceiling(float(L)/16.0)
else
    correctsize = L
end if

! determine the experimental and dictionary sizes in bytes
size_in_bytes_dict = Nd*correctsize*sizeof(correctsize)
size_in_bytes_expt = Ne*correctsize*sizeof(correctsize)
recordsize_correct = correctsize*4
patsz              = correctsize

! do a quick sanity check for the requested GPU memory 
call Message%printMessage(' --> Initializing OpenCL device')
CL = OpenCL_T()
Nres = Ne*Nd*4
call CL%query_platform_info(dinl%platid)
call CL%DI_memory_estimate(Nres, size_in_bytes_dict, size_in_bytes_expt, dinl%platid, dinl%devid)

if (trim(dinl%indexingmode).eq.'dynamic') then 
! override the point group number if this is an overlap master pattern
    if (MPDT%AveragedMP.eqv..TRUE.) then 
        pgnum = MPDT%newPGnumber
    end if

    !=====================================================
    ! make sure the minimum energy is set smaller than the maximum
    !=====================================================
    if (dinl%energymin.gt.dinl%energymax) then
        call Message%printMessage(' Minimum energy is larger than maximum energy; please correct input file')
        stop
    end if

    !=====================================================
    ! get the indices of the minimum and maximum energy
    !=====================================================
    Emin = nint((dinl%energymin - mcnl%Ehistmin)/mcnl%Ebinsize) +1
    if (Emin.lt.1)  Emin=1
    if (Emin.gt.MCDT%numEbins)  Emin=MCDT%numEbins

    Emax = nint((dinl%energymax - mcnl%Ehistmin)/mcnl%Ebinsize) + 1
    if (Emax .lt. 1) Emax = 1
    if (Emax .gt. MCDT%numEbins) Emax = MCDT%numEbins

    sz = shape(MPDT%mLPNH)
    numE = sz(3)

    ! intensity prefactor
    nel = float(mcnl%totnum_el) * float(mcnl%multiplier)
    emult = nAmpere * 1e-9 / nel  ! multiplicative factor to convert MC data to an equivalent incident beam of 1 nanoCoulomb
    ! intensity prefactor  (redefined by MDG, 3/23/18)
    prefactor = emult * dinl%beamcurrent * dinl%dwelltime * 1.0D-6
end if

!====================================
! init a bunch of parameters
!====================================
! binned pattern array
binx = dinl%numsx/dinl%binning
biny = dinl%numsy/dinl%binning
bindx = 1.0/float(dinl%binning)**2

! allocate the square-Lambert arrays
npy = mpnl%npx
if (trim(dinl%indexingmode).eq.'dynamic') then
  if (isECP.eqv..TRUE.) then 
    allocate(mLPNH2D(-mpnl%npx:mpnl%npx,-npy:npy))
    allocate(mLPSH2D(-mpnl%npx:mpnl%npx,-npy:npy))
    mLPNH2D = sum(MPDT%mLPNH,3)
    mLPSH2D = sum(MPDT%mLPSH,3)
  else
    allocate(mLPNH(-mpnl%npx:mpnl%npx,-npy:npy,MCDT%numEbins))
    allocate(mLPSH(-mpnl%npx:mpnl%npx,-npy:npy,MCDT%numEbins))
    allocate(accum_e_MC(MCDT%numEbins,dinl%numsx,dinl%numsy),stat=istat)
    accum_e_MC = det%accum_e_detector
    mLPNH = MPDT%mLPNH
    mLPSH = MPDT%mLPSH
  end if 
end if

!=====================================================
! SAMPLING OF RODRIGUES FUNDAMENTAL ZONE
!=====================================================
! if eulerfile is not defined, then we use the standard RFZ sampling;
! if it is defined, then we read the Eulerangle triplets from the file
! and generate the FZlist here... this can be useful to index patterns that
! have only a small misorientation range with respect to a known orientation,
! so that it is not necessary to scan all of orientation space.
if (trim(dinl%indexingmode).eq.'dynamic') then
    SO = so3_T(pgnum, zerolist='FZ') 

    if (trim(dinl%eulerfile).eq.'undefined') then
      call Message%printMessage(' Orientation space sampling mode set to RFZ')
      io_int(1) = pgnum
      io_int(2) = ncubochoric
      call Message%WriteValue(' Point group number and number of cubochoric sampling points : ',io_int,2,"(I4,',',I5)")

      call SO%sampleRFZ(ncubochoric)
      FZcnt = SO%getListCount('FZ')

      if (Clinked.eqv..TRUE.) then 
! generate the Euler dictionary array needed by the EMsoftWorkbench
        allocate(eudictarray(3*FZcnt))
        FZtmp => SO%getListHead('FZ')
        do ii = 1,FZcnt
          eu = FZtmp%rod%re()
          eudictarray((ii-1)*3+1:(ii-1)*3+3) = eu%e_copyd()
          FZtmp => FZtmp%next
        end do 
      end if
    else
    ! read the euler angle file and create the linked list
      call SO%getOrientationsfromFile(dinl%eulerfile)
      FZcnt = SO%getListCount('FZ')
      call Message%printMessage(' Orientation space sampling mode set to MIS')
      io_int(1) = pgnum
      io_int(2) = FZcnt
      call Message%WriteValue(' Point group number and number of sampling points : ',io_int,2,"(I4,',',I5)")
    end if

    ! allocate and fill FZarray for OpenMP parallelization
    allocate(FZarray(4,FZcnt),stat=istat)
    FZarray = 0.0

    FZtmp => SO%getListHead('FZ')
    do ii = 1,FZcnt
        FZarray(1:4,ii) = FZtmp%rod%r_copyd()
        FZtmp => FZtmp%next
    end do
    io_int(1) = FZcnt
    call Message%WriteValue(' Number of unique orientations sampled :        : ', io_int, 1, "(I8)")
! we can now delete the linked list since we have the FZarray
    call SO%delete_FZlist()
end if 


!================================
! INITIALIZATION OF OpenCL DEVICE
!================================

call CL%init_PDCCQ(platform, nump, dinl%platid, device, numd, dinl%devid, info, context, command_queue)

! read the cl source file
sourcefile = 'DictIndx.cl'
call CL%read_source_file(EMsoft, sourcefile, csource, slength)

! allocate device memory for experimental and dictionary patterns
cl_expt = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_expt, C_NULL_PTR, ierr)
call CL%error_check('DIdriver:clCreateBuffer', ierr)

cl_dict = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dict, C_NULL_PTR, ierr)
call CL%error_check('DIdriver:clCreateBuffer', ierr)

!================================
! the following lines were originally in the InnerProdGPU routine, but there is no need
! to execute them each time that routine is called so we move them here...
!================================
! create the program
pcnt = 1
psource = C_LOC(csource)
prog = clCreateProgramWithSource(context, pcnt, C_LOC(psource), C_LOC(slength), ierr)
call CL%error_check('InnerProdGPU:clCreateProgramWithSource', ierr)

! build the program
ierr = clBuildProgram(prog, numd, C_LOC(device), C_NULL_PTR, C_NULL_FUNPTR, C_NULL_PTR)

! get the compilation log
ierr2 = clGetProgramBuildInfo(prog, device(dinl%devid), CL_PROGRAM_BUILD_LOG, sizeof(source), C_LOC(source), cnum)
if (len(trim(source)) > 0) call Message%printMessage(trim(source(1:cnum)),frm='(A)')
call CL%error_check('InnerProdGPU:clBuildProgram', ierr)
call CL%error_check('InnerProdGPU:clGetProgramBuildInfo', ierr2)

call Message%printMessage(' Program Build Successful... Creating kernel')

! finally get the kernel and release the program
kernelname = 'InnerProd'
ckernelname = kernelname
ckernelname(10:10) = C_NULL_CHAR
kernel = clCreateKernel(prog, C_LOC(ckernelname), ierr)
call CL%error_check('InnerProdGPU:clCreateKernel', ierr)

ierr = clReleaseProgram(prog)
call CL%error_check('InnerProdGPU:clReleaseProgram', ierr)

! the remainder is done in the InnerProdGPU routine
!=========================================

!=========================================
! ALLOCATION AND INITIALIZATION OF ARRAYS
!=========================================
call Message%printMessage(' --> Allocating various arrays for indexing')

allocate(expt(Ne*correctsize),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for experimental patterns'
expt = 0.0

allocate(dict1(Nd*correctsize),dict2(Nd*correctsize),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for dictionary patterns'
dict1 = 0.0
dict2 = 0.0
dict => dict1
dicttranspose = 0.0

allocate(results(Ne*Nd),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for results'
results = 0.0

allocate(mask(binx,biny),masklin(L),stat=istat)
if (istat .ne. 0) stop 'Could not allocate arrays for masks'
mask = 1.0
masklin = 0.0

allocate(imageexpt(L),imageexptflt(correctsize),stat=istat)
allocate(tmpimageexpt(correctsize),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for reading experimental image patterns'
imageexpt = 0.0
imageexptflt = 0.0

allocate(meandict(correctsize),meanexpt(correctsize),imagedict(correctsize),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for mean dictionary and experimental patterns'
meandict = 0.0
meanexpt = 0.0

! allocate(pattern(dinl%numsx,dinl%numsy),binned(binx,biny),stat=istat)
allocate(pattern(binx,biny),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for EBSD pattern'
pattern = 0.0

allocate(resultarray(1:Nd),stat=istat)
if (istat .ne. 0) stop 'could not allocate result arrays'
resultarray = 0.0

allocate(indexarray(1:Nd),stat=istat)
if (istat .ne. 0) stop 'could not allocate index arrays'
indexarray = 0

allocate(indexlist(1:Nd*(ceiling(float(FZcnt)/float(Nd)))),stat=istat)
if (istat .ne. 0) stop 'could not allocate indexlist arrays'
do ii = 1,Nd*ceiling(float(FZcnt)/float(Nd))
    indexlist(ii) = ii
end do

allocate(resultmain(nnk,Ne*ceiling(float(totnumexpt)/float(Ne))),stat=istat)
if (istat .ne. 0) stop 'could not allocate main result array'
resultmain = -2.0

allocate(indexmain(nnk,Ne*ceiling(float(totnumexpt)/float(Ne))),stat=istat)
if (istat .ne. 0) stop 'could not allocate main index array'
indexmain = 0

allocate(resulttmp(2*nnk,Ne*ceiling(float(totnumexpt)/float(Ne))),stat=istat)
if (istat .ne. 0) stop 'could not allocate temporary result array'
resulttmp = -2.0

allocate(indextmp(2*nnk,Ne*ceiling(float(totnumexpt)/float(Ne))),stat=istat)
if (istat .ne. 0) stop 'could not allocate temporary index array'
indextmp = 0

allocate(eulerarray(1:3,Nd*ceiling(float(FZcnt)/float(Nd))),stat=istat)
if (istat .ne. 0) stop 'could not allocate euler array'
eulerarray = 0.0
if (trim(dinl%indexingmode).eq.'static') then
    eulerarray(1:3,1:FZcnt) = eulerarray2(1:3,1:FZcnt)
    deallocate(eulerarray2)
end if

allocate(exptIQ(totnumexpt), exptCI(totnumexpt), exptFit(totnumexpt), stat=istat)
if (istat .ne. 0) stop 'could not allocate exptIQ array'

allocate(rdata(binx,biny),fdata(binx,biny),stat=istat)
if (istat .ne. 0) stop 'could not allocate arrays for Hi-Pass filter'
rdata = 0.D0
fdata = 0.D0

!=====================================================
! determine loop variables to avoid having to duplicate 
! large sections of mostly identical code
!=====================================================
ratio = float(FZcnt)/float(Nd)
cratio = ceiling(ratio)
fratio = floor(ratio)

ratioE = float(totnumexpt)/float(Ne)
cratioE = ceiling(ratioE)
fratioE = floor(ratioE)

allocate(ppend(cratio),ppendE(cratioE))
ppend = (/ (Nd, i=1,cratio) /)
if (fratio.lt.cratio) then
  ppend(cratio) = MODULO(FZcnt,Nd)
end if

ppendE = (/ (Ne, i=1,cratioE) /)
if (fratioE.lt.cratioE) then
  ppendE(cratioE) = MODULO(totnumexpt,Ne)
end if

!=====================================================
! define the circular mask if necessary and convert to 1D vector
!=====================================================

if (trim(dinl%maskfile).ne.'undefined') then
! read the mask from file; the mask can be defined by a 2D array of 0 and 1 values
! that is stored in row form as strings, e.g.    
!    0000001110000000
!    0000011111000000
! ... etc
!
    f_exists = .FALSE.
    if (dinl%maskfile(1:1).ne.EMsoft%getConfigParameter('EMsoftnativedelimiter')) then
      fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(dinl%maskfile)
    else
      fname = trim(dinl%maskfile)
    end if 
    inquire(file=trim(fname), exist=f_exists)
    if (f_exists.eqv..TRUE.) then
      mask = 0.0
      open(unit=dataunit,file=trim(fname),status='old',form='formatted')
      do jj=biny,1,-1
        read(dataunit,"(A)") charline
        do ii=1,binx
          if (charline(ii:ii).eq.'1') mask(ii,jj) = 1.0
        end do
      end do
      close(unit=dataunit,status='keep')
    else
      call Message%printError('DIdriver',' maskfile '//trim(fname)//' does not exist')
    end if
else
    if (dinl%maskpattern.eq.'y') then
      do ii = 1,biny
          do jj = 1,binx
              if((ii-biny/2)**2 + (jj-binx/2)**2 .ge. dinl%maskradius**2) then
                  mask(jj,ii) = 0.0
              end if
          end do
      end do
    end if
end if

! convert the mask to a linear (1D) array
do ii = 1,biny
    do jj = 1,binx
        masklin((ii-1)*binx+jj) = mask(jj,ii)
    end do
end do

!=====================================================
! Preprocess all the experimental patterns and store
! them in a temporary file as vectors; also, create 
! an average dot product map to be stored in the h5ebsd output file
!=====================================================
call PreProcessPatterns(EMsoft, HDF, .FALSE., dinl, binx, biny, masklin, correctsize, totnumexpt, exptIQ=exptIQ)

!=====================================================
call Message%printMessage(' --> computing Average Dot Product map (ADP)')
call Message%printMessage(' ')

! re-open the temporary file
if (dinl%tmpfile(1:1).ne.EMsoft%getConfigParameter('EMsoftnativedelimiter')) then 
  fname = trim(EMsoft%generateFilePath('EMtmppathname'))//trim(dinl%tmpfile)
else
  fname = trim(dinl%tmpfile)
end if

open(unit=itmpexpt,file=trim(fname),&
     status='old',form='unformatted',access='direct',recl=recordsize_correct,iostat=ierr)

! use the getADPmap routine in the filters module
if (ROIselected.eqv..TRUE.) then
  allocate(dpmap(dinl%ROI(3)*dinl%ROI(4)))
  call getADPmap(itmpexpt, dinl%ROI(3)*dinl%ROI(4), L, dinl%ROI(3), dinl%ROI(4), dpmap)
else
  allocate(dpmap(totnumexpt))
  call getADPmap(itmpexpt, totnumexpt, L, dinl%ipf_wd, dinl%ipf_ht, dpmap)
end if

! we will leave the itmpexpt file open, since we'll be reading from it again...

!=====================================================
! MAIN COMPUTATIONAL LOOP (finally...)
!
! Some explanation is necessary here... the bulk of this code is 
! executed in OpenMP multithreaded mode, with nthreads threads.
! Thread 0 has a special role described below; threads 1 ... nthreads-1
! share the computation of the dictionary patterns, and wait for 
! thread 0 to finish, if necessary.
!
! Thread 0 takes the dictionary patterns computed by the other threads
! in the previous step in the dictionaryloop and sends them to the GPU,
! along with as many chunks of experimental data are to be handled (experimentalloop
! inside the thread 0 portion of the code); the experimental patterns 
! are then read from the temporary file (unit itmpexpt).  Once all dot
! products have been computed by the GPU, thread 0 will rank them largest
! to smallest and keep only the top nnk values along with their indices 
! into the array of Euler angle triplets.  If the other threads are still
! computing dictionary patterns, thread 0 will join them; otherwise
! thread 0 will immediately take the next batch of dictionary patterns 
! and start all over.
!
! The trick is for the user to determine the array chunk sizes so that 
! threads 1 ... nthreads-1 do not have to wait long for thread 0 to finish;
! this requires a bit of experimenting and observing the load on all the 
! system cores.  The load should always be approximately 100% x nthreads-1
! for an efficient execution.  The appropriate number of threads will depend
! on how powerful the GPU card is...
!=====================================================

call timer%makeTimeStamp() 
call timer%Time_tick(1)
call timer%Time_tick(2)

if (trim(dinl%indexingmode).eq.'dynamic') then
  call OMP_setNThreads(dinl%nthreads)
else
  call OMP_setNThreads(2)
end if

! define the jpar array of integer parameters
jpar(1) = dinl%binning
jpar(2) = dinl%numsx
jpar(3) = dinl%numsy
jpar(4) = mpnl%npx
jpar(5) = npy
jpar(6) = MCDT%numEbins
jpar(7) = numE

! do we need to allocate arrays for the cproc callback routine ?
if (Clinked.eqv..TRUE.) then 
  allocate(dparray(totnumexpt), indarray(totnumexpt))
! and get the C_LOC pointers to those arrays 
  dparr_cptr = C_LOC(dparray)
  indarr_cptr = C_LOC(indarray)
  euarr_cptr = C_LOC(eudictarray)
! and set the callback counters
  totn = cratio+1 
  dn = 1
  cn = 1 
  cancelled = .FALSE.
end if 

dictionaryloop: do ii = 1,cratio+1
    results = 0.0

! if ii is odd, then we use dict1 for the dictionary computation, and dict2 for the GPU
! (assuming ii>1); when ii is even we switch the two pointers 
    if (mod(ii,2).eq.1) then
      dict => dict1
      dict1 = 0.0
      T0dict => dict2   ! these are the patterns to be sent to the GPU
      if (verbose.eqv..TRUE.) call Message%WriteValue('','dict => dict1; T0dict => dict2')
    else
      dict => dict2
      dict2 = 0.0
      T0dict => dict1   ! these are the patterns to be sent to the GPU
      if (verbose.eqv..TRUE.) call Message%WriteValue('','dict => dict2; T0dict => dict1')
    end if

    if (verbose.eqv..TRUE.) then
      io_int(1) = ii
      call Message%WriteValue(' Dictionaryloop index = ',io_int,1)
    end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,iii,jj,ll,mm,pp,ierr,io_int, vlen, tock, ttime, dicttranspose, dictpatflt) &
!$OMP& PRIVATE(binned, ma, mi, patternintd, patterninteger, patternad, qu, ro, quat, imagedictflt,imagedictfltflip)

! allocate the local arrays that are used by each thread
    allocate(patterninteger(binx,biny))
    patterninteger = 0
    allocate(patternad(binx,biny),patternintd(binx,biny))
    patternad = 0.0
    patternintd = 0.0
    allocate(imagedictflt(correctsize),imagedictfltflip(correctsize))

    TID = OMP_GET_THREAD_NUM()

    if ((ii.eq.1).and.(TID.eq.0)) then 
      io_int(1) = OMP_GET_NUM_THREADS()
      call Message%WriteValue(' actual number of OpenMP threads  = ', io_int, 1)
    end if 

! the master thread should be the one working on the GPU computation
!$OMP MASTER
    if (ii.gt.1) then
      iii = ii-1        ! the index ii is already one ahead, since the GPU thread lags one cycle behind the others...
      if (verbose.eqv..TRUE.) then 
        if (associated(T0dict,dict1)) then 
          call Message%printMessage('   GPU thread is working on dict1')
        else
          call Message%printMessage('   GPU thread is working on dict2')
        end if
      end if

      allocate(dicttranspose(Nd*correctsize))
      dicttranspose = 0.0

      do ll = 1,correctsize
        do mm = 1,Nd
            dicttranspose((ll-1)*Nd+mm) = T0dict((mm-1)*correctsize+ll)
        end do
      end do

      ierr = clEnqueueWriteBuffer(command_queue, cl_dict, CL_TRUE, 0_8, size_in_bytes_dict, C_LOC(dicttranspose(1)), &
                                  0, C_NULL_PTR, C_NULL_PTR)
      call CL%error_check('DIdriver:clEnqueueWriteBuffer', ierr)

      mvres = 0.0

      experimentalloop: do jj = 1,cratioE

        expt = 0.0

        do pp = 1,ppendE(jj)   ! Ne or MODULO(totnumexpt,Ne)
          read(itmpexpt,rec=(jj-1)*Ne+pp) tmpimageexpt
          expt((pp-1)*correctsize+1:pp*correctsize) = tmpimageexpt
        end do

        ierr = clEnqueueWriteBuffer(command_queue, cl_expt, CL_TRUE, 0_8, size_in_bytes_expt, C_LOC(expt(1)), &
                                    0, C_NULL_PTR, C_NULL_PTR)
        call CL%error_check('DIdriver:clEnqueueWriteBuffer', ierr)

        call InnerProdGPU(CL,cl_expt,cl_dict,Ne,Nd,correctsize,results,numd,dinl%devid,kernel,context,command_queue)

        dp =  maxval(results)
        if (dp.gt.mvres) mvres = dp

! this might be simplified later for the remainder of the patterns
        do qq = 1,ppendE(jj)
            jjj = (jj-1)*Ne+qq
            resultarray(1:Nd) = results((qq-1)*Nd+1:qq*Nd)
            indexarray(1:Nd) = indexlist((iii-1)*Nd+1:iii*Nd)

            call SSORT(resultarray,indexarray,Nd,-2)
            resulttmp(nnk+1:2*nnk,jjj) = resultarray(1:nnk)
            indextmp(nnk+1:2*nnk,jjj) = indexarray(1:nnk)

            call SSORT(resulttmp(:,jjj),indextmp(:,jjj),2*nnk,-2)

            resultmain(1:nnk,jjj) = resulttmp(1:nnk,jjj)
            indexmain(1:nnk,jjj) = indextmp(1:nnk,jjj)
        end do

! handle the callback routines if requested 
        if (Clinked.eqv..TRUE.) then 
! has the cancel flag been set by the calling program ?
          if (cancel.ne.char(0)) cancelled = .TRUE.
! extract the first row from the indexmain and resultmain arrays, put them in 
! 1D arrays, and return the C-pointer to those arrays via the cproc callback routine 
          dparray(1:totnumexpt) = resultmain(1,1:totnumexpt) 
          indarray(1:totnumexpt) = indexmain(1,1:totnumexpt)
! and call the callback routine ... 
! callback arguments:  objAddress, loopCompleted, totalLoops, timeRemaining, dparray, indarray
          call proc(objAddress, FZcnt, euarr_cptr, dparr_cptr, indarr_cptr)
        end if

      end do experimentalloop

      deallocate(dicttranspose)

      io_real(1) = mvres
      io_real(2) = float(iii)/float(cratio)*100.0
      call Message%WriteValue('',io_real,2,"(' max. dot product = ',F10.6,';',F6.1,'% complete')")

      if (.not.Clinked) then
        if (mod(iii,10) .eq. 0) then
! do a remaining time estimate
! and print information
          if (iii.eq.10) then
              call timer%Time_tock(1)
              tock = timer%getInterval(1)
              ttime = float(tock) * float(cratio) / float(iii)
              tstop = ttime
              io_int(1:4) = (/iii,cratio, int(ttime/3600.0), int(mod(ttime,3600.0)/60.0)/)
              call Message%WriteValue('',io_int,4,"(' -> Completed cycle ',I5,' out of ',I5,'; est. total time ', &
                                      I4,' hrs',I3,' min')")
          else
              ttime = tstop * float(cratio-iii) / float(cratio)
              io_int(1:4) = (/iii,cratio, int(ttime/3600.0), int(mod(ttime,3600.0)/60.0)/)
              call Message%WriteValue('',io_int,4,"(' -> Completed cycle ',I5,' out of ',I5,'; est. remaining time ', &
                                      I4,' hrs',I3,' min')")
          end if
        end if
      end if
    else
       if (verbose.eqv..TRUE.) call Message%WriteValue('','        GPU thread is idling')
    end if  ! ii.gt.1

    if (Clinked.eqv..TRUE.) then 
! has the cancel flag been set by the calling program ?
      if (cancel.ne.char(0)) cancelled = .TRUE.
      ! get the timer value 
      if (iii.lt.5) then 
        ttime = 0.0
      else 
        if (iii.eq.5) then 
          call timer%Time_tock(1)
          tock = timer%getInterval(1)
          ttime = float(tock) * float(cratio) / float(iii)
          tstop = ttime
        else 
          ttime = tstop * float(cratio-iii) / float(cratio)
        end if 
      end if 
      call timeproc(objAddress, cn, totn, ttime) 
      cn = cn + dn
    end if
!$OMP END MASTER


! here we carry out the dictionary pattern computation, unless we are in the ii=cratio+1 step
    if (ii.lt.cratio+1) then
     if (verbose.eqv..TRUE.) then
       if (associated(dict,dict1)) then 
         io_int(1) = TID 
         call Message%WriteValue('    Thread ',io_int,1,"(I5,' is working on dict1')")
       else
         call Message%WriteValue('    Thread ',io_int,1,"(I5,' is working on dict2')")
       end if
     end if

     if (trim(dinl%indexingmode).eq.'dynamic') then
      allocate(binned(binx,biny))
!$OMP DO SCHEDULE(DYNAMIC)

      do pp = 1,ppend(ii)  !Nd or MODULO(FZcnt,Nd)
       if (cancelled.eqv..FALSE.) then
         binned = 0.0
         ro = r_T( rdinp = dble(FZarray(1:4,(ii-1)*Nd+pp)) )
         quat = ro%rq()
         qu = Quaternion_T( qd = quat%q_copyd() )

         if ( (isEBSD.eqv..TRUE.) .or. (isTKD.eqv..TRUE.) ) then 
           call EBSD%CalcEBSDPatternSingleFull(jpar,qu,accum_e_MC,mLPNH,mLPSH,det%rgx,&
                                               det%rgy,det%rgz,binned,Emin,Emax,mask,prefactor)
         else  ! ECP modality 
           call ECP%CalcECPatternSingle(ecpipar, qu, anglewf, mLPNH2D, mLPSH2D, kij, klist, binned, .FALSE.)
         end if 

         if (dinl%scalingmode .eq. 'gam') then
           binned = binned**dinl%gammavalue
         end if

! hi pass filtering
!      rdata = dble(binned)
!      fdata = HiPassFilter(rdata,dims,w)
!      binned = sngl(fdata)

! adaptive histogram equalization
         ma = maxval(binned)
         mi = minval(binned)
         
         patternintd = ((binned - mi)/ (ma-mi))
         patterninteger = nint(patternintd*255.0)
         patternad =  adhisteq(dinl%nregions,binx,biny,patterninteger)
         binned = float(patternad)

         imagedictflt = 0.0
         imagedictfltflip = 0.0
         do ll = 1,biny
           do mm = 1,binx
             imagedictflt((ll-1)*binx+mm) = binned(mm,ll)
           end do
         end do

! normalize and apply circular mask 
         imagedictflt(1:L) = imagedictflt(1:L) * masklin(1:L)
         vlen = vecnorm(imagedictflt(1:correctsize))
         if (vlen.ne.0.0) then
           imagedictflt(1:correctsize) = imagedictflt(1:correctsize)/vlen
         else
           imagedictflt(1:correctsize) = 0.0
         end if
         
         dict((pp-1)*correctsize+1:pp*correctsize) = imagedictflt(1:correctsize)
         ro = r_T( rdinp = dble(FZarray(1:4,(ii-1)*Nd+pp)) ) 
         eu = ro%re()
         eulerarray(1:3,(ii-1)*Nd+pp) = rtod * eu%e_copyd()
       end if
      end do
!$OMP END DO
      deallocate(binned)
    else  ! we are doing static indexing, so only 2 threads in total

! get a set of patterns from the precomputed dictionary file... 
! we'll use a hyperslab to read a block of preprocessed patterns from file 

      if (TID .ne. 0) then
! read data from the hyperslab
       dataset = SC_patterns
       dims2 = (/ correctsize, ppend(ii) /)
       offset2 = (/ 0, (ii-1)*Nd /)

       if(allocated(dictpatflt)) deallocate(dictpatflt)
       dictpatflt = HDF%readHyperslabFloatArray2D(dataset, offset2, dims2)
      
       do pp = 1,ppend(ii)  !Nd or MODULO(FZcnt,Nd)
         dict((pp-1)*correctsize+1:pp*correctsize) = dictpatflt(1:correctsize,pp)
       end do
     end if   
    end if

    if (verbose.eqv..TRUE.) then
       io_int(1) = TID
       call Message%WriteValue('',io_int,1,"('       Thread ',I2,' is done')")
    end if
   else
    if (verbose.eqv..TRUE.) then
       io_int(1) = TID
       call Message%WriteValue('',io_int,1,"('       Thread ',I2,' idling')")
    end if
   end if

   deallocate(patterninteger,patternad,patternintd,imagedictflt,imagedictfltflip)

! and we end the parallel section here (all threads will synchronize).
!$OMP END PARALLEL

if (cancelled.eqv..TRUE.) EXIT dictionaryloop

end do dictionaryloop

!-----
ierr = clReleaseMemObject(cl_dict)
call CL%error_check('DIdriver:clReleaseMemObject:cl_dict', ierr)

!-----
ierr = clReleaseMemObject(cl_expt)
call CL%error_check('DIdriver:clReleaseMemObject:cl_expt', ierr)

if (cancelled.eqv..FALSE.) then

  if (dinl%keeptmpfile.eq.'n') then
      close(itmpexpt,status='delete')
  else
      close(itmpexpt,status='keep')
  end if

! release the OpenCL kernel
  ierr = clReleaseKernel(kernel)
  call CL%error_check('DIdriver:clReleaseKernel', ierr)

  if (trim(dinl%indexingmode).eq.'static') then
    call HDF%pop(.TRUE.)
  end if

! perform some timing stuff
  call timer%Time_tock(2)
  tstop = timer%getInterval(2)
  io_real(1) = tstop
  call Message%WriteValue(' Indexing duration (system_clock, s)                : ',io_real,1,"(/,F14.3)")
  io_real(1) = float(totnumexpt)*float(FZcnt) / tstop
  call Message%WriteValue(' Number of pattern comparisons per second           : ',io_real,1,"(/,F14.3)")
  io_real(1) = float(totnumexpt) / tstop
  call Message%WriteValue(' Number of experimental patterns indexed per second : ',io_real,1,"(/,F14.3,/)")

! ===================
! MAIN OUTPUT SECTION
! ===================

! fill the ipar array with integer parameters that are needed to write the h5ebsd file
! (anything other than what is already in the dinl structure)
  ipar = 0
  ipar(1) = nnk
  ipar(2) = Ne*ceiling(float(totnumexpt)/float(Ne))
  ipar(3) = totnumexpt
  ipar(4) = Nd*ceiling(float(FZcnt)/float(Nd))
  ipar(5) = FZcnt
  ipar(6) = pgnum
  if (ROIselected.eqv..TRUE.) then
    ipar(7) = dinl%ROI(3)
    ipar(8) = dinl%ROI(4)
  else
    ipar(7) = dinl%ipf_wd
    ipar(8) = dinl%ipf_ht
  end if 

  allocate(OSMmap(jjend, iiiend))

  call timer%makeTimeStamp()
  tstre = timer%getTimeString()

  if (dinl%datafile.ne.'undefined') then 
    vendor = 'TSL'
    fname = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(dinl%datafile)
    call DIFT%setfilename(fname) 
    call DIFT%h5_writeFile(EMsoft, HDF, HDFnames, vendor, mcnl, xtalname, dstr, tstrb, tstre, ipar, resultmain, &
                           exptIQ, indexmain, eulerarray, dpmap, progname, nmldeffile, OSMmap)
    call Message%printMessage(' Data stored in h5 file : '//trim(dinl%datafile))
  end if

  VT = Vendor_T()
  call VT%set_Modality(MPFT%getModality())
  if (dinl%ctffile.ne.'undefined') then 
    fpar2(1) = mcnl%EkeV
    fpar2(2) = MCsig
    call VT%ctf_writeFile(EMsoft,cell,SG,dinl,ipar,fpar2,indexmain,eulerarray,resultmain, OSMmap, exptIQ)
    call Message%printMessage('Data stored in ctf file : '//trim(dinl%ctffile))
  end if
  
  if (dinl%angfile.ne.'undefined') then 
      fpar1(1) = WD
      call VT%ang_writeFile(EMsoft,cell,SG,dinl,ipar,fpar1,indexmain,eulerarray,resultmain,exptIQ)
      call Message%printMessage(' Data stored in ang file : '//trim(dinl%angfile))
  end if

! close the fortran HDF5 interface
  call closeFortranHDFInterface()

! if requested, we notify the user that this program has completed its run
  if (trim(EMsoft%getConfigParameter('Notify')).ne.'Off') then
    if (trim(dinl%Notify).eq.'On') then 
      NumLines = 3
      allocate(MessageLines(NumLines))
  
      call hostnm(c)
   
      MessageLines(1) = ' EMDI program has ended successfully'
      MessageLines(2) = ' Indexed data stored in '//trim(dinl%datafile)
      write (exectime,"(F15.0)") tstop  
      MessageLines(3) = ' Total execution time [s]: '//trim(exectime)
      TitleMessage = ' EMsoft on '//trim(c)
      i = PostMessage(EMsoft, MessageLines, NumLines, TitleMessage)
    end if
  end if
end if 
  
end associate

end subroutine DIdriver


!--------------------------------------------------------------------------
!
! SUBROUTINE:InnerProdGPU
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Perform the inner product computations for the dictionary approach
!
!> @param expt vector with list of observed patterns
!> @param dict vector with list of calculated patterns
!> @param Ne number of patterns in the expt vector
!> @param Nd number of patterns in the dict vector
!> @param L size of one single pattern
!> @param result result of the matrix multiplication
!> @param kernel opencl kernel pointer
!> @param context opencl context type
!> @param command_queue opencl command queue
!
!> @date 12/09/14  SS 1.0 original
!> @date 27/01/15  SS 1.1 modified to call the subroutine from mastersubroutine
!> @date 02/24/16 MDG 1.2 converted OpenCL calls to clfortran from fortrancl
!> @date 03/03/16 MDG 1.3 added C_NULL_CHAR to kernelname
!> @date 06/07/17 MDG 1.4 removed progoptions from Build Program call; caused some issues on Linux in Release mode
!> @date 11/13/17 MDG 2.0 moved several OpenCL init statements to main calling program
!--------------------------------------------------------------------------
recursive subroutine InnerProdGPU(CL,cl_expt,cl_dict,Ne,Nd,correctsize,results,numd,selnumd,kernel,context,command_queue)
!DEC$ ATTRIBUTES DLLEXPORT :: InnerProdGPU

use clfortran
use mod_CLsupport
use ISO_C_BINDING
use mod_io

IMPLICIT NONE

type(OpenCL_T),INTENT(INOUT)                        :: CL
integer(c_intptr_t),target,INTENT(INOUT)            :: cl_expt
!f2py intent(in,out) ::  cl_expt
integer(c_intptr_t),target,INTENT(INOUT)            :: cl_dict
!f2py intent(in,out) ::  cl_dict
integer(kind=4),INTENT(IN)                          :: Ne
integer(kind=4),INTENT(IN)                          :: Nd
real(kind=4),INTENT(OUT),target                     :: results(Ne*Nd)

integer(kind=4),INTENT(IN)                          :: correctsize
integer(kind=irg),INTENT(IN)                        :: numd, selnumd
integer(c_intptr_t),target,INTENT(INOUT)            :: context
!f2py intent(in,out) ::  context
integer(c_intptr_t),target,INTENT(INOUT)            :: kernel
!f2py intent(in,out) ::  kernel
integer(c_intptr_t),target,INTENT(INOUT)            :: command_queue
!f2py intent(in,out) ::  command_queue

integer(c_int32_t)                                  :: ierr, ierr2, pcnt
integer(c_intptr_t),target                          :: cl_result

real(kind=4)                                        :: dicttranspose(Nd*correctsize)
integer(kind=4),parameter                           :: iunit = 40
character(fnlen)                                    :: info ! info about the GPU
integer(kind=8),target                              :: globalsize(2),localsize(2)
integer(kind=4)                                     :: num,istat,i,j,ii,jj,kk, io_int(1)
integer(kind=4),target                              :: Wexp,Wdict
integer(kind=8)                                     :: size_in_bytes_expt,size_in_bytes_dict,size_in_bytes_result
integer(kind=irg)                                   :: irec

size_in_bytes_result = Ne*Nd*sizeof(results(1))
Wexp = correctsize
Wdict = Nd
localsize = (/16,16/)
globalsize = (/Ne,Nd/)

!=====================
! INITIALIZATION [mostly performed in the calling program]
!=====================

! create buffer
cl_result = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_result, C_NULL_PTR, ierr)
call CL%error_check('InnerProdGPU:clCreateBuffer', ierr)

! ---- 

! set kernel arguments
ierr =  clSetKernelArg(kernel, 0, sizeof(cl_expt), C_LOC(cl_expt))
call CL%error_check('InnerProdGPU:clSetKernelArg:cl_expt', ierr)

ierr = clSetKernelArg(kernel, 1, sizeof(cl_dict), C_LOC(cl_dict))
call CL%error_check('InnerProdGPU:clSetKernelArg:cl_dict', ierr)

ierr = clSetKernelArg(kernel, 2, sizeof(Wexp), C_LOC(Wexp))
call CL%error_check('InnerProdGPU:clSetKernelArg:Wexp', ierr)

ierr = clSetKernelArg(kernel, 3, sizeof(Wdict), C_LOC(Wdict))
call CL%error_check('InnerProdGPU:clSetKernelArg:Wdict', ierr)

ierr = clSetKernelArg(kernel, 4, sizeof(cl_result), C_LOC(cl_result))
call CL%error_check('InnerProdGPU:clSetKernelArg:cl_result', ierr)

!execute the kernel
ierr = clEnqueueNDRangeKernel(command_queue, kernel, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), &
                              0, C_NULL_PTR, C_NULL_PTR)
call CL%error_check('InnerProdGPU:clEnqueueNDRangeKernel', ierr)

! wait for the commands to finish
ierr = clFinish(command_queue)
call CL%error_check('InnerProdGPU:clFinish', ierr)

! read the resulting vector from device memory
ierr = clEnqueueReadBuffer(command_queue,cl_result,CL_TRUE,0_8,size_in_bytes_result,C_LOC(results(1)),0,C_NULL_PTR,C_NULL_PTR)
call CL%error_check('InnerProdGPU:clEnqueueReadBuffer', ierr)

ierr = clReleaseMemObject(cl_result)
call CL%error_check('InnerProdGPU:clReleaseMemObject:cl_result', ierr)

! ---
end subroutine InnerProdGPU
!--------------------------------------------------------------------------

end module mod_DI