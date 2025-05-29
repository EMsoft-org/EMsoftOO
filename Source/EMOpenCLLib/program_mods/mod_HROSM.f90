! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_HROSM
  !! author: MDG 
  !! version: 1.0 
  !! date: 05/21/25
  !!
  !! class definition for the EMHROSM program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMHROSM program
type, public :: HROSMNameListType
  integer(kind=irg)       :: nsamples   ! number of sampling points along radius of misorientation ball
  real(kind=sgl)          :: gangle     ! [deg] max grain misorientation angle for clustering
  real(kind=sgl)          :: misorang   ! [deg] misorientation ball radius for sampling
  character(fnlen)        :: dpfile     ! input dot product file
  character(fnlen)        :: OSMfile    ! output HDF5 file
  character(fnlen)        :: OSMtiff    ! new high resolution Orientation Similarity Map
  character(fnlen)        :: IPFmap     ! prefix of optional IPF maps
end type HROSMNameListType

! class definition
type, public :: HROSM_T
private 
  character(fnlen)        :: nmldeffile = 'EMHROSM.nml'
  type(HROSMNameListType) :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  ! procedure, pass(self) :: mygetNameList_
  procedure, pass(self) :: HROSM_

  ! generic, private :: mygetNameList => mygetNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: HROSM => HROSM_

end type HROSM_T

! the constructor routine for this class 
interface HROSM_T
  module procedure HROSM_constructor
end interface HROSM_T

contains

!--------------------------------------------------------------------------
type(HROSM_T) function HROSM_constructor( nmlfile ) result(HROSM)
!! author: MDG 
!! version: 1.0 
!! date: 05/21/25
!!
!! constructor for the HROSM_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call HROSM%readNameList(nmlfile)

end function HROSM_constructor

!--------------------------------------------------------------------------
subroutine HROSM_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 05/21/25
!!
!! destructor for the HROSM_T Class
 
IMPLICIT NONE

type(HROSM_T), INTENT(INOUT)  :: self 

call reportDestructor('HROSM_T')

end subroutine HROSM_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/21/25
!!
!! read the namelist from an nml file for the HROSM_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(HROSM_T), INTENT(INOUT)       :: self
character(fnlen),INTENT(IN)         :: nmlfile
 !! full path to namelist file
logical,OPTIONAL,INTENT(IN)         :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                      :: EMsoft 
type(IO_T)                          :: Message       
logical                             :: skipread = .FALSE.

integer(kind=irg)                   :: nsamples   ! number of sampling points along radius of misorientation ball
real(kind=sgl)                      :: gangle     ! [deg] max grain misorientation angle for clustering
real(kind=sgl)                      :: misorang   ! [deg] misorientation ball radius for sampling
character(fnlen)                    :: dpfile     ! input dot product file
character(fnlen)                    :: OSMfile    ! output HDF5 file
character(fnlen)                    :: OSMtiff    ! new high resolution Orientation Similarity Map
character(fnlen)                    :: IPFmap

! define the IO namelist to facilitate passing variables to the program.
namelist  / HROSMdata / nsamples, gangle, misorang, dpfile, OSMfile, OSMtiff, IPFmap 

nsamples = 20           ! number of sampling points along radius of misorientation ball
gangle = 5.0            ! [deg] max grain misorientation angle for clustering
misorang = 5.0          ! [deg] misorientation ball radius for sampling
dpfile = 'undefined'    ! input dot product file
OSMfile = 'undefined'   ! output HDF5 file
OSMtiff ='undefined'    ! new high resolution Orientation Similarity Map
IPFmap = 'undefined'    ! prefix for optional IPF maps

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 write (*,*) ' opening '//trim(nmlfile)
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=HROSMdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(dpfile).eq.'undefined') then
  call Message%printError('readNameList:',' dot product file name is undefined in '//nmlfile)
 end if

 if (trim(OSMfile).eq.'undefined') then
  call Message%printError('readNameList:',' OSM output file name is undefined in '//nmlfile)
 end if
end if

self%nml%nsamples = nsamples 
self%nml%gangle = gangle
self%nml%misorang = misorang
self%nml%dpfile = dpfile
self%nml%OSMfile = OSMfile
self%nml%OSMtiff = OSMtiff
self%nml%IPFmap = IPFmap

end subroutine readNameList_

! !--------------------------------------------------------------------------
! function mygetNameList_(self) result(nml)
! !DEC$ ATTRIBUTES DLLEXPORT :: mygetNameList_
! !! author: MDG 
! !! version: 1.0 
! !! date: 05/21/25
! !!
! !! pass the namelist for the HROSM_T Class to the calling program

! IMPLICIT NONE 

! class(HROSM_T), INTENT(INOUT)          :: self
! type(HROSMNameListType)                :: nml

! nml = self%nml

! end function mygetNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/21/25
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(HROSM_T), INTENT(INOUT)           :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 1, n_real = 2
integer(kind=irg)                       :: hdferr,  io_int(n_int)
real(kind=sgl)                          :: io_real(n_real)
character(20)                           :: intlist(n_int), reallist(n_real)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)

associate( enl => self%nml )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ enl%nsamples /)
intlist(1) = 'nsamples'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ enl%gangle, enl%misorang /)
reallist(1) = 'gangle'
reallist(2) = 'misorang'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the strings
dataset = 'dpfile'
line2(1) = trim(enl%dpfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create dpfile dataset', hdferr)

dataset = 'OSMfile'
line2(1) = trim(enl%OSMfile)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create OSMfile dataset', hdferr)

dataset = 'OSMtiff'
line2(1) = trim(enl%OSMtiff)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create OSMtiff dataset', hdferr)

dataset = 'IPFmap'
line2(1) = trim(enl%IPFmap)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create IPFmap dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine HROSM_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: HROSM_
!! author: MDG 
!! version: 1.0 
!! date: 05/21/25
!!
!! perform the computations

use mod_EMsoft
use mod_io
use mod_cluster
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_DIsupport
use mod_DIfiles
use mod_MCfiles
use mod_MPfiles
use mod_DI
use mod_io
use mod_image
use mod_IPF
use mod_IPFsupport
use mod_quaternions
use mod_rotations
use mod_so3
use mod_memory
use mod_timing
use ISO_C_BINDING
use mod_image
use mod_EBSD
use mod_crystallography
use mod_symmetry
use stringconstants

use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(HROSM_T), INTENT(INOUT)           :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(HDF_T)                             :: HDF
type(HDFnames_T)                        :: localHDFnames
type(IO_T)                              :: Message
type(Cell_T)                            :: cell
type(SpaceGroup_T)                      :: SG
type(EBSD_T)                            :: EBSD
type(DIfile_T)                          :: DIFT
type(MCfile_T)                          :: MCFT
type(MPfile_T)                          :: MPFT
type(DictionaryIndexingNameListType)    :: dinl
type(MCOpenCLNameListType)              :: mcnl
type(SEMmasterNameListType)             :: mpnl
type(Cluster_T)                         :: cluster
type(Quaternion_T)                      :: quat
type(QuaternionArray_T)                 :: sym, tmp, qAR
type(q_T)                               :: qu 
type(e_T)                               :: eu
type(r_T)                               :: ro  
type(so3_T)                             :: SO
type(memory_T)                          :: mem
type(Timing_T)                          :: timer
type(IPF_T)                             :: IPF 
type(IPFmap_T)                          :: IPFmap 


character(fnlen)                        :: DIfile, fname, xtalname, TIFF_filename, IPFmapfile, IPFmode
character(fnlen)                        :: dataname, datagroupname, groupname, attributename, dataset
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(2)                            :: listmode
integer(kind=irg)                       :: hdferr, io_int(2), nSamples, binx, biny, bindx, i, ir, ic, ROI(4), icnt, nt, &
                                           FZcnt, ii 
real(kind=sgl), allocatable             :: mainOSM(:,:), OSMmap(:,:), mainEuler(:,:,:), mainResult(:,:)  
real(kind=sgl)                          :: mi, ma
real(kind=sgl),allocatable              :: rodarray(:,:,:), maineu(:,:)
real(kind=sgl),allocatable              :: resultmain(:,:)
type(FZpointd),pointer                  :: FZlist, FZtmp

logical                                 :: verbose=.FALSE., f_exists
character(fnlen,kind=c_char)            :: HDF_FileVersion

! declare variables for use in object oriented image module
integer                                 :: iostat
character(len=128)                      :: iomsg
logical                                 :: isInteger
type(image_t)                           :: im
integer(int8)                           :: i8 (3,4)
integer(int8), allocatable              :: TIFF_image(:,:)


! outline of computations
! 1. read the dpfile to get all parameters, including refined orientations
! 2. use orientations to find grains via clustering algorithm in mod_cluster
! 3. compute the average orientation for each grain
! 4. loop over all grains:
!    a) generate sampling misorientation ball using nsamples value
!    b) generate a tmp file with pre-processed patterns only for this grain
!    c) perform DI calculation for this grain and the misorientation ball
!    d) compute the OSM for this grain
!    e) merge grain OSM values into main OSM
! 5. write all results to HDF5 OSMfile 
! 6. if requested, also produce a tiff merged OSM 

associate(osmnl=>self%nml, DIDT=>DIFT%DIDT, det=>EBSD%det, enl=>EBSD%nml)

timer = Timing_T()
tstrb = timer%getTimeString()

! 1. read the dpfile to get all parameters, including refined orientations
!    also read the Monte Carlo and Master pattern datasets
call openFortranHDFInterface()
HDF = HDF_T()
localHDFnames = HDFnames_T()

call localHDFnames%set_NMLfiles(SC_NMLfiles)
call localHDFnames%set_NMLfilename(SC_DictionaryIndexingNML)
call localHDFnames%set_NMLparameters(SC_NMLparameters)
call localHDFnames%set_NMLlist(SC_DictionaryIndexingNameListType)

DIfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(osmnl%dpfile)
call DIFT%readDotProductFile(EMsoft, HDF, localHDFnames, DIfile, hdferr, &
                             getRefinedEulerAngles = .TRUE.) 
dinl = DIFT%getNameList()

! 1a. read the Monte Carlo data file
call localHDFnames%set_ProgramData(SC_MCOpenCL)
call localHDFnames%set_NMLlist(SC_MCCLNameList)
call localHDFnames%set_NMLfilename(SC_MCOpenCLNML)
fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
call MCFT%setFileName(fname)
write (*,*) 'looking for ',trim(fname)
call MCFT%readMCfile(HDF, localHDFnames, getAccume=.TRUE.)
mcnl = MCFT%getnml()
xtalname = trim(mcnl%xtalname)
call Message%printMessage(' xtal file name : '//trim(xtalname))

! 1b. read the master pattern file
call localHDFnames%set_ProgramData(SC_EBSDmaster)
call localHDFnames%set_NMLlist(SC_EBSDmasterNameList)
call localHDFnames%set_NMLfilename(SC_EBSDmasterNML)
call localHDFnames%set_Variable(SC_MCOpenCL)
fname = EMsoft%generateFilePath('EMdatapathname',trim(dinl%masterfile))
call MPFT%setFileName(fname)
call MPFT%setModality('EBSD')
call MPFT%readMPfile(HDF, localHDFnames, mpnl, getmLPNH=.TRUE., getmLPSH=.TRUE.)

! 1c. we know that the master pattern file exists, and it also has all the
! crystallographic data in it, so we read that here instead of assuming
! that the actual .xtal file exists on this system ...
cell = Cell_T()
call cell%setFileName(mcnl%xtalname)
call cell%readDataHDF(SG, EMsoft)

! generate the detector arrays 
binx = dinl%exptnumsx/dinl%binning
biny = dinl%exptnumsy/dinl%binning
bindx = 1.0/float(dinl%binning)**2
! we also force the dictionary patterns to have this size 
dinl%numsx = binx
dinl%numsy = biny

mem = memory_T()
call mem%alloc(det%rgx, (/ dinl%numsx,dinl%numsy /), 'det%rgx') 
call mem%alloc(det%rgy, (/ dinl%numsx,dinl%numsy /), 'det%rgy') 
call mem%alloc(det%rgz, (/ dinl%numsx,dinl%numsy /), 'det%rgz') 
call mem%alloc(det%accum_e_detector, (/ MCFT%MCDT%numEbins,dinl%numsx,dinl%numsy /), 'det%accum_e_detector')
enl%numsx = dinl%numsx
enl%numsy = dinl%numsy
enl%xpc = dinl%xpc
enl%ypc = dinl%ypc
enl%delta = dinl%delta
enl%thetac = dinl%thetac
enl%L = dinl%L
enl%energymin = dinl%energymin
enl%energymax = dinl%energymax

call EBSD%GenerateDetector(MCFT, verbose)

! 2. use orientations to find grains via clustering algorithm in mod_cluster
! 3. this routine also does the orientation averaging using the von Mises-Fisher distribution...
cluster = Cluster_T( DIFT, osmnl%gangle )

io_int(1) = cluster%nGrains
call Message%WriteValue(' Number of grains found : ', io_int, 1)
call Message%printMessage(' Average grain orientations computed')
io_int(1) = count(cluster%kappa.eq.-1.D0)
call Message%WriteValue(' Number of non-converged average orientations : ', io_int, 1)

! save the grain ID map so we can take a look in IDL
! open(unit=dataunit,file='clustertest.data',status='unknown',form='unformatted')
! write (dataunit) cluster%grainID
! write (dataunit) real(cluster%avor)
! write (dataunit) real(cluster%kappa)
! close(unit=dataunit,status='keep')

! 4.  loop over all grains
! Since the DI step is parallel with GPU support, we need to do this grain by grain
! with each grain having a different number of patterns. We'll extract the relevant
! patterns from the main pattern file, pre-process them and store them in a tmp
! file.  We can read these in blocks but in many cases there will be a relatively 
! small number of patterns so we'll need to pad the array to maintain a multiple
! of 16. We'll need a slightly modified DIdriver routine to perform these runs; there
! is no need for this routine to produce the regular dp HDF5 file since that will
! be done by the present program.

allocate(mainOSM( cluster%ipf_wd, cluster%ipf_ht ) )
allocate(mainEuler( 3, cluster%ipf_wd, cluster%ipf_ht ) )
allocate(mainResult( cluster%ipf_wd, cluster%ipf_ht ) )
mainOSM = 0.0
mainEuler = 0.0
mainResult = 0.0
dinl%binning  = 1
nt = cluster%ipf_wd * cluster%ipf_ht

! do a portion of the array
! cluster%nGrains = 100
call Message%printMessage(' ')
call Message%printMessage(' ')

SO = so3_T( DIFT%DIDT%pgnum )
call SO%sample_isoCubeFilled(dble(osmnl%misorang), osmnl%nsamples)
listmode = 'CM'
io_int(1) = SO%getListCount(listmode)
call Message%WriteValue(' Starting indexing run; dictionary size : ', io_int,1)
call Message%printMessage(' ')
call SO%delete_FZlist('CM')

grainloop: do i=1,cluster%nGrains
  if (cluster%kappa(i).ne.-1.0) then 
    io_int = (/ i, cluster%nGrains /)
    call Message%WriteValue(' Indexing grain/total # grains ', io_int,2)
    io_int = (/ cluster%ROI(3,i), cluster%ROI(4,i) /)
    call Message%WriteValue(' OSM map size ', io_int,2)

  ! 4a. generate sampling misorientation ball around the averaged grain orientation;
  ! dictionary patterns will be computed for this list of orientations
  ! first generate the centered sampling cube; this is the same for all grains
    SO = so3_T( DIFT%DIDT%pgnum )
    call SO%sample_isoCubeFilled(dble(osmnl%misorang), osmnl%nsamples)
    listmode = 'CM'
    nSamples = SO%getListCount(listmode)
  ! then move the orientation ball to the averaged grain orientation
    qu = q_T( qdinp = cluster%avor(1:4,i) )
    ! write (*,*) ' avor, ROI = ', cluster%avor(1:4,i), cluster%ROI(1:4,i)
    ro = qu%qr()
    call SO%SampleIsoMisorientation(ro, dble(osmnl%misorang))

!     if (i.eq.2) then 
!       open(dataunit,file='grain2.txt',status='unknown',form='formatted')
!       FZcnt = SO%getListCount('CM')
! ! allocate and fill FZarray for OpenMP parallelization
!       FZtmp => SO%getListHead('CM')
!       do ii = 1,FZcnt
!           write(dataunit,"(4(F12.8,' '))") FZtmp%trod%r_copyd()
!           FZtmp => FZtmp%next
!       end do
!       close(dataunit,status='keep')
!     end if

  ! for each grain, we need to compute a bounding box that will be treated as the 
  ! standard ROI; using the dinl name list, we can modify the parameters for each grain
  ! along with the list of orientations to perform the DI run; the output is then 
  ! the best match orientations along with the list of N top-matches so that we can 
  ! compute an OSM for the ROI only, then copy those values into the overal OSM.
    dinl%ROI(1:4) = cluster%ROI(1:4, i)
    call OSMDIdriver(EMsoft, DIFT, MCFT, MPFT, dinl, mcnl, mpnl, cell, SG, EBSD, SO, &
                     OSMmap, resultmain, rodarray)

  ! copy the OSMmap parameters to the mainOSM array 
    do ic = 1, dinl%ipf_wd 
      do ir = 1, dinl%ipf_ht
        if (cluster%grainID(ic, ir).eq.i) then 
          mainOSM(ic,ir) = OSMmap(ic-cluster%ROI(1,i)+1, ir-cluster%ROI(2,i)+1)
          ro = r_T( rdinp = dble(rodarray(1:4,ic-cluster%ROI(1,i)+1, ir-cluster%ROI(2,i)+1)))
          eu = ro%re()
          mainEuler(1:3,ic,ir) = real( eu%e_copyd() )
          mainResult(ic,ir) = resultmain(ic-cluster%ROI(1,i)+1, ir-cluster%ROI(2,i)+1)
        end if 
      end do
    end do

  ! and delete the OSM array as well as the list of orientations
    deallocate(OSMmap, rodarray, resultmain)
    call SO%delete_FZlist('CM')
    call Message%printMessage(' ')
    ! if (i.eq.2) exit
  else
    io_int(1) = i 
    call Message%WriteValue(' skipping grain ', io_int,1)
    call Message%printMessage(' ')
  end if 
end do grainloop

timer = Timing_T()
dstr = timer%getDateString()
tstre = timer%getTimeString()

! used for initial debugging
! open(unit=dataunit, file ='OSM.data', status = 'unknown', form='unformatted')
! write(dataunit) mainOSM
! ! write(dataunit) indexmain
! close(unit=dataunit,status='keep')

! 5. write all results to an HDF5 file 

! get the filename; if it already exists, then delete it and create a new one
dataname = EMsoft%generateFilePath('EMdatapathname', osmnl%OSMfile)
inquire(file=trim(dataname), exist=f_exists)

if (f_exists) then
  open(unit=dataunit, file=trim(dataname), status='old',form='unformatted')
  close(unit=dataunit, status='delete')
end if

call localHDFnames%set_ProgramData(SC_HROSM)
call localHDFnames%set_NMLlist(SC_HROSMNameList)
call localHDFnames%set_NMLfilename(SC_HROSMNML)

! Create a new file using the default properties.
hdferr = HDF%createFile(dataname)

! write the EMheader to the file
datagroupname = trim(localHDFnames%get_ProgramData()) 
call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

! add the CrystalData group at the top level of the file
call cell%addXtalDataGroup(SG, EMsoft, HDF)

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(localHDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(localHDFnames%get_NMLfilename())    
hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
call HDF%pop()

! create a namelist group to write all the namelist files into
hdferr = HDF%createGroup(localHDFnames%get_NMLparameters())
call self%writeHDFNameList(HDF, localHDFnames)

! leave this group
call HDF%pop()

! then the remainder of the data in a EMData group
hdferr = HDF%createGroup(localHDFnames%get_EMData())

! here we add the data groupname and we attach to it a HDF_FileVersion attribute
hdferr = HDF%createGroup(datagroupname)
HDF_FileVersion = '4.0'
HDF_FileVersion = cstringify(HDF_FileVersion)
attributename = SC_HDFFileVersion
hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

dataset = 'nGrains'
  hdferr = HDF%writeDatasetInteger(dataset, cluster%nGrains)

dataset = 'grainID'
    hdferr = HDF%writeDatasetIntegerArray(dataset, cluster%grainID, cluster%ipf_wd, cluster%ipf_ht )

dataset = 'npixels'
    hdferr = HDF%writeDatasetIntegerArray(dataset, cluster%npixels, cluster%nGrains)

dataset = 'ROI'
    hdferr = HDF%writeDatasetIntegerArray(dataset, cluster%ROI, 4, cluster%nGrains)

dataset = 'avor'
    hdferr = HDF%writeDatasetDoubleArray(dataset, cluster%avor, 4, cluster%nGrains)

dataset = 'kappa'
    hdferr = HDF%writeDatasetDoubleArray(dataset, cluster%kappa, cluster%nGrains)

dataset = 'kam'
    hdferr = HDF%writeDatasetFloatArray(dataset, cluster%kam, cluster%ipf_wd, cluster%ipf_ht)

dataset = 'newOSM'
    hdferr = HDF%writeDatasetFloatArray(dataset, mainOSM, cluster%ipf_wd, cluster%ipf_ht)

dataset = 'newEuler'
    hdferr = HDF%writeDatasetFloatArray(dataset, mainEuler, 3, cluster%ipf_wd, cluster%ipf_ht)

dataset = 'newCI'
    hdferr = HDF%writeDatasetFloatArray(dataset, mainResult, cluster%ipf_wd, cluster%ipf_ht)

! =====================================================
! end of HDF_FileVersion = 4.0 write statements
! =====================================================

call HDF%popall()


! 6. if requested, also produce a tiff file with the mainOSM array
allocate(TIFF_image(dinl%ipf_wd,dinl%ipf_ht))
ma = maxval(mainOSM)
mi = minval(mainOSM)
TIFF_image = int(255 * (mainOSM - mi)/ (ma-mi))

im = image_t(TIFF_image)
if(im%empty()) call Message%printMessage("EMHROSM","failed to convert array to image")

! create the file
TIFF_filename = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(osmnl%OSMtiff)
call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
if(0.ne.iostat) then
  call Message%printMessage(" failed to write OSM map to file : "//iomsg)
else
  call Message%printMessage(' new OSM map written to '//trim(TIFF_filename))
end if

! 7. IPF maps if requested 
if (trim(osmnl%IPFmap).ne.'undefined') then 
! initialize the IPF map class; since we are not using an nmlfile argument here,
! we must manually initialize the parameters in this class
  IPF = IPF_T()
  allocate(maineu(3,nt))
  maineu = reshape(mainEuler,(/ 3, nt /))

  call tmp%QSym_Init(DIFT%DIDT%pgnum, sym)

! here we initialize the parameters of the IPF class; we will take a default file name 
! of IPFmapfile = 'IPFprefix_IPFZmap.tiff' with the current data path pre-pended.
  IPFmapfile = trim(osmnl%IPFmap)//'_IPFZmap.tiff'
  call IPF%set_nthreads(1)
  IPFmode = 'TSL'
  call IPF%set_IPFmode(IPFmode)
  qAR = QuaternionArray_T( n=nt, s='d' )
  do icnt=1,nt
    eu = e_T( edinp = dble(maineu(1:3,icnt))  )
    qu = eu%eq()
    quat = quaternion_T( qd = qu%q_copyd() )
    call qAR%insertQuatinArray( icnt, quat)
    end do 
! note the switch of x and y to get the same IPF map convention as DREAM.3D
  IPFmapfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(osmnl%IPFmap)//'_IPFXmap.tiff'
  call IPF%set_IPFfilename(IPFmapfile)
  call IPF%set_sampleDir( (/ 0, 1, 0 /) )
  call IPF%updateIPFmap(EMsoft, progname, cluster%ipf_wd, cluster%ipf_ht, DIFT%DIDT%pgnum, IPFmapfile, qAR, sym) 
  IPFmapfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(osmnl%IPFmap)//'_IPFYmap.tiff'
  call IPF%set_IPFfilename(IPFmapfile)
  call IPF%set_sampleDir( (/ 1, 0, 0 /) )
  call IPF%updateIPFmap(EMsoft, progname, cluster%ipf_wd, cluster%ipf_ht, DIFT%DIDT%pgnum, IPFmapfile, qAR, sym) 
  IPFmapfile = trim(EMsoft%generateFilePath('EMdatapathname'))//trim(osmnl%IPFmap)//'_IPFZmap.tiff'
  call IPF%set_IPFfilename(IPFmapfile)
  call IPF%set_sampleDir( (/ 0, 0, 1 /) )
  call IPF%updateIPFmap(EMsoft, progname, cluster%ipf_wd, cluster%ipf_ht, DIFT%DIDT%pgnum, IPFmapfile, qAR, sym) 
  call Message%printMessage(' IPF maps generated ')
end if

end associate

end subroutine HROSM_



end module mod_HROSM