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

module mod_PEDZA
    !! author: MDG 
    !! version: 1.0 
    !! date: 01/22/20
    !!
    !! class definition for the EMPEDZA program
  
  use mod_kinds
  use mod_global
  
  IMPLICIT NONE 
  
  ! namelist for the EMPEDZA program
  type, public :: PEDZANameListType
    integer(kind=irg)       :: k(3)
    integer(kind=irg)       :: fn(3)
    integer(kind=irg)       :: precsample
    integer(kind=irg)       :: precazimuthal
    integer(kind=irg)       :: npix
    integer(kind=irg)       :: nthreads
    real(kind=sgl)          :: voltage
    real(kind=sgl)          :: dmin
    real(kind=sgl)          :: precangle
    real(kind=sgl)          :: prechalfwidth
    real(kind=sgl)          :: thickness
    real(kind=sgl)          :: camlen
    character(5)            :: filemode
    character(fnlen)        :: xtalname
    character(fnlen)        :: outname
    character(fnlen)        :: BetheParametersFile
  end type PEDZANameListType
  
  ! class definition
  type, public :: PEDZA_T
  private 
    character(fnlen)          :: nmldeffile = 'EMPEDZA.nml'
    type(PEDZANameListType)   :: nml 
  
  contains
  private 
    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: PEDZA_
  
    generic, public :: getNameList => getNameList_
    generic, public :: writeHDFNameList => writeHDFNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: PEDZA => PEDZA_
  
  end type PEDZA_T
  
  ! the constructor routine for this class 
  interface PEDZA_T
    module procedure PEDZA_constructor
  end interface PEDZA_T
  
  contains
  
  !--------------------------------------------------------------------------
  type(PEDZA_T) function PEDZA_constructor( nmlfile ) result(PEDZA)
  !DEC$ ATTRIBUTES DLLEXPORT :: PEDZA_constructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! constructor for the PEDZA_T Class; reads the name list 
   
  IMPLICIT NONE
  
  character(fnlen), OPTIONAL   :: nmlfile 
  
  call PEDZA%readNameList(nmlfile)
  
  end function PEDZA_constructor
  
  !--------------------------------------------------------------------------
  subroutine PEDZA_destructor(self) 
  !DEC$ ATTRIBUTES DLLEXPORT :: PEDZA_destructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! destructor for the PEDZA_T Class
   
  IMPLICIT NONE
  
  type(PEDZA_T), INTENT(INOUT)  :: self 
  
  call reportDestructor('PEDZA_T')
  
  end subroutine PEDZA_destructor
  
  !--------------------------------------------------------------------------
  subroutine readNameList_(self, nmlfile, initonly)
  !DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! read the namelist from an self%nml file for the PEDZA_T Class 
  
  use mod_io 
  use mod_EMsoft
  
  IMPLICIT NONE 
  
  class(PEDZA_T), INTENT(INOUT)        :: self
  character(fnlen),INTENT(IN)          :: nmlfile
   !! full path to namelist file 
  logical,OPTIONAL,INTENT(IN)          :: initonly
   !! fill in the default values only; do not read the file
  
  type(EMsoft_T)                       :: EMsoft 
  type(IO_T)                           :: Message       
  logical                              :: skipread = .FALSE.
  
  integer(kind=irg)       :: k(3)
  integer(kind=irg)       :: fn(3)
  integer(kind=irg)       :: precsample
  integer(kind=irg)       :: precazimuthal
  integer(kind=irg)       :: npix
  integer(kind=irg)       :: nthreads
  real(kind=sgl)          :: voltage
  real(kind=sgl)          :: dmin
  real(kind=sgl)          :: precangle
  real(kind=sgl)          :: prechalfwidth
  real(kind=sgl)          :: thickness
  real(kind=sgl)          :: camlen
  character(5)            :: filemode
  character(fnlen)        :: xtalname
  character(fnlen)        :: outname
  
  ! define the IO namelist to facilitate passing variables to the program.
  namelist /EMPEDZA/ xtalname, voltage, k, fn, dmin, precangle, prechalfwidth, precsample, precazimuthal, &
                     thickness,  outname, npix, camlen, filemode, nthreads

  ! set the input parameters to default values (except for xtalname, which must be present)
  xtalname = 'undefined'          ! initial value to check that the keyword is present in the self%nml file
  voltage = 200.0                 ! acceleration voltage [kV]
  k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
  fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
  dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
  precangle = 10.472              ! beam precession angle [mrad]; default = 0.6 degrees
  prechalfwidth = 0.25            ! beam half width in the tilt direction [mrad]
  nthreads = 1                    ! number of OpenMP threads to start
  precsample = 10                 ! number of samples (concentric circles) in beam half width (total = 2*precsample + 1)
  precazimuthal = 360             ! number of azimuthal samples for each precession circle
  thickness = 10.0                ! sample thickness [nm]
  filemode = 'total'              ! 'total' mode or 'eachp'
  npix = 256                      ! output arrays will have size npix x npix
  outname = 'pedout.data'         ! output filename
  camlen = 1000.0                 ! camera length [mm]


  if (present(initonly)) then
    if (initonly) skipread = .TRUE.
  end if

  if (.not.skipread) then
  ! read the namelist file
  open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
  read(UNIT=dataunit,NML=EMPEDZA)
  close(UNIT=dataunit,STATUS='keep')

  ! check for required entries
  if (trim(xtalname).eq.'undefined') then
  call Message%printError('EMPEDZA:',' crystal structure file name is undefined in '//nmlfile)
  end if
  end if

  ! if we get here, then all appears to be ok, and we need to fill in the self%self%nml fields
  self%nml%xtalname = xtalname
  self%nml%voltage = voltage
  self%nml%k = k
  self%nml%fn = fn
  self%nml%dmin = dmin
  self%nml%precangle = precangle
  self%nml%prechalfwidth = prechalfwidth
  self%nml%precsample = precsample
  self%nml%precazimuthal = precazimuthal
  self%nml%thickness = thickness
  self%nml%filemode = filemode
  self%nml%npix = npix
  self%nml%nthreads = nthreads
  self%nml%outname = outname
  self%nml%camlen = camlen
  
  end subroutine readNameList_
  
  !--------------------------------------------------------------------------
  function getNameList_(self) result(nml)
  !DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! pass the namelist for the PEDZA_T Class to the calling program
  
  IMPLICIT NONE 
  
  class(PEDZA_T), INTENT(INOUT)          :: self
  type(PEDZANameListType)                :: nml

  nml = self%nml

  end function getNameList_
  
  !--------------------------------------------------------------------------
  recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
  !DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! write namelist to HDF file
  
  use mod_HDFsupport
  use mod_HDFnames
  use stringconstants 
  
  use ISO_C_BINDING
  
  IMPLICIT NONE
  
  class(PEDZA_T), INTENT(INOUT)           :: self 
  type(HDF_T), INTENT(INOUT)              :: HDF
  type(HDFnames_T), INTENT(INOUT)         :: HDFnames
  
  integer(kind=irg),parameter             :: n_int = 4, n_real = 6
  integer(kind=irg)                       :: hdferr, io_int(n_int)
  real(kind=sgl)                          :: io_real(n_real)
  character(20)                           :: intlist(n_int), reallist(n_real)
  character(fnlen)                        :: dataset, groupname
  character(fnlen,kind=c_char)            :: line2(1)  
  
  associate( pednl => self%nml )
  
    ! create the group for this namelist
  groupname = SC_PEDZANameList
  hdferr = HDF%createGroup(groupname)

  ! write all the single integers
  io_int = (/  pednl%precsample, pednl%precazimuthal, pednl%npix, pednl%nthreads /)
  intlist(1) = 'precsample'
  intlist(2) = 'precazimuthal'
  intlist(3) = 'npix'
  intlist(4) = 'nthreads'
  call HDF%writeNMLintegers(io_int, intlist, n_int)

  ! vectors
  dataset = SC_k
  hdferr = HDF%writeDatasetIntegerArray(dataset, pednl%k, 3)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDZANameList: unable to create k dataset', hdferr)

  dataset = SC_fn
  hdferr = HDF%writeDatasetIntegerArray(dataset, pednl%fn, 3)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDZANameList: unable to create fn dataset', hdferr)

  ! single reals
  io_real = (/ pednl%voltage, pednl%dmin, pednl%precangle, pednl%prechalfwidth, pednl%thickness, pednl%camlen /)
  reallist(1) = 'voltage'
  reallist(2) = 'dmin'
  reallist(3) = 'precangle'
  reallist(4) = 'prechalfwidth'
  reallist(5) = 'thickness'
  reallist(6) = 'camlen'
  call HDF%writeNMLreals(io_real, reallist, n_real)

  ! write all the strings
  dataset = SC_outname
  line2(1) = pednl%outname
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDZANameList: unable to create outname dataset', hdferr)

  dataset = SC_xtalname
  line2(1) = pednl%xtalname
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDZANameList: unable to create xtalname dataset', hdferr)

  dataset = SC_filemode
  line2(1) = pednl%filemode
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDZANameList: unable to create filemode dataset', hdferr)

  ! and pop this group off the stack
  call HDF%pop()

  end associate
  
  end subroutine writeHDFNameList_

  !--------------------------------------------------------------------------
  subroutine PEDZA_(self, EMsoft, progname, nmldeffile)
  !DEC$ ATTRIBUTES DLLEXPORT :: PEDZA_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! perform the computations
  
  use mod_EMsoft
  use mod_global
  use mod_crystallography
  use mod_HDFsupport
  use mod_HDFnames
  use mod_initializers
  use mod_gvectors
  use mod_kvectors
  use mod_io
  use mod_timing
  use mod_memory
  use mod_diffraction
  use mod_symmetry
  use mod_quaternions
  use mod_rotations
  use mod_so3
  use mod_math
  use mod_OMPsupport
  use stringconstants
  use HDF5

  IMPLICIT NONE 
  
  class(PEDZA_T), INTENT(INOUT)   :: self
  type(EMsoft_T), INTENT(INOUT)   :: EMsoft
  character(fnlen), INTENT(INOUT) :: progname 
  character(fnlen), INTENT(INOUT) :: nmldeffile

  type(Timing_T)                  :: timer
  type(SpaceGroup_T)              :: SG
  type(IO_T)                      :: Message
  type(Cell_T)                    :: Cell
  type(Diffraction_T)             :: Diff    
  type(gvectors_T)                :: reflist
  type(HDF_T)                     :: HDF
  type(HDFnames_T)                :: HDFnames
  type(kvectors_T)                :: kvec
  type(kvectorlist),pointer       :: ktmp
  type(memory_T)                  :: mem 
  
  real(kind=sgl)                  :: ktmax, io_real(3), bragg, thetac, sc, minten, pxy(2), galen, DM(2,2),DD,X(2),tstart,tstop, &
                                     frac, startthick, thick(1), klaue(2), thetam, kk(3), goffset,FN(3), kn, pmult
  integer(kind=irg)               :: ijmax,ga(3),gb(3),cnt, PX, numthick, ss, icnt, pgnum, ih, nunique, famnum, tickstart, &
                                     newcount,count_rate,count_max, io_int(6), i, j, isym, ir, skip, ghkl(3),hdferr, &
                                     npx, npy, numt, numk, ik, ip, jp, istat, dgn, nbeams, refcnt, gzero, TID, &
                                     ifamily, famhkl(3), inum, maxHOLZ, numksame, nns, nnw, nref, nt, NUMTHREADS
  character(3)                    :: method
  character(fnlen)                :: datafile, groupname, dataset, outname
  character(fnlen)                :: datagroupname

  real(kind=sgl),allocatable      :: disk(:,:),intarray(:),positions(:,:)
  integer(kind=irg),allocatable   :: familymult(:), familyhkl(:,:), whichHOLZ(:), gequiv(:,:), hklarray(:,:)
  real(kind=sgl),allocatable      :: inten(:,:)
  real(kind=dbl)                  :: s(3), mLambda
  logical                         :: verbose, silent, first, overwrite
  logical,allocatable             :: ksame(:)
  character(11)                   :: dstr
  character(15)                   :: tstrb
  character(15)                   :: tstre
  character(fnlen,kind=c_char)    :: line2(1)
  real(kind=sgl),allocatable      :: xx(:,:),yy(:,:),line(:),dot(:,:),pedpattern(:,:)
  character(len=1),allocatable    :: ped(:,:),pedpat(:,:)
  real(kind=sgl)                  :: rnmpp,Igmax,xc,yc,dx,dy,k(3),kp(3),Ig,ma,mi
  integer(kind=irg)               :: ww,tdp,sx,sy,nsize
  integer(HSIZE_T)                :: dims2(2)

  type(DynType),save              :: Dyn
  type(gnode),save                :: rlp
  type(reflisttype),pointer       :: firstw, rltmpa, khead
  type(BetheParameterType)        :: BetheParameters
  real(kind=sgl),allocatable      :: karray(:,:)
  integer(kind=irg),allocatable   :: kij(:,:)
  complex(kind=dbl),allocatable   :: DynMat(:,:)
  
  ! simplify the notation a little
  associate( pednl => self%nml )

  call openFortranHDFInterface()
  ! set the HDF group names for this program
  HDF = HDF_T()
  HDFnames = HDFnames_T()
      
  ! initialize the timing routines
  timer = Timing_T()
  tstrb = timer%getTimeString()
  call timer%Time_tick(1)

  mem = memory_T()

!===========================================================================================
! CRYSTALLOGRAPHY
!===========================================================================================
  verbose = .TRUE.

  call cell%setFileName(pednl%xtalname)
  
  call Diff%setrlpmethod('WK')
  call Diff%setV(dble(pednl%voltage))
 
  call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, pednl%dmin, verbose)

  ! set the foil normal 
  Dyn%FN = float(pednl%fn)
  call cell%NormVec(Dyn%FN, 'd')
  numt = 1
  goffset = 0.2
  thick(1) = pednl%thickness

  ! determine the point group number
  j=0
  do i=1,32
      if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
  end do
  isym = j  

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
  dgn = SG%GetPatternSymmetry(pednl%k,j,.TRUE.)
  pgnum = j
  isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers
  
  call cell%ShortestG(SG,pednl%k,ga,gb,isym)
  io_int(1:3)=ga(1:3)
  io_int(4:6)=gb(1:3)
  call Message%WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')')")
  
  kvec = kvectors_T()

  !=============================================
  !=============================================
  ! determine the list of contributing wave vectors
  call kvec%CalckvectorsPrecession(cell,Diff,dble(pednl%k),dble(ga),pednl%precangle,pednl%prechalfwidth,&
                                   pednl%precsample,pednl%precazimuthal,numk)
  
  !=============================================
  !=============================================
  ! set the scale parameter for a default camera length of 1000 mm.
  mLambda = Diff%getWaveLength()
  sc = mLambda * 1000.0 * 300.0 / 25.4  ! the absolute value does not matter and is derived from legacy Postscript code
  ! The original code used 300 dpi (hence 300/25.4) which was convenient for Postscript output; in the current case, we
  ! do not actually use the true value, but in the IDL visualization program, we scale the user defined camera length by
  ! 1000.0, and use this ratio to scale the diskoffset coordinates.  So, there's no absolute length scale, only a relative scale.

  ! next we create the output array 
  call mem%alloc(disk, (/ pednl%npix, pednl%npix /), 'disk', initval = 0.0, startdims = (/ -pednl%npix, -pednl%npix /))

  !===============================
  ! HDF5 I/O
  ! write out the data to the file

  ! ! Create a new file using the default properties.
  outname = EMsoft%generateFilePath('EMdatapathname',trim(pednl%outname))

  !=============================================
  ! create or update the HDF5 output file
  !=============================================
  call HDFnames%set_ProgramData(SC_PEDZA)
  call HDFnames%set_NMLlist(SC_PEDZANameList)
  call HDFnames%set_NMLfilename(SC_PEDZANML)

  ! Open an existing file or create a new file using the default properties.
  hdferr =  HDF%createFile(outname)

  ! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft,dstr, tstrb, tstre, progname, datagroupname)

    ! create a namelist group to write all the namelist files into
  groupname = SC_NMLfiles
  hdferr = HDF%createGroup(groupname)

  ! read the text file and write the array to the file
  dataset = SC_PEDkinNameList
  hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)

  ! leave this group
  call HDF%pop()

   ! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLparameters())
  call self%writeHDFNameList(HDF, HDFnames)

  ! leave this group
  call HDF%pop()

! then the remainder of the data in a EMData group
! [this is the list of variables needed for the IDL visualization program CBEDDisplay.pro and the EMCBEDpattern.f90 program]
  groupname = SC_EMData
  hdferr = HDF%createGroup(groupname)

  ! we need to write the image dimensions
  dataset = SC_numk
  hdferr = HDF%writeDatasetInteger(dataset, numk) 

  dataset = SC_ga
  hdferr = HDF%writeDatasetIntegerArray(dataset, ga, 3) 

  dataset = SC_gb
  hdferr = HDF%writeDatasetIntegerArray(dataset, gb, 3) 

  dataset = SC_lenga
  hdferr = HDF%writeDatasetFloat(dataset, cell%CalcLength(float(ga),'r')) 

  dataset = SC_lengb
  hdferr = HDF%writeDatasetFloat(dataset, cell%CalcLength(float(gb),'r')) 

  dataset = SC_angab
  hdferr = HDF%writeDatasetFloat(dataset, cell%CalcAngle(float(ga),float(gb),'r')) 

! we leave the file open until the program is done.
!=============================================
!=============================================

  DM(1,1) = cell%CalcDot(float(gb),float(gb),'c')
  DM(1,2) = -cell%CalcDot(float(ga),float(gb),'c')
  DM(2,1) = DM(1,2)
  DM(2,2) = cell%CalcDot(float(ga),float(ga),'c')
  DD = 1.0/(DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1))

  ! first convert the linked wave vector list to a set of regular arrays, so that we can 
! use OpenMP for this computation
  numk = kvec%get_numk()
  call mem%alloc(karray, (/ 4,numk /), 'karray', initval = 0.0) 
  call mem%alloc(kij, (/ 2,numk /), 'kij', initval = 0)
! point to the first beam direction
  ktmp => kvec%get_ListHead()
! and loop through the list, keeping k, kn, and i,j
  karray(1:3,1) = sngl(ktmp%k(1:3))
  karray(4,1) = sngl(ktmp%kn)
  kij(1:2,1) = (/ ktmp%i, ktmp%j /)
  do ik=2,numk
    ktmp => ktmp%next
    karray(1:3,ik) = sngl(ktmp%k(1:3))
    karray(4,ik) = sngl(ktmp%kn)
    kij(1:2,ik) = (/ ktmp%i, ktmp%j /)
  end do

  ! and remove the linked list
  call kvec%Delete_kvectorlist()

  !===========================================================================================
  !force dynamical matrix routine to read new Bethe parameters from file
  call Diff%SetBetheParameters(EMsoft, .FALSE., pednl%BetheParametersFile)

  ! the following call establishes the linked list of reflections and initializes their xg
  ! parameters to zero; these will then be filled in for the individual incident beam orientations
  kk(1:3) = float(pednl%k)
  FN = kk
  reflist = gvectors_T()
  call reflist%Initialize_ReflectionList(cell, SG, Diff, FN, kk, pednl%precangle, goffset,verbose)
  nref = reflist%get_nref()

  ! set the xg-field of this linked list to zero so that we can start accumulating intensities
  rltmpa => reflist%get_ListHead()
  rltmpa => rltmpa%next
  do icnt=1,nref
    rltmpa%xg = 0.0
    rltmpa=>rltmpa%next
  end do

  ! set up OpenMP parameters
  io_int(1)=numk
  call Message%WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")
  gzero = 1
  nt = 1
  
  verbose = .TRUE.

  ! set the foil normal 
  nullify(firstw)

  first = .TRUE.
  frac = 0.0
  
  ! loop over all beam orientations 
kvectorloop:  do ik = 1,numk
! generate the reflectionlist; this actually calls for a special routine, since the traditional ones
! are separate for the list generation and the subsequent application of the Bethe potentials; in the
! present case, we need to use an existing reflist, and in one pass determine the BethePotentials
! for a particular incident beam direction... we use a new routine called GetSubReflist() from gvectors module
   kk(1:3) = karray(1:3,ik)
   FN = kk
   call reflist%GetSubReflist(cell, Diff, FN, kk, firstw, nns, nnw, first)

   if (ik.eq.1) first = .FALSE.

 ! generate the dynamical matrix
    call mem%alloc(DynMat, (/ nns,nns /), 'DynMat', initval = CMPLX(0.D0,0.D0,8) )
    call reflist%GetDynMat(cell, Diff, firstw, DynMat, nns, nnw)

 ! allocate the intensity array to include both strong beams and weak beams (in that order)
    call mem%alloc(inten, (/ numt,nns /), 'inten', initval = 0.0)
 
 ! solve the dynamical eigenvalue equation for this beam direction
    kn = karray(4,ik)
    call reflist%CalcPEDint(cell,DynMat,kn,nns,nt,thick,inten)

 ! now we need to write these intensities to the appropriate slots in the reflist linked list... 
    rltmpa => reflist%Get_ListHead()
    rltmpa => rltmpa%next
    do icnt=1,nns
      rltmpa%xg = rltmpa%xg + inten(1,icnt)
      rltmpa=>rltmpa%nexts
    end do

! ! and remove the intensity array
    call mem%dealloc(DynMat, 'DynMat')
    call mem%dealloc(inten, 'inten')

! ! remove all the computed intensities for per pattern storage    
    if (pednl%filemode.eq.'eachp') then

! ! first we count how many reflections have none-zero intensity  
     refcnt = 0
     rltmpa => reflist%Get_ListHead()
     rltmpa => rltmpa%next
     do i=1,nref
       if (rltmpa%xg.ne.0.D0) refcnt = refcnt + 1
       rltmpa => rltmpa%next
     end do

! write refcnt
     write (dataunit) refcnt
     rltmpa => reflist%Get_ListHead()
     rltmpa => rltmpa%next
     do i=1,nref
       if (rltmpa%xg.ne.0.D0) then
! decompose this point w.r.t ga and gb
         X(1) = cell%CalcDot(float(rltmpa%hkl),float(ga),'c')
         X(2) = cell%CalcDot(float(rltmpa%hkl),float(gb),'c')
         X = matmul(DM,X) * DD
         write (dataunit) rltmpa%hkl
         write (dataunit) rltmpa%xg
         write (dataunit) X
       end if
       rltmpa => rltmpa%next
     end do
  
! and reset the intensities for the next run
     rltmpa => reflist%Get_ListHead()
     rltmpa => rltmpa%next
     do i=1,nref
       rltmpa%xg = 0.D0
       rltmpa => rltmpa%next
     end do
    end if

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call Message%WriteValue('       ', io_int, 1, "(1x,I3,'% completed')") 
    frac = frac + 0.05
   end if  

! reset the strong reflection pointers to null
   rltmpa => reflist%Get_ListHead()
   rltmpa => rltmpa%next
   do i=1,nref
     rltmpa%strong = .FALSE.
     rltmpa%weak = .FALSE.
     nullify(rltmpa%nextw)
     nullify(rltmpa%nexts)
     rltmpa=>rltmpa%next
   end do
   nullify(firstw)

  end do kvectorloop

  ! write the intensities as individual arrays to the HDF5 file
if (pednl%filemode.eq.'total') then

  ! first we count how many reflections have none-zero intensity  
    refcnt = 0
    rltmpa => reflist%Get_ListHead()
    rltmpa =>  rltmpa%next
     do i=1,nref
       if (rltmpa%xg.ne.0.D0) refcnt = refcnt + 1
       rltmpa => rltmpa%next
     end do
  
   ! write refcnt
   dataset = SC_refcnt
     hdferr = HDF%writeDatasetInteger(dataset, refcnt) 
  
  ! ! allocate arrays for output
     call mem%alloc(hklarray, (/ 3,refcnt /), 'hklarray', initval = 0) 
     call mem%alloc(intarray, (/ refcnt /), 'intarray', initval = 0.0)
     call mem%alloc(positions, (/ 2,refcnt /), 'positions', initval = 0.0)
  
     rltmpa => reflist%Get_ListHead()
     rltmpa =>  rltmpa%next
     icnt = 1
     do i=1,nref
       if (rltmpa%xg.ne.0.D0) then
    ! decompose this point w.r.t ga and gb
          X(1) = cell%CalcDot(float(rltmpa%hkl),float(ga),'c')
          X(2) = cell%CalcDot(float(rltmpa%hkl),float(gb),'c')
          X = matmul(DM,X) * DD
          hklarray(1:3,icnt) = rltmpa%hkl
          intarray(icnt) = rltmpa%xg
          positions(1:2,icnt) = X
          icnt = icnt+1
       end if
       rltmpa => rltmpa%next
     end do
  
  ! and write these arrays to the HDF5 file
  dataset = SC_hklarray
    hdferr = HDF%writeDatasetIntegerArray(dataset, hklarray, 3, refcnt) 
  
  dataset = SC_intarray
    hdferr = HDF%writeDatasetFloatArray(dataset, intarray, refcnt) 
  
  dataset = SC_positions
    hdferr = HDF%writeDatasetFloatArray(dataset, positions, 2, refcnt) 
  
  !=============================================
  !=============================================
  ! create the coordinate arrays for the Gaussian peaks that will represent the diffraction spots
    rnmpp = 1.0/0.075
    ww = 6
    tdp = 2*ww+1
    call mem%alloc(xx, (/ ww,ww /), 'xx', initval = 0.0, startdims = (/-ww, -ww/)) 
    call mem%alloc(yy, (/ ww,ww /), 'yy', initval = 0.0, startdims = (/-ww, -ww/)) 
    call mem%alloc(line, (/ ww /), 'line', initval = 0.0, startdims = (/ -ww /)) 
    call mem%alloc(dot, (/ ww,ww /), 'dot', initval = 0.0, startdims = (/-ww, -ww/)) 
    line = (/ (float(i),i=-ww,ww) /) * rnmpp
    xx = spread(line,dim=1,ncopies=2*ww+1)
    yy = transpose(xx)
  
  !=============================================
  !=============================================
  ! create the output array
    nsize = pednl%npix/2 + ww 
    call mem%alloc(pedpattern, (/ nsize, nsize /), 'pedpattern', initval = 0.0, startdims = (/-nsize,-nsize/))
    Igmax = maxval(intarray(2:refcnt))
    write (*,*) ' Maximum diffracted intensity : ',Igmax
    call mem%alloc(ped, (/ pednl%npix, pednl%npix /), 'ped', initval = char(0))
    call mem%alloc(pedpat, (/ nsize, nsize /), 'pedpat', initval = char(0), startdims = (/-nsize,-nsize/))
  
    k(1:3) = float(pednl%k)/mLambda
    rltmpa => reflist%Get_ListHead()
    rltmpa => rltmpa%next 
    do i = 2,refcnt
      rltmpa%sg = Diff%Calcsg(cell,float(hklarray(1:3,i)),k,FN)
      kp = k + float(hklarray(1:3,i)) + sngl(rltmpa%sg)
  
      xc = rnmpp * kp(1)
      yc = rnmpp * kp(2)

 ! and plot that spot as a small Gaussian in the pedpattern array, assuming it falls on the detector.

      if ((abs(xc).le.nsize-ww).and.(abs(yc).le.nsize-ww)) then
        sx = nint(xc)
        sy = nint(yc)
        dx = xc-sx
        dy = yc-sy
        Ig = intarray(i)
        dot = (Ig/Igmax)**0.2 * exp(-((xx-dx)**2+(yy-dy)**2)*0.0005)
        pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) = pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) + dot(-ww:ww,-ww:ww)
      end if
    end do
  
    ma = maxval(pedpattern)
    mi = minval(pedpattern)
    pedpattern = ((pedpattern - mi)/ (ma-mi))
    pedpat = char(nint(255.0*pedpattern))
    ped = ' '
  
  ! save the pedpattern to the ped array
    ped(1:pednl%npix,1:pednl%npix) = pedpat(-nsize+ww+1:nsize-ww,-nsize+ww+1:nsize-ww)
  
  dataset = SC_pedpattern
    dims2 =  (/ pednl%npix, pednl%npix /)
    hdferr = HDF%writeDatasetCharArray(dataset, ped, dims2)
  
    call HDF%pop()
  
    call timer%Time_tock(1)

  ! and update the end time
  groupname = SC_EMheader
    hdferr = HDF%openGroup(groupname)
  
  groupname = SC_PEDZA
    hdferr = HDF%openGroup(groupname)
  
  ! ! stop time /EMheader/StopTime 'character'
  ! dataset = SC_StopTime
  !   line2(1) = dstr//', '//tstre
  !   hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, overwrite)
  
    tstop = timer%getInterval(1)
  dataset = SC_Duration
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  
  ! ! close the datafile
    call HDF%popall()

    call Message%PrintMessage(' Data stored in '//pednl%outname,"(/A/)") 

    call mem%dealloc(hklarray, 'hklarray')
    call mem%dealloc(intarray, 'intarray')
    call mem%dealloc(positions, 'positions')
    call mem%dealloc(pedpattern, 'pedpattern')
    call mem%dealloc(pedpat, 'pedpat')
    call mem%dealloc(ped, 'ped')
    call mem%dealloc(line, 'line')
    call mem%dealloc(xx, 'xx')
    call mem%dealloc(yy, 'yy')
    call mem%dealloc(dot, 'dot')
  end if

  end associate

  call mem%dealloc(karray, 'karray')
  call mem%dealloc(kij, 'kij')

  end subroutine PEDZA_
  
  end module mod_PEDZA
