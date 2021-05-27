! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_PEDkin
    !! author: MDG 
    !! version: 1.0 
    !! date: 01/22/20
    !!
    !! class definition for the EMPEDkin program
  
  use mod_kinds
  use mod_global
  
  IMPLICIT NONE 
  
  ! namelist for the EMPEDkin program
  type, public :: PEDkinNameListType
    integer(kind=irg)       :: npix
    integer(kind=irg)       :: ncubochoric
    integer(kind=irg)       :: nthreads
    real(kind=sgl)          :: voltage
    real(kind=sgl)          :: thickness
    real(kind=sgl)          :: rnmpp
    real(kind=sgl)          :: dmin
    character(fnlen)        :: xtalname
    character(fnlen)        :: outname
    character(fnlen)        :: eulerfile
    character(4)            :: sampling
  end type PEDkinNameListType
  
  ! class definition
  type, public :: PEDkin_T
  private 
    character(fnlen)       :: nmldeffile = 'EMPEDkin.nml'
    type(PEDkinNameListType)  :: nml 
  
  contains
  private 
    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: PEDkin_
  
    generic, public :: getNameList => getNameList_
    generic, public :: writeHDFNameList => writeHDFNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: PEDkin => PEDkin_
  
  end type PEDkin_T
  
  ! the constructor routine for this class 
  interface PEDkin_T
    module procedure PEDkin_constructor
  end interface PEDkin_T
  
  contains
  
  !--------------------------------------------------------------------------
  type(PEDkin_T) function PEDkin_constructor( nmlfile ) result(PEDkin)
  !DEC$ ATTRIBUTES DLLEXPORT :: PEDkin_constructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! constructor for the PEDkin_T Class; reads the name list 
   
  IMPLICIT NONE
  
  character(fnlen), OPTIONAL   :: nmlfile 
  
  call PEDkin%readNameList(nmlfile)
  
  end function PEDkin_constructor
  
  !--------------------------------------------------------------------------
  subroutine PEDkin_destructor(self) 
  !DEC$ ATTRIBUTES DLLEXPORT :: PEDkin_destructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! destructor for the PEDkin_T Class
   
  IMPLICIT NONE
  
  type(PEDkin_T), INTENT(INOUT)  :: self 
  
  call reportDestructor('PEDkin_T')
  
  end subroutine PEDkin_destructor
  
  !--------------------------------------------------------------------------
  subroutine readNameList_(self, nmlfile, initonly)
  !DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! read the namelist from an nml file for the PEDkin_T Class 
  
  use mod_io 
  use mod_EMsoft
  
  IMPLICIT NONE 
  
  class(PEDkin_T), INTENT(INOUT)          :: self
  character(fnlen),INTENT(IN)          :: nmlfile
   !! full path to namelist file 
  logical,OPTIONAL,INTENT(IN)          :: initonly
   !! fill in the default values only; do not read the file
  
  type(EMsoft_T)                       :: EMsoft 
  type(IO_T)                           :: Message       
  logical                              :: skipread = .FALSE.
  
  integer(kind=irg)       :: npix
  integer(kind=irg)       :: ncubochoric
  integer(kind=irg)       :: nthreads
  real(kind=sgl)          :: voltage
  real(kind=sgl)          :: thickness
  real(kind=sgl)          :: rnmpp
  real(kind=sgl)          :: dmin
  character(fnlen)        :: xtalname
  character(fnlen)        :: outname
  character(fnlen)        :: eulerfile
  character(4)            :: sampling

! define the IO namelist to facilitate passing variables to the program.
  namelist /PEDkinNameList/ xtalname, voltage, npix, rnmpp, ncubochoric, nthreads, &
  thickness, outname , dmin, sampling, eulerfile

  ! set the input parameters to default values (except for xtalname, which must be present)
  xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
  voltage = 200.0              ! acceleration voltage [kV]
  nthreads = 1                    ! number of OpenMP threads to start
  thickness = 10.0                ! sample thickness [nm]
  npix = 256                      ! output arrays will have size npix x npix
  outname = 'pedout.data'         ! output filename
  dmin = 0.04                     ! smallest d-spacing [nm]
  ncubochoric = 100               ! number of samples along the cubochoric edge length
  rnmpp = 0.2                     ! nm^{-1} per pattern pixel
  sampling = 'dict'
  eulerfile = 'undefined'

  if (present(initonly)) then
    if (initonly) skipread = .TRUE.
  end if

  if (.not.skipread) then
    ! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=PEDkinNameList)
    close(UNIT=dataunit,STATUS='keep')

    ! check for required entries
    if (trim(xtalname).eq.'undefined') then
      call Message%printError('EMPED:',' crystal structure file name is undefined in '//nmlfile)
    end if

    if (sampling .eq. 'file' .and. trim(eulerfile) .eq. 'undefined') then
      call Message%printError('EMPED:',' euler angle file name is undefined in '//nmlfile)
    end if

  end if

  ! if we get here, then all appears to be ok, and we need to fill in the self%nml fields
  self%nml%xtalname = xtalname
  self%nml%voltage = voltage
  self%nml%thickness = thickness
  self%nml%dmin = dmin
  self%nml%npix = npix
  self%nml%nthreads = nthreads
  self%nml%outname = outname
  self%nml%rnmpp = rnmpp
  self%nml%ncubochoric = ncubochoric
  self%nml%sampling = sampling
  self%nml%eulerfile = eulerfile  
  
  end subroutine readNameList_
  
  !--------------------------------------------------------------------------
  function getNameList_(self) result(nml)
  !DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! pass the namelist for the PEDkin_T Class to the calling program
  
  IMPLICIT NONE 
  
  class(PEDkin_T), INTENT(INOUT)          :: self
  type(PEDkinNameListType)                :: nml
  
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
  
  class(PEDkin_T), INTENT(INOUT)          :: self 
  type(HDF_T), INTENT(INOUT)              :: HDF
  type(HDFnames_T), INTENT(INOUT)         :: HDFnames
  
  integer(kind=irg),parameter             :: n_int = 3, n_real = 4
  integer(kind=irg)                       :: hdferr,  io_int(n_int)
  real(kind=sgl)                          :: io_real(n_real)
  character(20)                           :: intlist(n_int), reallist(n_real)
  character(fnlen)                        :: dataset, sval(1),groupname
  character(fnlen,kind=c_char)            :: line2(1)
  
  associate( emnl => self%nml )
  
  ! create the group for this namelist
  groupname = trim(HDFnames%get_NMLlist())
  !groupname = SC_PEDkinNameList
  hdferr = HDF%createGroup(groupname)

  ! write all the single integers
  io_int = (/ emnl%npix, emnl%ncubochoric, emnl%nthreads /)
  intlist(1) = 'npix'
  intlist(2) = 'Ncubochoric'
  intlist(3) = 'nthreads'
  call HDF%writeNMLintegers(io_int, intlist, n_int)

  ! write all the single reals
  io_real = (/ emnl%voltage, emnl%thickness, emnl%dmin, emnl%rnmpp /)
  reallist(1) = 'voltage'
  reallist(2) = 'thickness'
  reallist(3) = 'dmin'
  reallist(4) = 'rnmpp'
  call HDF%writeNMLreals(io_real, reallist, n_real)

  ! write all the strings
  dataset = SC_outname
  line2(1) = emnl%outname
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDkinNameList: unable to create outname dataset', hdferr)

  dataset = SC_xtalname
  line2(1) = emnl%xtalname
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDkinNameList: unable to create xtalname dataset', hdferr)

  dataset = SC_eulerfile
  line2(1) = emnl%eulerfile
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDkinNameList: unable to create eulerfile dataset', hdferr)

  dataset = SC_sampling
  line2(1) = emnl%sampling
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
  if (hdferr.ne.0) call HDF%error_check('HDFwritePEDkinNameList: unable to create sampling dataset', hdferr)

  ! and pop this group off the stack
  call HDF%pop()

  end associate
  
  end subroutine writeHDFNameList_
  
  !--------------------------------------------------------------------------
  subroutine PEDkin_(self, EMsoft, progname, nmldeffile)
  !DEC$ ATTRIBUTES DLLEXPORT :: PEDkin_
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
  use mod_io
  use mod_timing
  use mod_diffraction
  use mod_symmetry
  use mod_quaternions
  use mod_rotations
  use mod_so3
  use mod_math
  use stringconstants
  use HDF5

  IMPLICIT NONE 
  
  class(PEDkin_T), INTENT(INOUT)          :: self
  type(EMsoft_T), INTENT(INOUT)           :: EMsoft
  character(fnlen), INTENT(INOUT)         :: progname
  character(fnlen), INTENT(INOUT)         :: nmldeffile

  type(Timing_T)                          :: timer
  type(SpaceGroup_T)                      :: SG
  type(IO_T)                              :: Message
  type(Cell_T)                            :: Cell
  type(Diffraction_T)                     :: Diff    
  type(gvectors_T)                        :: reflist
  type(HDF_T)                             :: HDF
  type(HDFnames_T)                        :: HDFnames
  type(so3_T)                             :: SO
  type(q_T)                               :: quat
  type(e_T)                               :: eu
  type(r_T)                               :: ro
  type(orientation_T)                     :: ot
  type(o_T)                               :: om

  integer(kind=irg)                       :: FZcnt, pgnum, ival, i, j, l
  type(FZpointd),pointer                  :: FZlist, FZtmp, FZtmp2
  real(kind=dbl)                          :: mLambda
  real(kind=sgl)                          :: la, dval, dmin, glen, gmax, io_real(3), k(3), sgmax, FN(3), xgmin, Ig, Igmax, & 
                                           maxint, w, ku(3), kp(3), rnmpp, dx, dy, tstart, tstop, x, y, ma, mi
  integer(kind=irg)                       :: gp(3), imh, imk, iml, nref, gg(3), ix, iy, iz, io_int(5), ww, nsize, tdp, sx, sy, hdferr, &
                                         ninbatch, nbatches, nremainder,ibatch,istat, gridtype, tickstart 
  integer(HSIZE_T)                        :: dims3(3), hdims(3), offset(3)
  logical                                 :: verbose, insert=.TRUE., overwrite=.TRUE., exists
  character(fnlen)                        :: groupname, dataset, outname
  character(fnlen)                        :: datagroupname
  character(11)                           :: dstr
  character(15)                           :: tstrb
  character(15)                           :: tstre


  character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
  character(fnlen,kind=c_char)                     :: line2(1)
  character(fnlen)                :: fname, atype

  real(kind=sgl),allocatable      :: pedpattern(:,:), xx(:,:), yy(:,:), line(:), dot(:,:), eulerarray(:,:), sglread(:)
  character(len=1),allocatable    :: pedp(:,:,:), pedpat(:,:)
  integer(kind=irg),allocatable   :: patinbatch(:)
  integer(kind=irg)               :: totnumberbatch, nextra
  type(DynType),save              :: Dyn
  type(gnode),save                :: rlp
  type(reflisttype),pointer       :: nexts, rltmpa

! type(HDFobjectStackType)          :: HDF_head


  sgmax = 0.50

  ! simplify the notation a little
  associate( cbednl => self%nml )

  call openFortranHDFInterface()
  ! set the HDF group names for this program
  HDF = HDF_T()
  HDFnames = HDFnames_T()
      
  ! initialize the timing routines
  timer = Timing_T()
  tstrb = timer%getTimeString()

!===========================================================================================
! CRYSTALLOGRAPHY
!===========================================================================================
  verbose = .TRUE.

  call cell%setFileName(cbednl%xtalname)
  
  call Diff%setrlpmethod('WK')
  call Diff%setV(dble(cbednl%voltage))
  
  call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, cbednl%dmin, verbose)
  
  ! determine the point group number
  j=0
  do i=1,32
      if (SGPG(i).le.SG%getSpaceGroupNumber()) j=i
  end do
  pgnum = j

  io_int(1) = pgnum
  call Message%WriteValue(' Point group number                       : ',io_int, 1, "(I3)")

!=============================================
!=============================================
! generation of all potential reflections inside a reciprocal space sphere
! computed from the camera length and the detector size ...

! first set the maximum |g| value that can possibly give rise to a diffracted beam on the detector (diagonal)
  gmax = sqrt(2.0) * float(cbednl%npix) * cbednl%rnmpp
  io_real(1) = gmax
  call Message%WriteValue(' Length of longest g-vector               : ', io_real, 1, "(F8.4)")


  reflist = gvectors_T()
  call reflist%Initialize_ReflectionList(cell, SG, Diff, gmax, Igmax, verbose)
  nref = reflist%get_nref()
  !=============================================
  !=============================================
  ! create the coordinate arrays for the Gaussian peaks that will represent the diffraction spots
  rnmpp = 1.0/cbednl%rnmpp
  ww = 6
  tdp = 2*ww+1
  allocate(xx(-ww:ww,-ww:ww), yy(-ww:ww,-ww:ww), line(-ww:ww), dot(-ww:ww,-ww:ww))
  line = (/ (float(i),i=-ww,ww) /) * rnmpp
  xx = spread(line,dim=1,ncopies=2*ww+1)
  yy = transpose(xx)


  !=============================================
  !=============================================
  ! create the output array
  nsize = cbednl%npix/2 + ww 
  allocate(pedpattern(-nsize:nsize,-nsize:nsize))
  maxint = Igmax 
  write (*,*) 'Maximum diffracted intensity : ',maxint


  
  !===============================
  ! HDF5 I/O
  ! write out the data to the file

  ! ! Create a new file using the default properties.
  outname = EMsoft%generateFilePath('EMdatapathname',trim(cbednl%outname))

  !=============================================
  ! create or update the HDF5 output file
  !=============================================
  call HDFnames%set_ProgramData(SC_PEDkin)
  call HDFnames%set_NMLlist(SC_PEDkinNameList)
  call HDFnames%set_NMLfilename(SC_PEDkinNML)

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


! then the remainder of the data in a EMData group
! [this is the list of variables needed for the IDL visualization program CBEDDisplay.pro and the EMCBEDpattern.f90 program]
  groupname = SC_EMData
  hdferr = HDF%createGroup(groupname)


!=====================================================
! SAMPLING OF RODRIGUES FUNDAMENTAL ZONE
!=====================================================
! if eulerfile is not defined, then we use the standard RFZ sampling;
! if it is defined, then we read the Eulerangle triplets from the file
! and generate the FZlist here... this can be useful to index patterns that
! have only a small misorientation range with respect to a known orientation,
! so that it is not necessary to scan all of orientation space.

 if (trim(cbednl%sampling).eq.'dict') then
  SO = so3_T(pgnum, zerolist='FZ')
  if (trim(cbednl%eulerfile).eq.'undefined') then
    call Message%printMessage(' Orientation space sampling mode set to RFZ')
    io_int(1) = pgnum
    io_int(2) = cbednl%ncubochoric
    call Message%WriteValue(' Point group number and number of cubochoric sampling points : ',io_int,2,"(I4,',',I5)")
    
    call SO%sampleRFZ(cbednl%ncubochoric)
    FZcnt = SO%getListCount('FZ')

    allocate(eulerarray(3,FZcnt))
    FZtmp => SO%getListHead('FZ')
    do i = 1,FZcnt
       eu = FZtmp%rod%re()
       eulerarray(1:3,i) = eu%e_copyd()
       FZtmp => FZtmp%next
    end do
 
  else
  ! read the euler angle file and create the linked list
    call SO%getOrientationsfromFile(cbednl%eulerfile)
    FZcnt = SO%getListCount('FZ')
    call Message%printMessage(' Orientation space sampling mode set to MIS')
    io_int(1) = pgnum
    io_int(2) = FZcnt
    call Message%WriteValue(' Point group number and number of sampling points : ',io_int,2,"(I4,',',I5)")
  end if
  io_int(1) = FZcnt

  call Message%WriteValue(' Number of unique orientations sampled :        : ', io_int, 1, "(I8)")
! we can now delete the linked list since we have the FZarray
 end if

! Euler angle array size
dataset = SC_FZcnt
hdferr = HDF%writeDatasetInteger(dataset, FZcnt)

! and write them to the HDF5 file
dataset = SC_EulerAngles
hdferr = HDF%writeDatasetFloatArray(dataset, eulerarray*180.0/sngl(cPi), 3, FZcnt)

!keep the datafile open for writing; we'll collect 1024 patterns in an array,
!and then use hyperslab writing to put them into the HDF5 file.

!=============================================
!=============================================
! and loop over all orientations

! allocate the hyperslab array that will hold the batches of PED patterns 
ninbatch = 1024
nbatches = floor(float(FZcnt)/float(ninbatch))
nremainder = mod(FZcnt,ninbatch)

nextra = 0
if(nremainder .gt. 0) nextra = 1
totnumberbatch = nbatches + nextra

allocate(patinbatch(totnumberbatch))
patinbatch = ninbatch
if(nextra .eq. 1) patinbatch(totnumberbatch) = nremainder
allocate(pedpat(-nsize:nsize,-nsize:nsize),stat=istat)

call Message%printMessage(' Starting main computation loop.')
call Message%printMessage(' ')

mLambda = Diff%getWaveLength()

batchloop: do i = 1, totnumberbatch       ! loop over all batches in the eulerarray
! allocate and set the individual output pattern to zero
    if(allocated(pedp)) deallocate(pedp)
    allocate(pedp(cbednl%npix,cbednl%npix,patinbatch(i)),stat=istat)
    pedp = ' '
    pedpattern = 0.0
    pedpat = ' '
    orientationloop: do l = 1,patinbatch(i)
        pedpattern = 0.0
        pedpat = ' '
        ival = (i-1)*patinbatch(i) + l

 ! convert the rodrigues vector to a passive rotation matrix.
        eu = e_T( edinp = dble(eulerarray(1:3,ival)) )
        om = eu%eo()
        

! multiplication with (0,0,1) produces the normalized beam direction in a
! cartesian reference frame; so now we can compute the excitation errors 
! for every reflection and keep only the ones that are sufficiently small
        k = (/ 0.0, 0.0, 1.0 /)
        ku = matmul(om%o_copyd(),k)
        FN = ku
        k = ku/sngl(mLambda)

! first we go through the entire reflection list and compute the excitation errors
! those points that satisfy the cutoff are linked via the nexts pointers
        rltmpa => reflist%Get_ListHead()
        rltmpa => rltmpa%next
        nexts => rltmpa
        do j=1,nref
          gg = rltmpa%hkl
          rltmpa%sg = Diff%Calcsg(cell,float(gg),k,FN)
! should we consider this point any further ? If so, add it to the strong reflection linked list
          if (abs(rltmpa%sg).le.sgmax) then 
            nexts%nexts => rltmpa
            nexts => rltmpa
          end if
          rltmpa => rltmpa%next
        end do

! then, for each point in the nexts list, we compute the components of k' = k+g+s
! and place them in the proper reference frame; we skip the incident beam since it is 
! meaningless in the kinematical approximation
        nexts => reflist%Get_ListHead()
      nexts => nexts%next%nexts
      do 
! determine the vector k'
        kp = k + float(nexts%hkl) + nexts%sg*ku
        kp = matmul(transpose(om%o_copyd()),kp)

! get the intensity for each point
        w = sngl(cPi)*nexts%sg*cbednl%thickness
        if (abs(w).lt.1.0e-6) then
          Ig = nexts%sangle  ! * (sngl(cPi)*cbednl%thickness/nexts%xg)**2
        else 
          Ig = nexts%sangle * (sin(w)/w)**2 ! * (sngl(cPi)*cbednl%thickness/nexts%xg)**2
        end if

! determine the spot coordinates on the detector
        x = rnmpp * kp(1)
        y = rnmpp * kp(2)

! and plot that spot as a small Gaussian in the pedpattern array, assuming it falls on the detector.
        if ((abs(x).le.nsize-ww).and.(abs(y).le.nsize-ww)) then
          sx = nint(x)
          sy = nint(y)
          dx = x-sx
          dy = y-sy
          dot = (Ig/Igmax)**0.2 * exp(-((xx-dx)**2+(yy-dy)**2)*0.0025)
          pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) = pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) + dot(-ww:ww,-ww:ww)
        end if

! and repeat this until the end of the list
        if (.not. associated(nexts%nexts)) EXIT
        nexts => nexts%nexts
      end do

      ma = maxval(pedpattern)
      mi = minval(pedpattern)
      pedpattern = ((pedpattern - mi)/ (ma-mi))
      pedpat = char(nint(255.0*pedpattern))

! save the pedpattern to the pedp array
      pedp(1:cbednl%npix,1:cbednl%npix,l) = pedpat(-nsize+ww+1:nsize-ww,-nsize+ww+1:nsize-ww)

! ! reset the nexts linked list and start over
      nexts => reflist%Get_ListHead()
      nexts => nexts%next
      rltmpa => nexts%nexts
      do 
        nullify(nexts%nexts)
        if (.not. associated(rltmpa%nexts)) EXIT
        nexts => rltmpa
        rltmpa => rltmpa%nexts
      end do
    end do orientationloop

    offset = (/ 0, 0, (i-1)*patinbatch(i) /)
    hdims = (/ cbednl%npix, cbednl%npix, FZcnt /)
    dims3 = (/ cbednl%npix, cbednl%npix, patinbatch(i) /)

    if (i.eq.1) then
        dataset = SC_PEDpatterns
        hdferr = HDF%writeHyperslabCharArray(dataset, pedp, hdims, offset, dims3)
    else
        dataset = SC_PEDpatterns
        hdferr = HDF%writeHyperslabCharArray(dataset, pedp, hdims, offset, dims3, insert)
    end if

    io_real(1) = float(sum(patinbatch(1:i)))*100/float(sum(patinbatch))
    call Message%WriteValue(' Completed ',io_real,1,'(F10.2," % ")')

 end do batchloop

 !call HDF%pop()

! ! and update the end time
! call timestamp(datestring=dstr, timestring=tstre)
! groupname = SC_EMheader
! hdferr = HDF_openGroup(groupname, HDF_head)

! groupname = SC_PEDkin
! hdferr = HDF_openGroup(groupname, HDF_head)

! ! stop time /EMheader/StopTime 'character'
! dataset = SC_StopTime
! line2(1) = dstr//', '//tstre
! hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, overwrite)

! tstop = Time_tock(tickstart)
! dataset = SC_Duration
! hdferr = HDF_writeDatasetFloat(dataset, tstop, HDF_head)

! leave this group and close the file
  call HDF%pop(.TRUE.)

  end associate

  end subroutine PEDkin_
  
  
  
  end module mod_PEDkin
