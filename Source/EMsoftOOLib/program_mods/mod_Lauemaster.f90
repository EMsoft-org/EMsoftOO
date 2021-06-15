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

module mod_Lauemaster
    !! author: MDG 
    !! version: 1.0 
    !! date: 01/22/20
    !!
    !! class definition for the EMLauemaster program
  
  use mod_kinds
  use mod_global
  
  IMPLICIT NONE 
  
  ! namelist for the EMLauemaster program
  type, public :: LauemasterNameListType
    integer(kind=irg)       :: npx
    integer(kind=irg)       :: patchw
    real(kind=sgl)          :: lambdamin
    real(kind=sgl)          :: lambdamax
    real(kind=dbl)          :: kappaVMF
    real(kind=dbl)          :: intfactor
    character(3)            :: outformat
    logical                 :: binarize
    character(fnlen)        :: SHT_folder
    character(fnlen)        :: SHT_formula
    character(fnlen)        :: SHT_name
    character(fnlen)        :: SHT_structuresymbol
    character(fnlen)        :: addtoKiltHub
    character(fnlen)        :: useDOI
    character(fnlen)        :: hdfname
    character(fnlen)        :: tiffname
    character(fnlen)        :: xtalname
  end type LauemasterNameListType
  
  ! class definition
  type, public :: Lauemaster_T
  private 
    character(fnlen)       :: nmldeffile = 'EMLauemaster.nml'
    type(LauemasterNameListType)  :: nml 
  
  contains
  private 
    procedure, pass(self) :: readNameList_
    procedure, pass(self) :: writeHDFNameList_
    procedure, pass(self) :: getNameList_
    procedure, pass(self) :: Lauemaster_
  
    generic, public :: getNameList => getNameList_
    generic, public :: writeHDFNameList => writeHDFNameList_
    generic, public :: readNameList => readNameList_
    generic, public :: Lauemaster => Lauemaster_
  
  end type Lauemaster_T
  
  ! the constructor routine for this class 
  interface Lauemaster_T
    module procedure Lauemaster_constructor
  end interface Lauemaster_T
  
  contains
  
  !--------------------------------------------------------------------------
  type(Lauemaster_T) function Lauemaster_constructor( nmlfile ) result(Lauemaster)
  !DEC$ ATTRIBUTES DLLEXPORT :: Lauemaster_constructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! constructor for the Lauemaster_T Class; reads the name list 
   
  IMPLICIT NONE
  
  character(fnlen), OPTIONAL   :: nmlfile 
  
  call Lauemaster%readNameList(nmlfile)
  
  end function Lauemaster_constructor
  
  !--------------------------------------------------------------------------
  subroutine Lauemaster_destructor(self) 
  !DEC$ ATTRIBUTES DLLEXPORT :: Lauemaster_destructor
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! destructor for the Lauemaster_T Class
   
  IMPLICIT NONE
  
  type(Lauemaster_T), INTENT(INOUT)  :: self 
  
  call reportDestructor('Lauemaster_T')
  
  end subroutine Lauemaster_destructor
  
  !--------------------------------------------------------------------------
  subroutine readNameList_(self, nmlfile, initonly)
  !DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! read the namelist from an nml file for the Lauemaster_T Class 
  
  use mod_io 
  use mod_EMsoft
  
  IMPLICIT NONE 
  
  class(Lauemaster_T), INTENT(INOUT)          :: self
  character(fnlen),INTENT(IN)          :: nmlfile
   !! full path to namelist file 
  logical,OPTIONAL,INTENT(IN)          :: initonly
   !! fill in the default values only; do not read the file
  
  type(EMsoft_T)                       :: EMsoft 
  type(IO_T)                           :: Message       
  logical                              :: skipread = .FALSE.
  
  integer(kind=irg)       :: npx
  integer(kind=irg)       :: patchw
  real(kind=sgl)          :: lambdamin
  real(kind=sgl)          :: lambdamax
  real(kind=dbl)          :: kappaVMF
  real(kind=dbl)          :: intfactor
  character(3)            :: outformat
  logical                 :: binarize
  character(fnlen)        :: SHT_folder
  character(fnlen)        :: SHT_formula
  character(fnlen)        :: SHT_name
  character(fnlen)        :: SHT_structuresymbol
  character(fnlen)        :: addtoKiltHub
  character(fnlen)        :: useDOI
  character(fnlen)        :: hdfname
  character(fnlen)        :: tiffname
  character(fnlen)        :: xtalname
  
  ! define the IO namelist to facilitate passing variables to the program.
  namelist  / LaueMasterData / npx, lambdamin, lambdamax, kappaVMF, hdfname, xtalname, &
            intfactor, tiffname, patchw, SHT_folder, SHT_formula, SHT_name, &
            SHT_structuresymbol, addtoKiltHub, useDOI, outformat, binarize

  npx = 500
  patchw = 5
  lambdamin = 0.10
  lambdamax = 0.16
  kappaVMF = 50000.D0
  intfactor = 0.0001D0
  outformat = 'LMP'
  SHT_folder = 'undefined'        ! folder to store SHT files, relative to EMDatapathname
  SHT_formula = 'undefined'       ! compound chemical formula, e.g., SiO2
  SHT_name = 'undefined'          ! compund name (e.g., forsterite)
  SHT_structuresymbol = 'undefined' ! StrukturBericht symbol (e.g., D0_22) or Pearson symbol (e.g., hP12), or ...
  addtoKiltHub = 'No'             ! file to be added to data base on kilthub.cmu.edu ?
  useDOI = 'undefined'            ! if no DOI is entered, then we use the Zenodo DOI for the .sht repository
  xtalname = 'undefined'
  hdfname = 'undefined'
  tiffname = 'undefined'
  binarize = .FALSE. 

  if (present(initonly)) then
    if (initonly) skipread = .TRUE.
  end if

  ! read the name list, depending on the class type
  if (.not.skipread) then
  ! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=LaueMasterData)
    close(UNIT=dataunit,STATUS='keep')
  
  ! check for required entries
    if (trim(xtalname).eq.'undefined') then
      call Message%printError('readNameList:',' crystal file name is undefined in '//nmlfile)
    end if

    if (outformat.eq.'SHT') then 
      ! for Legendre mode, the SHT_formula parameter MUST be present 
        if (trim(SHT_formula).eq.'undefined') then 
         call Message%printError('readNameList:',' SHT_formula must be defined in '//nmlfile)
        end if
     
        if (trim(SHT_folder).eq.'undefined') then 
         call Message%printError('readNameList:',' SHT_folder must be defined in '//nmlfile)
        end if
    end if 
  end if

  
  self%nml%npx = npx
  self%nml%patchw = patchw
  self%nml%lambdamin = lambdamin
  self%nml%lambdamax = lambdamax
  self%nml%kappaVMF = kappaVMF
  self%nml%intfactor = intfactor
  self%nml%xtalname = xtalname
  self%nml%outformat = outformat
  self%nml%hdfname = hdfname
  self%nml%tiffname = tiffname 
  self%nml%addtoKiltHub = addtoKiltHub
  self%nml%useDOI = useDOI
  self%nml%SHT_formula = SHT_formula
  self%nml%SHT_name = SHT_name
  self%nml%SHT_structuresymbol = SHT_structuresymbol
  self%nml%SHT_folder = trim(SHT_folder)
  self%nml%binarize = binarize

  end subroutine readNameList_
  
  !--------------------------------------------------------------------------
  function getNameList_(self) result(nml)
  !DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! pass the namelist for the Lauemaster_T Class to the calling program
  
  IMPLICIT NONE 
  
  class(Lauemaster_T), INTENT(INOUT)          :: self
  type(LauemasterNameListType)                :: nml
  
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
  
  class(Lauemaster_T), INTENT(INOUT)        :: self 
  type(HDF_T), INTENT(INOUT)              :: HDF
  type(HDFnames_T), INTENT(INOUT)         :: HDFnames
  
  integer(kind=irg),parameter             :: n_int = 11, n_real = 9
  integer(kind=irg)                       :: hdferr,  io_int(n_int)
  real(kind=sgl)                          :: io_real(n_real)
  character(20)                           :: intlist(n_int), reallist(n_real)
  character(fnlen)                        :: dataset, sval(1),groupname
  character(fnlen,kind=c_char)            :: line2(1)
  
  associate( mcnl => self%nml )
  
  end associate
  
  end subroutine writeHDFNameList_
  
  !--------------------------------------------------------------------------
  subroutine Lauemaster_(self, EMsoft, progname, nmldeffile)
  !DEC$ ATTRIBUTES DLLEXPORT :: Lauemaster_
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! perform the computations
  
  use mod_EMsoft
  use mod_initializers
  use mod_diffraction
  use mod_symmetry
  use mod_crystallography
  use mod_io
  use mod_gvectors
  use mod_kvectors
  use mod_math
  use mod_global
  use mod_timing
  use mod_Lambert
  use HDF5
  use mod_HDFsupport
  use mod_HDFnames
  use ISO_C_BINDING
  use omp_lib
  use mod_OMPsupport
  use stringconstants
  use mod_image
  use, intrinsic :: iso_fortran_env
  use mod_DSHT 
  use mod_fft_wrap
  use mod_LaueSupport

  IMPLICIT NONE 
  
  interface 
    recursive function writeShtFile (fn, iprm, fprm, doi, note, alm, aTy, aCd, vers, cprm) result(res) &
    bind(C, name ='writeShtFile_')

    use ISO_C_BINDING

    IMPLICIT NONE 

    character(c_char)             :: fn
    integer(c_int)                :: iprm(11) 
    real(c_float)                 :: fprm(25)
    character(c_char)             :: doi
    character(c_char)             :: note 
    real(C_DOUBLE_COMPLEX)        :: alm(2*(iprm(3)+3)*(iprm(3)+3))
    integer(c_int)                :: aTy(iprm(6))
    real(c_float)                 :: aCd(iprm(6),5)
    character(c_char)             :: vers 
    character(c_char)             :: cprm
    integer(c_int)                :: res
    end function writeShtFile
  end interface 


  class(Lauemaster_T), INTENT(INOUT)      :: self
  type(EMsoft_T), INTENT(INOUT)           :: EMsoft
  character(fnlen), INTENT(INOUT)         :: progname 
  character(fnlen),INTENT(IN)             :: nmldeffile

  type(Cell_T)                            :: cell
  type(DynType)                           :: Dyn
  type(Timing_T)                          :: timer
  type(Diffraction_T)                     :: Diff
  type(IO_T)                              :: Message
  type(Laue_T)                            :: reflist
  type(Lambert_T)                         :: L
  type(HDF_T)                             :: HDF
  type(SpaceGroup_T)                      :: SG
  type(kvectors_T)                        :: kvec
  type(HDFnames_T)                        :: HDFnames
  type(image_t)                           :: im, im2

  type(Laue_g_list), pointer :: rltmp
  logical 								                   :: verbose
  real(kind=sgl),allocatable                 :: mLPNH(:,:), mLPSH(:,:), masterSPNH(:,:), masterSPSH(:,:)
  integer(kind=irg)						               :: npx, npy, gcnt, ierr, nix, niy, nixp, niyp, i, j, w, istat, TIFF_nx, TIFF_ny, &
                                                hdferr, bw, d, ll, res, timestart, timestop, info, Lindex
  real(kind=sgl)							               :: xyz(3), kl(2), dx, dy, dxm, dym, Radius, mi, ma, tstart, tstop, sdev, mean
  real(kind=dbl)                             :: VMFscale, inten, p, LegendreLattitude
  character(fnlen)                           :: fname, hdfname, TIFF_filename, attributename, groupname, datagroupname, dataset, &
                                                 doiString, layout, SHTfile

  ! declare variables for use in object oriented image module
  integer                         :: iostat
  character(len=128)              :: iomsg
  logical                         :: isInteger, north, initLUT
  integer(int8), allocatable      :: TIFF_image(:,:)
  character(11)                   :: dstr
  character(15)                   :: tstrb
  character(15)                   :: tstre
  character(fnlen)                :: image_filename
  integer(int8)                   :: i8 (3,4), int8val
  integer(int8), allocatable      :: output_image(:,:)

  ! parameters for the .sht output file 
  integer(kind=irg),parameter     :: nipar=11, nfpar=25
  character(fnlen)                :: EMversion, cprm, note, notestring
  character(6)                    :: vstring
  character(8)                    :: vstring2
  character(fnlen)                :: revision
  integer(c_int32_t)              :: sgN      ! space group number [1,230]
  integer(c_int32_t)              :: sgS      ! space group setting [1,2]
  integer(c_int32_t)              :: numAt    ! number of atoms
  integer(c_int32_t),allocatable  :: aTy(:)   ! atom types (nAt atomic numbers)
  real(c_float),allocatable       :: aCd(:,:) ! atom coordinates, (nAt * 5 floats {x, y, z, occupancy, Debye-Waller in nm^2})
  real(c_float)                   :: lat(6)   ! lattice parameters {a, b, a, alpha, beta, gamma} (in nm / degree)
  real(c_float)                   :: fprm(nfpar) ! floating point parameters (float32 EMsoftED parameters in order)
  integer(c_int32_t)              :: iprm(nipar) ! integer parameters {# electrons, electron multiplier, numsx, npx, latgridtype}
  real(kind=dbl),allocatable      :: finalmLPNH(:,:), finalmLPSH(:,:), weights(:)
  real(kind=dbl),allocatable      :: LegendreArray(:), upd(:), diagonal(:)

  type(DiscreteSHT)               :: transformer 
  complex(kind=dbl), allocatable  :: almMaster(:,:)   ! spectra of master pattern to index against
  complex(kind=dbl), allocatable  :: almPat   (:,:)   ! work space to hold spectra of exerimental pattern
  real(kind=dbl),allocatable      :: alm(:)
  real(kind=sgl)                  :: dmin
  
  associate( lmnl => self%nml )

  ! basic explanation: this is a really simple and fast Laue master pattern; we compute all the plane normals 
  ! that fall inside the extended Ewald sphere volume.  For each we compute the kinematic intensity
  ! using the x-ray scattering factors.  Then we add a narrow Gaussian peak to the square Lambert projection
  ! (either Northern or Southern hemisphere) in the correct position, using spherical interpolation 
  ! (a von Mises-type distribution might be useful here ...).  Finally, standard output to an HDF5 file, or 
  ! output to an .sht file.

  ! lmnl components
  ! xtalname
  ! lambdamin
  ! lambdamax
  ! kappaVMF
  ! hdfname
  
  if (lmnl%outformat.eq.'SHT') then 
  npx = 193
    layout = 'legendre'
  else 
    npx = lmnl%npx
  end if
  

  call openFortranHDFInterface()
  HDF = HDF_T()

  ! set the HDF group names for this program
  HDFnames = HDFnames_T()

  ! initialize the timing routines
  timer = Timing_T()
  tstrb = timer%getTimeString()
  
  !=============================================
!=============================================
! crystallography section
  verbose = .TRUE.

  call cell%setFileName(lmnl%xtalname)

  dmin = 0.05
  call Diff%setV(dble(1.0))
  call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, useHDF=HDF)

  !=============================================
  !=============================================
  ! compute reflection list with kinematical intensities
  reflist = Laue_T()
  call  reflist%Init_Reflist(cell, SG, Diff,  gcnt, lmnl%lambdamin, lmnl%intfactor, verbose)

  !=============================================
  !=============================================
  ! populate the master pattern in square Lambert projection
  ! npx = lmnl%npx
    npy = npx
    allocate(mLPNH(-npx:npx,-npy:npy),stat=istat)
    allocate(mLPSH(-npx:npx,-npy:npy),stat=istat)
    mLPNH = 0.0
    mLPSH = 0.0

    
  !=============================================
  !=============================================
  ! precompute the Legendre array for the new lattitudinal grid values
    call Message%printMessage(' Computing Legendre lattitudinal grid values')
    allocate(diagonal(2*npx+1),upd(2*npx+1))
    diagonal = 0.D0
    upd = (/ (dble(i) / dsqrt(4.D0 * dble(i)**2 - 1.D0), i=1,2*npx+1) /)
    call dsterf(2*npx-1, diagonal, upd, info) 
  ! the eigenvalues are stored from smallest to largest and we need them in the opposite direction
    allocate(LegendreArray(0:2*npx))
    LegendreArray(0:2*npx) = diagonal(2*npx+1:1:-1)
  ! set the center eigenvalue to 0
    LegendreArray(npx) = 0.D0
    deallocate(diagonal, upd)

    ! the von Mises-Fisher distribution is defined by 
    !
    !   vmf(x;mu,kappa) = ( kappa / (4 pi sinh(kappa) ) ) exp[ kappa mu.x ]
    !
    ! and this is multiplied by the intensity;  since these can be really large numbers we 
    ! we will work with the large kappa expansion as well as logarithms; first of all, the 
    ! distribution becomes (for large kappa)
    !
    ! vmf(x;mu,kappa) = ( kappa exp[ kappa (mu.x-1) ]/  (2 pi)    (for large kappa)
    !
    ! multiplying this by the intensity I and taking the logarithm, we have
    !
    ! log(vmf) = (-1 + mu.x) kappa + Log(Inten) - Log(Pi) + Log(kappa) - Log(2) 
    ! 
    ! we'll take the constant part of this and call it VMFscale

    VMFscale = log(lmnl%kappaVMF) - log(2.D0) - log(cPi)

    ! set the size of the patches in the square Lambert space that we need to evaluate the VMF distribution for
    w = lmnl%patchw  ! this could become a part of the input namelist

    ! go through the entire reflection list
    rltmp => reflist%get_ListHead()

    do i=1,gcnt
  ! locate the nearest Lambert pixel (we need to make sure that the cartesian vector has unit length)
      call cell%NormVec(rltmp%xyz, 'c') 
  ! do we need to modify the direction cosines to coincide with the Legendre lattitudinal grid values?
      north = .TRUE.
      if (rltmp%xyz(3).lt.0) north=.FALSE.
      if (abs(rltmp%xyz(3)).ne.1.D0) then
        L = Lambert_T( xyz = sngl(rltmp%xyz))
        ierr = L%LambertSphereToSquare( kl ) 
        kl = kl * sngl(real(npx))
  ! here we need to be very careful to determine the index of the Legendre ring, NOT the Lambert ring !!!
        Lindex = npx 
        do while(LegendreArray(Lindex).lt.rltmp%xyz(3)) 
          Lindex = Lindex - 1
        end do  
        LegendreLattitude = LegendreArray( Lindex - 1)
  ! the factor p rescales the x and y components of kstar to maintain a unit vector
        p = sqrt((1.D0-LegendreLattitude**2)/(1.D0-rltmp%xyz(3)**2))
        rltmp%xyz = (/ p*rltmp%xyz(1), p*rltmp%xyz(2), LegendreLattitude /)
  ! rescale the coordinates in the Legendre square to be on the correct ring
        kl = kl * float(Lindex)/maxval(abs(kl))
      end if
        if (.not.north) rltmp%xyz(3) = -rltmp%xyz(3)
  ! and continue with the projection
        !call LambertgetInterpolation(sngl(rltmp%xyz), float(npx), npx, npy, nix, niy, nixp, niyp, dx, dy, dxm, dym)
  ! intensity with polarization correction
        inten = rltmp%sfs * rltmp%polar
        write(*,*) inten
        if (lmnl%binarize.eqv..TRUE.) inten = 1.0
  ! depending on the sign of xyz(3) we put this point in the Northern or Southern hemisphere, taking into account the
  ! special case of reflections along the equator which should appear in both hemisphere arrays.  The intensities are 
  ! computed on a small grid of w x w points on the Lambert projection, which are then interpolated from a "Gaussian" on
  ! the sphere. we use the von Mises-Fisher distribution with p=3
       call L%sampleVMF(sngl(rltmp%xyz), lmnl%kappaVMF, VMFscale, inten, npx, int(kl(1)), int(kl(2)), w, mLPNH, mLPSH, &
       LegendreArray)
  ! and go to the next point
      rltmp => rltmp%next
    end do 
    write(*,*) maxval(mLPNH)
  ! finally, make sure that the equator is copied into both arrays
    mLPSH(-npx,-npx:npx) = mLPNH(-npx,-npx:npx)
    mLPSH( npx,-npx:npx) = mLPNH( npx,-npx:npx)
    mLPSH(-npx:npx,-npx) = mLPNH(-npx:npx,-npx)
    mLPSH(-npx:npx, npx) = mLPNH(-npx:npx, npx)
  ! that completes the computation of the master pattern

  ! do we need to rebinarize?
  if (lmnl%binarize.eqv..TRUE.) then 
    where (mLPNH.gt.0.75) 
      mLPNH = 1.0
    end where 
    where (mLPSH.gt.0.75) 
      mLPSH = 1.0
    end where 
  end if 

  !=============================================
  !=============================================
  ! convert to stereographic projection
  allocate(masterSPNH(-npx:npx,-npy:npy))
  allocate(masterSPSH(-npx:npx,-npy:npy))
  masterSPNH = 0.0
  masterSPSH = 0.0


! get stereographic projections
  Radius = 1.0
  do i=-npx,npx 
    do j=-npx,npx 
      L = Lambert_T( xy = (/ float(i), float(j) /) / float(npx ))
      ierr = L%StereoGraphicInverse( xyz, Radius )
      xyz = xyz/vecnorm(xyz)
      if (ierr.ne.0) then 
        masterSPNH(i,j) = 0.0
        masterSPSH(i,j) = 0.0
      else
        masterSPNH(i,j) = InterpolateLambert(xyz, mLPNH, npx)
        masterSPSH(i,j) = InterpolateLambert(xyz, mLPSH, npx)
      end if
    end do
  end do

  !=============================================
  !=============================================
  ! we save an image for the Northern hemisphere in stereographic projection

  ! output the master pattern as a tiff file 

  fname = EMsoft%generateFilePath('EMdatapathname',trim(lmnl%tiffname))
  TIFF_filename = trim(fname)

  ! allocate memory for image
  TIFF_nx = 2*npx+1
  TIFF_ny = 2*npx+1
  allocate(TIFF_image(TIFF_nx,TIFF_ny))

  ! fill the image with whatever data you have (between 0 and 255)
  ma = maxval(masterSPNH)
  write(*,*) 'maximum intensity = ', ma 

  TIFF_image = int(255 * (masterSPNH/ma))

  ! set up the image_t structure
  im = image_t(TIFF_image)
  if(im%empty()) call Message%printMessage("ComputeLaueMasterPattern","failed to convert array to image")

  ! create the file
  call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
  if(0.ne.iostat) then
    call Message%printMessage("failed to write image to file : "//iomsg)
  else  
    call Message%printMessage('image written to '//trim(TIFF_filename))
  end if 
  deallocate(TIFF_image)

  !=============================================
  !=============================================
  ! save everything to HDF5 file
    
  if (lmnl%outformat.eq.'LMP') then   ! regular master pattern HDF5 output
    hdfname = EMsoft%generateFilePath('EMdatapathname',trim(lmnl%hdfname))
    !=============================================
    ! create or update the HDF5 output file
    !=============================================
    call HDFnames%set_ProgramData(SC_Lauemaster)
    call HDFnames%set_NMLlist(SC_LaueNameList)
    call HDFnames%set_NMLfilename(SC_LauemasterNML)
  
    ! Open an existing file or create a new file using the default properties.
    hdferr =  HDF%createFile(hdfname)

    ! write the EMheader to the file
    datagroupname = trim(HDFnames%get_ProgramData())
    call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)
  
  ! open or create a namelist group to write all the namelist files into
    groupname = SC_NMLfiles
    hdferr = HDF%createGroup(groupname)
  
  ! read the text file and write the array to the file
    dataset = SC_LauemasterNML
    hdferr = HDF%writeDatasetTextFile(dataset, nmldeffile)
  
  ! leave this group
    call HDF%pop()
  
  ! create a namelist group to write all the namelist files into
    groupname = SC_NMLparameters
    hdferr = HDF%createGroup(groupname)
  
    call self%writeHDFNameList(HDF, HDFnames)
  
  ! leave this group
    call HDF%pop()
  
  ! then the remainder of the data in a EMData group
    groupname = SC_EMData
    hdferr = HDF%createGroup(groupname)
  
  ! create the Lauemaster group and add a HDF_FileVersion attribbute to it 
    hdferr = HDF%createGroup(datagroupname)
    !HDF_FileVersion = '4.0'
    !attributename = SC_HDFFileVersion
    !hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)
  
  ! now start writing the ouput arrays...
    dataset = SC_mLPNH
    hdferr = HDF%writeDatasetFloatArray(dataset, mLPNH, 2*npx+1, 2*npx+1)
  
    dataset = SC_mLPSH
    hdferr = HDF%writeDatasetFloatArray(dataset, mLPSH, 2*npx+1, 2*npx+1)
  
    dataset = SC_masterSPNH
    hdferr = HDF%writeDatasetFloatArray(dataset, masterSPNH, 2*npx+1, 2*npx+1)
  
    dataset = SC_masterSPSH
    hdferr = HDF%writeDatasetFloatArray(dataset, masterSPSH, 2*npx+1, 2*npx+1)
  
  ! and close the file
    call HDF%pop(.TRUE.)
  ! close the Fortran interface
    call closeFortranHDFInterface()


  end if

  end associate

  end subroutine Lauemaster_
  
  
  
  end module mod_Lauemaster