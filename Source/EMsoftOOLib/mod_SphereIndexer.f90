!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!*                                                                     *
!* Copyright (c) 2019-2023, De Graef Group, Carnegie Mellon University *
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


module mod_SphereIndexer

use mod_kinds
use mod_imageOPs      ! rescale pattern
use mod_Detector      ! back project pattern
use mod_DSHT          ! spherical harmonic transform pattern
use mod_SHcorrelator  ! correlate with master pattern

  ! an indexing result
  type IdxRes
    real(kind=dbl) :: qu(4) ! (w, x, y, z)
    real(kind=dbl) :: xc    ! cross correlation at w, x, y, z
  end type IdxRes

  !@brief: type to hold data for back projection from detector to sphere
  type SphereIndexer
    type   (ImageRescaler      )              :: rescaler         ! rescale EBSD pattern
    type   (InverseEBSDDetector)              :: detector         ! back project rescaled pattern
    type   (DiscreteSHT        )              :: transformer      ! spherical harmonic transform back projected pattern
    type   (SphereCorrelator   )              :: correlator       ! correlate with master pattern
    
    complex(kind=dbl           ), allocatable :: almMaster(:,:)   ! spectra of master pattern to index against
    complex(kind=dbl           ), allocatable :: almPat   (:,:)   ! work space to hold spectra of exerimental pattern
    real   (kind=dbl           ), allocatable :: nh(:,:), sh(:,:) ! work space (spherical grid)
    real   (kind=dbl           ), allocatable :: patBuf(:,:)      ! work space (rescaled EBSD pattern)
    real   (kind=dbl           )              :: sclFct           ! scale factor for EBSD patterns prior to back projection
    logical                                   :: cirMsk           ! true/false to use/not use circular mask
  contains
    !@brief                      : initialize a spherical indexer
    !@param bw     [IN] integer  : indexing bandwidth
    !@param sig    [IN] real     : sample tilt (degrees)
    !@param L      [IN] integer  : sample -> detector distance (microns)
    !@param thetac [IN] real     : detector tilt (degrees)
    !@param delta  [IN] real     : detector pixel size (microns)
    !@param nx     [IN] integer  : detector width (pixels) 
    !@param ny     [IN] integer  : detector height (pixels)
    !@param mLPNH  [IN] real(:,:): north hemisphere of master pattern
    !@param mLPSH  [IN] real(:,:): south hemisphere of master pattern
    !@param cir    [IN] logical  : should a circluar mask be used
    procedure :: init      => SphereIndexer_Init

    !@brief: clean up a spherical indexer
    procedure :: destroy   => SphereIndexer_Destroy

    !@brief: clean up resources automatically
    final     ::              SphereIndexer_Finalize ! just calls destroy w/ polymorphism

    !@brief                                   : back project a pattern onto a spherical grid
    !@param pat[IN   ] real( 0:numsx, 0:numsy): ebsd pattern to back project to sphere
    !@param xpc[IN   ] real                   : x pattern center (pixels)
    !@param ypc[IN   ] real                   : y pattern center (pixels)
    !@param ref[IN   ] logical                : should newton's method refinement be used
    procedure :: index => SphereIndexer_Index ! this could be overloaded to enable back projection of different types
  end type SphereIndexer

contains


  !@brief               : initialize a spherical indexer
  !@param this          : object to initialize
  !@param bw     integer: indexing bandwidth
  !@param sig    real   : sample tilt (degrees)
  !@param L      integer: sample -> detector distance (microns)
  !@param thetac real   : detector tilt (degrees)
  !@param delta  real   : detector pixel size (microns)
  !@param nx     integer: detector width (pixels) 
  !@param ny     integer: detector height (pixels)
  !@param npx    integer: master pattern dimension (pixels)
  !@param mLPNH  real   : north hemisphere of master pattern
  !@param mLPSH  real   : south hemisphere of master pattern
  !@param cir    [IN] logical  : should a circluar mask be used
  recursive subroutine SphereIndexer_Init(this, bw, sig, L, thetac, delta, nx, ny, npx, mLPNH, mLPSH, cir)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereIndexer_Init
    use mod_io
    use mod_global
  implicit none
    class    (SphereIndexer),INTENT(INOUT)              :: this
    integer  (kind=irg     ),INTENT(IN   )              :: bw
    real     (kind=dbl     ),INTENT(IN   )              :: sig
    real     (kind=dbl     ),INTENT(IN   )              :: L
    real     (kind=dbl     ),INTENT(IN   )              :: thetac
    real     (kind=dbl     ),INTENT(IN   )              :: delta
    integer  (kind=irg     ),INTENT(IN   )              :: nx
    integer  (kind=irg     ),INTENT(IN   )              :: ny
    integer  (kind=irg     ),INTENT(IN   )              :: npx
    real     (kind=dbl     ),INTENT(IN   )              :: mLPNH(-npx:npx,-npx:npx)
    real     (kind=dbl     ),INTENT(IN   )              :: mLPSH(-npx:npx,-npx:npx)
    logical                 ,INTENT(IN   )              :: cir

    type(IO_T)                                          :: Message
    integer  (kind=irg     )                            :: d, dMP, nSph, nDet
    character(fnlen        )                            :: layout
    real     (kind=dbl     )                            :: omega
    ! clean up an existing object
    call this%destroy()

    ! save input
    this%cirMsk = cir

    ! determine spherical grid size from bandwidth
    d = bw / 2 + 1

    ! allocate spherical grid space
    allocate(this%nh(-d:d,-d:d))
    allocate(this%sh(-d:d,-d:d))
    
    ! build transformer
    layout = 'legendre'
    call this%transformer%init(d, bw, layout)

    ! interpolate master pattern onto spherical grid
    if(.not.all(shape(mLPNH).eq.shape(mLPSH))) then
      call Message%printError('SphereIndexer_Init', 'north / south hemisphere of master pattern are different shapes')
    endif
    if(ubound(mLPNH, 1).ne.ubound(mLPNH, 2)) then
      call Message%printError('SphereIndexer_Init', 'master pattern hemispheres must be square')
    endif
    if(lbound(mLPNH, 1).ne.lbound(mLPNH, 2)) then
      call Message%printError('SphereIndexer_Init', 'master pattern hemispheres must be square')
    endif
    if(abs(lbound(mLPNH, 1)).ne.ubound(mLPNH, 1)) then
      call Message%printError('SphereIndexer_Init', 'master pattern hemispheres must be symmetrically indexed about 0')
    endif
    dMP = ubound(mLPNH, 1)
    call LegendreInterp(dMP, mLPNH, mLPSH, d, this%nh, this%sh)

    ! compute spherical harmonic transform of master pattern
    allocate(this%almMaster(0:bw-1, 0:bw-1))
    allocate(this%almPat   (0:bw-1, 0:bw-1))
    call this%transformer%analyze(this%nh, this%sh, this%almMaster)

    ! build inverse detector using the original detector
    call this%detector%init(L, sig, thetac, d, delta, nx, ny)

    ! build a rescaler so images can be FFT rescaled before back projection
    omega = this%detector%solidAngle(0.D0, 0.D0) ! estimate fractional solid angle of EBSD detector
    nSph = nint(omega * (8 * d * d + 2) ) ! number of spherical pixels the EBSD detector covers
    nDet = nx * ny ! number of pixels on pattern
    this%sclFct = sqrt( real(nSph) * sqrt(2.D0) / real(nDet) ) ! this is the scale factor so that the pixel counts match
    call this%rescaler%init(nx, ny, this%sclFct)

    ! update inverse detector size
    call this%detector%resize(this%rescaler%wOut, this%rescaler%hOut)

    ! build correlator
    call this%correlator%init(bw)

    ! finally allocate space for rescaler output
    allocate(this%patBuf(this%rescaler%wOut,this%rescaler%hOut))

  end subroutine SphereIndexer_Init

  !@brief     : clean up an inverse EBSD detector type
  !@param this: SphereIndexer to clean up
  recursive subroutine SphereIndexer_Destroy(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereIndexer_Destroy
  implicit none
    class(SphereIndexer),INTENT(INOUT) :: this

    ! clean up working memory
    if (allocated(this%nh       ).eqv..TRUE.) deallocate(this%nh       )
    if (allocated(this%sh       ).eqv..TRUE.) deallocate(this%sh       )
    if (allocated(this%almMaster).eqv..TRUE.) deallocate(this%almMaster)
    if (allocated(this%almPat   ).eqv..TRUE.) deallocate(this%almPat   )
    if (allocated(this%patBuf   ).eqv..TRUE.) deallocate(this%patBuf   )

    ! type based members should clean themselves up automatically
    call this%transformer%destroy()
    call this%detector%destroy()
    call this%correlator%destroy()
    call this%rescaler%destroy()
  end subroutine SphereIndexer_Destroy

  !@brief     : clean up an inverse EBSD detector type
  !@param this: SphereIndexer to clean up
  recursive subroutine SphereIndexer_Finalize(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereIndexer_Finalize
  implicit none
    type (SphereIndexer),INTENT(INOUT) :: this
    call this%destroy()
  end subroutine SphereIndexer_Finalize

  !@brief     : clean up an inverse EBSD detector type
  !@param this: SphereIndexer to clean up
  recursive function SphereIndexer_Index(this, pat, xpc, ypc, ref) result(res)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereIndexer_Index
    use mod_Wigner
  implicit none
    class  (SphereIndexer),INTENT(INOUT) :: this
    real   (kind=dbl     ),INTENT(in   ) :: pat(:,:)
    real   (kind=dbl     ),INTENT(in   ) :: xpc
    real   (kind=dbl     ),INTENT(in   ) :: ypc
    type   (IdxRes       )               :: res

    real   (kind=dbl     )               :: eu(0:2), eps = 1.0D-2
    logical                              :: fMr, ref
    integer(kind=irg     )               :: fNf

    ! rescale pattern
    call this%rescaler%rescale(pat, this%patBuf, .true.) ! + remove DC background

    ! back project pattern to spherical grid
    this%nh = 0.D0
    this%sh = 0.D0
    call this%detector%unproject(this%patBuf, this%nh, this%sh, xpc * this%sclFct,  ypc * this%sclFct, this%cirMsk)

    ! compute SHT of pattern
    call this%transformer%analyze(this%nh, this%sh, this%almPat)

    ! cross correlate experimental and master pattern spectra
    fNf = 1
    fMr = .false.
    res%xc = this%correlator%correlate(this%almMaster, this%almPat, fMr, fNf, eu, ref, eps)

    ! convert result to quaternion
    res%qu = Wigner_zyz2qu(eu)


  end function SphereIndexer_Index


end module mod_SphereIndexer
