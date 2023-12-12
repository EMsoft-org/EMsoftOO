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


module mod_Detector

use mod_kinds

  !@brief: type to hold data for back projection from detector to sphere
  type InverseEBSDDetector
    ! all of these should be thread safe (read only)
    integer(kind=irg)             :: numsx, numsy ! EBSD detector size in pixels
    integer(kind=irg)             :: d            ! spherical grid half size - RQ(2,-d:d,-d:d)
    real   (kind=dbl),allocatable :: RQ(:,:,:)    ! position on the detector at spherical grid position i,j as [0,1] is 
                                                  ! (xpc/numsx, ypc/numsy) + RQ(:,i,j) + 0.5
  contains
    !@brief                    : initialize an inverse EBSD detector from geometry and spherical grid
    !@param L      [IN] integer: sample -> detector distance (microns)
    !@param sig    [IN] real   : sample tilt (degrees)
    !@param thetac [IN] real   : detector tilt (degrees)
    !@param d      [IN] integer: spherical grid half size - grid(-d:d)
    !@param delta  [IN] real   : detector pixel size (microns)
    !@param nx     [IN] integer: detector width (pixels) 
    !@param ny     [IN] integer: detector height (pixels)
    procedure :: init       => InverseEBSDDetector_Init

    !@brief                 : rescale a detector (~continuous binning)
    !@brief wNew[IN] integer: new detector width  in pixels (width  in microns is unchanged)
    !@brief hNew[IN] integer: new detector height in pixels (height in microns in unchanged)
    procedure :: resize     => InverseEBSDDetector_Resize

    !@brief: clean up an inverse EBSD detector type
    procedure :: destroy    => InverseEBSDDetector_Destroy

    !@brief: clean up resources automatically
    final     ::              InverseEBSDDetector_Finalize ! just calls destroy w/ polymorphism

    !@brief                                   : back project a pattern onto a spherical grid
    !@param pat[IN   ] real( 0:numsx, 0:numsy): ebsd pattern to back project to sphere
    !@param nh [INOUT] real(-d:d    ,-d:d    ): northern hemisphere to back project onto
    !@param sh [INOUT] real(-d:d    ,-d:d    ): northern hemisphere to back project onto
    !@param xpc[IN   ] real                   : x pattern center (pixels)
    !@param ypc[IN   ] real                   : y pattern center (pixels)
    !@param cir[IN   ] logical                : true/false to restrict back projection to circular mask
    procedure :: unproject  => InverseEBSDDetector_Unproject

    !@brief                : estimate the solid angle covered by this detector
    !@param xpc[IN   ] real: x pattern center (pixels)
    !@param ypc[IN   ] real: y pattern center (pixels)
    !@return               : fractional solid angle (i.e. [0,1] instead of [0,4*pi])
    procedure :: solidAngle => InverseEBSDDetector_SolidAngle
  end type InverseEBSDDetector

contains
  !@brief         : wrapper to initialize an inverse EBSD detector (used to be from namelists)
  !@param sinl    : spherical indexing namelist [removed]
  !@param mcnl    : monte carlo namelist        [removed]
  !@param L       : detector distance
  !@param sig     : sample tilt angle (degrees)
  !@param thetac  : detector tilt angle (degrees)
  !@param bw      : bandwidth parameter 
  !@param delta   : detector pixel size (microns)
  !@param patdims : detector dimensions (2D array)
  !@param detector: detector to initialize
  ! recursive subroutine GenerateInverseEBSDdetector(sinl, mcnl, detector)
  recursive subroutine GenerateInverseEBSDdetector(L, sig, thetac, bw, delta, patdims, detector)
  !DEC$ ATTRIBUTES DLLEXPORT :: GenerateInverseEBSDdetector
    ! use NameListTypedefs
    ! use sphtypedefs
  implicit none
    ! type   (SphInxNameListType ),INTENT(INOUT) :: sinl
    ! type   (MCCLNameListType   ),INTENT(INOUT) :: mcnl
    real(kind=dbl),INTENT(IN)                  :: L 
    real(kind=dbl),INTENT(IN)                  :: sig
    real(kind=dbl),INTENT(IN)                  :: thetac 
    integer(kind=irg),INTENT(IN)               :: bw
    real(kind=dbl),INTENT(IN)                  :: delta
    integer(kind=irg),INTENT(IN)               :: patdims(2)
    type   (InverseEBSDDetector),INTENT(INOUT) :: detector
    call detector%init(L, sig, thetac, bw + 1, delta, patdims(1), patdims(2)) ! initialize object
  end subroutine GenerateInverseEBSDdetector

  !@brief       : initialize an inverse EBSD detector from geometrry and spherical grid
  !@param this  : object to initialize
  !@param L     : sample -> detector distance (microns)
  !@param sig   : sample tilt (degrees)
  !@param thetac: detector tilt (degrees)
  !@param d     : spherical grid half size - grid(-d:d)
  !@param delta : detector pixel size (microns)
  !@param nx    : detector width (pixels) 
  !@param ny    : detector height (pixels)
  recursive subroutine InverseEBSDDetector_Init(this, L, sig, thetac, d, delta, nx, ny)
  !DEC$ ATTRIBUTES DLLEXPORT :: InverseEBSDDetector_Init
    use mod_DSHT ! SH_cosLats
    use mod_io
    use mod_global
  implicit none
    class    (InverseEBSDDetector),INTENT(INOUT) :: this
    real     (kind=dbl           ),INTENT(IN   ) :: L
    real     (kind=dbl           ),INTENT(IN   ) :: sig
    real     (kind=dbl           ),INTENT(IN   ) :: thetac
    integer  (kind=irg           ),INTENT(IN   ) :: d
    real     (kind=dbl           ),INTENT(IN   ) :: delta
    integer  (kind=irg           ),INTENT(IN   ) :: nx
    integer  (kind=irg           ),INTENT(IN   ) :: ny

    real     (kind=dbl           )               :: alp, ca, sa
    real     (kind=dbl           )               :: dc(3,-d:d,-d:d)
    character(fnlen              )               :: t = 'legendre'

    ! save needed constants
    this%numsx = nx
    this%numsy = ny
    this%d     = d

    ! allocate memory
    allocate(this%RQ(2,-d:d,-d:d)) ! Q(1) and R(2) for each i/j in spherical grid

    ! compute some intermediate values once
    dc = SH_dirCos(d, t) ! direction cosines
    alp = 0.5 * cPi - (sig - thetac) * cPi / 180.D0 ! angle between detector and sample plane in radians
    write (*,*) 'Inverse detector alpha angle = ', alp/dtor
    ca = cos(alp)
    sa = sin(alp)

    ! compute Q and R arrays from direction cosines (north hemisphere)
    this%RQ(2,:,:) = L / ( dc(1,:,:) * sa + dc(3,:,:) * ca ) ! this is 'd' in the c++ code

    this%RQ(1,:,:) =                         dc(2,:,:)   * this%RQ(2,:,:) ! R=c_y*d=c_y*L/(c_x*sin(alpha)+c_z*cos(alpha) ) 
    this%RQ(2,:,:) = ( sa * dc(3,:,:) - ca * dc(1,:,:) ) * this%RQ(2,:,:) ! Q=(sin(alpha)*c_z-cos(alpha)*c_x)*d

    ! we now have RQ in microns, convert to fractional detectors
    this%RQ(1,:,:) = this%RQ(1,:,:) / (delta * nx) ! convert to fractional detector
    this%RQ(2,:,:) = this%RQ(2,:,:) / (delta * ny) ! convert to fractional detector
  end subroutine InverseEBSDDetector_Init

  !@brief     : rescale a detector (~continuous binning)
  !@param this: object to resize
  !@param wNew: new detector width  in pixels (width  in microns is unchanged)
  !@param hNew: new detector height in pixels (height in microns in unchanged)
  recursive subroutine InverseEBSDDetector_Resize(this, wNew, hNew)
  !DEC$ ATTRIBUTES DLLEXPORT :: InverseEBSDDetector_Resize
  implicit none
    class    (InverseEBSDDetector),INTENT(INOUT) :: this
    integer  (kind=irg           ),INTENT(IN   ) :: wNew
    integer  (kind=irg           ),INTENT(IN   ) :: hNew
    this%numsx = wNew
    this%numsy = hNew
  end subroutine InverseEBSDDetector_Resize

  !@brief     : clean up an inverse EBSD detector type
  !@param this: InverseEBSDDetector to clean up
  recursive subroutine InverseEBSDDetector_Destroy(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: InverseEBSDDetector_Destroy
  implicit none
    class(InverseEBSDDetector),INTENT(INOUT) :: this
    if (allocated(this%RQ).eqv..TRUE.) deallocate(this%RQ)
  end subroutine InverseEBSDDetector_Destroy

  !@brief     : clean up an inverse EBSD detector type
  !@param this: InverseEBSDDetector to clean up
  recursive subroutine InverseEBSDDetector_Finalize(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: InverseEBSDDetector_Finalize
  implicit none
    type (InverseEBSDDetector),INTENT(INOUT) :: this
    call this%destroy()
  end subroutine InverseEBSDDetector_Finalize

  !@brief     : back project a pattern onto a spherical grid
  !@param this: initialized InverseEBSDDetector
  !@param pat : ebsd pattern to back project to sphere
  !@param nh  : northern hemisphere to back project onto
  !@param sh  : northern hemisphere to back project onto
  !@param xpc : x pattern center (pixels)
  !@param ypc : y pattern center (pixels)
  !@param cir : true/false to restrict back projection to circular mask
  !@param flip: optional parameter to flip the patterns upside down
  recursive subroutine InverseEBSDDetector_Unproject(this, pat, nh, sh, xpc, ypc, cir, flip)
  !DEC$ ATTRIBUTES DLLEXPORT :: InverseEBSDDetector_Unproject
  implicit none
    class  (InverseEBSDDetector),INTENT(IN   ) :: this
    real   (kind=dbl           ),INTENT(IN   ) :: pat(:,:) ! the code currently assumes that this is 0:numsx-1,0:numsy-1
    real   (kind=dbl           ),INTENT(INOUT) :: nh(-this%d:this%d, -this%d:this%d)
    real   (kind=dbl           ),INTENT(INOUT) :: sh(-this%d:this%d, -this%d:this%d)
    real   (kind=dbl           ),INTENT(IN   ) :: xpc
    real   (kind=dbl           ),INTENT(IN   ) :: ypc
    logical                     ,INTENT(IN   ) :: cir
    logical                     ,INTENT(IN   ),OPTIONAL :: flip

    real   (kind=dbl           )               :: x, y, fx, fy, cx, cy, v, mean, rx, ry, rr, hx, hy, hh, px, py, stdv
    integer(kind=irg           )               :: i, j, k, x0, y0, x1, y1, cnt
    logical                                    :: flipit = .FALSE.

    if (present(flip)) then
      if (flip.eqv..TRUE.) flipit = .TRUE.
    end if 

    ! start by filling sphere with 0
    nh = 0.D0
    sh = 0.D0

    ! compute radius of largest inscribed circle on detector
    hx = real(this%numsx) / 2.D0
    hy = real(this%numsy) / 2.D0
    hh = min(hx, hy)

    ! now loop over spherical grid computing position on detector and building mean of back projected balue
    cnt = 0
    mean = 0.D0
    do j = -this%d, this%d ! loop over first spherical grid direction
      do i = -this%d, this%d ! loop over second spherical grid direction
        ! compute position on EBSD detector for this spherical grid point (fractional detectors)
        x = xpc / this%numsx + this%RQ(1,i,j) ! [-0.5, 0.5]
        y = ypc / this%numsy + this%RQ(2,i,j) ! [-0.5, 0.5]

        ! check if inside circular mask if needed
        if(cir) then
          ! convert from [-0.5,0.5] in independant direction to [-1,1] in smaller direction
          rx = x * 2.0 * hx / hh
          ry = y * 2.0 * hy / hh
          rr = rx * rx + ry * ry
          if(rr.gt.1.D0) cycle ! skip points outside circular mask
        endif

        ! move origin from detector center to corner
        x = x + 0.5
        if(flipit.eqv..true.) y = - y ! flips pattern
        y = y + 0.5 ! [0, 1]
        ! if(.true.) y = 1.D0 - y ! flips pattern

        ! check if this spherical grid point projects onto the detector
        if(x.ge.0.D0.and.y.ge.0.D0.and.x.lt.1.D0.and.y.lt.1.D0) then ! on detector
          ! don't go out of bounds on the detector
          px = x * (this%numsx-1) + 1 ! continous pixel coordinate in [1,numsx]
          py = y * (this%numsy-1) + 1 ! continous pixel coordinate in [1,numsy]
          x0 = int(px) ! left bounding column
          y0 = int(py) ! left bounding column
          x1 = x0 + 1  ! right bounding column
          y1 = y0 + 1  ! top bounding row
          if(x1.gt.this%numsx) x1 = x0 ! handle exact edge
          if(y1.gt.this%numsy) y1 = y0 ! handle exact edge

          ! we're on the detector, bilinearly interpolate from the 4 bounding pixels
          fx = px - x0 ! fractional progress between pixels in x direction
          fy = py - y0 ! fractional progress between pixels in y direction
          cx = 1.D0 - fx  ! complement of fx
          cy = 1.D0 - fy  ! complement of fy

          v = cx * cy * pat(x0, y0) &
            + fx * cy * pat(x1, y0) &
            + cx * fy * pat(x0, y1) &
            + fx * fy * pat(x1, y1)

          ! save the value
          nh(i,j) = v
          mean = mean + v ! this should technically be weighted by the pixel area
          cnt  = cnt  + 1
        endif
      end do ! i
    end do ! j

    

    ! subtract mean and compute stdev
    ! this isn't quite the same as the C++ version which uses a proper area weighted calculation
    mean = mean / cnt
    stdv = 0.D0
    do j = -this%d, this%d ! loop over first spherical grid direction
      do i = -this%d, this%d ! loop over second spherical grid direction
        ! compute position on EBSD detector for this spherical grid point (fractional detectors)
        x = xpc / this%numsx + this%RQ(1,i,j) ! [-0.5, 0.5]
        y = ypc / this%numsy + this%RQ(2,i,j) ! [-0.5, 0.5]

        ! check if inside circular mask if needed
        if(cir) then
          ! convert from [-0.5,0.5] in independant direction to [-1,1] in smaller direction
          rx = x * 2.0 * hx / hh
          ry = y * 2.0 * hy / hh
          rr = rx * rx + ry * ry
          if(rr.gt.1.D0) cycle ! skip points outside circular mask
        endif

        ! move origin from detector center to corner
        x = x + 0.5
        y = y + 0.5

        ! check if this spherical grid point projects onto the detector
        if(x.ge.0.D0.and.y.ge.0.D0.and.x.lt.1.D0.and.y.lt.1.D0) then ! on detector
          nh(i,j) = nh(i,j) - mean ! subtract mean
          v = nh(i,j)
          stdv = stdv + v * v
        endif
      end do ! i
    end do ! j

    stdv = sqrt(stdv / cnt)
    nh = nh / stdv

  end subroutine InverseEBSDDetector_Unproject

  !@brief    : estimate the solid angle covered by this detector
  !@param xpc: x pattern center (pixels)
  !@param ypc: y pattern center (pixels)
  !@return   : fractional solid angle (i.e. [0,1] instead of [0,4*pi])
  recursive function InverseEBSDDetector_SolidAngle(this, xpc, ypc) result(omg)
  !DEC$ ATTRIBUTES DLLEXPORT :: InverseEBSDDetector_SolidAngle
  implicit none
    class  (InverseEBSDDetector),INTENT(IN   ) :: this
    real   (kind=dbl           ),INTENT(IN   ) :: xpc
    real   (kind=dbl           ),INTENT(IN   ) :: ypc
    real   (kind=dbl           )               :: omg
    integer                                    :: cnt, i, j, k, x0, y0
    real   (kind=dbl           )               :: x, y

    ! loop over spherical grid counting point that hit the detector
    cnt = 0
    do j = -this%d, this%d ! loop over first spherical grid direction
      do i = -this%d, this%d ! loop over second spherical grid direction

        ! compute position on EBSD detector for this spherical grid point (fractional detectors)
        x = xpc / this%numsx + this%RQ(1,i,j) + 0.5
        y = ypc / this%numsy + this%RQ(2,i,j) + 0.5

        ! check if this spherical grid point projects onto the detector
        if(x.ge.0.D0.and.y.ge.0.D0.and.y.le.1.D0.and.y.le.1.D0) then
          cnt = cnt + 1 ! inside detector
        endif
      enddo
    enddo

    ! solid angle is ~ fraction of pixels that hit the detector
    ! this is a pretty rough approximation for square legendre grids but should be good for square lambert
    omg = real(cnt) / real( 8 * this%d * this%d + 2 )
  end function InverseEBSDDetector_SolidAngle

end module mod_Detector
