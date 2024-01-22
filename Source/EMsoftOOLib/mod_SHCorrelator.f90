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
! EMsoft:SHCorrelator.f90
!--------------------------------------------------------------------------
!
! MODULE: SHCorrelator
!
!> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
!
!> @brief Spherical Harmonic cross-correlation routines 
!
!> @details 
!>  reference: Gutman, B., Wang, Y., Chan, T., Thompson, P. M., & Toga, A. W. (2008, October). Shape registration with 
!>     spherical cross correlation. In 2nd MICCAI workshop on mathematical foundations of computational anatomy (pp. 56-67).
!>  note: the reference is for the generic case (complex functions) but restricting to real valued functions allow for 
!>    savings from symmetry:
!>        \hat{f}(l, -m) = \hat{f}(l, m) * (-1)^m
!>  additionally since the decomposition of the Wigner D function to 2 Wigner d functions @ pi/2 introduces more symmetry:
!>        d^j_{-k,-m} = (-1)^(   k- m) d^j_{k,m}
!>        d^j_{ k,-m} = (-1)^(j+ k+2m) d^j_{k,m}
!>        d^j_{-k, m} = (-1)^(j+2k+3m) d^j_{k,m}
!>        d^j_{ m, k} = (-1)^(   k- m) d^j_{k,m}
!>  finally since the cross correlation of the real functions is also real there is another factor of 2 savings
!
!> @note: variable names are consistent with Gutman et. al. (see eq 12 for details) except 'j' is used in place of 'l'
! 
!> @date 01/29/19 MDG 1.0 original, based on Will Lenthe's C++ routines (sht_xcorr.hpp)
!--------------------------------------------------------------------------
module mod_SHcorrelator

use mod_kinds
use mod_fft_wrap

implicit none

private ! everything is private by default

public :: SphereXcConstants, SphereCorrelator

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                                    !!
  !!  these are the thread safe components required for a spherical XC  !!
  !!                                                                    !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type SphereXcConstants
    integer(kind=irg )             :: bw             ! maximum bandwidth to use (exclusive)
    integer(kind=sgl )             :: sl             ! side length of grid in euler space (2 * bandWidth - 1)
    integer(kind=irg )             :: bwP            ! bandwidth of padded sidelength
    integer(kind=sgl )             :: slP            ! sidelength padded to a fast FFT size
    real   (kind=dbl ),allocatable :: wigD (:,:,:  ) ! lookup table for wigner d functions
    type   (Real3DFFT)             :: plan           ! inverse 3D fft plan
  contains
    !@brief                : initialize a spherical cross correlation constants type
    !@param bw [IN] integer: bandwidth
    procedure :: init    => SphereXcConstants_Init

    !@brief: clean up a spherical cross correlation constants type
    procedure :: destroy => SphereXcConstants_Destroy

    !@brief: clean up resources automatically
    final     ::            SphereXcConstants_Finalize ! just calls destroy w/ polymorphism
  end type SphereXcConstants

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                                    !!
  !!       this is the main type to do spherical correlation with       !!
  !!                                                                    !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type SphereCorrelator
    ! these are the non-thread safe components required for a spherical cross correlation
    complex(kind=dbl         ),allocatable :: fm   (:,:    ) ! 2d lookup table to hold \hat{f}(j,m) * d^j_{m, k} for all j and m (for any given k)
    complex(kind=dbl         ),allocatable :: gn   (:      ) ! 1d lookup table to hold \bar{\hat{g}(j,n)} * d^j_{k,n} for all j (for any given k and n)
    ! type   (C_PTR            )             :: pFxc, pXc      ! c pointers for allocation of fxc and xc
    type   (FFTBuffer        )             :: pFxc, pXc      ! fft allocated arrays for pointers
    complex(C_DOUBLE_COMPLEX ),pointer     :: fxc  (:,:,:  ) ! fft of cross correlation (in half complex format)
    real   (C_DOUBLE         ),pointer     :: xc   (:,:,:  ) ! real space cross correlation (this is what we're after)
    real   (kind=dbl         ),allocatable :: dBeta(:,:,:,:) ! wigner (lowercase) d lookup table for arbitrary beta (for refinement)

    ! these are the thread safe (read only) components required for a spherical cross correlation
    type   (SphereXcConstants),allocatable :: xcLut
  contains
    !@brief                : initialize a spherical cross correlation calculator
    !@param bw [IN] integer: bandwidth
    procedure          :: init        => SphereCorrelator_Init

    !@brief: clean up a spherical cross correlation calculator
    procedure          :: destroy     => SphereCorrelator_Destroy

    !@brief: clean up resources automatically
    final              ::                SphereCorrelator_Finalize ! just calls destroy w/ polymorphism

    !@brief                                    : compute spherical cross correlation (result will be in this%xc)
    !@param flm [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the first function
    !@param gln [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the second function
    !@param fMr [IN   ] logical                : true/false if there is/isn't a mirror plane in the first function
    !@param fNf [IN   ] integer                : rotational symmetry about z axis in first function (1 for no rotational symmetry)
    !@param eu  [INOUT] real   (0:2)           : location to write rotation of maximum cross correlation as ZYZ euler angle
    !@param ref [IN   ] logical                : true/false to use/not use real space refinement
    !@param eps [IN   ] real                   : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
    !@return            real                   : maximum cross correlation
    procedure          :: correlate   => SphereCorrelator_Correlate ! this calls compute, findPeak, refinePeak underneath

    !@brief                                    : refine an orientation with newton's method
    !@param flm [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the first function
    !@param gln [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the second function
    !@param fMr [IN   ] logical                : true/false if there is/isn't a mirror plane in the first function
    !@param fNf [IN   ] integer                : rotational symmetry about z axis in first function (1 for no rotational symmetry)
    !@param eu  [INOUT] real   (0:2)           : orientation to refine
    !@param eps [IN   ] real                   : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
    !@return            real                   : maximum cross correlation
    procedure          :: refinePeak  => SphereCorrelator_RefinePeak ! not yet implemented

    ! start of private functions

    !@brief                                    : compute spherical cross correlation (result will be in this%xc)
    !@param flm [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the first function
    !@param gln [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the second function
    !@param fMr [IN   ] logical                : true/false if there is/isn't a mirror plane in the first function
    !@param fNf [IN   ] integer                : rotational symmetry about z axis in first function (1 for no rotational symmetry)
    procedure, private :: compute     => SphereCorrelator_Compute

    !@brief : get the maximum value and location of this%xc 
    !@return: index of maximum
    !@note  : this is trivial now but is implemented to support normalized correlation in the future
    procedure, private :: findPeak    => SphereCorrelator_FindPeak

    !@brief: compute the first and second derivatives of the cross correlation at a single rotation
    ! @param flm [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the first function
    ! @param gln [IN   ] cmplx(0:bw-1, 0:bw-1): spherical harmonic coefficients for the second function
    ! @param eu  [INOUT] real   (0:2)           : rotation to compute derivatives of cross correlation for as ZYZ euler angle
    ! @param jac [INOUT] real   (0:2)           : location to write jacobian of cross correlation {d/(d eu[0]), d/(d eu[1]), d/(d eu[2])}
    ! @param hes [INOUT] real   (0:8)           : location to write hessian (3x3 matrix as 9 component vector) of cross correlation hes_ij = d/(d eu[i]) * d/(d eu[j])
    ! @param mBW [IN   ] integer                : maximum bandwidth to use in calculation (must be <= bw)
    ! @param fMr [IN   ] logical                : true/false if there is/isn't a mirror plane in the first fuction
    ! @param fNf [IN   ] integer                : rotational symmetry about z axis in first function (1 for no rotational symmetry)
    ! @return    [IN   ] real                   : maximum cross correlation
    procedure, private :: derivatives => SphereCorrelator_Derivatives

  end type SphereCorrelator

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                     SphereXcConstants Members                      !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SphereXcConstants_Init
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief Initialize spherical cross correlation constants
  !
  !> @param this structure to initialize
  !> @param bw maximum bandwidth of spherical harmonic to use in correlation
  !
  !> @date 07/22/19 WCL 1.0 original
  !--------------------------------------------------------------------------
  recursive subroutine SphereXcConstants_Init(this, bw)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereXcConstants_Init
    use mod_Wigner
  implicit none
    class  (SphereXcConstants),INTENT(INOUT) :: this       ! structure to initialize
    integer                   ,INTENT(IN   ) :: bw         ! bandwidth to initialize with

    ! clean up an existing object
    call this%destroy()

    this%bw  = bw               ! save bandwidth
    this%sl  = bw * 2 - 1       ! save sidelength (-bw:bw)
    this%slP = this%sl          ! for now don't zero pad up (if there is a function to compute fast FFT sizes use it here)
    this%bwP = this%slP / 2 + 1 ! compute the padded bandwidth

    ! fill wigner d lookup table
    allocate(this%wigD (0:bw-1, 0:bw-1, 0:bw-1)) ! wigner d lookup table for beta = 0 (over all j, k, m)
    call Wigner_dTable1(bw, this%wigD)

    ! build fft plan (for now we'll just to the naive inverse fft)
    call this%plan%init(this%slP, FFT_PLAN_MEASURE)
  end subroutine SphereXcConstants_Init

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SphereXcConstants_Destroy
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief clean up spherical cross correlation constants
  !
  !> @param this structure to destroy
  !
  !> @date 07/22/19 WCL 1.0 original
  !--------------------------------------------------------------------------
  recursive subroutine SphereXcConstants_Destroy(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereXcConstants_Destroy
  implicit none
    class(SphereXcConstants),INTENT(INOUT) :: this ! structure to initialize
    if (allocated(this%wigD).eqv..TRUE.) deallocate(this%wigD) ! free wigner d lookup table
  end subroutine SphereXcConstants_Destroy

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SphereXcConstants_Finalize
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief clean up spherical cross correlation constants
  !
  !> @param this structure to destroy
  !
  !> @date 07/22/19 WCL 1.0 original
  !--------------------------------------------------------------------------
  recursive subroutine SphereXcConstants_Finalize(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereXcConstants_Finalize
  implicit none
    type(SphereXcConstants),INTENT(INOUT) :: this ! structure to initialize
    call this%destroy()
  end subroutine SphereXcConstants_Finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                       SHTCorrelator Members                        !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SphereCorrelator_Init
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief Initialize a spherical cross correlation calculator
  !
  !> @param this structure to initialize
  !> @param bw maximum bandwidth of spherical harmonic to use in correlation
  !
  !> @date 07/22/19 WCL 1.0 original
  !--------------------------------------------------------------------------
  recursive subroutine SphereCorrelator_Init(this, bw)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_Init
  implicit none
    class  (SphereCorrelator),INTENT(INOUT) :: this     ! structure to initialize
    integer                  ,INTENT(IN   ) :: bw       ! bandwidth to initialize with

    integer                                 :: bwP, slP ! since full name with xcLut is pretty long

    ! clean up an existing object
    call this%destroy()

    ! start by filling in read only constants
    allocate(this%xcLut)
    call this%xcLut%init(bw)

    ! allocate fortran work space
    allocate(this%fm   (0:bw-1, 0:bw-1               )) ! working space
    allocate(this%gn   (0:bw-1                       )) ! working space

    allocate(this%dBeta(0:1   ,0:bw-1, 0:bw-1, 0:bw-1)) ! wigner d lookup table for beta != 0 (for refinement)

    ! allocate fft arrays
    bwP = this%xcLut%bwP
    slP = this%xcLut%slP
    call this%pFxc%allocCplx(slP * slP * bwP)
    call this%pXc %allocReal(slP * slP * slP)

    ! wrap pointers for fortran access
    call c_f_pointer(this%pFxc%ptr, this%fxc, [bwP, slP, slP])
    call c_f_pointer(this%pXc %ptr, this%xc , [slP, slP, slP])
  end subroutine SphereCorrelator_Init

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SphereCorrelator_Destroy
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief cleanup a spherical cross correlation calculator 
  !
  !> @param this structure to destroy
  !
  !> @date 07/22/19 WCL 1.0 original
  !--------------------------------------------------------------------------
  recursive subroutine SphereCorrelator_Destroy(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_Destroy
  implicit none
    class(SphereCorrelator),INTENT(INOUT) :: this ! structure to destroy

    ! clean up lookup table
    if(allocated(this%xcLut)) then
      call this%xcLut%destroy()
      deallocate(this%xcLut)
    endif

    ! clean up fortran tables
    if (allocated(this%fm   ).eqv..TRUE.) deallocate(this%fm   )
    if (allocated(this%gn   ).eqv..TRUE.) deallocate(this%gn   )
    if (allocated(this%dBeta).eqv..TRUE.) deallocate(this%dBeta)

  end subroutine SphereCorrelator_Destroy

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SphereCorrelator_Finalize
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief clean up a spherical cross correlation calculator
  !
  !> @param this structure to destroy
  !
  !> @date 07/22/19 WCL 1.0 original
  !--------------------------------------------------------------------------
  recursive subroutine SphereCorrelator_Finalize(this)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_Finalize
  implicit none
    type(SphereCorrelator),INTENT(INOUT) :: this ! structure to destroy
    call this%destroy()
  end subroutine SphereCorrelator_Finalize

  !--------------------------------------------------------------------------
  !
  ! FUNCTION: SphereCorrelator_Correlate
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief   compute the cross correlation between two spherical functions
  !
  ! @param flm: spherical harmonic coefficients for the first function
  ! @param gln: spherical harmonic coefficients for the second function
  ! @param fMr: true/false if there is/isn't a mirror plane in the first function
  ! @param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
  ! @param eu : location to write rotation of maximum cross correlation as ZYZ euler angle
  ! @param ref: true/false to use/not use real space refinement
  ! @param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
  ! @return   : maximum cross correlation
  !
  !> @date 01/30/19 MDG 1.0 original, based on Will Lenthe's classes in sht_xcorr.hpp
  !--------------------------------------------------------------------------
  recursive function SphereCorrelator_Correlate(this, flm, gln, fMr, fNf, eu, ref, eps) result(peak)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_Correlate

    use mod_global 
    use mod_Wigner

  implicit none

    class  (SphereCorrelator),INTENT(INOUT) :: this
    complex(kind=dbl        ),INTENT(IN   ) :: flm(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    complex(kind=dbl        ),INTENT(IN   ) :: gln(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    logical                  ,INTENT(IN   ) :: fMr 
    integer(kind=irg        ),INTENT(IN   ) :: fNf
    real   (kind=dbl        ),INTENT(INOUT) :: eu(0:2)
    logical                  ,INTENT(IN   ) :: ref
    real   (kind=dbl        ),INTENT(IN   ) :: eps
    real   (kind=dbl        )               :: peak

    integer                                 :: indMax(3)                                ! location of maximum cross correlation
    real   (kind=dbl        )               :: nhood(-1:1,-1:1,-1:1), x(0:2), knm(0:2)  ! neighborhood around peak (for subpixel interpolation) + eu work space
    integer                                 :: k, n, m, kk, nn, mm                      ! periodic indices for neighborhood extraction
    logical                                 :: glide                                    ! did we cross the glide plane for neighborhood extraction

    ! compute grid of spherical cross correlation in this%xc
    call this%compute(flm, gln, fMr, fNf)

    ! find maximum cross correlation
    indMax = this%findPeak()

    ! extract 3x3 neighborhood around peak handling periodic boundary conditions
    do k = -1,1
        kk = indMax(3) + k ! compute neighbor index

        ! bring kk back to [1,bwP] if needed
        if(kk.lt.1) then
            kk = kk + this%xcLut%slP
            glide = .TRUE.
        elseif(kk.gt.this%xcLut%slP) then
            kk = kk - this%xcLut%slP
            glide = .TRUE.
        else
            glide = .FALSE.
        endif

        do n = -1,1
            nn = indMax(2) + n ! compute neighbor index

            ! bring nn back to [1,slP] if needed
            if(nn.lt.1) then
                nn = nn + this%xcLut%slP
            elseif(nn.gt.this%xcLut%slP) then
                nn = nn - this%xcLut%slP
            endif

            do m = -1,1
                mm = indMax(1) + m ! compute neighbor index
                ! bring mm back to [1,slP] if needed
                if(mm.lt.1) then
                    mm = mm + this%xcLut%slP
                elseif(mm.gt.this%xcLut%slP) then
                    mm = mm - this%xcLut%slP
                endif

                ! now extract neighbor
                nhood(k,n,m) = this%xc(mm,nn,kk)

            enddo ! m
        enddo ! n
    enddo ! k

    ! subpixel interpolate maxima
    x = 0.D0
    peak = SH_interpolateMaxima(nhood, x, 25)

    ! convert from fractional pixels to zyz euler angle
    knm = dble(indMax) - 1 ! convert from 1 to zero indexing
    eu(0) = ( (knm(0) + x(0) ) * 4 - this%xcLut%slP) * cPi / (this%xcLut%slP * 2) ! alpha
    eu(1) = ( (knm(2) + x(2) ) * 2 - this%xcLut%slP) * cPi / (this%xcLut%slP    ) ! beta
    eu(2) = ( (knm(1) + x(1) ) * 4 - this%xcLut%slP) * cPi / (this%xcLut%slP * 2) ! gamma

    ! do newton's method refinement if needed
    if(ref) then
      peak = this%refinePeak(flm, gln, fMr, fNf, eu, eps)
    endif

  end function SphereCorrelator_Correlate

  recursive function cholesky_solve(a, x, b, n) result(worked)
  !DEC$ ATTRIBUTES DLLEXPORT :: cholesky_solve
   implicit none
    real   (kind=dbl        ),INTENT(INOUT) :: a(0:n-1, 0:n-1) ! n by n symmetric matrix
    real   (kind=dbl        ),INTENT(INOUT) :: x(0:n-1)        ! n component vector (location to write result)
    real   (kind=dbl        ),INTENT(IN   ) :: b(0:n-1)        ! n component vector a * x = b
    integer(kind=irg        ),INTENT(IN   ) :: n               ! problem size
    logical                                 :: worked          ! return true on success, false otherwise
    logical                                 :: neg
    real                                    :: sum, d(0:n-1)
    integer(kind=irg        )               :: i, j, k
    real                                    :: eps = EPSILON(eps)

    ! sanity check dimensinos
    if(n.lt.1) then
      worked = .false.
      return
    endif

    ! check for negation
    if(a(0,0).lt.0.D0) then
      neg = .true.
    else
      neg = .false.
    endif

    ! start with in place decomposition
    do i = 0, n-1
      if( (a(i,i).lt.0.D0).neqv.neg) return ! not positive or negative definate
      do j = i, n-1
        sum = a(i,j) ! get A_{i,j} (or negative if needed)
        if(neg) sum = -sum
        do k = 0, i-1
          sum = sum - a(i,k) * a(j,k) ! accumulate A_{i,j} - \sum_{k=0}^{j-1} L_{i,k} * L^*_{j,k}
        enddo ! k
        if(i.eq.j) then ! we're computing L_{j,j} (save in d instead of overwriting elements of A)
          ! sum will always be real here since A_{i,i} is real and L_{i,k} * L^*_{j,k} -> L_{j,k} * L^*_{j,k} for i == j
          if(sum.lt.eps) then
            return ! need to use LDL instead (which is why we've kept the diagonal intact)
          endif
          d(i) = sqrt(sum) ! L_{j,j} = \sqrt{ A_{j,j} - \sum_{k=0}^{j-1} L_{j,k} * L_{j,k} }
        else ! we're computing L_{i,j}
          a(j,i) = sum / d(i) ! L_{i,j} = ( A_{i,j} - \sum_{k=0}^{j-1} L_{i,k} * L_{j,k} ) / L_{j,j}
        endif
      enddo ! j
    enddo ! i
    worked = .true. ! success

    ! now backsolve for x
    do i = 0, n-1
      x(i) = ( b(i) - dot_product(x(0:i-1), a(i,0:i-1) ) ) / d(i) ! solve L y = b for y with back substitution
    enddo
    do i = n-1, 0, -1 ! solve L^* x = y for x with forward substitution
      do j = n-1, i+1, -1
        x(i) = x(i) - a(j, i) * x(j)
      enddo
      x(i) = x(i) / d(i)
    enddo

    ! finally handle negation if needed
    if(neg) x = -x

  end function cholesky_solve

  !--------------------------------------------------------------------------
  !
  ! FUNCTION: SphereCorrelator_Correlate
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief   compute the cross correlation between two spherical functions
  !
  ! @param flm: spherical harmonic coefficients for the first function
  ! @param gln: spherical harmonic coefficients for the second function
  ! @param fMr: true/false if there is/isn't a mirror plane in the first function
  ! @param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
  ! @param eu : orientation to refine
  ! @param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
  ! @return   : maximum cross correlation
  !
  !> @date 07/22/19 WCL 1.0 original, based on Will Lenthe's classes in sht_xcorr.hpp
  !--------------------------------------------------------------------------
  recursive function SphereCorrelator_RefinePeak(this, flm, gln, fMr, fNf, eu, eps) result(peak)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_RefinePeak

    use mod_global 

  implicit none

    class  (SphereCorrelator),INTENT(INOUT) :: this
    complex(kind=dbl        ),INTENT(IN   ) :: flm(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    complex(kind=dbl        ),INTENT(IN   ) :: gln(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    logical                  ,INTENT(IN   ) :: fMr 
    integer(kind=irg        ),INTENT(IN   ) :: fNf
    real   (kind=dbl        ),INTENT(INOUT) :: eu(0:2)
    real   (kind=dbl        ),INTENT(IN   ) :: eps
    real   (kind=dbl        )               :: peak

    integer                                 :: iter, maxIter, info, n, nrhs, i
    real   (kind=dbl        )               :: eu0(0:2), absEps, euEps, jac(0:2), hes(0:8), step(0:2), prevMag2, mag2, det, &
                                               hespacked(6), eqB(3,1), dk, eval(0:2), d, p, p1, p2, r, q
    character(1)                            :: uplo
    logical                                 :: converged

    ! call FatalError('SphereCorrelator_RefinePeak', 'not yet implemented') ! I think the only thing missing is calling the cholesky decomposition

    ! initial setup
    eu0      = eu ! save a copy of the original orientation
    absEps   = eps * cPi * 2.D0 / this%xcLut%slP ! compute stopping criterion is factor of grid resolution

    euEps    = sqrt(EPSILON(1.D0)) ! epsilon for eu[1] being too close to 0 or pi
    maxIter  = 45 ! in tests of perfect rotations convergence generally takes at most 3 iterations (some edge cases near degeneracies are slower)
    prevMag2 = cPi * 6 / this%xcLut%slP ! first step better not be more than 1 pixel in each direction
    n = 3
    nrhs = 1
    uplo = 'U'

    ! newton's method loop
    do iter = 1,maxIter
      ! build hessian and derivatives
      peak = this%derivatives(flm, gln, eu, jac, hes, this%xcLut%bw, fMr, fNf)

      ! pack the Hessian matrix in upper diagonal form
      hespacked = -(/ hes(0), hes(1), hes(4), hes(2), hes(5), hes(8) /)
      eqB(1:3,1) = jac(0:2)
      call DPPSV(uplo, n, nrhs, hespacked, eqB, n, info) ! solve for step (this will fail for non definate matricies but thats good since we don't want a saddle point)
      ! eqB contains the solution if info = 0

      if(info.gt.0) then ! DPPSV failed (we're close to a singularity)
        if(hes(4).ne.hes(4)) then ! the beta derivatives are NaN
          ! if we're actually on the degeneracy - the beta^2 derivative is undefined, do the 1x1 sub problem
          step(0) = jac(0) / hes(0);
          step(1) = 0;
          step(2) = 0;
        else ! not degenerate, either close or we have some other type of numerical instability
          ! check how close to the degeneracy we are
          dk = abs( dmod(eu(1) / cPi, 1.D0) ) * this%xcLut%slP / 2 ! distance in pixels
          if(dK.lt.2.D0) then ! we're close to the degeneracy (this is a somewhat arbitrary cutoff)
            ! not degenerate close enough that DPPSV failed
            ! if we're just close to the degeneracy sovle 2x2 sub problem manually (3rd euler angle is false DoF)
            det = hes(0) * hes(4) - hes(1) * hes(1);
            ! if(det.lt.euEps) then ! 2x2 problem is also singular
            !     if(abs(det) .lt. euEps) call FatalError('SphereCorrelator_RefinePeak', &
            !       'singular matrix during reduced newton iteration') ! don't divide by 0 this should be extremely rare
            !     if(det .lt. euEps) call FatalError('SphereCorrelator_RefinePeak', &
            !       'converging to saddle during reduced newton iteration') ! it isn't clear how often this can happen (limiting the step size prevents occurrence in all testing)
            ! endif

            ! compute 2x2 step (dont update 3rd euler angle)
            step(0) = (jac(0) * hes(4) - jac(1) * hes(1)) / det;
            step(1) = (jac(1) * hes(0) - jac(0) * hes(1)) / det;
            step(2) = 0;

          else ! we're far from the degeneracy and converging towards a saddle point
            ! hopefully this just means we didn't start close enough to the maximum
            ! if we made it to here we're in relatively unexplored territory.
            ! I've only found a few edge cases where this occurs
            ! use ~saddle free newton: https://arxiv.org/abs/1406.2572
            ! compute the eigen docomposition of the hessian H such that H = Q * D * Q^T
            ! we know that D is mixed sign since cholesky decomposition failed
            ! let d be the most negative eigenvalue (min(D) since we are interested in maxima)
            ! compute step with H' = Q * I*d * Q^T ==> step = J / d ==> we only need eigenvalues
            ! there is a nice closed form solution for 3x3 matricies: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3×3_matrices
            p1 = hes(1) * hes(1) + hes(2) * hes(2) + hes(5) * hes(5)
            if(0.D0 == p1) then
              ! already diagonal
              eval(0) = hes(0)
              eval(1) = hes(4)
              eval(2) = hes(8)
            else
              ! compute prefactors
              q = (hes(0) + hes(4) + hes(8)) / 3 ! trace / 3
              p2 = p1 * 2
              p1 = hes(0) - q ! a11 - q
              p2 = p2 + p1 * p1;
              p1 = hes(4) - q ! a22 - q
              p2 = p2 + p1 * p1;
              p1 = hes(8) - q ! a33 - q
              p2 = p2 + p1 * p1; ! sumsq( diag(A) - q ) + 2 * p1
              p = sqrt(p2 / 6.D0)

              ! overwrite hes with hes - q * I
              hes(0) = hes(0) - q
              hes(4) = hes(4) - q
              hes(8) = hes(8) - q

              ! now compute det(hes)
              r = hes(0) * hes(4) * hes(8) - hes(0) * hes(5) * hes(7) &
                + hes(1) * hes(5) * hes(6) - hes(1) * hes(3) * hes(8) &
                + hes(2) * hes(3) * hes(7) - hes(2) * hes(4) * hes(6)
              r = r / (p * p * p * 2)
              if(r.le.-1.D0) then
                r = cPi / 3
              else if(r.ge.1.D0) then
                r = 0.D0
              else
                r = acos(r) / 3
              endif

              ! finally compute eigen values
              eval(0) = q + p * 2 * cos(r              ) ! this is the maximum eigen value
              eval(2) = q + p * 2 * cos(r + cPi * 2 / 3) ! this is the minimum eigen value
              eval(1) = q * 3 - eval(0) - eval(2) ! this is the intermediate value
            endif

            ! select scaling factor (most negative is overly conservative and leads to slow convergence)
            if(eval(1).lt.0) then
              d = ( eval(1) + eval(2) ) / 2 ! 2/3 negative eigen values, average
              ! d = eval(1)
            elseif(eval(2).lt.0) then
              d = eval(2) ! only 1 negative value
            else
              ! no negative eigen values
              eu = eu0
              peak = 0
              return
              ! call FatalError('SphereCorrelator_RefinePeak', 'too far from local maxima for refinement') ! this isn't impossible but should be unlikely since we're starting at a voxel maxima
            endif
            step = jac / d
          endif
        endif
      else ! info.gt.0
        ! DPPSV worked (normal update)
        step(0:2) = -eqB(1:3,1)
        mag2 = step(0)*step(0) + step(1)*step(1) + step(2)*step(2) ! compute step size
        ! if(mag2 > prevMag2) call FatalError('SphereCorrelator_RefinePeak', 'newton steps must always decrease in magnitude') ! so far only a observed near degeneracies
        ! prevMag2 = mag2 ! update step size
      endif

      ! now that we've computed the step, apply and check for convergence
      eu = eu - step ! apply step
      if( maxval( abs(step) ) .lt. absEps ) then ! check for convergence
        return ! we're done
      endif
    enddo

    ! we failed to converge if we made it this far
    eu = eu0
    peak = 0
    return
    ! call FatalError('SphereCorrelator_RefinePeak', 'failed to converge during cross correlation refinement')

  end function SphereCorrelator_RefinePeak

  !--------------------------------------------------------------------------
  !
  ! FUNCTION: SphereCorrelator_Compute
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief   compute the cross correlation between two spherical functions
  !
  ! @param flm: spherical harmonic coefficients for the first function
  ! @param gln: spherical harmonic coefficients for the second function
  ! @param fMr: true/false if there is/isn't a mirror plane in the first function
  ! @param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
  ! @param eu : location to write rotation of maximum cross correlation as ZYZ euler angle
  ! @param ref: true/false to use/not use real space refinement
  ! @param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
  ! @return   : maximum cross correlation
  !   
  ! from the original C++ code :
  ! //naive implementation (no symmetry) to make summation clear
  ! /*
  ! std::fill(fxc.begin(), fxc.end(), std::complex<Real>(0));
  ! const bool realFft = true;//true/false to use half sized fft format
  ! const size_t dm = realFft ? bw : sl;//length of fastest indexing dimension
  ! for(size_t ic = 0; ic < sl; ic++) {
  !     const int k = ic >= bw ? int(ic) - sl : ic;
  !     const size_t ak = std::abs(k);
  !     for(size_t ib = 0; ib < sl; ib++) {
  !         const int n = ib >= bw ? int(ib) - sl : ib;
  !         const size_t an = std::abs(n);
  !         const size_t maxKN = std::max(ak, an);
  !         for(size_t ia = 0; ia < dm; ia++) {
  !             const int m = ia >= bw ? int(ia) - sl : ia;
  !             const size_t am = std::abs(m);
  !             const size_t start = std::max<size_t>(am, maxKN);
  !             for(size_t j = start; j < bw; j++) {
  !                 const Real dlkm = wigD[ak * bw * bw + am * bw + j] * wigner::dSign(j, k, m);//wigner::d<Real>(j, k, m);
  !                 const Real dlnk = wigD[an * bw * bw + ak * bw + j] * wigner::dSign(j, n, k);//wigner::d<Real>(j, n, k);
  !                 const std::complex<Real>& vflm = flm[am * bw + j];//\hat{f}^l_{|m|}
  !                 const std::complex<Real>& vgln = gln[an * bw + j];//\hat{g}^l_{|n|}
  !                 const std::complex<Real> f = std::signbit(m) ? std::conj(vflm) * Real(0 == am % 2 ? 1 : -1) : vflm;//symmetry of real SHT coefficients
  !                 const std::complex<Real> g = std::signbit(n) ? std::conj(vgln) * Real(0 == an % 2 ? 1 : -1) : vgln;//symmetry of real SHT coefficients
  !                 fxc[ic * sl * dm + ib * dm + ia] += f * std::conj(g) * dlkm * dlnk;
  !             }
  !         }
  !     }
  ! }
  ! */
  !
  !> @date 01/30/19 MDG 1.0 original, based on Will Lenthe's classes in sht_xcorr.hpp
  !--------------------------------------------------------------------------
  recursive subroutine SphereCorrelator_Compute(this, flm, gln, fMr, fNf)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_Compute

    use mod_global
    use mod_Wigner

  implicit none

    class  (SphereCorrelator) ,INTENT(INOUT) :: this
    complex(kind=dbl         ),INTENT(IN   ) :: flm(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    complex(kind=dbl         ),INTENT(IN   ) :: gln(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    logical                   ,INTENT(IN   ) :: fMr 
    integer(kind=irg         ),INTENT(IN   ) :: fNf

    integer                                  :: mBw , flmFold, glnFold, slP, bwP           ! storage for useful values
    logical                                  :: fMir, gMir   , mirror                      ! storage for useful values
    complex(kind=dbl)                        :: cz = cmplx(0.D0, 0.D0)                     ! complex zero
    logical                                  :: nonZero(0:1,0:this%xcLut%bwP-1)            ! locations of non zero columns
    logical                                  :: bMirror, match0, match1, mir0, mir1, mFold ! intermediate variables for filling nonZero array
    logical                                  :: posN, posK, posKN                          ! check if there are conjugate values once per column
    logical                                  :: match, toggle                              ! inner most loop (m) checks
    integer                                  :: j, k, n, m, start, maxKN                   ! loop counters
    integer                                  :: k1, n1, m1, sl1                            ! convertion from 0 -> 1 indexing
    complex(kind=dbl)                        :: v, vnc, vp, vc                             ! working space

    ! the above quadruple loop is conceptually simple but has many redundant calculations
    ! the loop that follows is mathematically equivalent (for SHT of real functions only!!) but much faster:
    !  -use \hat{f}^l_{-m} = (-1)^m * \hat{f}^l_{m} for real valued functions (and the same for \hat{g}^l_{-n})
    !  -build the fft of a real valued cross correlation
    !  -precompute values of \hat{f}^l_{m} * d^l_{k,m}(\frac{\pi}{2}) (stored in fm)
    !  -precompute values of \hat{g}^l_{n} * d^l_{n,k}(\frac{\pi}{2}) (stored in gn)
    !  -eliminate redundant calculations from f * g and f * conj(g)

    ! save some useful values
    mBw     = this%xcLut%bw  ! maximum bandwidth to use in calculation (better be <= bw but could be smaller for slight speedup, effectively a top hat filter size)
    bwP     = this%xcLut%bwP ! zero padded bandwidth
    slP     = this%xcLut%slP ! zero padded side length
    sl1     = slP + 1        ! for 1 indexing
    flmFold = fNf            ! rotational symmetry of flm about z axis (flm[m*bw+j] == 0 if m % flmFold != 0)
    glnFold = 1              ! rotational symmetry of gln about z axis (gln[n*bw+j] == 0 if n % glnFold != 0)
    fMir    = fMr            ! true/false if (flm[m*bw+j] == 0 if (m+j) % 2 != 0)
    gMir    = .FALSE.        ! true/false if (gln[n*bw+j] == 0 if (n+j) % 2 != 0)
    mirror  = fMir.or.gMir   ! is there at least 1 mirror

    ! one value for each column for n % 2 == 0 and n % 2 == 1
    ! true if cross correlation at n, m is nonzero, false if it is a systemic zero
    ! using an integer array is much faster than bools becuase of the vector<bool> specialization
    nonZero(:,:) = .FALSE. ! initialize with false for padding columns
    bMirror = fMir.and.gMir
    do m = 0, mBW-1
        match0 = mod(m + 0, 2).eq.0         ! check if parity of m matches parity of even n
        match1 = mod(m + 1, 2).eq.0         ! check if parity of m matches parity of odd  n
        mir0   = bMirror.and.(.not.match0)  ! there are systemic zeros if both functions have a mirror plane and a parity mismatch
        mir1   = bMirror.and.(.not.match1)  ! there are systemic zeros if both functions have a mirror plane and a parity mismatch
        mFold  = mod(m, flmFold).ne.0       ! there are systemic zeros from rotational symmetry in flm
        nonZero(0,m) = .not.(mFold.or.mir0) ! systemic zeros in m for n%2 == 0
        nonZero(1,m) = .not.(mFold.or.mir1) ! systemic zeros in m for n%2 == 1
    enddo ! m

    ! loop over planes
    do k = 0, mBw-1
        k1 = k + 1 ! c_f_pointer doesn't currently support arbitrary lower bounds (fxc is 1 indexed and everything else is 0 indexed)

        ! precompute flm values * wigner d function
        do m = 0, mBw-1
            do j = max(m, k), mBw-1
                this%fm(j, m) = flm(j, m) * this%xcLut%wigD(j, m, k) ! f^j_m * wigner::d<Real>(j, k, m)
            enddo ! j
        enddo ! m (cols)

        ! loop over rows
        posK = k.gt.0 ! is there a negative k value
        do n = 0, bwP-1
            n1 = n + 1 ! c_f_pointer bounds adjustment
            if( ( mod(n, glnFold).eq.0 ).and.( n.lt.mBw ) ) then ! we haven't reached 0 padded rows and gln values are non-zero
                posN  = n.gt.0        ! is there a negative n value
                posKN = posK.and.posN ! are there negative n and k values
                maxKN = max(k, n)     ! compute maximum of k and n once

                ! precompute gln values * wigner d function
                do j = maxKN, mBw-1
                    this%gn(j) = conjg( gln(j, n) ) * this%xcLut%wigD(j, k, n) ! \hat{g}^j_n * wigner::d<Real>(j, n, k)
                enddo ! j

                ! loop over columns
                do m = 0, bwP-1
                    m1 = m + 1 ! c_f_pointer bounds adjustment
                    if( nonZero(mod(m, 2), m) ) then ! checking for systemic zeros from double mirror and parity mismatch here makes subsequent logic easy
                        v   = cz ! initialize dot product
                        vnc = cz ! initialize conjugate dot product
                        start = max(maxKN, m) ! max(k, n, m)
                        if(mirror) then ! there is a single mirror or double mirror with parity matching
                            if(fMir.and.( mod(start+m, 2).ne.0 )) start = start + 1 ! if fm[start] == 0 skip to next value (first value is zero)
                            if(gMir.and.( mod(start+n, 2).ne.0 )) start = start + 1 ! if gn[start] == 0 skip to next value (first value is zero) [for double mirrors no change here since parities match]
                            toggle = mod(start + m, 2).eq.0 ! we don't need to toggle since we're incrementing by 2 but we still may need to negate the result
                            do j = start, mBw-1, 2
                                ! do complex multiplication components by hand to eliminate duplicate flops from multiplying with conjugate
                                call SH_conjMult(this%fm(j, m), this%gn(j), vp, vc)
                                v   = v   + vp ! fm(j,m) *           gn[j]
                                vnc = vnc + vc ! fm(j,m) * std::conj(gn[j])
                            enddo ! j
                            if(toggle) vnc = -vnc ! negate sum instead of elements
                        else ! mirror
                            toggle = mod(start+m, 2).eq.0
                            do j = start, mBw-1 ! without a mierror we need to increment by 1 instead of 2
                                call SH_conjMult(this%fm(j, m), this%gn(j), vp, vc)
                                v   = v   + vp ! fm(j,m) *           gn[j]
                                if(toggle) then
                                    vnc = vnc + vc ! fm(j,m) * std::conj(gn[j])
                                else ! toggle
                                    vnc = vnc - vc ! fm(j,m) * std::conj(gn[j])
                                endif ! toggle
                                toggle = .not.toggle
                            enddo ! j
                        endif ! mirror
                        if(mod(k, 2).ne.0) vnc = -vnc ! correct for computing negative vnc depending on j/m parity

                        ! fill in symmetric values using symmetry from: wigner d function, sht of real signal, sht of real pattern
                        match = mod(m+n, 2).eq.0
                        if(posKN) then
                            this%fxc(m1,    n1,    k1) = v
                            this%fxc(m1,sl1-n ,sl1-k ) = vnc
                            if(match) then
                                this%fxc(m1,    n1,sl1-k ) =  v
                                this%fxc(m1,sl1-n ,    k1) =  vnc
                            else ! match
                                this%fxc(m1,    n1,sl1-k ) = -v
                                this%fxc(m1,sl1-n ,    k1) = -vnc
                            endif ! match
                        else ! posKN
                            this%fxc(m1,    n1,    k1) = v
                            if(match) then
                                if(posK) this%fxc(m1,    n1,sl1-k ) =  v
                                if(posN) this%fxc(m1,sl1-n ,    k1) =  vnc
                            else ! match
                                if(posK) this%fxc(m1,    n1,sl1-k ) = -v
                                if(posN) this%fxc(m1,sl1-n ,    k1) = -vnc
                            endif ! match
                        endif ! posKN
                    else ! nonzero
                        ! in systemic zero column
                        if(posKN) then
                            this%fxc(m1,    n1,    k1) = cz
                            this%fxc(m1,sl1-n ,sl1-k ) = cz
                            this%fxc(m1,    n1,sl1-k ) = cz
                            this%fxc(m1,sl1-n ,    k1) = cz
                        else ! posKN
                            this%fxc(m1,    n1,    k1) = cz
                            if(posK) this%fxc(m1,    n1,sl1-k ) = cz
                            if(posN) this%fxc(m1,sl1-n ,    k1) = cz
                        endif ! posKN
                    endif ! nonzero
                enddo ! m (cols)
            else ! nonzero gln values
                ! in systemic zero row
                this%fxc(:,    n1,    k1) = cz
                if(posKN) this%fxc(:,sl1-n ,sl1-k ) = cz
                if(posK ) this%fxc(:,    n1,sl1-k ) = cz
                if(posN ) this%fxc(:,sl1-n ,    k1) = cz
            endif ! nonzero gln values
        enddo ! n (rows)
    enddo ! k (slices)

    ! fill in zero pad slices if needed
    if(mBw < bwP) then
        this%fxc(:,:,  1+mBw:  1+bwP) = cz
        this%fxc(:,:,sl1-bwP:sl1-mBw) = cz
    endif ! mBw < bwP

    ! do inverse fft
    call this%xcLut%plan%inverse(this%fxc, this%xc)

  end subroutine SphereCorrelator_Compute

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SphereCorrelator_FindPeak
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief find the maximum cross correlation (in this%xc)
  !
  !> @return index of maximum
  !
  !> @date 07/22/19 WCL 1.0 original
  !--------------------------------------------------------------------------
  function SphereCorrelator_FindPeak(this) result(peak)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_FindPeak
  implicit none
    class  (SphereCorrelator),INTENT(IN) :: this
    integer                              :: peak(3)
    peak = maxloc(this%xc(:,:,1:this%xcLut%bwP)) ! get index of peak
  end function SphereCorrelator_FindPeak

  !--------------------------------------------------------------------------
  !
  ! FUNCTION: SphereCorrelator_Derivatives
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief   compute the first and second derivatives of the cross correlation at a single rotation
  !
  ! @param flm: spherical harmonic coefficients for the first function
  ! @param gln: spherical harmonic coefficients for the second function
  ! @param eu : rotation to compute derivatives of cross correlation for as ZYZ euler angle
  ! @param jac: location to write jacobian of cross correlation {d/(d eu[0]), d/(d eu[1]), d/(d eu[2])}
  ! @param hes: location to write hessian (3x3 matrix as 9 component vector) of cross correlation hes_ij = d/(d eu[i]) * d/(d eu[j])
  ! @param mBW: maximum bandwidth to use in calculation (must be <= bw)
  ! @param fMr: true/false if there is/isn't a mirror plane in the first fuction
  ! @param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
  ! @return   : maximum cross correlation
  !
  !> @date 03/18/19 MDG 1.0 original, based on Will Lenthe's classes in sht_xcorr.hpp
  !--------------------------------------------------------------------------
  recursive function SphereCorrelator_Derivatives(this, flm, gln, eu, jac, hes, mBW, fMr, fNf) result(corr)
  !DEC$ ATTRIBUTES DLLEXPORT :: SphereCorrelator_Derivatives

    use mod_global
    use mod_Wigner

  implicit none

    class  (SphereCorrelator),INTENT(INOUT) :: this
    complex(kind=dbl        ),INTENT(IN   ) :: flm(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    complex(kind=dbl        ),INTENT(IN   ) :: gln(0:this%xcLut%bw-1,0:this%xcLut%bw-1)
    real   (kind=dbl        ),INTENT(INOUT) :: eu(0:2)
    real   (kind=dbl        ),INTENT(INOUT) :: jac(0:2)
    real   (kind=dbl        ),INTENT(INOUT) :: hes(0:8)
    integer(kind=irg        ),INTENT(IN   ) :: mBW 
    logical                  ,INTENT(IN   ) :: fMr 
    integer(kind=irg        ),INTENT(IN   ) :: fNf
    real   (kind=dbl        )               :: corr 

    integer(kind=irg        )               :: flmFold, glnFold, dJ, m, start, j, n, k
    real   (kind=dbl        )               :: wrk(0:9), beta, tpi, sA, cA, sG, cG, t, csc, uA(0:2), tA(0:2), mm, mn, nn, sm, cm, &
                                               uG(0:2), tG(0:2), sn, cn, sign, sgn, coef2_0a, coef2_0b, coef2_1a, coef1_0PP, &
                                               coef1_0PN, coef2_0PP, coef2_0PN, coef2_1PP, coef2_1PN, d0P, d0N, &
                                               d0P_1, d0N_1, d0P_2, d0N_2, rjm, coef2_2, d1P, d1N, d2P, d2N, xc(0:9), xp(0:9)
    complex(kind=dbl        )               :: expAlpha, expGamma, agP, agN, vp, vc, vcPP, vcPP0, vcPP1, vpPN, vpPN0, vpPN1
    logical                                 :: nB, fMir, gMir, mirror, bmirror, mFold0, nFold0, match, mir0 

    tpi = 2.D0 * cPi 

    ! initialize terms with 0
    wrk = 0.D0 ! correlation, jacobian, hessian as 00, 11, 22, 01, 12, 20

    ! bring middle Euler angle to [-pi,pi] for Wigner calculations
    beta = mod(eu(1), tpi)
    if (beta.gt.cPi) then
        beta = beta - tpi
    else 
        if (beta.lt.-cPi) then
            beta = beta + tpi 
        end if 
    end if

    ! compute sin/cos of alpha/gamma once (multiple angles)
    sA = sin(eu(0)) ! sin(alpha)
    cA = cos(eu(0)) ! cos(alpha)
    sG = sin(eu(2)) ! sin(gamma)
    cG = cos(eu(2)) ! cos(gamma)

    ! precompute some values for on the fly Wigner (uppercase) D calculation
    t = cos(beta)
    nB = .TRUE.
    if (beta.gt.0.D0) nB = .FALSE.   ! std::signbit(beta)
    csc = 1.D0 / sqrt(1.D0 - t * t) ! csc(beta), cot(beta) is csc * t
    if (nb.eqv..TRUE.) csc = -csc
    call Wigner_dTable2(mBW, t, nB, this%dBeta) ! ::wigner::dTable(mBW, t, nB, dBeta.data());//compute wigner (lowercase) d(beta) once

    ! build symmetry information
    flmFold = fNf       ! rotational symmetry of flm about z axis (flm[m*bw+j] == 0 if m % flmFold != 0)
    glnFold = 1         ! rotational symmetry of gln about z axis (gln[n*bw+j] == 0 if n % glnFold != 0)
    fMir = fMr          ! true/false if (flm[m*bw+j] == 0 if (m+j) % 2 != 0)
    gMir = .FALSE.      ! true/false if (gln[n*bw+j] == 0 if (n+j) % 2 != 0)
    mirror = (fMir.or.gMir)     ! is there at least 1 mirror
    bMirror =(fMir.and.gMir)    ! do both functions have a mirror
    dJ = 2
    if (fMir.eqv..FALSE.) dJ = 1    ! dJ = fMir ? 2 : 1;

    !///////////////////////////////////////
    !/        loop over one order         //
    !///////////////////////////////////////
    uA = (/ 0.D0, sA * 2.D0, 1.D0 /)     ! recursion coefficients for chebyshev polynomial U_n(sin(alpha))
    tA = (/ 0.D0, cA       , 1.D0 /)     ! recursion coefficients for chebyshev polynomial T_n(cos(alpha))
    do  m = 0, mBW-1 
        !///////////////////////////////////////
        !/ efficiently compute exp(I m alpha) //
        !///////////////////////////////////////
        ! update chebyshev recursion and use to compute exp(I m alpha)
        if (m.lt.2) then ! use seed values for chebyshev recursion
            if (m.eq.0) then 
                uA(0) = 0.D0
                tA(0) = 1.D0
            else
                uA(0) = sA      ! sin(alpha * m)
                tA(0) = cA      ! cos(alpha * m)
            end if
        else  ! use chebyshev recursion
            ! compute chebyshev polynomials(m, alpha) and multiple angle sin/cos
            uA(0) = sA * uA(1) * 2.D0 - uA(2)      ! U_m(x) = 2 * x * U_{m-1}(x) - U_{m-2}(x)
            tA(0) = cA * tA(1) * 2.D0 - tA(2)      ! T_m(x) = 2 * x * T_{m-1}(x) - T_{m-2}(x)

            mm = -1.D0
            if (mod((m/2-1),2).eq.0) mm = 1.D0
            nn = -1.D0 
            if (mod((m-1)/2,2).eq.0) nn = 1.D0

            ! sm = m % 2 == 0 ? uA[1] * cA * (((m/2)-1) % 2 == 0 ? 1 : -1) : (uA[0] - sA * uA[1]) * (((m-1)/2) % 2 == 0 ? 1 : -1);//cos(alpha * m)
            if (mod(m,2).eq.0) then ! cos(alpha * m)            
                sm = uA(1) * cA * mm
            else 
                sm = (uA(0) - sA * uA(1)) * nn
            end if
            cm = tA(0)      ! cos(alpha * m)

            ! update recursion and store values of sin/cos(alpha * m)
            uA(2) = uA(1)       ! 
            uA(1) = uA(0)       ! update recursion coefficients for chebyshev polynomial of the first kind
            tA(2) = tA(1)       ! 
            tA(1) = tA(0)       ! update recursion coefficients for chebyshev polynomial of the second kind
            uA(0) = sm          ! store multiple sin value for subsequent access
            tA(0) = cm          ! store multiple cos value for subsequent access
        end if ! m.lt.2

        mFold0 = 0.ne.mod(m, flmFold)  ! there are systemic zeros from rotational symmetry in flm
        if (mFold0.eqv..TRUE.) CYCLE  ! continue;//flm[m * bw + j] == 0 so there is nothing to accumulate (but we still needed to update the multiple angle recursion)
        expAlpha = cmplx(tA(0), uA(0)) ! exp(I m alpha) = cos(m * alpha) + I sin(m * alpha)

        !///////////////////////////////////////
        !/       loop over other order        //
        !///////////////////////////////////////
        uG(0:2) = (/ 0.D0, sG * 2.D0, 1.D0 /)  ! recursion coefficients for chebyshev polynomial U_n(sin(gamma))
        tG(0:2) = (/ 0.D0, cG       , 1.D0 /)  ! recursion coefficients for chebyshev polynomial T_n(cos(gamma))

        do n = 0, mBW-1 
            !///////////////////////////////////////
            !/ efficiently compute exp(I n gamma) //
            !///////////////////////////////////////
            ! update chebyshev recursion and use to compute exp(I n gamma)          
            if (n.lt.2) then ! use seed values for chebyshev recursion
                if (n.eq.0) then 
                    uG(0) = 0.D0
                    tG(0) = 1.D0
                else
                    uG(0) = sG      ! sin(gamma * n)
                    tG(0) = cG      ! cos(gamma * n)
                end if
            else  ! use chebyshev recursion
                ! compute chebyshev polynomials(n, gamma) and multiple angle sin/cos
                uG(0) = sG * uG(1) * 2.D0 - uG(2)      ! U_m(x) = 2 * x * U_{m-1}(x) - U_{m-2}(x)
                tG(0) = cG * tG(1) * 2.D0 - tG(2)      ! T_m(x) = 2 * x * T_{m-1}(x) - T_{m-2}(x)

                mm = -1.D0
                if (mod((n/2-1),2).eq.0) mm = 1.D0
                nn = -1.D0 
                if (mod((n-1)/2,2).eq.0) nn = 1.D0

                ! sm = m % 2 == 0 ? uA[1] * cA * (((m/2)-1) % 2 == 0 ? 1 : -1) : (uA[0] - sA * uA[1]) * (((m-1)/2) % 2 == 0 ? 1 : -1);//cos(alpha * m)
                if (mod(n,2).eq.0) then ! cos(alpha * m)            
                    sn = uG(1) * cG * mm
                else 
                    sn = (uG(0) - sG * uG(1)) * nn
                end if
                cn = tG(0)      ! cos(alpha * m)

                ! update recursion and store values of sin/cos(alpha * m)
                uG(2) = uG(1)       ! 
                uG(1) = uG(0)       ! update recursion coefficients for chebyshev polynomial of the first kind
                tG(2) = tG(1)       ! 
                tG(1) = tG(0)       ! update recursion coefficients for chebyshev polynomial of the second kind
                uG(0) = sn          ! store multiple sin value for subsequent access
                tG(0) = cn          ! store multiple cos value for subsequent access
            end if ! n.lt.2

            nFold0 = 0.ne.mod(n, flmFold) ! there are systemic zeros from rotational symmetry in flm
            if (nFold0.eqv..TRUE.) CYCLE  ! flm[m * bw + j] == 0 so there is nothing to accumulate (but we still needed to update the multiple angle recursion)
            expGamma = cmplx(tG(0), uG(0)) ! exp(I n gamma) = cos(n * gamma) + I sin(n * gamma)

            ! handle the case of 2 mirror planes with different parity
            match = .FALSE. 
            if (mod(m + n, 2).eq.0) match = .TRUE. ! check if parity of m matches parity of n
            mir0 = (bMirror.and.(.not.match))      ! there are systemic zeros if both functions have a mirror plane and a parity mismatch
            if (mir0.eqv..TRUE.) CYCLE             ! flm * gln = 0 for any l (alternating between flm and gln being 0)

            !///////////////////////////////////////
            !/  compute degree independent terms  //
            !///////////////////////////////////////
            !/compute exp(I * +m * Alpha) * exp(I * +/-n * Gamma) for on the fly wigner (uppercase) D calculation
            call SH_conjMult(expAlpha, expGamma, agP, agN)
            sign = -1.D0
            if (mod(n+m,2).eq.0) sign = 1.D0
            sgn = -1.D0
            if (mod(n,2).eq.0) sgn = 1.D0
            agP = agP * sign
            agN = agN * sign * sgn

            ! compute some prefactors for calculating derivatives of d^j_{m,n}(beta)
            mm = m * m
            mn = m * n
            nn = n * n
            coef2_0a  =   t * t * mm + (nn - m)              
            coef2_0b  =   t * n * (1 - 2 * m)                
            coef2_1a  =   t *     (1 + 2 * m)                
            coef1_0PP = ( t * m      - n       ) * csc       
            coef1_0PN = ( t * m      + n       ) * csc       
            coef2_0PP = (coef2_0a    + coef2_0b) * csc * csc 
            coef2_0PN = (coef2_0a    - coef2_0b) * csc * csc 
            coef2_1PP = (coef2_1a    - 2 * n   ) * csc       
            coef2_1PN = (coef2_1a    + 2 * n   ) * csc       

            !///////////////////////////////////////
            !/   loop over degrees accumulating   //
            !///////////////////////////////////////
            start = max(m, n)
            if (fMir.and.(mod(start + m, 2).ne.0)) start = start + 1  ! if fm[start] == 0 skip to next value (first value is zero)
            if (gMir.and.(mod(start + n, 2).ne.0)) start = start + 1  ! if gn[start] == 0 skip to next value (first value is zero) [for double mirrors no change here since parities match]
            do j = start, mBW-1, dJ ! increment by 1 for no mirror planes 2 if any are present
                ! get wigner d^j_{m,+/-n} components
                d0P    =  this%dBeta(0,j,n,m)      ! d^j_{m  ,n}(     beta)
                d0N    =  this%dBeta(1,j,n,m)      ! d^j_{m  ,n}(pi - beta)
                if (m.ge.j) then 
                    d0P_1 = 0.D0 
                    d0N_1 = 0.D0
                else
                    d0P_1  = this%dBeta(0,j,n,m+1) ! d^j_{m+1,n}(     beta)
                    d0N_1  = this%dBeta(1,j,n,m+1) ! d^j_{m+1,n}(pi - beta)
                end if
                if ((m+1).ge.j) then 
                    d0P_2 = 0.D0 
                    d0N_2 = 0.D0
                else
                    d0P_2 = this%dBeta(0,j,n,m+2) ! d^j_{m+2,n}(     beta)
                    d0N_2 = this%dBeta(1,j,n,m+2) ! d^j_{m+2,n}(pi - beta)
                end if

                ! compute derivatives of d^j_{m,+/-n}(beta)
                rjm      = sqrt( dble( (j - m    ) * (j + m + 1) ) )
                if(j.eq.m) then
                  coef2_2 = 0.D0
                else
                  coef2_2  = sqrt( dble( (j - m - 1) * (j + m + 2) ) ) * rjm
                endif
                d1P = d0P * coef1_0PP - d0P_1 * rjm                                ! first  derivative of d^j_{+m,+n}(beta) w.r.t. beta
                d1N = d0N * coef1_0PN + d0N_1 * rjm                                ! first  derivative of d^j_{+m,-n}(beta) w.r.t. beta
                d2P = d0P * coef2_0PP - d0P_1 * rjm * coef2_1PP + d0P_2 * coef2_2  ! second derivative of d^j_{+m,+n}(beta) w.r.t. beta
                d2N = d0N * coef2_0PN + d0N_1 * rjm * coef2_1PN + d0N_2 * coef2_2  ! second derivative of d^j_{+m,-n}(beta) w.r.t. beta
                
                !compute f^l_m * g^l_n    \hat{f}^l_{+m} * hat{g}^l_{+n} and \hat{f}^l_{+m} * conj(hat{g}^l_{+n})
                call SH_conjMult(flm(j,m), gln(j,n), vp, vc)  ! do complex multiplication components by hand to eliminate duplicate flops from multiplying with conjugate
                if (mod(j+m,2).ne.0) vp = -vp

                ! compute components of cross correlation and partials
                vcPP  = vc   * agP  ! +n correlation term prefactor: \hat{f}^l_{+m} * conj(hat{g}^l_{+n}) * exp(I m alpha + I n gamma)
                vcPP0 = vcPP * d0P  ! +n correlation term: \hat{f}^l_{+m} * conj(hat{g}^l_{+n}) * D^l_{+m,+n}(alpha, beta, gamma)
                vcPP1 = vcPP * d1P  ! beta partial of +n correlation term
                vpPN  = vp   * agN  ! +n correlation term prefactor: \hat{f}^l_{+m} * conj(hat{g}^l_{-n}) * exp(I m alpha - I n gamma)
                vpPN0 = vpPN * d0N  ! -n correlation term: \hat{f}^l_{+m} * conj(hat{g}^l_{-n}) * D^l_{+m,-n}(alpha, beta, gamma)
                vpPN1 = vpPN * d1N  ! beta partial of -n correlation term

                ! compute contributions to cross correlation, jacobian, and hessian from +m,+n (the real part of contributions from -m,-n are the same for real functions)
                xc = (/ real(vcPP0), &                                                ! cross correlation
                        -aimag(vcPP0) * m ,  real (vcPP1)      , -aimag(vcPP0) * n , & ! alpha  , beta  , gamma   derivatives
                        -real (vcPP0) * mm,  real (vcPP ) * d2P, -real (vcPP0) * nn, & ! alpha^2, beta^2, gamma^2 derivatives
                        -aimag(vcPP1) * m , -aimag(vcPP1) * n  , -real (vcPP0) * mn /) ! alpha beta, beta gamma, gamma alpha derivatives

                ! compute contributions to cross correlation, jacobian, and hessian from +m,-n (the real part of contributions from -m,+n are the same for real functions)
                xp = (/ real(vpPN0), &                                                 ! cross correlation
                        -aimag(vpPN0) * m , real (vpPN1)       ,  aimag(vpPN0) * n , & ! alpha  , beta  , gamma   derivatives
                        -real (vpPN0) * mm, real (vpPN ) * d2N , -real (vpPN0) * nn, & ! alpha^2, beta^2, gamma^2 derivatives
                        -aimag(vpPN1) * m , aimag(vpPN1) * n   ,  real (vpPN0) * mn /) ! alpha beta, beta gamma, gamma alpha derivatives

                ! accumulate contributions
                wrk = wrk + xc                          ! std::transform(wrk, wrk + 10, xc, wrk, std::plus<Real>());//+m,+n
                if (n.gt.0) wrk = wrk + xp              ! std::transform(wrk, wrk + 10, xp, wrk, std::plus<Real>());//+m,-n
                if (m.gt.0) then 
                    wrk = wrk + xp                      ! std::transform(wrk, wrk + 10, xp, wrk, std::plus<Real>());//-m,+n
                    if (n.gt.0) wrk = wrk + xc          ! std::transform(wrk, wrk + 10, xc, wrk, std::plus<Real>());//-m,-n
                end if
            end do ! j
        end do ! n
    end do ! m

    ! copy result to outputs and return correlation
    jac = wrk(1:3)    ! d/dAlpha, d/dBeta ,d/dGamma

    ! hessian in row major order (for consistency w/ c++ code)
    hes(0) = wrk(4+0) ! d^2 / dAlpha^2
    hes(1) = wrk(4+3) ! d^2 / dAlpha dBeta
    hes(2) = wrk(4+5) ! d^2 / dAlpha dGamma
    hes(4) = wrk(4+1) ! d^2 / dBeta^2
    hes(5) = wrk(4+4) ! d^2 / dGamma dBeta
    hes(8) = wrk(4+2) ! d^2 / dGamma^2
    hes(3) = hes(1)   ! make symmetric
    hes(6) = hes(2)   ! make symmetric
    hes(7) = hes(5)   ! make symmetric

    ! save cross correlation
    corr = wrk(0)

  end function SphereCorrelator_Derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                          Helper Functions                          !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !--------------------------------------------------------------------------
  !
  ! FUNCTION: SH_interpolateMaxima
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief   interpolate subpixel peak location from a 3d voxel grid
  !
  ! @param p: neighborhood around peak
  ! @param x: location to store subpixel maxima location within neighborhood (x, y, z from -1->1)
  ! @param return: value of fit quadratic at maxima
  !
  !> @date 01/30/19 MDG 1.0 original, based on Will Lenthe's classes in sht_xcorr.hpp
  !--------------------------------------------------------------------------
  recursive function SH_interpolateMaxima(p, x, maxIter) result(vPeak)
  !DEC$ ATTRIBUTES DLLEXPORT :: SH_interpolateMaxima

  implicit none

    real   (kind=dbl),INTENT(IN   ) :: p(-1:1,-1:1,-1:1)
    real   (kind=dbl),INTENT(INOUT) :: x(0:2)
    integer(kind=irg),INTENT(IN   ) :: maxIter
    real   (kind=dbl)               :: vPeak

    integer(kind=irg)               :: i 
    real   (kind=dbl)               :: a000, a001, a002, a010, a020, a100, a200, a022, a011, a012, a021, a220, a110, a120, a122, &
                                       a210, a202, a101, a201, a102, a222, a211, a121, a112, a111, a212, a221, eps, xx, xy, &
                                       h00, h11, h22, h01, h12, h02, det, i00, i11, i22, i01, i12, i02, d0, d1, d2, step(0:2), &
                                       maxStep, yy, zz, zx, yz

    ! compute the 27 biquadratic coefficients, f(x,y,z) = a_{kji} x^i y^j z^k
    ! f(0, 0, 0) == a000
    a000 = p(0,0,0)

    ! f(1,0,0) = a000 + a100 + a200 && f(-1,0,0) = a000 - a100 + a200
    a001 = (p( 0, 0, 1) - p( 0, 0,-1)) / 2.D0
    a002 = (p( 0, 0, 1) + p( 0, 0,-1)) / 2.D0 - a000

    ! same relationships for y and z
    a010 = (p( 0, 1, 0) - p( 0,-1, 0)) / 2.D0
    a020 = (p( 0, 1, 0) + p( 0,-1, 0)) / 2.D0 - a000
    a100 = (p( 1, 0, 0) - p(-1, 0, 0)) / 2.D0
    a200 = (p( 1, 0, 0) + p(-1, 0, 0)) / 2.D0 - a000

    ! f( 1, 1,0) = a000 + a100 + a200 + a010 + a020 + a110 + a210 + a120 + a220
    ! f( 1,-1,0) = a000 + a100 + a200 - a010 + a020 - a110 - a210 + a120 + a220
    ! f(-1, 1,0) = a000 - a100 + a200 + a010 + a020 - a110 + a210 - a120 + a220
    ! f(-1,-1,0) = a000 - a100 + a200 - a010 + a020 + a110 - a210 - a120 + a220
    !  --> f( 1, 1,0) + f( 1,-1,0) + f(-1, 1,0) + f(-1,-1,0) = 4 * (a000 + a020 + a200 + a220)
    !  --> f( 1, 1,0) - f( 1,-1,0) - f(-1, 1,0) + f(-1,-1,0) = 4 * a110
    !  --> f( 1, 1,0) - f( 1,-1,0) + f(-1, 1,0) - f(-1,-1,0) = 4 * (a100 + a120)
    !  --> f( 1, 1,0) + f( 1,-1,0) - f(-1, 1,0) - f(-1,-1,0) = 4 * (a010 + a210)
    a022 = (p( 0, 1, 1) + p( 0, 1,-1) + p( 0,-1, 1) + p( 0,-1,-1)) / 4.D0 - a000 - a020 - a002
    a011 = (p( 0, 1, 1) - p( 0, 1,-1) - p( 0,-1, 1) + p( 0,-1,-1)) / 4.D0
    a012 = (p( 0, 1, 1) + p( 0, 1,-1) - p( 0,-1, 1) - p( 0,-1,-1)) / 4.D0 - a010
    a021 = (p( 0, 1, 1) - p( 0, 1,-1) + p( 0,-1, 1) - p( 0,-1,-1)) / 4.D0 - a001

    ! same relationships for yz and zx
    a220 = (p( 1, 1, 0) + p( 1,-1, 0) + p(-1, 1, 0) + p(-1,-1, 0)) / 4.D0 - a000 - a200 - a020
    a110 = (p( 1, 1, 0) - p( 1,-1, 0) - p(-1, 1, 0) + p(-1,-1, 0)) / 4.D0
    a120 = (p( 1, 1, 0) + p( 1,-1, 0) - p(-1, 1, 0) - p(-1,-1, 0)) / 4.D0 - a100
    a210 = (p( 1, 1, 0) - p( 1,-1, 0) + p(-1, 1, 0) - p(-1,-1, 0)) / 4.D0 - a010
    a202 = (p( 1, 0, 1) + p(-1, 0, 1) + p( 1, 0,-1) + p(-1, 0,-1)) / 4.D0 - a000 - a002 - a200
    a101 = (p( 1, 0, 1) - p(-1, 0, 1) - p( 1, 0,-1) + p(-1, 0,-1)) / 4.D0
    a201 = (p( 1, 0, 1) + p(-1, 0, 1) - p( 1, 0,-1) - p(-1, 0,-1)) / 4.D0 - a001
    a102 = (p( 1, 0, 1) - p(-1, 0, 1) + p( 1, 0,-1) - p(-1, 0,-1)) / 4.D0 - a100

    ! similar relationships for corners
    a222 = (p( 1, 1, 1) + p(-1,-1,-1) + p(-1, 1, 1) + p( 1,-1, 1) + p( 1, 1,-1) + p( 1,-1,-1) + p(-1, 1,-1) + p(-1,-1, 1)) &
           / 8.D0 - a000 - a200 - a020 - a002 - a022 - a202 - a220
    a211 = (p( 1, 1, 1) + p(-1,-1,-1) + p(-1, 1, 1) - p( 1,-1, 1) - p( 1, 1,-1) + p( 1,-1,-1) - p(-1, 1,-1) - p(-1,-1, 1)) &
           / 8.D0 - a011
    a121 = (p( 1, 1, 1) + p(-1,-1,-1) - p(-1, 1, 1) + p( 1,-1, 1) - p( 1, 1,-1) - p( 1,-1,-1) + p(-1, 1,-1) - p(-1,-1, 1)) &
           / 8.D0 - a101
    a112 = (p( 1, 1, 1) + p(-1,-1,-1) - p(-1, 1, 1) - p( 1,-1, 1) + p( 1, 1,-1) - p( 1,-1,-1) - p(-1, 1,-1) + p(-1,-1, 1)) &
           / 8.D0 - a110
    a111 = (p( 1, 1, 1) - p(-1,-1,-1) - p(-1, 1, 1) - p( 1,-1, 1) - p( 1, 1,-1) + p( 1,-1,-1) + p(-1, 1,-1) + p(-1,-1, 1)) &
           / 8.D0
    a122 = (p( 1, 1, 1) - p(-1,-1,-1) - p(-1, 1, 1) + p( 1,-1, 1) + p( 1, 1,-1) + p( 1,-1,-1) - p(-1, 1,-1) - p(-1,-1, 1)) &
           / 8.D0 - a100 - a120 - a102
    a212 = (p( 1, 1, 1) - p(-1,-1,-1) + p(-1, 1, 1) - p( 1,-1, 1) + p( 1, 1,-1) - p( 1,-1,-1) + p(-1, 1,-1) - p(-1,-1, 1)) &
           / 8.D0 - a010 - a012 - a210
    a221 = (p( 1, 1, 1) - p(-1,-1,-1) + p(-1, 1, 1) + p( 1,-1, 1) - p( 1, 1,-1) - p( 1,-1,-1) - p(-1, 1,-1) + p(-1,-1, 1)) &
           / 8.D0 - a001 - a201 - a021

    ! newton iterate to find maxima
    x = 0.D0  !  initial guess at maximum voxel (z,y,x)
    eps = sqrt(epsilon(1.D0))
    do i = 0, maxIter-1 
    ! compute components of hessian matrix
        xx = x(0) * x(0) 
        yy = x(1) * x(1) 
        zz = x(2) * x(2)
        xy = x(0) * x(1) 
        yz = x(1) * x(2) 
        zx = x(2) * x(0)
        h00 = (a200 + a210 * x(1) + a201 * x(2) + a220 * yy + a202 * zz + a211 * yz + a221 * yy * x(2) + a212 * x(1) * zz + &
               a222 * yy * zz) * 2.D0
        h11 = (a020 + a021 * x(2) + a120 * x(0) + a022 * zz + a220 * xx + a121 * zx + a122 * zz * x(0) + a221 * x(2) * xx + &
               a222 * zz * xx) * 2.D0
        h22 = (a002 + a102 * x(0) + a012 * x(1) + a202 * xx + a022 * yy + a112 * xy + a212 * xx * x(1) + a122 * x(0) * yy + &
               a222 * xx * yy) * 2.D0
        h01 = a110 + a111 * x(2) + a112 * zz + (a210 * x(0) + a120 * x(1) + a211 * zx + a121 * yz + a212 * x(0) * zz + &
              a122 * x(1) * zz + (a220 * xy + a221 * xy * x(2) + a222 * xy * zz) * 2) * 2.D0
        h12 = a011 + a111 * x(0) + a211 * xx + (a021 * x(1) + a012 * x(2) + a121 * xy + a112 * zx + a221 * x(1) * xx + &
              a212 * x(2) * xx + (a022 * yz + a122 * yz * x(0) + a222 * yz * xx) * 2) * 2.D0
        h02 = a101 + a111 * x(1) + a121 * yy + (a102 * x(2) + a201 * x(0) + a112 * yz + a211 * xy + a122 * x(2) * yy + &
              a221 * x(0) * yy + (a202 * zx + a212 * zx * x(1) + a222 * zx * yy) * 2) * 2.D0
        
    ! build inverse of hessian matrix
        det = h00 * h11 * h22 - h00 * h12 * h12 - h11 * h02 * h02 - h22 * h01 * h01 + h01 * h12 * h02 * 2.D0
        i00 = (h11 * h22 - h12 * h12) / det
        i11 = (h22 * h00 - h02 * h02) / det
        i22 = (h00 * h11 - h01 * h01) / det
        i01 = (h02 * h12 - h01 * h22) / det
        i12 = (h01 * h02 - h12 * h00) / det
        i02 = (h12 * h01 - h02 * h11) / det

    ! compute gradient
        d0 = a100 + a110 * x(1) + a101 * x(2) + a120 * yy + a102 * zz + a111 * yz + a121 * yy * x(2) + a112 * x(1) * zz + &
             a122 * yy * zz + x(0) * (a200 + a210 * x(1) + a201 * x(2) + a220 * yy + a202 * zz + a211 * yz + a221 * yy * x(2) + &
             a212 * x(1) * zz + a222 * yy * zz) * 2.D0
        d1 = a010 + a011 * x(2) + a110 * x(0) + a012 * zz + a210 * xx + a111 * zx + a112 * zz * x(0) + a211 * x(2) * xx + &
             a212 * zz * xx + x(1) * (a020 + a021 * x(2) + a120 * x(0) + a022 * zz + a220 * xx + a121 * zx + a122 * zz * x(0) + &
             a221 * x(2) * xx + a222 * zz * xx) * 2.D0
        d2 = a001 + a101 * x(0) + a011 * x(1) + a201 * xx + a021 * yy + a111 * xy + a211 * xx * x(1) + a121 * x(0) * yy + &
             a221 * xx * yy + x(2) * (a002 + a102 * x(0) + a012 * x(1) + a202 * xx + a022 * yy + a112 * xy + a212 * xx * x(1) + &
             a122 * x(0) * yy + a222 * xx * yy) * 2.D0

    ! update x
        step = (/ i00 * d0 + i01 * d1 + i02 * d2, i01 * d0 + i11 * d1 + i12 * d2, i02 * d0 + i12 * d1 + i22 * d2 /)
        x = x - step 

    ! check for convergence
    ! write (*,*) i, x, step 
        maxStep = maxval(abs(step)) 
        if (maxStep.lt.eps) EXIT
        if (i+1.eq.maxIter) x = 0.D0   ! don't interpolate if convergence wasn't reached
    end do

    ! compute interpolated value of maxima
    xx = x(0) * x(0) 
    yy = x(1) * x(1) 
    zz = x(2) * x(2)
    xy = x(0) * x(1) 
    yz = x(1) * x(2) 
    zx = x(2) * x(0)
    vPeak = a000                    + a111 * x(0) * x(1) * x(2) + a222 * xx   * yy   * zz + &
            a100 * x(0)             + a010 * x(1)               + a001 * x(2) + &
            a200 * xx               + a020 * yy                 + a002 * zz + &
            a110 * xy               + a011 * yz                 + a101 * zx + &
            a120 * x(0) * yy        + a012 * x(1) * zz          + a201 * x(2) * xx + &
            a210 * xx   * x(1)      + a021 * yy   * x(2)        + a102 * zz   * x(0) + &
            a220 * xx   * yy        + a022 * yy   * zz          + a202 * zz   * xx + &
            a112 * xy   * x(2)      + a211 * yz   * x(0)        + a121 * zx   * x(1) + &
            a122 * x(0) * yy   * zz + a212 * xx   * x(1) * zz   + a221 * xx   * yy   * x(2)

    ! zyx -> xyz
    xx = x(2)
    x(2) = x(0)
    x(0) = xx

  end function SH_interpolateMaxima


  !--------------------------------------------------------------------------
  !
  ! FUNCTION: SH_zyz2qu
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief   convert ZYZ euler angles to quaternion
  !
  ! @param eu : euler angles to convert to quaternion (Z, Y', Z'')
  ! @param pos: true/false to restrict rotation to (0,pi) (positive w)
  ! @result qu: location to write quaternion as w, x, y, z
  !   
  !> @date 01/30/19 MDG 1.0 original, based on Will Lenthe's classes in sht_xcorr.hpp
  !--------------------------------------------------------------------------
  recursive function SH_zyz2qu(eu, pos) result(qu)
  !DEC$ ATTRIBUTES DLLEXPORT :: SH_zyz2qu

    use mod_global

  implicit none

    real   (kind=dbl),INTENT(IN) :: eu(0:2)
    logical          ,INTENT(IN) :: pos 
    real   (kind=dbl)            :: qu(0:3) 

    real   (kind=dbl)            :: s, c, sigma, delta

    c = 0.5D0 * cos(eu(1))
    s = 0.5D0 * sin(eu(1))
    sigma = 0.5D0 * (eu(2) + eu(0))
    delta = 0.5D0 * (eu(2) - eu(0))

    ! this uses the epsijkd constant to define the 3D rotations convention (see constants.f90 module)
    qu(0) = c * cos(sigma)
    qu(1) = -epsijkd * s * sin(delta)
    qu(2) = -epsijkd * s * cos(delta)
    qu(3) = -epsijkd * s * sin(sigma)

    if ((pos.eqv..TRUE.).and.(qu(0).lt.0.D0)) qu = -qu 

  end function SH_zyz2qu

  !--------------------------------------------------------------------------
  !
  ! SUBROUTINE: SH_conjMult
  !
  !> @author Will Lenthe/Marc De Graef, Carnegie Mellon University
  !
  !> @brief   compute product of ab * cd and ad * conj(cd) without duplicate flops 
  !
  ! @param ab: first  complex number (a + bi)
  ! @param cd: second complex number (c + di)
  ! @param vp: location to write ab *      cd
  ! @param vc: location to write ab * conj(cd)
  ! @return  : agP, agN
  !
  !> @date 03/18/19 MDG 1.0 original, based on Will Lenthe's classes in sht_xcorr.hpp
  !--------------------------------------------------------------------------
  recursive subroutine SH_conjMult(ab, cd, vp, vc) 
  !DEC$ ATTRIBUTES DLLEXPORT :: SH_conjMult

  implicit none

    complex(kind=dbl),INTENT(IN   ) :: ab
    complex(kind=dbl),INTENT(IN   ) :: cd
    complex(kind=dbl),INTENT(INOUT) :: vp
    complex(kind=dbl),INTENT(INOUT) :: vc

    real(kind=dbl)                  :: rr, ri, ir, ii

    rr = real (ab) * real (cd)
    ri = real (ab) * aimag(cd)
    ir = aimag(ab) * real (cd)
    ii = aimag(ab) * aimag(cd)

    vp = cmplx(rr-ii, ir+ri) ! (a + bi) * (c + di)
    vc = cmplx(rr+ii, ir-ri) ! (a + bi) * (c - di)

  end subroutine SH_conjMult

end module mod_SHcorrelator
