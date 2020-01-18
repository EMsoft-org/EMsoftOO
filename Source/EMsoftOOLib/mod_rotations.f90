!--------------------------------------------------------------------------
!
! Function: genrot
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief generate a passive rotation representation, given the unit axis vector and the rotation angle
!
!> @param av 3-component vector (single precision)  
!> @param omega rotation angle (radians)
!>  
! 
!> @date 9/30/14   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function genrot(av,omega) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: genrot

use local
use constants
use error

IMPLICIT NONE

real(kind=sgl),INTENT(IN)       :: av(3)
real(kind=sgl),INTENT(IN)       :: omega

type(orientationtype)           :: res
real(kind=sgl)                  :: axang(4), s

! first make sure that the rotation angle is in the interval [0,pi]
if ((omega.lt.0.0).or.(omega.gt.sngl(cPi))) then 
  call FatalError('rotations:getrot','rotation angle must be in range [0,pi]')
  STOP
endif

! define the axis-angle pair with the correct sign of epsijk
axang(1:3) = -epsijk * av(1:3)
axang(4) = omega

! normalize the unit vector
s = sqrt(sum(av*av))
if (s.ne.0.0) then 
        axang(1:3) = axang(1:3)/s
else
        call FatalError('rotations:getrot','axis unit vector can not be [0,0,0]')
        STOP
endif

res = init_orientation(axang,'ax')

end function genrot

!--------------------------------------------------------------------------
!
! Function: genrot_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief generate a passive rotation representation, given the unit axis vector and the rotation angle
!
!> @param av 3-component vector (double precision)  
!> @param omega rotation angle (radians)
!>  
! 
!> @date 9/30/14   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function genrot_d(av,omega) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: genrot_d

use local
use constants
use error

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: av(3)
real(kind=dbl),INTENT(IN)       :: omega

type(orientationtyped)          :: res
real(kind=dbl)                  :: axang(4), s

! first make sure that the rotation angle is in the interval [0,pi]
if ((omega.lt.0.D0).or.(omega.gt.cPi)) then 
  call FatalError('rotations:getrot_d','rotation angle must be in range [0,pi]')
  STOP
endif

! define the axis-angle pair with the correct sign of epsijk
axang(1:3) = -epsijk * av(1:3)
axang(4) = omega

! normalize the unit vector
s = dsqrt(sum(av*av))
if (s.ne.0.D0) then 
        axang(1:3) = axang(1:3)/s
else
        call FatalError('rotations:getrot','axis unit vector can not be [0,0,0]')
        STOP
endif

res = init_orientation_d(axang,'ax')

end function genrot_d



!--------------------------------------------------------------------------
!
! FUNCTION: RodriguesProduct
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief multiply two Rodrigues vectors
!
!> @param roA first input Rodrigues vector components (single precision)
!> @param roB second input Rodrigues vector components (single precision)
!
!> @date 12/01/16   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function RodriguesProduct(roA,roB) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: RodriguesProduct

use local
use math

real(kind=sgl),INTENT(IN)       :: roA(3)
real(kind=sgl),INTENT(IN)       :: roB(3)

real(kind=sgl)                  :: res(3)

res = (roA + roB - cross3(roA,roB))/(1.0 - DOT_PRODUCT(roA,roB))

end function RodriguesProduct

!--------------------------------------------------------------------------
!
! FUNCTION: RodriguesProduct_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief multiply two Rodrigues vectors
!
!> @param roA first input Rodrigues vector components (single precision)
!> @param roB second input Rodrigues vector components (single precision)
!
!> @date 12/01/16   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function RodriguesProduct_d(roA,roB) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: RodriguesProduct_d

use local
use math

real(kind=dbl),INTENT(IN)       :: roA(3)
real(kind=dbl),INTENT(IN)       :: roB(3)

real(kind=dbl)                  :: res(3)

res = (roA + roB - cross3(roA,roB))/(1.D0 - DOT_PRODUCT(roA,roB))

end function RodriguesProduct_d

!--------------------------------------------------------------------------
!
! FUNCTION: RotVec_om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief rotate a vector using a rotation matrix, active or passive (single precision)
!
!> @details This routine provides a way for the user to transform a vector
!> and it returns the new vector components.  The user can use either a 
!> rotation matrix or a quaternion to define the transformation, and must
!> also specifiy whether an active or passive result is needed.  The quaternion
!> rotation is part of the quaternion.f90 file, in the routine quat_Lp or quat_Lpstar.
!
!> @param vec input vector components (single precision)
!> @param om orientation matrix (single precision)
!> @param ap active/passive switch
!>  
!
!> @date 8/18/14   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function RotVec_om(vec,om,ap) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: RotVec_om

use local

real(kind=sgl),INTENT(IN)       :: vec(3)
real(kind=sgl),INTENT(IN)       :: om(3,3)
character(1),INTENT(IN)         :: ap

real(kind=sgl)                  :: res(3)

if (ap.eq.'p') then
 res = matmul(om,vec)
else
 res = matmul(transpose(om),vec)
end if

end function RotVec_om

!--------------------------------------------------------------------------
!
! FUNCTION: RotVec_om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief rotate a vector using a rotation matrix, active or passive (double precision)
!
!> @details This routine provides a way for the user to transform a vector
!> and it returns the new vector components.  The user can use either a 
!> rotation matrix or a quaternion to define the transformation, and must
!> also specifiy whether an active or passive result is needed.  The quaternion
!> rotation is part of the quaternion.f90 file, in the routine quat_Lp or quat_Lpstar.
!
!> @param vec input vector components (double precision)
!> @param om orientation matrix (double precision)
!> @param ap active/passive switch
!>  
!
!> @date 8/18/14   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function RotVec_om_d(vec,om,ap) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: RotVec_om_d

use local

real(kind=dbl),INTENT(IN)       :: vec(3)
real(kind=dbl),INTENT(IN)       :: om(3,3)
character(1),INTENT(IN)         :: ap

real(kind=dbl)                  :: res(3)

if (ap.eq.'p') then
 res = matmul(om,vec)
else
 res = matmul(transpose(om),vec)
end if

end function RotVec_om_d

!--------------------------------------------------------------------------
!
! FUNCTION: RotTensor2_om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief rotate a second rank tensor using a rotation matrix, active or passive (single precision)
!
!> @param tensor input tensor components (single precision)
!> @param om orientation matrix (single precision)
!> @param ap active/passive switch
!>  
!
!> @date 8/18/14   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function RotTensor2_om(tensor,om,ap) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: RotTensor2_om

use local

real(kind=sgl),INTENT(IN)       :: tensor(3,3)
real(kind=sgl),INTENT(IN)       :: om(3,3)
character(1),INTENT(IN)         :: ap

real(kind=sgl)                  :: res(3,3)

if (ap.eq.'p') then
 res = matmul(matmul(om,tensor),transpose(om))
else
 res = matmul(matmul(transpose(om),tensor),om)
end if

end function RotTensor2_om


!--------------------------------------------------------------------------
!
! FUNCTION: RotTensor2_om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief rotate a second rank tensor using a rotation matrix, active or passive (double precision)
!
!> @param tensor input tensor components (double precision)
!> @param om orientation matrix (double precision)
!> @param ap active/passive switch
!>  
!
!> @date 8/18/14   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function RotTensor2_om_d(tensor,om,ap) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: RotTensor2_om_d

use local

real(kind=dbl),INTENT(IN)       :: tensor(3,3)
real(kind=dbl),INTENT(IN)       :: om(3,3)
character(1),INTENT(IN)         :: ap

real(kind=dbl)                  :: res(3,3)

if (ap.eq.'p') then
 res = matmul(matmul(om,tensor),transpose(om))
else
 res = matmul(matmul(transpose(om),tensor),om)
end if

end function RotTensor2_om_d


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! some basic averaging routines for quaternions
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! function: quat_average
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  computes the geometrical mean of a list of quaternions using the quaternion logarithm
! 
!> @param qlist quaternion list
!> @param numq number of quaternions in list
!> @param qstdev standard deviation quaternion
! 
!> @date 03/16/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function quat_average(qlist,numq,qstdev) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: quat_average

use local

real(kind=sgl),INTENT(IN)       :: qlist(4,numq)
integer(kind=irg)               :: numq
real(kind=sgl),INTENT(OUT)      :: qstdev(4)
real(kind=sgl)                  :: res(4)

integer(kind=irg)               :: i
real(kind=sgl)                  :: lsum(3), ax(4), qv, sqv, dfm(3)
real(kind=sgl)                  :: axanglist(3,numq)

! convert each one to a logarithm, which is really the unit rotation vector multiplied by half the rotation angle
do i=1,numq
! convert each quaternion to an axis angle pair
  ax = qu2ax(qlist(1:4,i))
  axanglist(1:3,i) =  ax(1:3)*0.5*ax(4)
end do

! compute the geometric mean
lsum = sum(axanglist,2)/float(numq)

! then get the standard deviation from the mean
dfm = 0.0
do i=1,numq
  dfm(1:3) = dfm(1:3) + (lsum - axanglist(1:3,i))**2
end do
dfm = sqrt(dfm/float(numq))

! then convert this average back to a quaternion via the exponentiation operation
qv = sqrt(sum(lsum*lsum))
sqv = sin(qv)/qv
res = (/ cos(qv), lsum(1)*sqv, lsum(2)*sqv, lsum(3)*sqv /)

! and do the same with the standard deviation quaternion
qv = sqrt(sum(dfm*dfm))
sqv = sin(qv)/qv
qstdev = (/ cos(qv), dfm(1)*sqv, dfm(2)*sqv, dfm(3)*sqv /)

end function quat_average

!--------------------------------------------------------------------------
!
! function: quat_average_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  computes the geometrical mean of a list of quaternions using the quaternion logarithm
! 
!> @param qlist quaternion list
!> @param numq number of quaternions in list
!> @param qstdev standard deviation quaternion
! 
!> @date 03/16/15 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function quat_average_d(qlist,numq,qstdev) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: quat_average_d

use local

real(kind=dbl),INTENT(IN)       :: qlist(4,numq)
integer(kind=irg)               :: numq
real(kind=dbl),INTENT(OUT)      :: qstdev(4)
real(kind=dbl)                  :: res(4)

integer(kind=irg)               :: i
real(kind=dbl)                  :: lsum(3), ax(4), qv, sqv, dfm(3)
real(kind=dbl)                  :: axanglist(3,numq)

! convert each one to a logarithm, which is really the unit rotation vector multiplied by half the rotation angle
do i=1,numq
! convert the quaternion to an axis angle pair
  ax = qu2ax_d(qlist(1:4,i))
  axanglist(1:3,i) =  ax(1:3)*0.5D0*ax(4)
end do

! compute the geometric mean
lsum = sum(axanglist,2)/dble(numq)

! then get the standard deviation from the mean
dfm = 0.D0
do i=1,numq
  dfm(1:3) = dfm(1:3) + (lsum - axanglist(1:3,i))**2
end do
dfm = dsqrt(dfm/dble(numq))

! then convert this average back to a quaternion via the exponentiation operation; check interval first ?
qv = dsqrt(sum(lsum*lsum))
sqv = dsin(qv)/qv
res = (/ dcos(qv), lsum(1)*sqv, lsum(2)*sqv, lsum(3)*sqv /)

! and do the same with the standard deviation quaternion
qv = dsqrt(sum(dfm*dfm))
sqv = dsin(qv)/qv
qstdev = (/ dcos(qv), dfm(1)*sqv, dfm(2)*sqv, dfm(3)*sqv /)

end function quat_average_d

