program tester 

use mod_global 
use mod_kinds
use mod_EMsoft
! use mod_symmetry
! use mod_crystallography
! use mod_QCsymmetry
! use mod_QCcrystallography
use mod_io
use mod_dirstats
! use mod_dualquaternions
use mod_quaternions
use mod_rotations
! use mod_octonions
! use mod_GBoctonions
! use mod_HDFsupport
! use HDF5
! use mod_vendors
! use mod_HallSG
! use mod_math
! use mod_platformsupport
! use mod_PGA3D
! use mod_PGA3Dsupport
use mod_axonometry
use mod_postscript


IMPLICIT NONE 

! type(PGA3D_T)             :: mv_plane, mv_line, mv, pt, mv_pp
! real(kind=dbl)            :: L, a, b, c, d, alpha, x,y,z, ord, sa, ca

! type(axonometry_T)          :: AXO 
! type(Postscript_T)          :: PS 
type(EMsoft_T)              :: EMsoft

type(c_T)                   :: cu 
type(q_T)                   :: qa, qb 

integer(kind=irg)           :: offsets(3,12), pnum, N 
integer(kind=irg)           :: nx, ny, i, j , k
real(kind=dbl)              :: delta 



! type(DirStat_T)             :: DSvmf, DSwat

! type(HDF_T)             :: HDF
! type(Vendor_T)          :: VT
! type(HallSG_T)          :: HSG 
! type(SpaceGroup_T)      :: SSG

! character(fnlen)        :: fname, groupname, inputtype, progname, progdesc, HDFstrings(10) 
! integer(kind=irg)       :: hdferr, itype, istat, ipf_wd, ipf_ht, sz(3), L, recordsize, &
!                            patsz, i, j, numsx, numsy, correctsize, s1, s2,HSGn, info, status
! real(kind=sgl),allocatable   :: exppatarray(:), tot(:), totold(:)
! real(kind=dbl),allocatable   :: SG(:,:,:)
! integer(HSIZE_T)        :: dims3(3), offset3(3)
! character(16)           :: HS
! real(kind=dbl),allocatable          :: SGdirec(:,:,:)
! real(kind=dbl)          :: z(11,11), fit(5), mp1, mp2, sig1, sig2 

! integer(C_INT32_T)  :: res

! type(Octonion_T)        :: o
! ! type(Quaternion_T)      :: qu1, qu2
! type(GBOctonion_T)      :: gb, newgb, sRL, sLR, sLL
! type(o_T)               :: Ra, Rb, Rpixs, tmpa, tmpb, ttt
! type(q_T)               :: qa, qb, qpixs, newa, newb
! type(Quaternion_T)      :: quat_pixs, quat_a, quat_b, quat_newa, quat_newb 
! real(kind=dbl)          :: mat_a(3,3), mat_b(3,3), mat_pixs(3,3), mat_inv(3,3), &
!                            mat_sigxs(3,3), mat_tmp(3,3)


! real(kind=dbl)          :: diffd
! real(kind=sgl)          :: diff

! ! threshold values 
! real(kind=dbl),parameter:: epsd = 1.0D-12
! real(kind=sgl),parameter:: eps  = 1.0E-7

! ! various parameters
! integer(kind=irg)       :: i, errcnt 

! type(Octonion_T) :: resdsum
! type(Octonion_T) :: resdsub
! type(Octonion_T) :: resdmult 
! type(Octonion_T) :: resdsmult 
! type(Octonion_T) :: resddiv 
! type(Octonion_T) :: resdconjg 
! type(Octonion_T) :: resdinv 
! type(Octonion_T) :: resdzero
! real(kind=dbl)   :: resdabs

! type(Octonion_T) :: resssum
! type(Octonion_T) :: resssub
! type(Octonion_T) :: ressmult 
! type(Octonion_T) :: resssmult 
! type(Octonion_T) :: ressdiv 
! type(Octonion_T) :: ressconjg 
! type(Octonion_T) :: ressinv 
! type(Octonion_T) :: resszero
! real(kind=sgl)   :: ressabs

! computation of average disorientations for fcc sampling of cubochoric space
! offset(1:3,1) = (/ 1, 1, 0 /)
! offset(1:3,2) = (/ 1,-1, 0 /)
! offset(1:3,3) = (/-1, 1, 0 /)
! offset(1:3,4) = (/-1,-1, 0 /)
! offset(1:3,5) = (/ 1, 0, 1 /)
! offset(1:3,6) = (/ 1, 0,-1 /)
! offset(1:3,7) = (/-1, 0, 1 /)
! offset(1:3,8) = (/-1, 0,-1 /)
! offset(1:3,9) = (/ 0, 1, 1 /)
! offset(1:3,10) = (/ 0, 1,-1 /)
! offset(1:3,11) = (/ 0,-1, 1 /)
! offset(1:3,12) = (/ 0,-1,-1 /)

! ! number of points and cubochoric grid spacing
! N = 20
! delta = cPi**(2.D0/3.D0) / dble(2*N)
! pnum = 4*(2*N)**3



















! DSvmf = DirStat_T(DStype='VMF')
! call DSvmf%UnitTests( 1, 100, 'VMF-logCp.txt')

! DSwat = DirStat_T(DStype='WAT')
! call DSwat%UnitTests( 1, 100, 'WAT-logCp.txt')








! ! simple test of the axonometry module
! progname = ' x '
! progdesc = ' y '
! axname = 'axotest.eps'
! EMsoft = EMsoft_T( progname, progdesc )
! PS = Postscript_T( progdesc, EMsoft, imanum = 1, dontask = .TRUE., psname = axname )
! AXO = axonometry_T( progdesc, axw = 6.5, xll = 3.5, yll = 3.0 )

! nx = 100
! ny = 150
! allocate( zz(nx,ny), x(nx), y(ny) )

! x = (/ (real(i), i=1,nx) /)/ real(nx) - 0.5
! y = (/ (real(i), i=1,ny) /)/ real(ny) - 0.5

! do i=1,nx
!     do j=1,ny
!         zz(i,j) = 15.0 * exp(- (x(i)**2+y(j)**2) * 100.0 )
!     end do 
! end do

! zz = cshift(zz, 15, 1)
! zz = zz - cshift(zz, -30, 1) 

! write (*,*) ' range = ', minval(zz), maxval(zz)

! g = 1.0
! call AXO%axonometry(PS,EMsoft,zz,nx,ny,g,axname)







! quick trial of the crystallographic and holomorphic chirality concepts
! for Grain Boundaries ... 

! quat_pixs = Quaternion_T( qd = (/ 0.D0, 1.D0, 0.D0, 0.D0 /) ) 
! qpixs = q_T( qdinp = (/ 0.D0, 1.D0, 0.D0, 0.D0 /) )
! o = Octonion_T( od = (/ 0.99513333D0,0.000000D0, 0.098537618D0,0.000000D0, &
!                         0.99513333D0,0.000000D0,-0.098537618D0,0.000000D0 /))
! gb = GBOctonion_T( oct = o )

! ! set up the quaternions and rotation matrices
! quat_a = gb%GBO_get_q(1)
! quat_b = gb%GBO_get_q(2)

! qa = q_T( qdinp = quat_a%get_quatd() )
! qb = q_T( qdinp = quat_b%get_quatd() )

! Ra = qa%qo()
! Rb = qb%qo()
! Rpixs = qpixs%qo()

! call Ra%o_print(' Ra : ')
! call Rb%o_print(' Rb : ')
! call Rpixs%o_print(' Rpixs : ')

! mat_a = Ra%o_copyd()
! mat_b = Rb%o_copyd()
! mat_pixs = Rpixs%o_copyd()
! mat_inv = 0.D0 
! mat_inv(1,1) = -1.0D0
! mat_inv(2,2) = -1.0D0
! mat_inv(3,3) = -1.0D0

! mat_sigxs = mat_inv 
! mat_sigxs(2,2) = 1.D0
! mat_sigxs(3,3) = 1.D0

! mat_pixs = -mat_sigxs

! write (*,*) ' Misorientation angle : ', 2.D0*acos( sum(quat_a%get_quatd()*quat_b%get_quatd()))/dtor
! call gb%oct_print(' input octonion : ')
! write (*,*) ' '

! ! if gb is in Sigma_RR, then what are the others ?
! write (*,*) ' Sigma_RR analysis '
! write (*,*) mat_pixs
! mat_tmp = matmul(mat_a,mat_pixs)
! tmpa = o_T( odinp = mat_tmp )
! mat_tmp = matmul(mat_b,mat_pixs)
! tmpb = o_T( odinp = mat_tmp )
! call tmpa%o_print(' Ra x Rpixs : ')
! call tmpb%o_print(' Rb x Rpixs : ')
! newa = tmpa%oq()
! newb = tmpb%oq()
! call newa%q_print(' newa :')
! call newb%q_print(' newb :')
! o = Octonion_T( od = (/ newa%q_copyd(), newb%q_copyd() /) )
! newgb = GBOctonion_T( oct = o )
! write (*,*) ' Misorientation angle : ', 2.D0*acos( sum(newa%q_copyd()*newb%q_copyd()))/dtor
! call newgb%oct_print(' starred quaternion: ')
! write (*,*) ' '

! ! Sigma_RL
! write (*,*) ' Sigma_RL analysis '
! mat_tmp = matmul(mat_inv, mat_b)
! tmpa = o_T( odinp = mat_a )
! tmpb = o_T( odinp = mat_tmp )
! call tmpb%o_print(' I x Rb : ')
! newa = tmpa%oq()
! newb = tmpb%oq()
! call newa%q_print(' newa :')
! call newb%q_print(' newb :')

! ! Sigma_LR
! write (*,*) ' Sigma_RL analysis '
! mat_tmp = matmul(mat_inv, mat_a)
! tmpa = o_T( odinp = mat_tmp )
! tmpb = o_T( odinp = mat_a )
! call tmpb%o_print(' I x Ra : ')
! newa = tmpa%oq()
! newb = tmpb%oq()
! call newa%q_print(' newa :')
! call newb%q_print(' newb :')

! ! Sigma_LL
! write (*,*) ' Sigma_RL analysis '
! mat_tmp = matmul(mat_inv, mat_a)
! tmpa = o_T( odinp = mat_tmp )
! mat_tmp = matmul(mat_inv, mat_b)
! tmpb = o_T( odinp = mat_tmp )
! call tmpa%o_print(' I x Ra : ')
! call tmpb%o_print(' I x Rb : ')
! newa = tmpa%oq()
! newb = tmpb%oq()
! call newa%q_print(' newa :')
! call newb%q_print(' newb :')

! !===================================================
! ! set the reference values (verified against values on <https://pypi.org/project/pyoctonion/#description>)

! ! correct answers for the Octonion_T class
! resdzero = Octonion_T()
! resdsum = Octonion_T( od = (/  2.0D0,   5.0D0,   8.0D0,  11.0D0,  14.0D0,   8.0D0,  11.0D0,  14.0D0 /) )
! resdsub = Octonion_T( od = (/  0.0D0,   -1.0D0,   -2.0D0,   -3.0D0,   -4.0D0,    4.0D0,    3.0D0,    2.0D0 /) )
! resdsmult = Octonion_T( od = (/ 1.41421356237310D0, 2.82842712474619D0, 4.24264068711929D0, 5.65685424949238D0, &
!                                 7.07106781186548D0, 8.48528137423857D0, 9.89949493661167D0,11.31370849898476D0 /) )
! resdmult = Octonion_T( od = (/-181.0D0,  -48.0D0,  -17.0D0,  -40.0D0,   83.0D0,    0.0D0,   35.0D0,    4.0D0 /) )
! resddiv = Octonion_T( od = (/  0.82805429864253D0,    0.23529411764706D0,    0.10407239819005D0,    0.21719457013575D0,  &
!                               -0.33031674208145D0,    0.05429864253394D0,   -0.09502262443439D0,    0.05429864253394D0 /) )
! resdabs = 14.282856857085701D0
! resdconjg = Octonion_T( od = (/  1.0D0,   -2.0D0,   -3.0D0,   -4.0D0,   -5.0D0,   -6.0D0,   -7.0D0,   -8.0D0 /) )
! resdinv = Octonion_T( od = (/  0.00490196078431D0,  -0.00980392156863D0,  -0.01470588235294D0,  -0.01960784313725D0,  &
!                               -0.02450980392157D0,  -0.02941176470588D0,  -0.03431372549020D0, -0.03921568627451D0 /) )

! resszero = Octonion_T()
! resssum = Octonion_T( o = (/  2.00,   5.00,   8.00,  11.00,  14.00,   8.00,  11.00,  14.00 /) )
! resssub = Octonion_T( o = (/  0.00,   -1.00,   -2.00,   -3.00,   -4.00,    4.00,    3.00,    2.00 /) )
! resssmult = Octonion_T( o = (/ 1.414214, 2.828427, 4.242640, 5.656854, 7.071068, 8.485281, 9.899495, 11.313708 /) )
! ressmult = Octonion_T( o = (/-181.00,  -48.00,  -17.00,  -40.00,   83.00,    0.00,   35.00,    4.00 /) )
! ressdiv = Octonion_T( o = (/  0.828054,   0.235294,   0.104072,   0.217195,  -0.330317,   0.054299,  -0.095023,   0.054299 /) )
! ressabs  = 14.2828569
! ressconjg = Octonion_T( o = (/  1.0,   -2.0,   -3.0,   -4.0,   -5.0,   -6.0,   -7.0,   -8.0 /) )
! ressinv = Octonion_T( o = (/  0.004902,  -0.009804,  -0.014706,  -0.019608,  -0.024510,  -0.029412,  -0.034314,  -0.039216 /) )

! ! initialize the error identifier to zero (should remain zero upon successful exit)
! res = 0

! !===================================================
! !=============Double Precision Tests================
! !===================================================
! ! call set_octonionprecision('d')
! ! call set_octonionGBmode(.FALSE.)

! !===================================================
! ! initialize zero quaternion 
! u = Octonion_T( od = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
! diffd = cabs(u)
! if (diffd.gt.epsd) then 
!   res = 1
!   write (*,"('double precision zero initialization test failed = ',D18.10)") diff
!   return
! end if

! if (.not.(u%octsequal(resdzero))) then 
!   res = 2
!   write (*,"('double precision zero comparison test failed = ',D18.10)") diff
!   return
! end if

! ! arithmetic tests 
! a = Octonion_T( od = (/ 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0 /) )
! b = Octonion_T( od = (/ 1.D0, 3.D0, 5.D0, 7.D0, 9.D0, 2.D0, 4.D0, 6.D0 /) )
! c = Octonion_T( od = (/ 8.D0, 7.D0, 6.D0, 5.D0, 4.D0, 3.D0, 2.D0, 1.D0 /) )

! d = a+b
! if (.not.(d%octsequal(resdsum))) then 
!   res = 3
!   write (*,"('double precision addition test failed = ')") 
!   return
! end if

! d = a-b
! if (.not.(d%octsequal(resdsub))) then 
!   res = 4
!   write (*,"('double precision subtraction test failed = ')")
!   return
! end if

! d = a*sqrt(2.D0)
! if (.not.(d%octsequal(resdsmult))) then 
!   res = 5
!   write (*,"('double precision scalar multiplication test failed = ')")
!   return
! end if

! d = a*b
! if (.not.(d%octsequal(resdmult))) then 
!   res = 6
!   write (*,"('double precision quaternion multiplication test failed = ')") 
!   return
! end if

! d = a/b
! if (.not.(d%octsequal(resddiv))) then 
!   res = 7
!   write (*,"('double precision division test failed = ')")
!   return
! end if

! d = conjg(a)
! if (.not.(d%octsequal(resdconjg))) then 
!   res = 8
!   write (*,"('double precision conjugation test failed = ')") 
!   return
! end if

! diffd = abs(cabs(a) - resdabs)
! if (diffd.gt.epsd) then 
!   res = 9
!   write (*,"('double precision norm test failed = ',D18.10)") diffd
!   return
! end if

! ! d = a%octinverse()
! ! if (.not.(d%octsequal(resdinv))) then 
! !   res = 10
! !   write (*,"('double precision inverse test failed = ')") 
! !   return
! ! end if

! ! from here on we work with unit quaternions
! call a%o_normalize()
! diffd = abs(cabs(a) - 1.D0)
! if (diffd.gt.epsd) then 
!   res = 11
!   write (*,"('double precision normalization test failed = ',D18.10)") diffd
!   return
! end if

! ! 
! ! call set_octonionGBmode(.TRUE.)
! ! call a%octnormalize()
! ! call a%oct_print(' GBOM normalization : ')
! ! diffd = abs(cabs(a) - 1.D0)
! ! if (diffd.gt.epsd) then 
! !   res = 12
! !   write (*,"('double precision normalization test failed = ',D18.10)") diffd
! !   return
! ! end if
! ! call set_octonionGBmode(.FALSE.)

! !===================================================
! !=============Single Precision Tests================
! !===================================================
! ! call set_octonionprecision('s')

! !===================================================
! ! initialize zero quaternion 
! u = Octonion_T( o = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /) )
! diffd = cabs(u)
! if (diffd.gt.eps) then 
!   res = 13
!   write (*,"('single precision zero initialization test failed = ',D18.10)") diff
!   return
! end if

! if (.not.(u%octsequal(resszero))) then 
!   res = 14
!   write (*,"('single precision zero comparison test failed = ',D18.10)") diff
!   return
! end if

! ! arithmetic tests 
! a = Octonion_T( o = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /) )
! b = Octonion_T( o = (/ 1.0, 3.0, 5.0, 7.0, 9.0, 2.0, 4.0, 6.0 /) )
! c = Octonion_T( o = (/ 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /) )
! d = Octonion_T( )

! d = a+b
! if (.not.(d%octsequal(resssum))) then 
!   res = 15
!   write (*,"('single precision addition test failed = ')") 
!   return
! end if

! d = a-b
! if (.not.(d%octsequal(resssub))) then 
!   res = 16
!   write (*,"('single precision subtraction test failed = ')")
!   return
! end if

! d = a*sqrt(2.D0)
! if (.not.(d%octsequal(resssmult))) then 
!   res = 17
!   write (*,"('single precision scalar multiplication test failed = ')")
!   return
! end if

! d = a*b
! if (.not.(d%octsequal(ressmult))) then 
!   res = 18
!   write (*,"('single precision quaternion multiplication test failed = ')") 
!   return
! end if

! d = a/b
! if (.not.(d%octsequal(ressdiv))) then 
!   res = 19
!   write (*,"('single precision division test failed = ')")
!   return
! end if

! d = conjg(a)
! if (.not.(d%octsequal(ressconjg))) then 
!   res = 20
!   write (*,"('single precision conjugation test failed = ')") 
!   return
! end if

! diff = abs(cabs(a) - ressabs)
! if (diff.gt.eps) then 
!   res = 21
!   write (*,"('single precision norm test failed = ',D18.10)") diffd
!   return
! end if

! ! d = a%octinverse()
! ! if (.not.(d%octsequal(ressinv))) then 
! !   res = 22
! !   write (*,"('single precision inverse test failed = ')") 
! !   return
! ! end if

! ! from here on we work with unit quaternions
! call a%o_normalize()
! diff = abs(cabs(a) - 1.0)
! if (diff.gt.eps) then 
!   res = 23
!   write (*,"('single precision normalization test failed = ',D18.10)") diffd
!   return
! end if

! !
! ! call set_octonionGBmode(.TRUE.)
! ! call a%octnormalize()
! ! call a%oct_print(' GBOM normalization : ')
! ! diff = abs(cabs(a) - 1.0)
! ! if (diff.gt.eps) then 
! !   res = 24
! !   write (*,"('single precision normalization test failed = ',D18.10)") diffd
! !   return
! ! end if
! ! call set_octonionGBmode(.FALSE.)


! write (*,*) ' if we get here then all tests are correctly performed '

! ! short GBoctonions test 
! qu1 = Quaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0 /) )
! qu2 = Quaternion_T( qd = (/ 0.D0, 0.D0, 1.D0, 0.D0 /) )

! gb = GBOctonion_T( qu1, qu2 )
! call gb%oct_print(' this should be a normalized octonion from two quaternions')









! L = 1000.D0 
! alpha = cPi*0.5D0 - 70.D0*dtor + 10.D0 * dtor
! sa = sin(alpha)
! ca = cos(alpha)
! a = L * sa
! b = 0.D0 
! c = L * ca
! d = L*L

! call PGA3D_initialize()

! mv_line = line(0.D0,0.D0,c)
! mv_line = mv_line%normalized()
! mv_plane = plane(a,b,c,d)
! mv_pp = plane(a,b,c,0.D0)
! call mv_plane%log('plane')
! call mv_line%log('line')

! mv = meet(mv_plane, mv_line)
! mv = mv%normalized()
! call mv%log('intersection')

! call getpoint(mv,x,y,z)
! write (*,*) ' point in detector frame : ',  y, ca*(x+a)-sa*(z+c), -sa*(x+a)-ca*(z+c)
! write (*,*) ' oriented distance to plane : ', ordisttoplane(mv,mv_plane)
! write (*,*) ' oriented distance to plane : ', ordisttoplane(mv,mv_pp)


! status = system_hostnm(fname)
! write (*,*) 'output of subroutine : ', trim(fname)

! info = system_hostnm(fname)
! write (*,*) 'output of function : ', trim(fname), info


! stop 


! mp1 = 5.D0 
! mp2 = 5.D0 
! sig1 = 0.5D0 
! sig2 = 0.5D0 

! z = reshape( (/  0.00001051,  0.00004493,  0.00012878,  0.00024743,  0.00031868,  0.00027512,  0.00015921,  0.00006176,  &
!   0.00001606,  0.00000280,  0.00000033,&
!   0.00010879,  0.00046518,  0.00133337,  0.00256190,  0.00329956,  0.00284861,  0.00164851,  0.00063949, &
!   0.00016629,  0.00002898,  0.00000339,&
!   0.00075503,  0.00322858,  0.00925424,  0.01778085,  0.02290057,  0.01977072,  0.01144144,  0.00443835, & 
!   0.00115410,  0.00020116,  0.00002350,&
!   0.00351266,  0.01502047,  0.04305395,  0.08272269,  0.10654143,  0.09198025,  0.05322956,  0.02064873, & 
!   0.00536928,  0.00093588,  0.00010935,&
!   0.01095445,  0.04684227,  0.13426640,  0.25797583,  0.33225605,  0.28684612,  0.16599969,  0.06439435, & 
!   0.01674443,  0.00291861,  0.00034101,&
!   0.02289956,  0.09792069,  0.28067513,  0.53928161,  0.69455954,  0.59963306,  0.34701149,  0.13461218, & 
!   0.03500314,  0.00610115,  0.00071285,&
!   0.03208825,  0.13721237,  0.39329890,  0.75567387,  0.97325867,  0.84024196,  0.48625340,  0.18862669, & 
!   0.04904851,  0.00854930,  0.00099889,&
!   0.03014026,  0.12888260,  0.36942284,  0.70979904,  0.91417489,  0.78923324,  0.45673433,  0.17717570, & 
!   0.04607092,  0.00803030,  0.00093825,&
!   0.01897712,  0.08114794,  0.23259852,  0.44690851,  0.57558902,  0.49692241,  0.28757218,  0.11155457, & 
!   0.02900748,  0.00505609,  0.00059075,&
!   0.00800932,  0.03424861,  0.09816855,  0.18861839,  0.24292819,  0.20972683,  0.12137026,  0.04708177, & 
!   0.01224265,  0.00213393,  0.00024933,&
!   0.00226591,  0.00968926,  0.02777282,  0.05336195,  0.06872671,  0.05933373,  0.03433681,  0.01331988, & 
!    0.00346356,  0.00060371,  0.00007054 /), (/11,11/) )

! z = transpose(z)

! fit = (/ z(6,6), mp1, mp2, sig1, sig2 /)
! call GaussianFit( 11, z, fit, info)

! write (*,*) 'fit = ', fit(1), fit(2) - 5.D0, fit(3)-5.D0, fit(4), fit(5) 



! stop

! HS = List_Hall_Symbols(62, HSGn)

! write (*,*) ' Hall Space Group number = ', HSGn, trim(HS)

! SSG = SpaceGroup_T( SGnumber = 62, useHall=.TRUE., HallSGnumber=HSGn )

! SGdirec = SSG%getSpaceGroupPGdirecMatrices()

! sz = shape(SGdirec)

! do i=1,sz(1) 
!   write(*,*) 'pg matrix ',i
!   do j=1,3
!     write (*,*) SGdirec(i,j,1:3)
!   end do 
! end do

! stop
! HSG = HallSG_T( HS )

! numsx = HSG%get_NHallgenerators()
! write (*,*) 'number of generators in '//trim(HS)//' : ', numsx
! allocate(SG(4,4,numsx))

! SG = HSG%get_Hall_SeitzGenerators()

! do i=1,numsx 
!   do j=1,4
!     write (*,*)  SG(j,1:4,i)
!   end do 
!   write (*,*) '-----'
! end do 


! stop
! ! HSG = HallSG_T( '-P 1', verbose=.TRUE. )

! ! HSG = HallSG_T( '-I 2xb', verbose=.TRUE. )
! ! HSG = HallSG_T( '-I 2xb (0 0 1)', verbose=.TRUE. )

! HSG = HallSG_T( '-P 31 2c', verbose=.TRUE. )
! HSG = HallSG_T( 'P 31 2c (0 0 1)', verbose=.TRUE. )

! stop


! HSG = HallSG_T( 'P 2 2ab -1ab', verbose=.TRUE. )
! HSG = HallSG_T( 'P 4ab 2ab -1ab', verbose=.TRUE. )
! HSG = HallSG_T( '-F 4 2 3', verbose=.TRUE. )
! HSG = HallSG_T( 'F 4d 2 3 -1cd', verbose=.TRUE. )

! stop
! stop
! progname = 'tester'
! progdesc = 'test program to read problematic HDF5 file'
! EMsoft = EMsoft_T( progname, progdesc)

! ! open the HDF interface
! call openFortranHDFInterface()
! HDF = HDF_T()
! fname = 'playarea/Oxford/Al-large.h5'
! fname = EMsoft%generateFilePath('EMdatapathname',trim(fname))

! inputtype = 'TSLHDF'
! VT = Vendor_T( inputtype )
! itype = VT%get_itype()
! call VT%set_filename(fname)

! ipf_wd = 750
! ipf_ht = 500
! numsx = 156 
! numsy = 128
! L = numsx*numsy 
! correctsize = 16*ceiling(float(L)/16.0)
! recordsize = correctsize*4
! patsz = correctsize

! HDFstrings = ''
! HDFstrings(1) = '1'
! HDFstrings(2) = 'EBSD'
! HDFstrings(3) = 'Data'
! HDFstrings(4) = 'Processed Patterns'
! ! open the pattern file
! istat = VT%openExpPatternFile(EMsoft, ipf_wd, L, recordsize, HDFstrings, HDF)

! allocate(exppatarray(patsz * ipf_wd), tot(ipf_wd), totold(ipf_wd),stat=istat)

! dims3 = (/ numsx, numsy, ipf_wd /)

! write (*,*) L, correctsize, recordsize, patsz 

! do i=250,275 ! ipf_ht
!   exppatarray = 0.0
!   offset3 = (/ 0, 0, (i-1)*ipf_wd /)
!   call VT%getExpPatternRow(i, ipf_wd, patsz, L, dims3, offset3, exppatarray, &
!                                      HDFstrings=HDFstrings, HDF=HDF)
!   write (*,*) 'row number = ', i
!   tot = 0.0
!   do j=1,ipf_wd
!     s1 = (j-1)*patsz+1
!     s2 = j*patsz
!     tot(j) = sum(exppatarray( s1:s2 ))
!   end do
!   if (i.eq.262) then 
!     write (*,*) tot(19960:patsz) 
!   else 
!     write (*,*) tot(19960:patsz) - totold(19960:patsz)
!   end if 
!   write (*,*) '--------------' 
!   totold = tot
! end do 

! call VT%closeExpPatternFile()





! type(SpaceGroup_T)          :: SG 
! type(QCSpaceGroup_T)        :: QCSG 
! type(IO_T)                  :: Message
! integer(kind=irg)           :: i, j, isg, SamplingType, isym

! type(DualQuaternion_T)      :: dq1, dq2, dqp, tmp
! real(kind=dbl)              :: v(3)

! type(Octonion_T)                :: a, b, c, d 
! real(kind=dbl)                  :: onorm
! real(kind=sgl)                  :: onorms

! a = Octonion_T( od = (/ 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0 /) )
! b = Octonion_T( od = (/ 1.D0, 3.D0, 5.D0, 7.D0, 9.D0, 2.D0, 4.D0, 6.D0 /) )
! c = Octonion_T( od = (/ 8.D0, 7.D0, 6.D0, 5.D0, 4.D0, 3.D0, 2.D0, 1.D0 /) )

! d = a+b 
! call d%oct_print('a+b = ')

! d = a-b 
! call d%oct_print('a-b = ')

! d = a*sqrt(2.D0)
! call d%oct_print('a*s = ')

! d = a*b
! call d%oct_print('a*b = ')

! d = a/b 
! call d%oct_print('a/b = ')

! onorm = cabs(a)
! write (*,*) 'a%norm = ', onorm

! d = conjg(a)
! call d%oct_print('conj(a) = ')

! d = a%octinverse()
! call d%oct_print('a%inv = ')

! a = Octonion_T( o = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /) )
! b = Octonion_T( o = (/ 1.0, 3.0, 5.0, 7.0, 9.0, 2.0, 4.0, 6.0 /) )
! c = Octonion_T( o = (/ 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /) )

! d = a+b 
! call d%oct_print('a+b = ')

! d = a-b 
! call d%oct_print('a-b = ')

! d = a*sqrt(2.D0)
! call d%oct_print('a*s = ')

! d = a*b 
! call d%oct_print('a*b = ')

! d = a/b 
! call d%oct_print('a/b = ')

! onorms = real(cabs(a))
! write (*,*) 'a%norm = ', onorm

! d = conjg(a)
! call d%oct_print('conj(a) = ')

! d = a%octinverse()
! call d%oct_print('a%inv = ')



! ! simple test of the dual quaternion package 
! dq1 = DualQuaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
! dq2 = DualQuaternion_T( qd = (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0 /) )

! dqp = dq1 + dq2 
! call dqp%dualquat_print()

! dqp = dq1 * dq2 
! call dqp%dualquat_print()

! dqp = dq1 / dq2 
! call dqp%dualquat_print()

! dqp = DualQuaternion_T( qd = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
! call dqp%generate_dualquat( cPi/2.D0, (/ 0.D0, 0.D0, 1.D0/), (/1.D0, 0.D0, 0.D0/), 5.D0 )
! ! dq2 = DualQuaternion_T( qd = (/ sqrt(3.D0), 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /) )
! ! dqp = dq2 * dqp
! call dqp%dualquat_print()
! v = (/ 10.D0, 0.D0, 0.D0 /)
! v = dqp%dualquat_Lp( v )
! write (*,*) v



! ! test for the 230 space groups to ensure that all symmetry codes are correct
! ! both for point groups and k-vector sampling 

! open(unit=dataunit, file='sym-new.txt', status='unknown',form='formatted')

! do isg=1,230 
!   call SG%setSpaceGroupNumber(isg)
!   call SG%setSpaceGrouptrigonal(.FALSE.)
!   if ((isg.ge.143).and.(isg.le.167)) call SG%setSpaceGrouptrigonal(.TRUE.)

!   j=0
!   do i=1,32
!     if (SGPG(i).le.isg) j=i
!   end do
!   isym = j
!   SamplingType = PGSamplingType(isym)

! ! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
!   if ((SamplingType.eq.-1).or.(isym.eq.14).or.(isym.eq.26)) then
!     SamplingType = SG%getHexvsRho(isym)
!   end if
!   write (dataunit,"(3I6)") isg, isym, SamplingType
! end do 

! close(unit=dataunit, status='keep')


! stop 

! do i=6,11 
!   if (i.ne.7) then 
!     QCSG = QCspacegroup_T( nD = 3, QCtype = 'Ico')
!     call QCSG%setSGnum(i)
! ! call QCcell_icosahedral%setMetricParametersQC()
!     call QCSG%GenerateQCSymmetry(dopg=.FALSE.) 
!     write (*,*) ' order of icosahedral space group ', i, ' equals ',QCSG%getnsym()
!   end if 
! end do

end program