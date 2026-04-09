! ###################################################################
! Copyright (c) 2016-2026, Marc De Graef Research Group/Carnegie Mellon University
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

module MODKRTest
  !! author: MDG 
  !! version: 1.0 
  !! date: 04/07/26
  !!
  !! test of the KR modules for remapping of a homochoric grid

use mod_kinds
use mod_global

contains 

subroutine MODKRExecuteTest(res) &
           bind(c, name='MODKRExecuteTest')    ! this routine is callable from a C/C++ program
!DEC$ ATTRIBUTES DLLEXPORT :: MODKRExecuteTest

use mod_KRsupport
use mod_KRcyclic
use mod_KRdihedral
use mod_KRtetrahedral
use mod_KRoctahedral 
use mod_KRicosahedral 

use,INTRINSIC :: ISO_C_BINDING
use, intrinsic :: iso_fortran_env, only: real64

IMPLICIT NONE

integer(C_INT32_T),INTENT(OUT)  :: res
real(kind=real64),parameter     :: epsilon = 1.0E-010_real64
integer(kind=irg)               :: i 
real(kind=real64)               :: test_in(3,6), diff
real(kind=real64)               :: h_in(3), h_out(3)
real(kind=real64)               :: C2_ref(3,6), C3_ref(3,6), C4_ref(3,6), C6_ref(3,6), &
                                   D2_ref(3,6), D3_ref(3,6), D4_ref(3,6), D6_ref(3,6), &
                                   T_ref(3,6), O_ref(3,6), I_ref(3,6)

!===================================================
! set the reference values (verified against Zach Varley's results)
test_in = reshape( &
    (/  0.000000000000000e+00_real64, 0.000000000000000e+00_real64, 0.000000000000000e+00_real64,&
        1.000000000000000e-01_real64, 2.000000000000000e-02_real64, 3.000000000000000e-02_real64,&
        1.500000000000000e-01_real64, 1.000000000000000e-01_real64, 5.000000000000000e-02_real64,&
        2.000000000000000e-01_real64, 5.000000000000000e-02_real64, 1.100000000000000e-01_real64,&
        2.500000000000000e-01_real64, 1.200000000000000e-01_real64, 1.000000000000000e-01_real64,&
        3.000000000000000e-01_real64, 1.500000000000000e-01_real64, 1.500000000000000e-01_real64 /), (/ 3,6 /) )

C2_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        9.572284703346512E-002_real64, 1.914456940669302E-002_real64, 1.547733986242521E-002_real64,&
        1.438171185505684E-001_real64, 9.587807903371225E-002_real64, 2.575094681938346E-002_real64,&
        1.885601240873353E-001_real64, 4.714003102183382E-002_real64, 5.804390515282085E-002_real64,&
        2.379749458288986E-001_real64, 1.142279739978713E-001_real64, 5.193666198720319E-002_real64,&
        2.839814125788652E-001_real64, 1.419907062894326E-001_real64, 7.854984979198544E-002_real64 /), (/ 3,6 /) )

C3_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        9.501675723171719E-002_real64, 1.900335144634344E-002_real64, 1.037099514075410E-002_real64,&
        1.427994687442143E-001_real64, 9.519964582947621E-002_real64, 1.725000722376921E-002_real64,&
        1.866073495920480E-001_real64, 4.665183739801201E-002_real64, 3.904847279075144E-002_real64,&
        2.359679583507703E-001_real64, 1.132646200083698E-001_real64, 3.484156766835086E-002_real64,&
        2.812748331509140E-001_real64, 1.406374165754570E-001_real64, 5.277139682589860E-002_real64 /), (/ 3,6 /) )

C4_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        9.477511596686281E-002_real64, 1.895502319337256E-002_real64, 7.791823574201853E-003_real64,&
        1.424513735444466E-001_real64, 9.496758236296440E-002_real64, 1.295877187115305E-002_real64,&
        1.859353177767037E-001_real64, 4.648382944417594E-002_real64, 2.937784909696136E-002_real64,&
        2.352798808136089E-001_real64, 1.129343427905323E-001_real64, 2.618718720565342E-002_real64,&
        2.803449991581404E-001_real64, 1.401724995790702E-001_real64, 3.968331165111177E-002_real64 /), (/ 3,6 /) )

C6_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        9.460420194595248E-002_real64, 1.892084038919049E-002_real64, 5.200958962064965E-003_real64,&
        1.422052145241237E-001_real64, 9.480347634941579E-002_real64, 8.649218353365556E-003_real64,&
        1.854588639669962E-001_real64, 4.636471599174904E-002_real64, 1.962859022119232E-002_real64,&
        2.347928315159212E-001_real64, 1.127005591276422E-001_real64, 1.748459803557552E-002_real64,&
        2.796862539033402E-001_real64, 1.398431269516701E-001_real64, 2.650512472789829E-002_real64 /), (/ 3,6 /) )

D2_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        5.819585839893087E-002_real64, 1.190379341523974E-002_real64, 2.081898711435249E-002_real64,&
        9.513918326546639E-002_real64, 6.431101244619619E-002_real64, 3.423171599031281E-002_real64,&
        1.213523131980589E-001_real64, 3.101486481918652E-002_real64, 7.629250874661539E-002_real64,&
        1.536222154691511E-001_real64, 7.512885854266711E-002_real64, 6.893845927131441E-002_real64,&
        1.875710276397318E-001_real64, 9.551340046790036E-002_real64, 1.033475884133597E-001_real64 /), (/ 3,6 /) )

D3_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        5.894220235498526E-002_real64, 1.199158106206740E-002_real64, 1.481861797601858E-002_real64,&
        9.288506947592473E-002_real64, 6.164500303625569E-002_real64, 2.440910640515023E-002_real64,&
        1.259963704733679E-001_real64, 3.200070924223777E-002_real64, 5.426306289857995E-002_real64,&
        1.566288650280200E-001_real64, 7.563783170441760E-002_real64, 4.885552693590505E-002_real64,&
        1.929490602109810E-001_real64, 9.694813705224181E-002_real64, 7.320786471542344E-002_real64 /), (/ 3,6 /) )

D4_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        5.926060128091191E-002_real64, 1.198296676404450E-002_real64, 1.135224543723497E-002_real64,&
        8.872788637935532E-002_real64, 5.888032700790971E-002_real64, 1.892039812383346E-002_real64,&
        1.279390722278293E-001_real64, 3.227712742574967E-002_real64, 4.156029475471136E-002_real64,&
        1.541003098483238E-001_real64, 7.373020136597484E-002_real64, 3.762666702789792E-002_real64,&
        1.890883805757348E-001_real64, 9.418647938250653E-002_real64, 5.648628251406934E-002_real64 /), (/ 3,6 /) )

D6_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        5.950768821692449E-002_real64, 1.194499196947368E-002_real64, 7.679671537911328E-003_real64,&
        8.797396828537178E-002_real64, 5.871307367796556E-002_real64, 1.283515488571061E-002_real64,&
        1.294301908823706E-001_real64, 3.239196582354556E-002_real64, 2.811121466003399E-002_real64,&
        1.500787308694375E-001_real64, 7.193122150062228E-002_real64, 2.566697574899849E-002_real64,&
        1.850443121369316E-001_real64, 9.241324150448472E-002_real64, 3.850725010117557E-002_real64 /), (/ 3,6 /) )

T_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        4.796231773301094E-002_real64, 7.590884211292534E-003_real64, 1.121230294654296E-002_real64,&
        6.446801955814822E-002_real64, 3.984031703634081E-002_real64, 1.826014678333339E-002_real64,&
        8.997072579817145E-002_real64, 1.816420440292857E-002_real64, 4.298143637148893E-002_real64,&
        1.091605885335856E-001_real64, 4.596800031159518E-002_real64, 3.734314292731306E-002_real64,&
        1.281033470975807E-001_real64, 5.655337415234034E-002_real64, 5.698184268163776E-002_real64 /), (/ 3,6 /) )

O_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        3.067028207841027E-002_real64, 7.953940577141768E-003_real64, 1.163382647065487E-002_real64,&
        5.247633948823253E-002_real64, 3.947260269110737E-002_real64, 2.051714832565471E-002_real64,&
        6.629697143322127E-002_real64, 2.026632782851863E-002_real64, 4.272725395501790E-002_real64,&
        8.336768354570702E-002_real64, 4.801174446848955E-002_real64, 4.050994627943932E-002_real64,&
        1.031434970228504E-001_real64, 6.118991373317085E-002_real64, 6.118991373317084E-002_real64 /), (/ 3,6 /) )

I_ref = reshape( &
     (/ 0.000000000000000E+000_real64, 0.000000000000000E+000_real64, 0.000000000000000E+000_real64,&
        2.428971010989204E-002_real64, 5.446101857316433E-003_real64, 6.995106347357991E-003_real64,&
        4.164554832139043E-002_real64, 2.810169094978036E-002_real64, 1.317827149345376E-002_real64,&
        4.796379509123005E-002_real64, 1.326505990022206E-002_real64, 2.633623077604095E-002_real64,&
        6.438780198565232E-002_real64, 3.297468312521545E-002_real64, 2.441115868175631E-002_real64,&
        7.702448383020780E-002_real64, 4.063095135596443E-002_real64, 3.699485246561798E-002_real64 /), (/ 3,6 /) )


res = 0
!===================================================

!===================================================
! in this test, we take six arbitrary homochoric points and apply
! all 11 Knothe--Rosenblatt (KR) rearrangements to them (mod_KRxxx modules)
!===================================================

!===================================================

! C2
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRcyclic(h_in, h_out, 2)
  diff = diff + maxval(abs(h_out(1:3)-C2_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 1
  write (*,"('Cyclic group 2 failed = ',D18.10)") diff
  return
end if

! C3
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRcyclic(h_in, h_out, 3)
  diff = diff + maxval(abs(h_out(1:3)-C3_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 2
  write (*,"('Cyclic group 3 failed = ',D18.10)") diff
  return
end if

! C4
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRcyclic(h_in, h_out, 4)
  diff = diff + maxval(abs(h_out(1:3)-C4_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 3
  write (*,"('Cyclic group 4 failed = ',D18.10)") diff
  return
end if

! C6
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRcyclic(h_in, h_out, 6)
  diff = diff + maxval(abs(h_out(1:3)-C6_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 4
  write (*,"('Cyclic group 6 failed = ',D18.10)") diff
  return
end if

! D2
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRdihedral(h_in, h_out, 2)
  diff = diff + maxval(abs(h_out(1:3)-D2_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 5
  write (*,"('Dihedral group 2 failed = ',D18.10)") diff
  return
end if

! D3
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRdihedral(h_in, h_out, 3)
  diff = diff + maxval(abs(h_out(1:3)-D3_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 6
  write (*,"('Dihedral group 3 failed = ',D18.10)") diff
  return
end if

! D4
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRdihedral(h_in, h_out, 4)
  diff = diff + maxval(abs(h_out(1:3)-D4_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 7
  write (*,"('Dihedral group 4 failed = ',D18.10)") diff
  return
end if

! D6
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRdihedral(h_in, h_out, 6)
  diff = diff + maxval(abs(h_out(1:3)-D6_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 8
  write (*,"('Dihedral group 6 failed = ',D18.10)") diff
  return
end if

! T
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRtetrahedral(h_in, h_out)
  diff = diff + maxval(abs(h_out(1:3)-T_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 9
  write (*,"('Tetrahedral group T failed = ',D18.10)") diff
  return
end if

! O
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRoctahedral(h_in, h_out)
  diff = diff + maxval(abs(h_out(1:3)-O_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 9
  write (*,"('Octahedral group O failed = ',D18.10)") diff
  return
end if

! I
diff = 0.0_real64
do i=1,6
  h_in = test_in(1:3,i)
  call KRicosahedral(h_in, h_out)
  diff = diff + maxval(abs(h_out(1:3)-I_ref(1:3,i)))
end do 
if (diff.gt.epsilon) then
  res = 9
  write (*,"('Icosahedral group I failed = ',D18.10)") diff
  return
end if

end subroutine MODKRExecuteTest

end module MODKRTest
