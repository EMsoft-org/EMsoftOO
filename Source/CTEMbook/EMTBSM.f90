! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef/Carnegie Mellon University
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

!--------------------------------------------------------------------------
! EMsoft:EMTBSM.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMTBSM
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Two-beam scattering matrix program (CTEM book)
!
!> @detail This is a very simplistic implementation to illustrate how the scattering matrix
!> approach can be implemented; this approach is also used in more sophisticated programs.
!
!> @date  4/28/01 MDG 1.0 original
!> @date  5/27/01 MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!> @date  2/15/24 MDG 4.0 rewrite for EMsoftOO
!--------------------------------------------------------------------------
program EMTBSM

use mod_kinds
use mod_global
use mod_EMsoft

IMPLICIT NONE

character(fnlen) :: progname = 'EMTBSM.f90'
character(fnlen) :: progdesc = 'CTEMbook: Two-beam computations: comparing scattering matrix and analytical'

type(EMsoft_T)   :: EMsoft

! initialize 
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 920 /) )

! do the computation
call CalcTBSM(EMsoft, progdesc)

end program

!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcTBSM
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Actual two-beam scattering matrix computation; this does not use the
!>        complex data type, but instead explicitly uses real and imaginary parts,
!>        just for fun...
!
!> @date   4/28/01  MDG 1.0 original
!> @date   5/27/01  MDG 2.0 f90
!> @date  4/18/13 MDG 3.0 rewrite
!> @date  2/16/24 MDG 4.0 complete rewrite for EMsoftOO
!--------------------------------------------------------------------------
subroutine CalcTBSM(EMsoft, progdesc)

use mod_kinds
use mod_global
use mod_EMsoft
use mod_io
use mod_crystallography
use mod_symmetry
use mod_diffraction
use mod_postscript
use mod_HDFsupport 
use mod_memory
use mod_timing
use mod_TB
use mod_misc

type(EMsoft_T),INTENT(INOUT)      :: EMsoft
character(fnlen),INTENT(IN)       :: progdesc

type(Cell_T)                      :: Cell 
type(Spacegroup_T)                :: SPG 
type(Diffraction_T)               :: Diff
type(IO_T)                        :: Message
type(Memory_T)                    :: mem
type(Timing_T)                    :: timer
type(Postscript_T)                :: PS
type(gnode)                       :: rlp

integer(kind=irg),parameter       :: ns=512, nt=512
integer(kind=irg)                 :: ind(3),io_int(1), g(3) 
real(kind=sgl)                    :: bg,xg,xgp,xgpz,wmax,tmax, &
                                     Ar(2,2), Ai(2,2), sg, dz, dsg, p(2),q(2), &
                                     BF(512,512,2),DF(512,512,2),It,Is,io_real(1)
real(kind=dbl),allocatable        :: phir(:,:),phii(:,:),SMr(:,:,:),SMi(:,:,:)
integer(kind=irg),allocatable     :: imaint(:,:)
character(fnlen)                  :: xtalname

call openFortranHDFInterface()

! first get the crystal data and microscope voltage
call Message%ReadValue(' Enter xtal file name : ', xtalname )
call cell%getCrystalData(xtalname, SPG, EMsoft)
call cell%calcPositions(SPG, 'v')

call Diff%GetVoltage(cell)
preg = 2.0 * sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18
call Diff%setrlpmethod('WK')

! get the reciprocal lattice vector
call Message%printMessage(' Enter diffraction vector :')
call GetIndex(SPG%getSpaceGrouphexset(), g, 'r')

! some more parameters
call Message%ReadValue(' Enter maximum value of w: ', io_real, 1)
wmax = io_real(1)
call Message%ReadValue(' Maximum foil thickness : ', io_real, 1)
tmax = io_real(1)

! normal aborption factor
call Diff%CalcUcg(cell, (/ 0, 0, 0 /))
rlp = Diff%getrlp()
xgpz= rlp%xgp

! extinction distance and anomalous absorption length
call Diff%CalcUcg(cell, g)
rlp = Diff%getrlp()
xgp = rlp%xgp
xg  = rlp%xg
imanum = 0
dz = tmax/float(nt-1)

! compute the scattering matrix SM for each orientation
dsg = 2.0*wmax/float(ns-1)
call Message%printMessage('Starting scattering matrix computation')
io_int(1) = 8*ns
call Message%WriteValue(' Number of computations per * = ', io_int, 1, "(I4)")

! allocate arrays
mem = Memory_T()
call mem%alloc(SMr, (/ 512,2,2 /), 'SMr')
call mem%alloc(SMi, (/ 512,2,2 /), 'SMi')
call mem%alloc(phir, (/ 512,2 /), 'phir')
call mem%alloc(phii, (/ 512,2 /), 'phii')
 
timer = Timing_T()
call timer%Time_tick()

 do i=1,ns
  sg = (-wmax+dsg*float(i))/xg
  call TBCalcSM(Ar,Ai,sg,dz,xg,xgp,xgpz,bg)
  SMr(i,1:2,1:2)=Ar(1:2,1:2)
  SMi(i,1:2,1:2)=Ai(1:2,1:2)
 end do

! compute intensities for the first row
  do i=1,ns
   BF(i,1,1)=SMr(i,1,1)**2+SMi(i,1,1)**2
   DF(i,1,1)=SMr(i,2,1)**2+SMi(i,2,1)**2

! the first column of SM is the actual initial wavefunction at z=dz
   phir(i,1)=SMr(i,1,1)
   phir(i,2)=SMr(i,2,1)
   phii(i,1)=SMi(i,1,1)
   phii(i,2)=SMi(i,2,1)
  end do

! loop to maximum thickness
  do l=2,nt
   do k=1,ns
    do i=1,2
     p(i)=0.D0
     q(i)=0.D0
     do j=1,2
      p(i)=p(i)+SMr(k,i,j)*phir(k,j)-SMi(k,i,j)*phii(k,j)
      q(i)=q(i)+SMi(k,i,j)*phir(k,j)+SMr(k,i,j)*phii(k,j)
     end do
    end do
    phir(k,1:2)=p(1:2)
    phii(k,1:2)=q(1:2)
   end do
   if (mod(l,8).eq.0) then
     call Message%printMessage('*',"(1A,$)")
   end if
   do j=1,ns
    BF(j,l,1) = phir(j,1)**2+phii(j,1)**2
    DF(j,l,1) = phir(j,2)**2+phii(j,2)**2
   end do
  end do
  call timer%Time_tock()
  call Message%printMessage('done')
  io_real(1) = timer%getInterval()
  call Message%WriteValue('  Total computation time [s] ', io_real, 1, "(F10.5)") 
  call Message%printMessage(' Starting direct analytical solution')

! next, redo the computation, but this time use TBCalcInten
  call timer%Time_tick()
  do i=1,ns
   sg = (-wmax+dsg*float(i))/xg
   do j=1,nt
    t = float(j)*dz
    call TBCalcInten(It,Is,sg,t,xg,xgp,xgpz,bg)
    BF(i,j,2) = It
    DF(i,j,2) = Is
   end do
   if (mod(i,8).eq.0) then
     call Message%printMessage('*',"(1A,$)")
   end if
  end do
  call timer%Time_tock()
  call Message%printMessage('done')
  io_real(1) = timer%getInterval()
  call Message%WriteValue('  Total computation time [s] ', io_real, 1, "(F10.5)") 

! compare the two computations
  bfdiff = 0.0
  dfdiff = 0.0
  do i=1,ns
   do j=1,nt
    bfdiff = bfdiff + (BF(i,j,1)-BF(i,j,2))**2
    dfdiff = dfdiff + (DF(i,j,1)-DF(i,j,2))**2
   end do
  end do
  bfdiff = bfdiff/float(ns)/float(nt)
  dfdiff = dfdiff/float(ns)/float(nt)
  io_real(1) = bfdiff
  call Message%WriteValue(' Average difference in BF ', io_real, 1, "(F10.4)")
  io_real(1) = dfdiff
  call Message%WriteValue(' Average difference in DF ', io_real, 1, "(F10.4)") 
  call Message%printMessage(' Images computed, preparing for PS output')

! open PostScript file
  call PS%openFile(progdesc, EMsoft)
  call PS%setpspage(0)
  call PS%newpage(.TRUE.,'Two Beam Rocking Curves')
  call PS%setfont(PSfonts(2),0.16)
  call PS%text(2.5,8.6,'Scattering Matrix Solution')
  call PS%text(2.5,2.1,'Direct Analytical Solution')
  call PS%setfont(PSfonts(2),0.10)
  call PS%text(0.5,1.8,'Input file')
  call PS%text(2.6,1.8,trim(cell%getFileName()) )
  call PS%text(0.5,1.6,'Active reflection')
  call PS%textint(2.5,1.6,' ',g(1))
  call PS%textint(2.6,1.6,' ',g(2))
  call PS%textint(2.7,1.6,' ',g(3))
  call PS%text(0.5,1.4,'Extinction distance  [nm]')
  call PS%textvar(2.5,1.4,' ',xg)
  call PS%text(0.5,1.2,'Anomalous absorption length [nm]')
  call PS%textvar(2.5,1.2,' ',xgp)
  call PS%text(0.5,1.0,'Normal absorption length [nm]')
  call PS%textvar(2.5,1.0,' ',xgpz)
  call PS%text(0.5,0.8,'Accelerating Voltage [V]')
  call PS%textvar(2.5,0.8,' ', sngl(Diff%getV()) )
  call PS%text(0.5,0.6,'Maximum w')
  call PS%textvar(2.5,0.6,' ',wmax)
  call PS%text(0.5,0.4,'Maximum thickness [nm]')
  call PS%textvar(2.5,0.4,' ',tmax)

! Bright Field image
  call mem%alloc(imaint, (/ ns, nt /), 'imaint')

  imanum = 0
  imaint(1:ns,1:nt)=int(255.0*BF(1:ns,1:nt,1))
  x0=0.2
  y0=5.5
  npx=ns
  npy=nt
  scl=3.0
  call PS%DumpImageDistort(imaint,imanum,x0,y0,npx,npy,scl,scl)

! Dark Field image
  imaint(1:ns,1:nt)=int(255.0*DF(1:ns,1:nt,1))
  x0=3.4
  y0=5.5
  npx=ns
  npy=nt
  scl=3.0
  call PS%DumpImageDistort(imaint,imanum,x0,y0,npx,npy,scl,scl)

! Bright Field image
  imaint(1:ns,1:nt)=int(255.0*BF(1:ns,1:nt,2))
  x0=0.2
  y0=2.3
  npx=ns
  npy=nt
  scl=3.0
  call PS%DumpImageDistort(imaint,imanum,x0,y0,npx,npy,scl,scl)

! Dark Field image
  imaint(1:ns,1:nt)=int(255.0*DF(1:ns,1:nt,2))
  x0=3.4
  y0=2.3
  npx=ns
  npy=nt
  scl=3.0
  call PS%DumpImageDistort(imaint,imanum,x0,y0,npx,npy,scl,scl)

! close Postscript file
  call PS%closeFile()

! deallocate arrays
call mem%dealloc(SMr, 'SMr')
call mem%dealloc(SMi, 'SMi')
call mem%dealloc(phir, 'phir')
call mem%dealloc(phii, 'phii')
call mem%dealloc(imaint, 'imaint')

call closeFortranHDFInterface()

end subroutine CalcTBSM
       

