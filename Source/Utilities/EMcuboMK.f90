! ###################################################################
! Copyright (c) 2015-2020, Marc De Graef Research Group/Carnegie Mellon University
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

program EMcuboMK
 !! author: MDG 
 !! version: 1.0 
 !! date: 01/25/20
 !!
 !! Generate MacKenzie histogram for a given rotational symmetry based on cubochoric sampling
 !!
 !! This is a simple utility program to generate a histogram (.csv) of a MacKenzie misorientation plot 
 !! for a given rotational symmetry group, based on a concentric sampling of cubochoric space.

use mod_kinds
use mod_global 
use mod_EMsoft
use mod_io
use mod_rotations
use mod_so3
use mod_symmetry 

IMPLICIT NONE 

character(fnlen)            :: progname = 'EMcuboMK.f90'
character(fnlen)            :: progdesc = 'Generate MacKenzie histogram for rotational symmetry based on cubochoric sampling'

type(EMsoft_T)                 :: EMsoft 
type(SpaceGroup_T)            :: SG 
type(so3_T)                    :: SO
type(IO_T)                    :: Message 
type(r_T)                    :: rod 
type(c_T)                    :: cu 

integer(kind=irg)            :: pgnum, FZorder, FZtype, nsteps, n, i, j, k, io_int(3), ntot, pgrotOrder
real(kind=dbl)                :: sedge, delta, x, y, z, tot, rho, ho2(3), ho1(3), vol, xx(4)
real(kind=dbl),allocatable    :: misor(:), histogram(:), e(:), mk(:)
logical                     :: b

! print the header information and handle command line arguments
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 906 /) )

! first ask for the rotational point group number
call SG%ListPointGroups()
call Message%ReadValue('Enter the desired point group number: ',io_int)
pgnum = io_int(1)
pgrotOrder = PGTHDorder(PGrot(pgnum))

! determine the Laue symmetry group and print it 
call Message%printMessage(' Corresponding pure rotational group is '//trim(PGTHD(PGrot(pgnum))), frm="(/A)")

! get the fundamental zone parameters
SO = so3_T(pgnum)
call SO%getFZtypeandorder(FZtype, FZorder)
io_int(1) = pgnum
io_int(2) = FZtype
io_int(3) = FZorder 
call Message%WriteValue(' pgnum, FZ type and order :',io_int,3)

! ask the user for the number of steps along the cubochoric semi edge length
call Message%printMessage(' ')
call Message%ReadValue('Enter the desired number of steps along the cubochoric semi-edge length: ',io_int)
nsteps = io_int(1)

! allocate output arrays 
allocate(misor(0:nsteps), histogram(0:nsteps), e(0:nsteps), mk(0:nsteps))
histogram = 0.D0

! set some auxiliary parameters
sedge = 0.5D0 * LPs%ap
delta = sedge / dble(nsteps)

! start the sampling, working in concentric cubes from the center outward
do n = 0, nsteps 
    if (n.eq.0) then 
        misor(0) = 0.D0
        e(0) = 0.D0
        histogram(0) = 0.D0 
    else ! loop over concentric cubes 
        e(n) = dble(n) * delta   ! semi-edge length of sub-cube
        ! get the angle first 
        if (n.lt.nsteps) then 
            cu = c_T( cdinp = (/ e(n), 0.D0, 0.D0 /) )
            rod = cu%cr()
            xx = rod%r_copyd()
            misor(n) = 2.D0 * atan(xx(4))
        else
            misor(n) = cPi
        end if
        ! cover the +/- x faces with a complete grid
        do j = -n, n
            y = dble(j) * delta
            do k = -n, n
                z = dble(k) * delta
                cu = c_T( cdinp = (/ e(n), y, z /) )
                rod = cu%cr()
                b = SO%IsinsideFZ(rod)
                if (b) histogram(n) = histogram(n) + 1.D0
                cu = c_T( cdinp = (/ -e(n), y, z/) )
                rod = cu%cr()
                b = SO%IsinsideFZ(rod)
                if (b) histogram(n) = histogram(n) + 1.D0
            end do 
        end do 
        ! then cover the +/- y faces with a grid, omitting the already covered points 
        do i = -n+1, n-1
            x = dble(i) * delta
            do k = -n, n
                z = dble(k) * delta
                cu = c_T( cdinp = (/ x, e(n), z /) )
                rod = cu%cr()
                b = SO%IsinsideFZ(rod)
                if (b) histogram(n) = histogram(n) + 1.D0
                cu = c_T( cdinp = (/ x, -e(n), z /) )
                rod = cu%cr()
                b = SO%IsinsideFZ(rod)
                if (b) histogram(n) = histogram(n) + 1.D0
            end do 
        end do 
        ! finally cover the +/- z faces with a grid, omitting the already covered points 
        do i = -n+1, n-1
            x = dble(i) * delta
            do j = -n+1, n-1
                y = dble(j) * delta
                cu = c_T( cdinp = (/ x, y, e(n) /) )
                rod = cu%cr()
                b = SO%IsinsideFZ(rod)
                if (b) histogram(n) = histogram(n) + 1.D0
                cu = c_T( cdinp = (/  x, y, -e(n) /) )
                rod = cu%cr()
                b = SO%IsinsideFZ(rod)
                if (b) histogram(n) = histogram(n) + 1.D0
            end do 
        end do 
        histogram(n) = histogram(n) * sin(misor(n)*0.5D0)**2 / dble( 2.D0 + 24.D0 * n**2 )
    end if
end do

call SO%getMacKenzieDistribution(nsteps, misor, mk)

! normalize 
histogram =  pgrotOrder * (4.D0 * cPi) * (cPi / 180.D0) * histogram / ( 2.D0 * cPi**2 )
misor = misor / dtor 

! and print the results to a text file
open(unit=20,file='EMcuboMK.csv',status='unknown',form='formatted')
write (20,"(A)") 'angle, sampled, theoretical'
write (20,"(I5,',',I5,',',I5)") nsteps, nsteps, nsteps
do n=0,nsteps
  write (20,"(F12.8,',',F12.8,',',F12.8)") misor(n), histogram(n), mk(n)
end do 
close(unit=20,status='keep')
call Message%printMessage('Results stored in EMcuboMK.csv file...', frm="(/A)")

end program EMcuboMK
