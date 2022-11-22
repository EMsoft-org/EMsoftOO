! ###################################################################
! Copyright (c) 2013-2022, Marc De Graef Research Group/Carnegie Mellon University
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

program EMHSGtest
  !! author: MDG 
  !! version: 1.0 
  !! date: 11/21/22
  !!
  !! output all Hall space group symbol Jones faithful strings to make sure 
  !! the Hall space group implementation is correct

use mod_EMsoft
use mod_global
use mod_io
use mod_symmetry
use mod_crystallography
use mod_HallSG

IMPLICIT NONE

type(IO_T)                  :: Message
type(SpaceGroup_T)          :: SG
      
character(3)                :: pos, fnum
integer(kind=irg)           :: p(4),ii,jj,i, io_int(1), Hnum, Hcnt, HSGnum, hi, hj
real(kind=sgl)              :: ppp, io_real(1)
real(kind=dbl),allocatable  :: sgdata(:,:,:)
character(fnlen)            :: mess
character(16)               :: HS

pos = 'xyz'
Hcnt = 1
do hi=1,230
  Hnum = get_Hall_nentries(hi)
  do hj=1,Hnum 
   HS = get_HallString(Hcnt)
write (*,*) 'generating cif file for '//trim(HS), hi, hj, Hnum, Hcnt, SG%getSpaceGroupMATnum()
  ! if (scan(HS,'R').eq.0) then 
   write (fnum,"(I3.3)") Hcnt
   SG = SpaceGroup_T( SGnumber = hi, useHall=.TRUE., HallSGnumber=Hcnt, setting=1)
   open (unit=10,file='HallGroups/H_'//fnum//'.cif',status='unknown',form='formatted')

   write(10,"(A)") '# keep this part the same for every file'
   write(10,"(A)") 'data_1'
   write(10,"(A)") '_cell_length_a    4'
   write(10,"(A)") '_cell_length_b    4'
   write(10,"(A)") '_cell_length_c    4'
   write(10,"(A)") '_cell_angle_alpha 90'
   write(10,"(A)") '_cell_angle_beta  90'
   write(10,"(A)") '_cell_angle_gamma 90'
   write(10,"(A)") 'loop_'
   write(10,"(A)") '_atom_site_label'
   write(10,"(A)") '_atom_site_fract_x'
   write(10,"(A)") '_atom_site_fract_y'
   write(10,"(A)") '_atom_site_fract_z'
   write(10,"(A)") 'H 0 0 0'
   write(10,"(A)") '# make this part differently for each'
   write(10,"(A)") '_space_group_name_Hall '''//trim(HS)//''''
   write(10,"(A)") 'loop_'
   write(10,"(A)") '_space_group_symop_operation_xyz'

! get the symmetry matrices 
   sgdata = SG%getSpaceGroupDataMatrices()

! loop over all symmetry matrices
   do i=1, SG%getSpaceGroupMATnum()
      write (10,"(A1)",advance="no") ''''
! loop over all rows
      do ii=1,3

    ! loop over colums (get the numbers)
         p(1:3)=int(sgdata(i,ii,1:3)) 
         ppp = sngl(sgdata(i,ii,4))

    ! print each entry 
    ! first the part containing x, y, and/or z
         do jj=1,3
          if (p(jj).ne.0) then
           mess(1:1)='+'
           if (p(jj).eq.-1) mess(1:1)='-'
           mess(2:2) = pos(jj:jj)
           write(10,"(A2)",advance="no") mess
          end if
         end do 

    ! if there is a translation component, print it
       if (ppp.ne.0.0) then
        mess(1:1)='+'
        if (ppp.lt.0.0) mess(1:1)='-'
        write(10,"(A1)",advance="no") mess
        write(10,"(f5.3)",advance="no") abs(ppp)
       end if

    ! print a comma, or close the brackets and do a newline
       if (ii.ne.3) then 
         write(10,"(A1)",advance="no") ','
        else
         write(10,"(A1)") ''''
       end if
      end do    
   end do
   close(unit=10,status='keep')
  ! end if 
   Hcnt = Hcnt + 1
  end do
end do

end program EMHSGtest
