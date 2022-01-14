! ###################################################################
! Copyright (c) 2014-2022, Marc De Graef Research Group/Carnegie Mellon University
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
! EMsoft:mod_defect.f90
!--------------------------------------------------------------------------
!
! MODULE: mod_defect
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides a routine to compute the displacement vector for an array of defects.
!
!> @date 04/29/11 MDG 1.0 original
!> @date 06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!> @date 11/13/13 MDG 2.1 fixed error with coordinate transformations (after very long bug search!)
!> @date 11/17/15 MDG 3.0 start of complete rewrite; this mod will now include a routine to read all defect info from a json file
!> @date 11/23/15 MDG 3.1 inserted all defect mods into this file instead of separate files
!> @date 12/08/15 MDG 3.2 added artificial distortion to inclusion field to mimic ellipsoidal shape (needs Eshelby for correct field)
!> @date 12/11/15 MDG 3.3 gave up on previous item and implemented full isotropic Eshelby ellipsoidal inclusion
!--------------------------------------------------------------------------

module mod_defect

use mod_kinds
use mod_global
use mod_quaternions

IMPLICIT NONE

type, public :: foiltype
  real(kind=dbl)                :: F(3), q(3),Fn(3),qn(3),brx,bry,brxy,cpx,cpy, &
                                   alP,alS,alR,beP,elmo(6,6),z0,zb,B(3),Bn(3),Bm(3)
  type(Quaternion_T)            :: a_fc, a_fm, a_mi, a_ic, a_mc, a_fi
  integer(kind=irg)             :: npix,npiy
  real(kind=sgl),allocatable    :: sg(:,:)
end type foiltype

type, public :: dislocationtype
  real(kind=dbl)                :: burg(3),burgd(3),u(3),un(3),g(3),gn(3),id,jd, zfrac, zu
  real(kind=dbl)                :: top(3), bottom(3)
  type(Quaternion_T)            :: a_dc, a_id, a_di, a_df
  complex(kind=dbl)             :: dismat(3,3),pa(3)
end type dislocationtype

type, public :: inclusiontype
        real(kind=sgl)          ::  xpos,ypos,zpos,radius,C
end type inclusiontype

type, public :: Einclusiontype
        real(kind=dbl)                  :: xyz(3), a123(3)
        real(kind=dbl)                  :: a1, a2, a3, principalaxes(3,3), permut(3,3), rotell(3,3), epsstar(3,3)
        real(kind=dbl)                  :: nu, omnu, pre, V, a12, a22, a32, asq(3), eta, ss1, svec(3), qs1, qvec1(3), &
                                           qvec2(3), Deltaij(3,3), kEl, preI1, preI3, s3, c1, c2, mith, math, thpre, &
                                           IIinside(3), IIJinside(3,3), xpos, ypos, zpos, ESV(6,6), EshelbyS(3,3,3,3)
        real(kind=dbl),allocatable      :: EFLUT(:), EELUT(:)
        integer(kind=irg)               :: nLUT
end type Einclusiontype

type, public :: stackingfaulttype
  real(kind=sgl)              :: lpu(3),tpu(3),lpb(3),lpbc(3),tpb(3),plane(3),sep,id,jd, &
                                lptop(3),lpbot(3),tptop(3),tpbot(3),thetan,a_if(3,3), &
                                lpr(3),tpr(3), Rdisp(3), poisson
  real(kind=sgl), allocatable :: zpos(:,:)
end type stackingfaulttype

type, public :: voidtype
        real(kind=sgl)       ::  xpos,ypos,zpos,radius
end type voidtype

type, public :: YDtype
  real(kind=dbl)             :: burg(3), burgd(3), u(3), un(3), g(3), gn(3), id, jd, zu, bs, be, bx, beta
  real(kind=dbl)             :: alpha, ca, sa, ta, cota,  top(3), bottom(3), sig
  type(Quaternion_T)         :: a_dc, a_id, a_di
end type YDtype

type, public :: apbtype
        real(kind=sgl)       ::  xpos,ypos,zpos,radius,w,Rdisp(3)
end type apbtype

type, public :: Defect_T
  integer(kind=irg)                        :: numvoids,numdisl,numYdisl,numsf,numinc,numEinc,numapb
  character(fnlen)                         :: foilname
  integer(kind=irg)                        :: Nmat,DF_g(3),DF_npix,DF_npiy,DF_nums,DF_numinclusion,DF_numvoid
  real(kind=sgl)                           :: DF_slice,DF_L,DF_gc(3),DF_gstar(3), DF_gf(3)
  type (foiltype)                          :: foil
  real(kind=sgl),allocatable               :: DF_foilsg(:,:),DF_R(:,:)
  type (dislocationtype), allocatable      :: DL(:)
  type (inclusiontype), allocatable        :: inclusions(:)
  type (Einclusiontype), allocatable       :: Einclusions(:)
  type (stackingfaulttype), allocatable    :: SF(:)
  type (voidtype), allocatable             :: voids(:)
  type (YDtype), allocatable               :: YD(:)
  type (apbtype), allocatable              :: apbs(:)

  contains
  private
  ! basic space group generating routines and related stuff
  procedure, pass(self) :: JSONreadDefectFile_
  procedure, pass(self) :: JSONreadFoilData_
  procedure, pass(self) :: InitializeDefects_
  procedure, pass(self) :: init_foil_data_
  procedure, pass(self) :: initialize_foil_geometry_
  procedure, pass(self) :: init_dislocation_data_
  procedure, pass(self) :: init_YSH_dislocation_data_
  procedure, pass(self) :: init_stacking_fault_data_
  procedure, pass(self) :: init_void_data_
  procedure, pass(self) :: init_inclusion_data_
  procedure, pass(self) :: makestackingfault_
  procedure, pass(self) :: makestackingfaultECCI_
  procedure, pass(self) :: makedislocation_
  procedure, pass(self) :: makeYSHdislocation_
  procedure, pass(self) :: YSHDisp_
  procedure, pass(self) :: Eshelby_disp_
  procedure, pass(self) :: calcR_

  generic, public :: JSONreadDefectFile => JSONreadDefectFile_
  generic, public :: JSONreadFoilData => JSONreadFoilData_
  generic, public :: InitializeDefects => InitializeDefects_
  generic, public :: init_foil_data => init_foil_data_
  generic, public :: initialize_foil_geometry => initialize_foil_geometry_
  generic, public :: init_dislocation_data => init_dislocation_data_
  generic, public :: init_YSH_dislocation_data => init_YSH_dislocation_data_
  generic, public :: init_stacking_fault_data => init_stacking_fault_data_
  generic, public :: init_void_data => init_void_data_
  generic, public :: init_inclusion_data => init_inclusion_data_
  generic, public :: makedislocation => makedislocation_
  generic, public :: makeYSHdislocation => makeYSHdislocation_
  generic, public :: makestackingfaultECCI => makestackingfaultECCI_
  generic, public :: makestackingfault => makestackingfault_
  generic, public :: YSHDisp => YSHDisp_
  generic, public :: Eshelby_disp => Eshelby_disp_
  generic, public :: calcR => calcR_

  end type Defect_T

  ! the constructor routine for this class
  interface Defect_T
    module procedure Defect_constructor
  end interface Defect_T

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in this class
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(Defect_T) function Defect_constructor() result(Def)
!DEC$ ATTRIBUTES DLLEXPORT :: Defect_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! constructor for the Defect Class

use mod_io

IMPLICIT NONE

end function Defect_constructor


!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadDefectFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief parse json file into defect namelist structures
!
!> @param cell unit cell pointer
!> @param jsonname input file name
!> @param defects defect structure, to be filled by this routine
!> @param error_cnt total number of errors encountered by json routines
!> @param verbose [optional] print a lot of output if present and true
!
!> @date 11/20/15 MDG 1.0 new routine
!> @date 12/08/15 MDG 1.1 added Einclusion defect type
!--------------------------------------------------------------------------
recursive subroutine JSONreadDefectFile_(self, EMsoft, cell, jsonname, error_cnt,verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: JSONreadDefectFile_

use ISO_C_BINDING
use mod_io
use mod_crystallography
use mod_JSONsupport
use mod_EMsoft

use, intrinsic :: iso_fortran_env, only: wp => real64

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)                        :: self
type(Cell_T),INTENT(IN)                               :: cell
type(IO_T)                                            :: Message
type(EMsoft_T), INTENT(INOUT)                         :: EMsoft
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt
logical,INTENT(IN),OPTIONAL                           :: verbose

type(json_file)                                       :: json    !the JSON structure read from the file:
type(json_value),pointer                              :: jval, child, child2, child3, child4
character(kind=jsonCK,len=:),allocatable                  :: name
integer(kind=irg)                                     :: i, j, jj, kk, v, io_int(3), jskip, ndis
integer(kind=irg)                                     :: vart,nc, nc2, nc3, nc4, nc5
logical                                               :: found
character(fnlen)                                      :: foilfilename, str, filename, dummystr
real(wp)                                              :: v4(4), v5(5), v6(6), v9(9), io_real(6)

v = 0
if (PRESENT(verbose)) then
  if (verbose) then
    v = 1
  end if
end if
dummystr = ''
! first of all, initialize json and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

filename = EMsoft%generateFilePath('EMdatapathname',trim(jsonname))

! populate the jval json_value structure
call json_parse(trim(filename), jval)
if (json_failed().eqv..TRUE.) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! get the top level file descriptor (should be the file name) and the number of its children nc (should really be 1)
  call json_info(jval,vart,nc,name)             ! jval name = filename
  if (v.eq.1) then
    call Message%WriteValue(' Defect file name : ',name,"(' ',A)")
  end if


! loop over the children
  do i=1,nc
    call json_get_child(jval,i,child)
    call json_info(child,vart,nc2,name)         ! child name = DefectDescriptors

! loop over those children, which are the actual defect descriptors and deal with each of them separately

! the foil date must be the first entry; if it is not, then loop until we find it
    found = .FALSE.
    nc2loop: do j=1,nc2
      call json_get_child(child,j,child2)
      call json_info(child2,vart,nc3,name)
      if (.not.(name.eq.'foil')) CYCLE nc2loop
! name = foil, so read the name of the foil descriptor file
      call json_get(child2,'foilfilename',child3,found)
      jskip = j
      if (found) then
        call json_get(child3, name)
        if (v.eq.1) call Message%WriteValue(' Foil file name = ',trim(name),"(' ',A)")
        self%foilname = trim(name)
      end if
    end do nc2loop

    if (.not.found) then
     call Message%printError('JSONreadDefectFile','JSON file does not contain a foilfilename entry')
    end if

 ! here we call the foil reading routine to first fill all the foil parameters
     call self%JSONreadFoilData(Emsoft, cell, error_cnt, verbose)
! then we need to get the total number of defects in the file, so that we can allocate
! the correct array sizes in the defects structure
    ndis = 0
    nc2loop2: do j=1,nc2
      if (j.eq.jskip) CYCLE nc2loop2
      call json_get_child(child,j,child2)
      call json_info(child2,vart,nc3,name)
      if (name.eq.'voids') then
        allocate(self%voids(nc3))
        self%numvoids = nc3
      end if
      if (name.eq.'inclusions') then
        allocate(self%inclusions(nc3))
        self%numinc = nc3
      end if
      if (name.eq.'Einclusions') then
        allocate(self%Einclusions(nc3))
        self%numEinc = nc3
      end if
      if (name.eq.'Ydislocations') then
        allocate(self%YD(nc3))
        self%numYdisl = nc3
      end if
      if (name.eq.'dislocations') then
        ndis = ndis + nc3
        self%numdisl = nc3
      end if
      if (name.eq.'stackingfaults') then
        allocate(self%SF(nc3))
        ndis = ndis + 2*nc3
        self%numsf = nc3
      end if
    end do nc2loop2
    if (ndis.gt.0) allocate(self%DL(ndis))

! now loop over all entries at the child level (Except for the foil data) and
! read the individual defect parameters; note that these are nested on level 4...
    ndis = 1
    nc2loop3: do j=1,nc2
      if (j.eq.jskip) CYCLE nc2loop3
      call json_get_child(child,j,child2)
      call json_info(child2,vart,nc3,name)

! dislocations
      if (name.eq.'dislocations') then
        do jj=1,nc3
         if (v.eq.1) then
           io_int(1) = jj
           call Message%WriteValue('   dislocation #  ',io_int,1,"(I4)")
         end if
         call json_get_child(child2,jj,child3)
         call json_info(child3,vart,nc4,name)
         do kk=1,nc4
          call json_get_child(child3,kk,child4)
          call json_info(child4,vart,nc5,name)
          if (name.eq.'id') then
            str = '        x-coordinate  = '
            self%DL(ndis)%id = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'jd') then
            str = '        y-coordinate  = '
            self%DL(ndis)%jd = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'zfrac') then
            str = '        zfrac         = '
            self%DL(ndis)%zfrac = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'u') then
            str = '        u             = '
            self%DL(ndis)%u = JSONgetDoubleVector(child4,nc5,str,v)
          end if
          if (name.eq.'bv') then
            str = '        bv            = '
            self%DL(ndis)%burg = JSONgetDoubleVector(child4,nc5,str,v)
          end if
         end do
         ndis = ndis + 1
        end do
        CYCLE nc2loop3
      end if

! Ydislocations
      if (name.eq.'Ydislocations') then
        do jj=1,nc3
         if (v.eq.1) then
           io_int(1) = jj
           call Message%WriteValue('   Ydislocation #  ',io_int,1,"(I4)")
         end if
         call json_get_child(child2,jj,child3)
         call json_info(child3,vart,nc4,name)
         do kk=1,nc4
          call json_get_child(child3,kk,child4)
          call json_info(child4,vart,nc5,name)
          if (name.eq.'id') then
            str = '        x-coordinate  = '
            self%YD(jj)%id = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'jd') then
            str = '        y-coordinate  = '
            self%YD(jj)%jd = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'poisson') then
            str = '        poisson       = '
            self%YD(jj)%sig = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'u') then
            str = '        u             = '
            self%YD(jj)%u = JSONgetDoubleVector(child4,nc5,str,v)
          end if
          if (name.eq.'bv') then
            str = '        bv            = '
            self%YD(jj)%burg = JSONgetDoubleVector(child4,nc5,str,v)
          end if
         end do
        end do
        CYCLE nc2loop3
      end if

! stacking faults
      if (name.eq.'stackingfaults') then
        do jj=1,nc3
         if (v.eq.1) then
           io_int(1) = jj
           call Message%WriteValue('   Stacking Fault #  ',io_int,1,"(I4)")
         end if
         call json_get_child(child2,jj,child3)
         call json_info(child3,vart,nc4,name)
         do kk=1,nc4
          call json_get_child(child3,kk,child4)
          call json_info(child4,vart,nc5,name)
          if (name.eq.'SFi') then
            str = '        x-coordinate  = '
            self%SF(jj)%id = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'SFj') then
            str = '        y-coordinate  = '
            self%SF(jj)%jd = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'SFsep') then
            str = '        SF separation = '
            self%SF(jj)%sep = JSONgetDouble(child4,str,v)
          end if
          if (name.eq.'SFplane') then
            str = '        SF plane      = '
            self%SF(jj)%plane = JSONgetDoubleVector(child4,nc5,str,v)
          end if
          if (name.eq.'SFlpu') then
            str = '        SF lpu        = '
            self%SF(jj)%lpu = JSONgetDoubleVector(child4,nc5,str,v)
          end if
          if (name.eq.'SFtpu') then
            str = '        SF tpu        = '
            self%SF(jj)%tpu = JSONgetDoubleVector(child4,nc5,str,v)
          end if
          if (name.eq.'SFlpb') then
            str = '        SF lpb        = '
            self%SF(jj)%lpb = JSONgetDoubleVector(child4,nc5,str,v)
          end if
          if (name.eq.'SFtpb') then
            str = '        SF tpb        = '
            self%SF(jj)%tpb = JSONgetDoubleVector(child4,nc5,str,v)
          end if
         end do
        end do
        CYCLE nc2loop3
      end if

! voids
      if (name.eq.'voids') then
        do jj=1,nc3
         call json_get_child(child2,jj,child3)
         call json_info(child3,vart,nc4,name)
         v4 = JSONgetDoubleVector(child3,nc4,dummystr,0)
         self%voids(jj)%xpos = v4(1)
         self%voids(jj)%ypos = v4(2)
         self%voids(jj)%zpos = v4(3)
         self%voids(jj)%radius = v4(4)
         if (v.eq.1) then
           io_real(1:4) = v4(1:4)
           io_int(1) = jj
           call Message%WriteValue(' void   ',io_int,1,"(I4)",advance="no")
           call Message%WriteValue('',io_real,4)
         end if
        end do
        CYCLE nc2loop3
      end if

! inclusions
      if (name.eq.'inclusions') then
        do jj=1,nc3
         call json_get_child(child2,jj,child3)
         call json_info(child3,vart,nc4,name)
         v5 = JSONgetDoubleVector(child3,nc4,dummystr,0)
         self%inclusions(jj)%xpos = v5(1)
         self%inclusions(jj)%ypos = v5(2)
         self%inclusions(jj)%zpos = v5(3)
         self%inclusions(jj)%radius = v5(4)
         self%inclusions(jj)%C = v5(5)
         if (v.eq.1) then
           io_real(1:5) = v5(1:5)
           io_int(1) = jj
           call Message%WriteValue(' inclusion    ',io_int,1,"(I4)",advance="no")
           call Message%WriteValue('',io_real,5)
         end if
        end do
        CYCLE nc2loop3
      end if

! Eshelby ellipsoidal inclusions (isotropic)
       if (name.eq.'Einclusions') then
         do jj=1,nc3
          if (v.eq.1) then
            io_int(1) = jj
            call Message%WriteValue('   Einclusion   #  ',io_int,1,"(I4)")
          end if
          call json_get_child(child2,jj,child3)
          call json_info(child3,vart,nc4,name)
          do kk=1,nc4
           call json_get_child(child3,kk,child4)
           call json_info(child4,vart,nc5,name)
           if (name.eq.'xyz') then
             str = '        xyz           = '
             self%Einclusions(jj)%xyz = JSONgetDoubleVector(child4,nc5,str,v)
           end if
           if (name.eq.'a123') then
             str = '        a123          = '
             self%Einclusions(jj)%a123 = JSONgetDoubleVector(child4,nc5,str,v)
           end if
           if (name.eq.'nu') then
             str = '        nu            = '
             self%Einclusions(jj)%nu = JSONgetDouble(child4,str,v)
           end if
           if (name.eq.'epsstarvoigt') then
             str = '        eps* (Voigt)  = '
             v6 = JSONgetDoubleVector(child4,nc5,str,v)
             self%Einclusions(jj)%epsstar(1,1) = v6(1)
             self%Einclusions(jj)%epsstar(2,2) = v6(2)
             self%Einclusions(jj)%epsstar(3,3) = v6(3)
             self%Einclusions(jj)%epsstar(1,2) = v6(6)*0.5D0
             self%Einclusions(jj)%epsstar(2,1) = v6(6)*0.5D0
             self%Einclusions(jj)%epsstar(1,3) = v6(5)*0.5D0
             self%Einclusions(jj)%epsstar(3,1) = v6(5)*0.5D0
             self%Einclusions(jj)%epsstar(2,3) = v6(4)*0.5D0
             self%Einclusions(jj)%epsstar(3,2) = v6(4)*0.5D0
           end if
           if (name.eq.'principalaxes') then
             str = '        principalaxes = '
             v9 = JSONgetDoubleVector(child4,nc5,str,v)
             self%Einclusions(jj)%principalaxes(1,1:3) = v9(1:3)
             self%Einclusions(jj)%principalaxes(2,1:3) = v9(4:6)
             self%Einclusions(jj)%principalaxes(3,1:3) = v9(7:9)
           end if
          end do
         end do
         CYCLE nc2loop3
       end if

! other defct types to be added here

    end do nc2loop3
  end do
end if

call JSON_failtest(error_cnt)

end subroutine JSONreadDefectFile_

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadFoilData
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief parse json foil file into defect namelist structures
!
!> @param cell unit cell pointer
!> @param defects defect structure, to be filled by this routine
!> @param error_cnt total number of errors encountered by json routines
!> @param verbose [optional] print a lot of output if present and true
!
!> @date 11/21/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
recursive subroutine JSONreadFoilData_(self, Emsoft, cell, error_cnt, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: JSONreadFoilData_

use ISO_C_BINDING
use mod_io
use mod_crystallography
use mod_EMsoft
use mod_JSONsupport

use, intrinsic :: iso_fortran_env, only: wp => real64

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)                        :: self
type(Cell_T),INTENT(IN)                               :: cell
type(IO_T)                                            :: Message
type(EMsoft_T), INTENT(INOUT)                         :: EMsoft
integer(kind=irg),INTENT(INOUT)                       :: error_cnt
logical,INTENT(IN),OPTIONAL                           :: verbose

type(json_value),pointer                              :: jval, child, child2, child3
type(json_value), pointer                             :: tmp_json_ptr
character(kind=jsonCK,len=:),allocatable                  :: name
integer(kind=irg)                                     :: v, i, j, jj, vart, nc, nc2, nc3, io_int(3)
real(kind=wp),dimension(:),allocatable                :: vec3
real(kind=wp)                                         :: val
real(kind=sgl)                                        :: io_real(6), x
logical                                               :: found
character(4),parameter                                :: row(6) = (/ 'row1', 'row2', 'row3', 'row4', 'row5', 'row6' /)
character(fnlen)                                      :: str, filename

v = 0
if (PRESENT(verbose)) then
  if (verbose) then
    v = 1
  end if
end if

! set the default values for all entries
self%foil%elmo = 0.0                         ! elastic moduli
self%foil%F = (/ 0.0,0.0,1.0 /)              ! foil normal in direct space Bravais reference frame
self%foil%q = (/ 1.0,0.0,0.0 /)              ! reciprocal space vector along primary tilt axis towards airlock
self%foil%alP = 0.0                          ! primary tilt angle in degrees
self%foil%alS = 0.0                          ! secondary tilt angle (for double tilt holder)
self%foil%alR = 0.0                          ! secondary tilt angle (for rotation tilt holder)
self%foil%beP = 0.0                          ! angle of primary tilt axis w.r.t. image bottom edge
self%foil%z0 = 100.0                         ! foil thickness in nm

! the following are not used currently, but need to be initialized properly
self%foil%brx = 0.0                          ! parameters to describe the foil shape as a quadratic surface
self%foil%bry = 0.0
self%foil%brxy = 0.0
self%foil%cpx = 0.0                          ! center of the foil quadratic surface within [-1,1] range in pixel coordinates
self%foil%cpy = 0.0

filename = EMsoft%generateFilePath('EMdatapathname',trim(self%foilname))

! json has alrady been initialized, so we should be ok directly reading the data from the file
call json_parse(trim(filename), jval)


if (json_failed().eqv..TRUE.) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
  call json_info(jval,vart,nc,name)             ! jval name = filename

! loop over the children (only 1)
  do i=1,nc
    call json_get_child(jval,i,child)
    call json_info(child,vart,nc2,name)         ! child name = FoilDescriptor

    nc2loop: do j=1,nc2
      call json_get_child(child,j,child2)
      call json_info(child2,vart,nc3,name)
! foil normal
        if (name.eq.'foilF') then
          str = '   Foil normal F = '
          self%foil%F = JSONgetDoubleVector(child2,nc3,str,v)
        end if
! foil q vector
        if (name.eq.'foilq') then
          str = '   Foil q-vector = '
          self%foil%q = JSONgetDoubleVector(child2,nc3,str,v)
        end if
! foil alP tilt
        if (name.eq.'foilalP') then
          str = '   Foil alP tilt = '
          self%foil%alP =JSONgetDouble(child2,str,v)
        end if
! foil alS tilt
        if (name.eq.'foilalS') then
          str = '   Foil alS tilt = '
          self%foil%alS =JSONgetDouble(child2,str,v)
        end if
! foil alR tilt
        if (name.eq.'foilalR') then
          str = '   Foil alR tilt = '
          self%foil%alR =JSONgetDouble(child2,str,v)
        end if
! foil thickness
        if (name.eq.'foilz0') then
          str = '   Foil thickness = '
          self%foil%z0 =JSONgetDouble(child2,str,v)
        end if
! foil shape parameters
        if (name.eq.'brx') then
          str = '   shape brx = '
          self%foil%brx =JSONgetDouble(child2,str,v)
        end if
! foil shape parameters
        if (name.eq.'bry') then
          str = '   shape bry = '
          self%foil%bry =JSONgetDouble(child2,str,v)
        end if
! foil shape parameters
        if (name.eq.'brxy') then
          str = '   shape brxy = '
          self%foil%brxy =JSONgetDouble(child2,str,v)
        end if
! foil shape parameters
        if (name.eq.'cpx') then
          str = '   shape cpx = '
          self%foil%cpx =JSONgetDouble(child2,str,v)
        end if
! foil shape parameters
        if (name.eq.'cpy') then
          str = '   shape cpy = '
          self%foil%cpy =JSONgetDouble(child2,str,v)
        end if
! foil elastic modulus tensor  (6x6 format)
        if (name.eq.'foilelmo') then
          str = ''
          do jj=1,6
           call json_get_child(child2,jj,child3)
           call json_info(child3,vart,nc3,name)
           if (name.eq.row(jj)) self%foil%elmo(jj,1:6) = JSONgetDoubleVector(child3,nc3,str,0)
          end do
          if (v.eq.1) then
            call Message%WriteValue('   Elastic moduli tensor (Voigt notation)','')
            do jj=1,6
              io_real(1:6) = self%foil%elmo(jj,1:6)
              call Message%WriteValue('',io_real,6)
            end do
          end if
        end if
    end do nc2loop
  end do
end if

call JSON_failtest(error_cnt)

! verify that the foil normal (in real space) and q (in reciprocal space) are orthogonal
! in other words, we do a cartesian dot product...
x = cell%CalcDot(self%foil%F,self%foil%q,'c')
if (abs(x).gt.0.005) then
  call Message%printMessage('Foil normal F must be orthogonal to q', frm = "(A)")
!  stop
end if

end subroutine JSONreadFoilData_

!--------------------------------------------------------------------------
!
! SUBROUTINE: InitializeDefects
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read defect information, and generate all data that can be precomputed for each defect
!
!> @date  11/22/15 MDG 1.0 original
!> @date  11/24/15 MDG 1.1 added Ydislocations, stacking faults, inclusions and voids
!--------------------------------------------------------------------------
recursive subroutine InitializeDefects_(self,EMsoft,cell,jsonname,npix,npiy,L,gf,error_cnt,verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: InitializeDefects_

use mod_io
use mod_JSONsupport
use mod_crystallography
use mod_EMsoft

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)          :: self
type(Cell_T)                            :: cell
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
type(IO_T)                              :: Message
character(fnlen),INTENT(IN)             :: jsonname
integer(kind=irg),INTENT(IN)            :: npix
integer(kind=irg),INTENT(IN)            :: npiy
real(kind=sgl),INTENT(IN)               :: L
real(kind=sgl),INTENT(IN)               :: gf(3)
integer(kind=irg),INTENT(INOUT)         :: error_cnt

logical,INTENT(IN),OPTIONAL             :: verbose

integer(kind=irg)                       :: v, i, io_int(1)

error_cnt = 0

v = 0
if (PRESENT(verbose)) then
  if (verbose) then
    v = 1
  end if
end if

self%numdisl = 0
self%numYdisl = 0
self%numsf = 0
self%numvoids = 0
self%numinc = 0
self%numEinc = 0

! first of all, we need to read all the defect data from the jsonname file, including the foil data
! note that the JSON file should have the .jsonc extension, isince it may have comments; those lines
! will be removed first
call self%JSONreadDefectFile(EMsoft, cell, jsonname, error_cnt, verbose)

if (v.eq.1) then
  call Message%printMessage('The following defects were initialized : ')
  io_int(1) = self%numdisl
  call Message%WriteValue('  Number of dislocations       : ',io_int,1)
  io_int(1) = self%numYdisl
  call Message%WriteValue('  Number of Yoffe dislocations : ',io_int,1)
  io_int(1) = self%numsf
  call Message%WriteValue('  Number of stacking faults    : ',io_int,1)
  io_int(1) = self%numinc
  call Message%WriteValue('  Number of inclusions         : ',io_int,1)
  io_int(1) = self%numEinc
  call Message%WriteValue('  Number of Eshelby inclusions : ',io_int,1)
  io_int(1) = self%numvoids
  call Message%WriteValue('  Number of voids              : ',io_int,1)
end if

! once we have this data, we need to initialize all other defect related parameters, including
! things like displacement field parameters etc...

! we begin with the foil itself
call self%init_foil_data(cell,npix,npiy,L,v)
if (v.eq.1) call Message%printMessage('========> completed foil generation')

! then we add the defects, starting with all the regular dislocations, if any
if (self%numdisl.ne.0) then
  call self%init_dislocation_data(cell,npix,npiy,gf,L,v)
  call Message%printMessage('========> completed dislocation generation')
end if

! then Ydislocations
if (self%numYdisl.ne.0) then
  call self%init_YSH_dislocation_data(cell,npix,npiy,gf,L,v)
  if (v.eq.1) call Message%printMessage('========> completed Ydislocation generation')
end if

! stacking faults
if (self%numsf.ne.0) then
  call self%init_stacking_fault_data(cell,L,npix,npiy,gf,v)
  if (v.eq.1) call Message%printMessage('========> completed SF generation')
end if

! inclusions
if (self%numinc.ne.0) then
  call self%init_inclusion_data(L,npix,npiy,v)
  if (v.eq.1) call Message%printMessage('========> completed inclusion generation')
end if

! ! Eshelby inclusions
! if (defects%numEinc.ne.0) then
!   do i=1,defects%numEinc
!     call InitializeEshelbyInclusion(cell,defects,i,v,L,npix,npiy)
!   end do
!   if (v.eq.1) call Message('========> completed E-inclusion generation')
! end if

! voids
if (self%numvoids.ne.0) then
  call self%init_void_data(L,npix,npiy,v)
  if (v.eq.1) call Message%printMessage('========> completed void generation')
end if

end subroutine InitializeDefects_


!--------------------------------------------------------------------------
!
! SUBROUTINE: init_dislocation_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  init dislocation namelist files
!
!> @param cell unit cell pointer
!> @param defects defect structure
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf
!> @param L
!> @param dinfo logical to trigger verbose output
!
!> @date 01/05/99 MDG 1.0 original
!> @date 05/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 06/04/13 MDG 3.0 rewrite
!> @date 06/09/14 MDG 4.0 added cell, DL argument
!> @date 11/22/15 MDG 4.1 old routine obsolete with Release 3.1; replaced by JsonreadDefectFile
!> @date 11/23/15 MDG 4.2 moved from dislocation.f90 to defectmodule.f90
!--------------------------------------------------------------------------
recursive subroutine init_dislocation_data_(self,cell,DF_npix,DF_npiy,DF_gf,L,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT :: init_dislocation_data_

use mod_io
use mod_crystallography

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)          :: self
type(Cell_T)                            :: cell
integer(kind=irg),INTENT(IN)            :: DF_npix, DF_npiy, dinfo
real(kind=sgl),INTENT(IN)               :: DF_gf(3), L

integer(kind=irg)                       :: i

! zfrac goes between -0.5 and +0.5, with -0.5 being the top surface and +0.5 the bottom
! this only really matters for dislocations that are parallel to the foil surfaces

! loop over all regular dislocations and initialize their displacement field parameters
   do i=1,self%numdisl  ! +2*defects%numsf   ! we do not deal with partials in stacking faults here ...

! center of dislocation inside the foil is transformed to foil coordinates [nm] with defects%DL(i)%kd=0 (center of foil) [verified 4/23/11]
! the point (0,0) is at the center of the image ... hence the factor of 0.5
    self%DL(i)%id = self%DL(i)%id * 0.5 * float(DF_npix) ! * L   scaling (zooming) is done later in the image reference frame...
    self%DL(i)%jd = self%DL(i)%jd * 0.5 * float(DF_npiy) ! * L
    self%DL(i)%g = DF_gf

! and pre-compute the dislocation displacement field parameters
    call self%makedislocation(cell,i,dinfo, L)
  end do

end subroutine init_dislocation_data_


!--------------------------------------------------------------------------
!
! SUBROUTINE: init_YSH_dislocation_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief init Yoffe dislocation
!
!> @param cell unit cell pointer
!> @param defects defects structure
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf
!> @param L
!> @param dinfo logical to trigger verbose output
!
!> @date  1/5/99  MDG 1.0 original
!> @date  5/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 03/25/13 MDG 3.0 updated IO
!> @date 11/21/13 MDG 3.1 verification
!> @date 06/10/14 MDG 4.0 added defects, cell and foil arguments
!> @date 11/23/15 MDG 4.1 made foil part of defects
!--------------------------------------------------------------------------
recursive subroutine init_YSH_dislocation_data_(self,cell,DF_npix,DF_npiy,DF_gf,L,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT :: init_YSH_dislocation_data_

use mod_io
use mod_crystallography

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)               :: self
type(Cell_T)                  :: cell

integer(kind=irg),INTENT(IN)    :: DF_npix, DF_npiy, dinfo
real(kind=sgl),INTENT(IN)       :: DF_gf(3), L

integer(kind=irg)               :: i
real(kind=sgl)                  :: id,jd,u(3),bv(3),poisson

! these are just the individual dislocations; the ones that belong to
! stacking faults are handled separately
do i=1,self%numYdisl
! top-of-the-foil intersection of dislocation line is transformed to foil coordinates [nm] with DL(i)%kd=0 (center of foil) [verified 4/23/11]
! the point (0,0) is at the center of the image ... hence the factor of 0.5
  self%YD(i)%id = self%YD(i)%id * 0.5 * float(DF_npix) ! * L   scaling (zooming) is done later in the image reference frame...
  self%YD(i)%jd = self%YD(i)%jd * 0.5 * float(DF_npiy) ! * L
  self%YD(i)%g = DF_gf

! and pre-compute the dislocation displacement field parameters
  call self%makeYSHdislocation(cell,i,dinfo, L)
end do

end subroutine init_YSH_dislocation_data_


!--------------------------------------------------------------------------
!
! SUBROUTINE: init_stacking_fault_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  read stacking fault namelist files
!
!> @param cell unit cell pointer
!> @param defects defect structure
!> @param DF_L
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf
!> @param dinfo logical to trigger verbose output
!> @param ECCI logical optional to indicate ECCI formatting rather than regular TEM
!
!> @date    1/5/99  MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!> @date   12/17/13 MDG 3.1 added ECCI mode
!> @date   06/09/14 MDG 4.0 added cell, defects arguments
!> @date   06/10/14 MDG 4.1 added foil argument
!--------------------------------------------------------------------------
recursive subroutine init_stacking_fault_data_(self,cell,DF_L,DF_npix,DF_npiy,DF_g,dinfo,ECCI)
!DEC$ ATTRIBUTES DLLEXPORT :: init_stacking_fault_data_

use mod_io
use mod_crystallography

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)           :: self
type(Cell_T)                            :: cell

integer(kind=irg),INTENT(IN)            :: DF_npix, DF_npiy, dinfo
real(kind=sgl),INTENT(IN)               :: DF_g(3)
real(kind=sgl),INTENT(IN)               :: DF_L
logical,INTENT(IN),OPTIONAL             :: ECCI

integer(kind=irg)                       :: i
real(kind=sgl)                          :: poisson


! read the namelist files for all of the stacking faults
 do i=1,self%numsf
!   SFR = (/ 0.0, 0.0, 0.0 /)
    poisson = 0.0
! transform the fault fractional coordinates to nm in the image reference frame
    self%SF(i)%id = self%SF(i)%id * 0.5 * float(DF_npix) ! * DF_L  (zooming is done later in the image reference frame)
    self%SF(i)%jd = self%SF(i)%jd * 0.5 * float(DF_npiy) ! * DF_L
    self%SF(i)%poisson = poisson

!   if (sum(abs(SFR)).eq.0.0) then
      self%SF(i)%Rdisp = self%SF(i)%lpb
!   else
!     self%SF(i)%Rdisp = SFR
!   end if
! initialize the stacking fault variables and both partial dislocations; this might depend
! on the imaging mode (TEM vs. ECCI); careful here, since the counting of dislocations has
! changed with respect to release 2.0 !!!
    if (present(ECCI)) then
     call self%makestackingfaultECCI(cell,i,DF_L,DF_npix,DF_npiy,DF_g,dinfo)
    self%numYdisl = self%numYdisl + 2
    else
      call self%makestackingfault(cell,i,DF_L,DF_npix,DF_npiy,DF_g,dinfo)
      self%numdisl = self%numdisl + 2
    end if
 end do

end subroutine init_stacking_fault_data_

!--------------------------------------------------------------------------
!
! SUBROUTINE: init_void_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  init void parameters
!
!> @param defects defects structure
!> @param DF_L column edge length
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param dinfo logical to trigger verbose output
!
!> @date 01/05/99 MDG 1.0 original
!> @date 05/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 03/25/13 MDG 3.0 updated IO
!> @date 06/09/14 MDG 4.0 added defects argument
!> @date 06/10/14 MDG 4.1 added foil argument
!> @date 11/23/15 MDG 4.2 removed foil and put it inside defects
!--------------------------------------------------------------------------
recursive subroutine init_void_data_(self,DF_L,DF_npix,DF_npiy,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT :: init_void_data_

use mod_io
use mod_quaternions

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)       :: self

integer(kind=irg),INTENT(IN)        :: dinfo,DF_npix,DF_npiy
real(kind=sgl),INTENT(IN)           :: DF_L

integer(kind=irg)                   :: i
real(kind=sgl)                      :: tmp(3)

type(Quaternion_T)                  :: a_fc_conjg

! read each subsequent line
do i=1,self%numvoids
  self%voids(i)%xpos = self%voids(i)%xpos * 0.5 * float(DF_npix) * DF_L
  self%voids(i)%ypos = self%voids(i)%ypos * 0.5 * float(DF_npiy) * DF_L
  self%voids(i)%zpos = self%voids(i)%zpos * self%foil%z0
! transform to the foil reference frame
  a_fc_conjg =  conjg(self%foil%a_fc)
  tmp = a_fc_conjg%quat_Lp(dble((/ self%voids(i)%xpos, self%voids(i)%ypos, self%voids(i)%zpos /)) )
  self%voids(i)%xpos = tmp(1)
  self%voids(i)%ypos = tmp(2)
  self%voids(i)%zpos = tmp(3)
  if (dinfo.eq.1) write (*,*) i,self%voids(i)%xpos,self%voids(i)%ypos,self%voids(i)%zpos,self%voids(i)%radius
end do

end subroutine init_void_data_

!--------------------------------------------------------------------------
!
! SUBROUTINE: init_inclusion_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  init inclusion parameters
!
!> @param defects defect structure
!> @param DF_L column edge length
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param dinfo logical to trigger verbose output
!
!> @date  01/05/99 MDG 1.0 original
!> @date  05/19/01 MDG 2.0 f90 version
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  03/25/13 MDG 3.0 updated IO
!> @date  06/09/14 MDG 4.0 added defects argument
!> @date  06/10/14 MDG 4.1 added foil argument
!> @date  11/23/15 MDG 4.2 made foil part of defects
!--------------------------------------------------------------------------
recursive subroutine init_inclusion_data_(self,DF_L,DF_npix,DF_npiy,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT :: init_inclusion_data_

use mod_io
use mod_quaternions

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)       :: self

integer(kind=irg),INTENT(IN)        :: dinfo,DF_npix,DF_npiy
real(kind=sgl),INTENT(IN)           :: DF_L

integer(kind=irg)                   :: i
real(kind=sgl)                      :: tmp(3)

type(Quaternion_T)                  :: a_fc_conjg
! read each subsequent line
do i=1,self%numinc
  self%inclusions(i)%xpos = self%inclusions(i)%xpos * 0.5 * float(DF_npix)*DF_L
  self%inclusions(i)%ypos = self%inclusions(i)%ypos * 0.5 * float(DF_npiy)*DF_L
  self%inclusions(i)%zpos = self%inclusions(i)%zpos * self%foil%z0         ! vertical fractional location in interval [-1,1]
  a_fc_conjg =  conjg(self%foil%a_fc)
  tmp = a_fc_conjg%quat_Lp(dble((/ self%inclusions(i)%xpos, self%inclusions(i)%ypos, &
        self%inclusions(i)%zpos /)) )
  self%inclusions(i)%xpos = tmp(1)
  self%inclusions(i)%ypos = tmp(2)
  self%inclusions(i)%zpos = tmp(3)
  if (dinfo.eq.1) write (*,*) i, self%inclusions(i)%xpos, self%inclusions(i)%ypos, self%inclusions(i)%zpos, &
                              self%inclusions(i)%radius, self%inclusions(i)%C
end do

end subroutine init_inclusion_data_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! here are the init_defect_data subroutines
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: init_foil_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  initializes the foil geometry data
!
!> @param cell unit cell pointer
!> @param defects defects structure
!> @param npix number of x image pixels
!> @param npiy number of y image pixels
!> @param L pixel size for column approximation
!> @param dinfo flag to print information
!
!> @date 11/22/15 MDG 1.0 new routine in Release 3.1; read portion replaced with JSONreadFoilData in JSONsupport
!--------------------------------------------------------------------------
recursive subroutine init_foil_data_(self,cell,npix,npiy,L,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT :: init_foil_data_

use mod_crystallography
use mod_io
use mod_global
use mod_rotations
use mod_kinds

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)  :: self
type(Cell_T)                    :: cell
type(IO_T)                      :: Message

integer(kind=irg),INTENT(IN)    :: npix, npiy, dinfo
real(kind=sgl),INTENT(IN)       :: L

real(kind=sgl)                  :: io_real(1)
real(kind=dbl)                  :: amat(3,3)
type(q_T)                       :: a_fm_qu
type(o_T)                       :: amat_om

! assign these values to the appropriate slots in foil%   [verified 4/23/11]
self%foil%alP = self%foil%alP*cPi/180.0            ! convert the tilt and rotation angles to radians
self%foil%alS = self%foil%alS*cPi/180.0
self%foil%alR = self%foil%alR*cPi/180.0
self%foil%npix = npix                        ! image size (this duplicates some values, but it's easier this way)
self%foil%npiy = npiy

! shape parameters
self%foil%cpx = self%foil%cpx * float(npix) * 0.5 * L  ! we'll define the foil shape center w.r.t. to the center of the image in [nm] coordinates
self%foil%cpy = self%foil%cpy * float(npiy) * 0.5 * L  !

! initialize a bunch of foil related quantities, using quaternions for all rotations
call self%initialize_foil_geometry(cell, dinfo)

! compute the projected thickness
a_fm_qu = q_T( qdinp =  self%foil%a_fm%get_quatd() )

amat_om = a_fm_qu%qo()
amat = amat_om%o_copyd()

self%foil%zb = self%foil%z0/amat(3,3)
if (dinfo.eq.1) then
  io_real(1)=self%foil%z0
  call Message%WriteValue('Nominal foil thickness = ', io_real, 1, "(F8.3)")
  io_real(1)=self%foil%zb
  call Message%WriteValue('Effective foil thickness = ', io_real, 1, "(F8.3/)")
end if

end subroutine init_foil_data_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! and here are the routines that actually compute all the defect parameters
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: initialize_foil_geometry
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  Initializes the foil geometry
!
!> @details This new implementation uses quaternions for all rotations.
!
!> @param cell unit cell pointer
!> @param dinfo
!
!> @date  1/ 5/99 MDG 1.0 original
!> @date  1/11/10 MDG 2.0 rewrite of beam direction part
!> @date  3/28/11 MDG 2.1 code verified
!> @date  4/23/11 MDG 2.2 redefined origin to be at center of image
!> @date  6/03/13 MDG 3.0 replaced rotation matrices by quaternions throughout
!> @date 10/30/13 MDG 3.1 complete debug of quaternion and rotation implementation
!> @date 06/09/14 MDG 4.0 added cell and foil as argument
!--------------------------------------------------------------------------
recursive subroutine initialize_foil_geometry_(self,cell,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT :: initialize_foil_geometry_

use mod_kinds
use mod_math
use mod_global
use mod_crystallography
use mod_symmetry
use mod_io
use mod_rotations
use mod_quaternions

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)  :: self
type(Cell_T)                    :: cell
integer(kind=sgl),INTENT(IN)    :: dinfo
type(orientation_T)             :: ot
type(IO_T)                      :: Message
type(o_T)                       :: a_fc_om
real(kind=dbl)                  :: ey(3), ex(3), tt, dx, dy
real(kind=sgl)                  :: io_real(3)
real(kind=dbl)                  :: cp,sp,cs,ss,cr,sr, ca, sa, a_fc_array(3,3)
integer(kind=irg)               :: i,j
!type(orientationtyped)          :: ot
character(10)                   :: pret

type(Quaternion_T)              :: qq1, qq2, qqf, quat, a_fm_conjg, a_fc_conjg
type(q_T)                       :: qu, a_fc_qu, a_mi_qu, a_fm_qu, a_mc_qu, a_ic_qu, a_fi_qu

associate( a_fm => self%foil%a_fm, a_mi => self%foil%a_mi, a_fc => self%foil%a_fc, &
          a_mc => self%foil%a_mc, a_ic => self%foil%a_ic, a_fi => self%foil%a_fi, &
          bry => self%foil%bry, cpy => self%foil%cpy, brxy => self%foil%brxy)

! determine the foil-to-microscope transformations [verified on 4/23/11, converted to quaternions on 6/4/13,
! verified 10/30/13, and again on 11/11/13 after some changes elsewhere]
  if (self%foil%alR.eq.0.D0) then
! the double tilt holder transformation a_fm; note quaternions, hence we need the half-angles !
! a_fm transforms a vector v FROM the microscope reference frame To the foil reference frame
! using the quat_rotate_vector routine.
    cp = dcos(self%foil%alP*0.5D0)
    sp = dsin(self%foil%alP*0.5D0)
    ca = dcos(self%foil%alP)
    sa = dsin(self%foil%alP)
    cs = dcos(self%foil%alS*0.5D0)
    ss = dsin(self%foil%alS*0.5D0)

    qq1 = Quaternion_T( qd = (/ cs, 0.D0, ss*ca, ss*sa /) )
    qq2 = Quaternion_T( qd = (/ cp, sp, 0.D0, 0.D0 /) )

    a_fm = conjg(qq1*qq2)
  else
! the rotation tilt holder transformation a_fm [verified on 4/23/11, converted to quaternions on 6/4/13,
! and again on 11/11/13 after changes elsewhere]
    cp = dcos(self%foil%alP*0.5D0)
    sp = dsin(self%foil%alP*0.5D0)
    cr = dcos(self%foil%alR*0.5D0)
    sr = dsin(self%foil%alR*0.5D0)
    ca = dcos(self%foil%alP)
    sa = dsin(self%foil%alP)

    qq1 = Quaternion_T( qd =  (/ cr,0.D0, -sr*sa, sr*ca /) )
    qq2 = Quaternion_T( qd = (/ cp, sp,0.D0,0.D0  /))

    a_fm = conjg(qq1*qq2)
  end if

  a_fc_qu = q_T( qdinp =  a_fm%get_quatd() )

  if (dinfo.eq.1) then
    pret = 'a_fm: '
    ot = orientation_T(a_fc_qu)
    call ot%print_orientation('d')
  end if

! a_mi (image to microscope) apart from a scale factor, these two are identical
! The EM book uses a beta rotation angle between the image and the microscope,
! but that is really not necessary because we already fix the image with respect to
! the microscope by defining q (the horizontal image direction) to point to the
! airlock. [verified 4/23/11, converted to quaternions on 6/4/13]
! So we'll keep this transformation equal to the identity at all times.

  a_mi = Quaternion_T(qd =  (/ 1.D0,0.D0,0.D0,0.D0 /))   ! identity quaternion
  a_mi_qu = q_T( qdinp =  a_mi%get_quatd() )
  if (dinfo.eq.1) then
    pret = 'a_mi: '
    ot = orientation_T(a_mi_qu)
    call ot%print_orientation('d')
  end if

! This allows us to get the beam direction, since we know the foil normal and tilt angles
! The beam direction is the inverse transform of the microscope e_z-axis to the foil reference frame [verified 11/12/13]
  a_fm_conjg = conjg(a_fm)

  self%foil%B = a_fm_conjg%quat_Lp((/ 0.D0,0.D0,-1.D0 /))
  self%foil%Bn = self%foil%B

  call cell%NormVec(self%foil%Bn,'c')
  if (dinfo.eq.1) then
    io_real(1:3) = self%foil%B(1:3)
    call Message%WriteValue('  Beam direction (foil reference frame) = ',io_real,3,"('[',F12.5,',',F12.5,',',F12.5,']')")
  end if

! transform both the foil normal and the q-vector to the crystal cartesian reference frame (eq. 8.8) [verified 4/23/11,
! and again on 11/12/13 afterchanges elsewhere]
  call cell%TransSpace(self%foil%F,self%foil%Fn,'d','c')
  call cell%TransSpace(self%foil%q,self%foil%qn,'r','c')
  call cell%NormVec(self%foil%Fn,'c')
  call cell%NormVec(self%foil%qn,'c')

! a_fc (crystal to foil)
  a_fc_array(3,1:3) = self%foil%Fn(1:3)
  a_fc_array(1,1:3) = self%foil%qn(1:3)
  call cell%CalcCross(self%foil%Fn,self%foil%qn,ey,'c','c',0)
  call cell%NormVec(ey,'c')
  a_fc_array(2,1:3) = ey(1:3)

  a_fc_om = o_T(odinp = a_fc_array)

  a_fc_qu = a_fc_om%oq()
  a_fc = Quaternion_T(qd = a_fc_qu%q_copyd())

  if (dinfo.eq.1) then
   pret = 'a_fc: '
   ot = orientation_T(a_fc_qu)
   call ot%print_orientation('d')
  end if

! a_mc (crystal to microscope)
  a_mc  = a_fm_conjg*a_fc
  a_mc_qu = q_T(qdinp =  a_mc%get_quatd())

  if (dinfo.eq.1) then
    pret = 'a_mc: '
    ot = orientation_T(a_mc_qu)
    call ot%print_orientation('d')
  end if


! a_ic (crystal to image)
  a_ic = conjg(a_mi)*a_mc
  a_ic_qu = q_T(qdinp =  a_ic%get_quatd())

  if (dinfo.eq.1) then
    pret = 'a_ic: '
    ot = orientation_T(a_ic_qu)
    call ot%print_orientation('d')
  end if

! a_fi (image to foil)
  a_fi = a_fc*conjg(a_ic)
  a_fi_qu = q_T(qdinp =  a_fi%get_quatd())

  if (dinfo.eq.1) then
    pret = 'a_fi: '
    ot = orientation_T(a_fi_qu)
    call ot%print_orientation('d')
  end if


! express the beam direction in the Bravais reference frame [verified 4/23/11, and again on 11/12/13
! after changes elsewhere]
  a_fc_conjg = conjg(a_fc)
  ex = a_fc_conjg%quat_Lp(dble(self%foil%Bn))

  call cell%TransSpace(ex,ey,'c','d')
  call cell%NormVec(ey,'c')
  if (dinfo.eq.1) then
    io_real(1:3) = ey(1:3)
    call Message%WriteValue('  Beam direction (crystal reference frame) = ', io_real, 3, "('[',F12.5,',',F12.5,',',F12.5,']'/)")
  end if

! define the foil shape (for now as an elliptic paraboloid z = brx * (x-xc)^2 + bry * (y-yc)^2)
if (.not.allocated(self%foil%sg)) allocate(self%foil%sg(self%foil%npix,self%foil%npiy))
! if the foil is not bent, then we set this array to zero, otherwise we compute the elliptical paraboloid
if ((self%foil%brx.eq.0.0).and.(bry.eq.0.0)) then
  if (dinfo.eq.1) then
    call Message%printMessage(' Initializing a flat foil ', frm = "(A)")
  end if
  self%foil%sg = 0.0
else
  dx = self%foil%npix*0.5
  dy = self%foil%npiy*0.5
  do i=1,self%foil%npix
   tt = self%foil%brx * (float(i)-dx-cpy)**2
    do j=1,self%foil%npiy
! initialize the foil shape function; we assume that the center of the elliptic paraboloid is at location (cpx,cpy)
! presumably, this surface could also be a saddle point if the brx and bry values have opposite sign ...
      self%foil%sg(i,j) = tt + bry * (float(j)-dy-cpy)**2+ 2.0*brxy * (float(j)-dy-cpy)*(float(i)-dx-cpy)
    end do
  end do
  if (dinfo.eq.1) then
    call Message%printMessage(' Initializing a bent foil ', frm = "(A)")
    io_real(1)=minval(self%foil%sg); io_real(2)=maxval(self%foil%sg);
    call Message%WriteValue('Range of local excitation error deviations : ', io_real, 2, "(F10.6,',',F10.6/)")
  end if
end if

end associate

end subroutine initialize_foil_geometry_


!--------------------------------------------------------------------------
!
! SUBROUTINE: makedislocation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  Compute the dismat displacement matrix for a given dislocation
!
!> @details This subroutine computes the matrix dismat that describes the displacement field
!> of a dislocation.  The routine needs the elastic moduli tensor, the transformation
!> matrix between the crystal and dislocation reference frames, and the dislocation
!> Burgers vector.  The routine computes the arrays dismat and pa, which should be used as follows:
!>
!> R_k = 2.0*real([ sum_a=1^3 (dismat(k,a)*log(Z_a)) ]),
!>
!> with Z_a = x_1 + pa(a)*x_2
!>
!>  [see CalcR subroutine for more information]
!>
!> We must also make sure that the x=0 plane of the defect reference frame contains the
!> incident beam direction, to avoid getting stacking-fault fringes in the wrong plane...
!> Actual stacking faults are added in using a different module (stacking_fault.f90).
!>
!> @param cell unit cell pointer
!> @param defects defect structure
!> @param inum
!> @param dinfo
!> @param DF_L column width
!
!> @date   1/ 5/99 MDG 1.0 original
!> @date   5/19/01 MDG 2.0 f90 version
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  06/04/13 MDG 3.0 rewrite
!> @date  10/30/13 MDG 3.1 debug of all rotation parts
!> @date  06/09/14 MDG 4.0 added cell, defects arguments
!> @date  06/10/14 MDG 4.1 added foil as argument
!> @date  11/21/15 MDG 4.2 moved foil into defects structure
!--------------------------------------------------------------------------
recursive subroutine makedislocation_(self,cell,inum,dinfo,DF_L)
!DEC$ ATTRIBUTES DLLEXPORT :: makedislocation_

use mod_math
use mod_global
use mod_crystallography
use mod_io
use mod_symmetry
use mod_rotations
use mod_quaternions

IMPLICIT NONE

class(Defect_T), INTENT(INOUT)          :: self
type(Cell_T)                            :: cell
!f2py intent(in,out) ::  defects
integer(kind=irg),INTENT(IN)            :: inum
integer(kind=irg),INTENT(IN)            :: dinfo
real(kind=sgl),INTENT(IN)               :: DF_L

type(foiltype)                          :: foil
real(kind=dbl)                          :: zz,zang,zmin
real(kind=sgl)                          :: ec(6,6),lec(6,6)
real(kind=dbl)                          :: a_dc(3,3),tmp(3),ex(3),ey(3)
real(kind=dbl)                          :: Bij(3,3),Hij(3,3)
complex(kind=dbl)                       :: a(0:6),b(0:6),c(0:6),d(0:6),e(0:6),ff(0:6),tt(5,0:6),s(0:6),roots(6), &
                                           zero,pasq(3),mat(3,3),aka(3,3),Lia(3,3),Mai(3,3),v(3),pas
integer(kind=irg)                       :: i,j,k,l,imin,ind(3),jnd(3)
type(Quaternion_T)                      :: a_fc_conjg
type(orientation_T)                     :: ot
type(o_T)                               :: a_dc_om
type(q_T)                               :: a_dc_qu
associate(foil => self%foil)

! convert line direction and g-vector to the Cartesian crystal reference frame
call cell%TransSpace(self%DL(inum)%u,self%DL(inum)%un,'d','c')
call cell%TransSpace(self%DL(inum)%g,self%DL(inum)%gn,'r','c')
! normalize both vectors
call cell%NormVec(self%DL(inum)%un,'c')
call cell%NormVec(self%DL(inum)%gn,'c')

! first find the length of the dislocation line inside the foil along the
! dislocation z-axis, which is the line direction; also, compute the intersection
! points of the line with the top and bottom surfaces  (all components in [nm])
zang = cell%CalcAngle(self%DL(inum)%u,dble(foil%F),'d')
zz = cos(zang)
if (abs(zz).gt.0.00001) then
  self%DL(inum)%zu = 0.5*foil%z0/zz
else
  self%DL(inum)%zu = 100000.0         ! this is when the dislocation is nearly parallel to the foil
end if

! transform the line direction to the foil reference frame
tmp =  foil%a_fc%quat_Lp(dble(self%DL(inum)%un) ) / DF_L

if (dinfo.eq.1) then
  write (*,*) 'transformed line direction ', tmp, zang, zz
end if

! determine the top and bottom intersection coordinates
if (zz.gt.0.0) then  ! u points to the top of the foil
    self%DL(inum)%top = (/ self%DL(inum)%id + tmp(1)*self%DL(inum)%zu, self%DL(inum)%jd + &
        tmp(2)*self%DL(inum)%zu, 0.5D0*foil%z0 /)
    self%DL(inum)%bottom = (/ self%DL(inum)%id - tmp(1)*self%DL(inum)%zu, self%DL(inum)%jd - &
        tmp(2)*self%DL(inum)%zu, -0.5D0*foil%z0 /)
else                 ! u points to the bottom of the foil
    self%DL(inum)%top = (/ self%DL(inum)%id - tmp(1)*self%DL(inum)%zu, self%DL(inum)%jd - &
        tmp(2)*self%DL(inum)%zu, -0.5D0*foil%z0 /)
    self%DL(inum)%bottom = (/ self%DL(inum)%id + tmp(1)*self%DL(inum)%zu, self%DL(inum)%jd + &
        tmp(2)*self%DL(inum)%zu, 0.5D0*foil%z0 /)
end if

if (dinfo.eq.1) then
  write (*,*) self%DL(inum)%id,self%DL(inum)%jd
  write (*,*) 'dislocation top intersection at ',self%DL(inum)%top
  write (*,*) 'dislocation bottom intersection at ',self%DL(inum)%bottom
end if

! a_dc (crystal to defect)  matrix corrected on 11/29/10 to put defect x-axis in the plane of u and B
if (dinfo.eq.1) then
  write (*,*) 'cartesian quantities'
  write (*,*) 'unit line direction = ',self%DL(inum)%un
  write (*,*) 'unit beam direction = ',foil%Bn
end if

! transform beam direction (currently in foil frame) to cartesian
a_fc_conjg = conjg(foil%a_fc)
tmp = a_fc_conjg%quat_Lp(dble(foil%Bn))
!tmp = quat_rotate_vector(conjg(foil%a_fc), (/ 0.0D0, 0.0D0, -1.0D0/) )
call cell%NormVec(tmp,'c')

! the defect z axis is the line direction and x is in the plane of u and B to avoid the intrinsic discontinuity (cut plane)
a_dc(3,1:3) = self%DL(inum)%un(1:3)
call cell%CalcCross(dble(self%DL(inum)%un),tmp,ex,'c','c',0)
call cell%NormVec(ex,'c')
a_dc(1,1:3) = ex(1:3)
call cell%CalcCross(dble(self%DL(inum)%un),ex,ey,'c','c',0)
call cell%NormVec(ey,'c')
a_dc(2,1:3) = ey(1:3)
a_dc_om = o_T(odinp = a_dc)
a_dc_qu = a_dc_om%oq()
self%DL(inum)%a_dc = Quaternion_T(qd = a_dc_qu%q_copyd())

!if (dinfo.eq.1) then
!  call PrintMatrixd('a_dc',a_dc)
!end if

! a_di (image to defect)
self%DL(inum)%a_di = self%DL(inum)%a_dc*conjg(foil%a_ic)
self%DL(inum)%a_id = conjg(self%DL(inum)%a_di)

!if (dinfo.eq.1) then
!  call print_orientation_d(init_orientation_d(self%DL(inum)%a_di,'qu'),'om','a_di:     ')
!  call print_orientation_d(init_orientation_d(self%DL(inum)%a_id,'qu'),'om','a_id:     ')
!end if

! finally, get the foil to defect transformation (used in defect module)
self%DL(inum)%a_df = self%DL(inum)%a_di*conjg(foil%a_fi)

! Burgers vector (in the defect reference frame !!!)
! first transform Burgers vector to crystal cartesian reference frame
call cell%TransSpace(dble(self%DL(inum)%burg),tmp,'d','c')
! then convert this to the defect reference frame
self%DL(inum)%burgd(1:3) = self%DL(inum)%a_dc%quat_Lp(dble(tmp))

if (dinfo.eq.1) then
  write (*,*) 'rotated burgers vector  = ', self%DL(inum)%burgd(1:3)
end if

! transform the elastic moduli
lec = foil%elmo

! transform lec to defect reference frame
!a_dc = qu2om(self%DL(inum)%a_dc)
a_dc_om = a_dc_qu%qo()

self%DL(inum)%a_dc = Quaternion_T(qd=a_dc_qu%q_copyd())


 call TransFourthRankTensor(a_dc,lec,ec)
if (dinfo.eq.1)  then
  write (*,*) 'Elasticity tensor in defect reference frame'
  do i=1,6
    write (*,"(6(F8.4,2x))") (ec(i,j),j=1,6)
  end do
  write (*,*) '----'
end if

! next, create the sextic polynomial
zero = cmplx(0.0,0.0,dbl)
a=zero; b=zero; c=zero; d=zero; e=zero; ff=zero
a(0:2) = (/ cmplx(ec(1,1),0.0,dbl), cmplx(ec(1,6)*2.0,0.0,dbl),     cmplx(ec(6,6),0.0,dbl) /)
b(0:2) = (/ cmplx(ec(6,6),0.0,dbl), cmplx(ec(2,6)*2.0,0.0,dbl),     cmplx(ec(2,2),0.0,dbl) /)
c(0:2) = (/ cmplx(ec(5,5),0.0,dbl), cmplx(ec(4,5)*2.0,0.0,dbl),     cmplx(ec(4,4),0.0,dbl) /)
d(0:2) = (/ cmplx(ec(5,6),0.0,dbl), cmplx(ec(4,6)+ec(2,5),0.0,dbl), cmplx(ec(2,4),0.0,dbl) /)
e(0:2) = (/ cmplx(ec(1,5),0.0,dbl), cmplx(ec(1,4)+ec(5,6),0.0,dbl), cmplx(ec(4,6),0.0,dbl) /)
ff(0:2) = (/ cmplx(ec(1,6),0.0,dbl), cmplx(ec(1,2)+ec(6,6),0.0,dbl), cmplx(ec(2,6),0.0,dbl) /)
tt = zero
s = zero

! matrix elements
do j=0,6
 do i=0,j
  tt(1,j) = tt(1,j) + a(j-i)*b(i)
  tt(2,j) = tt(2,j) + d(j-i)*e(i)
  tt(3,j) = tt(3,j) + a(j-i)*d(i)
  tt(4,j) = tt(4,j) + b(j-i)*e(i)
  tt(5,j) = tt(5,j) + c(j-i)*ff(i)
 end do
end do

! determinant leading to the sextic equation
do j=0,6
 do i=0,j
  s(j) = s(j) + tt(1,j-i)*c(i) + 2.0*tt(2,j-i)*ff(i) - tt(3,j-i)*d(i) - tt(4,j-i)*e(i) - tt(5,j-i)*ff(i)
 end do
end do

! get the complex root pairs
call zroots(s,roots)

! then, solve the equation for the vector A_k using the roots with positive imaginary part.
k=1
do j=1,6
  if (aimag(roots(j)).gt.0.0) then
    self%DL(inum)%pa(k) = roots(j)
    k=k+1
  end if
end do

! renumber them to avoid the symmetry degeneracy (see page 328 Head et al.)
v(1:3) = ec(5,5) + 2.0*self%DL(inum)%pa(1:3)*ec(4,5) + ec(4,4)*self%DL(inum)%pa(1:3)**2
zmin = 100.0
imin = 0

! where is the smallest value ?
do i=1,3
  if (abs(v(i)).lt.zmin) then
    imin=i
    zmin = abs(v(i))
  end if
end do

! is the 3rd one the smallest ? if not, then swap with the current 3rd one.
if (imin.ne.3) then
  pas = self%DL(inum)%pa(imin)
  self%DL(inum)%pa(imin)=self%DL(inum)%pa(3)
  self%DL(inum)%pa(3)=pas
end if
  pas = self%DL(inum)%pa(1)
  self%DL(inum)%pa(1)=self%DL(inum)%pa(2)
  self%DL(inum)%pa(2)=pas

! eliminate really small numbers
do i=1,3
  if (abs(aimag(self%DL(inum)%pa(i))).lt.1.0e-8)  self%DL(inum)%pa(i)=cmplx(real(self%DL(inum)%pa(i)),0.0,dbl)
  if (abs(real(self%DL(inum)%pa(i))).lt.1.0e-8)   self%DL(inum)%pa(i)=cmplx(0.0,aimag(self%DL(inum)%pa(i)),dbl)
end do
if (dinfo.eq.1) then
  write (*,*) ' sextic roots'
  do i=1,3
    write (*,*) self%DL(inum)%pa(i)
  end do
  write (*,*) '---'
end if

!  compute the A_ka vectors (see description on page 328 of Head et al.)
pasq = self%DL(inum)%pa**2
if (dinfo.eq.1) write (*,*) 'Aka vectors'
do k=1,3
  mat = zero
  mat(1,1) = ec(1,1)+2.D0*self%DL(inum)%pa(k)*ec(1,6)+ec(6,6)*pasq(k)
  mat(2,2) = ec(6,6)+2.D0*self%DL(inum)%pa(k)*ec(2,6)+ec(2,2)*pasq(k)
  mat(3,3) = ec(5,5)+2.D0*self%DL(inum)%pa(k)*ec(4,5)+ec(4,4)*pasq(k)
  mat(2,3) = ec(5,6)+self%DL(inum)%pa(k)*(ec(4,6)+ec(2,5))+ec(2,4)*pasq(k)
  mat(1,3) = ec(1,5)+self%DL(inum)%pa(k)*(ec(1,4)+ec(5,6))+ec(4,6)*pasq(k)
  mat(1,2) = ec(1,6)+self%DL(inum)%pa(k)*(ec(1,2)+ec(6,6))+ec(2,6)*pasq(k)
  if (k.eq.1) then
    aka(1,1) = mat(2,2)*mat(3,3)-mat(2,3)*mat(2,3)
    aka(1,2) = mat(1,3)*mat(2,3)-mat(1,2)*mat(3,3)
    aka(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
  end if
  if (k.eq.2) then
    aka(2,1) = mat(1,3)*mat(2,3)-mat(1,2)*mat(3,3)
    aka(2,2) = mat(1,1)*mat(3,3)-mat(1,3)*mat(1,3)
    aka(2,3) = mat(1,3)*mat(1,2)-mat(1,1)*mat(2,3)
  end if
  if (k.eq.3) then
    aka(3,1) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
    aka(3,2) = mat(1,3)*mat(1,2)-mat(1,1)*mat(2,3)
    aka(3,3) = mat(1,1)*mat(2,2)-mat(1,2)*mat(1,2)
  end if
  if (dinfo.eq.1) write (*,*) k,(aka(k,j),j=1,3)
end do
aka = transpose(aka)

! next, create the L_ialpha matrix
ind = (/ 6, 2, 4 /)
jnd = (/ 1, 6, 5 /)
Lia = zero
do i=1,3
 do j=1,3
  do k=1,3
   Lia(i,j) = Lia(i,j)+(ec(ind(i),jnd(k))+self%DL(inum)%pa(j)*ec(ind(i),ind(k)))*aka(k,j)
  end do
 end do
end do
! if (dinfo.eq.1)  call PrintMatrixcd('Lia ',Lia)

! and invert it
call cInvert(Lia,Mai)
! if (dinfo.eq.1)  call PrintMatrixcd('Mai ',Mai)

! compute Bij ( real matrix )
Bij = 0.D0
do i=1,3
 do j=1,3
  do k=1,3
   Bij(i,j) = Bij(i,j) - aimag(aka(i,k))*real(Mai(k,j)) - real(aka(i,k))*aimag(Mai(k,j))
  end do
 end do
end do
! if (dinfo.eq.1)  call PrintMatrixd('Bij ',Bij)

! and invert to get Hij
call mInvert(Bij,Hij,.FALSE.)
! if (dinfo.eq.1) call PrintMatrixd('Hij ',Hij)

! compute matrix (this is what actually gets to be used for the
! displacement field); needs to know the Burgers vector.
self%DL(inum)%dismat = zero
do k=1,3
 do l=1,3
   do i=1,3
    do j=1,3
     self%DL(inum)%dismat(k,l) = self%DL(inum)%dismat(k,l) + self%DL(inum)%burgd(i)*Hij(j,i)*Mai(l,j)*aka(k,l)
    end do
   end do
 end do
end do

! scale by 1/4pi
self%DL(inum)%dismat = self%DL(inum)%dismat*0.25D0/cPi
! if (dinfo.eq.1)  call PrintMatrixcd('dismat',self%DL(inum)%dismat)

! and return to calling routine
end associate

end subroutine makedislocation_


!--------------------------------------------------------------------------
!
! FUNCTION: YSHDisp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  compute the displacement field of an inclined dislocation intersecting the foil surface
!
!> @details compute the displacement field of an inclined dislocation intersecting the top surface of
!> the foil, taking into account surface relaxations for the isotropic elastic case (cubic only) ...
!>
!> equations are based on the Shaibani&Hazzledine 1981 paper, along with special limits for
!> the alpha->0 case, which were derived by MDG using Mathematica.
!
!> @paraqm defects defects structure
!> @param x dislocation x-coordinate
!> @param y dislocation y-coordinate
!> @param z dislocation z-coordinate
!> @param ii dislocation number
!
!> @todo There is a problem with dislocations normal to the foil surface, likely a typographical error
!> in the SH paper; this needs to be resolved further, which may require explicit repetition of all
!> analytical computations! Mathematica gives an infinite limit for the bx edge case when normal
!> to the foil surface.
!
!> @date    1/5/99  MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!> @date   11/21/13 MDG 3.1 verification
!> @date   06/09/14 MDG 4.0 added defects as argument
!--------------------------------------------------------------------------
recursive function YSHDisp_(self,x,y,z,ii) result(res)
!DEC$ ATTRIBUTES DLLEXPORT :: YSHDisp_

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)  :: self
real(kind=dbl),INTENT(IN)       :: x,y,z
integer(kind=irg),INTENT(IN)    :: ii

real(kind=dbl)                  :: eta, zeta, etap, zetap, r, oms, omts, xx, sgn, om, omp, AA, BB, BBp, th, &
                                 k, lam, alA, alB, u, v, w, ms, S, du, dv, dw, qe, me, De, qx, mx, Dx, rr, eps
real(kind=dbl)                  :: res(3)

! initialize the geometrical parameters
eta  = y*self%YD(ii)%ca - z*self%YD(ii)%sa
zeta = y*self%YD(ii)%sa + z*self%YD(ii)%ca
etap  = -y*self%YD(ii)%ca - z*self%YD(ii)%sa
zetap = y*self%YD(ii)%sa - z*self%YD(ii)%ca
r = sqrt(x**2+y**2+z**2)
oms = 1.D0-self%YD(ii)%sig
omts = 1.D0-2.D0*self%YD(ii)%sig

! cover the special case of negative x values (based on IDL tests)
xx = x
sgn = 1.D0
if (xx.lt.0.D0) then
  xx = dabs(x)
  sgn = -1.D0
else
  sgn = 1.D0
end if

! more parameters
om =  (datan2(y,xx)-datan2(eta,xx)+datan2(xx*r*self%YD(ii)%sa,eta*y+xx**2*self%YD(ii)%ca))
omp= (datan2(y,xx)-datan2(etap,xx)+datan2(xx*r*self%YD(ii)%sa,etap*y-xx**2*self%YD(ii)%ca))

AA = r-z
BB  = r-zeta
BBp = r-zetap
th = 2.D0*oms*(omp-om)
lam = omts*dlog(BBp/BB)
alA = dlog(AA)
alB = dlog(BB)


u = 0.D0
v = 0.D0
w = 0.D0
eps = 1.0D-6

! screw component first
if (abs(self%YD(ii)%bs).gt.eps) then
  ms = xx*sin(2.D0*self%YD(ii)%alpha)/r/BB
  S = self%YD(ii)%bs/(4.D0*cPi)
  if (self%YD(ii)%alpha.gt.0.01) then
    du = xx*ms+2.D0*eta*self%YD(ii)%ca**2/BB+2.D0*omts*self%YD(ii)%cota*(-1.D0+self%YD(ii)%ca+&
            self%YD(ii)%ca*alA-y*self%YD(ii)%sa/AA-alB)-sin(2.D0*self%YD(ii)%alpha)
    dv = y*ms-2.D0*xx*self%YD(ii)%ca/BB-self%YD(ii)%sa*(omp-om)+2.D0*omts*self%YD(ii)%cota* &
         (xx*self%YD(ii)%sa/AA-om*self%YD(ii)%ca)
    dw = z*ms+self%YD(ii)%ca*(omp-om)-2.D0*omts*om*self%YD(ii)%ca
  else
    du = 2.D0*y/(r-z)
    dv = -2.D0*xx*(r+z)/(xx**2+y**2)
    dw = cPi + datan2(y,xx) - datan2(-y,xx)
  end if
  u = u+du*S
  v = v-sgn*dv*S
  w = w+sgn*dw*S
end if

! then the edge component in the y-z plane
if (abs(self%YD(ii)%be).gt.eps) then
  qe = xx*(1.D0/BBp-1.D0/BB+2.D0*z*self%YD(ii)%ca/BB**2)
  me = -qe/r-4.D0*oms*xx*self%YD(ii)%ca**2/r/BB
  De = self%YD(ii)%be/(8.D0*cPi*oms)
  if (self%YD(ii)%alpha.gt.0.01) then
    k = 4.D0*oms*omts*self%YD(ii)%cota**2
    du = xx*me+lam+2.D0*self%YD(ii)%ca*(z+2.D0*oms*eta*self%YD(ii)%sa)/BB-4.D0*oms*self%YD(ii)%sa**2+&
           k*(1.D0-self%YD(ii)%ca-self%YD(ii)%ca*alA+y*self%YD(ii)%sa/AA+alB)
    dv = y*me+qe*self%YD(ii)%sa+th*self%YD(ii)%ca+k*(-xx*self%YD(ii)%sa/AA+om*self%YD(ii)%ca)
    dw = z*me+qe*self%YD(ii)%ca+th*self%YD(ii)%sa-2.D0*xx*self%YD(ii)%ca*(1.D0/BBp+omts/BB)+k*om*self%YD(ii)%sa
!    write (*,*) du,dv,dw
  else
    rr = xx**2+y**2
    du = 2.D0*z/(r-z)+4.D0*xx**2*(self%YD(ii)%sig*rr-r**2)/r/AA**2/(r+z)+2.D0*omts*oms*((xx**2+z*(z-r))/AA**2+alA)+&
         omts*dlog((r+z)/AA)
    dv = 4.D0*xx*y*(rr*self%YD(ii)%sig-r**2)/r/AA**2/(r+z)+2.D0*xx*y*(rr+2.D0*z*(r+z))*oms*omts/rr**2+&
            2.D0*oms*(cPi + datan2(y,xx) - datan2(-y,xx))
    dw = 4.D0*xx*rr*self%YD(ii)%sig*(z-2.D0*r*oms)/r/AA**2/(r+z)
  end if
  u = u+du*De
  v = v+sgn*dv*De
  w = w+sgn*dw*De
end if

! and finally the bx edge component
if (abs(self%YD(ii)%bx).gt.eps) then
  qx = etap/BBp-eta/BB-2.D0*z*eta*self%YD(ii)%ca/BB**2
  mx = -qx/r+2.D0*omts*y*self%YD(ii)%ca/r/BB
  Dx = self%YD(ii)%bx/(8.D0*cPi*oms)
  if (self%YD(ii)%alpha.gt.0.01) then
    k = 4.D0*oms*omts*self%YD(ii)%cota**2
    du = xx*mx+th+k*(xx*self%YD(ii)%ta/AA-om)
    dv = y*mx+qx*self%YD(ii)%sa-lam*self%YD(ii)%ca-2.D0*self%YD(ii)%ca*&
        (z*self%YD(ii)%ca+omts*y*self%YD(ii)%sa)/BB+ &
        k*(-1.D0+self%YD(ii)%ca-alA+y*self%YD(ii)%ta/AA+self%YD(ii)%ca*alB)
    dw = z*mx+qx*self%YD(ii)%ca-lam*self%YD(ii)%sa-2.D0*etap*self%YD(ii)%ca/BBp+&
        4.D0*self%YD(ii)%ca*(oms*y*self%YD(ii)%ca-omts*z*self%YD(ii)%sa)/BB+ &
        k*self%YD(ii)%ta*(self%YD(ii)%ca-alA+self%YD(ii)%ca*alB)+4.D0*oms*self%YD(ii)%ca*self%YD(ii)%cota
 else
    rr = xx**2+y**2
    du = -4.D0*xx*y*(rr*self%YD(ii)%sig-r**2)/r/AA**2/(r+z)-2.D0*xx*y*(rr+2.D0*z*(r+z))*oms*omts/rr**2+&
            2.D0*oms*(cPi + datan2(y,xx) - datan2(-y,xx))
    dv = 2.D0*z/(r-z)-4.D0*y**2*(self%YD(ii)%sig*rr-r**2)/r/AA**2/(r+z)+2.D0*omts*oms*(-1.D0+(z*(r-z)-y**2)/AA**2-alA)- &
            omts*dlog((r+z)/AA)
    dw = 0.D0     ! not sure if this limit is correct ... Mathematica gives a directedinfinity value for the limit, which might mean that the
    ! original YSH expression in the paper is incorrect for the w component ... this needs to be rederived and verified !!!
  end if

  u = u+sgn*du*Dx
  v = v+dv*Dx
  w = w+dw*Dx
end if

! and return the displacement components
res = (/ u,v,w /)

end function YSHDisp_

!--------------------------------------------------------------------------
!
! FUNCTION: makeYSHdislocation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief pre-compute geometrical parametersf or the Yoffe&Shaibani&Hazzledine (YSH)
!> surface-relaxed dislocation in an elastically isotropic matrix.
!
!> @details These parameters are then used in the CalcR routine.
!>
!> We implemented the YSH expressions instead of Yoffe's
!> since the former are more easily handled for numerical computations.
!>
!> SH have redefined the x-y-z reference frame used by Yoffe to fall along
!> the dislocation line itself.  As a result, the Burgers vector must be decomposed
!> into a screw component and two edge components, one in the plane of the
!> discontinuity, the other normal to that plane (which is by definition the x-axis).
!> Check the SH paper for more details.
!
!> @param cell unit cell pointer
!> @param defects defects structure
!> @param i dislocation number
!> @param dinfo triggers verbose output
!> @param L column edge length
!
!> @todo Convert IO to Write_Value calls
!
!> @date  1/5/99  MDG 1.0 original
!> @date  5/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 06/04/13 MDG 3.0 rewrite+added quaternions
!> @date 11/21/13 MDG 3.1 verification + rewrite of output handling
!> @date 06/09/14 MDG 4.0 added cell and defects arguments
!> @date 06/10/14 MDG 4.1 added foil argument
!--------------------------------------------------------------------------
recursive subroutine makeYSHdislocation_(self,cell,i,dinfo, L)
!DEC$ ATTRIBUTES DLLEXPORT :: makeYSHdislocation_

use mod_crystallography
use mod_io
use mod_rotations

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)       :: self
type(Cell_T)                        :: cell
type(IO_T)                          :: Message

integer(kind=irg),INTENT(IN)        :: i
integer(kind=irg),INTENT(IN)        :: dinfo
real(kind=sgl),INTENT(IN)           :: L

type(orientation_T)                 :: ot
type(o_T)                           :: a_di_om
type(q_T)                           :: a_di_qu

real(kind=dbl)                       :: alpha, beta, tu(3), tx(3), ty(3), te(3), tb(3), bl, fx(3), fy(3), fz(3), &
                                       dx, dy, a_di(3,3), io_real(3)

! first, determine the alpha angle between the
! negative z-axis, which is really the negative foil normal, and the line direction
! (make sure to reduce the angle to [0,90] interval).
! Each YSH dislocation must have a line direction that points INTO the foil !!!
alpha = cell%CalcAngle(self%foil%F,dble(self%YD(i)%u),'d')*180.0/cPi
if (alpha.ge.90.0) then
  alpha = 180.0-alpha
  if (dinfo.eq.1) then
    io_real(1:3) = self%foil%F(1:3)
    call Message%WriteValue('Foil normal = ', io_real, 3, "('[',3F5.1,']')")
    io_real(1:3) = self%YD(i)%u(1:3)
    call Message%WriteValue('line direction = ', io_real, 3, "('[',3F5.1,']')")
    io_real(1) = alpha
    call Message%WriteValue(' --> alpha angle = ', io_real, 1, "(F5.1)")
  end if
  alpha = alpha*cPi/180.0
else
  call Message%printError('makeYSHdislocation','YSH dislocations must have line directions pointing into the foil ! ')
end if

! normalize the line direction
call cell%TransSpace(self%YD(i)%u,tu,'d','c')
call cell%NormVec(tu,'c')

! consider the case of alpha=0 separately
if (alpha.gt.0.0) then
  call cell%TransSpace(self%foil%F,ty,'d','c')
  call cell%NormVec(ty,'c')                     !  F
  call cell%CalcCross(tu,ty,tx,'c','c',0)       ! x = u x F
  call cell%NormVec(tx,'c')
  call cell%CalcCross(tx,tu,te,'c','c',0)       ! e = x x u
  call cell%NormVec(te,'c')
  call cell%CalcCross(ty,tx,ty,'c','c',0)
  call cell%NormVec(ty,'c')
else
  tx = self%foil%qn
  call cell%CalcCross(tx,tu,te,'c','c',0)       ! e = x x u
  call cell%NormVec(te,'c')
end if
bl = cell%CalcLength(self%YD(i)%burg,'d')

if (dinfo.eq.1) then
  io_real(1:3) = tx(1:3)
  call Message%WriteValue(' tx = ',io_real, 3, "(3F5.1)")
  io_real(1:3) = te(1:3)
  call Message%WriteValue(' te = ',io_real, 3, "(3F5.1)")
  io_real(1:3) = tu(1:3)
  call Message%WriteValue(' tu = ',io_real, 3, "(3F5.1)")
  io_real(1:3) = ty(1:3)
  call Message%WriteValue(' ty = ',io_real, 3, "(3F5.1)")
  io_real(1) = bl
  call Message%WriteValue(' bl = ',io_real, 1, "(F8.3)")
end if

call cell%TransSpace(self%YD(i)%burg,tb,'d','c')
call cell%NormVec(tb,'c')
self%YD(i)%bx = bl * cell%CalcDot(tb,tx,'c')   ! edge component normal to cut plane
self%YD(i)%be = bl * cell%CalcDot(tb,te,'c')   ! edge component in cut plane
self%YD(i)%bs = bl * cell%CalcDot(tb,tu,'c')   ! screw component

if (dinfo.eq.1) then
  io_real(1:3) = (/ self%YD(i)%bx,self%YD(i)%be,self%YD(i)%bs /)
  call Message%WriteValue('Burgers vector components (bx,be,bs) ', io_real, 3, "(3F12.6)")
end if
! verified MDG 7/31/11


! we will also need to know the quaternion rotation between the dislocation reference frame
! and the foil reference frame, so that we can transform the foil coordinates to defect
! coordinates...  We need the angle beta between the defect x axis (tx) and the foil x axis,
! which is the first column of the foil%a_fc matrix ...   We must make sure that this angle
! is measured in a CCW sense.

! projection of defect x axis onto foil x and y axes
call cell%TransSpace(self%foil%q,fx,'d','c')
call cell%TransSpace(self%foil%F,fz,'d','c')
call cell%NormVec(fx,'c')
call cell%NormVec(fz,'c')
call cell%CalcCross(fz,fx,fy,'c','c',0)
dx = cell%CalcDot(tx,fx,'c')
dy = cell%CalcDot(tx,fy,'c')

! use the arctan function to get the angle with correct quadrant computation
self%YD(i)%beta = atan2(dy,dx) !+ cPi*0.5

if (dinfo.eq.1) then
  io_real(1) = dx
  call Message%WriteValue(' dx = ', io_real, 1, "(F8.3)")
  io_real(1) = dy
  call Message%WriteValue(' dy = ', io_real, 1, "(F8.3)")
  io_real(1) = self%YD(i)%beta
  call Message%WriteValue(' beta = ', io_real, 1, "(F8.3)")
end if

! convert to a quaternion
beta = self%YD(i)%beta
a_di(1,1:3) = (/ cos(beta), sin(beta), 0.D0 /)
a_di(2,1:3) = (/ -sin(beta), cos(beta), 0.D0 /)
a_di(3,1:3) = (/ 0.D0, 0.D0, 1.D0 /)

a_di_om = o_T(odinp = a_di)
a_di_qu = a_di_om%oq()
self%YD(i)%a_di = Quaternion_T(qd = a_di_qu%q_copyd())
self%YD(i)%a_id = conjg(self%YD(i)%a_di)
!
! if (dinfo.eq.1) then
!   write (*,*) 'beta = ',beta
!   write (*,*) self%YD(i)%a_di
!   write (*,*) self%YD(i)%a_id
! end if


! finally some geometrical parameters needed for the displacement field computation...
self%YD(i)%alpha =  alpha
self%YD(i)%ca = cos(alpha)
self%YD(i)%sa = sin(alpha)
self%YD(i)%ta = tan(alpha)
self%YD(i)%cota = 1.0/self%YD(i)%ta

! that's it! the rest is handled in the CalcR routine.

end subroutine makeYSHdislocation_

!--------------------------------------------------------------------------
!
! SUBROUTINE: makestackingfault
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute parameters for a stacking fault
!
!> @details  This subroutine computes the geometrical parameters for a
!> stacking fault.  It computes, among others, the coordinates of the
!> centers of the partial dislocations, the intersections of each dislocation
!> line with the top and bottom foil surface, and an array that indicates, for
!> each image pixel, whether or not the corresponding integration column
!> contains this fault; if it does not, the  value in the array is set to
!> -10000; if it does, then the value is equal to the point where the fault
!> plane intersects with the column, measured from the top surface.
!> In short, anything that is needed in the CalcR routine and can be computed
!> ahead of time, is computed here.  The routine also calls the makedislocation
!> routine to create the partials.
!
!> @param cell unit cell pointer
!> @param defects defect structure
!> @param inum
!> @param DF_L column edge length
!> @param nx
!> @param ny
!> @param DF_g
!> @param ndl
!> @param dinfo trigger for verbose output
!
!> @date   11/05/13 MDG 1.0 new attempt to replace faulty original routine
!> @date   11/13/13 MDG 1.1 traced error to problem with transformations in defectmodule
!> @date   11/13/13 MDG 1.2 changed SF normal transformation for zpos array computation (to be tested)
!> @date   06/09/14 MDG 2.0 added defects and cell as arguments
!> @date   06/10/14 MDG 2.1 added foil as argument
!--------------------------------------------------------------------------
recursive subroutine makestackingfault_(self,cell,inum,DF_L,nx,ny,DF_g,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT ::  makestackingfault_

use mod_math
use mod_rotations
use mod_crystallography
use mod_symmetry
use mod_quaternions

IMPLICIT NONE


class(Defect_T),INTENT(INOUT)   :: self
type(Cell_T)                    :: cell

real(kind=sgl),INTENT(IN)       :: DF_L
real(kind=sgl),INTENT(IN)       :: DF_g(3)
integer(kind=irg),INTENT(IN)    :: inum, nx, ny, dinfo

real(kind=sgl)                  :: fpn(3),am(4,4),midpoint(3), ex(3), ey(3),&
                                 lptopi(3),lpboti(3),tptopi(3),tpboti(3),det,A(4), xx(4), yy(4), tmp(3), &
                                 planenormal(3), rzero(3), unita(3)

integer(kind=irg)               :: i,j,info,ipiv,minx,maxx,miny,maxy, ndl

type(Quaternion_T)                      :: a_fc_conjg

ndl = self%numdisl

! we begin by computing the geometry in the foil reference frame, which is the cartesian frame
! for zero sample tilt;  sample tilts are applied once we known the partial dislocation geometry
call cell%TransSpace(self%SF(inum)%plane,tmp,'r','c')
call cell%NormVec(tmp,'c')
planenormal =  tmp


call cell%CalcCross(planenormal, (/ 0.0,0.0,1.0 /),  unita, 'c', 'c', 0)
call cell%NormVec(unita,'c')
fpn = planenormal
if (dinfo.eq.1) write (*,*) ' unita should have zero third component ',unita, fpn

! the fault plane goes through the point rzero in the foil center plane
rzero = (/ self%SF(inum)%id, self%SF(inum)%jd, 0.0 /)

! this leads to the location of the partial dislocation centers
self%SF(inum)%lpr = rzero + 0.5*self%SF(inum)%sep*unita / DF_L
self%SF(inum)%tpr = rzero - 0.5*self%SF(inum)%sep*unita / DF_L

if (dinfo.eq.1) write (*,*) 'lpr_i = ',self%SF(inum)%lpr(1:3)
if (dinfo.eq.1) write (*,*) 'tpr_i = ',self%SF(inum)%tpr(1:3)

! call makedislocation for each of the partials
  self%DL(ndl+1)%id = self%SF(inum)%lpr(1)
  self%DL(ndl+1)%jd = self%SF(inum)%lpr(2)
  self%DL(ndl+1)%u =  self%SF(inum)%lpu
  self%DL(ndl+1)%burg = self%SF(inum)%lpb
  self%DL(ndl+1)%g = DF_g
  call self%makedislocation(cell,ndl+1,dinfo,DF_L)
  if (dinfo.eq.1) write (*,*) 'Leading Partial Position ',self%DL(ndl+1)%id,self%DL(ndl+1)%jd

  self%DL(ndl+2)%id = self%SF(inum)%tpr(1)
  self%DL(ndl+2)%jd = self%SF(inum)%tpr(2)
  self%DL(ndl+2)%u = self%SF(inum)%tpu
  self%DL(ndl+2)%burg = self%SF(inum)%tpb
  self%DL(ndl+2)%g = DF_g
  call self%makedislocation(cell,ndl+2,dinfo,DF_L)
  if (dinfo.eq.1)  write (*,*) 'Trailing Partial Position ',self%DL(ndl+2)%id,self%DL(ndl+2)%jd

! copy the top and bottom dislocation intersections (computed in make_dislocation)
! into the corresponding variables of the SF record

self%SF(inum)%lpbot = self%DL(ndl+1)%bottom
self%SF(inum)%lptop = self%DL(ndl+1)%top
self%SF(inum)%tpbot = self%DL(ndl+2)%bottom
self%SF(inum)%tptop = self%DL(ndl+2)%top

! obviously, these four points need to lie in a single plane; at this point, we check that this is indeed the case
! by computing the volume of the tetrahedron formed by these four points; if the volume is zero, then the
! points are co-planar.  (Use LAPACK's LU-decomposition and compute the product of the diagonal elements of U)
am(1:4,1) = (/ self%SF(inum)%lptop(1:3),1.0 /)
am(1:4,2) = (/ self%SF(inum)%lpbot(1:3),1.0 /)
am(1:4,3) = (/ self%SF(inum)%tptop(1:3),1.0 /)
am(1:4,4) = (/ self%SF(inum)%tpbot(1:3),1.0 /)
call sgetrf(4,4,am,4,ipiv,info)
det = abs(am(1,1)*am(2,2)*am(3,3)*am(4,4))
if (dinfo.eq.1) write (*,*) 'determinant (should be zero) = ',det

! ok, next we need to figure out which image pixels lie on the projection of the stacking fault plane.
! We need to transform the corner points into the image reference frame !!!
lptopi = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%lptop))
lpboti = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%lpbot))
tptopi = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%tptop))
tpboti = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%tpbot))

if (dinfo.eq.1) then
  write (*,*) 'SF parameters :'
  write (*,*) lptopi,' <> ',lpboti
  write (*,*) tptopi,' <> ',tpboti
end if

! define the array that will contain the zpos values
allocate(self%SF(inum)%zpos(nx,ny))
self%SF(inum)%zpos = -10000.0    ! set all points to be outside the projected SF

! first determine the smaller box
minx = nint(min( lptopi(1),lpboti(1),tptopi(1),tpboti(1) )) -2
maxx = nint(max( lptopi(1),lpboti(1),tptopi(1),tpboti(1) )) +2
miny = nint(min( lptopi(2),lpboti(2),tptopi(2),tpboti(2) )) -2
maxy = nint(max( lptopi(2),lpboti(2),tptopi(2),tpboti(2) )) +2

! the fault edges may fall outside of the viewing frame (origin at the center !!!)
if (minx.lt.(-nx/2+1)) minx=-nx/2+1
if (maxx.gt.nx/2) maxx=nx/2
if (miny.lt.(-ny/2+1)) miny=-ny/2+1
if (maxy.gt.ny/2) maxy=ny/2

if (dinfo.eq.1) write (*,*) 'Integer fault box = ',minx,maxx,miny,maxy

! get the equation of the stacking fault plane in the image reference frame
! first the unit plane normal in image space
! we'll take two vectors: ex = from ltop to ttop; ey = from ltop to lbot
ex = tptopi - lptopi
ey = lpboti - lptopi
call cell%NormVec(ex,'c')
call cell%NormVec(ey,'c')
call cell%CalcCross(ex,ey,fpn,'c','c',0)

A(1:3) = fpn ! quat_LPstar( conjg(foil%a_fi), dble(fpn(1:3)) )
midpoint = 0.25*(lptopi+lpboti+tptopi+tpboti) ! quat_LPstar( conjg(foil%a_fi), dble(0.25*(lptopi+lpboti+tptopi+tpboti) ))
A(4) = sum(A(1:3)*midpoint)
if (dinfo.eq.1) write (*,*) 'fault plane parameters : ',A, midpoint

! rank the corner points so that the polygon is convex
! call rank_points(tpboti(1:2),lpboti(1:2),lptopi(1:2),tptopi(1:2),xx,yy)
xx = (/ lptopi(1), tptopi(1), tpboti(1),lpboti(1) /)
yy = (/ lptopi(2), tptopi(2), tpboti(2),lpboti(2) /)

! for all of the points inside this box:
do i=minx,maxx
  do j=miny,maxy
    if (point_inside_polygon( float(i), float(j), xx, yy ).gt.0) then
! the point lies inside the projected region,
! so we need the depth of the SF plane at this position, taking into account the
! proper coordinate transformation (depth must be expressed in image reference frame)
        self%SF(inum)%zpos(i+nx/2,j+ny/2) =  DF_L * ( A(4) - A(1)*float(i) - A(2)*float(j) )/A(3)
    end if
  end do
end do
if (dinfo.eq.1) write (*,*) 'fault plane pixels determined'

! let's also make sure that the SF displacement vector is translated to the
! cartesian reference frame, so that it can be used directly by the CalcR routine
self%SF(inum)%lpbc = matmul(cell%getdsm(),self%SF(inum)%Rdisp)

! that should do it for the stacking fault...  The rest
! takes place in the CalcR routine.
end subroutine makestackingfault_

!--------------------------------------------------------------------------
!
! SUBROUTINE: makestackingfaultECCI
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute parameters for a stacking fault in ECCI mode
!
!> @details  This subroutine computes the geometrical parameters for a
!> stacking fault.  It computes, among others, the coordinates of the surface
!> intersections of the partial dislocations, and an array that indicates, for
!> each image pixel, whether or not the corresponding integration column
!> contains this fault; if it does not, the  value in the array is set to
!> -10000; if it does, then the value is equal to the point where the fault
!> plane intersects with the column, measured from the top surface.
!> In short, anything that is needed in the CalcR routine and can be computed
!> ahead of time, is computed here.  The routine also calls the makedislocation
!> routine to create the partials.
!
!> @param cell unit cell pointer
!> @param defects defect structure
!> @param inum
!> @param DF_L column edge length
!> @param nx
!> @param ny
!> @param DF_g
!> @param ndl
!> @param dinfo trigger for verbose output
!
!> @date   11/05/13 MDG 1.0 new attempt to replace faulty original routine
!> @date   11/13/13 MDG 1.1 traced error to problem with transformations in defectmodule
!> @date   11/13/13 MDG 1.2 changed SF normal transformation for zpos array computation (to be tested)
!> @date   12/17/13 MDG 1.3 branch from original routine to deal with different ECCI geometry
!> @date   12/18/13 MDG 1.4 debug of stacking fault location array
!> @date   06/09/14 MDG 2.0 added cell, SF, YD arguments
!> @date   06/10/14 MDG 2.1 added foil argument
!--------------------------------------------------------------------------
recursive subroutine makestackingfaultECCI_(self,cell,inum,DF_L,nx,ny,DF_g,dinfo)
!DEC$ ATTRIBUTES DLLEXPORT :: makestackingfaultECCI_

use mod_math
use mod_rotations
use mod_crystallography
use mod_symmetry
use mod_quaternions

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)   :: self
type(Cell_T)                    :: cell

real(kind=sgl),INTENT(IN)       :: DF_L
real(kind=sgl),INTENT(IN)       :: DF_g(3)
integer(kind=irg),INTENT(IN)    :: inum, nx, ny, dinfo

real(kind=sgl)                  :: fpn(3),am(4,4),midpoint(3), ex(3), ey(3),&
                                 lptopi(3),lpboti(3),tptopi(3),tpboti(3),det,A(4), xx(4), yy(4), tmp(3), &
                                 planenormal(3), rzero(3), unita(3), lun(3), tun(3), zang, zu, zz, zpo

integer(kind=irg)               :: i,j,info,ipiv,minx,maxx,miny,maxy, ndl

type(Quaternion_T)                      :: a_fc_conjg

associate (ndl => self%numYdisl)

! we begin by computing the geometry in the foil reference frame, which is the cartesian frame
! for zero sample tilt;  sample tilts are applied once we known the partial dislocation geometry
call cell%TransSpace(self%SF(inum)%plane,tmp,'r','c')
call cell%NormVec(tmp,'c')
planenormal =  tmp

call cell%CalcCross(planenormal, (/ 0.0,0.0,1.0 /),  unita, 'c', 'c', 0)
call cell%NormVec(unita,'c')
fpn = planenormal
if (dinfo.eq.1) write (*,*) ' unita should have zero third component ',unita, fpn


! the fault plane goes through the point rzero in the foil top surface
rzero = (/ self%SF(inum)%id, self%SF(inum)%jd, sngl(self%foil%z0)*0.5 /)

! this leads to the location of the partial dislocation intersections, and these must
! be Yoffe dislocations !!!
self%SF(inum)%lpr = rzero + 0.5*self%SF(inum)%sep*unita / DF_L
self%SF(inum)%tpr = rzero - 0.5*self%SF(inum)%sep*unita / DF_L

if (dinfo.eq.1) write (*,*) 'lpr_i = ',self%SF(inum)%lpr(1:3)
if (dinfo.eq.1) write (*,*) 'tpr_i = ',self%SF(inum)%tpr(1:3)


! convert line directions to the Cartesian crystal reference frame
call cell%TransSpace(self%SF(inum)%lpu,lun,'d','c')
call cell%TransSpace(self%SF(inum)%tpu,tun,'d','c')
! normalize both vectors
call cell%NormVec(lun,'c')
call cell%NormVec(tun,'c')


! call makeYSHdislocation for each of the partials [must be done in main program]
!  if (.not. allocated(YD)) allocate(self%YD(3*maxself))

  self%YD(ndl+1)%id = self%SF(inum)%lpr(1)
  self%YD(ndl+1)%jd = self%SF(inum)%lpr(2)
  self%YD(ndl+1)%u(1:3) = dble(self%SF(inum)%lpu(1:3))
  self%YD(ndl+1)%burg(1:3) = dble(self%SF(inum)%lpb(1:3))
  self%YD(ndl+1)%g = DF_g
  self%YD(ndl+1)%sig = self%SF(inum)%poisson

  self%YD(ndl+2)%id = self%SF(inum)%tpr(1)
  self%YD(ndl+2)%jd = self%SF(inum)%tpr(2)
  self%YD(ndl+2)%u(1:3) = dble(self%SF(inum)%tpu(1:3))
  self%YD(ndl+2)%burg(1:3) = dble(self%SF(inum)%tpb(1:3))
  self%YD(ndl+2)%g = DF_g
  self%YD(ndl+2)%sig = self%SF(inum)%poisson

  call self%makeYSHdislocation(cell,ndl+1,dinfo,DF_L)
  if (dinfo.eq.1) write (*,*) 'Leading Partial Position ',self%YD(ndl+1)%id,self%YD(ndl+1)%jd
  call self%makeYSHdislocation(cell,ndl+2,dinfo,DF_L)
  if (dinfo.eq.1) write (*,*) 'Trailing Partial Position ',self%YD(ndl+2)%id,self%YD(ndl+2)%jd

! first find the length of the dislocation line inside the foil along the
! dislocation z-axis, which is the line direction; also, compute the intersection
! points of the line with the top and bottom surfaces  (all components in [nm])
zang = cell%CalcAngle(self%YD(ndl+1)%u,self%foil%F,'d')
zz = cos(zang)
if (abs(zz).gt.0.00001) then
  zu = abs(self%foil%z0/zz)
else
  zu = 100000.0         ! this is when the dislocation is nearly parallel to the foil
end if

! transform the line direction to the foil reference frame
a_fc_conjg = conjg(self%foil%a_fc)
tmp = a_fc_conjg%quat_Lp(dble(lun) ) / DF_L

! determine the top and bottom intersection coordinates
  self%YD(ndl+1)%top = (/ self%YD(ndl+1)%id, self%YD(ndl+1)%jd, self%foil%z0*0.5 /)
  if (zz.gt.0.0) then  ! u points to the top of the foil
    self%YD(ndl+1)%bottom = (/ self%YD(ndl+1)%id - tmp(1)*zu, self%YD(ndl+1)%jd - tmp(2)*zu, -0.5D0*self%foil%z0 /)
  else                 ! u points to the bottom of the foil
    self%YD(ndl+1)%bottom = (/ self%YD(ndl+1)%id + tmp(1)*zu, self%YD(ndl+1)%jd + tmp(2)*zu, -0.5D0*self%foil%z0 /)
  end if


! first find the length of the dislocation line inside the foil along the
! dislocation z-axis, which is the line direction; also, compute the intersection
! points of the line with the top and bottom surfaces  (all components in [nm])
zang = cell%CalcAngle(self%YD(ndl+2)%u,self%foil%F,'d')
zz = cos(zang)
if (abs(zz).gt.0.00001) then
  zu = abs(self%foil%z0/zz)
else
  zu = 100000.0         ! this is when the dislocation is nearly parallel to the foil
end if
! transform the line direction to the foil reference frame
  tmp = a_fc_conjg%quat_Lp(dble(tun) ) / DF_L

! determine the top and bottom intersection coordinates
  self%YD(ndl+2)%top = (/ self%YD(ndl+2)%id, self%YD(ndl+2)%jd, self%foil%z0*0.5 /)
  if (zz.gt.0.0) then  ! u points to the top of the foil
    self%YD(ndl+2)%bottom = (/ self%YD(ndl+2)%id - tmp(1)*zu, self%YD(ndl+2)%jd - tmp(2)*zu, -0.5D0*self%foil%z0 /)
  else                 ! u points to the bottom of the foil
    self%YD(ndl+2)%bottom = (/ self%YD(ndl+2)%id + tmp(1)*zu, self%YD(ndl+2)%jd + tmp(2)*zu, -0.5D0*self%foil%z0 /)
  end if


! copy the top and bottom dislocation intersections
! into the corresponding variables of the SF record

self%SF(inum)%lpbot = self%YD(ndl+1)%bottom
self%SF(inum)%lptop = self%YD(ndl+1)%top
self%SF(inum)%tpbot = self%YD(ndl+2)%bottom
self%SF(inum)%tptop = self%YD(ndl+2)%top

! obviously, these four points need to lie in a single plane; at this point, we check that this is indeed the case
! by computing the volume of the tetrahedron formed by these four points; if the volume is zero, then the
! points are co-planar.  (Use LAPACK's LU-decomposition and compute the product of the diagonal elements of U)
am(1:4,1) = (/ self%SF(inum)%lptop(1:3),1.0 /)
am(1:4,2) = (/ self%SF(inum)%lpbot(1:3),1.0 /)
am(1:4,3) = (/ self%SF(inum)%tptop(1:3),1.0 /)
am(1:4,4) = (/ self%SF(inum)%tpbot(1:3),1.0 /)
call sgetrf(4,4,am,4,ipiv,info)
det = abs(am(1,1)*am(2,2)*am(3,3)*am(4,4))
if (dinfo.eq.1) write (*,*) 'determinant (should be zero) = ',det

! ok, next we need to figure out which image pixels lie on the projection of the stacking fault plane.
! We need to transform the corner points into the image reference frame !!!

lptopi = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%lptop))
lpboti = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%lpbot))
tptopi = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%tptop))
tpboti = self%foil%a_fi%quat_Lp(dble(self%SF(inum)%tpbot))
if (dinfo.eq.1) then
  write (*,*) 'SF parameters :'
  write (*,*) lptopi,' <> ',lpboti
  write (*,*) tptopi,' <> ',tpboti
end if

! define the array that will contain the zpos values

!allocate(self%SF(inum)%zpos(nx,ny))


self%SF(inum)%zpos = -10000.0    ! set all points to be outside the projected SF

! first determine the smaller box
minx = nint(min( lptopi(1),lpboti(1),tptopi(1),tpboti(1) )) -2
maxx = nint(max( lptopi(1),lpboti(1),tptopi(1),tpboti(1) )) +2
miny = nint(min( lptopi(2),lpboti(2),tptopi(2),tpboti(2) )) -2
maxy = nint(max( lptopi(2),lpboti(2),tptopi(2),tpboti(2) )) +2

! the fault edges may fall outside of the viewing frame (origin at the center !!!)
if (minx.lt.(-nx/2+1)) minx=-nx/2+1
if (maxx.gt.nx/2) maxx=nx/2
if (miny.lt.(-ny/2+1)) miny=-ny/2+1
if (maxy.gt.ny/2) maxy=ny/2

if (dinfo.eq.1) write (*,*) 'Integer fault box = ',minx,maxx,miny,maxy

! get the equation of the stacking fault plane in the image reference frame
! first the unit plane normal in image space
! we'll take two vectors: ex = from ltop to ttop; ey = from ltop to lbot
ex = tptopi - lptopi
ey = lpboti - lptopi
call cell%NormVec(ex,'c')
call cell%NormVec(ey,'c')
call cell%CalcCross(ex,ey,fpn,'c','c',0)

A(1:3) = fpn ! quat_LPstar( conjg(foil%a_fi), dble(fpn(1:3)) )
midpoint = 0.25*(lptopi+lpboti+tptopi+tpboti) ! quat_LPstar( conjg(foil%a_fi), dble(0.25*(lptopi+lpboti+tptopi+tpboti) ))
A(4) = sum(A(1:3)*midpoint)
if (dinfo.eq.1) write (*,*) 'fault plane parameters : ',A, midpoint

! rank the corner points so that the polygon is convex
! call rank_points(tpboti(1:2),lpboti(1:2),lptopi(1:2),tptopi(1:2),xx,yy)
xx = (/ lptopi(1), tptopi(1), tpboti(1),lpboti(1) /)
yy = (/ lptopi(2), tptopi(2), tpboti(2),lpboti(2) /)

! for all of the points inside this box:
do i=minx,maxx
  do j=miny,maxy
    if (point_inside_polygon( float(i), float(j), xx, yy ).gt.0) then
! the point lies inside the projected region,
! so we need the depth of the SF plane at this position, taking into account the
! proper coordinate transformation (depth must be expressed in image reference frame)
        self%SF(inum)%zpos(i+nx/2,j+ny/2) = ( A(4) - A(1)*float(i) - A(2)*float(j) )/A(3)
    end if
  end do
end do
if (dinfo.eq.1) write (*,*) 'fault plane pixels determined'

! let's also make sure that the SF displacement vector is translated to the
! cartesian reference frame, so that it can be used directly by the CalcR routine
self%SF(inum)%lpbc = matmul(cell%getdsm(),self%SF(inum)%Rdisp)

! that should do it for the stacking fault...  The rest
! takes place in the CalcR routine.

end associate

end subroutine makestackingfaultECCI_

!--------------------------------------------------------------------------
!
! FUNCTION: Eshelby_disp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the actual displacement vector for an isotropic Eshelby ellipsoidal inclusion
!
!> We implemented the Eshelby expressions based on Mura's 1987 book.
!
!> @param defects defects structure
!> @param i Einclusion number
!> @param xyz coordinate triplet
!
!> @date 12/11/15 MDG 1.0 initial version based on trial IDL script
!> @date 12/13/15 MDG 1.1 corrections to some of the auxiliary expressions; sphere limit is now correct
!--------------------------------------------------------------------------
recursive function Eshelby_disp_(self, i, xyz) result(u)
!DEC$ ATTRIBUTES DLLEXPORT :: Eshelby_disp_
!
! implement the displacement field equations for an isotropic ellipsoidal inclusion
!

use mod_math

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)      :: self
integer(kind=irg),INTENT(IN)    :: i
real(kind=dbl),INTENT(IN)       :: xyz(3)
real(kind=dbl)                  :: u(3)

integer(kind=irg)               :: itheta, j, l, k
real(kind=dbl)                  :: Treps, r2(3), rsq, dd, v, w, modt, lambda, xx, thetaEl, dtheta, sigma, &
                                   II1, II3, II(3), IIJ(3,3), phici(3), psicjli(3,3,3), c, t1, t2, dt(3), Delta, EE, EF
complex(kind=dbl)               :: t, s, q

Treps = self%Einclusions(i)%epsstar(1,1)+self%Einclusions(i)%epsstar(2,2)+self%Einclusions(i)%epsstar(3,3)

r2 = xyz**2
rsq = sum(r2)
dd = dsqrt(xyz(1)**2/self%Einclusions(i)%a12+xyz(2)**2/self%Einclusions(i)%a22+xyz(3)**2/self%Einclusions(i)%a32)

!  first we need the lambda value for this point
lambda = 0.D0
if (dd.gt.1.D0) then
  s = cmplx( (rsq - self%Einclusions(i)%eta)**2 +sum(self%Einclusions(i)%svec*r2)-self%Einclusions(i)%ss1, 0.D0 )
  q = cmplx( -2.D0*rsq**3 - 3.D0*self%Einclusions(i)%eta*rsq**2 + rsq*(3.D0*self%Einclusions(i)%eta**2+ &
      sum(self%Einclusions(i)%qvec1*r2))-sum(self%Einclusions(i)%qvec2*r2)-self%Einclusions(i)%qs1, 0.D0 )
  t = (q + sqrt( q*q - 4.D0*s*s*s) )**(1.D0/3.D0)
  v = real(t)
  w = abs(aimag(t))
  modt = v**2+w**2
  lambda = (rsq-self%Einclusions(i)%eta)/3.D0 + (v+self%Einclusions(i)%s3*w)* &
           (2.D0*real(s)+self%Einclusions(i)%c1*modt)/self%Einclusions(i)%c2/modt
end if


! predefine a couple of parameters
dt = (/ self%Einclusions(i)%a12, self%Einclusions(i)%a22, self%Einclusions(i)%a32 /) + lambda
Delta = dsqrt(dt(1)*dt(2)*dt(3))


if (lambda.ne.0.D0) then
! next, we compute the elliptic integrals by interpolating the look-up tables
  thetaEl = dasin(dsqrt(self%Einclusions(i)%Deltaij(1,3)**2/dt(1)))
  xx = self%Einclusions(i)%thpre * (thetaEl-self%Einclusions(i)%mith)
  itheta = int(xx)+1
  dtheta = xx-int(xx)
  if (itheta.lt.self%Einclusions(i)%nLUT) then
    EF = (1.D0-dtheta)*self%Einclusions(i)%EFLUT(itheta) + dtheta*self%Einclusions(i)%EFLUT(itheta+1)
    EE = (1.D0-dtheta)*self%Einclusions(i)%EELUT(itheta) + dtheta*self%Einclusions(i)%EELUT(itheta+1)
  else
    EF = self%Einclusions(i)%EFLUT(itheta)
    EE = self%Einclusions(i)%EELUT(itheta)
  end if

! first order integrals
  II1 = self%Einclusions(i)%preI1 * ( EF - EE )
  II3 = self%Einclusions(i)%preI3 * ( self%Einclusions(i)%Deltaij(1,3)*dt(2)/Delta - EE )
  II = (/ II1, 3.D0*self%Einclusions(i)%V/Delta - II1 - II3, II3 /)
! second order integrals
  IIJ(1,2) = - (II(1) - II(2))/self%Einclusions(i)%Deltaij(1,2)**2
  IIJ(2,1) = IIJ(1,2)
  IIJ(1,3) = - (II(1) - II(3))/self%Einclusions(i)%Deltaij(1,3)**2
  IIJ(3,1) = IIJ(1,3)
  IIJ(2,3) = - (II(2) - II(3))/self%Einclusions(i)%Deltaij(2,3)**2
  IIJ(3,2) = IIJ(2,3)
  IIJ(1,1) = self%Einclusions(i)%V/dt(1)/Delta - (IIJ(1,2)+IIJ(1,3))/3.0
  IIJ(2,2) = self%Einclusions(i)%V/dt(2)/Delta - (IIJ(2,1)+IIJ(2,3))/3.0
  IIJ(3,3) = self%Einclusions(i)%V/dt(3)/Delta - (IIJ(3,1)+IIJ(3,2))/3.0
else
  II = self%Einclusions(i)%IIinside
  IIJ = self%Einclusions(i)%IIJinside
end if

! then we need to evaluate the large number of "tensor" components in the phi_{,j}
! and psi_{,jli} arrays, as defined in Muro's 1987 book; we've actually computed the
! derivatives ourselves similar to the computation of eq. (11.40)
phici = (/ -xyz(1)*II(1), -xyz(2)*II(2), -xyz(3)*II(3) /)

c = 0.D0
sigma = (xyz(1)/dt(1))**2 + (xyz(2)/dt(2))**2 + (xyz(3)/dt(3))**2
do j=1,3
 do l=1,3
  do k=1,3
   if (lambda.ne.0.D0) c = 3.D0*self%Einclusions(i)%V*lambda*xyz(k)*xyz(j)*xyz(l)/dt(k)/dt(j)/dt(l)/Delta/sigma
   t1 = kdelta(j,l) * xyz(k) * (II(k) - self%Einclusions(i)%asq(j)*IIJ(j,k))
   t2 = (kdelta(k,j)*xyz(l)+kdelta(k,l)*xyz(j))*(II(l)-self%Einclusions(i)%asq(j)*IIJ(j,l))
   psicjli(j,l,k) = -t1 - t2 + c
  end do
 end do
end do

do k=1,3
! middle term of Mura eq. (11.30)
  u(k) = -2.0*self%Einclusions(i)%nu*Treps*phici(k)
! last term
  do l=1,3
    u(k) = u(k) - 4.D0*(1.D0-self%Einclusions(i)%nu)*self%Einclusions(i)%epsstar(k,l)*phici(l)
  end do
! first term
  do j=1,3
    do l=1,3
      u(k) = u(k) + self%Einclusions(i)%epsstar(j,l) * psicjli(j,l,k)
    end do
  end do
end do
u = self%Einclusions(i)%pre * u

end function Eshelby_disp_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! and finally, here we compute the total displacements for an integration column
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcR
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief returns the total displacement vector for each slice in a column
!
!> @details Note that the end result MUST be expressed in the cartesian reference frame !
!
!> @param cell unit cell pointer
!> @param defects defect structure
!> @param i integer x coordinate
!> @param j integer y coordinate
!
!> @note This entire routine was thoroughly verified after the quaternion conversion !
!>
!> General comment for those who wish to add other defects...
!>
!> the general procedure to implement a defect displacement field is as follows:
!> - if the defect has its own reference frame, then transform (xpos,ypos,zpos) to
!>   that frame (see the dislocation section below for an example), and then do
!>   the computation of the displacement vector and express it in the cartesian frame.
!>
!> - if the defect uses the foil reference frame (e.g., voids, inclusions), then use tmpf
!>  as the current position vector.

!> @date  10/20/98 MDG 1.0 original
!> @date   5/22/01 MDG 2.0 f90
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  03/26/13 MDG 3.0 updated IO
!> @date  10/30/13 MDG 3.1 debug of coordinate rotations
!> @date  11/13/13 MDG 3.2 finally, the bug has been found!
!> @date  06/09/14 MDG 4.0 introduced defects argument and simplified routine
!> @date  06/10/14 MDG 4.1 added foil argument
!> @date  11/23/15 MDG 4.2 removed foil argument and placed it inside defects
!--------------------------------------------------------------------------
recursive subroutine CalcR_(self,cell,i,j)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcR_

use mod_crystallography
use mod_quaternions
use mod_rotations

IMPLICIT NONE

class(Defect_T),INTENT(INOUT)           :: self
type(Cell_T)                            :: cell

integer(kind=irg),INTENT(IN)            :: i,j

integer(kind=irg)                       :: k, islice, ii
real(kind=dbl)                          :: dis,xpos,ypos,zpos,sumR(3),thick,tmp(3),tmp2(3), &
                                           tmpf(3),u(3),zaamp,zaphase,zar,zai,zr(3),zi(3), &
                                           zt,fx,fy,fz,a_fm(3,3),ar    !,&
!                                nu,x,y,z,zn,t,pre,r1,r2,r3,th,rn

complex(kind=dbl)                       :: za(3)
complex(kind=sgl)                       :: zero
logical                                 :: isvoid

type(o_T)                               :: a_fm_om
type(q_T)                               :: a_fm_qu
type(Quaternion_T)                      :: a_dc_conj, a_di_conj, a_id_conj

! scale the image coordinates with respect to the origin at the center of the image
 xpos = float(i-self%DF_npix/2)*self%DF_L
 ypos = float(j-self%DF_npiy/2)*self%DF_L

! determine the starting point of the z-integration for the tilted foil
! this depends on the foil normal components which give the equation
! of the top foil plane as F . r = z0/2, from which we get zt...
 a_fm_qu = q_T( qdinp =  self%foil%a_fm%get_quatd() )
 a_fm_om = a_fm_qu%qo()

 a_fm =  a_fm_om%o_copyd()
 fx = a_fm(3,1)
 fy = a_fm(3,2)
 fz = a_fm(3,3)
 zt = self%foil%zb*0.5 - (fx*xpos + fy*ypos)/fz

! initialize some other variables
 thick = self%foil%zb
 zero = cmplx(0.0,0.0)

! loop over all slices (this is the main loop)
 sliceloop: do islice = 1,self%DF_nums

! zpos is the position down the column, starting at zt (in image coordinates)
    zpos = zt - float(islice)*self%DF_slice

! set the displacements to zero
    sumR = 0.0

! set the position in the foil reference frame
    tmpf = self%foil%a_fi%quat_Lp(dble( (/ xpos, ypos, zpos /)))


!------------
!----VOIDS---
!------------
! voids are easy to deal with; we simply return -10000 for each point tmpf that lies inside
! one of the voids; the calling routine then knows to use the void scattering matrix.
if (self%numvoids.ne.0) then
! are we inside a void ?
    isvoid = .FALSE.
    voidloop: do ii=1,self%numvoids
! subtract the void position from the current slice position to get the relative position vector
     tmp = tmpf -  (/ self%voids(ii)%xpos, self%voids(ii)%ypos, self%voids(ii)%zpos /)
     dis = cell%CalcLength(tmp,'c')
     if (dis.lt.self%voids(ii)%radius) then ! inside void
       isvoid = .TRUE.
       exit voidloop
     end if
    end do voidloop
! skip the rest of the computation for this slice if we are inside a void
    if (isvoid.eqv..TRUE.) then
      self%DF_R(islice,1) = -10000.0
      cycle sliceloop
    end if
 end if


! ok, if we get here, then we're not inside a void...

!------------------
!----CURVED FOIL---
!------------------
! first we take the foil shape into account using equations (8.28) and (8.29)
 sumR = sumR + float(islice)*self%DF_slice*self%foil%sg(i,j)*self%DF_gstar

!-----------------
!--DISLOCATIONS--
!-----------------
! let's put a few dislocations in ... (see section 8.4.2)
do ii=1,self%numdisl
! compute the difference vector between the current (xpos,ypos,zpos) in the foil reference frame
! and the defect center coordinate
  tmp2 =  tmpf - dble((/ self%DF_L*self%DL(ii)%id, self%DF_L*self%DL(ii)%jd, self%DL(ii)%zfrac*self%foil%z0 /))

! then convert the difference vector to the defect reference frame for this dislocation (we will only need the x and y coordinates)
  tmp = self%DL(ii)%a_df%quat_Lp(tmp2)


! compute x1 + p_alpha x2  (eq. 8.38)
  za(1:3) = tmp(1) + self%DL(ii)%pa(1:3)*tmp(2)
! compute the displacement vector u (eq. 8.38) [this expands the log of a complex number and takes the real part only,
! taking proper care of the branch cut]
   if (tmp(1).gt.0.0) then
   do k=1,3
    zar =  real(za(k))
    zai = aimag(za(k))
    zaamp = abs(za(k))
    zaphase = abs(zai/zar)
    zr(k) = log(zaamp)
    zi(k) = atan(zaphase)
    if (zar.le.0.0) then
      if (zai.lt.0.0) zi(k) = -cPi+zi(k)
      if (zai.eq.0.0) zi(k) = cPi
      if (zai.gt.0.0) zi(k) = cPi-zi(k)
    else
      if (zai.lt.0.0) zi(k) = -zi(k)
    end if
   end do
  else
   do k=1,3
    zar =  real(za(k))
    zai = aimag(za(k))
    zaamp = abs(za(k))
    zaphase = abs(zai/zar)
    zr(k) = log(zaamp)
    zi(k) = atan(zaphase)
    if (zar.le.0.0) then
      if (zai.gt.0.0) zi(k) = cPi-zi(k)
      if (zai.eq.0.0) zi(k) = cPi
      if (zai.lt.0.0) zi(k) = cPi+zi(k)
    else
      if (zai.lt.0.0) zi(k) = 2.0*cPi-zi(k)
      if (zai.eq.0.0) zi(k) = 0.0
    end if
   end do
  end if
  u = 2.0*real(matmul(self%DL(ii)%dismat,cmplx(zr,zi)))
! transform displacement vector u to the Cartesian crystal reference frame
  a_dc_conj = conjg(self%DL(ii)%a_dc)
  u = a_dc_conj%quat_Lp(dble(u) )
  sumR = sumR + u

end do
!
!-------------------------------------
!--SURFACE INTERSECTING DISLOCATIONS--
!-------------------------------------
! this part is mostly used for ECCI-type image simulations, not for EM or STEM,
! although it could probably be used there as well; we would need to extend it
! to incorporate both top and bottom foil surfaces

! do we have any dislocations with surface relaxations ?  YSH model
if (self%numYdisl.gt.0) then
   do ii=1,self%numYdisl
! first, figure out what the coordinates are in the YSH reference frame for this dislocation ...
! translate to the defect origin
     tmp =  tmpf -  (/ self%DF_L*self%YD(ii)%id, self%DF_L*self%YD(ii)%jd, self%foil%z0*0.5 /)

! rotate into the defect reference frame
     a_di_conj = conjg(self%YD(ii)%a_di)
     tmp = a_di_conj%quat_Lp(tmp)
! compute the displacement vector
!     u = sngl(YSHDisp(dble(tmp(2)),-dble(tmp(1)),dble(tmp(3)),ii))
     u = sngl(self%YSHDisp(dble(tmp(1)),dble(tmp(2)),dble(tmp(3)),ii))

! and rotate back to the image reference frame
     a_id_conj = conjg(self%YD(ii)%a_id)
     u = a_id_conj%quat_Lp(u)
     u = self%foil%a_ic%quat_Lp(u)

! that should do it !
     sumR = sumR + u

   end do
end if

!--------------------
!--STACKING FAULTS--
!--------------------
! stacking faults (this is easy because we've already done all the work in the stacking_fault module)
! all we need is the z-value at which the stacking fault plane is crossed in this particular image
! column; from that point on, we simply add the leading partial Burgers vector to the total displacement.
do ii=1,self%numsf
  if ((zpos.lt.self%SF(ii)%zpos(i,j)).and.(self%SF(ii)%zpos(i,j).ne.-10000.0)) then
    sumR = sumR + self%SF(ii)%lpbc
  end if
end do


! !--------------------
! !--LARGE INCLUSIONS--  currently commented out
! !--------------------
! ! Mader's expression for the displacement field of a large inclusion
! !   if (0.eq.1.) then
! !    nu = 0.25
! !    ce = 0.005
! !    rn = 25.0*DF_L
! !    x = (float(i-DF_npix/2)-0.5)*DF_L
! !    y = (float(j-DF_npiy/2)-0.5)*DF_L
! !    z = float(k)*DF_slice
! !    zn = 100.5*DF_slice
! !    t = DF_slice * DF_nums
! !    pre = (1.0+nu)/(3.0*(1.0-nu))*ce*rn**3
! !
! !    r1 = sqrt(x**2+y**2+(z-zn)**2)
! !    r2 = sqrt(x**2+y**2+(z+zn)**2)
! !    r3 = sqrt(x**2+y**2+(2.0*t-z-zn)**2)
! !
! !    if (((r1.eq.0.0).or.(r2.eq.0.0)).or.(r3.eq.0.0)) then
! !      return
! !    else
! !     dis = (1.0/r1**3+(3.0-4.0*nu)/r2**3-6.0*z*(z+zn)/r2**5+(3.0-4.0*nu)/r3**3-6.0*(t-z)*(2.0*t-z-zn)/r3**5)
! !     rx = x*dis
! !     ry = y*dis
! !     rz = (z-zn)/r1**3-(3.0-4.0*nu)*((z+zn)/r2**3+(2.0*t-z-zn)/r3**3)-6.0*z*(z+zn)**2/r2**5 + &
! !          2.0*z/r2**3+6.0*(t-z)*(2.0*t-z-zn)**2/r3**5-2.0*(t-z)/r3**3
! !
! !     sumR = pre*(/ rx, ry, rz /)
! !     return
! !    end if
! !   end if
!
!--------------------
!--SMALL INCLUSIONS--
!--------------------
! then the coherent precipitates, using the model in section 8.4.1
  if (self%numinc.gt.0) then
   do ii=1,self%numinc
! subtract the inclusion position from the current slice position to get the relative position vector
     tmp = tmpf - (/ self%inclusions(ii)%xpos, self%inclusions(ii)%ypos, self%inclusions(ii)%zpos /)
     dis = cell%CalcLength(tmp,'c')
     if (dis.ge.self%inclusions(ii)%radius) then ! outside particle
       tmp = tmp*(self%inclusions(ii)%radius/dis)**3
     end if
     sumR = sumR + self%inclusions(ii)%C*tmp
   end do
  end if

!---------------------
!--SMALL EINCLUSIONS--
!---------------------
! then the coherent ellipsoidally distorted precipitates, using Eshelby's model
   if (self%numEinc.gt.0) then
    do ii=1,self%numEinc
! subtract the inclusion position from the current slice position to get the relative position vector
      tmp = tmpf - (/ self%Einclusions(ii)%xpos, self%Einclusions(ii)%ypos, self%Einclusions(ii)%zpos /)
! and also get the position vector for the mirror image inclusion, to make sure we get a traction-free surface...
      tmp2 = tmpf - (/ self%Einclusions(ii)%xpos, self%Einclusions(ii)%ypos, -self%Einclusions(ii)%zpos /)
      u = self%Eshelby_disp(ii, tmp) + self%Eshelby_disp(ii, tmp2)
! we need to check the reference frame here !
      sumR = sumR + u
    end do
   end if

! TO BE IMPLEMENTED FOR RICHARD LESAR'S Discrete Dislocation Dynamics !
! finally any displacement fields defined by the user routine UserDisp
! sumR = sumR + UserDisp()


! TO BE IMPLEMENTED FOR YUNZHI WANG's Dislocation Simulations !
! finally any displacement fields defined by the user routine UserDisp
! sumR = sumR + UserDisp()

   self%DF_R(islice,1:3) = sumR(1:3)

  end do sliceloop ! main loop over the slices
end subroutine CalcR_


end module mod_defect
