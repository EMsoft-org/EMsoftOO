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

module mod_symmetry
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! all symmetry-related constants, derived types and methods
  !!
  !! routines to generate a space group based on the generator string; computation
  !! of orbits and families; computation of all atoms in a single or multiple unit cells
  !!
  !! original modification history
  !! @date  01/05/99 MDG 1.0 original
  !! @date  05/19/01 MDG 2.0 f90
  !! @date  11/27/01 MDG 2.1 added kind support
  !! @date  03/19/13 MDG 3.0 updated IO and such
  !! @date  01/10/14 MDG 4.0 SG is now part of the unitcell type
  !! @date  06/05/14 MDG 4.1 made cell an argument instead of global variable
  !! @date  09/05/16 MDG 4.2 added Wyckoff Position routines
  !!
  !! this module contains all the routines from the original symmetry.f90 module, as well as
  !! all symmetry-related constants and types from the typedefs.f90 module

use mod_kinds
use mod_global

IMPLICIT NONE
  private

!>  SYM_SGname all space group names
! TRICLINIC SPACE GROUPS
  character(11), public, dimension(237) :: SYM_SGname= (/" P  1      " ," P -1      ", & ! MONOCLINIC SPACE GROUPS
        " P 2       " ," P 21      " ," C 2       " ," P m       ", &
        " P c       " ," C m       " ," C c       " ," P 2/m     ", &
        " P 21/m    " ," C 2/m     " ," P 2/c     " ," P 21/c    ", &
        " C 2/c     ", &                                              ! ORTHORHOMBIC SPACE GROUPS
        " P 2 2 2   " ," P 2 2 21  " ," P 21 21 2 " ," P 21 21 21", &
        " C 2 2 21  " ," C 2 2 2   " ," F 2 2 2   " ," I 2 2 2   ", &
        " I 21 21 21" ," P m m 2   " ," P m c 21  " ," P c c 2   ", &
        " P m a 2   " ," P c a 21  " ," P n c 2   " ," P m n 21  ", &
        " P b a 2   " ," P n a 21  " ," P n n 2   " ," C m m 2   ", &
        " C m c 21  " ," C c c 2   " ," A m m 2   " ," A b m 2   ", &
        " A m a 2   " ," A b a 2   " ," F m m 2   " ," F d d 2   ", &
        " I m m 2   " ," I b a 2   " ," I m a 2   " ," P m m m   ", &
        " P n n n   " ," P c c m   " ," P b a n   " ," P m m a   ", &
        " P n n a   " ," P m n a   " ," P c c a   " ," P b a m   ", &
        " P c c n   " ," P b c m   " ," P n n m   " ," P m m n   ", &
        " P b c n   " ," P b c a   " ," P n m a   " ," C m c m   ", &
        " C m c a   " ," C m m m   " ," C c c m   " ," C m m a   ", &
        " C c c a   " ," F m m m   " ," F d d d   " ," I m m m   ", &
        " I b a m   " ," I b c a   " ," I m m a   ", &                ! TETRAGONAL SPACE GROUPS
        " P 4       " ," P 41      " ," P 42      " ," P 43      ", &
        " I 4       " ," I 41      " ," P -4      " ," I -4      ", &
        " P 4/m     " ," P 42/m    " ," P 4/n     " ," P 42/n    ", &
        " I 4/m     " ," I 41/a    " ," P 4 2 2   " ," P 4 21 2  ", &
        " P 41 2 2  " ," P 41 21 2 " ," P 42 2 2  " ," P 42 21 2 ", &
        " P 43 2 2  " ," P 43 21 2 " ," I 4 2 2   " ," I 41 2 2  ", &
        " P 4 m m   " ," P 4 b m   " ," P 42 c m  " ," P 42 n m  ", &
        " P 4 c c   " ," P 4 n c   " ," P 42 m c  " ," P 42 b c  ", &
        " I 4 m m   " ," I 4 c m   " ," I 41 m d  " ," I 41 c d  ", &
        " P -4 2 m  " ," P -4 2 c  " ," P -4 21 m " ," P -4 21 c ", &
        " P -4 m 2  " ," P -4 c 2  " ," P -4 b 2  " ," P -4 n 2  ", &
        " I -4 m 2  " ," I -4 c 2  " ," I -4 2 m  " ," I -4 2 d  ", &
        " P 4/m m m " ," P 4/m c c " ," P 4/n b m " ," P 4/n n c ", &
        " P 4/m b m " ," P 4/m n c " ," P 4/n m m " ," P 4/n c c ", &
        " P 42/m m c" ," P 42/m c m" ," P 42/n b c" ," P 42/n n m", &
        " P 42/m b c" ," P 42/m n m" ," P 42/n m c" ," P 42/n c m", &
        " I 4/m m m " ," I 4/m c m " ," I 41/a m d" ," I 41/a c d", & ! RHOMBOHEDRAL SPACE GROUPS
        " P 3       " ," P 31      " ," P 32      " ," R 3       ", &
        " P -3      " ," R -3      " ," P 3 1 2   " ," P 3 2 1   ", &
        " P 31 1 2  " ," P 31 2 1  " ," P 32 1 2  " ," P 32 2 1  ", &
        " R 3 2     " ," P 3 m 1   " ," P 3 1 m   " ," P 3 c 1   ", &
        " P 3 1 c   " ," R 3 m     " ," R 3 c     " ," P -3 1 m  ", &
        " P -3 1 c  " ," P -3 m 1  " ," P -3 c 1  " ," R -3 m    ", &
        " R -3 c    ", &                                              ! HEXAGONAL SPACE GROUPS
        " P 6       " ," P 61      " ," P 65      " ," P 62      ", &
        " P 64      " ," P 63      " ," P -6      " ," P 6/m     ", &
        " P 63/m    " ," P 6 2 2   " ," P 61 2 2  " ," P 65 2 2  ", &
        " P 62 2 2  " ," P 64 2 2  " ," P 63 2 2  " ," P 6 m m   ", &
        " P 6 c c   " ," P 63 c m  " ," P 63 m c  " ," P -6 m 2  ", &
        " P -6 c 2  " ," P -6 2 m  " ," P -6 2 c  " ," P 6/m m m ", &
        " P 6/m c c " ," P 63/m c m" ," P 63/m m c", &                ! CUBIC SPACE GROUPS
        " P 2 3     " ," F 2 3     " ," I 2 3     " ," P 21 3    ", &
        " I 21 3    " ," P m 3     " ," P n 3     " ," F m 3     ", &
        " F d 3     " ," I m 3     " ," P a 3     " ," I a 3     ", &
        " P 4 3 2   " ," P 42 3 2  " ," F 4 3 2   " ," F 41 3 2  ", &
        " I 4 3 2   " ," P 43 3 2  " ," P 41 3 2  " ," I 41 3 2  ", &
        " P -4 3 m  " ," F -4 3 m  " ," I -4 3 m  " ," P -4 3 n  ", &
        " F -4 3 c  " ," I -4 3 d  " ," P m 3 m   " ," P n 3 n   ", &
        " P m 3 n   " ," P n 3 m   " ," F m 3 m   " ," F m 3 c   ", &
        " F d 3 m   " ," F d 3 c   " ," I m 3 m   " ," I a 3 d   ", & ! TRIGONAL GROUPS RHOMBOHEDRAL SETTING
        " R 3   |146" ," R -3  |148" ," R 3 2 |155" ," R 3 m |160", &
        " R 3 c |161" ," R -3 m|166" ," R -3 c|167"/)
!DEC$ ATTRIBUTES DLLEXPORT :: SYM_SGname

!> extended Hermann-Mauguin symbols for the orthorhombic space groups in the following settings:
  character(8), public, dimension(6):: extendedOrthsettings = (/ &
    " a  b  c", " b  a -c", " c  a  b", "-c  b  a", " b  c  a", " a -c  b"  /)
!DEC$ ATTRIBUTES DLLEXPORT :: extendedOrthsettings

  character(11), public, dimension(6,59) :: extendedHMOrthsymbols = reshape( (/ &
    " P 2 2 2   ", " P 2 2 2   ", " P 2 2 2   ", " P 2 2 2   ", " P 2 2 2   ", " P 2 2 2   ", &
    " P 2 2 21  ", " P 2 2 21  ", " P 21 2 2  ", " P 21 2 2  ", " P 2 21 2  ", " P 2 21 2  ", &
    " P 21 21 2 ", " P 21 21 2 ", " P 2  21 21", " P 2  21 21", " P 21 2  21", " P 21 2  21", &
    " P 21 21 21", " P 21 21 21", " P 21 21 21", " P 21 21 21", " P 21 21 21", " P 21 21 21", &
    " C 2 2 21  ", " C 2 2 21  ", " A 21 2 2  ", " A 21 2 2  ", " B 2 21 2  ", " B 2 21 2  ", &
    " C 2 2 2   ", " C 2 2 2   ", " A 2 2 2   ", " A 2 2 2   ", " B 2 2 2   ", " B 2 2 2   ", &
    " F 2 2 2   ", " F 2 2 2   ", " F 2 2 2   ", " F 2 2 2   ", " F 2 2 2   ", " F 2 2 2   ", &
    " I 2 2 2   ", " I 2 2 2   ", " I 2 2 2   ", " I 2 2 2   ", " I 2 2 2   ", " I 2 2 2   ", &
    " I 21 21 21", " I 21 21 21", " I 21 21 21", " I 21 21 21", " I 21 21 21", " I 21 21 21", &
    " P m m 2   ", " P m m 2   ", " P 2 m m   ", " P 2 m m   ", " P m 2 m   ", " P m 2 m   ", &
    " P m c 21  ", " P c m 21  ", " P 21 m a  ", " P 21 a m  ", " P b 21 m  ", " P m 21 b  ", &
    " P c c 2   ", " P c c 2   ", " P 2 a a   ", " P 2 a a   ", " P b 2 b   ", " P b 2 b   ", &
    " P m a 2   ", " P b m 2   ", " P 2 m b   ", " P 2 c m   ", " P c 2 m   ", " P c 2 a   ", &
    " P c a 21  ", " P b c 21  ", " P 21 a b  ", " P 21 c a  ", " P c 21 b  ", " P b 21 a  ", &
    " P n c 2   ", " P c n 2   ", " P 2 n a   ", " P 2 a n   ", " P b 2 n   ", " P n 2 b   ", &
    " P m n 21  ", " P n m 21  ", " P 21 m n  ", " P 21 n m  ", " P n 21 m  ", " P m 21 n  ", &
    " P b a 2   ", " P b a 2   ", " P 2 c b   ", " P 2 c b   ", " P c 2 a   ", " P c 2 a   ", &
    " P n a 21  ", " P b n 21  ", " P 21 n b  ", " P 21 c n  ", " P c 21 n  ", " P n 21 a  ", &
    " P n n 2   ", " P n n 2   ", " P 2 n n   ", " P 2 n n   ", " P n 2 n   ", " P n 2 n   ", &
    " C m m 2   ", " C m m 2   ", " A 2 m m   ", " A 2 m m   ", " B m 2 m   ", " B m 2 m   ", &
    " C m c 21  ", " C c m 21  ", " A 21 m a  ", " A 21 a m  ", " B b 21 m  ", " B m 21 b  ", &
    " C c c 2   ", " C c c 2   ", " A 2 a a   ", " A 2 a a   ", " B b 2 b   ", " B b 2 b   ", &
    " A m m 2   ", " B m m 2   ", " B 2 m m   ", " C 2 m m   ", " C m 2 m   ", " A m 2 m   ", &
    " A b m 2   ", " B m a 2   ", " B 2 c m   ", " C 2 m b   ", " C m 2 a   ", " A c 2 m   ", &
    " A m a 2   ", " B b m 2   ", " B 2 m b   ", " C 2 c m   ", " C c 2 m   ", " A m 2 a   ", &
    " A b a 2   ", " B b a 2   ", " B 2 c b   ", " C 2 c b   ", " C c 2 a   ", " A c 2 a   ", &
    " F m m 2   ", " F m m 2   ", " F 2 m m   ", " F 2 m m   ", " F m 2 m   ", " F m 2 m   ", &
    " F d d 2   ", " F d d 2   ", " F 2 d d   ", " F 2 d d   ", " F d 2 d   ", " F d 2 d   ", &
    " I m m 2   ", " I m m 2   ", " I 2 m m   ", " I 2 m m   ", " I m 2 m   ", " I m 2 m   ", &
    " I b a 2   ", " I b a 2   ", " I 2 c b   ", " I 2 c b   ", " I c 2 a   ", " I c 2 a   ", &
    " I m a 2   ", " I b m 2   ", " I 2 m b   ", " I 2 c m   ", " I c 2 m   ", " I m 2 a   ", &
    " P m m m   ", " P m m m   ", " P m m m   ", " P m m m   ", " P m m m   ", " P m m m   ", &
    " P n n n   ", " P n n n   ", " P n n n   ", " P n n n   ", " P n n n   ", " P n n n   ", &
    " P c c m   ", " P c c m   ", " P m a a   ", " P m a a   ", " P b m b   ", " P b m b   ", &
    " P b a n   ", " P b a n   ", " P n c b   ", " P n c b   ", " P c n a   ", " P c n a   ", &
    " P m m a   ", " P m m b   ", " P b m m   ", " P c m m   ", " P m c m   ", " P m a m   ", &
    " P n n a   ", " P n n b   ", " P b n n   ", " P c n n   ", " P n c n   ", " P n a n   ", &
    " P m n a   ", " P n m b   ", " P b m n   ", " P c n m   ", " P n c m   ", " P m a n   ", &
    " P c c a   ", " P c c b   ", " P b a a   ", " P c a a   ", " P b c b   ", " P b a b   ", &
    " P b a m   ", " P b a m   ", " P m c b   ", " P m c b   ", " P c m a   ", " P c m a   ", &
    " P c c n   ", " P c c n   ", " P n a a   ", " P n a a   ", " P b n b   ", " P b n b   ", &
    " P b c m   ", " P c a m   ", " P m c a   ", " P m a b   ", " P b m a   ", " P c m b   ", &
    " P n n m   ", " P n n m   ", " P m n n   ", " P m n n   ", " P n m n   ", " P n m n   ", &
    " P m m n   ", " P m m n   ", " P n m m   ", " P n m m   ", " P m n m   ", " P m n m   ", &
    " P b c n   ", " P c a n   ", " P n c a   ", " P n a b   ", " P b n a   ", " P c n b   ", &
    " P b c a   ", " P c a b   ", " P b c a   ", " P c a b   ", " P b c a   ", " P c a b   ", &
    " P n m a   ", " P m n b   ", " P b n m   ", " P c m n   ", " P m c n   ", " P n a m   ", &
    " C m c m   ", " C c m m   ", " A m m a   ", " A m a m   ", " B b m m   ", " B m m b   ", &
    " C m c a   ", " C c m b   ", " A b m a   ", " A c a m   ", " B b c m   ", " B m a b   ", &
    " C m m m   ", " C m m m   ", " A m m m   ", " A m m m   ", " B m m m   ", " B m m m   ", &
    " C c c m   ", " C c c m   ", " A m a a   ", " A m m a   ", " B b m b   ", " B b m b   ", &
    " C m m a   ", " C m m b   ", " A b m m   ", " A c m m   ", " B m c m   ", " B m a m   ", &
    " C c c a   ", " C c c b   ", " A b a a   ", " A c a a   ", " B b c b   ", " B b a b   ", &
    " F m m m   ", " F m m m   ", " F m m m   ", " F m m m   ", " F m m m   ", " F m m m   ", &
    " F d d d   ", " F d d d   ", " F d d d   ", " F d d d   ", " F d d d   ", " F d d d   ", &
    " I m m m   ", " I m m m   ", " I m m m   ", " I m m m   ", " I m m m   ", " I m m m   ", &
    " I b a m   ", " I b a m   ", " I m c b   ", " I m c b   ", " I c m a   ", " I c m a   ", &
    " I b c a   ", " I c a b   ", " I b c a   ", " I c a b   ", " I b c a   ", " I c a b   ", &
    " I m m a   ", " I m m b   ", " I b m m   ", " I c m m   ", " I m c m   ", " I m a m   " /), (/6, 59/) )
!DEC$ ATTRIBUTES DLLEXPORT :: extendedHMOrthsymbols


!>  SYM_GL  encoded generator strings
  character(40), public, dimension(237) :: SYM_GL= (/  &
"000                                     ","100                                     ","01cOOO0                                 ", &
"01cODO0                                 ","02aDDOcOOO0                             ","01jOOO0                                 ", &
"01jOOD0                                 ","02aDDOjOOO0                             ","02aDDOjOOD0                             ", &
"11cOOO0                                 ","11cODO0                                 ","12aDDOcOOO0                             ", &
"11cOOD0                                 ","11cODD0                                 ","12aDDOcOOD0                             ", &
"02bOOOcOOO0                             ","02bOODcOOD0                             ","02bOOOcDDO0                             ", &
"02bDODcODD0                             ","03aDDObOODcOOD0                         ","03aDDObOOOcOOO0                         ", &
"04aODDaDODbOOOcOOO0                     ","03aDDDbOOOcOOO0                         ","03aDDDbDODcODD0                         ", &
"02bOOOjOOO0                             ","02bOODjOOD0                             ","02bOOOjOOD0                             ", &
"02bOOOjDOO0                             ","02bOODjDOO0                             ","02bOOOjODD0                             ", &
"02bDODjDOD0                             ","02bOOOjDDO0                             ","02bOODjDDO0                             ", &
"02bOOOjDDD0                             ","03aDDObOOOjOOO0                         ","03aDDObOODjOOD0                         ", &
"03aDDObOOOjOOD0                         ","03aODDbOOOjOOO0                         ","03aODDbOOOjODO0                         ", &
"03aODDbOOOjDOO0                         ","03aODDbOOOjDDO0                         ","04aODDaDODbOOOjOOO0                     ", &
"04aODDaDODbOOOjBBB0                     ","03aDDDbOOOjOOO0                         ","03aDDDbOOOjDDO0                         ", &
"03aDDDbOOOjDOO0                         ","12bOOOcOOO0                             ","03bOOOcOOOhDDD1BBB                      ", &
"12bOOOcOOD0                             ","03bOOOcOOOhDDO1BBO                      ","12bDOOcOOO0                             ", &
"12bDOOcDDD0                             ","12bDODcDOD0                             ","12bDOOcOOD0                             ", &
"12bOOOcDDO0                             ","12bDDOcODD0                             ","12bOODcODD0                             ", &
"12bOOOcDDD0                             ","03bOOOcDDOhDDO1BBO                      ","12bDDDcOOD0                             ", &
"12bDODcODD0                             ","12bDODcODO0                             ","13aDDObOODcOOD0                         ", &
"13aDDObODDcODD0                         ","13aDDObOOOcOOO0                         ","13aDDObOOOcOOD0                         ", &
"13aDDObODOcODO0                         ","04aDDObDDOcOOOhODD1OBB                  ","14aODDaDODbOOOcOOO0                     ", &
"05aODDaDODbOOOcOOOhBBB1ZZZ              ","13aDDDbOOOcOOO0                         ","13aDDDbOOOcDDO0                         ", &
"13aDDDbDODcODD0                         ","13aDDDbODOcODO0                         ","02bOOOgOOO0                             ", &
"02bOODgOOB0                             ","02bOOOgOOD0                             ","02bOODgOOF0                             ", &
"03aDDDbOOOgOOO0                         ","03aDDDbDDDgODB0                         ","02bOOOmOOO0                             ", &
"03aDDDbOOOmOOO0                         ","12bOOOgOOO0                             ","12bOOOgOOD0                             ", &
"03bOOOgDDOhDDO1YBO                      ","03bOOOgDDDhDDD1YYY                      ","13aDDDbOOOgOOO0                         ", &
"04aDDDbDDDgODBhODB1OYZ                  ","03bOOOgOOOcOOO0                         ","03bOOOgDDOcDDO0                         ", &
"03bOODgOOBcOOO0                         ","03bOODgDDBcDDB0                         ","03bOOOgOODcOOO0                         ", &
"03bOOOgDDDcDDD0                         ","03bOODgOOFcOOO0                         ","03bOODgDDFcDDF0                         ", &
"04aDDDbOOOgOOOcOOO0                     ","04aDDDbDDDgODBcDOF0                     ","03bOOOgOOOjOOO0                         ", &
"03bOOOgOOOjDDO0                         ","03bOOOgOODjOOD0                         ","03bOOOgDDDjDDD0                         ", &
"03bOOOgOOOjOOD0                         ","03bOOOgOOOjDDD0                         ","03bOOOgOODjOOO0                         ", &
"03bOOOgOODjDDO0                         ","04aDDDbOOOgOOOjOOO0                     ","04aDDDbOOOgOOOjOOD0                     ", &
"04aDDDbDDDgODBjOOO0                     ","04aDDDbDDDgODBjOOD0                     ","03bOOOmOOOcOOO0                         ", &
"03bOOOmOOOcOOD0                         ","03bOOOmOOOcDDO0                         ","03bOOOmOOOcDDD0                         ", &
"03bOOOmOOOjOOO0                         ","03bOOOmOOOjOOD0                         ","03bOOOmOOOjDDO0                         ", &
"03bOOOmOOOjDDD0                         ","04aDDDbOOOmOOOjOOO0                     ","04aDDDbOOOmOOOjOOD0                     ", &
"04aDDDbOOOmOOOcOOO0                     ","04aDDDbOOOmOOOcDOF0                     ","13bOOOgOOOcOOO0                         ", &
"13bOOOgOOOcOOD0                         ","04bOOOgOOOcOOOhDDO1YYO                  ","04bOOOgOOOcOOOhDDD1YYY                  ", &
"13bOOOgOOOcDDO0                         ","13bOOOgOOOcDDD0                         ","04bOOOgDDOcDDOhDDO1YBO                  ", &
"04bOOOgDDOcDDDhDDO1YBO                  ","13bOOOgOODcOOO0                         ","13bOOOgOODcOOD0                         ", &
"04bOOOgDDDcOODhDDD1YBY                  ","04bOOOgDDDcOOOhDDD1YBY                  ","13bOOOgOODcDDO0                         ", &
"13bOOOgDDDcDDD0                         ","04bOOOgDDDcDDDhDDD1YBY                  ","04bOOOgDDDcDDOhDDD1YBY                  ", &
"14aDDDbOOOgOOOcOOO0                     ","14aDDDbOOOgOOOcOOD0                     ","05aDDDbDDDgODBcDOFhODB1OBZ              ", &
"05aDDDbDDDgODBcDOBhODB1OBZ              ","01nOOO0                                 ","01nOOC0                                 ", &
"01nOOE0                                 ","02aECCnOOO0                             ","11nOOO0                                 ", &
"12aECCnOOO0                             ","02nOOOfOOO0                             ","02nOOOeOOO0                             ", &
"02nOOCfOOE0                             ","02nOOCeOOO0                             ","02nOOEfOOC0                             ", &
"02nOOEeOOO0                             ","03aECCnOOOeOOO0                         ","02nOOOkOOO0                             ", &
"02nOOOlOOO0                             ","02nOOOkOOD0                             ","02nOOOlOOD0                             ", &
"03aECCnOOOkOOO0                         ","03aECCnOOOkOOD0                         ","12nOOOfOOO0                             ", &
"12nOOOfOOD0                             ","12nOOOeOOO0                             ","12nOOOeOOD0                             ", &
"13aECCnOOOeOOO0                         ","13aECCnOOOeOOD0                         ","02nOOObOOO0                             ", &
"02nOOCbOOD0                             ","02nOOEbOOD0                             ","02nOOEbOOO0                             ", &
"02nOOCbOOO0                             ","02nOOObOOD0                             ","02nOOOiOOO0                             ", &
"12nOOObOOO0                             ","12nOOObOOD0                             ","03nOOObOOOeOOO0                         ", &
"03nOOCbOODeOOC0                         ","03nOOEbOODeOOE0                         ","03nOOEbOOOeOOE0                         ", &
"03nOOCbOOOeOOC0                         ","03nOOObOODeOOO0                         ","03nOOObOOOkOOO0                         ", &
"03nOOObOOOkOOD0                         ","03nOOObOODkOOD0                         ","03nOOObOODkOOO0                         ", &
"03nOOOiOOOkOOO0                         ","03nOOOiOODkOOD0                         ","03nOOOiOOOeOOO0                         ", &
"03nOOOiOODeOOO0                         ","13nOOObOOOeOOO0                         ","13nOOObOOOeOOD0                         ", &
"13nOOObOODeOOD0                         ","13nOOObOODeOOO0                         ","03bOOOcOOOdOOO0                         ", &
"05aODDaDODbOOOcOOOdOOO0                 ","04aDDDbOOOcOOOdOOO0                     ","03bDODcODDdOOO0                         ", &
"04aDDDbDODcODDdOOO0                     ","13bOOOcOOOdOOO0                         ","04bOOOcOOOdOOOhDDD1YYY                  ", &
"15aODDaDODbOOOcOOOdOOO0                 ","06aODDaDODbOOOcOOOdOOOhBBB1ZZZ          ","14aDDDbOOOcOOOdOOO0                     ", &
"13bDODcODDdOOO0                         ","14aDDDbDODcODDdOOO0                     ","04bOOOcOOOdOOOeOOO0                     ", &
"04bOOOcOOOdOOOeDDD0                     ","06aODDaDODbOOOcOOOdOOOeOOO0             ","06aODDaDODbODDcDDOdOOOeFBF0             ", &
"05aDDDbOOOcOOOdOOOeOOO0                 ","04bDODcODDdOOOeBFF0                     ","04bDODcODDdOOOeFBB0                     ", &
"05aDDDbDODcODDdOOOeFBB0                 ","04bOOOcOOOdOOOlOOO0                     ","06aODDaDODbOOOcOOOdOOOlOOO0             ", &
"05aDDDbOOOcOOOdOOOlOOO0                 ","04bOOOcOOOdOOOlDDD0                     ","06aODDaDODbOOOcOOOdOOOlDDD0             ", &
"05aDDDbDODcODDdOOOlBBB0                 ","14bOOOcOOOdOOOeOOO0                     ","05bOOOcOOOdOOOeOOOhDDD1YYY              ", &
"14bOOOcOOOdOOOeDDD0                     ","05bOOOcOOOdOOOeDDDhDDD1YYY              ","16aODDaDODbOOOcOOOdOOOeOOO0             ", &
"16aODDaDODbOOOcOOOdOOOeDDD0             ","07aODDaDODbODDcDDOdOOOeFBFhBBB1ZZZ      ","07aODDaDODbODDcDDOdOOOeFBFhFFF1XXX      ", &
"15aDDDbOOOcOOOdOOOeOOO0                 ","15aDDDbDODcODDdOOOeFBB0                 ","01dOOO0                                 ", &
"11dOOO0                                 ","02dOOOfOOO0                             ","02dOOOlOOO0                             ", &
"02dOOOlDDD0                             ","12dOOOfOOO0                             ","12dOOOfDDD0                             "/)
!DEC$ ATTRIBUTES DLLEXPORT :: SYM_GL

!> SGXsym contains the first space group of each crystal system
  integer(kind=irg), public, dimension(7) :: SGXsym = (/ 1, 3, 16, 75, 143, 168, 195 /)
!DEC$ ATTRIBUTES DLLEXPORT :: SGXsym

!>  SGPG contains the first space group # for a given point group
  integer(kind=irg), public, dimension(32):: SGPG =(/ &
                                          1,2,3,6,10,16,25,47,75,81,83,89,99,111,123,143, &
                                          147,149,156,162,168,174,175,177,183,187,191,195, &
                                          200,207,215,221/)
!DEC$ ATTRIBUTES DLLEXPORT :: SGPG

!>  SGsym contains the numbers of all the symmorphic space groups
  integer(kind=irg), public, dimension(73) :: SGsym =(/ &
                                            1,2,3,5,6,8,10,12,16,21,22,23,25,35,38,42,44,47, &
                                            65,69,71,75,79,81,82,83,87,89,97,99,107,111,115, &
                                            119,121,123,139,143,146,147,148,149,150,155,156, &
                                            157,160,162,164,166,168,174,175,177,183,187,189, &
                                            191,195,196,197,200,202,204,207,209,211,215,216, &
                                            217,221,225,229/)
!DEC$ ATTRIBUTES DLLEXPORT :: SGsym

!> SGsymnum contains the number of the symmorphic space group with the same point group symmetry
!>     this is necessary because sometimes the numbering of the space groups is not continuous in
!>     terms of the underlying point group ... e.g. 158 has the same point group as 156, not 157.
  integer(kind=irg), public, dimension(230) :: SGsymnum(230) = (/ &
                                           1,   2,   3,   3,   3,   6,   6,   6,   6,  10, &
                                          10,  10,  10,  10,  10,  16,  16,  16,  16,  16, &
                                          16,  16,  16,  16,  25,  25,  25,  25,  25,  25, &
                                          25,  25,  25,  25,  25,  25,  25,  25,  25,  25, &
                                          25,  25,  25,  25,  25,  25,  47,  47,  47,  47, &
                                          47,  47,  47,  47,  47,  47,  47,  47,  47,  47, &
                                          47,  47,  47,  47,  47,  47,  47,  47,  47,  47, &
                                          47,  47,  47,  47,  75,  75,  75,  75,  75,  75, &
                                          81,  81,  83,  83,  83,  83,  83,  83,  89,  89, &
                                          89,  89,  89,  89,  89,  89,  89,  89,  99,  99, &
                                          99,  99,  99,  99,  99,  99,  99,  99,  99,  99, &
                                         111, 111, 111, 111, 115, 115, 115, 115, 115, 115, &
                                         111, 111, 123, 123, 123, 123, 123, 123, 123, 123, &
                                         123, 123, 123, 123, 123, 123, 123, 123, 123, 123, &
                                         123, 123, 143, 143, 143, 143, 147, 147, 149, 150, &
                                         149, 150, 149, 150, 150, 156, 157, 156, 157, 156, &
                                         156, 162, 162, 164, 164, 164, 164, 168, 168, 168, &
                                         168, 168, 168, 174, 175, 175, 177, 177, 177, 177, &
                                         177, 177, 183, 183, 183, 183, 187, 187, 189, 189, &
                                         191, 191, 191, 191, 195, 195, 195, 195, 195, 200, &
                                         200, 200, 200, 200, 200, 200, 207, 207, 207, 207, &
                                         207, 207, 207, 207, 215, 215, 215, 215, 215, 215, &
                                         221, 221, 221, 221, 221, 221, 221, 221, 221, 221 /)
!DEC$ ATTRIBUTES DLLEXPORT :: SGsymnum

! space group order
  integer(kind=irg), public, dimension(230) :: SGorder(230) = (/ &
                                           1,   2,   2,   2,   4,   2,   2,   4,   4,   4, &
                                           4,   8,   4,   4,   8,   4,   4,   4,   4,   8, &
                                           8,  16,   8,   8,   4,   4,   4,   4,   4,   4, &
                                           4,   4,   4,   4,   8,   8,   8,   8,   8,   8, &
                                           8,  16,  16,   8,   8,   8,   8,   8,   8,   8, &
                                           8,   8,   8,   8,   8,   8,   8,   8,   8,   8, &
                                           8,   8,  16,  16,  16,  16,  16,  16,  32,  32, &
                                          16,  16,  16,  16,   4,   4,   4,   4,   8,   8, &
                                           4,   8,   8,   8,   8,   8,  16,  16,   8,   8, &
                                           8,   8,   8,   8,   8,   8,  16,  16,   8,   8, &
                                           8,   8,   8,   8,   8,   8,  16,  16,  16,  16, &
                                           8,   8,   8,   8,   8,   8,   8,   8,  16,  16, &
                                          16,  16,  16,  16,  16,  16,  16,  16,  16,  16, &
                                          16,  16,  16,  16,  16,  16,  16,  16,  32,  32, &
                                          32,  32,   3,   3,   3,   9,   6,  18,   6,   6, &
                                           6,   6,   6,   6,  18,   6,   6,   6,   6,  18, &
                                          18,  12,  12,  12,  12,  36,  36,   6,   6,   6, &
                                           6,   6,   6,   6,  12,  12,  12,  12,  12,  12, &
                                          12,  12,  12,  12,  12,  12,  12,  12,  12,  12, &
                                          24,  24,  24,  24,  12,  48,  24,  12,  24,  24, &
                                          24,  96,  96,  48,  24,  48,  24,  24,  96,  96, &
                                          48,  24,  24,  48,  24,  96,  48,  24,  96,  48, &
                                          48,  48,  48,  48, 192, 192, 192, 192,  96,  96 /)
!DEC$ ATTRIBUTES DLLEXPORT :: SGorder

! these parameters implement the diffraction group
! formalism described in the BESR paper.

!> 10 2D point group symbols in International Tables order
  character(10), public, dimension(0:11)  :: PGTWD = (/ &
                        ' none     ','    1     ','    2     ','    m     ','  2mm     ','    4     ', &
                        '  4mm     ','    3     ','   3m1    ','    6     ','  6mm     ','   31m    ' /)
!DEC$ ATTRIBUTES DLLEXPORT :: PGTWD

!> 10 2D point group orders in International Tables order
  integer(kind=irg), public, dimension(0:11) :: PGTWDorder = (/0,1,2,2,4,4,8,3,6,6,12,6/)
!DEC$ ATTRIBUTES DLLEXPORT :: PGTWDorder

!> inverse table for 2D point groups; this essentially implements the inverse of Table 4 in BESR paper for the Bright Field symmetry.
  integer(kind=irg), public, dimension(12,11) :: PGTWDinverse = reshape((/ &
                                   1,0,0,0,0,0,0,0,0,0,0,0,  1,2,0,0,0,0,0,0,0,0,0,0, &
                                   1,3,0,4,0,0,0,0,0,0,0,0,  1,3,0,5,0,0,0,0,0,0,0,0, &
                                   1,3,0,4,0,0,0,6,0,0,0,0,  1,0,7,0,0,0,0,0,0,0,0,0, &
                                   1,2,0,0,0,8,0,0,0,0,0,0,  1,3,0,0,0,9,0,0,0,0,0,0, &
                                   1,3,0,4,0,0,0,0,0,0,0,10, 1,3,7,4,0,0,0,0,0,0,0,0, &
                                   1,3,0,4,0,8,0,6,0,0,0,0 /), (/ 12,11 /))
!DEC$ ATTRIBUTES DLLEXPORT :: PGTWDinverse


!> 32 3D point group symbols in International Tables order; additional quasi-crystal rotational
!> groups are added at the end of the list
  character(5), public, dimension(36):: PGTHD =(/ '    1','   -1','    2','    m','  2/m','  222', &
                                                  '  mm2','  mmm','    4','   -4','  4/m','  422', &
                                                  '  4mm',' -42m','4/mmm','    3','   -3','   32', &
                                                  '   3m','  -3m','    6','   -6','  6/m','  622', &
                                                  '  6mm',' -6m2','6/mmm','   23','   m3','  432', &
                                                  ' -43m',' m-3m','  532','  822',' 1022',' 1222' /)
!DEC$ ATTRIBUTES DLLEXPORT :: PGTHD

!> 32 3D point group orders in International Tables order
  integer(kind=irg), public, dimension(32)  :: PGTHDorder = (/ 1, 2, 2, 2, 4, 4, 4, 8, 4, 8, &
                                                               8, 8, 8, 8,16, 3, 6, 6, 6,12, &
                                                               6,12,12,12,12,12,24,12,24,24, &
                                                              24,48 /)
!DEC$ ATTRIBUTES DLLEXPORT :: PGTHDorder

!> 3D point groups : purely rotational point groups corresponding to each point group
  integer(kind=irg), public, dimension(36)   :: PGrot = (/ &
                                                1,1,3,1,3,6,3,6,9,3,9,12,9,6,12,16,16, &
                                                18,16,18,21,16,21,24,21,18,24,28,28,30, &
                                                28,30,33,34,35,36/)
!DEC$ ATTRIBUTES DLLEXPORT :: PGrot

!> 3D point groups : Laue group number
  integer(kind=irg), public, dimension(36)   :: PGLaue =(/2,2,5,5,5,8,8,8,11,11,11,15,15,15,15,17,17, &
                                                         20,20,20,23,23,23,27,27,27,27,29,29,32,32,32, &
                                                         33,34,35,36/)
!DEC$ ATTRIBUTES DLLEXPORT :: PGLaue

!> 3D point groups : inverted Laue group number
  integer(kind=irg), public, dimension(36)   :: PGLaueinv = (/1,1,2,2,2,3,3,3,4,4,4,5,5,5,5,6,6, &
                                                              7,7,7,8,8,8,9,9,9,9,10,10,11,11,11, &
                                                              12,13,14,15/)
!DEC$ ATTRIBUTES DLLEXPORT :: PGLaueinv

!> 3D point groups mapped onto kvector sampling type (used for master pattern computations) [-1 for special cases]
  integer(kind=irg), public, dimension(36)   :: PGSamplingType = (/ &
                                                        1, 2, 3, 4, 5, 5, 5, 6, 5, 5, &
                                                        6, 6, 7,-1, 9,-1,-1,-1,-1,-1, &
                                                       15,12,17,16,18,-1,19, 3, 6, 6, &
                                                        8, 9, -1, -1, -1, -1 /)
!DEC$ ATTRIBUTES DLLEXPORT :: PGSamplingType

! There are 24 space groups with two origin choices.
! The symmetry of both sites is stored in the array
! sitesym; the space group numbers are stored
! in tworig.
! numbers of the space groups with two settings
integer(kind=irg),public, parameter  :: tworig(24)=(/48,50,59,68,70,85,86,88,125,126,129,130,133,134,137,138,&
                                                     141,142,201,203,222,224,227,228/)
!DEC$ ATTRIBUTES DLLEXPORT :: tworig

! site symmetry list
character(7),public, parameter :: sitesym(48) = (/ '222    ',' -1    ','222/n  ',' -1    ','mm2/n  ',' -1    ', &
                                                   '222    ',' -1    ','222    ',' -1    ','-4     ',' -1    ', &
                                                   '-4     ',' -1    ','-4     ',' -1    ','422    ','2/m    ', &
                                                   '422/n  ',' -1    ','-4m2   ','2/m    ','-4/ncn ',' -1    ', &
                                                   '-4121/c',' -1    ','-42m   ','2/m    ','-4m2/n ',' -1    ', &
                                                   '-4cg   ','2/m    ','-4m2   ','2/m    ','-4c21  ',' -1    ', &
                                                   '23     ',' -3    ','23     ',' -3    ','432    ',' -3    ', &
                                                   '-43m   ','-3m    ','-43m   ','-3m    ','23     ',' -3    '/)
!DEC$ ATTRIBUTES DLLEXPORT :: sitesym

!> 31 diffraction group symbols in BESR order
  character(5), public, dimension(31)  :: DG =(/ '    1','   1R','    2','   2R','  21R','   mR', &
                                                 '    m','  m1R','2mRmR','  2mm','2RmmR','2mm1R', &
                                                 '    4','   4R','  41R','4mRmR','  4mm','4RmmR', &
                                                 '4mm1R','    3','   6R','  3mR','   3m','6RmmR', &
                                                 '    6','  31R','  61R','6mRmR','  6mm',' 3m1R', &
                                                 '6mm1R' /)
!DEC$ ATTRIBUTES DLLEXPORT :: DG

!> 31 diffraction group orders in BESR order
  integer(kind=irg), public, dimension(31) :: DGorder =(/ 1, 2, 2, 2, 4, 2, 2, 4, 4, 4, 4, 8, &
                                                          4, 4, 8, 8, 8, 8,16, 3, 6, 6, 6,12, &
                                                          6, 6,12,12,12,12,24 /)
!DEC$ ATTRIBUTES DLLEXPORT :: DGorder

!> Bright Field planar point group for 31 diffraction groups (Table 2, column 2, BESR, with change in row ordering)
  integer(kind=irg), public, dimension(31) :: BFPG =(/ &
                                              1,2,2,1,2,3,3,4,4,4,3,4,5,5,5,6,6,6,6,7,7,8,8,8,9,9,9,10,10,10,10/)
!DEC$ ATTRIBUTES DLLEXPORT :: BFPG

!> Whole Pattern planar point group for 31 diffraction groups (Table 2, column 3, BESR, with change in row ordering)
  integer(kind=irg), public, dimension(31) :: WPPG =(/ &
                                              1,1,2,1,2,1,3,3,2,4,3,4,5,2,5,5,6,4,6,7,7,7,8,8,9,7,9,9,10,8,10/)
!DEC$ ATTRIBUTES DLLEXPORT :: WPPG

!> Dark Field planar point group for 31 diffraction groups (Table 2, column 4, BESR, with change in row ordering)
  integer(kind=irg), public, dimension(31) :: DFGN = (/ &
                                              1,2,1,1,2,1,1,2,1,1,1,2,1,1,2,1,1,1,2,1,1,1,1,1,1,2,2,1,1,2,2/)
!DEC$ ATTRIBUTES DLLEXPORT :: DFGN

!> Dark Field planar point group for 31 diffraction groups (Table 2, column 5, BESR, with change in row ordering)
  integer(kind=irg), public, dimension(31) :: DFSP = (/ &
                                              0,0,0,0,0,3,3,4,3,3,3,4,0,0,0,3,3,3,4,0,0,3,3,3,0,0,0,3,3,4,4/)
!DEC$ ATTRIBUTES DLLEXPORT :: DFSP

!> 10 projection diffraction groups in BESR order (Table 2, column 8, BESR, with change in row ordering)
  integer(kind=irg), public, dimension(31) :: PDG = (/ &
                                              2,2,5,5,5,8,8,8,12,12,12,12,15,15,15,19,19,19,19,26,27,30,30, &
                                              31,27,26,27,31,31,30,31/)
!DEC$ ATTRIBUTES DLLEXPORT :: PDG

!> short hand for .FALSE. logical parameter
  logical,parameter :: FF=.FALSE.

!> short hand for .TRUE. logical parameter
  logical,parameter :: TT=.TRUE.

!> Table 3 from BESR paper
  logical, public, dimension(32,31)  :: DGPG = reshape((/ &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,TT,FF,FF,TT, &
     FF,FF,FF,FF,TT,FF,FF,TT,FF,FF,TT,FF,FF,FF,TT,FF,FF,FF,FF,TT,FF,FF,TT,FF,FF,FF,TT,FF,TT,FF,FF,TT, &
     FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,TT,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,TT,FF,TT,FF,FF, &
     FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,TT,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,TT,FF,FF,FF,FF,TT,FF, &
     FF,FF,FF,TT,FF,FF,TT,FF,FF,FF,FF,FF,TT,TT,FF,FF,FF,FF,TT,FF,FF,TT,FF,FF,TT,TT,FF,FF,FF,FF,TT,FF, &
     FF,FF,TT,FF,FF,TT,TT,FF,TT,TT,FF,TT,TT,TT,FF,FF,FF,TT,FF,FF,TT,FF,FF,TT,TT,TT,FF,TT,FF,TT,TT,FF, &
     FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,TT,FF,FF,TT,FF,FF,TT,FF,FF,TT,FF,FF,FF,TT,FF,TT,FF,FF,TT,FF,FF,TT,FF,FF,FF,TT,FF,TT,FF,FF,TT, &
     FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,TT,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF, &
     TT,FF,TT,TT,FF,TT,TT,FF,TT,TT,FF,TT,TT,TT,FF,TT,FF,TT,TT,FF,TT,TT,FF,TT,TT,TT,FF,TT,FF,TT,TT,FF/), (/32,31/))
!DEC$ ATTRIBUTES DLLEXPORT :: DGPG

! the following arrays are used for the symmetry compression step in
! the spherical indexing (EMSphInx) package
  integer(kind=irg), public, dimension(230) :: SHT_ZRot = (/ &
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, &
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
                    2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, &
                    2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
                    4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
                    3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
                    6, 6, 3, 3, 3, 3, 6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, &
                    4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 /)
!DEC$ ATTRIBUTES DLLEXPORT :: SHT_ZRot

  integer(kind=irg), public, dimension(230) :: SHT_mirInv = (/ &
                    0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, &
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
                    3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, &
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                    0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
                    3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                    1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, &
                    0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 0, &
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
!DEC$ ATTRIBUTES DLLEXPORT :: SHT_mirInv

character(2), public, dimension(32) :: TSLsymtype = (/' 1',' 1',' 2',' 2',' 2','22','22','22', &
                                                       ' 4',' 4',' 4','42','42','42','42',' 3', &
                                                       ' 3','32','32','32',' 6',' 6',' 6','62', &
                                                       '62','62','62','23','23','43','43','43'/)
!DEC$ ATTRIBUTES DLLEXPORT :: TSLsymtype


!>  SYM_GL  encoded point group generator strings
  character(4), dimension(32) :: SYM_PGGL= (/  &
    '----', 'h---', 'c---', 'j---', 'ch--', 'bc--', 'bj--', 'bch-', &
    'g---', 'm---', 'gh--', 'cg--', 'gj--', 'cm--', 'cgh-', 'n---', &
    'hn--', 'en--', 'kn--', 'fhn-', 'bn--', 'in--', 'bhn-', 'ben-', &
    'bkn-', 'ikn-', 'benh', 'cd--', 'cdh-', 'dg--', 'dml-', 'dghl' /)

!DEC$ ATTRIBUTES DLLEXPORT :: SYM_PGGL

! here is the main symmetry class definition; this class must be able to function on its
! own, without knowledge of a particular crystal structure, but also in conjunction with the
! crystallography (Cell_T) class

  type, public :: SpaceGroup_T
    !! symmetry class definition
    private
      integer(kind=irg)     :: SGnumber = 0
       !! space group number
      integer(kind=irg)     :: PGnumber = 0
       !! space group number
      character(11)         :: SGname = '           '
       !! current space group name
      integer(kind=irg)     :: setting = 0
       !! space group setting (for those groups that have multiple settings)
      integer(kind=irg)     :: SGorder = 0
       !! space group order
      integer(kind=irg)     :: xtal_system = 0
       !! crystal system number
      logical               :: hexset = .FALSE.
       !! are we using a hexagonal setting?
      logical               :: symmorphic
       !! is this a symmorphic space group?
      integer(kind=irg)     :: GENnum = 0
       !! number of generator matrices
      integer(kind=irg)     :: MATnum = 0
       !! number of non-zero symmetry matrices
      integer(kind=irg)     :: NUMpt = 0
       !! number of point group operators
      logical               :: reduce = .TRUE.
       !! switch to enable/disable reduction to fundamental cell
      logical               :: trigonal = .FALSE.
       !! switch for hexagonal vs. rhombohedral settings
      logical               :: second = .FALSE.   ! we always use hexagonal setting for trigonal case
       !! switch for second setting of spacegroup (if any)
      logical               :: centrosym = .FALSE.
       !! switch for presence of centrosymmetry
      real(kind=dbl),allocatable :: data(:,:,:)   ! data(192,4,4)
       !! all symmetry matrices for a given spacegroup
      real(kind=dbl),allocatable,public :: direc(:,:,:)  ! direc(48,3,3)
       !! direct space point group matrices
      real(kind=dbl),allocatable :: recip(:,:,:)  ! recip(48,3,3)
       !! reciprocal space point group matrices
      real(kind=dbl)        :: c(4,4)
       !! dummy 4x4 matrix used for various computations
      logical, public       :: recip_pending = .TRUE.
       !! the reciprocal point group matrices require the metric tensors; their computation remains
       !! to be completed if recip_pending = .TRUE.; this can be completed by means of the
       !! fixRecipPG method which can be called the metric matrices are known. This is only
       !! relevant to the hexagonal/trigonal case.

    contains
    private
! basic space group generating routines and related stuff
      procedure, pass(self) :: getSpaceGroup_
      procedure, pass(self) :: GetSetting_
      procedure, pass(self) :: fillgen_
      procedure, pass(self) :: MakeGenerators_
      procedure, pass(self) :: matrixmult_
      procedure, pass(self) :: isitnew_
      procedure, pass(self) :: GenerateSymmetry_
      procedure, pass(self) :: ListPointGroups_
      procedure, pass(self) :: fixRecipPG_
      procedure, pass(self) :: resetSpaceGroup_
! routines to extract space group parameters
      procedure, pass(self) :: getSpaceGroupName_
      procedure, pass(self) :: getSpaceGroupOrder_
      procedure, pass(self) :: getSpaceGroupNumber_
      procedure, pass(self) :: getSpaceGroupGENnum_
      procedure, pass(self) :: getSpaceGroupMATnum_
      procedure, pass(self) :: getSpaceGroupNUMpt_
      procedure, pass(self) :: getSpaceGroupSetting_
      procedure, pass(self) :: getSpaceGroupCentro_
      procedure, pass(self) :: getSpaceGroupXtalSystem_
      procedure, pass(self) :: getSpaceGroupSymmorphic_
      procedure, pass(self) :: getSpaceGroupDataMatrices_
      procedure, pass(self) :: getSpaceGroupPGdirecMatrices_
      procedure, pass(self) :: getSpaceGroupPGrecipMatrices_
      procedure, pass(self) :: getSpaceGrouptrigonal_
      procedure, pass(self) :: getSpaceGroupsecond_
      procedure, pass(self) :: getSpaceGrouphexset_
! routines to set space group parameters
      procedure, pass(self) :: setSpaceGroupreduce_
      procedure, pass(self) :: setSpaceGrouphexset_
      procedure, pass(self) :: setSpaceGrouptrigonal_
      procedure, pass(self) :: setSpaceGroupsecond_
      procedure, pass(self) :: setSpaceGroupSetting_
      procedure, pass(self) :: setSpaceGroupNumber_
      procedure, pass(self) :: setSpaceGroupXtalSystem_
! general purpose routines that use symmetry
      procedure, pass(self) :: CalcFamily_
      procedure, pass(self) :: CalcOrbit_
      procedure, pass(self) :: CalcStar_
! routines to extract particular parameters
      procedure, pass(self) :: getPGnumber_
      procedure, pass(self) :: setPGnumber_
      procedure, pass(self) :: GetOrderinZone_
      procedure, pass(self) :: getLaueGroupNumber_
      procedure, pass(self) :: getHexvsRho_
! routines related to the Wyckoff positions
      procedure, pass(self) :: getmultiplicity_
      procedure, pass(self) :: getposition_
      procedure, pass(self) :: getWPstring_
      procedure, pass(self) :: printWyckoffPositions_
      procedure, pass(self) :: extractWyckoffposition_
      procedure, pass(self) :: interpretWyckoffletter_
! diffraction-related routines
      procedure, pass(self) :: IsGAllowed_
      procedure, pass(self) :: GetDiffractionGroup_
      procedure, pass(self) :: BFsymmetry_
      procedure, pass(self) :: GetPatternSymmetry_
! finally, define the destructor routine
      final :: SpaceGroup_destructor

! generic (public) function definitions and overloads
      generic, public :: getSpaceGroup => getSpaceGroup_
      generic, public :: resetSpaceGroup => resetSpaceGroup_
      generic, public :: getSpaceGroupName => getSpaceGroupName_
      generic, public :: getSpaceGroupGENnum => getSpaceGroupGENnum_
      generic, public :: getSpaceGroupMATnum => getSpaceGroupMATnum_
      generic, public :: getSpaceGroupNUMpt => getSpaceGroupNUMpt_
      generic, public :: getSpaceGroupOrder => getSpaceGroupOrder_
      generic, public :: getSpaceGroupNumber => getSpaceGroupNumber_
      generic, public :: getSpaceGroupSetting => getSpaceGroupSetting_
      generic, public :: getSpaceGroupCentro => getSpaceGroupCentro_
      generic, public :: getSpaceGroupXtalSystem => getSpaceGroupXtalSystem_
      generic, public :: getSpaceGroupSymmorphic => getSpaceGroupSymmorphic_
      generic, public :: getSpaceGroupDataMatrices => getSpaceGroupDataMatrices_
      generic, public :: getSpaceGroupPGdirecMatrices => getSpaceGroupPGdirecMatrices_
      generic, public :: getSpaceGroupPGrecipMatrices => getSpaceGroupPGrecipMatrices_
      generic, public :: getSpaceGrouphexset => getSpaceGrouphexset_
      generic, public :: getSpaceGrouptrigonal => getSpaceGrouptrigonal_
      generic, public :: getSpaceGroupsecond => getSpaceGroupsecond_
      generic, public :: setSpaceGroupreduce => setSpaceGroupreduce_
      generic, public :: setSpaceGrouphexset => setSpaceGrouphexset_
      generic, public :: setSpaceGrouptrigonal => setSpaceGrouptrigonal_
      generic, public :: setSpaceGroupsecond => setSpaceGroupsecond_
      generic, public :: setSpaceGroupSetting => setSpaceGroupSetting_
      generic, public :: setSpaceGroupNumber => setSpaceGroupNumber_
      generic, public :: setSpaceGroupXtalSystem => setSpaceGroupXtalSystem_
      generic, public :: GetSetting => GetSetting_
      generic, public :: GenerateSymmetry => GenerateSymmetry_
      generic, public :: ListPointGroups => ListPointGroups_
      generic, public :: fixRecipPG => fixRecipPG_

      generic, public :: CalcFamily => CalcFamily_
      generic, public :: CalcOrbit => CalcOrbit_
      generic, public :: CalcStar => CalcStar_

      generic, public :: getPGnumber => getPGnumber_
      generic, public :: setPGnumber => setPGnumber_
      generic, public :: GetOrderinZone => GetOrderinZone_
      generic, public :: getLaueGroupNumber => getLaueGroupNumber_
      generic, public :: getHexvsRho => getHexvsRho_
      generic, public :: getmultiplicity => getmultiplicity_
      generic, public :: getposition => getposition_
      generic, public :: getWPstring => getWPstring_
      generic, public :: printWyckoffPositions => printWyckoffPositions_
      generic, public :: extractWyckoffposition => extractWyckoffposition_
      generic, public :: interpretWyckoffletter => interpretWyckoffletter_
      generic, public :: IsGAllowed => IsGAllowed_
      generic, public :: GetDiffractionGroup => GetDiffractionGroup_
      generic, public :: BFsymmetry => BFsymmetry_
      generic, public :: GetPatternSymmetry => GetPatternSymmetry_

  end type SpaceGroup_T

! here is the main point symmetry class definition; this class fills in the direc arrays 
! in the Spacegroup_T class starting from the point group generators

  type, private :: PointGroup_T
    !! point group symmetry class definition (only generates the PG matrices, everything
    !! else is handled in the SpaceGroup_T class)
    private
      integer(kind=irg)           :: pgnum
      integer(kind=irg)           :: GENnum
      integer(kind=irg)           :: PGMATnum
      real(kind=dbl),allocatable  :: direc(:,:,:)  ! direc(48,3,3)
      real(kind=dbl)              :: c(3,3)

    contains
    private
! basic space group generating routines and related stuff
      procedure, pass(self) :: PGfillgen_
      procedure, pass(self) :: MakePGGenerators_
      procedure, pass(self) :: PGmatrixmult_
      procedure, pass(self) :: PGisitnew_
      procedure, pass(self) :: GeneratePGSymmetry_
! finally, define the destructor routine
      final :: PointGroup_destructor

! no generic (public) function definitions and overloads

  end type PointGroup_T

! the constructor routine for the SpaceGroup_T class
  interface SpaceGroup_T
    module procedure SpaceGroup_constructor
  end interface SpaceGroup_T

! the constructor routine for the PointGroup_T class
  interface PointGroup_T
    module procedure PointGroup_constructor
  end interface PointGroup_T

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in this class
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(SpaceGroup_T) function SpaceGroup_constructor( SGnumber, xtalSystem, setting ) result(SG)
!DEC$ ATTRIBUTES DLLEXPORT :: SpaceGroup_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! constructor for the SpaceGroup Class

use mod_io

IMPLICIT NONE

integer(kind=irg), intent(in), OPTIONAL :: SGnumber
integer(kind=irg), intent(in), OPTIONAL :: xtalSystem
integer(kind=irg), intent(in), OPTIONAL :: setting

type(IO_T)                              :: Message
integer(kind=irg)                       :: i, pgnum
integer(kind=irg),parameter             :: icv(7) = (/ 7, 6, 3, 2, 5, 4, 1 /)

! if there are no input parameters, then we start by asking for the crystal system
if ( .not.(present(SGnumber)) .and. .not.(present(xtalSystem)) .and. .not.(present(setting)) ) then
  call getXtalSystem_(SG)
  call getSpaceGroup_(SG)
  call getSetting_(SG)

else  ! at least one of the optional parameters are present

  if (present(SGnumber) .and. (.not.(present(xtalSystem))) ) then
! find the crystal system number from the space group number SGXsym = (/ 1, 3, 16, 75, 143, 168, 195 /)
    SG%SGnumber = SGnumber
    if (SGnumber.ge.SGXsym(7)) then 
        i = 7
    else
        i = 0
        do while(SGnumber.gt.SGXsym(i+1))
          i = i+1
        end do
    end if 
! convert this number to the EMsoft crystal system numbering scheme
    SG%xtal_system = icv(i)
  end if

  if ( (.not.(present(SGnumber))) .and. (present(xtalSystem)) ) then
    SG%xtal_system = xtalSystem
    call getSpaceGroup_(SG)
  end if

  if (present(SGnumber) .and. (present(xtalSystem)) ) then
    SG%SGnumber = SGnumber
    SG%xtal_system = xtalSystem
  end if

  if (present(setting)) then
    SG%setting = setting
  else
    call getSetting_(SG)
  end if

end if

! fill in the space group name  and order
SG%SGname = SYM_SGname(SG%SGnumber)
SG%SGorder = SGorder(SG%SGnumber)

! and convert the space group number into a point group number
pgnum = 0
do i=1,32
  if (SGPG(i).le.SG%SGnumber) pgnum = i
end do
call SG%setPGnumber(pgnum)

! allocate the arrays for symmetry operators
allocate( SG%data(SG%SGorder, 4, 4), SG%direc(PGTHDorder(pgnum),3,3), SG%recip(PGTHDorder(pgnum),3,3) )

SG%data = 0.D0
SG%direc = 0.D0
SG%recip = 0.D0

! generate all the symmetry operators as well as the corresponding point group symmetry
call GenerateSymmetry_(SG,.TRUE.)

end function SpaceGroup_constructor

!--------------------------------------------------------------------------
recursive subroutine SpaceGroup_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: SpaceGroup_destructor
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! clean up allocated arrays

IMPLICIT NONE

type(SpaceGroup_T), intent(inout)     :: self

call reportDestructor('SpaceGroup_T')

if (allocated(self%data)) deallocate(self%data)
if (allocated(self%direc)) deallocate(self%direc)
if (allocated(self%recip)) deallocate(self%recip)

end subroutine SpaceGroup_destructor

!--------------------------------------------------------------------------
type(PointGroup_T) function PointGroup_constructor( SG ) result(PG)
!DEC$ ATTRIBUTES DLLEXPORT :: PointGroup_constructor
  !! author: MDG
  !! version: 1.0
  !! date: 02/20/22
  !!
  !! constructor for the PointGroup Class

IMPLICIT NONE

class(SpaceGroup_T), INTENT(INOUT)  :: SG

integer(kind=irg)                   :: i, pgnum, SGnumber 

! get the space group number
SGnumber = SG%getSpaceGroupNumber()

! and convert the space group number into a point group number
pgnum = 0
do i=1,32
  if (SGPG(i).le.SGnumber) pgnum = i
end do
call SG%setPGnumber(pgnum)
PG%pgnum = pgnum

! allocate the array for symmetry operators
allocate( PG%direc(PGTHDorder(pgnum),3,3) )
PG%direc = 0.D0

! generate all the point group symmetry operators 
call PG%GeneratePGSymmetry_()

! and copy the direc array into the space group class 
if (allocated(SG%direc)) deallocate(SG%direc)
allocate(SG%direc(PGTHDorder(pgnum),3,3))
SG%direc = PG%direc
SG%NUMpt = PG%PGMATnum

! and clean up 
deallocate(PG%direc)

end function PointGroup_constructor

!--------------------------------------------------------------------------
recursive subroutine PointGroup_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: PointGroup_destructor
  !! author: MDG
  !! version: 1.0
  !! date: 02/20/22
  !!
  !! clean up allocated arrays

IMPLICIT NONE

type(PointGroup_T), intent(inout)     :: self

call reportDestructor('SpaceGroup_T')

if (allocated(self%direc)) deallocate(self%direc)

end subroutine PointGroup_destructor

!--------------------------------------------------------------------------
recursive subroutine resetSpaceGroup_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: resetSpaceGroup_
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! zero out an existing space group

IMPLICIT NONE

class(SpaceGroup_T), intent(inout)     :: self

if (allocated(self%data)) deallocate(self%data)
if (allocated(self%direc)) deallocate(self%direc)
if (allocated(self%recip)) deallocate(self%recip)

self%SGnumber = 0
self%SGname = ''
self%SGorder = 0
self%xtal_system = 0
self%NUMpt = 0
self%MATnum = 0
self%GENnum = 0
self%recip_pending = .TRUE.
self%second = .FALSE.
self%centrosym = .FALSE.

end subroutine resetSpaceGroup_

!--------------------------------------------------------------------------
recursive subroutine getXtalSystem_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getXtalSystem_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! ask the user for the crystal system

use mod_io

IMPLICIT NONE

class(SpaceGroup_T), intent(inout)     :: self

type(IO_T)                             :: Message
integer(kind=irg)                      :: io_int(1)

 call Message%printMessage(' Select the crystal system : ', frm = "(A)")
 call Message%printMessage('  1. Cubic ', frm = "(A)")
 call Message%printMessage('  2. Tetragonal ', frm = "(A)")
 call Message%printMessage('  3. Orthorhombic ', frm = "(A)")
 call Message%printMessage('  4. Hexagonal ', frm = "(A)")
 call Message%printMessage('  5. Trigonal ', frm = "(A)")
 call Message%printMessage('  6. Monoclinic ', frm = "(A)")
 call Message%printMessage('  7. Triclinic ', frm = "(A/)")
! call Message('  8. 2-D Quasi-Crystal', frm = "(A/)")
! call Message('  9. 3-D Quasi-Crystal', frm = "(A/)")

 call Message%printMessage(' Note about the trigonal system:', frm = "(A)")
 call Message%printMessage(' -------------------------------', frm = "(A)")
 call Message%printMessage(' Primitive trigonal crystals are defined with respect to a HEXAGONAL', frm = "(A)")
 call Message%printMessage(' reference frame.  Rhombohedral crystals can be referenced with', frm = "(A)")
 call Message%printMessage(' respect to a HEXAGONAL basis (first setting), or with respect to', frm = "(A)")
 call Message%printMessage(' a RHOMBOHEDRAL basis (second setting).  The default setting for ', frm = "(A)")
 call Message%printMessage(' trigonal symmetry is the hexagonal setting.  When you select', frm = "(A)")
 call Message%printMessage(' crystal system 5 above, you will be prompted for the setting. ', frm = "(A//)")
 call Message%ReadValue(' crystal system ---> ', io_int, 1)
 self%xtal_system = io_int(1)

end subroutine getXtalSystem_

!--------------------------------------------------------------------------
recursive subroutine getSpaceGroup_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroup_
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! asks the user for a space group number
  !!
  !! This routines lists all the relevant space groups for
  !! the present crystal system, and asks the user to pick one.
  !!
  !! the following space groups have a hexagonal and rhombohedral
  !! setting;  the generator matrices are different, and so are
  !! their entries in the SYM_GL array.\n
  !!  hexagonal setting 146: R 3           231\n
  !!  hexagonal setting 148: R -3          232\n
  !!  hexagonal setting 155: R 3 2         233\n
  !!  hexagonal setting 160: R 3 m         234\n
  !!  hexagonal setting 161: R 3 c         235\n
  !!  hexagonal setting 166: R -3 m        23\n
  !!  hexagonal setting 167: R -3 c        237\n
  !!
  !! TODO: automatic conversion from rhombohedral to hexagonal setting ...

use mod_io

IMPLICIT NONE

class(SpaceGroup_T), intent(inout)     :: self

type(IO_T)                             :: Message
integer(kind=irg)                      :: sgmin,sgmax,i,j,TRIG(7), io_int(1)
logical                                :: skip

 TRIG = (/ 146,148,155,160,161,166,167 /)
 skip = .FALSE.

 select case (self%xtal_system)
   case (1); sgmin = 195; sgmax = 230
   case (2); sgmin =  75; sgmax = 142
   case (3); sgmin =  16; sgmax =  74
   case (4); sgmin = 168; sgmax = 194
   case (5); if (self%second) then
               call Message%printMessage('The space groups below correspond to the ', frm = "(/A)")
               call Message%printMessage('second (rhombohedral) setting.', frm = "(A/)")
               call Message%printMessage('Please select one of these space groups.', frm = "(A/)")
               do i=1,7
                if ((mod(i,4).eq.0).or.(i.eq.7)) then
                  write (6,"(1x,i3,':',A11,5x)") TRIG(i),SYM_SGname(TRIG(i))
                else
                  write (6,"(1x,i3,':',A11,5x)",advance="no") TRIG(i),SYM_SGname(TRIG(i))
                end if
               end do
               call Message%printMessage(' -------------------------- ', frm = "(A)")
               call Message%ReadValue(' Enter space group number : ', io_int, 1)
               self%SGnumber = io_int(1)

  ! check for rhombohedral settings of rhombohedral space groups
               if (self%second) then
                 if (self%SGnumber.eq.146) self%SGnumber=231
                 if (self%SGnumber.eq.148) self%SGnumber=232
                 if (self%SGnumber.eq.155) self%SGnumber=233
                 if (self%SGnumber.eq.160) self%SGnumber=234
                 if (self%SGnumber.eq.161) self%SGnumber=235
                 if (self%SGnumber.eq.166) self%SGnumber=236
                 if (self%SGnumber.eq.167) self%SGnumber=237
               endif
               skip = .TRUE.
             else
              sgmin = 143
              sgmax = 167
             end if
   case (6); sgmin =   3; sgmax =  15
   case (7); sgmin =   1; sgmax =   2
 end select

! print out all the relevant space group names and numbers
 if (skip.eqv..FALSE.) then
  call Message%printMessage(' ', frm = "(/A/)")
  do i=sgmin,sgmax
   j=i-sgmin+1
   if ((mod(j,4).eq.0).or.(i.eq.sgmax)) then
    write (6,"(1x,i3,':',A11,5x)") i,SYM_SGname(i)
   else
    write (6,"(1x,i3,':',A11,5x)",advance="no") i,SYM_SGname(i)
   end if
  end do
  self%SGnumber = sgmin-1
  do while ((self%SGnumber.lt.sgmin).or.(self%SGnumber.gt.sgmax))
   call Message%printMessage(' -------------------------- ', frm = "(A)")
   call Message%ReadValue(' Enter space group number : ', io_int, 1)
   self%SGnumber = io_int(1)
   if ((self%SGnumber.lt.sgmin).or.(self%SGnumber.gt.sgmax)) then
    call Message%printMessage('Error in space group number ', frm = "(A)")
    call Message%printMessage('Crystal system / space group mismatch ', frm = "(A)")
   end if
  end do
 end if

end subroutine getSpaceGroup_

!--------------------------------------------------------------------------
recursive subroutine GetSetting_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GetSetting_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! space group first or second setting
  !!
  !! for the space groups with a second origin setting
  !! this routine asks which of the two settings to use

use mod_io

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT) :: self

integer(kind=irg)                 :: i, iset, isg, io_int(1)
type(IO_T)                        :: Message

 isg = 0
 do i=1,24
  if (tworig(i).eq.self%SGnumber) isg=i
 end do

 if (isg.ne.0) then
  iset = 0
  do while (iset.eq.0)
    call Message%printMessage(' ---------------------------------------------', frm = "(A)")
    call Message%printMessage(' This space group has two origin settings.', frm = "(A)")
    call Message%printMessage(' The first setting has site symmetry    : '//sitesym(2*isg-1), frm = "(A)")
    call Message%printMessage(' the second setting has site symmetry   : '//sitesym(2*isg), frm = "(A)")
    call Message%ReadValue(' Which setting do you wish to use (1/2) : ', io_int, 1)
    call Message%printMessage('---------------------------------------------', frm = "(A)")
    if ((io_int(1).eq.1).or.(io_int(1).eq.2)) then
      iset = io_int(1)
    else
      call Message%printMessage(' Value entered must be 1 or 2 !', frm = "(A)")
    end if
  end do
 else
  iset = 1   ! setting for space group with only one origin setting...
 end if

 self%setting = iset

end subroutine GetSetting_

!--------------------------------------------------------------------------
recursive subroutine fillgen_(self, t, isgn)
!DEC$ ATTRIBUTES DLLEXPORT :: fillgen_
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! (private) fills in a generator matrix based on an input 4-character code string

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)        :: self
character(1),INTENT(IN)                  :: t(4)
 !! 4-character input string
integer(kind=irg),INTENT(IN)             :: isgn
 !! indicates forward or reverse translation

integer(kind=irg)                        :: j
real(kind=dbl)                           :: sgn

! forward or reverse translation ?
 sgn=dble(isgn)

! first fill the array with zeroes and a 1 at 4,4
 self%c(1:4,1:4) = 0.0_dbl
 self%c(4,4) = 1.0_dbl

! then check for the particular matrix type
 select case (t(1))
  case ('a'); self%c(1,1) = 1.0_dbl; self%c(2,2) = 1.0_dbl; self%c(3,3) = 1.0_dbl
  case ('b'); self%c(1,1) =-1.0_dbl; self%c(2,2) =-1.0_dbl; self%c(3,3) = 1.0_dbl
  case ('c'); self%c(1,1) =-1.0_dbl; self%c(2,2) = 1.0_dbl; self%c(3,3) =-1.0_dbl
  case ('d'); self%c(1,3) = 1.0_dbl; self%c(2,1) = 1.0_dbl; self%c(3,2) = 1.0_dbl
  case ('e'); self%c(1,2) = 1.0_dbl; self%c(2,1) = 1.0_dbl; self%c(3,3) =-1.0_dbl
  case ('f'); self%c(1,2) =-1.0_dbl; self%c(2,1) =-1.0_dbl; self%c(3,3) =-1.0_dbl
  case ('g'); self%c(1,2) =-1.0_dbl; self%c(2,1) = 1.0_dbl; self%c(3,3) = 1.0_dbl
  case ('h'); self%c(1,1) =-1.0_dbl; self%c(2,2) =-1.0_dbl; self%c(3,3) =-1.0_dbl
  case ('i'); self%c(1,1) = 1.0_dbl; self%c(2,2) = 1.0_dbl; self%c(3,3) =-1.0_dbl
  case ('j'); self%c(1,1) = 1.0_dbl; self%c(2,2) =-1.0_dbl; self%c(3,3) = 1.0_dbl
  case ('k'); self%c(1,2) =-1.0_dbl; self%c(2,1) =-1.0_dbl; self%c(3,3) = 1.0_dbl
  case ('l'); self%c(1,2) = 1.0_dbl; self%c(2,1) = 1.0_dbl; self%c(3,3) = 1.0_dbl
  case ('m'); self%c(1,2) = 1.0_dbl; self%c(2,1) =-1.0_dbl; self%c(3,3) =-1.0_dbl
  case ('n'); self%c(1,2) =-1.0_dbl; self%c(2,1) = 1.0_dbl; self%c(2,2) =-1.0_dbl; self%c(3,3) = 1.0_dbl
 end select

! then fill in the translational component
 do j=2,4
  select case (t(j))
   case('A'); self%c(j-1,4) = sgn/6.0_dbl
   case('B'); self%c(j-1,4) = sgn/4.0_dbl
   case('C'); self%c(j-1,4) = sgn/3.0_dbl
   case('D'); self%c(j-1,4) = sgn/2.0_dbl
   case('E'); self%c(j-1,4) = sgn*2.0_dbl/3.0_dbl
   case('F'); self%c(j-1,4) = sgn*3.0_dbl/4.0_dbl
   case('G'); self%c(j-1,4) = sgn*5.0_dbl/6.0_dbl
   case('O'); self%c(j-1,4) = 0.0_dbl
   case('X'); self%c(j-1,4) = -sgn*3.0_dbl/8.0_dbl
   case('Y'); self%c(j-1,4) = -sgn/4.0_dbl
   case('Z'); self%c(j-1,4) = -sgn/8.0_dbl
  end select
 end do

end subroutine fillgen_

!--------------------------------------------------------------------------
recursive subroutine MakeGenerators_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: MakeGenerators_
  !! author: MDG
  !! version: 1.0
  !! date: 01/07/20
  !!
  !! (private) interprets the generator string and initializes all generator matrices

IMPLICIT NONE

class(SpaceGroup_T), INTENT(INOUT)  :: self

integer(kind=irg),parameter         :: QQ=48
integer(kind=irg)                   :: i,k,l,iset
real(kind=dbl)                      :: SYM_d(4,4),SYM_e(4,4)
real(kind=dbl), parameter           :: eps=0.0005_dbl
character(1)                        :: t(4)
character(40)                       :: genst

! initialize the encoded identity operator aOOO
 t = (/ 'a', 'O', 'O', 'O' /)
! compute its matrix
 call fillgen_(self,t,1)
! and put it first in the list of matrices
 self%data(1,:,:) = self%c(:,:)

! get the space group generator string
 genst = SYM_GL(self%SGnumber)

! initialize the number of generators
 self%GENnum = ichar(genst(2:2))-QQ

! create the generator matrices
 do i=2,2+self%GENnum - 1
   do k=1,4
     l=2+4*(i-2)+k
     t(k) = genst(l:l)
   end do
   call fillgen_(self,t,1)
   self%data(i,:,:) = self%c(:,:)
 end do

! this is where we are in the generator string
 i=2+4*self%GENnum+1

! if there is inversion symmetry, add the inversion to the generators
 if (genst(1:1).eq.'1') then
  self%centrosym=.TRUE.
  t = (/ 'h', 'O', 'O', 'O' /)
  call fillgen_(self,t,1)
  self%data(self%GENnum+2,:,:) = self%c(:,:)
  self%GENnum = self%GENnum+2
 else
  self%GENnum = self%GENnum+1
 end if

! now check for special origin conditions (choices 1 and 2)
 if (genst(i:i).ne.'0') then
  if (self%setting.eq.2) then
! second setting: apply translation transformation to generators
   t(1)='a'
   do k=2,4
    l=i+k-1
    t(k) = genst(l:l)
   end do
   do l=2,self%GENnum
! translate to first setting origin
    call fillgen_(self,t,-1)
    SYM_d(:,:)=self%data(l,:,:)
! apply generator
    SYM_e = matmul(SYM_d,self%c)
! translate back to second setting origin
    call fillgen_(self,t,1)
    SYM_d = matmul(self%c,SYM_e)
! reduce the translations to the fundamental unit cell
    do  k=1,3
     if (abs(SYM_d(k,4)).lt.eps) SYM_d(k,4)=0.0_dbl
     if (SYM_d(k,4).lt.0.0_dbl) then
      do while (SYM_d(k,4).lt.0.0_dbl)
       SYM_d(k,4)=SYM_d(k,4)+1.0_dbl
      end do
     end if
     if (SYM_d(k,4).ge.1.0_dbl) then
      do while (SYM_d(k,4).ge.1.0_dbl)
       SYM_d(k,4)=SYM_d(k,4)-1.0_dbl
      end do
     end if
     if (abs(SYM_d(k,4)-1.0_dbl).lt.eps) SYM_d(k,4)=0.0_dbl
    end do
! and store the result in the SYM_data array
    self%data(l,:,:)=SYM_d(:,:)
   end do
  end if ! if (setting.eq.2)
 end if ! if (genst(i:i).ne.'0')

end subroutine MakeGenerators_

!--------------------------------------------------------------------------
recursive subroutine matrixmult_(self, k1, k2)
!DEC$ ATTRIBUTES DLLEXPORT :: matrixmult_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! (private) multiplies two 4x4 symmetry matrices and brings
  !! the translation component back to the fundamental unit cell

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT) :: self

integer(kind=irg),INTENT(IN)      :: k1
 !! index of first 4x4 input matrix
integer(kind=irg),INTENT(IN)      :: k2
 !! index of second 4x4 input matrix
integer(kind=irg)                 :: i,j,k
 !! loop counters
real(kind=dbl),parameter          :: eps=0.0005_dbl
 !! truncation constant

 do i=1,4
  do j=1,4
   self%c(i,j) = 0.0_dbl
   do k=1,4
    self%c(i,j)=self%c(i,j)+self%data(k1,i,k)*self%data(k2,k,j)
   end do
  end do
 end do


! bring the translational part of the matrix back to
! the first unit cell and correct possible rounding errors
 do  k=1,3
  if (abs(self%c(k,4)).lt.eps) self%c(k,4)=0.0_dbl
  if (self%c(k,4).lt.0.0_dbl) then
   do while (self%c(k,4).lt.0.0_dbl)
    self%c(k,4)=self%c(k,4)+1.0_dbl
   end do
  end if
  if (self%c(k,4).gt.1.0_dbl) then
   do while (self%c(k,4).gt.1.0_dbl)
    self%c(k,4)=self%c(k,4)-1.0_dbl
   end do
  end if
  if (abs(self%c(k,4)-1.0_dbl).lt.eps) self%c(k,4)=0.0_dbl
 end do

end subroutine matrixmult_

!--------------------------------------------------------------------------
recursive function isitnew_(self,nsym) result(isnew)
!DEC$ ATTRIBUTES DLLEXPORT :: isitnew_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! (private) check whether or not this is a new operator by simply comparing it
  !! with all existing operators

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT) :: self

integer(kind=irg),INTENT(IN)      :: nsym
 !! index of matrix to be compared
logical                           :: isnew

integer(kind=irg)                 :: i,j,k,n
 !! loop counters
real(kind=dbl),parameter          :: eps=0.0005_dbl
 !! comparison threshold

 k=0
 n=0

 do while ((k.le.nsym).and.(n.ne.12))
  n=0
  k=k+1
  do i=1,3
   do j=1,4
    if (abs(self%c(i,j)- self%data(k,i,j)).lt.eps) n=n+1
   end do
  end do
 end do

 if (n.ne.12) then
  isnew=.TRUE.
 else
  isnew=.FALSE.
 end if

end function isitnew_

!--------------------------------------------------------------------------
recursive subroutine GenerateSymmetry_(self,dopg)
!DEC$ ATTRIBUTES DLLEXPORT :: GenerateSymmetry_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! compute all symmetry operators and store them in self%data.

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT) :: self
logical,INTENT(IN)                :: dopg
 !! logical to determine if point group matrices are to be computed as well

type(PointGroup_T)                :: PG 

integer(kind=irg)                 :: i,j,k,nsym,k1,k2,l1,l2       ! loop counters (mostly)
real(kind=dbl)                    :: q,sm                         ! auxiliary variables.

! create the space group generator matrices
 call MakeGenerators_(self)
 nsym = self%GENnum

! generate new elements from the squares of the generators
 do k=1,self%GENnum
  call matrixmult_(self,k,k)
  if (isitnew_(self,nsym).eqv..TRUE.) then
   nsym=nsym+1
   self%data(nsym,:,:) = self%c(:,:)
  end if
 end do

! generate the remainder of the factorgroup
 k1=1
 do while (k1.le.nsym)
  k2=k1+1
  do while (k2.le.nsym)
   call matrixmult_(self,k2,k1)
   if (isitnew_(self,nsym).eqv..TRUE.) then
    nsym=nsym+1
    self%data(nsym,:,:) = self%c(:,:)
    if (nsym.ge.192) then
     k2 = nsym
     k1 = nsym
    end if
   end if
   k2=k2+1
  end do
  k1=k1+1
 end do
 self%MATnum = nsym

! reduce the translation operators to the fundamental unit cell
 do i=1,self%MATnum
  do j=1,3
   self%data(i,j,4)=mod( self%data(i,j,4),1.0_dbl)
  end do
 end do

 self%recip_pending = .FALSE.

 if (dopg.eqv..TRUE.) then
! tag the point symmetry operators
! this is used to determine families of directions;
! for planes we must determine the transformed point symmetry
! operators SYM_recip() (this requires the metric tensors which we don't have at this point !!!)
  self%recip_pending = .TRUE.
  PG = PointGroup_T( self )
  self%NUMpt = PG%PGMATnum
! reciprocal space point group symmetry elements
! THIS REQUIRES THE CRYSTALLOGRAPHIC INFORMATION !!!!!!!!!!!
! we have a new routine below (fixrecipPG_) to cover this little bit...
! this method must be called after the metric tensors have been computed (calcMatrices)
! in the mod_crystallography module
 end if ! if (dopg.eq..TRUE.)
! this completes generation of the factor group

! and finally, let's determine whether or not this space group is symmorphic
 self%symmorphic = (minval(abs(SGsym - self%SGnumber)).eq.0)

end subroutine GenerateSymmetry_

!--------------------------------------------------------------------------
recursive subroutine fixRecipPG_(self, dmt, rmt)
!DEC$ ATTRIBUTES DLLEXPORT :: fixRecipPG_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! compute the reciprocal point symmetry matrices; they require the metric tensors.

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT) :: self
real(kind=dbl),INTENT(IN)         :: dmt(3,3)
 !! direct metric tensor
real(kind=dbl),INTENT(IN)         :: rmt(3,3)
 !! reciprocal metric tensor

integer(kind=irg)                 :: i, j, k, l1, l2
real(kind=dbl)                    :: q, sm

do i=1,self%NUMpt
! reciprocal space point group symmetry elements
  do j=1,3
   do k=1,3
    q=0.0_dbl
    do l1=1,3
     do l2=1,3
      q=q+dmt(j,l1)*self%direc(i,l1,l2)*rmt(l2,k)
     end do
    end do
    self%recip(i,j,k)=q
   end do
  end do
end do

self%recip_pending = .FALSE.

end subroutine fixRecipPG_

!--------------------------------------------------------------------------
recursive function getSpaceGroupName_(self) result(SGname)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupName_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group name string

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
character(11)                       :: SGname

SGname = self%SGname

end function getSpaceGroupName_

!--------------------------------------------------------------------------
recursive function getSpaceGroupOrder_(self) result(SGorder)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupOrder_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group order

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg)                   :: SGorder

SGorder = self%SGorder

end function getSpaceGroupOrder_

!--------------------------------------------------------------------------
recursive function getSpaceGroupNumber_(self) result(SGnumber)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupNumber_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group number

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg)                   :: SGnumber

SGnumber = self%SGnumber

end function getSpaceGroupNumber_

!--------------------------------------------------------------------------
recursive function getSpaceGroupSetting_(self) result(SGsetting)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupSetting_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group setting

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg)                   :: SGsetting

SGsetting = self%setting

end function getSpaceGroupSetting_

!--------------------------------------------------------------------------
recursive function getSpaceGroupGENnum_(self) result(g)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupGENnum_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! return the space group number of generators

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg)                   :: g

g = self%GENnum

end function getSpaceGroupGENnum_

!--------------------------------------------------------------------------
recursive function getSpaceGroupMATnum_(self) result(g)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupMATnum_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! return the space group number of symmetry matrices

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg)                   :: g

g = self%MATnum

end function getSpaceGroupMATnum_

!--------------------------------------------------------------------------
recursive function getSpaceGroupNUMpt_(self) result(g)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupNUMpt_
  !! author: MDG
  !! version: 1.0
  !! date: 01/13/20
  !!
  !! return the space group number of point symmetry matrices

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg)                   :: g

g = self%NUMpt

end function getSpaceGroupNUMpt_


!--------------------------------------------------------------------------
recursive function getSpaceGroupCentro_(self) result(SGcentrosym)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupCentro_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group centrosymmetry parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical                             :: SGcentrosym

SGcentrosym = self%centrosym

end function getSpaceGroupCentro_

!--------------------------------------------------------------------------
recursive function getSpaceGroupXtalSystem_(self) result(SGxtalsystem)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupXtalSystem_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group number

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg)                   :: SGxtalsystem

SGxtalsystem = self%xtal_system

end function getSpaceGroupXtalSystem_

!--------------------------------------------------------------------------
recursive function getSpaceGroupSymmorphic_(self) result(SGsymmorphic)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupSymmorphic_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group number

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical                             :: SGsymmorphic

SGsymmorphic = self%symmorphic

end function getSpaceGroupSymmorphic_

!--------------------------------------------------------------------------
recursive function getSpaceGroupDataMatrices_(self) result(SGdata)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupDataMatrices_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group data matrices

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
real(kind=dbl),allocatable          :: SGdata(:,:,:)

integer(kind=irg)                   :: sz(3)

sz = shape(self%data)

allocate(SGdata(sz(1),sz(2),sz(3)))

SGdata = self%data

end function getSpaceGroupDataMatrices_

!--------------------------------------------------------------------------
recursive function getSpaceGroupPGdirecMatrices_(self) result(SGdirec)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupPGdirecMatrices_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group direc matrices

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
real(kind=dbl),allocatable          :: SGdirec(:,:,:)

integer(kind=irg)                   :: sz(3)

sz = shape(self%direc)

allocate(SGdirec(sz(1),sz(2),sz(3)))

SGdirec = self%direc

end function getSpaceGroupPGdirecMatrices_

!--------------------------------------------------------------------------
recursive function getSpaceGroupPGrecipMatrices_(self) result(SGrecip)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupPGrecipMatrices_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! return the space group direc matrices

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
real(kind=dbl),allocatable          :: SGrecip(:,:,:)

integer(kind=irg)                   :: sz(3)

sz = shape(self%recip)

allocate(SGrecip(sz(1),sz(2),sz(3)))

SGrecip = self%recip

end function getSpaceGroupPGrecipMatrices_

!--------------------------------------------------------------------------
recursive function getSpaceGrouptrigonal_(self) result(trigonal)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGrouptrigonal_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! get the space group trigonal parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical                             :: trigonal

trigonal = self%trigonal

end function getSpaceGrouptrigonal_

!--------------------------------------------------------------------------
recursive subroutine setSpaceGroupSetting_(self, setting)
!DEC$ ATTRIBUTES DLLEXPORT :: setSpaceGroupSetting_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! set the space group reduce parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)        :: setting

self%setting = setting

end subroutine setSpaceGroupSetting_

!--------------------------------------------------------------------------
recursive subroutine setSpaceGroupNumber_(self, SGnum)
!DEC$ ATTRIBUTES DLLEXPORT :: setSpaceGroupNumber_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! set the space group number parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)        :: SGnum

self%SGnumber = SGnum

end subroutine setSpaceGroupNumber_

!--------------------------------------------------------------------------
recursive subroutine setSpaceGroupXtalSystem_(self, xs)
!DEC$ ATTRIBUTES DLLEXPORT :: setSpaceGroupXtalSystem_
  !! author: MDG
  !! version: 1.0
  !! date: 01/15/20
  !!
  !! set the space group crystal system parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
integer(kind=irg),INTENT(IN)        :: xs

self%xtal_system = xs

end subroutine setSpaceGroupXtalSystem_

!--------------------------------------------------------------------------
recursive subroutine setSpaceGroupreduce_(self, reduce)
!DEC$ ATTRIBUTES DLLEXPORT :: setSpaceGroupreduce_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! set the space group reduce parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical,INTENT(IN)                  :: reduce

self%reduce = reduce

end subroutine setSpaceGroupreduce_

!--------------------------------------------------------------------------
recursive function getSpaceGrouphexset_(self) result(hexset)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGrouphexset_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! set the space group hexset parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical                             :: hexset

hexset = self%hexset

end function getSpaceGrouphexset_

!--------------------------------------------------------------------------
recursive subroutine setSpaceGrouphexset_(self, hexset)
!DEC$ ATTRIBUTES DLLEXPORT :: setSpaceGrouphexset_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! set the space group hexset parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical,INTENT(IN)                  :: hexset

self%hexset = hexset

end subroutine setSpaceGrouphexset_

!--------------------------------------------------------------------------
recursive subroutine setSpaceGrouptrigonal_(self, trigonal)
!DEC$ ATTRIBUTES DLLEXPORT :: setSpaceGrouptrigonal_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! set the space group trigonal parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical,INTENT(IN)                  :: trigonal

self%trigonal = trigonal

end subroutine setSpaceGrouptrigonal_

!--------------------------------------------------------------------------
recursive subroutine setSpaceGroupsecond_(self, second)
!DEC$ ATTRIBUTES DLLEXPORT :: setSpaceGroupsecond_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! set the space group second parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical,INTENT(IN)                  :: second

self%second = second

end subroutine setSpaceGroupsecond_

!--------------------------------------------------------------------------
recursive function getSpaceGroupsecond_(self) result(second)
!DEC$ ATTRIBUTES DLLEXPORT :: getSpaceGroupsecond_
  !! author: MDG
  !! version: 1.0
  !! date: 01/11/20
  !!
  !! set the space group second parameter

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)   :: self
logical                             :: second

second = self%second

end function getSpaceGroupsecond_

!--------------------------------------------------------------------------
recursive subroutine ListPointGroups_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: ListPointGroups_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! list the crystallographic point groups in tabular form

use mod_io

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT) :: self

integer(kind=irg)                 :: i, j
type(IO_T)                        :: Message

call Message%printMessage('Crystallographic Point Groups')
call Message%printMessage('-----------------------------')

do i=1,32
 if (mod(i,8).eq.0) then
  write (6,"(1x,i3,':',A5,5x)") i,PGTHD(i)
 else
  write (6,"(1x,i3,':',A5,5x)",advance="no") i,PGTHD(i)
 end if
end do

end subroutine ListPointGroups_

!--------------------------------------------------------------------------
recursive subroutine CalcFamily_(self, ind, num, space, itmp)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcFamily_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! compute the indices of equivalent planes/directions and store them in the itmp array

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg),INTENT(OUT)           :: num
 !! number of equivalent entries generated
integer(kind=irg),INTENT(IN)            :: ind(3)
 !! input triplet
character(1),INTENT(IN)                 :: space
 !! 'd' or 'r'
integer(kind=irg),allocatable, INTENT(OUT) :: itmp(:,:)
 !! array used for family computations etc

integer(kind=irg)                       :: m,i,j
real(kind=sgl)                          :: h,k,l,ih,ik,il,idiff
logical                                 :: newpoint
real,parameter                          :: eps=0.0001_sgl

allocate(itmp(self%NUMpt, 3))

! first take the identity
 itmp = 0
 j=1
 itmp(j,1:3)=ind(1:3)
 h=float(ind(1))
 k=float(ind(2))
 l=float(ind(3))

! multiply with all point group elements
 do i=2,self%NUMpt
  if (space.eq.'d') then
   ih=self%direc(i,1,1)*h+self%direc(i,1,2)*k+self%direc(i,1,3)*l
   ik=self%direc(i,2,1)*h+self%direc(i,2,2)*k+self%direc(i,2,3)*l
   il=self%direc(i,3,1)*h+self%direc(i,3,2)*k+self%direc(i,3,3)*l
  else
   ih=self%recip(i,1,1)*h+self%recip(i,1,2)*k+self%recip(i,1,3)*l
   ik=self%recip(i,2,1)*h+self%recip(i,2,2)*k+self%recip(i,2,3)*l
   il=self%recip(i,3,1)*h+self%recip(i,3,2)*k+self%recip(i,3,3)*l
  end if

! is this a new point ?
  newpoint=.TRUE.
  do m=1,j+1
   idiff=(itmp(m,1)-ih)**2+(itmp(m,2)-ik)**2+(itmp(m,3)-il)**2
   if (idiff.lt.eps) newpoint=.FALSE.
  end do

  if (newpoint) then
   j=j+1
   itmp(j,1)=nint(ih)
   itmp(j,2)=nint(ik)
   itmp(j,3)=nint(il)
  endif

 end do
 num=j

end subroutine CalcFamily_

!--------------------------------------------------------------------------
recursive subroutine CalcOrbit_(self, site, n, ctmp)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcOrbit_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! compute the orbit of a point with coordinates site and return it in ctmp

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
real(kind=dbl),INTENT(IN)               :: site(3)
 !! input position
integer(kind=irg),INTENT(OUT)           :: n
 !! number of equivalent entries
real(kind=dbl),allocatable,INTENT(OUT)  :: ctmp(:,:)
 !! output array with orbit coordinates

real(kind=dbl)                          :: r(3),s(3),diff       ! auxiliary variables
real(kind=dbl), parameter               :: eps = 1.0D-4         ! comparison threshold
real(kind=dbl), parameter               :: eps2= 1.0D-6         ! comparison threshold
integer(kind=irg)                       :: i,j,k,mm             ! auxiliary variables
logical                                 :: new                  ! is this a new point ?

allocate(ctmp(self%MATnum, 3))

! get the atom coordinates
! and store them in the temporary array
 n = 1
 r = site
 ctmp(n,:)=r(:)

! get all the equivalent atom positions
 do i=2,self%MATnum
  do j=1,3
   s(j)=self%data(i,j,4)
   do k=1,3
    s(j)=s(j)+self%data(i,j,k)*r(k)
   end do
  end do

! sometimes the code below produces an incorrect answer when one of the fractional coordinates
! is slightly negative... we intercept such issues here ...
  do j=1,3
    if (abs(s(j)).lt.eps2) s(j) = 0.D0
  end do

! reduce to the fundamental unit cell if necessary
  if (self%reduce.eqv..TRUE.) then
   do j=1,3
    s(j) = mod(s(j)+100.0_dbl,1.0_dbl)
   end do
  end if

  do j=1,3
    if (abs(s(j)).lt.eps2) s(j) = 0.D0
  end do

! is this a new point ?
  new = .TRUE.
  do mm=1,n
   diff=0.0_dbl
   do j=1,3
    diff=diff+abs(ctmp(mm,j)-s(j))
   end do
   if (diff.lt.eps) then
     new = .FALSE.
   end if
  end do

! yes, it is a new point
  if (new.eqv..TRUE.) then
   n=n+1
   do j=1,3
    ctmp(n,j)=s(j)
   end do
  end if

 end do

end subroutine CalcOrbit_

!--------------------------------------------------------------------------
recursive subroutine CalcStar_(self, kk, n, stmp, space)
!DEC$ ATTRIBUTES DLLEXPORT :: CalcStar_
  !! author: MDG
  !! version: 1.0
  !! date: 01/09/20
  !!
  !! compute the star of a given reciprocal vector and outputs it in an array
  !!

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
real(kind=dbl),INTENT(IN)               :: kk(3)
 !! input vector
integer(kind=irg),INTENT(OUT)           :: n
 !! number of entries in equivalent vector array
real(kind=dbl),allocatable,INTENT(OUT)  :: stmp(:,:)
 !! output array with equivalent vectors
character(1),INTENT(IN)                 :: space
 !! 'd' or 'r'

integer(kind=irg)                       :: i,j,k,mm             ! various loop counters and such
real(kind=dbl)                          :: r(3),s(3),diff       ! auxiliary variables
real(kind=dbl),parameter                :: eps=1.0D-4           ! comparison threshold
logical                                 :: new                  ! logical (is this a new one?)

allocate(stmp(self%NUMpt,3))

 n=1
 r=kk
 stmp(n,1:3)=r(1:3)

! get all the equivalent reciprocal/direct space vectors
 do i=2,self%NUMpt
  do j=1,3
   s(j)=0.0_dbl
   do k=1,3
    if (space.eq.'r') then
     s(j)=s(j)+self%recip(i,j,k)*r(k)
    else
     s(j)=s(j)+self%direc(i,j,k)*r(k)
    end if
   end do
  end do

! is this a new point ?
  new = .TRUE.
  do mm=1,n
   diff=0.0_dbl
   do j=1,3
    diff=diff+abs(stmp(mm,j)-s(j))
   end do
   if (diff.le.eps) then
     new  = .FALSE.
   endif
  end do

! yes, it is a new point
  if (new.eqv..TRUE.) then
   n=n+1
   stmp(n,1:3)=s(1:3)
  end if

 end do

end subroutine CalcStar_

!--------------------------------------------------------------------------
function getPGnumber_(self) result(pgnum)
!DEC$ ATTRIBUTES DLLEXPORT :: getPGnumber_
  !! author: MDG
  !! version: 1.0
  !! date: 01/24/20
  !!
  !! convert the space group number to a point group number

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg)                       :: pgnum

integer(kind=irg)                       :: i

! find the rotational symmetry group number
if (self%SGnumber.ge.221) then
  i = 32
else
  i=0
  do while (SGPG(i+1).le.self%SGnumber)
    i = i+1
  end do
end if

self%PGnumber = i
pgnum = i

end function getPGnumber_

!--------------------------------------------------------------------------
subroutine setPGnumber_(self, pgnum)
!DEC$ ATTRIBUTES DLLEXPORT :: setPGnumber_
  !! author: MDG
  !! version: 1.0
  !! date: 01/25/20
  !!
  !! set the point group number

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)            :: pgnum

self%PGnumber = pgnum

end subroutine setPGnumber_

!--------------------------------------------------------------------------
recursive subroutine GetOrderinZone_(self, k, il, num, jcnt, itmp)
!DEC$ ATTRIBUTES DLLEXPORT :: GetOrderinZone_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! determine the order of the subfamily of a reciprocal
  !! lattice family belonging to a zone.

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
real(kind=sgl),INTENT(IN)               :: k(3)
 !! input vector (zone axis)
integer(kind=irg),INTENT(OUT)           :: il(48)
 !! output index list
integer(kind=irg),INTENT(IN)            :: num
 !! number of  entries to test
integer(kind=irg),INTENT(OUT)           :: jcnt
 !! number of entries in output
integer(kind=irg),INTENT(IN)            :: itmp(48,3)
 !! array used for family computations etc

integer(kind=irg)                       :: i            ! loop counter
real(kind=sgl)                          :: gn(3)        ! auxiliary variable
real(kind=sgl),parameter                :: eps=1.0E-5   ! threshold value

 jcnt = 0
 do i=1,num
  gn(1:3) = float(itmp(i,1:3))
  if (abs(sum(gn*k)).lt.eps) then
    jcnt = jcnt+1
    il(jcnt) = i
  end if
 end do

end subroutine GetOrderinZone_

!--------------------------------------------------------------------------
recursive function getLaueGroupNumber_(self) result(LGN)
!DEC$ ATTRIBUTES DLLEXPORT :: getLaueGroupNumber_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! get the Laue point group number for a given space group number

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg)                       :: LGN, i

! Oxford software uses the following symmetry conversion table:
! LG_Triclinic = 1,
! LG_Monoclinic = 2,
! LG_Orthorhombic = 3,
! LG_Tetragonal_Low = 4,
! LG_Tetragonal_High = 5,
! LG_Trigonal_Low = 6,
! LG_Trigonal_High = 7,
! LG_Hexagonal_Low = 8,
! LG_Hexagonal_High = 9,
! LG_Cubic_Low = 10,
! LG_Cubic_High = 11,
! UnknownSymmetry = 12    -> this value is not used in EMsoft
! this function returns one of the above numbers for a given space group number

! find the rotational symmetry group number
if (self%SGnumber.ge.221) then
  i = 32
else
  i=0
  do while (SGPG(i+1).le.self%SGnumber)
    i = i+1
  end do
end if

LGN = PGLaueinv(i)

end function getLaueGroupNumber_

!--------------------------------------------------------------------------
recursive function getHexvsRho_(self, pgnum) result(stnum)
!DEC$ ATTRIBUTES DLLEXPORT :: getHexvsRho_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! convert a 3D point group number to a SamplingType for trigonal symmetry

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)            :: pgnum
integer(kind=irg)                       :: stnum, sg

sg = self%SGnumber
! Is this a trigonal group
if (self%trigonal.eqv..TRUE.) then ! yes, it is

! go through all the trigonal space groups from 143 to 167 and set the correct sampling type number stnum
! point group 3
  if ((sg.ge.143).and.(sg.le.146)) stnum = 10

! point group bar3 [corrected on 7/31/18 by MDG]
  if (sg.eq.147) stnum = 12
  if (sg.eq.148) stnum = 12

! point group 32
  if ((sg.eq.149).or.(sg.eq.151).or.(sg.eq.153)) stnum = 12
  if ((sg.eq.150).or.(sg.eq.152).or.(sg.eq.154)) stnum = 12
  if (sg.eq.155) stnum = 12

! point group 3m
  if ((sg.eq.156).or.(sg.eq.158)) stnum = 14
  if ((sg.eq.157).or.(sg.eq.159)) stnum = 15
  if ((sg.eq.160).or.(sg.eq.161)) stnum = 14

! point group bar3m
  if ((sg.eq.162).or.(sg.eq.163)) stnum = 17
  if ((sg.eq.164).or.(sg.eq.165)) stnum = 16
  if ((sg.eq.166).or.(sg.eq.167)) stnum = 16
else
! this must be either point group 14 or 26, each with two settings
  if (pgnum.eq.14) then
    if ((sg.ge.115).and.(sg.le.120)) then
      stnum = 6
    else
      stnum = 8
    end if
  end if
  if (pgnum.eq.26) then
    if ((sg.eq.187).or.(sg.eq.188)) then
      stnum = 16
    else
      stnum = 17
    end if
 end if
end if

end function getHexvsRho_

!--------------------------------------------------------------------------
pure recursive function getmultiplicity_(self, t) result(stmult)
!DEC$ ATTRIBUTES DLLEXPORT :: getmultiplicity_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! convert multiplicity letter into a number string  (for Wyckoff positions)

IMPLICIT NONE

class(SpaceGroup_T),INTENT(IN)       :: self
character(1),INTENT(IN)              :: t
character(3)                         :: stmult

 select case (t(1:1))
  case ('A');  stmult = '1'
  case ('B');  stmult = '2'
  case ('C');  stmult = '3'
  case ('D');  stmult = '4'
  case ('E');  stmult = '6'
  case ('F');  stmult = '8'
  case ('G');  stmult = '9'
  case ('H');  stmult = '12'
  case ('I');  stmult = '16'
  case ('J');  stmult = '18'
  case ('K');  stmult = '24'
  case ('L');  stmult = '32'
  case ('M');  stmult = '36'
  case ('N');  stmult = '48'
  case ('O');  stmult = '64'
  case ('P');  stmult = '96'
  case ('Q');  stmult = '192'
 end select

end function getmultiplicity_

!--------------------------------------------------------------------------
pure recursive function getposition_(self, t) result(st)
!DEC$ ATTRIBUTES DLLEXPORT :: getposition_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! convert multiplicity letter into a number string (for Wyckoff positions)

IMPLICIT NONE

class(SpaceGroup_T),INTENT(IN)       :: self
character(3),INTENT(IN)              :: t
character(fnlen)                     :: st

integer(kind=irg)                    :: i

st = ' ( '
do i=1,3
 select case (t(i:i))
  case ('a');  st = trim(st)//' 0'
  case ('b');  st = trim(st)//' 1/2'
  case ('c');  st = trim(st)//' 1/4'
  case ('d');  st = trim(st)//' 3/4'
  case ('e');  st = trim(st)//' 1/6'
  case ('f');  st = trim(st)//' 1/3'
  case ('g');  st = trim(st)//' 2/3'
  case ('h');  st = trim(st)//' 5/6'
  case ('i');  st = trim(st)//' 1/8'
  case ('j');  st = trim(st)//' 3/8'
  case ('k');  st = trim(st)//' 5/8'
  case ('l');  st = trim(st)//' 7/8'
  case ('m');  st = trim(st)//' -x'
  case ('n');  st = trim(st)//' -y'
  case ('o');  st = trim(st)//' y+1/4'
  case ('p');  st = trim(st)//' x+1/4'
  case ('q');  st = trim(st)//' -y+1/4'
  case ('r');  st = trim(st)//' y+1/2'
  case ('s');  st = trim(st)//' -y+1/2'
  case ('t');  st = trim(st)//' x+1/4'
  case ('u');  st = trim(st)//' x+1/2'
  case ('v');  st = trim(st)//' -x+1/2'
  case ('w');  st = trim(st)//' 2x'
  case ('x');  st = trim(st)//' x'
  case ('y');  st = trim(st)//' y'
  case ('z');  st = trim(st)//' z'
 end select
 if (i.ne.3) then
   st = trim(st)//','
 else
   st = trim(st)//' )'
 end if
end do

end function getposition_

!--------------------------------------------------------------------------
recursive subroutine getWPstring_(self, wpstring, WPfile)
!DEC$ ATTRIBUTES DLLEXPORT :: getWPstring_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! retrieve the Wyckoff positions encoder string
  !!
  !! in the calling program, use:
  !!    WPfile = EMsoft%getConfigParameter('WyckoffPositionsfilename')

use mod_io

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
character(fnlen),INTENT(INOUT)          :: wpstring
 !! output string
character(fnlen),INTENT(IN)             :: WPfile
 !! full path to Wyckoff code file

type(IO_T)                              :: Message

character(fnlen)                        :: line
logical                                 :: fexists
integer(kind=irg)                       :: ios, i, j

inquire(file=trim(WPfile),exist=fexists)
if (.not.fexists) then
  call Message%printError('getWPstring_','Wyckoff Positions file not found')
end if

open(UNIT=dataunit,FILE=trim(WPfile), STATUS='old', FORM='formatted', ACCESS='sequential')
wpstring = ''
do i = 1, self%SGnumber
 line = ''
 read(dataunit,'(I2,A)',iostat=ios) j, line
 if (ios.ne.0) then
  exit
 end if
end do
CLOSE(UNIT=dataunit, STATUS='keep')

! wp contains the string of encoded special positions
wpstring = trim(line)

end subroutine getWPstring_

!--------------------------------------------------------------------------
recursive function interpretWyckoffletter_(self, t, x, y, z) result(st)
!DEC$ ATTRIBUTES DLLEXPORT :: interpretWyckoffletter_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! convert Wyckff letter into numerical value

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
character(1),INTENT(IN)                 :: t
real(kind=sgl),INTENT(IN)               :: x, y, z

real(kind=sgl)                          :: st

select case (t(1:1))
  case ('a');  st = 0.0
  case ('b');  st = 1.0/2.0
  case ('c');  st = 1.0/4.0
  case ('d');  st = 3.0/4.0
  case ('e');  st = 1.0/6.0
  case ('f');  st = 1.0/3.0
  case ('g');  st = 2.0/3.0
  case ('h');  st = 5.0/6.0
  case ('i');  st = 1.0/8.0
  case ('j');  st = 3.0/8.0
  case ('k');  st = 5.0/8.0
  case ('l');  st = 7.0/8.0
  case ('m');  st = -x
  case ('n');  st = -y
  case ('o');  st = y+1.0/4.0
  case ('p');  st = x+1.0/4.0
  case ('q');  st = -y+1.0/4.0
  case ('r');  st = y+1.0/2.0
  case ('s');  st = -y+1.0/2.0
  case ('t');  st = x+1.0/4.0
  case ('u');  st = x+1.0/2.0
  case ('v');  st = -x+1.0/2.0
  case ('w');  st = 2.0*x
  case ('x');  st = x
  case ('y');  st = y
  case ('z');  st = z
end select

end function interpretWyckoffletter_

!--------------------------------------------------------------------------
recursive subroutine printWyckoffPositions_(self, wpstring, WPfile, WyckoffList)
!DEC$ ATTRIBUTES DLLEXPORT :: printWyckoffPositions_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! Print a list of Wyckoff positions for the current space group
  !!
  !! in the calling program, use:
  !!    WPfile = EMsoft%getConfigParameter('WyckoffPositionsfilename')

use mod_io

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
character(fnlen),INTENT(INOUT)          :: wpstring
character(fnlen),INTENT(IN)             :: WPfile
character(6),INTENT(OUT),OPTIONAL       :: WyckoffList(27)

type(IO_T)                              :: Message
character(fnlen)                        :: line, wpf, mess, gpstring, spstring
logical                                 :: fexists
integer(kind=irg)                       :: ios, i, j, numsp, ipos, io_int(1)
character(3)                            :: gpmul, spmul, pos
character(5)                            :: splett

call getWPstring_(self, wpstring, WPfile)

! number of special positions encoded in string
numsp = (len(trim(wpstring))-1)/4

io_int(1) = self%SGnumber
call Message%WriteValue('Wyckoff positions for space group ',io_int, 1, "(I4)")
mess = '-------------------------------------'
call Message%printMessage(mess)

! next, loop through all the sets of 4 characters and interpret them
if (numsp.ne.0) then
  do i=1,numsp
    ipos = (i-1)*4+1
! get the multiplicity
    spmul = getmultiplicity_(self, wpstring(ipos:ipos))
! and interpret the position symbols
    do j=1,3
      pos(j:j) = wpstring(ipos+j:ipos+j)
    end do
    spstring = getposition_(self, pos)
! and finally the Wyckoff letter
    splett = char(96+i)
! and put them all together
    mess = trim(spmul)//trim(splett)//trim(spstring)
    call Message%printMessage(mess)
! should we store this in the optional WyckoffList variable?
    if (PRESENT(WyckoffList)) then
      WyckoffList(i) =  trim(spmul)//trim(splett)
    end if
  end do
end if

! print the general position
if (numsp.eq.26) then
  splett = 'alpha'
else
  splett = char(96+numsp+1)
end if
ipos = numsp*4+1
gpmul = getmultiplicity_(self, wpstring(ipos:ipos))
mess = trim(gpmul)//trim(splett)//' ( x, y, z )'
! should we store this in the optional WyckoffList variable?
if (PRESENT(WyckoffList)) then
  WyckoffList(numsp+1) =  trim(gpmul)//trim(splett)
end if
call Message%printMessage(mess)

end subroutine printWyckoffPositions_

!--------------------------------------------------------------------------
recursive subroutine extractWyckoffposition_(self, Wyckoffpos, pt)
!DEC$ ATTRIBUTES DLLEXPORT :: extractWyckoffposition_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! Extract the coordinates by interpreting the Wyckoff encoded string

use mod_io

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
character(3),INTENT(IN)                 :: Wyckoffpos
real(kind=sgl),INTENT(OUT)              :: pt(3)

type(IO_T)                              :: Message
integer(kind=irg)                       :: i, accum, findx, findy, findz
real(kind=sgl)                          :: io_real(3), xval, yval, zval

! first we need to figure out which coordinate values we need to ask for
accum = 0
findx = scan(Wyckoffpos,'mptuvwx')
if (findx.ne.0) accum = accum + 1
findy = scan(Wyckoffpos,'noqrsy')
if (findy.ne.0) accum = accum + 2
findz = scan(Wyckoffpos,'z')
if (findz.ne.0) accum = accum + 4

! then we do a case statement to ask for the correct values:
select case (accum)
   case (0)

   case (1)
              call Message%ReadValue(' Enter x value : ',io_real,1)
              xval = io_real(1)
   case (2)
              call Message%ReadValue(' Enter y value : ',io_real,1)
              yval = io_real(1)
   case (3)
              call Message%ReadValue(' Enter x, y values : ',io_real,2)
              xval = io_real(1)
              yval = io_real(2)
   case (4)
              call Message%ReadValue(' Enter z value : ',io_real,1)
              zval = io_real(1)
   case (5)
              call Message%ReadValue(' Enter x, z values : ',io_real,2)
              xval = io_real(1)
              zval = io_real(2)
   case (6)
              call Message%ReadValue(' Enter y, z values : ',io_real,2)
              yval = io_real(1)
              zval = io_real(2)
   case (7)
              call Message%ReadValue(' Enter x, y, z values : ',io_real,3)
              xval = io_real(1)
              yval = io_real(2)
              zval = io_real(3)
   case default
              call Message%printError('extractWyckoffposition','Unknown Wyckoff symbol')
end select

do i=1,3
  pt(i) = interpretWyckoffletter_(self, Wyckoffpos(i:i), xval, yval, zval)
end do

end subroutine extractWyckoffposition_

!--------------------------------------------------------------------------
recursive function IsGAllowed_(self, g) result(allowed)
!DEC$ ATTRIBUTES DLLEXPORT :: IsGAllowed_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! determine whether or not a given reflection is absent due to
  !! lattice centering operations.

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)            :: g(3)
 !! input reciprocal lattice vector
logical                                 :: allowed

integer(kind=irg)                       :: seo          ! auxiliary variable
character(1)                            :: lc           ! first letter of space group name

! Determine whether or not this vector is
! actually allowed by the lattice centering
 lc(1:1) =  self%SGname(2:2)
 allowed = .TRUE.
 select case (lc)
  case ('P'); ! all reflections allowed for a primitive lattice
  case ('F'); seo = sum(mod(g+100,2)); if ((seo.eq.1).or.(seo.eq.2)) allowed = .FALSE.
  case ('I'); seo = mod(sum(g)+100,2); if (seo.eq.1) allowed = .FALSE.
  case ('A'); seo = mod(g(2)+g(3)+100,2); if (seo.eq.1) allowed = .FALSE.
  case ('B'); seo = mod(g(1)+g(3)+100,2); if (seo.eq.1) allowed = .FALSE.
  case ('C'); seo = mod(g(1)+g(2)+100,2); if (seo.eq.1) allowed = .FALSE.
  case ('R'); if (self%hexset) then
               seo = mod(-g(1)+g(2)+g(3)+90,3); if (seo.ne.0) allowed = .FALSE.
              endif ! otherwise reflections are all allowed
 end select

end function IsGAllowed_

!--------------------------------------------------------------------------
recursive subroutine BFsymmetry_(self, uvw, j, isym, ir)
!DEC$ ATTRIBUTES DLLEXPORT :: BFsymmetry_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! routine to determine the symmetry of a bright field bend center [uses Table 4 in BESR paper]
  !!
  !! we must intercept the special cases where more than one symmetry of the same order can
  !! occur in a given Laue group;  e.g.  Laue group 2/m has two bright field symmetries
  !! of order 2:  2 and m
  !!
  !! This happens for the following Laue groups:
  !!   2/m     [010] -> 2    [u0w] -> m
  !!   -3m     [11.0] -> 2   [u-u.w] -> m
  !!
  !! The normal conversion from the reduced order ir to the actual Bright Field
  !! symmetry uses the PGTWDinverse array to determine the 2D symmetry

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)            :: uvw(3)
 !! zone axis indices
integer(kind=irg),INTENT(IN)            :: j
 !! index into Laue group list
integer(kind=irg),INTENT(OUT)           :: isym
 !! keeps track of special cases
integer(kind=irg),INTENT(OUT)           :: ir
 !! index of point group

integer(kind=irg)                       :: orderPG, Lauenum, ng         ! auxiliary variables
real(kind=dbl),allocatable              :: kstar(:,:)                   ! star variable


 orderPG = self%NUMpt
 Lauenum = PGLaueinv(j)
 call CalcStar_(self, dble(uvw), ng, kstar, 'd')
 ir = orderPG/ng

! take care of special cases
 isym = PGTWDinverse(ir,Lauenum)
 if ((Lauenum.eq.2).and.(ir.eq.2)) then   ! this deals with Laue Group 2/m
  if (uvw(2).eq.0) isym=3
 end if
 if ((Lauenum.eq.7).and.(ir.eq.2)) then   ! and this covers -3m
  if (uvw(1).eq.-uvw(2)) isym=3
 end if

end subroutine BFsymmetry_

!--------------------------------------------------------------------------
recursive function GetPatternSymmetry_(self, uvw, pgnum, verbose) result(dgn)
!DEC$ ATTRIBUTES DLLEXPORT :: GetPatternSymmetry_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! Determine the diffraction group number and optionally produce output

use mod_io

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)            :: uvw(3)
 !! zone axis indices
integer(kind=irg),INTENT(IN)            :: pgnum
 !! point group number
logical,INTENT(IN),OPTIONAL             :: verbose
 !! print output or not

type(IO_T)                              :: Message
integer(kind=irg)                       :: dgn          ! output diffraction group number
integer(kind=irg)                       :: io_int(3)

! get the diffraction group number
 dgn = GetDiffractionGroup_(self, uvw, pgnum)

! and print some general information if needed.
 if (present(verbose)) then
   call Message%printMessage(' ',"(A)")
   io_int(1:3) = uvw(1:3)
   call Message%WriteValue(' Zone Axis : ',io_int,3,"('[',I3,I3,I3,']')")
   io_int(1) = pgnum
   call Message%printMessage(' Crystal point group           : '//PGTHD(pgnum), frm = "(/A)")
   call Message%WriteValue(' Crystal point group number    : ', io_int, 1, "(I3)")
   call Message%printMessage(' Laue group                    : '// PGTHD(PGLaue(pgnum)), frm = "(A)")
   io_int(1) = dgn
   call Message%printMessage(' Diffraction group             : '//DG(dgn), frm = "(A)")
   call Message%WriteValue(' Diffraction group number      : ', io_int, 1, "(I3)")
   call Message%printMessage(' Projection diffraction group  : '// DG(PDG(dgn)), frm = "(A/)")

   call Message%printMessage(' Bright Field symmetry         : '//PGTWD(BFPG(dgn)), frm = "(A)")
   call Message%printMessage(' Whole Pattern symmetry        : '//PGTWD(WPPG(dgn)), frm = "(A)")
   call Message%printMessage(' Dark Field general symmetry   : '//PGTWD(DFGN(dgn)), frm = "(A)")
   call Message%printMessage(' Dark Field special symmetry   : '//PGTWD(DFSP(dgn)), frm = "(A/)")
 end if

end function GetPatternSymmetry_

!--------------------------------------------------------------------------
recursive function GetDiffractionGroup_(self, uvw, pgn) result(dgn)
!DEC$ ATTRIBUTES DLLEXPORT :: GetDiffractionGroup_
  !! author: MDG
  !! version: 1.0
  !! date: 01/10/20
  !!
  !! get the diffraction group number for this zone axis orientation
  !!
  !! This implements Table 4 in BESR76 or Table 7.2 in the EM book
  !! This is a bit tricky and we want to do this as efficiently as possible...
  !! In this first version, let's compute the order of the family of directions
  !! and based on that and the PG number determine what the diffraction group is.
  !! The data in this routine was painstakingly entered after staring at both the
  !! BESR table and the International Tables for Crystallography, point group section,
  !! Table 10.2.2.

use mod_io

IMPLICIT NONE

class(SpaceGroup_T),INTENT(INOUT)       :: self
integer(kind=irg),INTENT(IN)            :: uvw(3)
 !! zone axis indices
integer(kind=irg),INTENT(IN)            :: pgn
 !! point group number

type(IO_T)                              :: Message
real(kind=dbl),allocatable              :: kstar(:,:)                   ! star variable
integer(kind=irg)                       :: ng                           ! number of members in star
integer(kind=irg)                       :: auvw(3), mina, nz, i, s      ! abs(uvw), min(auvw), number of zeroes
integer(kind=irg)                       :: dgn                          ! (output) diffraction group number
logical                                 :: found

! compute the order of the family of directions (in direct space)
 call CalcStar_(self, dble(uvw), ng, kstar, 'd')

! determine some parameters that might be useful in deciding the correct diffraction group symmetry
 auvw = iabs(uvw)
 mina = minval(auvw)
 nz = 0   ! how many zeroes are there in the index symbol ?
 do i=1,3
   if (uvw(i).eq.0) nz = nz+1
 end do

! very long nested case statement to cover each of the 32 point groups
! we'll put them in reverse order since structures in materials
! science are more likely to have a higher symmetry...
!
! These statements were carefully entered, but it is certainly possible that
! there are errors remaining... Please report any potential errors to the author.
select case (pgn)
 case (32) ! m -3 m
     select case (ng)
         case (6)
             dgn = 19 ! 4mm1R
         case (8)
             dgn = 24 ! 6RmmR
         case (12)
             dgn = 12 ! 2mm1R
         case (24)
             dgn = 11 ! 2RmmR
         case (48)
             dgn = 4  ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (31) ! -4 3 m
     select case (ng)
         case (4)
             dgn = 23 ! 3m
         case (6)
             dgn = 18 !4RmmR
         case (12)
             if (mina.eq.0) then
                dgn = 8 ! m1R
            else
                dgn = 7 ! m
            end if
         case (24)
             if (mina.eq.0) then
                dgn = 6 ! mR
            else
                dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (30) ! 432
     select case (ng)
         case (6)
             dgn = 16 ! 4mRmR
         case (8)
             dgn = 22 ! 3mR
         case (12)
             dgn = 9  ! 2mRmR
         case (24)
              if (mina.eq.0) then
                  dgn = 6 ! mR
              else
           if ( (auvw(1).eq.auvw(2)).or.(auvw(1).eq.auvw(3)).or.(auvw(2).eq.auvw(3)) ) then
                dgn = 6 ! mR
              else
                dgn = 1 ! 1
              end if
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 !------------
 case (29) ! m3
     select case (ng)
         case (6)
             dgn = 12 ! 2mm1R
         case (8)
             dgn = 21 ! 6R
         case (12)
             dgn = 11 ! 2RmmR
         case (24)
             dgn = 4  ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (28) ! 23
     select case (ng)
         case (4)
             dgn = 20 ! 3
         case (6)
             dgn = 9  ! 2mRmR
        case (12)
              if (mina.eq.0) then
                  dgn = 6 ! mR
              else
                 dgn = 1 ! 1
             end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 !------------
 case (27) ! 6/mmm
     select case (ng)
         case (2)
             dgn = 31 ! 6mm1R
         case (6)
             dgn = 12 ! 2mm1R
         case (12)
             dgn = 11 ! 2RmmR
        case (24)
             dgn = 4  ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (26) ! -6m2
     select case (ng)
         case (2)
             dgn = 30 ! 3m1R
         case (3)
             dgn = 10 ! 2mm
         case (6)
             if (auvw(3).ne.0) then
                dgn = 7 ! m
            else  ! check if [11.0] is part of kstar
                found = .FALSE.
                do i=1,ng
                    s = sum( int(kstar(i,1:3)) - (/ 1, 1, 0 /) )
                    if (s.eq.0) found=.TRUE.
                end do
                if (found) then
                    dgn = 8 ! m1R
                else
                    dgn = 7 ! m
                end if
            end if
        case (12)
             found = .FALSE.
            do i=1,ng
                s = kstar(i,1)-kstar(i,2)
                if (s.eq.0) found=.TRUE.
            end do
            if (found) then
                dgn = 6  ! mR
            else
                dgn = 1  ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (25) ! 6mm
     select case (ng)
         case (1)
             dgn = 29 ! 6mm
         case (6)
             if (auvw(3).eq.0) then
                dgn = 8 ! m1R
            else
                 dgn = 7 ! m
            end if
         case (12)
             if (auvw(3).eq.0) then
                dgn = 6 ! mR
            else
                 dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (24) ! 622
     select case (ng)
         case (2)
             dgn = 28 ! 6mRmR
         case (6)
            dgn = 9  ! 2mRmR
         case (12)   ! there are three cases that produce mR, and one that produces 1
             found = .FALSE.
             if (auvw(3).eq.0) found=.TRUE.
             do i=1,ng  ! look for [hh.l]
                s = kstar(i,1)-kstar(i,2)
                if (s.eq.0) found=.TRUE.
            end do
             do i=1,ng  ! look for [h -h.l]
                s = kstar(i,1)+kstar(i,2)
                if (s.eq.0) found=.TRUE.
            end do
             if (found) then
                dgn = 6 ! mR
            else
                 dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 !------------
 case (23) ! 6/m
     select case (ng)
         case (2)
             dgn = 27 ! 61R
         case (6)
             dgn = 11 ! 2RmmR
         case (12)
            dgn = 4  ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (22) ! -6
     select case (ng)
         case (2)
             dgn = 26 ! 31R
         case (3)
             dgn = 7  ! m
         case (6)
            dgn = 1  ! 1
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (21) ! -6
     select case (ng)
         case (1)
             dgn = 25 ! 6
         case (6)
            if (auvw(3).eq.0) then
                dgn = 6  ! mR
            else
                dgn = 1  ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 !------------
 case (20) ! -3m
     select case (ng)
         case (2)
             dgn = 24 ! 6RmmR
         case (6)
            if (auvw(3).eq.0) then
                found = .FALSE.
                do i=1,ng
                    s = kstar(i,1)+kstar(i,2)
                    if (s.eq.0) found=.TRUE.
                end do
                if (found) then
                    dgn = 11  ! 2RmmR
                else
                    dgn = 5  ! 21R
                end if
            else
                dgn = 11  ! 2RmmR
            end if
          case (12)
            dgn = 4  ! 2R
        case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (19) ! -3m1
     select case (ng)
         case (1)
             dgn = 23 ! 3m
         case (3)
            if (auvw(3).eq.0) then
                dgn = 2  ! 1R
            else
                dgn = 7  ! m
            end if
          case (6)
            dgn = 1  ! 1
        case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (18) ! 321
     select case (ng)
         case (2)
             dgn = 22 ! 3mR
         case (3)
            dgn = 3  ! 2
          case (6)
            found = .FALSE.
            do i=1,ng
                s = kstar(i,1)+kstar(i,2)
                if (s.eq.0) found=.TRUE.
            end do
            if (found) then
                dgn = 6  ! mR
            else
                dgn = 1  ! 1
            end if
        case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 !------------
 case (17) ! -3
     select case (ng)
         case (2)
             dgn = 21 ! 6R
         case (6)
             dgn = 4  ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (16) ! 3
     select case (ng)
         case (1)
             dgn = 20 ! 3
         case (3)
             dgn = 1  ! 1
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 !------------
 case (15) ! 4/mmm
     select case (ng)
         case (2)
             dgn = 19 ! 4mm1R
         case (4)
             dgn = 12 ! 2mm1R
         case (8)
            dgn = 11 ! 2RmmR
         case (16)
             dgn = 4  ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (14) ! 4/mmm
     select case (ng)
         case (2)
             dgn = 18 ! 4RmmR
         case (4)
             if (nz.eq.2) then
                dgn = 9 ! 2mRmR
            else
                if (nz.eq.1) then
                    dgn = 8 ! m1R
                else
                    dgn = 7 ! m
                end if
            end if
         case (8)
            if (nz.eq.1) then
                dgn = 6 ! mR
            else
                dgn = 1  ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
     end select
 case (13) ! 4mmm
     select case (ng)
         case (1)
             dgn = 17 ! 4mm
         case (4)
             if ( (nz.eq.2).or. ( (nz.eq.1).and.(auvw(3).eq.0))  ) then
                dgn = 8 ! m1R
            else
                dgn = 7 ! m
            end if
         case (8)
            if (nz.eq.1) then
                dgn = 6 ! mR
            else
                dgn = 1  ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
 case (12) ! 422
     select case (ng)
         case (2)
             dgn = 16 ! 4mRmR
         case (4)
            dgn = 9  ! 2mRmR
         case (8)
             auvw = auvw - auvw(1)
            if ( (auvw(2).ne.0).and.(auvw(3).ne.0) ) then
                dgn = 1 ! 1
            else
                dgn = 6 ! mR
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
!------------
 case (11) ! 4/m
     select case (ng)
         case (2)
             dgn = 15 ! 41R
         case (4)
            dgn = 11 ! 2RmmR
         case (8)
            dgn = 4  ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
 case (10) ! -4
     select case (ng)
         case (2)
             dgn = 14 ! 4R
         case (4)
            if (auvw(3).eq.0) then
                dgn = 6 ! mR
            else
                dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
 case (9) ! 4
     select case (ng)
         case (1)
             dgn = 13 ! 4
         case (4)
            if (auvw(3).eq.0) then
                dgn = 6 ! mR
            else
                dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
!------------
 case (8) ! mmm
     select case (ng)
         case (2)
             dgn = 12 ! 2mm1R
         case (4)
             dgn = 11 ! 2RmmR
         case (8)
            dgn = 4 ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
 case (7) ! mm2
     select case (ng)
         case (1)
             dgn = 10 ! 2mm
         case (2)
             if (nz.eq.1) then
                 dgn = 8 ! m1R
             else
                dgn = 7 ! m
            end if
         case (4)
             if (auvw(3).eq.0) then
                dgn = 6 ! mR
            else
                dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
 case (6) ! 222
     select case (ng)
         case (2)
             dgn = 9 ! 2mRmR
         case (4)
             if (nz.eq.1) then
                dgn = 6 ! mR
            else
                dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
!------------
 case (5) ! 2/m
     select case (ng)
         case (2)
            if (sum(iabs( uvw - (/0, 1, 0/) )).eq.0) then
                dgn = 5  ! 21R
            else
                dgn = 11 ! 2RmmR
            end if
         case (4)
             dgn = 4 ! 2R
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
 case (4) ! m
     select case (ng)
         case (1)
            dgn = 7 ! m
         case (2)
            if (sum(iabs( uvw - (/0, 1, 0/) )).eq.0) then
                dgn = 2 ! 1R
            else
                 dgn = 1 ! 1
            end if
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
 case (3) ! 2
     select case (ng)
          case (1)
            if (sum(iabs( uvw - (/0, 1, 0/) )).eq.0) then
                dgn = 3 ! 2
            else
                dgn = 6 ! mR
            end if
         case (2)
            dgn = 1 ! 1
         case default
            call Message%printMessage(' -> incorrect number of equivalent directions in point group '//PGTHD(pgn), frm = "(A)")
    end select
!------------
 case (2) ! -1
     dgn = 4 ! 2R
 case (1) ! 1
     dgn = 1 ! 1
 case default
    call Message%printMessage('unknown point group')
    dgn = 0
end select

end function GetDiffractionGroup_

!--------------------------------------------------------------------------
recursive subroutine PGfillgen_(self, t)
!DEC$ ATTRIBUTES DLLEXPORT :: PGfillgen_
  !! author: MDG
  !! version: 1.0
  !! date: 02/20/22
  !!
  !! (private) fills in a generator matrix based on a single character code 

IMPLICIT NONE

class(PointGroup_T),INTENT(INOUT)       :: self
character(1),INTENT(IN)                 :: t

! first fill the array with zeroes 
 self%c = 0.D0

! then check for the particular matrix type
 select case (t)
  case ('a'); self%c(1,1) = 1.D0; self%c(2,2) = 1.D0; self%c(3,3) = 1.D0
  case ('b'); self%c(1,1) =-1.D0; self%c(2,2) =-1.D0; self%c(3,3) = 1.D0
  case ('c'); self%c(1,1) =-1.D0; self%c(2,2) = 1.D0; self%c(3,3) =-1.D0
  case ('d'); self%c(1,3) = 1.D0; self%c(2,1) = 1.D0; self%c(3,2) = 1.D0
  case ('e'); self%c(1,2) = 1.D0; self%c(2,1) = 1.D0; self%c(3,3) =-1.D0
  case ('f'); self%c(1,2) =-1.D0; self%c(2,1) =-1.D0; self%c(3,3) =-1.D0
  case ('g'); self%c(1,2) =-1.D0; self%c(2,1) = 1.D0; self%c(3,3) = 1.D0
  case ('h'); self%c(1,1) =-1.D0; self%c(2,2) =-1.D0; self%c(3,3) =-1.D0
  case ('i'); self%c(1,1) = 1.D0; self%c(2,2) = 1.D0; self%c(3,3) =-1.D0
  case ('j'); self%c(1,1) = 1.D0; self%c(2,2) =-1.D0; self%c(3,3) = 1.D0
  case ('k'); self%c(1,2) =-1.D0; self%c(2,1) =-1.D0; self%c(3,3) = 1.D0
  case ('l'); self%c(1,2) = 1.D0; self%c(2,1) = 1.D0; self%c(3,3) = 1.D0
  case ('m'); self%c(1,2) = 1.D0; self%c(2,1) =-1.D0; self%c(3,3) =-1.D0
  case ('n'); self%c(1,2) =-1.D0; self%c(2,1) = 1.D0; self%c(2,2) =-1.D0; self%c(3,3) = 1.D0
 end select

end subroutine PGfillgen_

!--------------------------------------------------------------------------
recursive subroutine MakePGGenerators_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: MakePGGenerators_
  !! author: MDG
  !! version: 1.0
  !! date: 02/20/22
  !!
  !! (private) interprets the generator string and initializes all generator matrices

IMPLICIT NONE

class(PointGroup_T), INTENT(INOUT)  :: self

integer(kind=irg)                   :: i
character(1)                        :: t
character(4)                        :: genst

! initialize the encoded identity operator 
t = 'a'
! compute its matrix
call PGfillgen_(self,t)
! and put it first in the list of matrices
self%direc(1,:,:) = self%c(:,:)

! get the point group generator string
genst = SYM_PGGL(self%pgnum)

! initialize the number of generators
self%GENnum = 0
do i=1,4
 if (genst(i:i).ne.'-') self%GENnum = self%GENnum+1
end do 

! create the generator matrices
if (self%GENnum.ne.0) then 
  do i=1,self%GENnum
   t = genst(i:i)
   call PGfillgen_(self,t)
   self%direc(i+1,:,:) = self%c(:,:)
  end do
end if 

! we need to include the identity as a generator
self%GENnum = self%GENnum+1

end subroutine MakePGGenerators_

!--------------------------------------------------------------------------
recursive subroutine PGmatrixmult_(self, k1, k2)
!DEC$ ATTRIBUTES DLLEXPORT :: PGmatrixmult_
  !! author: MDG
  !! version: 1.0
  !! date: 02/20/22
  !!
  !! (private) multiplies two 3x3 symmetry matrices 

IMPLICIT NONE

class(PointGroup_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)      :: k1
integer(kind=irg),INTENT(IN)      :: k2

self%c = matmul(self%direc(k1,:,:), self%direc(k2,:,:))

end subroutine PGmatrixmult_

!--------------------------------------------------------------------------
recursive function PGisitnew_(self,nsym) result(isnew)
!DEC$ ATTRIBUTES DLLEXPORT :: PGisitnew_
  !! author: MDG
  !! version: 1.0
  !! date: 02/20/22
  !!
  !! (private) check whether or not this is a new operator by simply comparing it
  !! with all existing operators

IMPLICIT NONE

class(PointGroup_T),INTENT(INOUT) :: self
integer(kind=irg),INTENT(IN)      :: nsym
logical                           :: isnew

integer(kind=irg)                 :: k,n
 !! loop counters
real(kind=dbl),parameter          :: eps=0.0005_dbl
 !! comparison threshold

k=0
n=0
isnew = .TRUE.

kloop: do while (k.le.nsym)
  k=k+1
  if (sum(abs(self%c - self%direc(k,:,:) )).lt.eps) then 
    isnew = .FALSE.
    EXIT kloop
  end if 
end do kloop

end function PGisitnew_

!--------------------------------------------------------------------------
recursive subroutine GeneratePGSymmetry_(self)
!DEC$ ATTRIBUTES DLLEXPORT :: GeneratePGSymmetry_
  !! author: MDG
  !! version: 1.0
  !! date: 02/20/22
  !!
  !! compute all symmetry operators and store them in self%direc.

IMPLICIT NONE

class(PointGroup_T),INTENT(INOUT) :: self

integer(kind=irg)                 :: i,j,k,nsym,k1,k2,l1,l2       ! loop counters (mostly)
real(kind=dbl)                    :: q,sm                         ! auxiliary variables.

! create the space group generator matrices
call self%MakePGGenerators_()
nsym = self%GENnum

! generate new elements from the squares of the generators
do k=1,self%GENnum
  call self%PGmatrixmult_(k,k)
  if (self%PGisitnew_(nsym).eqv..TRUE.) then
   nsym=nsym+1
   self%direc(nsym,:,:) = self%c(:,:)
  end if
end do

! generate the remainder of the factorgroup
 k1=1
 do while (k1.le.nsym)
  k2=k1+1
  do while (k2.le.nsym)
   call PGmatrixmult_(self,k2,k1)
   if (PGisitnew_(self,nsym).eqv..TRUE.) then
    nsym=nsym+1
    self%direc(nsym,:,:) = self%c(:,:)
    if (nsym.ge.48) then
     k2 = nsym
     k1 = nsym
    end if
   end if
   k2=k2+1
  end do
  k1=k1+1
 end do
 self%PGMATnum = nsym

end subroutine GeneratePGSymmetry_

end module mod_symmetry
