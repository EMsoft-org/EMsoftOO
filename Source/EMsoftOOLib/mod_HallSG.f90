! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_HallSG
  !! author: MDG 
  !! version: 1.0 
  !! date: 11/17/22
  !!
  !! class definition for the Hall Space Group notation
  !!
  !! This class is used as an alternative to the regular space group notation
  !! for crystal structure entry.  The class can be instantiated with the -H
  !! option to the EMmkxtal program, which will then allow the user to either 
  !! select one of the pre-defined Hall symbols or enter a custom symbol.
  !! The get_Hall_SeitzGenerators routine can be used to copy the Seitz matrices
  !! into the list of generator matrices from the mod_symmetry class which then 
  !! generates the entire set of symmetry operators.

use mod_kinds
use mod_global

IMPLICIT NONE 

public :: List_Hall_Symbols, get_HallString, get_Hall_nentries 

integer(kind=irg), parameter, private  :: Hall_SGstart(230) = (/ & 
   1,   2,   3,   6,   9,  18,  21,  30,  39,  57, &
  60,  63,  72,  81,  90, 108, 109, 112, 115, 116, &
 119, 122, 123, 124, 125, 128, 134, 137, 143, 149, &
 155, 161, 164, 170, 173, 176, 182, 185, 191, 197, &
 203, 209, 212, 215, 218, 221, 227, 228, 230, 233, &
 239, 245, 251, 257, 263, 266, 269, 275, 278, 284, &
 290, 292, 298, 304, 310, 313, 316, 322, 334, 335, &
 337, 338, 341, 343, 349, 350, 351, 352, 353, 354, &
 355, 356, 357, 358, 359, 361, 363, 364, 366, 367, &
 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, &
 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, &
 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, &
 398, 399, 400, 401, 402, 404, 406, 407, 408, 410, &
 412, 413, 414, 416, 418, 419, 420, 422, 424, 425, &
 426, 428, 430, 431, 432, 433, 435, 436, 438, 439, &
 440, 441, 442, 443, 444, 446, 447, 448, 449, 450, &
 452, 454, 455, 456, 457, 458, 460, 462, 463, 464, &
 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, &
 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, &
 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, &
 495, 497, 498, 500, 501, 502, 503, 504, 505, 506, &
 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, &
 517, 518, 520, 521, 523, 524, 525, 527, 529, 530 /)
 
integer(kind=irg), parameter, private  :: Hall_SGnentries(230) = (/ & 
  1,  1,  3,  3,  9,  3,  9,  9, 18,  3, &
  3,  9,  9,  9, 18,  1,  3,  3,  1,  3, &
  3,  1,  1,  1,  3,  6,  3,  6,  6,  6, &
  6,  3,  6,  3,  3,  6,  3,  6,  6,  6, &
  6,  3,  3,  3,  3,  6,  1,  2,  3,  6, &
  6,  6,  6,  6,  3,  3,  6,  3,  6,  6, &
  2,  6,  6,  6,  3,  3,  6, 12,  1,  2, &
  1,  3,  2,  6,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  2,  2,  1,  2,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  2,  2,  1,  1,  2,  2, &
  1,  1,  2,  2,  1,  1,  2,  2,  1,  1, &
  2,  2,  1,  1,  1,  2,  1,  2,  1,  1, &
  1,  1,  1,  1,  2,  1,  1,  1,  1,  2, &
  2,  1,  1,  1,  1,  2,  2,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  2,  1,  2,  1,  1,  1,  1,  1,  1,  1, &
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, &
  1,  2,  1,  2,  1,  1,  2,  2,  1,  1 /)
 
character(8), parameter, private  :: Hall_SGlabels(530) = (/ & 
' 1      ',' 2      ',' 3:b    ',' 3:c    ',' 3:a    ',' 4:b    ',' 4:c    ',' 4:a    ',' 5:b1   ',' 5:b2   ', &
' 5:b3   ',' 5:c1   ',' 5:c2   ',' 5:c3   ',' 5:a1   ',' 5:a2   ',' 5:a3   ',' 6:b    ',' 6:c    ',' 6:a    ', &
' 7:b1   ',' 7:b2   ',' 7:b3   ',' 7:c1   ',' 7:c2   ',' 7:c3   ',' 7:a1   ',' 7:a2   ',' 7:a3   ',' 8:b1   ', &
' 8:b2   ',' 8:b3   ',' 8:c1   ',' 8:c2   ',' 8:c3   ',' 8:a1   ',' 8:a2   ',' 8:a3   ',' 9:b1   ',' 9:b2   ', &
' 9:b3   ',' 9:-b1  ',' 9:-b2  ',' 9:-b3  ',' 9:c1   ',' 9:c2   ',' 9:c3   ',' 9:-c1  ',' 9:-c2  ',' 9:-c3  ', &
' 9:a1   ',' 9:a2   ',' 9:a3   ',' 9:-a1  ',' 9:-a2  ',' 9:-a3  ','10:b    ','10:c    ','10:a    ','11:b    ', &
'11:c    ','11:a    ','12:b1   ','12:b2   ','12:b3   ','12:c1   ','12:c2   ','12:c3   ','12:a1   ','12:a2   ', &
'12:a3   ','13:b1   ','13:b2   ','13:b3   ','13:c1   ','13:c2   ','13:c3   ','13:a1   ','13:a2   ','13:a3   ', &
'14:b1   ','14:b2   ','14:b3   ','14:c1   ','14:c2   ','14:c3   ','14:a1   ','14:a2   ','14:a3   ','15:b1   ', &
'15:b2   ','15:b3   ','15:-b1  ','15:-b2  ','15:-b3  ','15:c1   ','15:c2   ','15:c3   ','15:-c1  ','15:-c2  ', &
'15:-c3  ','15:a1   ','15:a2   ','15:a3   ','15:-a1  ','15:-a2  ','15:-a3  ','16      ','17      ','17:cab  ', &
'17:bca  ','18      ','18:cab  ','18:bca  ','19      ','20      ','20:cab  ','20:bca  ','21      ','21:cab  ', &
'21:bca  ','22      ','23      ','24      ','25      ','25:cab  ','25:bca  ','26      ','26:ba-c ','26:cab  ', &
'26:-cba ','26:bca  ','26:a-cb ','27      ','27:cab  ','27:bca  ','28      ','28:ba-c ','28:cab  ','28:-cba ', &
'28:bca  ','28:a-cb ','29      ','29:ba-c ','29:cab  ','29:-cba ','29:bca  ','29:a-cb ','30      ','30:ba-c ', &
'30:cab  ','30:-cba ','30:bca  ','30:a-cb ','31      ','31:ba-c ','31:cab  ','31:-cba ','31:bca  ','31:a-cb ', &
'32      ','32:cab  ','32:bca  ','33      ','33:ba-c ','33:cab  ','33:-cba ','33:bca  ','33:a-cb ','34      ', &
'34:cab  ','34:bca  ','35      ','35:cab  ','35:bca  ','36      ','36:ba-c ','36:cab  ','36:-cba ','36:bca  ', &
'36:a-cb ','37      ','37:cab  ','37:bca  ','38      ','38:ba-c ','38:cab  ','38:-cba ','38:bca  ','38:a-cb ', &
'39      ','39:ba-c ','39:cab  ','39:-cba ','39:bca  ','39:a-cb ','40      ','40:ba-c ','40:cab  ','40:-cba ', &
'40:bca  ','40:a-cb ','41      ','41:ba-c ','41:cab  ','41:-cba ','41:bca  ','41:a-cb ','42      ','42:cab  ', &
'42:bca  ','43      ','43:cab  ','43:bca  ','44      ','44:cab  ','44:bca  ','45      ','45:cab  ','45:bca  ', &
'46      ','46:ba-c ','46:cab  ','46:-cba ','46:bca  ','46:a-cb ','47      ','48:1    ','48:2    ','49      ', &
'49:cab  ','49:bca  ','50:1    ','50:2    ','50:1cab ','50:2cab ','50:1bca ','50:2bca ','51      ','51:ba-c ', &
'51:cab  ','51:-cba ','51:bca  ','51:a-cb ','52      ','52:ba-c ','52:cab  ','52:-cba ','52:bca  ','52:a-cb ', &
'53      ','53:ba-c ','53:cab  ','53:-cba ','53:bca  ','53:a-cb ','54      ','54:ba-c ','54:cab  ','54:-cba ', &
'54:bca  ','54:a-cb ','55      ','55:cab  ','55:bca  ','56      ','56:cab  ','56:bca  ','57      ','57:ba-c ', &
'57:cab  ','57:-cba ','57:bca  ','57:a-cb ','58      ','58:cab  ','58:bca  ','59:1    ','59:2    ','59:1cab ', &
'59:2cab ','59:1bca ','59:2bca ','60      ','60:ba-c ','60:cab  ','60:-cba ','60:bca  ','60:a-cb ','61      ', &
'61:ba-c ','62      ','62:ba-c ','62:cab  ','62:-cba ','62:bca  ','62:a-cb ','63      ','63:ba-c ','63:cab  ', &
'63:-cba ','63:bca  ','63:a-cb ','64      ','64:ba-c ','64:cab  ','64:-cba ','64:bca  ','64:a-cb ','65      ', &
'65:cab  ','65:bca  ','66      ','66:cab  ','66:bca  ','67      ','67:ba-c ','67:cab  ','67:-cba ','67:bca  ', &
'67:a-cb ','68:1    ','68:2    ','68:1ba-c','68:2ba-c','68:1cab ','68:2cab ','68:1-cba','68:2-cba','68:1bca ', &
'68:2bca ','68:1a-cb','68:2a-cb','69      ','70:1    ','70:2    ','71      ','72      ','72:cab  ','72:bca  ', &
'73      ','73:ba-c ','74      ','74:ba-c ','74:cab  ','74:-cba ','74:bca  ','74:a-cb ','75      ','76      ', &
'77      ','78      ','79      ','80      ','81      ','82      ','83      ','84      ','85:1    ','85:2    ', &
'86:1    ','86:2    ','87      ','88:1    ','88:2    ','89      ','90      ','91      ','92      ','93      ', &
'94      ','95      ','96      ','97      ','98      ','99      ','100     ','101     ','102     ','103     ', &
'104     ','105     ','106     ','107     ','108     ','109     ','110     ','111     ','112     ','113     ', &
'114     ','115     ','116     ','117     ','118     ','119     ','120     ','121     ','122     ','123     ', &
'124     ','125:1   ','125:2   ','126:1   ','126:2   ','127     ','128     ','129:1   ','129:2   ','130:1   ', &
'130:2   ','131     ','132     ','133:1   ','133:2   ','134:1   ','134:2   ','135     ','136     ','137:1   ', &
'137:2   ','138:1   ','138:2   ','139     ','140     ','141:1   ','141:2   ','142:1   ','142:2   ','143     ', &
'144     ','145     ','146:H   ','146:R   ','147     ','148:H   ','148:R   ','149     ','150     ','151     ', &
'152     ','153     ','154     ','155:H   ','155:R   ','156     ','157     ','158     ','159     ','160:H   ', &
'160:R   ','161:H   ','161:R   ','162     ','163     ','164     ','165     ','166:H   ','166:R   ','167:H   ', &
'167:R   ','168     ','169     ','170     ','171     ','172     ','173     ','174     ','175     ','176     ', &
'177     ','178     ','179     ','180     ','181     ','182     ','183     ','184     ','185     ','186     ', &
'187     ','188     ','189     ','190     ','191     ','192     ','193     ','194     ','195     ','196     ', &
'197     ','198     ','199     ','200     ','201:1   ','201:2   ','202     ','203:1   ','203:2   ','204     ', &
'205     ','206     ','207     ','208     ','209     ','210     ','211     ','212     ','213     ','214     ', &
'215     ','216     ','217     ','218     ','219     ','220     ','221     ','222:1   ','222:2   ','223     ', &
'224:1   ','224:2   ','225     ','226     ','227:1   ','227:2   ','228:1   ','228:2   ','229     ','230     '/)
 
character(12), parameter, private  :: Hall_Intlabels(530) = (/ & 
'P 1         ','P -1        ','P 1 2 1     ','P 1 1 2     ','P 2 1 1     ', &
'P 1 21 1    ','P 1 1 21    ','P 21 1 1    ','C 1 2 1     ','A 1 2 1     ', &
'I 1 2 1     ','A 1 1 2     ','B 1 1 2     ','I 1 1 2     ','B 2 1 1     ', &
'C 2 1 1     ','I 2 1 1     ','P 1 m 1     ','P 1 1 m     ','P m 1 1     ', &
'P 1 c 1     ','P 1 n 1     ','P 1 a 1     ','P 1 1 a     ','P 1 1 n     ', &
'P 1 1 b     ','P b 1 1     ','P n 1 1     ','P c 1 1     ','C 1 m 1     ', &
'A 1 m 1     ','I 1 m 1     ','A 1 1 m     ','B 1 1 m     ','I 1 1 m     ', &
'B m 1 1     ','C m 1 1     ','I m 1 1     ','C 1 c 1     ','A 1 n 1     ', &
'I 1 a 1     ','A 1 a 1     ','C 1 n 1     ','I 1 c 1     ','A 1 1 a     ', &
'B 1 1 n     ','I 1 1 b     ','B 1 1 b     ','A 1 1 n     ','I 1 1 a     ', &
'B b 1 1     ','C n 1 1     ','I c 1 1     ','C c 1 1     ','B n 1 1     ', &
'I b 1 1     ','P 1 2/m 1   ','P 1 1 2/m   ','P 2/m 1 1   ','P 1 21/m 1  ', &
'P 1 1 21/m  ','P 21/m 1 1  ','C 1 2/m 1   ','A 1 2/m 1   ','I 1 2/m 1   ', &
'A 1 1 2/m   ','B 1 1 2/m   ','I 1 1 2/m   ','B 2/m 1 1   ','C 2/m 1 1   ', &
'I 2/m 1 1   ','P 1 2/c 1   ','P 1 2/n 1   ','P 1 2/a 1   ','P 1 1 2/a   ', &
'P 1 1 2/n   ','P 1 1 2/b   ','P 2/b 1 1   ','P 2/n 1 1   ','P 2/c 1 1   ', &
'P 1 21/c 1  ','P 1 21/n 1  ','P 1 21/a 1  ','P 1 1 21/a  ','P 1 1 21/n  ', &
'P 1 1 21/b  ','P 21/b 1 1  ','P 21/n 1 1  ','P 21/c 1 1  ','C 1 2/c 1   ', &
'A 1 2/n 1   ','I 1 2/a 1   ','A 1 2/a 1   ','C 1 2/n 1   ','I 1 2/c 1   ', &
'A 1 1 2/a   ','B 1 1 2/n   ','I 1 1 2/b   ','B 1 1 2/b   ','A 1 1 2/n   ', &
'I 1 1 2/a   ','B 2/b 1 1   ','C 2/n 1 1   ','I 2/c 1 1   ','C 2/c 1 1   ', &
'B 2/n 1 1   ','I 2/b 1 1   ','P 2 2 2     ','P 2 2 21    ','P 21 2 2    ', &
'P 2 21 2    ','P 21 21 2   ','P 2 21 21   ','P 21 2 21   ','P 21 21 21  ', &
'C 2 2 21    ','A 21 2 2    ','B 2 21 2    ','C 2 2 2     ','A 2 2 2     ', &
'B 2 2 2     ','F 2 2 2     ','I 2 2 2     ','I 21 21 21  ','P m m 2     ', &
'P 2 m m     ','P m 2 m     ','P m c 21    ','P c m 21    ','P 21 m a    ', &
'P 21 a m    ','P b 21 m    ','P m 21 b    ','P c c 2     ','P 2 a a     ', &
'P b 2 b     ','P m a 2     ','P b m 2     ','P 2 m b     ','P 2 c m     ', &
'P c 2 m     ','P m 2 a     ','P c a 21    ','P b c 21    ','P 21 a b    ', &
'P 21 c a    ','P c 21 b    ','P b 21 a    ','P n c 2     ','P c n 2     ', &
'P 2 n a     ','P 2 a n     ','P b 2 n     ','P n 2 b     ','P m n 21    ', &
'P n m 21    ','P 21 m n    ','P 21 n m    ','P n 21 m    ','P m 21 n    ', &
'P b a 2     ','P 2 c b     ','P c 2 a     ','P n a 21    ','P b n 21    ', &
'P 21 n b    ','P 21 c n    ','P c 21 n    ','P n 21 a    ','P n n 2     ', &
'P 2 n n     ','P n 2 n     ','C m m 2     ','A 2 m m     ','B m 2 m     ', &
'C m c 21    ','C c m 21    ','A 21 m a    ','A 21 a m    ','B b 21 m    ', &
'B m 21 b    ','C c c 2     ','A 2 a a     ','B b 2 b     ','A m m 2     ', &
'B m m 2     ','B 2 m m     ','C 2 m m     ','C m 2 m     ','A m 2 m     ', &
'A b m 2     ','B m a 2     ','B 2 c m     ','C 2 m b     ','C m 2 a     ', &
'A c 2 m     ','A m a 2     ','B b m 2     ','B 2 m b     ','C 2 c m     ', &
'C c 2 m     ','A m 2 a     ','A b a 2     ','B b a 2     ','B 2 c b     ', &
'C 2 c b     ','C c 2 a     ','A c 2 a     ','F m m 2     ','F 2 m m     ', &
'F m 2 m     ','F d d 2     ','F 2 d d     ','F d 2 d     ','I m m 2     ', &
'I 2 m m     ','I m 2 m     ','I b a 2     ','I 2 c b     ','I c 2 a     ', &
'I m a 2     ','I b m 2     ','I 2 m b     ','I 2 c m     ','I c 2 m     ', &
'I m 2 a     ','P m m m     ','P n n n:1   ','P n n n:2   ','P c c m     ', &
'P m a a     ','P b m b     ','P b a n:1   ','P b a n:2   ','P n c b:1   ', &
'P n c b:2   ','P c n a:1   ','P c n a:2   ','P m m a     ','P m m b     ', &
'P b m m     ','P c m m     ','P m c m     ','P m a m     ','P n n a     ', &
'P n n b     ','P b n n     ','P c n n     ','P n c n     ','P n a n     ', &
'P m n a     ','P n m b     ','P b m n     ','P c n m     ','P n c m     ', &
'P m a n     ','P c c a     ','P c c b     ','P b a a     ','P c a a     ', &
'P b c b     ','P b a b     ','P b a m     ','P m c b     ','P c m a     ', &
'P c c n     ','P n a a     ','P b n b     ','P b c m     ','P c a m     ', &
'P m c a     ','P m a b     ','P b m a     ','P c m b     ','P n n m     ', &
'P m n n     ','P n m n     ','P m m n:1   ','P m m n:2   ','P n m m:1   ', &
'P n m m:2   ','P m n m:1   ','P m n m:2   ','P b c n     ','P c a n     ', &
'P n c a     ','P n a b     ','P b n a     ','P c n b     ','P b c a     ', &
'P c a b     ','P n m a     ','P m n b     ','P b n m     ','P c m n     ', &
'P m c n     ','P n a m     ','C m c m     ','C c m m     ','A m m a     ', &
'A m a m     ','B b m m     ','B m m b     ','C m c a     ','C c m b     ', &
'A b m a     ','A c a m     ','B b c m     ','B m a b     ','C m m m     ', &
'A m m m     ','B m m m     ','C c c m     ','A m a a     ','B b m b     ', &
'C m m a     ','C m m b     ','A b m m     ','A c m m     ','B m c m     ', &
'B m a m     ','C c c a:1   ','C c c a:2   ','C c c b:1   ','C c c b:2   ', &
'A b a a:1   ','A b a a:2   ','A c a a:1   ','A c a a:2   ','B b c b:1   ', &
'B b c b:2   ','B b a b:1   ','B b a b:2   ','F m m m     ','F d d d:1   ', &
'F d d d:2   ','I m m m     ','I b a m     ','I m c b     ','I c m a     ', &
'I b c a     ','I c a b     ','I m m a     ','I m m b     ','I b m m     ', &
'I c m m     ','I m c m     ','I m a m     ','P 4         ','P 41        ', &
'P 42        ','P 43        ','I 4         ','I 41        ','P -4        ', &
'I -4        ','P 4/m       ','P 42/m      ','P 4/n:1     ','P 4/n:2     ', &
'P 42/n:1    ','P 42/n:2    ','I 4/m       ','I 41/a:1    ','I 41/a:2    ', &
'P 4 2 2     ','P 42 1 2    ','P 41 2 2    ','P 41 21 2   ','P 42 2 2    ', &
'P 42 21 2   ','P 43 2 2    ','P 43 21 2   ','I 4 2 2     ','I 41 2 2    ', &
'P 4 m m     ','P 4 b m     ','P 42 c m    ','P 42 n m    ','P 4 c c     ', &
'P 4 n c     ','P 42 m c    ','P 42 b c    ','I 4 m m     ','I 4 c m     ', &
'I 41 m d    ','I 41 c d    ','P -4 2 m    ','P -4 2 c    ','P -4 21 m   ', &
'P -4 21 c   ','P -4 m 2    ','P -4 c 2    ','P -4 b 2    ','P -4 n 2    ', &
'I -4 m 2    ','I -4 c 2    ','I -4 2 m    ','I -4 2 d    ','P 4/m m m   ', &
'P 4/m c c   ','P 4/n b m:1 ','P 4/n b m:2 ','P 4/n n c:1 ','P 4/n n c:2 ', &
'P 4/m b m   ','P 4/m n c   ','P 4/n m m:1 ','P 4/n m m:2 ','P 4/n c c:1 ', &
'P 4/n c c:2 ','P 42/m m c  ','P 42/m c m  ','P 42/n b c:1','P 42/n b c:2', &
'P 42/n n m:1','P 42/n n m:2','P 42/m b c  ','P 42/m n m  ','P 42/n m c:1', &
'P 42/n m c:2','P 42/n c m:1','P 42/n c m:2','I 4/m m m   ','I 4/m c m   ', &
'I 41/a m d:1','I 41/a m d:2','I 41/a c d:1','I 41/a c d:2','P 3         ', &
'P 31        ','P 32        ','R 3:H       ','R 3:R       ','P -3        ', &
'R -3:H      ','R -3:R      ','P 3 1 2     ','P 3 2 1     ','P 31 1 2    ', &
'P 31 2 1    ','P 32 1 2    ','P 32 2 1    ','R 32:H      ','R 32:R      ', &
'P 3 m 1     ','P 3 1 m     ','P 3 c 1     ','P 3 1 c     ','R 3 m:H     ', &
'R 3 m:R     ','R 3 c:H     ','R 3 c:R     ','P -3 1 m    ','P -3 1 c    ', &
'P -3 m 1    ','P -3 c 1    ','R -3 m:H    ','R -3 m:R    ','R -3 c:H    ', &
'R -3 c:R    ','P 6         ','P 61        ','P 65        ','P 62        ', &
'P 64        ','P 63        ','P -6        ','P 6/m       ','P 63/m      ', &
'P 6 2 2     ','P 61 2 2    ','P 65 2 2    ','P 62 2 2    ','P 64 2 2    ', &
'P 63 2 2    ','P 6 m m     ','P 6 c c     ','P 63 c m    ','P 63 m c    ', &
'P -6 m 2    ','P -6 c 2    ','P -6 2 m    ','P -6 2 c    ','P 6/m m m   ', &
'P 6/m c c   ','P 63/m c m  ','P 63/m m c  ','P 2 3       ','F 2 3       ', &
'I 2 3       ','P 21 3      ','I 21 3      ','P m -3      ','P n -3:1    ', &
'P n -3:2    ','F m -3      ','F d -3:1    ','F d -3:2    ','I m -3      ', &
'P a -3      ','I a -3      ','P 4 3 2     ','P 42 3 2    ','F 4 3 2     ', &
'F 41 3 2    ','I 4 3 2     ','P 43 3 2    ','P 41 3 2    ','I 41 3 2    ', &
'P -4 3 m    ','F -4 3 m    ','I -4 3 m    ','P -4 3 n    ','F -4 3 c    ', &
'I -4 3 d    ','P m -3 m    ','P n -3 n:1  ','P n -3 n:2  ','P m -3 n    ', &
'P n -3 m:1  ','P n -3 m:2  ','F m -3 m    ','F m -3 c    ','F d -3 m:1  ', &
'F d -3 m:2  ','F d -3 c:1  ','F d -3 c:2  ','I m -3 m    ','I a -3 d    '/)
 
character(16), parameter, private  :: Hall_labels(530) = (/ & 
'P 1             ','-P 1            ','P 2y            ','P 2             ','P 2x            ', &
'P 2yb           ','P 2c            ','P 2xa           ','C 2y            ','A 2y            ', &
'I 2y            ','A 2             ','B 2             ','I 2             ','B 2x            ', &
'C 2x            ','I 2x            ','P -2y           ','P -2            ','P -2x           ', &
'P -2yc          ','P -2yac         ','P -2ya          ','P -2a           ','P -2ab          ', &
'P -2b           ','P -2xb          ','P -2xbc         ','P -2xc          ','C -2y           ', &
'A -2y           ','I -2y           ','A -2            ','B -2            ','I -2            ', &
'B -2x           ','C -2x           ','I -2x           ','C -2yc          ','A -2yac         ', &
'I -2ya          ','A -2ya          ','C -2ybc         ','I -2yc          ','A -2a           ', &
'B -2bc          ','I -2b           ','B -2b           ','A -2ac          ','I -2a           ', &
'B -2xb          ','C -2xbc         ','I -2xc          ','C -2xc          ','B -2xbc         ', &
'I -2xb          ','-P 2y           ','-P 2            ','-P 2x           ','-P 2yb          ', &
'-P 2c           ','-P 2xa          ','-C 2y           ','-A 2y           ','-I 2y           ', &
'-A 2            ','-B 2            ','-I 2            ','-B 2x           ','-C 2x           ', &
'-I 2x           ','-P 2yc          ','-P 2yac         ','-P 2ya          ','-P 2a           ', &
'-P 2ab          ','-P 2b           ','-P 2xb          ','-P 2xbc         ','-P 2xc          ', &
'-P 2ybc         ','-P 2yn          ','-P 2yab         ','-P 2ac          ','-P 2n           ', &
'-P 2bc          ','-P 2xab         ','-P 2xn          ','-P 2xac         ','-C 2yc          ', &
'-A 2yac         ','-I 2ya          ','-A 2ya          ','-C 2ybc         ','-I 2yc          ', &
'-A 2a           ','-B 2bc          ','-I 2b           ','-B 2b           ','-A 2ac          ', &
'-I 2a           ','-B 2xb          ','-C 2xbc         ','-I 2xc          ','-C 2xc          ', &
'-B 2xbc         ','-I 2xb          ','P 2 2           ','P 2c 2          ','P 2a 2a         ', &
'P 2 2b          ','P 2 2ab         ','P 2bc 2         ','P 2ac 2ac       ','P 2ac 2ab       ', &
'C 2c 2          ','A 2a 2a         ','B 2 2b          ','C 2 2           ','A 2 2           ', &
'B 2 2           ','F 2 2           ','I 2 2           ','I 2b 2c         ','P 2 -2          ', &
'P -2 2          ','P -2 -2         ','P 2c -2         ','P 2c -2c        ','P -2a 2a        ', &
'P -2 2a         ','P -2 -2b        ','P -2b -2        ','P 2 -2c         ','P -2a 2         ', &
'P -2b -2b       ','P 2 -2a         ','P 2 -2b         ','P -2b 2         ','P -2c 2         ', &
'P -2c -2c       ','P -2a -2a       ','P 2c -2ac       ','P 2c -2b        ','P -2b 2a        ', &
'P -2ac 2a       ','P -2bc -2c      ','P -2a -2ab      ','P 2 -2bc        ','P 2 -2ac        ', &
'P -2ac 2        ','P -2ab 2        ','P -2ab -2ab     ','P -2bc -2bc     ','P 2ac -2        ', &
'P 2bc -2bc      ','P -2ab 2ab      ','P -2 2ac        ','P -2 -2bc       ','P -2ab -2       ', &
'P 2 -2ab        ','P -2bc 2        ','P -2ac -2ac     ','P 2c -2n        ','P 2c -2ab       ', &
'P -2bc 2a       ','P -2n 2a        ','P -2n -2ac      ','P -2ac -2n      ','P 2 -2n         ', &
'P -2n 2         ','P -2n -2n       ','C 2 -2          ','A -2 2          ','B -2 -2         ', &
'C 2c -2         ','C 2c -2c        ','A -2a 2a        ','A -2 2a         ','B -2 -2b        ', &
'B -2b -2        ','C 2 -2c         ','A -2a 2         ','B -2b -2b       ','A 2 -2          ', &
'B 2 -2          ','B -2 2          ','C -2 2          ','C -2 -2         ','A -2 -2         ', &
'A 2 -2c         ','B 2 -2c         ','B -2c 2         ','C -2b 2         ','C -2b -2b       ', &
'A -2c -2c       ','A 2 -2a         ','B 2 -2b         ','B -2b 2         ','C -2c 2         ', &
'C -2c -2c       ','A -2a -2a       ','A 2 -2ac        ','B 2 -2bc        ','B -2bc 2        ', &
'C -2bc 2        ','C -2bc -2bc     ','A -2ac -2ac     ','F 2 -2          ','F -2 2          ', &
'F -2 -2         ','F 2 -2d         ','F -2d 2         ','F -2d -2d       ','I 2 -2          ', &
'I -2 2          ','I -2 -2         ','I 2 -2c         ','I -2a 2         ','I -2b -2b       ', &
'I 2 -2a         ','I 2 -2b         ','I -2b 2         ','I -2c 2         ','I -2c -2c       ', &
'I -2a -2a       ','-P 2 2          ','P 2 2 -1n       ','-P 2ab 2bc      ','-P 2 2c         ', &
'-P 2a 2         ','-P 2b 2b        ','P 2 2 -1ab      ','-P 2ab 2b       ','P 2 2 -1bc      ', &
'-P 2b 2bc       ','P 2 2 -1ac      ','-P 2a 2c        ','-P 2a 2a        ','-P 2b 2         ', &
'-P 2 2b         ','-P 2c 2c        ','-P 2c 2         ','-P 2 2a         ','-P 2a 2bc       ', &
'-P 2b 2n        ','-P 2n 2b        ','-P 2ab 2c       ','-P 2ab 2n       ','-P 2n 2bc       ', &
'-P 2ac 2        ','-P 2bc 2bc      ','-P 2ab 2ab      ','-P 2 2ac        ','-P 2 2bc        ', &
'-P 2ab 2        ','-P 2a 2ac       ','-P 2b 2c        ','-P 2a 2b        ','-P 2ac 2c       ', &
'-P 2bc 2b       ','-P 2b 2ab       ','-P 2 2ab        ','-P 2bc 2        ','-P 2ac 2ac      ', &
'-P 2ab 2ac      ','-P 2ac 2bc      ','-P 2bc 2ab      ','-P 2c 2b        ','-P 2c 2ac       ', &
'-P 2ac 2a       ','-P 2b 2a        ','-P 2a 2ab       ','-P 2bc 2c       ','-P 2 2n         ', &
'-P 2n 2         ','-P 2n 2n        ','P 2 2ab -1ab    ','-P 2ab 2a       ','P 2bc 2 -1bc    ', &
'-P 2c 2bc       ','P 2ac 2ac -1ac  ','-P 2c 2a        ','-P 2n 2ab       ','-P 2n 2c        ', &
'-P 2a 2n        ','-P 2bc 2n       ','-P 2ac 2b       ','-P 2b 2ac       ','-P 2ac 2ab      ', &
'-P 2bc 2ac      ','-P 2ac 2n       ','-P 2bc 2a       ','-P 2c 2ab       ','-P 2n 2ac       ', &
'-P 2n 2a        ','-P 2c 2n        ','-C 2c 2         ','-C 2c 2c        ','-A 2a 2a        ', &
'-A 2 2a         ','-B 2 2b         ','-B 2b 2         ','-C 2bc 2        ','-C 2bc 2bc      ', &
'-A 2ac 2ac      ','-A 2 2ac        ','-B 2 2bc        ','-B 2bc 2        ','-C 2 2          ', &
'-A 2 2          ','-B 2 2          ','-C 2 2c         ','-A 2a 2         ','-B 2b 2b        ', &
'-C 2b 2         ','-C 2b 2b        ','-A 2c 2c        ','-A 2 2c         ','-B 2 2c         ', &
'-B 2c 2         ','C 2 2 -1bc      ','-C 2b 2bc       ','C 2 2 -1bc      ','-C 2b 2c        ', &
'A 2 2 -1ac      ','-A 2a 2c        ','A 2 2 -1ac      ','-A 2ac 2c       ','B 2 2 -1bc      ', &
'-B 2bc 2b       ','B 2 2 -1bc      ','-B 2b 2bc       ','-F 2 2          ','F 2 2 -1d       ', &
'-F 2uv 2vw      ','-I 2 2          ','-I 2 2c         ','-I 2a 2         ','-I 2b 2b        ', &
'-I 2b 2c        ','-I 2a 2b        ','-I 2b 2         ','-I 2a 2a        ','-I 2c 2c        ', &
'-I 2 2b         ','-I 2 2a         ','-I 2c 2         ','P 4             ','P 4w            ', &
'P 4c            ','P 4cw           ','I 4             ','I 4bw           ','P -4            ', &
'I -4            ','-P 4            ','-P 4c           ','P 4ab -1ab      ','-P 4a           ', &
'P 4n -1n        ','-P 4bc          ','-I 4            ','I 4bw -1bw      ','-I 4ad          ', &
'P 4 2           ','P 4ab 2ab       ','P 4w 2c         ','P 4abw 2nw      ','P 4c 2          ', &
'P 4n 2n         ','P 4cw 2c        ','P 4nw 2abw      ','I 4 2           ','I 4bw 2bw       ', &
'P 4 -2          ','P 4 -2ab        ','P 4c -2c        ','P 4n -2n        ','P 4 -2c         ', &
'P 4 -2n         ','P 4c -2         ','P 4c -2ab       ','I 4 -2          ','I 4 -2c         ', &
'I 4bw -2        ','I 4bw -2c       ','P -4 2          ','P -4 2c         ','P -4 2ab        ', &
'P -4 2n         ','P -4 -2         ','P -4 -2c        ','P -4 -2ab       ','P -4 -2n        ', &
'I -4 -2         ','I -4 -2c        ','I -4 2          ','I -4 2bw        ','-P 4 2          ', &
'-P 4 2c         ','P 4 2 -1ab      ','-P 4a 2b        ','P 4 2 -1n       ','-P 4a 2bc       ', &
'-P 4 2ab        ','-P 4 2n         ','P 4ab 2ab -1ab  ','-P 4a 2a        ','P 4ab 2n -1ab   ', &
'-P 4a 2ac       ','-P 4c 2         ','-P 4c 2c        ','P 4n 2c -1n     ','-P 4ac 2b       ', &
'P 4n 2 -1n      ','-P 4ac 2bc      ','-P 4c 2ab       ','-P 4n 2n        ','P 4n 2n -1n     ', &
'-P 4ac 2a       ','P 4n 2ab -1n    ','-P 4ac 2ac      ','-I 4 2          ','-I 4 2c         ', &
'I 4bw 2bw -1bw  ','-I 4bd 2        ','I 4bw 2aw -1bw  ','-I 4bd 2c       ','P 3             ', &
'P 31            ','P 32            ','R 3             ','P 3*            ','-P 3            ', &
'-R 3            ','-P 3*           ','P 3 2           ','P 3 2"          ','P 31 2c (0 0 1) ', &
'P 31 2"         ','P 32 2c (0 0 -1)','P 32 2"         ','R 3 2"          ','P 3* 2          ', &
'P 3 -2"         ','P 3 -2          ','P 3 -2"c        ','P 3 -2c         ','R 3 -2"         ', &
'P 3* -2         ','R 3 -2"c        ','P 3* -2n        ','-P 3 2          ','-P 3 2c         ', &
'-P 3 2"         ','-P 3 2"c        ','-R 3 2"         ','-P 3* 2         ','-R 3 2"c        ', &
'-P 3* 2n        ','P 6             ','P 61            ','P 65            ','P 62            ', &
'P 64            ','P 6c            ','P -6            ','-P 6            ','-P 6c           ', &
'P 6 2           ','P 61 2 (0 0 -1) ','P 65 2 (0 0 1)  ','P 62 2c (0 0 1) ','P 64 2c (0 0 -1)', &
'P 6c 2c         ','P 6 -2          ','P 6 -2c         ','P 6c -2         ','P 6c -2c        ', &
'P -6 2          ','P -6c 2         ','P -6 -2         ','P -6c -2c       ','-P 6 2          ', &
'-P 6 2c         ','-P 6c 2         ','-P 6c 2c        ','P 2 2 3         ','F 2 2 3         ', &
'I 2 2 3         ','P 2ac 2ab 3     ','I 2b 2c 3       ','-P 2 2 3        ','P 2 2 3 -1n     ', &
'-P 2ab 2bc 3    ','-F 2 2 3        ','F 2 2 3 -1d     ','-F 2uv 2vw 3    ','-I 2 2 3        ', &
'-P 2ac 2ab 3    ','-I 2b 2c 3      ','P 4 2 3         ','P 4n 2 3        ','F 4 2 3         ', &
'F 4d 2 3        ','I 4 2 3         ','P 4acd 2ab 3    ','P 4bd 2ab 3     ','I 4bd 2c 3      ', &
'P -4 2 3        ','F -4 2 3        ','I -4 2 3        ','P -4n 2 3       ','F -4c 2 3       ', &
'I -4bd 2c 3     ','-P 4 2 3        ','P 4 2 3 -1n     ','-P 4a 2bc 3     ','-P 4n 2 3       ', &
'P 4n 2 3 -1n    ','-P 4bc 2bc 3    ','-F 4 2 3        ','-F 4c 2 3       ','F 4d 2 3 -1d    ', &
'-F 4vw 2vw 3    ','F 4d 2 3 -1cd   ','-F 4cvw 2vw 3   ','-I 4 2 3        ','-I 4bd 2c 3     '/)
 

integer(kind=irg), parameter, public :: HallmatrixID(530) = (/ &
                0,  0,  0,  1,  2,  0,  1,  2,  0,  3,  4,  1,  5,  6,  2,  7,  8,  0,  1,  2, & 
                0,  3,  4,  1,  5,  6,  2,  7,  8,  0,  3,  4,  1,  5,  6,  2,  7,  8,  0,  3, &
                4,  3,  0,  4,  1,  5,  6,  5,  1,  6,  2,  7,  8,  7,  2,  8,  0,  1,  2,  0, &
                1,  2,  0,  3,  4,  1,  5,  6,  2,  7,  8,  0,  3,  4,  1,  5,  6,  2,  7,  8, &
                0,  3,  4,  1,  5,  6,  2,  7,  8,  0,  3,  4,  3,  0,  4,  1,  5,  6,  5,  1, &
                6,  2,  7,  8,  7,  2,  8,  0,  0,  1,  2,  0,  1,  2,  0,  0,  1,  2,  0,  1, &
                2,  0,  0,  0,  0,  1,  2,  0,  9,  1, 10,  2, 11,  0,  1,  2,  0,  9,  1, 10, &
                2, 11,  0,  9,  1, 10,  2, 11,  0,  9,  1, 10,  2, 11,  0,  9,  1, 10,  2, 11, &
                0,  1,  2,  0,  9,  1, 10,  2, 11,  0,  1,  2,  0,  1,  2,  0,  9,  1, 10,  2, & 
               11,  0,  1,  2,  0,  9,  1, 10,  2, 11,  0,  9,  1, 10,  2, 11,  0,  9,  1, 10, & 
                2, 11,  0,  9,  1, 10,  2, 11,  0,  1,  2,  0,  1,  2,  0,  1,  2,  0,  1,  2, & 
                0,  9,  1, 10,  2, 11,  0,  0,  0,  0,  1,  2,  0,  0,  1,  1,  2,  2,  0,  9, & 
                1, 10,  2, 11,  0,  9,  1, 10,  2, 11,  0,  9,  1, 10,  2, 11,  0,  9,  1, 10, &  
                2, 11,  0,  1,  2,  0,  1,  2,  0,  9,  1, 10,  2, 11,  0,  1,  2,  0,  0,  1, &
                1,  2,  2,  0,  9,  1, 10,  2, 11,  0,  9,  0,  9,  1, 10,  2, 11,  0,  9,  1, & 
               10,  2, 11,  0,  9,  1, 10,  2, 11,  0,  1,  2,  0,  1,  2,  0,  0,  1,  1,  2, &  
                2,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  0,  0,  0,  0,  0,  1,  2, &  
                0,  0,  0,  0,  1,  1,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &  
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &  
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &  
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &  
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 12,  0,  0, 12,  0,  0,  0, &  
                0,  0,  0,  0, 12,  0,  0,  0,  0,  0, 12,  0, 12,  0,  0,  0,  0,  0, 12,  0, & 
               12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &  
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &  
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)

! class definition
type, public :: HallSG_T
private 
  character(30)                     :: Hall_Symbol
  integer(kind=irg)                 :: Hall_SGnumber
  integer(kind=irg)                 :: m          ! total number of generators
  integer(kind=irg)                 :: mc         ! total number of centering generators
  integer(kind=irg)                 :: mr         ! total number of rotation generators
  integer(kind=irg)                 :: seq(4)     ! sequence of rotation orders (important for interpretation)
  character(1)                      :: axis(4)    ! axis of each rotational generator
  real(kind=dbl),allocatable,public :: SeitzGenerators(:,:,:)
  integer(kind=irg)                 :: spaces(10)
  real(kind=dbl),public             :: kvec_transform(3,3)  ! transformation matrix to apply to kvectors in master patterns

contains
private 
  procedure, pass(self) :: build_Hall_SeitzGenerators_
  procedure, pass(self) :: analyze_Hall_symbol_
  procedure, pass(self) :: get_Hall_SeitzGenerators_
  procedure, pass(self) :: print_Hall_Seitzgenerators_
  procedure, pass(self) :: get_Hall_Ngenerators_
  procedure, pass(self) :: get_NHallgenerators_
  procedure, pass(self) :: get_kvectortransform_
  procedure, pass(self) :: get_HallSGlabel_
  procedure, pass(self) :: set_HallSGnumber_

  generic, public :: get_Hall_SeitzGenerators => get_Hall_SeitzGenerators_
  generic, public :: get_NHallgenerators => get_NHallgenerators_
  generic, public :: get_HallSGlabel => get_HallSGlabel_
  generic, public :: set_HallSGnumber => set_HallSGnumber_

end type HallSG_T

! the constructor routine for this class 
interface HallSG_T
  module procedure HallSG_constructor
end interface HallSG_T

contains

!--------------------------------------------------------------------------
type(HallSG_T) function HallSG_constructor( HS, verbose ) result(HallSG)
!DEC$ ATTRIBUTES DLLEXPORT :: HallSG_constructor
!! author: MDG 
!! version: 1.0 
!! date: 11/17/22
!!
!! constructor for the HallSG_T Class; sets the symbol and generates the Seitz matrices
 
IMPLICIT NONE

character(*), INTENT(IN)        :: HS
logical,INTENT(IN),OPTIONAL     :: verbose

integer(kind=irg)               :: i, j 

! set the Hall space group symbol
HallSG%Hall_Symbol = trim(adjustl(HS))

! get the Hall space group symbol number, given the string 
j = 1
do while (j.lt.531) 
  if (HallSG%Hall_Symbol.eq.trim(Hall_labels(j))) HallSG%Hall_SGnumber = j
  j = j+1
end do 

! first we get the number of generators and the number of Seitz matrices
if (present(verbose)) then 
  if (verbose.eqv..TRUE.) call HallSG%get_Hall_Ngenerators_( verbose ) 
else
  call HallSG%get_Hall_Ngenerators_() 
end if 

! then we allocate the necessary arrays of Seitz matrices
allocate(HallSG%SeitzGenerators(4,4,HallSG%m))

! initialize all of them to be the 4x4 identity matrix
HallSG%SeitzGenerators = 0.D0
do i=1,HallSG%m
  do j=1,4 
    HallSG%SeitzGenerators(j,j,i) = 1.D0 
  end do
end do 

! and extract the generators; then we can use the existing routines in mod_symmetry to 
! find all symmetry operators
if (present(verbose)) then 
  if (verbose.eqv..TRUE.) call HallSG%build_Hall_SeitzGenerators_( verbose ) 
else
  call HallSG%build_Hall_SeitzGenerators_()
end if 

! and finally determine the transformation matrix that weill be needed for all 
! master pattern computations, to bring the irreducible part of the Kikuchi sphere 
! in the right orientation...  This can be derived from the Hall_SGlabels.
call HallSG%get_kvectortransform_()

end function HallSG_constructor

!--------------------------------------------------------------------------
subroutine HallSG_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: HallSG_destructor
!! author: MDG 
!! version: 1.0 
!! date: 11/17/22
!!
!! destructor for the HallSG_T Class
 
IMPLICIT NONE

type(HallSG_T), INTENT(INOUT)  :: self 

call reportDestructor('HallSG_T')

end subroutine HallSG_destructor

!--------------------------------------------------------------------------
function get_Hall_SeitzGenerators_(self) result(Hallgenerators)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Hall_SeitzGenerators_
!! author: MDG 
!! version: 1.0 
!! date: 11/20/22
!!
!! return an allocated array with the generator Seitz matrices
!! calling routine must allocate the properly dimensioned array

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
real(kind=dbl)                        :: Hallgenerators(4,4,self%m)

Hallgenerators = self%SeitzGenerators

end function get_Hall_SeitzGenerators_

!--------------------------------------------------------------------------
function get_NHallgenerators_(self) result( N )
!DEC$ ATTRIBUTES DLLEXPORT :: get_NHallgenerators_
!! author: MDG 
!! version: 1.0 
!! date: 11/20/22
!!
!! return the number of generators

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
integer(kind=irg)                     :: N 

N = self%m

end function get_NHallgenerators_

!--------------------------------------------------------------------------
subroutine build_Hall_SeitzGenerators_(self, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: build_Hall_SeitzGenerators_
!! author: MDG 
!! version: 1.0 
!! date: 11/18/22
!!
!! build all the generator Seitz matrices

use mod_io 

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
logical,INTENT(IN),OPTIONAL           :: verbose

type(IO_T)                            :: Message

character(30)                         :: hst, hst2, hst3, hst4 
character(1)                          :: sp = char(32)
integer(kind=irg)                     :: cgen, newgen, m, i, j
real(kind=dbl)                        :: t1(3), vmat(4,4), vmatn(4,4)
logical                               :: done 

! if there is a translation vector at the end of the Hall symbol, 
! remove it after setting the vector; we'll apply this to all the 
! generators after we create all of them.
hst = trim(self%Hall_Symbol)
t1 = 0.D0
if (scan(hst,'(').ne.0) then ! yes, there is a translation vector
  m = index(trim(hst),'(')
  if (index(trim(hst),'(0 0 1)').ne.0) t1 = (/ 0.D0, 0.D0, 1.D0/12.D0 /)
  if (index(trim(hst),'(0 0 -1)').ne.0) t1 = (/ 0.D0, 0.D0,-1.D0/12.D0 /)
  hst2 = trim(hst(1:m-1))
else
  hst2 = trim(hst)
end if 

! the first generator is always the identity, so we don't need to do anything;
! we'll start with the second one 
cgen = 2 
! extract the lattice centering information 
m = index(trim(hst2), sp)
hst3 = hst2(1:m-1)
if (scan(trim(hst3),'P').eq.0) then 
  if (scan(trim(hst3),'A').ne.0) then 
      self%SeitzGenerators(1:3,4,cgen) = (/ 0.D0, 0.5D0, 0.5D0 /)
      cgen = 3
  else if (scan(trim(hst3),'B').ne.0) then
      self%SeitzGenerators(1:3,4,cgen) = (/ 0.5D0, 0.D0, 0.5D0 /)
      cgen = 3
  else if (scan(trim(hst3),'C').ne.0) then
      self%SeitzGenerators(1:3,4,cgen) = (/ 0.5D0, 0.5D0, 0.D0 /)
      cgen = 3
  else if (scan(trim(hst3),'I').ne.0) then
      self%SeitzGenerators(1:3,4,cgen) = (/ 0.5D0, 0.5D0, 0.5D0 /)
      cgen = 3
  else if (scan(trim(hst3),'R').ne.0) then
      self%SeitzGenerators(1:3,4,cgen) = (/ 2.D0, 1.D0, 1.D0 /) / 3.D0
      cgen = 3
      self%SeitzGenerators(1:3,4,cgen) = (/ 1.D0, 2.D0, 2.D0 /) / 3.D0
      cgen = 4
  else if (scan(trim(hst3),'S').ne.0) then
      self%SeitzGenerators(1:3,4,cgen) = (/ 1.D0, 1.D0, 2.D0 /) / 3.D0
      cgen = 3
      self%SeitzGenerators(1:3,4,cgen) = (/ 2.D0, 2.D0, 1.D0 /) / 3.D0
      cgen = 4
  else if (scan(trim(hst3),'T').ne.0) then
      self%SeitzGenerators(1:3,4,cgen) = (/ 1.D0, 2.D0, 1.D0 /) / 3.D0
      cgen = 3
      self%SeitzGenerators(1:3,4,cgen) = (/ 2.D0, 1.D0, 2.D0 /) / 3.D0
      cgen = 4
  else if (scan(trim(hst3),'F').ne.0) then
      self%SeitzGenerators(1:3,4,cgen) = (/ 0.D0, 0.5D0, 0.5D0 /)
      cgen = 3
      self%SeitzGenerators(1:3,4,cgen) = (/ 0.5D0, 0.D0, 0.5D0 /)
      cgen = 4
      self%SeitzGenerators(1:3,4,cgen) = (/ 0.5D0, 0.5D0, 0.D0 /)
      cgen = 5
  else
    call Message%printError('build_Hall_SeitzGenerators_', 'unknown lattice symbol '//trim(hst3))
  end if 
! else 
!   ! we deal with the primitive space groups 
!   if (scan(trim(hst3),'-').ne.0) then
!     self%SeitzGenerators(:,:,2) = self%SeitzGenerators(:,:,1)
!     self%SeitzGenerators(1:3,1:3,2) = -self%SeitzGenerators(1:3,1:3,2)
!     cgen = 3   
!   end if 
end if 

! is there an inversion in the lattice symbol ?
if (scan(trim(hst3),'-').ne.0) then
  do i=1,cgen-1 
    newgen = cgen-1+i 
    self%SeitzGenerators(:,:,newgen) = self%SeitzGenerators(:,:,i)
    self%SeitzGenerators(1:3,1:3,newgen) = -self%SeitzGenerators(1:3,1:3,newgen)
  end do
  cgen = 2*(cgen-1)+1
end if

! do the remaining rotational symbols 
do i=cgen, self%m 
  j = i-cgen+1 
  hst3 = trim(adjustl(self%Hall_Symbol(self%spaces(j):self%spaces(j+1))))
  call self%analyze_Hall_symbol_(hst3, cgen-1+j, j)
end do 

! finally, if |t1|.ne.0, apply the origin transformation to the generators 
if (sum(abs(t1)).ne.0.D0) then 
  vmat = 0.D0 
  vmat(1,1) = 1.D0 
  vmat(2,2) = 1.D0 
  vmat(3,3) = 1.D0 
  vmat(4,4) = 1.D0 
  vmatn = vmat 
  vmat(1:3,4) = t1(1:3)
  vmatn(1:3,4) = -t1(1:3)
  do i=1,self%m
    self%SeitzGenerators(:,:,i) = matmul(vmat, matmul(self%SeitzGenerators(:,:,i), vmatn))
! bring the translation part inside the unit cell 
    do j=1,3 
      if (self%SeitzGenerators(j,4,i).lt.0.D0) self%SeitzGenerators(j,4,i) = 10.D0 + self%SeitzGenerators(j,4,i)
      if (self%SeitzGenerators(j,4,i).gt.1.D0) self%SeitzGenerators(j,4,i) = mod(self%SeitzGenerators(j,4,i),1.D0)
    end do
  end do  
end if

if (present(verbose)) then 
  if (verbose.eqv..TRUE.) call self%print_Hall_Seitzgenerators_(self%m)
end if 

end subroutine build_Hall_SeitzGenerators_

!--------------------------------------------------------------------------
subroutine analyze_Hall_symbol_(self, hst, cgen, jseq)
!DEC$ ATTRIBUTES DLLEXPORT :: analyze_Hall_symbol_
!! author: MDG 
!! version: 1.0 
!! date: 11/18/22
!!
!! print the current set of generators

use mod_io 

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
character(*), INTENT(INOUT)           :: hst 
integer(kind=irg),INTENT(IN)          :: cgen 
integer(kind=irg),INTENT(IN)          :: jseq

type(IO_T)                            :: Message
logical                               :: minus 
integer(kind=irg)                     :: ora, ax, screw 
real(kind=dbl)                        :: avec(3)
real(kind=dbl)                        :: pitch(8) = (/ 1.D0/3.D0, 2.D0/3.D0, 1.D0/4.D0, 3.D0/4.D0, &
                                                       1.D0/6.D0, 1.D0/3.D0, 2.D0/3.D0, 5.D0/6.D0 /)
character(20)                         :: stmp

! is there a minus sign ?
minus = .FALSE.
if (scan(hst,'-').eq.1) then 
  minus = .TRUE. 
! remove this symbol
  stmp = hst(2:)
  hst = adjustl(trim(stmp))
end if 

! what is the order of the rotation axis (ora)
if (scan(hst,'1').eq.1) then
  ora=1
else if (scan(hst,'2').eq.1) then
  ora=2
else if (scan(hst,'3').eq.1) then
  ora=3
else if (scan(hst,'4').eq.1) then
  ora=4
else if (scan(hst,'6').eq.1) then
  ora=6
else
  call Message%printError('analyze_Hall_symbol_',' invalid rotation axis order')
end if
self%seq(jseq) = ora ! we keep the sequence of rotation operators to assist in the 
                     ! interpretation of the simplified Hall symbols

! is this potentially a screw axis ?
screw = 0
avec = 0.D0
if (index(hst,'31').ne.0) then
  screw = 1
else if (index(hst,'32').ne.0) then
  screw = 2
else if (index(hst,'41').ne.0) then
  screw = 3
else if (index(hst,'43').ne.0) then
  screw = 4
else if (index(hst,'61').ne.0) then
  screw = 5
else if (index(hst,'62').ne.0) then
  screw = 6
else if (index(hst,'64').ne.0) then
  screw = 7
else if (index(hst,'65').ne.0) then
  screw = 8
end if

! remove this symbol
if (screw.eq.0) then 
  stmp = hst(2:)
else
  stmp = hst(3:)
end if
hst = trim(stmp)

! what is the rotation axis orientation ?
! the ordering is     [a, b, c, a+b, a+b+c] according to the simplified Hall symbol convention
! they correspond to  [x, y, z,  " ,   *  ] in the generator string
self%axis(jseq) = '0'
if (screw.ne.0) avec(3) = pitch(screw)
if (index(hst,'x').ne.0) then
  self%axis(jseq) = 'x'
  if (screw.ne.0) avec(1) = pitch(screw)
else if (index(hst,'y').ne.0) then
  self%axis(jseq) = 'y'
  if (screw.ne.0) avec(2) = pitch(screw)
else if (index(hst,'z').ne.0) then
  self%axis(jseq) = 'z'
  if (screw.ne.0) avec(3) = pitch(screw)
else if (index(hst,'''').ne.0) then 
  self%axis(jseq) = ''''
else if (index(hst,'"').ne.0) then
  self%axis(jseq) = '"'
else if (index(hst,'*').ne.0) then
  self%axis(jseq) = '*'
end if

! if there is no explicit axis symbol present, maybe it is an implicit one
! so we need to apply the ordering logic for the simplified Hall notation
! see https://scripts.iucr.org/cgi-bin/paper?a19707 for all the details 

! 1. the first rotation has an axis direction of c
! 2. the second rotation (if N is 2) has an axis direction of
!       a     if preceded by an N of 2 or 4
!       a-b   if preceded by an N of 3 or 6
! 3. the third rotation (N is always 3) has an axis direction of a+b+c

! apply rule 1 (in particular necessary for the rhombohedral/trigonal case)
if ( (jseq.eq.1) .and. (self%axis(jseq).eq.'0') ) self%axis(jseq) = 'z'

! apply rule 2.
if ( (jseq.gt.1) .and. (self%axis(jseq).eq.'0') ) then 
  if ( self%seq(jseq).eq.2) then 
    if ( (self%seq(jseq-1).eq.2) .or. (self%seq(jseq-1).eq.4)  ) then 
      self%axis(jseq) = 'x'
    else if ( (self%seq(jseq-1).eq.3) .or. (self%seq(jseq-1).eq.6) ) then
      self%axis(jseq) = ''''
    end if 
  end if 
end if 

! apply rule 3. 
if ( (self%seq(jseq).eq.3) .and. (self%axis(jseq).eq.'0') .and. (screw.eq.0) ) self%axis(jseq) = '*'

! next we need to 
! if there is no rotation axis orientation symbol present (i.e., ax=0), then we 
! have the c-axis by default
if ( (self%axis(jseq).eq.'0') .or. (self%axis(jseq).eq.'z') ) then 
  select case(ora)  ! c-axis case
    case(1) 
      ! nothing to be done since the generator is already initialized to the identity matrix
    case(2)
      self%SeitzGenerators(:,:,cgen) = reshape( (/-1.D0, 0.D0, 0.D0, 0.D0,-1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), (/3,3/) )
    case(3)
      self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 1.D0, 0.D0,-1.D0,-1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), (/3,3/) )
    case(4)
      self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 1.D0, 0.D0,-1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), (/3,3/) )
    case(6)
      self%SeitzGenerators(:,:,cgen) = reshape( (/ 1.D0, 1.D0, 0.D0,-1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), (/3,3/) )
  end select 
else if (self%axis(jseq).eq.'x') then ! a-axis case
  select case(ora)
    case(1) 
      ! nothing to be done since the generator is already initialized to the identity matrix
    case(2)
      self%SeitzGenerators(:,:,cgen) = reshape( (/1.D0, 0.D0, 0.D0, 0.D0,-1.D0, 0.D0, 0.D0, 0.D0,-1.D0 /), (/3,3/) )
    case(3)
      self%SeitzGenerators(:,:,cgen) = reshape( (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,-1.D0,-1.D0 /), (/3,3/) )
    case(4)
      self%SeitzGenerators(:,:,cgen) = reshape( (/1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0,-1.D0, 0.D0 /), (/3,3/) )
    case(6)
      self%SeitzGenerators(:,:,cgen) = reshape( (/1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 0.D0,-1.D0, 0.D0 /), (/3,3/) )
  end select 
else if (self%axis(jseq).eq.'y') then ! b-axis case
  select case(ora)
    case(1) 
      ! nothing to be done since the generator is already initialized to the identity matrix
    case(2)
      self%SeitzGenerators(:,:,cgen) = reshape( (/-1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0,-1.D0 /), (/3,3/) )
    case(3)
      self%SeitzGenerators(:,:,cgen) = reshape( (/-1.D0, 0.D0,-1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 0.D0 /), (/3,3/) )
    case(4)
      self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 0.D0,-1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 0.D0 /), (/3,3/) )
    case(6)
      self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 0.D0,-1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 1.D0 /), (/3,3/) )
  end select 
else if (self%axis(jseq).eq.'''') then 
  if (self%axis(jseq-1).eq.'x') then 
    self%SeitzGenerators(:,:,cgen) = reshape( (/-1.D0, 0.D0, 0.D0, 0.D0, 0.D0,-1.D0, 0.D0,-1.D0, 0.D0 /), (/3,3/) )
  else if (self%axis(jseq-1).eq.'y') then 
    self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 0.D0,-1.D0, 0.D0,-1.D0, 0.D0,-1.D0, 0.D0, 0.D0 /), (/3,3/) )
  else if ( (self%axis(jseq-1).eq.'0') .or. (self%axis(jseq-1).eq.'z') .or. (self%axis(jseq-1).eq.'*')  ) then 
    self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0,-1.D0, 0.D0,-1.D0, 0.D0, 0.D0, 0.D0, 0.D0,-1.D0 /), (/3,3/) )
  end if  
else if (self%axis(jseq).eq.'"') then 
  if (self%axis(jseq-1).eq.'x') then 
    self%SeitzGenerators(:,:,cgen) = reshape( (/-1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0 /), (/3,3/) )
  else if (self%axis(jseq-1).eq.'y') then 
    self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 0.D0, 1.D0, 0.D0,-1.D0, 0.D0, 1.D0, 0.D0, 0.D0 /), (/3,3/) )
  else if ( (self%axis(jseq-1).eq.'0') .or. (self%axis(jseq-1).eq.'z') ) then 
    self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0,-1.D0 /), (/3,3/) )
  end if  
else if (self%axis(jseq).eq.'*') then 
    self%SeitzGenerators(:,:,cgen) = reshape( (/ 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 1.D0, 0.D0, 0.D0 /), (/3,3/) )
end if 

! do we need to apply the inversion to this generator ?
if (minus.eqv..TRUE.) then 
  self%SeitzGenerators(1:3,1:3,cgen) = -self%SeitzGenerators(1:3,1:3,cgen)
end if 

! finally, determine the translation component unless we have a screw axis
! in which case we already have avec properly determined
if (screw.eq.0) then 
  avec = 0.D0
  if (index(hst,'a').ne.0) avec = avec + (/ 0.5D0, 0.D0, 0.D0 /)
  if (index(hst,'b').ne.0) avec = avec + (/ 0.D0, 0.5D0, 0.D0 /)
  if (index(hst,'c').ne.0) avec = avec + (/ 0.D0, 0.D0, 0.5D0 /)
  if (index(hst,'n').ne.0) avec = avec + (/ 0.5D0, 0.5D0, 0.5D0 /)
  if (index(hst,'u').ne.0) avec = avec + (/ 0.25D0, 0.D0, 0.D0 /)
  if (index(hst,'v').ne.0) avec = avec + (/ 0.D0, 0.25D0, 0.D0 /)
  if (index(hst,'w').ne.0) avec = avec + (/ 0.D0, 0.D0, 0.25D0 /)
  if (index(hst,'d').ne.0) avec = avec + (/ 0.25D0, 0.25D0, 0.25D0 /)
end if

self%SeitzGenerators(1:3,4,cgen) = avec(1:3)

end subroutine analyze_Hall_symbol_

!--------------------------------------------------------------------------
subroutine print_Hall_Seitzgenerators_(self, cnt)
!DEC$ ATTRIBUTES DLLEXPORT :: print_Hall_Seitzgenerators_
!! author: MDG 
!! version: 1.0 
!! date: 11/18/22
!!
!! print the current set of generators

use mod_io 

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
integer(kind=irg),INTENT(IN)          :: cnt

type(IO_T)                            :: Message

real(kind=dbl)                        :: io_d(4)
integer(kind=irg)                     :: i,j,io_int(1)

do i=1,cnt 
  io_int(1) = i
  call Message%WriteValue('------------- ',io_int,1)
  do j=1,4 
    io_d(1:4) = self%SeitzGenerators(j,1:4,i)
    call Message%WriteValue('',io_d,4,"(3(F10.6,','),F10.6)")
  end do 
end do

end subroutine print_Hall_Seitzgenerators_

!--------------------------------------------------------------------------
function get_Hall_nentries(SGnum) result(m)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Hall_nentries
!! author: MDG 
!! version: 1.0 
!! date: 11/21/22
!!
!! get the number of Hall groups for a given spade group

IMPLICIT NONE 

integer(kind=irg),INTENT(IN)          :: SGnum
integer(kind=irg)                     :: m 

m = Hall_SGnentries(SGnum)

end function get_Hall_nentries

!--------------------------------------------------------------------------
subroutine get_Hall_Ngenerators_(self, verbose)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Hall_Ngenerators_
!! author: MDG 
!! version: 1.0 
!! date: 11/17/22
!!
!! get the number of generators 

use mod_io 

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
logical,INTENT(IN),OPTIONAL           :: verbose

type(IO_T)                            :: Message
character(30)                         :: hst, hst2, hst3, hst4 
character(1)                          :: sp = char(32)
integer(kind=irg)                     :: m, gcnt, tcnt, io_int(1)

hst = ''
hst2 = ''

if (present(verbose)) then 
  if (verbose.eqv..TRUE.) call Message%printMessage(' Analyzing Hall symbol '//trim(self%Hall_Symbol))
endif 

! if there is a translation vector at the end of the Hall symbol, 
! remove it before determining the number of generators 
hst = trim(self%Hall_Symbol)
if (scan(hst,'(').ne.0) then ! yes, there is a translation vector
  m = index(trim(hst),'(')
  hst2 = hst(1:m-1)
else
  hst2 = hst
end if 
hst4 = hst2

! next, count the number of spaces in the symbol
! each of those corresponds to a generator
m = index(trim(hst4), sp)
gcnt = 1
self%spaces(1) = m
do while(m.gt.0)
  hst3 = ''
  hst3 = hst4(m+1:)
  if (trim(hst3) == '1') then  ! intercept P 1 and -P 1
    m = -1 
    exit 
  end if
  if (len(trim(hst3)).eq.0) exit  
  m = index(trim(hst3), sp)
  if (m.gt.0) then 
    gcnt = gcnt+1
    self%spaces(gcnt) = self%spaces(gcnt-1)+m
  end if 
  hst4 = hst3
end do 
self%spaces(gcnt+1) = len(trim(hst2))
if (m.eq.-1) then 
  gcnt = 0
  self%spaces(2:) = 0 
end if 

! then we analyze the first sub string for the centering symbol and a minus sign
m = index(hst2,' ')
hst = hst2(1:m-1)
if (scan(hst,'P').ne.0) tcnt = 1
if ( (scan(hst,'A').ne.0).or.(scan(hst,'B').ne.0).or.(scan(hst,'C').ne.0).or.(scan(hst,'I').ne.0) ) tcnt = 2
if ( (scan(hst,'R').ne.0).or.(scan(hst,'S').ne.0).or.(scan(hst,'T').ne.0) ) tcnt = 3
if (scan(hst,'F').ne.0) tcnt = 4
if (scan(hst,'-').ne.0) tcnt = 2*tcnt

! and set the total number of generators
self%m = tcnt + gcnt
self%mc = tcnt  ! keep separate track of the centering-related generators
self%mr = gcnt  ! and the rotational generators

if (present(verbose)) then 
  if (verbose.eqv..TRUE.) then
    io_int(1) = self%m 
    call Message%WriteValue('Number of generators found in '//trim(self%Hall_Symbol)//' : ',io_int,1)
  end if 
end if 

end subroutine get_Hall_Ngenerators_

!--------------------------------------------------------------------------
function List_Hall_Symbols( SGnumber, HallSGnumber, verbose ) result(HS)
!DEC$ ATTRIBUTES DLLEXPORT :: List_Hall_Symbols
!! author: MDG 
!! version: 1.0 
!! date: 11/20/22
!!
!! List all the standard Hall symbols for a given space group number 

use mod_io 

IMPLICIT NONE 

integer(kind=irg),INTENT(IN)    :: SGnumber 
integer(kind=irg),INTENT(OUT)   :: HallSGnumber 
logical,INTENT(IN),OPTIONAL     :: verbose
character(16)                   :: HS

type(IO_T)                      :: Message 

integer(kind=irg)               :: i, rd_int(1) 
character(2)                    :: n 

call Message%printMessage(' ------------------------------------')
call Message%printMessage(' This program will use the Hall Space Group symbols;')
call Message%printMessage(' please select the correct Hall Space Group entry from')
call Message%printMessage(' the list below')
call Message%printMessage(' ------------------------------------')
call Message%printMessage('   SG label Int. Symbol  Hall Symbol ')

do i=Hall_SGstart(SGnumber), Hall_SGstart(SGnumber) + Hall_SGnentries(SGnumber) - 1 
  write (n,"(I2)") i-Hall_SGstart(SGnumber)+1
  call Message%printMessage( n//' : '//Hall_SGlabels(i)//Hall_Intlabels(i)//Hall_labels(i) )
end do 

call Message%printMessage('')
call Message%ReadValue('Enter Hall space group selection : ', rd_int, 1 )

! and return the Hall Symbol for the selected space group 
HallSGnumber = Hall_SGstart(SGnumber) + rd_int(1) - 1 
HS = Hall_labels( HallSGnumber )

if (present(verbose)) then 
  if (verbose.eqv..TRUE.) then
    call Message%printMessage('')
    call Message%printMessage(' This program will use the following Hall symbol : '//trim(HS))
  end if 
end if 

end function List_Hall_Symbols

!--------------------------------------------------------------------------
function get_HallString( Hnumber, SGshort, SGsym ) result(HS)
!DEC$ ATTRIBUTES DLLEXPORT :: get_HallString
!! author: MDG 
!! version: 1.0 
!! date: 11/20/22
!!
!! return the standard Hall symbol for a given Hall space group number 

IMPLICIT NONE 

integer(kind=irg),INTENT(IN)      :: Hnumber 
character(8),INTENT(OUT),OPTIONAL :: SGshort
character(12),INTENT(OUT),OPTIONAL:: SGsym
character(16)                     :: HS

HS = Hall_labels( Hnumber )

if (present(SGshort)) SGshort = Hall_SGlabels( Hnumber )
if (present(SGsym)) SGsym = Hall_Intlabels( Hnumber )

end function get_HallString

!--------------------------------------------------------------------------
subroutine set_HallSGnumber_( self, Hnum ) 
!DEC$ ATTRIBUTES DLLEXPORT :: set_HallSGnumber_
!! author: MDG 
!! version: 1.0 
!! date: 11/23/22
!!
!!set the Hall space group number 

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
integer(kind=irg),INTENT(IN)          :: Hnum 

self%Hall_SGnumber = Hnum

end subroutine set_HallSGnumber_

!--------------------------------------------------------------------------
function get_HallSGlabel_( self, Hnum ) result(SGlabel)
!DEC$ ATTRIBUTES DLLEXPORT :: get_HallSGlabel_
!! author: MDG 
!! version: 1.0 
!! date: 11/23/22
!!
!!set the Hall space group number 

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self
integer(kind=irg),INTENT(IN)          :: Hnum 
character(8)                          :: SGlabel

SGlabel = trim(Hall_SGlabels( Hnum ))

end function get_HallSGlabel_

!--------------------------------------------------------------------------
subroutine get_kvectortransform_( self ) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_kvectortransform_
!! author: MDG 
!! version: 1.0 
!! date: 11/23/22
!!
!! determine the transformation matrix for k-vector sampling in an arbitrary
!! Hall space group...
!! 
!! In the older EMsoft code, the orthorhombic settings were encoded by:
!!  character(8), public, dimension(6):: extendedOrthsettings = (/ &
!!    " a  b  c", " b  a -c", " c  a  b", "-c  b  a", " b  c  a", " a -c  b"  /)
!! In addition to the identity, there are only 12 unique transformation matrices 
!! that take the standard setting into any of the other settings... the integer 
!! in the HallmatrixID array indicates which matrix should be used for each of the 
!! 530 Hall space group symbols.

IMPLICIT NONE 

class(HallSG_T), INTENT(INOUT)        :: self

real(kind=dbl), parameter             :: m = -1.D0, z = 0.D0, p = 1.D0, s = 1.D0/3.D0, t = 2.D0/3.D0

! initialize the matrix to the identity matrix
self%kvec_transform = reshape( (/ p, z, z, z, p, z, z, z, p /), (/3,3/) )

select case (HallmatrixID( self%Hall_SGnumber ))
  case(1)
    self%kvec_transform = reshape( (/ z, z, p, p, z, z, z, p, z /), (/3,3/) )
  case(2)
    self%kvec_transform = reshape( (/ z, p, z, z, z, p, p, z, z /), (/3,3/) )
  case(3)
    self%kvec_transform = reshape( (/ m, z, m, z, p, z, p, z, z /), (/3,3/) )
  case(4)
    self%kvec_transform = reshape( (/ z, z, p, z, p, z, m, z, m /), (/3,3/) )
  case(5)
    self%kvec_transform = reshape( (/ p, z, z, m, z, m, z, p, z /), (/3,3/) )
  case(6)
    self%kvec_transform = reshape( (/ m, z, m, z, z, p, z, p, z /), (/3,3/) )
  case(7)
    self%kvec_transform = reshape( (/ z, p, z, p, z, z, m, z, m /), (/3,3/) )
  case(8)
    self%kvec_transform = reshape( (/ z, p, z, m, z, m, z, z, p /), (/3,3/) )
  case(9)
    self%kvec_transform = reshape( (/ z, p, z, p, z, z, z, z, m /), (/3,3/) )
  case(10)
    self%kvec_transform = reshape( (/ z, z, m, z, p, z, p, z, z /), (/3,3/) )
  case(11)
    self%kvec_transform = reshape( (/ p, z, z, z, z, m, z, p, z /), (/3,3/) )
  case(12)
    self%kvec_transform = reshape( (/ t, s, s,-s, s, s,-s,-t, s /), (/3,3/) )
end select 

! self%kvec_transform = transpose(self%kvec_transform)

end subroutine get_kvectortransform_

end module mod_HallSG