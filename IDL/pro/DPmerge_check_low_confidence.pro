;
; Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
; All rights reserved.
;
; Redistribution and use in source and binary forms, with or without modification, are 
; permitted provided that the following conditions are met:
;
;     - Redistributions of source code must retain the above copyright notice, this list 
;        of conditions and the following disclaimer.
;     - Redistributions in binary form must reproduce the above copyright notice, this 
;        list of conditions and the following disclaimer in the documentation and/or 
;        other materials provided with the distribution.
;     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
;        of its contributors may be used to endorse or promote products derived from 
;        this software without specific prior written permission.
;
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
; IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
; ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
; LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
; SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
; CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
; OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
; USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
; ###################################################################
;--------------------------------------------------------------------------
; EMsoft:DPmerge_check_low_confidence.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_check_low_confidence.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief special event handler for all the CW_BGROUP calls, since CW_BGROUP does not support event_pro
;
;> @date 04/07/23 MDG 1.0 first version
;--------------------------------------------------------------------------
function DPmerge_check_low_confidence, x, y, nn

;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common triangleparms, E, Z, xoff, yoff, alpha

islow = replicate(0,nn)

dy = y
dx = x

; line 1
Ax = alpha * E
Ay = 0.0
Bx = 0.5 * alpha * E
By = alpha * E * z
p1 = ( (Bx-Ax)*(dy-Ay)-(By-Ay)*(dx-Ax) ) 
p1 = p1 lt 0.0

; line 2
Ax = (1.0-alpha) * E * 0.5
Ay = (1.0-alpha) * E * z
Bx = E - 0.5 * (1.0-alpha) * E
By = (1.0-alpha) * E * z
p2 = dy lt By

; line 3
Ax = E - alpha * E * 0.5
Ay = alpha * E * z
Bx = (1.0-alpha) * E
By = 0.0
p3 = ( (Bx-Ax)*(dy-Ay)-(By-Ay)*(dx-Ax) ) 
p3 = p3 lt 0.0

q = where ( (p1 eq 1) and (p2 eq 1) and (p3 eq 1), cnt )
if (cnt gt 0) then islow[q] = 1

return, islow
end

