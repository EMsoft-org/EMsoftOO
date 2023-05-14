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
; EMsoft:DPmerge_annotate_triangle.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_annotate_triangle.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief special event handler for all the CW_BGROUP calls, since CW_BGROUP does not support event_pro
;
;> @date 04/07/23 MDG 1.0 first version
;--------------------------------------------------------------------------
pro DPmerge_annotate_triangle, dummy

;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common triangleparms, E, Z, xoff, yoff, alpha


for i=1,9 do begin
  i1 = fix( E*(1.0-0.1*i-0.1*(10-i)*0.5) )
  j1 = fix( z*E*0.1*(10-i) )
  i2 = E/2-(i1 - E/2)
  j2 = j1
  xyouts,xoff+i1-20,yoff+j1-5,string(10*i,format="(I2.2)"),charsize=1.2,/dev, color=0
  xyouts,xoff+i2+5,yoff+j2-5,string(100-10*i,format="(I2.2)"),charsize=1.2,/dev, color=0
  xyouts,xoff+fix(E*0.1*i)-8,yoff-15,string(10*i,format="(I2.2)"),charsize=1.2,/dev, color=0
endfor

xyouts,xoff-8,yoff-8,string(1,format="(I1.1)"),charsize=1.5,charthick=3,/dev, color=0
xyouts,xoff+fix(E)+8,yoff-8,string(2,format="(I1.1)"),charsize=1.5,charthick=3,/dev, color=0
xyouts,xoff+fix(E)/2-5,yoff+z*E,string(3,format="(I1.1)"),charsize=1.5,charthick=3,/dev, color=0


end
