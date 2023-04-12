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
; EMsoft:DPmerge_ctriangle.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_ctriangle.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief special event handler for all the CW_BGROUP calls, since CW_BGROUP does not support event_pro
;
;> @date 04/07/23 MDG 1.0 first version
;--------------------------------------------------------------------------
function DPmerge_ctriangle, dummy

;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common triangleparms, E, Z, xoff, yoff, alpha

cimage = replicate(255B, 3, fix(E+2*xoff), fix((E+2*xoff)*z))

white = [255B, 255B, 255B]
for i=0,fix(E) do begin
  for j=0,fix(E*z) do begin
    vB = float(j)/E/z
    vR = 1.0-float(i)/E-vB*0.5
    if (((vR+vB) le 1.0) and (vR ge 0.0)) then begin
      m = max( [vR, 1.0-vR-vB, vB] )
      point = byte( 255.0 * [vR, 1.0-vR-vB, vB] / m )
;     if (check_low_confidence([i],[j],1) eq 3) then begin
;       cimage[0:2,xoff+i,yoff+j] = white
;     end else begin
        cimage[0:2,xoff+i,yoff+j] = point
;     endelse
    endif
  endfor
endfor

tvscl,cimage,true=1

; indicate the confidence line endpoints

for i=1,9 do begin
  i1 = fix(E*0.1*i)
  j1 = 0
  j2 = fix( z*E*0.1*i )
  i2 = fix( E*(1.0-0.1*(10-i)-0.1*i*0.5) )
  plots,xoff+[i1,i2],yoff+[j1,j2],/dev,color='000000'x, thick=1.5
endfor

for i=0,9 do begin
  i1 = fix( E*(1.0-0.1*i-0.1*(10-i)*0.5) )
  j1 = fix( z*E*0.1*(10-i) )
  i2 = E/2-(i1 - E/2)
  j2 = j1
  plots,xoff+[i1,i2],yoff+[j1,j2],/dev,color='000000'x, thick=1.5
endfor

for i = 1,9 do begin
  i1 = fix(E*0.1*(10-i))
  j1 = 0
  j2 = fix( z*E*0.1*i )
  i2 = fix( E*(1.0-0.1*(10-i)-0.1*i*0.5) )
  i2 = E/2-(i2 - E/2)
  plots,xoff+[i1,i2],yoff+[j1,j2],/dev,color='000000'x, thick=1.5
endfor

if (alpha ne 0.0) then begin
  al = alpha

  i1 = fix(E*al)
  j1 = 0
  i2 = fix( E*al/2.0 )
  j2 = fix( z*E*al )
  plots,xoff+[i1,i2],yoff+[j1,j2],/dev,color='FFFFFF'x, thick=2.0

  i1 = fix( E*(1.0-al)*0.5 )
  j1 = fix( z*E*(1.0-al) )
  i2 = E-i1
  j2 = j1
  plots,xoff+[i1,i2],yoff+[j1,j2],/dev,color='FFFFFF'x, thick=2.0

  i1 = fix(E*(1.0-al))
  j1 = 0
  j2 = fix( z*E*al )
  i2 = fix( E*al*0.5 )
  i2 = E-i2
  plots,xoff+[i1,i2],yoff+[j1,j2],/dev,color='FFFFFF'x, thick=2.0
  empty
endif

cimage = tvrd(0,0,660,571,true=1)

;wdelete,10

return, cimage
end
