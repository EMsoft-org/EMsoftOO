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
; EMsoft:DPmerge_cstrip.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_cstrip.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief create a single strip for a binary phase confidence map
;
;> @date 04/10/23 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
function DPmerge_cstrip,alpha=alpha

E = 600.0
H = 60
;

xoff = 30
yoff = 30

cimage = replicate(255B, 3, fix(E+2*xoff), fix((H+2*yoff)))

white = [255B, 255B, 255B]
for i=0,fix(E) do begin
    vR = 1.0-float(i)/E
    m = max( [vR, 1.0-vR, 0.0] )
    point = byte( 255.0 * [vR, 1.0-vR, 0.0] / m )
    for j=1,H do cimage[0:2,xoff+i,yoff+j-1] = point[0:2]
endfor

for i=1,9 do for j=1,H do cimage[0:2,xoff+i*60,yoff+j-1] = [0B,0B,0B]

if (arg_present(alpha)) then begin
  lpos = 300.0 - 3.0*alpha
  for j=1,H do cimage[0:2,xoff+lpos,yoff+j-1] = [255B,255B,255B]
  for j=1,H do cimage[0:2,xoff+600-lpos,yoff+j-1] = [255B,255B,255B]
endif

return,cimage
end
