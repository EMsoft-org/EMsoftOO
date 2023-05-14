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
; EMsoft:DPmerge_binary.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_binary.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief phase confidence map computation for the binary case 
;
;> @date 04/10/23 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro DPmerge_binary, numnearmatches, cstrip, phasemap=phasemap, dpmap=dpmap, clev=clev

;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata

common DIdata_common, DIdata, w, h 

if (not arg_present(clev)) then clev = 0.0
;
; routine to create a binary confidence level map
;
; takes two h5 dot product files as input, along
; with a number of near matches to consider; returns
; a color binary confidence level map and, optionally,
; a color phase map.  The clevlines option adds white 
; confidence level line to the plot; omitted if absent.
;
; This routine makes all the drawings in the Z buffer device

numx = DPmergedata.ipf_wd
numy = DPmergedata.ipf_ht
nn = numx*numy
numm = numnearmatches

Atdp = DIdata[0].TDP
Btdp = DIdata[1].TDP

if arg_present(phasemap) then begin
  phasemap = bytarr(3,numx,numy)
endif

if arg_present(dpmap) then begin 
  dpmap = fltarr(nn)
; and fill the dotproductmap with the highest values for the two phases
  for i=0,nn-1 do begin
    m = [Atdp[0,i],Btdp[0,i]]
    dpmap[i] = max(m)
  endfor
  dpmap = congrid(reverse(reform(dpmap, numx, numy),2),2*numx,2*numy)
endif

; truncate the array and scale the dot products by the weight factors
Atdp = Atdp[0:numm-1,0:nn-1] * DPmergedata.wf[0]
Btdp = Btdp[0:numm-1,0:nn-1] * DPmergedata.wf[1]

MA = fltarr(nn)
MB = fltarr(nn)

; figure out the intensity scaling via color saturation
mindp = min([Atdp, Btdp], max=maxdp)

Atdp = (Atdp-mindp)/(maxdp-mindp)
Btdp = (Btdp-mindp)/(maxdp-mindp)

line = [replicate(1,numm),replicate(2,numm)]

for i=0,nn-1 do begin
  dp = [Atdp[*,i],Btdp[*,i]]
  s = reverse(sort(dp))
  newline = line[s]
  newline = newline[0:numm-1]
  dp = dp[s]
  denominator = total(dp[0:numm-1])
  q = where(newline eq 1,cnt)
  if (cnt gt 0.0) then MA[i] = total(dp[q])/denominator else MA[i] = 0.0
  q = where(newline eq 2,cnt)
  if (cnt gt 0.0) then MB[i] = total(dp[q])/denominator else MB[i] = 0.0
endfor

MA= reform(MA,numx,numy)
MB = reform(MB,numx,numy)

if arg_present(phasemap) then begin
  CImap = bytarr(3,numx,numy)
  for i=0,numx-1 do begin
    for j=0,numy-1 do begin
      m = max([MA[i,j],MB[i,j]])
      CImap[0:2,i,j] = byte(255.0*[MA[i,j],MB[i,j],0.0] / m )
    endfor
  endfor
endif

; next generate the map in the Z-buffer
window,10,xsi=660,ysi=120,retain=2,/pixmap
erase

E = 600
H = 60

xoff = 30
yoff = 30 + H/2

if arg_present(clev) then begin
  bg = DPmerge_cstrip(alpha=clev)
end else begin
  bg = DPmerge_cstrip()
end

bgdim = size(bg,/dimensions)
tvscl,bg, true=1

; transform the data points to linear coordinates
x = reform(mB,nn)*E

; generate a histogram for each of the potential points inside the triangle
np = 31

gridx = fltarr(np)
k = 0
dx = E/float(np-1)
for j=0,np-1 do begin
    gridx[k] = j*dx
    k += 1
endfor 

histn = histogram(x,min=0.0,max=E,binsize=dx)

h = 2.5*fix(alog10(histn+1.0))

nt = 360.
th = findgen(nt)*!dtor
ct = cos(th)
st = sin(th)

for k=0,np-1 do if (h[k] gt 0.001) then polyfill,xoff+gridx[k]+h[k]*ct,yoff+h[k]*st,/dev,color=0

cstrip = tvrd(0,0,bgdim[1],bgdim[2],TRUE=1)

if (clev ne 0.0) then begin
  lpos = 300.0 - 3.0*clev
  q = where((x ge lpos) and (x le (E-lpos)), cnt)
  if (cnt gt 0) then begin
    for i=0,2 do begin
      slice = reform(CImap[i,*,*])
      slice[q] = 0B
      CImap[i,0:*,0:*] = slice
    endfor
  endif
endif

if (arg_present(phasemap)) then begin
  bigCImap = CImap
  phasemap = CImap
endif

if (DPmergedata.palette eq 1) then begin
  slice = reform(cstrip[0,0:*,0:*])
  cstrip[2,0:*,0:*] = slice
  if arg_present(phasemap) then begin
    slice = reform(bigCImap[0,0:*,0:*])
    bigCImap[2,0:*,0:*] = slice
    phasemap = bigCImap
  endif
endif

; and clean up the pixmap window
wdelete,10

end

