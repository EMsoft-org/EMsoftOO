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
; EMsoft:DPmerge_ternary.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_ternary.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief phase confidence map computation for the binary case 
;
;> @date 04/12/23 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro DPmerge_ternary, numnearmatches, cmap, phasemap=phasemap, dpmap=dpmap, clevlines=clevlines, histval=histval

;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata

common DIdata_common, DIdata, w, h 
common triangleparms, E, Z, xoff, yoff, alpha

alpha = 0.0

if (arg_present(clevlines)) then alpha = (1.0-clevlines/40.0+1.5) * 0.2

wset,DPmergedata.CIdrawID
erase, color=255

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
Ctdp = DIdata[2].TDP

if arg_present(phasemap) then begin
  phasemap = bytarr(3,numx,numy)
endif

if arg_present(dpmap) then begin 
  dpmap = fltarr(nn)
; and fill the dotproductmap with the highest values for the two phases
  for i=0,nn-1 do begin
    m = [Atdp[0,i],Btdp[0,i],Ctdp[0,i]]
    dpmap[i] = max(m)
  endfor
  dpmap = congrid(reverse(reform(dpmap, numx, numy),2),2*numx,2*numy)
endif

; truncate the array and scale the dot products by the weight factors
Atdp = Atdp[0:numm-1,0:nn-1] * DPmergedata.wf[0]
Btdp = Btdp[0:numm-1,0:nn-1] * DPmergedata.wf[1]
Ctdp = Ctdp[0:numm-1,0:nn-1] * DPmergedata.wf[2]

MA = fltarr(nn)
MB = fltarr(nn)
MC = fltarr(nn)
CImap = bytarr(3,numx, numy)

; figure out the intensity scaling via color saturation
mindp = min([Atdp, Btdp, Ctdp], max=maxdp)

Atdp = (Atdp-mindp)/(maxdp-mindp)
Btdp = (Btdp-mindp)/(maxdp-mindp)
Ctdp = (Ctdp-mindp)/(maxdp-mindp)

avdp = reform((Atdp[0,*] + Btdp[0,*] + Ctdp[0,*]) / 3.0)

line = [replicate(1,numm),replicate(2,numm),replicate(3,numm)]

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
  q = where(newline eq 3,cnt)
  if (cnt gt 0.0) then MC[i] = total(dp[q])/denominator else MC[i] = 0.0
endfor

MA = reform(MA,numx,numy)
MB = reform(MB,numx,numy)
MC = reform(MC,numx,numy)
avdp = reform(avdp,numx,numy)

if arg_present(phasemap) then begin
  CImap = bytarr(3,numx,numy)
  for i=0,numx-1 do begin
    for j=0,numy-1 do begin
      m = max([MA[i,j],MB[i,j],MC[i,j]])
      CImap[0:2,i,j] = byte(255.0*[MA[i,j],MB[i,j],MC[i,j]] / m )
    endfor
  endfor
endif

; next generate the points in a triangular color map
E = 600
z = sqrt(3.0)*0.5
xoff = 30
yoff = 30

device,decomposed=0
window,10,xsi=660,ysi=571,retain=2,/pixmap
erase

; generate the color triangle
bg = DPmerge_ctriangle()
bgdim = size(bg,/dimensions)

; transform the data points to the triangular coordinates 
y = reform(MC,nn) * E * z
x = -E*(reform(MA,nn)+reform(MC,nn)*0.5-1.0)

; generate a histogram for each of the potential points inside the triangle
nummm = 30
np = (nummm+2)*(nummm+1)/2
; first define the coordinates of the grid points
gridx = fltarr(np)
gridy = fltarr(np)
k = 0
dx = E/float(nummm)
dy = dx*z
for j=0,nummm do begin
  for i=0,nummm-j do begin
    gridx[k] = j*dx*0.5+i*dx
    gridy[k] = j*dy
    k += 1
  endfor
endfor 

avdp = reform(avdp,nn)
avdp -= min(avdp)
avdp /= max(avdp)

; then fill the histogram bins based on Euclidean distance to the grid centers
histn = fltarr(np)
dthr = dx*0.5
for i=0,np-1 do begin
  d = sqrt( (x - gridx[i])^2 + (y - gridy[i])^2 )
  for j=0,nn-1 do if (d[j] lt dthr) then histn[i] += 1.0 
endfor

if (arg_present(histvalues)) then histvalues = histn

; rescale the histogram to a logarithmic scale
h = fix(alog10(histn+1.0))*2.5

; draw filled circles at all histogram bins that have more than 10 counts in them
nt = 360.
th = findgen(nt)*!dtor
ct = cos(th)
st = sin(th)
for k=0,np-1 do if (h[k] gt 0.001) then polyfill,xoff+gridx[k]+h[k]*ct,yoff+gridy[k]+h[k]*st,/dev,color=0

; draw the circle set for the legend
gx = 500
gy = reverse(500 - findgen(max(h)/2.5+1)*30)
for k=fix(max(h)/2.5),1,-1 do polyfill,xoff+gx+k*2.5*ct,yoff+gy[k]+k*2.5*st,/dev,color=0
empty

; finally, read the color image into an array that is returned to the calling program
cmap = tvrd(0,0,bgdim[1],bgdim[2],TRUE=1)

; delete the drawing window since it is no longer needed.
wdelete,10

; next, go through the entire set of point coordinates and for those that fall
; inside the region below the confidence level, turn those pixels to white in
; the CImap...

;islow = check_low_confidence(reform(MA,nn), reform(MC,nn), nn)
islow = DPmerge_check_low_confidence(x, y, nn)
q = where(islow eq 1, cnt)

if (cnt gt 0) then begin
  for i=0,2 do begin
    slice = reform(CImap[i,*,*])
    slice[q] = 255B
    CImap[i,0:*,0:*] = slice
  endfor
endif

; scale up the phasemap
if (arg_present(phasemap)) then begin
  bigCImap = bytarr(3,2*numx,2*numy)
  for i=0,2 do begin
    slice = reform(CImap[i,*,*])
    slice = reverse(congrid(slice,2*numx,2*numy),2)
    bigCImap[i,0:*,0:*] = slice
  endfor
  phasemap = bigCImap
end

; do we need to convert the output to a color scheme more suitable for people with deuteranopia ?
if (DPmergedata.palette eq 1) then begin  
  slice = reform(cmap[0,0:*,0:*])
  cmap[2,0:*,0:*] = slice
  if (arg_present(phasemap)) then begin
    slice = reform(bigCImap[0,0:*,0:*])
    bigCImap[2,0:*,0:*] = slice
    phasemap = bigCImap
  endif
endif


end

