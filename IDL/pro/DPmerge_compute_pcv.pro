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
; EMsoft:DPmerge_compute_pcv.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_compute_pcv.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Perform the phase confidence vector computation + display results for 2 or 3 phases
;
;> @date 04/10/23 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro DPmerge_compute_pcv,dummy
 
;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata

common DIdata_common, DIdata, w, h 

if (DPmergedata.Nphases eq 2) then begin 
    np = 15
    cmap = bytarr(3,660, 41 * np)
    for i=2,2+np-1 do begin
      DPmerge_binary,2*i-2,cstrip, clev=1
      cstrip = cstrip[0:2,0:*,60-20:60+20]
      cmap[0:2,0:*,(i-2)*41:(i-1)*41-1] = cstrip
    endfor
    wset,DPmergedata.CIdrawID
    tvscl,cmap,true=1

    for i=0,np-1 do xyouts,3,18 + i*41,string(2*i+2,format="(I2.2)"),/dev, color=0

    Core_Print,'Binary phase confidence index maps'
endif 

if (DPmergedata.Nphases eq 3) then begin 
  clev = 1
  DPmerge_ternary,DPmergedata.Mval,cmap,phasemap=pmap, clevlines=clev, dpmap=dpmap
  wset,DPmergedata.CIdrawID
  tv,cmap,true=1
  Core_Print,'Ternary phase confidence index maps'
endif


end