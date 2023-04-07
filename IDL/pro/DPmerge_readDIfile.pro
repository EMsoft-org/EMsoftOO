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
; EMsoft:DPmerge_readDIfile.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_readDIfile.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Read necessary information from a Dictionary Indexing output file 
;
;> @date 04/07/23 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro DPmerge_readDIfile, h5name, sel, status

common DIdata_common, DIdata, w, h
common DPmerge_data_common, DPmergedata

status = 0

; read all the necessary data from the DI h5 file
file_id = H5F_OPEN(h5name)
; get the ROI dimensions
group_id = H5G_OPEN(file_id,'NMLparameters/EBSDIndexingNameListType')
; group_id = H5G_OPEN(file_id,'NMLparameters/')
dset_id = H5D_OPEN(group_id,'ipf_wd')
w = H5D_READ(dset_id)
w = w[0]
H5D_close,dset_id
dset_id = H5D_OPEN(group_id,'ipf_ht')
h = H5D_READ(dset_id)
h = h[0]
H5D_close,dset_id
dset_id = H5D_OPEN(group_id,'nnk')
DPmergedata.maxM = H5D_READ(dset_id)
H5D_close,dset_id
H5G_close,group_id

; get the dimensions of the TopDotProduct array first so we can declare the DIdata structure
group_id = H5G_OPEN(file_id,'Scan 1/EBSD/Data')
dset_id = H5D_OPEN(group_id,'TopDotProductList')
tdp = H5D_READ(dset_id)
H5D_close,dset_id
sz = size(tdp, /dimensions)
DPmergedata.maxM = sz[0]

entry = { DIarrays, ADP:bytarr(w,h), CI:bytarr(w,h), IQ:bytarr(w,h), OSM:bytarr(w,h), TDP:fltarr(sz[0],sz[1]), ROIdims:lonarr(2) }
entry.ROIdims = [w,h]
if (sel eq 0) then begin 
    DIdata = replicate(entry, DPmergedata.Nphases)
    DPmergedata.ipf_wd = w 
    DPmergedata.ipf_ht = h 
end else begin
    if ( ( w ne DPmergedata.ipf_wd ) or ( h ne DPmergedata.ipf_ht) ) then begin
        status = -5
        goto, skip
    endif 
endelse

DIdata[sel] = entry 
DIdata[sel].TDP = tdp

dset_id = H5D_OPEN(group_id,'AvDotProductMap')
DIdata[sel].ADP = H5D_READ(dset_id)
H5D_close,dset_id
dset_id = H5D_OPEN(group_id,'IQMap')
DIdata[sel].IQ = H5D_READ(dset_id)
H5D_close,dset_id
dset_id = H5D_OPEN(group_id,'OSM')
DIdata[sel].OSM = bytscl(reform(H5D_READ(dset_id),w,h))
H5D_close,dset_id
dset_id = H5D_OPEN(group_id,'CIMap')
DIdata[sel].CI = H5D_READ(dset_id)
H5D_close,dset_id

skip:
H5G_close,group_id
H5F_close,file_id

end