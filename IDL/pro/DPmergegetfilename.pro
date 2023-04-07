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
; EMsoft:DPmergegetfilename.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmergegetfilename.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Display an interface and ask user to select an input DI file
;
;> @date 04/07/23 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro DPmergegetfilename,validfile
 
;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata

common DIdata_common, DIdata, w, h 

  if (total(size(DIdata)) ne 0) then DELVAR, DIdata

  validfile = 0

  s = ''
  cd,current = s
  DPmergedata.homefolder = s
  if (DPmergedata.DProot eq 'undefined') then begin
    DPmergedata.DProot = DPmergedata.homefolder
  end 

  rootpath = DPmergedata.DProot

  for ii=0,DPmergedata.Nphases-1 do begin
    res=dialog_pickfile(title='Select DI file '+string(ii+1,format="(I1)"),path=rootpath,filter='*.h5')
    if (res eq '') then begin
	    Core_Print,'No selection made'
	    goto, skip
    end
	  validfile = 1
  	finfo = file_info(res)
    DPmergedata.filesize[ii] = finfo.size
; find the last folder separator
    spos = strpos(res,'/',/reverse_search)
  	dpos = strpos(res,'.',/reverse_search)
  	plen = strlen(res)
  	DPmergedata.DIpathname[ii] = strmid(res,0,spos)
  	DPmergedata.DIfilename[ii] = strmid(res,spos+1)
  	DPmergedata.DIsuffix[ii] = strmid(res,dpos+1)
    DPmergedata.DImergeroot[ii] = DPmergedata.DIpathname[ii]

  	WIDGET_CONTROL, SET_VALUE=DPmergedata.DIfilename[ii], DPmergewidget_s.dpfile[ii]

  	Core_Print,' full path '+res
  	Core_Print,' path '+DPmergedata.DIpathname[ii]
  	Core_Print,' data file '+DPmergedata.DIfilename[ii]
  	Core_Print,' suffix '+DPmergedata.DIsuffix[ii]
    Core_Print,' ---> Loading data'

    DPmerge_readDIfile, res, ii, status

    if (ii gt 0) then begin   ; we need to make sure that the dimensions of ROIs in all files are identical
                              ; if they are not identical, we abort the data loading here and displace an error message
      if (status eq -5) then begin
        Core_Print,'The ROI dimensions of this data file are different from those of the previous file !!!'
        Core_Print,'  Aborting remaining file reading'
        DELVAR, DIdata 
        goto, skip
      endif
    endif
    Core_Print,'      Done' ,/blank

  endfor

; if we get here, then we need to activate a bunch of widgets 
  Core_Print,' ROI dimensions '+string(DPmergedata.ipf_wd,format="(I4)")+' by '+string(DPmergedata.ipf_ht,format="(I4)")

  for ii=0,DPmergedata.Nphases-1 do begin
    WIDGET_CONTROL, DPmergewidget_s.ADP[ii], sensitive = 1
    WIDGET_CONTROL, DPmergewidget_s.CI[ii], sensitive = 1
    WIDGET_CONTROL, DPmergewidget_s.IQ[ii], sensitive = 1
    WIDGET_CONTROL, DPmergewidget_s.OSM[ii], sensitive = 1
  endfor

; at this point, we should compute the standard phase map, using the highest dot product for each phase
; and compute the phase confidence vector f.  Once that is done, then we can start playing around with 
; the value of M to be used for the computation of f.
DPmerge_compute_pcv

skip:
end

