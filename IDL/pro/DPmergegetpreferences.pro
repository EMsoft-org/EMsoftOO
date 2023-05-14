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
; EMsoft:DPmergegetpreferences.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmergegetpreferences.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief read the preferences file and initialize all relevant widgets
;
;> @date 04/06/33 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro DPmergegetpreferences,noprint=noprint
 
;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata

; does the preferences file exist ?
rs = file_test(DPmergedata.appdir+DPmergedata.prefname)

if (rs eq 1) then begin
  s = ''
  i = 0
  openr,1,DPmergedata.appdir+DPmergedata.prefname
  readf,1,i
  DPmergedata.nprefs = i

; next, do a little loop and read name:value pairs
  for i=0,DPmergedata.nprefs-1 do begin
    readf,1,s
    spos = strpos(s,'::')
    nm = strmid(s,0,spos)
    val = strmid(s,spos+2)
    case nm of 
; root folder
	 'DProot': begin 
              DPmergedata.DImergeroot[0]=val
              DPmergedata.DProot = val
             end   

; various parameters

; window locations
  	'xlocation': DPmergedata.xlocation = float(val)
  	'ylocation': DPmergedata.ylocation = float(val)
  	'xlocationcontrol': DPmergedata.xlocationcontrol = float(val)
  	'ylocationcontrol': DPmergedata.ylocationcontrol = float(val)
    'xlocationdisplay': DPmergedata.xlocationdisplay = float(val)
    'ylocationdisplay': DPmergedata.ylocationdisplay = float(val)
    'xlocationCIdisplay': DPmergedata.xlocationCIdisplay = float(val)
    'ylocationCIdisplay': DPmergedata.ylocationCIdisplay = float(val)

    else: MESSAGE,'unknown option for preferences file'
    endcase
  endfor

  close,1
end else begin
  s = ''
  cd,current=s
  DPmergedata.DProot=s
; prefs file does not exist yet, so let's create it with default values
  if not keyword_set(noprint) then Core_Print,'Creating preferences file '+DPmergedata.appdir+DPmergedata.prefname
  if keyword_set(noprint) then DPmergewritepreferences,/noprint else DPmergewritepreferences
endelse

end

