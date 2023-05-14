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
; EMsoft:DPmergewritepreferences.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmergewritepreferences.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief write the DPmerge preferences file
;
;> @date 04/06/23 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro DPmergewritepreferences,noprint=noprint
 
;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata

; prefs file
  openw,1,DPmergedata.appdir+DPmergedata.prefname
  nprefs = 9
  DPmergedata.nprefs = nprefs
  printf,1,nprefs
  printf,1,'DProot::'+DPmergedata.DImergeroot[0]

; fixed parameters


; window locations
  printf,1,'xlocation::'+string(DPmergedata.xlocation,format="(F6.1)")
  printf,1,'ylocation::'+string(DPmergedata.ylocation,format="(F6.1)")
  printf,1,'xlocationcontrol::'+string(DPmergedata.xlocationcontrol,format="(F6.1)")
  printf,1,'ylocationcontrol::'+string(DPmergedata.ylocationcontrol,format="(F6.1)")
  printf,1,'xlocationdisplay::'+string(DPmergedata.xlocationdisplay,format="(F6.1)")
  printf,1,'ylocationdisplay::'+string(DPmergedata.ylocationdisplay,format="(F6.1)")
  printf,1,'xlocationCIdisplay::'+string(DPmergedata.xlocationCIdisplay,format="(F6.1)")
  printf,1,'ylocationCIdisplay::'+string(DPmergedata.ylocationCIdisplay,format="(F6.1)")
; and close the file
  close,1

  if not keyword_set(noprint) then Core_Print,'The preferences file '+DPmergedata.appdir+DPmergedata.prefname+' was successfully saved '

end

