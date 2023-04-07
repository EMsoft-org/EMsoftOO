
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
; EMsoft:Core_PhaseLine.pro
;--------------------------------------------------------------------------
;
; PROGRAM: Core_PhaseLine.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief generate a widget line for DPmerge.pro
;
;> @date 04/06/23 MDG 1.0 first attempt at a user-friendly interface
;--------------------------------------------------------------------------
function Core_PhaseLine, block2, LBL

common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata

line1 = WIDGET_BASE(block2, XSIZE=1200, /ROW)
item1 = WIDGET_LABEL(line1, VALUE=LBL, font=fontstr, /ALIGN_LEFT)

; then insert the 5 buttons
for i=0,4 do begin 
    uval = strtrim(LBL,1) + string(i+1,format="(I1)")
    if (LBL eq 'ADP') then begin
        DPmergewidget_s.ADP[i]= WIDGET_BUTTON(line1, $
                                UVALUE=uval, $
                                VALUE=string(i+1,format="(I1)"), $
                                FONT=fontstrlarge, $
                                EVENT_PRO='DPmerge_event', $
                                SENSITIVE=0)
    endif
    if (LBL eq ' CI') then begin
        DPmergewidget_s.CI[i]= WIDGET_BUTTON(line1, $
                                UVALUE=uval, $
                                VALUE=string(i+1,format="(I1)"), $
                                FONT=fontstrlarge, $
                                EVENT_PRO='DPmerge_event', $
                                SENSITIVE=0)
    endif
    if (LBL eq ' IQ') then begin
        DPmergewidget_s.IQ[i]= WIDGET_BUTTON(line1, $
                                UVALUE=uval, $
                                VALUE=string(i+1,format="(I1)"), $
                                FONT=fontstrlarge, $
                                EVENT_PRO='DPmerge_event', $
                                SENSITIVE=0)
    endif
    if (LBL eq 'OSM') then begin
        DPmergewidget_s.OSM[i]= WIDGET_BUTTON(line1, $
                                UVALUE=uval, $
                                VALUE=string(i+1,format="(I1)"), $
                                FONT=fontstrlarge, $
                                EVENT_PRO='DPmerge_event', $
                                SENSITIVE=0)
    endif

end

return,0
end
