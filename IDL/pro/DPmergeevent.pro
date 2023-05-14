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
; EMsoft:DPmergeevent.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmergeevent.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief special event handler for all the CW_BGROUP calls, since CW_BGROUP does not support event_pro
;
;> @date 04/07/23 MDG 1.0 first version
;--------------------------------------------------------------------------
function DPmergeevent, event

;------------------------------------------------------------
; common blocks
common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

CASE eventval OF

        ; 'DISPLAYOPTION' : begin   ; this comes from the DPmerge_control widget...
        ;         DPmergedata.displayoption = Core_WidgetChoiceEvent( DPmergewidget_s.displayoption, 'Set option to ',/value)
        ;         DPmerge_showpattern
        ; endcase

        'PALETTE' : begin
                DPmergedata.palette = Core_WidgetChoiceEvent( DPmergewidget_s.palette, 'Set palette to ',/value)
                DPmerge_compute_pcv
                clev = 1.0
                DPmerge_binary,DPmergedata.Mval, cstrip, phasemap=pmap
                wset,DPmergedata.drawID
                tvscl,pmap,true=1
        endcase

        'IMAGEFORMAT' : begin
                DPmergedata.imageformat = Core_WidgetChoiceEvent( DPmergewidget_s.imageformat, 'Set option to ',/value)
        endcase

else: MESSAGE, "DPmergeevent: Event User Value Not Found"

endcase

return,eventval
end 
