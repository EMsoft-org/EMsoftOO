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
; EMsoft:DPmerge_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main event handler for DPmerge.pro program
;
;> @date 04/07/23 MDG 1.0 first attempt at a user-friendly interface
;--------------------------------------------------------------------------
pro DPmerge_event,event

common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata
common CommonCore, status, logmode, logunit
common DIdata_common, DIdata, w, h


if (event.id eq DPmergewidget_s.base) then begin
  DPmergedata.xlocation = event.x
  DPmergedata.ylocation = event.y-25
end else begin

  WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value
  
  CASE eventval OF
        'DIFILE': begin
; ask the user to select the data file
		DPmergegetfilename,validfile

; then start the actual display widget and put the experimental pattern in it...
; better yet, we can check what the dimensions of the current drawing area are in 
; case the new experimental pattern has different dimensions...
                if (XRegistered("DPmerge_displaybase") NE 0) then begin
		  WIDGET_CONTROL, DPmergewidget_s.displaybase, /DESTROY
                endif 
                DPmerge_display

        endcase

        'NPHASES': begin
                DPmergedata.hipasscutoff = Core_WidgetEvent( DPmergewidget_s.hipasscutoff,  'hipass cut off set to ', '(F9.2)', /flt)
                WIDGET_CONTROL, DPmergewidget_s.mkjson, sensitive=1
        endcase

 
 	'QUIT': begin
		DPmergewritepreferences
; do a general cleanup of potentially open widgets
                if (XRegistered("DPmerge_control") NE 0) then begin
		  WIDGET_CONTROL, DPmergewidget_s.controlbase, /DESTROY
                endif
                if (XRegistered("DPmerge_display") NE 0) then begin
		  WIDGET_CONTROL, DPmergewidget_s.displaybase, /DESTROY
                endif
 		Core_Print,'Quitting program',/blank
                wait,1.0
		WIDGET_CONTROL, DPmergewidget_s.base, /DESTROY
		!EXCEPT=1
	endcase

; here are the display map options 
        'ADP1': begin 
                  erase 
                  tv,DIdata[0].ADP
        endcase

        'ADP2': begin 
                  erase 
                  tv,DIdata[1].ADP
        endcase

        'ADP3': begin 
                  erase 
                  tv,DIdata[2].ADP
        endcase

        'ADP4': begin 
                  erase 
                  tv,DIdata[3].ADP
        endcase

        'ADP5': begin 
                  erase 
                  tv,DIdata[4].ADP
        endcase

        'OSM1': begin 
                  erase 
                  tv,DIdata[0].OSM
        endcase

        'OSM2': begin 
                  erase 
                  tv,DIdata[1].OSM
        endcase

        'OSM3': begin 
                  erase 
                  tv,DIdata[2].OSM
        endcase

        'OSM4': begin 
                  erase 
                  tv,DIdata[3].OSM
        endcase

        'OSM5': begin 
                  erase 
                  tv,DIdata[4].OSM
        endcase

        'IQ1': begin 
                  erase 
                  tv,DIdata[0].IQ
        endcase

        'IQ2': begin 
                  erase 
                  tv,DIdata[1].IQ
        endcase

        'IQ3': begin 
                  erase 
                  tv,DIdata[2].IQ
        endcase

        'IQ4': begin 
                  erase 
                  tv,DIdata[3].IQ
        endcase

        'IQ5': begin 
                  erase 
                  tv,DIdata[4].IQ
        endcase

        'CI1': begin 
                  erase 
                  tv,DIdata[0].CI
        endcase

        'CI2': begin 
                  erase 
                  tv,DIdata[1].CI
        endcase

        'CI3': begin 
                  erase 
                  tv,DIdata[2].CI
        endcase

        'CI4': begin 
                  erase 
                  tv,DIdata[3].CI
        endcase

        'CI5': begin 
                  erase 
                  tv,DIdata[4].CI
        endcase

  else: MESSAGE, "DPmerge_event: Event User Step "+eventval+" Not Found"

  endcase


endelse

end
