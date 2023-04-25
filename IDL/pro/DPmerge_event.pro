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

; then start the actual display widget 
                if (XRegistered("DPmerge_displaybase") NE 0) then begin
		              WIDGET_CONTROL, DPmergewidget_s.displaybase, /DESTROY
                endif 
                DPmerge_display
        endcase

        'NPHASES': begin
                DPmergedata.Nphases= Core_WidgetEvent( DPmergewidget_s.Nphases,  ' Number of phases set to ', '(F9.2)', /flt)
                if (DPmergedata.Nphases gt 5) then begin 
                        Core_Print,' Maximum number of phases is 5; resetting value to 5'
                        DPmerge.Nphases = 5 
                endif 
        endcase

        'MVAL': begin
                DPmergedata.Mval = Core_WidgetEvent( DPmergewidget_s.Mval,  ' Mval set to ', '(F9.2)')
                if (DPmergedata.Nphases eq 2) then begin
                        clev = 1.0
                        DPmerge_binary,DPmergedata.Mval, cstrip, phasemap=pmap
                        wset,DPmergedata.drawID
                        tvscl,pmap,true=1
                end else if (DPmergedata.Nphases eq 3) then begin
                        clev = 1
                        DPmerge_ternary,DPmergedata.Mval,cmap,phasemap=pmap, clevlines=clev, dpmap=dpmap
                        wset,DPmergedata.drawID
                        tv,pmap,true=1
                        wset,DPmergedata.CIdrawID
                        tv,cmap,true=1
                        DPmerge_annotate_triangle
                end

        endcase

        'PALETTE': begin
                DPmergedata.palette= Core_WidgetEvent( DPmergewidget_s.palette,  ' Color palette set to ', '(I2)')
        endcase

        'RESET': begin   ; reset all fields to the program starting status
                       if (XRegistered("DPmerge_display") NE 0) then begin
                          WIDGET_CONTROL, DPmergewidget_s.displaybase, /DESTROY
                        endif
                        if (XRegistered("DPmerge_CIdisplay") NE 0) then begin
                          WIDGET_CONTROL, DPmergewidget_s.CIdisplaybase, /DESTROY
                        endif

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
                        if (XRegistered("DPmerge_CIdisplay") NE 0) then begin
                          WIDGET_CONTROL, DPmergewidget_s.CIdisplaybase, /DESTROY
                        endif
             		Core_Print,'Quitting program',/blank
                          wait,1.0
            		WIDGET_CONTROL, DPmergewidget_s.base, /DESTROY
            		!EXCEPT=1
        endcase

; dot product weight factors

        'WF1': begin 
                DPmergedata.wf[0]= Core_WidgetEvent( DPmergewidget_s.wf1,  ' Weight factor 1 set to ', '(F8.5)')
                DPmerge_compute_pcv
                if (DPmergedata.Nphases eq 2) then begin
                        clev = 1.0
                        DPmerge_binary,DPmergedata.Mval, cstrip, phasemap=pmap
                        wset,DPmergedata.drawID
                        tvscl,pmap,true=1
                end else if (DPmergedata.Nphases eq 3) then begin
                        clev = 1
                        DPmerge_ternary,DPmergedata.Mval,cmap,phasemap=pmap, clevlines=clev, dpmap=dpmap
                        wset,DPmergedata.drawID
                        tv,pmap,true=1
                        wset,DPmergedata.CIdrawID
                        tv,cmap,true=1
                        DPmerge_annotate_triangle
                end
        endcase

        'WF2': begin 
                DPmergedata.wf[1]= Core_WidgetEvent( DPmergewidget_s.wf2,  ' Weight factor 2 set to ', '(F8.5)')
                DPmerge_compute_pcv
                if (DPmergedata.Nphases eq 2) then begin
                        clev = 1.0
                        DPmerge_binary,DPmergedata.Mval, cstrip, phasemap=pmap
                        wset,DPmergedata.drawID
                        tvscl,pmap,true=1
                end else if (DPmergedata.Nphases eq 3) then begin
                        clev = 1
                        DPmerge_ternary,DPmergedata.Mval,cmap,phasemap=pmap, clevlines=clev, dpmap=dpmap
                        wset,DPmergedata.drawID
                        tv,pmap,true=1
                        wset,DPmergedata.CIdrawID
                        tv,cmap,true=1
                        DPmerge_annotate_triangle
                end
        endcase

        'WF3': begin 
                DPmergedata.wf[2]= Core_WidgetEvent( DPmergewidget_s.wf3,  ' Weight factor 3 set to ', '(F8.5)')
                DPmerge_compute_pcv
                if (DPmergedata.Nphases eq 2) then begin
                        clev = 1.0
                        DPmerge_binary,DPmergedata.Mval, cstrip, phasemap=pmap
                        wset,DPmergedata.drawID
                        tvscl,pmap,true=1
                end else if (DPmergedata.Nphases eq 3) then begin
                        clev = 1
                        DPmerge_ternary,DPmergedata.Mval,cmap,phasemap=pmap, clevlines=clev, dpmap=dpmap
                        wset,DPmergedata.drawID
                        tv,pmap,true=1
                        wset,DPmergedata.CIdrawID
                        tv,cmap,true=1
                        DPmerge_annotate_triangle
                end
        endcase

        'WF4': begin 
                DPmergedata.wf[3]= Core_WidgetEvent( DPmergewidget_s.wf4,  ' Weight factor 4 set to ', '(F8.5)')
                DPmerge_compute_pcv
                if (DPmergedata.Nphases eq 2) then begin
                        clev = 1.0
                        DPmerge_binary,DPmergedata.Mval, cstrip, phasemap=pmap
                        wset,DPmergedata.drawID
                        tvscl,pmap,true=1
                end else if (DPmergedata.Nphases eq 3) then begin
                        clev = 1
                        DPmerge_ternary,DPmergedata.Mval,cmap,phasemap=pmap, clevlines=clev, dpmap=dpmap
                        wset,DPmergedata.drawID
                        tv,pmap,true=1
                        wset,DPmergedata.CIdrawID
                        tv,cmap,true=1
                        DPmerge_annotate_triangle
                end
        endcase

        'WF5': begin 
                DPmergedata.wf[4]= Core_WidgetEvent( DPmergewidget_s.wf5,  ' Weight factor 5 set to ', '(F8.5)')
                DPmerge_compute_pcv
                if (DPmergedata.Nphases eq 2) then begin
                        clev = 1.0
                        DPmerge_binary,DPmergedata.Mval, cstrip, phasemap=pmap
                        wset,DPmergedata.drawID
                        tvscl,pmap,true=1
                end else if (DPmergedata.Nphases eq 3) then begin
                        clev = 1
                        DPmerge_ternary,DPmergedata.Mval,cmap,phasemap=pmap, clevlines=clev, dpmap=dpmap
                        wset,DPmergedata.drawID
                        tv,pmap,true=1
                        wset,DPmergedata.CIdrawID
                        tv,cmap,true=1
                        DPmerge_annotate_triangle
                end
        endcase

; here are the display map options 
        'ADP1': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[0].ADP
        endcase

        'ADP2': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[1].ADP
        endcase

        'ADP3': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[2].ADP
        endcase

        'ADP4': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[3].ADP
        endcase

        'ADP5': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[4].ADP
        endcase

        'OSM1': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[0].OSM
        endcase

        'OSM2': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[1].OSM
        endcase

        'OSM3': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[2].OSM
        endcase

        'OSM4': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[3].OSM
        endcase

        'OSM5': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[4].OSM
        endcase

        'IQ1': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[0].IQ
        endcase

        'IQ2': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[1].IQ
        endcase

        'IQ3': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[2].IQ
        endcase

        'IQ4': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[3].IQ
        endcase

        'IQ5': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[4].IQ
        endcase

        'CI1': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[0].CI
        endcase

        'CI2': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[1].CI
        endcase

        'CI3': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[2].CI
        endcase

        'CI4': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[3].CI
        endcase

        'CI5': begin 
                  wset,DPmergedata.drawID
                  erase 
                  tv,DIdata[4].CI
        endcase

  else: MESSAGE, "DPmerge_event: Event User Step "+eventval+" Not Found"

  endcase


endelse

end
