@Core_WidgetEvent
@Core_WidgetChoiceEvent
@Core_WText
@Core_WTextE
@Core_PhaseLine
@Core_Print
@Core_getenv
@Core_quatmult
@Core_quat_Lp
@Core_eu2qu
@Core_qu2eu
@Core_LambertSphereToSquare
@Core_histnd
@Core_mind
@DPmergeevent
@DPmerge_event
@DPmergegetpreferences
@DPmergewritepreferences
@DPmergegetfilename
@DPmerge_readDIfile
@DPmerge_display
@DPmerge_display_event
@DPmerge_CIdisplay
@DPmerge_CIdisplay_event
@DPmerge_binary
@DPmerge_cstrip
@DPmerge_compute_pcv
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
; EMsoft:DPmerge.pro
;--------------------------------------------------------------------------
;
; PROGRAM: DPmerge.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Interactive program to merge two or more DI dot product files
;
;> @details 04/06/23 a new GUI to merge two or more DI dot product files and 
;> generate phase maps with orientation and phase confidence metrics
;> Code is bsed on the following paper:
;> F. Ram and M. De Graef. "Phase Differentiation by Electron Backscatter Diffraction 
;> using the Dictionary Indexing Approach". Acta Materialia. 144, 352-364 (2018).
;> 
;>
;> @date 04/06/23 MDG 1.0 first attempt at a user-friendly interface
;--------------------------------------------------------------------------
pro DPmerge,dummy

common DPmerge_widget_common, DPmergewidget_s
common DPmerge_data_common, DPmergedata
common CommonCore, status, logmode, logunit
common DIdata_common, DIdata, w, h

; before we do anything, we make sure that the location of the app_user_dir is set 
appdir = app_user_dir('EMsoft','EMsoftPackage','VMapps','Virtual Machine Apps',['This folder is used by vitual machine apps within EMsoft'],1)

!EXCEPT = 0
logmode = 0
logunit = 10

; widget structure
DPmergewidget_s = {widgetstruct, $
            base:long(0), $                     
            mainstop:long(0), $                     
        	controlbase:long(0), $                 
            displaybase:long(0), $               
            CIdisplaybase:long(0), $               
            displayoption:long(0), $
            cancelbutton:long(0), $
            cancelwidget:long(0), $
            logodraw:long(0), $
            logodrawid:long(0), $
            imageformat:long(0), $
            draw:long(0), $
            drawid:long(0), $
            CIdraw:long(0), $
            CIdrawid:long(0), $
            status:long(0), $
            reset:long(0), $
            Nphases:long(0), $
            dpfile:lonarr(5), $
            palette:long(0), $
            wf1:long(0), $
            wf2:long(0), $
            wf3:long(0), $
            wf4:long(0), $
            wf5:long(0), $
            min:long(0), $
            max:long(0), $
            Mval:long(0), $
            DIloadfile:long(0), $
            ADP:lonarr(5), $
            IQ:lonarr(5), $
            CI:lonarr(5), $
            OSM:lonarr(5), $
            test:long(0) }

; data structure
DPmergedata = {DPmergedatastruct, $
                Nphases:fix(2), $
                Mval:fix(1), $
                maxM:fix(0), $
                wf:replicate(1.0,5), $
                dpfiles:strarr(5), $
                ipf_wd:fix(0), $
                ipf_ht:fix(0), $
                CI_wd:fix(660), $
                CI_ht:fix(615), $
                scrdimx:fix(0), $
                scrdimy:fix(0), $
                xlocation:fix(0), $
                ylocation:fix(0), $
                xlocationcontrol:fix(0), $
                ylocationcontrol:fix(0), $
                xlocationdisplay:fix(0), $
                ylocationdisplay:fix(0), $
                xlocationCIdisplay:fix(0), $
                ylocationCIdisplay:fix(0), $
                navstepsize:float(0), $
                navselector:long(0), $
                imageformat:long(0), $
                displayoption:long(0), $
                drawID:long(0), $
                CIdrawID:long(0), $
                palette:fix(0), $
                DImergeroot:replicate('',5), $
                DProot:'undefined', $
                filesize:intarr(5), $
                DIpathname:replicate('',5), $
                DIfilename:replicate('',5), $
                DIsuffix:replicate('',5), $
                pathname:'', $
                suffix:'', $
                compute:long(0), $
                homefolder:'', $
                nprefs:long(0), $
                appdir: appdir, $               ; location of the user application folder
                prefname: 'DPmergegui.prefs', $    ; filename of preferences file (will be located inside data.appdir)
                foldersep: '/', $               ; folder separator character ('/' for OS X and Linux, '\' for Windows)
                test:long(0) }

; set the foldersep string
if ( (!version.os ne 'darwin') and (!version.os ne 'linux') ) then DPmergedata.foldersep = '\'
DPmergedata.appdir = DPmergedata.appdir+DPmergedata.foldersep

;------------------------------------------------------------
; get the display window size to 80% of the current screen size (but be careful with double screens ... )
; We'll need to guess whether or not the user has a double screen: if the aspect ratio is larger than 16/9,
; then there are likely two screens, so we need to limit ourselves to just the first one...
; This should really become a core function that we can call from all programs.
device,decomposed = 0
device, GET_SCREEN_SIZE = scr

sar = float(scr[0])/float(scr[1])
if (sar gt (1.1*16.0/9.0)) then begin
	scr[0] = scr[0]/2
end
DPmergedata.scrdimy = scr[1] * 0.8
DPmergedata.scrdimx = scr[0]
DPmergedata.xlocation = DPmergedata.scrdimx / 8.0
DPmergedata.ylocation = DPmergedata.scrdimx / 8.0

;------------------------------------------------------------
; does the preferences file exist ?  If not, create it, otherwise read it
; this should also fill in some of the default values for the refinable parameters and the stepsizes and such
DPmergegetpreferences,/noprint

;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
; a few font strings (this will need to be redone for Windows systems)
fontstr='-adobe-new century schoolbook-bold-r-normal--14-100-100-100-p-87-iso8859-1'
fontstrlarge='-adobe-new century schoolbook-medium-r-normal--20-140-100-100-p-103-iso8859-1'
fontstrsmall='-adobe-new century schoolbook-medium-r-normal--14-100-100-100-p-82-iso8859-1'

;------------------------------------------------------------
; create the top level widget
DPmergewidget_s.base = WIDGET_BASE(TITLE='Dictionary Indexing dpmerge Program', $
                        /COLUMN, $
                        XSIZE=1220, $
                        /ALIGN_LEFT, $
                        /TLB_MOVE_EVENTS, $
			            EVENT_PRO='DPmerge_event', $
                        XOFFSET=DPmergedata.xlocation, $
                        YOFFSET=DPmergedata.ylocation)

block0 = WIDGET_BASE(DPmergewidget_s.base, $
			XSIZE=1220, $
			/ALIGN_CENTER, $
			/ROW)

;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------
; create the two main columns for the top row
; block 1 is the left column, with the logo 
block1 = WIDGET_BASE(block0, $
			/FRAME, $
			XSIZE=610, $
			/ALIGN_CENTER, $
			/ROW)

DPmergewidget_s.logodraw = WIDGET_DRAW(block1, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/FRAME, $
			/ALIGN_CENTER, $
			XSIZE=600, $
			YSIZE=200)

;------------------------------------------------------------
; block 2 is the right column, with the input file widgets
; we're asking for the master pattern file name, which will load everything needed.
; then we display here the energyfile name, and perhaps a few other items
block2 = WIDGET_BASE(block0, $
			XSIZE=610, $
			/FRAME, $
			/COLUMN)

;----------  Number of dot product files to consider <= 5
item2 = WIDGET_BASE(block2, /ROW, XSIZE=650, /ALIGN_LEFT)
DPmergewidget_s.Nphases= Core_WTextE(item2,'Number of phases to consider', fontstr, 250, 25, 10, 1, string(DPmergedata.Nphases,format="(I4)"),'NPHASES','DPmerge_event')

DPmergewidget_s.DIloadfile = WIDGET_BUTTON(item2, $
                                UVALUE='DIFILE', $
                                VALUE='Load DI files', $
                                EVENT_PRO='DPmerge_event', $
                                SENSITIVE=1, $
                                /FRAME)

DPmergewidget_s.reset = WIDGET_BUTTON(item2, $
                                UVALUE='RESET', $
                                VALUE='Reset', $
                                EVENT_PRO='DPmerge_event', $
                                SENSITIVE=1, $
                                /FRAME)

DPmergewidget_s.mainstop = WIDGET_BUTTON(item2, $
                                UVALUE='QUIT', $
                                VALUE='Quit', $
                                EVENT_PRO='DPmerge_event', $
                                SENSITIVE=1, $
                                /FRAME)

;---------- experimental pattern file name
file1 = WIDGET_BASE(block2, /ROW, XSIZE=600, /ALIGN_CENTER)
DPmergewidget_s.dpfile[0] = Core_WText(file1,'DI file 1', fontstr, 75, 25, 80, 1, DPmergedata.dpfiles[0])
file1 = WIDGET_BASE(block2, /ROW, XSIZE=600, /ALIGN_CENTER)
DPmergewidget_s.dpfile[1] = Core_WText(file1,'DI file 2', fontstr, 75, 25, 80, 1, DPmergedata.dpfiles[1])
file1 = WIDGET_BASE(block2, /ROW, XSIZE=600, /ALIGN_CENTER)
DPmergewidget_s.dpfile[2] = Core_WText(file1,'DI file 3', fontstr, 75, 25, 80, 1, DPmergedata.dpfiles[2])
file1 = WIDGET_BASE(block2, /ROW, XSIZE=600, /ALIGN_CENTER)
DPmergewidget_s.dpfile[3] = Core_WText(file1,'DI file 4', fontstr, 75, 25, 80, 1, DPmergedata.dpfiles[3])
file1 = WIDGET_BASE(block2, /ROW, XSIZE=600, /ALIGN_CENTER)
DPmergewidget_s.dpfile[4] = Core_WText(file1,'DI file 5', fontstr, 75, 25, 80, 1, DPmergedata.dpfiles[4])


;------------------------------------------------------------
;------------------------------------------------------------
;------------------------------------------------------------

block1 = WIDGET_BASE(DPmergewidget_s.base, $
			XSIZE=1200, $
			/ALIGN_CENTER, $
			/ROW)

block2 = WIDGET_BASE(block1, $
			/FRAME, $
			XSIZE=800, $
			/ALIGN_CENTER, $
			/COLUMN)

;------------------------------------------------------------
tmp1 = WIDGET_LABEL(block2, VALUE='Adjustable Dot Product Weight Factors', font=fontstrlarge, /ALIGN_LEFT)


line1 = WIDGET_BASE(block2, XSIZE=1200, /ROW)

;----------  weight factors for each phase 
item2 = WIDGET_BASE(line1, /ROW, XSIZE=1000, /ALIGN_LEFT)
DPmergewidget_s.wf1= Core_WTextE(item2,'Phase 1', fontstr, 70, 25, 10, 1, string(DPmergedata.wf[0],format="(F9.2)"),'WF1','DPmerge_event')
DPmergewidget_s.wf2= Core_WTextE(item2,'Phase 2', fontstr, 70, 25, 10, 1, string(DPmergedata.wf[1],format="(F9.2)"),'WF2','DPmerge_event')
DPmergewidget_s.wf3= Core_WTextE(item2,'Phase 3', fontstr, 70, 25, 10, 1, string(DPmergedata.wf[2],format="(F9.2)"),'WF3','DPmerge_event')
DPmergewidget_s.wf4= Core_WTextE(item2,'Phase 4', fontstr, 70, 25, 10, 1, string(DPmergedata.wf[3],format="(F9.2)"),'WF4','DPmerge_event')
DPmergewidget_s.wf5= Core_WTextE(item2,'Phase 5', fontstr, 70, 25, 10, 1, string(DPmergedata.wf[4],format="(F9.2)"),'WF5','DPmerge_event')

line2 = WIDGET_BASE(block2, XSIZE=700, /ROW)

;----------  M value 
item1 = WIDGET_BASE(line2, /ROW, XSIZE=350, /ALIGN_LEFT)
DPmergewidget_s.Mval = Core_WTextE(item1,'M value', fontstr, 75, 25, 10, 1, string(DPmergedata.Mval,format="(I4)"),'MVAL','DPmerge_event')

item2 = WIDGET_BASE(line2, /ROW, XSIZE=350, /ALIGN_LEFT)
vals = ['Standard','Deuteranopia']
DPmergewidget_s.palette= CW_BGROUP(item2, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT='Color Palette', $
                        EVENT_FUNC ='DPmergeevent', $
                        UVALUE='PALETTE', $
                        SET_VALUE=DPmergedata.palette)

;----------  display controls
block2 = WIDGET_BASE(block1, $
			/FRAME, $
			XSIZE=450, $
			/ALIGN_LEFT, $
			/COLUMN)

tmp1 = WIDGET_LABEL(block2, VALUE='Display Controls', font=fontstrlarge, /ALIGN_LEFT)

; and create each row of fitting widgets
ret = Core_PhaseLine(block2, 'ADP')
ret = Core_PhaseLine(block2, ' IQ')
ret = Core_PhaseLine(block2, ' CI')
ret = Core_PhaseLine(block2, 'OSM')

block1 = WIDGET_BASE(DPmergewidget_s.base, $
            XSIZE=1200, $
            /ALIGN_CENTER, $
            /ROW)

; then we have the program message window
DPmergewidget_s.status= WIDGET_TEXT(block1, $
			XSIZE=195, $
			YSIZE=10, $
			/SCROLL, $
			VALUE=' ',$
			/ALIGN_LEFT)

; the following is needed by the Core_Print routine
status = DPmergewidget_s.status 

;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,DPmergewidget_s.base,/REALIZE

; realize the draw widgets
WIDGET_CONTROL, DPmergewidget_s.logodraw, GET_VALUE=drawID
DPmergewidget_s.logodrawID = drawID
;
read_jpeg,'Resources/EMsoftVBFFlogo.jpeg',logo
wset,DPmergewidget_s.logodrawID
tvscl,logo,true=1

; and hand over control to the xmanager
XMANAGER,"DPmerge",DPmergewidget_s.base,/NO_BLOCK

end ; program
