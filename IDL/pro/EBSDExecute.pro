;
; Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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
; EMsoft:EBSDExecute.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDExecute.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main routine for creation of nml file and execution of EMEBSD code
;
;> @date 05/22/14 MDG 1.0 first version
;> @date 04/14/15 MDG 1.1 added HDF5 support
;> @date 11/04/15 MDG 2.0 modified routine for call_external computation of EBSD pattern
;> @date 04/14/21 MDG 3.0 updated for EMsoftOO
;--------------------------------------------------------------------------
pro EBSDExecute, status, single=single

; the keyword /single indicates that only one pattern should be computed

;------------------------------------------------------------
; common blocks
common SEM_widget_common, SEMwidget_s
common SEM_data_common, SEMdata
common EBSDpatterns, pattern, image, finalpattern
common EBSD_anglearrays, euler, quaternions
common EBSDmasks, circularmask
common EBSD_rawdata, accum_e, accum_z, mLPNH, mLPSH
common getenv_common, librarylocation

status = 1

if (SEMdata.EMsoftpathname eq 'path_unknown') then SEMdata.EMsoftpathname = Core_getenv()

; check whether the mask needs to be recomputed or not
s = size(circularmask)
dbin = 2^SEMdata.detbinning
sm = min( [SEMdata.detnumsx/dbin, SEMdata.detnumsy/dbin] )
if (s[0] ne sm) then begin
  d = shift(dist(sm),sm/2,sm/2)
  d[where(d le sm/2)] = 1.0
  d[where(d gt sm/2)] = 0.0
  circularmask = fltarr(SEMdata.detnumsx/dbin, SEMdata.detnumsy/dbin)
  if (sm eq SEMdata.detnumsx/dbin) then begin
    dm = (SEMdata.detnumsy/dbin - sm)/2
    circularmask[0,dm] = d
  end else begin
    dm = (SEMdata.detnumsx/dbin - sm)/2
    circularmask[dm,0] = d
  end
endif

ipar = replicate(0L,80)    ; parameter defined in getEBSDPatternsWrapper in mod_wrappers.f90
ipar[0] = SEMdata.mcnsx
ipar[8] = SEMdata.numset
ipar[11] = SEMdata.mcenergynumbin
ipar[16] = SEMdata.mpimx
ipar[18] = SEMdata.detnumsx
ipar[19] = SEMdata.detnumsy
ipar[20] = SEMdata.numangles
ipar[21] = SEMdata.detbinning
ipar[22] = SEMdata.detnumsx/SEMdata.detbinning
ipar[23] = SEMdata.detnumsx/SEMdata.detbinning
ipar[24] = 1
if (SEMdata.angletype eq 'qu') then ipar[24] = 0 
ipar[25] = 0

fpar = replicate(0.0,80)   ; parameter defined in getEBSDPatternsWrapper in mod_wrappers.f90 
fpar[0] = SEMdata.mcvangle
fpar[1] = SEMdata.detomega
fpar[14] = SEMdata.detxpc + float(SEMdata.detxs) * SEMdata.detxss / SEMdata.detdelta
fpar[15] = SEMdata.detypc - float(SEMdata.detys) * SEMdata.detyss / SEMdata.detdelta
fpar[16] = SEMdata.detdelta
fpar[17] = SEMdata.dettheta
fpar[18] = SEMdata.detL
fpar[19] = SEMdata.detbeamcurrent
fpar[20] = SEMdata.detdwelltime
fpar[21] = SEMdata.gammavalue

callname = 'getEBSDPatternsWrapper'

if keyword_set(single) then begin
; and here is the (single) quaternion that represents the Euler angle triplet
  quats  = Core_eu2qu( [SEMdata.detphi1, SEMdata.detphi, SEMdata.detphi2] )
  quats = reform(quats,4,1)

  ipar[20] = 1L

; initialize the simulated pattern array
  EBSDpattern = fltarr(SEMdata.detnumsx,SEMdata.detnumsy)
  EBSDpattern = reform(EBSDpattern,SEMdata.detnumsx,SEMdata.detnumsy,1)

; call the EMsoft wrapper routine from EMdymod.f90
  if (!version.os eq 'darwin') then begin
    res = call_external(librarylocation+'/libEMsoftLib.dylib', callname, $
                        ipar, fpar, EBSDpattern, quats, accum_e, mLPNH, mLPSH, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)
  endif

  if (!version.os eq 'Win32') then begin
    res = call_external(librarylocation+'/EMsoftLib.dll', callname, $
                        ipar, fpar, EBSDpattern, quats, accum_e, mLPNH, mLPSH, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)
  endif

  if (!version.os eq 'linux') then begin
    res = call_external(librarylocation+'/libEMsoftLib.so', callname, $
                        ipar, fpar, EBSDpattern, quats, accum_e, mLPNH, mLPSH, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)
  endif

  if (res ne 1.0) then begin
    Core_print,'getEBSDPatternsWrapper return code = '+string(res,format="(F4.1)")
    status = 0
  end 

; IDL takes the image origin at the bottom left corner, so we need to flip this pattern vertically 
  pattern = reverse(reform(EBSDpattern),2)
  pattern = reform(pattern,SEMdata.detnumsx,SEMdata.detnumsy,1)

end else begin ; computation of multiple EBSDpatterns

  if (SEMdata.numangles gt 500) then begin
    Core_Print,'',/blank
    Core_Print,'You are computing more than 500 EBSPs; this will take a while...'
    Core_Print,'The program will not provide any further updates until the run has been completed.'
    Core_Print,'',/blank
  endif

  pattern = fltarr(SEMdata.detnumsx,SEMdata.detnumsy,SEMdata.numangles)

; initialize the simulated pattern array
  EBSDpattern = fltarr(SEMdata.detnumsx,SEMdata.detnumsy,SEMdata.numangles)

; call the EMsoft wrapper routine from EMdymod.f90
  if (!version.os eq 'darwin') then begin
    res = call_external(librarylocation+'/libEMsoftLib.dylib', callname, $
                        ipar, fpar, EBSDpattern, quaternions, accum_e, mLPNH, mLPSH, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)
  endif

  if (!version.os eq 'Win32') then begin
    res = call_external(librarylocation+'/EMsoftLib.dll', callname, $
                        ipar, fpar, EBSDpattern, quaternions, accum_e, mLPNH, mLPSH, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)
  endif

  if (!version.os eq 'linux') then begin
    res = call_external(librarylocation+'/libEMsoftLib.so', callname, $
                        ipar, fpar, EBSDpattern, quaternions, accum_e, mLPNH, mLPSH, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)
  endif

  if (res ne 1.0) then begin
    Core_print,'getEBSDPatternsWrapper return code = '+string(res,format="(F4.1)")
    status = 0
  end 

; IDL takes the image origin at the bottom left corner, so we need to flip this pattern vertically 
  for ii=0,SEMdata.numangles-1 do begin
    slice = reverse(reform(EBSDpattern[0:*,0:*,ii]),2)
    pattern[0:*,0:*,ii] = slice[0:*,0:*]
  endfor
endelse


end

