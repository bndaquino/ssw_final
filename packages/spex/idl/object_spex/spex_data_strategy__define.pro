;+
; Time-stamp: <Tue Jul 26 2005 17:05:25 csillag tounesol.gsfc.nasa.gov>
; Modifications:
;  15-Jul-2004, Kim.  Added spex_interval_filter
;  24-Nov-2004, Sandhia.  Added a new field called spex_data_name to the argument
;                           list for read_data (which is instrument-specific).
;                           It will return the name of the instrument.
;	21-Jul-2005, Kim.  Added /this_class in setunits and getunits for speed
;	17-Aug-2005, Kim.  Added spex_deconvolved and spex_pseudo_livetime to calling args for
;		read_data, preview, get_spectra, and setinfoparams
;	Aug-2005 - Kim, Andre.  Lots of changes to accommodate image cube spectra in process,
;	    setinfoparams, setspectra methods.
;	26-May-2006, Kim.  Added spex_data_pos to args.
;	23-Jun-2006, Kim.  Added get_source_pos method
;	30-Jul-2006, Kim.  Got rid of preview method.  Initialize spex_tband to time range
;	  divided into 5 intervals.  Added GET method so we could apply spex_ut_offset to
;	  spex_ut_edges whenever we get them.
;	2-Aug-2006, Kim.  In process method, read file only if force or filename changed.
;	  Otherwise, apply timeshift and save spex_ut_edges, and set_last_update flag
;	8-Aug-2006, Kim.  If data has fewer than 10 time or energy intervals, set
;	  spex_tband and spex_eband to the raw intervals.
;	30-Aug-2006, Kim. Corrected test for eband inside edg in setinfoparams method
;	20-Nov-2006, Kim. In process, if spectrum undefined after reading, throw exception
; 16-Jun-2008,  Kim. spex_specfile is now allowed to be an array so self.specfile must now be 
;   a pointer (was a scalar string). Also, use same_data to compare old and new values.
; 29-Jul-2008, Kim. In setinfoparams, changed kev to keV
; 12-Aug-2009, Kim.  Added cleanup method
; 28-Nov-2009, Kim.  Added self.data_sel, and check if it's changed in process
; 19-Feb-2010, Kim.  In process, if more than one input file, do loop here (previously did it in 
;   messenger_specfile read_data), also use accum_time to save a subset of what was in file. Also,
;   default eband bins are now defined by dividing en range into 5 logarithmically spaced bins.
; 07-Jul-2010, Kim. In process, if spex_accum_time is invalid for input file, set to 0.,0.
; 14-Feb-2011, Kim. In setinfoparams, for fermi_gbm, don't just create spex_eband by dividing energy
;  range into 5 bins, hard-code them to bins that R. Schwartz recommended.
; 11-Aug-2011, Kim. Added spex_def_eband arg to read and setinfoparams calls (so default bands are
;  now set in the data-specific read routines), and use them if eband not set otherwise. (at some point,
;  should probably take eband and tband out of strategy info and put in data info)
; 17-Aug-2011, Kim. When reading data from multiple files, only append part of array from new file where
;  times are later than old file's times.  Fermi gbm has overlapping data in files.
; 12-Feb-2014, Kim. In process, when looping through files, if one file is bad or had no good data,
;  keep going to next file. Previously aborted.
; 25-Nov-2015, Kim. Fixed error when computing spex_eband when spex_def_eband is not defined 
;  (was using log(edg) instead of edg)
; 08-Sep-2016, Kim. In setinfoparams, reset spex_eband when min or max of spex_eband is outside range of new edges
; 19-Jan-2023, Kim. Added get_time_plot_args method, just returns -1. But if any data class (e.g. minxss) has its own
;   get_time_plot_args method, it will be called instead of this one to allow plot args specific to that class.
;---------------------------------------------------------------------------

function spex_data_strategy::init, source = source, _extra = _extra

ret = self->framework::init( source = source, $
                             control = spex_data_strategy_control(), $
                             info = spex_data_strategy_info() )

self->set, _extra = _extra
self.specfile = ptr_new('')
self.data_sel = ''

return, ret

end

;------------------------------------------------------------------------------

pro spex_data_strategy::cleanup
free_var, self.specfile
self->framework::cleanup
end
 
;------------------------------------------------------------------------------
; the method is necessary so that, e.g. hessi_multi_image strategy can have its own
; preview of the input file without doing a read_data (which asks user for ROI and
; constructs spectra in ROI).

;pro spex_data_strategy::preview, $
;                      spectrum,  errors,  livetime,  $
;                      spex_respinfo, spex_file_time,  spex_ut_edges,  spex_ct_edges,  $
;                      spex_area, spex_title, spex_detectors,  $
;                      spex_interval_filter, spex_units, spex_data_name, $
;                      spex_deconvolved, spex_pseudo_livetime, spex_data_pos, $
;                      err_code=err_code
;
;self -> read_data, $
;  spectrum,  errors,  livetime,  $
;  spex_respinfo, spex_file_time,  spex_ut_edges,  spex_ct_edges,  $
;  spex_area, spex_title, spex_detectors,  $
;  spex_interval_filter, spex_units, spex_data_name, $
;  spex_deconvolved, spex_pseudo_livetime, spex_data_pos, $
;  err_code=err_code
;
;end

;---------------------------------------------------------------------------------

pro spex_data_strategy::setinfoparams, $
    spex_respinfo, spex_file_time,  spex_ut_edges,  spex_ct_edges,  $
    spex_area, spex_title, spex_detectors,  $
    spex_interval_filter, spex_units, spex_data_name, $
    spex_deconvolved, spex_pseudo_livetime, spex_data_pos, spex_def_eband=spex_def_eband

units_str = {spex_units}
units_str.data_name = spex_data_name
units_str.data_type = 'counts'
units_str.data = spex_units
units_str.drm = 'counts photons!u-1!n'
units_str.area = 'cm!u-2!n'
units_str.energy = 'keV!u-1!n'
units_str.time = 's!u-1!n'
units_str.ct_edges = 'keV'
units_str.ut_edges = 'UTC'
units_str.ph_edges = 'keV'

; if default eband not passed in, set it to 5 even log bins
edg = minmax(spex_ct_edges)
if ~exist(spex_def_eband) then begin
  logedg = alog10(edg)
  spex_def_eband = get_edges(10.^( logedg[0] + indgen(6) * (logedg[1]-logedg[0]) / 5. ), /edges_2)
endif

; set the info parameters --- originally for spex_hessi_specfile
self -> Set, spex_respinfo = spex_respinfo, $
  spex_file_time = spex_file_time, $
  spex_ut_edges = spex_ut_edges, $
  spex_ct_edges = spex_ct_edges, $
  spex_area = spex_area, $
  spex_title = spex_title, $
  spex_file_units = units_str, $
  spex_detectors = spex_detectors, $
  spex_interval_filter = spex_interval_filter, $
  spex_deconvolved = spex_deconvolved, $
  spex_pseudo_livetime = spex_pseudo_livetime, $
  spex_data_pos = spex_data_pos, $
  spex_def_eband = spex_def_eband

; set the info parameters --- originally for spex_data

; If eband is already defined, only keep the bands that are within new edge boundaries.
; If none qualify, or it wasn't defined, then 
;   If <= 10 energy bins, set eband to energy edges
;   Use spex_def_eband default ebands for this instrument

eband = self -> get(/spex_eband)
if eband[0] ne -1 then begin
  ; If spex_eband was already set, reset them only if any of them lie outside the range of the new edges.
  ; Want to preserve user's values if they've set them explicitly, when it seems reasonable.
  ; Prevously was leaving them if any overlapped with new edges (and only keeping the overlapping ones),
  ; but this cause a problem with the konus _1 and _2 files.
  if min(eband) lt min(edg) or max(eband) gt max(edg) then eband = -1 
;    eband = eband > min(edg) < max(edg)
;    q = where ((eband[1,*]-eband[0,*]) gt 0., count)
;    eband = count gt 0 ? eband[*,q] : -1
endif
if eband[0] eq -1 then if n_elements(spex_ct_edges[0,*]) le 10 then eband = spex_ct_edges
if eband[0] eq -1 then eband = spex_def_eband
self -> set, spex_eband= eband

if n_elements(spex_ut_edges[0,*]) le 10 then begin
	self -> set, spex_tband = spex_ut_edges
endif else begin
	tim = minmax(spex_ut_edges)
	self -> set, spex_tband = get_edges(tim[0] + indgen(5) * (tim[1]-tim[0]) / 4., /edges_2)
endelse

end

;---------------------------------------------------------------------------------

;function spex_data_strategy::get, spex_ut_edges=spex_ut_edges, _extra=_extra, found=found, not_found=not_found
;
;ret = self->framework::get(spex_ut_edges=spex_ut_edges, _extra=_extra, found=found, not_found=not_found)
;if keyword_set(spex_ut_edges) then begin
;	ut_offset = self->framework::get(/spex_ut_offset)
;	if is_struct(ret) then ret.spex_ut_edges = ret.spex_ut_edges + ut_offset else $
;		ret = ret + ut_offset
;endif
;
;return, ret
;end

;---------------------------------------------------------------------------------

pro spex_data_strategy::process, force=force, _extra = _extra

IF Keyword_Set( _EXTRA ) THEN self->Set, _EXTRA = _extra

dprint,'in spex_data_strategy::process'

err_code = 0
err_msg = ''

if keyword_set(force) or $
 not same_data(self->get( /spex_specfile ), *(self.specfile)) or $
 not same_data(self->get( /spex_data_sel ), self.data_sel) or $
 not same_data(self->get(/spex_accum_time), self.accum_time) then begin
	;print, 'in spex_data_strategy::process for new input file(s)'
	
	nfiles = n_elements(self->get(/spex_specfile))
	
	for i=0,nfiles-1 do begin
	  dprint, 'in spex_data_strategy::process, reading data file... 
  	self->read_data, file_index=i, spectrumi,  errorsi,  livetimei,  $
  	  spex_respinfo, spex_file_time,  spex_ut_edgesi,  spex_ct_edges,  $
  	  spex_area, spex_title, spex_detectors,  $
  	  spex_interval_filter, spex_units, spex_data_name, $
  	  spex_deconvolved, spex_pseudo_livetime, spex_data_pos, spex_def_eband=spex_def_eband, $
  	  err_code=err_code
  
  	if err_code then begin
  	  message,'Spectrum file is invalid.', /cont
  	  goto, nextfile
  	endif
  	  
  	if n_elements( spectrumi ) eq 0 then begin
  	  message,'No spectrum data for file.', /cont
  	  goto, nextfile
  	endif
  	
  	if ~exist(spectrum) then begin
  	  spectrum = spectrumi
      errors = errorsi
      livetime = livetimei
      spex_ut_edges = spex_ut_edgesi
  	endif else begin
  	  q = where (spex_ut_edgesi[0,*] - last_item(spex_ut_edges) ge -.0005, count)
  	  if count gt 0 then begin
  	    spectrum = [[spectrum], [spectrumi[*,q[0]:*]]]
        errors = [[errors], [errorsi[*,q[0]:*]]]
        livetime =  [[livetime], [livetimei[*,q[0]:*]]]
        spex_ut_edges = [[spex_ut_edges], [spex_ut_edgesi[*,q[0]:*]]]
      endif
    endelse
    spex_file_time = minmax(spex_ut_edges)
    
    nextfile:
  endfor
 
  if n_elements( spectrumi ) eq 0 then message,'No spectrum data in any of the input files.'
  
  accum_time = self->get(/spex_accum_time)
  t_msg = ''
  if valid_time_range(accum_time) then begin
    self.accum_time = accum_time
    q = where_within (spex_ut_edges, accum_time, count)
    if count gt 0 then begin
      spectrum = spectrum[*,q]
      errors = errors[*,q]
      livetime = livetime[*,q]
      spex_ut_edges = spex_ut_edges[*,q]
    endif else t_msg = 'No times found in your requested accum_time: ' + format_intervals(accum_time,/ut) + ' Reading entire file.'
  endif else if ~same_data(accum_time,[0.d,0.d]) then t_msg = 'Invalid accum_time, Reading entire file.'
    
  if t_msg ne '' then begin
    self->set, spex_accum_time = [0.d,0.d]
    message, t_msg, /cont
  endif
  
	self->setinfoparams, spex_respinfo, spex_file_time,  spex_ut_edges,  spex_ct_edges,  $
	    spex_area, spex_title, spex_detectors,  $
	    spex_interval_filter, spex_units, spex_data_name, $
	    spex_deconvolved, spex_pseudo_livetime, spex_data_pos, spex_def_eband=spex_def_eband

	self.ut_offset = 0.

	self->setspectra, spectrum, errors, livetime

	*(self.specfile) = self->get( /spex_specfile )
	self.data_sel = self->get( /spex_data_sel )
	
endif

new_ut_offset = self->get(/spex_ut_offset)
if new_ut_offset ne self.ut_offset then begin
	diff = new_ut_offset - self.ut_offset
	;print, 'UT offset.  Old, new, diff', self.ut_offset, new_ut_offset, diff
	self -> set, spex_ut_edges = self->get(/spex_ut_edges) + diff
	self.ut_offset = new_ut_offset

	print, 'Warning !!!!!!!!'
	print,'Shifting time array for data (spex_ut_edges).  If you have set '
	print,'  other time intervals (bk, fit) they will not be changed, but the '
	print,'  data in those intervals is now different and will be reaccumulated.'

	self->framework::set_last_update	;need to do this explicitly since not calling setdata

endif

end

;---------------------------------------------------------------------------------

pro spex_data_strategy::setspectra, spectrum, errors, livetime

; Any error bars or data points with zero livetime are set to zero.
zero_index =  where( livetime eq 0,  count )
if count gt 0 then begin
  errors[zero_index] =  0
  spectrum[zero_index] =  0
endif

; kim, 4/4/05 reform because for single energy, imagecube returns [1,n], spectrum returns [n]
; so we'll save as [n]
final_data =  {data: reform(spectrum), $
	edata: reform(errors), $
	ltime: reform(livetime)}

units_str = self->get( /spex_file_units )

; want to always store data in counts units.
; first have to set current units into this object's units structure
final_data = self -> convert_units(final_data, 'counts', units_str)

; two sets of units_str to set
self -> setunits, units_str, /orig
self -> setunits, units_str

;stop
self -> framework::setdata, final_data

end

;--------------------------------------------------------------------------

function spex_data_strategy::getunits, orig = orig

  if keyword_set(orig) then return, self -> get(/spex_data_origunits, /this_class)
  return, self -> get(/spex_data_units, /this_class)

end

;--------------------------------------------------------------------------

pro spex_data_strategy::setunits, units, orig = orig

if keyword_set(orig) then self -> set, spex_data_origunits=units, /this_class else $
    self -> set, spex_data_units=units, /this_class

end

;-------------------------------------------------------------------------------

function spex_data_strategy::get_source_pos
return, [-9999.,-9999.]
end


;-------------------------------------------------------------------------------

function spex_data_strategy::get_time_plot_args  ; kim added this function
return, -1
end

;-------------------------------------------------------------------------------

pro spex_data_strategy__define

dummy = {spex_data_strategy, $
         specfile: ptr_new(), $
         data_sel: '', $
         accum_time: [0.d, 0.d], $
         ut_offset: 0., $
         inherits spex_gen }

end
