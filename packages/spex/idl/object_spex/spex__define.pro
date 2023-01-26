;+
; Name: SPEX__DEFINE
;
; Purpose: Main OSPEX object
;
; Category: OSPEX
;
; Written: 2003, Kim Tolbert
; Modifications:
;   22-Jun-2004, Kim.  Replaced spex_interval_range with spex_intervals_tofit
;   15-Jul-2004, Kim.  Added _extra to fitsummary and setupsummary so now they can
;     pass that to text_output (for /print and /file_text for printing and sending to file)
;   16-Jul-2004, Kim.  Moved preview stuff out of here and into spex_data and spex_drm - just
;     call those from this preview
;   20-Jul-2004, Kim.  Added filter info to fitsummary and setupsummary output.
;   09-Aug-2004, Sandhia,  Added a method setParams to allow user to set parameters for OSPEX manually.
;   26-Aug-2004, Sandhia,  Modified savefit to write fit parameters to a FITS file or IDL save
;                          file based on a keyword value.
;   08-Sep-2004, Sandhia,  implemented code to write the spex_summ structure to a fits file in
;                             the format compatible with a RATE FITS file.  The new fits file will
;                             now have three extensions - primary, rate and ebounds.
;   09-Sep-2004, Sandhia,  restorefit (from FITS file option) will now call spex_read_fit_results routine
;                             to copy the FITS results in OSPEX_summ structure.
;   12-Sep-2004, Kim.  writescript wasn't writing bk_time_interval if bk_sep=1. Fixed.
;   17-Sep-2004, Kim.  Added writing bk_rate and bk_error to FITS file
;   19-Sep-2004, Kim.  Added 4 background things to plot in plot_summ, and
;     2 bk things to calc in calc_summ.
;   20-Sep-2004, Kim.  Added an extension to fit results FITS file for OSPEX control params
;   24-Sep-2004, Sandhia, Filled mjdref, timezero, tstarti, tstopi, tstartf and tstopf, poisserr,
;                              timesys, clockcor, observer and version fields of
;                              rate_struct before calling mk_rate_hdr.  Deleted code to write these
;                              parameters to fits file from add_rate_keywords routine.
;   5-Oct-2004, Kim.  Changed writescript to write a procedure instead of a script that is @'ed.
;   5-Oct-2004, Kim.  Added runscript method to set params from a script with option to init first.
;   5-Oct-2004, Kim.  Added adjust_intervals method to force fit intervals to data boundaries
;   20-Oct-2004, Sandhia,  Separated the code to handle writing fits files (data or
;                              fit results) in a new file spex__fitswrite
;   2-Dec-2004, Kim.  Added chisq_full item to calculate in calc_summ method.
;   11-Apr-2005, Kim.  In intervals, use plot_time,/data instead of plot, class=spex_data, and
;     added /show_filter on time plot
;   23-Jun-2005, Kim. Changed file filter for pickfile to include .dat, .fits for data, .rmf for srm
;	10-Aug-2005, Kim.  Added area to fitresults display, added current date to default name for
;	  script, and added some comments at the beginning of the script written to explain how to
;	  use it.  And added a ' READING...' while reading the input file.
;	Sep-2005, Kim.  Use xmessage for preview and fitsheader instead of prstr
;	30-Sep-2005, Kim.  Added roi, roi_config, and save_roi methods - just so you can call
;	  those spex_image methods from a spex object.  Just passes through.
;	  Also, in writescript, if obj is using an image cube, then temporary obj used to
;	  get default params must have strategy set to an image also, or params aren't there.
;	7-Oct-2005, Kim.  Added list_roi method, so can call spex_image::list_roi from a spex obj.
;	9-Feb-2006, Kim.  Added set method - allows us to trap when a new input file is set, so we
;	  can reset key parameters to defaults.
;	  Added albedo info to setupsummary
;	  get_script_params method renamed to get_param_names, and added new_input keyword (also added
;	    an 'N' flag to param list in ospex_script_params.txt to indicate which params should be
;	    reset when input changes, and now read that file with rd_tfile to parse it)
;	  clearall method renamed to init_params, and added new_input keyword.
;	  Added is_image_input method to test if input file is an image cube and use it in the roi methods.
;	15-Mar-2006, Kim.  Added epsilon=1.e-5 to call to find_contig_ranges
;	7-Apr-2006, Kim.  Added call to func_comp_kw in fitsummary, added list_function method,
;	  and call it in setupsummary.
;	9-May-2006, Kim.  Disabled option to write fit results in save (genx) file
;   23-Jun-2006, Kim. Set default for source_xy whenever a new input file is set (in set method)
;	2-Jul-2006, Kim.  Added tband to intervals method.
;	17-Jul-2006, Kim.  To enable SOXS data:
;		In set_file, added .out as acceptable file name extension
;		In fitsheader, check if file is a FITS file before calling mrdfits
;	3-Aug-2006, Kim.  Added nodialog keyword to restorefits method
;	18-Sep-2006, Kim.  Added bk_ratio info to setupsummary output
;	19-Oct-2006, ABShah Change *.out to *.LES for SOXS instrument.
;   12-Apr-2007, Kim.  Added call to spex_old_defaults
;	31-May-2007, Kim.  Added find_bad_intervals wrapper for spex_fitint::find_bad_intervals
;	13-Jun-2007, Kim. For Messenger data, added .csv extension in set_file
; 9-May-2008, Kim. In plot_summ, added panel_id tag to plot_params so it can be separate from id (title)
; 16-Jun-2008, Kim. Spex_specfile is now allowed to be an array of files (for MESSENGER).  In cases where
;   we're just testing the file, changed to file[0]. Change equality test to same_data call.
; 24-Jun-2008, Kim. Default script name now lowercase.
; 25-Jun-2008, Kim. Added write_model_textfile method.
; 27-Aug-2008, Kim. In plot_summ, allow summing plot (dim1_enab_sum=1)
; 24-Nov-2008, Kim. In plot_summ, added interval # and energy range n to plot items.  In calc_summ, added
;   enrange to calc items.
; 3-Dec-2008, Kim. In fitsummary, added parameter number column at far left.
; 4-Dec-2008, Kim. In plot_summ, added plot_free field to table, and now plot bar showing free mask for params that
;     have plot_free set to 1.  In write_model_textfile, use energies from calc_func_components structure.
; 5-Dec-2008, Kim. In calc_summ and plot_summ, added thick_energy.  Also, changed order of table in plot_summ.
; 10-Dec-2008, Kim. Call brmthick_power with eps=1.e-6
; 11-Dec-2008, Kim. Added plot_summ_map method, and added ff_map and enrange_map to plot_summ choices
; 12-Dec-2008, Kim. Added plaw_atbreak to calc_summ and plot_summ, and show func component name on free/fixed map
; 17-Dec-2008, Kim. Changed plaw_atbreak to plaw_at50
; 16-Jan-2009, Kim.  Major changes to calc_summ,plot_summ, plot_summ_tables
; 19-Feb-2009, Kim.  Minor change to chisq_full calc in calc_summ
; 18-Mar-2009, KIm.  Added show_free option in plot_summ
;  2-Jul-2009, Kim.  Added thick2 for Non-thermal Electron Energy Flux option in plot_summ_tables and calc_summ methods
; 17-Aug-2009, Kim.  In plot_summ_map, destroy spex_map object. memory leak.
; 26-Sep-2009, Kim.  In calc_summ, added doc, made it so old item names translate to new ones, just use new ones.
; 28-Oct-2009, Kim.  In set_file, added .rsp and .pha extensions for fermi gbm files.
; 28-Nov-2009, Kim.  In set_file, added wda* to filters for pickfile.  Changed filter from
;   extensions to full names.  Don't check if file chosen has one of allowed extensions, since now
;   'wda*' is allowed.  If file isn't right type, reader will report error. Also in set_file, added
;   /dialog on call to set specfile, for popup for data type selection.
; 03-Feb-2010, Kim. In fitsheader, use is_fits to test if file is fits format. Previously looked for .fits in filename.
; 18-Jun-2010, Kim. In plot_summ_tables and calc_summ, fixed error.  Normalized residuals are what's stored in
;   spex_summ, raw residuals are calculated from that times the errors.  I had it backward. (so norm_resid as something
;   to calculate is now raw_resid.)
; 30-Jun-2010, Kim.  In set_file, added rsp2 extension for gbm drm files.
; 07-Jul-2010, Kim.  Added set_time method.  Use locate_file to find file.
; 25-Jan-2011, Kim. In set_file added '.fit' for spec file types
; 9-Jun-2011, Kim. In plot_summ_tables and calc_summ, added thick_integ_flux for thick2_vnorm component,
;   and made thick_energy work for thick2_vnorm as well as other thicks.  Fixed error in calc_summ for
;   flux_at_e when single time interval, but two energies.  Also in plot_summ, if ask for time plot,
;   but single time, print message and values but no plot.
; 19-Oct-2011, Kim. Multiple filters to dialog_pickfile don't work on unix unless they are in the form
;   of a string with no space, separated by ';'.  Previously used string array (works on Window, not Unix).
; 10-Nov-2011, Kim. Added plot_log to plot_summ_tables vtime structure, and set default for log scale to
;   1 for all normalization parameters. Use plot_log flag in plot_summ method.
; 25-Sep-2012, Kim. In plot_summ_tables, added units for conv, resid, and raw_resid. Added check in calc_summ
;   for whether any spex_summ values have been stored yet.  Also when calculating chisq_full, make sure chisq
;   isn't 0.
; 08-Feb-2013, For thick2_rc, added same special options as for thick2_vnorm
; 01-Jul-2013, In SET, check env var OSPEX_NEWFILE_PARAMS_INIT do decide whether to init params when input
;   file is changed.  If not set or has a value of 1, do init params. Normally safer to initialize, but if user
;   knows new files is compatible with old (e.g. same time interval, different detector), then OK.
;   setenv,'OSPEX_NEWFILE_PARAMS_INIT=0'  to turn off initializing params.  User must remember to turn back on.
; 28-Aug-2013, Kim. In fitsummary, added a column showing ndex number within each component, with -0- for first.
; 27-Sep-2013, Kim. In setupsummary, changed output for background setup for changes in profile methods.
; 18-Mar-2014, Kim. In calc_summ, for flux_at_e, make bin .02 keV wide around user's requested energy (was 1 keV)
;   Also, change label for flux_at_e option. In plot_summ, convert flux units to use positioning commands (!u, !n).
; 19-Mar-2014, Kim. Added ct_flux_at_e option (Time Profile for Count Flux at E keV for Model). Just added for All components.
; 11-Jul-2014, Kim. Added items 'edf' and 'therm_energy_plasma' to values to calculate in calc_summ, and added them to
;   tables in plot_summ_tables. In plot_summ, if only one interval for edf then call fit object's plot_edf method to get
;   nice EDF plot, otherwise overlay edf for each interval on one plot. Use func keyword on calls to fint_from_fref.
; 14-Jul-2014, Kim. In calc_summ, call calc_nontherm_electron_energy_flux instead of fint_from_fref and brmthick_power,
;   and check if strat keyword input for 'thick_energy' item.
; 19-Sep-2014, Kim. In calc_summ, added 'integ_photon_flux' item.  Also added it to vtime table in plot_summ_tables.
; 13-Jan-2015, Kim. Added restore_roi here, so users can call from ospex top obj, not from the data strategy object
; 30-Oct-2017, Kim. In plot_summ_tables method, added plot_bk flag to ven structure to indicate which item name should be used
;   to plot background for each item (if left blank, don't plot background for this item).  In plot_summ method, added
;   overlay_bk keyword to plot background (using new ven tag) on data plots (counts or photon, rate or flux).
; 14-Nov-2017, Kim. In plot_summ, in call to do_plot, changed keyword from plotman to use_plotman to match change in spex_gen
; 09-Jul-2018, Kim. Added textfile method
; 06-Dec-2019, Kim. Added dct_flux_at_e item to calc_summ and plot_summ_tables, so we can plot data count flux vs
;   time with error bars in 'view spex fit results'.
; 06-Apr-2022, Kim. When calling therm_energy_plasma, don't multiply em by 1.d49.
; 19-Jan-2023, Kim. Added .sav to types of files to show in dialog when browsing for input file (for minxss)
;-
;---------------------------------------------------------------------------

function spex::init, source=source, no_gui=no_gui, _extra=_extra

  if not since_version('5.6') then begin
    print,'SORRY - you must have IDL Version 5.6 or later to run OSPEX.  Aborting.'
    return,0
  endif

  if not keyword_set( source ) then source = obj_new( 'spex_fit', _extra=_extra )

  ret = self->framework::init( source = source, $
    ;control={spex_control}, $
    ;info = {spex_info}, $
    _extra = _extra )

  spex_old_defaults, self, _extra=_extra

  if not keyword_set(no_gui) then self -> gui, _extra=_extra

  return, ret
end

;--------------------------------------------------------------------------

pro spex::set, _extra=_extra

  old_specfile = self->get(/spex_specfile)

  self -> framework::set, _extra=_extra

  ;if old_specfile eq '' then return

  if not same_data(old_specfile, self->get(/spex_specfile)) then begin
    ; if env var OSPEX_NEWFILE_PARAMS_INIT isn't set at all, or is set to 1, initialize parameters when changing to new input file
    ; if OSPEX_NEWFILE_PARAMS_INIT is set to 0, don't init ospex params
    env = getenv('OSPEX_NEWFILE_PARAMS_INIT')
    params_init = env eq '' ? 1 : fix(env)
    if old_specfile[0] ne '' and params_init then self -> init_params, /new_input, param_names=param_names
    self -> set, spex_source_xy = (self->get(/obj,class='spex_data')) -> get_source_pos()
  endif

end

;--------------------------------------------------------------------------
; Helper routine for setting srm or spectrum file.  If dialog set, pops up widget dialog to browse to file.
; Tries to read file, and if there's a problem, doesn't set that file into ospex.
; When setting spectrum file, checks to see if should also set drm file.

pro spex::set_file, file, accum_time=accum_time, $
  srm=srm, $
  drm=drm, $  ; same as srm
  dialog=dialog, $
  FOUND_FILE=found_file, $
  TITLE=title, $
  NOINTERACTIVE=nointeractive, $
  QUIET=quiet, $
  _REF_EXTRA=_ref_extra

  quiet = keyword_set(quiet)
  checkvar, dialog, 0

  checkvar, srm, drm
  checkvar, title, srm ? 'Please Select SRM file' : 'Please Select Spectrum or Image File'

  ; on unix this array of filters doesn't work, must be string, separated with ;, no spaces
  ;filter = srm ? ['*.fits', '*.rmf', '*.dat', '*.rsp', '*.rsp2'] : $
  ;  ['*.fit', '*.fits', '*.dat', '*.les', '*.csv', '*.pha', 'wda*']
  filter = srm ? '*.fits;*.rmf;*.dat;*.rsp;*.rsp*,*.sav' : $
  '*.fit;*.fits;*.dat;*.les;*.csv;*.pha;wda*'
  ;exts = srm ? ['fits', 'rmf', 'dat', 'rsp'] : ['fits', 'dat', 'les', 'csv', 'pha']
  multiple_files = srm ? 0 : 1

  found_file = self->locate_file(file, status=status)
  exists = status
  ;found_file = loc_file( file, COUNT=count )
  ;exists = count GT 0

  msg = ''
  if dialog or (~status and ~keyword_set(nointeractive))  then begin

    found_file = ssw_pickfile( file=file[0], $
      exists=exists, $
      title=title, $
      filter=filter, $
      multiple_files=multiple_files, $ ;added 16-jun-2008
      _extra=_ref_extra )

    case 1 of
      found_file[0] eq '': msg = 'ERROR: No file selected for reading.'
      exists eq 0: msg = 'ERROR: Selected file not found - ' + found_file[0]
      else: msg = ''
    endcase

  endif else begin
    if not exists then msg = 'ERROR: Specified file not found - ' + file
  endelse

  if msg ne '' then begin
    message, msg, /cont
    return
  endif

  if srm then begin
    oldfile = self -> get(/spex_drmfile)
    self -> set, spex_drmfile = found_file[0]
    dummy = self -> getdata(class='spex_drm')
    if dummy[0] eq -1 then self -> set, spex_drmfile=oldfile
  endif else begin
    oldfile = self -> get(/spex_specfile)
    if not spex_get_nointeractive() then xmessage,['', '      Reading...       ', ''], wbase=wxmessage
    self -> set, spex_specfile = found_file, /dialog
    if exist(accum_time) then self -> set, spex_accum_time=accum_time
    dummy = self -> getdata(class='spex_data')
    if xalive(wxmessage) then widget_control, wxmessage, /destroy
    if size(/tname, dummy) ne 'STRUCT' then self -> set, spex_specfile=oldfile else begin
      spex_respinfo = self -> get(/spex_respinfo)
      drmfile = is_string(spex_respinfo[0]) ? loc_file(spex_respinfo) : ''
      self -> set, spex_drmfile=drmfile
    endelse
  endelse

end

;--------------------------------------------------------------------------
; Helper routine when setting spex_accum_time.  After reading input file for new
; time, looks for correct response file to set (for gbm, depends on time as well as
; input file)
pro spex::set_time, time_range, status=status

  status = 0
  if ~valid_time_range(time_range) then begin
    message, /cont, 'Invalid time range. '
    return
  endif

  file_time = self -> get(/spex_file_time)
  if total(file_time) eq 0. then begin
    message, /cont, 'You must first select an input file.'
    return
  endif

  if has_overlap(time_range, file_time) then begin
    status = 1
    self -> set, spex_accum_time = time_range
    dummy = self -> getdata(class='spex_data')
    if is_struct(dummy) then begin
      spex_respinfo = self -> get(/spex_respinfo)
      drmfile = is_string(spex_respinfo[0]) ? spex_respinfo : ''
      self -> set, spex_drmfile=drmfile
    endif
  endif else message,'Requested time not contained in file or no file set.  Not using.', /cont

end

;--------------------------------------------------------------------------

function spex::getunits, class_name=class_name, _extra=_extra

  if keyword_set(class_name) then begin
    obj = self -> get(/obj, class_name=class_name)
    return, obj -> getunits(_extra=_extra)
  endif

  return, -1
end

;--------------------------------------------------------------------------

function spex::getdata, $
  class_name=class_name, $
  _extra=_extra

  @spex_insert_catch

  checkvar, class_name, 'spex_data'

  data = self -> framework::getdata(class_name=class_name, _extra=_extra)

  return, data

end

;--------------------------------------------------------------------------

pro spex::cleanup

  plotman_obj = self -> get(/spex_plotman_obj)
  if obj_valid(plotman_obj) then obj_destroy, plotman_obj

  self -> framework::cleanup

end

;--------------------------------------------------------------------------
; loop only applies to defining bk_time when bk_sep is set - loops through the
; separate energy bands

pro spex::intervals, $
  eband=eband, $
  tband=tband, $
  erange=erange, $
  bk_time=bk_time, $
  loop=loop, $
  bk_eband=bk_eband, $
  fit=fit, $
  _extra=_extra

  case 1 of

    keyword_set(eband) : (self -> get(/obj,class_name='spex_data')) -> intervals, /energy, _extra=_extra

    keyword_set(tband) : (self -> get(/obj,class_name='spex_data')) -> intervals, _extra=_extra

    keyword_set(erange) : (self -> get(/obj,class_name='spex_fitrange')) -> intervals, _extra=_extra

    keyword_set(bk_eband) : (self -> get(/obj,class_name='spex_bkint')) -> intervals, /energy, _extra=_extra

    keyword_set(bk_time) : begin
      bk_sep = self -> get(/spex_bk_sep)
      if keyword_set(loop) and bk_sep then begin
        bands = self -> get(/spex_bk_eband)
        if bands[0] ne -1 then begin
          for ib = 0,n_elements(bands)/2-1 do begin
            ;self -> plotman, class_name='spex_data', interval=bands[*,ib]
            self -> plot_time, /data, interval=bands[*,ib], /show_filter, _extra=_extra
            self -> intervals, /bk_time, this_band=ib, _extra=_extra
          endfor
        endif else message, /cont, 'No background energy bands set.  Aborting.'
      endif else (self -> get(/obj,class_name='spex_bkint')) -> intervals, _extra=_extra
    end

    keyword_set(fit) : (self -> get(/obj,class_name='spex_fitint')) -> intervals, _extra=_extra

    else: (self -> get(/obj,class_name='spex_fitint')) -> intervals, _extra=_extra
  endcase

end

;------------------------------------------------------------------------------

pro spex::adjust_intervals
  intervals = self->get(/spex_fit_time_interval)
  index = self -> edges2index (newedges=newedges, intervals=intervals, /do_time, got_int=got_int)
  if got_int then self -> set, spex_fit_time_interval=newedges
end

;------------------------------------------------------------------------------

pro spex::xfit_comp, _extra=_extra
  (self -> get(/obj,class='spex_fit')) -> xfit_comp, _extra=_extra
end

;------------------------------------------------------------------------------

pro spex::fitsummary, out=out, _extra=_extra

  struct = self -> get(/spex_summ)

  npad = 15

  if not self->valid_summ_params() then out = 'No Fit Parameters stored yet.' else begin
    times = format_intervals(struct.spex_summ_time_interval, /ut, /prefix)
    out = ['Current Fit Results     ' + strmid (anytim(!stime, /trunc, /vms), 0, 17), $
      '', $
      'Fit Function: ' + struct.spex_summ_fit_function, $
      'Detectors Used: ' + self -> get(/spex_detectors) + $
      '    Area: ' + trim(self->get(/spex_summ_area)), $
      '', $
      strpad('', 6) + $
      strpad('Fit Params', npad) + $
      strpad('Sigma', npad) + $
      strpad('Start Param',npad) + $
      strpad('Minimum',npad) + $
      strpad('Maximum', npad) + $
      strpad('Free', npad) ]

    ntime = n_elements(times)
    nparams = n_elements(struct.spex_summ_params[*,0])

    ind = fit_function_query(struct.spex_summ_fit_function, /param_index)
    diff = ind[*,1] - ind[*,0] + 1
    for i=0,n_elements(diff)-1 do ii = append_arr(ii, indgen(diff[i]))
    index = strpad(trim(ii),3) + ' '
    q0 = where(ii eq 0)
    index[q0] = ' -0-'

    for i=0,ntime-1 do begin
      if not struct.spex_summ_fit_done[i] then goto, next
      eindex = where (struct.spex_summ_emask[*,i])
      if eindex[0] eq -1 then goto, next
      eranges = find_contig_ranges(struct.spex_summ_energy[*,eindex], epsilon=1.e-5)
      out =[out, $
        '', $
        times[i] + ', Filter: ' + trim(struct.spex_summ_filter[i]) + $
        '  Energy Range: ' + arr2str(format_intervals(eranges, format='(f9.2)')), $
        '  Chisq=' + trim(struct.spex_summ_chisq[i],'(f8.2)') + $
        '  MaxIter=' + trim(struct.spex_summ_maxiter[i]) + $
        '  #Iter=' + trim(struct.spex_summ_niter[i]) + $
        '  Uncert=' + trim(struct.spex_summ_uncert[i], '(f4.2)') + $
        '  Stop Msg=' + struct.spex_summ_stop_msg[i] ]
      if fit_comp_kw(struct, time_int=i, summary_string=ss) then out = [out, '  '+ss]

      for ip = 0, nparams-1 do begin
        out = [out, $
          strpad(trim(ip),2) + $
          index[ip] + $
          strpad(trim(struct.spex_summ_params[ip,i]), npad) + $
          strpad(trim(struct.spex_summ_sigmas[ip,i]), npad) + $
          strpad(trim(struct.spex_summ_starting_params[ip,i]), npad) + $
          strpad(trim(struct.spex_summ_minima[ip,i]), npad) + $
          strpad(trim(struct.spex_summ_maxima[ip,i]), npad) + $
          strpad(trim(struct.spex_summ_free_mask[ip,i]+0), npad) ]
      endfor
      next:
    endfor
  endelse

  ;prstr, out

  ;xmessage, out, font='fixedsys', xsize=max(strlen(out)), ysize=n_elements(out)+5, $
  ;   title='Fit Results'

  text_output, out, title='Fit Results', _extra=_extra
end

;------------------------------------------------------------------------------
; Write the model data (computed fit function components separate and combined) into a text file.
; All arguments to calc_func_components are allowed - see calc_func_components for explanation.
; In addition to those:
; out - return string array of output
; filename - file name to write to

pro spex::write_model_textfile, $
  this_interval=this_interval, $
  spex_units=spex_units, $
  photons=photons, $
  out=out, $
  filename=filename, $
  _extra=_extra

  default, spex_units, 'counts'
  default, photons, 0

  s = self -> calc_func_components (this_interval=this_interval, $
    spex_units=spex_units, photons=photons,  _extra=_extra)
  e = s.ct_energy

  if is_struct(s) and n_elements(e) gt 1 then begin

    if size(filename, /tname) ne 'STRING' then begin
      atim = strlowcase( trim( str_replace(anytim(!stime, /vms, /date), '-', '_') ))
      filename = 'ospex_model_data_'+atim+'.pro'
      filename = dialog_pickfile (path=curdir(), filter='*.txt', $
        file=filename, $
        title = 'Select output file',  $
        group=group, $
        get_path=path)
    endif

    if filename ne '' then begin

      out = ''
      npad = 15
      ncols = n_elements(s.id)
      nint = n_elements(this_interval)
      ctsorph = ['count', 'photon']

      for jint = 0,nint-1 do begin
        int = this_interval[jint]
        if exist(this_interval) then out = [out, '', $
          'Interval ' + trim(int) + '  Units are ' + spex_units + ' in ' + ctsorph[photons] + ' space']
        out = [out, 'Column 0: Energy Edge Low (keV)','Column 1: Energy Edge High (keV)']

        out = [out,  $
          'Column ' + trim(indgen(ncols)+2) + ': ' + s.id]
        data = strpad(trim(e[0,*]),npad) + strpad(trim(e[1,*]),npad)
        for i=0,ncols-1 do data = data+strpad(trim(s.yvals[*,i,jint]),npad)
        out = [out, data]
      endfor
      text_output, out, file_text=filename
    endif else begin
      print,'No output file selected.  Aborting.'
    endelse

  endif else begin
    print,'No model data to write.  Aborting.'
  endelse

end

;------------------------------------------------------------------------------

; spec and drm are just flags 0/1 saying which file to preview.
; gets the file name out of the parameters.
pro spex::preview, spec=spec, drm=drm, nomore=nomore, out=out

  error=0
  ;catch,error
  if error then begin
    catch, /cancel
    ;print,'in spex::preview catch handler'
    return
  endif

  spec = keyword_set(spec)
  drm = keyword_set(drm)
  if not (spec or drm) then begin
    spec=1 & drm=1
  endif

  out = ''

  if spec then begin
    (self -> get(/obj,class='spex_data')) -> preview, out=out_data, /nomore
    out = [out, out_data]
  endif

  if drm then begin
    if n_elements(out) gt 1 then out=[out,'','']
    (self -> get(/obj,class='spex_drm')) -> preview, out=out_drm, /nomore
    out = [out, out_drm]
  endif

  if not keyword_set(nomore) and n_elements(out) gt 1 then xmessage, out  ;prstr, strjustify(out,/box)

end

;------------------------------------------------------------------------------

; spec and drm are just flags 0/1 saying which file to preview.
; gets the file name out of the parameters.
pro spex::fitsheader, spec=spec, drm=drm, nomore=nomore, out=out

  error=0
  ;catch,error
  if error then begin
    catch, /cancel
    ;print,'in spex::fitsheader catch handler'
    return
  endif

  spec = keyword_set(spec)
  drm = keyword_set(drm)
  if not (spec or drm) then begin
    spec=1 & drm=1
  endif

  out = ''

  if spec then begin
    specfile = self -> get(/spex_specfile)
    if is_string(specfile) then begin
      if is_fits(specfile) then data = mrdfits(specfile, 0, out1, /silent) else $
        out1 = 'Not a FITS file.  No header information.'
      out = [out, 'Spectrum or Image File: ' + specfile, '', out1, '', '']
    endif else out = [out, 'No spectrum or image file selected.']
  endif

  if drm then begin
    srmfile = self -> get(/spex_drmfile)
    if is_string(srmfile) then begin
      if is_fits(srmfile) then data = mrdfits(srmfile, 0, out2, /silent) else $
        out2 = 'Not a FITS file.  No header information.'
      out = [out, 'SRM File: ' + srmfile, '', out2]
    endif else out = [out, 'No srm file selected.']
  endif

  if not keyword_set(nomore) and n_elements(out) gt 1 then xmessage, out ; prstr, strjustify(out,/box)

end

;------------------------------------------------------------------------------

pro spex::setupsummary, short=short, out=out, _extra=_extra

  enab = ['Disabled','Enabled']
  filetime = self -> get(/spex_file_time)
  file_units = self -> get(/spex_data_origunits)

  out = ['Current OSPEX setup     ' + strmid (anytim(!stime, /trunc, /vms), 0, 17), $
    '' ]

  self->preview, /nomore, out=out_preview
  out = [out, out_preview]

  albedo_correct = self->get(/spex_albedo_correct)
  if albedo_correct then begin
    anis = self -> get(/spex_anisotropy)
    angle = self -> get(/spex_source_angle)
    xy = self -> get(/spex_source_xy)
    al_out = ['Albedo Correction: Enabled     Anistropy: ' + trim(anis), $
      'Source Angle: ' + trim(angle) + '    Source Position (X,Y): ' + arr2str(trim(xy, '(f8.2)')) ]
  endif else al_out = 'Albedo Correction: Disabled'

  out = [out, '', al_out]


  bk_sep = self -> get(/spex_bk_sep)
  ;bk_ratio = self -> get(/spex_bk_ratio)
  ;bk_sm_width = bk_ratio ? '   Profile Half Smoothing Width (#points): ' + trim(self -> get(/spex_bk_sm_width)) : ''
  bk_sm_width = '   Profile Half Smoothing Width (#points): ' + trim(self -> get(/spex_bk_sm_width))
  bk_eband = self -> get(/spex_bk_eband)
  a_bands = format_intervals(bk_eband) + ' keV'
  neband = bk_sep ? n_elements(bk_eband)/2 : 1
  a_neband = bk_sep ? '   Number of bands: ' + trim(neband) : ''
  bk_out = ['Separate Background in Energy Bands: ' + enab[bk_sep] + a_neband, bk_sm_width]
  ;if bk_sep then bk_out = [bk_out, 'Use Ratio to High Band: ' + enab[bk_ratio] + bk_sm_width]

  for i=0,neband-1 do begin
    bk_time = self -> get(this_band=i, /this_time)
    bk_order = self -> get(this_band=i, /this_order)
    bk_out = [bk_out, $
      'Background Energy Band: ' + a_bands[i] + '  Method: ' + trim(bk_order) ]
    if bk_time[0] ne -1 then begin
      a_times = format_intervals(bk_time, /ut)
      for j=0,n_elements(a_times)-1 do bk_out=[bk_out, '     ' + a_times[j] ]
    endif else bk_out=[bk_out,'     No times defined for this band.']
  endfor
  bk_out = [bk_out, 'Note: Method = 0/1/2/3/4/5/6 means 0Poly,1Poly,2Poly,3Poly,Exp,High E Profile,This E Profile']


  out = [out, '', bk_out]

  fit_int = self -> get(/spex_fit_time_interval)
  filter = self -> get(/spex_fitint_filter)
  if fit_int[0] eq -1 then fit_out = 'No fit time intervals defined' else begin
    fit_out = 'Number of fit time intervals: ' + trim(n_elements(fit_int[0,*]))
    if not keyword_set(short) then $
      fit_out = [fit_out, format_intervals(fit_int, /ut, /prefix) + '   Filter: ' + trim(filter)]
  endelse
  out = [out, '', fit_out]

  fit_function = self -> get(/fit_function)
  loop_mode = (['Automatic', 'Manual on First Interval', 'Manual on All Intervals'])[self->get(/spex_fit_manual)]
  out = [out, $
    '', $
    'Start method: ' + self->get(/spex_fit_start_method) + $
    '  Loop Mode through intervals: ' + loop_mode, $
    '', $
    'Fit Function: ' + fit_function]

  self->list_function, /simple, /nolist, out=list_func
  out = [out, list_func]

  erange = self -> get(/spex_erange)
  a_erange = same_data(erange,[0.,0.]) ? 'All' : format_intervals(erange, format='(f9.2)')
  if a_erange[0] eq 'None' then a_erange = 'All'
  int_tofit = self -> get(/spex_intervals_tofit)
  if int_tofit[0] eq -1 then a_int_tofit = 'All' else begin
    d1 = find_contig(int_tofit, d2, ind)
    if size(ind,/n_dim) gt 1 then ind = tranpose(ind)
    a_int_tofit = format_intervals(int_tofit[ind])
  endelse

  out = [out, $
    '', $
    'Energy range to fit over: ' + arr2str(a_erange) + ' keV', $
    'Intervals to fit: ' + a_int_tofit ]

  text_output, out, title='Setup Summary', _extra=_extra

end

;--------------------------------------------------------------------------

pro spex::dofit, _extra=_extra

  fit = self -> getdata(class_name='spex_fit', _extra=_extra, /force)

end

;--------------------------------------------------------------------------

; Function to calculate quantities from the spex_summ data.  Currently available values for item:
;
;'data_count_rate' - obs-back data in count rate
;'data_count_flux'  - obs-back data in count flux
;'data_photon_rate' - obs-back data in photon rate
;'data_photon_flux' - obs-back data in photon flux  (default)
;'model_count_rate' - model in count rate
;'model_count_flux' - model in count flux
;'model_photon_rate' - model in photon rate
;'model_photon_flux' - model in photon flux
;'background_count_rate' - background count rate
;'background_count_flux' - background count flux
;'background_photon_rate' - background photon rate
;'background_photon_flux' - background photon flux
;'thick_energy' - non-thermal Electron Energy Flux
;'thick_integ_flux' - integrated electron flux (for thick2_vnorm or thick2_rc - same as a[0] for thick)
;'flux_at_e' - photon flux time profile at E keV (for model)
;'ct_flux_at_e' - count flux time profile at E keV (from model)
;'dct_flux_at_e' - count flux time profile at E keV (from data)
;'integ_photon_flux' - integrated photon flux over energies in func_energy (for model)
;'chisq_full' - full chi-squared
;'raw_resid' - raw residuals
;'enrange' - energy range
;'edf' - electron distribution function (electron spectrum)
;'therm_energy_plasma' - thermal energy of plasma

; enrange - nth energy range limit (n selected by index keyword) for all intervals
;
; index is index for items that require an index (e.g. enrange)
; func_energy - scalar or array of energies to calculate flux_at_e at
;
; strat - strategy name of function component for items that need one (strategy is the component name followed by #n where
;  n identifies which component when >1 of a particular component, usually 0 (e.g. line#0, but if func is 'line+line', then
;  the two strategies are 'line#0' and 'line#1')
;
; dim1_id contains output labels for multiple traces

function spex::calc_summ, summ=summ, $
  item=item_in, $
  this_interval=this_interval, $
  errors=errors, $
  index=index, $
  func_energy=func_energy, $
  strat=strat, $
  dim1_id=dim1_id, $
  err_msg=err_msg

  err_msg = ''

  checkvar, summ, self -> get(/spex_summ)
  checkvar, item_in, 'ph_flux'
  checkvar, this_interval, 0
  checkvar, index, -1

  ; if spex_summ_time_interval is still a pointer, then nothing stored
  if ptr_chk(summ.spex_summ_time_interval) then begin
    err_msg = 'There are no fit results stored yet.  Aborting'
    message, /cont, err_msg
    return, -1
  endif

  ; make old item names still work by translating to new names
  old = ['ct_flux',         'ph_flux',          'bk_ct_flux',            'bk_ph_flux',             'ct_model']
  new = ['data_count_flux', 'data_photon_flux', 'background_count_flux', 'background_photon_flux', 'model_count_flux']
  q = where (item_in eq old, c)
  if c gt 0 then item = new[q[0]] else item = item_in

  item = strlowcase(item)
  if stregex(item, '^model|^data|^background',/bool) then begin
    tmp = str2arr(item, '_')
    item = tmp[0]
    photons = tmp[1] eq 'photon'
    units = tmp[2]
  endif

  done = summ.spex_summ_fit_done
  ntimes = n_elements(done)

  ninterval = n_elements(this_interval)

  ewidth = get_edge_products(summ.spex_summ_energy, /width)
  nen = n_elements(ewidth)
  if ninterval gt 1 then ewidth = rebin (ewidth, nen, ninterval)

  ; get time and energy bin widths, expand array of widths to n,intervals so we can multiply by data
  twidth = get_edge_products(summ.spex_summ_time_interval[*,this_interval], /width)
  twidth = rebin (twidth,  ninterval, nen)

  area = summ.spex_summ_area
  conv = summ.spex_summ_conv[*,this_interval]

  tags = strlowcase(tag_names(summ))
  tags = ssw_strsplit(tags, 'spex_summ_', /tail) ; remove spex_summ_ from all tag names

  ; first check if item is just a field in spex_summ structure to return
  if is_member(item, tags) then begin
    q = (where (tags eq item) )[0]
    ; field is an array dimensioned [nen, ntimes], then extract data for this_interval
    ; otherwise it's across all intervals
    if same_data(size(summ.(q),/dim), long([nen, ntimes])) then begin
      val = (summ.(q))[*,this_interval]
    endif else begin
      val = reform( index eq -1 ?  summ.(q) : (summ.(q))[index,*] )
      if item eq 'params' then errors = reform(summ.spex_summ_sigmas[index,*])
    endelse
    goto, getout

  endif

  case item of

    'model': begin
      checkvar, units, 'flux'
      val = summ.spex_summ_ph_model[*,this_interval]
      if not photons then val = val * conv
      case units of
        ;          'counts': val = val * twidth * area * ewidth
        'rate': val = val * area * ewidth
        'flux':
        else: err_msg = 'Invalid units choice = ' + units + ' for item ' + item
      end
    end

    'data': begin
      checkvar, units, 'rate'
      val = summ.spex_summ_ct_rate[*,this_interval]
      errors = summ.spex_summ_ct_error[*,this_interval]
      if photons then begin
        val = f_div(val, conv)
        errors = f_div(errors, conv)
      endif
      case units of
        ;          'counts': begin
        ;             val = val * twidth
        ;             errors = errors * twidth
        ;             end
        'rate':
        'flux': begin
          val = val / area / ewidth
          errors = errors / area / ewidth
        end
        else: err_msg = 'Invalid units choice = ' + units + ' for item ' + item
      end
    end

    'background': begin
      checkvar, units, 'rate'
      val = summ.spex_summ_bk_rate[*,this_interval]
      errors = summ.spex_summ_bk_error[*,this_interval]
      if photons then begin
        val = f_div(val, conv)
        errors = f_div(errors, conv)
      endif
      case units of
        ;          'counts': begin
        ;             val = val * twidth
        ;             errors = errors * twidth
        ;             end
        'rate':
        'flux': begin
          val = val / area / ewidth
          errors = errors / area / ewidth
        end
        else: err_msg = 'Invalid units choice = ' + units + ' for item ' + item
      end
    end
    ;
    ;    'ct_flux': begin
    ;       val = summ.spex_summ_ct_rate[*,this_interval] / area / ewidth
    ;       errors = summ.spex_summ_ct_error[*,this_interval] / area / ewidth
    ;       end
    ;    'ph_flux': begin
    ;       val = f_div (summ.spex_summ_ct_rate[*,this_interval] / area / ewidth, conv)
    ;       errors = f_div (summ.spex_summ_ct_error[*,this_interval] / area / ewidth, conv)
    ;       end
    ;    'bk_ct_flux': begin
    ;       val = summ.spex_summ_bk_rate[*,this_interval] / area / ewidth
    ;       errors = summ.spex_summ_bk_error[*,this_interval] / area / ewidth
    ;       end
    ;    'bk_ph_flux': begin
    ;       val = f_div (summ.spex_summ_bk_rate[*,this_interval] / area / ewidth, conv)
    ;       errors = f_div (summ.spex_summ_bk_error[*,this_interval] / area / ewidth, conv)
    ;       end
    ;
    ;    'ct_model': val = summ.spex_summ_ph_model[*,this_interval] * conv

    'thick_energy': begin   ; this is for thick or thick2 or thick2_vnorm or thick2_rc

      if exist(strat) then begin
        ind = fit_function_query(summ.spex_summ_fit_function, comp=strat, /param_index)
        comp = ssw_strsplit(strat+'#', '#')
        status = check_func_electron(comp, vnorm=vnorm, athin=athin)
        if ~status || athin then ind = -1
      endif else begin
        comp_thick=['thick', 'thick2', 'thick2_vnorm', 'thick2_rc', 'photon_thick', 'thick_nui']
        ind=-1
        for i = 0,n_elements(comp_thick)-1 do begin
          ind = fit_function_query(summ.spex_summ_fit_function, comp=comp_thick[i], /param_index)
          if ind[0] ne -1 then begin
            comp = comp_thick[i]
            break
          endif
        endfor
      endelse

      if ind[0] eq -1 then begin
        val = -1
        message,/cont,'No specified component (or no thick component at all ) in fit function.  Not calculating thick_energy.'
      endif else begin
        tp = summ.spex_summ_params[ind[0]:ind[1],*]
        val = fltarr(ntimes)
        for i=0,ntimes -1 do begin
          if done[i] then val[i] = calc_nontherm_electron_energy_flux(tp[*,i], func=comp)
          ;           if done[i] then begin
          ;             if comp eq 'thick2_vnorm' or comp eq 'thick2_rc' then begin
          ;               ; calculate integrated electron flux from flux at ref energy for these two functions
          ;               tp[0,i] = Fint_from_Fref(tp[0,i],tp[1,i],tp[2,i],tp[3,i],tp[4,i],tp[5,i],tp[6,i]), func=comp)
          ;             endif
          ;             val[i] = brmthick_power(tp[0,i],tp[1,i],tp[2,i],tp[3,i],tp[4,i],tp[5,i],eps=1.e-6)
          ;           endif
        endfor
      endelse
    end

    'thick_integ_flux': begin   ; this is for thick2_vnorm or thick2_rc (corresponds to value of a[0] for thick2)

      if exist(strat) then begin
        good_func = stregex(strat,'thick2_vnorm.*|thick2_rc.*', /bool)
        ind = good_func ? fit_function_query(summ.spex_summ_fit_function, comp=strat, /param_index) : -1
        comp = ssw_strsplit(strat+'#', '#')
      endif else begin
        comp = 'thick2_vnorm'
        ind = fit_function_query(summ.spex_summ_fit_function, comp=comp, /param_index)
        if ind[0] eq -1 then begin
          comp = 'thick2_rc'
          ind = fit_function_query(summ.spex_summ_fit_function, comp=comp, /param_index)
        endif
      endelse

      if ind[0] ne -1 then begin
        tp = summ.spex_summ_params[ind[0]:ind[1],*]
        val = fltarr(ntimes)
        for i=0,ntimes -1 do $
          if done[i] then val[i] = Fint_from_Fref(tp[0,i],tp[1,i],tp[2,i],tp[3,i],tp[4,i],tp[5,i],tp[6,i], func=comp)
      endif else err_msg = 'No specified component (or no thick2_vnorm or thick2_rc component at all). Aborting.'
    end

    'flux_at_e': begin
      checkvar, func_energy, 50.
      checkvar, strat, ''

      for ie = 0,n_elements(func_energy)-1 do begin
        if not is_number(func_energy[ie]) or func_energy[ie] lt 0 then begin
          message, 'Invalid energy : ' + trim(func_energy[ie]) + '. Changing to 50.',/cont
          func_energy[ie] = 50.
        endif
      endfor
      en = func_energy
      ; we'll convert energies to 2xn, but some functions interpret one 2 element array as 2 mean
      ; energies instead of one lo/hi bin.  So if only one energy, append another and then get rid
      ; of it afterward
      if n_elements(func_energy) eq 1 then begin
        one_energy = 1
        en = en + [0,20.] ; append second energy 20 keV higher - arbitrary. will delete later
      endif else one_energy = 0

      nen = n_elements(en)
      val = reform(fltarr(ntimes,nen))

      en = transpose( [[en-.01], [en+.01]] ) ; need 2xn, make .02 keV wide

      p = summ.spex_summ_params
      obj = obj_new('fit_function')
      obj -> set, fit_function=summ.spex_summ_fit_function, fit_xvals = en
      for i=0,ntimes-1 do begin
        obj -> set, fit_comp_params = p[*,i]
        if done[i] then begin
          ; if single time, may still have multiple energies in en so val will be [nen]
          if ntimes eq 1 then val[*] = obj -> getdata(strategy_name = strat) else $
            val[i,*] = obj -> getdata(strategy_name = strat)
        endif
      endfor
      if one_energy then val = val[*,0]

      ;           ind = fit_function_query(summ.spex_summ_fit_function, comp=comp, /param_index)
      ;           if ind[0] ne -1 then begin
      ;              p = summ.spex_summ_params[ind[0]:ind[1],*]
      ;           	  for i=-0,ntimes-1 do $
      ;              	if done[i] then val[i,*] = call_function('f_'+comp, en, p[*,i])
      ;           endif else begin
      ;              val = -1
      ;              message,'Invalid function component name ' + comp + ' Aborting.', /cont
      ;           endelse

      dim1_id = trim(func_energy,'(f5.1)') + ' keV'
      return, val
    end

    'ct_flux_at_e': begin
      checkvar, func_energy, 50.
      q = where (func_energy gt min(summ.spex_summ_energy) and func_energy le max(summ.spex_summ_energy), count)
      if count eq 0 then begin
        err_msg = 'Your requested energies are outside the valid range of energies.'
      endif else begin
        func_energy = func_energy[q]
        ind = value_locate(summ.spex_summ_energy[0,*], func_energy)
        ct_model_flux = summ.spex_summ_ph_model * summ.spex_summ_conv
        dim1_id = trim(func_energy,'(f5.1)') + ' keV'
        val = reform(transpose(ct_model_flux[ind,*]))
      endelse
    end

    'dct_flux_at_e': begin
      checkvar, func_energy, 50.
      units = 'flux'
      val = summ.spex_summ_ct_rate
      errors = summ.spex_summ_ct_error
      nt = n_elements(val[0,*])
      ew = ewidth # (intarr(nt)+1)
      val = val / area / ew
      errors = errors / area / ew
      q = where (func_energy gt min(summ.spex_summ_energy) and func_energy le max(summ.spex_summ_energy), count)
      if count eq 0 then begin
        err_msg = 'Your requested energies are outside the valid range of energies.'
      endif else begin
        func_energy = func_energy[q]
        ind = value_locate(summ.spex_summ_energy[0,*], func_energy)
        dim1_id = trim(func_energy,'(f5.1)') + ' keV'
        val = reform(transpose(val[ind,*]))
        errors = reform(transpose(errors[ind,*]))
      endelse
    end

    'integ_photon_flux': begin
      checkvar, func_energy, [10.,50.]
      f_e = float(func_energy)  ; protect input arg
      if n_elements(f_e) eq 1 then f_e = [f_e[0], 50.]
      f_e = f_e(sort(f_e[0:1])) ; in case there are more than 2, and sort if out of order
      ee = summ.spex_summ_energy
      ii = where(ee[0,*] ge f_e[0] and ee[1,*] le f_e[1]) > 0
      ew = get_edges(ee, /width)
      val = total(summ.spex_summ_ph_model[ii,*] * (ew[ii]#(intarr(ntimes)+1)), 1)
      eused = minmax(ee[*,ii])
      dim1_id = 'Integral from ' + arr2str(trim(eused,'(f5.1)'), ' to ') + ' keV'
    end

    'chisq_full': begin
      val = fltarr(ntimes)
      for i = 0,ntimes-1 do begin
        if done[i] and summ.spex_summ_chisq[i] ne 0.  then begin
          q = where (summ.spex_summ_emask[*,i], ndata )
          q = where (summ.spex_summ_free_mask[*,i], nfree )
          val[i] = summ.spex_summ_chisq[i] * (ndata - nfree -1)
        endif
      endfor
    end

    'raw_resid': val = summ.spex_summ_resid[*,this_interval] * summ.spex_summ_ct_error[*,this_interval]

    'enrange': begin
      emask = summ.spex_summ_emask
      ninterval = n_elements(emask[0,*])
      val = fltarr(ninterval)
      for i=0,ninterval-1 do begin
        q = where (emask[*,i], count)
        if count gt 0 then begin
          erange = find_contig_ranges( summ.spex_summ_energy[*,q], epsilon=1.e-5 )
          ; treat erange as 1-d array, and get 'index' element
          val[i] = index lt n_elements(erange) ? erange[index] : 0.
        endif
      endfor
    end

    'edf': begin
      func = summ.spex_summ_fit_function
      comp_arr = str2arr(func, '+')
      q = where(strpos(comp_arr, 'thick') ne -1  or  strpos(comp_arr, 'thin') ne -1, count)
      if count eq 0 then begin
        err_msg = 'Your function has no component to plot the electron distribution for. Aborting.'
      endif else begin
        if count gt 1 then $
          message, /cont, 'You have >1 component to plot edf for. Using first one.'
        fcomp = comp_arr[q[0]]
        if check_func_electron(fcomp) then begin
          pind = fit_function_query(func, comp=fcomp, /param_elem)
          par = summ.spex_summ_params[*,this_interval]
          val = fltarr(nen, ninterval)
          for i=0,ninterval-1 do begin
            par = summ.spex_summ_params[pind,this_interval[i]]
            val[*,i] = calc_electron_distribution(energy=summ.spex_summ_energy, par=par, func=fcomp)
          endfor
        endif else err_msg = 'Can not plot electron distribution for ' + fcomp + '. Aborting'
      endelse
    end

    'therm_energy_plasma': begin
      if exist(strat) then begin
        good_func = stregex(strat,'vth.*|vth_abun*.*', /bool)
        ind = good_func ? fit_function_query(summ.spex_summ_fit_function, comp=strat, /param_index) : -1
      endif else begin
        ind = fit_function_query(summ.spex_summ_fit_function, comp='vth', /param_index)
        if ind[0] eq -1 then ind = fit_function_query(summ.spex_summ_fit_function, comp='vth_abun', /param_index)
      endelse
      if ind[0] ne -1 then begin
        p = summ.spex_summ_params[ind[0]:ind[1],*]
        val = reform(therm_energy_plasma(p[0,*], p[1,*], /t_in_kev))
      endif else err_msg = 'Can only calculate thick_integ_flux for thick2_vnorm or thick2_rc component. Aborting.'
    end

    else: err_msg = 'Unrecognized calc_summ item choice:' + item

  endcase

  if err_msg ne '' then begin
    message, /cont, err_msg
    return, -1
  endif

  getout:
  return, n_elements(val) eq 1 ? val[0] : val
end

;--------------------------------------------------------------------------
; Return tables of time, energy, and map items available for plotting from summary data
; item should be tag name (without spex_summ_) if it's a field in summ structure, otherwise it
; should be one of the allowed items in calc_summ.
; Returns 4 structures:
;  vtime is for items that are a function of time
;  ven is for items that are a function of energy
;  vmap is for items that will be drawn as a map
;  vall contains a summary of the other three -  a list of item names, and a type field indicating
;    which list the item is in: 0,1,2 = time,energy,map item
;
; For each structure:
;   item - name of item to plot
;   label - description of item.  Used for GUI selection, and plot label
; Additional fields for vtime structure:
;   comp - name of function component
;   strat - name of function strategy (even if duplicate components, strat is unique, e.g. line#0)
;   index - index into full array of params for this param (i.e. if function is vth+bpow, the first
;     bpow parameter would have index 3, since 0,1,2 are for vth)
;   plot_free - flag indicating this item should have a bar indicating free/fixed on plot
; Additional fields for ven structure:
;   plot_func - flag indicating whether overlay of function on this item makes sense
;


pro spex::plot_summ_tables, vtime=vtime, ven=ven, vmap=vmap, vall=vall

  summ = self -> get(/spex_summ)
  func_comps = str2arr(summ.spex_summ_fit_function,'+')
  fitobj = obj_new('fit_function',fit_function=summ.spex_summ_fit_function)
  func_strats = fitobj -> get(/compman_strategy)
  obj_destroy, fitobj

  ; First parameter for all function components except these should be plotted on log scale
  comp_nolog = ['albedo', 'drm_mod', 'gain_mod', 'pileup_mod']

  tmp = {item: '',  label: '', comp: '', strat: '', index: -1, plot_free: 0, plot_log: 0}

  nparams = n_elements(summ.spex_summ_params[*,0])
  vtime = replicate(tmp, nparams+40)
  ind = 0
  param_ind = 0
  for i = 0,n_elements(func_comps)-1 do begin
    str = fit_model_components(func_comps[i], /struct)
    for j=0,str.nparams-1 do begin
      vtime[ind].item = 'params'
      ;    vtime[ind].uvalue = 'a'+trim(param_ind)
      vtime[ind].label = 'a['+trim(j)+'] ' + str.param_desc[j]
      vtime[ind].comp = func_comps[i]
      vtime[ind].strat = func_strats[i]
      vtime[ind].index = param_ind
      param_ind = param_ind + 1
      vtime[ind].plot_free = 1
      ; 0'th param should be plotted on log scale (except for comps in comp_nolog array)
      if j eq 0 and ~is_member(func_comps[i], comp_nolog) then vtime[ind].plot_log = 1
      ind = ind + 1
    endfor
    if is_member(func_comps[i], ['thick', 'thick2', 'thick2_vnorm', 'thick2_rc']) then begin
      vtime[ind].item = 'thick_energy'
      vtime[ind].label = 'Non-thermal Electron Energy Flux, erg s^(-1)'
      vtime[ind].comp = func_comps[i]
      vtime[ind].strat = func_strats[i]
      ind = ind + 1
    endif
    if is_member(func_comps[i], ['thick2_vnorm', 'thick2_rc']) then begin
      vtime[ind].item = 'thick_integ_flux'
      vtime[ind].label = 'Integrated Electron Flux, 10^35 electron s^(-1)'
      vtime[ind].comp = func_comps[i]
      vtime[ind].strat = func_strats[i]
      ind = ind + 1
    endif
    if is_member(func_comps[i], ['vth','vth_abun']) then begin
      vtime[ind].item = 'therm_energy_plasma'
      vtime[ind].label = 'Thermal Energy of Plasma, ergs'
      vtime[ind].comp = func_comps[i]
      vtime[ind].strat = func_strats[i]
      ind = ind + 1
    endif
    vtime[ind].item = 'flux_at_e'
    vtime[ind].label = 'Photon Flux at E keV, photons s^(-1) cm^(-2) keV^(-1)'
    vtime[ind].comp = func_comps[i]
    vtime[ind].strat = func_strats[i]
    vtime[ind].plot_log = 1
    ind = ind + 1
  endfor
  vtime[ind].item = 'chisq'
  vtime[ind].label = 'Reduced Chi-Squared'
  ind = ind + 1
  vtime[ind].item = 'chisq_full'
  vtime[ind].label = 'Full Chi-Squared'
  ind = ind + 1
  vtime[ind].item = 'flux_at_e'
  vtime[ind].label = 'Photon Flux at E keV, photons s^(-1) cm^(-2) keV^(-1)'
  vtime[ind].plot_log = 1
  ind = ind + 1
  vtime[ind].item = 'ct_flux_at_e'
  vtime[ind].label = 'Model Count Flux at E keV, counts s^(-1) cm^(-2) keV^(-1)'
  vtime[ind].plot_log = 1
  ind = ind + 1
  vtime[ind].item = 'dct_flux_at_e'
  vtime[ind].label = 'Data Count Flux at E keV, counts s^(-1) cm^(-2) keV^(-1)'
  vtime[ind].plot_log = 1
  ind = ind + 1
  vtime[ind].item = 'integ_photon_flux'
  vtime[ind].label = 'Integrated Photon Flux, photons s^(-1) cm^(-2)'
  vtime[ind].plot_log = 1
  vtime = vtime[0:ind]
  nt = ind+1


  nm = 2
  vmap = replicate({vmap, item: '', label: ''}, nm)
  vmap[0] = [{vmap, 'ff_map', 'Free/Fixed Parameter Map'}, $
    {vmap, 'enrange_map', 'Energy Range Map'}]

  nen=16
  ven = replicate( {ven, item: '', label: '',  plot_func: 0, plot_bk: ''}, nen)
  ven[0] = [ {ven, 'model_photon_flux', 'Model Photon Flux', 0, 'background_photon_flux'}, $
    {ven, 'model_photon_rate', 'Model Photon Rate', 0, 'background_photon_rate'}, $
    {ven, 'model_count_flux', 'Model Count Flux', 0, 'background_count_flux'}, $
    {ven, 'model_count_rate', 'Model Count Rate', 0, 'background_count_rate'}, $
    {ven, 'data_photon_flux', 'Data Photon Flux', 1, 'background_photon_flux'}, $
    {ven, 'data_photon_rate', 'Data Photon Rate', 1, 'background_photon_rate'}, $
    {ven, 'data_count_flux', 'Data Count Flux', 1, 'background_count_flux'}, $
    {ven, 'data_count_rate', 'Data Count Rate', 1, 'background_count_rate'}, $
    {ven, 'background_photon_flux', 'Background Photon Flux', 0, ''}, $
    {ven, 'background_photon_rate', 'Background Photon Rate', 0, ''}, $
    {ven, 'background_count_flux', 'Background Count Flux', 0, ''}, $
    {ven, 'background_count_rate', 'Background Count Rate', 0, ''}, $
    {ven, 'conv', 'Conversion Factors (counts/photon)', 0, ''}, $
    {ven, 'resid', 'Normalized Residuals (sigma)', 0, ''}, $
    {ven, 'raw_resid', 'Raw Residuals (counts/sec)', 0, ''}, $
    {ven, 'edf', 'Electron Spectrum', 0, ''} ]

  vall = {item: [vtime.item,ven.item,vmap.item], type: [intarr(nt), 1+intarr(nen),2+intarr(nm)]}
  s = get_uniq(vall.item, sindex)
  vall = {item: vall.item[sindex], type: vall.type[sindex]}

end

;--------------------------------------------------------------------------

; plot_summ plots items from the spex_summ structure of stored fit results
; xplot can be 'time', 'interval' or 'energy', or any of the items in plot_summ_tables structure
; yplot = If xplot is time or interval, any of the items in plot_summ_tables vtime or vmap structure.
;         If xplot is energy, any of the items in plot_summ_tables vtime structure.
;         If xplot isn't one of those three, then xplot and yplot must both be items in same structure (vtime or ven)
;
; xindex, yindex  - index into full spex_summ_params array for parameter to plot, only apply when
;   xplot, yplot is an item requiring an index (e.g. params)
; xstrat, ystrat - function component strategy to plot (e.g. if 2 lines, first line is line#0, 2nd is
;   line#1), only applies for xplot,yplot for a function component
; xlabel, ylabel - plot labels, if passed, will override what plot_summ sets
; show_error - if set, overlay errors on plot where applicable
; show_free - if set, show free/fixed bar on time plots of parameters
; show_erange - if set, mark range of fitted energies where applicable
; overlay_model - if set, overlay model where applicable
; this_interval - scalar or array of interval numbers to plot (for energy plots)
; psym, xlog, ylog - if passed, will override what plot_summ sets.
; no_plotman - If set, plot in regular window, not plotman window
; Can also pass in most plot keywords in _extra

pro spex::plot_summ, $
  xplot=xplot, xindex=xindex, xstrat=xstrat, $
  yplot=yplot, yindex=yindex, ystrat=ystrat, $
  xlabel=xlabel, ylabel=ylabel, $
  show_errors=show_errors, $
  show_erange=show_erange, $
  show_free=show_free, $
  overlay_model=overlay_model, $
  overlay_bk=overlay_bk, $
  this_interval=this_interval, $
  psym=psym, xlog=xlog, ylog=ylog, $
  no_plotman=no_plotman, $
  err_msg=err_msg, $
  _extra=_extra

  if not self->valid_summ_params() then begin
    err_msg = 'There are no fit results stored yet.  Aborting.'
    goto, error_exit
  endif

  self -> plot_summ_tables, vtime=vtime, ven=ven, vmap=vmap, vall=vall

  use_plotman = (keyword_set(no_plotman) eq 0)

  checkvar, xplot,'time'
  checkvar, xindex, -1
  checkvar, xstrat, ''
  checkvar, yplot, 'chisq'
  checkvar, yindex, -1
  checkvar, ystrat, ''
  checkvar, overlay_model, 0
  checkvar, overlay_bk, 0

  checkvar, this_interval, self->get(/spex_interval_index)

  err_msg = ''

  xplot = strlowcase(xplot)
  yplot = strlowcase(yplot)

  summ = self -> get(/spex_summ)

  overlay_obj = -1

  bad = where (summ.spex_summ_fit_done eq 0., nbad)  ; for time plots, insert NaN's if interval wasn't fit

  if xplot eq 'time' or xplot eq 'interval' then begin

    ; if requested a map, call plot_summ_map and return
    yitem = where (vmap.item eq yplot, count)
    if count gt 0 then begin
      self -> plot_summ_map, xplot, vmap[yitem[0]].item, summ, err_msg, no_plotman=no_plotman, _extra=_extra
      return
    endif

    if xplot eq 'interval' then begin
      nx = n_elements(summ.spex_summ_fit_done)
      xdata = findgen(nx) ; use findgen, not indgen, or NaNs will be 0s
      if nbad gt 0 then xdata[bad] = !values.f_nan
      plot_type = 'xyplot'
      checkvar, xlabel, 'Interval #'
      checkvar,psym,-2
    endif else begin
      xdata = summ.spex_summ_time_interval
      nx = n_elements(xdata)/2  ; since for time intervals, have start,end of each
      utbase = min(xdata)
      xdata = xdata - utbase
      if nbad gt 0 then xdata[*,bad] = !values.f_nan
      plot_type = 'utplot'
      checkvar, xlabel, ''
      checkvar,psym,10
    endelse

    ; prepend blank line so if we plot free bar, there's room
    label = ['', 'Fit Function: ' + summ.spex_summ_fit_function]

    iy = where (vtime.item eq yplot and vtime.index eq yindex and vtime.strat eq ystrat, count)
    if count eq 0 then goto, error_exit

    iy = iy[0]
    checkvar, ylabel, vtime[iy].strat + ' ' + vtime[iy].label
    title = ylabel
    ystrat = vtime[iy].strat

    ydata = self->calc_summ(summ=summ, item=vtime[iy].item, index=yindex, strat=ystrat, error=yerr, _extra=_extra, dim1_id=dim1_id)
    if ydata[0] eq -1 then goto, error_exit
    if nx eq 1 then begin
      print,'Can not plot a single time interval. Values are '
      if xplot eq 'interval' then print, trim(xdata), + '    ' + arr2str(ydata) else $
        print, format_intervals(xdata+utbase, /ut) + '   ' + arr2str(ydata)
      goto,error_exit
    endif

    if nbad gt 0 then ydata[bad] = !values.f_nan

    checkvar,xlog,0
    ;  checkvar,ylog, (strpos(vtime[iy].item,'flux_at_e') ne -1)  ; linear for everything other than flux_at_e plots
    checkvar,ylog, vtime[iy].plot_log

    if keyword_set(show_free) and vtime[iy].plot_free then begin
      ; if plotting bar showing which intervals this param was free for, do it
      ; by setting up the addplot capability in xyplot.  xyplot will call
      ; procedure named addplot_name with arguments in addplot_arg after drawing plot.
      free = summ.spex_summ_free_mask[yindex,*]
      if nbad gt 0 then free[bad] = !values.f_nan
      find_changes, free, index, state, index2d=index2d, count=count
      qf = where(state eq 1, c)
      if c gt 0 then begin
        addplot_name = 'mark_intervals'
        ; sbar, ebar are start, end time or interval when free mask was 1
        if size(xdata,/n_dim) eq 2 then begin
          sbar = xdata[0,index2d[0,qf]]
          ebar = xdata[1,index2d[1,qf]]
        endif else begin
          sbar = xdata[index2d[0,qf]]
          ebar = xdata[index2d[1,qf]]
        endelse
        valid_plotman = 0
        if not keyword_set(no_plotman) then begin
          plotman_obj= self -> get_plotman_obj(valid=valid_plotman)
          ; if using plotman get colors, otherwise will use different linestyles
          if valid_plotman then begin
            cn = plotman_obj -> get(/color_names)
            barcolor = cn.green
          endif
        endif
        if not valid_plotman then begin
          tvlct,rsave,gsave,bsave,/get
          linecolors
          barcolor = 7 ; green
        endif
        addplot_arg = {sbar:sbar, ebar:ebar, label:'Free', bottom:0, utplot:0, color: barcolor}
      endif
    endif

  endif else begin      ; end of x equals time or interval branch

    if xplot eq 'energy' then begin
      iy = where (ven.item eq yplot, count)
      if count eq 0 then goto, error_exit

      checkvar, xlabel, 'Energy (keV)'
      overlay_obj = -1

      xdata = summ.spex_summ_energy
      yname = ven[iy].item
      ; if electron distribution function is selected, and only one interval is selected, then use
      ; plot_edf method to get nice edf plot.  Otherwise do the usual call to calc_summ for all the
      ; intervals requested and plot the edf for all intervals on one plot
      if yname eq 'edf' and n_elements(this_interval) eq 1 then begin
        (self->get(/obj,class='spex_fit'))->plot_edf, interval=this_interval
        return
      endif
      checkvar, ylabel, ven[iy].label
      title=ylabel

      ydata = self->calc_summ(summ=summ, item=yname, this_interval=this_interval, errors=errors)
      if exist(errors) then yerr = errors
      if ydata[0] eq -1 then goto, error_exit

      if overlay_model and ven[iy].plot_func then begin
        units = strpos(yname, 'rate') eq -1 ? 'flux' : 'rate'
        photons = strpos(yname,'ph') eq -1 ? 0 : 1
        overlay_obj = self->make_fitplot_obj ( spex_units=units, $
          photons=photons, /overlay, this_interval=this_interval, /use_fitted)
      endif

      if keyword_set(show_erange) then begin
        erange = self -> get_erange (/use_fitted, this_interval=this_interval)
        if erange[0] ne -1 then begin
          addplot_name = 'oplot_xmarkers'
          addplot_arg = {intervals: erange}
          if n_elements(this_interval) gt 1 then $
            message,'Note: energy intervals for first selected interval are displayed.', /cont
        endif
      endif
      plot_type = 'xyplot'
      dim1_id = trim(this_interval) +'. ' + format_intervals(summ.spex_summ_time_interval[*,this_interval], /ut)
      checkvar,xlog,1
      checkvar,ylog, (strpos(yname,'resid') eq -1)  ; log for everything other than residual plots
      checkvar,psym,10

      if overlay_bk and ~(ven[iy].plot_bk eq '') then begin
        ydata_bk = self->calc_summ(summ=summ, item=ven[iy].plot_bk, this_interval=this_interval, errors=errors)
        if max(ydata_bk) ne 0. then begin
          ydata = [[ydata],[ydata_bk]]
          if exist(yerr) && exist(errors) then yerr = [[yerr], [errors]]
          dim1_id = [dim1_id, dim1_id + ' Bk']
        endif
      endif

    endif else begin   ; end of x equals energy branch

      ; one item vs another branch is below
      qx = where (vall.item eq xplot,countx) & qy = where (vall.item eq yplot, county)
      if countx eq 0 or county eq 0 then goto, error_exit  ; x and y item valid?
      if n_elements(get_uniq(vall.type[qy])) ne 1 then goto, error_exit ; all y's the same type?
      if vall.type[qx] ne vall.type[qy] then goto, error_exit  ; x and y type the same?
      plot_type = 'xyplot'
      checkvar,psym,0
      case vall.type[qx] of
        0: begin  ;time items
          ix = (where (vtime.item eq xplot and vtime.index eq xindex and vtime.strat eq xstrat, count))[0]
          checkvar, xlabel, vtime[ix].strat + ' ' + vtime[ix].label
          xstrat = vtime[ix].strat
          xdata = self->calc_summ(summ=summ, item=vtime[ix].item, index=xindex, strat=xstrat, _extra=_extra)
          iy = (where (vtime.item eq yplot and vtime.index eq yindex and vtime.strat eq ystrat, count))[0]
          checkvar, ylabel, vtime[iy].strat + ' ' + vtime[iy].label
          ystrat = vtime[iy].strat
          title = ylabel
          ydata = self->calc_summ(summ=summ, item=vtime[iy].item, index=yindex, strat=ystrat, error=yerr, _extra=_extra, dim1_id=dim1_id)
          if xdata[0] eq -1 or ydata[0] eq -1 then goto, error_exit
          if nbad gt 0 then xdata[bad] = !values.f_nan
          if nbad gt 0 then ydata[bad] = !values.f_nan
          checkvar, xlog, vtime[ix].plot_log
          checkvar, ylog, vtime[iy].plot_log
        end
        1: begin  ;energy items
          ix = (where (ven.item eq xplot, count))[0]
          checkvar, xlabel, ven[ix].label
          xdata = self->calc_summ(summ=summ, item=xplot, this_interval=this_interval, _extra=_extra)
          iy = (where (ven.item eq yplot, count))[0]
          checkvar, ylabel, ven[iy].label
          title = ylabel
          ydata = self->calc_summ(summ=summ, item=yplot, error=yerr, this_interval=this_interval, _extra=_extra)
          if xdata[0] eq -1 or ydata[0] eq -1 then goto, error_exit
        end
        else: goto, error_exit
      endcase
    endelse
  endelse

  checkvar, addplot_name, ''
  checkvar, addplot_arg, -1
  checkvar, utbase, 0

  ; convert units in text string to units with superscripts for plot labels
  units_pos = strpos(xlabel,'s^(-1) cm^(-2) keV^(-1)')
  if units_pos gt -1 then xlabel = strmid(xlabel,0,units_pos) + 's!u-1!n cm!u-2!n keV!u-1'
  units_pos = strpos(ylabel,'s^(-1) cm^(-2) keV^(-1)')
  if units_pos gt -1 then ylabel = strmid(ylabel,0,units_pos) + 's!u-1!n cm!u-2!n keV!u-1'

  plot_params = {plot_type:plot_type, $
    xdata: xdata, $
    utbase: anytim(utbase,/vms), $
    ydata: ydata , $
    id:'', $
    panel_id: title, $
    data_unit:ylabel, $
    xtitle:xlabel, $
    dim1_enab_sum:1, $
    overlay_obj: overlay_obj, $
    addplot_name: addplot_name, $
    addplot_arg: addplot_arg }
  if exist(dim1_id) then plot_params=add_tag (plot_params, dim1_id, 'dim1_id')
  if keyword_set(show_errors) and exist(yerr) then plot_params=add_tag (plot_params, yerr, 'edata')
  if exist(label) then plot_params=add_tag(plot_params, label, 'label')

  self -> do_plot, use_plotman=use_plotman, plot_params=plot_params, psym=psym, xlog=xlog, ylog=ylog, _extra=_extra
  if exist(rsave) then tvlct,rsave,gsave,bsave
  return


  error_exit:
  prstr,/nomore, 'No data to plot'

end

;--------------------------------------------------------------------------

; Plot summary results in map form.
; xplot must be interval
; yplot must be an item in vmap structure.
; summ - spex_summ structure
; no_plotman - If set, plot in regular window, not plotman window
; Doesn't work for xplot=time yet.  Tried to use specplot obj (see commented out code below)
; but didn't work right.

pro spex::plot_summ_map, xplot, yplot, summ, err_msg, no_plotman=no_plotman, _extra=_extra

  done = summ.spex_summ_fit_done
  qnofit = where (done eq 0, count_nofit)
  time0 = anytim(/vms,summ.spex_summ_time_interval[0])

  if xplot eq 'time' then $
    message,'X axis = time not implemented yet.  Plotting Interval # on X axis.', /cont

  case yplot of
    'ff_map': begin
      free=summ.spex_summ_free_mask
      ff=free + 2b
      if count_nofit gt 0 then ff[*,qnofit] = 1b
      ff = transpose(ff)
      sz = size(ff,/dim)
      obj = obj_new('spex_map')
      obj -> setmap,data=ff,dx=1,dy=1,xc=sz[0]/2.,yc=sz[1]/2.,$
        xunits='Interval Index',yunits='Parameter Index',$
        id='Free/Fixed Mask', time=time0
      obj -> set,legend=['Not fit', 'Fixed', 'Free'], /xgrid, /ygrid
      panel_desc = 'Free/Fixed Mask'
      func = summ.spex_summ_fit_function
      func_arr = str2arr(func,'+')
      ind = fit_function_query(func,/param_index)
      yvals = ind[*,0] + .2
      obj -> set, label_arr = func_arr, xlabpos = -1, ylabpos = yvals
    end

    'enrange_map': begin
      nint = n_elements(summ.spex_summ_emask[0,*])
      en = summ.spex_summ_energy
      enmm = minmax(en)
      ; divide energy into .1 keV bins so they'll be smaller than any actual bins
      envals = enmm[0] + findgen((enmm[1]-enmm[0])/.1) * .1
      enmap = bytarr(n_elements(envals),nint) + 2b
      for int=0,nint-1 do begin
        if done[int] eq 0 then continue
        eindex = where(summ.spex_summ_emask[*,int])
        eranges = find_contig_ranges(en[*,eindex], epsilon=1.e-5)
        ind = value_locate(envals, eranges)
        for j=0,n_elements(eranges[0,*])-1 do enmap[ind[0,j]:ind[1,j], int] = 3b
      endfor
      if count_nofit gt 0 then enmap[*,qnofit] = 1b
      enmap = transpose(enmap)
      obj = obj_new('spex_map')
      obj -> setmap, data=enmap, dx=1,dy=.1, xc=nint/2., yc=average(enmm), $
        xunits='Interval Index', yunits='keV', $
        id='Energy Range', time=time0
      obj -> set,legend=['Not fit', 'Energy not used', 'Energy used'], /xgrid
      panel_desc = 'Energy Range'

      ;    xdata = findgen(n_elements(ff[*,0]))
      ;    obj2 = obj_new('specplot',xdata,ff)
      ;    obj2 -> set,dim1_vals=.5+findgen(n_elements(ff[0,*])),/no_ut,drange=[0.,6.],utbase=0.
      ;    obj->plot
      ;    plotman_obj -> new_panel,input=obj2, plot_type='specplot', cbar=0,  legend_color=0, desc='Free/Fixed Mask'
    end

  endcase

  valid_plotman = 0
  if not keyword_set(no_plotman) then begin
    plotman_obj= self -> get_plotman_obj(valid=valid_plotman)
    if valid_plotman then plotman_obj -> new_panel, input=obj, plot_type='image', $
      square_scale=0, cbar=0, legend_color=0, legend_loc=2, desc=panel_desc, _extra=_extra
  endif
  if not valid_plotman then obj->plot, _extra=_extra
  destroy, obj

  ;times = o->get(/spex_summ_time_interval)
  ;utbase = min(times)
  ;xdata=average(anytim(times),1)-utbase
  ;obj=obj_new('specplot',xdata,ff)
  ;obj->set,dim1_vals=.5+findgen(n_elements(ff[0,*])),utbase=utbase,drange=[0.,6.]
  ;obj->plot
  ;p->new_panel,input=obj, plot_type='specplot', cbar=0,  legend_color=0, desc='Free/Fixed Mask'


end

;--------------------------------------------------------------------------
; function to return the names of control and admin parameters
; if new_input is set, only return those params that have an N after their name in the text file

function spex::get_param_names, new_input=new_input

  filename = 'spex_params_for_script.txt'
  file = concat_dir('SSW_OSPEX', filename)
  file  = findfile(file, count=count)
  if count eq 0 then begin
    message,'Can not find file ' + filename + '.  Aborting.  Is SSW_OSPEX defined?',/cont
    retall
  endif

  ;params = rd_ascii(file, error=error)
  ;if error then return, ''
  params = rd_tfile(file, 2, delim='#', nocomment=';')
  if params[0] eq '' then return, ''

  if keyword_set(new_input) then begin
    q = where (strpos(params[1,*], 'N') ne -1, count)
    if count gt 0 then params = params[0,q]
  endif else params = params[0,*]

  params = strlowcase(reform(params))

  return, params
end

;--------------------------------------------------------------------------

; Procedure to allow use to set parameters manually.
;
pro spex::setParams


  plotman_obj = self->get_plotman_obj(/nocreate,valid=valid)
  if valid then group=plotman_obj->get(/plot_base)

  control = self->get(/control)

  status = xset_struct("OSPEX", group, control, _extra)

  if status then self->set,_extra=_extra
end

;--------------------------------------------------------------------------
; Function is_image_input returns 1 if the input source is an image cube.
; Also returns image object in keyword if requested.
function spex::is_image_input, img_obj=img_obj
  img_obj = (self->get(/obj,class='spex_data')) -> getstrategy()
  return, is_class( img_obj, 'spex_image')
end


;--------------------------------------------------------------------------

; Procedure to write a script to set up an ospex session using the parameter setup
; in the current session of ospex.  Script will be a procedure named whatever the file
; name is, and will have obj=obj as an argument.  If obj is an OSPEX object, it will
; be used, otherwise it will be created first.
;
; Keywords:
; outfile - name of script file to write.  If not set, user will be prompted for name
; restorefit - If set, then save fit results in a file, and include the restorefit
;   command in the script to restore the results from that file
; fit_outfile - if restorefit is set, then this is name of fit results output file.  If
;   not supplied, user will be prompted for name
; gui - if set, command in script that creates ospex object will start with gui
;

pro spex::writescript, outfile=outfile, $
  restorefit=restorefit, $
  fit_outfile=fit_outfile, $
  gui=gui, sav=sav

  params = self -> get_param_names()
  if params[0] eq '' then return

  if size(outfile, /tname) ne 'STRING' then begin
    atim = strlowcase( trim( str_replace(anytim(!stime, /vms, /date), '-', '_') ))
    outfile = 'ospex_script_'+atim+'.pro'
    outfile = dialog_pickfile (path=curdir(), filter='*.pro', $
      file=outfile, $
      title = 'Select output file',  $
      group=group, $
      get_path=path)
  endif

  if outfile ne '' then begin

    if keyword_set(restorefit) then $
      self -> savefit, outfile=fit_outfile, sav=sav else $
      fit_outfile = ''

    ; First get all default parameters for OSPEX by creating a dummy object.
    ; This is really kludgy, but if spex_data strategy is a spex_image class, then we
    ; won't get the spex_image default params unless we set that strategy in the dummy object.

    def_obj = ospex(/no_gui)
    if self->is_image_input() then $
      (def_obj->get(/obj,class='spex_data')) -> setstrategy,'SPEX_HESSI_IMAGE'
    defp = def_obj->get()
    deftags = strlowcase(tag_names(defp))
    obj_destroy, def_obj

    curp = self->get()
    curtags = strlowcase(tag_names(curp))

    gui = keyword_set(gui) ? '()' : '(/no_gui)'
    break_file, outfile, dum,dum, pro_name
    out = ['; OSPEX script created ' + systime(0) + ' by OSPEX writescript method.', $
    ';', $
    ';  Call this script with the keyword argument, obj=obj to return the ', $
    ';  OSPEX object reference for use at the command line as well as in the GUI.', $
    ';  For example: ', $
    ';     ' + pro_name + ', obj=obj', $
    ';', $
    ';  Note that this script simply sets parameters in the OSPEX object as they', $
    ';  were when you wrote the script, and optionally restores fit results.', $
    ';  To make OSPEX do anything in this script, you need to add some action commands.', $
    ';  For instance, the command', $
    ';     obj -> dofit, /all', $
    ';  would tell OSPEX to do fits in all your fit time intervals.', $
    ';  See the OSPEX methods section in the OSPEX documentation at ', $
    ';  http://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_explanation.htm', $
    ';  for a complete list of methods and their arguments.', $
    ';', $
    'pro ' + pro_name + ', obj=obj', $
      "if not is_class(obj,'SPEX',/quiet) then obj = ospex"+gui]

    bk_sep = self->get(/spex_bk_sep)

    for i = 0,n_elements(params)-1 do begin
      qc = where(params[i] eq curtags, countc)
      qd = where (params[i] eq deftags, countd)
      if countc gt 0 and countd gt 0 then begin
        ic = qc[0]
        id = qd[0]
        if not same_data(curp.(ic), defp.(id)) then begin
          val = curp.(ic)

          ; We need absolute times to be anytim ascii format (not seconds) or
          ;   we will lose precision.
          ; This is a kludge because we're assuming that all time params have 'time' in name
          ;   and are double.  If < 1.e6, it's a relative time.
          if size(val, /tname) eq 'DOUBLE' and $
            (strpos(params[i], 'time') ne -1 or strpos(params[i],'tband') ne -1) then begin
            if val[0] gt 1.e6 then val = anytim(val,/vms)
          endif

          ; If separate background bands, the bk times are pointer arrays - have to set per band, but
          ; otherwise, set normally
          if params[i] eq 'spex_bk_time_interval' or $
            params[i] eq 'spex_bk_eband' or $
            params[i] eq 'spex_bk_order' then begin
            bk_eband = self -> get(/spex_bk_eband)
            if not bk_sep then bk_eband = bk_eband[*,0]
            nband = (bk_eband[0] eq -1) ? 1 : n_elements(bk_eband)/2
            case params[i] of
              'spex_bk_time_interval': begin
                for ib=0,nband-1 do begin
                  time = self -> get(this_band=ib, /this_time)
                  if time[0] gt 1.e6 then time = anytim(time, /vms)
                  str = bk_sep ? 'this_band=' + val2string(ib) + ', this_time=' + val2string(time) : $
                    'spex_bk_time_interval=' + val2string(time)
                  out = [out, 'obj-> set, ' + str]
                endfor
              end
              'spex_bk_order': begin
                order = self -> get(/spex_bk_order)
                if not bk_sep then order = order[0]
                out = [out, 'obj-> set, spex_bk_order=' + val2string(order)]
              end
              'spex_bk_eband': if bk_sep then out = [out, 'obj-> set, spex_bk_eband=' + val2string(bk_eband)]
            endcase
          endif else begin
            sval = val2string(val)
            ;print,params[i] + '= ', +sval
            if sval[0] ne 'BAD' then out = [out, 'obj-> set, ' + params[i] + '= ' + sval ]
          endelse
        endif
      endif
    endfor

    roimsg = ''
    if self->is_image_input() then begin
      self->save_roi, spex_roi_outfile, done=done
      if done then begin
        out = [out, "obj -> set, spex_roi_infile='" + spex_roi_outfile + "'" ]
        roimsg = ' Script includes setting ROI input file command.'
      endif else roimsg = ' No ROI save file created; no spex_roi_infile in script.'

    endif

    if fit_outfile eq '' then begin
      fitmsg = '.  No fit results saved; no restorefit in script.'

    endif else begin
      out = [out, "obj -> restorefit, file='" + fit_outfile + "'"]
      fitmsg = ' Script includes restorefit command.'
    endelse

    out = [out, 'end']

    ; Wrap long lines on ', '
    out = wrap_txt(out, delim=', ', length=90)


    wrt_ascii, out, outfile, err_msg=err_msg
    if err_msg eq '' then print,'Wrote commands in script file ' + outfile + fitmsg + roimsg

  endif else print,'No output file name selected. No file written.'

end

;--------------------------------------------------------------------------

pro spex::runscript, file=file, init=init

  if size(file, /tname) ne 'STRING' then begin
    file = dialog_pickfile (path=curdir(), filter=['*.pro'], $
      file=file, $
      title = 'Select OSPEX script file to run',  $
      group=group, $
      get_path=path)
  endif

  if file ne '' then begin
    if keyword_set(init) then begin
      self -> init_params
    endif
    break_file, file, disk, dir, pro_name
    cd, disk+dir, current=savedir
    resolve_routine, pro_name
    cd, savedir
    call_procedure, pro_name, obj=self
    message, 'Ran procedure ' + pro_name + ' to set OSPEX parameters', /info
  endif else message,'No script file selected.  Script not run.', /info

end

;--------------------------------------------------------------------------

pro spex::savefit, outfile=outfile, sav=sav

  ; Check sav keyword to determine the format for the fit results file
  if (keyword_set( sav )) then begin
    self->savefit_sav, outfile=outfile
  endif else begin
    self->fitswrite, outfile=outfile, /fit_results
  endelse

end

;--------------------------------------------------------------------------

pro spex::savefit_sav, outfile=outfile

  message, 'Fit results save output file are no longer supported. Use FITS format.', /cont
  return

  ;fit_results = self -> get(/spex_summ)
  ;
  ;IF (n_elements(fit_results.SPEX_SUMM_TIME_INTERVAL) gt 1) THEN BEGIN
  ;
  ;    if size(outfile, /tname) ne 'STRING' then begin
  ;        atim = trim( str_replace(anytim(!stime, /vms, /date), '-', '_') )
  ;        outfile = 'ospex_results_'+atim+'.geny'
  ;        outfile = dialog_pickfile (path=curdir(), filter='*.geny', $
  ;           file=outfile, $
  ;           title = 'Select output save file name',  $
  ;           group=group, $
  ;           get_path=path)
  ;    endif
  ;
  ;    if outfile ne '' then begin
  ;        savegenx, fit_results, file=outfile, /over
  ;        message, 'Fit results saved in file ' + outfile, /info
  ;    endif else message, 'No output file name selected.  Not saving fit results.', /info
  ;
  ;endif  ELSE BEGIN
  ;    message, 'No Fit Results to save.  Aborting.', /info
  ;    outfile = ' '
  ;ENDELSE

end

;--------------------------------------------------------------------------

; restorefit method reads a FITS or geny file and restores the fits results (all spex_summ...
; parameters) into the object.
pro spex::restorefit, file=file, nodialog=nodialog

  if size(file, /tname) ne 'STRING' then begin
    if not keyword_set(nodialog) then begin
      file = dialog_pickfile (path=curdir(), filter='*.fits', $
        file=file, $
        title = 'Select geny or FITS file to restore fit results from',  $
        group=group, $
        get_path=path)
    endif else file = ''
  endif

  if file ne '' then begin
    fit_results = spex_read_fit_results(file)
    if size(/tname, fit_results) ne 'STRUCT' then begin
      message, 'Aborting restorefit.', /cont
      return
    endif

    ;initialize all spex_summ params to nothing (in case structure from file doesn't have all tags)
    self -> clearsumm

    ; set structure from file into spex_summ parameters
    self -> set, _extra = fit_results

    msg = 'Fit results restored from ' + file + ' and set into spex_summ... params in ospex object.'
    if spex_get_nointeractive() or keyword_set(nodialog) then print,msg else $
      a=dialog_message(msg, /info)
  endif else message,'No input file selected.  No fit results restored.', /info

end

;--------------------------------------------------------------------------
; init_params method resets all params (in spex_params_for_script.txt) to default values
; If new_input keyword is set, reset only those parameters that should be reinitialized
;   when a new input file is selected.
; param_names is an output keyword - returns the names of the parameters changed

pro spex::init_params, new_input=new_input, param_names=param_names

  param_names = self -> get_param_names(new_input=new_input)
  if param_names[0] eq '' then return

  ; First get all default parameters for OSPEX by creating a dummy object.
  ; This is really kludgy, but if spex_data strategy is a spex_image class, then we
  ; won't get the spex_image default params unless we set that strategy in the dummy object.

  temp_obj = ospex(/no_gui)
  if self->is_image_input() then $
    (temp_obj->get(/obj,class='spex_data')) -> setstrategy,'SPEX_HESSI_IMAGE'
  def_params = temp_obj -> get()
  obj_destroy, temp_obj

  clear_params = str_subset(def_params, param_names)

  ; temporarily set to nointeractive so spex_data::set won't prompt for file names
  ; for spex_specfile and spex_drmfile when we set them to blanks (the defaults)
  nointer = spex_get_nointeractive()
  setenv,'OSPEX_NOINTERACTIVE=1'
  self -> set, _extra=clear_params
  setenv,'OSPEX_NOINTERACTIVE=' + trim(nointer)

  if keyword_set(new_input) then begin
    ;message,'New input source.  Reinitializing key parameters to default values.', /info
  endif else begin
    self -> clearsumm
    message, 'All parameters set to default values.', /info
  endelse
end

;--------------------------------------------------------------------------

pro spex::clearsumm
  summ = self -> get(/spex_summ)
  null_summ = {spex_fit_info}
  self -> set, _extra = str_subset(null_summ, tag_names(summ))
end

;--------------------------------------------------------------------------

pro spex::roi, _extra=_extra
  if self ->is_image_input(img_obj=img_obj) then img_obj -> roi, spex_obj=self, _extra=_extra
end

;--------------------------------------------------------------------------

pro spex::roi_config, _extra=_extra
  if self ->is_image_input(img_obj=img_obj) then img_obj -> roi_config, _extra=_extra
end

;--------------------------------------------------------------------------

pro spex::save_roi, spex_roi_outfile, _ref_extra=_extra
  if self ->is_image_input(img_obj=img_obj) then img_obj -> save_roi, spex_roi_outfile, _extra=_extra
end

;--------------------------------------------------------------------------

pro spex::restore_roi, spex_roi_infile, _ref_extra=_extra
  if self ->is_image_input(img_obj=img_obj) then img_obj -> restore_roi, spex_roi_infile, _extra=_extra
end

;--------------------------------------------------------------------------

pro spex::list_roi, spex_roi_file
  if self ->is_image_input(img_obj=img_obj) then img_obj -> list_roi, spex_roi_file
end

;--------------------------------------------------------------------------

pro spex::panel_show, _extra=_extra
  if self ->is_image_input(img_obj=img_obj) then img_obj -> panel_show, _extra=_extra
end

;--------------------------------------------------------------------------

pro spex::list_function, out=out, _extra=_extra
  (self->get(/obj,class='fit_comp_manager'))->list_function,out=out, _extra=_extra
end

;--------------------------------------------------------------------------

pro spex::find_bad_intervals, bad=bad, _extra=_extra
  (self -> get(/obj,class='spex_fitint')) -> find_bad_intervals, bad=bad, _extra=_extra
end

;--------------------------------------------------------------------------
; Method texfile to write output text file of raw data either in raw energy bins or spex_eband bins in units of spex_units.
; binned_eband - if set, use spex_eband binning, otherwise raw
; filename - output file name.  If not set, and /nodialog is set, will be autonamed spex_instr_yyyymmdd_hhmm.txt.
; nodialog - if set, don't bring up pickfile dialog
; spex_units - 'count', 'rate', 'flux' option
; Added 09-Jul-2018, Kim.
;

pro spex::textfile, binned_eband=binned_eband, filename=filename, nodialog=nodialog, spex_units=spex_units, _extra=_extra

  checkvar, spex_units, 'counts'
  binned_eband = keyword_set(binned_eband)

  data = self->getdata(class='spex_data', /force, spex_units=spex_units)

  times = self->getaxis(/ut, /mean)
  ntimes = n_elements(times)

  units_str = self->get(/spex_data_units)
  data_name = strlowcase(units_str.data_name)

  if ~is_string(filename) then begin
    filename = 'spex_' + data_name + '_' + time2file(anytim(times[0])) + '.txt'
    if ~keyword_set(nodialog) then begin
      filename = dialog_pickfile (path=curdir(), filter='*.txt', $
        file=filename, title = 'Select output text file name',  get_path=path)
      if filename eq '' then begin
        print, 'No output file selected.  Aborting.'
        return
      endif
    endif
  endif

  if binned_eband then begin
    en = self->get(/spex_eband)
    nen = n_elements(en[0,*])
    ydata = (self->get(/obj,class='spex_data')) -> bin_data (data=data, intervals=en, er=z_err)
  endif else begin
    en = self->get(/spex_ct_edges)
    nen = n_elements(en[0,*])
    ydata = data.data
  endelse

  sen = str_replace(format_intervals(en,format='(g15.3)'),' to ', '-')
  format = '(a26,' + trim(nen+1) + 'a15)'
  head = [strupcase(data_name) + ' data for time interval: ' + anytim(times[0],/vms) + ' to ' + anytim(times[ntimes-1],/vms), $
    'Current time: ' + !stime, $
    'Data units: ' + spex_units, $
    'Energies are ' + (binned_eband ? 'binned into spex_eband bins' : 'raw energy bins') + '(' + trim(nen) + ' bins)', $
    '', $
    string('Time at center of bin', sen, format=format), $
    string('', reproduce('kev',nen), format=format), $
    '']
  out = strarr(ntimes)
  format = '(a26,' + trim(nen) + 'g15.5)'
  for i = 0,ntimes-1 do out[i] = string( anytim(times[i], /vms), ydata[*,i], format=format)

  out = [head, out]
  prstr, file=filename, out
  print, 'Text file written: ' + filename
end

;--------------------------------------------------------------------------
pro spex__define

  dummy = {spex, $
    inherits spex_gen }
end

;--------------------------------------------------------------------------
