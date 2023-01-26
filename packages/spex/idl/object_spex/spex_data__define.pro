;+
;
; Modifications:
;   19-Jul-2004, Kim.  Moved preview method out of spex__define to here
;   11-Aug-2004, Kim.  In plot_setup, check if # energy bins is 1, and if so,
;     don't try to plot spectrum or spectrogram.
;   16-Sep-2004, Kim.  Added show_err option in plot_setup method.
;   18-Nov-2004, Sandhia.  Added code to manage XSM-specific data.
;   17_Aug-2005, Kim.  Added spex_deconvolved, spex_pseudo_livetime to preview args.
;   13-Feb-2006, Kim.  In set method, handle invalid FITS files - use status keyword in
;     getpar call and print error message for bad or not found files
;     In plot_setup, call get_plot_title to get title
;	26-May-2006, Kim.  Added spex_data_pos to args.
;	23-Jun-2006, Kim.  Added get_source_pos method
;	30-Jun-2006, Kim.  Added tband (in set and intervals pre and post hook methods)
;	Jun 2006,    Kim.  Added SPEX__SPECFILE strategy, Changed preview to do a getdata
;	  and get the info from the obj (previously tried to keep preview separate, not putting
;	  the data in the obj, but it's not worth it).  Made dim1_sum=0 the default for time plots.
;	13-Sep-2006, Kim.  In plot_setup, call spex_apply_drm for errors too if photons
;	31-Oct-2006 A.B. Shah - New SOXS .les identified by 'soxs' or 'corr' in file name
;	13-Jun-2007, Kim. Added SPEX_MESSENGER_SPECFILE strategy
;	25-Mar-2008, Kim. In plot_setup, name xlog, ylog in args, and use checkvar to set them
;	28-Mar-2008, Kim. In set, check just file part of name for 'fits', not entire path
; 29-Apr-2008, Kim. use uppercase when checking for xrs in file name for MESSENGER data
; 16-Jun-2008,  Kim. spex_specfile is now allowed to be an array.  Use file[0] in most cases
;   where it is treated as a scalar.
; 31-Mar-2009, Kim.  soxs files now identified by .les extension
; 04-Aug-2009, Kim.  In intervals_pre_hook, ensure that spex_eband is float
; 17-Sep-2009, Kim. Added setunits method
; 28-Oct-2009, Kim. Added fermi gbm  strategy
; 27-Nov-2009, Kim. Added yohkoh wbs strategy.  Also, in set, added dialog keyword for popup widget
;   for data type selection for yohkoh and soxs (have multiple data types in one file)
; 18-Feb-2010, Kim. Get data type classes from function spex_datatypes, instead of hardcoded. Next -
;   should make setting the data strategy table-driven.
;   Also added spex_accum_time to set, so can change to sec.
; 07-Jul-2010, Kim.  Added select_type method.  For files that contain more than one type of data,
;   this selects which type, either through a widget selection tool, or by looking at spex_data_sel param.
;   Call new method for SOXS and YOHKOH_WBS data.  Also, call locate_file instead of just assuming
;   spex_specfile is in cur dir (will implement spex_data_dir later)
; 25-Jan-2011, Kim. Added Fermi LAT.  Restructured set method to use is_fits function and
;   get_fits_instr, and use them for determining data type.
; 23-Feb-2011, Kim. In set, call getpar and getheader with /silent
; 31-Oct-2011, Kim. In set, sort spex_specfile before setting so they will be in time order
; 7-Nov-2011, Kim. When only have one file, the sort in previous change turned spex_specfile into a vector, which
;   caused other problems for data types (like RHESSI) that are expected to have a single file. Only do sort for
;   multiple files.
;   04-Oct-2012, Kim.  Added is_image_input method
; 19-Feb-2013, Kim.  Added SMM HXRBS.
; 10-Oct-2014, Kim. If input file not recognized, set strategy to spex_any_specfile
; 24-Nov-2014, Kim. In set, trap spex_file_reader - to set, change strategy to spex_any_specfile, set it, then change strategy back.
;  Moved setting spex_any_specfile to outside of 'if fits thenelse' block. And remind user to set spex_file_reader param
; 22-Jan-2015, Kim. For spectrogram, change y title (data_unit) to 'keV'
; 25-Jun-2015, Kim. For XSM data, some FITS files don't have instrument in header, so look for ext name with 'xsm' in it
; 10-Nov-2015, Kim. Previously set weighted_sum on for anything other than counts plots.  That's wrong. For time plots,
;  weighted_sum should only be set for flux plots.  For spectra, weighted_sum should be set for rate or flux plots.  For
;  spectrograms, shouldn't be set (don't allow summing anyway).  (So biggest error was on rate time profiles - combining was
;  weighting by energy band widths, instead of just summing.)
; 11-Mar-2016, Kim. In plot_setup method, for non-spectrogram plots, y units weren't being set right for photon plots
;  (always showed counts...)
; 28-Jul-2016, Added Konus 
; 27-Sep-2016, Richard.Schwartz@nasa.gov, added STIX image cube in set method
; 29-nov-2017, Richard.Schwartz@nasa.gov, added added a bit more for STIX RATE files but still using HESSI where
; we should have STIX named procs
; 21-Mar-2019, Kim. Added SMM GRS input file option
; 19-Jan-2023, Kim. In plot_setup method, added call to get_time_plot_args, and if any added them to plot_params struct
;-
;---------------------------------------------------------------------------

pro spex_data_test

  o =  spex_data()
  ;o->set,  spex_spec =  '../ospex/hsi_spectrum_20031028_110633.fits'
  o->set, spex_specfile = 'hsi_spectrum_20020220_105002.fits'
  sp =  o->getdata()
  help,  sp,  /struct
  info =  o -> get( /info )
  help,  info,  /struct
  obj_destroy, o
  
  o = ospex(/no_gui)
  
  o->set, spex_specfile = 'hsi_spectrum_20020220_105002.fits'
  o->set, spex_drmfile = 'hsi_srm_20020220_105002.fits'
  ;o->set, spex_specfile = 'hsi_spectrum_20020421_004040.fits'
  ;o->set, spex_drmfile = 'hsi_srm_20020421_004040.fits'
  
  o-> set, spex_bk_time_int = [['20-Feb-2002 10:56:23.040', '20-Feb-2002 10:57:02.009'], $
    ['20-Feb-2002 11:22:13.179', '20-Feb-2002 11:22:47.820'] ]
    
  o->set, spex_eband=get_edge_products([3,22,43,100,240],/edges_2)
  
  o->set, spex_fit_time_inte = [ ['20-Feb-2002 11:06:03.259', '20-Feb-2002 11:06:11.919'], $
    ['20-Feb-2002 11:06:11.919', '20-Feb-2002 11:06:24.909'], $
    ['20-Feb-2002 11:06:24.909', '20-Feb-2002 11:06:33.570'] ]
    
    
  o->set, spex_erange=[19,190]
  
  o->set, fit_function='vth+bpow'
  
  o->set, fit_comp_param=[1.0e-005,1., .5, 3., 45., 4.5]
  o->set, fit_comp_free = [0,0,1,1,1,1]
  o->set, fit_comp_min = [1.e-20, .5, 1.e-10, 1.7, 10., 1.7]
  
  o->set, spex_fit_auto=1
  o->set, spex_autoplot_enable=1
  
  o->gui
  o->dofit, /all
  
  o->fitsummary
  
  obj_destroy, o
  
  print,  'test with the spex gui'
  o = ospex( /no )
  o->set,  spex_specfile = 'hsi_spectrum_20020220_105002.fits'
  o->set,  spex_drmfile =  'hsi_srm_20020220_105002.fits'
  o->gui
  obj_destroy,  o
  
  o =  spex_data()
  ;o->set,  multi_fits = 'test_cube_small.fits'
  o->set,  spex_specfile = 'test_cube_small.fits'
  data =  o->getdata()
  help,  data
  o->set,  spex_specfile = 'hsi_spectrum_20020220_105002.fits'
  o->set,  spex_drmfile =  'hsi_srm_20020220_105002.fits'
  data2=  o->getdata()
  help,  data2
  obj_destroy,  o
  
  o =  ospex()
  ;o->set,  spex_specfile = 'test_cube_small.fits'
  ;o->set,  multi_fits = 'test_cube_small.fits'
  o->set, spex_specfile = 'hsi_imagecube_20040421_1500_14x2.fits'
  
  o->gui
  obj_destroy,  o
  
  o =  ospex()
  o->set,  spex_specfile = 'cc_60s_ispec_1311.fits'
  o->xfit
  
  o->gui
  o =  spex_data()
  o->set,  spex_specfile = 'cc_60s_ispec_1311.fits'
  data =  o->getdata()
  help,  data
  
  
end

;---------------------------------------------------------------------------

function spex_data::init, source = source, _extra = _extra

  data_classes = spex_datatypes(/class)
  ret = self -> Spex_Gen_Strategy_Holder::INIT( data_classes, $
    SOURCE = source,  $
    info = {spex_data_info} )
    
  self -> set, spex_data_dir = curdir()
  self -> set, _extra = _extra
  
  return, ret
  
end

;--------------------------------------------------------------------------

; added dialog option for popping up widget for selecting data type within file for
; files that have multiple data types (yohkoh, soxs)

pro spex_data::set,  spex_data_source=spex_data_source, $
  spex_specfile = spex_specfile, $
  spex_tband=spex_tband, $
  spex_accum_time=spex_accum_time, $
  dialog=dialog, $
  spex_file_reader=spex_file_reader, $
  _extra = _extra, done = done, not_found = not_found
  
  if exist(spex_data_source) then begin
    self -> setstrategy, spex_data_source
    (self->getstrategy()) -> initialize
  endif
  
  if exist(spex_accum_time) then self->spex_gen_strategy_holder::set, spex_accum_time=anytim(spex_accum_time)
  
  if exist(spex_tband) then  self->spex_gen_strategy_holder::set, spex_tband=anytim(spex_tband)
  
  if is_string(spex_file_reader) then begin
    curr_strategy = self->getstrategy()
    curr_strat_class = obj_class(curr_strategy)
    self->setstrategy, 'SPEX_ANY_SPECFILE'
    (self->getstrategy())->set,spex_file_reader=spex_file_reader
    self->setstrategy, curr_strat_class
  endif
  
  if exist( spex_specfile ) then begin
    spex_specfile = self -> locate_file (spex_specfile, status=status)
    if status eq 0 then goto, set_extra
    
    ; If input file is being set, read enough of it to determine type
    ; of file so we can set the spex_data strategy.
    
    ; some of data types can only handle a scalar, so if only one element in spex_specfile, make scalar
    if n_elements(spex_specfile) eq 1 then spex_specfile = spex_specfile[0]
    
    is_fits = is_fits(spex_specfile[0])
    strat_set = 0
    err_msg = ''
    
    if is_fits then begin
    
      ; for FITS files, use fitsread object to get header, and instrument name, and first extension info
      reader = fitsread( filename =  spex_specfile[0] )
      status =  reader -> getstatus( )
      
      if status eq 0 and not spex_get_nointeractive() then begin
        spex_specfile =  reader -> select_file( )
        status =  reader -> getstatus( )
      endif
      
      if status eq 0 then err_msg = 'Input file not found: ' + spex_specfile[0] else begin
      
        if stregex(file_basename(spex_specfile[0]),'^KW.*pha',/bool) then instr = 'konus' else begin
        
          header = reader -> getheader(/silent)
          instr = strtrim(strlowcase(get_fits_instr(header)))
          if instr eq 'rhessi' then instr = 'hessi'
          if instr eq 'hard x-ray burst spectrometer (hxrbs)' then instr = 'hxrbs'
          
          if instr ne 'hxrbs' then begin
            reader -> set,  extension =  1
            extname =  reader -> getpar( 'EXTNAME', status=status, /silent )
            if status eq 0 then begin
              err_msg = 'ERROR reading input file: ' + spex_specfile[0]
              goto, done_strategy
            endif
            if instr eq '' and strpos(strlowcase(extname), 'xsm') ne -1 then instr='xsm'
          endif
        endelse
        
        case instr of
          'hessi': begin
          
            ; It's either an image cube or a spectrum file
            
            if (size(extname, /type) EQ 7) then begin
              if extname EQ 'HESSI Image' or extname eq 'CONTROL PARAMETERS' THEN begin
                self -> setstrategy, 'SPEX_HESSI_IMAGE'
                strat_set = 1
                ; this is ugly but we'll take care of this later
                self -> set,  im_input_fits = spex_specfile[0]
              endif else begin
                if strtrim(extname) EQ 'RATE' then begin
                  self->setstrategy, 'SPEX_HESSI_SPECFILE'
                  strat_set = 1
                endif
              endelse
            endif
          end
          'stix': begin
            ; It's either an image cube or a spectrum file

            is_cube = stregex( /boo, /fold, fxpar( header,'filetype' ),'cube')
            if is_cube THEN begin
              self -> setstrategy, 'SPEX_STIX_IMAGE'
              strat_set = 1
              ; this is ugly but we'll take care of this later
              self -> set,  im_input_fits = spex_specfile[0]
            endif else begin
              ;The stix spectral file reader object goes here in the future
                                if strtrim(extname) EQ 'RATE' then begin
                                  self->setstrategy, 'SPEX_HESSI_SPECFILE'
                                  strat_set = 1
                                  endif

            endelse
          end
    
          'xsm': begin
            self->setstrategy, 'SPEX_XSM_SPECFILE'
            strat_set = 1
          end
          
          'gbm': begin
            self -> setstrategy, 'SPEX_FERMI_GBM_SPECFILE'
            strat_set = 1
          end
          
          'lat': begin
            self -> setstrategy, 'SPEX_FERMI_LAT_SPECFILE'
            strat_set = 1
          end
          
          'hxrbs': begin
            self->setstrategy, 'SPEX_SMM_HXRBS_SPECFILE'
            strat_set = 1
          end
          
          'grs': begin
            self->setstrategy, 'SPEX_SMM_GRS_SPECFILE'
            strat_set = 1
          end
          
          'konus': begin
            self->setstrategy, 'SPEX_KONUS_SPECFILE'
            strat_set = 1
          end
          
          'minxss-1': begin
            self->setstrategy, 'SPEX_MINXSS_SPECFILE'
            strat_set = 1
          end
          
          'minxss-2': begin
            self->setstrategy, 'SPEX_MINXSS_SPECFILE'
            strat_set = 1
          end
          
          else: begin
          
            ; To check for XSM data, either instr='xsm' or the extension contains a field called
            ;      INTEGRATION_TIME or TIMEDEL.
            ttype3 = reader->getpar('TTYPE3', /silent)
            if instr eq 'xsm' or ((strtrim(ttype3) EQ 'INTEGRATION_TIME') OR (strtrim(ttype3) EQ 'TIMEDEL')) then begin
              ;NOTE: XSM warning can be removed when XSM is happy with calibration.
              msg = ['  ', $
                'WARNING:  XSM calibration is not complete. ', $
                'You may proceed only if you agree to check your results with the XSM Authorities!', $
                '']
              answer = 'y'
              print,msg
              if ~spex_get_nointeractive() then answer = xanswer ([msg,' Proceed? '], /str, default=1, /suppress)
              if answer eq 'y' then begin
                self->setstrategy, 'SPEX_XSM_SPECFILE'
                strat_set = 1
              endif
            endif
          end
          
        endcase
      endelse
    endif else begin   ; end of is_fits branch
    
      ; If input file isn't a FITS file, and has '.les' in the name, then it's SOXS ; 3/31/09 (was 'soxs')
      if strpos(spex_specfile[0], '.les') ne -1 then begin
        self -> setstrategy, 'SPEX_SOXS_SPECFILE'
        strat_set = 1
        spex_data_sel = self -> select_type( ['Si', 'CZT'], dialog=dialog)
      endif
      
      ; If input file isn't a FITS file, and has XRS in the name, then it's MESSENGER
      if strpos(strupcase(spex_specfile[0]), 'XRS') ne -1 then begin
        self -> setstrategy, 'SPEX_MESSENGER_SPECFILE'
        strat_set = 1
      endif
      
      ; If input file isn't a FITS file, and has wda in the name, then it's YOHKOH WBS
      if strpos(strlowcase(spex_specfile[0]), 'wda') ne -1 then begin
        self -> setstrategy, 'SPEX_YOHKOH_WBS_SPECFILE'
        strat_set = 1
        spex_data_sel = self -> select_type( ['GRS1', 'GRS2', 'HXS'], dialog=dialog)
      endif
      
      ; If input file isn't a FITS file, and has minxss in the name, then it's MINXSS
      if strpos(strupcase(spex_specfile[0]), 'MINXSS') ne -1 then begin
        self -> setstrategy, 'SPEX_MINXSS_SPECFILE'
        strat_set = 1
      endif
      
      ; If input file isn't a FITS file, and has daxss in the name, then it's DAXSS
      if strpos(strupcase(spex_specfile[0]), 'DAXSS') ne -1 then begin
        self -> setstrategy, 'SPEX_DAXSS_SPECFILE'
        strat_set = 1
      endif
      
    endelse
    
    if ~strat_set then begin
      ; If input file wasn't recognized by any of tests above, then set strategy to SPEX_ANY_SPECFILE
      ; and see if it can figure out what to do with it
      self->setstrategy, 'SPEX_ANY_SPECFILE'
      strat_set = 1
      message, /info, 'Input file is not recognized.  Be sure to set spex_file_reader parameter to base name of routine to read data and drm files.'
    endif
    
    if strat_set then begin
      ; sort file names so they will be in time order
      if n_elements(spex_specfile) gt 1 then begin
        s=sort(spex_specfile)
        spex_specfile = spex_specfile[s]
      endif
      self->spex_gen_strategy_holder::set, spex_specfile = spex_specfile, spex_data_sel=spex_data_sel, $
        done = done, not_found = not_found
    endif ;else begin
    ;       if err_msg eq '' then err_msg = 'Input file is not a recognized data file ' + spex_specfile[0]
    ;	   message, /info, err_msg
    ;    endelse
    
    done_strategy:
    if exist(reader) then obj_destroy,  reader
    
  endif
  
  set_extra:
  if keyword_set(_extra ) then self->spex_gen_strategy_holder::set, $
    _extra = _extra, done = done, not_found = not_found
    
end

;--------------------------------------------------------------------------
; if dialog is set, pop up widget to select data type from 'types'. If dialog
; not set, and current value of spex_data_sel isn't one of 'types', then return
; first element of 'types'

function spex_data::select_type, types, dialog=dialog

  if keyword_set(dialog) then begin
    choice = xchoice('Select data type', '   ' + types + '   ' )
    spex_data_sel = types[choice]
  endif else begin
    data_sel = self->get(/spex_data_sel)
    if not is_member(strupcase(data_sel), types) then begin
      spex_data_sel = types[0]
      message,'Setting data type selection to ' + spex_data_sel + '.', /continue
    endif
  endelse
  
  return, spex_data_sel
  
end

;--------------------------------------------------------------------------

pro spex_data::preview, out=out, nomore=nomore

  ;this_strat = self -> getstrategy()
  ;
  ;this_strat -> preview, spectrum,  errors,  livetime,  $
  ;    spex_respinfo, spex_file_time,  spex_ut_edges,  spex_ct_edges,  $
  ;    spex_area, spex_title, spex_detectors,  $
  ;    spex_interval_filter, spex_units, spex_data_name, $
  ;    spex_deconvolved, spex_pseudo_livetime, spex_data_pos, $
  ;    err_code=err_code
  
  data = self->getdata()
  specfile = arr2str(self -> get(/spex_specfile))
  
  if not is_struct(data) then begin
  
    out = ['', 'No Spectrum or Image Input Data available.  ' + specfile]
    
  endif else begin
  
    units_str = self -> get(/spex_data_units)
    spex_data_name = units_str.data_name
    spex_ct_edges = self->get(/spex_ct_edges)
    spex_ut_edges = self->get(/spex_ut_edges)
    spex_file_time = self -> get(/spex_file_time)
    respinfo = self -> get(/spex_respinfo)
    spex_area = self -> get(/spex_area)
    spex_detectors = self -> get(/spex_detectors)
    
    if specfile eq '' then specfile = 'None'
    n_eband = n_elements (spex_ct_edges[0,*])
    emm = format_intervals (minmax (spex_ct_edges))
    n_time = n_elements (spex_ut_edges[0,*])
    tr = format_intervals([spex_file_time[0], spex_file_time[1]], /ut)
    if is_number(respinfo[0]) then begin
      n = n_elements(respinfo)
      respinfo = arr2str(trim(respinfo[0:(n-1)<10],'(g8.2)'))
      if n gt 10 then respinfo=respinfo + '...'
    endif
    if not is_string(respinfo) then respinfo = 'None'
    
    out = [ $
      'Spectrum or Image File Summary', $
      'File name: ' + specfile, $
      'Data type: ' + spex_data_name, $
      '# Time Bins: ' + trim(n_time) + '  Time range: ' + tr, $
      '# Energy Bins: ' + trim(n_eband) + '  Energy range: ' + emm, $
      'Area: ' + trim(spex_area), $
      'Detectors Used: ' + spex_detectors, $
      'Response Info: ' + respinfo]
  endelse
  
end

;--------------------------------------------------------------------------
; Plot method plots raw (directly from file) counts, count rate, or count flux
;
; intervals - if set to -1, use raw intervals, otherwise see below
; pl_time - if set, plot is vs time (default)
;   If intervals aren't set, uses spex_ebands and plots them separately
;   Can pass in a set of energy bands to plot, and it plots them separately
; pl_energy - if set, plots spectrum
;   If intervals aren't set, uses raw time intervals, and sums them
;   Can pass in a set of time intervals to plot, and it sums them
;   In plotman, can still access the separate intervals
; pl_spec - if set, plots spectrogram
;   Intervals not used here - doesn't bin in energy
; dim1_sum - sum over 1st dimension.  Default for spectra is 0, for time plots is 1


function spex_data::plot_setup, $
  intervals = intervals, $
  pl_time = pl_time, $
  pl_energy = pl_energy, $
  pl_spec = pl_spec, $
  photons = photons, $
  drm_eff = drm_eff, $
  plot_params = plot_params, $
  dim1_sum = dim1_sum, $
  show_err = show_err, $
  xlog = xlog, ylog = ylog, $
  _extra = _extra
  
  error_catch = 0
  if spex_get_debug() eq 0 then catch, error_catch
  if error_catch ne 0 then begin
    catch, /cancel
    message, !error_state.msg, /cont
    message,'Error making plot.', /cont
    return,  0
  endif
  
  photons = keyword_set(photons)
  if photons and not keyword_set(drm_eff) then message,'No efficiency factors provided.  Cannot display photons.'
  
  pl_time = keyword_set(pl_time)
  pl_energy = keyword_set(pl_energy)
  pl_spec = keyword_set(pl_spec)
  if not (pl_time or pl_energy or pl_spec) then pl_time = 1
  
  data = self -> getdata(_extra=_extra)
  if not is_struct(data) then message, 'No data accumulated.'
  if data.data[0] eq -1 then message, 'No data accumulated.'
  
  label = 'Detectors: ' + self -> get(/spex_detectors)
  
  title = self -> get_plot_title(photons=photons)
  
  units_str = self -> getunits()
  data_unit = units_str.data
  
  energies = self -> getaxis(/ct_energy, /edges_2)
  nenergies = n_elements(energies[0,*])
  
  if pl_energy then begin
  
    if nenergies eq 1 then message,'Only one energy bin. Can not plot spectrum.'
    checkvar, intervals, -1  ; if no time bins passed in, don't bin
    z = self -> bin_data (data=data, intervals=intervals, er=z_err, /do_time, newedges=newedges)
    if z[0] eq -1 then return,  0
    
    xdata = energies
    ydata = z
    plot_type='xyplot'
    title = title + ' vs Energy'
    xtitle = 'Energy (keV)'
    dim1_id = format_intervals(newedges, /ut) + ' (Data with Bk)'
    dim1_vals = anytim(newedges, /vms)
    checkvar, dim1_sum, 1
    dim1_unit = ''
    checkvar, xlog, 1
    checkvar, ylog, 1
    weighted_sum = units_str.data_type ne 'counts'
    
  endif else begin
  
    times = self -> getaxis(/ut, /edges_2)
    dim1_unit = units_str.ct_edges
    checkvar, xlog, 0
    checkvar, ylog, 1
    xtitle=''
    
    if pl_time then begin
    
      checkvar, intervals, self->get(/spex_eband)  ; if no energy bins passed in, use eband
      z = self -> bin_data (data=data, intervals=intervals, er=z_err, newedges=newedges)
      if z[0] eq -1 then return,  0
      dim1_vals = newedges
      utbase = min(times)
      xdata = anytim(times) - utbase
      if n_elements(xdata) le 2 then begin
        message, /cont, 'Error - only one time interval.  Aborting plot.'
        return, 0
      endif
      ydata = nenergies gt 1 ? transpose(z) : z
      z_err = nenergies gt 1 ? transpose(z_err) : z_err
      plot_type = 'utplot'
      title = title + ' vs Time'
      dim1_id = format_intervals(newedges, format='(f9.1)') + ' keV (Data with Bk)'
      checkvar,dim1_sum, 1
      weighted_sum = units_str.data_type eq 'flux'
      more_plot_args = (self->getstrategy())->get_time_plot_args()  ;kim added 19-jan-2023
      
    endif else begin   ; spectrogram
      if nenergies eq 1 then message,'Only one energy bin. Can not plot spectrogram.'
      intervals = -1    ; don't bin in energy for spectrogram
      z = self -> bin_data (data=data, intervals=intervals, er=z_err, newedges=newedges)
      if z[0] eq -1 then return,  0
      utbase = min(times)
      xdata = average(anytim(times),1) - utbase
      if n_elements(xdata) le 1 then begin
        message, /cont, 'Error - only one time interval.  Aborting plot.'
        return, 0
      endif
      ydata = transpose(z)
      z_err = transpose(z_err)
      plot_type = 'specplot'
      title = title + ' Spectrogram'
      data_unit = 'keV'
      dim1_vals = average(newedges,1)
      dim1_id = ''
      checkvar,dim1_sum, 0
      weighted_sum = 0
    endelse
    
  endelse
  
  if photons then begin
    spex_apply_eff, ydata, drm_eff, units_str=units_str
    spex_apply_eff, z_err, drm_eff
    if ~pl_spec then data_unit = units_str.data  ; replace units count... with photons...
  endif
  
  checkvar, utbase, 0
  plot_params = { $
    plot_type: plot_type, $
    xdata: xdata, $
    utbase: anytim(utbase,/vms), $
    ydata: ydata, $
    id: title, $
    label: label, $
    data_unit: data_unit, $
    dim1_id: dim1_id, $
    dim1_vals: dim1_vals, $
    dim1_sum: dim1_sum, $
    dim1_unit: dim1_unit, $
    dim1_enab_sum: 1, $
    weighted_sum: weighted_sum, $
    xlog: xlog, $
    ylog: ylog, $
    xtitle: xtitle $
  }
  
  if is_struct(more_plot_args) then plot_params = join_struct(plot_params, more_plot_args)  ; kim added 19-jan-2023
  
  if keyword_set(show_err)  and not pl_spec then plot_params = add_tag (plot_params, z_err, 'edata')
  return, 1
end

;--------------------------------------------------------------------------

function spex_data::getunits, orig=orig

  strategy =  self -> getstrategy()
  return, strategy -> getunits( orig =  orig )
  
end

;--------------------------------------------------------------------------

pro spex_data::setunits, units, orig=orig

  if keyword_set(orig) then self -> set, spex_data_origunits=units else $
    self -> set, spex_data_units=units
    
end

;------------------------------------------------------------------------------

pro spex_data::intervals_pre_hook, $
  full_options=full_options, $
  intervals=intervals, $
  valid_range=valid_range, $
  title=title, $
  energy=energy, $
  type=type, $
  spex_units=spex_units, $
  abort=abort
  
  checkvar, spex_units, 'flux'
  
  ; make sure data has been retrieved for current settings.  May be a better way?  !!!!
  data = self -> getdata()
  if not is_struct(data) then abort = 1
  
  if abort then return
  
  checkvar, energy, 0
  
  if energy then begin
    intervals = self -> get(/spex_eband) * 1. ; make float, so won't be interpreted as integer times
    
    title='Select Energy Bands for Displaying Data'
    type='Pre-bin'
    
    if not keyword_set(full_options) then $
      if not self->valid_plot(/xyplot) then self -> plotman, /pl_energy, spex_units=spex_units
      
    valid_range = minmax (self->get(/spex_ct_edges))
  endif else begin
  
    intervals = self -> get(/spex_tband)
    
    title='Select Time Bands for Displaying Data'
    type='Pre-bin'
    
    if not keyword_set(full_options) then $
      if not self->valid_plot(/utplot) then self -> plotman, /pl_time, spex_units=spex_units
      
    valid_range = minmax (self->get(/spex_ut_edges))
    
  endelse
  
end

;------------------------------------------------------------------------------

pro spex_data::intervals_post_hook, energy=energy, _extra=_extra, bins

  if keyword_set(energy) then self -> set, spex_eband=bins else self -> set, spex_tband=bins
  
end

;--------------------------------------------------------------------------

function spex_data::get_source_pos

  strategy =  self -> getstrategy()
  return, strategy -> get_source_pos()
  
end

;-------------------------------------------------------------------------------

function spex_data::is_image_input
  strat_obj = self -> getstrategy()
  return, is_class( strat_obj, 'spex_image')
end

;--------------------------------------------------------------------------

pro spex_data__define

  dummy = {spex_data, $
    inherits spex_gen_strategy_holder }
    
end