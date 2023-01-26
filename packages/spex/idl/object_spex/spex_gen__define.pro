
;+
; Name: SPEX_GEN__DEFINE
;
; Purpose: This abstract class is inherited by some of the OSPEX objects to
;	handle general functions.
;
; Category: OSPEX
;
; Written: 2003, Kim Tolbert
; Modifications:
;   15-Jul-2004, Kim.  bad_interval keyword added to edges2index and bindata
;	16-Jul-2004, Kim.  Got rid of extra stuff in plot_setup (like overlaying
;	  background and fits - should use the gen plot_spectrum method for that.
;	19-Jul-2004, Kim.  Added filter specification in convert_fitfunc_units and get_drm_eff
;	   and added get_fitint_filter method.
;	20-Jul-2004, Kim.  Added show_filter option to plot_time
;	10-Aug-2004, Kim.  Added use_fitted keyword to get_fitint_filter
;	11-Aug-2004, Kim.  In bin_data, don't transpose data array if binning over energy and there's
;	   only one energy bin.
;	6-Oct-2004, Kim.  In valid_summ_params, use spex_fit_time_used instead of
;	   spex_fit_time_interval, and check max of absolute value of diff, not max of diff
;	10-Oct-2004, Kim.  In plot_spectrum, even though we're not actually using spex_fitint obj, do
;	   getdata on it so that info params will be set. And use spex_fit_time_used intervals.
;	17-Nov-2004, Kim.  In plot_spectrum, when overlaying fit, was showing last erange used
;	   instead of getting the erange from the spex_summ params when it could.  If multiple
;	   intervals are plotted, gets erange from first interval selected.
;	11-Apr-2005, Kim.  In plot_time, plot filter bars by adding to plotman properties,so
;		they'll still appear on zooms, prints, etc.
;	7-Jun-2005, Kim.  Added fit_color to make_fitplot_obj method so we can pass in colors for the
;		fit component overlays even if we're not using plotman
;	28-Jun-2005, Kim.  Added show_chi keyword to make_fitplot_obj method
;	17-Aug-2005, Kim. In convert_units, use spex_deconvolved to choose 'photons' or 'counts'. Added
;		allint keyword to plot_spectrum method.  Added units_widget method.
;	25-Aug-2005, Kim.  Added class= to calls to get wherever possible to speed things up/  Also
;		use f_div in do_binning in case iwidth is 0.  Also, plot_spectrum and plot_time, added
;		region numbers to labels for imagecube data.
;	12-Oct-2005, Kim.  In units_widget, put droplist in a base so can store uvalues in base -
;		otherwise crashes in 5.6
;	9-Feb-2006, Kim.  Added full srm indication to label in plot_spectrum if /photon set.
;	  Added get_plot_title function to consolidate what plot for each data type was doing.
;	  Also, remove the /replace flag on new_panel call so we keep all panels even if they have same name
;	15-Mar-2006, Kim.  Added epsilon=1.e-5 to call to find_contig_ranges
;	7-Apr-2006, Kim. Changes for addition of spectrum,model keywords:
;	  In get_drm_eff, make fit func obj and pass that to get_diagonal instead of passing
;	  func and params.  Added make_func_obj method.  In make_fitplot_obj, call make_func_obj
;     to create and set up fit function, instead of doing it in make_fitplot_obj,
;     and added short summary of keyword values to plot label for fit functions
;	27-Jun-2006, Kim.  Added func_obj to call to apply_drm (in convert_fitfunc_units) to
;	  accommodate new drm_mod function for RHESSI (to fine-tune detector params in DRM)
;	2-Jul-2006, Kim.  Added option to plot_spectrum method to bin in spex_tband time intervals
;	30-Jul-2006, Kim.  If want allint, but there are > 100, don't plot them separately - too many - bin
;	  into a single bin containing from min to max of all times.
;	1-Aug-2006, Kim.  In valid_summ_params, corrected bug - use interval[i] not i in loop
;	1-Sep-2006, Kim.  Define linestyles for fitplot obj even if using plotman (for b/w)
;	11-Jul-2007, Kim. edges2index was wrong - was using start/end to find which bin to
;	  put interval in and then adjusting for no overlap, but that didn't work unless working
;     on all bins at once, i.e. if just binning into 3-6kev, and next call binned into 6-12kev,
;     they might overlap. Now use mid points.  Also this will take care of input bins that overlap.
;	4-Nov-2007, Kim. In do_plot, set _extra in call to setdefaults.
;	8-Feb-2008, Kim.  In make_fitplot_obj:If pileup_mod is part of function, don't apply it for
;	  separate component - just for the composite function.  NOTE:  if change order of plotting,
;	  need to change this - right now, composite function is plotted first.
;	5-Mar-2008, Kim. In make_fitplot_obj method, added this_interval keyword in call to make_func_obj
;	25-Mar-2008, Kim.
;	1. In plot_spectrum, if we're not using fitted valued, print label on plot 'No fit done.'
;	2. In plot_spectrum, call check_times if need to use spex_summ info.
;	3. In valid_summ_params, use new spex_allow_diff_energy param to decide whether to check if energies
;	   values in spex_summ and current match.
;	27-Mar-2008, Kim.  In make_fitplot_obj, removed show_chi arg.  If chisq is passed in or retrieved in
;	  call to make_func_obj, then inclue chisq in label for plot. Otherwise, put 'Not fit' in label. In
;	  make_func_obj, don't set chisq unless getting vals from saved fit results.
; 8-Apr-2008, Kim.  In do_plot, use plot_params.panel_id (new tag) for panel id if it exists.  Previously
;   used plot_params.id (which is title), but in some cases want that to be blank.
; 23-Jun-2008, Kim.  Added calc_func_components method to calculate and return function components
;   in any units in photons or counts, separately or combined.  Code was part of make_fitplot_obj, which
;   now calls calc_func_components
; 29-Jul-2008,, Kim.  In make_fitplot_obj, added fit_legend_loc keyword. Also if dim1_color passed in,
;   use it even if using plotman.
; 15-Aug-2008, Kim. In plot_spectrum, call empty to ensure last erange line is drawn.
; 6-Nov-2008, Kim.  In plot_spectrum, make it so energy range lines persist (by using more_commands in plotman)
; 4-Dec-2008, Kim.  In calc_func_components, use spex_summ_energy if use_fitted is set (previously crashed here
;   if hadn't read spectrum file since didn't have ct_energy).  Add energy bins to output structure.
;   In make_fitplot_obj, use energy edges from structure returned by calc_func_components.
;   In valid_summ_params, added more info to question about continuing if energy bins different
; 16-Jan-2009,, Kim.  Added get_erange method. In plot_spectrum, use add_plot and oplot_xmarkers to show energies.
; 30-Mar-2009, Kim. Added obj and get_plot_obj keywords to do_plot, plot_spectrum, and plot_time.  (get_plot_obj is
;   passed through _extra to do_plot).  If get_plot_obj is set, just pass back plot object, don't 
;   plot it.  This feature is used in show_synop.
; 10-Jul-2009, Kim.  Now allow spex_roi_use to multi-valued, so check spex_roi_use[0] in plot routines
; 04-Aug-2009, Kim.  In do_plot, destroy obj after plotting.  Memory leak. 
; 17-Aug-2009, Kim.  In convert_fitfunc_units and get_drm_eff, destroy fit obj.  memory leak.
; 17-Sep-2009, Kim. In convert_units, return structure (previously returned .data, and error in keyword arg)
;   Also in calc_count, calc_rate, and calc_flux, convert obsdata, eobsdata, bkdata, and ebkdata (in addition
;   to primary data and error) if they're in the structure.  (Added them to fitint data structure)
; 26-Sep-2009, Kim. In valid_summ_params, added summ to args.  If passed use it, otherwise get(/spex_summ)
; 02-Oct-2009, Kim. In calc_count, calc_rate, calc_flux, don't do bk if bk[0] = -1
; 28-Oct-2009, Kim. In valid_drm, use compare_float instead of same_data
; 28-Nov-2009, Kim. In plot_time and plot_spectrum, added a catch (can't use spex_insert_catch
;   because it's not a function (can't return,-1)
; 9-Dec-2009, Kim. In bin_data and do_binning, added ave_err and ave_corr to  allow different 
;  handling of data where error should be averaged (e.g. bk data computed from fit).  ave_err
;  option is normally set in spex_bk::bin_data method depending on value of spex_bk_poisson_error
; 01-Jun-2010, Kim. Added ph_edges keyword to calc_func_components
; 07-Jul-2010, Kim. Added locate_file and get_fitint_time methods.
; 26-Aug-2011, Kim. In calc_count and calc_rate, use bkltime not data.ltime for bk since they could be different
; 17-Oct-2011, Kim. In get_drm_eff, pass time for current index in loop, not full time array, to get_diagonal
; 09-Nov-2011, Kim. In do_plot, call obj_new with utbase keyword
; 26-Nov-2011, Kim. In calc_func_components, call convert with disable_albedo set for all components 
;   except 0 (sum of all comp). Previously pileup was disabled for all except 0, but albedo corr. was
;   being applied to separate components as well as composite.  Also, now compute and return pileup
;   and albedo contribution (by taking differences) so they can be shown on plots.  Also made yvals double
;   (it matters now because taking differences)
; 15-Dec-2011, Kim. In do_plot, if utbase isn't in plot_params struct, add it with val=0.
; 09-Aug-2012, Kim. In make_fitplot_obj, added more colors to dim1_color array
; 14-Nov-2012, Kim. Added plot_spectrum_resid method to automatically stack spectrum and resid plots, 70%,30%
; 30-Jul-2013, Kim. In bin_data, for energy binning case with nenergy>1, wasn't passing back the binned ltime.
; 17-Sep-2013, Kim. In bin_data, when binning over energy, don't bin ltime, just reproduce the ltime to correct dim.
; 09-Apr-2014, Kim. In calc_func_components, call convert_fitfunc_units with this_strat so apply_drm knows when it's
;   working on the entire function (this_strat=''), or just a component. (Needed for adding in components that don't
;   go through DRM. )
; 16-Apr-2014, Kim. In calc_func_components, added obj=self to fit_comp_kw call to enable retrieving source_xy and source_angle
; 26-Oct-2014, Kim. In do_plot, combine _extra with plot_params since xyplot obj now handles
;  x,y, and main title.  Also don't need to do anything special with xtitle now.
; 29-Sep-2015, Kim. In bin_data, when binning by energy, need to transpose ltime so it's [nenergy_bin,ntime] like result and eresult
; 11-Apr-2017, Kim. In plot_spectrum_resid previously called plot_spectrum with legend_loc=0 - removed that
; 03-Oct-2017, Kim. Added spex_image_spectrum_source info to label on plots
; 13-Nov-2017, Kim. Brian accidently found that passing plotman_obj to plot_spectrum or plot_time caused
;   an error that closed the entire ospex gui's plotman (because it gets passed through _extra, it goes into plotman's
;   plot_args structure and then gets destroyed in update_panel method. Note: plotman was changed to not allow saving
;   plotman in plot_args).  Made changes to allow plotman_obj keyword because it adds flexibility - now can send a
;   plot to a different plotman instance instead of one in ospex gui. Changes include calling out 
;   plotman_obj=plotman_obj keyword in call and changing plotman keyword in do_plot to use_plotman (to avoid
;   keyword duplicate error). Note, if plotman_obj is not passed in (and no_plotman isn't set), uses ospex gui plotman.
; 24-Aug-2022, Kim. In calc_func_components method, removed this_strat=strat_arr[i] keyword in call to get full spectrum without
;   pileup_mod correction. That keyword was preventing apply_drm from adding in the ...nodrm component.
; 19-Jan-2023, Kim. In plot_time method, added check for histogram in plot_params and calling setdefaults with psym (for minxss)
;
;-
;---------------------------------------------------------------------------


; function to check if current plotman window is a valid plot window of the specified
; type ('xyplot', 'utplot', 'image', 'spec' passed in _extra)

function spex_gen::valid_plot, _extra=_extra, nocreate=nocreate
plotman_obj = self -> get_plotman_obj(valid=valid, nocreate=nocreate)
return, valid ? plotman_obj -> valid_window(_extra=_extra) : 0
end

;--------------------------------------------------------------------------

; function to return plotman object, and whether or not it's valid to use.

function spex_gen::get_plotman_obj, valid=valid, nocreate=nocreate

checkvar, nocreate, 0

plotman_obj = self -> get(/spex_plotman_obj)

; if plotman object doesn't already exist, create it, unless $OSPEX_NOINTERACTIVE or nocreate is set

if not is_class(plotman_obj, 'plotman') then begin
	if not (spex_get_nointeractive() or nocreate) then begin
		plotman_obj = plotman(/multi)
		self -> set, spex_plotman_obj=plotman_obj
	endif
endif

valid = is_class(plotman_obj, 'plotman')

return, plotman_obj
end

;--------------------------------------------------------------------------

pro spex_gen::plot, use_plotman=use_plotman, _extra=_extra

status = self -> plot_setup (plot_params=plot_params, _extra=_extra)

if status then self -> do_plot, plot_params=plot_params, use_plotman=use_plotman, _extra=_extra

end

;--------------------------------------------------------------------------
; function to draw plot, either in plotman window (if plotman was specified and
; plotman object is valid) or a simple plot window.
;
; the /replace keyword on new_panel call means replace panel of same description (so
; we don't generate a gazillion panels).

; if obj keyword is present, just pass out the xy or utplot object.  Don't draw the plot.

pro spex_gen::do_plot, plot_params=plot_params, use_plotman=use_plotman,  $
  get_plot_obj=get_plot_obj, obj=obj, plotman_obj=plotman_obj, _extra=_extra

error_catch = 0
if spex_get_debug() eq 0 then catch, error_catch
if error_catch ne 0 then begin
	catch, /cancel
	message, !error_state.msg, /cont
	message,'Error making plot.', /cont
	return
endif

if keyword_set(use_plotman) then if ~is_class(plotman_obj, 'plotman') then plotman_obj = self -> get_plotman_obj()
valid = is_class(plotman_obj, 'plotman')

plot_type = plot_params.plot_type
; if plot_params struct doesn't have a utbase tag, add with val=0. so we can call obj_new below with it
if ~tag_exist(plot_params,'utbase') then plot_params = add_tag(plot_params,0.d0,'utbase')

obj = (tag_exist(plot_params,'edata')) ? $
	obj_new(plot_type, plot_params.xdata, plot_params.ydata, plot_params.edata, utbase=plot_params.utbase) : $
	obj_new(plot_type, plot_params.xdata, plot_params.ydata, utbase=plot_params.utbase)

; take these tags out of plot_params so they won't get passed to set routine.
plot_params = rem_tag ( plot_params, ['plot_type', 'xdata', 'ydata', 'edata'] )

plot_params = join_struct(_extra, plot_params)

obj -> set, _extra=plot_params

; This is necessary because the plot default is histogram, and for some data (minxss) want to override that
hist_tag = get_tag_value(plot_params, /histogram, err=err)
if err eq 0 && hist_tag eq 0 then psym_extra = 0 ; psym_extra won't have a value if plot_params.histogram didn't exist or isn't 0

if not keyword_set(get_plot_obj) then begin

  if valid then begin

	  status = plotman_obj -> setdefaults (input=obj, plot_type=plot_type, _extra=_extra, psym = psym_extra)
  ;	plotman_obj -> set, _extra=_extra
  ; if plot_params.xtitle ne '' then plotman_obj -> set, xtitle=plot_params.xtitle
	  panel_id = tag_exist(plot_params,'panel_id') ? plot_params.panel_id : plot_params.id
	  plotman_obj -> new_panel, panel_id;, /replace

  endif else obj->plot, _extra=_extra
  
  obj_destroy,obj  ; added 4-aug-2009
  
endif

end

;--------------------------------------------------------------------------

pro spex_gen::plotman, class_name=class_name, _extra=_extra

if keyword_set(class_name) then $
	(self->get(/obj,class_name=class_name))->plot, /use_plotman, _extra=_extra $
else $
	self -> plot, /use_plotman, _extra=_extra

end

;--------------------------------------------------------------------------

; function to return strategy name associated with a component name.  If there's more than
; one component with comp_name, returns first strategy name found.

; Either pass in the comp_name in keyword comp_name, or pass in component name or strategy
; in _extra, via e.g. /vth

function spex_gen::name2strat, comp_name=comp_name, strategy_name=strategy_name, _extra=_extra

; this seems silly, but it allows us to call this routine without having to figure out what's
; in _extra
if exist(strategy_name) then return, strategy_name

compman_name = self -> get(/compman_name, class='fit_comp_manager')
compman_strategy = self -> get(/compman_strategy, class='fit_comp_manager')

;compman_name = self -> framework::get(/compman_name)
;compman_strategy = self -> framework::get(/compman_strategy)

ind = -1

if keyword_set(comp_name) then begin
	ind = where (comp_name eq compman_name)
endif else begin
	if keyword_set(_extra) then $
		ind = where_arr( strlowcase(compman_name), strlowcase(tag_names(_extra)) )
endelse

return, ind[0] eq -1 ? '' : compman_strategy[ind[0]]
end

;--------------------------------------------------------------------------

function spex_gen::strat2name, strategy_name=strategy_name
return,  ssw_strsplit(strategy_name, '#', /head)
end

;--------------------------------------------------------------------------
; function to return first/last indices in an array of parameters that apply to
; that strategy component (e.g. for 'vth+bpow' if strategy for bpow is selected,
; then the returned values would be [2,5] since vth uses indices [0,1]
function spex_gen::strat2paramindex, strategy_name

;compman_name = self -> framework::get(/compman_name)
compman_name = self -> get(/compman_name, class='fit_comp_manager')
;compman_strategy = self -> framework::get(/compman_strategy)
compman_strategy = self -> get(/compman_strategy, class='fit_comp_manager')

ind = where (strategy_name eq compman_strategy, count)

if count ne 0 then begin
	endi = -1
	for i = 0,ind do begin
		starti = endi + 1
		endi = starti + fit_comp_defaults(compman_name[i], /nparam) - 1
	endfor
	return, [starti,endi]
endif else return, [-1,-1]

end

;--------------------------------------------------------------------------

; this function is necessary because in _extra we might have this_component='xxx' or
; /vth, etc for which function to get params for, but since these aren't parameters
; to actually retrieve, framework will return a structure even when only one
; parameter is requested.  Tried changing framework, but that broke other things
; because they (e.g. spectrogram /binning) expect this behavior.
; so with this function if just getting a single param, it won't return a structure.

function spex_gen::getfit,  _extra=_extra

spex_fit = self -> get(/obj, class='spex_fit')
if size(spex_fit,/tname) ne 'OBJREF' then ret = self -> get(_extra=_extra) else $
	ret = spex_fit -> get(_extra=_extra)
if size(ret,/tname) eq 'STRUCT' then if n_tags(ret) eq 1 then return, ret.(0)
return, ret
end

;--------------------------------------------------------------------------

function spex_gen::getalg, _extra=_extra
fitalg = self -> get(/obj,class_name='spex_fitalg')
if size(fitalg,/tname) ne 'OBJREF' then return, -1
this_alg = fitalg -> getstrategy()
return, this_alg -> getalg(_extra=_extra)
end


;--------------------------------------------------------------------------

function spex_gen::getaxis, ct_energy=ct_energy, ph_energy=ph_energy, ut=ut, _extra=_extra

case 1 of
	keyword_set(ct_energy): edges = self -> get(/spex_ct_edges)
	keyword_set(ph_energy): edges = self -> get(/spex_drm_ph_edges, class='spex_drm')
	else: edges = self -> get(/spex_ut_edges)
endcase

if keyword_set(_extra) then $
	return, get_edges(edges, _extra=_extra) else $
	return, get_edges(edges, /mean)

end

;---------------------------------------------------------------------------

; generic intervals method.  if /full_options is set, brings up widget that
; lets you set intervals in about 20 different ways.  Otherwise, just
; brings up the graphical interval selector.
; intervals_pre_hook is defined in each concrete class to do the non-generic
; stuff beforehand - like setting the initial interval values and plotting
; the right type of plot.
; intervals_post_hook is also defined in each concrete class to save the
; intervals in the correct parameter.

pro spex_gen::intervals, full_options=full_options, energy=energy, _extra=_extra

if spex_get_nointeractive(/message) then return

abort = 0

self -> intervals_pre_hook, full_options=full_options, $
	intervals=intervals, $
	valid_range=valid_range, $
	title=title, $
	type=type, $
	energy=energy, $
	abort=abort, $
	 _extra=_extra

if abort then return

bins = -1

plotman_obj = self -> get_plotman_obj()

if keyword_set(full_options) then begin
	bins = xsel_intervals( energy=energy, $
	                       input_intervals=intervals, $
	                       plotman_obj=plotman_obj, $
	                       group=plotman_obj->get(/plot_base), $
	                       valid_range=valid_range, $
	                       title=title, $
	                       type=type, $
	                       /show_start, $
	                       /force, $
	                       _extra=_ref_extra )

endif else begin

    plotman_obj -> intervals, title=title, $
                            type=type, $
                            intervals=intervals, $  
                            /show_start, $
                            /force, $
                            _extra=_ref_extra

    bins = plotman_obj->get( /intervals )


endelse

if bins[0] ne -99 then self -> intervals_post_hook, energy=energy, _extra=_extra, bins

end

;--------------------------------------------------------------------------
; units are counts unless spex_units is specified (choices are 'flux','counts','rate')
; Note that spex_units does not get passed to process method.
; Added original_units. This is so that this problem can be avoided:
;  When a user calls getdata with a requested spex_units, setunits is called to record
;  the units returned.  But internally, the spex process methods call getdata, and we don't want
;  to affect the spex_xxx_units then, or it will change what user did.  So internally in
;  spex, we'll call this with /original_units, and not change spex_xxx_units, but use
;  getunits(/orig).

function spex_gen::getdata, $
	spex_units=spex_units, $
	original_units=original_units, $
	class_name=class_name, $
	_extra=_extra

@spex_insert_catch

if keyword_set(class_name) then $
	return, self -> framework::getdata(class_name=class_name, $
		spex_units=spex_units, $
		original_units=original_units, $
		_extra=_extra)

data = self -> framework::getdata(_extra=_extra)

if not is_struct(data) or keyword_set(original_units) then return, data

checkvar, spex_units, 'counts'
units_str = self -> getunits(/orig)

if data.data[0] ne -1 and (spex_units ne units_str.data_type) then $
	data = self -> convert_units (data, spex_units, units_str)

;photons = 1
;if keyword_set(photons) then eff = (self->get(/obj,class='spex_drm')) -> get_diagonal()
;data.data = data.data / eff
;data.edata = data.edata / eff ; !!! check

self -> setunits, units_str
return, data

end

;---------------------------------------------------------------------------
; Function to convert between any of the possible units, counts, rate or flux
; data structure is returned, also units_str is changed
; data is the data structure (.data, .edata, .ltime)
; data_unit - string of units to convert to ('counts', 'rate', or 'flux')

function spex_gen::convert_units, data, data_unit, units_str

; allow fuzzy units selection (anything close is OK)
units_choices = ['flux', 'rate', 'counts']
index = str_fuzzy (['fl', 'ra', 'co'], data_unit)
data_unit = index eq -1 ? 'flux' : units_choices[index]

ph_or_c = self->get(/spex_deconvolved) ? 'photons' : 'counts'



case data_unit of
	'counts':  begin
		new_data = self -> calc_count(data=data, er=er, units_str=units_str)
		units_str.data = ph_or_c
		units_str.data_type = 'counts'
		end
	'rate': begin
		new_data = self -> calc_rate(data=data, er=er, units_str=units_str)
		units_str.data = arr2str([ph_or_c,units_str.time], ' ')
		units_str.data_type = 'rate'
		end
	'flux': begin
		new_data = self -> calc_flux(data=data, er=er, units_str=units_str)
		units_str.data = arr2str([ph_or_c,units_str.time,units_str.area,units_str.energy], ' ')
		units_str.data_type = 'flux'
		end
endcase

;data.(0) = new_data
;data.(1) = er

return, new_data
end
;---------------------------------------------------------------------------

function spex_gen::calc_count, data=data, eresult=eresult, units_str=units_str

@spex_insert_catch

obs = tag_exist(data, 'obsdata')
bk =  tag_exist(data, 'bkdata')
if bk then if data.bkdata[0] eq -1 then bk = 0
 
if units_str.data_type eq 'counts' then begin
	return, data
endif else begin
	data = self -> calc_rate (data=data, eresult=eresult, units_str=units_str)
	data.(0) = data.(0) * data.(2)
	data.(1) = data.(1) * data.(2)
	if obs then begin
	  data.obsdata = data.obsdata * data.(2)
	  data.eobsdata = data.eobsdata * data.(2)
	endif
	if bk then begin
    data.bkdata = data.bkdata * data.bkltime
    data.ebkdata = data.ebkdata * data.bkltime
  endif	 
endelse

return, data
end

;---------------------------------------------------------------------------

function spex_gen::calc_rate, data=data, eresult=eresult, units_str=units_str

@spex_insert_catch

obs = tag_exist(data, 'obsdata')
bk = tag_exist(data, 'bkdata')
if bk then if data.bkdata[0] eq -1 then bk = 0

case units_str.data_type of

	'counts': begin
		data.(0) = f_div(data.(0), data.(2))
		data.(1) = f_div(data.(1), data.(2))
		if obs then begin
		  data.obsdata = f_div(data.obsdata, data.(2))
		  data.eobsdata = f_div(data.eobsdata, data.(2))
		end
    if bk then begin
      data.bkdata = f_div(data.bkdata, data.bkltime)
      data.ebkdata = f_div(data.ebkdata, data.bkltime)
    end
    end		
    
	'rate': 
	
	'flux': begin

		ewidth = self -> getaxis(/ct_energy, /width)
		ntime = n_elements(data.(0)[0,*])
		if ntime gt 1 then ewidth = ewidth # make_array(ntime, /float, value=1.)
		area = self -> get(/spex_area)
		data.(0) = data.(0) * ewidth * area
		data.(1) = data.(1) * ewidth * area
    if obs then begin
      data.obsdata = data.obsdata * ewidth * area
      data.eobsdata = data.eobsdata * ewidth * area
    end
    if bk then begin
      data.bkdata = data.bkdata * ewidth * area
      data.ebkdata = data.ebkdata * ewidth * area
    end       
		end
endcase

return, data
end

;---------------------------------------------------------------------------

function spex_gen::calc_flux, data=data, eresult=eresult, units_str=units_str

@spex_insert_catch

obs = tag_exist(data, 'obsdata')
bk =  tag_exist(data, 'bkdata')
if bk then if data.bkdata[0] eq -1 then bk = 0

if units_str.data_type ne 'flux' then begin 

	data = self -> calc_rate (data=data, eresult=eresult, units_str=units_str)

	ewidth = self -> getaxis(/ct_energy, /width)
	ntime = n_elements(self -> getaxis(/ut, /mean))
	if ntime gt 1 then ewidth = ewidth # make_array(ntime, /float, value=1.)
	area = self -> get(/spex_area)
	data.(0) = f_div (data.(0), ewidth) / area
	data.(1) = f_div (data.(1), ewidth) / area
  if obs then begin
    data.obsdata = f_div(data.obsdata, ewidth) / area
    data.eobsdata = f_div(data.eobsdata, ewidth) / area
  end
  if bk then begin
    data.bkdata = f_div(data.bkdata, ewidth) / area
    data.ebkdata = f_div(data.ebkdata, ewidth) / area
  end
         
endif

return, data

end

;---------------------------------------------------------------------------
; generic binning - binning is over the second index in data
; index - start/end indices for each bin
; data, edata - input data and error to bin
; width, mult - bin width and multiplying factors
; ave_err - if set, average errors
; ave_corr - if ave_err is set, use ave_corr to correct data before averaging (e.g.
;   want to average in rate, so if input is counts, then ave_corr will be the livetimes)
;   
; Output of function is binned data
; Output keyword:
; eresult - binned error
; 
; Added ave_err, ave_corr 9-Dec-2009 to allow different handling of
;   data where error should be averaged (e.g. bk data computed from fit)
function do_binning, index, data, edata, width=width, mult=mult, $
  eresult=eresult, ave_err=ave_err, ave_corr=ave_corr

checkvar, mult, 1.
checkvar, width, 1.
checkvar, ave_corr, 1.
checkvar, ave_err, 0

do_er = arg_present(eresult)

nbin = n_elements(index[0,*])
ndata = n_elements(data[*,0])

result = fltarr (ndata, nbin)
if do_er then eresult = fltarr (ndata, nbin)

use_width = n_elements(width) gt 1
width_is_2d = size(width, /n_dim) eq 2

for i = 0, nbin-1 do begin

	i1 = index[0,i]  & i2 = index[1,i]

	if i1 eq i2 then begin
		result[*,i] = data[*,i1]
		if do_er then eresult[*,i] = edata[*,i1]
	endif else begin
		if use_width then iwidth = width_is_2d ? ftotal(width[*,i1:i2], 2) : ftotal (width[i1:i2]) $
			else iwidth = width[0]
		factor = n_elements(mult) gt 1 ? mult[*,i1:i2] : mult[0]
		result [*,i] = f_div (ftotal (factor * data[*,i1:i2], 2), iwidth)
		if do_er then begin
		  if ave_err then begin
        ; if ave_err set, then just average error. Use ave_corr because averaging
        ; in time can only be done on rates or flux, and averaging in energy can only be done
        ; on flux, ave_corr is used to change to needed units before averaging, then undo.
        ; if ave_corr is 1., can't index it so put what we'll use into ave_corr_use
        ; (Note: averaging is the same as saying relative error on combined bin (as rate or flux)
        ; has to be the same as relative error on a single bin)
        ave_corr_use = same_size(edata,ave_corr) ? ave_corr[*,i1:i2] : ave_corr		  
		    eresult[*,i] = average( edata[*,i1:i2] / ave_corr_use, 2) * total(ave_corr_use,2)
		  endif else begin
		    eresult [*,i] = f_div (sqrt (ftotal ((factor * edata[*,i1:i2])^2, 2)), iwidth)
		  endelse
		endif

	endelse
endfor

return, result
end

;---------------------------------------------------------------------------
; INPUT keywords:
; data input is a structure whose tags are:
;  1. data (nenergy, ntime)
;  2. error in data (nenergy, ntime)
;  3. livetime (nenergy, ntime)
; Units of data structure can be counts, rate, or flux.
; intervals - time or energy intervals to bin into
; units_str - units structure of input data structure. (if not provided, does getunits()
; do_time_bin - if set, bin into time bins. Default is energy bins.
; ave_err - if set, combine errors by averaging (e.g. for fitted background)
; 
; OUTPUT of function is binned data (i.e. just data.data binned)
; 
; OUTPUT keywords:
; eresult - binned errors corresponding to binned data
; ltime - binned livetime corresponding to binned data
; newedges - edges of bins used (may not be exactly what user asked for - force to data edges)
; index - array of start/end indices to use for each bin
; bad_interval - indices of intervals for which no data was found
; 
; bad_interval added 7/14/04
; ave_err added 9-Dec-2009
; transpose ltime in energy binning - 29-sep-2015

function spex_gen::bin_data, data=data, intervals=intervals, units_str=units_str, $
	do_time_bin=do_time_bin, eresult=eresult, ltime=ltime, newedges=newedges, index=index, $
	bad_interval=bad_interval, ave_err=ave_err

@spex_insert_catch

if not keyword_set(units_str) then units_str= self -> getunits()

index = self -> edges2index (intervals=intervals, $
	newedges=newedges, $
	do_time_bin=do_time_bin, $
	got_int=got_int, $
	bad_interval=bad_interval)

if got_int then begin

	if (size(index,/dim))[0] ne 2 then begin
		msg = 'No data found in requested bins.'
		@spex_message
	endif

	if keyword_set(do_time_bin) then begin

		; binning over times
		; if not counts, then have to multiply by livetime before binning
		; if counts, and error is being averaged, have to pass ltime (in ave_corr), so can average rate.
		width = units_str.data_type eq 'counts' ? 1. : data.(2)
		mult = units_str.data_type eq 'counts' ? 1. : data.(2)
		ave_corr = units_str.data eq 'counts' ? data.(2) : 1.
		result = do_binning (index, data.(0), data.(1), width=width, mult=mult, $
		  eresult=eresult, ave_err=ave_err, ave_corr=ave_corr)
		ltime = do_binning (index, data.(2))

	endif else begin

		; binning over energies.  transpose since do_binning bins over second index
		; if flux, then construct array of energy widths (ntime,nenergy) to multiply by before binning
		; if not flux, and error is being averaged, have to pass energy width (in ave_corr), so can average flux.
		nenergy = n_elements(self -> getaxis(/ct_energy, /mean))
		en_width = self -> getaxis(/ct_energy, /width)
		width = units_str.data_type eq 'flux' ? en_width : 1.
		ntime = n_elements(data.(0)[0,*])
		mult = units_str.data_type eq 'flux' ? make_array(ntime,/float,val=1.) # width : 1.
		ave_corr = units_str.data_type ne 'flux' ? make_array(ntime,/float,val=1.) # en_width : 1.

		if nenergy gt 1 then begin
			result = do_binning (index, transpose(data.(0)), transpose(data.(1)), width=width, mult=mult, $
			  eresult=eresult, ave_err=ave_err, ave_corr=ave_corr)
			result = transpose(result)
			eresult = transpose(eresult)
			nbin = n_elements(index[0,*])
;     ltime = do_binning(index,transpose(data.(2)))			
			; since binning in energy, and live time is same for all energies for each time bin, don't bin, just reproduce nbin times
			ltime = reform((data.(2))[0,*])
			ltime = transpose(reproduce(ltime, nbin))
		endif else begin
			result = do_binning (index, data.(0), data.(1), width=width, mult=mult, $
			  eresult=eresult, ave_err=ave_err, ave_corr=ave_corr)
;     ltime = do_binning(index,data.(2))			  
			ltime = reform(data.(2))
		endelse
	endelse

endif else begin
	result = data.(0)
	eresult = data.(1)
	ltime = data.(2)
endelse

return, result

end

;---------------------------------------------------------------------------

;function spex_gen::bin_ltime, data, intervals=intervals, newedges=newedges
;
;checkvar, intervals, -1
;got_int = (size(intervals, /dim))[0] eq 2
;
;if got_int then begin
;
;	ind = self -> edges2index (intervals=intervals, newedges=newedges, do_time_bin=do_time_bin)
;
;	if (size(ind,/dim))[0] ne 2 then begin
;		msg = 'No data found in requested bins.'
;		@spex_message
;	endif
;
;	result = do_binning (ind, data.(0))
;endif else begin
;	result = data.(0)
;	newedges = self->getaxis(/ut, /edges_2)
;endelse
;
;return, result
;
;end

;---------------------------------------------------------------------------
; function to return indices ([2,n] - start/end indices) of data bins in each interval.
;
; bad_interval added 7/14/04 to return indices of intervals for which no data was found

function spex_gen::edges2index, $
	newedges=newedges, $
	intervals=intervals, $
	do_time_bin=do_time_bin, $
	got_int=got_int, $
	bad_interval=bad_interval

@spex_insert_catch

checkvar, intervals, -1
got_int = (size(intervals, /dim))[0] eq 2

bad_interval = -1

if keyword_set(do_time_bin) then begin
	edges = self -> getaxis(/ut, /edges_2)
	int = got_int ? anytim(intervals,/sec) : edges
	newedges = edges
	; if e.g. energies instead of times entered, end up with one value
	if (size(int, /dim))[0] ne 2 then begin
		message,'Invalid time intervals. Aborting.', /cont
		got_int = 0
		return, -1
	endif
endif else begin
	edges = self -> getaxis(/ct_energy, /edges_2)
	int = got_int ? intervals : edges
	newedges=edges
endelse

n_int = n_elements(int[0,*])
index = lonarr(2,n_int)

; got_int is true if we're really binning into intervals. Otherwise, just return indgen.
if got_int then begin

; If large times, make them relative to min time.    maybe not necessary.
;	if int[0,0] gt 1.e6 then begin
;		tref = min ([int[*],edges[*]])
;		int = int - tref
;		edges = edges - tref
;	endif

	; for each bin, find indices such that midpoint of edges is within bin start/end

	mids = get_edges(edges, /mean)
	for i=0,n_int-1 do begin
		q = where (mids ge int[0,i] and mids lt int[1,i], count)
		case count of
			0: index[*,i] = [-1,-1]
			1: index[*,i] = [q,q]
			else: index[*,i] = minmax(q)
		endcase
	endfor

	;print,'in edges2index: index= ', index

;	; for each bin, find edges such that edges start is <= bin end and edges end is >= bin start
;	for i = 0, n_int-1 do begin
;		delta = edges[0,0] gt 1.e6 ? 2.e-3 : 0.
;		q = where ( edges[0,*] lt int[1,i]-delta and edges[1,*] gt int[0,i]+delta, count )
;		case count of
;			0:  index[*,i] = [-1,-1]
;			1:  index[*,i] = [q,q]
;			else: index[*,i] = minmax(q)
;		endcase
;
;		; Don't want to use the same raw bin in two output bins
;		; If not on first bin, check whether start of current bin equals end of previous bin
;		; If so, then adjust current bin to start one later up, or previous bin to start one earlier,
;		; depending on which is closer to the raw edge
;		if i gt 0 then begin
;			if index[0,i] eq index[1,i-1] and index[0,i] ne -1 then begin
;				edge = edges[*,index[0,i]]
;				if abs(int[0,i] - edge[1]) lt abs(int[1,i-1] - edge[0]) then begin
;					index[0,i] = index[0,i] + 1  ; adjust current bin up by one
;				endif else begin
;					index[1,i-1] = index[1,i-1] - 1	; adjust previous bin down by one
;					; check that after adjustment it's still a valid bin
;					if index[0,i-1] gt index[1,i-1] then index[*,i-1] = [-1,-1]
;				endelse
;			endif
;		endif
;
;		if index[0,i] gt index[1,i] then index[*,i] = [-1,-1]  ; if start > end of bin, invalid bin
;	endfor

	; only use valid bins (no -1's)
	q = where (index[0,*] ne -1 and index[1,*] ne -1, count, complement=bad_interval)
	if count gt 0 then begin
		index = index[*,q]
		newedges = dblarr(2,n_elements(index[0,*]) )
		newedges[0,*] = edges[0,index[0,*]]
		newedges[1,*] = edges[1,index[1,*]]
	endif else begin
		index = -1
		newedges = -1
		;got_int = 0
	endelse

endif else begin
	; no binning required
	index[0,*] = indgen(n_int)
	index[1,*] = indgen(n_int)
endelse

return, index
end

;---------------------------------------------------------------------------

; convert_fitfunc_units method converts the calculated fit function array to counts if requested,
; (or stays in photons) and returns units of counts, rate, or flux.
; yvals is array of fit function values
; fit function units always start as photons/cm^2/sec/keV
; If output in counts is requested, yvals must be at photon energy bins to start with,
;   and will be returned in count energy bins.
; If output in photons is requested, yvals should be at count energy bins (since
;   converting to counts or rate uses ewidth and ltime (for counts), and they're at ct edges
; Output options:
;   photons = 0 or 1 for count or photon output
;   spex_units = 'counts', 'rate',  or 'flux'
; Default is photons=0, spex_units='counts'
; (spex_units='counts' and type= 'photons' means photon counts)
; Need interval number to get correct filter for SRM, and some conversions need width of time interval.
;   Use spex_interval_index unless this_interval is passed in


function spex_gen::convert_fitfunc_units, $
	yvals, $
	photons=photons, $
	spex_units=spex_units, $
	this_interval=this_interval, $
	use_fitted=use_fitted, $
	_extra=_extra

checkvar, photons, 0
checkvar, use_fitted, 0

; allow fuzzy units selection (anything close is OK)
units_choices = ['flux', 'rate', 'counts']
index = str_fuzzy (['fl', 'ra', 'co'], spex_units)
spex_units = index eq -1 ? 'counts' : units_choices[index]

checkvar, this_interval, self->get(/spex_interval_index)

filter = self -> get_fitint_filter(this_interval=this_interval, use_fitted=use_fitted)
time =   self -> get_fitint_time  (this_interval=this_interval, use_fitted=use_fitted)


if use_fitted then begin
	area = self -> get(/spex_summ_area, class='spex_fit')
	e_width = get_edge_products(self -> get(/spex_summ_energy,class='spex_fit'), /width)
endif else begin
	area = self -> get(/spex_drm_area, class='spex_drm')
	e_width = self -> getaxis(/ct_energy, /width)
endelse

; yvals starts as photons in flux units, so for that case don't do anything.

case photons of
	0: begin
		valid_drm = self -> valid_drm(use_fitted=use_fitted)
		if valid_drm eq 1 then begin
			; yvals should be at photon energy bins, apply_drm returns it in count energy bins
			func_obj = self->make_func_obj (use_fitted=use_fitted, this_interval=this_interval)
			data = (self->get(/obj,class='spex_drm')) -> $
				apply_drm(yvals, this_filter=filter, this_time=time, func_obj=func_obj, _extra=_extra)
			; apply_drm returns count/sec, so if rate requested, this is result
			destroy, func_obj
			case spex_units of
				'counts': begin
					fitdata = self -> framework::getdata(class='spex_fitint')
					if not is_struct(fitdata) then return, -1
					return,data * fitdata.ltime[*,this_interval]
					end
				'rate': return, data
				'flux': return, f_div(data, area * e_width)
			endcase

		endif else begin
			if use_fitted then begin
				; converts to count flux (not rate as above)
				data = yvals * (self -> get(/spex_summ_conv, class='spex_fit'))[*,this_interval]
				case spex_units of
					'counts': return, -1
					'rate': return, data * area * e_width
					'flux': return, data
				endcase
			endif else return, -1
		endelse
		end

	1: begin
		case spex_units of
			'counts': begin
				fitdata = self -> framework::getdata(class='spex_fitint')
				if not is_struct(fitdata) then return, -1
				return, yvals * area * e_width * fitdata.ltime[*,this_interval]
				end
			'rate':   return, yvals * area * e_width
			'flux':   return, yvals
		endcase
		end

endcase

end

;---------------------------------------------------------------------------

function spex_gen::get_erange, use_fitted=use_fitted, this_interval=this_interval

if use_fitted then begin
  eindex = where ( (self -> get (/spex_summ_emask,class='spex_fit'))[*,this_interval[0]], neind )
  if neind eq 0 then erange = -1 else $
    erange = find_contig_ranges( (self -> get(/spex_summ_energy,class='spex_fit'))[*,eindex], epsilon=1.e-5 )
endif else erange = self -> get(/spex_erange)

return, erange
end

;---------------------------------------------------------------------------
; get filter for this_interval fitint interval(s).  If this_interval not passed in
; use spex_interval_index.

function spex_gen::get_fitint_filter, this_interval=this_interval, use_fitted=use_fitted

@spex_insert_catch

checkvar, this_interval, self->get(/spex_interval_index)

if keyword_set(use_fitted) then begin
	;if /use_fitted, get fitint filters out of summary parameters.  If wasn't there,
	; will be null pointer, so set to -1
	fitint_filter = self -> get(/spex_summ_filter, class='spex_fit')
	if size(fitint_filter,/tname) eq 'POINTER' then fitint_filter=-1
	; make sure filters correspond to current definition of fit intervals by
	; doing a getdata.
endif else begin
	tmp = self -> framework::getdata(class='spex_fitint')

	fitint_filter = self -> get(/spex_fitint_filter)
endelse

filter = fitint_filter[0] eq -1 ? rebin([-1], n_elements(this_interval)) : fitint_filter[this_interval]

return, n_elements(filter) eq 1 ? filter[0] : filter
end

;---------------------------------------------------------------------------
; get time for 'this_interval' fitint interval(s).  If this_interval not passed in
; use spex_interval_index.

function spex_gen::get_fitint_time, this_interval=this_interval, use_fitted=use_fitted

@spex_insert_catch

checkvar, this_interval, self->get(/spex_interval_index)

if keyword_set(use_fitted) then begin
  ;if /use_fitted, get fitint times out of summary parameters.  If wasn't there,
  ; will be null pointer, so set to -1
  fitint_time = self -> get(/spex_summ_time_interval, class='spex_fit')
  if size(fitint_time,/tname) eq 'POINTER' then fitint_time=-1.d0
endif else begin
  ; make sure time correspond to current definition of fit intervals by
  ; doing a getdata.
  tmp = self -> framework::getdata(class='spex_fitint')

  fitint_time = self -> get(/spex_fit_time_interval )
endelse

time = fitint_time[0] eq -1.d0 ? rebin([-1.d0], 2, n_elements(this_interval)) : fitint_time[*,this_interval]

return, n_elements(time) eq 2 ? time[*,0] : time
end

;---------------------------------------------------------------------------

; Get efficiency factors used to convert from counts to photons.
;
; this_interval - scalar or vector of fit intervals to get eff. fact. for.
;   Default is spex_interval_index.
; use_fitted - if set, use fit parameters from spex_summ if available.  Default is 1.
;   If not set, use current fit_function and fit_comp_params.
; (NOTE: interval is important for getting the SRM for the right filter state or time, fit parameters
;   are important in computing diagonals.)


function spex_gen::get_drm_eff, $
	this_interval=this_interval, $
	use_fitted=use_fitted

checkvar, use_fitted, 1
this_use_fitted = use_fitted

checkvar, this_interval, self->get(/spex_interval_index)

valid_drm = self -> valid_drm(use_fitted=use_fitted, /message)
if valid_drm lt 0 then return, -1
if valid_drm eq 2 then this_use_fitted = 0

filter = self -> get_fitint_filter(this_interval=this_interval)
time =   self -> get_fitint_time  (this_interval=this_interval)

;func_def = self -> get(/fit_function, class='fit_function')
;params_def = self -> get(/fit_comp_params, class='fit_comp_manager')	;???? want class=?

drm_obj = self -> get(/obj, class='spex_drm')

valid = this_use_fitted ? $
	self -> valid_summ_params(interval=this_interval, /check_match) : $
	rebin([1], n_elements(this_interval))

for i = 0,n_elements(this_interval)-1 do begin
	int = this_interval[i]
	fitobj = self->make_func_obj (use_fitted=(this_use_fitted and valid[i]), $
			this_interval=int, status=status)
	eff = drm_obj -> get_diagonal(fitobj, this_filter=filter[i], this_time=time[*,i])
	drm_eff = exist(drm_eff) ? [[drm_eff], [eff]] : eff
	destroy, fitobj
endfor

;if exist(this_interval) then begin
;	if this_use_fitted then begin
;;		summ = self -> get(/spex_summ, class='spex_fit')		;????? want class=?
;		valid = self -> valid_summ_params(interval=this_interval, /check_match)
;	endif
;	for i = 0, n_elements(this_interval)-1 do begin
;;		int = this_interval[i]
;		this_filter = filter[i]
;		fitobj = self->make_func_obj (use_fitted=(this_use_fitted and valid[i]), $
;			this_interval=this_interval[i], status=status)
;;		params = params_def
;;		func = func_def
;;		if this_use_fitted then if valid[i] then begin
;;			params = summ.spex_summ_params[*,int]
;;			func = summ.spex_summ_fit_function
;;		endif
;;		eff = drm_obj -> get_diagonal(func, params, this_filter=this_filter)
;		eff = drm_obj -> get_diagonal(fitobj, this_filter=this_filter)
;		drm_eff = exist(drm_eff) ? [[drm_eff], [eff]] : eff
;	endfor
;;endif else drm_eff = drm_obj -> get_diagonal(func_def, params_def, this_filter=filter[0])
;endif else drm_eff = drm_obj -> get_diagonal(fitobj, this_filter=filter[0])

return, drm_eff

end

;---------------------------------------------------------------------------

; check whether data has been saved in the spex_summ... structure.  If interval is set,
; then make sure that interval has data in it.
; if check_match is set, then check that times and ct energies match.
; If interval is an array, return an array of 1's and 0's, otherwise scalar 1 or 0.

function spex_gen::valid_summ_params, summ=summ, interval=interval, check_match=check_match, message=message

checkvar, check_match, 0
checkvar, message, 0

if not keyword_set(summ) then summ = self -> get(/spex_summ)		;????? want class=?

; initialize return value to 0, or an array of 0's for each interval
retval = exist(interval) ? interval*0 : 0

if size(summ, /tname) ne 'STRUCT' then return, retval

if size(summ.spex_summ_fit_done, /tname) eq 'POINTER' then return, retval

check_energy = self -> get(/spex_allow_diff_energy) eq 0

if check_match then begin
	fit_time = self -> get(/spex_fit_time_used)
	fit_range = minmax(fit_time)
	summ_range = minmax(summ.spex_summ_time_interval)
	if (fit_range[1] lt summ_range[0]) or (fit_range[0] gt summ_range[1]) then return, retval*0
	ct_energy = self -> getaxis(/ct_energy, /edges_2)
endif

if exist(interval) then begin
	for i=0,n_elements(interval)-1 do begin
		good = 1
		int = interval[i]
		if int ge 0 and int lt n_elements(summ.spex_summ_fit_done) then begin
			if not summ.spex_summ_fit_done[int] then good=0
		endif else good=0
		if good and check_match then begin
			if n_elements(fit_time[0,*])-1 ge int then begin	; corrected from i to int 21-mar-2008
				if max(abs(fit_time[*,int] - summ.spex_summ_time_interval[*,int])) gt 1. then good=0
				if good and check_energy then begin
					if not same_data(ct_energy, summ.spex_summ_energy) then good=0
					if not good and not spex_get_nointeractive() then begin
					  ncurr = n_elements(ct_energy[0,*]) & nsumm = n_elements(summ.spex_summ_energy[0,*])
					  mmc = minmax(ct_energy) & mms = minmax(summ.spex_summ_energy)
						msg = ['Current energy edges are not the same as the energy edges in the spex_summ', $
							'structure so the saved fit results may not be valid for the current input data.  ', $
							'', $
							'Current # bins   = ' + trim(ncurr) + '  Min,max = ' + arr2str(trim(mmc)) + ' keV', $
							'spex_summ # bins = ' + trim(ncurr) + '  Min,max = ' + arr2str(trim(mms)) + ' keV', $
							'', $
							'If they are slightly different (e.g. RHESSI native energy bins ', $ 
							'are slightly different for each detector), the spex_summ parameters ', $
							'are probably valid to use.  If they are really different, you may have ', $
							'made an error in file selection.', $
							'', $
							'Do you want to proceed? ', $
							'', $
							'If you answer YES, the spex_allow_diff_energy flag will be ', $
							"set to 1 and you won't be asked again.  The stored fit results will be used.", $
							"If you answer NO, the stored fit results won't be used - the current fit parameters will be used."]

						yesno = dialog_message(msg, /question)
						if yesno eq 'Yes' then begin
							self -> set,/spex_allow_diff_energy
							good = 1
						endif
					endif

				endif
			endif

		endif
		if not good and message then $
			message,/cont,'Data times or energy for interval ' + $
			trim(int) + ' do not match time or energy in saved fit results. Not using saved.'
		retval[i] = good
	endfor
endif else retval = 1

return, n_elements(retval) eq 1 ? retval[0] : retval
end

;---------------------------------------------------------------------------

function spex_gen::get_plot_title, photons=photons

units_str = self -> getunits()
yes_photons = (keyword_set(photons) or stregex(units_str.data,'photon', /boolean, /fold_case) )
space = yes_photons ? ' Photon' : ' Count'
case units_str.data_type of
    'counts': title = 'SPEX ' + units_str.data_name + space + 's'
    'rate': title = 'SPEX ' + units_str.data_name + space + ' Rate'
    else: title = 'SPEX ' + units_str.data_name + space + ' Flux'
endcase
return, title
end

;---------------------------------------------------------------------------

function spex_gen::plot_setup, class_name=class_name, plot_params=plot_params, _extra=_extra

; if object doesn't have a plot_setup method, the code below will call this function
; again, so just return an error.
if get_caller() eq 'SPEX_GEN::PLOT_SETUP' then begin
	message, 'This class does not have a plot method: '+obj_class(self),/cont
	return, 0
endif

if keyword_set(class_name) then begin
	obj = self -> get(/obj, class_name=class_name)
	return, obj -> plot_setup (plot_params=plot_params, _extra=_extra)

endif else return, self -> plot_setup (plot_params=plot_params, _extra=_extra)

end

;---------------------------------------------------------------------------

;pro spex_gen::plot_data, $
;	bksub=bksub, $
;	overlay_bk=overlay_bk, $
;	photons=photons, $
;	_extra=_extra
;
;fit_times = self -> get(/spex_fit_time_interval)
;
;if keyword_set(photons) then drm_eff = (self -> get(/obj, class='spex_drm'))-> get_diagonal()
;
;class = keyword_set(bksub) ? 'spex_bksub' : 'spex_data'
;status = self -> plot_setup (class=class, plot_params=plot_params, /pl_energy, intervals=fit_times, photons=photons, $
;	drm_eff=drm_eff, _extra=_extra)
;
;if status then begin
;
;	if keyword_set(overlay_bk) then begin
;		status = self -> plot_setup (class='spex_bk', plot_params=plot_params_bk, /pl_energy, intervals=fit_times, photons=photons, $
;			drm_eff=drm_eff, _extra=_extra)
;		if status then begin
;			plot_params = rep_tag_value(plot_params, [ [plot_params.ydata], [plot_params_bk.ydata]], 'ydata')
;			plot_params = rep_tag_value(plot_params, [ plot_params.dim1_id, plot_params_bk.dim1_id], 'dim1_id')
;			plot_params = rep_tag_value(plot_params, [ [plot_params.dim1_vals], [plot_params_bk.dim1_vals]], 'dim1_vals')
;			plot_params.dim1_sum = 0
;		endif
;	endif
;
;	self -> do_plot, /plotman, plot_params=plot_params, _extra=_extra
;
;endif
;
;end

;---------------------------------------------------------------------------

;_extra can include spex_units, comb_func, sep_func, show_err
;
; bksub - if set, use background-subtracted data, otherwise use data including background for primary plot
; overlay_bksub - if set, overlay background-subtracted data
; overlay_back - if set, overlay background
; photons - if set, plot is in photons
; show_fit - if set, show fit function
; use_fitted - if set, then if show_fit is set, show fit function from spex_summ fitted results.
;	Otherwise, show current fit parameters stored in fit_comp structure.  Default is 1.
; origint - If set, use time intervals in original input file.  Otherwise use spex_fit_time_interval.
; tband - If set, use spex_tband time intervals binning for plot (default is the first one, unless allint
;   or this_interval is set).  If tband is set, but there are no tband intervals, reverts to origint.
; this_interval - scalar or vector of interval indices to plot.  If not specified, then
;	for origint=1, use interval 0, otherwise use current spex_interval_index.
; allint - If set, override this_interval, and select all intervals (either original or fit ints)
; no_unique_label - If set, don't label plot with Interval #.  Since the label will then
;	be generic, panels will be replaced instead of preserved for every plot.
; no_plotman - If set, plot in regular window, not plotman window
; if obj keyword is present, just pass out the xy or utplot object.  Don't draw the plot.
;
; Default is to plot on fitint time boundaries, for the current value of spex_interval_index
;
; Since we allow time binning on original file interval boundaries, fit interval boundaries,
; or tband boundaries,
; we'll use the bksub object here for background-subtracted data and bin it to the fit intervals
; if necessary, instead of the fitint object since the fitint object is already binned
; to fit interval boundaries.

pro spex_gen::plot_spectrum, $
	bksub=bksub, $
	overlay_bksub=overlay_bksub, $
	overlay_back=overlay_back, $
	photons=photons, $
	show_fit=show_fit, $
	use_fitted=use_fitted, $
	origint=origint, $
	this_interval=this_interval, $
	allint=allint, $
	tband=tband, $
	no_unique_label=no_unique_label, $
	no_plotman=no_plotman, $
	obj=obj, $
	_extra=_extra

error_catch = 0
if spex_get_debug() eq 0 then catch, error_catch
if error_catch ne 0 then begin

  catch, /cancel

  ; was exception generated by a MESSAGE call?
  if !error_state.name eq 'IDL_M_USER_ERR' then print,!error_state.msg else begin $
    help, /last_message, out=out
    prstr,out, /nomore
  endelse
  return
endif

checkvar, bksub, 0
checkvar, overlay_bksub, 0
if bksub then overlay_bksub=0
checkvar, overlay_back, 0
checkvar, photons, 0
checkvar, show_fit, 0
checkvar, use_fitted, 1
checkvar, origint, 0
checkvar, tband, 0
use_plotman = keyword_set(no_plotman) eq 0
phot = photons	; use local variable, so we don't change returned value

tmp = self -> framework::getdata(class='spex_data') ; just to make sure times, etc are read in
if not is_struct(tmp) then return

int_times = -1

; if tband is selected, use spex_tband, otherwise original
if tband then begin
	checkvar, this_interval, 0
	int_times = self -> get(/spex_tband)
	if int_times[0] eq -1 then origint=1
	show_fit = 0
	phot = 0
endif

if origint then begin
	checkvar, this_interval, 0
	int_times = self -> getaxis(/ut, /edges_2)
	show_fit = 0	; added 20-jul-04 - no filter info if not using fitint intervals
	phot = 0	; added 20-jul-04
endif

if int_times[0] eq -1 then begin
	tmp = self -> getdata(class='spex_fitint') ; just to make sure spex_fit_time_used is defined
	int_times = self -> get(/spex_fit_time_used)
	checkvar, this_interval, self->get(/spex_interval_index)

	if int_times[0] eq -1 then begin
		int_times = self -> get(/spex_file_time)
		this_interval=0
	endif
endif

if keyword_set(allint) then begin
	nt = n_elements(int_times)/2
	if nt gt 100 then begin
		this_interval = 0
		int_times = minmax(int_times)
		print,'Too many intervals ('+trim(nt)+') to plot individually.  Set spex_tband to group time bands.'
	endif else this_interval = findgen(nt)
endif

this_interval = fix(get_uniq(this_interval > 0 < (n_elements(int_times)/2 - 1)))
int_times = int_times[*,this_interval]

; first make sure we can use the options (photons, use_fitted, show_fit) requested.

if phot or show_fit then begin

	if use_fitted then begin
		; if using fitted info, first check if time intervals need to be adjusted in spex_summ
		; to match current intervals, otherwise interval # is wrong.
		spex_summ = (self -> get(/obj,class='spex_fit')) -> check_times()
		valid_summ = self -> valid_summ_params(interval=this_interval, /check_match, /message)
		if min(valid_summ) eq 0 then use_fitted = 0
	endif

	valid_drm = self -> valid_drm (use_fitted=use_fitted, /message)

	if phot then begin
		drm_eff = -1
		if valid_drm gt 0 then $
			drm_eff = self -> get_drm_eff(this_interval=this_interval, use_fitted=use_fitted)
		if drm_eff[0] eq -1 then begin
			phot = 0
			message,/cont, 'Not plotting in photons because drm not valid.'
		endif
	endif

	if show_fit and not phot and valid_drm lt 0 then begin
		show_fit = 0
		message,/cont, 'Not plotting fit data in counts because drm not valid.'
	endif

	if show_fit and use_fitted then begin
		valid_summ = self -> valid_summ_params(interval=this_interval, /check_match)
		if min(valid_summ) eq 0 then begin
			show_fit = 0
			message,/cont, "Not plotting saved fit data because it doesn't match current data."
		endif
	endif

endif

class = bksub ? 'spex_bksub' : 'spex_data'
status = self -> plot_setup (class=class, plot_params=plot_params, /pl_energy, intervals=int_times, $
	photons=phot, drm_eff=drm_eff, dim1_sum=0, _extra=_extra)

if status then begin

	if keyword_set(overlay_back) then begin
		status = self -> plot_setup (class='spex_bk', plot_params=plot_params2, /pl_energy, $
			intervals=int_times, photons=phot, drm_eff=drm_eff, _extra=_extra)
		if status then plot_params = spex_merge_plots (plot_params, plot_params2)
	endif

	if keyword_set(overlay_bksub) then begin
		status = self -> plot_setup (class='spex_bksub', plot_params=plot_params2, /pl_energy, intervals=int_times, photons=phot, $
			drm_eff=drm_eff, _extra=_extra)
		if status then plot_params = spex_merge_plots (plot_params, plot_params2)
	endif

	if show_fit then begin

		fitplot_obj = self -> make_fitplot_obj (/overlay, photons=phot, $
			this_interval=this_interval, use_fitted=use_fitted, $
			no_plotman=no_plotman, _extra=_extra)
		plot_params = add_tag (plot_params, fitplot_obj, 'overlay_obj')
		plot_params.id = keyword_set(no_unique_label) ? plot_params.id + ' with Fit Function' : $
			plot_params.id + ' with Fit Function, Interval ' + trim(this_interval[0])

	endif

	plot_params = add_tag (plot_params, 2, 'legend_loc')	;put legend in upper right corner

	; if we're doing an image cube (roi_use ne -1) then add regions label to label
	spex_roi_use = self->get(/spex_roi_use)
	if spex_roi_use[0] ne -1 then begin
		roi = (self->get(/spex_roi_integrate) eq 1) ? '  All Regions' : $
			'  Region ' + arr2str(trim(spex_roi_use))
		if phot and (self->get(/spex_image_full_srm)) then roi = [roi, 'Full SRM used']
		spec_source = self->get(/spex_image_spectrum_source)
		if spec_source ne '' then roi = [roi, spec_source + ' Spectrum']
		new_label = tag_exist(plot_params,'label') ? [roi,plot_params.label] : roi
		plot_params = rep_tag_value(plot_params, new_label, 'label')
	endif

;	If drawing lines for energy ranges, for plotman, set up the commands  in cmds and add as more_commands
;	parameter in plotman.  If not using plotman, draw the lines after making the plot.
	if show_fit then begin
		erange = self -> get_erange (use_fitted=use_fitted, this_interval=this_interval)
    if erange[0] ne -1 then begin	
       plot_params = add_tag (plot_params, 'oplot_xmarkers', 'addplot_name')
       plot_params = add_tag (plot_params, {intervals: erange}, 'addplot_arg')
    endif
	endif

	self -> do_plot, use_plotman=use_plotman, plot_params=plot_params, obj=obj, _extra=_extra
endif

end

;---------------------------------------------------------------------------

; plot_spectrum_resid method plots the spectrum for the interval(s) selected with whatever options are in _extra,
; and plots the residuals for the same interval(s), and then overplots (stacks) the residuals on the spectrum plot
; with a 70%, 30% split
; Arguments:
; this_interval - interval number(s) to plot
; _extra - any of the keywords for plot_spectrum, e.g. /bksub, spex_units=' ', /show_fit, /show_err, /photons, etc.

pro spex_gen::plot_spectrum_resid, this_interval=interval,  _extra=_extra

checkvar, interval, 0

pobj = self->get_plotman_obj()

(self -> get(/obj,class='spex_fit')) -> plot_resid, interval=interval, _extra=_extra, /main_window, $
  legend_loc=0, ytitle= ' ', status=status
if ~status then return

resid_panel = pobj->get(/current_panel_desc)

; Removed legend_loc=0 from plot_spectrum call, request from Brian 11-Apr-2017
self->plot_spectrum,  this_interval=interval, _extra=_extra

pobj->select
overlay_panel = pobj->get(/overlay_panel)
overlay_panel[1] = resid_panel
overlay_ysize = pobj->get(/overlay_ysize)
overlay_ysize[0] = 70
overlay_ysize[1] = 30
pobj->set,overlay_panel=overlay_panel, overlay_ysize=overlay_ysize, /overlay_squish
pobj->plot

end

;---------------------------------------------------------------------------

; General time profile plotting method.
; Up to 3 types of data selected are overlaid. If no data type is selected, then
; 'data' (data with background) is used.  Grouped by spex_ebands energy bands unless
; origint or bkint is set.
; data - If set, then show data (with background)
; bksub - If set, then show background-subtracted data
; bkdata - If set, then show background
; origint - If set, use original energy intervals (see this_interval)
; bkint - If set, use spex_bk_eband energy intervals (see this_interval)
; this_interval - scalar or vector, plot only those energy intervals
; this_band - same as this_interval (added for consistency with get and set methods)
; no_plotman - If set, plot to regular plot window.
; if obj keyword is present, just pass out the xy or utplot object.  Don't draw the plot.
; show_filter - If set, display filter states at top of plot
; _extra can contain spex_units

; Examples:
; o -> plot_time, /data, /bkdata, spex_units='flux', this_interval=2
; o -> plot_time, /bksub, spex_units='rate', /origint, this_interval=indgen(10)
; o -> plot_time, /data, /bkdata, /bkint, this_interval=1

pro spex_gen::plot_time, $
	data=data, $
	bksub=bksub, $
	bkdata=bkdata, $
	origint=origint, $
	bkint=bkint, $
	this_interval=this_interval, $
	this_band=this_band, $  ; same as this_interval
	no_plotman=no_plotman, $
	show_filter=show_filter, $
	obj=obj, $
	_extra=_extra

error_catch = 0
if spex_get_debug() eq 0 then catch, error_catch
if error_catch ne 0 then begin

  catch, /cancel

  ; was exception generated by a MESSAGE call?
  if !error_state.name eq 'IDL_M_USER_ERR' then print,!error_state.msg else begin $
    help, /last_message, out=out
    prstr,out, /nomore
  endelse
  return
endif

checkvar, origint, 0
checkvar, bkint, 0
use_plotman = keyword_set(no_plotman) eq 0

if exist(this_band) then this_interval=this_band

checkvar, data, 0
checkvar, bksub, 0
checkvar, bkdata, 0
if data+bksub+bkdata eq 0 then data = 1

if data then class = 'spex_data' else begin
	if bksub then class = 'spex_bksub' else class = 'spex_bk'
endelse

data = self -> framework::getdata(class='spex_data') ; just to make sure times, energies are read in
if not is_struct(data) then return

if origint then begin
	checkvar, this_interval, 0
	energy_int = self -> getaxis(/ct_energy, /edges_2)
endif else begin
	if bkint then begin
		energy_int = self->get(/spex_bk_sep) ? self -> get(/spex_bk_eband) : -1
	endif else energy_int = self -> get(/spex_eband)
	if energy_int[0] eq -1 then begin
		energy_int = minmax(self -> getaxis(/ct_energy, /edges_2))
		this_interval=0
	endif else checkvar, this_interval, indgen(n_elements(energy_int)/2)
endelse

this_interval = get_uniq(this_interval > 0 < (n_elements(energy_int)/2 - 1))
energy_int = energy_int[*,this_interval]

; call spex_gen::plot_setup instead of self->plot_setup so that if class is passed it will get used
; 11-apr-2005, Kim.
status = self -> spex_gen::plot_setup (class=class, plot_params=plot_params, intervals=energy_int, $
	photons=phot, drm_eff=drm_eff, dim1_sum=0, _extra=_extra)

if status then begin

	if bksub and class ne 'spex_bksub' then begin
		status = self -> plot_setup (class='spex_bksub', plot_params=plot_params2, $
			intervals=energy_int, photons=phot, drm_eff=drm_eff, _extra=_extra)
		if status then plot_params = spex_merge_plots (plot_params, plot_params2)
	endif

	if bkdata and class ne 'spex_bk' then begin
		status = self -> plot_setup (class='spex_bk', plot_params=plot_params2, $
			intervals=energy_int, photons=phot, drm_eff=drm_eff, _extra=_extra)
		if status then plot_params = spex_merge_plots (plot_params, plot_params2)
	endif

	; if we're doing an image cube (roi_use ne -1) then add regions label to label
	spex_roi_use = self->get(/spex_roi_use)
	if spex_roi_use[0] ne -1 then begin
		roi = (self->get(/spex_roi_integrate) eq 1) ? '  All Regions' : $
			'  Region ' + arr2str(trim(spex_roi_use))
		new_label = tag_exist(plot_params,'label') ? [roi,plot_params.label] : roi
		plot_params = rep_tag_value(plot_params, new_label, 'label')
	endif

	if keyword_set(show_filter) then begin
		str = self -> get_filter_bars()
		if is_struct(str) then begin
			;if there is filter info, move label down, and add plot routine name and structure
			; to plot_params structure.
			nuniq = n_elements(get_uniq(str.label))
			plot_params = rep_tag_value (plot_params, [replicate('',nuniq), [plot_params.label]], 'label')
			plot_params = add_tag(plot_params, 'spex_draw_bars', 'addplot_name')
			plot_params = add_tag(plot_params, str, 'addplot_arg')
		endif
	endif

	self -> do_plot, use_plotman=use_plotman, plot_params=plot_params, obj=obj, _extra=_extra

; changed drawing the filter bars to a routine, spex_draw_bars, that is automatically called
; by xyplot obj
;	if keyword_set(show_filter) then begin
;		if use_plotman then begin
;			p = self -> get_plotman_obj(valid=valid)
;			if valid then p->select else use_plotman=0
;		endif
;
;		self -> filter_bars, str, /plot
;		if use_plotman then p -> unselect
;	endif

endif
end

;---------------------------------------------------------------------------

; check that drm is valid by checking if :
; 1. drm has count energy edges
; 2. drm count energy edges equals current spex_ct_edges (or if use_fitted set, spex_summ_energy)
; 3. drm detectors equals data detectors (or is use_fitted set, spex_summ_detused)
;
; returns 1 if everything OK.  returns 2 if use_fitted was set, and fitted energies don't match, but
;	current edges do match (means we can't use this drm for the saved fitted data)
; returns 0 otherwise
;
; if message is set, prints the reason for a return of other than 1.

function spex_gen::valid_drm, use_fitted=use_fitted, message=message

checkvar, message, 0
checkvar, use_fitted, 0

valid = 1

drm = self -> framework::getdata(class='spex_drm') ; just to make sure it's read in and initialized
if drm[0] eq -1 then code = 0 else begin

	drm_ct_edges = self->get(/spex_drm_ct_edges)
	if drm_ct_edges[0] eq -1 then code = 1 else begin

		if use_fitted then begin
			fit_edges = self -> get(/spex_summ_energy, class='spex_fit')
		;	dets = self -> get(/spex_summ_detused)
			if compare_float(fit_edges, drm_ct_edges) then return, 1
			code = 2
			valid = 2
		endif else code = 3
		ct_edges = self->get(/spex_ct_edges)
		dets = self -> get(/spex_detectors)
		if compare_float(ct_edges, drm_ct_edges) then begin
			if use_fitted and message then $
				message,/cont,'DRM is not valid for saved fit data, but is valid for current data.'
			return, valid
		endif
	endelse
endelse

if message then begin
	msg = ['No data in drm.', $
			'No count energy edges in drm.', $
			'Saved Fit energy edges and drm energy edges do not match.', $
			'Current energy edges and drm energy edges do not match.' ]
	message, /cont, 'DRM is not valid because:  ' + msg[code]
endif

return, 0
end

;---------------------------------------------------------------------------
; calc_func_components - Function to calculate and return function components in any units, 
; in photon or count space, separately or combined.  This code was part of make_fitplot_obj, 
; which now calls this function (23-Jun-2008)
; Arguments: 
; photons - if set, function is in photon space (otherwise count space)
; spex_units - 'counts','rate','flux' units for function (default is 'counts')
; all_func - if set, return all component functions separately and combined
; comb_func - if set just return combined components
; sep_func - if set, just return separate components
; Note: if /all_func isn't set, then looks in _extra for name of function to return, like /vth.  If
;   no component is requested in _extra, then returns all separately and combined unless
;   /comb_func or /sep_func is set.
; use_fitted - if set, use the fit parameters saved in spex_summ structure for calculating function.  If 
;   not set, then use current parameters in the fit_comp structure.  Default is 1 if this_interval is set.
; chisq - returns value of chisq from fit results (or in some cases, when fit hasn't been saved yet, chisq is
;   already set and isn't changed here)
; this_interval - scalar or array of interval #'s to use (only makes senses if use_fitted is set)
; ph_edges - if set, calculates function on photon edges.  Only works for /photons and spex_units='flux'
;
; Returns a structure with 
;   yvals - function values, dimensioned (nenergy, nfunc_component, ninterval)
;   id - string label for each function, dimensioned nfunc_component (includes param values if nintervals=1)
;   strat_arr - string of function component strategy name, dimensioned nfunc_component
;   Suppose if we have 77 energies, function ='vth+bpow', 3 time intervals, and we call
;    Example 1:
;     a=o->calc_func_components(spex_units='flux', /all_func, this_int=[0,1],/photon, chisq=ch) , we get:
;     YVALS           FLOAT     Array[77, 3, 2]
;     ID              STRING    Array[3]
;     STRAT_ARR       STRING    Array[3]
;     The function components are returned for the combined function, and then for each separate component.
;     so a.id = ['vth+bpow', 'vth', 'bpow']
;     and a.yvals[*,0,0] are function values for the combined function, 0th interval
;           a.yvals[*,1,0] are function values for vth component, 0th interval, etc.
;    Example 2:
;     a=o->calc_func_components(spex_units='flux', /all_func, this_int=1,/photon, chisq=ch)
;     YVALS           FLOAT     Array[77, 3] ; 2nd dimension is for the combined (vth+bpow) and separate vth and bpow components
;     ID              STRING    Array[3]  ; ['vth+bpow', 'vth 1.00e-005,1.00,1.00  full mewe 3.92e-005', 'bpow 0.857,3.27,55.3,4.36'] 
;     STRAT_ARR       STRING    Array[3]  ; ['', 'vth#0', 'bpow#0']
;
function spex_gen::calc_func_components, $
  photons=photons, $
  spex_units=spex_units, $
  all_func=all_func, $
  comb_func=comb_func, $
  sep_func=sep_func, $
  use_fitted=use_fitted, $
  chisq=chisq, $
  this_interval=this_interval, $
  ph_edges=ph_edges, $
  _extra=_extra
  
photons = keyword_set(photons)
checkvar, use_fitted, exist(this_interval)
checkvar, spex_units, 'counts'
all_func = keyword_set(all_func)
checkvar, this_interval, self->get(/spex_interval_index)
ph_edges = keyword_set(ph_edges)

; don't need to check for valid_drm if we want photons
valid_drm = photons ? 0 : self->valid_drm(use_fitted=use_fitted)
if use_fitted then begin
  valid = self -> valid_summ_params(interval=this_interval)
  if total(valid) eq 0 then return, -1
  if valid_drm ne 1 and not photons then begin
    all_func=1
    comb_func=1
  endif
  energy = self->get(/spex_summ_energy)
endif else begin
  if not valid_drm then if not photons then return, -1
  valid=rebin([1], n_elements(this_interval))
  energy = self->getaxis(/ct_energy, /edges_2)
endelse

; if want fit function in photons, we won't apply drm, can use count edges here.
; if want fit function in counts, need to apply drm, so need to use drm_ph_edges here
if ph_edges then begin
  ; If user requested use_ph_edges, then can't return 'counts' units, have to return 'photons'.
  ; If returning photons, in ph bins, set energy to the ph edges
  if not photons or spex_units ne 'flux' then begin
    message,'Can only return func in photon bins in photon space flux units.', /cont
    return, -1
  endif 
  energy = self->getaxis(/ph_energy, /edges_2)  ; this is now really photon energy bins
  use_ph_edges = 1
endif else use_ph_edges = (valid_drm eq 1 and not photons)

; make new fit_function object so we don't change the one in the object chain.
; if use_fitted, make a fresh fit_function object.  If not, clone current one, in case
; strategy_name is an argument here.  For example,
;    if user had deleted bpow#0 and asks for strategy_name=bpow#1, in a fresh object, the first
;    bpow would be called bpow#0 and there would be no bpow#1.
; first create a dummy fitobj just so we can get function component information out of it.  Then
; while looping through intervals, create fitobj for the time intervals requested
; use this_interval in this call even though not needed because otherwise make_func_obj will try
; to use spex_interval_index which may not be valid for current spex_summ params. 5-mar-2008
fitobj = self->make_func_obj(use_fitted=use_fitted, use_ph_edges=use_ph_edges, $
        this_interval=this_interval[0], status=status)
if not status then return, -1
fit_function = fitobj->get(/fit_function)

; strat_arr will be name of fit components for each trace.  '' means all combined.  If doing all components
; separately and combined, prepend '' to list of separate component names, so combined will be first.
strat_arr = ''
if not all_func then begin
  if keyword_set(_extra) then strat_arr = fitobj -> name2strat(_extra=_extra)
endif
if strat_arr[0] eq '' then begin  ; no component requested in _extra keywords, so do all
  all_func = 1
  if keyword_set(comb_func) then strat_arr = '' else begin
    compman_strat = fitobj -> get(/compman_strategy)
    if keyword_set(sep_func) then strat_arr = compman_strat else $
      strat_arr = n_elements(compman_strat) eq 1 ? compman_strat[0] : ['', compman_strat]
  endelse
endif
obj_destroy,fitobj


; Set up y array to hold fit function values for each energy bin, for each function component, for each interval.
; (number of energy bins will always be number of ct bins)
; id will hold name of each function component for label.  id will contain the parameter values if only one
; interval is requested (and if separate and combined is requested, param values are only on the id's for separate)

n_plot_energy = n_elements(energy[0,*])
n_comp = n_elements(strat_arr)  ; number of components to plot
ninterval = n_elements(this_interval)
yvals = dblarr(n_plot_energy, n_comp, ninterval)
id = strarr(n_comp)

for int = 0,ninterval-1 do begin
  if valid[int] then begin
    fitobj = self->make_func_obj(use_fitted=use_fitted, use_ph_edges=use_ph_edges, $
      this_interval=this_interval[int], chisq=chisq, status=status)

    for i = 0,n_comp-1 do begin
      
      ; for albedo and pileup components, previously returned 0s (since they're pseudo functions). Now
      ; we calculate the contribution from albedo and pileup by taking difference in spectrum with and 
      ; without that correction, so we can plot those contributions.
      ; 
      ; Don't apply albedo or pileup correction to the other individual components - only the composite spectrum
      ; NOTE:  if change order of plotting,need to change this - right now, composite function is the 0th index in yvals.
      ; 
      ; Note that albedo can be enabled either as a function component, or through the spex_albedo_correct
      ; parameter.  In either case, only apply albedo to composite spectrum.  Only calculate albedo contribution
      ; for plotting if albedo was included as a function component.
      
      albedo = strpos(strat_arr[i], 'albedo') ne -1  ; will be 1 when we're doing the albedo component
      pileup = strpos(strat_arr[i], 'pileup_mod') ne -1  ; will be 1 when we're doing the pileup component
 
      if ph_edges then begin
      
        ; must be photon output, and must be in 'flux' units
        case 1 of
        albedo: begin
          ; doing separate albedo component
          yv = strat_arr[0] eq '' ? yvals[*,0,int] : fitobj -> getdata(strategy_name='')
          alb_params = fitobj -> get(/fit_comp_params,comp_name='albedo')  
          anis=alb_params[0]
          theta = self -> get(/spex_source_angle)
          yvfull = (self->get(/obj,class='spex_drm')) -> albedo(yv, theta=theta, anisotropy=anis, energy=energy)
          ynew = yvfull - yv 
          end
        else: ynew = fitobj -> getdata(strategy_name=strat_arr[i])
        endcase
        
      endif else begin
      
        ; can be photons or counts, and units can be counts,rate,or flux - convert_fitfunc_units will take care of both.
        ; Note that pileup and albedo correction are normally applied in apply_drm
        case 1 of
        pileup: begin
          ; doing separate pileup component
          ; convert yv flux to requested units and space with and without pileup correction, and take difference 
          yv = fitobj -> getdata(strategy_name='')
          yfull = strat_arr[0] eq '' ? yvals[*,0,int] : self -> convert_fitfunc_units (yv, photons=photons, $
            spex_units=spex_units, this_interval=this_interval[int], use_fitted=use_fitted, this_strat=strat_arr[i])
          yfull_nop = self -> convert_fitfunc_units (yv, photons=photons, spex_units=spex_units, $
            this_interval=this_interval[int], use_fitted=use_fitted, /disable_pileup_mod);, this_strat=strat_arr[i]) ;removed 8/24/2022, kim
          ynew = yfull - yfull_nop
          end
        albedo: begin
          ; doing separate albedo component. 
          yv = fitobj -> getdata(strategy_name='')
          if photons then begin
            ; for photon output, call albedo method directly since need to supply the energy edges (yv here is on
            ; count edges, but default for albedo matrix is photon edges). Take difference of flux with albedo
            ; correction (yvfull) to without (yv) 
            alb_params = fitobj -> get(/fit_comp_params,comp_name='albedo')  
            anis=alb_params[0]
            theta = self -> get(/spex_source_angle)
            yvfull = (self->get(/obj,class='spex_drm')) -> albedo(yv, theta=theta, anisotropy=anis, energy=energy)
            yv_alb = yvfull - yv
            ynew = self -> convert_fitfunc_units(yv_alb, photons=photons, spex_units=spex_units, $
              this_interval=this_interval[int], use_fitted=use_fitted, this_strat=strat_arr[i])
          endif else begin
            ; if not photon output, then yv was calculated on photon edges.  Can let convert_fitfunc_untis call
            ; apply_drm as normal, and it will call albedo method.  Make sure pileup correction is not done on 
            ; both the spectrum with and without albedo            
            yfull_nop = self -> convert_fitfunc_units (yv, photons=photons, spex_units=spex_units, $
              this_interval=this_interval[int], use_fitted=use_fitted, /disable_pileup_mod, this_strat=strat_arr[i])
            yfull_nop_noa = self -> convert_fitfunc_units (yv, photons=photons, spex_units=spex_units, $
              this_interval=this_interval[int], use_fitted=use_fitted, /disable_pileup_mod, /disable_albedo, this_strat=strat_arr[i])
            ynew = yfull_nop - yfull_nop_noa
          endelse
          end
        else:begin
          yv = fitobj -> getdata(strategy_name=strat_arr[i])
          ynew = self -> convert_fitfunc_units (yv, photons=photons, spex_units=spex_units, $
            this_interval=this_interval[int], use_fitted=use_fitted, $
            disable_pileup_mod=(i gt 0), disable_albedo=(i gt 0), this_strat=strat_arr[i])
          end
        endcase
      endelse 
      
      yvals[*,i,int] = ynew[0] eq -1 ? !values.f_nan : ynew
      comp_params = fitobj -> getfit(/fit_comp_params, strategy_name=strat_arr[i])
      ; if only one interval, then can put params and kw plot label in label
      str_params = ninterval eq 1 ? ' ' + arr2str (trim(comp_params, '(g9.3)'), ',') : ''
      pl = ''
      ; only get keyword plot labels if doing a single interval
      if ninterval eq 1 then kw = $
        fit_comp_kw(fitobj->getfit(/fit_comp, strategy_name=strat_arr[i]), plot_label=pl, obj=self)
      ;strat_arr[i] = '' is the sum of all components.  If there's only one, show the parameters
      ; next to it, otherwise, show the parameters next to the separate component they apply to.
      if strat_arr[i] eq '' then begin
        id[i] = n_comp eq 1 ? fit_function + str_params + pl: fit_function
      endif else begin
        id[i] = fitobj->strat2name(strategy_name=strat_arr[i]) + str_params + ' ' + pl
      endelse
    endfor  ; end of component loop over i
    obj_destroy,fitobj
  endif
endfor ; end of interval loop over int

;obj_destroy, fitobj
if total(yvals) eq 0. then return, -1 ; no data to plot
return, {yvals: yvals, id: id, strat_arr: strat_arr, ct_energy: energy}
end
  
;---------------------------------------------------------------------------

; make_func_obj method
;
; Purpose: spex_gen method to create a fit function object from either the active
; parameters/energies or the from the parameters/energies stored in the
; spex_summ structure.
;
; Keywords:
;  use_fitted - If set, use spex_summ structure values. Otherwise, use
;    current values in fit_comp... parameters.  Default=0
;  this_interval - If use_fitted is set, use params from this interval
;  use_ph_edges - If set, use photon energies instead of count energies. Default=0
;  chisq - chisq for interval if a fit was done is returned here
;  status - 0/1 means failure / success
;
; Added: March 2006
; 27-mar-2008 - don't return chisq unless using fitted results
;
;----------------------------------------------------------------------------------------

function spex_gen::make_func_obj, $
	use_fitted=use_fitted, this_interval=this_interval, $
	use_ph_edges=use_ph_edges, $
	chisq=chisq, $
	status=status

checkvar, use_fitted, 0
checkvar, this_interval, self->get(/spex_interval_index)
checkvar, use_ph_edges, 0

status = 0

if use_fitted then begin

	; if plotting fitted function, then get energies, fit params and chisq out of the
	; spex_summ structures.
	summ = self->get(/spex_summ, class='spex_fit')
	ct_energy = summ.spex_summ_energy
	fit_function = summ.spex_summ_fit_function
	params = summ.spex_summ_params
	fitobj = obj_new('fit_function', fit_function=fit_function)
	fitobj -> set, fit_comp_params=params[*,this_interval], $
		fit_comp_spectrum=summ.spex_summ_func_spectrum[*,this_interval], $
		fit_comp_model=summ.spex_summ_func_model[*,this_interval]
	chisq = summ.spex_summ_chisq[this_interval]

endif else begin

	ct_energy = self -> getaxis(/ct_energy, /edges_2)
	params = self -> get(/fit_comp_params)
	fitobj = obj_clone(self->get(/obj,class='fit_function'))
;	chisq = self ->getalg(/chisq)

endelse

fit_xvals = use_ph_edges ? self->getaxis(/ph_energy, /edges_2) : ct_energy
n_energy = n_elements(fit_xvals)/2
if n_energy le 1 then return, -1  ; either energy is a null pointer, or has bad energy values

fitobj -> set, fit_xvals=fit_xvals

status=1
return, fitobj
end

;---------------------------------------------------------------------------

; Make a plot object for displaying the fit function components.
; fit_color - array of colors indices to use for each fit component plotted
; chisq - if chisq passed in (as it is when plotting in spex_fit__xfit_comp while fitting) or
;   if chisq is returend by calc_func_components, then chisq is written in label on plot.
;   Otherwise, 'No Fit Done' is in label.
; this_interval - scalar or array of interval #'s to use (only makes senses if use_fitted is set)
; overlay - if set, we want the fitplot_obj that is returned to be used as an overlay in another xy plot obj
; no_plotman - if set, don't use plotman for plotting
; 
; Look in calc_func_components for other keywords to pass through _extra
; When called from spex_fit__xfit_comp widget to plot a single component, strategy_name is in _extra.
;

function spex_gen::make_fitplot_obj, $
	fit_color=dim1_color, $
	fit_legend_loc = fit_legend_loc, $
	chisq=chisq, $
	this_interval=this_interval, $
	overlay=overlay, $
	no_plotman=no_plotman, $
	_extra=_extra

yf = self -> calc_func_components (this_interval=this_interval, chisq=chisq, _extra=_extra)
if not is_struct(yf) then return, -1

ct_energy = yf.ct_energy
n_plot_energy = n_elements(ct_energy[0,*])
n_comp = n_elements(yf.strat_arr)	; number of components to plot
ninterval = n_elements(this_interval)

; Create fit plot object.  If just for overlay, only set dim1 info, and move label to bottom
; Otherwise, set title, etc. in fit plot object

yvals = reform(yf.yvals, n_plot_energy, n_comp*ninterval)
fitplot_obj = obj_new('xyplot', ct_energy, yvals)	;use mean of ct_energy if don't want histogram plot

if exist(chisq) then begin
	if ninterval eq 1 then begin
		label_chisq = 'Chi-square = ' + trim(chisq, '(f7.2)')
		label = 'Fit Interval ' + trim(this_interval) + '   ' + label_chisq
	endif
endif else label = 'No fit done.'

if keyword_set(overlay) then begin
	valid = 0
	if not keyword_set(no_plotman) then plotman_obj= self -> get_plotman_obj(valid=valid)
	; if using plotman get colors, otherwise will use different linestyles
	if valid then begin
		cn = plotman_obj -> get(/color_names)
		if not keyword_set(dim1_color) then begin 
		   dim1_color = [cn.red, cn.green, cn.yellow, cn.pink, cn.cyan, cn.violet, cn.olive, cn.maroon, cn.lime, cn.orange, cn.green]
		   dim1_color = [dim1_color, dim1_color]  ; make array longer than we'd ever need it, otherwise last color repeats
		   dim1_color = dim1_color[indgen(n_comp)]
		endif
	endif
	dim1_linestyles = indgen(6)+1
	fitplot_obj -> set, $
		dim1_id = yf.id, $
		dim1_color=dim1_color, $
		dim1_linestyles=dim1_linestyles, $
		label=label, $
		legend_loc=fcheck(fit_legend_loc,3)
endif else begin
	fitplot_obj -> set, $
		id='Fit Function', $
		dim1_id=yf.id, $
		dim1_enabsum=0, $
		label=label, $
		/xlog, /ylog
endelse

return, fitplot_obj
end

;---------------------------------------------------------------------------
;
pro spex_gen::units_widget, wunits, parent=parent, uvalue=uvalue, $
	units=units, getval=getval, update=update

if keyword_set(getval) then begin
	parent = widget_info(wunits, /parent)
	widget_control, parent, get_uvalue=unit_vals
	units = trim(unit_vals[widget_info(wunits, /droplist_select)])
	return
endif

allow_counts = not self->get(/spex_pseudo_livetime)
units_new = allow_counts ? ['Counts', 'Rate', 'Flux'] : ['Rate', 'Flux']
short = allow_counts ? ['co', 'ra', 'fl'] : ['ra', 'fl']

if not keyword_set(update) then begin
	checkvar, uvalue, 'none'
	; in 5.6 can't do a get_value on a droplist widget, so no way
	; to get the items in the list.  So put the droplist in an invisible
	; base, and put the list of units in the uvalue of the base.
	wbase = widget_base (parent, /row, uvalue=units_new)
	wunits = widget_droplist (wbase, title='Plot Units: ', $
		value=trim(units_new,'(a9)'), uvalue=uvalue)
	if keyword_set(units) then begin
		index = str_fuzzy (short, units)
		if index eq -1 then index=1
	endif else index = n_elements(units_new)-1
	widget_control, wunits, set_droplist_select=index
endif else begin
	parent = widget_info(wunits, /parent)
	widget_control, parent, get_uvalue=units_old
	index = widget_info(wunits, /droplist_select)
	if n_elements(units_old) lt n_elements(units_new) then index = index + 1
	if n_elements(units_old) gt n_elements(units_new) then index = (index - 1) > 0
	widget_control, wunits, set_value=trim(units_new,'(a9)'), set_droplist_select=index
	widget_control, parent, set_uvalue=units_new
endelse

end

;---------------------------------------------------------------------------

function spex_gen::locate_file, file, status=status
status = 1
   
if file[0] eq '' then begin
  status = 0
  return, ''
endif   
; make sure file exists.  If no path provided, try current directory first, then spex_data_dir

if file_test(file[0], /read) then return, file

path = file_dirname(file[0])
dir = self->get(/spex_data_dir)

if path eq '.' and dir ne '' then begin  
  file = concat_dir(dir, file)
  if file_test(file[0], /read) then return, file
  msg = 'Also checked in spex_data_dir = ' + dir
endif else msg = ''

message, 'File(s) not found or unreadable. ' + msg + ' .   NOT SETTING FILE(S).', /info
status = 0
return, ''

end

;---------------------------------------------------------------------------

pro spex_gen__define

dummy = {spex_gen, $
         INHERITS framework}
end
