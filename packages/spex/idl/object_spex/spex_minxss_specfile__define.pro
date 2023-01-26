;+
;
; NAME:
;   spex_minxss_specfile__define
;
; PURPOSE:
;   Provides read_data method for both MinXSS-1 and MinXSS-2.
;   NOTE: This procedure is based on spex_hessi_specfile__define.
;
;
; CATEGORY:
;       SPECTRAL FITTING, SPEX - minxss
;
; CALLING SEQUENCE:
;
; CALLS:
;   read_messenger_4_ospex
;
; INPUTS:
;
; OUTPUTS:
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; PROCEDURE:
;
; Written 19-Jan-2022, Brendan D'Aquino
; 
; Modification History
; 19-Jan-2023, Kim. Added get_time_plot_args method to allow setting default for histogram option to 0
;------------------------------------------------------------------------------

; Reformats mission data for use in OSPEX.
; Input: a data structure, expected to be formatted as the level 1 data on the MinXSS data products site,
;   but with inner the .time structure unbundled (see self->unbundle_time).
; Output: reformatted data to be used in OSPEX. Performs basic tasks like:
;   calculating energy bin edges from the data's energy bin centers,
;   calculating UT time edges from data's TAI-formatted times,
;   and packing other relevant information into an array of structures, with one structure per integration period.
function spex_minxss_specfile::format_to_ospex, minxss_data_structure

  ; minxss_data_structure is the x123 field of the .sav data product
  ; So it's an array of structures, where each structure includes a spectrum and metadata
  ;   for a certain integration period

  ; UTC and GPS count seconds from Jan 6 1980 (with GPS ignoring leap seconds)
  ; Seems like UT counts from Jan 1 1979 instead
  ; So we construct an array of integration time edges in GPS, then convert to UTC, and then
  ;   add to that the number of seconds between Jan 1 1979 and Jan 6 1980 to finally get UT time edges (which OSPEX needs)
  ; ut_edges is a [2,n_times] array, where [0,m] is the start of the mth time and [1,m] the end
  seconds_between_19790101_19800106 = 60.*60.*24.*(365+5)
  tai_edges = transpose([[minxss_data_structure.tai], $
    [minxss_data_structure.tai + minxss_data_structure.integration_time]])
  ut_edges = anytim(tai_edges, /sec, fiducial='tai')

  ; I don't know how the previous code handled energy edges. It seems like energy_bin_center_array was passed
  ; in as an array of edges (not of centers) that could just be returned as is. Instead, I'm going to calculate
  ; the edges based on the energy centers I find in the data product. I'm assuming that there are no gaps
  ; between energy bins, so that the width of each energy bin is equal to the distance between adjacent centers,
  ; and also that the bin width is the same across all bins.
  ; energy_edges is a [2, n_energy] array where [0,n] is the start of the nth bin and [1,n] the end
  energy_bin_width = minxss_data_structure[0].energy[1] - minxss_data_structure[0].energy[0]
  energy_edges = transpose([[minxss_data_structure[0].energy - energy_bin_width/2], $
    [minxss_data_structure[0].energy + energy_bin_width/2]])

  n_times = n_elements(ut_edges[0,*])  ; number of integration periods

  minxss_x123_ospex_structure_0 = {total_counts: minxss_data_structure[0].spectrum_total_counts, $
    uncertainty_total_counts: minxss_data_structure[0].spectrum_total_counts_precision, $
    slow_counts: minxss_data_structure[0].x123_slow_count, $
    integration_time: minxss_data_structure[0].integration_time, $
    count_rate: minxss_data_structure[0].spectrum_cps, $
    uncertainty_count_rate: minxss_data_structure[0].spectrum_cps_precision, $
    energy_bins: energy_edges, $
    time_HUMAN: minxss_data_structure[0].HUMAN, $
    start_time: ut_edges[0,0], $
    end_time: ut_edges[1,0], $
    minxss_version_flag: minxss_data_structure[0].flight_model}

  minxss_x123_ospex_structure_temp = replicate(minxss_x123_ospex_structure_0, n_times)

  for m = 0, n_times - 1 do begin
    minxss_x123_ospex_structure_temp[m].total_counts = minxss_data_structure[m].spectrum_total_counts
    minxss_x123_ospex_structure_temp[m].uncertainty_total_counts = minxss_data_structure[m].spectrum_total_counts_precision
    minxss_x123_ospex_structure_temp[m].slow_counts = minxss_data_structure[m].x123_slow_count
    minxss_x123_ospex_structure_temp[m].integration_time = minxss_data_structure[m].integration_time
    minxss_x123_ospex_structure_temp[m].count_rate = minxss_data_structure[m].spectrum_cps
    minxss_x123_ospex_structure_temp[m].uncertainty_count_rate = minxss_data_structure[m].spectrum_cps_precision
    minxss_x123_ospex_structure_temp[m].energy_bins = energy_edges
    minxss_x123_ospex_structure_temp[m].time_HUMAN = minxss_data_structure[m].HUMAN
    minxss_x123_ospex_structure_temp[m].start_time = ut_edges[0,m]
    minxss_x123_ospex_structure_temp[m].end_time = ut_edges[1,m]
    minxss_x123_ospex_structure_temp[m].minxss_version_flag = minxss_data_structure[m].flight_model
  endfor

  output_minxss_x123_ospex_structure = minxss_x123_ospex_structure_temp
  
  return, output_minxss_x123_ospex_structure

end

; The .sav mission data has time bundled into its own .time structure, but the FITS version instead has the
; time formats at the top level. This method copies the fields of the .time structure into the top level of
; each structure in the .sav data, creating a uniform formatting between .sav and FITS files that allows
; the code in self->format_to_ospex to be used on either file type.  
function spex_minxss_specfile::unbundle_time, sav_data

  time_formats = {iso: sav_data[0].time.iso, $
    human: sav_data[0].time.human, $
    yyyymmdd: sav_data[0].time.yyyymmdd, $
    yyyydoy: sav_data[0].time.yyyydoy, $
    hhmmss: sav_data[0].time.hhmmss, $
    sod: sav_data[0].time.sod, $
    fod: sav_data[0].time.fod, $
    jd: sav_data[0].time.jd, $
    tai: sav_data[0].time.tai, $
    spacecraftgpsformat: sav_data[0].time.spacecraftgpsformat}
    
  template = create_struct(sav_data[0], time_formats)
  unbundled = replicate(template, n_elements(sav_data))
  struct_assign, sav_data, unbundled
  
  for i = 0, n_elements(sav_data) - 1 do begin
    unbundled[i].iso = sav_data[i].time.iso
    unbundled[i].human = sav_data[i].time.human
    unbundled[i].yyyymmdd = sav_data[i].time.yyyymmdd
    unbundled[i].yyyydoy = sav_data[i].time.yyyydoy
    unbundled[i].hhmmss = sav_data[i].time.hhmmss
    unbundled[i].sod = sav_data[i].time.sod
    unbundled[i].fod = sav_data[i].time.fod
    unbundled[i].jd = sav_data[i].time.jd
    unbundled[i].tai = sav_data[i].time.tai
    unbundled[i].spacecraftgpsformat = sav_data[i].time.spacecraftgpsformat
  endfor
  
  return, unbundled

end

; Gets mission data from the user-set specfile (.sav or FITS) in an OSPEX-friendly format.
function spex_minxss_specfile::get_minxss_data

  path = self->get(/spex_specfile)
  
  if strpos(path, '.sav') ne -1 then begin
    restore, path
    unbundled = self->unbundle_time(minxsslevel1.x123)
    data = self->format_to_ospex(unbundled)
    return, data
  endif else if strpos(path, '.fits') ne -1 then begin
    minxsslevel1 = mrdfits(path, 1)
    data = self->format_to_ospex(minxsslevel1)
    return, data
  endif else begin
    message, 'MinXSS specfile must be either .fits or .sav, but given path was ' + path
  endelse

end

; Gets the MinXSS-1 or MinXSS-2 DRM from SSW.
function spex_minxss_specfile::get_minxss_drm
  ; How will OSPEX get this path? Where does the path go when I call o->set, spex_drmfile?
  path = '/Users/bdaquino/MinXSS_OSPEX/Code/cmoore/drm_complete_minxss_x123_fm1_all_ospex_n_26_no_electrons.fits'
  drm = mrdfits(path, 1)
  
  min_energy_kev = 0.3
  max_energy_kev = 25.
  index_in_range = WHERE((drm.edges_out[0,*] gt min_energy_kev) $ 
    and (drm.edges_out[0,*] lt max_energy_kev) $
    and (drm.edges_in[0,*] gt min_energy_kev) $
    and (drm.edges_in[0,*] lt max_energy_kev), n_index_in_range)
  drm_index_min = min(index_in_range)
  drm_index_max = max(index_in_range)
  
  respinfo = {drm: drm.repsonse_matrix[drm_index_min:drm_index_max, drm_index_min:drm_index_max], $
    edges_in: drm.edges_in[*,index_in_range], $
    edges_out: drm.edges_out[*,index_in_range]}
  
  to_return = {area: drm.aperture_area_cm, $
    respinfo: respinfo}
    
  return, to_return

end

; Selects from a given array of bins ([2,n]) the ones between a min and max value
; and returns the indexes of the bins kept.
function spex_minxss_specfile::trim_bins, bins, min_val, max_val, idx=idx

  idx_in_range = where((bins[0,*] gt min_val) and (bins[0,*] lt max_val))
  idx = idx_in_range
  return, bins[*,idx_in_range]

end

; Determines whether the specfile has data for MinXSS-1 or MinXSS-2.
function spex_minxss_specfile::get_instrument

  path = self->get(/spex_specfile)
  is_fits = strpos(path, '.fits') ne -1
  
  if is_fits then begin
    fits_read, '~/data/minxss1/my_m5.fits', dat, hdr, /header_only, /pdu
    instrument_idx = where(strpos(hdr, 'INSTRUME') ne -1)
    instrument_card = hdr[instrument_idx]
    name = strmid(instrument_card, strpos(instrument_card, 'MinXSS'), strlen('MinXSS-X'))
    return, name
  endif
  
  if (strpos(strlowcase(path), 'minxss1') ne -1) or (strpos(strlowcase(path), 'minxss-1') ne -1) then begin
    return, 'MinXSS-1'
  endif

  if (strpos(strlowcase(path), 'minxss2') ne -1) or (strpos(strlowcase(path), 'minxss-2') ne -1) then begin
    return, 'MinXSS-2'
  endif

end

;------------------------------------------------------------------------------

pro spex_minxss_specfile::read_data, file_index=file_index, $
  spectrum,  errors,  livetime,  $
  spex_respinfo, spex_file_time,  spex_ut_edges,  spex_ct_edges,  $
  spex_area, spex_title, spex_detectors,  $
  spex_interval_filter, spex_units, spex_data_name, $
  spex_deconvolved, spex_pseudo_livetime, spex_data_pos, spex_def_eband=spex_def_eband, $
  err_code=err_code, _extra=_extra
  
  data = self->get_minxss_data()
  drm  = self->get_minxss_drm()
  instrument = self->get_instrument()
  
  n_times = n_elements(data)
  
  min_energy_kev = 0.3
  max_energy_kev = 25.02  ; for MinXSS-1: arbitrary +0.02 to match dimensions with DRM (831)
  trimmed_ebins = self->trim_bins(data[0].energy_bins, min_energy_kev, max_energy_kev, idx=ebin_idx)
  n_energy = n_elements(trimmed_ebins[0,*])
  
  ltimes = dblarr(n_energy, n_times)
  for k = 0, n_times - 1 do begin
    ltimes[*,k] = data[k].integration_time
  endfor
  
  spectrum = data.total_counts[ebin_idx, *]                         ; dimensions: [n_energy, n_times]
  errors = data.uncertainty_total_counts[ebin_idx, *]               ; dimensions: [n_energy, n_times]
  livetime = ltimes                                                 ; dimensions: [n_energy, n_times]
  spex_file_time = [data[0].start_time, data[-1].end_time]          ; Start and end times of the entire file
  spex_ut_edges = transpose([[data.start_time], [data.end_time]])   ; dimensions: [2, n_times], format: UT

  spex_respinfo = drm.respinfo                                      ; structure with {matrix: [n_out, n_in], edges_in: [2, n_in], edges_out: [2,n_out]}
  spex_ct_edges = trimmed_ebins                                     ; dimensions: [2, n_energy]
  spex_area = drm.area                                              ; units: cm^2
  spex_title = instrument + ' Spectrum'
  spex_detectors = 'X123'
  spex_units = 'counts'
  spex_interval_filter = -1
  spex_data_name = instrument
  spex_deconvolved = 0
  spex_pseudo_livetime = 1
  spex_data_pos = [0.,0.]

end

function spex_minxss_specfile::get_time_plot_args ; kim added this function 19-jan-2023
  return, {histogram: 0}
end

;------------------------------------------------------------------------------
pro spex_minxss_specfile__define

  self = {spex_minxss_specfile, $
    INHERITS spex_data_strategy }

END
