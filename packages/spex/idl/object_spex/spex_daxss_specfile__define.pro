; This file contains spex_daxss_specfile's read_data method and a few helpers

; Based on minxss_x123_level1_make_ospex_structure
; Just reformats the .sav data to be easier to use for OSPEX
function spex_daxss_specfile::sav_to_ospex, daxss_data_structure

  ; daxss_data_structure is the x123 field of the .sav data product
  ; So it's an array of structures, where each structure includes a spectrum and metadata
  ;   for a certain integration period (maybe ~30 seconds each)

  ; UTC and GPS count seconds from Jan 6 1980 (with GPS ignoring leap seconds)
  ; Seems like UT counts from Jan 1 1979 instead
  ; So we construct an array of integration time edges in GPS, then convert to UTC, and then
  ;   add to that the number of seconds between Jan 1 1979 and Jan 6 1980 to finally get UT time edges (which OSPEX needs)
  ; ut_edges is a [2,n_times] array, where [0,m] is the start of the mth time and [1,m] the end
  seconds_between_19790101_19800106 = 60.*60.*24.*(365+5)
  utc_edges = gps2utc(transpose([[daxss_data_structure.time_gps], $
    [daxss_data_structure.time_gps + daxss_data_structure.integration_time]]))
  ut_edges = utc_edges + seconds_between_19790101_19800106

  ; I don't know how the previous code handled energy edges. It seems like energy_bin_center_array was passed
  ; in as an array of edges (not of centers) that could just be returned as is. Instead, I'm going to calculate
  ; the edges based on the energy centers I find in the data product. I'm assuming that there are no gaps
  ; between energy bins, so that the width of each energy bin is equal to the distance between adjacent centers,
  ; and also that the bin width is the same across all bins.
  ; energy_edges is a [2, n_energy] array where [0,n] is the start of the nth bin and [1,n] the end
  energy_bin_width = daxss_data_structure[0].energy[1] - daxss_data_structure[0].energy[0]
  energy_edges = transpose([[daxss_data_structure[0].energy - energy_bin_width/2], $
    [daxss_data_structure[0].energy + energy_bin_width/2]])

  ; In the below code, I deleted the index_valid_drm parts. I'll just return all the data that the .sav file has.
  ; Deleted the section for average_data
  ; Deleted the field called ut_edges. I don't think I need to put the time edges in each and every structure
  ;   since each structure will tell its own start and end times anyway.
  ; I don't have access to the detector area here. Maybe I can get it from the DRM?
  ; Removed the n_times gt 1 guard on the for loop. I think it's okay to run the for loop even just once?
  ;   I don't think there's any issue with using replicate(x,1)
  ; Removed all the other time fields (except human, just for convenience)

  n_times = n_elements(ut_edges[0,*])  ; number of integration periods

  daxss_x123_ospex_structure_0 = {total_counts: daxss_data_structure[0].spectrum_cps * daxss_data_structure[0].integration_time, $
    uncertainty_total_counts: daxss_data_structure[0].spectrum_cps_accuracy * daxss_data_structure[0].integration_time, $
    integration_time: daxss_data_structure[0].integration_time, $
    energy_bins: energy_edges, $
    time_ISO: daxss_data_structure[0].time_ISO, $
    start_time: ut_edges[0,0], $
    end_time: ut_edges[1,0]}

  daxss_x123_ospex_structure_temp = replicate(daxss_x123_ospex_structure_0, n_times)

  for m = 0, n_times - 1 do begin
    daxss_x123_ospex_structure_temp[m].total_counts = daxss_data_structure[m].spectrum_cps * daxss_data_structure[m].integration_time
    daxss_x123_ospex_structure_temp[m].uncertainty_total_counts = daxss_data_structure[m].spectrum_cps_precision * daxss_data_structure[m].integration_time
    daxss_x123_ospex_structure_temp[m].integration_time = daxss_data_structure[m].integration_time
    daxss_x123_ospex_structure_temp[m].energy_bins = energy_edges
    daxss_x123_ospex_structure_temp[m].time_ISO = daxss_data_structure[m].time_ISO
    daxss_x123_ospex_structure_temp[m].start_time = ut_edges[0,m]
    daxss_x123_ospex_structure_temp[m].end_time = ut_edges[1,m]
  endfor
  
  return, daxss_x123_ospex_structure_temp 

end


function spex_daxss_specfile::get_daxss_data

  path = self->get(/spex_specfile)
  
  restore, path
  data = self->sav_to_ospex(daxss_level1_data)
  return, data

end


function spex_daxss_specfile::get_daxss_drm
  ; How will OSPEX get this path? Where does the path go when I call o->set, spex_drmfile?
;  path = '~/data/daxss/minxss_fm3_ARF.fits'
;  drm = mrdfits(path, 1)
;  
;  min_energy_kev = 0.3
;  max_energy_kev = 25.
;  index_drm_range = WHERE((drm.edges_out[0,*] gt min_energy_kev) and (drm.edges_out[0,*] lt max_energy_kev) and (drm.edges_in[0,*] gt min_energy_kev) and (drm.edges_in[0,*] lt max_energy_kev), n_index_drm_range)
;  drm_index_min = min(index_drm_range)
;  drm_index_max = max(index_drm_range)
;  
;  
;  
;  trimmed_in = self->trim_bins(drm.edges_in, idx=idx_in)
;  trimmed_out = self->trim_bins(drm.edges_out, idx=idx_out)
;  trimmed_matrix = drm.repsonse_matrix[idx_out,*]
;  trimmed_matrix = trimmed_matrix[*,idx_in]
;  
;  respinfo = {drm: drm.repsonse_matrix[drm_index_min:drm_index_max, drm_index_min:drm_index_max], $
;    edges_in: drm.edges_in[*,index_drm_range], $
;    edges_out: drm.edges_out[*,index_drm_range]}
  
  dummy_respinfo = {drm: identity(1001), $
    edges_in: findgen(1001)/1001*(25-3) + 3, $
    edges_out: findgen(1001)/1001*(25-3) + 3}
    
  dummy_area = 1.
  
  to_return = {area: dummy_area, $
    respinfo: dummy_respinfo}

  return, to_return
  
end

; Selects from a given array of bins ([2,n]) the ones between a min and max value
; and returns the indexes of the bins kept
function spex_daxss_specfile::trim_bins, bins, idx=idx

  min_energy_kev = 0.30
  max_energy_kev = 25.00
;  idx_in_range = where((bins[0,*] gt min_energy_kev) and (bins[1,*] lt max_energy_kev))
  idx_in_range = where((bins[0,*] gt min_energy_kev) and (bins[0,*] lt max_energy_kev+0.02))  ; arbitrary to match dimensions with drm (831)
  idx = idx_in_range
  return, bins[*,idx_in_range]

end

;------------------------------------------------------------------------------

pro spex_daxss_specfile::read_data, file_index=file_index, $
  spectrum,  errors,  livetime,  $
  spex_respinfo, spex_file_time,  spex_ut_edges,  spex_ct_edges,  $
  spex_area, spex_title, spex_detectors,  $
  spex_interval_filter, spex_units, spex_data_name, $
  spex_deconvolved, spex_pseudo_livetime, spex_data_pos, spex_def_eband=spex_def_eband, $
  err_code=err_code, _extra=_extra
  
  data = self->get_daxss_data()
  drm  = self->get_daxss_drm()
  
  n_times = n_elements(data)
  
  trimmed_ebins = self->trim_bins(data[0].energy_bins, idx=ebin_idx)
  n_energy = n_elements(trimmed_ebins[0,*])
  
  ltimes = dblarr(n_energy, n_times)
  for k = 0, n_times - 1 do begin
    ; Same livetime for each energy bin
    ltimes[*,k] = data[k].integration_time
  endfor
  stop
  spectrum = data.total_counts[ebin_idx, *]                         ; dimensions: [n_energy, n_times]
  errors = data.uncertainty_total_counts[ebin_idx, *]               ; dimensions: [n_energy, n_times]
  livetime = ltimes                                                 ; dimensions: [n_energy, n_times]
  spex_file_time = [data[0].start_time, data[-1].end_time]          ; Are these supposed to be the start and end times of the entire file?
                                                                    ;   HESSI's don't look like a number of seconds since Jan 1 1979 -- what are the units?
                                                                    ;   Also this is just computed later using minmax(spex_ut_edges) (line 230 spex_data_strategy)
  spex_ut_edges = transpose([[data.start_time], [data.end_time]])   ; dimensions: [2, n_times], units: seconds since Jan 1 1979

  spex_respinfo = drm.respinfo                                      ; structure with {matrix: [n_out, n_in], edges_in: [2, n_in], edges_out: [2,n_out]}
  spex_ct_edges = trimmed_ebins                                     ; dimensions: [2, n_energy]
  spex_area = drm.area                                              ; units: cm^2
  spex_title = 'DAXSS Spectrum'
  spex_detectors = 'X123'
  spex_units = 'counts'
;  spex_interval_filter = ''
  spex_interval_filter = -1
  spex_data_name = 'DAXSS'
  spex_deconvolved = 0
  spex_pseudo_livetime = 1
  spex_data_pos = [0.,0.]

end

;------------------------------------------------------------------------------
pro spex_daxss_specfile__define

  self = {spex_daxss_specfile, $
    INHERITS spex_data_strategy }

END
