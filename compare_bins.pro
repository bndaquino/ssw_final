; This file is not part of SSW. It just demonstrates how to calculate bin edges
; given an array of bin centers if the bin centers are not evenly spaced.
; The resulting edges will line up, leaving no overlapping bins.

pro compare_bins

  restore, '~/data/minxss1/my_m5.sav'
  data_centers = minxsslevel1.x123[0].energy
  bin_widths = data_centers[1:-1] - data_centers
  
  data_ebins = dblarr(2, n_elements(data_centers))
  last_center_idx = n_elements(data_centers) - 1
  data_ebins[0,0] = data_centers[0] - bin_widths[0]/2
  data_ebins[1,0] = data_centers[0] + bin_widths[0]/2
  data_ebins[0,last_center_idx] = data_centers[last_center_idx] - bin_widths[-1]/2
  data_ebins[1,last_center_idx] = data_centers[last_center_idx] + bin_widths[-1]/2
  for i = 1, n_elements(bin_widths) - 1 do begin
    ; bins_widths[i] is the bin width between data_centers[i] and data_centers[i+1]
    rw = bin_widths[i]
    lw = bin_widths[i-1]
    ; data_ebins[0,i] is the left bin of the bucket that is centered on data_centers[i] 
    data_ebins[0,i] = data_centers[i] - lw/2
    data_ebins[1,i] = data_centers[i] + rw/2
  endfor
  
  drm = mrdfits('~/MinXSS_OSPEX/Code/cmoore/drm_complete_minxss_x123_fm1_all_ospex_n_26_no_electrons.fits', 1)
  drm_in = drm.edges_in
  drm_out = drm.edges_out
  
  
  ; Trimming
  
  min_kev = 0.3
  max_kev = 25.

  data_ebins_trimmed = data_ebins[*,where((data_ebins[0,*] ge min_kev) and (data_ebins[1,*] le max_kev))]
  
  drm_in_trimmed = drm_in[*,where((drm_in[0,*] ge min_kev) and (drm_in[1,*] le max_kev))]
  drm_out_trimmed = drm_out[*,where((drm_out[0,*] ge min_kev) and (drm_out[1,*] le max_kev))]
  
  stop

end