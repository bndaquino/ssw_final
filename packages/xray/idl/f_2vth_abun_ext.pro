;+
;Name:
;   F_2VTH_ABUN_EXT
;PURPOSE:
;   This function returns the sum of two f_vth functions.  f_vth returns the 
;   optically thin thermal bremsstrahlung radiation function
;   as differential spectrum seen at Earth in units of photon/(cm2 s keV).
;   Input to f_2vth_abun includes two pairs of emission measure and temperature, but only one
;   set of abundance values.
;   Same as f_2vth except relative abundances for all low-fip elements are controlled by
;   separate parameters.
;
;CATEGORY:
;   SPECTRA, XRAYS (>1 keV)
;INPUTS:
;   E  energy vector in keV
;   apar[0]  em_49, emission measure units of 10^49
;   apar[1]  KT, plasma temperature in keV
;   apar[2]  em_49, emission measure units of 10^49
;   apar[3]  KT, plasma temperature in keV
;   apar[4]  Relative abundance for Fe
;   apar[5]  Relative abundance for Ca
;   apar[6]  Relative abundance for S
;   apar[7]  Relative abundance for Mg
;   apar[8]  Relative abundance for Si
;   apar[9]  Relative abundance for Ar
;   apar[10]  Relative abundance for He, C, N, O, F, Ne, Na, Al, K
;   apar[11]  Relative abundance for Ni
;            Abundances are relative to coronal abundance for Chianti
;            Relative to solar abundance for Mewe
;           (unless user selects a different abundance table manually)
;           
;KEYWORD INPUTS:
;   LINES - Only return lines
;   CONTINUUM - Only return continuum (same as setting NOLINE keyword)
;   NOLINE - Only return continuum
;   CHIANTI - If set, use chianti_kev
;   MEWE - If set, use mewe_kev
;   REL_ABUN - 2x2 array giving Fe, Ni abundance [ 26,x],[28,x] ],  If rel_abun keyword not used,
;     the value of x is taken from apar(2) to make 2x5 array for Fe, Ni, Si, and Ca.
;	  S is also included but at half the deviation from nominal as the others.  If that's not there either, x is 1.
;
;   Defaults are to return full spectrum (lines+continuum), chianti
;
;CALLS:
;   mk_contiguous, BREM_49, MEWE_KEV, CHIANTI_KEV
;
;Common Blocks:
;   None
;
;Method:
;   If energy array starts at gt 8 keV, then noline is set to calc pure
;     free-free continuum from either chianti_kev or mewe_kev.
;   If edges weren't passed as 2xn array, use Brem_49 function.
;
; History:
;   Kim, 10-Oct-2013
;  22-Nov-2022, Christopher S. Moore. Added changes to accommodate MinXSS-1, MinXSS-2, DAXSS, and other soft X-ray spectrometers;                                                                                   ;-

function f_2vth_abun_ext, e, apar,	_extra=_extra

em1 = apar[0]
t1 = apar[1]
em2 = apar[2]
t2 = apar[3]
abun = apar[4:*]

return, f_vth_abun(e, [em1,t1,abun], _extra=_extra)  +  f_vth_abun(e, [em2,t2,abun], _extra=_extra)
end
