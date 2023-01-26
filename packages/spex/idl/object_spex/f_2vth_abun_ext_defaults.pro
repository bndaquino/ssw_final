;+
; NAME:
;	F_2VTH_ABUN_DEFAULTS_EXT
;
; PURPOSE: Function to return default values for
;   parameters, minimum and maximum range, free parameter mask,
;   spectrum and model to use, and all spectrum and model options when
;   fitting to f_2vth_abun_ext function.
;
; CALLING SEQUENCE: defaults = f_2vth_abun_defaults_ext()
;
; INPUTS:
;	None
; OUTPUTS:
;	Structure containing default values
;
; MODIFICATION HISTORY:
; Kim Tolbert, 10-Oct-2013
; Christopher S. More, 30-Jun-2017  Added defaults for 8th parameter (Ni abundances)
;
;-
;------------------------------------------------------------------------------

FUNCTION F_2VTH_ABUN_EXT_DEFAULTS

defaults = { $
  fit_comp_params:           [1e0,   .2,    1e0,   .2,     1.,  1.,  1.,  1., 1., 1.,  1., 1.], $
  fit_comp_minima:           [1e-20, 1e-1, 1e-20, 1e-1, .01, .01, .01, .01, .01, .01, .01, .01], $
  fit_comp_maxima:           [1e20,  4.,   1e20,  4.,   10., 10., 10., 10., 10., 10., 10., 10.], $
  fit_comp_free_mask:        [1b,    1b,   1b,    1b,    1b,  0b,  0b,  0b, 0b, 0b,  0b, 0b], $

  fit_comp_spectrum:         'full', $
  fit_comp_model:            'chianti', $

  fc_spectrum_options:		['full','continuum','line'], $
  fc_model_options:			['chianti', 'mewe'] $
}

RETURN, defaults

END
