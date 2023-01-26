;+
; Project     : HESSI
;
; Name        : SPEX_DATATYPES
;
; Purpose     : Function to return the data types available in OSPEX, as well as some information about them.
;
; Category    : spex
;
; Syntax      : datatypes = spex_datatypes()
;
; Input:
;   input - optional. class name or object of one of these classes to return information for
;
; Keywords    :
;   class - if set, just return array of classes
;
; Output      : Array of structures containing:
;    class - object class name
;    name  - data type nume
;    web   - byte, 0 or 1 means not/ is available via http over web
;    drmfile - byte, 0 or 1 means doesn't / does require drm file
;    dets  - string array of detectors available for data type
;    life  - string containing lifetime duration of data type
;
; Written: Kim Tolbert February 2010
; Modifications:
; 07-Jul-2010, Kim. Added iput arg, and comment tag for returned structure
; 25-Jan-2011, Kim. Added Fermi LAT
; 19-Feb-2013, Kim. Added SMM HXRBS
; 18-Dec-2013, Kim. Set LAT web byte on to enable getting over web
; 10-Oct-2014, Kim. Added SPEX_ANY_SPECFILE
; 14-Jul-2015, Kim. Set SMM HXRBS web byte on to enable getting over web, and updated soxs and messenger lifetimes
; 28-Jul-2016, Kim. Added konus, and sort structure by name field.
; 04-Aug-2016, Kim. Corrected Konus lifetime
; 08-Sep-2016, Kim. Changed konus currently available data to 2012-2016.
; 27-Sep-2016, Richard.Schwartz@nasa.gov, added SPEX_STIX_IMAGE for STIX image cube
; 24-Aug-2017, Kim. Changed konus availability to 2012 - present.
; 21-Mar-2019, Kim. Added SMM GRS
; 31-Oct-2022, Brendan added MinXSS (s now has 15 elements)
;-
function spex_datatypes, input, class=class

  checkvar, class, 0

  ; class name, short name, files on web, use drm file, detectors, lifetime, comment

  struct = {class: '', name: '', web: 0b, drmfile: 0b, dets: strarr(30), life: '', comment: ''}
  tags = tag_names(struct)

  s = replicate(struct,16)

  dets = strarr(30)
  ;soxs_dets = ['Si', 'CZT', strarr(28)]
  gbm_dets = ['b0','b1','n'+trim(indgen(10)), 'na', 'nb', strarr(16)]

  drm_file = 1

  life = ''
  life_soxs = '1-Jul-2003 - 2-May-2011'
  life_mess = '~2004 - 17-Sep-2014, spotty coverage'
  life_gbm = '12-Aug-2008 - present'
  life_lat = '12-Aug-2008 - present'
  life_hxrbs = '19-Feb-1980 - 21-Nov-1989. Flares.'
  life_konus = '12-Nov-1994 - present. Flares.'
  life_grs = '19-Feb-1980 - 21-Nov-1989. Flares.'

  comment = ''
  comment_soxs = 'CZT and SI data are in same file.'
  comment_wbs = 'GRS1,GRS2,HXS are in same file.'
  comment_konus = '*_1.pha and *_1.rmf are for 20-1250 keV, *_2.pha and *_2.rmf are for 0.25-15.8 MeV'


  s[0] = create_struct(tags, 'SPEX_HESSI_SPECFILE', 'HESSI', 0b, 1b, dets, life, comment)
  s[1] = create_struct(tags, 'SPEX_HESSI_IMAGE', 'HESSI IMAGE', 0b, 0b, dets, life, comment)
  s[2] = create_struct(tags, 'SPEX_USER_DATA', 'USER',  0b, 2b, dets, life, comment)
  s[3] = create_struct(tags, 'SPEX_XSM_SPECFILE', 'XSM', 0b, 1b, dets, life, comment)
  s[4] = create_struct(tags, 'SPEX_SOXS_SPECFILE', 'SOXS', 1b, 0b, dets, life_soxs, comment_soxs)
  s[5] = create_struct(tags, 'SPEX_MESSENGER_SPECFILE', 'MESSENGER', 1b, 0b, dets, life_mess, comment)
  s[6] = create_struct(tags, 'SPEX_FERMI_GBM_SPECFILE', 'FERMI_GBM', 1b, 1b, gbm_dets, life_gbm, comment)
  s[7] = create_struct(tags, 'SPEX_YOHKOH_WBS_SPECFILE', 'YOHKOH_WBS', 0b, 0b, dets, life, comment_wbs)
  s[8] = create_struct(tags, 'SPEX_FERMI_LAT_SPECFILE', 'FERMI_LAT', 1b, 1b, dets, life_lat, comment)
  s[9] = create_struct(tags, 'SPEX_SMM_HXRBS_SPECFILE', 'SMM_HXRBS', 1b, 0b, dets, life_hxrbs, comment)
  s[10]= create_struct(tags, 'SPEX_KONUS_SPECFILE', 'KONUS', 1b, 1b, dets, life_konus, comment_konus)
  s[11]= create_struct(tags, 'SPEX_ANY_SPECFILE', 'ANY', 0b, 2b, dets, life, comment)
  s[12]= create_struct(tags, 'SPEX_STIX_IMAGE', 'STIX IMAGE', 0b, 0b, dets, life, comment)
  s[13] = create_struct(tags, 'SPEX_SMM_GRS_SPECFILE', 'SMM_GRS', 1b, 0b, dets, life_grs, comment)
  s[14] = create_struct(tags, 'SPEX_MINXSS_SPECFILE', 'MINXSS', 0b, 0b, dets, life, comment)
  s[15] = create_struct(tags, 'SPEX_DAXSS_SPECFILE', 'DAXSS', 0b, 0b, dets, life, comment)

  ; Sort alphabetically by name field
  sorted = sort(s.name)
  s = s[sorted]

  if exist(input) then begin
    findclass = obj_valid(input) ? obj_class(input) : input
    q = where(stregex(s.class, strupcase(findclass), /bool) eq 1, count)
    if count gt 0 then return, class ? s[q[0]].class : s[q[0]] else return, -1
  endif

  return, class ? s.class : s
end