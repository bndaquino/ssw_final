# Summary of changes

Major changes:
- Added spex_minxss_specfile__define (currently using the same one for MinXSS-1 and MinXSS-2)
- Added spex_daxss_specfile__define
- Modified spex_data__define to set strategies for MinXSS (FITS, sav) and DAXSS (just sav)
- Modified spex_datatypes to include strategies for MinXSS and DAXSS
- Modified spex_xtime_data_info.txt to include mission dates

Other changes (maybe to be reverted):
- Added if fc.fit_comp_function eq 'vth_abun' then stop to fit_comp__define? See line 99.
- Modified spex_fit__xfit_comp for loop on line 242 for debugging -- should change this back

Todo:
- Add Chris's comments to spex_xtime_data_info.txt
- Try implementing data download from sockets
	- send Crisel Messenger example
- Consider refactoring MinXSS specfiles into separate ones for MinXSS-1 and MinXSS-2, finding a way to reuse common code
- Add test files and test data to the repository for testing on other machines
- Add real DAXSS and MinXSS-2 DRMs
- Add new Chianti files (currently only tracking /packages/spex/idl -- I think these are in a different directory)
- Find a good place to add gps2utc and its dependencies (yd2jd and tai_utc)
	- instead, removing dependency on gps2utc by using anytim (built-in to SSW)
- Add compare_bins's way of computing bins (even though the naive way seems to be working)
- For the future: Kim will tell us how to set variable DRM paths from GUI/command line