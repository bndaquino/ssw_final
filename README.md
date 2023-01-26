# Summary of changes

Major changes:
- Added spex_minxss_specfile__define (currently using the same one for MinXSS-1 and MinXSS-2)
- Added spex_daxss_specfile__define (needs testing)
- Modified spex_data__define to set strategies for MinXSS (FITS, sav) and DAXSS (just sav)
- Modified spex_datatypes to include strategies for MinXSS and DAXSS
- Modified spex_xtime_data_info.txt to include mission dates
- Kim modified spex_gen__define and spex_data_strategy__define to disable histogramming by default for MinXSS
- Kim modified spex__define to enable .sav file browsing by default in the GUI (needs to send us the updated version)
- Removed dependency on gps2utc

Todo:
- Add real DAXSS and MinXSS-2 DRMs
- For the future: Kim will tell us how to set variable DRM paths from GUI/command line
- Try implementing data download from sockets
- Consider refactoring MinXSS specfiles into separate ones for MinXSS-1 and MinXSS-2, finding a way to reuse common code
- Add test files and test data to the repository for testing on other machines
- Add compare_bins's way of computing energy bins to spex_minxss_specfile's format_to_ospex method (even though the naive way seems to be working)

Note:
- compare_bins included. Not part of SSW, but demonstrates a better way to calculate bin edges from bin centers
- data/ directory is not part of SSW, but includes some testing data to use (FITS and .sav)
