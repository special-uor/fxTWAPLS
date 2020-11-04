## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

In this version, I have:

* Removed example for `hex_logo` that was attempting "to write to the user library".
* Made utilitarian functions for internal use only, declutering the documentation.
* Remove the external datasets, as requested for our data provider.
* Due to the previous point, examples can't be run, thus I have added the `\dontrun` flag.
