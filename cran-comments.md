## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Resubmission
This is a resubmission. In this version I have:

* Add local binding for local variables

* Excluded the README.Rmd from build

* Removed package name from the DESCRIPTION file

## Resubmission (2020-10-06)
* Removing `\dontrun` blocks for examples and replacing for `\donttest` (for long 
examples)

* Limitting the number of cores to 2 for examples

* Describing acronyms in DESCRIPTION file

* Removing `<<-` operator to `<-`

* Updating to version `0.0.2`
