## Test environments
* Ubuntu 20.04 LTS with gcc 10.1.0
  R version 4.2.1
* Ubuntu 20.04 LTS with gcc 10.1.0
  R version 4.2.1 using Valgrind
* Ubuntu 20.04 LTS with gcc 10.1.0
  R devel 2022-06-09 r82474 with ASAN and UBSAN
* GitHub actions on windows-latest (release), macOS-latest (release), 
  ubuntu-20.04 (release), ubuntu-20.04 (old-release), and ubuntu-20.04 (devel)
* win-builder (devel, oldrelease, and release)
  
## R CMD check results
There were no WARNINGs or ERRORs.

There is a NOTE about the package size in some cases.

This is a resubmission. The following has been changed:

 - The Description field has been updated. The redundant "Provides functions 
   to" has been removed and the adverb "possibly" has been removed.
 - The LICENSE.md file has been removed.
 - There is currently no reference describing the method but a paper that is
   not on arXiv is submitted. A reference to this paper will be added later.
 - \value sections have been added to plot_marker() and plot_surv().
