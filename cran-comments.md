## New version submission

CRAN archived VAJointSurv on 2026-02-15 because SimSurvNMarker was unavailable.
SimSurvNMarker 0.1.4 is now on CRAN.

Version 0.1.2 replaces the native pipe in the plot_surv() example.
This resolves the R 4.1 syntax NOTE while keeping the declared R minimum.

## Local checks

On Debian 13 (x86_64), R 4.5.0 with SimSurvNMarker 0.1.4:

    R CMD build VAJointSurv
    R CMD check --as-cran VAJointSurv_0.1.2.tar.gz

The check completed with zero errors, zero warnings, and one NOTE.
The incoming feasibility check identified the package as a new submission
because CRAN had archived the previous release.
