## New version submission

CRAN archived VAJointSurv on 2026-02-15 because SimSurvNMarker was unavailable.
SimSurvNMarker 0.1.4 is now on CRAN, and VAJointSurv 0.1.2 removes the R 4.1
pipe from the plot_surv() example.

## Check environments and results

- Debian 13 (x86_64), R 4.5.0 with SimSurvNMarker 0.1.4, `--as-cran`.
  Zero errors and zero warnings. The incoming feasibility NOTE identifies
  the archived package as a new submission.
- macOS Sequoia 1.52 with a M1 chip with R 4.6.2.
- devtools::check_win_devel(), devtools::check_win_release(), and
  devtools::check_win_oldrelease(). There is only a similar note.
- With `--no-manual --as-cran`, checks passed on Ubuntu 26.04.1 LTS
  (x86_64, R 4.5.3, R 4.6.2, and R-devel r90971), macOS 26.7.2
  (arm64, R 4.6.2), and Windows Server 2022 (x86_64, R 4.6.2).
- Ubuntu 26.04.1 LTS (x86_64), R-devel r90971 and r90985, Clang 23.1.6,
  `--no-manual --no-stop-on-test-error`. Separate AddressSanitizer and
  UndefinedBehaviorSanitizer checks passed, including examples and tests.
- Fedora Linux 44 (x86_64), R-devel r90971, GCC 16.2.1 with address and
  undefined behavior sanitizers, `--no-manual --no-build-vignettes --as-cran`.
  R CMD check returned OK, including examples, donttest examples, and tests. 
  The sanitizer reported loads and stores with insufficient space for a double 
  in the bundled NLopt source file bobyqa.c and not in this package.
- Fedora Linux 44 (x86_64), R-devel r90971 and r90985, GCC 16.2.1, Valgrind 3.27.1,
  `--no-manual --use-valgrind --no-stop-on-test-error`. Examples and tests
  passed.
