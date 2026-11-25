## New version submission

CRAN archived VAJointSurv on 2026-02-15 because SimSurvNMarker was unavailable.
SimSurvNMarker 0.1.4 is now on CRAN, and VAJointSurv 0.1.2 removes the R 4.1
pipe from the plot_surv() example.

## Check environments and results

- Debian 13 (x86_64), R 4.5.0 with SimSurvNMarker 0.1.4, `--as-cran`.
  Zero errors and zero warnings. The incoming feasibility NOTE identifies
  the archived package as a new submission. A second NOTE reported that
  the checker could not verify the current time.
- With `--no-manual --as-cran`, checks passed on Ubuntu 26.04.1 LTS
  (x86_64, R 4.5.3, R 4.6.2, and R-devel r90971), macOS 26.7.2
  (arm64, R 4.6.2), and Windows Server 2022 (x86_64, R 4.6.2).
- Ubuntu 26.04.1 LTS (x86_64), R-devel r90971, Clang 23.1.6,
  `--no-manual --no-stop-on-test-error`. Separate AddressSanitizer and
  UndefinedBehaviorSanitizer checks passed, including examples and tests.
- Fedora Linux 44 (x86_64), R-devel r90959, GCC 16.2.1 with
  AddressSanitizer. Package installation stopped during lazy loading because
  the nloptr binary could not load libRlapack.so. The examples and tests
  did not run.
- Fedora Linux 44 (x86_64), R-devel r90971, GCC 16.2.1, Valgrind 3.27.1,
  `--no-manual --use-valgrind --no-stop-on-test-error`. Examples and tests
  passed. Memcheck reported one possibly lost allocation in each process
  (400 bytes during examples and 416 bytes during tests). Both traces
  point to thread-local storage allocated at OpenMP worker creation.
  Memcheck reported zero bytes definitely or indirectly lost.
