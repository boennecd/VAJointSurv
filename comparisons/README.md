# Comparisons with Other Packages
This directory contains examples where this package is compared with the JM, 
JMbayes, and joineRML package. The 
[`comp-simple.R`](comp-simple.R) contains a simple example
where the baseline hazard is from a Weibull model and all other time-varying 
effects are polynomials. The output is shown in the 
[`comp-simple.Rout`](comp-simple.Rout) file and the plots are in 
[`comp-simple.pdf`](comp-simple.pdf).

The 
[`comp.R`](comp.R) contains a similar example but where the
time-varying effects are splines. The output is shown in the 
[`comp.Rout`](comp.Rout) file  and the plots are in 
[`comp.pdf`](comp.pdf).

In both cases, the main output is at the end of the files. The output is created 
by calling the following commands on Ubuntu

```bash
R CMD BATCH --no-save --no-restore comp-simple.R
mv Rplots.pdf comp-simple.pdf
R CMD BATCH --no-save --no-restore comp.R
mv Rplots.pdf comp.pdf
```
