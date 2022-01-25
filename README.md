Joint Survival and Marker Models
================

[![R-CMD-check](https://github.com/boennecd/VAJointSurv/workflows/R-CMD-check/badge.svg)](https://github.com/boennecd/VAJointSurv/actions)

This package provides means of estimating joint models for a mixture of
different types of markers and different types of survival outcomes. The
models are estimated with Gaussian variational approximations (GVAs).
Description of the supported models and multiple examples are given in
this document. A presentation of the package is available
[here](https://rpubs.com/boennecd/MEB-VAJointSurv) which includes a
simulation study to investigate bias of the parameter estimates and
coverage of approximate likelihood ratio based confidence intervals.

A comparison with the JM, JMbayes, and joineRML is provided on Github at
[github.com/boennecd/VAJointSurv/tree/main/comparisons](https://github.com/boennecd/VAJointSurv/tree/main/comparisons).

## Installation

The package can be installed from Github by calling:

``` r
remotes::install_github("boennecd/VAJointSurv", build_vignettes = TRUE)
```

The package currently requires the development version of psqn. You can
install the development version by calling:

``` r
remotes::install_github("boennecd/psqn", build_vignettes = TRUE)
```

## The Model

We will start by covering the model, then cover some examples, and end
with details about the implementation.

### The Markers

We assume that there are
![L](https://render.githubusercontent.com/render/math?math=L "L")
observed markers at each time point. We let
![\\vec Y\_{ij}](https://render.githubusercontent.com/render/math?math=%5Cvec%20Y_%7Bij%7D "\vec Y_{ij}")
denote the
![j](https://render.githubusercontent.com/render/math?math=j "j")th
observed markers for individual
![i](https://render.githubusercontent.com/render/math?math=i "i") and
let
![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D "s_{ij}")
be the observation time. The model for the markers is

![\\begin{align\*}
\\vec Y\_{ij}&= \\vec\\mu\_i(s\_{ij}, \\vec U\_i) + \\vec\\epsilon\_{ij} \\\\
\\epsilon\_{ij} &\\sim N^{(L)}(\\vec 0, \\Sigma) \\\\
\\mu\_{i1}(s, \\vec U\_i) &= \\vec x\_{i1}^\\top\\vec\\gamma\_1 + \\vec g\_1(s)^\\top\\vec\\beta\_1 + 
  \\vec m\_1(s)^\\top\\vec U\_{i1} \\\\
\\vdots &\\hphantom{=}\\vdots\\\\\\
\\mu\_{iL}(s, \\vec U\_i) &= \\vec x\_{iL}^\\top\\vec\\gamma\_L + \\vec g\_L(s)^\\top\\vec\\beta\_L + 
  \\vec m\_L(s)^\\top\\vec U\_{iL} \\\\
\\vec U\_i  &= \\begin{pmatrix}
    \\vec U\_{i1} \\\\ \\vdots \\\\ \\vec U\_{iL}
  \\end{pmatrix}\\sim N^{(R)}(\\vec0, \\Psi) 
\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%0A%5Cvec%20Y_%7Bij%7D%26%3D%20%5Cvec%5Cmu_i%28s_%7Bij%7D%2C%20%5Cvec%20U_i%29%20%2B%20%5Cvec%5Cepsilon_%7Bij%7D%20%5C%5C%0A%5Cepsilon_%7Bij%7D%20%26%5Csim%20N%5E%7B%28L%29%7D%28%5Cvec%200%2C%20%5CSigma%29%20%5C%5C%0A%5Cmu_%7Bi1%7D%28s%2C%20%5Cvec%20U_i%29%20%26%3D%20%5Cvec%20x_%7Bi1%7D%5E%5Ctop%5Cvec%5Cgamma_1%20%2B%20%5Cvec%20g_1%28s%29%5E%5Ctop%5Cvec%5Cbeta_1%20%2B%20%0A%20%20%5Cvec%20m_1%28s%29%5E%5Ctop%5Cvec%20U_%7Bi1%7D%20%5C%5C%0A%5Cvdots%20%26%5Chphantom%7B%3D%7D%5Cvdots%5C%5C%5C%0A%5Cmu_%7BiL%7D%28s%2C%20%5Cvec%20U_i%29%20%26%3D%20%5Cvec%20x_%7BiL%7D%5E%5Ctop%5Cvec%5Cgamma_L%20%2B%20%5Cvec%20g_L%28s%29%5E%5Ctop%5Cvec%5Cbeta_L%20%2B%20%0A%20%20%5Cvec%20m_L%28s%29%5E%5Ctop%5Cvec%20U_%7BiL%7D%20%5C%5C%0A%5Cvec%20U_i%20%20%26%3D%20%5Cbegin%7Bpmatrix%7D%0A%20%20%20%20%5Cvec%20U_%7Bi1%7D%20%5C%5C%20%5Cvdots%20%5C%5C%20%5Cvec%20U_%7BiL%7D%0A%20%20%5Cend%7Bpmatrix%7D%5Csim%20N%5E%7B%28R%29%7D%28%5Cvec0%2C%20%5CPsi%29%20%0A%5Cend%7Balign%2A%7D "\begin{align*}
\vec Y_{ij}&= \vec\mu_i(s_{ij}, \vec U_i) + \vec\epsilon_{ij} \\
\epsilon_{ij} &\sim N^{(L)}(\vec 0, \Sigma) \\
\mu_{i1}(s, \vec U_i) &= \vec x_{i1}^\top\vec\gamma_1 + \vec g_1(s)^\top\vec\beta_1 + 
  \vec m_1(s)^\top\vec U_{i1} \\
\vdots &\hphantom{=}\vdots\\\
\mu_{iL}(s, \vec U_i) &= \vec x_{iL}^\top\vec\gamma_L + \vec g_L(s)^\top\vec\beta_L + 
  \vec m_L(s)^\top\vec U_{iL} \\
\vec U_i  &= \begin{pmatrix}
    \vec U_{i1} \\ \vdots \\ \vec U_{iL}
  \end{pmatrix}\sim N^{(R)}(\vec0, \Psi) 
\end{align*}")

where
![\\vec\\epsilon\_{ij}](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cepsilon_%7Bij%7D "\vec\epsilon_{ij}")
is an error term of the observation,
![\\vec U\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%20U_i "\vec U_i")
is the ![R](https://render.githubusercontent.com/render/math?math=R "R")
dimensional unobserved random effect for individual
![i](https://render.githubusercontent.com/render/math?math=i "i"), and
the
![\\vec g\_k](https://render.githubusercontent.com/render/math?math=%5Cvec%20g_k "\vec g_k")s
and
![\\vec m\_k](https://render.githubusercontent.com/render/math?math=%5Cvec%20m_k "\vec m_k")s
are basis expansions. We can write the model more compactly by letting

![G(s) = \\begin{pmatrix} 
  \\vec g\_1(s)^\\top & 0^\\top & \\cdots & \\vec 0^\\top \\\\
  \\vec 0^\\top & \\vec g\_2(s)^\\top & \\ddots & \\vdots \\\\
  \\vdots & \\ddots & \\ddots & \\vec 0^\\top \\\\
  \\vec 0^\\top & \\cdots & \\vec 0^\\top & \\vec g\_L(s)^\\top \\end{pmatrix}](https://render.githubusercontent.com/render/math?math=G%28s%29%20%3D%20%5Cbegin%7Bpmatrix%7D%20%0A%20%20%5Cvec%20g_1%28s%29%5E%5Ctop%20%26%200%5E%5Ctop%20%26%20%5Ccdots%20%26%20%5Cvec%200%5E%5Ctop%20%5C%5C%0A%20%20%5Cvec%200%5E%5Ctop%20%26%20%5Cvec%20g_2%28s%29%5E%5Ctop%20%26%20%5Cddots%20%26%20%5Cvdots%20%5C%5C%0A%20%20%5Cvdots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cvec%200%5E%5Ctop%20%5C%5C%0A%20%20%5Cvec%200%5E%5Ctop%20%26%20%5Ccdots%20%26%20%5Cvec%200%5E%5Ctop%20%26%20%5Cvec%20g_L%28s%29%5E%5Ctop%20%5Cend%7Bpmatrix%7D "G(s) = \begin{pmatrix} 
  \vec g_1(s)^\top & 0^\top & \cdots & \vec 0^\top \\
  \vec 0^\top & \vec g_2(s)^\top & \ddots & \vdots \\
  \vdots & \ddots & \ddots & \vec 0^\top \\
  \vec 0^\top & \cdots & \vec 0^\top & \vec g_L(s)^\top \end{pmatrix}")

and defining
![M(s)](https://render.githubusercontent.com/render/math?math=M%28s%29 "M(s)")
and
![X\_i](https://render.githubusercontent.com/render/math?math=X_i "X_i")
similarly. Furthermore, let
![\\vec\\gamma = (\\vec\\gamma\_1^\\top,\\dots,\\vec\\gamma\_L^\\top)](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cgamma%20%3D%20%28%5Cvec%5Cgamma_1%5E%5Ctop%2C%5Cdots%2C%5Cvec%5Cgamma_L%5E%5Ctop%29 "\vec\gamma = (\vec\gamma_1^\top,\dots,\vec\gamma_L^\top)")
and define
![\\vec\\beta](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cbeta "\vec\beta")
similarly. Then the conditional mean vector at time
![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D "s_{ij}")
is

![\\vec\\mu\_i(s\_{ij}, \\vec U\_i) = X\_i\\vec\\gamma + G(s\_{ij})\\vec\\beta + M(s\_{ij})\\vec U\_i.](https://render.githubusercontent.com/render/math?math=%5Cvec%5Cmu_i%28s_%7Bij%7D%2C%20%5Cvec%20U_i%29%20%3D%20X_i%5Cvec%5Cgamma%20%2B%20G%28s_%7Bij%7D%29%5Cvec%5Cbeta%20%2B%20M%28s_%7Bij%7D%29%5Cvec%20U_i. "\vec\mu_i(s_{ij}, \vec U_i) = X_i\vec\gamma + G(s_{ij})\vec\beta + M(s_{ij})\vec U_i.")

### The Survival Outcomes

We assume that there are
![H](https://render.githubusercontent.com/render/math?math=H "H")
different types of survival outcomes. The conditional hazards of the
survival outcomes are

![\\begin{align\*}
h\_{i1}(t\\mid \\vec U\_i, \\vec\\xi\_i) &= \\exp\\left( 
  \\vec z\_{i1}^\\top\\vec \\delta\_1 + 
  \\omega\_1^\\top \\vec b\_1(t) + 
  \\vec\\alpha\_1^\\top M(t)\\vec U\_i + \\xi\_{i1}
  \\right) \\\\
\\vdots &\\hphantom{=}\\vdots \\\\
h\_{iH}(t\\mid \\vec U\_i,\\vec \\xi\_i) &= \\exp\\left(  
  \\vec z\_{iL}^\\top\\vec \\delta\_H + 
  \\omega\_H^\\top\\vec b\_H(t) + 
  \\vec\\alpha\_H^\\top M(t)\\vec U\_i + \\xi\_{iH}
  \\right) \\\\
\\vec\\xi\_i = \\begin{pmatrix}\\xi\_{i1} \\\\ \\vdots \\\\ \\xi\_{iH}\\end{pmatrix}
  &\\sim N^{(H)}(\\vec 0, \\Xi).
\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%0Ah_%7Bi1%7D%28t%5Cmid%20%5Cvec%20U_i%2C%20%5Cvec%5Cxi_i%29%20%26%3D%20%5Cexp%5Cleft%28%20%0A%20%20%5Cvec%20z_%7Bi1%7D%5E%5Ctop%5Cvec%20%5Cdelta_1%20%2B%20%0A%20%20%5Comega_1%5E%5Ctop%20%5Cvec%20b_1%28t%29%20%2B%20%0A%20%20%5Cvec%5Calpha_1%5E%5Ctop%20M%28t%29%5Cvec%20U_i%20%2B%20%5Cxi_%7Bi1%7D%0A%20%20%5Cright%29%20%5C%5C%0A%5Cvdots%20%26%5Chphantom%7B%3D%7D%5Cvdots%20%5C%5C%0Ah_%7BiH%7D%28t%5Cmid%20%5Cvec%20U_i%2C%5Cvec%20%5Cxi_i%29%20%26%3D%20%5Cexp%5Cleft%28%20%20%0A%20%20%5Cvec%20z_%7BiL%7D%5E%5Ctop%5Cvec%20%5Cdelta_H%20%2B%20%0A%20%20%5Comega_H%5E%5Ctop%5Cvec%20b_H%28t%29%20%2B%20%0A%20%20%5Cvec%5Calpha_H%5E%5Ctop%20M%28t%29%5Cvec%20U_i%20%2B%20%5Cxi_%7BiH%7D%0A%20%20%5Cright%29%20%5C%5C%0A%5Cvec%5Cxi_i%20%3D%20%5Cbegin%7Bpmatrix%7D%5Cxi_%7Bi1%7D%20%5C%5C%20%5Cvdots%20%5C%5C%20%5Cxi_%7BiH%7D%5Cend%7Bpmatrix%7D%0A%20%20%26%5Csim%20N%5E%7B%28H%29%7D%28%5Cvec%200%2C%20%5CXi%29.%0A%5Cend%7Balign%2A%7D "\begin{align*}
h_{i1}(t\mid \vec U_i, \vec\xi_i) &= \exp\left( 
  \vec z_{i1}^\top\vec \delta_1 + 
  \omega_1^\top \vec b_1(t) + 
  \vec\alpha_1^\top M(t)\vec U_i + \xi_{i1}
  \right) \\
\vdots &\hphantom{=}\vdots \\
h_{iH}(t\mid \vec U_i,\vec \xi_i) &= \exp\left(  
  \vec z_{iL}^\top\vec \delta_H + 
  \omega_H^\top\vec b_H(t) + 
  \vec\alpha_H^\top M(t)\vec U_i + \xi_{iH}
  \right) \\
\vec\xi_i = \begin{pmatrix}\xi_{i1} \\ \vdots \\ \xi_{iH}\end{pmatrix}
  &\sim N^{(H)}(\vec 0, \Xi).
\end{align*}")

The
![\\vec\\alpha](https://render.githubusercontent.com/render/math?math=%5Cvec%5Calpha "\vec\alpha")s
are association parameters which makes the markers and the survival
outcomes marginally dependent. The
![\\exp\\xi\_{ih}](https://render.githubusercontent.com/render/math?math=%5Cexp%5Cxi_%7Bih%7D "\exp\xi_{ih}")s
are frailty effects which makes the survival outcomes marginally
dependent even if
![\\vec\\alpha = \\vec 0](https://render.githubusercontent.com/render/math?math=%5Cvec%5Calpha%20%3D%20%5Cvec%200 "\vec\alpha = \vec 0").
By default, these are not included as it is assumed that events are
often terminal in which case identifying the frailty may be hard and
evidence of a frailty is both evidence of non-proportional hazards and
heterogeneity.

The observation process, the
![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D "s_{ij}")s,
can be modeled as one of
![H](https://render.githubusercontent.com/render/math?math=H "H")
different types of survival process. This can be done by adding each
![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D "s_{ij}")
as a left-truncated outcome of the given type of process except for
![s\_{i1}](https://render.githubusercontent.com/render/math?math=s_%7Bi1%7D "s_{i1}")
which is not left-truncated.

#### More General Survival Sub-model

The conditional hazards shown above only allows the hazard to be
associated with the population deviation for the mean marker through the
current value

![M(t)\\vec U\_i = \\begin{pmatrix}
    \\vec m\_1(t)^\\top\\vec U\_{i1} \\\\
    \\vdots \\\\
    \\vec m\_L(t)^\\top\\vec U\_{i1}
  \\end{pmatrix}.](https://render.githubusercontent.com/render/math?math=M%28t%29%5Cvec%20U_i%20%3D%20%5Cbegin%7Bpmatrix%7D%0A%20%20%20%20%5Cvec%20m_1%28t%29%5E%5Ctop%5Cvec%20U_%7Bi1%7D%20%5C%5C%0A%20%20%20%20%5Cvdots%20%5C%5C%0A%20%20%20%20%5Cvec%20m_L%28t%29%5E%5Ctop%5Cvec%20U_%7Bi1%7D%0A%20%20%5Cend%7Bpmatrix%7D. "M(t)\vec U_i = \begin{pmatrix}
    \vec m_1(t)^\top\vec U_{i1} \\
    \vdots \\
    \vec m_L(t)^\top\vec U_{i1}
  \end{pmatrix}.")

However, the researcher may be interested in using the cumulative value,
![\\int^t \\vec m\_k(s)^\\top ds\\vec U\_{ik}](https://render.githubusercontent.com/render/math?math=%5Cint%5Et%20%5Cvec%20m_k%28s%29%5E%5Ctop%20ds%5Cvec%20U_%7Bik%7D "\int^t \vec m_k(s)^\top ds\vec U_{ik}"),
the derivative,
![\\partial\\vec m\_k(t)^\\top\\vec U\_{ik}/\\partial t](https://render.githubusercontent.com/render/math?math=%5Cpartial%5Cvec%20m_k%28t%29%5E%5Ctop%5Cvec%20U_%7Bik%7D%2F%5Cpartial%20t "\partial\vec m_k(t)^\top\vec U_{ik}/\partial t"),
or a combination of these and the current value. This is supported in
the package through the `ders` argument of `surv_term`. As an example,
if we let
![\\vec d\_k(t) = \\int^t \\vec m\_k(s) ds](https://render.githubusercontent.com/render/math?math=%5Cvec%20d_k%28t%29%20%3D%20%5Cint%5Et%20%5Cvec%20m_k%28s%29%20ds "\vec d_k(t) = \int^t \vec m_k(s) ds"),
![\\vec r\_k(t) = \\partial\\vec m\_k(t)/\\partial t](https://render.githubusercontent.com/render/math?math=%5Cvec%20r_k%28t%29%20%3D%20%5Cpartial%5Cvec%20m_k%28t%29%2F%5Cpartial%20t "\vec r_k(t) = \partial\vec m_k(t)/\partial t"),
and
![L = 3](https://render.githubusercontent.com/render/math?math=L%20%3D%203 "L = 3")
then it is possible to replace the
![\\vec\\alpha\_1^\\top M(t)\\vec U\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%5Calpha_1%5E%5Ctop%20M%28t%29%5Cvec%20U_i "\vec\alpha_1^\top M(t)\vec U_i")
in the hazard with

![\\vec\\alpha\_1^\\top\\begin{pmatrix}
  \\vec d\_1(t)^\\top \\vec U\_{i1} \\\\
  \\vec m\_1(t)^\\top \\vec U\_{i1} \\\\
  \\vec m\_2(t)^\\top \\vec U\_{i2} \\\\
  \\vec d\_3(t)^\\top \\vec U\_{i3} \\\\
  \\vec r\_3(t)^\\top \\vec U\_{i3} \\\\
\\end{pmatrix}](https://render.githubusercontent.com/render/math?math=%5Cvec%5Calpha_1%5E%5Ctop%5Cbegin%7Bpmatrix%7D%0A%20%20%5Cvec%20d_1%28t%29%5E%5Ctop%20%5Cvec%20U_%7Bi1%7D%20%5C%5C%0A%20%20%5Cvec%20m_1%28t%29%5E%5Ctop%20%5Cvec%20U_%7Bi1%7D%20%5C%5C%0A%20%20%5Cvec%20m_2%28t%29%5E%5Ctop%20%5Cvec%20U_%7Bi2%7D%20%5C%5C%0A%20%20%5Cvec%20d_3%28t%29%5E%5Ctop%20%5Cvec%20U_%7Bi3%7D%20%5C%5C%0A%20%20%5Cvec%20r_3%28t%29%5E%5Ctop%20%5Cvec%20U_%7Bi3%7D%20%5C%5C%0A%5Cend%7Bpmatrix%7D "\vec\alpha_1^\top\begin{pmatrix}
  \vec d_1(t)^\top \vec U_{i1} \\
  \vec m_1(t)^\top \vec U_{i1} \\
  \vec m_2(t)^\top \vec U_{i2} \\
  \vec d_3(t)^\top \vec U_{i3} \\
  \vec r_3(t)^\top \vec U_{i3} \\
\end{pmatrix}")

That is, a cumulative and current value for the first marker, only the
current value for the second marker, and the cumulative and derivative
for the third marker.

### Model Support

The package can handle

-   Mixed type of basis expansions in
    ![G](https://render.githubusercontent.com/render/math?math=G "G"),
    ![M](https://render.githubusercontent.com/render/math?math=M "M"),
    and
    ![b\_1,\\dots,b\_H](https://render.githubusercontent.com/render/math?math=b_1%2C%5Cdots%2Cb_H "b_1,\dots,b_H").
-   A maximum
    ![L = 31](https://render.githubusercontent.com/render/math?math=L%20%3D%2031 "L = 31")
    different markers, and an arbitrary number survival outcome types,
    ![H](https://render.githubusercontent.com/render/math?math=H "H")
    (though, identification of parameters quickly becomes an issue).
-   Only observing a subset of the
    ![L](https://render.githubusercontent.com/render/math?math=L "L")
    markers at a given point in time
    ![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D "s_{ij}")
    (i.e. some may be missing).
-   Any number of observed markers and survival outcomes for each
    individual.
-   Left-truncation and right-censoring.

The following is not supported

-   Interval-censoring and left-censoring are not supported although
    this is not a complicated extension.
-   There cannot be multiple observed markers of the same type at each
    point
    ![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D "s_{ij}")
    for each individual
    ![i](https://render.githubusercontent.com/render/math?math=i "i").
-   There cannot be time-varying effects for a covariate. This is not
    too complicated to implement but not yet implemented.
-   Competing events, and multi-state models are not supported.

## Examples

We show a few examples of using this package and illustrate its API in
this section. The examples includes:

-   Only one marker in the [Univariate Marker
    Example](#univariate-marker-example) section where we can compare
    with lme4. The section also includes some remarks about the
    implemented expansions.
-   Two markers in the [Two Markers](#two-markers) section.
-   A recurrent event in the [Recurrent Event](#recurrent-event)
    section. It is also shown how to change the quadrature rule that is
    used to approximate the approximate expected cumulative hazard.
-   An example with two markers and a recurrent event is given in the
    [Two Markers and a Recurrent
    Event](#two-markers-and-a-recurrent-event) section. The section also
    includes an example of how to construct approximate profile
    likelihood based confidence intervals and how to obtain Wald type
    tests and confidence intervals from the observed information matrix.
-   A similar example is provided but without the frailty in the [Two
    Markers and a Recurrent Event Without
    Frailty](#two-markers-and-a-recurrent-event-without-frailty)
    section. Notice that this is the default with `surv_term`. It is
    also shown how one can get the variational parameters for each
    individual and change the starting values of the covariance
    matrices.
-   An example with two markers, the observation process, and a terminal
    event is given in the [Two Markers, the Observation Time Process,
    and a Terminal
    Event](#two-markers-the-observation-time-process-and-a-terminal-event)
    section. The observation time process and the markers are not
    marginally independent in the true model. Thus, they need to be
    modeled together or the dependence accounted for in some other way.
    Examples and comments regarding caching of the expansions are also
    provided.
-   A similar example is provided in the [Two Markers, the Observation
    Time Process, and a Terminal Event with Delayed
    Entry](#two-markers-the-observation-time-process-and-a-terminal-event-with-delayed-entry)
    section but where there is delayed entry.
-   A similar example is provided in the [Two Markers, the Observation
    Time Process, a Terminal Event and Mixed
    Dependencies](#two-markers-the-observation-time-process-a-terminal-event-and-mixed-dependencies)
    section but where a mixture of the cumulative, the present value,
    and the derivative is used in the sub-survival models.

We start first with some general remarks.

### General Remarks

There is not intercept by default in the splines. Usually, you almost
always want one in the both the random effect and fixed effect part. The
intercept in included by default in `formula` (unless you do
`. ~ ... - 1` and there is not factor variable). Thus, you do not want
an intercept in the `time_fixef` in both `marker_term` and `surv_term`.
On the other hand, you do need yourself to include the intercept in the
`time_rng` argument of `marker_term`. Thus, a typical call may look like

``` r
mark <- marker_term(
  y ~ X1 + X2, # there is an intercept by default
  # `intercept` is FALSE by default 
  time_fixef = poly_term(time, degree = 3), 
  # we need to add the intercept
  time_rng = poly_term(time, degree = 1, intercept = TRUE))

s_term <- surv_term(
  Surv(tstop, event) ~ X1, # there is an intercept by default
  # `intercept` is FALSE by default
  time_fixef = poly_term(tstop, degree = 1))
```

You may run into computational issues if you have a very flexible random
effect model. This is quite typical with mixed models in a frequentist
setting also with other estimation methods. One outcome of a too
flexible random effect model is slow convergence or a `NaN` output
because some of the scale matrices do not have full rank.

### Univariate Marker Example

We make a comparison below with one marker only with the lme4 package.
We do this as we can check that

-   We get almost the same (the lower bound can be equal to the log
    marginal likelihood in this case).
-   The lower bound is less than or equal to the log marginal
    likelihood.

lme4 should be used in this case as the linear model is tractable and
this is used by the lme4 package. The computation time comparison is not
fair as we use four threads with the methods in this package.

``` r
# settings for the simulation
library(splines)
g_func <- function(x)
  ns(x, knots = c(.5, 1.5), Boundary.knots = c(0, 2))
m_func <- function(x)
  ns(x, knots = 1, Boundary.knots = c(0, 2), intercept = TRUE)

fixef_vary_marker <- c(1.4, -1.2, -2.1) # beta
fixef_marker <- c(-.5, 1) # gamma

# Psi 
vcov_vary <- structure(
  c(0.18, 0.05, -0.05, 0.05, 0.34, -0.25, -0.05, -0.25, 0.24), 
  .Dim = c(3L, 3L))
vcov_marker <- matrix(.6^2, 1) # Sigma

# plot the true mean curve along with 95% confidence pointwise quantiles
library(VAJointSurv)
par(mar = c(5, 5, 1, 1))
plot_marker(
  time_fixef = ns_term(
    knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
  time_rng = ns_term(
    knots = 1, Boundary.knots = c(0, 2), intercept = TRUE), 
  fixef_vary = fixef_vary_marker, x_range = c(0, 2), vcov_vary = vcov_vary)
```

![](man/figures/README-univariate-1.png)<!-- -->

``` r
library(mvtnorm)
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  sig_sqrt <- sqrt(vcov_marker)
  dat <- lapply(1:n_ids, function(id){
    # sample the number of outcomes and the fixed effect covariates
    n_obs <- sample.int(8L, 1L)
    obs_time <- sort(runif(n_obs, 0, 2))
    X <- cbind(1, rnorm(n_obs))
    colnames(X) <- paste0("X", 1:NCOL(X) - 1L)
    
    # sample the outcomes
    U <- drop(rmvnorm(1, sigma = vcov_vary))
    eta <- X %*% fixef_marker + drop(g_func(obs_time)) %*% fixef_vary_marker + 
      drop(m_func(obs_time)) %*% U
    y <- drop(eta) + rnorm(n_obs, sd = sig_sqrt)
    
    cbind(y = y, X = X[, -1, drop = FALSE], time = obs_time, id = id)
  })
  
  # combine the data and return
  out <- as.data.frame(do.call(rbind, dat))
  out$id <- as.integer(out$id)
  out[sample.int(NROW(out)), ] # the order does not matter
}

# sample a moderately large data set
set.seed(2)
dat <- sim_dat(2000L)

# the number of observed outcomes
nrow(dat)
#> [1] 8992

# example of one individual's outcomes
subset(dat, id == 1)
#>         y       X1   time id
#> 5 -3.6206 -0.23970 1.8877  1
#> 2  0.1270 -0.08025 1.1467  1
#> 3 -0.3545  0.13242 1.4047  1
#> 1 -1.1089 -1.13038 0.3361  1
#> 4 -0.2127  0.70795 1.8869  1

# fit the model with lme4
library(lme4)
#> Loading required package: Matrix
system.time(
  fit <- lmer(y ~ X1 + g_func(time) + (m_func(time) - 1| id), dat, 
              control = lmerControl(optimizer = "bobyqa"),
              # to compare with the lower bound from this package
              REML = FALSE))
#>    user  system elapsed 
#>   0.864   0.004   0.868

# the maximum log likelihood
print(logLik(fit), digits = 8)
#> 'log Lik.' -9161.7723 (df=12)

# estimate the model with this package. Get the object we need for the 
# optimization
system.time(comp_obj <- joint_ms_ptr(
  markers = marker_term(
    y ~ X1, id = id, dat, 
    time_fixef = ns_term(time, knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
    time_rng = ns_term(time, knots = 1, Boundary.knots = c(0, 2), 
                       intercept = TRUE)),
  max_threads = 4L))
#>    user  system elapsed 
#>   0.015   0.000   0.016

# get the starting values
system.time(start_val <- joint_ms_start_val(comp_obj))
#>    user  system elapsed 
#>   2.675   0.019   1.176

# lower bound at the starting values
print(-attr(start_val, "value"), digits = 8)
#> [1] -9161.7803

# check that the gradient is correct
f <- function(x){
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(start_val, 12 + 2 * 9)), 
          head(joint_ms_lb_gr(comp_obj, start_val), 12 + 2 * 9))
#> [1] "Mean relative difference: 6.563e-06"

# find the maximum lower bound estimate
system.time(opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 1000L,
                                    cg_tol = .2, c2 = .1))
#>    user  system elapsed 
#>   1.207   0.000   0.304

# check the gradient norm. We may need to reduce the convergence tolerance if 
# this is not small. In can also be a sign of convergence issues
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.665
opt_out$info # convergence code (0 == 'OK')
#> [1] 0
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -9161.7763
print(logLik(fit), digits = 8) # maximum likelihood from lme4
#> 'log Lik.' -9161.7723 (df=12)

# compare the estimates with the actual values. Start with the fixed effects
fmt_ests <- joint_ms_format(comp_obj, opt_out$par)
rbind(lme4 = fixef(fit),
      `This package` = do.call(c, fmt_ests$markers[[1]]), 
      Actual = c(fixef_marker, fixef_vary_marker))
#>              (Intercept)    X1 g_func(time)1 g_func(time)2 g_func(time)3
#> lme4             -0.5206 1.007         1.372        -1.104        -2.089
#> This package     -0.5206 1.007         1.372        -1.104        -2.089
#> Actual           -0.5000 1.000         1.400        -1.200        -2.100

# then the random effect covariance matrix
vcov_est <- VarCorr(fit)[["id"]]
attributes(vcov_est)[c("stddev", "correlation")] <- NULL
vcov_est # lme4
#>               m_func(time)1 m_func(time)2 m_func(time)3
#> m_func(time)1       0.16228       0.09634      -0.04534
#> m_func(time)2       0.09634       0.29559      -0.26821
#> m_func(time)3      -0.04534      -0.26821       0.29652
fmt_ests$vcov$vcov_vary # this package
#>          [,1]     [,2]     [,3]
#> [1,]  0.16578  0.09519 -0.04554
#> [2,]  0.09519  0.29598 -0.26831
#> [3,] -0.04554 -0.26831  0.29718
vcov_vary # the actual values
#>       [,1]  [,2]  [,3]
#> [1,]  0.18  0.05 -0.05
#> [2,]  0.05  0.34 -0.25
#> [3,] -0.05 -0.25  0.24

# the error term variance
c(lme4 = attr(VarCorr(fit), "sc")^2, 
  `This package` = fmt_ests$vcov$vcov_marker, 
  Actual = vcov_marker)
#>         lme4 This package       Actual 
#>       0.3667       0.3664       0.3600
```

#### Note on Basis Expansions

This is a technical section which you may skip. Special basis expansions
for e.g. `poly` and `ns` are implemented in R as `poly_term` and
`ns_term`. The reason is that the integration shown in the [Survival
Outcomes](#survival-outcomes) section has to be done many times. Thus,
we have implemented all the expansions in C++ to be used in C++ to
reduce the cost of the evaluations.

The `_term` functions provide the R interface to the C++ functions. Each
return a list with elements like their R counterparts and an `eval`
element. The latter element can be used to evaluate the basis, the
derivatives of the basis if the `der` argument is greater than zero, or
the integral

![\\vec d(u) = \\int\_l^u \\vec m(s) ds](https://render.githubusercontent.com/render/math?math=%5Cvec%20d%28u%29%20%3D%20%5Cint_l%5Eu%20%5Cvec%20m%28s%29%20ds "\vec d(u) = \int_l^u \vec m(s) ds")

if the `der` argument is equal to minus one. The lower limit is passed
with the `lower_limit` argument to the `eval` function. We illustrate
this below first with a raw polynomial. We also show that the
derivatives and integral are correct.

``` r
# raw polynomial with an without an intercept
pl <- poly_term(degree = 3, raw = TRUE, intercept = TRUE)
pl_no_inter <- poly_term(degree = 3, raw = TRUE, intercept = FALSE)
t(pl         $eval(1:2))
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    1    1    1
#> [2,]    1    2    4    8
t(pl_no_inter$eval(1:2))
#>      [,1] [,2] [,3]
#> [1,]    1    1    1
#> [2,]    2    4    8
outer(1:2, 0:3, `^`) # R version
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    1    1    1
#> [2,]    1    2    4    8

# derivatives
library(numDeriv)
t(pl$eval(2:3, der = 1))
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    4   12
#> [2,]    0    1    6   27
t(pl_no_inter$eval(2:3, der = 1))
#>      [,1] [,2] [,3]
#> [1,]    1    4   12
#> [2,]    1    6   27
t(pl$eval(2:3, der = 2)) # second derivative
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    2   12
#> [2,]    0    0    2   18
t(pl_no_inter$eval(2:3, der = 2)) # second derivative
#>      [,1] [,2] [,3]
#> [1,]    0    2   12
#> [2,]    0    2   18
# trivial to do compute by hand but we do it with numDeriv anyway
t(sapply(2:3, function(x) jacobian(function(z) outer(z, 0:3, `^`), x)))
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    4   12
#> [2,]    0    1    6   27

# integral
t(pl$eval(3:4, der = -1, lower_limit = 0))
#>      [,1] [,2]  [,3]  [,4]
#> [1,]    3  4.5  9.00 20.25
#> [2,]    4  8.0 21.33 64.00
t(pl_no_inter$eval(3:4, der = -1, lower_limit = 0))
#>      [,1]  [,2]  [,3]
#> [1,]  4.5  9.00 20.25
#> [2,]  8.0 21.33 64.00
# we could do this by hand but we write a function which we can reuse
do_int <- function(xs, f, lower_limit){
  n_terms <- length(f(xs[1]))
  g <- function(x, i) f(x)[, i]
  outer(xs, 1:n_terms, Vectorize(
    function(x, i) integrate(g, lower_limit, x, i = i)$value))
}
do_int(3:4, function(x) outer(x, 0:3, `^`), 0)
#>      [,1] [,2]  [,3]  [,4]
#> [1,]    3  4.5  9.00 20.25
#> [2,]    4  8.0 21.33 64.00
```

The functionality also works for orthogonal polynomials.

``` r
# equally works with orthogonal polynomials
xx <- seq(0, 4, length.out = 10)
pl <- poly_term(xx, degree = 3, raw = FALSE, intercept = TRUE)
t(pl$eval(2:3))
#>      [,1]      [,2]    [,3]       [,4]
#> [1,]    1 5.500e-17 -0.3590  2.391e-16
#> [2,]    1 2.477e-01 -0.1387 -3.881e-01
predict(poly(xx, degree = 3), 2:3)
#>              1       2          3
#> [1,] 5.500e-17 -0.3590  2.391e-16
#> [2,] 2.477e-01 -0.1387 -3.881e-01

# derivatives
t(pl$eval(2:3, der = 1))
#>      [,1]   [,2]      [,3]     [,4]
#> [1,]    0 0.2477 1.110e-16 -0.59310
#> [2,]    0 0.2477 4.406e-01  0.02176
t(sapply(2:3, function(x) 
  jacobian(function(z) predict(poly(xx, degree = 3), z), x)))
#>        [,1]   [,2]     [,3]
#> [1,] 0.2477 0.0000 -0.59310
#> [2,] 0.2477 0.4406  0.02176

# integral 
t(pl$eval(3:4, der = -1, lower_limit = 0))
#>      [,1]       [,2]    [,3]      [,4]
#> [1,]    3 -3.716e-01 -0.4162 1.211e-01
#> [2,]    4  2.220e-16 -0.2611 1.776e-15
do_int(3:4, function(x) predict(poly(xx, degree = 3), x), 0)
#>            [,1]    [,2]      [,3]
#> [1,] -3.716e-01 -0.4162 1.211e-01
#> [2,]  2.248e-16 -0.2611 5.231e-16
```

The derivatives and integral are also available for B-splines and
natural cubic splines through the `ns_term` and `bs_term` functions.

``` r
# B-splines
library(splines)
xx <- 0:10
b <- bs_term(xx, df = 4, Boundary.knots = range(xx))
t(b$eval(4.33))
#>        [,1]   [,2]   [,3] [,4]
#> [1,] 0.3598 0.4755 0.1624    0
c(predict(bs(xx, df = 4, Boundary.knots = range(xx)), 4.33))
#> [1] 0.3598 0.4755 0.1624 0.0000

# derivatives
t(b$eval(4.33, der = 1))
#>         [,1]    [,2]   [,3] [,4]
#> [1,] -0.1713 0.06963 0.1125    0
f <- function(z) suppressWarnings(
  predict(bs(xx, df = 4, Boundary.knots = range(xx)), z))
t(jacobian(f, 4.33))
#>         [,1]    [,2]   [,3] [,4]
#> [1,] -0.1713 0.06963 0.1125    0

# integrals
all.equal( t(b$eval(12, der = -1, lower_limit = 11)), do_int(12, f, 11), 1e-6)
#> [1] TRUE
all.equal( t(b$eval(11, der = -1, lower_limit =  0)), do_int(11, f,  0), 1e-6)
#> [1] TRUE
all.equal(-t(b$eval( 0, der = -1, lower_limit = 11)), do_int(11, f,  0), 1e-6)
#> [1] TRUE
all.equal( t(b$eval( 5, der = -1, lower_limit =  1)), do_int( 5, f,  1), 1e-6)
#> [1] TRUE
all.equal( t(b$eval(-1, der = -1, lower_limit =  2)), do_int(-1, f,  2), 1e-6)
#> [1] TRUE

# natural cubic spline
xx <- 0:10
n <- ns_term(xx, df = 3, Boundary.knots = range(xx))
t(n$eval(4.33))
#>        [,1]   [,2]   [,3]
#> [1,] 0.1723 0.5256 -0.346
c(predict(ns(xx, df = 3, Boundary.knots = range(xx)), 4.33))
#> [1]  0.1723  0.5256 -0.3460

# derivatives
t(n$eval(4.33, der = 1))
#>        [,1]     [,2]    [,3]
#> [1,] 0.2186 -0.05737 0.05165
f <- function(z) predict(ns(xx, df = 3, Boundary.knots = range(xx)), z)
t(jacobian(f, 4.33))
#>        [,1]     [,2]    [,3]
#> [1,] 0.2186 -0.05737 0.05165

# integrals
all.equal( t(n$eval(12, der = -1, lower_limit = 11)), do_int(12, f, 11), 1e-6)
#> [1] TRUE
all.equal( t(n$eval(11, der = -1, lower_limit =  0)), do_int(11, f,  0), 1e-6)
#> [1] TRUE
all.equal(-t(n$eval( 0, der = -1, lower_limit = 11)), do_int(11, f,  0), 1e-6)
#> [1] TRUE
all.equal( t(n$eval( 7, der = -1, lower_limit = -2)), do_int( 7, f, -2), 1e-6)
#> [1] TRUE
all.equal( t(n$eval(-1, der = -1, lower_limit = -3)), do_int(-1, f, -3), 1e-6)
#> [1] TRUE
```

Caution: the B-splines are only tested with `degree = 3`!

#### Log Transformations

The expansions can also be used on the log scale. This is mainly
implemented to allow the log of the baseline hazard to be parameterized
in terms of an expansions in log time to have the Weibull model as a
special case. A few examples are given below.

``` r
# raw log polynomial with an without an intercept
pl <- poly_term(degree = 3, raw = TRUE, intercept = TRUE, use_log = TRUE)
pl_no_inter <- poly_term(degree = 3, raw = TRUE, intercept = FALSE, 
                         use_log = TRUE)
t(pl         $eval(1:2))
#>      [,1]   [,2]   [,3]  [,4]
#> [1,]    1 0.0000 0.0000 0.000
#> [2,]    1 0.6931 0.4805 0.333
t(pl_no_inter$eval(1:2))
#>        [,1]   [,2]  [,3]
#> [1,] 0.0000 0.0000 0.000
#> [2,] 0.6931 0.4805 0.333
outer(log(1:2), 0:3, `^`) # R version
#>      [,1]   [,2]   [,3]  [,4]
#> [1,]    1 0.0000 0.0000 0.000
#> [2,]    1 0.6931 0.4805 0.333

# derivatives
library(numDeriv)
t(pl$eval(2:3, der = 1))
#>      [,1]   [,2]   [,3]   [,4]
#> [1,]    0 0.5000 0.6931 0.7207
#> [2,]    0 0.3333 0.7324 1.2069
t(pl_no_inter$eval(2:3, der = 1))
#>        [,1]   [,2]   [,3]
#> [1,] 0.5000 0.6931 0.7207
#> [2,] 0.3333 0.7324 1.2069
# trivial to do compute by hand but we do it with numDeriv anyway
t(sapply(2:3, function(x) jacobian(function(z) outer(log(z), 0:3, `^`), x)))
#>      [,1]   [,2]   [,3]   [,4]
#> [1,]    0 0.5000 0.6931 0.7207
#> [2,]    0 0.3333 0.7324 1.2069
```

``` r
# B-splines
library(splines)
xx <- 1:10
b <- bs_term(xx, df = 4, use_log = TRUE)
t(b$eval(4.33))
#>        [,1]   [,2]   [,3] [,4]
#> [1,] 0.1736 0.4746 0.3491    0
c(predict(bs(log(xx), df = 4), log(4.33)))
#> [1] 0.1736 0.4746 0.3491 0.0000

# derivatives
t(b$eval(4.33, der = 1))
#>         [,1]     [,2]  [,3] [,4]
#> [1,] -0.1223 -0.03495 0.165    0
f <- function(z) suppressWarnings(predict(bs(log(xx), df = 4), log(z)))
t(jacobian(f, 4.33))
#>         [,1]     [,2]  [,3] [,4]
#> [1,] -0.1223 -0.03495 0.165    0

# natural cubic spline
xx <- 1:10
n <- ns_term(xx, df = 3, use_log = TRUE)
t(n$eval(4.33))
#>        [,1]   [,2]   [,3]
#> [1,] 0.3855 0.4276 -0.307
c(predict(ns(log(xx), df = 3), log(4.33)))
#> [1]  0.3855  0.4276 -0.3070

# derivatives
t(n$eval(4.33, der = 1))
#>        [,1]     [,2]    [,3]
#> [1,] 0.2536 -0.09559 0.07548
f <- function(z) predict(ns(log(xx), df = 3), log(z))
t(jacobian(f, 4.33))
#>        [,1]     [,2]    [,3]
#> [1,] 0.2536 -0.09559 0.07548
```

### Two Markers

We provide an example with two markers instead of one. We observe one of
the two markers or both at each observation time. There is also a closed
form solution in this case. Thus, this is not the optimal way of doing
this and it is only shown as an example and as a sanity check.

``` r
# settings for the simulation
library(splines)
g_funcs <- list(
  function(x)
    ns(x, knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
  function(x)
    # a raw polynomial
    outer(x, 1:2, `^`))
m_funcs <- list(
  function(x)
    ns(x, knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
  function(x)
    # a raw polynomial
    outer(x, 0:1, `^`))

fixef_vary_marker <- list(c(1.4, -1.2, -2.1), c(.5, .67)) # beta
fixef_marker <- list(c(-.5, 1), .25) # gamma

# Psi
vcov_vary <- structure(
  c(0.35, 0.02, -0.05, 0.01, 0.02, 0.12, -0.06, -0.01, -0.05, -0.06, 0.32, 0.09, 0.01, -0.01, 0.09, 0.12),
  .Dim = c(4L, 4L))
vcov_marker <- matrix(c(.6^2, .1, .1, .4^2), 2) # Sigma

# plot the markers' mean curve
library(VAJointSurv)
par(mar = c(5, 5, 1, 1))
plot_marker(
  time_fixef = ns_term(
    knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
  time_rng = ns_term(
    knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
  fixef_vary = fixef_vary_marker[[1]], x_range = c(0, 2), 
  vcov_vary = vcov_vary[1:2, 1:2], ylab = "Marker 1")
```

![](man/figures/README-two_markers-1.png)<!-- -->

``` r
plot_marker(
  time_fixef = poly_term(degree = 2, raw = TRUE),
  time_rng = poly_term(degree = 1, raw = TRUE, intercept = TRUE), 
  fixef_vary = fixef_vary_marker[[2]], x_range = c(0, 2), 
  vcov_vary = vcov_vary[3:4, 3:4], ylab = "Marker 2")
```

![](man/figures/README-two_markers-2.png)<!-- -->

``` r
library(mvtnorm)
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  dat <- lapply(1:n_ids, function(id){
    # sample the number of outcomes and the fixed effect covariates
    n_obs <- sample.int(8L, 1L)
    obs_time <- sort(runif(n_obs, 0, 2))
    X1 <- cbind(1, rnorm(n_obs))
    X2 <- matrix(1, n_obs)
    colnames(X1) <- paste0("X1_", 1:NCOL(X1) - 1L)
    colnames(X2) <- paste0("X2_", 1:NCOL(X2) - 1L)
    X <- list(X1, X2)
    
    # sample the outcomes
    U <- drop(rmvnorm(1, sigma = vcov_vary))
    eta <- sapply(1:2, function(i)
      X[[i]] %*% fixef_marker[[i]] + 
      drop(g_funcs[[i]](obs_time)) %*% fixef_vary_marker[[i]] + 
      drop(m_funcs[[i]](obs_time)) %*% U[1:2 + (i == 2) * 2])
    
    ys <- eta + rmvnorm(n_obs, sigma = vcov_marker)
    colnames(ys) <- paste0("Y", 1:2)
    
    # mask some observations
    do_mask <- sample.int(3L, n_obs, replace = TRUE) 
    ys[do_mask == 2, 1] <- NA
    ys[do_mask == 3, 2] <- NA
    
    X <- do.call(cbind, lapply(X, function(x) x[, -1, drop = FALSE]))
    cbind(ys, X, time = obs_time, id = id)
  })
  
  # combine the data and return
  out <- as.data.frame(do.call(rbind, dat))
  out$id <- as.integer(out$id)
  out[sample.int(NROW(out)), ] # the order does not matter
}

# sample a moderately large data set
set.seed(2)
dat <- sim_dat(1000L)

# example of the data for one individual
subset(dat, id == 1)
#>       Y1     Y2     X1_1   time id
#> 3     NA 3.1915  0.13242 1.4047  1
#> 5     NA 4.3373 -0.23970 1.8877  1
#> 4     NA 4.8558  0.70795 1.8869  1
#> 2  1.399 1.7030 -0.08025 1.1467  1
#> 1 -1.286 0.4379 -1.13038 0.3361  1

# estimate the model with this package. Get the object we need for the 
# optimization
marker_1 <- marker_term(
    Y1 ~ X1_1, id = id, subset(dat, !is.na(Y1)), 
    time_fixef = ns_term(time, knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
    time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 2), 
                       intercept = TRUE))
marker_2 <- marker_term(
    Y2 ~ 1, id = id, subset(dat, !is.na(Y2)), 
    time_fixef = poly_term(time, degree = 2, raw = TRUE),
    time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

comp_obj <- joint_ms_ptr(markers = list(marker_1, marker_2), 
                         max_threads = 4L)
rm(marker_1, marker_2)

# get the starting values
system.time(start_val <- joint_ms_start_val(comp_obj))
#>    user  system elapsed 
#>   2.175   0.004   0.680

# lower bound at the starting values
print(-attr(start_val, "value"), digits = 8)
#> [1] -5729.7904

# check that the gradient is correct
f <- function(x){
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(start_val, 22 + 2 * 14)), 
          head(joint_ms_lb_gr(comp_obj, start_val), 22 + 2 * 14))
#> [1] "Mean relative difference: 5.618e-05"

# find the maximum lower bound estimate
system.time(opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 1000L, 
                                    pre_method = 3L, cg_tol = .2, c2 = .1, 
                                    gr_tol = .1))
#>    user  system elapsed 
#>   0.326   0.000   0.083

# we set gr_tol in the call so this is the convergence criterion for the 
# gradient
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.09769
opt_out$info # convergence code (0 == 'OK')
#> [1] 0
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -5729.7904
opt_out$counts
#> function gradient     n_cg 
#>       42       30       41

# find the maximum lower bound with lbfgs
library(lbfgsb3c)
system.time(lbfgs_res <- lbfgsb3c(
  start_val, function(x) joint_ms_lb(comp_obj, x),
  function(x) joint_ms_lb_gr(comp_obj, x), 
  control = list(factr = 1e-8, maxit = 10000L)))
#>    user  system elapsed 
#>  44.238   0.031  11.224
lbfgs_res$convergence # convergence code (0 == 'OK')
#> [1] 0
print(-lbfgs_res$value, digits = 8)  # maximum lower bound value
#> [1] -5729.7903
lbfgs_res$counts # may have hit maxit!
#> [1] 5064 5064

# compare the estimates with the actual values. Start with the fixed effects
fmt_ests <- joint_ms_format(comp_obj, opt_out$par)
fmt_ests_lbfgs <- joint_ms_format(comp_obj, lbfgs_res$par)

# the parameters for the first marker
fmt_ests$markers[[1]] 
#> $fixef
#> [1] -0.4967  1.0249
#> 
#> $fixef_vary
#> [1]  1.544 -1.237 -2.114
fmt_ests_lbfgs$markers[[1]] # with lbfgs
#> $fixef
#> [1] -0.4966  1.0249
#> 
#> $fixef_vary
#> [1]  1.544 -1.238 -2.114

fixef_marker[[1]] # true values
#> [1] -0.5  1.0
fixef_vary_marker[[1]] # true values
#> [1]  1.4 -1.2 -2.1

# the parameters for the second marker
fmt_ests$markers[[2]]
#> $fixef
#> [1] 0.2552
#> 
#> $fixef_vary
#> [1] 0.5381 0.6592
fmt_ests_lbfgs$markers[[2]] # with lbfgs
#> $fixef
#> [1] 0.2552
#> 
#> $fixef_vary
#> [1] 0.5382 0.6591
fixef_marker[[2]] # true values
#> [1] 0.25
fixef_vary_marker[[2]] # true values
#> [1] 0.50 0.67

# the parameters for covariance matrix of the random effects
fmt_ests$vcov$vcov_vary
#>            [,1]     [,2]      [,3]      [,4]
#> [1,]  0.4155040  0.02347 -0.006174 0.0004588
#> [2,]  0.0234666  0.19181 -0.087433 0.0167084
#> [3,] -0.0061743 -0.08743  0.331594 0.0688998
#> [4,]  0.0004588  0.01671  0.068900 0.0948096
fmt_ests_lbfgs$vcov$vcov_vary # with lbfgs
#>            [,1]     [,2]      [,3]      [,4]
#> [1,]  0.4155512  0.02345 -0.006135 0.0004166
#> [2,]  0.0234543  0.19150 -0.087508 0.0167449
#> [3,] -0.0061348 -0.08751  0.331640 0.0688572
#> [4,]  0.0004166  0.01674  0.068857 0.0948669
vcov_vary # the true values
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.35  0.02 -0.05  0.01
#> [2,]  0.02  0.12 -0.06 -0.01
#> [3,] -0.05 -0.06  0.32  0.09
#> [4,]  0.01 -0.01  0.09  0.12

# the parameters for the error term variance
fmt_ests$vcov$vcov_marker
#>        [,1]   [,2]
#> [1,] 0.3575 0.1003
#> [2,] 0.1003 0.1601
fmt_ests_lbfgs$vcov$vcov_marker # with lbfgs
#>        [,1]   [,2]
#> [1,] 0.3575 0.1003
#> [2,] 0.1003 0.1601
vcov_marker
#>      [,1] [,2]
#> [1,] 0.36 0.10
#> [2,] 0.10 0.16
```

### Recurrent Event

We simulate from a model with a recurrent event in this section. This is
a more relevant example as there is not a closed form solution of the
marginal likelihood.

``` r
# the survival parameters
fixef_surv <- c(-.5, .25) # delta
fixef_vary_surv <- c(.5, .1, -.2, .11) # omega
vcov_surv <- matrix(.2^2, 1) # Xi

library(splines)
b_func <- function(x)
  bs(x, knots = 1, Boundary.knots = c(0, 2))

# sample a few survival curves and plot them 
library(SimSurvNMarker)
# time points where we evaluate the conditional survival functions
tis <- seq(0, 2, length.out = 50)
set.seed(1)
# compute the conditional survival functions disregarding the fixed effect
surv_vals <- replicate(100, {
  err <- rnorm(1, sd = sqrt(vcov_surv))
  eval_surv_base_fun(
    ti = tis, omega = fixef_vary_surv, b_func = b_func, 
    gl_dat = get_gl_rule(100), delta = fixef_surv[1] + err)
})

par(mar = c(5, 5, 1, 1))
matplot(tis, surv_vals, type = "l", col = gray(0, .1), bty = "l", lty = 1,
        ylim = c(0, 1), xaxs = "i", yaxs = "i", xlab = "Time", 
        ylab = "Conditional survival function")
```

![](man/figures/README-only_recurrent-1.png)<!-- -->

``` r
# compute the average to get an estimate of the marginal survival curve
plot(tis, rowMeans(surv_vals), type = "l",  bty = "l", ylim = c(0, 1), 
     xaxs = "i", yaxs = "i", xlab = "Time", 
     ylab = "Marginal survival function")
```

![](man/figures/README-only_recurrent-2.png)<!-- -->

``` r
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  gl_dat <- get_gl_rule(100L)
  dat <- lapply(1:n_ids, function(id){
    # simulate the survival outcome
    rng_surv <- rnorm(1, sd = sqrt(vcov_surv))
    Z <- c(1, runif(1, -1, 1))
    log_haz_offset <- sum(Z * fixef_surv) + rng_surv

    # the conditional survival function
    surv_func <- function(ti)
      eval_surv_base_fun(
        ti = ti, omega = fixef_vary_surv, b_func = b_func, gl_dat = gl_dat,
        delta = log_haz_offset)

    # simulate the recurrent events one at a time
    max_sample <- 10L
    Z <- matrix(rep(Z, each = max_sample), max_sample)
    event <- y <- lf_trunc <- rep(NA_real_, max_sample)
    lf_trunc_i <- 0
    rngs <- runif(max_sample)
    left_trunc_surv <- 1
    for(i in 1:max_sample){
      # sample a random uniform variable and invert the survival function
      rng_i <- rngs[i]
      lf_trunc[i] <- lf_trunc_i

      root_func <- function(x) rng_i - surv_func(x) / left_trunc_surv
      if(root_func(2) < 0){
        # the observation is right-censored and we can exit
        y[i] <- 2
        event[i] <- 0
        break
      }

      # we need to invert the survival function to find the observation time
      root <- uniroot(root_func, c(lf_trunc_i, 2), tol = 1e-6)
      lf_trunc_i <- y[i] <- root$root
      event[i] <- 1
      left_trunc_surv <- surv_func(y[i])
    }

    # keep the needed data and return
    colnames(Z) <- paste0("Z_", 1:NCOL(Z) - 1L)
    surv_data <- cbind(lf_trunc = lf_trunc[1:i], y = y[1:i], event = event[1:i],
                       Z[1:i, -1, drop = FALSE], id = id)

    list(surv_data = surv_data)
  })

  # combine the data and return
  surv_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "surv_data")))
  surv_data$id <- as.integer(surv_data$id)
  # the order does not matter
  surv_data <- surv_data[sample.int(NROW(surv_data)), ]

  list(surv_data = surv_data)
}

# sample a moderately large data set
set.seed(3)
dat <- sim_dat(1000L)

# we show a few properties of the data below
NROW(dat$surv_data) # number of survival outcomes
#> [1] 2388
sum(dat$surv_data$event) # number of observed events
#> [1] 1388

# example of the data for three individuals
subset(dat$surv_data, id %in% 1:3)
#>   lf_trunc      y event     Z_1 id
#> 5   1.6383 2.0000     0  0.1208  3
#> 3   0.0000 0.3526     1  0.1208  3
#> 4   0.3526 1.6383     1  0.1208  3
#> 1   0.0000 2.0000     0 -0.2301  1
#> 2   0.0000 2.0000     0  0.6594  2

# distribution of observed events per individual
proportions(table(table(dat$surv_data$id) - 1L))
#> 
#>     0     1     2     3     4     5     6 
#> 0.266 0.343 0.220 0.103 0.050 0.012 0.006

# estimate the model with this package. Get the object we need for the
# optimization
library(survival)
library(VAJointSurv)
surv_obj <- surv_term(
  Surv(lf_trunc, y, event) ~ Z_1, id = id, dat$surv_data,
  time_fixef = bs_term(y, knots = 1, Boundary.knots = c(0, 2)), 
  with_frailty = TRUE)

comp_obj <- joint_ms_ptr(markers = list(),
                         survival_terms = surv_obj, max_threads = 4L)
rm(surv_obj)

# get the starting values
system.time(start_val <- joint_ms_start_val(comp_obj))
#>    user  system elapsed 
#>   0.241   0.000   0.157

# lower bound at the starting values
print(-attr(start_val, "value"), digits = 8)
#> [1] -1864.1515

# check that the gradient is correct
f <- function(x){
  start_val <- comp_obj$start_val
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(comp_obj$start_val, 7 + 2 * 2)), 
          head(joint_ms_lb_gr(comp_obj, comp_obj$start_val), 7 + 2 * 2), 
          tolerance = 1e-6)
#> [1] TRUE

# find the maximum lower bound estimate
system.time(opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 1000L, 
                                    pre_method = 3L, cg_tol = .2, c2 = .1))
#>    user  system elapsed 
#>   0.776   0.000   0.195

# check the gradient norm. We may need to reduce the convergence tolerance if 
# this is not small. In can also be a sign of convergence issues
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.1184
opt_out$info # convergence code (0 == 'OK')
#> [1] 0
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -1863.2795
opt_out$counts
#> function gradient     n_cg 
#>      156      101      165

# check the estimates
fmt_est <- joint_ms_format(comp_obj, opt_out$par)

rbind(estimate = fmt_est$survival[[1]]$fixef, truth = fixef_surv)
#>                        
#> estimate -0.5058 0.2582
#> truth    -0.5000 0.2500
rbind(estimate = fmt_est$survival[[1]]$fixef_vary, truth = fixef_vary_surv)
#>                                       
#> estimate 0.3965 0.2539 -0.3217 0.07835
#> truth    0.5000 0.1000 -0.2000 0.11000
c(estimate = sqrt(fmt_est$vcov$vcov_surv), truth = sqrt(vcov_surv))
#> estimate    truth 
#>    0.228    0.200
```

#### Note on Quadrature Rule

A quadrature rule is used to integrate the expected cumulative hazard.
The user can pass any quadrature rule with
![n](https://render.githubusercontent.com/render/math?math=n "n") nodes
and weights each denoted by
![(n\_i, w\_i)](https://render.githubusercontent.com/render/math?math=%28n_i%2C%20w_i%29 "(n_i, w_i)")
such that

![\\int\_0^1f(x) dx \\approx \\sum\_{i = 1}^nw\_i f(n\_i).](https://render.githubusercontent.com/render/math?math=%5Cint_0%5E1f%28x%29%20dx%20%5Capprox%20%5Csum_%7Bi%20%3D%201%7D%5Enw_i%20f%28n_i%29. "\int_0^1f(x) dx \approx \sum_{i = 1}^nw_i f(n_i).")

By default Gauss–Legendre quadrature is used. A very simple alternative
is to use the midpoint rule. This is illustrated below.

``` r
# computes the midpoint rule
mid_rule <- function(n)
  list(node = seq(0, 1, length.out = n), weight = rep(1 / n, n))

# compare the midpoint with different number of nodes
for(n in 2^(3:12)){
  res <- joint_ms_lb(comp_obj, opt_out$par, n_threads = 4L, 
                     quad_rule = mid_rule(n))
  cat(sprintf("# nodes, lower bound  %5d %12.4f\n", n, res))
}
#> # nodes, lower bound      8    1859.8265
#> # nodes, lower bound     16    1861.5170
#> # nodes, lower bound     32    1862.3858
#> # nodes, lower bound     64    1862.8294
#> # nodes, lower bound    128    1863.0536
#> # nodes, lower bound    256    1863.1663
#> # nodes, lower bound    512    1863.2228
#> # nodes, lower bound   1024    1863.2511
#> # nodes, lower bound   2048    1863.2652
#> # nodes, lower bound   4096    1863.2723

# the maximum lower bound value we got with Gauss–Legendre quadrature
print(-opt_out$value, digits = 8)
#> [1] -1863.2795
```

We can compare this with using Gauss–Legendre quadrature with different
number of nodes.

``` r
# get the rules
library(gaussquad)
#> Loading required package: polynom
#> Loading required package: orthopolynom
rules <- legendre.quadrature.rules(2^8, normalized = TRUE)

# rescale the rule
rules <- lapply(rules, function(x)
  list(node = .5 * (x[, "x"] + 1), weight = .5 * x[, "w"]))

# use different number of nodes
for(n in 2^(2:8)){
  res <- joint_ms_lb(comp_obj, opt_out$par, n_threads = 4L, 
                     quad_rule = rules[[n]])
  cat(sprintf("# nodes, lower bound  %5d %12.8f\n", n, res))
}
#> # nodes, lower bound      4 1863.11441967
#> # nodes, lower bound      8 1863.26694563
#> # nodes, lower bound     16 1863.27861371
#> # nodes, lower bound     32 1863.27932305
#> # nodes, lower bound     64 1863.27936970
#> # nodes, lower bound    128 1863.27937266
#> # nodes, lower bound    256 1863.27937285

# the result we got
print(-opt_out$value, digits = 14)
#> [1] -1863.2795232469
length(comp_obj$quad_rule$node) # the number of nodes we used
#> [1] 25
```

Evaluating the approximate expected cumulative hazard is by far the most
computationally expensive part of the likelihood. Therefore, it may be
nice to use fewer nodes e.g. when finding starting values.

### Two Markers and a Recurrent Event

We simulate from a model with two different types of markers and a
recurrent event in this section.

A Weibull model is selected for the baseline hazard by taking a
polynomial of degree one and using the log of time. As an illustration,
we select a more flexible expansion in the baseline hazard in the model
we estimate.

``` r
# settings for the simulation
library(splines)
g_funcs <- list(
  function(x)
    ns(x, knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
  function(x)
    # a raw polynomial
    outer(x, 1:2, `^`))
m_funcs <- list(
  function(x)
    ns(x, knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
  function(x)
    # a raw polynomial
    outer(x, 0:1, `^`))

fixef_vary_marker <- list(c(1.4, -1.2, -2.1), c(.5, .67)) # beta
fixef_marker <- list(c(-.5, 1), .25) # gamma

# Psi
vcov_vary <- structure(
  c(0.35, 0.02, -0.05, 0.01, 0.02, 0.12, -0.06, -0.01, -0.05, -0.06, 0.32, 0.09, 0.01, -0.01, 0.09, 0.12),
  .Dim = c(4L, 4L))
vcov_marker <- matrix(c(.6^2, .1, .1, .4^2), 2) # Sigma

# plot the markers' mean curve
library(VAJointSurv)
par(mar = c(5, 5, 1, 1))
plot_marker(
  time_fixef = ns_term(
    knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
  time_rng = ns_term(
    knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
  fixef_vary = fixef_vary_marker[[1]], x_range = c(0, 2), 
  vcov_vary = vcov_vary[1:2, 1:2], ylab = "Marker 1")
```

![](man/figures/README-recurrent_and_marker-1.png)<!-- -->

``` r
plot_marker(
  time_fixef = poly_term(degree = 2, raw = TRUE),
  time_rng = poly_term(degree = 1, raw = TRUE, intercept = TRUE), 
  fixef_vary = fixef_vary_marker[[2]], x_range = c(0, 2), 
  vcov_vary = vcov_vary[3:4, 3:4], ylab = "Marker 2")
```

![](man/figures/README-recurrent_and_marker-2.png)<!-- -->

``` r
# the survival parameters
fixef_surv <- c(-.5, .25) # delta
associations <- c(-.8, .7) # alpha
fixef_vary_surv <- c(.5) # omega
vcov_surv <- matrix(.2^2, 1) # Xi

b_term <- poly_term(degree = 1L, use_log = TRUE, raw = TRUE)
b_func <- function(x)
  t(b_term$eval(x))

# plot the hazard with pointwise quantiles
plot_surv(
  time_fixef = b_term, 
  time_rng = list(
    ns_term(knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
    poly_term(degree = 1, raw = TRUE, intercept = TRUE)), 
  x_range = c(0, 2), vcov_vary = vcov_vary, frailty_var = vcov_surv,
  ps = c(.25, .5, .75), log_hazard_shift = fixef_surv[1], 
  fixef_vary = fixef_vary_surv, associations = associations)
```

![](man/figures/README-recurrent_and_marker-3.png)<!-- -->

``` r
# without the markers
plot_surv(
  time_fixef = b_term, 
  time_rng = list(
    ns_term(knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
    poly_term(degree = 1, raw = TRUE, intercept = TRUE)), 
  x_range = c(0, 2), vcov_vary = matrix(0, 4, 4), frailty_var = vcov_surv,
  ps = c(.25, .5, .75), log_hazard_shift = fixef_surv[1], 
  fixef_vary = fixef_vary_surv, associations = associations)
```

![](man/figures/README-recurrent_and_marker-4.png)<!-- -->

``` r
# compute a few survival curves and plot them 
library(SimSurvNMarker)
 # time points where we evaluate the conditional survival functions
tis <- seq(0, 2, length.out = 50)
set.seed(1)
Us <- mvtnorm::rmvnorm(250, sigma = vcov_vary) # the random effects
# compute the conditional survival functions disregarding the fixed effect
surv_vals <- apply(Us, 1, function(U){
  expansion <- function(x)
    cbind(b_func(x), m_funcs[[1]](x) %*% U[1:2], 
          m_funcs[[2]](x) %*% U[3:4])
  
  err <- rnorm(1, sd = sqrt(vcov_surv))
  
  eval_surv_base_fun(
    ti = tis, omega = c(fixef_vary_surv, associations), b_func = expansion, 
    gl_dat = get_gl_rule(100), delta = fixef_surv[1] + err)
})

matplot(tis, surv_vals, type = "l", col = gray(0, .1), bty = "l", lty = 1,
        ylim = c(0, 1), xaxs = "i", yaxs = "i", xlab = "Time", 
        ylab = "Conditional survival function")
```

![](man/figures/README-recurrent_and_marker-5.png)<!-- -->

``` r
# compute the average to get an estimate of the marginal survival curve
plot(tis, rowMeans(surv_vals), type = "l",  bty = "l", ylim = c(0, 1), 
     xaxs = "i", yaxs = "i", xlab = "Time", 
     ylab = "Marginal survival function")
```

![](man/figures/README-recurrent_and_marker-6.png)<!-- -->

``` r
library(mvtnorm)
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  gl_dat <- get_gl_rule(100L)
  dat <- lapply(1:n_ids, function(id){
    # sample the number of outcomes and the fixed effect covariates
    n_obs <- sample.int(8L, 1L)
    obs_time <- sort(runif(n_obs, 0, 2))
    X1 <- cbind(1, rnorm(n_obs))
    X2 <- matrix(1, n_obs)
    colnames(X1) <- paste0("X1_", 1:NCOL(X1) - 1L)
    colnames(X2) <- paste0("X2_", 1:NCOL(X2) - 1L)
    X <- list(X1, X2)

    # sample the outcomes
    U <- drop(rmvnorm(1, sigma = vcov_vary))
    eta <- sapply(1:2, function(i)
      X[[i]] %*% fixef_marker[[i]] +
      drop(g_funcs[[i]](obs_time)) %*% fixef_vary_marker[[i]] +
      drop(m_funcs[[i]](obs_time)) %*% U[1:2 + (i == 2) * 2])

    ys <- eta + rmvnorm(n_obs, sigma = vcov_marker)
    colnames(ys) <- paste0("Y", 1:2)

    # mask some observations
    do_mask <- sample.int(3L, n_obs, replace = TRUE)
    ys[do_mask == 2, 1] <- NA
    ys[do_mask == 3, 2] <- NA

    X <- do.call(cbind, lapply(X, function(x) x[, -1, drop = FALSE]))
    marker_data <- cbind(ys, X, time = obs_time, id = id)
    
    # clean up 
    rm(list = setdiff(ls(), c("U", "marker_data", "id")))

    # simulate the survival outcome
    rng_surv <- rnorm(1, sd = sqrt(vcov_surv))
    Z <- c(1, runif(1, -1, 1))
    log_haz_offset <- sum(Z * fixef_surv) + rng_surv

    expansion <- function(x)
      cbind(b_func(x), m_funcs[[1]](x) %*% U[1:2],
            m_funcs[[2]](x) %*% U[3:4])

    # the conditional survival function
    surv_func <- function(ti)
      eval_surv_base_fun(
        ti = ti, omega = c(fixef_vary_surv, associations), b_func = expansion,
        gl_dat = gl_dat, delta = log_haz_offset)

    # simulate the recurrent events one at a time
    max_sample <- 10L
    Z <- matrix(rep(Z, each = max_sample), max_sample)
    event <- y <- lf_trunc <- rep(NA_real_, max_sample)
    lf_trunc_i <- 0
    left_trunc_surv <- 1
    for(i in 1:max_sample){
      # sample a random uniform variable and invert the survival function
      rng_i <- runif(1)
      lf_trunc[i] <- lf_trunc_i

      root_func <- function(x) rng_i - surv_func(x) / left_trunc_surv
      if(root_func(2) < 0){
        # the observation is right-censored and we can exit
        y[i] <- 2
        event[i] <- 0
        break
      }

      # we need to invert the survival function to find the observation time
      root <- uniroot(root_func, c(lf_trunc_i, 2), tol = 1e-6)
      lf_trunc_i <- y[i] <- root$root
      event[i] <- 1
      left_trunc_surv <- surv_func(y[i])
    }

    # keep the needed data and return
    colnames(Z) <- paste0("Z_", 1:NCOL(Z) - 1L)

    surv_data <- cbind(lf_trunc = lf_trunc[1:i], y = y[1:i], event = event[1:i],
                       Z[1:i, -1, drop = FALSE], id = id)

    list(marker_data = marker_data, surv_data = surv_data)
  })

  # combine the data and return
  marker_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "marker_data")))
  marker_data$id <- as.integer(marker_data$id)
  # the order does not matter
  marker_data <- marker_data[sample.int(NROW(marker_data)), ]

  surv_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "surv_data")))
  surv_data$id <- as.integer(surv_data$id)
  # the order does not matter
  surv_data <- surv_data[sample.int(NROW(surv_data)), ]

  list(marker_data = marker_data, surv_data = surv_data)
}

# sample a moderately large data set
set.seed(2)
dat <- sim_dat(1000L)

# we show a few properties of the data below
NROW(dat$marker_data) # number of observed marker
#> [1] 4463
NROW(dat$surv_data) # number of survival outcomes
#> [1] 2495
sum(dat$surv_data$event) # number of observed events
#> [1] 1500

# the data for one individual
subset(dat$marker_data, id == 1)
#>       Y1     Y2     X1_1   time id
#> 5     NA 4.3373 -0.23970 1.8877  1
#> 4     NA 4.8558  0.70795 1.8869  1
#> 3     NA 3.1915  0.13242 1.4047  1
#> 2  1.399 1.7030 -0.08025 1.1467  1
#> 1 -1.286 0.4379 -1.13038 0.3361  1
subset(dat$surv_data, id == 1)
#>   lf_trunc     y event     Z_1 id
#> 2    1.091 1.205     1 -0.9706  1
#> 3    1.205 2.000     0 -0.9706  1
#> 1    0.000 1.091     1 -0.9706  1

# distribution of observed events per individual
proportions(table(table(dat$surv_data$id) - 1L))
#> 
#>     0     1     2     3     4     5     6     7     8     9 
#> 0.320 0.320 0.154 0.090 0.057 0.031 0.008 0.005 0.004 0.011

# estimate the model with this package. Get the object we need for the
# optimization
marker_1 <- marker_term(
    Y1 ~ X1_1, id = id, subset(dat$marker_data, !is.na(Y1)),
    time_fixef = ns_term(time, knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
    time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 2),
                       intercept = TRUE))
marker_2 <- marker_term(
    Y2 ~ 1, id = id, subset(dat$marker_data, !is.na(Y2)),
    time_fixef = poly_term(time, degree = 2, raw = TRUE),
    time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

surv_obj <- surv_term(
  Surv(lf_trunc, y, event) ~ Z_1, id = id, dat$surv_data,
  # we select a more flexible model for the baseline hazard
  time_fixef = poly_term(y, degree = 3L, raw = TRUE, use_log = TRUE),
  with_frailty = TRUE)

comp_obj <- joint_ms_ptr(markers = list(marker_1, marker_2),
                         survival_terms = surv_obj, max_threads = 4L)
rm(marker_1, marker_2, surv_obj)

# get the starting values
system.time(start_val <- joint_ms_start_val(comp_obj, gr_tol = .1))
#>    user  system elapsed 
#>   3.847   0.004   1.233

# lower bound at the starting values
print(-attr(start_val, "value"), digits = 8)
#> [1] -7542.5634

# check that the gradient is correct
f <- function(x){
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(start_val, 29 + 2 * 20)), 
          head(joint_ms_lb_gr(comp_obj, start_val), 29 + 2 * 20))
#> [1] "Mean relative difference: 6.73e-08"

# find the maximum lower bound estimate
system.time(opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 1000L, 
                                    pre_method = 3L, cg_tol = .2, c2 = .1, 
                                    gr_tol = .1))
#>    user  system elapsed 
#>  14.258   0.000   3.566

# we set gr_tol in the call so this is the convergence criterion for the 
# gradient
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.09921
opt_out$info # convergence code (0 == 'OK')
#> [1] 0
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -7257.0219
opt_out$counts
#> function gradient     n_cg 
#>      490      354      970

# find the maximum lower bound with lbfgs
library(lbfgsb3c)
system.time(lbfgs_res <- lbfgsb3c(
  start_val, function(x) joint_ms_lb(comp_obj, x),
  function(x) joint_ms_lb_gr(comp_obj, x), 
  control = list(factr = 1e-8, maxit = 2000L)))
#>    user  system elapsed 
#>  63.695   0.019  16.038
lbfgs_res$convergence # convergence code (0 == 'OK')
#> [1] 1
print(-lbfgs_res$value, digits = 8)  # maximum lower bound value
#> [1] -7257.4153
lbfgs_res$counts # may have hit maxit!
#> [1] 2000 2000

# compare the estimates with the actual values. Start with the fixed effects
fmt_ests <- joint_ms_format(comp_obj, opt_out$par)

# the parameters for the first marker
fmt_ests$markers[[1]] 
#> $fixef
#>                 
#> -0.4923  1.0171 
#> 
#> $fixef_vary
#>                      
#>  1.403 -1.191 -2.031

fixef_marker[[1]] # true values
#> [1] -0.5  1.0
fixef_vary_marker[[1]] # true values
#> [1]  1.4 -1.2 -2.1

# the parameters for the second marker
fmt_ests$markers[[2]]
#> $fixef
#>       
#> 0.244 
#> 
#> $fixef_vary
#>               
#> 0.6133 0.6280
fixef_marker[[2]] # true values
#> [1] 0.25
fixef_vary_marker[[2]] # true values
#> [1] 0.50 0.67

# the fixed effects for the survival outcome and the association parameters
fmt_ests$survival[[1]]
#> $fixef
#>                 
#> -0.5070  0.2613 
#> 
#> $fixef_vary
#>                         
#> 0.40382 0.07573 0.03291 
#> 
#> $associations
#>                 
#> -0.7299  0.7049
fixef_surv
#> [1] -0.50  0.25
fixef_vary_surv
#> [1] 0.5
associations
#> [1] -0.8  0.7

# the parameters for covariance matrix of the random effects
fmt_ests$vcov$vcov_vary
#>           [,1]      [,2]     [,3]      [,4]
#> [1,]  0.316939  0.008933 -0.04181  0.012407
#> [2,]  0.008933  0.167065 -0.05168 -0.004297
#> [3,] -0.041814 -0.051679  0.30518  0.090884
#> [4,]  0.012407 -0.004297  0.09088  0.139207
vcov_vary # the true values
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.35  0.02 -0.05  0.01
#> [2,]  0.02  0.12 -0.06 -0.01
#> [3,] -0.05 -0.06  0.32  0.09
#> [4,]  0.01 -0.01  0.09  0.12

# the parameters for the error term covariance matrix
fmt_ests$vcov$vcov_marker
#>         [,1]    [,2]
#> [1,] 0.37498 0.09843
#> [2,] 0.09843 0.16335
vcov_marker
#>      [,1] [,2]
#> [1,] 0.36 0.10
#> [2,] 0.10 0.16

# the parameters for the frailty covariance matrix
fmt_ests$vcov$vcov_surv
#>         [,1]
#> [1,] 0.03787
vcov_surv
#>      [,1]
#> [1,] 0.04
```

#### Observed Information Matrix and Approximate Wald Intervals

The package supplies an approximation of the observed information matrix
and the full Hessian of all the parameters. This is illustrated below.

``` r
# compute the Hessian
system.time(hess <- joint_ms_hess(comp_obj, par = opt_out$par))
#>    user  system elapsed 
#>   7.801   0.058   7.857
dim(hess$hessian_all) # the full matrix!
#> [1] 20029 20029

# compare parts it with those from numerical differentiation from R
n_comp <- 150L
hess_num <- numDeriv::jacobian(
  function(x){
    par <- opt_out$par
    par[1:n_comp] <- x
    joint_ms_lb_gr(comp_obj, par)[1:n_comp]
  }, head(opt_out$par, n_comp))

# did we get the same?
all.equal(hess_num, as.matrix(hess$hessian_all[1:n_comp, 1:n_comp]), 
          check.attributes = FALSE)
#> [1] TRUE

# compute the covariance matrix from the approximate observed information matrix
# the part only for the model parameters accounting for the variational 
# parameters
dim(hess$hessian)
#> [1] 29 29
obs_mat <- -hess$hessian # the observed information matrix
SEs <- sqrt(diag(solve(-obs_mat))) # approximate standard errors

# Wald type confidence intervals for the association parameters
idx_assoc <- comp_obj$indices$survival[[1]]$associations
matrix(opt_out$par[idx_assoc] + outer(SEs[idx_assoc], c(-1, 1) * 1.96), 2,
       dimnames = list(c("alpha_1", "alpha_2"), 
                       sprintf("%.2f pct", c(2.5, 97.5))))
#>         2.50 pct 97.50 pct
#> alpha_1  -1.0613   -0.3985
#> alpha_2   0.6314    0.7785

# we can do this for all parameters. We illustrate this by showing the
# estimates along with standard errors
do_comb <- function(est, se){
  if(is.list(est))
    setNames(mapply(do_comb, est, se, SIMPLIFY = FALSE), names(est))
  else 
    rbind(Estimates = est, SE = se)
}

par_fmt <- joint_ms_format(comp_obj, opt_out$par)
par_fmt_SE <- joint_ms_format(comp_obj, SEs)
est_n_se <- 
  do_comb(par_fmt[c("markers", "survival")], par_fmt_SE[c("markers", "survival")])

# show the result
est_n_se$markers[[1]]
#> $fixef
#>                           
#> Estimates -0.49227 1.01709
#> SE         0.04692 0.01205
#> 
#> $fixef_vary
#>                                  
#> Estimates 1.40268 -1.1909 -2.0313
#> SE        0.05532  0.1179  0.0451
fixef_marker[[1]]
#> [1] -0.5  1.0
fixef_vary_marker[[1]]
#> [1]  1.4 -1.2 -2.1

est_n_se$markers[[2]]
#> $fixef
#>                  
#> Estimates 0.24398
#> SE        0.03147
#> 
#> $fixef_vary
#>                         
#> Estimates 0.6133 0.62797
#> SE        0.0616 0.02978
fixef_marker[[2]]
#> [1] 0.25
fixef_vary_marker[[2]]
#> [1] 0.50 0.67

est_n_se$survival[[1]]
#> $fixef
#>                           
#> Estimates -0.50701 0.26125
#> SE         0.04753 0.05066
#> 
#> $fixef_vary
#>                                  
#> Estimates 0.40382 0.07573 0.03291
#> SE        0.05953 0.08393 0.02596
#> 
#> $associations
#>                          
#> Estimates -0.7299 0.70493
#> SE         0.1691 0.03752
fixef_surv
#> [1] -0.50  0.25
fixef_vary_surv
#> [1] 0.5
associations
#> [1] -0.8  0.7
```

We can compute standard error estimates for the covariance matrices’
parameters but this requires an application of the delta method because
of the parameterization of the covariance matrices.

#### Approximate Likelihood Ratio based Confidence Intervals

Assuming that the lower bound is tight, we can construct approximate
likelihood ratio based confidence intervals using the lower bound. We
show how to do this below with the `mask` argument of `joint_ms_opt`.

``` r
# fixed input 
level <- .95 # confidence level
which_fix <- 14L
# see the indices element for which element is fixed. It is an association 
# parameter in this case
comp_obj$indices$survival[[1]]$associations # confidence interval for 
#> [1] 14 15

# assumed plausible values of the parameter. The joint_ms_profile function
# shown later finds these automatically
params <- opt_out$par[which_fix]
params <- seq(params - .5, params + .5, length.out = 15)

# find the maximum lower bound values
lbs_max <- sapply(params, function(x){
  par <- opt_out$par
  par[which_fix] <- x
  res <- joint_ms_opt(comp_obj, par = par, max_it = 1000L, 
                      pre_method = 3L, cg_tol = .2, c2 = .1, 
                      # -1 needed in the psqn package (zero-based indices)
                      mask = which_fix - 1L)
  
  # return the maximum lower bound if the method converged
  if(!res$convergence)
    stop("did not converge")
  -res$value
})
```

``` r
# find the critical value and the approximate confidence interval
z_vals <- sqrt(pmax(-opt_out$value - lbs_max, 0) * 2)
z_vals <- ifelse(params < opt_out$par[which_fix], -1, 1) * z_vals

alpha <- 1 - level
pvs <- c(alpha / 2, 1 - alpha/2)
conf_int <- setNames(approx(z_vals, params, xout = qnorm(pvs))$y,
                     sprintf("%.2f pct.", 100 * pvs))
conf_int # the approximate confidence interval
#>  2.50 pct. 97.50 pct. 
#>    -1.0786    -0.4066

# plot the approximate log profile likelihood and highlight the critical value
par(mar = c(5, 5, 1, 1))
plot(params, lbs_max, pch = 16, bty = "l", xlab = "Association parameter", 
     ylab = "Approximate log profile likelihood")
grid()
smooth_est <- smooth.spline(params, lbs_max)
lines(predict(smooth_est, seq(min(params), max(params), length.out = 100)))

# mark the confidence interval
abline(h = -opt_out$value - qchisq(level, 1) / 2, lty = 2)
abline(v = conf_int, lty = 2)
```

![](man/figures/README-res_manual_pl-1.png)<!-- -->

The `joint_ms_profile` uses similar steps to the above to find
approximate profile likelihood based confidence intervals. An example is
shown below.

``` r
# construct the approximate likelihood ratio based confidence interval
system.time(
  prof_conf <- joint_ms_profile(
    comp_obj, opt_out = opt_out, which_prof = which_fix, delta = .25, 
    level = level, 
    max_it = 1000L, pre_method = 3L, cg_tol = .2, c2 = .1))
#> 
#> Finding the upper limit of the approximate profile likelihood curve
#> LogLike: -7258.1642 at        -0.479927
#> LogLike: -7261.6330 at        -0.229927
#> LogLike: -7259.3810 at        -0.371898. Lb, target, ub: -7259.3810, -7258.9426, -7258.1642
#> LogLike: -7258.8961 at        -0.410612. Lb, target, ub: -7259.3810, -7258.9426, -7258.8961
#> 
#> Finding the lower limit of the approximate profile likelihood curve
#> LogLike: -7258.0398 at        -0.979927
#> LogLike: -7260.7686 at        -1.229927
#> LogLike: -7259.3109 at        -1.108780. Lb, target, ub: -7259.3109, -7258.9426, -7258.0398
#> LogLike: -7258.8517 at        -1.067638. Lb, target, ub: -7259.3109, -7258.9426, -7258.8517
#> LogLike: -7257.0219 at        -0.729927
#>    user  system elapsed 
#>  92.235   0.036  23.155
```

``` r
prof_conf$confs # the approximate confidence interval
#>  2.50 pct. 97.50 pct. 
#>    -1.0761    -0.4067

# plot the approximate log profile likelihood and highlight the critical value
par(mar = c(5, 5, 1, 1))

with(prof_conf, {
  plot(xs, p_log_Lik, pch = 16, bty = "l", 
     xlab = "Association parameter", ylab = "Approximate log profile likelihood")
  grid()
  smooth_est <- smooth.spline(xs, p_log_Lik)
  lines(predict(smooth_est, seq(min(xs), max(xs), length.out = 100)))
  
  abline(v = confs, lty = 2)
})
```

![](man/figures/README-res_joint_ms_profile-1.png)<!-- -->

The Hessian can be used to possibly decrease the computation time by
passing it to the `hess` argument as illustrated below. Starting values
along the profile likelihood curve is then found using an approximate
normal distribution.

``` r
# construct the approximate likelihood ratio based confidence interval
system.time(
  prof_conf_fast <- joint_ms_profile(
    comp_obj, opt_out = opt_out, which_prof = which_fix, 
    # we can use the standard error
    delta = 2.5 * sqrt(diag(solve(hess$hessian)))[which_fix], 
    level = level, 
    max_it = 1000L, pre_method = 3L, cg_tol = .2, c2 = .1,
    hess = hess))
#> Computing the vector to find the starting values along the profile likelihood curve
#> 
#> Finding the upper limit of the approximate profile likelihood curve
#> LogLike: -7260.3168 at        -0.307236
#> LogLike: -7257.0219 at        -0.729927
#> LogLike: -7257.8975 at        -0.510901. Lb, target, ub: -7260.3168, -7258.9426, -7257.8975
#> LogLike: -7258.9267 at        -0.407876. Lb, target, ub: -7260.3168, -7258.9426, -7258.9267
#> LogLike: -7259.0789 at        -0.395498. Lb, target, ub: -7259.0789, -7258.9426, -7258.9267
#> 
#> Finding the lower limit of the approximate profile likelihood curve
#> LogLike: -7259.8185 at        -1.152618
#> LogLike: -7257.0219 at        -0.729927
#> LogLike: -7258.1151 at        -0.987975. Lb, target, ub: -7259.8185, -7258.9426, -7258.1151
#> LogLike: -7258.9017 at        -1.072598. Lb, target, ub: -7259.8185, -7258.9426, -7258.9017
#> LogLike: -7259.0338 at        -1.084992. Lb, target, ub: -7259.0338, -7258.9426, -7258.9017
#> LogLike: -7257.0219 at        -0.729927
#>    user  system elapsed 
#>  67.824   1.241  27.067
```

``` r
# we got the same
prof_conf$confs 
#>  2.50 pct. 97.50 pct. 
#>    -1.0761    -0.4067
prof_conf_fast$confs
#>  2.50 pct. 97.50 pct. 
#>    -1.0765    -0.4066
```

### Two Markers and a Recurrent Event Without Frailty

In this section, we estimate a model like before but without the
frailty.

``` r
# set the frailty variance to ~0
vcov_surv <- matrix(1e-6^2, 1)

# plot the hazard with pointwise quantiles
par(mar = c(5, 5, 1, 1))
plot_surv(
  time_fixef = b_term, 
  time_rng = list(
    ns_term(knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
    poly_term(degree = 1, raw = TRUE, intercept = TRUE)), 
  x_range = c(0, 2), vcov_vary = vcov_vary, frailty_var = vcov_surv,
  ps = c(.25, .5, .75), log_hazard_shift = fixef_surv[1], 
  fixef_vary = fixef_vary_surv, associations = associations)
```

![](man/figures/README-no_frailty_recurrent_and_marker-1.png)<!-- -->

``` r
# sample a moderately large data set
set.seed(2)
dat <- sim_dat(1000L)

# we show a few properties of the data below
NROW(dat$marker_data) # number of observed marker
#> [1] 4480
NROW(dat$surv_data) # number of survival outcomes
#> [1] 2566
sum(dat$surv_data$event) # number of observed events
#> [1] 1571

# the data for one individual
subset(dat$marker_data, id == 1)
#>       Y1     Y2     X1_1   time id
#> 5     NA 4.3373 -0.23970 1.8877  1
#> 4     NA 4.8558  0.70795 1.8869  1
#> 3     NA 3.1915  0.13242 1.4047  1
#> 2  1.399 1.7030 -0.08025 1.1467  1
#> 1 -1.286 0.4379 -1.13038 0.3361  1
subset(dat$surv_data, id == 1)
#>   lf_trunc     y event     Z_1 id
#> 2    1.206 1.330     1 -0.9706  1
#> 3    1.330 2.000     0 -0.9706  1
#> 1    0.000 1.206     1 -0.9706  1

# distribution of observed events per individual
proportions(table(table(dat$surv_data$id) - 1L))
#> 
#>     0     1     2     3     4     5     6     7     8     9 
#> 0.317 0.279 0.184 0.102 0.057 0.029 0.010 0.007 0.004 0.011

# estimate the model with this package. Get the object we need for the
# optimization
marker_1 <- marker_term(
    Y1 ~ X1_1, id = id, subset(dat$marker_data, !is.na(Y1)),
    time_fixef = ns_term(time, knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
    time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 2),
                       intercept = TRUE))
marker_2 <- marker_term(
    Y2 ~ 1, id = id, subset(dat$marker_data, !is.na(Y2)),
    time_fixef = poly_term(time, degree = 2, raw = TRUE),
    time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

surv_obj <- surv_term(
  Surv(lf_trunc, y, event) ~ Z_1, id = id, dat$surv_data,
  # we select a more flexible model for the baseline hazard
  time_fixef = poly_term(y, degree = 3L, raw = TRUE, use_log = TRUE))

comp_obj <- joint_ms_ptr(markers = list(marker_1, marker_2),
                         survival_terms = surv_obj, max_threads = 4L)
rm(marker_1, marker_2, surv_obj)

# get the starting values
system.time(start_val <- joint_ms_start_val(comp_obj))
#>    user  system elapsed 
#>   2.565   0.000   0.742

# lower bound at the starting values
print(-attr(start_val, "value"), digits = 8)
#> [1] -7423.5688

# check that the gradient is correct
f <- function(x){
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(start_val, 28 + 2 * 14)), 
          head(joint_ms_lb_gr(comp_obj, start_val), 28 + 2 * 14))
#> [1] "Mean relative difference: 6.211e-08"

# find the maximum lower bound estimate
system.time(opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 1000L, 
                                    pre_method = 3L, cg_tol = .2, c2 = .1, 
                                    rel_eps = 1e-10))
#>    user  system elapsed 
#>   9.849   0.000   2.465

# check the gradient norm. We may need to reduce the convergence tolerance if 
# this is not small. In can also be a sign of convergence issues
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.02053
opt_out$info # convergence code (0 == 'OK')
#> [1] 0
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -7093.4009
opt_out$counts
#> function gradient     n_cg 
#>      404      302      675

# compare the estimates with the actual values. Start with the fixed effects
fmt_ests <- joint_ms_format(comp_obj, opt_out$par)

# the parameters for the first marker
fmt_ests$markers[[1]] 
#> $fixef
#> [1] -0.5361  1.0058
#> 
#> $fixef_vary
#> [1]  1.439 -1.134 -2.067

fixef_marker[[1]] # true values
#> [1] -0.5  1.0
fixef_vary_marker[[1]] # true values
#> [1]  1.4 -1.2 -2.1

# the parameters for the second marker
fmt_ests$markers[[2]]
#> $fixef
#> [1] 0.2043
#> 
#> $fixef_vary
#> [1] 0.6436 0.6173
fixef_marker[[2]] # true values
#> [1] 0.25
fixef_vary_marker[[2]] # true values
#> [1] 0.50 0.67

# the fixed effects for the survival outcome and the association parameters
fmt_ests$survival[[1]]
#> $fixef
#> [1] -0.4367  0.2499
#> 
#> $fixef_vary
#> [1]  0.426448 -0.007011  0.010087
#> 
#> $associations
#> [1] -0.8853  0.7342
fixef_surv
#> [1] -0.50  0.25
fixef_vary_surv
#> [1] 0.5
associations
#> [1] -0.8  0.7

# the parameters for covariance matrix of the random effects
fmt_ests$vcov$vcov_vary
#>          [,1]     [,2]     [,3]     [,4]
#> [1,]  0.30520  0.01911 -0.06209  0.02126
#> [2,]  0.01911  0.08416 -0.03360 -0.02355
#> [3,] -0.06209 -0.03360  0.30456  0.08503
#> [4,]  0.02126 -0.02355  0.08503  0.14097
vcov_vary # the true values
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.35  0.02 -0.05  0.01
#> [2,]  0.02  0.12 -0.06 -0.01
#> [3,] -0.05 -0.06  0.32  0.09
#> [4,]  0.01 -0.01  0.09  0.12

# the parameters for the error term covariance matrix
fmt_ests$vcov$vcov_marker
#>         [,1]    [,2]
#> [1,] 0.36068 0.09824
#> [2,] 0.09824 0.15530
vcov_marker
#>      [,1] [,2]
#> [1,] 0.36 0.10
#> [2,] 0.10 0.16

# the parameters for the frailty covariance matrix
fmt_ests$vcov$vcov_surv # NULL as there is not any
#> NULL
```

#### Using the Estimated Variational Parameters

The `joint_ms_va_par` function can be used to extract the estimated
variational parameters for each individual. These can be used to get
marginal estimates as the variational approximation for each individual
is an approximation of the conditional distribution of the random
effects given the observed data.

``` r
# get the variational parameters
va_par <- joint_ms_va_par(comp_obj, opt_out$par)

# it is a list with the estimated variational parameters for each individual
va_par[c("157", "345", "878")]
#> $`157`
#> $`157`$mean
#> [1] -0.47333 -0.02992  0.85784  0.37737
#> 
#> $`157`$vcov
#>           [,1]       [,2]      [,3]       [,4]
#> [1,]  0.196325 -0.0094636 -0.011033  0.0255751
#> [2,] -0.009464  0.0674062  0.003957  0.0006269
#> [3,] -0.011033  0.0039570  0.098740 -0.0459766
#> [4,]  0.025575  0.0006269 -0.045977  0.0412375
#> 
#> 
#> $`345`
#> $`345`$mean
#> [1] -0.03226 -0.11443  0.92850  0.42210
#> 
#> $`345`$vcov
#>           [,1]      [,2]      [,3]      [,4]
#> [1,]  0.184690  0.015114 -0.005131  0.032541
#> [2,]  0.015114  0.072407 -0.001721  0.001451
#> [3,] -0.005131 -0.001721  0.038315 -0.029652
#> [4,]  0.032541  0.001451 -0.029652  0.064874
#> 
#> 
#> $`878`
#> $`878`$mean
#> [1] 0.23917 0.03412 0.25500 0.14915
#> 
#> $`878`$vcov
#>           [,1]       [,2]      [,3]       [,4]
#> [1,]  0.189760  2.004e-02 -0.004531  1.121e-02
#> [2,]  0.020041  7.358e-02 -0.004654 -4.778e-05
#> [3,] -0.004531 -4.654e-03  0.055993 -2.594e-02
#> [4,]  0.011205 -4.778e-05 -0.025936  4.042e-02

# they can be used to approximate marginal measures such as the mean curve for 
# the individual along with 95% pointwise confidence intervals as shown below 
# for where the true curve is. You can skip the code below and go the plots after
plot_mean_curve <- function(who){
  # get the data for the individual
  m_dat <- subset(dat$marker_data, id == who)
  
  # plot the estimated mean population curve with pointwise confidence intervals
  par(mfcol = c(1, 2), bty = "l")
  fmt_ests <- joint_ms_format(comp_obj, opt_out$par)
  plot_marker(
    time_fixef = ns_term(
      knots = c(.5, 1.5), Boundary.knots = c(0, 2)),
    time_rng = ns_term(
      knots = numeric(), Boundary.knots = c(0, 2), intercept = TRUE), 
    fixef_vary = fmt_ests$markers[[1]]$fixef_vary, x_range = c(0, 2), 
    vcov_vary = fmt_ests$vcov$vcov_vary[1:2, 1:2], ylab = "Marker 1")
  
  # add the points for the individual
  y_hat <- m_dat$Y1 - cbind(1, m_dat$X1_1) %*% fmt_ests$markers[[1]]$fixef
  points(m_dat$time, y_hat, pch = 16)
  
  # plot the estimated mean curve with pointwise confidence intervals
  va_par_who <- va_par[[who]]
  xs <- seq(0, 2, length.out = 100)
  mea <- g_funcs[[1]](xs) %*% fmt_ests$markers[[1]]$fixef_vary + 
    m_funcs[[1]](xs) %*% va_par_who$mean[1:2]
  lines(xs, mea, lty = 2)
  se <- apply(m_funcs[[1]](xs), 1, function(x)
    sqrt(x %*% va_par_who$vcov[1:2, 1:2] %*% x))
  polygon(c(xs, rev(xs)), c(mea + 1.96 * se, rev(mea - 1.95 * se)), border = NA,
          col = gray(0, .1))
  
  # plot the estimated mean population curve with pointwise confidence intervals
  plot_marker(
    time_fixef = poly_term(degree = 2, raw = TRUE),
    time_rng = poly_term(degree = 1, raw = TRUE, intercept = TRUE), 
    fixef_vary = fmt_ests$markers[[2]]$fixef_vary, x_range = c(0, 2), 
    vcov_vary = fmt_ests$vcov$vcov_vary[3:4, 3:4], ylab = "Marker 2")
  
  # add the points for the individual
  y_hat <- m_dat$Y2 - fmt_ests$markers[[2]]$fixef
  points(m_dat$time, y_hat, pch = 16)
  
  # plot the estimated mean curve with confidence intervals
  xs <- seq(0, 2, length.out = 100)
  mea <- g_funcs[[2]](xs) %*% fmt_ests$markers[[2]]$fixef_vary + 
    m_funcs[[2]](xs) %*% va_par_who$mean[2:3]
  lines(xs, mea, lty = 2)
  se <- apply(m_funcs[[2]](xs), 1, function(x)
    sqrt(x %*% va_par_who$vcov[3:4, 3:4] %*% x))
  polygon(c(xs, rev(xs)), c(mea + 1.96 * se, rev(mea - 1.95 * se)), border = NA,
          col = gray(0, .1))
}

# plot the curves for three individuals
plot_mean_curve(157)
```

![](man/figures/README-use_VA_par-1.png)<!-- -->

``` r
plot_mean_curve(235)
```

![](man/figures/README-use_VA_par-2.png)<!-- -->

``` r
plot_mean_curve(878)
```

![](man/figures/README-use_VA_par-3.png)<!-- -->

### Two Markers, the Observation Time Process, and a Terminal Event

We simulate and fit a model in this section where we have two markers, a
recurrent event which is the observation times, and a terminal event.

``` r
# settings for the simulation
library(splines)
g_funcs <- list(
  function(x)
    ns(x, knots = c(3.33, 6.67), Boundary.knots = c(0, 10)),
  function(x)
    # a raw polynomial
    outer(x, 1:2, `^`))
m_funcs <- list(
  function(x)
    ns(x, knots = numeric(), Boundary.knots = c(0, 10), intercept = TRUE),
  function(x)
    # a raw polynomial
    outer(x, 0:1, `^`))

fixef_vary_marker <- list(c(1.4, 1.2, -2.1), c(.5, -.02)) # beta
fixef_marker <- list(c(-.5, 2), 1) # gamma

# Psi
vcov_vary <- structure(c(0.35, 0.08, -0.05, 0.01, 0.08, 1.92, -0.24, -0.04,
                   -0.05, -0.24, 0.32, 0.09, 0.01, -0.04, 0.09, 0.12),
                 .Dim = c(4L, 4L))
vcov_marker <- matrix(c(.6^2, .1, .1, .4^2), 2)

# plot the markers' mean curve
library(VAJointSurv)
par(mar = c(5, 5, 1, 1))
plot_marker(
  time_fixef = ns_term(
    knots = c(3.33, 6.67), Boundary.knots = c(0, 10), intercept = FALSE),
  time_rng = ns_term(
    knots = numeric(), Boundary.knots = c(0, 10), intercept = TRUE),
  fixef_vary = fixef_vary_marker[[1]], x_range = c(0, 10),
  vcov_vary = vcov_vary[1:2, 1:2], ylab = "Marker 1")
```

![](man/figures/README-obs_process_markers_and_recurrent-1.png)<!-- -->

``` r
plot_marker(
  time_fixef = poly_term(degree = 2, raw = TRUE),
  time_rng = poly_term(degree = 1, raw = TRUE, intercept = TRUE),
  fixef_vary = fixef_vary_marker[[2]], x_range = c(0, 10),
  vcov_vary = vcov_vary[3:4, 3:4], ylab = "Marker 2")
```

![](man/figures/README-obs_process_markers_and_recurrent-2.png)<!-- -->

``` r
# the survival parameters
vcov_surv <- matrix(c(.2^2, .15^2, .15^2, .25^2), 2) # Xi 

fixef_surv <- list(c(-1, .25), .2)
associations <- list(c(.6, -.4), c(-.7, .2))
fixef_vary_surv <- list(c(.5, .1, -.2, .11),
                          c(-1, -.25))

b_funcs <- list(
  function(x) bs(x, knots = 5, Boundary.knots = c(0, 10)),
  function(x) ns(x, knots = 5, Boundary.knots = c(0, 10)))

# plot the log hazard with the 25%, 50% and 75% quantiles
library(SimSurvNMarker)

plot_surv(
  time_fixef = bs_term(knots = 5, Boundary.knots = c(0, 10)),
  time_rng = list(
    ns_term(knots = numeric(), Boundary.knots = c(0, 10), intercept = TRUE),
    poly_term(degree = 1, raw = TRUE, intercept = TRUE)),
  x_range = c(0, 10), fixef_vary = fixef_vary_surv[[1]],
  vcov_vary = vcov_vary, frailty_var = vcov_surv[1, 1], ps = c(.25, .5, .75),
  associations = associations[[1]], log_hazard_shift = fixef_surv[[1]][1],
  ylab = "Terminal event")
```

![](man/figures/README-obs_process_markers_and_recurrent-3.png)<!-- -->

``` r
plot_surv(
  time_fixef = ns_term(knots = 5, Boundary.knots = c(0, 10)),
  time_rng = list(
    ns_term(knots = numeric(), Boundary.knots = c(0, 10), intercept = TRUE),
    poly_term(degree = 1, raw = TRUE, intercept = TRUE)),
  x_range = c(0, 10), fixef_vary = fixef_vary_surv[[2]],
  vcov_vary = vcov_vary, frailty_var = vcov_surv[2, 2], ps = c(.25, .5, .75),
  associations = associations[[2]], log_hazard_shift = fixef_surv[[2]][1],
  ylab = "Observation process")
```

![](man/figures/README-obs_process_markers_and_recurrent-4.png)<!-- -->

``` r
library(mvtnorm)
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  gl_dat <- get_gl_rule(100L)
  dat <- lapply(1:n_ids, function(id){
    # sample the terminal event time and the censoring time
    cens <- min(rexp(1, rate = 1/10), 10)
    U <- drop(rmvnorm(1, sigma = vcov_vary))

    frailties <- drop(rmvnorm(1, sigma = vcov_surv))
    Z1 <- c(1, runif(1, -1, 1))
    log_haz_offset <- sum(Z1 * fixef_surv[[1]]) + frailties[1]

    # assign the conditional hazard function
    expansion <- function(x, b_func)
      cbind(b_func(x), m_funcs[[1]](x) %*% U[1:2],
            m_funcs[[2]](x) %*% U[3:4])
    surv_func <- function(ti, fixef_vary_surv, associations, b_func){
      formals(expansion)$b_func <- b_func
      eval_surv_base_fun(
        ti = ti, omega = c(fixef_vary_surv, associations), b_func = expansion,
        gl_dat = gl_dat, delta = log_haz_offset)
    }

    # sample the survival time
    rng <- runif(1)
    root_func <- function(x, rng)
      rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[1]],
                      associations = associations[[1]], b_func = b_funcs[[1]])

    if(root_func(cens, rng) < 0){
      # the observation is censored
      y_terminal <- cens
      event <- 0
    } else {
      # find the event time
      root <- uniroot(root_func, c(0, cens), tol = 1e-6, rng = rng)
      y_terminal <- root$root
      event <- 1

    }

    terminal_outcome <- cbind(y = y_terminal, event = event, Z1 = Z1[2],
                              id = id)

    # clean up
    rm(list = setdiff(ls(), c("y_terminal", "terminal_outcome", "expansion",
                              "surv_func", "frailties", "U", "id")))

    # simulate the observation times
    Z2 <- 1
    log_haz_offset <- sum(Z2 * fixef_surv[[2]]) + frailties[2]

    root_func <- function(x, left_trunc_surv, rng)
      rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[2]],
                      associations = associations[[2]], b_func = b_funcs[[2]]) /
      left_trunc_surv

    max_sample <- 1000L
    left_trunc_surv <- 1
    Z2 <- matrix(rep(Z2, each = max_sample), max_sample)
    event <- y <- lf_trunc <- rep(NA_real_, max_sample)
    lf_trunc_i <- 0
    for(i in 1:max_sample){
      # sample a random uniform variable and invert the survival function
      rng_i <- runif(1)
      lf_trunc[i] <- lf_trunc_i

      if(root_func(y_terminal, left_trunc_surv, rng_i) < 0){
        # the observation is right-censored and we can exit
        y[i] <- y_terminal
        event[i] <- 0
        break
      }

      # we need to invert the survival function to find the observation time
      root <- uniroot(root_func, c(lf_trunc_i, y_terminal), tol = 1e-6,
                      left_trunc_surv = left_trunc_surv, rng = rng_i)
      lf_trunc_i <- y[i] <- root$root
      event[i] <- 1
      left_trunc_surv <- surv_func(
        y[i], fixef_vary_surv = fixef_vary_surv[[2]], associations = associations[[2]],
        b_func = b_funcs[[2]])
    }

    colnames(Z2) <- paste0("Z", 1:NCOL(Z2) - 1L)
    obs_process <- cbind(lf_trunc = lf_trunc[1:i], y = y[1:i],
                         event = event[1:i], Z2[1:i, -1, drop = FALSE],
                         id = id)

    # clean up
    rm(list = setdiff(ls(), c("terminal_outcome", "U", "id",
                              "obs_process")))

    # sample the number of outcomes and the fixed effect covariates
    obs_time <- c(0, obs_process[obs_process[, "event"] == 1, "y"])
    n_obs <- length(obs_time)
    X1 <- cbind(1, rnorm(n_obs))
    X2 <- matrix(1, n_obs)
    colnames(X1) <- paste0("X1_", 1:NCOL(X1) - 1L)
    colnames(X2) <- paste0("X2_", 1:NCOL(X2) - 1L)
    X <- list(X1, X2)

    # sample the outcomes
    eta <- sapply(1:2, function(i)
      X[[i]] %*% fixef_marker[[i]] +
        drop(g_funcs[[i]](obs_time)) %*% fixef_vary_marker[[i]] +
        drop(m_funcs[[i]](obs_time)) %*% U[1:2 + (i == 2) * 2])

    ys <- eta + rmvnorm(n_obs, sigma = vcov_marker)
    colnames(ys) <- paste0("Y", 1:2)

    # mask some observations
    do_mask <- sample.int(3L, n_obs, replace = TRUE)
    ys[do_mask == 2, 1] <- NA
    ys[do_mask == 3, 2] <- NA

    X <- do.call(cbind, lapply(X, function(x) x[, -1, drop = FALSE]))
    marker_data <- cbind(ys, X, time = obs_time, id = id)

    return(list(marker_data = marker_data, obs_process = obs_process,
                terminal_outcome = terminal_outcome))
  })

  # combine the data and return
  marker_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "marker_data")))
  marker_data$id <- as.integer(marker_data$id)
  # the order does not matter
  marker_data <- marker_data[sample.int(NROW(marker_data)), ]

  obs_process <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "obs_process")))
  obs_process$id <- as.integer(obs_process$id)
  # the order does not matter
  obs_process <- obs_process[sample.int(NROW(obs_process)), ]

  terminal_outcome <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "terminal_outcome")))
  terminal_outcome$id <- as.integer(terminal_outcome$id)
  # the order does not matter
  terminal_outcome <- terminal_outcome[sample.int(NROW(terminal_outcome)), ]

  list(marker_data = marker_data, obs_process = obs_process,
       terminal_outcome = terminal_outcome)
}

# sample a moderately large data set
set.seed(2)
dat <- sim_dat(1000L)

# we show a few properties of the data below
mean(dat$terminal_outcome$event) # mean event rate
#> [1] 0.782
sum(dat$obs_process$event) # number of observed markers less the individuals
#> [1] 2421
NROW(dat$marker_data) # number of observed markers less the individuals
#> [1] 3421

# distribution of observed marker per individual
proportions(table(table(dat$obs_process$id)))
#> 
#>     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23 
#> 0.353 0.194 0.135 0.076 0.069 0.045 0.030 0.024 0.023 0.010 0.005 0.006 0.008 0.003 0.001 0.004 0.001 0.003 0.001 0.001 0.001 0.003 0.001 
#>    24    25    26 
#> 0.001 0.001 0.001

# show data for one individual
subset(dat$marker_data, id == 1)
#>          Y1    Y2     X1_1   time id
#> X.3      NA 4.041 -0.02815 2.4032  1
#> X.2      NA 2.164 -1.03429 0.2988  1
#> X.1 -0.3812 2.646 -0.28571 0.2603  1
#> X   -2.6308 1.584 -1.43968 0.0000  1
subset(dat$obs_process, id == 1)
#>   lf_trunc      y event id
#> 2   0.2603 0.2988     1  1
#> 1   0.0000 0.2603     1  1
#> 4   2.4032 2.8900     0  1
#> 3   0.2988 2.4032     1  1
subset(dat$terminal_outcome, id == 1)
#>      y event      Z1 id
#> 1 2.89     1 -0.6384  1

# estimate the model with this package. Get the object we need for the
# optimization
marker_1 <- marker_term(
  Y1 ~ X1_1, id = id, subset(dat$marker_data, !is.na(Y1)),
  time_fixef = ns_term(time, knots = c(3.33, 6.67), Boundary.knots = c(0, 10)),
  time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 10),
                     intercept = TRUE))
marker_2 <- marker_term(
  Y2 ~ 1, id = id, subset(dat$marker_data, !is.na(Y2)),
  time_fixef = poly_term(time, degree = 2, raw = TRUE),
  time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

library(survival)
surv_terminal <- surv_term(
  Surv(y, event) ~ Z1, id = id, dat$terminal_outcome,
  time_fixef = bs_term(y, knots = 5, Boundary.knots = c(0, 10)),
  with_frailty = TRUE)
surv_obs <- surv_term(
  Surv(lf_trunc, y, event) ~ 1, id = id, dat$obs_process,
  time_fixef = ns_term(y, knots = 5, Boundary.knots = c(0, 10)),
  with_frailty = TRUE)

comp_obj <- joint_ms_ptr(markers = list(marker_1, marker_2),
                         survival_terms = list(surv_terminal, surv_obs),
                         max_threads = 4L)
rm(marker_1, marker_2, surv_terminal, surv_obs)

# get the starting values
system.time(start_val <- joint_ms_start_val(comp_obj, gr_tol = .1))
#>    user  system elapsed 
#>  26.787   0.080   7.168

# lower bound at the starting values
print(-attr(start_val, "value"), digits = 8)
#> [1] -7885.9054

# check that the gradient is correct
f <- function(x){
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(start_val, 37 + 2 * 27)),
          head(joint_ms_lb_gr(comp_obj, start_val), 37 + 2 * 27), 
          tolerance = 1e-6)
#> [1] TRUE

# find the maximum lower bound estimate
system.time(opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 2000L,
                                    pre_method = 3L, cg_tol = .2, c2 = .1,
                                    gr_tol = .1))
#>    user  system elapsed 
#> 329.727   0.148  82.567

# we set gr_tol in the call so this is the convergence criterion for the 
# gradient
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.09267
opt_out$info # convergence code (0 == 'OK')
#> [1] 0
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -7600.1698
opt_out$counts
#> function gradient     n_cg 
#>     4530     3231    32080

# compare the estimates with the actual values. Start with the fixed effects
fmt_ests <- joint_ms_format(comp_obj, opt_out$par)

# the parameters for the first marker
fmt_ests$markers[[1]]
#> $fixef
#>                 
#> -0.5136  1.9734 
#> 
#> $fixef_vary
#>                      
#>  1.400  1.560 -1.606

fixef_marker[[1]] # true values
#> [1] -0.5  2.0
fixef_vary_marker[[1]] # true values
#> [1]  1.4  1.2 -2.1

# the parameters for the second marker
fmt_ests$markers[[2]]
#> $fixef
#>       
#> 1.036 
#> 
#> $fixef_vary
#>                   
#>  0.48416 -0.02007
fixef_marker[[2]] # true values
#> [1] 1
fixef_vary_marker[[2]] # true values
#> [1]  0.50 -0.02

# the fixed effects for the survival outcome and the association parameters
# for the terminal event
fmt_ests$survival[[1]]
#> $fixef
#>                 
#> -0.8571  0.2568 
#> 
#> $fixef_vary
#>                                     
#> -0.07581  1.17619 -1.33199  0.31682 
#> 
#> $associations
#>                 
#>  0.5673 -0.3621
fixef_surv[[1]]
#> [1] -1.00  0.25
fixef_vary_surv[[1]]
#> [1]  0.50  0.10 -0.20  0.11
associations[[1]]
#> [1]  0.6 -0.4

# same for the observation process
fmt_ests$survival[[2]]
#> $fixef
#>        
#> 0.2476 
#> 
#> $fixef_vary
#>                 
#> -1.3554 -0.3508 
#> 
#> $associations
#>                 
#> -0.7371  0.2080
fixef_surv[[2]]
#> [1] 0.2
fixef_vary_surv[[2]]
#> [1] -1.00 -0.25
associations[[2]]
#> [1] -0.7  0.2

# the parameters for covariance matrix of the random effects
fmt_ests$vcov$vcov_vary
#>          [,1]      [,2]     [,3]     [,4]
#> [1,]  0.29424 -0.145972  0.05297 0.001100
#> [2,] -0.14597  1.480889 -0.12433 0.007764
#> [3,]  0.05297 -0.124332  0.27871 0.091859
#> [4,]  0.00110  0.007764  0.09186 0.116753
vcov_vary # the true values
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.35  0.08 -0.05  0.01
#> [2,]  0.08  1.92 -0.24 -0.04
#> [3,] -0.05 -0.24  0.32  0.09
#> [4,]  0.01 -0.04  0.09  0.12

# the parameters for the error term covariance matrix
fmt_ests$vcov$vcov_marker
#>        [,1]   [,2]
#> [1,] 0.3890 0.1104
#> [2,] 0.1104 0.1688
vcov_marker
#>      [,1] [,2]
#> [1,] 0.36 0.10
#> [2,] 0.10 0.16

# the parameters for the frailty covariance matrix
fmt_ests$vcov$vcov_surv
#>           [,1]     [,2]
#> [1,]  0.007361 -0.01350
#> [2,] -0.013496  0.02766
vcov_surv
#>        [,1]   [,2]
#> [1,] 0.0400 0.0225
#> [2,] 0.0225 0.0625
```

#### Caching Expansions

Some basis expansions like `bs_term` and `ns_term` take relatively long
time to evaluate in the approximation of the approximate expected
cumulative hazard. Thus, it may be advantageous to save the expansions
if the same quadrature rule is used. This is done by setting the
`cache_expansions` argument to true. The pros of doing this is that the
expensive basis expansions are only evaluated once which may decrees the
computation time. The cons are

-   We no longer use the CPU cache efficiently with present hardware.
    Thus, you may not see great advantages of using many threads and the
    performance  
    may even be worse. This is more likely to be an issue with larger
    data set.
-   It requires a lot more memory which may be an issue larger problems.

We illustrate this by showing the computation time where we change the
number of threads we use.

``` r
# with caching
w_caching <- bench::mark(
  `w/ caching 1 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 1, cache_expansions = TRUE),
  `w/ caching 2 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 2, cache_expansions = TRUE),
  `w/ caching 3 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 3, cache_expansions = TRUE),
  `w/ caching 4 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 4, cache_expansions = TRUE))
w_caching[, c("expression", "median")]
#> # A tibble: 4 × 2
#>   expression            median
#>   <bch:expr>          <bch:tm>
#> 1 w/ caching 1 thread   7.92ms
#> 2 w/ caching 2 thread   4.04ms
#> 3 w/ caching 3 thread   2.85ms
#> 4 w/ caching 4 thread   2.27ms

# difference between one and four threads
with(w_caching, median[4] / median[1]) 
#> [1] 287ms

# w/o caching
wo_caching <- bench::mark(
  `w/o caching 1 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 1, cache_expansions = FALSE),
  `w/o caching 2 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 2, cache_expansions = FALSE),
  `w/o caching 3 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 3, cache_expansions = FALSE),
  `w/o caching 4 thread` = 
    joint_ms_lb(comp_obj, opt_out$par, n_threads = 4, cache_expansions = FALSE))
wo_caching[, c("expression", "median")]
#> # A tibble: 4 × 2
#>   expression             median
#>   <bch:expr>           <bch:tm>
#> 1 w/o caching 1 thread  20.65ms
#> 2 w/o caching 2 thread  10.62ms
#> 3 w/o caching 3 thread   8.88ms
#> 4 w/o caching 4 thread   6.56ms

# difference between one and four threads
with(wo_caching, median[4] / median[1]) 
#> [1] 318ms
```

### Two Markers, the Observation Time Process, and a Terminal Event with Delayed Entry

We alter the previous example in this section by adding delayed entries.
In this case, the likelihood has to be altered to account for the
delayed entry as described by Crowther et al. (2016) and Berg and
Drepper (2016).

To show the likelihood we use, assume that there is only one type of
survival outcome,
![H = 1](https://render.githubusercontent.com/render/math?math=H%20%3D%201 "H = 1").
The log marginal likelihood of individual
![i](https://render.githubusercontent.com/render/math?math=i "i") where
we do not properly account for the delayed entry is

![\\log E\\left( f\_i(\\vec y\_i\\mid \\vec U\_i)h\_{i1}(t\_{i1}\\mid \\vec U\_i)^{d\_{i1}} \\frac{S\_{i1}(t\_{i1}\\mid \\vec U\_i)}{S\_{i1}(v\_{i1}\\mid \\vec U\_i)}\\right)](https://render.githubusercontent.com/render/math?math=%5Clog%20E%5Cleft%28%20f_i%28%5Cvec%20y_i%5Cmid%20%5Cvec%20U_i%29h_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5E%7Bd_%7Bi1%7D%7D%20%5Cfrac%7BS_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%7D%7BS_%7Bi1%7D%28v_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%7D%5Cright%29 "\log E\left( f_i(\vec y_i\mid \vec U_i)h_{i1}(t_{i1}\mid \vec U_i)^{d_{i1}} \frac{S_{i1}(t_{i1}\mid \vec U_i)}{S_{i1}(v_{i1}\mid \vec U_i)}\right)")

where
![\\vec y\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%20y_i "\vec y_i")
is observed markers of individual
![i](https://render.githubusercontent.com/render/math?math=i "i"),
![f\_i](https://render.githubusercontent.com/render/math?math=f_i "f_i")
is the conditional density for the markers of individual
![i](https://render.githubusercontent.com/render/math?math=i "i") given
the random effects
![\\vec U\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%20U_i "\vec U_i"),
![t\_{i1}](https://render.githubusercontent.com/render/math?math=t_%7Bi1%7D "t_{i1}")
is the observed time,
![d\_{i1}](https://render.githubusercontent.com/render/math?math=d_%7Bi1%7D "d_{i1}")
is an event indicator,
![h\_{i1}](https://render.githubusercontent.com/render/math?math=h_%7Bi1%7D "h_{i1}")
is the conditional hazard,
![S\_{i1}](https://render.githubusercontent.com/render/math?math=S_%7Bi1%7D "S_{i1}")
is the conditional survival function, and
![v\_{i1}\\in \[0,\\infty)](https://render.githubusercontent.com/render/math?math=v_%7Bi1%7D%5Cin%20%5B0%2C%5Cinfty%29 "v_{i1}\in [0,\infty)")
is the left-truncation time. We approximate the log marginal likelihood
term with a lower bound given by

![\\log E\\left(  f\_i(\\vec y\_i\\mid \\vec U\_i)h\_{i1}(t\_{i1}\\mid \\vec U\_i)^{d\_{i1}}  \\frac{S\_{i1}(t\_{i1}\\mid \\vec U\_i)}{S\_{i1}(v\_{i1}\\mid \\vec U\_i)}\\right) \\geq  E\_{Q\_{i1}}\\left(\\log\\left(  f\_i(\\vec y\_i\\mid \\vec U\_i)h\_{i1}(t\_{i1}\\mid \\vec U\_i)^{d\_{i1}}  \\frac{S\_{i1}(t\_{i1}\\mid \\vec U\_i)v(\\vec U\_i)}  {S\_{i1}(v\_{i1}\\mid \\vec U\_i)q\_{\\vec\\theta\_{i1}}(\\vec U\_i)}  \\right)\\right)](https://render.githubusercontent.com/render/math?math=%5Clog%20E%5Cleft%28%20%20f_i%28%5Cvec%20y_i%5Cmid%20%5Cvec%20U_i%29h_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5E%7Bd_%7Bi1%7D%7D%20%20%5Cfrac%7BS_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%7D%7BS_%7Bi1%7D%28v_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%7D%5Cright%29%20%5Cgeq%20%20E_%7BQ_%7Bi1%7D%7D%5Cleft%28%5Clog%5Cleft%28%20%20f_i%28%5Cvec%20y_i%5Cmid%20%5Cvec%20U_i%29h_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5E%7Bd_%7Bi1%7D%7D%20%20%5Cfrac%7BS_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29v%28%5Cvec%20U_i%29%7D%20%20%7BS_%7Bi1%7D%28v_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29q_%7B%5Cvec%5Ctheta_%7Bi1%7D%7D%28%5Cvec%20U_i%29%7D%20%20%5Cright%29%5Cright%29 "\log E\left(  f_i(\vec y_i\mid \vec U_i)h_{i1}(t_{i1}\mid \vec U_i)^{d_{i1}}  \frac{S_{i1}(t_{i1}\mid \vec U_i)}{S_{i1}(v_{i1}\mid \vec U_i)}\right) \geq  E_{Q_{i1}}\left(\log\left(  f_i(\vec y_i\mid \vec U_i)h_{i1}(t_{i1}\mid \vec U_i)^{d_{i1}}  \frac{S_{i1}(t_{i1}\mid \vec U_i)v(\vec U_i)}  {S_{i1}(v_{i1}\mid \vec U_i)q_{\vec\theta_{i1}}(\vec U_i)}  \right)\right)")

where the expectation on the right-hand side is using the density
![q\_{\\vec\\theta\_{i1}}](https://render.githubusercontent.com/render/math?math=q_%7B%5Cvec%5Ctheta_%7Bi1%7D%7D "q_{\vec\theta_{i1}}")
and ![v](https://render.githubusercontent.com/render/math?math=v "v") is
the unconditional density of the random effects. Since this is a lower
bound, we can jointly optimize the lower bound over the
![\\vec\\theta\_{i1}](https://render.githubusercontent.com/render/math?math=%5Cvec%5Ctheta_%7Bi1%7D "\vec\theta_{i1}")s
and the model parameters.

The correct likelihood with delayed entry is

![\\log E\\left(  f\_i(\\vec y\_i\\mid \\vec U\_i)h\_{i1}(t\_{i1}\\mid \\vec U\_i)^{d\_{i1}}  S\_{i1}(t\_{i1}\\mid \\vec U\_i)\\right)  - \\log E\\left(S\_{i1}(v\_{i1}\\mid \\vec U\_i)\\right).](https://render.githubusercontent.com/render/math?math=%5Clog%20E%5Cleft%28%20%20f_i%28%5Cvec%20y_i%5Cmid%20%5Cvec%20U_i%29h_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5E%7Bd_%7Bi1%7D%7D%20%20S_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5Cright%29%20%20-%20%5Clog%20E%5Cleft%28S_%7Bi1%7D%28v_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5Cright%29. "\log E\left(  f_i(\vec y_i\mid \vec U_i)h_{i1}(t_{i1}\mid \vec U_i)^{d_{i1}}  S_{i1}(t_{i1}\mid \vec U_i)\right)  - \log E\left(S_{i1}(v_{i1}\mid \vec U_i)\right).")

Of course, we can replace this with an approximation in the form of

![E\_{Q\_{i1}}\\left(\\log\\left(  f\_i(\\vec y\_i\\mid \\vec U\_i)h\_{i1}(t\_{i1}\\mid \\vec U\_i)^{d\_{i1}}  S\_{i1}(t\_{i1}\\mid \\vec U\_i)\\frac{v(\\vec U\_i)}{q\_{\\theta\_{i1}}(\\vec U\_i)}  \\right)\\right)  - E\_{Q\_{i2}}\\left(\\log\\left(S\_{i1}(v\_{i1}\\mid \\vec U\_i)  \\frac{v(\\vec U\_i)}{q\_{\\theta\_{i2}}(\\vec U\_i)}\\right)\\right).](https://render.githubusercontent.com/render/math?math=E_%7BQ_%7Bi1%7D%7D%5Cleft%28%5Clog%5Cleft%28%20%20f_i%28%5Cvec%20y_i%5Cmid%20%5Cvec%20U_i%29h_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5E%7Bd_%7Bi1%7D%7D%20%20S_%7Bi1%7D%28t_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5Cfrac%7Bv%28%5Cvec%20U_i%29%7D%7Bq_%7B%5Ctheta_%7Bi1%7D%7D%28%5Cvec%20U_i%29%7D%20%20%5Cright%29%5Cright%29%20%20-%20E_%7BQ_%7Bi2%7D%7D%5Cleft%28%5Clog%5Cleft%28S_%7Bi1%7D%28v_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%20%20%5Cfrac%7Bv%28%5Cvec%20U_i%29%7D%7Bq_%7B%5Ctheta_%7Bi2%7D%7D%28%5Cvec%20U_i%29%7D%5Cright%29%5Cright%29. "E_{Q_{i1}}\left(\log\left(  f_i(\vec y_i\mid \vec U_i)h_{i1}(t_{i1}\mid \vec U_i)^{d_{i1}}  S_{i1}(t_{i1}\mid \vec U_i)\frac{v(\vec U_i)}{q_{\theta_{i1}}(\vec U_i)}  \right)\right)  - E_{Q_{i2}}\left(\log\left(S_{i1}(v_{i1}\mid \vec U_i)  \frac{v(\vec U_i)}{q_{\theta_{i2}}(\vec U_i)}\right)\right).")

This has some implications:

1.  We still need to maximize over the
    ![\\vec\\theta\_{i1}](https://render.githubusercontent.com/render/math?math=%5Cvec%5Ctheta_%7Bi1%7D "\vec\theta_{i1}")s
    and the model parameters. However, we have to minimize over the
    ![\\vec\\theta\_{i2}](https://render.githubusercontent.com/render/math?math=%5Cvec%5Ctheta_%7Bi2%7D "\vec\theta_{i2}")s
    as they provide an upper bound on the log marginal likelihood. Thus,
    we have a maxmin problem rather than a pure maximization problem.
2.  The approximation is no longer guaranteed to be a lower bound on the
    log marginal likelihood because we are adding a lower bound on one
    term with an upper bound on the other term.
3.  The approximation of
    ![\\log E\\left(S\_{i1}(v\_{i1}\\mid \\vec U\_i)\\right)](https://render.githubusercontent.com/render/math?math=%5Clog%20E%5Cleft%28S_%7Bi1%7D%28v_%7Bi1%7D%5Cmid%20%5Cvec%20U_i%29%5Cright%29 "\log E\left(S_{i1}(v_{i1}\mid \vec U_i)\right)")
    may be poor as we do not have the multivariate normal density from
    the markers.

Currently, we compute the terms in the marginal likelihood due to
delayed entry with adaptive Gauss-Hermite quadrature. This can be quite
slow if there are many random effects that needs to marginalized out.
However, it seems that quite few quadrature nodes can be used to
approximate the terms from delayed entry when the marginal survival
probabilities are close to one (and the variance of the random effects
is not extremely large). However, more quadrature nodes seems to be
needed for the gradient of marginal likelihood terms due to delayed
entry. Thus, you may experience odd behavior (very small improvements in
the likelihood approximation at each step) of the optimizer if more
quadrature nodes are needed.

We do believe that the maxmin problem mentioned in 1. can be solved
efficiently but our current “brute force” attempts where not fast.

Delayed entry is handled by the `delayed` argument of `surv_term` as
illustrated in the example below.

``` r
library(mvtnorm)
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  gl_dat <- get_gl_rule(100L)
  dat <- lapply(1:n_ids, function(id){
    # sample the delayed entry time
    delayed_entry <- pmax(runif(1, -1, 5), 0)
    
    # sample the censoring time
    cens <- -Inf
    while(cens < delayed_entry)
      cens <- min(rexp(1, rate = 1/10), 10)
    
    # sample the terminal event time and the random effects
    Z1 <- c(1, runif(1, -1, 1))
    
    y_terminal <- -Inf 
    while(y_terminal < delayed_entry){
      U <- drop(rmvnorm(1, sigma = vcov_vary))
      frailties <- drop(rmvnorm(1, sigma = vcov_surv))
      log_haz_offset <- sum(Z1 * fixef_surv[[1]]) + frailties[1]
  
      # assign the conditional hazard function
      expansion <- function(x, b_func)
        cbind(b_func(x), m_funcs[[1]](x) %*% U[1:2],
              m_funcs[[2]](x) %*% U[3:4])
      surv_func <- function(ti, fixef_vary_surv, associations, b_func){
        formals(expansion)$b_func <- b_func
        eval_surv_base_fun(
          ti = ti, omega = c(fixef_vary_surv, associations), b_func = expansion,
          gl_dat = gl_dat, delta = log_haz_offset)
      }
  
      # sample the survival time
      rng <- runif(1)
      root_func <- function(x, rng)
        rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[1]],
                        associations = associations[[1]], b_func = b_funcs[[1]])
  
      if(root_func(cens, rng) < 0){
        # the observation is censored
        y_terminal <- cens
        event <- 0
      } else {
        # find the event time
        root <- uniroot(root_func, c(0, cens), tol = 1e-6, rng = rng)
        y_terminal <- root$root
        event <- 1
  
      }
    }

    terminal_outcome <- cbind(y = y_terminal, event = event, Z1 = Z1[2],
                              id = id, delayed_entry = delayed_entry)

    # clean up
    rm(list = setdiff(ls(), c(
      "y_terminal", "terminal_outcome", "expansion", "surv_func", "frailties", 
      "U", "id", "delayed_entry")))

    # simulate the observation times
    Z2 <- 1
    log_haz_offset <- sum(Z2 * fixef_surv[[2]]) + frailties[2]

    root_func <- function(x, left_trunc_surv, rng)
      rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[2]],
                      associations = associations[[2]], b_func = b_funcs[[2]]) /
      left_trunc_surv

    max_sample <- 1000L
    left_trunc_surv <- 1
    Z2 <- matrix(rep(Z2, each = max_sample), max_sample)
    event <- y <- lf_trunc <- rep(NA_real_, max_sample)
    lf_trunc_i <- 0
    for(i in 1:max_sample){
      # sample a random uniform variable and invert the survival function
      rng_i <- runif(1)
      lf_trunc[i] <- lf_trunc_i

      if(root_func(y_terminal, left_trunc_surv, rng_i) < 0){
        # the observation is right-censored and we can exit
        y[i] <- y_terminal
        event[i] <- 0
        break
      }

      # we need to invert the survival function to find the observation time
      root <- uniroot(root_func, c(lf_trunc_i, y_terminal), tol = 1e-6,
                      left_trunc_surv = left_trunc_surv, rng = rng_i)
      lf_trunc_i <- y[i] <- root$root
      event[i] <- 1
      left_trunc_surv <- surv_func(
        y[i], fixef_vary_surv = fixef_vary_surv[[2]], associations = associations[[2]],
        b_func = b_funcs[[2]])
    }

    colnames(Z2) <- paste0("Z", 1:NCOL(Z2) - 1L)
    obs_process <- cbind(lf_trunc = lf_trunc[1:i], y = y[1:i],
                         event = event[1:i], Z2[1:i, -1, drop = FALSE],
                         id = id)
    
    # account for the delayed entry 
    obs_process[, "lf_trunc"] <- pmax(delayed_entry, obs_process[, "lf_trunc"])
    obs_process <- obs_process[
      obs_process[, "y"] > delayed_entry, , drop = FALSE]

    # clean up
    rm(list = setdiff(ls(), c("terminal_outcome", "U", "id",
                              "obs_process", "delayed_entry")))

    # sample the number of outcomes and the fixed effect covariates
    obs_time <- c(delayed_entry, obs_process[obs_process[, "event"] == 1, "y"])
    
    n_obs <- length(obs_time)
    X1 <- cbind(1, rnorm(n_obs))
    X2 <- matrix(1, n_obs)
    colnames(X1) <- paste0("X1_", 1:NCOL(X1) - 1L)
    colnames(X2) <- paste0("X2_", 1:NCOL(X2) - 1L)
    X <- list(X1, X2)

    # sample the outcomes
    eta <- sapply(1:2, function(i)
      X[[i]] %*% fixef_marker[[i]] +
        drop(g_funcs[[i]](obs_time)) %*% fixef_vary_marker[[i]] +
        drop(m_funcs[[i]](obs_time)) %*% U[1:2 + (i == 2) * 2])

    ys <- eta + rmvnorm(n_obs, sigma = vcov_marker)
    colnames(ys) <- paste0("Y", 1:2)

    # mask some observations
    do_mask <- sample.int(3L, n_obs, replace = TRUE)
    ys[do_mask == 2, 1] <- NA
    ys[do_mask == 3, 2] <- NA

    X <- do.call(cbind, lapply(X, function(x) x[, -1, drop = FALSE]))
    marker_data <- cbind(ys, X, time = obs_time, id = id)

    return(list(marker_data = marker_data, obs_process = obs_process,
                terminal_outcome = terminal_outcome))
  })

  # combine the data and return
  marker_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "marker_data")))
  marker_data$id <- as.integer(marker_data$id)
  # the order does not matter
  marker_data <- marker_data[sample.int(NROW(marker_data)), ]

  obs_process <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "obs_process")))
  obs_process$id <- as.integer(obs_process$id)
  # the order does not matter
  obs_process <- obs_process[sample.int(NROW(obs_process)), ]

  terminal_outcome <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "terminal_outcome")))
  terminal_outcome$id <- as.integer(terminal_outcome$id)
  # the order does not matter
  terminal_outcome <- terminal_outcome[sample.int(NROW(terminal_outcome)), ]

  list(marker_data = marker_data, obs_process = obs_process,
       terminal_outcome = terminal_outcome)
}

# sample a moderate sized data set
set.seed(2)
dat <- sim_dat(1000L)

# we show a few properties of the data below
mean(dat$terminal_outcome$event) # mean event rate
#> [1] 0.662
sum(dat$obs_process$event) # number of observed markers less the individuals
#> [1] 2729
NROW(dat$marker_data) # number of observed markers less the individuals
#> [1] 3729
# fraction of individuals with a delayed entry time
mean(dat$terminal_outcome$delayed_entry > 0) 
#> [1] 0.845
 
# distribution of observed marker per individual
proportions(table(table(dat$obs_process$id)))
#> 
#>     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23 
#> 0.343 0.186 0.140 0.085 0.047 0.040 0.046 0.021 0.016 0.016 0.010 0.011 0.006 0.004 0.005 0.004 0.003 0.001 0.001 0.003 0.001 0.003 0.001 
#>    24    25    27    30    33    43 
#> 0.001 0.001 0.001 0.002 0.001 0.001

# show data for one individual
subset(dat$marker_data, id == 1)
#>          Y1    Y2    X1_1   time id
#> X        NA 2.029  0.4177 0.1093  1
#> X.2 -2.0826 2.277 -0.3927 0.8605  1
#> X.1 -0.9811    NA  0.9818 0.5344  1
subset(dat$obs_process, id == 1)
#>   lf_trunc      y event id
#> 2   0.5344 0.8605     1  1
#> 1   0.1093 0.5344     1  1
#> 3   0.8605 1.1949     0  1
subset(dat$terminal_outcome, id == 1)
#>       y event     Z1 id delayed_entry
#> 1 1.195     1 0.1467  1        0.1093

# estimate the model with this package. Get the object we need for the
# optimization while NOT account and accounting for the delayed entry
marker_1 <- marker_term(
  Y1 ~ X1_1, id = id, subset(dat$marker_data, !is.na(Y1)),
  time_fixef = ns_term(time, knots = c(3.33, 6.67), Boundary.knots = c(0, 10)),
  time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 10),
                     intercept = TRUE))
marker_2 <- marker_term(
  Y2 ~ 1, id = id, subset(dat$marker_data, !is.na(Y2)),
  time_fixef = poly_term(time, degree = 2, raw = TRUE),
  time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

# the wrong way 
library(survival)
surv_terminal_wrong <- surv_term(
  Surv(delayed_entry, y, event) ~ Z1, id = id, dat$terminal_outcome,
  time_fixef = bs_term(y, knots = 5, Boundary.knots = c(0, 10)),
  with_frailty = TRUE)

# the right way 
surv_terminal <- surv_term(
  Surv(delayed_entry, y, event) ~ Z1, id = id, dat$terminal_outcome,
  time_fixef = bs_term(y, knots = 5, Boundary.knots = c(0, 10)),
  with_frailty = TRUE, 
  # some have delayed entry
  delayed = delayed_entry > 0)

surv_obs <- surv_term(
  Surv(lf_trunc, y, event) ~ 1, id = id, dat$obs_process,
  time_fixef = ns_term(y, knots = 5, Boundary.knots = c(0, 10)),
  with_frailty = TRUE)

# the wrong way 
comp_obj_wrong <- joint_ms_ptr(
  markers = list(marker_1, marker_2),
  survival_terms = list(surv_terminal_wrong, surv_obs),
  max_threads = 4L)

# the right way
comp_obj <- joint_ms_ptr(
  markers = list(marker_1, marker_2),
  survival_terms = list(surv_terminal, surv_obs),
  max_threads = 4L)

# get the starting values
system.time(start_val_wrong <- joint_ms_start_val(comp_obj_wrong, gr_tol = .1))
#>    user  system elapsed 
#>  14.103   0.008   3.836
system.time(start_val <- joint_ms_start_val(comp_obj, gr_tol = .1))
#>    user  system elapsed 
#>  22.766   0.031   5.978

# lower bound at the starting values
print(-attr(start_val_wrong, "value"), digits = 8)
#> [1] -8725.0019
print(-attr(start_val, "value"), digits = 8)
#> [1] -8725.0413

# check that the gradient is correct
f <- function(x){
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(start_val, 37 + 2 * 27)),
          head(joint_ms_lb_gr(comp_obj, start_val), 37 + 2 * 27),
          tolerance = 1e-6)
#> [1] TRUE

# find the maximum lower bound estimate
system.time(opt_out_wrong <- joint_ms_opt(
  comp_obj_wrong, par = start_val_wrong, max_it = 2000L, pre_method = 3L, 
  cg_tol = .2, c2 = .1, gr_tol = .1))
#>    user  system elapsed 
#>  76.672   0.012  19.173

system.time(opt_out <- joint_ms_opt(
  comp_obj, par = opt_out_wrong$par, max_it = 2000L, pre_method = 3L, 
  # the function is more expensive. Thus, it makes more sense to get a better
  # approximation in the conjugate gradient step
  cg_tol = .1, c2 = .1, 
  # it takes much longer to get the same precision
  gr_tol = 1))
#>     user   system  elapsed 
#> 1837.573    0.486  460.513

# we set gr_tol in the call so this is the convergence criterion for the
# gradient
sqrt(sum(joint_ms_lb_gr(comp_obj_wrong, opt_out_wrong$par)^2))
#> [1] 0.09974
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.9794

opt_out_wrong$info # convergence code (0 == 'OK')
#> [1] 0
opt_out$info # convergence code (0 == 'OK')
#> [1] 0

print(-opt_out_wrong$value, digits = 8) # maximum lower bound value
#> [1] -8351.736
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -8333.1739

opt_out_wrong$counts
#> function gradient     n_cg 
#>     1456     1015     3276
opt_out$counts
#> function gradient     n_cg 
#>     2288     1571    14890

# compare the estimates with the actual values. Start with the fixed effects
fmt_ests_wrong <- joint_ms_format(comp_obj_wrong, opt_out_wrong$par)
fmt_ests <- joint_ms_format(comp_obj, opt_out$par)

# the parameters for the first marker
mapply(rbind, SIMPLIFY = FALSE,
       wrong = fmt_ests_wrong$markers[[1]], 
       right = fmt_ests      $markers[[1]])
#> $fixef
#>                    
#> wrong -0.5421 2.004
#> right -0.4757 2.006
#> 
#> $fixef_vary
#>                         
#> wrong 1.679 1.369 -1.675
#> right 1.540 1.235 -2.052

fixef_marker[[1]] # true values
#> [1] -0.5  2.0
fixef_vary_marker[[1]] # true values
#> [1]  1.4  1.2 -2.1

# the parameters for the second marker
mapply(rbind, SIMPLIFY = FALSE,
       wrong = fmt_ests_wrong$markers[[2]], 
       right = fmt_ests      $markers[[2]])
#> $fixef
#>            
#> wrong 1.021
#> right 1.001
#> 
#> $fixef_vary
#>                      
#> wrong 0.6660 -0.02026
#> right 0.5442 -0.01898

fixef_marker[[2]] # true values
#> [1] 1
fixef_vary_marker[[2]] # true values
#> [1]  0.50 -0.02

# the fixed effects for the survival outcome and the association parameters
# for the terminal event
mapply(rbind, SIMPLIFY = FALSE,
       wrong = fmt_ests_wrong$survival[[1]], 
       right = fmt_ests      $survival[[1]])
#> $fixef
#>                     
#> wrong -1.0630 0.2639
#> right -0.9873 0.2905
#> 
#> $fixef_vary
#>                                    
#> wrong 0.3087 -0.2337 -0.6420 0.1967
#> right 0.3493  0.1313 -0.4404 0.3969
#> 
#> $associations
#>                     
#> wrong 0.8757 -0.4008
#> right 0.7787 -0.4141

fixef_surv[[1]]
#> [1] -1.00  0.25
fixef_vary_surv[[1]]
#> [1]  0.50  0.10 -0.20  0.11
associations[[1]]
#> [1]  0.6 -0.4

# same for the observation process
mapply(rbind, SIMPLIFY = FALSE,
       wrong = fmt_ests_wrong$survival[[2]], 
       right = fmt_ests      $survival[[2]])
#> $fixef
#>             
#> wrong 0.3081
#> right 0.2538
#> 
#> $fixef_vary
#>                     
#> wrong -1.089 -0.3681
#> right -1.222 -0.2858
#> 
#> $associations
#>                     
#> wrong -0.6652 0.2161
#> right -0.6633 0.2156

fixef_surv[[2]]
#> [1] 0.2
fixef_vary_surv[[2]]
#> [1] -1.00 -0.25
associations[[2]]
#> [1] -0.7  0.2

# the parameters for covariance matrix of the random effects
fmt_ests_wrong$vcov$vcov_vary
#>          [,1]     [,2]     [,3]     [,4]
#> [1,]  0.45734  0.16471 -0.05149  0.03625
#> [2,]  0.16471  1.76579 -0.28499 -0.06962
#> [3,] -0.05149 -0.28499  0.37750  0.07415
#> [4,]  0.03625 -0.06962  0.07415  0.09781
fmt_ests      $vcov$vcov_vary
#>          [,1]     [,2]     [,3]     [,4]
#> [1,]  0.47654  0.14431 -0.03266  0.01979
#> [2,]  0.14431  1.83347 -0.27366 -0.03223
#> [3,] -0.03266 -0.27366  0.36566  0.08044
#> [4,]  0.01979 -0.03223  0.08044  0.10748
vcov_vary # the true values
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.35  0.08 -0.05  0.01
#> [2,]  0.08  1.92 -0.24 -0.04
#> [3,] -0.05 -0.24  0.32  0.09
#> [4,]  0.01 -0.04  0.09  0.12

norm(fmt_ests_wrong$vcov$vcov_vary - vcov_vary, "F")
#> [1] 0.2473
norm(fmt_ests      $vcov$vcov_vary - vcov_vary, "F")
#> [1] 0.1933

# the parameters for the error term covariance matrix
fmt_ests_wrong$vcov$vcov_marker
#>         [,1]    [,2]
#> [1,] 0.35539 0.08921
#> [2,] 0.08921 0.15723
fmt_ests      $vcov$vcov_marker
#>         [,1]    [,2]
#> [1,] 0.35467 0.08864
#> [2,] 0.08864 0.15697
vcov_marker
#>      [,1] [,2]
#> [1,] 0.36 0.10
#> [2,] 0.10 0.16

norm(fmt_ests_wrong$vcov$vcov_marker - vcov_marker, "F")
#> [1] 0.01618
norm(fmt_ests      $vcov$vcov_marker - vcov_marker, "F")
#> [1] 0.01719

# the parameters for the frailty covariance matrix
fmt_ests_wrong$vcov$vcov_surv
#>         [,1]    [,2]
#> [1,] 0.01633 0.02719
#> [2,] 0.02719 0.05935
fmt_ests      $vcov$vcov_surv
#>         [,1]    [,2]
#> [1,] 0.01389 0.02352
#> [2,] 0.02352 0.05506
vcov_surv
#>        [,1]   [,2]
#> [1,] 0.0400 0.0225
#> [2,] 0.0225 0.0625

norm(fmt_ests_wrong$vcov$vcov_surv - vcov_surv, "F")
#> [1] 0.02478
norm(fmt_ests      $vcov$vcov_surv - vcov_surv, "F")
#> [1] 0.02719
```

### Two Markers, the Observation Time Process, a Terminal Event, and Mixed Dependencies

We simulate and fit a model in this section where we have two markers, a
recurrent event which is the observation times, and a terminal event
like in the previous section. Unlike in the previous section, we let the
two time-to-event outcomes depend on the integral or derivative of the
deviation of the mean marker and the current value.

``` r
# settings for the simulation
library(VAJointSurv)
g1 <- ns_term(knots = c(3.33, 6.67), Boundary.knots = c(0, 10))
g2 <- poly_term(degree = 2, raw = TRUE)
g_funcs <- list(g1$eval, g2$eval)

m1 <- ns_term(knots = numeric(), Boundary.knots = c(0, 10), intercept = TRUE)
m2 <- poly_term(degree = 1, intercept = TRUE, raw = TRUE)
m_funcs <- list(m1$eval, m2$eval)

fixef_vary_marker <- list(c(1.4, 1.2, -2.1), c(.5, -.02)) # beta
fixef_marker <- list(c(-.5, 2), 1) # gamma

# Psi
vcov_vary <- structure(c(0.35, 0.08, -0.05, 0.01, 0.08, 1.92, -0.24, -0.04,
                   -0.05, -0.24, 0.32, 0.09, 0.01, -0.04, 0.09, 0.12),
                 .Dim = c(4L, 4L))
vcov_marker <- matrix(c(.6^2, .1, .1, .4^2), 2)

# plot the markers' mean curve
par(mar = c(5, 5, 1, 1))
plot_marker(
  time_fixef = g1, time_rng = m1,
  fixef_vary = fixef_vary_marker[[1]], x_range = c(0, 10),
  vcov_vary = vcov_vary[1:2, 1:2], ylab = "Marker 1")
```

![](man/figures/README-obs_mixed_process_markers_and_recurrent-1.png)<!-- -->

``` r
plot_marker(
  time_fixef = g2, time_rng = m2,
  fixef_vary = fixef_vary_marker[[2]], x_range = c(0, 10),
  vcov_vary = vcov_vary[3:4, 3:4], ylab = "Marker 2")
```

![](man/figures/README-obs_mixed_process_markers_and_recurrent-2.png)<!-- -->

``` r
# the survival parameters
vcov_surv <- matrix(c(.2^2, .15^2, .15^2, .25^2), 2) # Xi 

fixef_surv <- list(c(-1, .25), .2)
# needs more parameters now
associations <- list(c(-.4, .6, -.4), c(-.7, .2, -1)) 
fixef_vary_surv <- list(c(.5, .1, -.2, .11),
                          c(-1, -.25))

# specify the dependence with the random effect from the markers
ders <- list(
  list(c(-1L, 0L), # cumulative and present value 
       0L),        # present value
  list(0L,         # present value
       c(0L, 1L))) # present value and derivative

b_funcs <- list(
  function(x) bs(x, knots = 5, Boundary.knots = c(0, 10)),
  function(x) ns(x, knots = 5, Boundary.knots = c(0, 10)))

# plot the log hazard with the 25%, 50% and 75% quantiles
library(SimSurvNMarker)

plot_surv(
  time_fixef = bs_term(knots = 5, Boundary.knots = c(0, 10)),
  time_rng = list(m1, m2), ders = ders[[1]],
  x_range = c(0, 10), fixef_vary = fixef_vary_surv[[1]],
  vcov_vary = vcov_vary, frailty_var = vcov_surv[1, 1], ps = c(.25, .5, .75),
  associations = associations[[1]], log_hazard_shift = fixef_surv[[1]][1],
  ylab = "Terminal event")
```

![](man/figures/README-obs_mixed_process_markers_and_recurrent-3.png)<!-- -->

``` r
plot_surv(
  time_fixef = ns_term(knots = 5, Boundary.knots = c(0, 10)),
  time_rng = list(m1, m2), ders = ders[[2]],
  x_range = c(0, 10), fixef_vary = fixef_vary_surv[[2]],
  vcov_vary = vcov_vary, frailty_var = vcov_surv[2, 2], ps = c(.25, .5, .75),
  associations = associations[[2]], log_hazard_shift = fixef_surv[[2]][1],
  ylab = "Observation process")
```

![](man/figures/README-obs_mixed_process_markers_and_recurrent-4.png)<!-- -->

``` r
library(mvtnorm)
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  gl_dat <- get_gl_rule(100L)
  dat <- lapply(1:n_ids, function(id){
    # sample the terminal event time and the censoring time
    cens <- min(rexp(1, rate = 1/10), 10)
    U <- drop(rmvnorm(1, sigma = vcov_vary))

    frailties <- drop(rmvnorm(1, sigma = vcov_surv))
    Z1 <- c(1, runif(1, -1, 1))
    log_haz_offset <- sum(Z1 * fixef_surv[[1]]) + frailties[1]

    # assign the conditional hazard function
    expansion <- function(x, b_func)
      cbind(b_func(x), 
            t(U[1:2] %*% m_funcs[[1]](x, ders[[1]][[1]][1])),
            t(U[1:2] %*% m_funcs[[1]](x, ders[[1]][[1]][2])),
            t(U[3:4] %*% m_funcs[[2]](x, ders[[1]][[2]][1])))
    surv_func <- function(ti, fixef_vary_surv, associations, b_func){
      formals(expansion)$b_func <- b_func
      eval_surv_base_fun(
        ti = ti, omega = c(fixef_vary_surv, associations), b_func = expansion,
        gl_dat = gl_dat, delta = log_haz_offset)
    }

    # sample the survival time
    rng <- runif(1)
    root_func <- function(x, rng)
      rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[1]],
                      associations = associations[[1]], b_func = b_funcs[[1]])

    if(root_func(cens, rng) < 0){
      # the observation is censored
      y_terminal <- cens
      event <- 0
    } else {
      # find the event time
      root <- uniroot(root_func, c(0, cens), tol = 1e-6, rng = rng)
      y_terminal <- root$root
      event <- 1

    }

    terminal_outcome <- cbind(y = y_terminal, event = event, Z1 = Z1[2],
                              id = id)

    # clean up
    rm(list = setdiff(ls(), c("y_terminal", "terminal_outcome", "expansion",
                              "surv_func", "frailties", "U", "id")))

    # simulate the observation times
    Z2 <- 1
    log_haz_offset <- sum(Z2 * fixef_surv[[2]]) + frailties[2]

    expansion <- function(x, b_func)
      cbind(b_func(x), 
            t(U[1:2] %*% m_funcs[[1]](x, ders[[2]][[1]][1])),
            t(U[3:4] %*% m_funcs[[2]](x, ders[[2]][[2]][1])),
            t(U[3:4] %*% m_funcs[[2]](x, ders[[2]][[2]][2])))
    root_func <- function(x, left_trunc_surv, rng)
      rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[2]],
                      associations = associations[[2]], b_func = b_funcs[[2]]) /
      left_trunc_surv

    max_sample <- 1000L
    left_trunc_surv <- 1
    Z2 <- matrix(rep(Z2, each = max_sample), max_sample)
    event <- y <- lf_trunc <- rep(NA_real_, max_sample)
    lf_trunc_i <- 0
    for(i in 1:max_sample){
      # sample a random uniform variable and invert the survival function
      rng_i <- runif(1)
      lf_trunc[i] <- lf_trunc_i

      if(root_func(y_terminal, left_trunc_surv, rng_i) < 0){
        # the observation is right-censored and we can exit
        y[i] <- y_terminal
        event[i] <- 0
        break
      }

      # we need to invert the survival function to find the observation time
      root <- uniroot(root_func, c(lf_trunc_i, y_terminal), tol = 1e-6,
                      left_trunc_surv = left_trunc_surv, rng = rng_i)
      lf_trunc_i <- y[i] <- root$root
      event[i] <- 1
      left_trunc_surv <- surv_func(
        y[i], fixef_vary_surv = fixef_vary_surv[[2]], associations = associations[[2]],
        b_func = b_funcs[[2]])
    }

    colnames(Z2) <- paste0("Z", 1:NCOL(Z2) - 1L)
    obs_process <- cbind(lf_trunc = lf_trunc[1:i], y = y[1:i],
                         event = event[1:i], Z2[1:i, -1, drop = FALSE],
                         id = id)

    # clean up
    rm(list = setdiff(ls(), c("terminal_outcome", "U", "id",
                              "obs_process")))

    # sample the number of outcomes and the fixed effect covariates
    obs_time <- c(0, obs_process[obs_process[, "event"] == 1, "y"])
    n_obs <- length(obs_time)
    X1 <- cbind(1, rnorm(n_obs))
    X2 <- matrix(1, n_obs)
    colnames(X1) <- paste0("X1_", 1:NCOL(X1) - 1L)
    colnames(X2) <- paste0("X2_", 1:NCOL(X2) - 1L)
    X <- list(X1, X2)

    # sample the outcomes
    eta <- sapply(1:2, function(i)
      X[[i]] %*% fixef_marker[[i]] +
        drop(t(g_funcs[[i]](obs_time))) %*% fixef_vary_marker[[i]] +
        drop(t(m_funcs[[i]](obs_time))) %*% U[1:2 + (i == 2) * 2])

    ys <- eta + rmvnorm(n_obs, sigma = vcov_marker)
    colnames(ys) <- paste0("Y", 1:2)

    # mask some observations
    do_mask <- sample.int(3L, n_obs, replace = TRUE)
    ys[do_mask == 2, 1] <- NA
    ys[do_mask == 3, 2] <- NA

    X <- do.call(cbind, lapply(X, function(x) x[, -1, drop = FALSE]))
    marker_data <- cbind(ys, X, time = obs_time, id = id)

    return(list(marker_data = marker_data, obs_process = obs_process,
                terminal_outcome = terminal_outcome))
  })

  # combine the data and return
  marker_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "marker_data")))
  marker_data$id <- as.integer(marker_data$id)
  # the order does not matter
  marker_data <- marker_data[sample.int(NROW(marker_data)), ]

  obs_process <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "obs_process")))
  obs_process$id <- as.integer(obs_process$id)
  # the order does not matter
  obs_process <- obs_process[sample.int(NROW(obs_process)), ]

  terminal_outcome <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "terminal_outcome")))
  terminal_outcome$id <- as.integer(terminal_outcome$id)
  # the order does not matter
  terminal_outcome <- terminal_outcome[sample.int(NROW(terminal_outcome)), ]

  list(marker_data = marker_data, obs_process = obs_process,
       terminal_outcome = terminal_outcome)
}

# sample a moderately large data set
set.seed(2)
dat <- sim_dat(1000L)

# we show a few properties of the data below
mean(dat$terminal_outcome$event) # mean event rate
#> [1] 0.77
sum(dat$obs_process$event) # number of observed markers less the individuals
#> [1] 2425
NROW(dat$marker_data) # number of observed markers less the individuals
#> [1] 3425

# distribution of observed marker per individual
proportions(table(table(dat$obs_process$id)))
#> 
#>     1     2     3     4     5     6     7     8     9    10    11    12    13    14    16    17    27    29    31    38 
#> 0.319 0.190 0.141 0.109 0.078 0.041 0.033 0.027 0.013 0.012 0.015 0.003 0.004 0.005 0.003 0.003 0.001 0.001 0.001 0.001

# show data for one individual
subset(dat$marker_data, id == 1)
#>          Y1    Y2     X1_1   time id
#> X.6      NA 8.382 -0.56843 7.8724  1
#> X.3  0.7605    NA -0.97556 3.6495  1
#> X.7 -3.5564 9.325 -1.03442 8.4949  1
#> X.2 -0.5523    NA -1.11525 0.4897  1
#> X.9 -3.9194 9.798 -0.98139 9.5827  1
#> X.1      NA 1.726 -0.36632 0.4271  1
#> X        NA 1.987 -0.02815 0.0000  1
#> X.4  4.1065 6.374  1.12102 4.9005  1
#> X.8 -1.5241 9.405 -0.53302 8.6015  1
#> X.5  0.3933    NA  0.32444 7.5675  1
subset(dat$obs_process, id == 1)
#>    lf_trunc       y event id
#> 8    8.4949  8.6015     1  1
#> 9    8.6015  9.5827     1  1
#> 2    0.4271  0.4897     1  1
#> 1    0.0000  0.4271     1  1
#> 7    7.8724  8.4949     1  1
#> 10   9.5827 10.0000     0  1
#> 5    4.9005  7.5675     1  1
#> 6    7.5675  7.8724     1  1
#> 4    3.6495  4.9005     1  1
#> 3    0.4897  3.6495     1  1
subset(dat$terminal_outcome, id == 1)
#>    y event      Z1 id
#> 1 10     0 -0.6384  1

# estimate the model with this package. Get the object we need for the
# optimization
marker_1 <- marker_term(
  Y1 ~ X1_1, id = id, subset(dat$marker_data, !is.na(Y1)),
  time_fixef = ns_term(time, knots = c(3.33, 6.67), Boundary.knots = c(0, 10)),
  time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 10),
                     intercept = TRUE))
marker_2 <- marker_term(
  Y2 ~ 1, id = id, subset(dat$marker_data, !is.na(Y2)),
  time_fixef = poly_term(time, degree = 2, raw = TRUE),
  time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

library(survival)
surv_terminal <- surv_term(
  Surv(y, event) ~ Z1, id = id, dat$terminal_outcome,
  time_fixef = bs_term(y, knots = 5, Boundary.knots = c(0, 10)), 
  ders = ders[[1]], with_frailty = TRUE)
surv_obs <- surv_term(
  Surv(lf_trunc, y, event) ~ 1, id = id, dat$obs_process,
  time_fixef = ns_term(y, knots = 5, Boundary.knots = c(0, 10)), 
  ders = ders[[2]], with_frailty = TRUE)

comp_obj <- joint_ms_ptr(markers = list(marker_1, marker_2),
                         survival_terms = list(surv_terminal, surv_obs),
                         max_threads = 4L)
rm(marker_1, marker_2, surv_terminal, surv_obs)

# get the starting values
system.time(start_val <- joint_ms_start_val(comp_obj, gr_tol = .1))
#> Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge with max|grad| = 0.00249365 (tol =
#> 0.002, component 1)
#>    user  system elapsed 
#>  18.689   0.012   4.958

# lower bound at the starting values
print(-attr(start_val, "value"), digits = 8)
#> [1] -7950.3559

# check that the gradient is correct
f <- function(x){
  start_val[seq_along(x)] <- x
  joint_ms_lb(comp_obj, start_val)
}

all.equal(numDeriv::grad(f, head(start_val, 39 + 2 * 27)),
          head(joint_ms_lb_gr(comp_obj, start_val), 39 + 2 * 27), 
          tolerance = 1e-6)
#> [1] TRUE

# find the maximum lower bound estimate
system.time(opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 1000L,
                                    pre_method = 3L, cg_tol = .2, c2 = .1,
                                    gr_tol = .1))
#>    user  system elapsed 
#> 207.052   0.064  51.780

# we set gr_tol in the call so this is the convergence criterion for the 
# gradient
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2))
#> [1] 0.09917
opt_out$info # convergence code (0 == 'OK')
#> [1] 0
print(-opt_out$value, digits = 8) # maximum lower bound value
#> [1] -7596.2563
opt_out$counts
#> function gradient     n_cg 
#>     3570     2531    15124

# compare the estimates with the actual values. Start with the fixed effects
fmt_ests <- joint_ms_format(comp_obj, opt_out$par)

# the parameters for the first marker
fmt_ests$markers[[1]]
#> $fixef
#>                 
#> -0.4302  2.0078 
#> 
#> $fixef_vary
#>                      
#>  1.355  1.218 -1.801

fixef_marker[[1]] # true values
#> [1] -0.5  2.0
fixef_vary_marker[[1]] # true values
#> [1]  1.4  1.2 -2.1

# the parameters for the second marker
fmt_ests$markers[[2]]
#> $fixef
#>        
#> 0.9901 
#> 
#> $fixef_vary
#>                   
#>  0.53128 -0.02203
fixef_marker[[2]] # true values
#> [1] 1
fixef_vary_marker[[2]] # true values
#> [1]  0.50 -0.02

# the fixed effects for the survival outcome and the association parameters
# for the terminal event
fmt_ests$survival[[1]]
#> $fixef
#>                 
#> -1.0845  0.2985 
#> 
#> $fixef_vary
#>                                     
#>  0.73891  0.07372 -0.11053 -0.43122 
#> 
#> $associations
#>                         
#> -0.4091  0.5455 -0.4309
fixef_surv[[1]]
#> [1] -1.00  0.25
fixef_vary_surv[[1]]
#> [1]  0.50  0.10 -0.20  0.11
associations[[1]]
#> [1] -0.4  0.6 -0.4

# same for the observation process
fmt_ests$survival[[2]]
#> $fixef
#>        
#> 0.2204 
#> 
#> $fixef_vary
#>                 
#> -1.2803 -0.2924 
#> 
#> $associations
#>                         
#> -0.7104  0.1854 -0.7895
fixef_surv[[2]]
#> [1] 0.2
fixef_vary_surv[[2]]
#> [1] -1.00 -0.25
associations[[2]]
#> [1] -0.7  0.2 -1.0

# the parameters for covariance matrix of the random effects
fmt_ests$vcov$vcov_vary
#>          [,1]    [,2]     [,3]     [,4]
#> [1,]  0.34378  0.2413 -0.07396 -0.03169
#> [2,]  0.24130  2.1631 -0.26015 -0.11240
#> [3,] -0.07396 -0.2602  0.27222  0.10383
#> [4,] -0.03169 -0.1124  0.10383  0.12940
vcov_vary # the true values
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  0.35  0.08 -0.05  0.01
#> [2,]  0.08  1.92 -0.24 -0.04
#> [3,] -0.05 -0.24  0.32  0.09
#> [4,]  0.01 -0.04  0.09  0.12

# the parameters for the error term covariance matrix
fmt_ests$vcov$vcov_marker
#>        [,1]   [,2]
#> [1,] 0.3332 0.0930
#> [2,] 0.0930 0.1621
vcov_marker
#>      [,1] [,2]
#> [1,] 0.36 0.10
#> [2,] 0.10 0.16

# the parameters for the frailty covariance matrix
fmt_ests$vcov$vcov_surv
#>          [,1]     [,2]
#> [1,]  0.06048 -0.04688
#> [2,] -0.04688  0.03692
vcov_surv
#>        [,1]   [,2]
#> [1,] 0.0400 0.0225
#> [2,] 0.0225 0.0625
```

## Technical Details

We provide a few technical details in this section. The concatenated
coefficient is given by

![\\vec\\theta = \\begin{pmatrix}
  \\vec\\gamma \\\\ \\vec\\beta\\\\
  \\vec\\omega\_1 \\\\ \\vec\\delta\_1 \\\\ \\vec\\alpha\_1 \\\\
  \\vdots \\\\
  \\vec\\omega\_H \\\\ \\vec\\delta\_H \\\\ \\vec\\alpha\_H \\\\
    \\text{vec}(\\Sigma) \\\\ \\text{vec}(\\Psi) \\\\ \\text{vec}(\\Xi)
\\end{pmatrix}](https://render.githubusercontent.com/render/math?math=%5Cvec%5Ctheta%20%3D%20%5Cbegin%7Bpmatrix%7D%0A%20%20%5Cvec%5Cgamma%20%5C%5C%20%5Cvec%5Cbeta%5C%5C%0A%20%20%5Cvec%5Comega_1%20%5C%5C%20%5Cvec%5Cdelta_1%20%5C%5C%20%5Cvec%5Calpha_1%20%5C%5C%0A%20%20%5Cvdots%20%5C%5C%0A%20%20%5Cvec%5Comega_H%20%5C%5C%20%5Cvec%5Cdelta_H%20%5C%5C%20%5Cvec%5Calpha_H%20%5C%5C%0A%20%20%20%20%5Ctext%7Bvec%7D%28%5CSigma%29%20%5C%5C%20%5Ctext%7Bvec%7D%28%5CPsi%29%20%5C%5C%20%5Ctext%7Bvec%7D%28%5CXi%29%0A%5Cend%7Bpmatrix%7D "\vec\theta = \begin{pmatrix}
  \vec\gamma \\ \vec\beta\\
  \vec\omega_1 \\ \vec\delta_1 \\ \vec\alpha_1 \\
  \vdots \\
  \vec\omega_H \\ \vec\delta_H \\ \vec\alpha_H \\
    \text{vec}(\Sigma) \\ \text{vec}(\Psi) \\ \text{vec}(\Xi)
\end{pmatrix}")

where
![\\text{vec}(\\cdot)](https://render.githubusercontent.com/render/math?math=%5Ctext%7Bvec%7D%28%5Ccdot%29 "\text{vec}(\cdot)")
is the vectorization function that stacks the column of a matrix on top
of each other.
![\\vec\\theta](https://render.githubusercontent.com/render/math?math=%5Cvec%5Ctheta "\vec\theta")
is with the parameters from the variational approximation. We work with
log Cholesky decompositions of the covariance matrices in practice. The
`joint_ms_format` is a utility functions which returns a list with each
of parameter separately.

There are different types of terms in the lower bound in the GVA. We
cover the types below. The lower bound has many parameters even with a
moderate amount of parameters. However, we quickly optimize the lower
bound using the [psqn package](https://github.com/boennecd/psqn).

## Kullback–Leibler Divergence Term

In the GVA, we assume that the conditional distribution of
![(\\vec U\_i^\\top, \\vec\\xi\_i^\\top)^\\top](https://render.githubusercontent.com/render/math?math=%28%5Cvec%20U_i%5E%5Ctop%2C%20%5Cvec%5Cxi_i%5E%5Ctop%29%5E%5Ctop "(\vec U_i^\top, \vec\xi_i^\top)^\top")
is a normal distribution with mean
![\\vec\\zeta\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%5Czeta_i "\vec\zeta_i")
and covariance matrix
![\\Omega\_i](https://render.githubusercontent.com/render/math?math=%5COmega_i "\Omega_i").
One of the terms in the lower bound of the GVA is the Kullback–Leibler
(KL) divergence term between the unconditional distribution of
![(\\vec U\_i^\\top, \\vec\\xi\_i^\\top)^\\top](https://render.githubusercontent.com/render/math?math=%28%5Cvec%20U_i%5E%5Ctop%2C%20%5Cvec%5Cxi_i%5E%5Ctop%29%5E%5Ctop "(\vec U_i^\top, \vec\xi_i^\top)^\top")
and the assumed conditional distribution. This term is given by

![\\begin{multline\*}
\\frac{1}{2}\\Big(\\log\\lvert\\Omega\_i\\rvert - \\log\\lvert\\Psi\\rvert - \\log\\lvert\\Xi\\rvert -\\vec\\zeta\_{i,1:R}^\\top\\Psi^{-1}\\vec\\zeta\_{i,1:R}-\\vec\\zeta\_{i,(-1:R)}^\\top\\Xi^{-1}\\vec\\zeta\_{i,(-1:R)} \\\\- \\text{tr}\\Omega\_{i,1:R,1:R}\\Psi^{-1}- \\text{tr}\\Omega\_{i,(-1:R),(-1:R)}\\Xi^{-1}+ R + H\\Big)\\end{multline\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bmultline%2A%7D%0A%5Cfrac%7B1%7D%7B2%7D%5CBig%28%5Clog%5Clvert%5COmega_i%5Crvert%20-%20%5Clog%5Clvert%5CPsi%5Crvert%20-%20%5Clog%5Clvert%5CXi%5Crvert%20-%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%5E%5Ctop%5CPsi%5E%7B-1%7D%5Cvec%5Czeta_%7Bi%2C1%3AR%7D-%5Cvec%5Czeta_%7Bi%2C%28-1%3AR%29%7D%5E%5Ctop%5CXi%5E%7B-1%7D%5Cvec%5Czeta_%7Bi%2C%28-1%3AR%29%7D%20%5C%5C-%20%5Ctext%7Btr%7D%5COmega_%7Bi%2C1%3AR%2C1%3AR%7D%5CPsi%5E%7B-1%7D-%20%5Ctext%7Btr%7D%5COmega_%7Bi%2C%28-1%3AR%29%2C%28-1%3AR%29%7D%5CXi%5E%7B-1%7D%2B%20R%20%2B%20H%5CBig%29%5Cend%7Bmultline%2A%7D "\begin{multline*}
\frac{1}{2}\Big(\log\lvert\Omega_i\rvert - \log\lvert\Psi\rvert - \log\lvert\Xi\rvert -\vec\zeta_{i,1:R}^\top\Psi^{-1}\vec\zeta_{i,1:R}-\vec\zeta_{i,(-1:R)}^\top\Xi^{-1}\vec\zeta_{i,(-1:R)} \\- \text{tr}\Omega_{i,1:R,1:R}\Psi^{-1}- \text{tr}\Omega_{i,(-1:R),(-1:R)}\Xi^{-1}+ R + H\Big)\end{multline*}")

where
![\\text{tr}(\\cdot)](https://render.githubusercontent.com/render/math?math=%5Ctext%7Btr%7D%28%5Ccdot%29 "\text{tr}(\cdot)")
returns the trace of a matrix. Thus, the derivatives w.r.t.
![\\Omega\_i](https://render.githubusercontent.com/render/math?math=%5COmega_i "\Omega_i"),
![\\Psi](https://render.githubusercontent.com/render/math?math=%5CPsi "\Psi"),
![\\Xi](https://render.githubusercontent.com/render/math?math=%5CXi "\Xi"),
![\\vec\\zeta\_i](https://render.githubusercontent.com/render/math?math=%5Cvec%5Czeta_i "\vec\zeta_i")
are respectively

![\\frac{1}{2}\\left(\\Omega\_i^{-1} - \\begin{pmatrix} \\Psi^{-1} & 0 \\\\ 0 & \\Xi^{-1}\\end{pmatrix}\\right)](https://render.githubusercontent.com/render/math?math=%5Cfrac%7B1%7D%7B2%7D%5Cleft%28%5COmega_i%5E%7B-1%7D%20-%20%5Cbegin%7Bpmatrix%7D%20%5CPsi%5E%7B-1%7D%20%26%200%20%5C%5C%200%20%26%20%5CXi%5E%7B-1%7D%5Cend%7Bpmatrix%7D%5Cright%29 "\frac{1}{2}\left(\Omega_i^{-1} - \begin{pmatrix} \Psi^{-1} & 0 \\ 0 & \Xi^{-1}\end{pmatrix}\right)")

![\\frac{1}{2}\\Psi^{-1}(\\vec\\zeta\_{i,1:R}\\vec\\zeta\_{i,1:R}^\\top + \\Omega\_{i,1:R,1:R} - \\Psi)\\Psi^{-1}](https://render.githubusercontent.com/render/math?math=%5Cfrac%7B1%7D%7B2%7D%5CPsi%5E%7B-1%7D%28%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%5E%5Ctop%20%2B%20%5COmega_%7Bi%2C1%3AR%2C1%3AR%7D%20-%20%5CPsi%29%5CPsi%5E%7B-1%7D "\frac{1}{2}\Psi^{-1}(\vec\zeta_{i,1:R}\vec\zeta_{i,1:R}^\top + \Omega_{i,1:R,1:R} - \Psi)\Psi^{-1}")

![\\frac{1}{2}\\Xi^{-1}(\\vec\\zeta\_{i,(-1:R)}\\vec\\zeta\_{i,(-1:R)}^\\top+ \\Omega\_{i,(-1:R),(-1:R)} - \\Xi)\\Xi^{-1}.](https://render.githubusercontent.com/render/math?math=%5Cfrac%7B1%7D%7B2%7D%5CXi%5E%7B-1%7D%28%5Cvec%5Czeta_%7Bi%2C%28-1%3AR%29%7D%5Cvec%5Czeta_%7Bi%2C%28-1%3AR%29%7D%5E%5Ctop%2B%20%5COmega_%7Bi%2C%28-1%3AR%29%2C%28-1%3AR%29%7D%20-%20%5CXi%29%5CXi%5E%7B-1%7D. "\frac{1}{2}\Xi^{-1}(\vec\zeta_{i,(-1:R)}\vec\zeta_{i,(-1:R)}^\top+ \Omega_{i,(-1:R),(-1:R)} - \Xi)\Xi^{-1}.")

![-\\begin{pmatrix}
  \\Psi^{-1} & 0 \\\\
  0 & \\Xi^{-1} 
\\end{pmatrix}\\vec\\zeta\_i](https://render.githubusercontent.com/render/math?math=-%5Cbegin%7Bpmatrix%7D%0A%20%20%5CPsi%5E%7B-1%7D%20%26%200%20%5C%5C%0A%20%200%20%26%20%5CXi%5E%7B-1%7D%20%0A%5Cend%7Bpmatrix%7D%5Cvec%5Czeta_i "-\begin{pmatrix}
  \Psi^{-1} & 0 \\
  0 & \Xi^{-1} 
\end{pmatrix}\vec\zeta_i")

## Marker Terms

The log conditional density term of observation
![j](https://render.githubusercontent.com/render/math?math=j "j") of
individual
![i](https://render.githubusercontent.com/render/math?math=i "i") at
time
![s\_{ij}](https://render.githubusercontent.com/render/math?math=s_%7Bij%7D "s_{ij}")
is

![\\int \\log\\left(\\phi^{(L)}\\left(\\vec y\_{ij}; X\_i\\vec\\gamma + G(s\_{ij})\\vec\\beta + M(s\_{ij})\\vec w, \\Sigma\\right)\\right)\\phi^{(R)}\\left(\\vec w; \\vec\\zeta\_{i,1:R}, \\Omega\_{i,1:R,1:R}\\right) d\\vec w](https://render.githubusercontent.com/render/math?math=%5Cint%20%5Clog%5Cleft%28%5Cphi%5E%7B%28L%29%7D%5Cleft%28%5Cvec%20y_%7Bij%7D%3B%20X_i%5Cvec%5Cgamma%20%2B%20G%28s_%7Bij%7D%29%5Cvec%5Cbeta%20%2B%20M%28s_%7Bij%7D%29%5Cvec%20w%2C%20%5CSigma%5Cright%29%5Cright%29%5Cphi%5E%7B%28R%29%7D%5Cleft%28%5Cvec%20w%3B%20%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%2C%20%5COmega_%7Bi%2C1%3AR%2C1%3AR%7D%5Cright%29%20d%5Cvec%20w "\int \log\left(\phi^{(L)}\left(\vec y_{ij}; X_i\vec\gamma + G(s_{ij})\vec\beta + M(s_{ij})\vec w, \Sigma\right)\right)\phi^{(R)}\left(\vec w; \vec\zeta_{i,1:R}, \Omega_{i,1:R,1:R}\right) d\vec w")

which gives

![\\begin{multline\*}-\\frac{L}{2}\\log 2\\pi  -\\frac {1}{2} \\log\\lvert\\Sigma \\rvert - \\frac{1}{2}\\left(\\vec y\_{ij} - X\_i\\vec\\gamma - G(s\_{ij})\\vec\\beta - M(s\_{ij})\\vec\\zeta\_{i,1:R}\\right)^\\top \\\\ \\Sigma^{-1}\\left(\\vec y\_{ij} - X\_i\\vec\\gamma - G(s\_{ij})\\vec\\beta - M(s\_{ij})\\vec\\zeta\_{i,1:R}\\right) - \\frac{1}{2}\\text{tr}\\Sigma^{-1}M(s\_{ij})\\Omega\_{i,1:R,1:R}M(s\_{ij})^\\top. \\end{multline\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bmultline%2A%7D-%5Cfrac%7BL%7D%7B2%7D%5Clog%202%5Cpi%20%20-%5Cfrac%20%7B1%7D%7B2%7D%20%5Clog%5Clvert%5CSigma%20%5Crvert%20-%20%5Cfrac%7B1%7D%7B2%7D%5Cleft%28%5Cvec%20y_%7Bij%7D%20-%20X_i%5Cvec%5Cgamma%20-%20G%28s_%7Bij%7D%29%5Cvec%5Cbeta%20-%20M%28s_%7Bij%7D%29%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%5Cright%29%5E%5Ctop%20%5C%5C%20%5CSigma%5E%7B-1%7D%5Cleft%28%5Cvec%20y_%7Bij%7D%20-%20X_i%5Cvec%5Cgamma%20-%20G%28s_%7Bij%7D%29%5Cvec%5Cbeta%20-%20M%28s_%7Bij%7D%29%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%5Cright%29%20-%20%5Cfrac%7B1%7D%7B2%7D%5Ctext%7Btr%7D%5CSigma%5E%7B-1%7DM%28s_%7Bij%7D%29%5COmega_%7Bi%2C1%3AR%2C1%3AR%7DM%28s_%7Bij%7D%29%5E%5Ctop.%20%5Cend%7Bmultline%2A%7D "\begin{multline*}-\frac{L}{2}\log 2\pi  -\frac {1}{2} \log\lvert\Sigma \rvert - \frac{1}{2}\left(\vec y_{ij} - X_i\vec\gamma - G(s_{ij})\vec\beta - M(s_{ij})\vec\zeta_{i,1:R}\right)^\top \\ \Sigma^{-1}\left(\vec y_{ij} - X_i\vec\gamma - G(s_{ij})\vec\beta - M(s_{ij})\vec\zeta_{i,1:R}\right) - \frac{1}{2}\text{tr}\Sigma^{-1}M(s_{ij})\Omega_{i,1:R,1:R}M(s_{ij})^\top. \end{multline*}")

## Survival Outcomes

Let
![T\_{i1}, \\dots, T\_{iL}](https://render.githubusercontent.com/render/math?math=T_%7Bi1%7D%2C%20%5Cdots%2C%20T_%7BiL%7D "T_{i1}, \dots, T_{iL}")
be the minimum of the censoring time or the observed event time and
![D\_{i1}, \\dots, D\_{iL}](https://render.githubusercontent.com/render/math?math=D_%7Bi1%7D%2C%20%5Cdots%2C%20D_%7BiL%7D "D_{i1}, \dots, D_{iL}")
be event indicators. The the lower bound terms for the
![k](https://render.githubusercontent.com/render/math?math=k "k")th time
to event is

![\\int \\left(d\_{ik}\\log h\_{ik}(t\_{ik}\\mid \\vec w\_{1:R}, w\_{R + k})-\\int\_0^{t\_{ik}}  h\_{ik}(s\\mid \\vec w\_{1:R}, w\_{R + k})ds \\right)\\phi^{(R + H)}(\\vec w; \\vec\\zeta\_i, \\Omega\_i)d\\vec w.](https://render.githubusercontent.com/render/math?math=%5Cint%20%5Cleft%28d_%7Bik%7D%5Clog%20h_%7Bik%7D%28t_%7Bik%7D%5Cmid%20%5Cvec%20w_%7B1%3AR%7D%2C%20w_%7BR%20%2B%20k%7D%29-%5Cint_0%5E%7Bt_%7Bik%7D%7D%20%20h_%7Bik%7D%28s%5Cmid%20%5Cvec%20w_%7B1%3AR%7D%2C%20w_%7BR%20%2B%20k%7D%29ds%20%5Cright%29%5Cphi%5E%7B%28R%20%2B%20H%29%7D%28%5Cvec%20w%3B%20%5Cvec%5Czeta_i%2C%20%5COmega_i%29d%5Cvec%20w. "\int \left(d_{ik}\log h_{ik}(t_{ik}\mid \vec w_{1:R}, w_{R + k})-\int_0^{t_{ik}}  h_{ik}(s\mid \vec w_{1:R}, w_{R + k})ds \right)\phi^{(R + H)}(\vec w; \vec\zeta_i, \Omega_i)d\vec w.")

Left-truncation is easily handled by replacing the zero in the lower
limit by the left-truncation point. The log hazard term gives the
following lower bound terms

![d\_{ik}\\left(\\vec z\_{ik}^\\top\\vec\\delta\_k + \\omega\_k^\\top \\vec b\_k(t\_{ik}) +  \\vec\\alpha\_k^\\top M(t\_{ik})\\vec\\zeta\_{i,1:R} + \\zeta\_{i,R + k}\\right).](https://render.githubusercontent.com/render/math?math=d_%7Bik%7D%5Cleft%28%5Cvec%20z_%7Bik%7D%5E%5Ctop%5Cvec%5Cdelta_k%20%2B%20%5Comega_k%5E%5Ctop%20%5Cvec%20b_k%28t_%7Bik%7D%29%20%2B%20%20%5Cvec%5Calpha_k%5E%5Ctop%20M%28t_%7Bik%7D%29%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%20%2B%20%5Czeta_%7Bi%2CR%20%2B%20k%7D%5Cright%29. "d_{ik}\left(\vec z_{ik}^\top\vec\delta_k + \omega_k^\top \vec b_k(t_{ik}) +  \vec\alpha_k^\top M(t_{ik})\vec\zeta_{i,1:R} + \zeta_{i,R + k}\right).")

As for the expected cumulative hazard term, assume that we can
interchange the order of integration. Then the lower bound term is

![-\\exp(\\vec z\_{ik}^\\top\\vec\\delta\_k + \\zeta\_{R +k})\\int\_0^{t\_{ij}} \\exp\\left(\\omega\_k^\\top \\vec b\_k(s) +  \\vec\\alpha\_k^\\top M(s)\\vec\\zeta\_{i,1:R} +\\frac{1}{2}  (\\vec\\alpha\_k^\\top, 1) O\_{ik}(s)\\begin{pmatrix} \\vec\\alpha\_k \\\\ 1 \\end{pmatrix}\\right) ds](https://render.githubusercontent.com/render/math?math=-%5Cexp%28%5Cvec%20z_%7Bik%7D%5E%5Ctop%5Cvec%5Cdelta_k%20%2B%20%5Czeta_%7BR%20%2Bk%7D%29%5Cint_0%5E%7Bt_%7Bij%7D%7D%20%5Cexp%5Cleft%28%5Comega_k%5E%5Ctop%20%5Cvec%20b_k%28s%29%20%2B%20%20%5Cvec%5Calpha_k%5E%5Ctop%20M%28s%29%5Cvec%5Czeta_%7Bi%2C1%3AR%7D%20%2B%5Cfrac%7B1%7D%7B2%7D%20%20%28%5Cvec%5Calpha_k%5E%5Ctop%2C%201%29%20O_%7Bik%7D%28s%29%5Cbegin%7Bpmatrix%7D%20%5Cvec%5Calpha_k%20%5C%5C%201%20%5Cend%7Bpmatrix%7D%5Cright%29%20ds "-\exp(\vec z_{ik}^\top\vec\delta_k + \zeta_{R +k})\int_0^{t_{ij}} \exp\left(\omega_k^\top \vec b_k(s) +  \vec\alpha_k^\top M(s)\vec\zeta_{i,1:R} +\frac{1}{2}  (\vec\alpha_k^\top, 1) O_{ik}(s)\begin{pmatrix} \vec\alpha_k \\ 1 \end{pmatrix}\right) ds")

where

![O\_{ik}(s) = \\begin{pmatrix}M(s) & 0 \\\\ 0 & 1\\end{pmatrix}
  \\Omega\_{(1:R, R + k), (1:R, R + k)}
    \\begin{pmatrix}M(s)^\\top & 0 \\\\ 0 & 1\\end{pmatrix}.](https://render.githubusercontent.com/render/math?math=O_%7Bik%7D%28s%29%20%3D%20%5Cbegin%7Bpmatrix%7DM%28s%29%20%26%200%20%5C%5C%200%20%26%201%5Cend%7Bpmatrix%7D%0A%20%20%5COmega_%7B%281%3AR%2C%20R%20%2B%20k%29%2C%20%281%3AR%2C%20R%20%2B%20k%29%7D%0A%20%20%20%20%5Cbegin%7Bpmatrix%7DM%28s%29%5E%5Ctop%20%26%200%20%5C%5C%200%20%26%201%5Cend%7Bpmatrix%7D. "O_{ik}(s) = \begin{pmatrix}M(s) & 0 \\ 0 & 1\end{pmatrix}
  \Omega_{(1:R, R + k), (1:R, R + k)}
    \begin{pmatrix}M(s)^\top & 0 \\ 0 & 1\end{pmatrix}.")

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Gerard16" class="csl-entry">

Berg, Gerard J. van den, and Bettina Drepper. 2016. “Inference for
Shared-Frailty Survival Models with Left-Truncated Data.” *Econometric
Reviews* 35 (6): 1075–98.
<https://doi.org/10.1080/07474938.2014.975640>.

</div>

<div id="ref-Crowther16" class="csl-entry">

Crowther, Michael J., Therese M.-L Andersson, Paul C. Lambert, Keith R.
Abrams, and Keith Humphreys. 2016. “Joint Modelling of Longitudinal and
Survival Data: Incorporating Delayed Entry and an Assessment of Model
Misspecification.” *Statistics in Medicine* 35 (7): 1193–1209.
https://doi.org/<https://doi.org/10.1002/sim.6779>.

</div>

</div>
