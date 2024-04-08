
# SIHR

The package **SIHR**(**S**tatistics **I**nference for
**H**igh-dimensional **R**egression) facilitates statistical inference
in high-dimensional generalized linear models (GLM) with continuous and
binary outcomes. It offers tools to construct confidence intervals and
to conduct hypothesis testing for low-dimensional objectives
(e.g. $e_1^\intercal \beta$) in both one-sample and two-sample
regression regimes.

## Background

In numerous scenarios, regression problems often involve a number of
dimensions $p$ surpassing the sample size $n$. Traditional estimators,
derived through penalized maximum likelihood methods like Lasso and
Ridge, are unsuitable for statistical inference due to significant
estimation bias introduced by the penalty term. Our package leverages
debiasing methods, addressing a broad spectrum of inference challenges
in high-dimensional GLMs, applicable to both continuous and binary
outcomes.

To demonstrate the effect of bias correction, we conducted 250
simulation rounds. The data generation process is defined with
$n=p=200$: for $1\leq i\leq n$, 

```math
X_i \sim \mathcal{N}({\bf 0}_p, {\bf I}_p),\quad Y_i = X_i^\intercal \beta + \epsilon_i, \quad \beta = ({\bf 0.5}_{5}, \frac{1}{5}, \frac{2}{5}, \frac{3}{5}, \frac{4}{5}, \frac{5}{5}, {\bf 0}_{p-10}),\quad \epsilon_i\sim \mathcal{N}(0,1).
```

Our objective was to conduct inference on the first coefficient
$e_1^\intercal \beta = 0.5$. However, as illustrated in the left
subfigure, Lasso estimators from the **glmnet** package exhibit
considerable bias. Conversely, our **SIHR** estimators, shown in the
right subfigure, are unbiased and thus ready for valid statistical
inference.

<div class="figure" style="text-align: center" align="center">

<img src="./images/SIHR-bias.png" alt="Illustration of Bias Correction" width="80%" />
<p class="caption">
Illustration of Bias Correction
</p>

</div>

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zywang0701/SIHR")
```

## Getting Started

The package consists of 5 main functions in one or two samples regime,
with different low-dimensional objectives, respectively.

<div class="figure" style="text-align: center" align="center">

<img src="./images/Outline-SIHR.png" alt="Outline of SIHR" width="80%" />
<p class="caption">
Outline of SIHR
</p>

</div>

To illustrate the usage of these functions, we apply SIHR to simulated
data in the example.
