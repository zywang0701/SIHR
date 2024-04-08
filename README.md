
# SIHR

This package facilitates statistical inference in high-dimensional
generalized linear models (GLM) with continuous and binary outcomes. It
offers tools to construct confidence intervals and to conduct hypothesis
testing for low-dimensional objectives (e.g. $e_1^\intercal \beta$) in
both one-sample and two-sample regression scenarios.

## Background

In numerous scenarios, regression problems often involve a number of
dimensions $p$ surpassing the sample size $n$. Traditional estimators,
derived through penalized maximum likelihood methods like Lasso and
Ridge, are unsuitable for statistical inference due to significant
estimation bias introduced by the penalty term. Our package leverages
debiasing methods, addressing a broad spectrum of inference challenges
in high-dimensional GLMs, applicable to both continuous and binary
outcomes.

To demonstrate the effect of bias correction, we conducted 500
simulation rounds. The data generation process is defined with
$n=p=200$: $$
X_i \sim \mathcal{N}({\bf 0}_p, {\bf I}_p),\quad Y_i = X_i^\intercal \beta + \epsilon_i, \quad \beta = ({\bf 0.5}_{5}, \frac{1}{5}, \frac{2}{5}, \frac{3}{5}, \frac{4}{5}, \frac{5}{5}, {\bf 0}_{p-10}),\quad \epsilon_i\sim \mathcal{N}(0,1).
$$ Our objective was to conduct inference on the first coefficient
$e_1^\intercal \beta = 0.5$. However, as illustrated in the left
subfigure, Lasso estimators from the *glmnet* package exhibit
considerable bias. Conversely, our *SIHR* estimators, shown in the right
subfigure, are unbiased and thus ready for valid statistical inference.

<div class="figure" style="text-align: center">

<embed src="./images/effect-debias.pdf" title="Illustration of Bias Correction" width="80%" type="application/pdf" />
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

<!-- ## Results -->
<!-- We consider the high-dimensional GLM: for $1\leq i\leq n$, -->
<!-- \begin{equation} -->
<!--     \mathbb{E}(y_i \mid X_{i\cdot}) = f(X_{i\cdot}^\intercal \beta),\quad \textrm{with}\; -->
<!--     f(z) = \begin{cases} -->
<!--         z & \quad \textrm{for linear model;}\\ -->
<!--         \exp{(z)} / \left[1 + \exp{(z)} \right] & \quad \textrm{for logistic model;} \\ -->
<!--         %\Phi(z) & \quad \textrm{for probit model} -->
<!--     \end{cases} -->
<!--     \label{eq: glm} -->
<!-- \end{equation} -->
<!-- where $\beta \in \mathbb{R}^p$ denotes the high-dimensional regression vector. In addition to the one-sample setting, we examine the statistical inference methods for the two-sample regression models. Particularly, we generalize the regression model in \eqref{eq: glm} and consider: -->
<!-- \begin{equation} -->
<!--     \mathbb{E}(y_i^{(k)} \mid X_{i\cdot}^{(k)}) = f(X_{i\cdot}^{(k)\intercal} \beta^{(k)}) \quad \textrm{with}\; k=1,2 \; \textrm{and}\; 1\leq i\leq n_k, -->
<!--     \label{eq: two sample glm} -->
<!-- \end{equation} -->
<!-- where $f(\cdot)$ is the pre-specified link function defined as \eqref{eq: glm}. -->
<!-- This package consists of five main functions \texttt{LF}, \texttt{QF}, \texttt{CATE}, \texttt{InnProd}, and \texttt{Dist} implementing the statistical inferences for five different quantities, under the one-sample model \eqref{eq: glm} or two-sample model \eqref{eq: two sample glm}. -->
<!-- \begin{enumerate} -->
<!--     \item \texttt{LF}, abbreviated for linear functional, implements the inference approach for $\xnew^\intercal \beta$, with $\xnew \in \RR^p$ denoting a loading vector. -->
<!--     With $\xnew = e_j$ as a special case, \texttt{LF} infers the regression coefficient $\beta_j$. -->
<!--     \item \texttt{QF}, abbreviated for quadratic functional, makes inferences for $\beta^{\intercal} A \beta$. $A$ is either a pre-specified submatrix or the unknown covariance matrix $\mathbb{E}(X_{i\cdot}X_{i\cdot}^\intercal)$. -->
<!--     \item \texttt{CATE}, abbreviated for conditional average treatment effect, is to make inference for $f(\xnew^\intercal \beta^{(2)}) - f(\xnew^\intercal \beta^{(1)})$. This difference measures the discrepancy between conditional means, closely related to the conditional average treatment effect for the new observation with covariates $\xnew$. -->
<!--     \item \texttt{InnProd}, abbreviated for inner products, implements the statistical inference for $\beta^{(1)\intercal} A \beta^{(2)}$. The inner products measure the similarity between the high-dimensional vectors $\beta^{(1)}$ and $\beta^{(2)}$, which is useful in capturing the genetic relatedness in the GWAS applications. -->
<!--     \item \texttt{Dist}, short-handed for distance, makes inferences for the weighted distances  $\gamma^\intercal A \gamma$ with $\gamma = \beta^{(2)} - \beta^{(1)}$. The distance measure is useful in comparing different high-dimensional regression vectors. -->
<!-- \end{enumerate} -->
