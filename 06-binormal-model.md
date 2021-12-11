# Binormal model {#binormal-model}










## TBA How much finished {#binormal-model-how-much-finished}
70%




## TBA Introduction {#binormal-model-introduction}
The equal variance binormal model was described in TBA Chapter 02. The ratings method of acquiring ROC data and calculation of operating points was discussed in TBA Chapter 04. It was shown, TBA Fig. \@ref(fig:ratings-paradigm-eq-var-fit-a), that for a clinical dataset the unequal-variance binormal model visually fitted the data better than the equal-variance binormal model, although how the unequal variance fit was obtained was not discussed. This chapter deals with details of the unequal-variance binormal model, often abbreviated to **binormal model**, establishes necessary notation, and derives expressions for sensitivity, specificity and the area under the predicted ROC curve). 

The binormal model describes univariate datasets, in which there is *one ROC rating per case*, as in a single observer interpreting cases, one at a time, in a single modality. By convention the qualifier "univariate" is often omitted. In TBA Chapter 21 a bivariate model will be described where each case yields two ratings, as in a single observer interpreting cases in two modalities, or the homologous problem of two observers interpreting cases in a single modality. 

The main aim of this chapter is to demystify statistical curve fitting. With the passing of Dorfman, Metz and Swensson, parametric modeling is being neglected. Researchers are instead focusing on non-parametric analysis using the empirical AUC. While useful and practical, empirical AUC yields almost no insight into what is limiting performance. Taking the mystery out of curve fitting will allow the reader to appreciate later chapters that describe more complex fitting methods, which yield important insights into factors limiting performance. 

Here is the organization of this chapter. It starts with a description of the binormal model and how it accommodates data binning. An important point, on which there is much confusion, on the invariance of the binormal model to arbitrary monotone transformations of the ratings is explicated with an example. Expressions for sensitivity and specificity are derived. Two notations used to characterize the binormal model are explained. Expressions for the pdfs of the binormal model are derived. A simple linear fitting method is illustrated: this used to be the only recourse a researcher had before Dorfman and Alf's seminal publication [@RN1081]. The maximum likelihood method for estimating parameters of the binormal model is detailed. Validation of the fitting method is described, i.e., how can one be confident that the fitting method, which makes normality and other assumptions, is valid for a dataset arising from an unknown distribution. The Appendix has a detailed derivation, originally published in a terse paper [@thompson1989statistical] on the partial-area under the ROC curve. The partial-area is defined by the area under the binormal ROC curve from $FPF = 0$ to $FPF = c$, where $0 \leq c \leq 1$. As a special case $c = 1$ yields the total area under the binormal ROC. 

## Binormal model {#binormal-model-definition}

### The basic model

The unequal-variance binormal model (henceforth abbreviated to binormal model; when I mean equal variances, it will be made explicit) is defined by (capital letters indicate random variables and their lower-case counterparts are realized values):

\begin{equation} 
Z_{k_tt}\sim N\left ( \mu_t,\sigma_{t}^{2} \right );t=1,2
(\#eq:binormal-model-z-samples-1)
\end{equation}

where 

\begin{equation} 
\left.
\begin{aligned}
\mu_1=&0\\
\mu_2=&\mu\\
\sigma_{1}^{2}=&1\\
\sigma_{2}^{2}=&\sigma^{2}
\end{aligned}
\right \}
(\#eq:binormal-model-z-samples-2)
\end{equation}


Eqn. \@ref(eq:binormal-model-z-samples-1) states that the z-samples for non-diseased cases are distributed as a $N(0,1)$  distribution, i.e., the unit normal distribution, while the z-samples for diseased cases are distributed as a  $N(\mu,\sigma^2)$ distribution, i.e., a normal distribution with mean $\mu$  and variance $\sigma^2$. *This is a 2-parameter model of the z-samples, not counting additional threshold parameters needed for data binning.* ^[A more complicated version of this model allows the mean of the non-diseased distribution to be non-zero and its variance different from unity. The 4-parameter model is no more general than the 2-parameter model. The reason is that one is free to transform the decision variable, and associated thresholds, by applying arbitrary monotonic increasing function transformation, which do not change the ordering of the ratings and hence do not change the ROC curve. So if the mean of the noise distribution were non-zero, subtracting this value from all Z-samples would shift the effective mean of the non-diseased distribution to zero (the shifted Z-values are monotonically related to the original values) and the mean of the shifted diseased distribution becomes $\mu_2-\mu_1$. Next, one scales or divides (division by a positive number is also a monotonic transformation) all the Z-samples by $\sigma_1$, resulting in the scaled non-diseased distribution having unit variance, and the scaled diseased distribution has mean $\frac{\mu_2-\mu_1}{\sigma_1}$  and variance $(\frac{\sigma_2}{\sigma_1})^2$. Therefore, if one starts with 4 parameters then one can, by simple shifting and scaling operations, reduce the model to 2 parameters, as in Eqn. \@ref(eq:binormal-model-z-samples-1). [The author has seen a publication on Bayesian ROC estimation using the four-parameter model.]]



### Additional parameters for binned data
In an R-rating ROC study the observed ratings $r$ take on integer values, 1 through $R$, it being understood that higher ratings correspond to greater confidence for disease. Defining dummy cutoffs $\zeta_0 = -\infty$ and  $\zeta_R = +\infty$, the binning rule for a case with realized z-sample z is (Chapter \@ref(ratings-paradigm), Eqn. \@ref(eq:ratings-paradigm-binningRule)):

\begin{equation} 
\text{if} \left (\zeta_{r-1} \le z \le \zeta_r  \right )\Rightarrow \text {rating} = r
(\#eq:binormal-modelZBinning)
\end{equation}


<img src="06-binormal-model_files/figure-html/unnamed-chunk-3-1.png" width="672" />

In the unequal-variance binormal model, the variance $\sigma^2$ of the z-samples for diseased cases is allowed to be different from unity. Most ROC datasets are consistent with  $\sigma > 1$. The above figure, generated with  $\mu = 1.5, \sigma = 1.5, \zeta_1 = -2, \zeta_2 = -0.5, \zeta_3 = 1, \zeta_4 = 2.5$, illustrates how realized z-samples are converted to ratings, i.e., application of the binning rule \@ref(eq:binormal-modelZBinning). For example, a case with  z-sample equal to -2.5 would be rated "1", and one with  z-sample equal to -1 would be rated "2", cases with z-samples greater than 2.5 would be rated "5", etc.

### Sensitivity and specificity
Let $Z_t$ denote the random z-sample for truth state $t$ ($t$ = 1 for non-diseased and $t$ = 2 for diseased cases).  Since the distribution of z-samples from disease-free cases is $N(0,1)$, the expression for specificity, Chapter "Modeling Binary Paradigm", Eqn. 3.13, applies. It is reproduced below: 

\begin{equation} 
Sp\left ( \zeta \right )=P\left ( Z_1 < \zeta \right )=\Phi\left ( \zeta \right )
(\#eq:binormal-modelSp)
\end{equation}

To obtain an expression for sensitivity, consider that for truth state $t = 2$, the random variable $\frac{Z_2-\mu}{\sigma}$  is distributed as $N(0,1)$: 

\begin{equation*} 
\frac{Z_2-\mu}{\sigma}\sim N\left ( 0,1 \right )
\end{equation*}

Sensitivity is $P\left ( Z_2 > \zeta \right )$, which implies, because $\sigma$ is positive (subtract $\mu$ from both sides of the “greater than” symbol and divide by $\sigma$):

\begin{equation} 
Se\left ( \zeta | \mu, \sigma \right )= P\left ( Z_2 > \zeta \right )=P\left ( \frac{Z_2-\mu}{\sigma} > \frac{\zeta-\mu}{\sigma} \right )
(\#eq:binormal-modelSe)
\end{equation}

The right-hand-side can be rewritten as follows:

\begin{equation*} 
Se\left ( \zeta | \mu, \sigma \right )= 1 - P\left ( \frac{Z_2-\mu}{\sigma} \leq  \frac{\zeta-\mu}{\sigma} \right )\\
=1-\Phi\left (  \frac{\zeta-\mu}{\sigma}\right )=\Phi\left (  \frac{\mu-\zeta}{\sigma}\right )
\end{equation*}

Summarizing, the formulae for the specificity and sensitivity for the binormal model are: 

\begin{equation} 
Sp\left ( \zeta \right ) = \Phi\left ( \zeta \right )\\
Se\left ( \zeta | \mu, \sigma \right ) = \Phi\left (  \frac{\mu-\zeta}{\sigma}\right )
(\#eq:binormal-modelSeSp)
\end{equation}

The coordinates of the operating point defined by $\zeta$ are given by:

\begin{equation} 
FPF\left ( \zeta \right ) = 1 - Sp\left ( \zeta \right ) = 1 - \Phi\left ( \zeta \right ) = \Phi\left ( -\zeta \right )
(\#eq:binormal-model-fpf)
\end{equation}

\begin{equation} 
TPF\left ( \zeta | \mu, \sigma \right ) = \Phi\left ( \frac{\mu-\zeta}{\sigma} \right )
(\#eq:binormal-model-tpf)
\end{equation}

These expressions allow calculation of the operating point for any $\zeta$. An equation for a curve is usually expressed as $y=f(x)$. An expression of this form for the ROC curve, i.e., the y-coordinate (TPF) expressed as a function of the x-coordinate (FPF), follows upon inversion of the expression for FPF, Eqn.  \@ref(eq:binormal-model-fpf):

\begin{equation} 
\zeta = -\Phi^{-1}\left ( FPF \right )
(\#eq:binormal-modelZeta)
\end{equation}

Substitution of Eqn. \@ref(eq:binormal-modelZeta) in Eqn. \@ref(eq:binormal-model-tpf) yields:

\begin{equation} 
TPF = \Phi\left ( \frac{\mu + \Phi^{-1}\left (FPF  \right )}{\sigma} \right )
(\#eq:binormal-model-roc-curve1)
\end{equation}

This equation gives the dependence of TPF on FPF, i.e., the equation for the ROC curve. It will be put into standard notation next.  


### Binormal model in conventional notation
The following notation is widely used in the literature: 

\begin{equation} 
a=\frac{\mu}{\sigma};b=\frac{1}{\sigma}
(\#eq:binormal-modelabParameters)
\end{equation}

The reason for the $(a,b)$  instead of the  $(\mu,\sigma)$ notation is that Dorfman and Alf assumed, in their seminal paper [@RN1081], that the diseased distribution (signal distribution in signal detection theory) had unit variance, and the non-diseased distribution (noise) had standard deviation $b$ ($b > 0$) or variance $b^2$, and that the separation of the two distributions was $a$, see figure below. In this example: $a = 1.11$ and $b = 0.556$, corresponding to $\mu = 2$ and $\sigma = 1.8$. Dorfman and Alf's fundamental contribution, namely estimating these parameters from ratings data, to be described below, led to the widespread usage of the  $(a,b)$ parameters estimated by their software (RSCORE), and its newer variants (e.g., RSCORE–II, ROCFIT and ROCKIT). 

By dividing the z-samples by $b$, the variance of the distribution labeled "Noise" becomes unity, its mean stays at zero, and the variance of the distribution labeled "Signal" becomes $1/b$, and its mean becomes $a/b$, as shown below. It illustrates that the inverses of Eqn. \@ref(eq:binormal-modelabParameters) are:

\begin{equation} 
\mu=\frac{a}{b};\sigma=\frac{1}{b}
(\#eq:binormal-modelabParametersInv)
\end{equation}

Eqns. \@ref(eq:binormal-modelabParameters) and \@ref(eq:binormal-modelabParametersInv) allow conversion from one notation to another.




```r
grid.arrange(p1,p2,ncol=2)
```

<div class="figure">
<img src="06-binormal-model_files/figure-html/binormal-model-ab2-mu-sigma-1.png" alt="Plot A shows the definitions of the (a,b) parameters of the binormal model. In plot B the x-axis has been rescaled so that the noise distribution has unit variance, thereby illustrations between (a,b) and the ($\mu,\sigma$) parameters." width="672" />
<p class="caption">(\#fig:binormal-model-ab2-mu-sigma)Plot A shows the definitions of the (a,b) parameters of the binormal model. In plot B the x-axis has been rescaled so that the noise distribution has unit variance, thereby illustrations between (a,b) and the ($\mu,\sigma$) parameters.</p>
</div>


## Binormal ROC curve {#binormal-model-roc-curve}
Using the $(a,b)$ notation, Eqn. \@ref(eq:binormal-model-roc-curve1) for the ROC curve reduces to:

\begin{equation} 
TPF = \Phi\left ( a+ b \Phi^{-1}\left (FPF  \right ) \right )
(\#eq:binormal-model-roc-curve)
\end{equation}

Since $\Phi^{-1}(FPF)$  is an increasing function of its argument $FPF$, and $b > 0$, the argument of the  $\Phi$ function is an increasing function of $FPF$. Since $\Phi$  is a monotonically increasing function of its argument, $TPF$ is a monotonically increasing function of $FPF$. This is true regardless of the sign of $a$. If $FPF = 0$, then $\Phi^{-1}(0) = -\infty$  and $TPF = 0$. If $FPF = 1$, then $\Phi^{-1}(1) = +\infty$ and $TPF = 1$. Regardless of the value of $a$, as long as $b \ge 0$, the ROC curve starts at (0,0) and increases monotonically ending at (1,1).

From Eqn. \@ref(eq:binormal-model-fpf) and Eqn. \@ref(eq:binormal-model-tpf), the expressions for $FPF$ and $TPF$ in terms of model parameters $(a,b)$ are:

\begin{equation}
\left.
\begin{aligned}
FPF\left ( \zeta \right ) &= \Phi\left ( -\zeta \right )\\
TPF &= \Phi\left ( a - b \zeta \right )
\end{aligned}
\right \}
(\#eq:binormal-model-op-point-ab)
\end{equation}

## Scalar threshold-independent measure {#binormal-model-scalar-measure}

Sensitivity-specificity is a dual (two-valued) measure of performance. Using a dual measure it is difficult to unambiguously compare two systems since one cannot separate the effect of reporting threshold from the measures. For example, if sensitivity is higher for one system but specificity is higher for another, this could be due to different thresholds. Sensitivity and specificity depend on the threshold. As the threshold changes, sensitivitya nd specificity are both affected in opposite directions. Desirable is a scalard measure of performance that takes this variation into account and does not depend on any specific threshold. 

Generally accepted measures are the partial-area $A_{z;c}$ under the ROC, Eqn. \@ref(eq:binormal-model-partial-area-final), the full-area $A_z$ under the ROC, Eqn. \@ref(eq:binormal-model-ab-2az), and the $d'$ index Eqn. \@ref(eq:binormal-model-ab-2dprime). 

Before deriving analytical expressions for these measures let us further examine the premise that sensitivity-specificity is undesirable because it is a 2D measure. A trivial way to convert it to a scalar measure is to sum the two values: high sensitivity and high specificity are both desirable, so a high value of their sum is certainly also desirable. In fact this is the basis for the Youden index, defined as sensitivity plust specificity minus one [@youden1950index]. (Subtracting one makes the Youden index range from 0 to 1.) However, this index varies with the position of the operating point on the ROC curve. (The operating point at which it is maximum is often thought of as the optimal operating point on the ROC curve.)

To emphasize, we desire a scalar measure that is threshold independent.



### Partial AUC {#binormal-model-partial-auc}

While this is a scalar measure, it does depend on choice of operating point. It is included here as it yields, as a special case, a scalar measure that does not depend on choice of operating point. The details are in Section \@ref(binormal-model-appendix-1-partial-auc), which derives the formula for the partial-area under the unequal-variance binormal model. The final result is:

\begin{equation}
A_{z;c} = \int_{z_2=-\infty}^{\Phi^{-1}\left ( c \right )}   \int_{z_1=-\infty}^{\frac{a}{\sqrt{1+b^2}}} \phi\left ( z_1,z_2;\rho \right ) dz_1dz_2
(\#eq:binormal-model-partial-area-final)
\end{equation}

The threshold $\zeta_1$ corresponding to $FPF = c$ is given by:

\begin{equation} 
\zeta_1 = - \Phi^{-1}   \left ( c \right )
(\#eq:binormal-model-zeta-c-relation)
\end{equation}

$A_{z;c}$ is the area under the partial ROC curve extending from $FPF = 0$ to $FPF = c$ and $\phi\left ( z_1,z_2;\rho \right )$ is the standard bivariate normal distribution, where the correlation coefficient $\rho$ of the distribution is defined by:


\begin{equation}
\rho = - \frac{b}{\sqrt{1+b^2}}
(\#eq:binormal-model-rho-final)
\end{equation}

The bivariate 2D integral can be evaluated numerically. The following code illustrates calculation of the partial-area measure using the function `pmvnorm` in `R` package `mvtnorm`. The following parameter values were used: $a = 2$, $b = 1$ and $\zeta_1 = 1.5$. (The parameter $b$ was deliberately chosen equal to one so that we do not have to worry about improper ROC curves.)




```{.r .numberLines}
a <- 2;b <- 1; zeta1 <- 1.5
A_z <- pnorm(a/sqrt(1+b^2))
opPtx <- pnorm(-zeta1)
opPty <- pnorm(a - b * zeta1)
rho <- -b/sqrt(1+b^2)
Lower1 <- -Inf
Upper1 <- qnorm(opPtx)
Lower2 <- -Inf
Upper2 <- a/sqrt(1+b^2)
sigma <- rbind(c(1, rho), c(rho, 1))
A_zc <- as.numeric(pmvnorm(
  c(Lower1, Lower2), 
  c(Upper1, Upper2), 
  sigma = sigma))
```


The partial-area measure is $A_{z;c}$ = 0.0352195. The corresponding full-area measure is $A_z$ = 0.9213504. $A_{z;c}$ is small because the reporting threshold is high. However, $A_{z;c}$ should not be confused with true performance of the observer, as shown in Section \@ref(binormal-model-partial-true).



### Full AUC {#binormal-model-full-auc}

A special case of this formula is the area under the full ROC curve, shown below using both parameterizations of the binormal model:

\begin{equation} 
A_z=\Phi\left ( \frac{a}{\sqrt{1+b^2}} \right )=\Phi\left ( \frac{\mu}{\sqrt{1+\sigma^2}} \right )
(\#eq:binormal-model-ab-2az)
\end{equation} 

The binormal fitted AUC increases as $a$ increases or as $b$ decreases. Equivalently, it increases as $\mu$ increases or as $\sigma$ decreases. In the example just given, the full AUC is $A_z$ = 0.9213504.


### The d' measure {#binormal-model-d-prime}

The $d'$ parameter is defined as the separation of two unit-variance normal distributions yielding the same AUC as that predicted by the $(a,b)$ parameter binormal model. It is defined by:

\begin{equation} 
d'=\sqrt{2}\Phi^{-1}\left ( A_z \right )
(\#eq:binormal-model-ab-2dprime)
\end{equation} 




The d' index corresponding to the above binormal parameters is 2. The transformation from an index that ranges from 0.5 to 1 to one that ranges from 0 to infinity can be viewed as desirable. The d' index can be regarded as a perceptual signal-to-noise-ratio. 




## Partial AUC vs. true performance {#binormal-model-partial-true}

* A *partial-area observer* such as in Section \@ref(binormal-model-partial-auc) rates cases as follows: for the sub-set of cases defined by $z \ge \zeta_1$ the observer reports *explicit* ratings exactly equal to the observed z-samples (or some monotonic transformation of the z-samples). For the remaining cases the observer assigns a *fixed value rating that is smaller than $\zeta_1$* (the exact value does not matter; these cases are said to be assigned *implicit* ratings). 

* In contrast, the *full-area observer* reports explcit ratings *for all cases*. 

*To measure true performance of the partial-area observer one must, of course, include all cases.* The ROC curve extends continously from the origin to the solid dot *plus the area under the dotted line* extending from the solid dot to (1,1). True performance, the area under the continuous section plus that under the straight line extension, is denoted $A_{z;c,TRUE}$ and is defined by:


\begin{equation} 
A_{z;c,\text{TRUE}} = A_{z;c} + \frac {\left ( 1 - FPF \right ) \left ( 1 + TPF \right )}{2}
(\#eq:binormal-model-partial-area-true-performance)
\end{equation}

In other words one adds to $A_{z;c}$ the area of the trapezoid with bases each equal to $(1 - FPF)$ and opposing sides equal to $TPF$ and unity. 

Since the partial-area observer does not preserve ordering information, *true performance of a partial-area observer is smaller than performance $A_z$ of a full-area observer*. 

\begin{equation} 
A_{z;c,\text{TRUE}} \le A_{z}
(\#eq:binormal-model-true-performance-az-inequality)
\end{equation}



True performance is illustrated with the following simulation 2AFC study. The `Wilcoxon` function, defined next, can be thought of as the mathematical equivalent of a 2AFC study, conducted with all possible pairings of non-diseased and diseased cases. For each pairing, if the z-sample of the diseased case exceeds that of the non-diseased case one adds unity to a zero-initialized counter; if it is smaller one does nothing; if they are equal one adds 0.5; and finally one divides by the number of comparisons.



```{.r .numberLines}
Wilcoxon <- function (zk1, zk2)
{
  K1 = length(zk1)
  K2 = length(zk2)
  W <- 0
    for (k1 in 1:K1) {
      W <- W + sum(zk1[k1] < zk2)
      W <- W + 0.5 * sum(zk1[k1] == zk2)
    }
    W <- W/K1/K2
  return (W)
}
```



The following code saves 10,000 pairs of ratings in two arrays: `z[1,]` and `z[2,]`. The first array corresponds to non-diseased cases and the second to diseased cases. Note the usage, at lines 3-4, of the $a,b$ values to define the two distributions. The array `zc`, initially a copy of `z`, is selectively binned by setting, lines 6-7, all ratings less than $\zeta_1$ to -100. The ordering information for these z-samples is lost. 




```{.r .numberLines}
nPairs <- 10000
z <- array(dim = c(2, nPairs))
z[1,] <- rnorm(nPairs, sd = b)
z[2,] <- rnorm(nPairs, mean = a, sd = 1)
zc <- z
zc[1,z[1,] < zeta1] <- -100 # ratings of partial area observer
zc[2,z[2,] < zeta1] <- -100 # do:
```



The following code prints the predicted and observed full areas under the ROCs followed by the predicted and observed true performances. With this many cases sampling variability is small and the predicted and observed values are close. 




```
#> A_z predicted =  0.9213504 
#> A_z observed =  0.9193161
#> A_z{c;true} predicted =  0.8244498 
#> A_z{c;true} observed =  0.825695
```


Note that:

* $A_{z;c,\text{TRUE}} < A_z$, because ordering information is lost for all cases with z-samples less than $\zeta_1$. 
* $A_{z;c,\text{TRUE}} >> A_{z;c}$, because of the large contribution from the area under the straight line, left poanel Fig. \@ref(fig:binormal-model-threshold-dependence-2). 


## Illustrative plots {#binormal-model-illustrative-plots}

In the ROC plots below the partial-area observer curve is shown as a continuous line extending from the origin to the limiting point *plus* a dotted line extending from the limiting point to (1,1). The continuous section is determined by cumulating cases with z-samples $z \ge \zeta_1$ while the (1,1) point is determined by cumulating all cases. 


The ROC curve for both types of observers is shown in the left panel of \@ref(fig:binormal-model-threshold-dependence-2) for the following parameters: $a = 2$, $b = 1$ and $\zeta_1 = 1.5$; $\zeta_1$ corresponds to $c \equiv FPF =  \Phi(-\zeta_1)$ = 0.0668072 and $TPF =  \Phi(a - b\zeta_1)$ = 0.6914625. In other words the limiting point coordinates are (0.067, 0.691), shown in the plot by the solid dot. Partial AUC $A_{z;c}$ equals 0.0352195. The full-area ROC curve, shown by the complete solid curve, extends from (0,0) to (1,1), the area under which is $A_z$ = 0.9213504. 















<div class="figure">
<img src="06-binormal-model_files/figure-html/binormal-model-threshold-dependence-2-1.png" alt="Left panel: binormal ROC curve corresponding to a = 2 and b = 1. The dot is the operating point corresponding to $\zeta_1 = 1.5$. The continuous curve extending from the origin to (1,1) represents the full ROC. Note that in the region above the dot the continuus curve is above the dotted line, meaning true performance of an observer who only rates a sub-set of cases is less than performance of an observer who rates all cases. Right panel: variation of true performance with FPF; at FPF = 0 the plot starts at ordinate equal to 0.5 and levels out at FPF = 1 at AUC = $A_z = 0.921.$" width="672" />
<p class="caption">(\#fig:binormal-model-threshold-dependence-2)Left panel: binormal ROC curve corresponding to a = 2 and b = 1. The dot is the operating point corresponding to $\zeta_1 = 1.5$. The continuous curve extending from the origin to (1,1) represents the full ROC. Note that in the region above the dot the continuus curve is above the dotted line, meaning true performance of an observer who only rates a sub-set of cases is less than performance of an observer who rates all cases. Right panel: variation of true performance with FPF; at FPF = 0 the plot starts at ordinate equal to 0.5 and levels out at FPF = 1 at AUC = $A_z = 0.921.$</p>
</div>


As FPF increases true-performance increases. the right panel of Fig. \@ref(fig:binormal-model-threshold-dependence-2) shows the variation of true performance $A_{z;c,\text{TRUE}}$ with FPF. The curve starts from (0, 0.5) and ends at (1.000, 0.921). For low values of FPF the curve is very steep while for FPF > 0.25 the curve levels out, approaching the maximum value defined by $A_z$ = 0.9213504. True performance is maximized at $\zeta_1 = -\infty$.



Fig. \@ref(fig:binormal-model-threshold-dependence-3), left panel, corresponding to $a = 1$, $b = 0.2$ and $\zeta_1 = 1.5$, shows an improper ROC curve. The dashed line is well above the continuous curve and true performance is maximized at a finite value of $\zeta_1$, corresponding to $FPF =  0.153$, see right panel. This is an invalid conclusion since an improper ROC curve is a fitting artifact of the binormal model easily avoided by using modern curve-fitting methods (eg., PROPROC, CBM or RSM). TBA However, since the wAFROC has an operating characteristic with an improper-like feature but which is not a fitting artifact, this example serves a purpose, elaborated on in TBA Chapter (optim-op-point), where it is shown that by maximizing the area under the wAFROC one can find the optimal threshold of an algorithmic observer.








```
#> true performance max occurs at FPF =  0.1525
```


<div class="figure">
<img src="06-binormal-model_files/figure-html/binormal-model-threshold-dependence-3-1.png" alt="The left panel shows the visibly improper ROC curve for a = 1 and b = 0.2. The solid line is below the dotted line. The right panel shows the variation of true performance pAUC_TRUE with FPF. True performance is maximized at FPF = 0.153. Since improper ROC fits are fitting artifacts, this example does not negate the previous finding that true performance for a proper ROC curve is maximized by setting the threshold to report all cases, i.e., FPF = 1." width="672" />
<p class="caption">(\#fig:binormal-model-threshold-dependence-3)The left panel shows the visibly improper ROC curve for a = 1 and b = 0.2. The solid line is below the dotted line. The right panel shows the variation of true performance pAUC_TRUE with FPF. True performance is maximized at FPF = 0.153. Since improper ROC fits are fitting artifacts, this example does not negate the previous finding that true performance for a proper ROC curve is maximized by setting the threshold to report all cases, i.e., FPF = 1.</p>
</div>


## Geometrical argument {#binormal-model-geometrical-argument}

Defining geometrical features of a proper ROC are: 

* As one moves up the curve the slope decreases monotonically; 
* At each point the slope is greater than that of the straight line connecting the point to (1,1); 
* The curve ends at (1,1). 

The geometry ensures that true performance for a proper ROC is maximized at $\zeta_1 = -\infty$, i.e., at FPF = 1, as in Fig. \@ref(fig:binormal-model-threshold-dependence-2), right panel.


## Optimal operating point on ROC {#binormal-model-optimal-op-pt} 
We have seen that optimal ROC AUC is achieved by setting $\zeta_1 = -\infty$, i.e., by reporting all cases as diseased. Of course, from clinical considerations, this is nonsense. Consider screening mammography, where typically for every 1000 cases only 5 are malignant. Recalling everybody would incur huge costs from having to rule out cancer in 995 actually non-diseased patients. Of course the 5 malignant cancers would be confirmed at the follow-up diagnostic mammography examination. But one can clearly see that the benefit of correctly detecting the 5 malignancies is far outweighed by the 995 unnecessary recalls. And if one is going to recall everybody, why perform the initial screeing mammography exam?

So what is going on? The problem is that AUC measures classification performance in a 2AFC task. A screening examination is not a 2AFC task: the radiologist is not presented two cases, one non-diseased and one diseased, and asked to pick the diseased patient. Rather, the radiologist is shown images of a single patient, and the object is to maximize the detection rate while minimizing false positives. 

To address this optimization task one needs to know the costs and benefits of the four decision outcomes in the binary paradigm: true and false positives, and true and false negatives. This has been addressed in [@metz1978rocmethodology]. Here is the reasonaing. Let 

* $C_0$ denote the overhead cost of performing the imaging examination, 
* $C_{\text{TP}}$ denote the cost of a true positive decision (a benefit can be expressed as a negative cost), 
* $C_{\text{FN}}$ denote the cost of a false negative decision,
* $C_{\text{FP}}$ denote the cost of a false positive decision, and 
* $C_{\text{TN}}$ denote the cost (or negative benefit) of a true negative decision.    


It is shown [@metz1978rocmethodology] that the average cost of the examination is:


\begin{equation} 
\overline{C} = C_0 + C_{\text{TP}} P(\text{TP}) 
+ C_{\text{TN}} P(\text{TN})
+ C_{\text{FP}} P(\text{FP})
+ C_{\text{FN}} P(\text{FN})
(\#eq:binormal-model-cost)
\end{equation}

In this equation $P(\text{TP})$ is the probability of a TP-event, etc. These probabilities are related to disease prevalence $P(+)$ and the operating point by:

\begin{equation} 
\left.
\begin{aligned}
P(\text{TP}) =& P(+) \text{TPF} \\
P(\text{TN}) =& (1-P(+)) (1-\text{FPF}) \\
P(\text{FP}) =& (1-P(+)) \text{FPF} \\
P(\text{FN}) =& P(+) (1-\text{TPF}) 
\end{aligned}
\right \}
(\#eq:binormal-model-cost2)
\end{equation}


With these substitutions one gets for the average cost:


\begin{equation} 
\overline{C} = C_0 + C_{\text{TN}} P(-) + C_{\text{FN}}  P(+) \\
+ (C_{\text{TP}} - C_{\text{FN}}) P(+) \text{TPF} \\
+ (C_{\text{FP}} - C_{\text{TN}}) P(-) \text{FPF} \\
(\#eq:binormal-model-cost3)
\end{equation}

Equating the derivative of the average cost to zero, to minimize the average cost, one gets:


\begin{equation} 
\frac{d(\text{TPF})}{d(\text{FPF})} = \frac{C_{\text{FP}} - C_{\text{TN}}}{C_{\text{FN}} - C_{\text{TP}}}\frac{P(-)}{P(+)}
(\#eq:binormal-model-cost4)
\end{equation}

This defines the slope $\frac{d(\text{TPF})}{d(\text{FPF})}$ of the ROC at the optimal operating point, i.e., the point that minimizes the average cost of the examination. Note that $P(-) = 1 - P(+)$. 

* If disease prevalence is high, then the optimal operating point is where the slope of the ROC is low, which is near the  upper-right corner. With mostly diseased cases it makes sense to set the operating point at high sensitivity and low specificity. Conversely, with low prevalence, one should set the operating point at low sensitivity and high specificity.

* For a given disease prevalence, if the cost of a FP decision is high (or if the benefit of a TN is high - recall that a benefit is the same as a negative cost), then the optimal operating point is where the slope of the ROC is high, which is near the lower-left corner. One sets the operating point at low sensitivity and high specificity.

* For a given disease prevalence, if the cost of a FN decision is high (or if the benefit of a TP is high), then the optimal operating point is where the slope of the ROC is low, which is near the  upper-right corner. One sets the operating point at high sensitivity and low specificity.


The costs and benefits are often difficult to quantify. If one assumes that the right hand side of Eqn. \@ref(eq:binormal-model-cost4) equals unity (e.g., the four costs / benefits are equal and disease-prevalence is 50%) then the optimal operating point is defined by that point on the ROC curve where the slope is unity, which is the point of nearest approach of the curve to the upper-left corner. This corresponds to maximizing the Youden index [@youden1950index], defined as the sum of sensitivity and specificity minus one. This is demonstrated in the following code.


```r
a <- 2;b <- 1
z <- seq(-3,5.5,0.05)
FPF <- pnorm(-z)
TPF <- pnorm(a - b*z)
Youden <- TPF + (1 - FPF) - 1
curve <- data.frame(FPF = FPF, TPF = TPF, YOU = Youden)
dist <- sqrt(FPF^2 + (1 - TPF)^2)
p1 <- ggplot2::ggplot(curve, aes(x = FPF, y = TPF)) + 
  geom_line() +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1))
p2 <- ggplot2::ggplot(curve, aes(x = FPF, y = YOU)) + 
  geom_line() +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1))
indxDist <- which(dist == min(dist))
indxYoud <- which(Youden == max(Youden))
if (indxDist != indxYoud) stop("The two indices are different") else {
  cat("Op Pt corresponding to max Youden and min distance is: \nFPF = ", 
      FPF[indxDist], 
      "\nTPF = ", 
      TPF[indxDist])
}
#> Op Pt corresponding to max Youden and min distance is: 
#> FPF =  0.1586553 
#> TPF =  0.8413447
```


<div class="figure">
<img src="06-binormal-model_files/figure-html/binormal-model-youden-max-1.png" alt="Left panel: binormal ROC curve corresponding to a = 2 and b = 1. Right panel: variation of Youden index with FPF; the plot shows a maximum at FPF = 0.1586553; this corresponds to the nearest approch of the ROC curve to the upper-left corner." width="672" />
<p class="caption">(\#fig:binormal-model-youden-max)Left panel: binormal ROC curve corresponding to a = 2 and b = 1. Right panel: variation of Youden index with FPF; the plot shows a maximum at FPF = 0.1586553; this corresponds to the nearest approch of the ROC curve to the upper-left corner.</p>
</div>


## Discussion{#binormal-model-discussion}
The binormal model is historically very important and the contribution by Dorfman and Alf [@RN1081] was seminal. Prior to their work, there was no valid way of estimating AUC from observed ratings counts. Their work and a key paper [@RN1487] accelerated research using ROC methods. The number of publications using their algorithm, and the more modern versions developed by Metz and colleagues, is probably well in excess of 500. Because of its key role, I have endeavored to take out some of the mystery about how the binormal model parameters are estimated. In particular, a common misunderstanding that the binormal model assumptions are violated by real datasets, when in fact it is quite robust to apparent deviations from normality, is addressed. 

A good understanding of this chapter should enable the reader to better understand alternative ROC models, discussed later.

It has been stated that the `b`-parameter of the binormal model is generally observed to be less than one, consistent with the diseased distribution being wider than the non-diseased one. The ROC literature is largely silent on the reason for this finding. One reason, namely location uncertainty, is presented in Chapter "Predictions of the RSM", where RSM stands for Radiological Search Model. Basically, if the location of the lesion is unknown, then z-samples from diseased cases can be of two types, samples from the correct lesion location, or samples from other non-lesion locations. The resulting mixture distribution will then appear to have larger variance than the corresponding samples from non-diseased cases. This type of mixing need not be restricted to location uncertainty. Even is location is known, if the lesions are non-homogenous (e.g., they contain a range of contrasts) then a similar mixture-distribution induced broadening is expected. The contaminated binormal model (CBM) - see Chapter TBA - also predicts that the diseased distribution is wider than the non-diseased one.

The fact that the `b`-parameter is less than unity implies that the predicted ROC curve is improper, meaning its slope is not monotone decreasing as the operating point moves up the curve. The result is that a portion of the curve, near (1,1) that crosses the chance-diagonal and hooks upward approaching (1,1) with infinite slope. Ways of fitting proper ROC curves are described in Chapter "Other proper ROC models". Usually the hook is not readily visible, which has been used as an excuse to ignore the problem. For example, in Fig. 6.4, one would have to "zoom-in" on the upper right corner to see it, but the reader should make no mistake about it, the hook is there as  . 

A recent example is Fig. 1 in the publication resulting from the Digital Mammographic Imaging Screening Trial (DMIST) clinical trial [@RN1784] involving 49,528 asymptomatic women from 33 clinical sites and involving 153 radiologists, where each of the film modality ROC plots crosses the chance diagonal and hooks upwards to (1,1), which as is known, results anytime $b <1$.

The unphysical nature of the hook (predicting worse than chance-level performance for supposedly expert readers) is not the only reason for seeking alternate ROC models. The binormal model is susceptible to degeneracy problems. If the dataset does not provide any interior operating points (i.e., all observed points lie on the axes defined by FPF = 0 or TPF = 1) then the model fits these points with b = 0. The resulting straight-line segment fits do not make physical sense. These problems are addressed by the contaminated binormal model16 to be discussed in Chapter "Other proper ROC models". The first paper in the series has particularly readable accounts of data degeneracy.

To this day the binormal model is widely used to fit ROC datasets. In spite of its limitations, the binormal model has been very useful in bringing a level of quantification to this field that did not exist prior to [@RN1081].

## Appendix I: Density functions {#binormal-model-pdf-curves-appendix-1}
According to Eqn. \@ref(eq:binormal-model-z-samples-1) the probability that a z-sample is smaller than a specified threshold  $\zeta$, i.e., the CDF function, is:

\begin{equation*} 
P\left ( Z \le \zeta \mid  Z\sim N\left ( 0,1 \right ) \right ) = 1-FPF\left ( \zeta \right ) = \Phi \left ( \zeta  \right )
\end{equation*}

\begin{equation*} 
P\left ( Z \le \zeta \mid  Z\sim N\left ( \mu,\sigma^2 \right ) \right ) = 1-TPF\left ( \zeta \right ) = \Phi \left ( \frac{\zeta - \mu}{\sigma}  \right )
\end{equation*}

Since the *pdf* is the derivative of the corresponding CDF function, it follows that (the subscripts N and D denote non-diseased and diseased cases, respectively):

\begin{equation*} 
pdf_N\left ( \zeta \right ) = \frac{\partial \Phi\left ( \zeta \right )}{\partial \zeta} = \phi\left ( \zeta \right ) \equiv \frac{1}{\sqrt{2 \pi}}\exp\left ( -\frac{\zeta^2}{2} \right )
\end{equation*}

\begin{equation*} 
pdf_D\left ( \zeta \right ) = \frac{\partial \Phi\left ( \frac{\zeta - \mu}{\sigma} \right )}{\partial \zeta} = \frac{1}{\sigma} \phi\left ( \frac{\zeta - \mu}{\sigma} \right ) \equiv \frac{1}{\sqrt{2 \pi}\sigma}\exp\left ( -\frac{\left (\zeta-\mu  \right )^2}{2\sigma} \right )
\end{equation*}

The second equation can be written in $(a,b)$  notation as:

\begin{equation*} 
pdf_D\left ( \zeta \right ) = b\phi\left ( b\zeta-a \right ) = \frac{b}{\sqrt{2 \pi}}\exp\left ( -\frac{\left (b\zeta - a \right )^2}{2} \right )
\end{equation*}


## Appendix II: Area under binormal ROC {#binormal-model-appendix-2}
### General case (partial-area) {#binormal-model-appendix-1-partial-auc}
This section is based on [@thompson1989statistical]. In what follows, FPF is abbreviates to x and TPF to y. Then the equation for the ROC curve is \@ref(eq:binormal-model-roc-curve):

\begin{equation} 
y=\Phi\left ( a + b \Phi^{-1} \left ( x \right ) \right )
(\#eq:binormal-model-xy)
\end{equation}

The partial-area under the ROC curve from $x = 0$ to $x = c$, where $0 ≤ c ≤ 1$, is given by: 

\begin{equation} 
A_{z;c} = \int_{0}^{c} y dx =  \int_{0}^{c} dx \Phi\left ( a + b \Phi^{-1} \left ( x \right ) \right ) 
(\#eq:binormal-model-partial-area-1)
\end{equation}

Define change of variable:

\begin{equation} 
x = \Phi\left ( x_1 \right )
(\#eq:binormal-model-change-variable-1)
\end{equation}

which implies:

\begin{equation}
\left. 
\begin{aligned}
x_1 = &  \Phi^{-1}\left ( x \right ) \\
dx = &  dx_1 \phi\left ( x_1 \right )
\end{aligned}
\right \}
(\#eq:binormal-model-dx-dx1)
\end{equation}

This yields:

\begin{equation}
\left. 
\begin{aligned}
A_{z;c} = &  \int_{0}^{c} dx \Phi \left ( a + b x_1 \right )\\
 = &  \int_{-\infty}^{\Phi^{-1}\left ( c \right )} dx_1 \phi\left ( x_1 \right )\Phi \left ( a + b x_1 \right )
\end{aligned}
\right \}
(\#eq:binormal-model-partial-area-2)
\end{equation}



The right hand side of Eqn. \@ref(eq:binormal-model-partial-area-2) can be expressed as an integral over the bivariate normal distribution as follows. From the definition of the $\Phi$ function the above integral can be written as the following double integral:


\begin{equation}
A_{z;c} = \int_{x_1=-\infty}^{\Phi^{-1}\left ( c \right )} dx_1 \phi\left ( x_1 \right ) \int_{x_2=-\infty}^{a+bx_1} \phi\left ( x_2 \right ) dx_2
(\#eq:binormal-model-partial-area-double-integral-1)
\end{equation}


Change variables from $(x_1, x_2)$ to $(z_1, z_2)$ as follows: 


\begin{equation}
\left. 
\begin{aligned}
z_2= &  x_1 \\
z_1= &  \left ( x_2 - b x_1 \right ) f
\end{aligned}
\right \}
(\#eq:binormal-model-change-variable-2)
\end{equation}


Here $f$ is a quantity to be determined, which will allow us to complete the transformation to the desired bivariate integral. The second equation above can be written as:


\begin{equation}
x_2= \frac{z_1}{f} + bx_1 = \frac{z_1}{f}  + bz_2
(\#eq:binormal-model-dx2-z1-f)
\end{equation}


The Jacobian [@stein1992calculus] of the transformation is


\begin{equation}
J = \left (\begin{matrix}
0 & 1\\ 
\frac{1}{f} & b
\end{matrix}
  \right )(\#eq:binormal-model-jacobian)
\end{equation}


The magnitude of the determinant of J is $1/f$.


From a theorem in calculus [@stein1992calculus], the double integral over $(x_1,x_2)$ can be expressed in terms of a double integral over $(z_1,z_2)$ as follows:


\begin{equation}
A_{z;c} = \frac{1}{f}\int_{z_2=-\infty}^{\Phi^{-1}\left ( c \right )} dz_2 \phi\left ( z_2 \right ) \int_{z_1=-\infty}^{z_{1}^{UL}} \phi\left ( \frac{z_1}{f} + b z_2\right ) dz_1
(\#eq:binormal-model-partial-area-double-integral-2)
\end{equation}


The upper limit of the inner integral can be calculated as follows. Using the second equation in Eqn. \@ref(eq:binormal-model-change-variable-2):


\begin{equation}
z_{1}^{UL}=\left ( x_{2}^{UL} - b x_1\right ) f = \left ( a+b x_1-bx_1 \right ) f = af
(\#eq:binormal-model-upper-limit-simplified)
\end{equation}


Eqn. \@ref(eq:binormal-model-partial-area-double-integral-2) simplifies to: 


\begin{equation}
A_{z;c} = \frac{1}{f}\int_{z_2=-\infty}^{\Phi^{-1}\left ( c \right )} dz_2 \phi\left ( z_2 \right ) \int_{z_1=-\infty}^{af} \phi\left ( \frac{z_1}{f} + b z_2\right ) dz_1
(\#eq:binormal-model-partial-area-double-integral-3)
\end{equation}


Perform a change of variable from $f$ to a correlation-like quantity $\rho$ defined by:


\begin{equation}
f= \sqrt{1-\rho^2}
(\#eq:binormal-model-partial-area-rho2)
\end{equation}


Define $\rho$ in terms of the b-parameter as follows:


\begin{equation}
b\sqrt{1-\rho^2} = - \rho
(\#eq:binormal-model-partial-b-rho)
\end{equation}


This implies that $\rho$ is given by:


\begin{equation}
\rho = - \frac{b}{\sqrt{1+b^2}}
(\#eq:binormal-model-rho-final1)
\end{equation}


The argument of the right-most $\phi$ function in Eqn. \@ref(eq:binormal-model-partial-area-double-integral-3) simplifies as follows:


\begin{equation}
\frac{z_1}{f} + bz_2 = \frac{z_1 + b z_2 \sqrt{1-\rho^2}}{\sqrt{1-\rho^2}} = \frac{z_1 - \rho z_2}{\sqrt{1-\rho^2}}
(\#eq:binormal-model-partial-simplification)
\end{equation}


The expression for the partial-area under the ROC reduces to:


\begin{equation}
A_{z;c} = \frac{1}{\sqrt{1-\rho^2}}\int_{z_2=-\infty}^{\Phi^{-1}\left ( c \right )} dz_2 \phi\left ( z_2 \right ) \int_{z_1=-\infty}^{a\sqrt{1-\rho^2}} \phi\left ( \frac{z_1 - \rho z_2}{\sqrt{1-\rho^2}}\right ) dz_1
(\#eq:binormal-model-partial-area-double-integral-4)
\end{equation}


Eqn. \@ref(eq:binormal-model-partial-b-rho) implies:


\begin{equation}
1 - \rho^2 = \frac{1}{1+b^2}\\
(\#eq:binormal-model-simplification-4)
\end{equation}



Therefore,


\begin{equation}
A_{z;c} = \frac{1}{\sqrt{1-\rho^2}}\int_{z_2=-\infty}^{\Phi^{-1}\left ( c \right )} dz_2 \phi\left ( z_2 \right ) \int_{z_1=-\infty}^{\frac{a}{\sqrt{1+b^2}}} \phi\left ( \frac{z_1 - \rho z_2}{\sqrt{1-\rho^2}}\right ) dz_1
(\#eq:binormal-model-partial-area-double-integral-5)
\end{equation}


The standard bivariate normal distribution with correlation coefficient $\rho$ is defined by:


\begin{equation}
\phi\left ( z_1,z_2;\rho \right ) = \frac{1}{2 \pi \sqrt{1-\rho^2}} \exp  \left ( -\frac{z_1^2 - 2 \rho z_1 z_2 + z_2^2}{2(1-\rho^2)}  \right )
(\#eq:binormal-model-definition-standard-bivariate-distribution)
\end{equation}


The standard normal distribution is defined by:

\begin{equation}
\phi\left ( z\right ) = \frac{1}{\sqrt{2 \pi }} \exp  \left ( -\frac{z^2}{2}  \right )
(\#eq:binormal-model-definition-standard-normal-distribution)
\end{equation}

It can be shown using these definitions that:



\begin{equation}
\phi\left ( z_1,z_2;\rho \right ) = \frac{1}{\sqrt{1-\rho^2}}\phi\left ( z_2 \right )\phi\left ( \frac{z_1-\rho z_2}{\sqrt{1-\rho^2}} \right )
(\#eq:binormal-model-standard-bivariate-distribution-2)
\end{equation}


Using this form the expression for the partial-area is:


\begin{equation}
A_{z;c} = \int_{z_2=-\infty}^{\Phi^{-1}\left ( c \right )}   \int_{z_1=-\infty}^{\frac{a}{\sqrt{1+b^2}}} \phi\left ( z_1,z_2;\rho \right ) dz_1dz_2
(\#eq:binormal-model-partial-area-final1)
\end{equation}



### Special case (total area) {#binormal-model-appendix-1-total-auc}

Since $c$ is the upper limit of FPF, setting $c = 1$ yields the total area under the binormal ROC curve ^[Since the integral over $z_2$ is over the entire range it integrates out to unity leaving the one-dimensional density function $\phi\left ( z_1 \right )$ inside the integral. The last step follows from the definition of the $\Phi$ function.]:


\begin{equation}
\left. 
\begin{aligned}
A_{z} = & \int_{z_2=-\infty}^{\infty}   \int_{z_1=-\infty}^{\frac{a}{\sqrt{1+b^2}}} \phi\left ( z_1,z_2;\rho \right ) dz_1dz_2 \\
= &  \int_{z_1=-\infty}^{\frac{a}{\sqrt{1+b^2}}} \phi\left ( z_1 \right ) dz_1 \\
= & \Phi\left ( \frac{a}{\sqrt{1+b^2}} \right )
\end{aligned}
\right \}
(\#eq:binormal-model-partial-area-special-case)
\end{equation}


An equivalent forms for the total area under the unequal variance binormal ROC curve is:


\begin{equation}
\left. 
\begin{aligned}
A_{z} =  & \Phi\left ( \frac{a}{\sqrt{1+b^2}} \right ) \\
= & \Phi\left ( \frac{\frac{a}{b}}{\sqrt{1+\frac{1}{b^2}}} \right ) \\
= & \Phi\left ( \frac{\mu}{\sqrt{1 + \sigma^2}} \right )
\end{aligned}
\right \}
(\#eq:binormal-model-partial-area-special-case-2)
\end{equation}

## Appendix III: Invariance property of pdfs {#binormal-model-invariance-property-appendix-3}
The binormal model is not as restrictive as might appear at first sight. Any monotone increasing transformation $Y=f(Z)$ applied to the observed z-samples, and the associated thresholds, will yield the same observed data, e.g., Table \@ref(tab:ratings-paradigmExampleTable). This is because such a transformation leaves the ordering of the ratings unaltered and hence results in the same operating points. While the distributions for $Y$ will not be binormal (i.e., two independent normal distributions), one can safely "pretend" that one is still dealing with an underlying binormal model. An alternative way of stating this is that any pair of distributions is allowed as long as they are reducible to a binormal model form by a monotonic increasing transformation of  Y: e.g., $Z=f^{-1}$. [If $f$ is a monotone increasing function of its argument, so is  $f^{-1}$}.]  For this reason, the term “pair of latent underlying normal distributions” is sometimes used to describe the binormal model. The robustness of the binormal model has been investigated [@RN1216; @RN100]. The referenced paper by Dorfman et al has an excellent discussion of the robustness of the binormal model.

The robustness of the binormal model, i.e., the flexibility allowed by the infinite choices of monotonic increasing functions, application of each of which leaves the ordering of the data unaltered, is widely misunderstood. The non-Gaussian appearance of histograms of ratings in ROC studies can lead one to incorrect conclusions that the binormal model is inapplicable to these datasets. To quote a reviewer of one of my recent papers:

> I have had multiple encounters with statisticians who do not understand this difference.... They show me histograms of data, and tell me that the data is obviously not normal, therefore the binormal model should not be used. 

The reviewer is correct. The misconception is illustrated next.





```r
# shows that monotone transformations have no effect on 
# AUC even though the pdfs look non-gaussian
# common misconception about ROC analysis
fArray <- c(0.1,0.5,0.9)
seedArray <- c(10,11,12)
for (row in 1:3) { 
  f <- fArray[row]
  seed <- seedArray[row]
  set.seed(seed) 
  # numbers of cases simulated
  K1 <- 900
  K2 <- 1000
  mu1 <- 30
  sigma1 <- 7
  mu2 <- 55
  sigma2 <- 7 
  # Simulate true gaussian ratings using above parameter values
  z1 <- rnorm(K1,mean = mu1,sd = sigma1)
  z1[z1>100] <- 100;z1[z1<0] <- 0 # constrain to 0 to 100
  z2 <- rnorm(K2,mean = mu2,sd = sigma2)
  z2[z2>100] <- 100;z2[z2<0] <- 0 # constrain to 0 to 100
  # calculate AUC for true Gaussian ratings
  AUC1 <- TrapezoidalArea(z1, z2)
  Gaussians <- c(z1, z2)
  # display histograms of true Gaussian ratings, A1, A2 or A3
  x <- data.frame(x=Gaussians) #  line 27
  x <-  
    ggplot(data = x, mapping = aes(x = x)) +
    geom_histogram(binwidth = 5, color = "black", fill="grey") + 
    xlab(label = "Original Rating") + 
    ggtitle(label = paste0("A", row, ": ", "Gaussians"))
  print(x)
  z <- seq(0.0, 100, 0.1)
  # transform the latent Gaussians to true Gaussians
  transformation <- 
    data.frame(
      x = z, 
      z =  Y(z,mu1,mu2,sigma1,sigma2,f))
  # display transformation functions, B1, B2 or B3
  x <- 
    ggplot(mapping = aes(x = x, y = z)) + 
    geom_line(data = transformation, size = 1) +
    xlab(label = "Original Rating") +
    ylab(label = "Transformed Rating") + 
    ggtitle(label = paste0("B", row, ": ","Monotone Transformation"))
  print(x)
  y <- Y(c(z1, z2),mu1,mu2,sigma1,sigma2,f)
  y1 <- y[1:K1];y2 <- y[(K1+1):(K1+K2)]
  # calculate AUC for transformed ratings
  AUC2 <- TrapezoidalArea( y1, y2)
  # display histograms of latent Gaussian ratings, C1, C2 or C3
  x <- data.frame(x=y)
  x <-  ggplot(data = x, mapping = aes(x = x)) +
    geom_histogram(binwidth = 5, color = "black", fill="grey") +
    xlab(label = "Transformed Rating") + 
    ggtitle(label = paste0("C", row, ": ", "Latent Gaussians"))
  print(x)
# print AUCs, note they are identical (for each row)  
options(digits = 9)
  cat("row =", row, ", seed =", seed, ", f =", f, 
      "\nAUC of actual Gaussians =", AUC1, 
      "\nAUC of latent Gaussians =", AUC2, "\n")
}
#> row = 1 , seed = 10 , f = 0.1 
#> AUC of actual Gaussians = 0.99308 
#> AUC of latent Gaussians = 0.99308
#> row = 2 , seed = 11 , f = 0.5 
#> AUC of actual Gaussians = 0.993668889 
#> AUC of latent Gaussians = 0.993668889
#> row = 3 , seed = 12 , f = 0.9 
#> AUC of actual Gaussians = 0.995041111 
#> AUC of latent Gaussians = 0.995041111
```

<img src="06-binormal-model_files/figure-html/unnamed-chunk-15-1.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-2.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-3.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-4.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-5.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-6.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-7.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-8.png" width="33%" /><img src="06-binormal-model_files/figure-html/unnamed-chunk-15-9.png" width="33%" />

**Figure captions (A1 - C3):** Illustrating the invariance of ROC analysis to arbitrary monotone transformations of the ratings. Each row contains 3 plots: labeled 1, 2 and 3. Each column contains 3 plots labeled A, B and C. So, for example, plot C2 refers to the second row and third column. The for-loop generates the plot one row at a time. Each of the latent Gaussian plots C1, C2 and C3 appears not binormal. However, using the inverse of the monotone transformations shown B1, B2 and B3, they can be transformed to the binormal model histograms A1, A2 and A3. Plot A1 shows the histogram of simulated ratings from a binormal model. Two peaks, one at 30 and the other at 55 are evident (by design, all ratings in this figure are in the range 0 to 100). Plot B1 shows the monotone transformation for $f = 0.1$. Plot C1 shows the histogram of the transformed rating. The choice of $f$ leads to a transformed rating histogram that is peaked near the high end of the rating scale. For A1 and C1 the corresponding AUCs are identical (0.993080000). Plot A2 is for a different seed value, plot B2 is the transformation for $f = 0.5$ and now the transformed histogram is almost flat, plot C2. For plots A2 and C2 the corresponding AUCs are identical (0.993668889). Plot A3 is for a different seed value, B3 is the transformation for $f = 0.9$ and the transformed histogram C3 is peaked near the low end of the transformed rating scale. For plots A3 and (C3) the corresponding AUCs are identical (0.995041111).

The idea is to simulate continuous ratings data in the range 0 to 100 from a binormal model. $K_1 = 900$ non-diseased cases are sampled from a Gaussian centered at $\mu_1$ = 30 and standard deviation $\sigma_1 = 7$. $K_2 = 1000$ diseased cases are sampled from a Gaussian centered at $\mu_2$ = 55 and standard deviation $\sigma_2$ = 7. The variable $f$, which is in the range (0,1), controls the shape of the transformed distribution. If $f$ is small, the transformed distribution will be peaked towards 0 and if $f$ is unity, it will be peaked at 100. If $f$ equals 0.5, the transformed distribution is flat. Insight into the reason for this transformation is in [@RN300], Chapter 7: it has to do with transformations of random variables. The transformation function, $Y(Z)$, implements:

\begin{equation} 
Y\left ( Z \right )=\left [ \left ( 1-f \right )\Phi\left ( \frac{Z-\mu_1}{\sigma_1} \right )+f\Phi\left ( \frac{Z-\mu_2}{\sigma_2} \right ) \right ]100
(\#eq:binormal-modelDemoMisconception)
\end{equation}

The multiplication by 100 ensures that the transformed variable is in the range 0 to 100 (if not, it is code-constrained to be). The code realizes the random samples, calculates the empirical AUC, displays the histogram of the true binormal samples, plots the transformation function, calculates the empirical AUC using the transformed samples, and plots the histogram of the transformed samples (the latent binormal). 

* B1 shows the transformation for $f = 0.1$. The steep initial rise of the curve has the effect of flattening the histogram of the transformed ratings at the low end of the rating scale, C1. Conversely, the flat nature of the curve near upper end of the rating range has the effect of causing the histogram of the transformed variable to peak in that range. 
* B2 shows the transformation for $f = 0.5$. This time the latent rating histogram, C2, is almost flat over the entire range, definitely not visually binormal. 
* B3 shows the transformation for $f = 0.9$. This time the transformed rating histogram, C3, is peaked at the low end of the transformed rating scale.
* The output lists the values of the seed variable and the value of the shape parameter $f$. *For each value of seed and the shape parameter, the AUCs of the actual Gaussians and the transformed variables are identical*. 
* The values of the parameters were chosen to best illustrate the true binormal nature of the plots A2 and A3. This has the effect of making the AUCs close to unity. 

The histograms in C1, C2 and C3 appear to be non-Gaussian. The corresponding non-diseased and diseased ratings will fail tests of normality. [Showing this is left as an exercise for the reader.] Nevertheless, they are latent Gaussians in the sense that the inverses of the transformations shown in B1, B2 and B3 will yield histograms that are strictly binormal, i.e., A1, A2 and A3. By appropriate changes to the monotone transformation function, the histograms shown in C1, C2 and C3 can be made to resemble a wide variety of shapes, for example, quasi-bimodal (don't confuse bimodal with binormal) histograms.] 

**Visual examination of the shape of the histograms of ratings, or standard tests for normality, yield little, if any, insight into whether the underlying binormal model assumptions are being violated.**



## Appendix IV: Fitting an ROC curve {#binormal-model-curve-fitting-appendix-4}
### JAVA fitted ROC curve
This section, described in the physical book, has been abbreviated to a [relevant website](http://www.rad.jhmi.edu/jeng/javarad/roc/JROCFITi.html).

### Simplistic straight line fit to the ROC curve
To be described next is a method for fitting data such as in Table \@ref(tab:ratings-paradigmExampleTable) to the binormal model, i.e., determining the parameters $(a,b)$ and the thresholds $\zeta_r , \quad r = 1, 2, ..., R-1$, to best fit, in some to-be-defined sense, the observed cell counts. The most common method uses an algorithm called maximum likelihood. But before getting to that, I describe the least-square method, which is conceptually simpler, but not really applicable, as will be explained shortly.

#### Least-squares estimation
By applying the function $\Phi^{-1}$  to both sides of Eqn. \@ref(eq:binormal-model-roc-curve1), one gets (the "inverse" function cancels the "forward" function on the right hand side):

\begin{equation*} 
\Phi^{-1}\left ( TPF \right ) = a + b \Phi^{-1}\left ( FPF \right )
\end{equation*}

This suggests that a plot of $y = \Phi^{-1}\left ( TPF \right )$  vs. $x=\Phi^{-1}\left ( FPF \right )$ is expected to follow a straight line with slope $b$ and intercept $a$. Fitting a straight line to such data is generally performed by the method of least-squares, a capability present in most software packages and spreadsheets. Alternatively, one can simply visually draw the best straight line that fits the points, memorably referred to [@RN300] as "chi-by-eye". This was the way parameters of the binormal model were estimated prior to Dorfman and Alf's work [@RN1081]. The least-squares method is a quantitative way of accomplishing the same aim. If $\left ( x_t,y_t \right )$ are the data points, one constructs $S$, the sum of the squared deviations of the observed ordinates from the predicted values (since $R$ is the number of ratings bins, the summation runs over the $R-1$ operating points):

\begin{equation*} 
S  = \sum_{i=1}^{R-1}\left ( y_i - \left ( a + bx_i \right ) \right )^2
\end{equation*} 

The idea is to minimize S with respect to the parameters $(a,b)$. One approach is to differentiate this with respect to $a$ and $b$ and equate each resulting derivate expression to zero. This yields two equations in two unknowns, which are solved for $a$ and $b$. If the reader has never done this before, one should go through these steps at least once, but it would be smarter in future to use software that does all this. In `R` the least-squares fitting function is `lm(y~x)`, which in its simplest form fits a linear model `lm(y~x)`  using the method of least-squares (in case you are wondering `lm` stands for linear model, a whole branch of statistics in itself; in this example one is using its simplest capability).


```r
# ML estimates of a and b (from Eng JAVA program)
# a <- 1.3204; b <- 0.6075 
# # these are not used in program; just here for comparison

FPF <- c(0.017, 0.050, 0.183, 0.5)  
# this is from Table 6.11, last two rows
TPF <- c(0.440, 0.680, 0.780, 0.900)
# ...do...

PhiInvFPF <- qnorm(FPF)
# apply the PHI_INV function
PhiInvTPF <- qnorm(TPF)
# ... do ... 

fit <- lm(PhiInvTPF~PhiInvFPF)
print(fit)
#> 
#> Call:
#> lm(formula = PhiInvTPF ~ PhiInvFPF)
#> 
#> Coefficients:
#> (Intercept)    PhiInvFPF  
#>    1.328844     0.630746
```


<div class="figure">
<img src="06-binormal-model_files/figure-html/binormal-modelStraightLineFit-1.png" alt="The straight line fit method of estimating parameters of the fitting model." width="672" />
<p class="caption">(\#fig:binormal-modelStraightLineFit)The straight line fit method of estimating parameters of the fitting model.</p>
</div>

Fig. \@ref(fig:binormal-modelStraightLineFit) shows operating points from Table \@ref(tab:ratings-paradigmExampleTable), transformed by the $\Phi^{-1}$ function; the slope of the line is the least-squares estimate of the $b$ parameter and the intercept is the corresponding $a$ parameter of the binormal model.

The last line contains the least squares estimated values, $a$ = 1.3288 and $b$ = 0.6307. The corresponding maximum likelihood estimates of these parameters, as yielded by the Eng web code, Appendix B, are listed in line 4 of the main program: $a$ = 1.3204 and $b$ = 0.6075. The estimates appear to be close, particularly the estimate of $a$ , but there are a few things wrong with the least-squares approach. First, the method of least squares assumes that the data points are independent. Because of the manner in which they are constructed, namely by cumulating points, the independence assumption is not valid for ROC operating points. Cumulating the 4 and 5 responses constrains the resulting operating point to be above and to the right of the point obtained by cumulating the 5 responses only, so the data points are definitely not independent. Similarly, cumulating the 3, 4 and 5 responses constrains the resulting operating point to be above and to the right of the point obtained by cumulating the 4 and 5 responses, and so on. The second problem is the linear least-squares method assumes there is no error in measuring x; the only source of error that is accounted for is in the y-coordinate. In fact, both coordinates of an ROC operating point are subject to sampling error. Third, disregard of error in the x-direction is further implicit in the estimates of the thresholds, which according to Eqn. (6.2.19), is given by:

\begin{equation*} 
\zeta_r = - \Phi^{-1}\left ( FPF_r \right )
\end{equation*} 

These are "rigid" estimates that assume no error in the FPF values. As was shown in Chapter \@ref(binary-task), 95% confidence intervals apply to these estimates.

A historical note: prior to computers and easy access to statistical functions the analyst had to use a special plotting paper, termed “double probability paper”, that converted probabilities into x and y distances using the inverse function.

### Maximum likelihood estimation (MLE)
The approach taken by Dorfman and Alf was to maximize the likelihood function instead of S. The likelihood function is the probability of the observed data given a set of parameter values, i.e.,

\begin{equation*} 
\text {L} \equiv P\left ( data \mid \text {parameters} \right )
\end{equation*} 

Generally "data" is suppressed, so likelihood is a function of the parameters; but "data" is always implicit. With reference to Fig. 6.1, the probability of a non-diseased case yielding a count in the 2nd bin equals the area under the curve labeled "Noise" bounded by the vertical lines at $\zeta_1$ and $\zeta_2$. In general, the probability of a non-diseased case yielding a count in the $r^\text{th}$ bin equals the area under the curve labeled "Noise" bounded by the vertical lines at  $\zeta_{r-1}$ and $\zeta_r$.  Since the area to the left of a threshold is the CDF corresponding to that threshold, the required probability is $\Phi\left ( \zeta_r \right ) - \Phi\left ( \zeta_{r-1} \right )$; we are simply subtracting two expressions for specificity, Eqn. (6.2.5). 

\begin{equation*} 
\text {count in non-diseased bin } r = \Phi\left ( \zeta_r \right ) - \Phi\left ( \zeta_{r-1} \right )
\end{equation*} 

Similarly, the probability of a diseased case yielding a count in the rth bin equals the area under the curve labeled "Signal" bounded by the vertical lines at  $\zeta_{r-1}$ and $\zeta_r$. The area under the diseased distribution to the left of threshold $\zeta_r$  is the $1 - TPF$ at that threshold:

\begin{equation*} 
1 - \Phi\left ( \frac{\mu-\zeta_r}{\sigma} \right ) = \Phi\left ( \frac{\zeta_r - \mu}{\sigma} \right )
\end{equation*} 

The area between the two thresholds is:

\begin{align*} 
P\left ( \text{count in diseased bin }r \right ) &= \Phi\left ( \frac{\zeta_r - \mu}{\sigma} \right ) - \Phi\left ( \frac{\zeta_{r-1} - \mu}{\sigma} \right ) \\
&= \Phi\left ( b\zeta_r-a \right ) - \Phi\left ( b\zeta_{r-1}-a \right )
\end{align*}

Let $K_{1r}$  denote the number of non-diseased cases in the rth bin, and $K_{2r}$  denotes the number of diseased cases in the rth bin. Consider the number of counts $K_{1r}$ in non-diseased case bin $r$. Since the probability of each count is  $\Phi\left ( \zeta_{r+1}  \right ) - \Phi\left ( \zeta_r  \right )$, the probability of the observed number of counts, assuming the counts are independent, is  ${\left(\Phi\left ( \zeta_{r+1}  \right ) - \Phi\left ( \zeta_r  \right )  \right )}^{K_{1r}}$. Similarly, the probability of observing  counts in diseased case bin $r$ is ${\left (\Phi\left ( b\zeta_{r+1}-a  \right ) - \Phi\left ( b\zeta_r-a  \right )  \right )}^{K_{2r}}$, subject to the same independence assumption. The probability of simultaneously observing $K_{1r}$  counts in non-diseased case bin r and $K_{2r}$  counts in diseased case bin $r$ is the product of these individual probabilities (again, an independence assumption is being used): 

\begin{equation*} 
\left (\Phi\left ( \zeta_{r+1}  \right ) - \Phi\left ( \zeta_r  \right )  \right )^{K_{1r}} \left (\Phi\left ( b\zeta_{r+1}-a  \right ) - \Phi\left ( b\zeta_r-a  \right )  \right )^{K_{2r}}
\end{equation*} 

Similar expressions apply for all integer values of $r$ ranging from $1,2,...,R$. Therefore the probability of observing the entire data set is the product of expressions like Eqn. (6.4.5), over all values of $r$:

\begin{equation} 
\prod_{r=1}^{R}\left [\left (\Phi\left ( \zeta_{r+1}  \right ) - \Phi\left ( \zeta_r  \right )  \right )^{K_{1r}} \left (\Phi\left ( b\zeta_{r+1}-a  \right ) - \Phi\left ( b\zeta_r-a  \right )  \right )^{K_{2r}}  \right ]
(\#eq:binormal-modelProductProb)
\end{equation} 

We are almost there. A specific combination of $K_{11},K_{12},...,K_{1R}$ counts from $K_1$ non-diseased cases and counts $K_{21},K_{22},...,K_{2R}$ from $K_2$ diseased cases can occur the following number of times (given by the multinomial factor shown below):

\begin{equation} 
\frac{K_1!}{\prod_{r=1}^{R}K_{1r}!}\frac{K_2!}{\prod_{r=1}^{R}K_{2r}!}
(\#eq:binormal-modelCombFactor)
\end{equation} 

The likelihood function is the product of Eqn. \@ref(eq:binormal-modelProductProb) and Eqn.  \@ref(eq:binormal-modelCombFactor):

\begin{equation} 
\begin{split}
L\left ( a,b,\overrightarrow{\zeta} \right ) &= \left (\frac{K_1!}{\prod_{r=1}^{R}K_{1r}!}\frac{K_2!}{\prod_{r=1}^{R}K_{2r}!}  \right ) \times \\
&\quad\prod_{r=1}^{R}\left [\left (\Phi\left ( \zeta_{r+1}  \right ) - \Phi\left ( \zeta_r  \right )  \right )^{K_{1r}} \left (\Phi\left ( b\zeta_{r+1}-a  \right ) - \Phi\left ( b\zeta_r-a  \right )  \right )^{K_{2r}}  \right ]
\end{split}
(\#eq:binormal-modelLikelihood)
\end{equation}

The left hand side of Eqn. \@ref(eq:binormal-modelLikelihood) shows explicitly the dependence of the likelihood function on the parameters of the model, namely $a,b,\overrightarrow{\zeta}$, where the vector of thresholds $\overrightarrow{\zeta}$ is a compact notation for the set of thresholds $\zeta_1,\zeta_2,...,\zeta_R$, (note that since $\zeta_0 = -\infty$, and $\zeta_R = +\infty$, only $R-1$ free threshold parameters are involved, and the total number of free parameters in the model is $R+1$). For example, for a 5-rating ROC study, the total number of free parameters is 6, i.e., $a$, $b$ and 4 thresholds $\zeta_1,\zeta_2,\zeta_3,\zeta_4$.

Eqn.  \@ref(eq:binormal-modelLikelihood) is forbidding but here comes a simplification. The difference of probabilities such as $\Phi\left ( \zeta_r  \right )-\Phi\left ( \zeta_{r-1}  \right )$ is guaranteed to be positive and less than one [the $\Phi$ function is a probability, i.e., in the range 0 to 1, and since $\zeta_r$ is greater than $\zeta_{r-1}$, the difference is positive and less than one]. When the difference is raised to the power of $K_{1r}$ (a non-negative integer) a very small number can result. Multiplication of all these small numbers may result in an even smaller number, which may be too small to be represented as a floating-point value, especially as the number of counts increases. To prevent this we resort to a trick. Instead of maximizing the likelihood function $L\left ( a,b,\overrightarrow{\zeta} \right )$ we choose to maximize the logarithm of the likelihood function (the base of the logarithm is immaterial). The logarithm of the likelihood function is: 

\begin{equation} 
LL\left ( a,b,\overrightarrow{\zeta} \right )=\log \left ( L\left ( a,b,\overrightarrow{\zeta} \right ) \right )
(\#eq:binormal-modelLogLikelihood)
\end{equation} 

Since the logarithm is a monotonically increasing function of its argument, maximizing the logarithm of the likelihood function is equivalent to maximizing the likelihood function. Taking the logarithm converts the product symbols in Eqn. (6.4.8) to summations, so instead of multiplying small numbers one is adding them, thereby avoiding underflow errors. Another simplification is that one can ignore the logarithm of the multinomial factor involving the factorials, because these do not depend on the parameters of the model. Putting all this together, we get the following expression for the logarithm of the likelihood function:

\begin{equation} 
\begin{split}
LL\left ( a,b,\overrightarrow{\zeta} \right ) \propto& \sum_{r=1}^{R} K_{1r}\log \left ( \Phi\left ( \zeta_{r+1} \right ) - \Phi\left ( \zeta_r \right ) \right ) \\
&+ \sum_{r=1}^{R} K_{2r}\log \left ( \Phi\left (b \zeta_{r+1} - a \right ) - \Phi\left ( b \zeta_r - a \right ) \right ) 
\end{split}
(\#eq:binormal-modelLL)
\end{equation} 

The left hand side of Eqn. \@ref(eq:binormal-modelLL) is a function of the model parameters $a,b,\overrightarrow{\zeta}$ and the observed data, the latter being the counts contained in the vectors $\overrightarrow{K_1}$ and $\overrightarrow{K_2}$, where the vector notation is used as a compact form for the counts $K_{11},K_{12},...,K_{1R}$ and $K_{21},K_{22},...,K_{2R}$, respectively. The right hand side of Eqn. \@ref(eq:binormal-modelLL) is monotonically related to the probability of observing the data given the model parameters $a,b,\overrightarrow{\zeta}$. If the choice of model parameters is poor, then the probability of observing the data will be small and log likelihood will be small. With a better choice of model parameters the probability and log likelihood will increase. With optimal choice of model parameters the probability and log likelihood will be maximized, and the corresponding optimal values of the model parameters are called maximum likelihood estimates (MLEs). These are the estimates produced by the programs RSCORE and ROCFIT. 

### Code implementing MLE


```r

# ML estimates of a and b (from Eng JAVA program)
# a <- 1.3204; b <- 0.6075 
# these are not used in program; just there for comparison

K1t <- c(30, 19, 8, 2, 1)
K2t <- c(5,  6, 5, 12, 22)
dataset <- Df2RJafrocDataset(K1t, K2t, InputIsCountsTable = TRUE)
retFit <- FitBinormalRoc(dataset)
retFit[1:5]
#> $a
#> [1] 1.32045261
#> 
#> $b
#> [1] 0.607492932
#> 
#> $zetas
#>      zetaFwd1      zetaFwd2      zetaFwd3      zetaFwd4 
#> 0.00768054675 0.89627306763 1.51564784976 2.39672209865 
#> 
#> $AUC
#> [1] 0.870452157
#> 
#> $StdAUC
#>              [,1]
#> [1,] 0.0379042262
print(retFit$fittedPlot)
```

<img src="06-binormal-model_files/figure-html/unnamed-chunk-17-1.png" width="672" />

Note the usage of the `RJafroc` package [@R-RJafroc]. Specifically, the function `FitBinormalRoc`. The ratings table is converted to an `RJafroc` dataset object, followed by application of the fitting function. The results, contained in `retFit` should be compared to those obtained from the [website implementation of ROCFIT](http://www.rad.jhmi.edu/jeng/javarad/roc/JROCFITi.html).

## Appendix V: Validating fitting model {#binormal-model-curve-fitting-validation-appendix-5}
The above ROC curve is a good visual fit to the observed operating points. Quantification of the validity of the fitting model is accomplished by calculating the Pearson goodness-of-fit test [@RN2656], also known as the chi-square test, which uses the statistic defined by [@RN1492]:

\begin{equation} 
C^2=\sum_{t=1}^{2}\sum_{r=1}^{R}\frac{\left (K_{tr}-\left \langle K_{tr} \right \rangle  \right )^2}{\left \langle K_{tr} \right \rangle}\\
K_{tr} \geq 5
(\#eq:binormal-modelGoodnessFit)
\end{equation} 

The expected values are given by:

\begin{equation}
\begin{split}
\left \langle K_{1r} \right \rangle &=K_1\left ( \Phi\left ( \zeta_{r+1} \right ) - \Phi\left ( \zeta_r \right )  \right ) \\
\left \langle K_{2r} \right \rangle &=K_2\left ( \Phi\left ( a\zeta_{r+1}-b \right ) - \Phi\left ( a\zeta_r - b\right )  \right )
\end{split}
(\#eq:binormal-modelGoodnessFitExpVals)
\end{equation} 

These expressions should make sense: the difference between the two CDF functions is the probability of a count in the specified bin, and multiplication by the total number of relevant cases should yield the expected counts (a non-integer). 

It can be shown that under the null hypothesis that the assumed probability distribution functions for the counts equals the true probability distributions, i.e., the model is valid, the statistic $C^2$ is distributed as:

\begin{equation} 
C^2\sim \chi_{df}^{2}
(\#eq:binormal-modelGoodnessFitDistr)
\end{equation} 

Here $C^2\sim \chi_{df}^{2}$  is the chi-square distribution with degrees of freedom  *df* defined by:

\begin{equation} 
df=\left ( R-1 \right )+\left ( R-1 \right )-\left (2+ R-1 \right )=\left ( R-3 \right )
(\#eq:binormal-modelGoodnessFitdf)
\end{equation} 

The right hand side of the above equation has been written in an expansive form to illustrate the general rule: for $R$ non-diseased cells in the ratings table, the degree of freedom is $R-1$: this is because when all but one cells are specified, the last is determined, because they must sum to $K_1$ . Similarly, the degree of freedom for the diseased cells is also $R-1$. Last, we need to subtract the number of free parameters in the model, which is $(2+R-1)$, i.e., the  $a,b$ parameters and the $R-1$ thresholds. It is evident that if $R = 3$ then $df = 0$. In this situation, there are only two non-trivial operating points and the straight-line fit shown will pass through both of them. With two basic parameters, fitting two points is trivial, and goodness of fit cannot be calculated.

Under the null hypothesis (i.e., model is valid) $C^2$ is distributed as $\chi_{df}^{2}$. Therefore, one computes the probability that this statistic is larger than the observed value, called the *p-value*. If this probability is very small, that means that the deviations of the observed values of the cell counts from the expected values are so large that it is unlikely that the model is correct. The degree of unlikeliness is quantified by the p-value. Poor fits lead to small p values.

At the 5% significance level, one concludes that the fit is not good if $p < 0.05$. In practice one occasionally accepts smaller values of $p$, $p > 0.001$ before completely abandoning a model. It is known that adoption of a stricter criterion, e.g., $p > 0.05$, can occasionally lead to rejection of a retrospectively valid model [@RN300].

### Estimating the covariance matrix
TBA See book chapter 6.4.3. This is implemented in `RJafroc.`

### Estimating the variance of Az
TBA See book chapter 6.4.4. This is implemented in `RJafroc`.


## References {#binormal-model-references}



