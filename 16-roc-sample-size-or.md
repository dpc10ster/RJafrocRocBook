# Sample size estimation for ROC studies  OR method {#roc-sample-size-or}






## TBA How much finished {#roc-sample-size-or-how-much-finished}
70%




## Introduction {#roc-sample-size-or-introduction}

 

## Statistical Power {#StatPower2}


\begin{equation}
Power = 1 - \beta
(\#eq:DefinitionStatPower1)
\end{equation}


### Sample size estimation for random-reader random-cases
For convenience the OR model is repeated below with the case-set index suppressed:

\begin{equation}
Y_{n(ijk)}=\mu+\tau_i+R_j+C_k+(\tau R)_{ij}+(\tau C)_{ik}+(RC)_{jk}+(\tau RC)_{ijk}+\epsilon_ {n(ijk)}
(\#eq:ORModelSsOR)
\end{equation}

As usual, the treatment effects $\tau_i$  are subject to the constraint that they sum to zero. The observed effect size (a random variable) is defined by:

\begin{equation}
d=\theta_{1\bullet}-\theta_{2\bullet}
(\#eq:EffectSize1)
\end{equation}

It is a realization of a random variable, so one has some leeway in the choice of anticipated effect size. In the significance-testing procedure described in TBA Chapter 09 interest was in the distribution of the F-statistic when the NH is true. For sample size estimation, one needs to know the distribution of the statistic when the NH is false. It was shown that then the observed F-statistic TBA Eqn. (9.35) is distributed as a non-central F-distribution  $F_{ndf,ddf,\Delta}$ with non-centrality parameter $\Delta$: 

\begin{equation}
F_{DBM|AH} \sim F_{ndf,ddf,\Delta}
(\#eq:FDBMSampling1)
\end{equation}

The non-centrality parameter   was defined, Eqn. TBA (9.34), by:

\begin{equation}
\Delta=\frac{JK\sigma_{Y;\tau}^2}{\left ( \sigma_{Y;\epsilon}^2 + \sigma_{Y;\tau RC}^2 \right )+K\sigma_{Y;\tau R}^2+J\sigma_{Y;\tau C}^2}
(\#eq:DefDelta1)
\end{equation}

To minimize confusion, this equation has been rewritten here using the subscript $Y$ to explicitly denote pseudo-value derived quantities (in TBA Chapter 09 this subscript was suppressed. 

The estimate of $\sigma_{Y;\tau C}^2$ can turn out to bee negative. To avoid a negative denominator, Hillis suggests the following modification:

\begin{equation}
\Delta=\frac{JK\sigma_{Y;\tau}^2}{\left ( \sigma_{Y;\epsilon}^2 + \sigma_{Y;\tau RC}^2 \right )+K\sigma_{Y;\tau R}^2+\max \left (J\sigma_{Y;\tau C}^2 ,0 \right )}
(\#eq:DefDeltaHillis1)
\end{equation}

This expression depends on three variance components, $(\sigma_{Y;\epsilon}^2 + \sigma_{Y;\tau RC}^2)$ - the two terms are inseparable - $\sigma_{Y;\tau R}^2$ and $\sigma_{Y;\tau C}^2$. The $ddf$ term appearing in TBA Eqn. (11.4) was defined by TBA Eqn. (9.24) - this quantity does not change between NH and AH:

\begin{equation}
ddf_H=\frac{\left [MSTR+\max(MSTR-MSTRC,0)  \right ]^2}{\frac{[MSTR]^2}{(I-1)(J-1)}}
(\#eq:ddfH1)
\end{equation}

The mean squares in this expression can be expressed in terms of the three variance-components appearing in TBA Eqn. (11.6). Hillis and Berbaum [@RN1476] have derived these expression and they will not be repeated here (Eqn. 4 in the cited reference). RJafroc implements a function to calculate the mean squares, `UtilMeanSquares()`, which allows ddf to be calculated using Eqn. TBA (11.7). The sample size functions in this package need only the three variance-components (the formula for $ddf_H$ is implemented internally). 

For two treatments, since the individual treatment effects must be the negatives of each other (because they sum to zero), it is easily shown that:

\begin{equation}
\sigma_{Y;\tau}^2=\frac{d^2}{2}
(\#eq:sigma2Tau1)
\end{equation}
 

### Dependence of statistical power on estimates of model parameters {#roc-sample-size-or-dependence-of-stats-power}
Examination of the expression for  , Eqn. (11.5), shows that statistical power increases if:

* The numerator is large. This occurs if: (a) the anticipated effect-size $d$ is large. Since effect-size enters as the *square*, TBA Eqn. (11.8), it is has a particularly strong effect; (b) If $J \times K$ is large. Both of these results should be obvious, as a large effect size and a large sample size should result in increased probability of rejecting the NH. 
* The denominator is small. The first term in the denominator is  $\left ( \sigma_{Y;\epsilon}^2 + \sigma_{Y;\tau RC}^2 \right )$. These two terms cannot be separated. This is the residual variability of the jackknife pseudovalues. It should make sense that the smaller the variability, the larger is the non-centrality parameter and the statistical power. 
* The next term in the denominator is $K\sigma_{Y;\tau R}^2$, the treatment-reader variance component multiplied by the total number of cases. The reader variance $\sigma_{Y;R}^2$ has no effect on statistical power, because it has an equal effect on both treatments and cancels out in the difference. Instead, it is the treatment-reader variance $\sigma_{Y;R}^2$  that contributes "noise" tending to confound the estimate of the effect-size. 
* The variance components estimated by the ANOVA procedure are realizations of random variables and as such subject to noise (there actually exists a beast such as variance of a variance). The presence of the $K$ term, usually large, can amplify the effect of noise in the estimate of $\sigma_{Y;R}^2$, making the sample size estimation procedure less accurate.
* The final term in the denominator is  $J\sigma_{Y;\tau C}^2$. The variance $\sigma_{Y;C}^2$ has no impact on statistical power, as it cancels out in the difference. The treatment-case variance component introduces "noise" into the estimate of the effect size, thereby decreasing power. Since it is multiplied by J, the number of readers, and typically $J<<K$, the error amplification effect on accuracy of the sample size estimate is not as bad as with the treatment-reader variance component.
* Accuracy of sample size estimation, essentially estimating confidence intervals for statistical power, is addressed in [@RN2027].

### Formulae for random-reader random-case (RRRC) sample size estimation {#roc-sample-size-or-RRRC-sample-size-estimation}


### Significance testing {#roc-sample-size-or-sig-testing}

### p-value and confidence interval {#roc-sample-size-or-pvalue-ci}

### Comparing DBM to Obuchowski and Rockette for single-reader multiple-treatments {#roc-sample-size-or-CompareDBM2OR}
Having performed a pilot study and planning to perform a pivotal study, sample size estimation follows the following procedure, which assumes that both reader and case are treated as random factors. Different formulae, described later, apply when either reader or case is treated as a fixed factor.

* Perform OR analysis on the pilot data. This yields the observed effect size as well as estimates of all relevant variance components and mean squares appearing in TBA Eqn. (11.5) and Eqn. (11.7).
* This is the difficult but critical part: make an educated guess regarding the effect-size, $d$, that one is interested in "detecting" (i.e., hoping to reject the NH with probability $1-\beta$). The author prefers the term "anticipated" effect-size to "true" effect-size (the latter implies knowledge of the true difference between the modalities which, as noted earlier, would obviate the need for a pivotal study). 
* Two scenarios are considered below. In the first scenario, the effect-size is assumed equal to that observed in the pilot study, i.e., $d = d_{obs}$. 
* In the second, so-called "best-case" scenario, one assumes that the anticipate value of $d$ is the observed value plus two-sigma of the confidence interval, in the correct direction, of course, i.e., $d=\left | d_{obs} \right |+2\sigma$. Here $\sigma$ is one-fourth the width of the 95% confidence interval for $d_{obs}$. Anticipating more than $2\sigma$  greater than the observed effect-size would be overly optimistic. The width of the CI implies that chances are less than 2.5% that the anticipated value is at or beyond the overly optimistic value. These points will become clearer when example datasets are analyzed below.
*	Calculate statistical power using the distribution implied by Eqn. (11.4), to calculate the probability that a random value of the relevant F-statistic will exceed the critical value, as in §11.3.2.
* If power is below the desired or "target" power, one tries successively larger value of $J$ and / or $K$ until the target power is reached. 


## Formulae for fixed-reader random-case (FRRC) sample size estimation {#roc-sample-size-or-FRRC-sample-size-estimation}
It was shown in TBA §9.8.2 that for fixed-reader analysis the non-centrality parameter is defined by:

\begin{equation}
\Delta=\frac{JK\sigma_{Y;\tau}^2}{\sigma_{Y;\epsilon}^2+\sigma_{Y;\tau RC}^2+J\sigma_{Y;\tau C}^2}
(\#eq:DeltaFRRC1)
\end{equation}

The sampling distribution of the F-statistic under the AH is:

\begin{equation}
F_{AH|R}\equiv \frac{MST}{MSTC}\sim F_{I-1,(I-1)(K-1),\Delta}
(\#eq:SamplingFFRRC1)
\end{equation}

### Formulae for random-reader fixed-case (RRFC) sample size estimation {#roc-sample-size-or-RRFC-sample-size-estimation}
It is shown in TBA §9.9 that for fixed-case analysis the non-centrality parameter is defined by:

\begin{equation}
\Delta=\frac{JK\sigma_{Y;\tau}^2}{\sigma_{Y;\epsilon}^2+\sigma_{Y;\tau RC}^2+K\sigma_{Y;\tau R}^2}
(\#eq:DeltaFRRFC1)
\end{equation}

Under the AH, the test statistic is distributed as a non-central F-distribution as follows:

\begin{equation}
F_{AH|C}\equiv \frac{MST}{MSTR}\sim F_{I-1,(I-1)(J-1),\Delta}
(\#eq:SamplingFRRFC1)
\end{equation}

### Example 1
In the first example the Van Dyke dataset is regarded as a pilot study. Two implementations are shown, a direct application of the relevant formulae, including usage of the mean squares, which in principle can be calculated from the three variance-components. This is then compared to the `RJafroc` implementation. 

Shown first is the "open" implementation. 


```r
alpha <- 0.05;cat("alpha = ", alpha, "\n")
#> alpha =  0.05
rocData <- dataset02 # select Van Dyke dataset
retDbm <- StSignificanceTesting(dataset = rocData, FOM = "Wilcoxon", method = "DBM") 
varYTR <- retDbm$ANOVA$VarCom["VarTR","Estimates"]
varYTC <- retDbm$ANOVA$VarCom["VarTC","Estimates"]
varYEps <- retDbm$ANOVA$VarCom["VarErr","Estimates"]
effectSize <- retDbm$FOMs$trtMeanDiffs["trt0-trt1","Estimate"]
cat("effect size = ", effectSize, "\n")
#> effect size =  -0.043800322

#RRRC
J <- 10; K <- 163
ncp <- (0.5*J*K*(effectSize)^2)/(K*varYTR+max(J*varYTC,0)+varYEps)
MS <- UtilMeanSquares(rocData, FOM = "Wilcoxon", method = "DBM")
ddf <- (MS$msTR+max(MS$msTC-MS$msTRC,0))^2/(MS$msTR^2)*(J-1)
FCrit <- qf(1 - alpha, 1, ddf)
Power <- 1-pf(FCrit, 1, ddf, ncp = ncp)
data.frame("J"= J,  "K" = K, "FCrit" = FCrit, "ddf" = ddf, "ncp" = ncp, "RRRCPower" = Power)
#>    J   K     FCrit       ddf       ncp  RRRCPower
#> 1 10 163 4.1270572 34.334268 8.1269825 0.79111255

#FRRC
J <- 10; K <- 133
ncp <- (0.5*J*K*(effectSize)^2)/(max(J*varYTC,0)+varYEps)
ddf <- (K-1)
FCrit <- qf(1 - alpha, 1, ddf)
Power <- 1-pf(FCrit, 1, ddf, ncp = ncp)
data.frame("J"= J,  "K" = K, "FCrit" = FCrit, "ddf" = ddf, "ncp" = ncp, "RRRCPower" = Power)
#>    J   K    FCrit ddf       ncp  RRRCPower
#> 1 10 133 3.912875 132 7.9873835 0.80111671

#RRFC
J <- 10; K <- 53
ncp <- (0.5*J*K*(effectSize)^2)/(K*varYTR+varYEps)
ddf <- (J-1)
FCrit <- qf(1 - alpha, 1, ddf)
Power <- 1-pf(FCrit, 1, ddf, ncp = ncp)
data.frame("J"= J,  "K" = K, "FCrit" = FCrit, "ddf" = ddf, "ncp" = ncp, "RRRCPower" = Power)
#>    J  K    FCrit ddf       ncp  RRRCPower
#> 1 10 53 5.117355   9 10.048716 0.80496663
```

For 10 readers, the numbers of cases needed for 80% power is largest (163) for RRRC and least for RRFC (53). For all three analyses, the expectation of 80% power is met - the numbers of cases and readers were chosen to achieve close to 80% statistical power. Intermediate quantities such as the critical value of the F-statistic, `ddf` and `ncp` are shown. The reader should confirm that the code does in fact implement the relevant formulae. Shown next is the `RJafroc` implementation. The relevant file is mainSsDbm.R, a listing of which follows: 

### Fixed-reader random-case (FRRC) analysis {#roc-sample-size-or-FRRCAnalysis}

### Random-reader fixed-case (RRFC) analysis {#roc-sample-size-or-RRFCAnalysis}

### Single-treatment multiple-reader analysis {#roc-sample-size-or-STMRAnalysis}

## Discussion/Summary/3

## References {#roc-sample-size-or-references}

