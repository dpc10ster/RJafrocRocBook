# Significance Testing using the DBM Method {#dbm-analysis-significance-testing}








## TBA How much finished {#dbm-analysis-significance-testing-how-much-finished}
60%





## The DBM sampling model
DBM = Dorfman Berbaum Metz

The figure-of-merit has three indices:  

* A treatment index $i$, where $i$ runs from 1 to $I$, where $I$ is the total number of treatments.  
* A reader index $j$, where $j$ runs from 1 to $J$, where $J$ is the total number of readers.  
* The case-sample index $\{c\}$, where  $\{1\}$  i.e., $c$ = 1, denotes a set of cases, $K_1$ non-diseased and $K_2$ diseased, interpreted by all readers in all treatments, and other integer values of $c$ correspond to other independent sets of cases that, although not in fact interpreted by the readers, could potentially be “interpreted” using resampling methods such as the bootstrap or the jackknife. 

The approach [@RN204] taken by DBM was to use the jackknife resampling method to calculate FOM pseudovalues ${Y'}_{ijk}$ defined by (the reason for the prime will become clear shortly):

\begin{equation}
Y'_{ijk}=K\theta_{ij}-(K-1)\theta_{ij(k)}
(\#eq:pseudoValPrime)
\end{equation}

Here $\theta_{ij}$ is the estimate of the figure-of-merit for reader $j$ interpreting all cases in treatment $i$ and $\theta_{ij(k)}$ is the corresponding figure of merit with case $k$ *deleted* from the analysis. To keep the notation compact the case-sample index $\{1\}$ on every figure of merit symbol is suppressed. 

Recall from book Chapter 07 that the jackknife is a way of teasing out the case-dependence: the left hand side of Equation  \@ref(eq:pseudoValPrime) has a case index $k$, with $k$ running from 1 to $K$, where $K$ is the total number of cases:  $K=K_1+K_2$. 

Hillis et al [@RN1866] proposed a centering transformation on the pseudovalues (he terms it "normalized" pseudovalues, but to me "centering" is a more accurate and descriptive term - *Normalize: (In mathematics) multiply (a series, function, or item of data) by a factor that makes the norm or some associated quantity such as an integral equal to a desired value (usually 1). New Oxford American Dictionary, 2016*):

\begin{equation}
Y_{ijk}=Y'_{ijk}+\left (\theta_{ij} - Y'_{ij\bullet}  \right )
(\#eq:pValCentered)
\end{equation}

**Note: the bullet symbol denotes an average over the corresponding index.**

The effect of this transformation is that the average of the centered pseudovalues over the case index is identical to the corresponding estimate of the figure of merit:

\begin{equation}
Y_{ij\bullet}=Y'_{ij\bullet}+\left (\theta_{ij} - Y'_{ij\bullet}  \right )=\theta_{ij}
(\#eq:EffectOfCentering)
\end{equation}

This has the advantage that all confidence intervals are properly centered. The transformation is unnecessary if one uses the Wilcoxon as the figure-of-merit, as the pseudovalues calculated using the Wilcoxon as the figure of merit are "naturally" centered, i.e.,

$\theta_{ij} - Y'_{ij\bullet} = 0$

*It is understood that, unless explicitly stated otherwise, all calculations from now on will use centered pseudovalues.*

Consider $N$ replications of a MRMC study, where a replication means repetition of the study with the same treatments, readers and case-set $\{C=1\}$. For $N$ replications per treatment-reader-case combination, the DBM model for the pseudovalues is ($n$ is the replication index, usually $n$ = 1, but kept here for now):

\begin{equation}
Y_{n(ijk)}  = \mu + \tau_i+ R_j + C_k + (\tau R)_{ij}+ (\tau C)_{ik}+ (R C)_{jk} + (\tau RC)_{ijk}+ \epsilon_{n(ijk)} 
(\#eq:DefDBMModel)
\end{equation}

The term $\mu$ is a constant. By definition, the treatment effect $\tau_i$  is subject to the constraint:

\begin{equation}
\sum_{i=1}^{I}\tau_i=0\Rightarrow \tau_\bullet=0
(\#eq:constraintTau)
\end{equation}

This constraint ensures that $\mu$ has the interpretation of the average of the pseudovalues over treatments, readers and cases. 

The (nesting) notation for the replication index, i.e., $n(ijk)$, implies $n$ observations for treatment-reader-case combination $ijk$. With no replications ($N$ = 1) it is convenient to omit the n-symbol.

The parameter $\tau_i$ is estimated as follows: 

\begin{equation}
Y_{ijk} \equiv Y_{1(ijk)}\\
\tau_i = Y_{i \bullet \bullet} -Y_{\bullet \bullet \bullet} 
(\#eq:estimatingTau)
\end{equation}

*The basic assumption of the DBM model is that the pseudovalues can be regarded as independent and identically distributed observations. That being the case, the pseudovalues can be analyzed by standard ANOVA techniques.* Since pseduovalues are computed from a common dataset, this assumption is, non-intuitive. However, for the special case of Wilcoxon figure of merit, it is justified.

### Explanation of terms in the model 
The right hand side of Eqn. \@ref(eq:pseudoValPrime) consists of one fixed and 7 random effects. The current analysis assumes readers and cases as random factors (RRRC), so by definition $R_j$ and $C_k$ are random effects, and moreover, any term that includes a random factor is a random effect; for example, $(\tau R)_{ij}$ is a random effect because it includes the $R$ factor. Here is a list of the random terms: 

\begin{equation}
R_j, C_k, (\tau R)_{ij}, (\tau C)_{ik}, (RC)_{jk},  (\tau RC)_{ijk},  \epsilon_{ijk}
(\#eq:DBMRandomTerms)
\end{equation}


**Assumption:** Each of the random effects is modeled as a random sample from mutually independent zero-mean normal distributions with variances as specified below:

\begin{align}
\left.\begin{array}{rll}
{R_j}&\sim& N\left ( 0,\sigma_{R}^{2} \right ) \\[0.5em]
{C_k}&\sim& N\left ( 0,\sigma_{C}^{2} \right ) \\[0.5em]
{(\tau R)_{ij}}&\sim& N\left ( 0,\sigma_{\tau R}^{2} \right ) \\[0.5em]
{(\tau C)_{ik}}&\sim& N\left ( 0,\sigma_{\tau C}^{2} \right ) \\[0.5em]
{(RC)_{jk}}&\sim& N\left ( 0,\sigma_{RC}^{2} \right ) \\[0.5em]
{(\tau RC)_{ijk}}&\sim& N\left ( 0,\sigma_{\tau RC}^{2} \right ) \\[0.5em]
\epsilon_{ijk} &\sim& N\left ( 0,\sigma_{\epsilon}^{2} \right )
\end{array}\right\}
(\#eq:samplingOfDbmTerms)
\end{align}

Equation \@ref(eq:samplingOfDbmTerms) defines the meanings of the variance components appearing in Equation \@ref(eq:DBMRandomTerms). One could have placed a $Y$ subscript (or superscript) on each of the variances, as they describe fluctuations of the pseudovalues, not FOM values. However, this tends to clutter the notation. So here is the convention:

**Unless explicitly stated otherwise, all variance symbols in this chapter refer to pseudovalues. ** 
Another convention: $(\tau R)_{ij}$ is *not* the product of the treatment and reader factors, rather it is a single factor, namely the treatment-reader factor with $IJ$ levels, subscripted by the index $ij$ and similarly for the other product-like terms in Equation \@ref(eq:samplingOfDbmTerms).

### Meanings of variance components in the DBM model (**TBA this section can be improved**)
The variances defined in \@ref(eq:samplingOfDbmTerms) are collectively termed *variance components*. Specifically, they are jackknife pseudovalue variance components, to be distinguished from figure of merit (FOM) variance components to be introduced in TBA Chapter 10. They are in order: $\sigma_{R}^{2} ,\sigma_{C}^{2} \sigma_{\tau R}^{2},\sigma_{\tau C}^{2},\sigma_{RC}^{2}, \sigma_{\tau RC}^{2},\sigma_{\epsilon}^{2}$. They have the following meanings.

*	The term $\sigma_{R}^{2}$ is the variance of readers that is independent of treatment or case, which are modeled separately. It is not to be confused with the terms $\sigma_{br+wr}^{2}$ and $\sigma_{cs+wr}^{2}$ used in §9.3, which describe the variability of $\theta$ measured under specified conditions. [A jackknife pseudovalue is a weighted difference of FOM like quantities, TBA \@ref(eq:pseudoValPrime). Its meaning will be explored later. For now, *a pseudovalue variance is distinct from a FOM variance*.]
*	The term $\sigma_{C}^{2}$ is the variance of cases that is independent of treatment or reader.
*	The term $\sigma_{\tau R}^{2}$ is the treatment-dependent variance of readers that was excluded in the definition of $\sigma_{R}^{2}$. If one were to sample readers and treatments for the same case-set, the net variance would be  $\sigma_{R}^{2}+\sigma_{\tau R}^{2}+\sigma_{\epsilon}^{2}$. 
*	The term $\sigma_{\tau C}^{2}$ is the treatment-dependent variance of cases that was excluded in the definition of $\sigma_{C}^{2}$. So, if one were to sample cases and treatments for the same readers, the net variance would be $\sigma_{C}^{2}+\sigma_{\tau C}^{2}+\sigma_{\epsilon}^{2}$. 
*	The term $\sigma_{RC}^{2}$ is the treatment-independent variance of readers and cases that were excluded in the definitions of $\sigma_{R}^{2}$ and $\sigma_{C}^{2}$. So, if one were to sample readers and cases for the same treatment, the net variance would be  $\sigma_{R}^{2}+\sigma_{C}^{2}+\sigma_{RC}^{2}+\sigma_{\epsilon}^{2}$. 
*	The term $\sigma_{\tau RC}^{2}$ is the variance of treatments, readers and cases that were excluded in the definitions of all the preceding terms in TBA \@ref(eq:pseudoValPrime). So, if one were to sample treatments, readers and cases the net variance would be $\sigma_{R}^{2}+\sigma_{C}^{2}+\sigma_{\tau C}^{2}+\sigma_{RC}^{2}+\sigma_{\tau RC}^{2}+\sigma_{\epsilon}^{2}$. 
*	The last term, $\sigma_{\epsilon}^{2}$  describes the variance arising from different replications of the study using the same treatments, readers and cases. Measuring this variance requires repeating the study several ($N$) times with the same treatments, readers and cases, and computing the variance of $Y_{n(ijk)}$ , where the additional $n$-index refers to true replications, $n$ = 1, 2, ..., $N$. 

\begin{equation}
\sigma_{\epsilon}^{2}=\frac{1}{IJK}\sum_{i=1}^{I}\sum_{j=1}^{J}\sum_{k=1}^{k}\frac{1}{N-1}\sum_{n=1}^{N}\left ( Y_{n(ijk)} - Y_{\bullet (ijk)} \right )^2
(\#eq:EstimatingEpsilon)
\end{equation}

The right hand side of TBA \@ref(eq:pseudoValPrime) is the variance of $Y_{n(ijk)}$, for specific $ijk$, with respect to the replication index $n$, averaged over all $ijk$. In practice $N$ = 1 (i.e., there are no replications) and this variance cannot be estimated (it would imply dividing by zero). It has the meaning of *reader inconsistency*, usually termed *within-reader* variability. As will be shown later, the presence of this inestimable term does not limit ones ability to perform significance testing on the treatment effect without having to replicate the whole study, as implied in earlier work [@RN1450].

An equation like TBA \@ref(eq:pseudoValPrime) is termed a *linear model* with the left hand side, the pseudovalue "observations", modeled by a sum of fixed and random terms. Specifically it is a *mixed model*, because the right hand side has both fixed and random effects. Statistical methods have been developed for analysis of such linear models. One estimates the terms on the right hand side of TBA \@ref(eq:pseudoValPrime), it being understood that for the random effects, one estimates the variances of the zero-mean normal distributions, TBA \@ref(eq:pseudoValPrime)Eqn. (9.7), from which the samples are obtained (by assumption).

Estimating the fixed effects is trivial. The term $\mu$ is estimated by averaging the left hand side of TBA \@ref(eq:pseudoValPrime)Eqn. (9.4) over all three indices (since $N$ = 1): $\mu=Y_{\bullet \bullet \bullet}$

Because of the way the treatment effect is defined, TBA \@ref(eq:pseudoValPrime) Eqn. (9.5), averaging, which involves summing, over the treatment-index $i$, yields zero, and all of the remaining random terms yield zero upon averaging, because they are individually sampled from zero-mean normal distributions. To estimate the treatment effect one takes the difference $\tau_i=Y_{\bullet \bullet \bullet}-\mu$.

It can be easily seen that the reader and case averaged difference between two different treatments $i$ and $i'$  is estimated by $\tau_i-\tau_{i'} = Y_{i \bullet \bullet} - Y_{i' \bullet \bullet}$.

Estimating the strengths of the random terms is a little more complicated. It involves methods adapted from least squares, or maximum likelihood, and more esoteric ways. I do not feel comfortable going into these methods. Instead, results are presented and arguments are made to make them plausible. The starting point is definitions of quantities called **mean squares** and their expected values.

### Definitions of mean-squares 
Again, to be clear, one chould put a $Y$ subscript (or superscript) on each of the following definitions, but that would make the notation unnecessarily cumbersome. 

*In this chapter, all mean-square quantities are calculated using pseudovalues, not figure-of-merit values. The presence of three subscripts on Y should make this clear. Also the replication index and the nesting notation are suppressed. The notation is abbreviated so MST is the mean square corresponding to the treatment effect, etc.*

The definitions of the mean-squares below match those (where provided) in [@RN1476, page 1261]. 

\begin{align}
\left.\begin{array}{rll}
\text{MST}&=&\frac{JK\sum_{i=1}^{I}\left ( Y_{i \bullet \bullet} - Y_{ \bullet \bullet \bullet} \right )^2}{I-1}\\[0.5em]
\text{MSR}&=&\frac{IK\sum_{j=1}^{J}\left ( Y_{\bullet j \bullet} - Y_{ \bullet \bullet \bullet} \right )^2}{J-1}\\[0.5em] 
\text{MS(C)}&=&\frac{IJ\sum_{k=1}^{K}\left ( Y_{\bullet \bullet k} - Y_{ \bullet \bullet \bullet} \right )^2}{K-1}\\[0.5em] 
\text{MSTR}&=&\frac{K\sum_{i=1}^{I}\sum_{j=1}^{J}\left ( Y_{i j \bullet} - Y_{i \bullet \bullet} - Y_{\bullet j \bullet} + Y_{ \bullet \bullet \bullet} \right )^2}{(I-1)(J-1)}\\[0.5em] 
\text{MSTC}&=&\frac{J\sum_{i=1}^{I}\sum_{k=1}^{K}\left ( Y_{i \bullet k} - Y_{i \bullet \bullet} - Y_{\bullet \bullet k} + Y_{ \bullet \bullet \bullet} \right )^2}{(I-1)(K-1)}\\[0.5em] 
\text{MSRC}&=&\frac{I\sum_{j=1}^{J}\sum_{k=1}^{K}\left ( Y_{\bullet j k} - Y_{\bullet j \bullet} - Y_{\bullet \bullet k} + Y_{ \bullet \bullet \bullet} \right )^2}{(J-1)(K-1)}\\[0.5em] 
\text{MSTRC}&=&\frac{\sum_{i=1}^{I}\sum_{j=1}^{J}\sum_{k=1}^{K}\left ( Y_{i j k} - Y_{i j \bullet} - Y_{i \bullet k} - Y_{\bullet j k} + Y_{i \bullet \bullet} + Y_{\bullet j \bullet} + Y_{\bullet \bullet k} - Y_{ \bullet \bullet \bullet} \right )^2}{(I-1)(J-1)K-1)}
\end{array}\right\}
(\#eq:MeanSquares)
\end{align}

Note the absence of $MSE$, corresponding to the $\epsilon$ term on the right hand side of \@ref(eq:MeanSquares). With only one observation per treatment-reader-case combination, MSE cannot be estimated; it effectively gets absorbed into the $MSTRC$ term. 

## Expected values of mean squares 
> "In our original formulation [2], expected mean squares for the ANOVA were derived from a restricted parameterization in which mixed-factor interactions sum to zero over indexes of fixed effects. In the restricted parameterization, the mixed effects are correlated, parameters are sometimes awkward to define [17], and extension to unbalanced designs is dubious [17, 18]. In this article, we recommend the unrestricted parameterization. The restricted and unrestricted parameterizations are special cases of a general model by Scheffe [19] that allows an arbitrary covariance structure among experimental units within a level of a random factor. Tables 1 and 2 show the ANOVA tables with expected mean squares for the unrestricted formulation."
>
> --- [@RN2079]

The *observed* mean squares defined in Equation \@ref(eq:MeanSquares) can be calculated directly from the *observed* pseudovalues. The next step in the analysis is to obtain expressions for their *expected* values in terms of the variances defined in \@ref(eq:MeanSquares). Assuming no replications, i.e., $N$ = 1, the expected mean squares are as follows, Table Table \@ref(tab:ExpValMs); understanding how this table is derived, would lead me outside my expertise and the scope of this book; suffice to say that these are *unconstrained* estimates (as summarized in the quotation above) which are different from the *constrained* estimates appearing in the original DBM publication [@RN204]. 

Source          df               E(MS)
-------         ----             ------------------------------------------------------------
T               (I-1)            $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$ + $J\sigma_{\tau C}^{2}$ + $JK\sigma_{\tau}^{2}$ 
R               (J-1)            $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ + $IK\sigma_{R}^{2}$ + $K\sigma_{\tau R}^{2}$
C               (K-1)            $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ + $IJ\sigma_{C}^{2}$ + $J\sigma_{\tau C}^{2}$
TR              (I-1)(J-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$  
TC              (I-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $J\sigma_{\tau C}^{2}$
RC              (J-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ 
TRC             (I-1)(J-1)(K-1)  $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ 
$\epsilon$      $N-1=0$          $\sigma_{\epsilon}^{2}$

Table: (\#tab:ExpValMs) Unconstrained expected values of mean-squares, as in [@RN2079]

* In Table \@ref(tab:ExpValMs) the following notation is used as a shorthand:

\begin{equation}
\sigma_{\tau}^{2}=\frac{1}{I-1}\sum_{i=1}^{I}\left ( Y_{i \bullet \bullet} - Y_{\bullet \bullet \bullet} \right )^2
(\#eq:defnVarTau)
\end{equation}

Since treatment is a fixed effect, the variance symbol $\sigma_{\tau}^{2}$, which is used for notational consistency in Table \@ref(tab:ExpValMs), could cause confusion. The right hand side "looks like" a variance, indeed one that could be calculated for just two treatments but, of course, random sampling from a *distribution of treatments* is not the intent of the notation. 

## Random-reader random-case (RRRC) analysis {#dbm-analysis-significance-testing-RRRC-analysis}
Both readers and cases are regarded as random factors. The expected mean squares in Table Table \@ref(tab:ExpValMs) are variance-like quantities; specifically, they are weighted linear combinations of the variances appearing in \@ref(eq:samplingOfDbmTerms). For single factors the column headed "degrees of freedom" ($df$) is one less than the number of levels of the corresponding factor; estimating a variance requires first estimating the mean, which imposes a constraint, thereby decreasing $df$ by one. For interaction terms, $df$ is the product of the degrees of freedom for the individual factors. As an example, the term $(\tau RC)_{ijk}$ contains three individual factors, and therefore $df = (I-1)(J-1)(K-1)$. The number of degrees of freedom can be thought of as the amount of information available in estimating a mean square. As a special case, with no replications, the $\epsilon$ term has zero $df$ as $N-1 = 0$. With only one observation $Y_{1(ijk)}$ there is no information to estimate the variance corresponding to the $\epsilon$  term. To estimate this term one needs to replicate the study several times – each time the same readers interpret the same cases in all treatments – a very boring task for the reader and totally unnecessary from the researcher's point of view.

### Calculation of mean squares: an example
* We choose `dataset02` to illustrate calculation of mean squares for pseudovalues. This is referred to in the book as the "VD" dataset [@RN1993]. It consists of 114 cases, 45 of which are diseased, interpreted in two treatments by five radiologists using the ROC paradigm. 

* The first line computes the pseudovalues using the `RJafroc` function `UtilPseudoValues()`, and the second line extracts the numbers of treatments, readers and cases. The following lines calculate, using Equation \@ref(eq:MeanSquares) the mean-squares. After displaying the results of the calculation, the results are compared to those calculated by the `RJafroc` function `UtilMeanSquares()`.


```r
Y <- UtilPseudoValues(dataset02, FOM = "Wilcoxon")$jkPseudoValues

I <- dim(Y)[1];J <- dim(Y)[2];K <- dim(Y)[3]

msT <- 0
for (i in 1:I)  {
  msT <- msT + (mean(Y[i, , ]) - mean(Y))^2
}
msT <- msT * J * K/(I - 1)

msR <- 0
for (j in 1:J) {
  msR <- msR + (mean(Y[, j, ]) - mean(Y))^2
}
msR <- msR * I * K/(J - 1)

msC <- 0
for (k in 1:K) {
  msC <- msC + (mean(Y[, , k]) - mean(Y))^2
}
msC <- msC * I * J/(K - 1)

msTR <- 0
for (i in 1:I) {
  for (j in 1:J) {
    msTR <- msTR + 
      (mean(Y[i, j, ]) - mean(Y[i, , ]) - mean(Y[, j, ]) + mean(Y))^2
  }
}
msTR <- msTR * K/((I - 1) * (J - 1))

msTC <- 0
for (i in 1:I) {
  for (k in 1:K) {
    msTC <- msTC + 
      (mean(Y[i, , k]) - mean(Y[i, , ]) - mean(Y[, , k]) + mean(Y))^2
  }
  msTC <- msTC * J/((I - 1) * (K - 1))
} 

msTC <- 0
for (i in 1:I) {
  for (k in 1:K) { # OK
    msTC <- msTC + 
      (mean(Y[i, , k]) - mean(Y[i, , ]) - mean(Y[, , k]) + mean(Y))^2
  }
}
msTC <- msTC * J/((I - 1) * (K - 1))

msRC <- 0
for (j in 1:J) {
  for (k in 1:K) {
    msRC <- msRC + 
      (mean(Y[, j, k]) - mean(Y[, j, ]) - mean(Y[, , k]) + mean(Y))^2
  }
}
msRC <- msRC * I/((J - 1) * (K - 1))

msTRC <- 0
for (i in 1:I) {
  for (j in 1:J) {
    for (k in 1:K) {
      msTRC <- msTRC + (Y[i, j, k] - mean(Y[i, j, ]) - 
                          mean(Y[i, , k]) - mean(Y[, j, k]) + 
                          mean(Y[i, , ]) + mean(Y[, j, ]) + 
                          mean(Y[, , k]) - mean(Y))^2
    }
  }
}
msTRC <- msTRC/((I - 1) * (J - 1) * (K - 1))

data.frame("msT" = msT, "msR" = msR, "msC" = msC, 
           "msTR" = msTR, "msTC" = msTC, 
           "msRC" = msRC, "msTRC" = msTRC)
#>         msT       msR       msC       msTR       msTC       msRC     msTRC
#> 1 0.5467634 0.4373268 0.3968699 0.06281749 0.09984808 0.06450106 0.0399716

as.data.frame(UtilMeanSquares(dataset02)[1:7])
#>         msT       msR       msC       msTR       msTC       msRC     msTRC
#> 1 0.5467634 0.4373268 0.3968699 0.06281749 0.09984808 0.06450106 0.0399716
```

### Significance testing {#dbm-analysis-significance-testing-sig-testing}
If the NH of no treatment effect is true, i.e., if $\sigma_{\tau}^{2}$ = 0, then according to Table \@ref(tab:ExpValMs) the following holds (the last term in the row labeled $T$ in Table \@ref(tab:ExpValMs) drops out):

\begin{equation}
E\left ( MST\mid NH \right ) = \sigma_{\epsilon}^{2} + \sigma_{\tau RC}^{2} + K\sigma_{\tau R}^{2} + J\sigma_{\tau C}^{2}
(\#eq:ExpMST)
\end{equation}

Also, the following linear combination is equal to $E\left ( MST\mid NH \right )$:

\begin{align}
\begin{split}
&E\left ( MSTR \right ) + E\left ( MSTC \right )  - E\left ( MSTRC \right ) \\
&= \left (\sigma_{\epsilon}^{2} + \sigma_{\tau RC}^{2} + K\sigma_{\tau R}^{2} \right ) + \left (\sigma_{\epsilon}^{2} + \sigma_{\tau RC}^{2} + J\sigma_{\tau C}^{2} \right ) -\left (\sigma_{\epsilon}^{2} + \sigma_{\tau RC}^{2}  \right ) \\
&= \sigma_{\epsilon}^{2} + \sigma_{\tau RC}^{2} + J \sigma_{\tau C}^{2} +  K\sigma_{\tau R}^{2} \\
&= E\left ( MST\mid NH \right )
\end{split}
(\#eq:linearComb)
\end{align}

Therefore, under the NH, the ratio: 

\begin{equation}
\frac{E\left ( MST\mid NH \right )}{E\left ( MSTR \right ) + E\left ( MSTC \right )  - E\left ( MSTRC \right )} = 1
(\#eq:ConstrFRatio)
\end{equation}

In practice, one does not know the expected values – that would require averaging each of these quantities, regarded as random variables, over their respective distributions. Therefore, one defines the following statistic, denoted  $F_{DBM}$, using the observed values of the mean squares, calculated almost trivially as in the previous example, using their definitions in Equation \@ref(eq:MeanSquares):

\begin{equation}
F_{DBM} = \frac{MST}{MSTR + MSTC - MSTRC}
(\#eq:DefFStatRRRC)
\end{equation}

$F_{DBM}$ is a realization of a random variable. A non-zero treatment effect, i.e.,  $\sigma_{\tau}^{2} > 0$, will cause the ratio to be larger than one, because $E\left ( MST \right)$  will be larger, see row labeled $T$ in Table \@ref(tab:ExpValMs). Therefore values of $F_{DBM} > 1$  will tend to reject the NH. Drawing on a theorem from statistics [@RN1492], under the NH the ratio of two independent mean squares is distributed as a (central) F-statistic with degrees of freedom corresponding to those of the mean squares forming the numerator and denominator of the ratio (Theorem 12.2.5 in “An Introduction to Mathematical Statistics and Its Applications”). To perform  hypothesis testing one needs the distribution,  under the NH, of the statistic defined by Eqn. \@ref(eq:DefFStatRRRC). This is completely analogous to Chapter 08 where knowledge of the distribution of AUC under the NH enabled testing the null hypothesis that the observed value of AUC equals a pre-specified value. 

Under the NH, $F_{DBM|NH}$ is distributed according to the F-distribution characterized by two numbers: 

* A numerator degrees of freedom ($\text{ndf}$) – determined by the degrees of freedom of the numerator, $MST$, of the ratio comprising the F-statistic, i.e., $I – 1$, and 
* A denominator degrees of freedom ($\text{ddf}$) - determined by the degrees of freedom of the denominator, $MSTR + MSTC - MSTRC$, of the ratio comprising the F-statistic, to be described in the next section. 

Summarizing,

\begin{align}
\left.\begin{array}{rll}
F_{DBM|NH} \sim F_{\text{ndf},\text{ddf}} \\
\text{ndf}=I-1
\end{array}\right\}
(\#eq:NhDistrDBMDefFStatRRRC)
\end{align}

The next topic is estimating $ddf$.

### The Satterthwaite approximation
The denominator of the F-ratio is $MSTR+MSTC-MSTRC$. This is not a *simple* mean square (I am using terminology in the Satterthwaite papers - he means any mean square defined by equations such as in Equation \@ref(eq:MeanSquares)). Rather it is a *linear combination of mean squares* (with coefficients 1, 1 and -1), and the resulting value could even be negative leading to a negative $F_{DBM|NH}$, which is an illegal value for a sample from an F-distribution (a ratio of two variances). In 1941 Satterthwaite [@RN2359; @RN2360] proposed an approximate degree of freedom for a linear combination of simple mean square quantities. TBA Online Appendix 9.A explains the approximation in more detail. The end result is that the mean square quantity described in Equation \@ref(eq:DefFStatRRRC) has an approximate degree of freedom defined by (this is called the *Satterthwaite's approximation*):

\begin{equation}
ddf_{Sat}=\frac{\left ( MSTR + MSTC - MSTRC \right )^2}{\left ( \frac{MSTR^2}{(I-1)(J-1)} + \frac{MSTC^2}{(I-1)(K-1)} + \frac{MSTRC^2}{(I-1)(J-1)(K-1)}  \right )}
(\#eq:ddfSatter)
\end{equation}

The subscript $Sat$ is for Satterthwaite. From Equation \@ref(eq:ddfSatter) it should be fairly obvious that in general $ddf_{Sat}$ is not an integer. To accommodate possible negative estimates of the denominator of Equation \@ref(eq:ddfSatter), the original DBM method [@RN204] proposed, depending on the signs of $\sigma_{\tau R}^2$ and $\sigma_{\tau C}^2$, four expressions for the F-statistic and corresponding expressions for $ddf$. Rather than repeat them here, since they have been superseded by the method described below, the interested reader is referred to Eqn. 6 and Eqn. 7 in Reference [@RN1866].

Instead Hillis [@RN1865] proposed the following statistic for testing the null hypothesis: 

\begin{equation}
F_{DBM} = \frac{MST}{MSTR + \max \left (MSTC - MSTRC, 0  \right )}
(\#eq:FStatHillis)
\end{equation}

Now the denominator cannot be negative. One can think of the F-statistic $F_{DBM}$ as a signal-to-noise ratio like quantity, with the difference that both numerator and denominator are variance like quantities. If the "variance" represented by the treatment effect is larger than the variance of the noise tending to mask the treatment effect, then $F_{DBM}$  tends to be large, which makes the observed treatment "variance" stand out more clearly compared to the noise, and the NH is more likely to be rejected. Hillis in [@RN1772] has shown that the left hand side of Equation \@ref(eq:FStatHillis) is distributed as an F-statistic with $\text{ndf} = I-1$ and denominator degrees of freedom $ddf_H$ defined by:

\begin{equation}
ddf_H =\frac{\left ( MSTR + \max \left (MSTC - MSTRC,0  \right ) \right )^2}{\text{MSTR}^2}(I-1)(J-1)
(\#eq:DefddfH)
\end{equation}

Summarizing,

\begin{equation}
F_{DBM} \sim F_{\text{ndf},\text{ddf}_H} \\
\text{ndf}=I-1
(\#eq:NhDistrDBMDefFStatRRRCHillis)
\end{equation}

Instead of 4 rules, as in the original DBM method, the Hillis modification involves just one rule, summarized by Equations \@ref(eq:DefddfH) through \@ref(eq:NhDistrDBMDefFStatRRRCHillis). Moreover, the F-statistic is constrained to non-negative values. Using simulation testing [@RN1866] he has been shown that the modified DBM method has better null hypothesis behavior than the original DBM method. The latter tended to be too conservative, typically yielding Type I error rates smaller than the expected 5% for $\alpha$ = 0.05. 

### Decision rules, p-value and confidence intervals
The *critical* value of the F-distribution, denoted $F_{1-\alpha,\text{ndf},\text{ddf}_H}$, is defined such that fraction $1-\alpha$ of the distribution lies to the left of the critical value, in other words it is the $1-\alpha$ *quantile* of the F-distribution: 

\begin{equation}
\Pr\left ( F\leq F_{1-\alpha,\text{ndf},\text{ddf}_H} \mid F\sim F_{\text{ndf},\text{ddf}_H}\right ) = 1 - \alpha
(\#eq:critValFStat)
\end{equation}

The critical value $F_{1-\alpha,\text{ndf},\text{ddf}_H}$ increases as $\alpha$ decreases. The value of $\alpha$, generally chosen to be 0.05, termed the *nominal* $\alpha$, is fixed. The decision rule is that if $F_{DBM} > F_{1-\alpha, \text{ndf}, \text{ddf}_H}$ one rejects the NH and otherwise one does not. It follows, from the definition of $F_{DBM}$, Equation \@ref(eq:FStatHillis), that rejection of the NH is more likely to occur if:

* $F_{DBM}$ is large, which occurs if $MST$ is large, meaning the treatment effect is large
* $MSTR + \max \left (MSTC - MSTRC,0  \right )$ is small, see comments following TBA \@ref(eq:pseudoValPrime) Eqn. (9.23). 
*	$\alpha$ is large: for then $F_{1-\alpha,\text{ndf},\text{ddf}_H}$ decreases and is more likely to be exceeded by the observed value of $F_{DBM}$. 
*	\text{ndf} is large: the more the number of treatment pairings, the greater the chance that at least one pairing will reject the NH. This is one reason sample size calculations are rarely conducted for more than 2-treatments. 
* $\text{ddf}_H$	is large: this causes the critical value to decrease, see below, and is more likely to be exceeded by the observed value of $F_{DBM}$. 

#### p-value of the F-test 
**The p-value of the test is the probability, under the NH, that an equal or larger value of the F-statistic than observed $F_{DBM}$ could occur by chance.** In other words, it is the area under the (central) F-distribution $F_{\text{ndf},\text{ddf}}$ that lies to the right of the observed value of $F_{DBM}$:

\begin{equation}
p=\Pr\left ( F > F_{DBM} \mid F \sim F_{\text{ndf},\text{ddf}_H} \right )
(\#eq:pValueRRRC)
\end{equation}

#### Confidence intervals for inter-treatment FOM differences
If $p < \alpha$ then the NH that all treatments are identical is rejected at significance level $\alpha$. That informs the researcher that there exists at least one treatment-pair that has a difference significantly different from zero. To identify which pair(s) are different, one calculates confidence intervals for each paired difference. Hillis in [@RN1772] has shown that the $(1-\alpha)$ confidence interval for $Y_{i \bullet \bullet} - Y_{i' \bullet \bullet}$ is given by:

\begin{equation}
CI_{1-\alpha}=\left ( Y_{i \bullet \bullet} - Y_{i' \bullet \bullet} \right ) \pm t_{\alpha/2;\text{ddf}_H} \sqrt{\frac{2}{JK}\left ( MSTR + \max\left ( MSTC-MSTRC,0 \right ) \right )}
(\#eq:CIRRRCDBM)
\end{equation}

Here $t_{\alpha/2;\text{ddf}_H}$ is that value such that $\alpha/2$  of the *central t-distribution* with  $\text{ddf}_H$ degrees of freedom is contained in the upper tail of the distribution: 

\begin{equation}
\Pr\left ( T>t_{\alpha/2;\text{ddf}_H} \right )=\alpha/2
(\#eq:tDistTailProb)
\end{equation}

Since centered pseudovalues were used:   

\begin{equation}
\left ( Y_{i \bullet \bullet} - Y_{i' \bullet \bullet} \right )=\left ( \theta_{i \bullet } - \theta_{i' \bullet} \right )
\end{equation}

Therefore, Equation \@ref(eq:CIRRRCDBM) can be rewritten:

\begin{equation}
CI_{1-\alpha}=\left ( \theta_{i \bullet} - \theta_{i' \bullet} \right ) \pm t_{\alpha/2;\text{ddf}_H} \sqrt{\frac{2}{JK}\left ( MSTR + \max\left ( MSTC-MSTRC,0 \right ) \right )}
(\#eq:confIntervalRRRC)
\end{equation}

For two treatments any of the following equivalent rules could be adopted to reject the NH: 

* $F_{DBM} > F_{1-\alpha,\text{ndf},\text{ddf}_H}$
* $p < \alpha$
* $CI_{1-alpha}$ excludes zero

For more than two treatments the first two rules are equivalent and if a significant difference is found using either of them, then one can use the confidence intervals to determine which treatment pair differences are significantly different from zero. The first F-test is called the *overall F-test* and the subsequent tests the *treatment-pair t-tests*. One only conducts treatment pair t-tests if the overall F-test yields a significant result.

#### Code illustrating the F-statistic, ddf and p-value for RRRC analysis, Van Dyke data
Line 1 defines $\alpha$. Line 2 forms a data frame from previously calculated mean-squares. Line 3 calculates the denominator appearing in Equation \@ref(eq:FStatHillis). Line 4 computes the observed value of $F_{DBM}$, namely the ratio of the numerator and denominator in Equation \@ref(eq:FStatHillis). Line 5 sets $\text{ndf}$ to $I - 1$. Line 6 computes $\text{ddf}_H$. Line 7 computes the critical value of the F-distribution $F_{crit}\equiv F_{\text{ndf},\text{ddf}_H}$. Line 8 calculates the p-value, using the definition Equation \@ref(eq:pValueRRRC). Line 9 prints out the just calculated quantities. The next line uses the `RJafroc` function `StSignificanceTesting()` and the 2nd last line prints out corresponding `RJafroc`-computed quantities. Note the correspondences between the values just computed and those provide by `RJafroc`. Note that the FOM difference is not significant at the 5% level of significance as $p > \alpha$. The last line shows that $F_{DBM}$ does not exceed $F_{crit}$. The two rules are equivalent.


```r
alpha <- 0.05
retMS <- data.frame("msT" = msT, "msR" = msR, "msC" = msC, 
                    "msTR" = msTR, "msTC" = msTC, 
                    "msRC" = msRC, "msTRC" = msTRC)
F_DBM_den <- retMS$msTR+max(retMS$msTC - retMS$msTRC,0) 
F_DBM <- retMS$msT / F_DBM_den 
ndf <- (I-1)
ddf_H <- (F_DBM_den^2/retMS$msTR^2)*(I-1)*(J-1)
FCrit <- qf(1 - alpha, ndf, ddf_H)
pValueH <- 1 - pf(F_DBM, ndf, ddf_H)
data.frame("F_DBM" = F_DBM, "ddf_H"= ddf_H, "pValueH" = pValueH) # Line 9
#>      F_DBM    ddf_H    pValueH
#> 1 4.456319 15.25967 0.05166569
retRJafroc <- StSignificanceTesting(dataset02, 
                                    FOM = "Wilcoxon",
                                    method = "DBM")
data.frame("F_DBM" = retRJafroc$RRRC$FTests$FStat[1], 
           "ddf_H"= retRJafroc$RRRC$FTests$DF[2], 
           "pValueH" = retRJafroc$RRRC$FTests$p[1])
#>       F_DBM     ddf_H     pValueH
#> 1 4.4563187 15.259675 0.051665686
F_DBM > FCrit
#> [1] FALSE
```

#### Code illustrating the inter-treatment confidence interval for RRRC analysis, Van Dyke data
Line 1 computes the FOM matrix using function `UtilFigureOfMerit`. The next 9 lines compute the treatment FOM differences. The next line `nDiffs` (for "number of differences") evaluates to 1, as with two treatments, there is only one difference. The next line initializes `CI_DIFF_FOM_RRRC`, which stands for "confidence intervals, FOM differences, for RRRC analysis". The next 8 lines evaluate, using Equation \@ref(eq:confIntervalRRRC), and prints the lower value, the mid-point and the upper value of the confidence interval. Finally, these values are compared to those yielded by `RJafroc`. The FOM difference is not significant, whether viewed from the point of view of the F-statistic not exceeding the critical value, the observed p-value being larger than alpha or the 95% CI for the FOM difference including zero. 

```r
theta <- as.matrix(UtilFigureOfMerit(dataset02, FOM = "Wilcoxon"))
theta_i_dot <- array(dim = I)
for (i in 1:I) theta_i_dot[i] <- mean(theta[i,])
trtDiff <- array(dim = c(I,I))
for (i1 in 1:(I-1)) {    
  for (i2 in (i1+1):I) {
    trtDiff[i1,i2] <- theta_i_dot[i1]- theta_i_dot[i2]    
  }
}
trtDiff <- trtDiff[!is.na(trtDiff)]
nDiffs <- I*(I-1)/2
CI_DIFF_FOM_RRRC <- array(dim = c(nDiffs, 3))
for (i in 1 : nDiffs) {
  CI_DIFF_FOM_RRRC[i,1] <- qt(alpha/2,df = ddf_H)*sqrt(2*F_DBM_den/J/K) + trtDiff[i]
  CI_DIFF_FOM_RRRC[i,2] <- trtDiff[i]
  CI_DIFF_FOM_RRRC[i,3] <- qt(1-alpha/2,df = ddf_H)*sqrt(2*F_DBM_den/J/K) + trtDiff[i]
  print(data.frame("Lower" = CI_DIFF_FOM_RRRC[i,1], 
                   "Mid" = CI_DIFF_FOM_RRRC[i,2], 
                   "Upper" = CI_DIFF_FOM_RRRC[i,3]))
}
#>          Lower          Mid         Upper
#> 1 -0.087959499 -0.043800322 0.00035885444
data.frame("Lower" = retRJafroc$RRRC$ciDiffTrt[1,"CILower"], 
           "Mid" = retRJafroc$RRRC$ciDiffTrt[1,"Estimate"], 
           "Upper" = retRJafroc$RRRC$ciDiffTrt[1,"CIUpper"])
#>          Lower          Mid         Upper
#> 1 -0.087959499 -0.043800322 0.00035885444
```


## Sample size estimation for random-reader random-case generalization
### The non-centrality parameter
In the significance-testing procedure just described, the relevant distribution was that of the F-statistic when the NH is true, Equation \@ref(eq:NhDistrDBMDefFStatRRRCHillis). *For sample size estimation, one needs to know the distribution of the statistic when the NH is false.* In the latter condition (i.e., the AH) the observed F-statistic, defined by Equation \@ref(eq:DefFStatRRRC), is distributed as a *non-central* F-distribution  $F_{\text{ndf},{\text{ddf}}_H,\Delta}$ with *non-centrality parameter* $\Delta$: 

\begin{equation}
F_{DBM|AH} \sim F_{\text{ndf},ddf_H,\Delta}
(\#eq:AhDistrDBMFStatRRRC)
\end{equation}

The non-centrality parameter $\Delta$ is defined, compare [@RN1476] Eqn. 6, by:

\begin{equation*}
\Delta=\frac{JK\sigma_{\tau}^2}{\sigma_{\epsilon}^2+K\sigma_{\tau R}^2+J\sigma_{\tau C}^2}
\end{equation*}

The parameters $\sigma_{\tau}^2$, $\sigma_{\tau R}^2$ and $\sigma_{\tau C}^2$ appearing in this equation are identical to three of the six variances describing the DBM model, Equation \@ref(eq:DefDBMModel). The estimates of $\sigma_{\tau R}^2$ and/or $\sigma_{\tau C}^2$ can turn out to be negative (if either of these parameters is close to zero, an estimate from a small pilot study can be negative). To avoid a possibly negative denominator, [@RN1476] suggest the following modifications (see sentence following Eqn. 4 in cited paper):

\begin{equation}
\Delta=\frac{JK\sigma_{\tau}^2}{\sigma_{\epsilon}^2+\max(K\sigma_{\tau R}^2,0)+\max(J\sigma_{\tau C}^2,0)}
(\#eq:DefDelta)
\end{equation}

The observed effect size $d$, a realization of a random variable, is defined by (the bullet represents an average over the reader index):

\begin{equation}
d=\theta_{1\bullet}-\theta_{2\bullet}
(\#eq:EffectSize)
\end{equation}

For two treatments, since the individual treatment effects must be the negatives of each other (because they sum to zero, see \@ref(eq:constraintTau)), it follows that:

\begin{equation}
\sigma_{\tau}^2=\frac{d^2}{2}
(\#eq:sigma2Tau)
\end{equation}
 
Therefore, for two treatments the numerator of the expression for $\Delta$ is $JKd^2/2$. Dividing numerator and denominator of Equation \@ref(eq:DefDelta) by $K$, one gets the final expression for $\Delta$, as coded in `RJafroc`, namely:

\begin{equation}
\Delta=\frac{Jd^2/2}{\max(\sigma_{\tau R}^2,0)+(\sigma_{\epsilon}^2+\max(J\sigma_{\tau C}^2,0))/K}
(\#eq:DeltaCoded)
\end{equation}

The variances, $\sigma_{\tau}^2$, $\sigma_{\tau R}^2$ and $\sigma_{\tau C}^2$, appearing in Equation \@ref(eq:DeltaCoded), can be calculated from the observed mean squares using the following equations, see [@RN1476] Eqn. 4, 

\begin{equation}
\left.\begin{array}{rl}
\sigma_{\epsilon}^2&={\text{MSTRC}}^*\\[1em] 
\sigma_{\tau R}^2&=\displaystyle\frac{{\text{MSTR}}^*-{\text{MSTRC}}^*}{K^*}\\[1em]
\sigma_{\tau C}^2&=\displaystyle\frac{{\text{MSTC}}^*-{\text{MSTRC}}^*}{J^*}
\end{array}\right\}
(\#eq:HBEqn4Pilot)
\end{equation}

* Here the asterisk is used to (consistently) denote quantities, including the mean squares, pertaining to the *pilot* study. 
* In particular, $J^*$ and $K^*$ denote the numbers of readers and cases, respectively, *in the pilot study*, while $J$ and $K$, appearing elsewhere, for example in Equation \@ref(eq:DeltaCoded), are the corresponding numbers for the *planned or pivotal study*. 
* The three variances, determined from the pilot study via Equation \@ref(eq:HBEqn4Pilot), are assumed to apply unchanged to the pivotal study (as they are sample-size independent parameters of the DBM model). 

### The denominator degrees of freedom
* (The numerator degrees of freedom of the non-central $F$ distribution is always unity.) It remains to calculate the appropriate denominator degrees of freedom for the pivotal study. This is denoted $df_2$, to distinguish it from $ddf_H$, where the latter applies to the pilot study as in Equation \@ref(eq:DefddfH). 
* The starting point is Equation \@ref(eq:DefddfH) with the left hand side replaced by $df_2$, and with the emphasis that *all quantities appearing in it apply to the pivotal study*. 
* The mean squares appearing in Equation \@ref(eq:DefddfH) can be related to the variances by an equation analogous to Equation \@ref(eq:HBEqn4Pilot), except that, again, all quantities in it apply to the *pivotal* study (note the absence of asterisks): 

\begin{equation}
\left.\begin{array}{rl}
\sigma_{\epsilon}^2&=MSTRC\\[1em]
\sigma_{\tau R}^2&=\displaystyle\frac{MSTR-MSTRC}{K}\\[1em]
\sigma_{\tau C}^2&=\displaystyle\frac{MSTC-MSTRC}{J}
\end{array}\right\}
(\#eq:HBEqn4Pivotal)
\end{equation}

Substituting from Equation \@ref(eq:HBEqn4Pivotal) into Equation \@ref(eq:DefddfH) with the left hand side replaced by $df_2$, and dividing numerator and denominator by $K^2$, one has the final expression as coded in `RJafroc`:

\begin{equation}
df_2 = \frac{(\max(\sigma_{\tau R}^2,0)+(\max(J\sigma_{\tau C}^2,0)+\sigma_{\epsilon}^2)/K)^2}{(\max(\sigma_{\tau R}^2,0)+\sigma_{\epsilon}^2/K)^2}(J-1)
(\#eq:ddfCoded)
\end{equation}

### Example of sample size estimation, RRRC generalization
The Van Dyke dataset is regarded as a pilot study. In the first block of code function `StSignificanceTesting()` is used to get the DBM variances (i.e., `VarTR` = $\sigma_{\tau R}^2$, etc.) and the effect size $d$. 


```r
rocData <- dataset02 # select Van Dyke dataset
retDbm <- StSignificanceTesting(dataset = rocData,
                                FOM = "Wilcoxon",
                                method = "DBM") 
VarTR <- retDbm$ANOVA$VarCom["VarTR","Estimates"]
VarTC <- retDbm$ANOVA$VarCom["VarTC","Estimates"]
VarErr <- retDbm$ANOVA$VarCom["VarErr","Estimates"]
d <- retDbm$FOMs$trtMeanDiffs["trt0-trt1","Estimate"]
```

The observed effect size is -0.04380032. The sign is negative as the reader-averaged second modality has greater FOM than the first. The next code block shows implementation of the RRRC formulae just presented. The values of $J$ and $K$ were preselected to achieve 80% power, as verified from the final line of the output. 


```r
#RRRC
J <- 10; K <- 163
den <- max(VarTR, 0) + (VarErr + J * max(VarTC, 0)) / K 
deltaRRRC <- (d^2 * J/2) / den
df2 <- den^2 * (J - 1) / (max(VarTR, 0) + VarErr / K)^2
fvalueRRRC <- qf(1 - alpha, 1, df2)
Power <- 1-pf(fvalueRRRC, 1, df2, ncp = deltaRRRC)
data.frame("J"= J,  "K" = K, "fvalueRRRC" = fvalueRRRC, "df2" = df2, "deltaRRRC" = deltaRRRC, "PowerRRRC" = Power)
#>    J   K fvalueRRRC       df2 deltaRRRC  PowerRRRC
#> 1 10 163  3.9930236 63.137871 8.1269825 0.80156249
```

## Significance testing and sample size estimation for fixed-reader random-case generalization
The extension to FRRC generalization is as follows. One sets $\sigma_R^2 = 0$ and $\sigma_{\tau R}^2 = 0$ in the DBM model \@ref(eq:DefDBMModel). The F-statistic for testing the NH and its distribution under the NH is:

\begin{equation}
F=\frac{\text{MST}}{\text{MSTC}} \sim F_{I-1,(I-1)(K-1)}
(\#eq:DefFStatFRRC-DBM1)
\end{equation}

The NH is rejected if the observed value of $F$ exceeds the critical value defined by $F_{\alpha, I-1,(I-1)(K-1)}$. For two modalities the denominator degrees of freedom is $df_2 = K-1$. The expression for the non-centrality parameter follows from \@ref(eq:DeltaCoded) upon setting $\sigma_{\tau R}^2 = 0$.

\begin{equation}
\Delta=\frac{Jd^2/2}{(\sigma_{\epsilon}^2+\max(J\sigma_{\tau C}^2,0))/K}
(\#eq:DefDeltaFRRC)
\end{equation}

These equations are coded in the following code-chunk:


```r
#FRRC
# set VarTC = 0 in RRRC formulae
J <- 10; K <- 133
den <- (VarErr + J * max(VarTC, 0)) / K
deltaFRRC <- (d^2 * J/2) / den
df2FRRC <- K - 1
fvalueFRRC <- qf(1 - alpha, 1, df2FRRC)
powerFRRC <- pf(fvalueFRRC, 1, df2FRRC, ncp = deltaFRRC, FALSE)
data.frame("J"= J,  "K" = K, "fvalueFRRC" = fvalueFRRC, "df2" = df2FRRC, "deltaFRRC" = deltaFRRC, "powerFRRC" = powerFRRC)
#>    J   K fvalueFRRC df2 deltaFRRC  powerFRRC
#> 1 10 133   3.912875 132 7.9873835 0.80111671
```

## Significance testing and sample size estimation for random-reader fixed-case generalization
The extension to RRFC generalization is as follows. One sets $\sigma_C^2 = 0$ and $\sigma_{\tau C}^2 = 0$ in the DBM model \@ref(eq:DefDBMModel). The F-statistic for testing the NH and its distribution under the NH is:

\begin{equation}
F=\frac{\text{MST}}{\text{MSTR}} \sim F_{I-1,(I-1)(J-1)}
(\#eq:DefFStatFRRC-DBM2)
\end{equation}

The NH is rejected if the observed value of $F$ exceeds the critical value defined by $F_{\alpha, I-1,(I-1)(J-1)}$. For two modalities the denominator degrees of freedom is $df_2 = J-1$. The expression for the non-centrality parameter follows from \@ref(eq:DeltaCoded) upon setting $\sigma_{\tau C}^2 = 0$.

\begin{equation}
\Delta=\frac{Jd^2/2}{\max(\sigma_{\tau R}^2,0)+\sigma_{\epsilon}^2/K}
(\#eq:DefDeltaRRFC)
\end{equation}

These equations are coded in the following code-chunk:


```r
#RRFC
# set VarTR = 0 in RRRC formulae
J <- 10; K <- 53
den <- max(VarTR, 0) + VarErr/K
deltaRRFC <- (d^2 * J/2) / den
df2RRFC <- J - 1
fvalueRRFC <- qf(1 - alpha, 1, df2RRFC)
powerRRFC <- pf(fvalueRRFC, 1, df2RRFC, ncp = deltaRRFC, FALSE)
data.frame("J"= J,  "K" = K, "fvalueRRFC" = fvalueRRFC, "df2" = df2RRFC, "deltaRRFC" = deltaRRFC, "powerRRFC" = powerRRFC)
#>    J  K fvalueRRFC df2 deltaRRFC  powerRRFC
#> 1 10 53   5.117355   9 10.048716 0.80496663
```

It is evident that for this dataset, for 10 readers, the numbers of cases needed for 80% power is largest (163) for RRRC and least for RRFC (53). For all three analyses, the expectation of 80% power is met - the numbers of cases and readers were deliberately chosen to achieve close to 80% statistical power.  

## Summary TBA {#dbm-analysis-significance-testing-summary}
This chapter has detailed analysis of MRMC ROC data using the DBM method. A reason for the level of detail is that almost all of the material carries over to other data collection paradigms, and a thorough understanding of the relatively simple ROC paradigm data is helpful to understanding the more complex ones. 

DBM has been used in several hundred ROC studies (Prof. Kevin Berbaum, private communication ca. 2010). While the method allows generalization of a study finding, e.g., rejection of the NH, to the population of readers and cases, I believe this is sometimes taken too literally. If a study is done at a single hospital, then the radiologists tend to be more homogenous as compared to sampling radiologists from different hospitals. This is because close interactions between radiologists at a hospital tend to homogenize reading styles and performance. A similar issue applies to patient characteristics, which are also expected to vary more between different geographical locations than within a given location served by the hospital. This means is that single hospital study based p-values may tend to be biased downwards, declaring differences that may not be replicable if a wider sampling "net" were used using the same sample size. The price paid for a wider sampling net is that one must use more readers and cases to achieve the same sensitivity to genuine treatment effects, i.e., statistical power (i.e., there is no "free-lunch").

A third MRMC ROC method, due to Clarkson, Kupinski and Barrett19,20, implemented in open-source JAVA software by Gallas and colleagues22,44 (http://didsr.github.io/iMRMC/) is available on the web. Clarkson et al19,20 provide a probabilistic rationale for the DBM model, provided the figure of merit is the empirical $AUC$. The method is elegant but it is only applicable as long as one is using the empirical AUC as the figure of merit (FOM) for quantifying observer performance. In contrast the DBM approach outlined in this chapter, and the approach outlined in the following chapter, are applicable to any scalar FOM. Broader applicability ensures that significance-testing methods described in this, and the following chapter, apply to other ROC FOMs, such as binormal model or other fitted AUCs, and more importantly, to other observer performance paradigms, such as free-response ROC paradigm. An advantage of the Clarkson et al. approach is that it predicts truth-state dependence of the variance components. One knows from modeling ROC data that diseased cases tend to have greater variance than non-diseased ones, and there is no reason to suspect that similar differences do not exist between the variance components. 

Testing validity of an analysis method via simulation testing is only as good as the simulator used to generate the datasets, and this is where current research is at a bottleneck. The simulator plays a central role in ROC analysis. In my opinion this is not widely appreciated. In contrast, simulators are taken very seriously in other disciplines, such as cosmology, high-energy physics and weather forecasting. The simulator used to validate3 DBM is that proposed by Roe and Metz39 in 1997. This simulator has several shortcomings. (a) It assumes that the ratings are distributed like an equal-variance binormal model, which is not true for most clinical datasets (recall that the b-parameter of the binormal model is usually less than one). Work extending this simulator to unequal variance has been published3. (b) It does not take into account that some lesions are not visible, which is the basis of the contaminated binormal model (CBM). A CBM model based simulator would use equal variance distributions with the difference that the distribution for diseased cases would be a mixture distribution with two peaks. The radiological search model (RSM) of free-response data, Chapter 16 &17 also implies a mixture distribution for diseased cases, and it goes farther, as it predicts some cases yield no z-samples, which means they will always be rated in the lowest bin no matter how low the reporting threshold. Both CBM and RSM account for truth dependence by accounting for the underlying perceptual process. (c) The Roe-Metz simulator is out dated; the parameter values are based on datasets then available (prior to 1997). Medical imaging technology has changed substantially in the intervening decades. d Finally, the methodology used to arrive at the proposed parameter values is not clearly described. Needed is a more realistic simulator, incorporating knowledge from alternative ROC models and paradigms that is calibrated, by a clearly defined method, to current datasets. 

Since ROC studies in medical imaging have serious health-care related consequences, no method should be used unless it has been thoroughly validated. Much work still remains to be done in proper simulator design, on which validation is dependent.

## Things for me to think about
### Expected values of mean squares 

Assuming no replications the expected mean squares are as follows, Table Table \@ref(tab:ExpValMs); understanding how this table is derived, would lead me outside my expertise and the scope of this book; suffice to say that these are *unconstrained* estimates (as summarized in the quotation above) which are different from the *constrained* estimates appearing in the original DBM publication [@RN204], Table 9.2; the differences between these two types of estimates is summarized in [@RN2079]. For reference, Table 9.3 is the table published in the most recent paper that I am aware of [@RN2508]. All three tables are different! **In this chapter I will stick to Table Table \@ref(tab:ExpValMs) for the subsequent development.** 

Source          df               E(MS)
-------         ----             ------------------------------------------------------------
T               (I-1)            $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$ + $J\sigma_{\tau C}^{2}$ + $JK\sigma_{\tau}^{2}$ 
R               (J-1)            $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ + $IK\sigma_{R}^{2}$ + $K\sigma_{\tau R}^{2}$
C               (K-1)            $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ + $IJ\sigma_{C}^{2}$ + $J\sigma_{\tau C}^{2}$
TR              (I-1)(J-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$  
TC              (I-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $J\sigma_{\tau C}^{2}$
RC              (J-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ 
TRC             (I-1)(J-1)(K-1)  $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ 
$\epsilon$      $N-1=0$          $\sigma_{\epsilon}^{2}$

Table: Table 9.1 Unconstrained expected values of mean-squares, as in [@RN2079]

Source          df               E(MS)
-------         ----             ------------------------------------------------------------
T               (I-1)            $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$ + $J\sigma_{\tau C}^{2}$ + $JK\sigma_{\tau}^{2}$ 
R               (J-1)            $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ + $IK\sigma_{R}^{2}$
C               (K-1)            $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ + $IJ\sigma_{C}^{2}$
TR              (I-1)(J-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$  
TC              (I-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $J\sigma_{\tau C}^{2}$
RC              (J-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $I\sigma_{RC}^{2}$ 
TRC             (I-1)(J-1)(K-1)  $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ 
$\epsilon$      0                $\sigma_{\epsilon}^{2}$

Table: Table 9.2 Constrained expected values of mean-squares, as in [@RN204] 

Source          df               E(MS)
-------         ----             ------------------------------------------------------------
T               (I-1)            $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$ + $J\sigma_{\tau C}^{2}$ + $JK\sigma_{\tau}^{2}$ 
R               (J-1)            $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $I\sigma_{RC}^{2}$ + $IK\sigma_{R}^{2}$ + $K\sigma_{\tau R}^{2}$
C               (K-1)            $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $I\sigma_{RC}^{2}$ + $IJ\sigma_{C}^{2}$ + $J\sigma_{\tau C}^{2}$
TR              (I-1)(J-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $K\sigma_{\tau R}^{2}$  
TC              (I-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $J\sigma_{\tau C}^{2}$
RC              (J-1)(K-1)       $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ + $I\sigma_{RC}^{2}$ 
TRC             (I-1)(J-1)(K-1)  $\sigma_{\epsilon}^{2}$ + $\sigma_{\tau RC}^{2}$ 
$\epsilon$      0                $\sigma_{\epsilon}^{2}$

Table: Table 9.3 As in Hillis "marginal-means ANOVA paper" [@RN2508]

## References {#dbm-analysis-significance-testing-references}

