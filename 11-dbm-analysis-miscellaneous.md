# DBM method special cases {#dbm-analysis-special-cases}
Special cases of DBM analysis are described here, namely fixed-reader random-case (FRRC), sub-special case of which is Single-reader multiple-treatment analysis, and random-reader fixed-case (RRFC). 


## TBA How much finished {#dbm-analysis-special-cases-fits-how-much-finished}
30%





## Fixed-reader random-case (FRRC) analysis {#FRRCDBMAnalysis}
The model is the same as in Eqn. \@ref(eq:DefDBMModel) except one sets $\sigma_{R}^{2}$ = $\sigma_{\tau R}^{2}$ = 0 in Table \@ref(tab:ExpValMs). The appropriate test statistic is: 

\begin{equation}
\frac{E\left ( MST \right )}{E\left ( MSTC \right )} = \frac{\sigma_{\epsilon}^{2}+\sigma_{\tau RC}^{2}+J\sigma_{\tau C}^{2}+JK\sigma_{\tau}^{2}}{\sigma_{\epsilon}^{2}+\sigma_{\tau RC}^{2}+J\sigma_{\tau C}^{2}}
\end{equation}

Under the null hypothesis $\sigma_{\tau}^{2} = 0$:

\begin{equation}
\frac{E\left ( MST \right )}{E\left ( MSTC \right )} = 1
\end{equation}

The F-statistic is (replacing *expected* with *observed* values):

\begin{equation}
F_{DBM|R}=\frac{MST}{MSTC}
(\#eq:FStatFRRC-DBM)
\end{equation}

The observed value $F_{DBM|R}$ (the Roe-Metz notation [@RN1124] is used which indicates that the factor appearing to the right of the vertical bar is regarded as fixed) is distributed as an F-statistic with $\text{ndf}$ = $I – 1$ and $ddf = (I-1)(K-1)$; the degrees of freedom follow from the rows labeled $T$ and $TC$ in TBA Table Table \@ref(tab:ExpValMs). Therefore, the distribution of the observed value is (no Satterthwaite approximation needed this time as both numerator and denominator are simple mean-squares):

\begin{equation}
F_{DBM|R} \sim F_{I-1,(I-1)(K-1)}
(\#eq:SamplingFStatFRRC)
\end{equation}

The null hypothesis is rejected if the observed value of the F- statistic exceeds the critical value:

\begin{equation}
F_{DBM|R} > F_{1-\alpha,I-1,(I-1)(K-1)}
(\#eq:NhRejectRuleFRRC)
\end{equation}

The p-value of the test is the probability that a random sample from the F-distribution TBA \@ref(eq:pseudoValPrime) Eqn. (9.39), exceeds the observed value:

\begin{equation}
p=\Pr\left ( F> F_{DBM|R} \mid F \sim F_{I-1,(I-1)(K-1)} \right )
(\#eq:pFRRC)
\end{equation}

The $(1-\alpha)$  confidence interval for the inter-treatment reader-averaged difference FOM is given by:

\begin{equation}
CI_{1-\alpha}=\left ( \theta_{i \bullet} - \theta_{i' \bullet} \right ) \pm t_{\alpha/2,(I-1)(K-1)}\sqrt{2\frac{MST}{JK}}
(\#eq:confIntervalFRRC)
\end{equation}

### Single-reader multiple-treatment analysis {#FRRCSingleReaderDBMAnalysis}
With a single reader interpreting cases in two or more treatments, the reader factor must necessarily be regarded as fixed. The preceding analysis is applicable. One simply puts $J = 1$ in the equations above. 

#### Example 5: Code illustrating p-values for FRRC analysis, Van Dyke data

```r
alpha <- 0.05
retMS <- UtilMeanSquares(dataset02)
I <- length(dataset02$ratings$NL[,1,1,1])
J <- length(dataset02$ratings$NL[1,,1,1])
K <- length(dataset02$ratings$NL[1,1,,1])
FDbmFR <- retMS$msT / retMS$msTC
ndf <- (I-1); ddf <- (I-1)*(K-1)
pValue <- 1 - pf(FDbmFR, ndf, ddf)

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

std_DIFF_FOM_FRRC <- sqrt(2*retMS$msTC/J/K)
nDiffs <- I*(I-1)/2
CI_DIFF_FOM_FRRC <- array(dim = c(nDiffs, 3))
for (i in 1 : nDiffs) {
  CI_DIFF_FOM_FRRC[i,1] <- qt(alpha/2,df = ddf)*std_DIFF_FOM_FRRC + trtDiff[i]
  CI_DIFF_FOM_FRRC[i,2] <- trtDiff[i]
  CI_DIFF_FOM_FRRC[i,3] <- qt(1-alpha/2,df = ddf)*std_DIFF_FOM_FRRC + trtDiff[i]
  print(data.frame("pValue" = pValue, 
                   "Lower" = CI_DIFF_FOM_FRRC[i,1], 
                   "Mid" = CI_DIFF_FOM_FRRC[i,2], 
                   "Upper" = CI_DIFF_FOM_FRRC[i,3]))
}
#>       pValue       Lower         Mid        Upper
#> 1 0.02103497 -0.08088303 -0.04380032 -0.006717613

retRJafroc <- StSignificanceTesting(dataset02, FOM = "Wilcoxon", method = "DBM")

data.frame("pValue" = retRJafroc$FRRC$FTests$p[1],
           "Lower" = retRJafroc$FRRC$ciDiffTrt[1,"CILower"], 
           "Mid" = retRJafroc$FRRC$ciDiffTrt[1,"Estimate"], 
           "Upper" = retRJafroc$FRRC$ciDiffTrt[1,"CIUpper"])
#>        pValue        Lower          Mid         Upper
#> 1 0.021034969 -0.080883031 -0.043800322 -0.0067176131
```

As one might expect, if one "freezes" reader variability, the FOM difference becomes significant, whether viewed from the point of view of the F-statistic exceeding the critical value, the observed p-value being smaller than alpha or the 95% CI for the difference FOM not including zero. 

## Random-reader fixed-case (RRFC) analysis {#dbm-analysis-special-cases-RRFCAnalysis}
The model is the same as in TBA \@ref(eq:pseudoValPrime) Eqn. (9.4) except one puts $\sigma_C^2 = \sigma_{\tau C}^2 =0$ in Table Table \@ref(tab:ExpValMs). It follows that: 

\begin{equation}
\frac{E(MST)}{E(MSTR)}=\frac{\sigma_\epsilon^2+\sigma_{\tau RC}^2+K\sigma_{\tau R}^2+JK\sigma_{\tau}^2}{\sigma_\epsilon^2+\sigma_{\tau RC}^2+K\sigma_{\tau R}^2}
\end{equation}

Under the null hypothesis $\sigma_\tau^2 = 0$:

\begin{equation}
\frac{E(MST)}{E(MSTR)}=1
\end{equation}

Therefore, one defines the F-statistic (replacing expected values with observed values) by:

\begin{equation}
F_{DBM|C} \sim \frac{MST}{MSTR}
(\#eq:FStatRRFC-Misc)
\end{equation}

The observed value $F_{DBM|C}$ is distributed as an F-statistic with $ndf = I – 1$ and $ddf = (I-1)(J-1)$, see rows labeled $T$ and $TR$ in Table Table \@ref(tab:ExpValMs).

\begin{equation}
F_{DBM|C} \sim F_{I-1,(I-1)(J-1))}
(\#eq:SamplingFStatRRFC)
\end{equation}

The null hypothesis is rejected if the observed value of the F statistic exceeds the critical value:

\begin{equation}
F_{DBM|C} > F_{1-\alpha,I-1,(I-1)(J-1))}
(\#eq:NhRejectRuleRRFC)
\end{equation}

The p-value of the test is the probability that a random sample from the distribution exceeds the observed value:

\begin{equation}
p=\Pr\left ( F>F_{DBM|C} \mid F \sim F_{I-1,(I-1)(J-1)} \right )
(\#eq:pRRFC)
\end{equation}

The confidence interval for inter-treatment differences is given by (TBA check this):

\begin{equation}
CI_{1-\alpha}=\left ( \theta_{i \bullet} - \theta_{i' \bullet} \right ) \pm t_{\alpha/2,(I-1)(J-1)}\sqrt{2\frac{MSTR}{JK}}
(\#eq:confIntervalRRFC)
\end{equation}

#### Example 6: Code illustrating analysis for RRFC analysis, Van Dyke data


```r
FDbmFC <- retMS$msT / retMS$msTR
ndf <- (I-1)
ddf <- (I-1)*(J-1)
pValue <- 1 - pf(FDbmFC, ndf, ddf)

nDiffs <- I*(I-1)/2
CI_DIFF_FOM_RRFC <- array(dim = c(nDiffs, 3))
for (i in 1 : nDiffs) {
  CI_DIFF_FOM_RRFC[i,1] <- qt(alpha/2,df = ddf)*sqrt(2*retMS$msTR/J/K) + trtDiff[i]
  CI_DIFF_FOM_RRFC[i,2] <- trtDiff[i]
  CI_DIFF_FOM_RRFC[i,3] <- qt(1-alpha/2,df = ddf)*sqrt(2*retMS$msTR/J/K) + trtDiff[i]
  print(data.frame("pValue" = pValue, 
                   "Lower" = CI_DIFF_FOM_RRFC[i,1], 
                   "Mid" = CI_DIFF_FOM_RRFC[i,2], 
                   "Upper" = CI_DIFF_FOM_RRFC[i,3]))
}
#>        pValue        Lower          Mid         Upper
#> 1 0.041958752 -0.085020224 -0.043800322 -0.0025804202
data.frame("pValue" = retRJafroc$RRFC$FTests$p[1],
           "Lower" = retRJafroc$RRFC$ciDiffTrt[1,"CILower"], 
           "Mid" = retRJafroc$RRFC$ciDiffTrt[1,"Estimate"], 
           "Upper" = retRJafroc$RRFC$ciDiffTrt[1,"CIUpper"])
#>        pValue        Lower          Mid         Upper
#> 1 0.041958752 -0.085020224 -0.043800322 -0.0025804202
```


## References {#dbm-analysis-special-cases-references}

