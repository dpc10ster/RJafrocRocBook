# Obuchowski Rockette Applications {#or-applications} 







## TBA How much finished {#or-applications-how-much-finished}
80%


## Introduction {#or-applications-introduction}  

This chapter illustrates Obuchowski-Rockette analysis with several examples. The first example is a full-blown "hand-calculation" for `dataset02`, showing explicit implementations of formulae presented in the previous chapter. The second example shows application of the `RJafroc` package function `StSignificanceTesting()` to the same dataset: this function encapsulates all formulae and accomplishes all analyses with one function call. The third example shows application of the `StSignificanceTesting()` function to an ROC dataset derived from the Federica Zanca dataset [@RN1882], which has five modalities and four readers. This illustrates multiple treatment pairings (in contrast, `dataset02` has only one treatment pairing). The fourth example shows application of `StSignificanceTesting()` to `dataset04`, which is an **FROC** dataset (in contrast to the previous examples, which employed **ROC** datasets). It illustrates the key difference involved in FROC analysis, namely the choice of figure of merit. The final example again uses `dataset04`, i.e., FROC data, *but this time we use DBM analysis*. Since DBM analysis is pseudovalue based, and the figure of merit is not the empirical AUC under the ROC, one may expect to see differences from the previously presented OR analysis on the same dataset.

Each analysis involves the following steps: 

* Calculate the figure of merit; 
* Calculate the variance-covariance matrix and mean-squares;
* Calculate the NH statistic, p-value and confidence interval(s).
* For each analysis, three sub-analyses are shown: 
    + random-reader random-case (RRRC),
    + fixed-reader random-case (FRRC), and
    + random-reader fixed-case (RRFC).

## Hand calculation {#or-applications-dataset02-hand}

Dataset `dataset02` is well-know in the literature [@RN1993] as it has been widely used to illustrate advances in ROC methodology. The following code extract the numbers of modalities, readers and cases for `dataset02` and defines strings `modalityID`, `readerID` and `diffTRName` that are needed for the hand-calculations.


```r
I <- length(dataset02$ratings$NL[,1,1,1])
J <- length(dataset02$ratings$NL[1,,1,1])
K <- length(dataset02$ratings$NL[1,1,,1])
modalityID <- dataset02$descriptions$modalityID
readerID <- dataset02$descriptions$readerID
diffTRName <- array(dim = choose(I, 2))
ii <- 1
for (i in 1:I) {
  if (i == I) 
    break
  for (ip in (i + 1):I) {
    diffTRName[ii] <- 
      paste0("trt", modalityID[i], 
             sep = "-", "trt", modalityID[ip])
    ii <- ii + 1
  }
}
```

The dataset consists of I = 2 treatments,  J = 5 readers and  K = 114 cases.

### Random-Reader Random-Case (RRRC) analysis {#or-applications-RRRC-dataset02-hand}
* The first step is to calculate the figures of merit using `UtilFigureOfMerit()`. 
* Note that the `FOM` argument has to be explicitly specified as there is no default.


```r
foms <- UtilFigureOfMerit(dataset02, FOM = "Wilcoxon")
print(foms, digits = 4)
#>        rdr0   rdr1   rdr2   rdr3   rdr4
#> trt0 0.9196 0.8588 0.9039 0.9731 0.8298
#> trt1 0.9478 0.9053 0.9217 0.9994 0.9300
```

* For example, for the first treatment, `"trt0"`, the second reader `"rdr1"` figure of merit is 0.8587762.
* The next step is to calculate the variance-covariance matrix and the mean-squares.
* The function `UtilORVarComponentsFactorial()` returns these quantities, which are saved to `vc`. 
* The `Factorial` in the function name is because this code applies to the factorial design. A different function is used for a split-plot design.


```r
vc <- UtilORVarComponentsFactorial(
  dataset02, FOM = "Wilcoxon", covEstMethod = "jackknife")
print(vc, digits = 4)
#> $TRanova
#>          SS DF       MS
#> T  0.004796  1 0.004796
#> R  0.015345  4 0.003836
#> TR 0.002204  4 0.000551
#> 
#> $VarCom
#>       Estimates   Rhos
#> VarR  0.0015350     NA
#> VarTR 0.0002004     NA
#> Cov1  0.0003466 0.4320
#> Cov2  0.0003441 0.4289
#> Cov3  0.0002390 0.2979
#> Var   0.0008023     NA
#> 
#> $IndividualTrt
#>      DF msREachTrt varEachTrt cov2EachTrt
#> trt0  4   0.003083  0.0010141   0.0004840
#> trt1  4   0.001305  0.0005905   0.0002042
#> 
#> $IndividualRdr
#>      DF msTEachRdr varEachRdr cov1EachRdr
#> rdr0  1  0.0003971  0.0006989   3.735e-04
#> rdr1  1  0.0010829  0.0011061   7.602e-04
#> rdr2  1  0.0001597  0.0008423   3.553e-04
#> rdr3  1  0.0003445  0.0001506   1.083e-06
#> rdr4  1  0.0050161  0.0012136   2.430e-04
```

* The next step is the calculate the NH testing statistic. 
* The relevant equation is Eqn. \@ref(eq:F-ORH-RRRC). 
* `vc` contains the values needed in this equation, as follows:
    + MS(T) is in `vc$TRanova["T", "MS"]`, whose value is 0.0047962. 
    + MS(TR) is in `vc$TRanova["TR", "MS"]`, whose value is 5.5103062\times 10^{-4}. 
    + `Cov2` is in `vc$VarCom["Cov2", "Estimates"]`, whose value is 3.4407483\times 10^{-4}. 
    + `Cov3` is in `vc$VarCom["Cov3", "Estimates"]`, whose value is 2.3902837\times 10^{-4}. 

Applying Eqn. \@ref(eq:F-ORH-RRRC) one gets (`den` is the denominator on the right hand side of the referenced equation) and F_ORH_RRRC is the value of the F-statistic:


```r
den <- vc$TRanova["TR", "MS"] + 
  J* max(vc$VarCom["Cov2", "Estimates"] - 
           vc$VarCom["Cov3", "Estimates"],0)
F_ORH_RRRC <- vc$TRanova["T", "MS"]/den
print(F_ORH_RRRC, digits = 4)
#> [1] 4.456
```

* The F-statistic has numerator degrees of freedom $\text{ndf} = I - 1$ and denominator degrees of freedom, `ddf`, to be calculated next.
* From the previous chapter, `ddf` is calculated using Eqn. \@ref(eq:ddfH-RRRC)). The numerator of `ddf` is identical to `den^2`, where `den` was calculated in the preceding code block. The implementation follows:


```r
ddf <- den^2*(I-1)*(J-1)/(vc$TRanova["TR", "MS"])^2
print(ddf, digits = 4)
#> [1] 15.26
```

* The next step is calculation of the p-value for rejecting the NH
* The relevant equation is Eqn. \@ref(eq:pValueOR-RRRC) whose implementation follows: 


```r
p <- 1 - pf(F_ORH_RRRC, I - 1, ddf)
print(p, digits = 4)
#> [1] 0.05167
```

* The difference is not significant at $\alpha$ = 0.05. 
* The next step is to calculate confidence intervals.
* Since `I` = 2, their is only one paired difference in reader-averaged FOMs, namely, the first treatment minus the second.


```r
trtMeans <- rowMeans(foms)
trtMeanDiffs <- trtMeans[1] - trtMeans[2]
names(trtMeanDiffs) <- "trt0-trt1"
print(trtMeans, digits = 4)
#>   trt0   trt1 
#> 0.8970 0.9408
print(trtMeanDiffs, digits = 4)
#> trt0-trt1 
#>   -0.0438
```

* `trtMeans`contains the reader-averaged figures of merit for each treatment.
* `trtMeanDiffs`contains the reader-averaged difference figure of merit.
* From the previous chapter, the $(1-\alpha)$ confidence interval for $\theta_{1 \bullet} - \theta_{2 \bullet}$ is given by Eqn. \@ref(eq:CI-DiffFomRRRC), in which equation the expression inside the square-root symbol is `2/J*den`. 
* $\alpha$, the significance level of the test, is set to 0.05. 
* The implementation follows:


```r
alpha <- 0.05
stdErr <- sqrt(2/J*den)
t_crit <- abs(qt(alpha/2, ddf))
CI_RRRC <- c(trtMeanDiffs - t_crit*stdErr, 
             trtMeanDiffs + t_crit*stdErr)
names(CI_RRRC) <- c("Lower", "Upper")
print(CI_RRRC, digits = 4)
#>      Lower      Upper 
#> -0.0879595  0.0003589
```

The confidence interval includes zero, which confirms the F-statistic finding that the reader-averaged FOM difference between treatments is not significant. 

Calculated next is the confidence interval for the reader-averaged FOM for each treatment, i.e. $CI_{1-\alpha,RRRC,\theta_{i \bullet}}$. The relevant equations are Eqn. \@ref(eq:CI-RRRC-df-IndvlTrt) and Eqn. \@ref(eq:CI-RRRC-IndvlTrt). The implementation follows:


```r
df_i <- array(dim = I)
den_i <- array(dim = I)
stdErr_i <- array(dim = I)
ci <- array(dim = c(I, 2))
CI_RRRC_IndvlTrt <- data.frame()
for (i in 1:I) {
  den_i[i] <- vc$IndividualTrt[i, "msREachTrt"] + 
    J * max(vc$IndividualTrt[i, "cov2EachTrt"], 0)
  df_i[i] <- 
    (den_i[i])^2/(vc$IndividualTrt[i, "msREachTrt"])^2 * (J - 1)
  stdErr_i[i] <- sqrt(den_i[i]/J)
  ci[i,] <- 
    c(trtMeans[i] + qt(alpha/2, df_i[i]) * stdErr_i[i], 
      trtMeans[i] + qt(1-alpha/2, df_i[i]) * stdErr_i[i])
  rowName <- paste0("trt", modalityID[i])
  CI_RRRC_IndvlTrt <- rbind(
    CI_RRRC_IndvlTrt, 
    data.frame(Estimate = trtMeans[i], 
               StdErr = stdErr_i[i],
               DFi = df_i[i],
               CILower = ci[i,1],
               CIUpper = ci[i,2],
               Cov2i = vc$IndividualTrt[i,"cov2EachTrt"],
               row.names = rowName,
               stringsAsFactors = FALSE))
}
print(CI_RRRC_IndvlTrt, digits = 4)
#>      Estimate  StdErr   DFi CILower CIUpper     Cov2i
#> trt0   0.8970 0.03317 12.74  0.8252  0.9689 0.0004840
#> trt1   0.9408 0.02157 12.71  0.8941  0.9875 0.0002042
```


### Fixed-Reader Random-Case (FRRC) analysis {#or-applications-FRRC-dataset02-hand}
* The chi-square statistic is calculated using Eqn. \@ref(eq:DefFStatFRRC-OR) and Eqn. \@ref(eq:ChisqStatFRRC-OR). 
* The needed quantities are in `vc`. 
* For example, MS(T) is in vc$TRanova["T", "MS"], see above. Likewise for `Cov2` and `Cov3`.
* The remaining needed quantities are:
+ `Var` is in `vc$VarCom["Var", "Estimates"]`, whose value is 8.0228827\times 10^{-4}. 
+ `Cov1` is in `vc$VarCom["Cov1", "Estimates"]`, whose value is 3.4661371\times 10^{-4}. 
* The degree of freedom is $I-1$.
* The implementation follows:


```r
den_FRRC <- vc$VarCom["Var","Estimates"] - 
  vc$VarCom["Cov1","Estimates"] + 
  (J - 1) * max(vc$VarCom["Cov2","Estimates"] - 
                  vc$VarCom["Cov3","Estimates"] ,0)
chisqVal <- (I-1)*vc$TRanova["T","MS"]/den_FRRC
p <- 1 - pchisq(chisqVal, I - 1)
FTests <- data.frame(MS = c(vc$TRanova["T", "MS"], den_FRRC),
                     Chisq = c(chisqVal,NA),
                     DF = c(I - 1, NA),
                     p = c(p,NA),
                     row.names = c("Treatment", "Error"),
                     stringsAsFactors = FALSE)
print(FTests, digits = 4)
#>                  MS Chisq DF       p
#> Treatment 0.0047962 5.476  1 0.01928
#> Error     0.0008759    NA NA      NA
```

* Since p < 0.05, one has a significant finding. 
* Freezing reader variability shows a significant difference between the treatments. 
* The downside is that the conclusion applies only to the readers used in the study.
* The next step is to calculate the confidence interval for the reader-averaged FOM difference, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* The relevant equation is Eqn. \@ref(eq:CIDiffFomFRRC-OR), whose implementation follows.


```r
stdErr <- sqrt(2 * den_FRRC/J)
zStat <- vector()
PrGTz <- vector()
CI <- array(dim = c(choose(I,2),2))
for (i in 1:choose(I,2)) {
  zStat[i] <- trtMeanDiffs[i]/stdErr
  PrGTz[i] <- 2 * pnorm(abs(zStat[i]), lower.tail = FALSE)
  CI[i, ] <- c(trtMeanDiffs[i] + qnorm(alpha/2) * stdErr, 
               trtMeanDiffs[i] + qnorm(1-alpha/2) * stdErr)
}
ciDiffTrtFRRC <- data.frame(Estimate = trtMeanDiffs, 
                            StdErr = rep(stdErr, choose(I, 2)),
                            z = zStat, 
                            PrGTz = PrGTz, 
                            CILower = CI[,1],
                            CIUpper = CI[,2], 
                            row.names = diffTRName,
                            stringsAsFactors = FALSE)
print(ciDiffTrtFRRC, digits = 4)
#>           Estimate  StdErr     z   PrGTz  CILower   CIUpper
#> trt0-trt1  -0.0438 0.01872 -2.34 0.01928 -0.08049 -0.007115
```

* Consistent with the chi-square statistic significant finding, one finds that the treatment difference confidence interval does not include zero.
* The next step is to calculate the confidence interval for the reader-averaged figures of merit for each treatment, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet}}$.
* The relevant formula is in Eqn. \@ref(eq:CIIndTrtFomFRRC-OR), whose implementation follows:


```r
stdErr <- vector()
df <- vector()
CI <- array(dim = c(I,2))
ciAvgRdrEachTrt <- data.frame()
for (i in 1:I) {
  df[i] <- K - 1
  stdErr[i] <- 
    sqrt((vc$IndividualTrt[i,"varEachTrt"] + 
            (J-1)*max(vc$IndividualTrt[i,"cov2EachTrt"],0))/J)
  CI[i, ] <- c(trtMeans[i] + qnorm(alpha/2) * stdErr[i],
               trtMeans[i] + qnorm(1-alpha/2) * stdErr[i])
  rowName <- paste0("trt", modalityID[i])
  ciAvgRdrEachTrt <- 
    rbind(ciAvgRdrEachTrt, 
          data.frame(Estimate = trtMeans[i], 
                     StdErr = stdErr[i],
                     DF = df[i],
                     CILower = CI[i,1],
                     CIUpper = CI[i,2],
                     row.names = rowName,
                     stringsAsFactors = FALSE))
}
print(ciAvgRdrEachTrt, digits = 4)
#>      Estimate  StdErr  DF CILower CIUpper
#> trt0   0.8970 0.02429 113  0.8494  0.9446
#> trt1   0.9408 0.01678 113  0.9080  0.9737
```
* Finally, one calculates confidence intervals for the FOM differences for individual readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i j} - \theta_{i' j}}$. 
* The relevant formula is in Eqn. \@ref(eq:CIIndRdrDiffFomFRRC-OR), whose implementation follows:


```r
trtMeanDiffs1 <- array(dim = c(J, choose(I, 2)))
Reader <- array(dim = c(J, choose(I, 2)))
stdErr <- array(dim = c(J, choose(I, 2)))
zStat <- array(dim = c(J, choose(I, 2)))
trDiffNames <- array(dim = c(J, choose(I, 2)))
PrGTz <- array(dim = c(J, choose(I, 2)))
CIReader <- array(dim = c(J, choose(I, 2),2))
ciDiffTrtEachRdr <- data.frame()
for (j in 1:J) {
  Reader[j,] <- rep(readerID[j], choose(I, 2))
  stdErr[j,] <- 
    sqrt(
      2 * 
        (vc$IndividualRdr[j,"varEachRdr"] - 
           vc$IndividualRdr[j,"cov1EachRdr"]))
  pair <- 1
  for (i in 1:I) {
    if (i == I) break
    for (ip in (i + 1):I) {
      trtMeanDiffs1[j, pair] <- foms[i, j] - foms[ip, j]
      trDiffNames[j,pair] <- diffTRName[pair]
      zStat[j,pair] <- trtMeanDiffs1[j,pair]/stdErr[j,pair]
      PrGTz[j,pair] <- 
        2 * pnorm(abs(zStat[j,pair]), lower.tail = FALSE)
      CIReader[j, pair,] <- 
        c(trtMeanDiffs1[j,pair] + 
            qnorm(alpha/2) * stdErr[j,pair], 
          trtMeanDiffs1[j,pair] + 
            qnorm(1-alpha/2) * stdErr[j,pair])
      rowName <- 
        paste0("rdr", Reader[j,pair], "::", trDiffNames[j, pair])
      ciDiffTrtEachRdr <- rbind(
        ciDiffTrtEachRdr, 
        data.frame(Estimate = trtMeanDiffs1[j, pair], 
                   StdErr = stdErr[j,pair], 
                   z = zStat[j, pair], 
                   PrGTz = PrGTz[j, pair], 
                   CILower = CIReader[j, pair,1],
                   CIUpper = CIReader[j, pair,2],
                   row.names = rowName,
                   stringsAsFactors = FALSE))
      pair <- pair + 1
    }
  }
}
print(ciDiffTrtEachRdr, digits = 3)
#>                 Estimate StdErr      z  PrGTz CILower  CIUpper
#> rdr0::trt0-trt1  -0.0282 0.0255 -1.105 0.2693 -0.0782  0.02182
#> rdr1::trt0-trt1  -0.0465 0.0263 -1.769 0.0768 -0.0981  0.00501
#> rdr2::trt0-trt1  -0.0179 0.0312 -0.573 0.5668 -0.0790  0.04330
#> rdr3::trt0-trt1  -0.0262 0.0173 -1.518 0.1290 -0.0601  0.00764
#> rdr4::trt0-trt1  -0.1002 0.0441 -2.273 0.0230 -0.1865 -0.01381
```

The notation in the first column shows the reader and the treatment pairing. For example, `rdr1::trt0-trt1` means the FOM difference for reader `rdr1`. Only the fifth reader, i.e., `rdr4`, shows a significant difference between the treatments: the p-value is 0.023001 and the confidence interval also does not include zero. The large FOM difference for this reader, -0.100161, was enough to result in a significant finding for FRRC analysis. The FOM differences for the other readers are about a factor of 2.1522491 or more smaller than that for this reader.

### Random-Reader Fixed-Case (RRFC) analysis {#or-applications-RRFC-dataset02-hand}
The F-statistic is shown in Eqn. \@ref(eq:DefFStatRRFC). This time `ndf` = $I-1$ and `ddf` = $(I-1) \times (J-1)$, the values proposed in the Obuchowski-Rockette paper. The implementation follows:


```r
den <- vc$TRanova["TR","MS"]
f <- vc$TRanova["T","MS"]/den
ddf <- ((I - 1) * (J - 1))
p <- 1 - pf(f, I - 1, ddf)
FTests_RRFC <- 
  data.frame(DF = c(I-1,(I-1)*(J-1)), 
             MS = c(vc$TRanova["T","MS"],vc$TRanova["TR","MS"]), 
             F = c(f,NA),  p = c(p,NA), 
             row.names = c("T","TR"), 
             stringsAsFactors = FALSE)
print(FTests_RRFC, digits = 4)
#>    DF       MS     F       p
#> T   1 0.004796 8.704 0.04196
#> TR  4 0.000551    NA      NA
```

Freezing case variability also results in a significant finding, but the conclusion is only applicable to the specific case set used in the study. Next one calculates confidence intervals for the reader-averaged FOM differences, the relevant formula is in Eqn. \@ref(eq:CIDiffFomRRFC), whose implementation follows.


```r
stdErr <- sqrt(2 * den/J)
tStat <- vector()
PrGTt <- vector()
CI <- array(dim = c(choose(I,2), 2))
for (i in 1:choose(I,2)) {
  tStat[i] <- trtMeanDiffs[i]/stdErr
  PrGTt[i] <- 2 * 
    pt(abs(tStat[i]), ddf, lower.tail = FALSE)
  CI[i, ] <- c(trtMeanDiffs[i] + qt(alpha/2, ddf) * stdErr, 
               trtMeanDiffs[i] + qt(1-alpha/2, ddf) * stdErr)
}
ciDiffTrt_RRFC <- 
  data.frame(Estimate = trtMeanDiffs, 
             StdErr = rep(stdErr, choose(I, 2)), 
             DF = rep(ddf, choose(I, 2)), 
             t = tStat, 
             PrGTt = PrGTt, 
             CILower = CI[,1],
             CIUpper = CI[,2],
             row.names = diffTRName, 
             stringsAsFactors = FALSE)

print(ciDiffTrt_RRFC, digits = 4)
#>           Estimate  StdErr DF     t   PrGTt  CILower  CIUpper
#> trt0-trt1  -0.0438 0.01485  4 -2.95 0.04196 -0.08502 -0.00258
```
* As expected because the overall F-test showed significance, the confidence interval does not include zero (the p-value is identical to that found by the F-test). 
* This completes the hand calculations.

## RJafroc: dataset02 {#or-applications-dataset02-RJafroc}

The second example shows application of the `RJafroc` package function `StSignificanceTesting()` to `dataset02`. This function encapsulates all formulae discussed previously and accomplishes the analyses with a single function call. It returns an object, denoted `st1` below, that contains all results of the analysis. It is a `list` with the following components:

* `FOMs`, this in turn is a `list` containing the following data frames: 
    + `foms`, the individual treatment-reader figures of merit, i.e., $\theta_{i j}$, 
    + `trtMeans`, the treatment figures of merit averaged over readers, i.e., $\theta_{i \bullet}$,
    + `trtMeanDiffs`, the inter-treatment figures of merit differences averaged over readers, i.e., $\theta_{i \bullet}-\theta_{i' \bullet}$.

* `ANOVA`, a `list` containing the following data frames: 
    + `TRanova`, the treatment-reader ANOVA table,
    + `VarCom`, Obuchowski-Rockette variance-covariances and correlations,
    + `IndividualTrt`, the mean-squares, `Var` and `Cov2` calculated over individual treatments,
    + `IndividualRdr`, the mean-squares, `Var` and `Cov1` calculated over individual readers.

* `RRRC`, a `list` containing the following data frames: 
    + `FTests`, the results of the F-test,
    + `ciDiffTrt`, the confidence intervals for inter-treatment FOM differences, averaged over readers, denoted $CI_{1-\alpha,RRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$ in the previous chapter,
    + `ciAvgRdrEachTrt`, the confidence intervals for individual treatment FOMs, averaged over readers, denoted $CI_{1-\alpha,RRRC,\theta_{i \bullet}}$ in the previous chapter.

* `FRRC`, a `list` containing the following data frames: 
    + `FTests`, the results of the F-tests, which in this case specializes to chi-square tests,
    + `ciDiffTrt`, the confidence intervals for inter-treatment FOM differences, averaged over readers, denoted $CI_{1-\alpha,FRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$ in the previous chapter,
    + `ciAvgRdrEachTrt`, the confidence intervals for individual treatment FOMs, averaged over readers, denoted $CI_{1-\alpha,FRRC,\theta_{i \bullet}}$ in the previous chapter,
    + `ciDiffTrtEachRdr`, the confidence intervals for inter-treatment FOM differences for individual readers, denoted $CI_{1-\alpha,FRRC,\theta_{ij} - \theta_{i'j}}$ in the previous chapter,
    + `IndividualRdrVarCov1`, the individual reader variance-covariances and means squares.

* `RRFC`, a `list` containing the following data frames: 
    + `FTests`, the results of the F-tests, which in this case specializes to chi-square tests,
    + `ciDiffTrt`, the confidence intervals for inter-treatment FOM differences, averaged over readers, denoted $CI_{1-\alpha,RRFC,\theta_{i \bullet} - \theta_{i' \bullet}}$ in the previous chapter,
    + `ciAvgRdrEachTrt`, the confidence intervals for indvidual treatment FOMs, averaged over readers, denoted $CI_{1-\alpha,RRFC,\theta_{i \bullet}}$ in the previous chapter.

In the interest of clarity, in the first example using the `RJafroc` package the components of the returned object `st1` are listed separately and described explicitly. In the interest of brevity, in subsequent examples the object is listed in its entirety.

Online help on the `StSignificanceTesting()` function is available:


```r
?`StSignificanceTesting`
```

The lower right `RStudio` panel contains the online description. Click on the small up-and-right pointing arrow icon to expand this to a new window. 

### Random-Reader Random-Case (RRRC) analysis {#or-applications-RRRC-dataset02-RJafroc}
* Since `analysisOption` is not explicitly specified in the following code, the function `StSignificanceTesting` performs all three analyses: `RRRC`, `FRRC` and `RRFC`.
* Likewise, the significance level of the test, also an argument, `alpha`, defaults to 0.05. 
* The code below applies `StSignificanceTesting()` and saves the returned object to `st1`. 
* The first member of this object, a  `list` named `FOMs`, is then displayed. 
* `FOMs` contains three data frames: 
    + `FOMS$foms`, the figures of merit for each treatment and reader, 
    + `FOMS$trtMeans`, the figures of merit for each treatment averaged over readers, and 
    + `FOMS$trtMeanDiffs`, the inter-treatment difference figures of merit averaged over readers. The difference is always the first treatment minus the second, etc., in this example, `trt0` minus `trt1`.


```r
st1 <- StSignificanceTesting(dataset02, FOM = "Wilcoxon", method = "OR")
print(st1$FOMs, digits = 4)
#> $foms
#>        rdr0   rdr1   rdr2   rdr3   rdr4
#> trt0 0.9196 0.8588 0.9039 0.9731 0.8298
#> trt1 0.9478 0.9053 0.9217 0.9994 0.9300
#> 
#> $trtMeans
#>      Estimate
#> trt0   0.8970
#> trt1   0.9408
#> 
#> $trtMeanDiffs
#>           Estimate
#> trt0-trt1  -0.0438
```

* Displayed next are the variance components and mean-squares contained in the `ANOVA` `list`. 
    * `ANOVA$TRanova` contains the treatment-reader ANOVA table, i.e. the sum of squares, the degrees of freedom and the mean-squares, for treatment, reader and treatment-reader factors, i.e., `T`, `R` and `TR`.
    * `ANOVA$VarCom` contains the OR variance components and the correlations.
    * `ANOVA$IndividualTrt` contains the quantities necessary for individual treatment analyses.
    * `ANOVA$IndividualRdr` contains the quantities necessary for individual reader analyses.


```r
print(st1$ANOVA, digits = 4)
#> $TRanova
#>          SS DF       MS
#> T  0.004796  1 0.004796
#> R  0.015345  4 0.003836
#> TR 0.002204  4 0.000551
#> 
#> $VarCom
#>       Estimates   Rhos
#> VarR  0.0015350     NA
#> VarTR 0.0002004     NA
#> Cov1  0.0003466 0.4320
#> Cov2  0.0003441 0.4289
#> Cov3  0.0002390 0.2979
#> Var   0.0008023     NA
#> 
#> $IndividualTrt
#>      DF msREachTrt varEachTrt cov2EachTrt
#> trt0  4   0.003083  0.0010141   0.0004840
#> trt1  4   0.001305  0.0005905   0.0002042
#> 
#> $IndividualRdr
#>      DF msTEachRdr varEachRdr cov1EachRdr
#> rdr0  1  0.0003971  0.0006989   3.735e-04
#> rdr1  1  0.0010829  0.0011061   7.602e-04
#> rdr2  1  0.0001597  0.0008423   3.553e-04
#> rdr3  1  0.0003445  0.0001506   1.083e-06
#> rdr4  1  0.0050161  0.0012136   2.430e-04
```

* Displayed next are the results of the RRRC significance test, contained in `st1$RRRC`.


```r
print(st1$RRRC$FTests, digits = 4)
#>              DF       MS FStat       p
#> Treatment  1.00 0.004796 4.456 0.05167
#> Error     15.26 0.001076    NA      NA
```

* `st1$RRRC$FTests` contains the results of the F-tests: the degrees of freedom, the mean-squares, the observed value of the F-statistic and the p-value for rejecting the NH, listed separately, where applicable, for the treatment and error terms. 
* For example, the treatment mean squares is `st1$RRRC$FTests["Treatment", "MS"]` whose value is 0.00479617.


```r
print(st1$RRRC$ciDiffTrt, digits = 3)
#>           Estimate StdErr   DF     t  PrGTt CILower  CIUpper
#> trt0-trt1  -0.0438 0.0207 15.3 -2.11 0.0517  -0.088 0.000359
```

* `st1$RRRC$ciDiffTrt` contains the results of the confidence intervals for the inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.


```r
print(st1$RRRC$ciAvgRdrEachTrt, digits = 4)
#>      Estimate  StdErr    DF CILower CIUpper      Cov2
#> trt0   0.8970 0.03317 12.74  0.8252  0.9689 0.0004840
#> trt1   0.9408 0.02157 12.71  0.8941  0.9875 0.0002042
```

* `st1$RRRC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet}}$.

### Fixed-Reader Random-Case (FRRC) analysis {#or-applications-FRRC-dataset02-RJafroc}

* Displayed next are the results of FRRC analysis, contained in `st1$FRRC`.
* `st1$FRRC$FTests` contains the results of the F-tests: the degrees of freedom, the mean-squares, the observed value of the F-statistic and the p-value for rejecting the NH, listed separately, where applicable, for the treatment and error terms. 
* For example, the treatment mean squares is `st1$FRRC$FTests["Treatment", "MS"]` whose value is 0.00479617.


```r
print(st1$FRRC$FTests, digits = 4)
#>                  MS Chisq DF       p
#> Treatment 0.0047962 5.476  1 0.01928
#> Error     0.0008759    NA NA      NA
```

* Note that this time the output lists a chi-square distribution observed value, 5.47595324, with degree of freedom $df = I -1 = 1$.
* The listed mean-squares and the p-value agree with the previously performed hand calculations.
* For FRRC analysis the value of the chi-square statistic is significant and the p-value is smaller than $\alpha$.


```r
print(st1$FRRC$ciDiffTrt, digits = 4)
#>           Estimate  StdErr     z   PrGTz  CILower   CIUpper
#> trt0-trt1  -0.0438 0.01872 -2.34 0.01928 -0.08049 -0.007115
```

* `st1$FRRC$ciDiffTrt` contains confidence intervals for inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* The confidence interval excludes zero, and the p-value, listed under `PrGTz` (for probability greater than `z`) is smaller than 0.05.
* One could be using the t-distribution with infinite degrees of freedom, but this is identical to the normal distribution. Hence the listed value is a `z` statistic, i.e., `z = -0.043800322/0.018717483` = -2.34007543.


```r
print(st1$FRRC$ciAvgRdrEachTrt, digits = 4)
#>      Estimate  StdErr  DF CILower CIUpper
#> trt0   0.8970 0.02429 113  0.8494  0.9446
#> trt1   0.9408 0.01678 113  0.9080  0.9737
```

* `st1$FRRC$st1$FRRC$ciAvgRdrEachTrt` contains confidence intervals for individual treatment FOMs, averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet}}$.



```r
print(st1$FRRC$ciDiffTrtEachRdr, digits = 3)
#>                 Estimate StdErr      z  PrGTz CILower  CIUpper
#> rdr0::trt0-trt1  -0.0282 0.0255 -1.105 0.2693 -0.0782  0.02182
#> rdr1::trt0-trt1  -0.0465 0.0263 -1.769 0.0768 -0.0981  0.00501
#> rdr2::trt0-trt1  -0.0179 0.0312 -0.573 0.5668 -0.0790  0.04330
#> rdr3::trt0-trt1  -0.0262 0.0173 -1.518 0.1290 -0.0601  0.00764
#> rdr4::trt0-trt1  -0.1002 0.0441 -2.273 0.0230 -0.1865 -0.01381
```

* `st1$FRRC$st1$FRRC$ciDiffTrtEachRdr` contains confidence intervals for inter-treatment difference FOMs, for each reader, i.e., $CI_{1-\alpha,FRRC,\theta_{i j} - \theta_{i' j}}$.

### Random-Reader Fixed-Case (RRFC) analysis {#or-applications-RRFC-dataset02-RJafroc}


```r
print(st1$RRFC$FTests, digits = 4)
#>    DF       MS     F       p
#> T   1 0.004796 8.704 0.04196
#> TR  4 0.000551    NA      NA
```

* `st1$RRFC$FTests` contains results of the F-test: the degrees of freedom, the mean-squares, the observed value of the F-statistic and the p-value for rejecting the NH, listed separately, where applicable, for the treatment and treatment-reader terms. The latter is also termed the "error term". 
* For example, the treatment-reader mean squares is `st1$RRFC$FTests["TR", "MS"]` whose value is 5.51030622\times 10^{-4}.


```r
print(st1$RRFC$ciDiffTrt, digits = 4)
#>           Estimate  StdErr DF     t   PrGTt  CILower  CIUpper
#> trt0-trt1  -0.0438 0.01485  4 -2.95 0.04196 -0.08502 -0.00258
```

* `st1$RRFC$ciDiffTrt` contains confidence intervals for the inter-treatment paired difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet} - \theta_{i' \bullet}}$.



```r
print(st1$RRFC$ciAvgRdrEachTrt, digits = 4)
#>      Estimate  StdErr DF CILower CIUpper
#> Trt0   0.8970 0.02483  4  0.8281  0.9660
#> Trt1   0.9408 0.01615  4  0.8960  0.9857
```

* `st1$RRFC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet}}$.

## RJafroc: dataset04 {#or-applications-dataset04-RJafroc}
* The third example uses the Federica Zanca dataset [@RN1882], i.e., `dataset04`, which has five modalities and four readers. 
* It illustrates the situation when multiple treatment pairings are involved. In contrast, the previous example had only one treatment pairing.
* Since this is an FROC dataset, in order to keep it comparable with the previous example, one converts it to an inferred-ROC dataset.
* The function `DfFroc2Roc(dataset04)` converts, using the highest-rating, the FROC dataset to an inferred-ROC dataset.
* The results are contained in `st2`. 
* As noted earlier, this time the object is listed in its entirety.


```r
ds <- DfFroc2Roc(dataset04) # convert to ROC
I <- length(ds$ratings$NL[,1,1,1])
J <- length(ds$ratings$NL[1,,1,1])
cat("I = ", I, ", J = ", J, "\n")
#> I =  5 , J =  4
st2 <- StSignificanceTesting(ds, FOM = "Wilcoxon", method = "OR")
print(st2, digits = 3)
#> $FOMs
#> $FOMs$foms
#>       rdr1  rdr2  rdr3  rdr4
#> trt1 0.904 0.798 0.812 0.866
#> trt2 0.864 0.845 0.821 0.872
#> trt3 0.813 0.816 0.753 0.857
#> trt4 0.902 0.832 0.789 0.880
#> trt5 0.841 0.773 0.771 0.848
#> 
#> $FOMs$trtMeans
#>      Estimate
#> trt1    0.845
#> trt2    0.850
#> trt3    0.810
#> trt4    0.851
#> trt5    0.808
#> 
#> $FOMs$trtMeanDiffs
#>            Estimate
#> trt1-trt2 -0.005100
#> trt1-trt3  0.035325
#> trt1-trt4 -0.005412
#> trt1-trt5  0.036775
#> trt2-trt3  0.040425
#> trt2-trt4 -0.000312
#> trt2-trt5  0.041875
#> trt3-trt4 -0.040737
#> trt3-trt5  0.001450
#> trt4-trt5  0.042187
#> 
#> 
#> $ANOVA
#> $ANOVA$TRanova
#>         SS DF       MS
#> T  0.00759  4 0.001897
#> R  0.02188  3 0.007294
#> TR 0.00555 12 0.000462
#> 
#> $ANOVA$VarCom
#>       Estimates  Rhos
#> VarR   1.28e-03    NA
#> VarTR -1.09e-05    NA
#> Cov1   2.95e-04 0.374
#> Cov2   2.33e-04 0.296
#> Cov3   2.12e-04 0.269
#> Var    7.89e-04    NA
#> 
#> $ANOVA$IndividualTrt
#>      DF msREachTrt varEachTrt cov2EachTrt
#> trt1  3   0.002422   0.000711    0.000211
#> trt2  3   0.000523   0.000751    0.000266
#> trt3  3   0.001855   0.000876    0.000246
#> trt4  3   0.002578   0.000727    0.000220
#> trt5  3   0.001766   0.000882    0.000222
#> 
#> $ANOVA$IndividualRdr
#>      DF msTEachRdr varEachRdr cov1EachRdr
#> rdr1  4   0.001551   0.000689    0.000215
#> rdr2  4   0.000794   0.000824    0.000346
#> rdr3  4   0.000786   0.001009    0.000354
#> rdr4  4   0.000153   0.000635    0.000265
#> 
#> 
#> $RRRC
#> $RRRC$FTests
#>             DF       MS FStat      p
#> Treatment  4.0 0.001897  3.47 0.0305
#> Error     16.8 0.000547    NA     NA
#> 
#> $RRRC$ciDiffTrt
#>            Estimate StdErr   DF       t  PrGTt   CILower  CIUpper
#> trt1-trt2 -0.005100 0.0165 16.8 -0.3084 0.7616 -0.040021  0.02982
#> trt1-trt3  0.035325 0.0165 16.8  2.1361 0.0477  0.000404  0.07025
#> trt1-trt4 -0.005412 0.0165 16.8 -0.3273 0.7475 -0.040334  0.02951
#> trt1-trt5  0.036775 0.0165 16.8  2.2238 0.0402  0.001854  0.07170
#> trt2-trt3  0.040425 0.0165 16.8  2.4445 0.0258  0.005504  0.07535
#> trt2-trt4 -0.000312 0.0165 16.8 -0.0189 0.9851 -0.035234  0.03461
#> trt2-trt5  0.041875 0.0165 16.8  2.5322 0.0216  0.006954  0.07680
#> trt3-trt4 -0.040737 0.0165 16.8 -2.4634 0.0249 -0.075659 -0.00582
#> trt3-trt5  0.001450 0.0165 16.8  0.0877 0.9312 -0.033471  0.03637
#> trt4-trt5  0.042187 0.0165 16.8  2.5511 0.0208  0.007266  0.07711
#> 
#> $RRRC$ciAvgRdrEachTrt
#>      Estimate StdErr    DF CILower CIUpper     Cov2
#> trt1    0.845 0.0286  5.46   0.774   0.917 0.000211
#> trt2    0.850 0.0199 27.72   0.809   0.891 0.000266
#> trt3    0.810 0.0266  7.04   0.747   0.873 0.000246
#> trt4    0.851 0.0294  5.40   0.777   0.925 0.000220
#> trt5    0.808 0.0258  6.78   0.747   0.870 0.000222
#> 
#> 
#> $FRRC
#> $FRRC$FTests
#>                 MS Chisq DF       p
#> Treatment 0.001897  13.6  4 0.00868
#> Error     0.000558    NA NA      NA
#> 
#> $FRRC$ciDiffTrt
#>            Estimate StdErr       z  PrGTz  CILower CIUpper
#> trt1-trt2 -0.005100 0.0167 -0.3054 0.7601 -0.03783  0.0276
#> trt1-trt3  0.035325 0.0167  2.1151 0.0344  0.00259  0.0681
#> trt1-trt4 -0.005412 0.0167 -0.3241 0.7459 -0.03815  0.0273
#> trt1-trt5  0.036775 0.0167  2.2019 0.0277  0.00404  0.0695
#> trt2-trt3  0.040425 0.0167  2.4204 0.0155  0.00769  0.0732
#> trt2-trt4 -0.000312 0.0167 -0.0187 0.9851 -0.03305  0.0324
#> trt2-trt5  0.041875 0.0167  2.5073 0.0122  0.00914  0.0746
#> trt3-trt4 -0.040737 0.0167 -2.4392 0.0147 -0.07347 -0.0080
#> trt3-trt5  0.001450 0.0167  0.0868 0.9308 -0.03128  0.0342
#> trt4-trt5  0.042187 0.0167  2.5260 0.0115  0.00945  0.0749
#> 
#> $FRRC$ciAvgRdrEachTrt
#>      Estimate StdErr  DF CILower CIUpper
#> trt1    0.845 0.0183 199   0.809   0.881
#> trt2    0.850 0.0197 199   0.812   0.889
#> trt3    0.810 0.0201 199   0.770   0.849
#> trt4    0.851 0.0186 199   0.814   0.887
#> trt5    0.808 0.0197 199   0.770   0.847
#> 
#> $FRRC$ciDiffTrtEachRdr
#>                 Estimate StdErr       z   PrGTz  CILower CIUpper
#> rdr1::trt1-trt2  0.04000 0.0308  1.2989 0.19400 -0.02036  0.1004
#> rdr1::trt1-trt3  0.09130 0.0308  2.9646 0.00303  0.03094  0.1517
#> rdr1::trt1-trt4  0.00190 0.0308  0.0617 0.95081 -0.05846  0.0623
#> rdr1::trt1-trt5  0.06285 0.0308  2.0408 0.04127  0.00249  0.1232
#> rdr1::trt2-trt3  0.05130 0.0308  1.6658 0.09576 -0.00906  0.1117
#> rdr1::trt2-trt4 -0.03810 0.0308 -1.2372 0.21603 -0.09846  0.0223
#> rdr1::trt2-trt5  0.02285 0.0308  0.7420 0.45811 -0.03751  0.0832
#> rdr1::trt3-trt4 -0.08940 0.0308 -2.9029 0.00370 -0.14976 -0.0290
#> rdr1::trt3-trt5 -0.02845 0.0308 -0.9238 0.35559 -0.08881  0.0319
#> rdr1::trt4-trt5  0.06095 0.0308  1.9791 0.04780  0.00059  0.1213
#> rdr2::trt1-trt2 -0.04650 0.0309 -1.5039 0.13260 -0.10710  0.0141
#> rdr2::trt1-trt3 -0.01815 0.0309 -0.5870 0.55719 -0.07875  0.0424
#> rdr2::trt1-trt4 -0.03330 0.0309 -1.0770 0.28147 -0.09390  0.0273
#> rdr2::trt1-trt5  0.02520 0.0309  0.8150 0.41505 -0.03540  0.0858
#> rdr2::trt2-trt3  0.02835 0.0309  0.9169 0.35918 -0.03225  0.0889
#> rdr2::trt2-trt4  0.01320 0.0309  0.4269 0.66943 -0.04740  0.0738
#> rdr2::trt2-trt5  0.07170 0.0309  2.3190 0.02040  0.01110  0.1323
#> rdr2::trt3-trt4 -0.01515 0.0309 -0.4900 0.62414 -0.07575  0.0454
#> rdr2::trt3-trt5  0.04335 0.0309  1.4021 0.16090 -0.01725  0.1039
#> rdr2::trt4-trt5  0.05850 0.0309  1.8921 0.05848 -0.00210  0.1191
#> rdr3::trt1-trt2 -0.00875 0.0362 -0.2418 0.80896 -0.07969  0.0622
#> rdr3::trt1-trt3  0.05900 0.0362  1.6302 0.10307 -0.01194  0.1299
#> rdr3::trt1-trt4  0.02310 0.0362  0.6383 0.52331 -0.04784  0.0940
#> rdr3::trt1-trt5  0.04060 0.0362  1.1218 0.26196 -0.03034  0.1115
#> rdr3::trt2-trt3  0.06775 0.0362  1.8719 0.06122 -0.00319  0.1387
#> rdr3::trt2-trt4  0.03185 0.0362  0.8800 0.37885 -0.03909  0.1028
#> rdr3::trt2-trt5  0.04935 0.0362  1.3635 0.17271 -0.02159  0.1203
#> rdr3::trt3-trt4 -0.03590 0.0362 -0.9919 0.32124 -0.10684  0.0350
#> rdr3::trt3-trt5 -0.01840 0.0362 -0.5084 0.61118 -0.08934  0.0525
#> rdr3::trt4-trt5  0.01750 0.0362  0.4835 0.62872 -0.05344  0.0884
#> rdr4::trt1-trt2 -0.00515 0.0272 -0.1893 0.84987 -0.05848  0.0482
#> rdr4::trt1-trt3  0.00915 0.0272  0.3363 0.73664 -0.04418  0.0625
#> rdr4::trt1-trt4 -0.01335 0.0272 -0.4907 0.62366 -0.06668  0.0400
#> rdr4::trt1-trt5  0.01845 0.0272  0.6781 0.49770 -0.03488  0.0718
#> rdr4::trt2-trt3  0.01430 0.0272  0.5256 0.59918 -0.03903  0.0676
#> rdr4::trt2-trt4 -0.00820 0.0272 -0.3014 0.76312 -0.06153  0.0451
#> rdr4::trt2-trt5  0.02360 0.0272  0.8674 0.38572 -0.02973  0.0769
#> rdr4::trt3-trt4 -0.02250 0.0272 -0.8270 0.40825 -0.07583  0.0308
#> rdr4::trt3-trt5  0.00930 0.0272  0.3418 0.73249 -0.04403  0.0626
#> rdr4::trt4-trt5  0.03180 0.0272  1.1688 0.24249 -0.02153  0.0851
#> 
#> $FRRC$IndividualRdrVarCov1
#>      varEachRdr cov1EachRdr
#> rdr1   0.000689    0.000215
#> rdr2   0.000824    0.000346
#> rdr3   0.001009    0.000354
#> rdr4   0.000635    0.000265
#> 
#> 
#> $RRFC
#> $RRFC$FTests
#>    DF       MS   F      p
#> T   4 0.001897 4.1 0.0253
#> TR 12 0.000462  NA     NA
#> 
#> $RRFC$ciDiffTrt
#>            Estimate StdErr DF       t  PrGTt  CILower  CIUpper
#> trt1-trt2 -0.005100 0.0152 12 -0.3355 0.7431 -0.03822  0.02802
#> trt1-trt3  0.035325 0.0152 12  2.3237 0.0385  0.00220  0.06845
#> trt1-trt4 -0.005412 0.0152 12 -0.3560 0.7280 -0.03854  0.02771
#> trt1-trt5  0.036775 0.0152 12  2.4191 0.0324  0.00365  0.06990
#> trt2-trt3  0.040425 0.0152 12  2.6592 0.0208  0.00730  0.07355
#> trt2-trt4 -0.000312 0.0152 12 -0.0206 0.9839 -0.03344  0.03281
#> trt2-trt5  0.041875 0.0152 12  2.7546 0.0175  0.00875  0.07500
#> trt3-trt4 -0.040737 0.0152 12 -2.6797 0.0200 -0.07386 -0.00761
#> trt3-trt5  0.001450 0.0152 12  0.0954 0.9256 -0.03167  0.03457
#> trt4-trt5  0.042187 0.0152 12  2.7751 0.0168  0.00906  0.07531
#> 
#> $RRFC$ciAvgRdrEachTrt
#>      Estimate StdErr DF CILower CIUpper
#> Trt1    0.845 0.0246  3   0.767   0.923
#> Trt2    0.850 0.0114  3   0.814   0.887
#> Trt3    0.810 0.0215  3   0.741   0.878
#> Trt4    0.851 0.0254  3   0.770   0.931
#> Trt5    0.808 0.0210  3   0.742   0.875
```

### Random-Reader Random-Case (RRRC) analysis {#or-applications-RRRC-dataset04}

* `st2$RRRC$FTests` contains the results of the F-test.
* In this example `ndf` = 4 because there are I = 5 treatments. Since the p-value is less than 0.05, at least one treatment-pairing FOM difference is significantly different from zero.

* `st2$RRRC$ciDiffTrt` contains the confidence intervals for the inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* With I = 5 treatments there are 10 distinct treatment-pairings. 
* Looking at the `PrGTt` (for probability greater than `t`) column, one finds six pairings that are significant: `trt1-trt3`, `trt1-trt5`, etc. The smallest p-value is for the `trt4-trt5` pairing. 

* `st2$RRRC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet}}$.
* Looking at the `Estimate` column one confirms that `trt5` has the smallest FOM while `trt4` has the highest.

### Fixed-Reader Random-Case (FRRC) analysis {#or-applications-FRRC-dataset04}

* `st2$FRRC$FTests` contains results of the F-tests, which in this situation is actually a chi-square test of the NH.
* Again, `ndf` = 4 because there are I = 5 treatments. Since the p-value is less than 0.05, at least one treatment-pairing FOM difference is significantly different from zero.

* `st2$FRRC$ciDiffTrt` contains confidence intervals for the inter-treatment paired difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* With I = 5 treatments there are 10 distinct treatment-pairings. 
* Looking at the `PrGTt` column, one finds six pairings that are significant: `trt1-trt3`, `trt1-trt5`, etc. The smallest p-value is for the `trt4-trt5` pairing. 

* `st2$FRRC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet}}$.
* The `Estimate` column confirms that `trt5` has the smallest FOM while `trt4` has the highest.

### Random-Reader Fixed-Case (RRFC) analysis {#or-applications-RRFC-dataset04}

* `st2$RRFC$FTests` contains the results of the F-test of the NH.
* Again, `ndf` = 4 because there are I = 5 treatments. Since the p-value is less than 0.05, at least one treatment-pairing FOM difference is significantly different from zero.

* `st2$RRFC$ciDiffTrt` contains confidence intervals for the inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* With I = 5 treatments there are 10 distinct treatment-pairings. 
* The `PrGTt` column shows that six pairings are significant: `trt1-trt3`, `trt1-trt5`, etc. The smallest p-value is for the `trt4-trt5` pairing. 

* `st2$RRFC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet}}$.
* The `Estimate` column confirms that `trt5` has the smallest FOM while `trt4` has the highest (the `Estimates` column is identical for RRRC, FRRC and RRFC analyses).

## RJafroc: dataset04, FROC {#or-applications-dataset04-FROC-RJafroc}
* The fourth example uses `dataset04`, but this time we use the FROC data, specifically, we do not convert it to inferred-ROC. 
* Since this is an FROC dataset, one needs to use an FROC figure of merit. 
* In this example the weighted AFROC figure of merit `FOM = "wAFROC"` is specified. This is the recommended figure of merit when both normal and abnormal cases are present in the dataset.
* If the dataset does not contain normal cases, then the weighted AFROC1 figure of merit `FOM = "wAFROC1"` should be specified. 
* The results are contained in `st3`. 
* As noted earlier, this time the object is listed in its entirety.


```r
ds <- dataset04 # do NOT convert to ROC
FOM <- "wAFROC" 
st3 <- StSignificanceTesting(ds, FOM = FOM, method = "OR")
print(st3, digits = 3)
#> $FOMs
#> $FOMs$foms
#>       rdr1  rdr3  rdr4  rdr5
#> trt1 0.779 0.725 0.704 0.805
#> trt2 0.787 0.727 0.723 0.804
#> trt3 0.730 0.716 0.672 0.773
#> trt4 0.810 0.743 0.694 0.829
#> trt5 0.749 0.682 0.655 0.771
#> 
#> $FOMs$trtMeans
#>      Estimate
#> trt1    0.753
#> trt2    0.760
#> trt3    0.723
#> trt4    0.769
#> trt5    0.714
#> 
#> $FOMs$trtMeanDiffs
#>           Estimate
#> trt1-trt2 -0.00686
#> trt1-trt3  0.03061
#> trt1-trt4 -0.01604
#> trt1-trt5  0.03884
#> trt2-trt3  0.03747
#> trt2-trt4 -0.00918
#> trt2-trt5  0.04570
#> trt3-trt4 -0.04665
#> trt3-trt5  0.00823
#> trt4-trt5  0.05488
#> 
#> 
#> $ANOVA
#> $ANOVA$TRanova
#>         SS DF      MS
#> T  0.00927  4 0.00232
#> R  0.03540  3 0.01180
#> TR 0.00204 12 0.00017
#> 
#> $ANOVA$VarCom
#>       Estimates  Rhos
#> VarR   0.002209    NA
#> VarTR -0.000305    NA
#> Cov1   0.000422 0.455
#> Cov2   0.000336 0.362
#> Cov3   0.000304 0.328
#> Var    0.000928    NA
#> 
#> $ANOVA$IndividualTrt
#>      DF msREachTrt varEachTrt cov2EachTrt
#> trt1  3    0.00221   0.000877    0.000333
#> trt2  3    0.00171   0.000939    0.000380
#> trt3  3    0.00171   0.000970    0.000297
#> trt4  3    0.00386   0.000859    0.000311
#> trt5  3    0.00298   0.000995    0.000359
#> 
#> $ANOVA$IndividualRdr
#>      DF msTEachRdr varEachRdr cov1EachRdr
#> rdr1  4   0.001014   0.000883    0.000412
#> rdr3  4   0.000509   0.000897    0.000436
#> rdr4  4   0.000698   0.001171    0.000495
#> rdr5  4   0.000604   0.000762    0.000345
#> 
#> 
#> $RRRC
#> $RRRC$FTests
#>             DF       MS FStat        p
#> Treatment  4.0 0.002317   7.8 0.000117
#> Error     36.8 0.000297    NA       NA
#> 
#> $RRRC$ciDiffTrt
#>           Estimate StdErr   DF      t    PrGTt  CILower  CIUpper
#> trt1-trt2 -0.00686 0.0122 36.8 -0.563 5.77e-01 -0.03155  0.01784
#> trt1-trt3  0.03061 0.0122 36.8  2.512 1.65e-02  0.00592  0.05531
#> trt1-trt4 -0.01604 0.0122 36.8 -1.316 1.96e-01 -0.04073  0.00866
#> trt1-trt5  0.03884 0.0122 36.8  3.188 2.92e-03  0.01415  0.06354
#> trt2-trt3  0.03747 0.0122 36.8  3.075 3.96e-03  0.01278  0.06217
#> trt2-trt4 -0.00918 0.0122 36.8 -0.753 4.56e-01 -0.03387  0.01552
#> trt2-trt5  0.04570 0.0122 36.8  3.750 6.07e-04  0.02100  0.07040
#> trt3-trt4 -0.04665 0.0122 36.8 -3.828 4.85e-04 -0.07135 -0.02195
#> trt3-trt5  0.00823 0.0122 36.8  0.675 5.04e-01 -0.01647  0.03292
#> trt4-trt5  0.05488 0.0122 36.8  4.504 6.52e-05  0.03018  0.07957
#> 
#> $RRRC$ciAvgRdrEachTrt
#>      Estimate StdErr    DF CILower CIUpper     Cov2
#> trt1    0.753 0.0298  7.71   0.684   0.822 0.000333
#> trt2    0.760 0.0284 10.69   0.697   0.823 0.000380
#> trt3    0.723 0.0269  8.62   0.661   0.784 0.000297
#> trt4    0.769 0.0357  5.24   0.679   0.860 0.000311
#> trt5    0.714 0.0333  6.59   0.635   0.794 0.000359
#> 
#> 
#> $FRRC
#> $FRRC$FTests
#>                 MS Chisq DF       p
#> Treatment 0.002317  15.4  4 0.00393
#> Error     0.000602    NA NA      NA
#> 
#> $FRRC$ciDiffTrt
#>           Estimate StdErr      z   PrGTz  CILower CIUpper
#> trt1-trt2 -0.00686 0.0173 -0.395 0.69260 -0.04085  0.0271
#> trt1-trt3  0.03061 0.0173  1.765 0.07753 -0.00338  0.0646
#> trt1-trt4 -0.01604 0.0173 -0.925 0.35518 -0.05003  0.0180
#> trt1-trt5  0.03884 0.0173  2.240 0.02511  0.00485  0.0728
#> trt2-trt3  0.03747 0.0173  2.161 0.03073  0.00348  0.0715
#> trt2-trt4 -0.00918 0.0173 -0.529 0.59662 -0.04317  0.0248
#> trt2-trt5  0.04570 0.0173  2.635 0.00841  0.01171  0.0797
#> trt3-trt4 -0.04665 0.0173 -2.690 0.00715 -0.08064 -0.0127
#> trt3-trt5  0.00823 0.0173  0.474 0.63515 -0.02576  0.0422
#> trt4-trt5  0.05488 0.0173  3.164 0.00155  0.02089  0.0889
#> 
#> $FRRC$ciAvgRdrEachTrt
#>      Estimate StdErr  DF CILower CIUpper
#> trt1    0.753 0.0217 199   0.711   0.796
#> trt2    0.760 0.0228 199   0.715   0.805
#> trt3    0.723 0.0216 199   0.680   0.765
#> trt4    0.769 0.0212 199   0.728   0.811
#> trt5    0.714 0.0228 199   0.670   0.759
#> 
#> $FRRC$ciDiffTrtEachRdr
#>                 Estimate StdErr       z   PrGTz  CILower   CIUpper
#> rdr1::trt1-trt2 -0.00773 0.0307 -0.2520 0.80105 -0.06788  0.052416
#> rdr1::trt1-trt3  0.04957 0.0307  1.6154 0.10622 -0.01057  0.109724
#> rdr1::trt1-trt4 -0.03087 0.0307 -1.0058 0.31451 -0.09102  0.029282
#> rdr1::trt1-trt5  0.03047 0.0307  0.9928 0.32083 -0.02968  0.090616
#> rdr1::trt2-trt3  0.05731 0.0307  1.8674 0.06185 -0.00284  0.117457
#> rdr1::trt2-trt4 -0.02313 0.0307 -0.7538 0.45097 -0.08328  0.037016
#> rdr1::trt2-trt5  0.03820 0.0307  1.2448 0.21322 -0.02195  0.098349
#> rdr1::trt3-trt4 -0.08044 0.0307 -2.6212 0.00876 -0.14059 -0.020293
#> rdr1::trt3-trt5 -0.01911 0.0307 -0.6226 0.53352 -0.07926  0.041041
#> rdr1::trt4-trt5  0.06133 0.0307  1.9986 0.04566  0.00118  0.121482
#> rdr3::trt1-trt2 -0.00201 0.0304 -0.0661 0.94726 -0.06152  0.057504
#> rdr3::trt1-trt3  0.00913 0.0304  0.3008 0.76357 -0.05038  0.068646
#> rdr3::trt1-trt4 -0.01822 0.0304 -0.6002 0.54836 -0.07774  0.041287
#> rdr3::trt1-trt5  0.04262 0.0304  1.4035 0.16046 -0.01690  0.102129
#> rdr3::trt2-trt3  0.01114 0.0304  0.3669 0.71367 -0.04837  0.070654
#> rdr3::trt2-trt4 -0.01622 0.0304 -0.5341 0.59329 -0.07573  0.043296
#> rdr3::trt2-trt5  0.04462 0.0304  1.4697 0.14165 -0.01489  0.104137
#> rdr3::trt3-trt4 -0.02736 0.0304 -0.9010 0.36758 -0.08687  0.032154
#> rdr3::trt3-trt5  0.03348 0.0304  1.1027 0.27014 -0.02603  0.092996
#> rdr3::trt4-trt5  0.06084 0.0304  2.0037 0.04510  0.00133  0.120354
#> rdr4::trt1-trt2 -0.01899 0.0368 -0.5166 0.60543 -0.09104  0.053061
#> rdr4::trt1-trt3  0.03132 0.0368  0.8519 0.39429 -0.04074  0.103370
#> rdr4::trt1-trt4  0.00927 0.0368  0.2521 0.80099 -0.06279  0.081320
#> rdr4::trt1-trt5  0.04845 0.0368  1.3179 0.18753 -0.02360  0.120503
#> rdr4::trt2-trt3  0.05031 0.0368  1.3685 0.17116 -0.02174  0.122361
#> rdr4::trt2-trt4  0.02826 0.0368  0.7687 0.44209 -0.04379  0.100311
#> rdr4::trt2-trt5  0.06744 0.0368  1.8345 0.06658 -0.00461  0.139495
#> rdr4::trt3-trt4 -0.02205 0.0368 -0.5998 0.54864 -0.09410  0.050003
#> rdr4::trt3-trt5  0.01713 0.0368  0.4661 0.64118 -0.05492  0.089186
#> rdr4::trt4-trt5  0.03918 0.0368  1.0659 0.28649 -0.03287  0.111236
#> rdr5::trt1-trt2  0.00131 0.0289  0.0453 0.96385 -0.05526  0.057881
#> rdr5::trt1-trt3  0.03243 0.0289  1.1237 0.26116 -0.02414  0.089006
#> rdr5::trt1-trt4 -0.02432 0.0289 -0.8425 0.39953 -0.08089  0.032256
#> rdr5::trt1-trt5  0.03384 0.0289  1.1724 0.24102 -0.02273  0.090414
#> rdr5::trt2-trt3  0.03112 0.0289  1.0783 0.28089 -0.02545  0.087698
#> rdr5::trt2-trt4 -0.02563 0.0289 -0.8878 0.37466 -0.08220  0.030948
#> rdr5::trt2-trt5  0.03253 0.0289  1.1271 0.25969 -0.02404  0.089106
#> rdr5::trt3-trt4 -0.05675 0.0289 -1.9661 0.04929 -0.11332 -0.000177
#> rdr5::trt3-trt5  0.00141 0.0289  0.0488 0.96109 -0.05516  0.057981
#> rdr5::trt4-trt5  0.05816 0.0289  2.0149 0.04391  0.00159  0.114731
#> 
#> $FRRC$IndividualRdrVarCov1
#>      varEachRdr cov1EachRdr
#> rdr1   0.000883    0.000412
#> rdr3   0.000897    0.000436
#> rdr4   0.001171    0.000495
#> rdr5   0.000762    0.000345
#> 
#> 
#> $RRFC
#> $RRFC$FTests
#>    DF      MS    F        p
#> T   4 0.00232 13.7 0.000202
#> TR 12 0.00017   NA       NA
#> 
#> $RRFC$ciDiffTrt
#>           Estimate  StdErr DF      t    PrGTt CILower  CIUpper
#> trt1-trt2 -0.00686 0.00921 12 -0.745 4.71e-01 -0.0269  0.01321
#> trt1-trt3  0.03061 0.00921 12  3.324 6.06e-03  0.0106  0.05068
#> trt1-trt4 -0.01604 0.00921 12 -1.741 1.07e-01 -0.0361  0.00403
#> trt1-trt5  0.03884 0.00921 12  4.218 1.19e-03  0.0188  0.05891
#> trt2-trt3  0.03747 0.00921 12  4.069 1.56e-03  0.0174  0.05754
#> trt2-trt4 -0.00918 0.00921 12 -0.997 3.39e-01 -0.0292  0.01089
#> trt2-trt5  0.04570 0.00921 12  4.963 3.29e-04  0.0256  0.06576
#> trt3-trt4 -0.04665 0.00921 12 -5.066 2.77e-04 -0.0667 -0.02659
#> trt3-trt5  0.00823 0.00921 12  0.894 3.89e-01 -0.0118  0.02829
#> trt4-trt5  0.05488 0.00921 12  5.959 6.62e-05  0.0348  0.07494
#> 
#> $RRFC$ciAvgRdrEachTrt
#>      Estimate StdErr DF CILower CIUpper
#> Trt1    0.753 0.0235  3   0.678   0.828
#> Trt2    0.760 0.0207  3   0.694   0.826
#> Trt3    0.723 0.0207  3   0.657   0.788
#> Trt4    0.769 0.0311  3   0.670   0.868
#> Trt5    0.714 0.0273  3   0.627   0.801
```

### Random-Reader Random-Case (RRRC) analysis {#or-applications-RRRC-dataset04-FROC}

* `st3$RRRC$FTests` contains the results of the F-tests.
* The p-value is much smaller than that obtained after converting to an ROC dataset. Specifically, for FROC analysis, the p-value is 1.17105004\times 10^{-4} while that for ROC analysis is 0.03054456. The F-statistic and the `ddf` are both larger for FROC analysis, both of of which result in increased probability of rejecting the NH, i.e., FROC analysis has greater power than ROC analysis.
* The increased power of FROC analysis has been confirmed in simulation studies [@RN1331].

* `st3$RRRC$ciDiffTrt` contains the confidence intervals for the inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* With I = 5 treatments there are 10 distinct treatment-pairings. 
* Looking at the `PrGTt` (for probability greater than `t`) column, one finds six pairings that are significant: `trt1-trt3`, `trt1-trt5`, etc. The smallest p-value is for the `trt4-trt5` pairing. The findings are consistent with the prior ROC analysis, the difference being the smaller p-values. 

* `st3$RRRC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet}}$.
* Looking at the `Estimate` column one confirms that `trt5` has the smallest FOM while `trt4` has the highest (the `Estimates` column is identical for RRRC, FRRC and RRFC analyses).

* `st3$RRRC$st1$RRRC$ciDiffTrtEachRdr` contains confidence intervals for inter-treatment difference FOMs, for each reader, i.e., $CI_{1-\alpha,RRRC,\theta_{i j} - \theta_{i' j}}$.

### Fixed-Reader Random-Case (FRRC) analysis {#or-applications-FRRC-dataset04-FROC}

* `st3$FRRC$FTests` contains results of the F-test of the NH.
* Again, `ndf` = 4 because there are I = 5 treatments. Since the p-value is less than 0.05, at least one treatment-pairing FOM difference is significantly different from zero.

* `st3$FRRC$ciDiffTrt` contains the confidence intervals for the inter-treatment paired difference FOMs averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* With I = 5 treatments there are 10 distinct treatment-pairings. 
* Looking at the `PrGTt` (for probability greater than `t`) column, one finds six pairings that are significant: `trt1-trt3`, `trt1-trt5`, etc. The smallest p-value is for the `trt4-trt5` pairing. The findings are consistent with the prior ROC analysis, the difference being the smaller p-values. 

* `st3$FRRC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet}}$.
* Looking at the `Estimate` column one confirms that `trt5` has the smallest FOM while `trt4` has the highest.

* `st3$FRRC$st1$FRRC$ciDiffTrtEachRdr` contains confidence intervals for inter-treatment difference FOMs, for each reader, i.e., $CI_{1-\alpha,FRRC,\theta_{i j} - \theta_{i' j}}$.

### Random-Reader Fixed-Case (RRFC) analysis {#or-applications-RRFC-dataset04-FROC}

* `st3$RRFC$FTests` contains results of the F-test of the NH.
* Again, `ndf` = 4 because there are I = 5 treatments. Since the p-value is less than 0.05, at least one treatment-pairing FOM difference is significantly different from zero.

* `st3$RRFC$ciDiffTrt` contains confidence intervals for the inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet} - \theta_{i' \bullet}}$.

* `st3$RRFC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet}}$.
* The `Estimate` column confirms that `trt5` has the smallest FOM while `trt4` has the highest (the `Estimates` column is identical for RRRC, FRRC and RRFC analyses).

## RJafroc: dataset04, FROC/DBM {#or-applications-dataset04-FROC-DBM-RJafroc}
* The fourth example again uses `dataset04`, i.e., FROC data, *but this time using DBM analysis*.
* The key difference below is in the call to `StSignificanceTesting()` function, where we set `method = "DBM"`.
* Since DBM analysis is pseudovalue based, and the figure of merit is not the empirical AUC under the ROC, one expects to see differences from the previously presented OR analysis, contained in `st3`.


```r
st4 <- StSignificanceTesting(ds, FOM = FOM, method = "DBM") 
# Note: using DBM analysis
print(st4, digits = 3)
#> $FOMs
#> $FOMs$foms
#>       rdr1  rdr3  rdr4  rdr5
#> trt1 0.779 0.725 0.704 0.805
#> trt2 0.787 0.727 0.723 0.804
#> trt3 0.730 0.716 0.672 0.773
#> trt4 0.810 0.743 0.694 0.829
#> trt5 0.749 0.682 0.655 0.771
#> 
#> $FOMs$trtMeans
#>      Estimate
#> trt1    0.753
#> trt2    0.760
#> trt3    0.723
#> trt4    0.769
#> trt5    0.714
#> 
#> $FOMs$trtMeanDiffs
#>           Estimate
#> trt1-trt2 -0.00686
#> trt1-trt3  0.03061
#> trt1-trt4 -0.01604
#> trt1-trt5  0.03884
#> trt2-trt3  0.03747
#> trt2-trt4 -0.00918
#> trt2-trt5  0.04570
#> trt3-trt4 -0.04665
#> trt3-trt5  0.00823
#> trt4-trt5  0.05488
#> 
#> 
#> $ANOVA
#> $ANOVA$TRCanova
#>            SS   DF     MS
#> T       1.853    4 0.4633
#> R       7.081    3 2.3603
#> C     289.602  199 1.4553
#> TR      0.407   12 0.0339
#> TC     95.772  796 0.1203
#> RC    126.902  597 0.2126
#> TRC   226.479 2388 0.0948
#> Total 748.096 3999     NA
#> 
#> $ANOVA$VarCom
#>        Estimates
#> VarR    0.002209
#> VarC    0.060862
#> VarTR  -0.000305
#> VarTC   0.006369
#> VarRC   0.023545
#> VarErr  0.094841
#> 
#> $ANOVA$IndividualTrt
#>       DF  Trt1  Trt2  Trt3  Trt4  Trt5
#> msR    3 0.442 0.343 0.342 0.772 0.597
#> msC  199 0.375 0.416 0.372 0.358 0.415
#> msRC 597 0.109 0.112 0.134 0.110 0.127
#> 
#> $ANOVA$IndividualRdr
#>       DF   rdr1   rdr3  rdr4   rdr5
#> msT    4 0.2027 0.1019 0.140 0.1208
#> msC  199 0.5064 0.5278 0.630 0.4285
#> msTC 796 0.0942 0.0922 0.135 0.0833
#> 
#> 
#> $RRRC
#> $RRRC$FTests
#>             DF     MS FStat        p
#> Treatment  4.0 0.4633   7.8 0.000117
#> Error     36.8 0.0594    NA       NA
#> 
#> $RRRC$ciDiffTrt
#>           Estimate StdErr   DF      t    PrGTt  CILower  CIUpper
#> trt1-trt2 -0.00686 0.0122 36.8 -0.563 5.77e-01 -0.03155  0.01784
#> trt1-trt3  0.03061 0.0122 36.8  2.512 1.65e-02  0.00592  0.05531
#> trt1-trt4 -0.01604 0.0122 36.8 -1.316 1.96e-01 -0.04073  0.00866
#> trt1-trt5  0.03884 0.0122 36.8  3.188 2.92e-03  0.01415  0.06354
#> trt2-trt3  0.03747 0.0122 36.8  3.075 3.96e-03  0.01278  0.06217
#> trt2-trt4 -0.00918 0.0122 36.8 -0.753 4.56e-01 -0.03387  0.01552
#> trt2-trt5  0.04570 0.0122 36.8  3.750 6.07e-04  0.02100  0.07040
#> trt3-trt4 -0.04665 0.0122 36.8 -3.828 4.85e-04 -0.07135 -0.02195
#> trt3-trt5  0.00823 0.0122 36.8  0.675 5.04e-01 -0.01647  0.03292
#> trt4-trt5  0.05488 0.0122 36.8  4.504 6.52e-05  0.03018  0.07957
#> 
#> $RRRC$ciAvgRdrEachTrt
#>      Estimate StdErr    DF CILower CIUpper
#> trt1    0.753 0.0298  7.71   0.684   0.822
#> trt2    0.760 0.0284 10.69   0.697   0.823
#> trt3    0.723 0.0269  8.62   0.661   0.784
#> trt4    0.769 0.0357  5.24   0.679   0.860
#> trt5    0.714 0.0333  6.59   0.635   0.794
#> 
#> 
#> $FRRC
#> $FRRC$FTests
#>            DF    MS FStat       p
#> Treatment   4 0.463  3.85 0.00416
#> Error     796 0.120    NA      NA
#> 
#> $FRRC$ciDiffTrt
#>           Estimate StdErr  DF      t   PrGTt  CILower CIUpper
#> trt1-trt2 -0.00686 0.0173 796 -0.395 0.69271 -0.04090  0.0272
#> trt1-trt3  0.03061 0.0173 796  1.765 0.07791 -0.00343  0.0647
#> trt1-trt4 -0.01604 0.0173 796 -0.925 0.35546 -0.05008  0.0180
#> trt1-trt5  0.03884 0.0173 796  2.240 0.02539  0.00480  0.0729
#> trt2-trt3  0.03747 0.0173 796  2.161 0.03103  0.00343  0.0715
#> trt2-trt4 -0.00918 0.0173 796 -0.529 0.59677 -0.04322  0.0249
#> trt2-trt5  0.04570 0.0173 796  2.635 0.00858  0.01166  0.0797
#> trt3-trt4 -0.04665 0.0173 796 -2.690 0.00730 -0.08069 -0.0126
#> trt3-trt5  0.00823 0.0173 796  0.474 0.63528 -0.02581  0.0423
#> trt4-trt5  0.05488 0.0173 796  3.164 0.00161  0.02084  0.0889
#> 
#> $FRRC$ciAvgRdrEachTrt
#>      Estimate StdErr  DF CILower CIUpper
#> trt1    0.753 0.0217 199   0.711   0.796
#> trt2    0.760 0.0228 199   0.715   0.805
#> trt3    0.723 0.0216 199   0.680   0.765
#> trt4    0.769 0.0212 199   0.728   0.811
#> trt5    0.714 0.0228 199   0.669   0.759
#> 
#> $FRRC$ciDiffTrtEachRdr
#>                 Estimate StdErr  DF       t   PrGTt   CILower   CIUpper
#> rdr1::trt1-trt2 -0.00773 0.0307 199 -0.2520 0.80131 -0.068250  0.052784
#> rdr1::trt1-trt3  0.04957 0.0307 199  1.6154 0.10781 -0.010942  0.110092
#> rdr1::trt1-trt4 -0.03087 0.0307 199 -1.0058 0.31573 -0.091384  0.029650
#> rdr1::trt1-trt5  0.03047 0.0307 199  0.9928 0.32203 -0.030050  0.090984
#> rdr1::trt2-trt3  0.05731 0.0307 199  1.8674 0.06332 -0.003209  0.117825
#> rdr1::trt2-trt4 -0.02313 0.0307 199 -0.7538 0.45186 -0.083650  0.037384
#> rdr1::trt2-trt5  0.03820 0.0307 199  1.2448 0.21469 -0.022317  0.098717
#> rdr1::trt3-trt4 -0.08044 0.0307 199 -2.6212 0.00944 -0.140959 -0.019925
#> rdr1::trt3-trt5 -0.01911 0.0307 199 -0.6226 0.53423 -0.079625  0.041409
#> rdr1::trt4-trt5  0.06133 0.0307 199  1.9986 0.04702  0.000816  0.121850
#> rdr3::trt1-trt2 -0.00201 0.0304 199 -0.0661 0.94733 -0.061885  0.057868
#> rdr3::trt1-trt3  0.00913 0.0304 199  0.3008 0.76389 -0.050743  0.069010
#> rdr3::trt1-trt4 -0.01822 0.0304 199 -0.6002 0.54904 -0.078102  0.041652
#> rdr3::trt1-trt5  0.04262 0.0304 199  1.4035 0.16202 -0.017260  0.102493
#> rdr3::trt2-trt3  0.01114 0.0304 199  0.3669 0.71406 -0.048735  0.071018
#> rdr3::trt2-trt4 -0.01622 0.0304 199 -0.5341 0.59389 -0.076093  0.043660
#> rdr3::trt2-trt5  0.04462 0.0304 199  1.4697 0.14323 -0.015252  0.104502
#> rdr3::trt3-trt4 -0.02736 0.0304 199 -0.9010 0.36867 -0.087235  0.032518
#> rdr3::trt3-trt5  0.03348 0.0304 199  1.1027 0.27148 -0.026393  0.093360
#> rdr3::trt4-trt5  0.06084 0.0304 199  2.0037 0.04645  0.000965  0.120718
#> rdr4::trt1-trt2 -0.01899 0.0368 199 -0.5166 0.60600 -0.091485  0.053502
#> rdr4::trt1-trt3  0.03132 0.0368 199  0.8519 0.39531 -0.041177  0.103810
#> rdr4::trt1-trt4  0.00927 0.0368 199  0.2521 0.80125 -0.063227  0.081760
#> rdr4::trt1-trt5  0.04845 0.0368 199  1.3179 0.18904 -0.024044  0.120944
#> rdr4::trt2-trt3  0.05031 0.0368 199  1.3685 0.17271 -0.022185  0.122802
#> rdr4::trt2-trt4  0.02826 0.0368 199  0.7687 0.44300 -0.044235  0.100752
#> rdr4::trt2-trt5  0.06744 0.0368 199  1.8345 0.06807 -0.005052  0.139935
#> rdr4::trt3-trt4 -0.02205 0.0368 199 -0.5998 0.54932 -0.094544  0.050444
#> rdr4::trt3-trt5  0.01713 0.0368 199  0.4661 0.64168 -0.055360  0.089627
#> rdr4::trt4-trt5  0.03918 0.0368 199  1.0659 0.28778 -0.033310  0.111677
#> rdr5::trt1-trt2  0.00131 0.0289 199  0.0453 0.96389 -0.055610  0.058227
#> rdr5::trt1-trt3  0.03243 0.0289 199  1.1237 0.26251 -0.024485  0.089352
#> rdr5::trt1-trt4 -0.02432 0.0289 199 -0.8425 0.40055 -0.081235  0.032602
#> rdr5::trt1-trt5  0.03384 0.0289 199  1.1724 0.24242 -0.023077  0.090760
#> rdr5::trt2-trt3  0.03112 0.0289 199  1.0783 0.28219 -0.025794  0.088044
#> rdr5::trt2-trt4 -0.02563 0.0289 199 -0.8878 0.37573 -0.082544  0.031294
#> rdr5::trt2-trt5  0.03253 0.0289 199  1.1271 0.26105 -0.024385  0.089452
#> rdr5::trt3-trt4 -0.05675 0.0289 199 -1.9661 0.05068 -0.113669  0.000169
#> rdr5::trt3-trt5  0.00141 0.0289 199  0.0488 0.96113 -0.055510  0.058327
#> rdr5::trt4-trt5  0.05816 0.0289 199  2.0149 0.04526  0.001240  0.115077
#> 
#> 
#> $RRFC
#> $RRFC$FTests
#>           DF     MS FStat        p
#> Treatment  4 0.4633  13.7 0.000202
#> Error     12 0.0339    NA       NA
#> 
#> $RRFC$ciDiffTrt
#>           Estimate  StdErr DF      t    PrGTt CILower  CIUpper
#> trt1-trt2 -0.00686 0.00921 12 -0.745 4.71e-01 -0.0269  0.01321
#> trt1-trt3  0.03061 0.00921 12  3.324 6.06e-03  0.0106  0.05068
#> trt1-trt4 -0.01604 0.00921 12 -1.741 1.07e-01 -0.0361  0.00403
#> trt1-trt5  0.03884 0.00921 12  4.218 1.19e-03  0.0188  0.05891
#> trt2-trt3  0.03747 0.00921 12  4.069 1.56e-03  0.0174  0.05754
#> trt2-trt4 -0.00918 0.00921 12 -0.997 3.39e-01 -0.0292  0.01089
#> trt2-trt5  0.04570 0.00921 12  4.963 3.29e-04  0.0256  0.06576
#> trt3-trt4 -0.04665 0.00921 12 -5.066 2.77e-04 -0.0667 -0.02659
#> trt3-trt5  0.00823 0.00921 12  0.894 3.89e-01 -0.0118  0.02829
#> trt4-trt5  0.05488 0.00921 12  5.959 6.62e-05  0.0348  0.07494
#> 
#> $RRFC$ciAvgRdrEachTrt
#>      Estimate StdErr DF CILower CIUpper
#> trt1    0.753 0.0235  3   0.678   0.828
#> trt2    0.760 0.0207  3   0.694   0.826
#> trt3    0.723 0.0207  3   0.657   0.788
#> trt4    0.769 0.0311  3   0.670   0.868
#> trt5    0.714 0.0273  3   0.627   0.801
```

### Random-Reader Random-Case (RRRC) analysis {#or-applications-RRRC-dataset04-FROC-DBM}

* `st4$RRRC$FTests` contains the results of the F-test of the NH.

* `st4$RRRC$ciDiffTrt` contains the confidence intervals for the inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.

* `st4$RRRC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRRC,\theta_{i \bullet}}$.

### Fixed-Reader Random-Case (FRRC) analysis {#or-applications-FRRC-dataset04-FROC-DBM}

* `st4$FRRC$FTests` contains results of the F-test of the NH, which is actually a chi-square statistic.

* `st4$FRRC$ciDiffTrt` contains confidence intervals for the inter-treatment difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet} - \theta_{i' \bullet}}$.
* With I = 5 treatments there are 10 distinct treatment-pairings. 
* Looking at the `PrGTt` (for probability greater than `t`) column, one finds six pairings that are significant: `trt1-trt3`, `trt1-trt5`, etc. The smallest p-value is for the `trt4-trt5` pairing. The findings are consistent with the prior ROC analysis, the difference being the smaller p-values. 

* `st4$FRRC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,FRRC,\theta_{i \bullet}}$.

* `st4$FRRC$ciDiffTrtEachRdr` contains confidence intervals for inter-treatment difference FOMs, for each reader, i.e., $CI_{1-\alpha,FRRC,\theta_{i j} - \theta_{i' j}}$.

### Random-Reader Fixed-Case (RRFC) analysis {#or-applications-RRFC-dataset04-FROC-DBM}

* `st4$RRFC$FTests` contains the results of the F-test of the NH.

* `st4$RRFC$ciDiffTrt` contains confidence intervals for the inter-treatment paired difference FOMs, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet} - \theta_{i' \bullet}}$.

* `st4$RRFC$ciAvgRdrEachTrt` contains confidence intervals for each treatment, averaged over readers, i.e., $CI_{1-\alpha,RRFC,\theta_{i \bullet}}$.
* The `Estimate` column confirms that `trt5` has the smallest FOM while `trt4` has the highest (the `Estimates` column is identical for RRRC, FRRC and RRFC analyses).


## Summary{#or-applications-Summary}
## Discussion{#or-applications-Discussion}
## Tentative {#ToMullOver1-tentative}


```r
ds1 <- dataset04 # do NOT convert to ROC
# comment/uncomment following code to disable/enable unequal weights
# K2 <- length(ds1$ratings$LL[1,1,,1])
# weights <- array(dim = c(K2, max(ds1$lesions$perCase)))
# perCase <- ds1$lesions$perCase
# for (k2 in 1:K2) {
#   sum <- 0
#   for (el in 1:perCase[k2]) {
#     weights[k2,el] <- 1/el
#     sum <- sum + 1/el
#   }
#   weights[k2,1:perCase[k2]] <- weights[k2,1:perCase[k2]] / sum
# }
# ds1$lesions$weights <- weights
ds <- ds1
FOM <- "wAFROC" # also try wAFROC1, MaxLLF and MaxNLF
st5 <- StSignificanceTesting(ds, FOM = FOM, method = "OR")
print(st5, digits = 4)
```

A comparison was run between results of OR and DBM for the FROC dataset. Except for `FRRC`, where differences are expected (because `ddf` in the former is $\infty$, while that in the later is $(I-1)\times(J-1))$, the results for the p-values were identical. This was true for the following FOMs: `wAFROC`, with equal and unequal weights, and `MaxLLF`. The confidence intervals (again, excluding `FRRC`) were identical for `FOM` = `wAFROC`. Slight differences were observed for `FOM` = `MaxLLF`.  

## References {#or-applications-references}

