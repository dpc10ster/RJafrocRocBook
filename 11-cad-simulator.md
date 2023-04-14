# CAD simulator {#cad-simulator}



## Introduction

Single-modality multiple-reader decision variable simulator

CAD and a group of radiologists interpreting the same cases is a single-modality multiple-reader paradigm. It is necessary to calibrate the simulator to an actual clinical dataset, i.e., assign realistic values to the parameters of the simulation model. A calibrated simulator will yield samples not identical, but statistically identical to the original dataset. The criteria used to assess statistical identity, is described later. Traditional Roe-Metz simulator, involving univariate sampling, is described in the online appendix material. This chapter describes a new method, involving multivariate sampling, were used.


### Calibrating the simulator to a specific dataset

The multivariate simulator was calibrated to the dataset used to illustrate the analysis8 in Chapter 22, i.e., the Hupse-Karssemeijer dataset. The study used the location receiver operating characteristic (LROC) paradigm. For the purpose of the current chapter, the indicated location was ignored, i.e., only the ROC rating is used. 

#### Calibrating parameters to individual radiologists

For each radiologist  ,  , one fits the clinical data using CBM , yielding  , and a 2 x 2 covariance matrix  . In the study in question J = 9. Threshold related covariance terms are excluded from this and the following description. Since each radiologist is an equally valid sample from the population, for parsimony, the values were averaged over  , obtaining  , where the dot symbol represents an average over the replaced index. For clarity, these are henceforth denoted  , where the "1R" subscript denotes single radiologist-based estimates, i.e., 

 	. 	(23.5)

In Eqn. (23.5),   is the CBM separation parameter when readers are analyzed individually using the CBM method and the parameters averaged. Similarly, the parameter   is the average CBM disease-visibility parameter. These are the calibrated CBM parameters for the individual readers.

#### Calibrating parameters to paired radiologists 

Reader index 0 is reserved for CAD. For each pairing of different radiologists in the clinical dataset,  , one fits the data using CORCBM1 yielding  , and the 6 x 6 covariance matrix  . The first index on each   is the case-truth index  . Since for each pairing  there exists a pairing  , upon averaging over all pairings the results must be symmetric with respect to   and  , yielding    and the average covariance matrix  . This is cast in more transparent notation as follows: 

 	. 	(23.6)

The "2R" and "RR" subscripts denote that the paired-radiologist fitting-method was used to estimate each parameter. Since different algorithms are used, one does not expect estimate   to equal  .  

#### Reconciling single radiologist and paired radiologist parameters

To simultaneously simulate all readers according to a multivariate distribution, see below, one needs a unique value for the radiologist separation parameter. This is accomplished by averaging the two estimates:

 	. 	(23.7)

A similar averaging is applied to the CBM visibility parameter:

 	.	(23.8)

Therefore, the final calibrated CBM parameters for paired radiologists are:

 	.	(23.9)

#### Calibrating parameters to CAD

Like the calibration of individual radiologists, the individual performance of CAD (regarded as reader "0") can be estimated using CBM model, yielding: 

 	. 	(23.10)

#### CAD-Radiologist pairings

Upon pairing CAD with radiologist  , CORCBM yields  and  . Here   represents the estimate of the CBM separation parameter when CAD is paired with radiologist   and likewise for the other model parameters. As before, one averages over the radiologist index, yielding    and  , where the subscript "CR" denotes a CAD vs. radiologist pairing. To reconcile parameters, the CAD based estimates of the non-correlation parameters are averaged with the corresponding CAD vs. radiologist values, yielding:

 	.	(23.11)

This can be written in more transparent notation as follows (C = CAD; CR = CAD vs. Radiologist): 

 	. 	(23.12)

Therefore, the final calibrated parameters for CAD paired with radiologist data are: 

 	. 	(23.13)

Eqn. (23.5), Eqn. (23.9) and Eqn. (23.13) complete the calibration process. Here   is the separation parameter for CAD while   is the separation parameter for the average radiologist (when measured paired with CAD). Likewise,   is the disease-visibility parameter for CAD while   is the corresponding parameter for the average radiologist (when measured paired with CAD). Under the null hypothesis the CAD and the radiologist parameters must be identical:

 	.	(23.14)

### Simulating data using the calibrated simulator

Let   denote the number of radiologists, non-diseased and diseased cases, respectively, to be simulated, and CAD is reader "0". For each pair of different readers  ,  , where the range of either   or   is 0 to J, each sample from the multivariate normal distribution yields a row-vector of 6 – parameters (C = CAD, R = radiologist):

 	.	(23.15)

The two parameters, characterizing the first (i.e., "X") arm of the pairing, are   and the corresponding ones characterizing the second (i.e., "Y") arm are  . The last two parameters   are the decision variable correlations for non-diseased and diseased cases, respectively, assuming, for diseased cases, that the disease is visible in both X- and Y-arms to the radiologist (intermediate visibility conditions will be accommodated later). 

•	In Eqn. (23.15), the upper line on the right hand side expresses the sampling when one of the readers is CAD,   or  . In this situation the other reader is the radiologist and one samples  . The first two samples,  , correspond to CAD and the next two correspond to the radiologist paired with CAD. 
•	In Eqn. (23.15), the lower line on the right hand side expresses the sampling when none of the readers is CAD,   and  . In this situation both readers are radiologist and one samples  . The first two samples,  , correspond to one radiologist and the next two correspond to another radiologist.
•	Because of the nature of random sampling, a parameter value for pairing   was not necessarily identical to that for  ; the different values were replaced by their average.

An average (and therefore more stable) value   for reader ,  , is obtained by averaging   over  , when   or averaging   over  , when  . This can be expressed by: 

 	.	(23.16) 

A similar expression applies to average value  : 

 	.	(23.17) 

One may wonder why the correlations do not have to be likewise averaged. This is because the correlation is between the two arms of the pairing, while the other parameters are specific to one or the other arm of the pairing.

#### Example with J = 3

Following is an example of generating Z-samples for CAD and three radiologists, J = 3 (four readers total). 

##### Non-diseased cases

For non-diseased cases, one need only construct the covariance matrix  . The subscript “1” indicates its truth state, and “J + 1” indicates its dimension along one "edge" of the matrix. All diagonal elements are set to unity.  The off-diagonal elements are filled as follows:

  	. 	(23.18)

The covariance matrix is shown below. The correlations were sampled in Eqn. (23.15). One is using the non-diseased case correlation   as the intent is to sample new non-diseased case z-samples.

 	. 	(23.19)

##### Diseased cases

For diseased cases, define the visibility-condition   (row) vector indexed by  . Each row contains   elements, each of which indicate whether disease was visible to the reader specified by the column number in the row-vector. If the disease was visible, the corresponding element is one and otherwise it is zero. For example, if the   element of   is one, then the disease was visible to reader  ,  . 

For J = 3, there are   visibility conditions, for which the visibility-condition vectors are shown in Table 23.1 under  . For example,  , means that in visibility-condition  , disease was visible to CAD (the first "one" in the array) and to all radiologists. Likewise,   means disease was invisible to CAD but was visible to all radiologists;   means disease was invisible to CAD and the first radiologist, but was visible to the other radiologists and finally,   means the disease was invisible to CAD or any of the radiologists.

 
Table 23.1: Shown in this table, for J = 3, are 16 visibility conditions. The visibility-condition row-vectors are shown in column 2. the corresponding   vector is in column= 2. The covariance matrices are listed in column 3 and the last column lists the probability of occurrence of each condition. All matrices are symmetric.
 

The corresponding   vectors follow from the second column in Table 23.1. For example, since  , it follows that . Since the disease was visible, the first value is the separation parameter for CAD, i.e.,  , the second value is the separation parameter for the first radiologist, i.e.,  , etc. The needed   were calculated via Eqn. (23.16). As another example, from  , it follows that . Since the disease was invisible to CAD and the first radiologist, the first two values are zeroes and the third value is the separation parameter for the second radiologist, i.e.,  , and the last value is the separation parameter for the third radiologist, i.e.,  . 

In Table 23.1, the covariance matrix   has ones along the diagonal (this is true of all covariance matrices in Table 23.1), the value at row 1 (CAD,   = 0) and column 2 (radiologist 1,   = 1) is  . One is using   because the disease was visible to both readers implied by row 1 and column 2, and the "01" indexing follows from the readers specified by the row and column indices. As another example, the value at row 2 (radiologist 1,   = 1) and column 3 (radiologist 2,   = 2) is  . One is using   because the disease was visible to both readers implied by the specified row and column indices. The "12" indexing follows from the readers implied by the specified row and column indices.

As a more complicated example, consider  . Since  , the disease was invisible to CAD or the first radiologist. This implies that the correlation at row 1 (  = 0) and column 2 (  = 1) must be  , subscripted with the appropriate reader indices, i.e., "01". The value at row 1 and column 3 is the average of   and  , i.e.,  , subscripted with the appropriate reader indices,  . As a final example, the value at row 3 and column 4 is  , because disease was visible to the second and third radiologists (  = 2,   = 3), and these are the appropriate indices. 

The last column of Table 23.1 lists the probabilities associate with each visibility-condition. For example, for visibility-condition   = 1,   because disease was visible to CAD, the probability of which is  , it was visible to radiologist 1, the probability of which is  , it was visible to radiologist 2, the probability of which is  and it was visible to radiologist 3, the probability of which is  . Because the radiologists are independent, the probability of the visibility-condition   = 1 is the product of the component probabilities, which is  . Similarly, the probability of observing visibility-condition   = 2 is  , because disease was invisible to CAD, the probability of which is  , etc. One can confirm that the probabilities listed in the last column of Table 23.1 sum to unity. 

To determine the number of diseased cases in each visibility-condition one samples the multinomial distribution with trial size   and cell probabilities specified by  .

##### Illustration using R code

The software implementing the method described above is implemented in a file, in the software directory, named mainCadVsRadCalibValidate.R. 

1.	Highlight lines 1 – 26 and click Run. Echoed lines are not shown below (warnings are ignorable). Enter q <- 1 in console window. Highlight lines 28 – 29 and click Run. Highlight lines 35– 91 and click Run. Highlight simuMu, line 85, and click Run. Highlight simuAlpha, line 86, and click Run. Highlight rhoNorMatr, line 88, and click Run. Highlight rhoAbn2Matr, line 89, and click Run. Highlight rhoNorMatr, line 88, and click Run. Highlight rhoAbn2Matr, line 89, and click Run. This should yield all values shown below under Part A. Highlight lines 96 – 111 and click Run. 

2.	Enter l <- 16 in console window. Highlight lines 113 – 131 and click Run. Highlight conditionArray, line 96, and click Run. Highlight probVector, line 101, and click Run. Highlight muTemp, line 114, and click Run. Highlight sigmaTemp, line 116, and click Run. This should yield all values shown below under Part B. This corresponds to c = 1 in Table 23.1. None of the elements of muTemp are zero, consistent with the fact that in this condition disease was visible to CAD and all radiologists.

3.	Enter l <- 14 and repeat #2. This corresponds to c = 2 in Table 23.1. The first element of muTemp is zero, consistent with the fact that in this condition disease was invisible to CAD.

4.	Enter l <- 13 and repeat #2. This corresponds to c = 3 in Table 23.1. The first two elements of muTemp are zero, consistent with the fact that in this condition disease was invisible to CAD and the first radiologist.

5.	Enter l <- 1 and repeat #2. This corresponds to c = 16 in Table 23.1. Note that all elements of muTemp are zero, consistent with the fact that in this condition disease was invisible to CAD and all radiologists.
23.3.2.1.3.1: Code output (partial)

Part A
seed = 1 , binning= TRUE , simuJ= 4 , K1= 120 , K2= 80
> simuMu
[1] 2.011985 2.061436 2.035396 2.092945
> simuAlpha
[1] 0.8879874 0.8728608 0.8820709 0.8531668
> rhoNorMatr
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.2430901 0.1093731 0.1802754
[2,] 0.2430901 1.0000000 0.5145082 0.5683869
[3,] 0.1093731 0.5145082 1.0000000 0.4721771
[4,] 0.1802754 0.5683869 0.4721771 1.0000000
> rhoAbn2Matr
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.7539593 0.4511450 0.7447770
[2,] 0.7539593 1.0000000 0.7295952 0.7648653
[3,] 0.4511450 0.7295952 1.0000000 0.7859280
[4,] 0.7447770 0.7648653 0.7859280 1.0000000


CAD and a group of radiologists interpreting the same cases is a single-modality multiple-reader paradigm. It is necessary to calibrate the simulator to an actual clinical dataset, i.e., assign realistic values to the parameters of the simulation model. A calibrated simulator will yield samples not identical, but statistically identical to the original dataset. The criteria used to assess statistical identity, is described later. Traditional Roe-Metz simulator, involving univariate sampling, is described in the online appendix material. This chapter describes a new method, involving multivariate sampling, were used. 


### Calibrating the simulator to a specific dataset

The multivariate simulator was calibrated to the dataset used to illustrate the analysis8 in Chapter 22, i.e., the Hupse-Karssemeijer dataset. The study used the location receiver operating characteristic (LROC) paradigm9. For the purpose of the current demonstration, the indicated location was ignored, i.e., only the ROC rating is used. 


#### Calibrating parameters to individual radiologists

For each radiologist  ,  , one fits the clinical data using CBM , yielding  , and a 2 x 2 covariance matrix  . In the study in question J = 9. Threshold related covariance terms are excluded from this and the following description. Since each radiologist is an equally valid sample from the population, for parsimony, the values were averaged over  , obtaining  , where the dot symbol represents an average over the replaced index. For clarity, these are henceforth denoted  , where the "1R" subscript denotes single radiologist-based estimates, i.e., 

 	. 	(23.5)

In Eqn. (23.5),   is the CBM separation parameter when readers are analyzed individually using the CBM method and the parameters averaged. Similarly, the parameter   is the average CBM disease-visibility parameter. These are the calibrated CBM parameters for the individual readers.


#### Calibrating parameters to paired radiologists 

Reader index 0 is reserved for CAD. For each pairing of different radiologists in the clinical dataset,  , one fits the data using CORCBM1 yielding  , and the 6 x 6 covariance matrix  . The first index on each   is the case-truth index  . Since for each pairing  there exists a pairing  , upon averaging over all pairings the results must be symmetric with respect to   and  , yielding    and the average covariance matrix  . This is cast in more transparent notation as follows: 

 	. 	(23.6)

The "2R" and "RR" subscripts denote that the paired-radiologist fitting-method was used to estimate each parameter. Since different algorithms are used, one does not expect estimate   to equal  .  


#### Reconciling single radiologist and paired radiologist parameters

To simultaneously simulate all readers according to a multivariate distribution, see below, one needs a unique value for the radiologist separation parameter. This is accomplished by averaging the two estimates:

 	. 	(23.7)

A similar averaging is applied to the CBM visibility parameter:

 	.	(23.8)

Therefore, the final calibrated CBM parameters for paired radiologists are:

 	.	(23.9)

#### Calibrating parameters to CAD

Like the calibration of individual radiologists, the individual performance of CAD (regarded as reader "0") can be estimated using CBM model, yielding: 

 	. 	(23.10)



#### CAD-Radiologist pairings

Upon pairing CAD with radiologist  , CORCBM yields  and  . Here   represents the estimate of the CBM separation parameter when CAD is paired with radiologist   and likewise for the other model parameters. As before, one averages over the radiologist index, yielding    and  , where the subscript "CR" denotes a CAD vs. radiologist pairing. To reconcile parameters, the CAD based estimates of the non-correlation parameters are averaged with the corresponding CAD vs. radiologist values, yielding:

 	.	(23.11)

This can be written in more transparent notation as follows (C = CAD; CR = CAD vs. Radiologist): 

 	. 	(23.12)

Therefore, the final calibrated parameters for CAD paired with radiologist data are: 

 	. 	(23.13)

Eqn. (23.5), Eqn. (23.9) and Eqn. (23.13) complete the calibration process. Here   is the separation parameter for CAD while   is the separation parameter for the average radiologist (when measured paired with CAD). Likewise,   is the disease-visibility parameter for CAD while   is the corresponding parameter for the average radiologist (when measured paired with CAD). Under the null hypothesis the CAD and the radiologist parameters must be identical:

 	.	(23.14)


### Simulating data using the calibrated simulator

Let   denote the number of radiologists, non-diseased and diseased cases, respectively, to be simulated, and CAD is reader "0". For each pair of different readers  ,  , where the range of either   or   is 0 to J, each sample from the multivariate normal distribution yields a row-vector of 6 – parameters (C = CAD, R = radiologist):

 	.	(23.15)

The two parameters, characterizing the first (i.e., "X") arm of the pairing, are   and the corresponding ones characterizing the second (i.e., "Y") arm are  . The last two parameters   are the decision variable correlations for non-diseased and diseased cases, respectively, assuming, for diseased cases, that the disease is visible in both X- and Y-arms to the radiologist (intermediate visibility conditions will be accommodated later). 

•	In Eqn. (23.15), the upper line on the right hand side expresses the sampling when one of the readers is CAD,   or  . In this situation the other reader is the radiologist and one samples  . The first two samples,  , correspond to CAD and the next two correspond to the radiologist paired with CAD. 
•	In Eqn. (23.15), the lower line on the right hand side expresses the sampling when none of the readers is CAD,   and  . In this situation both readers are radiologist and one samples  . The first two samples,  , correspond to one radiologist and the next two correspond to another radiologist.
•	Because of the nature of random sampling, a parameter value for pairing   was not necessarily identical to that for  ; the different values were replaced by their average.

An average (and therefore more stable) value   for reader ,  , is obtained by averaging   over  , when   or averaging   over  , when  . This can be expressed by: 

 	.	(23.16) 

A similar expression applies to average value  : 

 	.	(23.17) 

One may wonder why the correlations do not have to be likewise averaged. This is because the correlation is between the two arms of the pairing, while the other parameters are specific to one or the other arm of the pairing.


#### Example with J = 3

Following is an example of generating Z-samples for CAD and three radiologists, J = 3 (four readers total). 


##### Non-diseased cases

For non-diseased cases, one need only construct the covariance matrix  . The subscript “1” indicates its truth state, and “J + 1” indicates its dimension along one "edge" of the matrix. All diagonal elements are set to unity.  The off-diagonal elements are filled as follows:

  	. 	(23.18)

The covariance matrix is shown below. The correlations were sampled in Eqn. (23.15). One is using the non-diseased case correlation   as the intent is to sample new non-diseased case z-samples.

 	. 	(23.19)


##### Diseased cases

For diseased cases, define the visibility-condition   (row) vector indexed by  . Each row contains   elements, each of which indicate whether disease was visible to the reader specified by the column number in the row-vector. If the disease was visible, the corresponding element is one and otherwise it is zero. For example, if the   element of   is one, then the disease was visible to reader  ,  . 

For J = 3, there are   visibility conditions, for which the visibility-condition vectors are shown in Table 23.1 under  . For example,  , means that in visibility-condition  , disease was visible to CAD (the first "one" in the array) and to all radiologists. Likewise,   means disease was invisible to CAD but was visible to all radiologists;   means disease was invisible to CAD and the first radiologist, but was visible to the other radiologists and finally,   means the disease was invisible to CAD or any of the radiologists.

Table 23.1: Shown in this table, for J = 3, are 16 visibility conditions. The visibility-condition row-vectors are shown in column 2. the corresponding   vector is in column= 2. The covariance matrices are listed in column 3 and the last column lists the probability of occurrence of each condition. All matrices are symmetric.


The corresponding   vectors follow from the second column in Table 23.1. For example, since  , it follows that . Since the disease was visible, the first value is the separation parameter for CAD, i.e.,  , the second value is the separation parameter for the first radiologist, i.e.,  , etc. The needed   were calculated via Eqn. (23.16). As another example, from  , it follows that . Since the disease was invisible to CAD and the first radiologist, the first two values are zeroes and the third value is the separation parameter for the second radiologist, i.e.,  , and the last value is the separation parameter for the third radiologist, i.e.,  . 

In Table 23.1, the covariance matrix   has ones along the diagonal (this is true of all covariance matrices in Table 23.1), the value at row 1 (CAD,   = 0) and column 2 (radiologist 1,   = 1) is  . One is using   because the disease was visible to both readers implied by row 1 and column 2, and the "01" indexing follows from the readers specified by the row and column indices. As another example, the value at row 2 (radiologist 1,   = 1) and column 3 (radiologist 2,   = 2) is  . One is using   because the disease was visible to both readers implied by the specified row and column indices. The "12" indexing follows from the readers implied by the specified row and column indices.

As a more complicated example, consider  . Since  , the disease was invisible to CAD or the first radiologist. This implies that the correlation at row 1 (  = 0) and column 2 (  = 1) must be  , subscripted with the appropriate reader indices, i.e., "01". The value at row 1 and column 3 is the average of   and  , i.e.,  , subscripted with the appropriate reader indices,  . As a final example, the value at row 3 and column 4 is  , because disease was visible to the second and third radiologists (  = 2,   = 3), and these are the appropriate indices. 

The last column of Table 23.1 lists the probabilities associate with each visibility-condition. For example, for visibility-condition   = 1,   because disease was visible to CAD, the probability of which is  , it was visible to radiologist 1, the probability of which is  , it was visible to radiologist 2, the probability of which is  and it was visible to radiologist 3, the probability of which is  . Because the radiologists are independent, the probability of the visibility-condition   = 1 is the product of the component probabilities, which is  . Similarly, the probability of observing visibility-condition   = 2 is  , because disease was invisible to CAD, the probability of which is  , etc. One can confirm that the probabilities listed in the last column of Table 23.1 sum to unity. 

To determine the number of diseased cases in each visibility-condition one samples the multinomial distribution with trial size   and cell probabilities specified by  . 


##### Illustration using R code

The software implementing the method described above is implemented in a file, in the software directory, named mainCadVsRadCalibValidate.R. 

1.	Highlight lines 1 – 26 and click Run. Echoed lines are not shown below (warnings are ignorable). Enter q <- 1 in console window. Highlight lines 28 – 29 and click Run. Highlight lines 35– 91 and click Run. Highlight simuMu, line 85, and click Run. Highlight simuAlpha, line 86, and click Run. Highlight rhoNorMatr, line 88, and click Run. Highlight rhoAbn2Matr, line 89, and click Run. Highlight rhoNorMatr, line 88, and click Run. Highlight rhoAbn2Matr, line 89, and click Run. This should yield all values shown below under Part A. Highlight lines 96 – 111 and click Run. 
2.	Enter l <- 16 in console window. Highlight lines 113 – 131 and click Run. Highlight conditionArray, line 96, and click Run. Highlight probVector, line 101, and click Run. Highlight muTemp, line 114, and click Run. Highlight sigmaTemp, line 116, and click Run. This should yield all values shown below under Part B. This corresponds to c = 1 in Table 23.1. None of the elements of muTemp are zero, consistent with the fact that in this condition disease was visible to CAD and all radiologists.
3.	Enter l <- 14 and repeat #2. This corresponds to c = 2 in Table 23.1. The first element of muTemp is zero, consistent with the fact that in this condition disease was invisible to CAD.
4.	Enter l <- 13 and repeat #2. This corresponds to c = 3 in Table 23.1. The first two elements of muTemp are zero, consistent with the fact that in this condition disease was invisible to CAD and the first radiologist.
5.	Enter l <- 1 and repeat #2. This corresponds to c = 16 in Table 23.1. Note that all elements of muTemp are zero, consistent with the fact that in this condition disease was invisible to CAD and all radiologists.


23.3.2.1.3.1: Code output (partial)



#### General case

Let   denote column   of the row vector (consisting of zeroes and ones) specified by  . The rule for calculating the elements of the mean vector for visibility-condition   is:

 	.	(23.20)

The rule for calculating the covariance matrix for visibility-condition   is:

  	.	(23.21)

The rule for calculating the probability vector   of the different visibility-conditions is:

 	.	(23.22)

#### Using the simulator

The   non-diseased ratings are generated as follows:

  	.	(23.23)

For diseased cases, one samples the multinomial distribution  times with cell probabilities as specified in . This yields the number of diseased cases in each visibility-condition  . Let   denote the number of diseased cases in visibility-condition  . One generates  samples from the multivariate normal distribution:

  	.	(23.24)

This is repeated for all values of  . This completes the simulation of continuous ratings for a single-modality and (J +1) readers ROC dataset, of which the first reader is CAD. Since the original dataset was binned into 6 bins, the simulated datasets were likewise binned into 6 bins.
23.4: Calibration, validation of simulator and testing its NH behavior 
Table 23.2, Table 23.3 and Table 23.4, summarize the results of the calibration process.

## Calibration, validation of simulator and testing its NH behavior 

Table 23.2, Table 23.3 and Table 23.4, summarize the results of the calibration process.


### Calibration of the simulator 

Table 23.2: This table shows, for the clinical dataset, the calibrated values for the parameter vectors needed in Eqn. (23.15).
	 

Table 23.3: This table shows, for the clinical dataset, the calibrated values for   needed in Eqn. (23.15).
	 
 
 
 
Table 23.3: This table shows, for the clinical dataset, the calibrated values for   needed in Eqn. (23.15).
	 
 

### Validation of simulator and testing its NH behavior 

The code was run with different values of  , as indicated in Table 23.5. The seed variable was set to NULL, line 28, which generates random seeds . Data binning to 6 bins was used, line 26 and line 138 - 140. When the total number of cases is different from that in the clinical dataset, the values of   and   corresponding to simulated dataset  , need to be multiplied by  , line 143. Here NEW refers to the simulated datasets and ORG refers to the original dataset. Confidence intervals for   and   were obtained by bootstrapping readers and cases 2000 times.

Table 23.5: This table summarizes results of validation of the simulation method and results of NH testing. For the original dataset   = 0.00033 (0.000216, 0.000573) an   = 0.00087 (0.000636, 0.00109). In the table   = average of 2000 values of  ,   = average of 2000 values of  , where s is the simulation index  . Instances where the 95% confidence interval for the original dataset did not include the corresponding simulation averaged estimate are indicated with an asterisk. The NH rejection rates were within that expected for 2000 simulations, i.e., they are all in the range (0.04, 0.06).

Results of the evaluation, summarized in Table 23.5 show that, with one exception, as indicated by the asterisk, the estimates of   and   are contained within the 95% confidence interval of the corresponding values for the original data. The estimates of   (average of rows 1 – 7 yields 0.00081) are close to that for the original dataset:   = 0.00087 (0.000636, 0.00109). The estimates of   (corresponding average = 0.00023) are smaller by about 30% than that for the original dataset:   = 0.00033 (0.000216, 0.000573). Row 8 in Table 23.5 is identical to row 7, except that the data was not binned. The variance is smaller, suggesting that binning introduces additional noise into the data, which seems intuitively reasonable to the author. The tendency for the covariance term to be underestimated could be due to over-simplification of the model in our attempt at parsimony. If parameter averaging were not performed, one expects greater variability. 


## Discussion

This chapter describes a method for designing a ratings simulator that is statistically matched ("calibrated") to the single-modality multiple-reader ROC Hupse-Karssemeijer dataset. Showing that it yields, upon analysis with the ORH method, a figure-of-merit variance structure that is consistent with that of the original dataset validates the method. Furthermore, when the NH condition was imposed, the analysis method described in Chapter 22 rejects the NH consistent with the nominal  , thereby validation the analysis method. 

It was noted in previous chapters that the Roe and Metz (RM) simulator is outdated and moreover, there does not exist a systematic way of estimating its parameters. The online appendix to this chapter details the calibration of the Roe and Metz model to the Hupse-Karssemeijer dataset (it is specific to CAD vs. radiologists). When calibrated and the null hypothesis imposed, simulations yield the correct rejection rate, consistent with 5%. The method yielded   = 0.00048, consistent with the original data, but  = 0.00133, which is outside the 95% CI of the original data estimate. 

The method described in this chapter is currently being extended to multiple modalities. It will then be used to test the NH behavior of DBMH and ORH analyses for simulators calibrated to different clinical datasets.

