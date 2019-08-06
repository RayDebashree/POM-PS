### Description
POM-PS (proportional odds model adjusted for propensity score) is a novel statistical approach for testing association of one or more secondary traits with a single genetic marker within the framework of a case-control design. The R function `pomps` implements this association test. For details of this statistical method, please refer/cite:

Ray, D. and Basu, S. "[A Novel Association Test for Multiple Secondary Phenotypes from a
Case-Control GWAS](http://onlinelibrary.wiley.com/wol1/doi/10.1002/gepi.22045/full)". *Genetic Epidemiology*, 41:413-426, 2017. [PMID: 28393390](https://www.ncbi.nlm.nih.gov/pubmed/28393390)

**Key Words:** Case-control design; Cross-phenotype association; GWAS; Joint modeling; Multiple traits;
Multivariate analysis; Propensity Score; Proportional Odds Model; Secondary traits; Stratification Score


### Requirements
R (>= 2.14.0), MASS, survey


### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/POM-PS/blob/master/pomps_v1.15.R?raw=TRUE")
```


### Changes
Version 1.15 - February 28, 2017
> First public release of the software.


### Usage

#### Simple example
```{r}
pomps(Y, X, D=NULL, COV=NULL, PS=NULL, method="POM-PS", add.D.as.COV=FALSE, msg.mute=FALSE, no.format.check=FALSE, ...)
```
#### Arguments
| Input | Description |
| ---: | --- |
| `Y` | The `nxK` phenotype matrix, where `n` is the number of individuals and `K` is the number of secondary phenotypes. The joint association of all `K` phenotypes with the single marker will be tested. `Y` needs to be in R data frame format. |
| `X` | The `nx1` column matrix for the single genetic marker, where `n` is the number of individuals. `X` needs to be in R data frame format. |
| `D` | The `nx1` column matrix for the case-control status (primary phenotype), where `n` is the number of individuals. `D` needs to be in R data frame format. |
| `COV` | The `nxq` matrix of covariates that need to be adjusted in the model. `q` is the number of such covariates. `COV` needs to be in R data frame format. The default value is NULL, i.e., it is assumed there is no covariate in the model. The propensity score, although adjusted as a covariate, is input separately and not included in `COV`. |
| `PS` | The `nx1` column matrix for the propensity score (estimated probability of being a case), where `n` is the number of individuals. It can be calculated using the `getPS()` function included in this software. `PS` needs to be in R data frame format. |
| `method` | The method to be used for testing genetic association of traits. Default value is `POM-PS`. The other possible value is `POM`, which tests the genetic association of traits in a proportional odds framework without any adjustment of the propensity score. `POM` may be used when the traits are from a random sample. |
| `add.D.as.COV` | Default value is FALSE. If TRUE, the case-control status `D` is adjusted as a covariate. |
| `msg.mute` | Default value is FALSE, which allows messages to print when the analysis is in progress. |
| `no.format.check` | Default value is FALSE, which allows the code to check if the input parameters conform to the required format. |
| `...` | Additional arguments to be passed to `polr()` (used when proportional odds model framework is used) or `glm()` (used when logistic model framework is used). |

#### Value
| Output | Description |
| ---: | --- |
| `coef` | The regression coefficients in the model. |
| `stat` | The likelihood ratio test (LRT) statistic for testing the null hypothesis of no genetic association. |
| `df` | The degrees of freedom of the asymptotic null distribution. |
| `pvalue` | The p-value from testing the null hypothesis. |
| `n.obs` | Number of individuals (with complete observations) used for testing association. Individuals with missing observations in `Y`, `X,` `D`, `PS` or `COV` are removed. |
| `error.msg` | Saves any error message that arises during the analysis. If no error is encountered, message "OK" is returned. |


### A Working Example
```
source("pomps_v1.15.R")
# simulate 2 phenotypes on 1000 individuals
Y<-mvrnorm(n=1000, mu=c(0,0), Sigma=matrix(c(1,0.2,0.2,1),2,2))

# simulate 3 covariates for disease status
Sig<-matrix(c(1,0.1,-0.45,0.1,1,0.6,-0.45,0.6,1),3,3)
Z<-mvrnorm(n=1000, mu=c(0,0,0), Sigma=Sig)

# simulate error for disease status
E<-rt(n=1000, df=1,ncp=0)

# simulate disease status
D<-rep(0,1000)
D[which(rowSums(cbind(Y,Z,E))>=0)]<-1

# simulate a single marker for 1000 individuals
X<-matrix(rbinom(n=1000, size=2, prob=0.2), ncol=1)	# additive model

# required data-frame formats
Y<-as.data.frame(Y)
D<-as.data.frame(D)
Z<-as.data.frame(Z)
X<-as.data.frame(X)

# unique column names for the data-frames
colnames(Y)<-paste("Y",1:2,sep="")
colnames(Z)<-paste("Z",1:3,sep="")
colnames(X)<-"X"
colnames(D)<-"D"

## apply POM to test association
out1<-pomps(Y=Y, X=X, COV=NULL, method="POM")
# POM test statistic and p-value
t1<-out1$stat
p1<-out1$pvalue

## apply POM-PS to test association
S=getPS(Y=Y, D=D, covars=Z)
out2<-pomps(Y=Y, X=X, COV=NULL, PS=S)
# POM-PS test statistic and p-value
t2<-out2$stat
p2<-out2$pvalue

```

### Notes
1. The proportional odds framework of POM-PS allows one or more secondary phenotypes, which may be binary and/or continuous.

2. The method POM-PS and its software is designed for unrelated individuals. If you have two cohorts with overlapping samples and you want to analyze the combined sample, it is desirable to exclude the overlapping individuals, and any related individuals. 

3. The genotype X takes values 0, 1 or 2. For a given sample, if X has only two possible values, a logistic regression (`glm`) is used instead of proportional odds regression (`polr`).

4. POM-PS is a single-variant association test, and is not expected to work well for rare variants (i.e., when the variant allele-frequency is very low).

5. Apart from function `pomps()`, this software has other relevant functions such as `getPS()` (for calculating propensity score of an individual given secondary phenotypes, case-control status and other covariates - excluding genetic variants) and `pom.ipw()` (for implementing an inverse-probability-weighted proportional odds regression model). Details of these functions are provided in the attached manual.
    * Remember: The estimated propensity score using the `getPS()` function needs to be obtained only once for a given GWAS.
    * Caution: While `pomps()` removes missing observations (if any) from the input data, `getPS()` is not programmed to handle missing observations (NA values) and hence we need to provide complete (non-missing) data as input for `getPS()`.
    * Caution: `pom.ipw()` often encounters convergence issues; use with caution.

6. POM-PS currently analyzes secondary phenotypes when the primary trait is binary. For primary trait with three or more ordinal categories, one may employ a proportional odds model for calculating the propensity score.
