###------------- Code for POM-PS and POM-IPSW for testing association of multiple secondary phenotypes with a single SNP -------------------
#
# Cite: "A Novel Association Test for Multiple Secondary Phenotypes from a Case-Control GWAS".
#
# POM-PS:   A proportional odds model for the effect of a genetic marker on multiple secondary traits with adjustment
#           for the propensity score (the conditional probability of being a case)
# POM-IPSW: A proportional odds model for the effect of a genetic marker on multiple secondary traits with weights
#           equal to the inverse of propensity scores

##--------------------------------------------- Version 1.15 (dated February 28, 2017) --------------------------------------------------
# Corresponding Authors: Debashree Ray <debashr@umich.edu> ; Saonli Basu <saonli@umn.edu>

############################################
library(MASS)
library(survey)

message("===============================")
message("     POM-PS v1.15 is loaded")
message("===============================")
message("If you use this software, please cite:")
message("Ray et al.(2017) A Novel Association Test for Multiple Secondary Phenotypes")
message("    from a Case-Control GWAS. Genetic Epidemiology, 41. DOI:10.1002/gepi.22045")
message("-------------------------------")


############################################
#-----------------------------------------------Begin: functions to avoid some survey package issues with environment
update.svyglm <- function(object, formula) {
    env <- environment(terms(object))
    if (is.null(env)) {
        env <- parent.frame()
    }
   exp <- update.default(object, formula, evaluate=FALSE)
   eval(exp, env)
}

update.svyolr <- function(object, formula) {
    env <- environment(terms(object))
    if (is.null(env)) {
        env <- parent.frame()
    }
   exp <- update.default(object, formula, evaluate=FALSE)	
   eval(exp, env)
}

model.frame.svyolr <- function(formula, ...) {
    env <- environment(terms(formula))
    if (is.null(env)) {
        env <- parent.frame()
    }
    mcall <- match.call(svyolr, formula$call)
    design <- eval(mcall$design, envir=env)
    formula <- eval(mcall$formula, envir=env)
    mf <- model.frame(formula, model.frame(design))
    w <- weights(design, type = "sampling")
    if (is.null(naa <- attr(mf, "na.action"))) 
        mf[["(weights)"]] <- w
    else mf[["(weights)"]] <- w[-naa]
    mf
}
#-----------------------------------------------End: functions to avoid some survey package issues with environment

#---------------- function for identifying number of parameters for which starting values are required in optim
# this function may come in handy when there is any error in optim and starting values need to be changed/user-specified
.get.start.length<-function(formula, design)
{
    m <- model.frame(formula, model.frame(design), na.action = na.pass)[,1] # first column is response (factor) in this model matrix
    start.length <- length(attr(terms(formula), "term.labels")) + (length(levels(m))-1)
    return(start.length)
}

.get.start.length.polr<-function(formula, X)
{
   start.length <- length(attr(terms(formula), "term.labels")) + (length(unique(X))-1)
   return(start.length)
}


#----------------- function for calculating propensity score P[D|Y,covars]
getPS<-function(Y, D, covars=NULL, no.format.check=FALSE, ...)
{
   if(!no.format.check) {
     # check formats and dimensions of inputs
     if(class(Y)!="data.frame" | class(D)!="data.frame" | (!is.null(covars) & class(covars)!="data.frame") )
       	stop("Inputs Y(phenotype), D(case-control status) and covars(covariates, if provided) must be input in data frame format.")
     n<-nrow(Y)
     if(nrow(D)!=n)
        stop("Sample sizes (no. of rows) in input data frames do not match.")
     if(ncol(D)!=1) stop("Only 1 column for D (case-control status) is allowed.")
     if(length(unique(D[,1]))!=2) stop('D must have only 2 possible values: 0 and 1.')
   }
   if(!is.null(covars)) {
       	if(nrow(covars)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.")
       	traitcov<-cbind(Y,covars)
   }else traitcov<-Y
   PS<-predict(glm(as.matrix(D)~as.matrix(traitcov), family="binomial", maxit=1e5, ...),type="response")
   return(as.data.frame(PS))
}

#---------------- function for calculating propensity-score weights (probabilities)
.PSW<-function(i,D,PS) return( D[i,1]*PS[i,1] + (1-D[i,1])*(1-PS[i,1]) )         # inverse will be taken by svydesign()
.getPSW<-function(D, PS)
{
  D<-as.matrix(D)
  PS<-as.matrix(PS)
  n<-nrow(D)
  if(nrow(PS)!=n) stop('Sample sizes do not match!')
  W<-sapply(1:n, .PSW, D, PS)
  return(as.data.frame(W))
}

#---------------- function for format and dimension checks of inputs in pom.ipw
.format.check<-function(Y, X, D=NULL, COV=NULL, PS=NULL, weights=NULL, method="POM-IPSW", test.method="LRT", msg.mute=FALSE)
{
   # check formats
   if(class(Y)!="data.frame" | class(X)!="data.frame")
       	stop("Inputs Y(phenotype), X(SNP genotype) must be in data frame format.")
   if(!is.null(D) & class(D)!="data.frame") stop("Input D(case-control status, if provided) must be in data frame format.")

   if(is.null(PS)){
       if(method=="POM-IPSW")
          stop('Propensity Scores not provided for method="POM-IPSW".')
   }else {
      if(class(PS)!="data.frame") stop("Input PS(propensity score), if provided, must be in data frame format.")
   }
   if((method=="POM-IPSW" & is.null(PS)) | (method=="POM-IPSW" & is.null(D)))
      stop('Need to provide both propensity scores PS and disease status D for method="POM-IPSW".')
   if(is.null(weights)){
      if(method=="POM-IPW") stop('Input "weights" (inverse probability weights) not provided for method="POM-IPW"')
   }else {
      if(class(weights)=="data.frame") stop('Input "weights" (if provided) must NOT be in data frame format. Input a vector instead.')
   }
   if(!is.null(COV) & class(COV)!="data.frame") stop("Input COV(covariates), if provided, must be in data frame format.")
   if(method!="POM-IPSW" & method!="POM-IPW") stop('For method, choose either "POM-IPSW" or "POM-IPW".')
   if(!is.null(PS) & method=="POM-IPSW" & !isTRUE(msg.mute))
       	message('Caution: POM-IPSW method may have inflated type I error when testing genetic associations of multiple secondary traits. POM-PS is recommended.')
   if(test.method!="Wald" & test.method!="LRT")
       	stop('Testing method not recognized. Input either "Wald" (faster) or "LRT" (robust; default).')

   # check dimensions
   n<-nrow(Y)
   if(nrow(X)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.")
   if(!is.null(D)) { if(nrow(D)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.") }
   if(!is.null(PS)) { if(nrow(PS)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.") }
   if(!is.null(weights)) { if(length(weights)!=n) stop("Length of weight vector (inverse probability weights) do not match no. of rows in other input data frames.") }
   if(!is.null(COV)) { if(nrow(COV)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.") }
   if(ncol(X)!=1) stop("Data frame X (SNP genotype) should have a single column.")
   if(!is.null(PS)) { if(ncol(PS)!=1) stop("Data frame PS (propensity score) should have a single column.") }
}

#---------------- function for format and dimension checks of inputs in pomps
.format.check.pomps<-function(Y, X, D=NULL, COV=NULL, PS=NULL, method="POM-PS", msg.mute=FALSE)
{
   # check formats
   if(class(Y)!="data.frame" | class(X)!="data.frame")
       	stop("Inputs Y(phenotype), X(SNP genotype) must be in data frame format.")
   if(!is.null(D) & class(D)!="data.frame") stop("Input D(case-control status, if provided) must be in data frame format.")

   if(is.null(PS)){
       if(method=="POM-PS")
          stop('Propensity Scores not provided for method="POM-PS".')
   }else {
      if(class(PS)!="data.frame") stop("Input PS(propensity score), if provided, must be in data frame format.")
   }
   if(!is.null(COV) & class(COV)!="data.frame") stop("Input COV(covariates), if provided, must be in data frame format.")

   # check dimensions
   n<-nrow(Y)
   if(nrow(X)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.")
   if(!is.null(D)) { if(nrow(D)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.") }
   if(!is.null(PS)) { if(nrow(PS)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.") }
   if(!is.null(COV)) { if(nrow(COV)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.") }
   if(ncol(X)!=1) stop("Data frame X (SNP genotype) should have a single column.")
   if(!is.null(PS)) { if(ncol(PS)!=1) stop("Data frame PS (propensity score) should have a single column.") }
}


#--------------------------------------- Main function for POM-PS (uses polr) --------------------------------
pomps<-function(Y, X, D=NULL, COV=NULL, PS=NULL, method="POM-PS", add.D.as.COV=FALSE, msg.mute=FALSE, no.format.check=FALSE, ...)
{
   #--------------------------- CHECKS ----------------------------
   if(!no.format.check) .format.check.pomps(Y, X, D, COV, PS, method, msg.mute)

   # check names
   n<-nrow(Y)   # sample size
   q<-0         # no. of covariates
   dataf<-cbind(Y,X)
        Yname<-colnames(Y) ;    Xname<-colnames(X) ;    Dname<-NULL ;    PSname<-NULL ;  COVname<-NULL
   if(!is.null(D)) { dataf<-cbind(dataf, D) ; Dname<-colnames(D) }
   if(!is.null(PS)) { dataf<-cbind(dataf, PS); PSname<-colnames(PS) }
   if(!is.null(COV)) { dataf<-cbind(dataf, COV) ; COVname<-colnames(COV) ; q<-ncol(COV) }
   if( length(unique(colnames(dataf)))!=ncol(dataf) )
        stop("One or more data frame inputs have common column names. Please provide distinct column names across all input data frames.")

   # other checks
   dataf<-na.omit(dataf)
   nobs<-nrow(dataf)
   if(nobs!=n & !isTRUE(msg.mute)) message('Removing samples with missing observations...')
   k<-ncol(Y)
   if(length(unique(dataf[,k+1]))==1)
     stop('All individuals have the same genotype!')
   if(!(is.null(D)))
        if(length(unique(dataf[,k+2]))>2)
           stop('Something wrong with case-control status D. D must take only two possible values: 0 or 1.')
   if(length(unique(dataf[,k+1]))>2) {
     if(!is.factor(dataf[,k+1])) dataf[,k+1]<-as.factor(dataf[,k+1])
   }else {
     if(is.factor(dataf[,k+1])) dataf[,k+1]<-as.numeric(as.character(dataf[,k+1]))
   }

   #------------------------------ FORMULAE -------------------------
   if(method=="POM-PS") {
     # formula for adjusting PS as a covariate
     preds.ps<-paste(c(Yname,PSname,COVname),collapse="+")
     formula.ps<-as.formula(paste(Xname, "~", preds.ps, sep=""))
	# for null model
	preds.ps0<-paste(c(PSname,COVname),collapse="+")
	formula.ps0<-as.formula(paste(Xname, "~", preds.ps0, sep=""))
   }else {
     # formula for no PS adjustment as covariate
     preds.nops<-paste(c(Yname,COVname),collapse="+")
     formula.nops<-as.formula(paste( Xname, "~", preds.nops, sep="" ))
	# for null model
	if(!is.null(COV)) {
	   preds.nops0<-paste(c(COVname),collapse="+")
	   formula.nops0<-as.formula(paste( Xname, "~", preds.nops0, sep="" ))
	}else {
	   formula.nops0<-as.formula(paste( Xname, "~ 1", sep="" ))
	}
   }

   # add D as covariate if asked
   if(isTRUE(add.D.as.COV)) {
        if(is.null(D)) stop('D was not provided; can not add D as covariate. Input D or use add.D.as.COV=FALSE.')
	if(method=="POM-PS"){
          formula.ps<-update(formula.ps, as.formula(paste("~.+",Dname,sep="")))
          formula.ps0<-update(formula.ps0, as.formula(paste("~.+",Dname,sep="")))
	}else {
          formula.nops<-update(formula.nops, as.formula(paste("~.+",Dname,sep="")))
          formula.nops0<-update(formula.nops0, as.formula(paste("~.+",Dname,sep="")))
	}
   }

   #------------------------------MODEL FIT-------------------------
   if(length(unique(dataf[,k+1]))==2)
   {
     if(!isTRUE(msg.mute)) message('Warning: Only 2 possible values for genotype X; fitting logistic model instead of proportional odds model.')
     if(method=="POM-PS" & !is.null(PS)) {
        fit<-try(glm(formula.ps, family="binomial", data=dataf, ...), silent=TRUE)
        fit0<-try(glm(formula.ps0, family="binomial", data=dataf, ...), silent=TRUE)
     }else {
        fit<-try(glm(formula.nops, family="binomial", data=dataf, ...), silent=TRUE)
        fit0<-try(glm(formula.nops0, family="binomial", data=dataf, ...), silent=TRUE)
     }
   }else
   {
     if(method=="POM-PS" & !is.null(PS)) {
        fit<-try(polr(formula.ps, method="logistic", data=dataf, Hess=TRUE, ...), silent=TRUE)
        fit0<-try(polr(formula.ps0, method="logistic", data=dataf, Hess=TRUE, ...), silent=TRUE)
     }else {
        fit<-try(polr(formula.nops, method="logistic", data=dataf, Hess=TRUE, ...), silent=TRUE)
        fit0<-try(polr(formula.nops0, method="logistic", data=dataf, Hess=TRUE, ...), silent=TRUE)
     }
   }

   # Output
   if(inherits(fit, "try-error")){
     # Get length of 'start' (i.e., no. of parameters to estimate)
     if(method=="POM-PS" & !is.null(PS)) {
        start.length<-.get.start.length.polr(formula.ps, dataf[,k+1])
     }else { 
        start.length<-.get.start.length.polr(formula.nops, dataf[,k+1])
     }
     error.msg<-paste(fit[1],"-- For",Xname,"may need to change full model starting values parameter 'start' of length",start.length)
     if(!isTRUE(msg.mute)) message(error.msg)
     betas<-rep(NA, start.length)
     stat<-df<-pval<-NA
   }else if(inherits(fit0, "try-error")){
     # Get length of 'start' (i.e., no. of parameters to estimate)
     if(method=="POM-PS" & !is.null(PS)) {
        start.length<-.get.start.length.polr(formula.ps0, dataf[,k+1])
     }else {
        start.length<-.get.start.length.polr(formula.nops0, dataf[,k+1])
     }
     error.msg<-paste(fit0[1],"-- For",Xname,"may need to change null model starting values parameter 'start' of length",start.length)
     if(!isTRUE(msg.mute)) message(error.msg)
     betas<-summary(fit)$"coefficients"[,1]
     stat<-df<-pval<-NA
   }else {
     betas<-summary(fit)$"coefficients"[,1]

     #------------------------------TESTING-------------------------
     stat<-c(2*(logLik(fit)-logLik(fit0)))
     if("polr" %in% class(fit)){
	df<-fit$"edf"-fit0$"edf"
     }else df<-fit$"rank"-fit0$"rank"
     pval<-pchisq(q=stat, df=df, ncp=0, lower.tail=FALSE)
     error.msg<-"OK"
   }
   return(list(coef=betas, stat=stat, df=df, pvalue=pval, n.obs=nobs, error.msg=error.msg))
}


#-------------------------------------- Main function for POM-IP(S)W (uses svyolr) ---------------------------------
pom.ipw<-function(Y, X, D=NULL, COV=NULL, PS=NULL, weights=NULL, method="POM-IPSW", test.method="LRT", msg.mute=FALSE, no.format.check=FALSE, ...)
{
   #--------------------------- CHECKS ----------------------------
   if(!no.format.check) .format.check(Y, X, D, COV, PS, weights, method, test.method, msg.mute)

   # check names
   n<-nrow(Y)	# sample size
   q<-0         # no. of covariates
   ID<-1:n
   dataf<-cbind(ID,Y,X)
       	Yname<-colnames(Y) ;    Xname<-colnames(X) ;    Dname<-NULL ;    PSname<-NULL ;  COVname<-NULL
   if(!is.null(D)) { dataf<-cbind(dataf, D) ; Dname<-colnames(D) }
   if(!is.null(PS)) { dataf<-cbind(dataf, PS); PSname<-colnames(PS) }
   if(!is.null(COV)) { dataf<-cbind(dataf, COV) ; COVname<-colnames(COV) ; q<-ncol(COV) }
   if(!is.null(weights)) dataf<-cbind(dataf, weights) 
   if( length(unique(colnames(dataf)))!=ncol(dataf) )
       	stop("One or more data frame inputs have common column names. Please provide distinct column names across all input data frames.")

   # other checks
   dataf<-na.omit(dataf)
   nobs<-nrow(dataf)
   if(nobs!=n & !isTRUE(msg.mute)) message('Removing samples with missing observations...')
   k<-ncol(Y)
   if(length(unique(dataf[,k+2]))==1)
     stop('All individuals have the same genotype!')
   if(!(is.null(D))) 
	if(length(unique(dataf[,k+3]))>2)
       	   stop('Something wrong with case-control status D. D must take only two possible values: 0 or 1.')
   if(length(unique(dataf[,k+2]))>2) {
     if(!is.factor(dataf[,k+2])) dataf[,k+2]<-as.factor(dataf[,k+2])
   }else {
     if(is.factor(dataf[,k+2])) dataf[,k+2]<-as.numeric(as.character(dataf[,k+2]))
   }

   #------------------------------DESIGNS & FORMULAE-------------------------
   # Design: Independent sampling design (with replacement) with equal weight for each sample
     des<-svydesign(id=~ID, weights=NULL, data=dataf, probs=NULL, strata=NULL, fpc=NULL)

   # Design for IPSW
     if(method=="POM-IPSW" & !is.null(PS))
     {
	# calculate PSW
       	W<-.getPSW(D=dataf[,k+3], PS=dataf[,k+4])
       	dataf<-cbind(dataf,W)
          dataf<-dataf[which(dataf$W>1e-3),]            ## fix numerically zero values of inverse propensity scores
          dataf<-na.omit(dataf)
          nobs<-nrow(dataf)
       	des.wtd<-svydesign(id=~ID, weights=NULL, data=dataf, probs=~W, strata=NULL, fpc=NULL)
     }

   # Design for IPW
   if(method=="POM-IPW" & !is.null(weights))
     des.wtd2<-svydesign(id=~ID, weights=~weights, data=dataf, probs=NULL, strata=NULL, fpc=NULL)

   # formula for covariates
   preds.nops<-paste(c(Yname,COVname),collapse="+")
   formula.nops<-as.formula(paste( Xname, "~", preds.nops, sep="" ))

   #------------------------------MODEL FIT-------------------------
   if(length(unique(dataf[,k+2]))==2)
   {
     if(!isTRUE(msg.mute)) message('Warning: Only 2 possible values for genotype X; fitting logistic model instead of proportional odds model.')
     if(method=="POM-IPSW" & !is.null(PS)) {
       	fit<-svyglm(formula.nops, design=des.wtd, ...)
     }else if(method=="POM-IPW" & !is.null(weights)) {
       	fit<-svyglm(formula.nops, design=des.wtd2, ...)
     }else {
       	stop("POM-IPSW works only when PS is provided. POM-IPW works only when weights are provided.")
     }
   }else
   {
     if(method=="POM-IPSW" & !is.null(PS)) {
        fit<-svyolr(formula.nops, design=des.wtd, ...)
     }else if(method=="POM-IPW" & !is.null(weights)) {
       	fit<-svyolr(formula.nops, design=des.wtd2, ...)
     }else {
        stop("POM-IPSW works only when PS is provided. POM-IPW works only when weights are provided.")
     }
   }

   # Output
   if(inherits(fit, "try-error")){
     # Get length of 'start' (i.e., no. of parameters to estimate)
     if(method=="POM-IPSW" & !is.null(PS)) {
	start.length<-.get.start.length(formula.nops, des.wtd)
     }else { #if(method=="POM-IPW" & !is.null(weights)) 
	start.length<-.get.start.length(formula.nops, des.wtd2)
     }
     error.msg<-paste(fit[1],"-- For",Xname,"may need to change optim() parameter 'start' of length",start.length)
     if(!isTRUE(msg.mute)) message(error.msg)
     betas<-rep(NA, start.length)
     stat<-df<-pval<-test.method.used<-NA
   }else{
     betas<-coef(fit)
     #------------------------------TESTING-------------------------
     formula.test<-as.formula(paste("~",paste(Yname,collapse="+"),sep=""))
     # Wald is faster but LRT is robust and hence recommended
	# "LRT method will not work if the model had starting values supplied for the regression coefficients."
	# Wald will be used instead.
     if(test.method=="Wald" | (test.method=="LRT" & hasArg(start))) {
	if(test.method=="LRT" & hasArg(start) & !isTRUE(msg.mute))
	   message(paste('Warning: test.method="LRT" does not work if user specifies "start" parameter; Instead "Wald" will be used for',Xname))
        out<-regTermTest(model=fit, test.terms=formula.test, method="Wald", df=Inf)  #df=Inf input produces chi-sq test instead of F-test
        stat<-out$chisq[1]
	test.method.used="Wald"
       	df<-out$df
     	pval<-out$p[1]
     	error.msg<-"OK"
     }else {
        out<-try(regTermTest(model=fit, test.terms=formula.test, method="LRT"), silent=TRUE)
	if(inherits(out, "try-error")) {
	  error.msg<-paste(out[1],"-- Error in saddlepoint approximation of LRT for",Xname)
	  if(!isTRUE(msg.mute)) message(error.msg)
	  stat<-df<-pval<-NA
	}else {
	  stat<-out$chisq/mean(out$lambda)
       	  df<-out$df
     	  pval<-out$p[1]
     	  error.msg<-"OK"
	}
	test.method.used="LRT"	  
     }

   }

   return(list(coef=betas, stat=stat, df=df, pvalue=pval, test.method=test.method.used, n.obs=nobs, error.msg=error.msg))
}


