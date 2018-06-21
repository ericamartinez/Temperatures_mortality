
################################################################################
# FUNCTION FOR COMPUTING ATTRIBUTABLE MEASURES FROM DLNM
#    Checking the effectiveness of the plan
#    Have the AF decreased considering the risk in P2 and temperatures in P2?
#    We compare the number of deaths attributed to cold and heat considering:
#         Risk P1 + Temperatures P2
#         Risk P2 + Temperatures P2 
################################################################################

################################################################################
#   - data: DATA FOR PROVINCE i
#   - coef1 AND vcov1: COEF AND VCOV FOR BLUP1 (BEFORE PERIOD) PROVINCE i
#   - coef2 AND vcov2: COEF AND VCOV FOR BLUP2 (AFTER PERIOD) PROVINCE i
#   - type: EITHER "an" OR "af" FOR ATTRIBUTABLE NUMBER OR FRACTION
#   - dir: EITHER "back" OR "forw" FOR BACKWARD OR FORWARD PERSPECTIVES
#   - tot: IF TRUE, THE TOTAL ATTRIBUTABLE RISK IS COMPUTED
#   - cen: THE REFERENCE VALUE USED AS COUNTERFACTUAL SCENARIO
#   - range: THE RANGE OF EXPOSURE. IF NULL, THE WHOLE RANGE IS USED
#   - sim: IF SIMULATION SAMPLES SHOULD BE RETURNED. ONLY FOR tot=TRUE
#   - nsim: NUMBER OF SIMULATION SAMPLES
################################################################################


attrdl_effectiveness <- function (data,coef1,vcov1,coef2,vcov2,type="an",dir="back",tot=TRUE,
                                  cen,range=NULL,sim=FALSE,nsim=5000){
  
  x <- data[data$bef_aft==1,]$tempmax_compl #temperatures P2
  cases <- data[data$bef_aft==1,]$adeath #deaths P2
  
  dat1=data[data$bef_aft==0,] #Before period (P1)
  dat2=data[data$bef_aft==1,] #After period (P2)
  
  data=dat2
  
  # DEFINE THE CROSSBASIS
  argvar <- list(fun="bs",degree=2,knots=quantile(data$tempmax_compl,varper/100,na.rm=T))
  argvar1 <- list(fun="bs",degree=2,knots=quantile(dat1$tempmax_compl,varper/100,na.rm=T))
  
  # We define 3 internal knots placed at 10,75,90th percentiles
  cb <- crossbasis(data$tempmax_compl, lag=lag, argvar=argvar, arglag=arglag)
  cb1 <- crossbasis(dat1$tempmax_compl, lag=lag, argvar=argvar1, arglag=arglag)
  
  
  name <- deparse(substitute(cb))
  
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) 
    x[x<range[1]| x>range[2]] <- cen
  
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"
  lag <- attr(cb,"lag")
  if(NCOL(x)==1L) {
    at <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2])) else 
      matrix(rep(x,diff(lag)+1),length(x))
    # if cen is a vector, uncomment the following line
    # at2 <- if(dir=="back") tsModel:::Lag(cen,seq(lag[1],lag[2])) else 
    #   matrix(rep(cen,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(at <- x)!=diff(lag)+1) 
      stop("dimension of 'x' not compatible with 'basis'")
  }
  
  # cases: TRANFORM IN MEAN OF FUTURE CASES IF dir="forw"
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    cases <- rowMeans(cases)
  } else {
    if(dir=="forw") cases <- rowMeans(tsModel:::Lag(cases,-seq(lag[1],lag[2])))
  }
  
  # IF REDUCED ESTIMATES ARE PROVIDED
  typebasis <- ifelse(length(coef2)!=ncol(cb),"one","cb")
  
  ################################################################################
  #
  # PREPARE THE ARGUMENTS FOR THE BASIS TRANSFORMATION
  predvar <- if(typebasis=="one") x else seq(nrow(at))
  predlag <- if(typebasis=="one") 0 else dlnm:::seqlag(lag)
  #  
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON typebasis)
  lagvec <- rep(predlag, each = length(predvar))
  
  basisvar.no.centered <- do.call("onebasis", c(list(x = at), attr(cb, "argvar")))
  basisvar.no.centered1 <- do.call("onebasis", c(list(x = at), attr(cb1, "argvar")))
  
  basislag <- do.call("onebasis", c(list(x = lagvec), attr(cb, "arglag")))
  
  basiscen <- do.call("onebasis", c(list(x = cen), attr(cb, "argvar")))
  basiscen1 <- do.call("onebasis", c(list(x = cen), attr(cb1, "argvar")))
  
  ## If we want to center using another temperature vector:
  # Use the following basiscen if you want to center to a vector
  #basiscen <- do.call("onebasis", c(list(x = at2), attr(cb, "argvar")))
  
  basisvar <- scale(basisvar.no.centered, center = basiscen, scale = FALSE)
  basisvar1 <- scale(basisvar.no.centered1, center = basiscen1, scale = FALSE)
  ## note that basisvar is the difference between the real values and the value temp=cen (everyday)
  
  Xpred <- tensor.prod.model.matrix(list(basisvar, basislag))
  Xpred1 <- tensor.prod.model.matrix(list(basisvar1, basislag))
  
  Xpredall <- 0
  Xpredall1 <- 0
  for (j in seq(length(predlag))) {
    ind <- seq(length(predvar))+length(predvar)*(j-1)
    Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
    Xpredall1 <- Xpredall1 + Xpred1[ind,,drop=FALSE]
  }
  
  ################################################################################
  #
  # COMPUTE AF AND AN 
  # Risk P2 + temperatures P2
  af2 <- 1-exp(-rowSums(as.matrix(Xpredall%*%coef2)))
  an2 <- af2*cases
  
  # Risk P1 + temperatures P2
  af1 <- 1-exp(-rowSums(as.matrix(Xpredall1%*%coef1)))
  an1 <- af1*cases
  
  #
  # TOTAL (ACCOUNTING FOR MISSING)
  isna <- is.na(an2)
  af2 <- sum(an2[!isna])/sum(cases[!isna])
  an2 <- af2*sum(cases,na.rm=T)
  
  isna1 <- is.na(an1)
  af1 <- sum(an1[!isna1])/sum(cases[!isna1])
  an1 <- af1*sum(cases,na.rm=T)
  
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  
  # NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
  if (sim){
    # SAMPLE COEF
    k <- length(coef2)
    k1 <- length(coef1)
    eigen <- eigen(vcov2)
    eigen1 <- eigen(vcov1)
    X <- matrix(rnorm(length(coef2)*nsim),nsim)
    X1 <- matrix(rnorm(length(coef1)*nsim),nsim)
    coefsim <- coef2 + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    coefsim1 <- coef1 + eigen1$vectors %*% diag(sqrt(eigen1$values),k1) %*% t(X1)
    
    # RUN THE LOOP
    afsim <- apply(coefsim,2, function(coefi) {
      ani <- (1-exp(-drop(Xpredall%*%coefi)))*cases
      sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
    })
    ansim <- afsim*sum(cases,na.rm=T)
    
    afsim1 <- apply(coefsim1,2, function(coefi) {
      ani1 <- (1-exp(-drop(Xpredall1%*%coefi)))*cases
      sum(ani1[!is.na(ani1)])/sum(cases[!is.na(ani1)])
    })
    ansim1 <- afsim1*sum(cases,na.rm=T)
  }
  #
  ################################################################################
  #
  res2 <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an2 else af2    
  }
  #
  
  res1 <- if(sim) {
    if(type=="an") ansim1 else afsim1
  } else {
    if(type=="an") an1 else af1    
  }
  
  res <- c(res1,res2)
  return (res)

}



