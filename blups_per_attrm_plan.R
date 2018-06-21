################################################################################
################################################################################
# "Temporal changes in temperature-related mortality in Spain and effect of 
#        the implementation of a Heat Health Prevention Plan"
#   
#   ISGlobal  
#   June 2018
#   
#
################################################################################

blups_plan <- function(provinces,provincies_n,dlist,varper,lag,arglag,name_var,
                  dfseas,dftrend){
  
  
  ################################################################################
  # CREATE THE OBJECTS TO STORE THE RESULTS
  
  # COEFFICIENTS AND COVARIANCE FOR OVERALL CUMULATIVE SUMMARY
  coef1 <- coef2 <- matrix(NA,length(provinces), 17,
                           dimnames= list(provincies_n))
  # The matrix has 15 columns
  vcov1 <- vcov2 <- vector("list",length(provinces))
  names(vcov1) <- names(vcov2) <- provincies_n
  
 
  ################################################################################
  # RUN THE LOOP
  
  # LOOP
  time <- proc.time()[3]
  for(i in seq(length(dlist))) {
      
    # PRINT
    cat(i,"")
    
    # EXTRACT THE DATA
    data <- dlist[[i]]
    dat1=data[data$bef_aft==0,] #Before period (P1)
    dat2=data[data$bef_aft==1,] #After period (P2)
    
    data=dat2
    
    # DEFINE THE CROSSBASIS
    
    argvar <- list(fun="bs",degree=2,knots=quantile(data$tempmax_compl,varper/100,na.rm=T))
    argvar1 <- list(fun="bs",degree=2,knots=quantile(dat1$tempmax_compl,varper/100,na.rm=T))
    
    # We define 3 internal knots placed at 10,75,90th percentiles
    cb <- crossbasis(data$tempmax_compl, lag=lag, argvar=argvar, arglag=arglag)
    cb1 <- crossbasis(dat1$tempmax_compl, lag=lag, argvar=argvar1, arglag=arglag)
    # cb1: crossbasis with temperatures P1  
    
    
    # RUN THE MODEL AND OBTAIN PREDICTIONS
    # NB: NO CENTERING NEEDED HERE, AS THIS DOES NOT AFFECT COEF-VCOV
    # We include day of the week, holidays and hospitalizations due to influenza
    
    # model1: considering risk P1 + temperatures P1data$t=1:dim(data)[1]
    dat1$t <- 1:dim(dat1)[1]
    
    outcome1 <- dat1[[name_var]]
    
    model1 <- glm(outcome1 ~ hw_plan + cb1 + dow + hday + phday + total_influenza_h +
                    ns (doy,df=dfseas): factor(yyyy) + 
                    ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
                  dat1, family=quasipoisson, na.action="na.exclude")
    
    data$t=1:dim(data)[1]
    outcome <- data[[name_var]]
    
    # model: considering risk P2 + temperatures P2
    model2 <- glm(outcome ~ hw_plan + cb + dow + hday + phday + total_influenza_h +
                    ns (doy,df=dfseas): factor(yyyy) + 
                    ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
                  data, family=quasipoisson, na.action="na.exclude")
    
    # name <- deparse(substitute(cb))
    # cond <- paste0(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
    
    cond <- "cb[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}"
    
    
    model.class1 <- class(model1)
    coefP1 <- dlnm:::getcoef(model1,model.class1)
    ind1 <- grep(cond,names(coefP1))
    ind1 <- c(2,ind1)
    # Coefficients cb + hw_plan
    coefP1 <- coefP1[ind1]
    vcovP1 <- dlnm:::getvcov(model1,model.class1)[ind1,ind1,drop=FALSE]
    
    
    model.class2 <- class(model2)
    coefP2 <- dlnm:::getcoef(model2,model.class2)
    ind2 <- grep(cond,names(coefP2))
    ind2 <- c(2,ind2)
    # Coefficients cb + hw_plan
    coefP2 <- coefP2[ind2]
    vcovP2 <- dlnm:::getvcov(model2,model.class2)[ind2,ind2,drop=FALSE]
     
    
    # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
    coef1[i,] <- coefP1
    vcov1[[i]] <- vcovP1
    
    coef2[i,] <- coefP2
    vcov2[[i]] <- vcovP2
    
    
  }
  proc.time()[3]-time
  
 
  ################################################################################
  # META-ANALYSIS
  # We pool the estimated location-specific overall cumulative exposure-response associations 
    
  # RUN THE MODELS FOR BEFORE AND AFTER PERIODS
  mv1 <- mvmeta(coef1~1,vcov1,provincies_n,control=list(showiter=T,maxiter=10000))
  #summary(mv1)
  mv2 <- mvmeta(coef2~1,vcov2,provincies_n,control=list(showiter=T,maxiter=10000))
  #summary(mv2)
  
   
  ################################################################################
  # OBTAIN BLUPS (best linear unbiased predictions)
  # For random effects models, predictions are the sum of the average results 
  # of the fixed part of the model plus the predicted random effects variances
  # We can draw the exposure-response curves in each city
  
  blupP1 <- blup(mv1,vcov=T)
  blupP2 <- blup(mv2,vcov=T)
  
  # We separate coefficents of hw_plan and cb (before and after periods)
  coef_blupP1_hw_plan <- lapply(blupP1, function(x) x$blup[[1]])
  coef_blupP1_cb <- lapply(blupP1, function(x) x$blup[2:17])
  
  vcov_blupP1_hw_plan <- lapply(blupP1, function(x) x$vcov[1,1])
  vcov_blupP1_cb <- lapply(blupP1, function(x) x$vcov[c(2:17),c(2:17)])
  
  coef_blupP2_hw_plan <- lapply(blupP2, function(x) x$blup[[1]])
  coef_blupP2_cb <- lapply(blupP2, function(x) x$blup[2:17])
  
  vcov_blupP2_hw_plan <- lapply(blupP2, function(x) x$vcov[1,1])
  vcov_blupP2_cb <- lapply(blupP2, function(x) x$vcov[c(2:17),c(2:17)])
  
  
  return(list(coef_blupP1_hw_plan,coef_blupP1_cb,coef_blupP2_hw_plan, coef_blupP2_cb,
              vcov_blupP1_hw_plan,vcov_blupP1_cb,vcov_blupP2_hw_plan,vcov_blupP2_cb,
              provinces,provincies_n,dlist))
  
}

