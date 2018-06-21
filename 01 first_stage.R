
################################################################################
# "Temporal changes in temperature-related mortality in Spain and effect of 
#        the implementation of a Heat Health Prevention Plan"
#   
#   ISGlobal  
#   June 2018
#   
#
################################################################################


################################################################################
# FIRST-STAGE ANALYSIS: RUN THE MODEL IN EACH PROVINCE, REDUCE AND SAVE
################################################################################

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS

# COEFFICIENTS AND COVARIANCE FOR OVERALL CUMULATIVE SUMMARY
coef <- coef1 <- coef2 <- coefint <- matrix(NA,length(provinces), length(varper)+vardegree,
               dimnames= list(provincies_n))

# The matrix has 5 columns

vcov <- vcov1 <- vcov2 <- vcovint <- vector("list",length(provinces))
names(vcov) <- names(vcov1) <- names(vcov2) <- names(vcovint) <- provincies_n



################################################################################
# RUN THE LOOP

# LOOP
time <- proc.time()[3]
for(i in seq(length(dlist))) {
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]

  # BEFORE PERIOD (1993-2002)
  data1 <- subset(data, data$yyyy<=2002)

  data1$t=1:dim(data1)[1]

  argvar1 <- list(fun="bs",degree=2,knots=quantile(data1$tempmax_compl,varper/100,na.rm=T))
  
  cb1 <- crossbasis(data1$tempmax_compl, lag=lag, 
                    argvar=argvar1,arglag=arglag)
  
  model1 <- glm(adeath ~ cb1 + dow + hday + phday + total_influenza_h +
                  ns (doy,df=dfseas): factor(yyyy) + 
                  ns(t, df=round(length(unique(data1$yyyy))/dftrend/10)), 
                data1, family=quasipoisson, na.action="na.exclude")
  
  # AFTER PERIOD (2004-2013)
  data2 <- subset(data, data$yyyy>=2004)
  data2$t=1:dim(data2)[1]

  argvar2 <- list(fun="bs",degree=2,knots=quantile(data2$tempmax_compl,varper/100,na.rm=T))

  cb2 <- crossbasis(data2$tempmax_compl, lag=lag, 
                    argvar=argvar2, arglag=arglag)

  model2 <- glm(adeath ~ cb2 + dow + hday + phday + total_influenza_h +
                  ns (doy,df=dfseas): factor(yyyy) + 
                  ns(t, df=round(length(unique(data2$yyyy))/dftrend/10)), 
                data2, family=quasipoisson, na.action="na.exclude")
  
  # Predicted effects: extract from the model those parameters corresponding to 
  # cross-basis variables through functions coef and vcov
  pred1 <- crosspred(cb1,model1)
  pred2 <- crosspred(cb2,model2)
  
     
  # REDUCTION TO OVERALL CUMULATIVE
  # Sum the effects of all lags in order to eliminate one dimension of the association. 
  #Sum (acumulate) the risk during the lag period

  red1 <- crossreduce(cb1, model1)
  red2 <- crossreduce(cb2, model2)
  
  # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL

  coef1[i,] <- coef(red1)
  vcov1[[i]] <- vcov(red1)
  
  coef2[i,] <- coef(red2)
  vcov2[[i]] <- vcov(red2)

 
}
proc.time()[3]-time

#


