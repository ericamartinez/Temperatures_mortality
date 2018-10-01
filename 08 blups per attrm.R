################################################################################
# "Impact of ambient temperatures on mortality in Spain (1993-2013)"
#   
#   SEPARATED ANALYSIS FOR PERIOD 1 AND PERIOD 2 +  BLUPS
#   
#   MODEL1 <- risk P1 + temperatures P1
#   MODEL -> risk P2 + temperatures P2
#
#   IMPORTANT: execute this code before attributable measures (06 attrm_new)
#   ISGlobal  
#   May 2017
#   
#
################################################################################



################################################################################
# FIRST-STAGE ANALYSIS: RUN THE MODEL IN EACH PROVINCE, REDUCE AND SAVE
################################################################################

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS

# COEFFICIENTS AND COVARIANCE FOR OVERALL CUMULATIVE SUMMARY
coef1 <- coef2 <- matrix(NA,length(provinces), 25,
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
  #cb <- crossbasis(data$tempmax_compl_censored, lag=lag, argvar=argvar, arglag=arglag)
  
  cb1 <- crossbasis(dat1$tempmax_compl, lag=lag, argvar=argvar1, arglag=arglag)
  # cb1: crossbasis with temperatures P1  
  
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  # NB: NO CENTERING NEEDED HERE, AS THIS DOES NOT AFFECT COEF-VCOV
  # We include day of the week, holidays and hospitalizations due to influenza
 
   # model1: considering risk P1 + temperatures P1data$t=1:dim(data)[1]
  dat1$t <- 1:dim(dat1)[1]
  
  model1 <- glm(adeath ~ cb1 + dow + hday + phday + total_influenza_h +
                  ns (doy,df=dfseas): factor(yyyy) + 
                  ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
                dat1, family=quasipoisson, na.action="na.exclude")
  
  data$t=1:dim(data)[1]
  # model: considering risk P2 + temperatures P2
  model2 <- glm(adeath ~ cb + dow + hday + phday + total_influenza_h +
                 ns (doy,df=dfseas): factor(yyyy) + 
                 ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
               data, family=quasipoisson, na.action="na.exclude")
  
  name <- deparse(substitute(cb))
  cond <- paste0(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
  
  model.class1 <- class(model1)
  coefP1 <- dlnm:::getcoef(model1,model.class1)
  ind1 <- grep(cond,names(coefP1))
  coefP1 <- coefP1[ind1]
  vcovP1 <- dlnm:::getvcov(model1,model.class1)[ind1,ind1,drop=FALSE]
  
  model.class2 <- class(model2)
  coefP2 <- dlnm:::getcoef(model2,model.class2)
  ind2 <- grep(cond,names(coefP2))
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
# (combina els efectes de cada provincia i obtenim una estimaci? global a nivell de pa?s)

# RUN THE MODELS FOR BEFORE AND AFTER PERIODS
mv1 <- mvmeta(coef1~1,vcov1,provincies_n,control=list(showiter=T,maxiter=1000))
#summary(mv1)
mv2 <- mvmeta(coef2~1,vcov2,provincies_n,control=list(showiter=T,maxiter=1000))
#summary(mv2)


################################################################################
# OBTAIN BLUPS (best linear unbiased predictions)
# For random effects models, predictions are the sum of the average results 
# of the fixed part of the model plus the predicted random effects variances
# We can draw the exposure-response curves in each city

blupP1 <- blup(mv1,vcov=T)
blupP2 <- blup(mv2,vcov=T)

save(blupP1,blupP2,file="T:\\Xavi\\ETEC\\Erica\\mortality\\separate models\\blupP.RData")
save(blupP1,blupP2,file="H:\\doctorat\\Mortality\\02_Stata\\01_do\\R code\\Definitiu\\blupP.RData")


