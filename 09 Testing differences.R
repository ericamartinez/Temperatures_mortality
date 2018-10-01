

################################################################################
# "Temporal changes in temperature-related mortality in Spain and effect of 
#        the implementation of a Heat Health Prevention Plan"
#   
#   ISGlobal  
#   June 2018
#   
#
################################################################################

load("H:\\doctorat\\Mortality\\02_Stata\\03_data\\tempDEATHS.Rdata")

# LOAD THE PACKAGES
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel)

# Day of the year (1:365)
tempDEATHS$doy <- as.numeric(strftime(tempDEATHS$date, format = "%j"))

# We exclude 2003
tempDEATHS <- subset(tempDEATHS, tempDEATHS$yyyy!=2003)

tempDEATHS$dow <- as.factor(tempDEATHS$dow) 

###Pick holidays (summer and Christmas)
tempDEATHS$phday <- 0

# First three weeks of August
tempDEATHS[tempDEATHS$dd < 22 & tempDEATHS$mm == 8, 'phday'] <- 1
# From the 23rd of December until the 6th of January
tempDEATHS[(tempDEATHS$dd > 22 & tempDEATHS$mm == 12) | 
             (tempDEATHS$dd < 7 & tempDEATHS$mm == 1), 'phday'] <- 1

# Binary variables for before and after periods
tempDEATHS$int1 <- 0
tempDEATHS[tempDEATHS$yyyy >= 2004, 'int1'] <- 1

tempDEATHS$int2 <- 0
tempDEATHS[tempDEATHS$yyyy < 2004, 'int2'] <- 1

provincies_n <- list("Alava", "Albacete", "Alicante", "Almeria", "Avila", 
                     "Badajoz", "Illes Balears", "Barcelona", "Burgos", "Caceres", "Cadiz", 
                     "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", 
                     "Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen", "Leon", 
                     "Lleida", "La Rioja", "Lugo", "Madrid", "Malaga", "Murcia", "Navarra", 
                     "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
                     "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", 
                     "Tarragona", "Teruel", "Toledo", "Valencia", "Valladolid", "Vizcaya", 
                     "Zamora", "Zaragoza")


# ARRANGE THE DATA AS A LIST OF DATA SETS
provinces <- as.character(unique(tempDEATHS$province)) # My provinces
dlist <- lapply(provinces,function(x) tempDEATHS[tempDEATHS$province==x,]) 
# Create a list with 50 provinces 
#(agafa el data frame i el converteix a llista de tants elements com provincies)
names(dlist) <- provincies_n

# PARAMETERS FOR THE EXPOSURE-RESPONSE FUNCTION
# 3 internal knots placed at the 10th, 75th and 90th percentiles of 
# location-specific temperature distribution
varper <- c(10,75,90)
vardegree <- 2

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3 #Number of knots for lag model
arglag<- list(knots=logknots(lag,lagnk))

# DEEGRES OF FREEDOM
dfseas <- 8  #Seasonality
dftrend <- 1 #Long-term trend



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


########################################################################################################
# SECOND-STAGE ANALYSIS: MULTIVARIATE META-ANALYSIS OF THE REDUCED COEF AND THEN COMPUTATION OF BLUP
########################################################################################################

# CREATE AVERAGE TEMPERATURE AND RANGE AS META-PREDICTORS
avgtmean <- sapply(dlist,function(x) mean(x$tempmax_compl,na.rm=T))
rangetmean <- sapply(dlist,function(x) diff(range(x$tempmax_compl,na.rm=T)))

################################################################################
# META-ANALYSIS
# We pool the estimated location-specific overall cumulative exposure-response associations 


# RUN THE MODELS FOR BEFORE AND AFTER PERIODS
mv1 <- mvmeta(coef1~1,vcov1,provincies_n,control=list(showiter=T))
summary(mv1)
mv2 <- mvmeta(coef2~1,vcov2,provincies_n,control=list(showiter=T))
summary(mv2)


################################################################################
# OBTAIN BLUPS (best linear unbiased predictions)
# For random effects models, predictions are the sum of the average results 
# of the fixed part of the model plus the predicted random effects variances
# We can draw the exposure-response curves in each city

blup1 <- blup(mv1,vcov=T)
blup2 <- blup(mv2,vcov=T)


################################################################################
# RE-CENTERING

# GENERATE THE MATRIX FOR STORING THE RESULTS
minperccity <- mintempcity <- rep(NA,length(dlist))
names(mintempcity) <- names(minperccity) <- provincies_n

minperccity_1 <- minperccity_2 <- minperccity
mintempcity_1 <- mintempcity_2 <- mintempcity

predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
tmeancountry <- rowMeans(sapply(dlist,function(x) quantile(jitter(x$tempmax_compl),
                                                           predper/100,na.rm=T)))

# DEFINE MINIMUM MORTALITY VALUES: EXCLUDE LOW AND VERY HOT TEMPERATURE
for(i in seq(length(dlist))) {
  data <- dlist[[i]]
  predvar <- quantile(data$tempmax_compl,1:99/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  
  # We define 3 knots at the 10th, 75th and 90th percentiles
  argvar <- list(x=predvar,fun="bs",
                 degree=2, knots=quantile(data$tempmax_compl, varper/100, na.rm=T),
                 Bound=range(data$tempmax_compl,na.rm=T))
  
  bvar <- do.call(onebasis,argvar)
  minperccity[i] <- (1:99)[which.min((bvar%*%blup1[[i]]$blup))]
  mintempcity[i] <- quantile(data$tempmax_compl,minperccity[i]/100,na.rm=T)
}

# COUNTRY-SPECIFIC POINTS OF MINIMUM MORTALITY
(minperccountry <- median(minperccity))

#

## BEFORE PERIOD ##
for(i in seq(length(dlist))) {
  data <- dlist[[i]]
  data1 <- subset(data, data$yyyy<=2002)
  
  predvar <- quantile(data1$tempmax_compl,1:99/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  
  # We define 3 knots at the 10th, 75th and 90th percentiles
  argvar <- list(x=predvar,fun="bs",
                 degree=2, knots=quantile(data1$tempmax_compl, varper/100, na.rm=T),
                 Bound=range(data1$tempmax_compl,na.rm=T))
  
  bvar <- do.call(onebasis,argvar)
  minperccity_1[i] <- (1:99)[which.min((bvar%*%blup1[[i]]$blup))]
  mintempcity_1[i] <- quantile(data1$tempmax_compl,minperccity_1[i]/100,na.rm=T)
}

# COUNTRY-SPECIFIC POINTS OF MINIMUM MORTALITY
(minperccountry_1 <- median(minperccity_1))

## AFTER PERIOD ##
for(i in seq(length(dlist))) {
  data <- dlist[[i]]
  data2 <- subset(data, data$yyyy>=2004)
  
  predvar <- quantile(data2$tempmax_compl,1:99/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  
  # We define 3 knots at the 10th, 75th and 90th percentiles
  argvar <- list(x=predvar,fun="bs",
                 degree=2, knots=quantile(data2$tempmax_compl, varper/100, na.rm=T),
                 Bound=range(data2$tempmax_compl,na.rm=T))
  
  bvar <- do.call(onebasis,argvar)
  minperccity_2[i] <- (1:99)[which.min((bvar%*%blup2[[i]]$blup))]
  mintempcity_2[i] <- quantile(data2$tempmax_compl,minperccity_2[i]/100,na.rm=T)
}

# COUNTRY-SPECIFIC POINTS OF MINIMUM MORTALITY
(minperccountry_2 <- median(minperccity_2))

#



################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

################################################################################

prov <- c("Alava", "Albacete", "Alicante", "Almeria", "Avila", "Badajoz", "Illes Balears", "Barcelona", "Burgos",
          "Caceres", "Cadiz", "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", "Granada",
          "Guadalajara", "Guipuzkoa", "Huelva", "Huesca", "Jaen", "Leon", "Lleida", "La Rioja", "Lugo", "Madrid",
          "Malaga", "Murcia", "Navarra", "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
          "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", "Tarragona", "Teruel", "Toledo",
          "Valencia", "Valladolid", "Vizcaia", "Zamora", "Zaragoza")

datanew <- data.frame(avgtmean=mean(tapply(avgtmean,prov,mean)),
                      rangetmean=mean(tapply(rangetmean,prov,mean)))

# PREDICT THE POOLED COEFFICIENTS FOR EACH MODEL
#mvpred <- predict(mv,datanew,vcov=T,format="list")

mvpred1 <- predict(mv1,datanew,vcov=T,format="list")
mvpred2 <- predict(mv2,datanew,vcov=T,format="list")


# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
set.seed(13041975)
tmeancountry <- rowMeans(sapply(dlist,function(x) quantile(jitter(x$tempmax_compl),
                                                           predper/100,na.rm=T)))

# DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATION
bvar <- onebasis(tmeancountry,fun="bs",degree=2, 
                 knots=tmeancountry[paste(varper,".0%",sep="")])

cenindcountry1 <- which.min(bvar%*%mvpred1$fit)
cenindcountry2 <- which.min(bvar%*%mvpred2$fit)

# DEFINE CENTERING PERCENTILE FOR COUNTRY

cenpercountry1 <- pmin(pmax(predper[cenindcountry1],10),90)
cenpercountry2 <- pmin(pmax(predper[cenindcountry2],10),90)

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

# OBTAIN THE CENTERED PREDICTIONS

cen1 <- tmeancountry[paste(cenpercountry1,".0%",sep="")]
cen2 <- tmeancountry[paste(cenpercountry2,".0%",sep="")]

cp1 <- crosspred(bvar,coef=mvpred1$fit,vcov=mvpred1$vcov,model.link="log",
                 at=tmeancountry,cen=cen1)
cp2 <- crosspred(bvar,coef=mvpred2$fit,vcov=mvpred2$vcov,model.link="log",
                 at=tmeancountry,cen=cen2)
#


################################################################################
# EFFECTS BY COUNTRY (TABLE 2 IN THE MANUSCRIPT)
# NB: NOT IDENTICAL TO MANUSCRIPT, AS BASED ON UK ONLY

# FUNCTION FOR MULTIVARIATE WALD TEST (FIRST VERSION, BASED ON COEF-VCOV)
fwald1 <- function(b1,V1,b2=NULL,V2=NULL) {
  invVp <- if(is.null(b2)) solve(V1) else solve(V1+V2)
  b <- if(is.null(b2)) b1 else b1-b2
  stat <- t(b)%*%invVp%*%(b)
  df <- length(b1)
  pchisq(stat,df,lower.tail=FALSE)
}

# MULTIVARIATE TESTS FOR A NULL INTERACTION (P-VALUE)
fwald1(coef(cp1),vcov(cp1),coef(cp2),vcov(cp2))



