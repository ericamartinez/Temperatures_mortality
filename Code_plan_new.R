
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



################################################################################
# PREPARE THE DATA
################################################################################

# LOAD THE PACKAGES
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel);library(mgcv)
library(foreach); library(doSNOW); 
library(Epi); library(metafor)


load("tempDEATHS.Rdata")


############### ACTIVATION PLAN INDICATOR ######################################
# Create indicator activation of the Plan according to the thresholds established
# for each province

tempDEATHS$prov <- as.integer(tempDEATHS$province)

library(zoo)
tempDEATHS$hw_indicator <- 0

prov_thresholds <- data.frame(
  prov=c(seq(1,50,by=1)),
  threshold_max=c(34,36,35,35,33,40,35,30.5,33,38,33,33,35,41,33,32,34,39,35,36,37,34,39,
                  33,37,36,33,36.5,36,38,36,37,33,36,33,33,35,33,35,34,41,34,33,35,38,34,
                  36,37,35,38),
  threshold_min=c(20,20,23,24,22,21,22,22,20,23,24,23,22,22,20,21,20,23,21,22,22,20,25,20,
                  21,22,20,21,23,22,22,21,20,21,23,22,20,23,22,20,22,20,22,20,22,23,21,21,
                  22,21)
)

for(i in 1:nrow(prov_thresholds)){
  p <- tempDEATHS$prov == i
  hw_threshold <- tempDEATHS[p,]$tempmax_compl > prov_thresholds[i, 'threshold_max'] & tempDEATHS[p,]$tempmin_compl > prov_thresholds[i, 'threshold_min']
  hw_window <- rollapply(hw_threshold, width=5, sum, partial=TRUE)
  hw_indicator <- rollapply(hw_window, width=5, max, partial=TRUE)
  tempDEATHS[p,]$hw_indicator <- hw_indicator  
}




# Create heat wave indicator according to the Plan (0: no activation / 1: activation)
tempDEATHS$hw_plan <- 0
tempDEATHS[tempDEATHS$hw_indicator > 0, 'hw_plan'] <- 1


# Day of the year (1:365)
tempDEATHS$doy <- as.numeric(strftime(tempDEATHS$date, format = "%j"))
tempDEATHS$dow <- as.factor(tempDEATHS$dow) 

# We exclude 2003
tempDEATHS <- subset(tempDEATHS, tempDEATHS$yyyy!=2003)

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

provincies_n_total <- list("Alava", "Albacete", "Alicante", "Almeria", "Avila", 
                     "Badajoz", "Illes Balears", "Barcelona", "Burgos", "Caceres", "Cadiz", 
                     "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", 
                     "Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen", "Leon", 
                     "Lleida", "La Rioja", "Lugo", "Madrid", "Malaga", "Murcia", "Navarra", 
                     "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
                     "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", 
                     "Tarragona", "Teruel", "Toledo", "Valencia", "Valladolid", "Vizcaya", 
                     "Zamora", "Zaragoza")


# Select months when the Plan is activated
tempDEATHS <- subset(tempDEATHS, (tempDEATHS$mm>=6 & tempDEATHS$mm<=8) | (tempDEATHS$mm==9 & tempDEATHS$dd<=15))


# ARRANGE THE DATA AS A LIST OF DATA SETS
provinces_total <- as.character(unique(tempDEATHS$province)) # My provinces
dlist_total <- lapply(provinces_total,function(x) tempDEATHS[tempDEATHS$province==x,]) 
# Create a list with 50 provinces 
names(dlist_total) <- provincies_n_total

# PARAMETERS FOR THE EXPOSURE-RESPONSE FUNCTION
# 2 internal knots placed at 75th and 90th percentiles of 
# location-specific temperature distribution
varper <- c(75,90)
vardegree <- 2

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 2 #Number of knots for lag model
arglag<- list(knots=logknots(lag,lagnk))

# DEEGRES OF FREEDOM
dfseas <- 4  #Seasonality
dftrend <- 1 #Long-term trend


# CREATE THE MATRICES TO STORE THE RESULTS
main.eff <- added.eff <- matrix(NA,length(provincies_n_total),2, 
                                dimnames=list(provincies_n_total,paste("hw",c("est","sd"),sep=".")))
main.eff_before <- added.eff_before <- main.eff_after <- added.eff_after <- main.eff

# COEFFICIENTS AND COVARIANCE FOR OVERALL CUMULATIVE SUMMARY
coef <- coef1 <- coef2 <- matrix(NA,length(provinces_total), length(varper)+vardegree,
                                            dimnames= list(provincies_n_total))

vcov <- vcov1 <- vcov2 <- vector("list",length(provinces_total))
names(vcov1) <- names(vcov2) <- provincies_n_total


for(i in seq(length(dlist_total))) {
  
  # FIRST STAGE
  
  data <- dlist_total[[i]]
  percentiles <- quantile(data$tempmax_compl,c(75,90,92.5,95,97.5)/100,na.rm=T)
  range <- round(range(data$tempmax_compl,na.rm=T),0)
   
  
  # RUN THE MODELS
  # BEFORE
  data1 <- subset(data, data$yyyy < 2003)
  if (sum(data1$hw_plan)>0){
    # DEFINE THE CROSSBASIS
    argvar <- list(fun="bs",degree=2,knots=quantile(data1$tempmax_compl,varper/100,na.rm=T))
    
    # We define 2 internal knots placed at 75,90th percentiles
    cb1 <- crossbasis(data1$tempmax_compl, lag=lag, argvar=argvar, 
                     arglag=arglag)
  
    data1$t=1:dim(data1)[1]
  
    model1 <- glm(adeath ~ hw_plan + cb1 + dow + hday + phday + total_influenza_h +
                  ns (doy,df=dfseas): factor(yyyy) + 
                  ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
                data1, family=quasipoisson, na.action="na.exclude")
    
    red1 <- crossreduce(cb1,model1)
    coef1[i,] <- coef(red1)
    vcov1[[i]] <- vcov(red1)
  
    # AFTER
    data2 <- subset(data, data$yyyy > 2003)
  
    # DEFINE THE CROSSBASIS
    argvar <- list(fun="bs",degree=2,knots=quantile(data2$tempmax_compl,varper/100,na.rm=T))
  
    # We define 2 internal knots placed at 75,90th percentiles
    cb2 <- crossbasis(data2$tempmax_compl, lag=lag, argvar=argvar, 
                     arglag=arglag)
  
    data2$t=1:dim(data2)[1]
  
    model2 <- glm(adeath ~ hw_plan + cb2 + dow + hday + phday + total_influenza_h +
                  ns (doy,df=dfseas): factor(yyyy) + 
                  ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
                data2, family=quasipoisson, na.action="na.exclude")
    
    red2 <- crossreduce(cb2,model2)
    coef2[i,] <- coef(red2)
    vcov2[[i]] <- vcov(red2)
  }
  
  # SAVE MAIN EFFECT - BEFORE PERIOD
  if(sum(data1$hw_plan)>0) {
    tmedian <- median(data$tempmax_compl[data$hw_plan==1],na.rm=T)
    
    pred.before <- crosspred(cb1,model1,at=c((range[1]+1):(range[2]-1),tmedian),
                      cen=percentiles[1])
    pred.after <- crosspred(cb2,model2,at=c((range[1]+1):(range[2]-1),tmedian),
                                    cen=percentiles[1])
    main.eff_before[i,c(1,2)] <- cbind(pred.before$allfit,
                                      pred.before$allse)[as.character(tmedian),]
    main.eff_after[i,c(1,2)] <- cbind(pred.after$allfit,
                                              pred.after$allse)[as.character(tmedian),]
    # SAVE ADDED EFFECT
    added.eff_before[i,c(1,2)] <- ci.lin(model1)["hw_plan",1:2]
    added.eff_after[i,c(1,2)] <- ci.lin(model2)["hw_plan",1:2]
  } else {main.eff_before[i,c(1,2)] <- c(NA,NA)
          main.eff_after[i,c(1,2)] <- c(NA,NA)}
  
}
  


########################################################################################################
# SECOND-STAGE ANALYSIS: UNIVARIATE META-ANALYSIS
########################################################################################################


# META-ANALYSIS
# We pool the estimated location-specific overall cumulative exposure-response associations 

main.eff_before <- main.eff_before[-c(3),]
main.eff_after <- main.eff_after[-c(3),]

# RUN THE MODELS FOR BEFORE AND AFTER PERIODS
pool.main.before <- rma.uni(yi=main.eff_before[,1], sei=main.eff_before[,2])
summary(pool.main.before)
pool.main.after <- rma.uni(yi=main.eff_after[,1], sei=main.eff_after[,2])
summary(pool.main.after)

pool.added.before <- rma.uni(yi=added.eff_before[,1], sei=added.eff_before[,2])
summary(pool.added.before)
pool.added.after <- rma.uni(yi=added.eff_after[,1], sei=added.eff_after[,2])
summary(pool.added.after)


##################################################################

label <- paste("hw_plan")
table1.before <- table1.after <- matrix(NA,1,7,dimnames=list(label,
                                      c("N comm","Est.main","95%CI.main","P-het.added","Est.added",
                                        "95%CI.added","P-het.added")))

# FILL TABLE1
table1.before[1,] <- c(sum(!is.na(added.eff_before[,1])),
                  round(exp(pool.main.before$b)*100-100,1),
                  paste(round(exp(pool.main.before$b-1.96*pool.main.before$se)*100-100,1),"to",
                        round(exp(pool.main.before$b+1.96*pool.main.before$se)*100-100,1)),
                  round(pool.main.before$QEp,3),
                  round(exp(pool.added.before$b)*100-100,1),
                  paste(round(exp(pool.added.before$b-1.96*pool.added.before$se)*100-100,1),"to",
                        round(exp(pool.added.before$b+1.96*pool.added.before$se)*100-100,1)),
                  round(pool.added.before$QEp,3))
  
table1.after[1,] <- c(sum(!is.na(added.eff_after[,1])),
                  round(exp(pool.main.after$b)*100-100,1),
                  paste(round(exp(pool.main.after$b-1.96*pool.main.after$se)*100-100,1),"to",
                        round(exp(pool.main.after$b+1.96*pool.main.after$se)*100-100,1)),
                  round(pool.main.after$QEp,3),
                  round(exp(pool.added.after$b)*100-100,1),
                  paste(round(exp(pool.added.after$b-1.96*pool.added.after$se)*100-100,1),"to",
                        round(exp(pool.added.after$b+1.96*pool.added.after$se)*100-100,1)),
                  round(pool.added.after$QEp,3))


table1.before
table1.after



########################################################################################################
# MULTIVARIATE META-ANALYSIS OF THE REDUCED COEF AND THEN COMPUTATION OF BLUP
########################################################################################################

prov_valides <- !is.na(coef1[,1])
dlist_valides <- dlist_total[prov_valides]
provincies_n_valides <- provincies_n_total[prov_valides]
coef1_valides <- coef1[prov_valides,]
coef2_valides <- coef2[prov_valides,]
vcov1_valides <- vcov1[prov_valides]
vcov2_valides <- vcov2[prov_valides]
provinces_valides <- provinces_total[prov_valides]
coef_addeff_1 <- added.eff_before[prov_valides,1]
coef_addeff_2 <- added.eff_after[prov_valides,1]


source("blups_per_attrm_plan.R")

### BLUPS FOR ATTRIBUTABLE MEASURES (WITHOUT ADDED EFFECTS)###
blupsattr<- blups(provinces_valides,provincies_n_valides,dlist_valides,varper,lag,arglag,"adeath",dfseas,dftrend)
blupP1 <- blupsattr[[1]]; blupP2 <- blupsattr[[2]]


### BLUPS FOR ATTRIBUTABLE MEASURES (WITH ADDED EFFECTS)###
blupsattr_plan<- blups_plan(provinces_valides,provincies_n_valides,dlist_valides,varper,lag,arglag,"adeath",dfseas,dftrend)
coef_blupP1_hw_plan <- blupsattr_plan[[1]]; coef_blupP1_cb <- blupsattr_plan[[2]]
coef_blupP2_hw_plan <- blupsattr_plan[[3]]; coef_blupP2_cb <- blupsattr_plan[[4]]
vcov_blupP1_cb <- blupsattr_plan[[6]]; vcov_blupP2_cb <- blupsattr_plan[[8]]
vcov_blupP1_hw_plan <- blupsattr_plan[[5]]
vcov_blupP2_hw_plan <- blupsattr_plan[[7]]

###########################################################################
###########################################################################
## ATTRIBUTABLE MEASURES

source("attrdl_effectiveness.R")
source("attrdl_effectiveness_added_effect.R")

# CREATE THE VECTORS TO STORE THE TOTAL MORTALITY (ACCOUNTING FOR MISSING)
totdeath <- rep(NA,length(provincies_n_valides))
names(totdeath) <- provincies_n_valides

totdeath_1 <- totdeath_2 <- totdeath

# CREATE THE MATRIX TO STORE THE ATTRIBUTABLE INJURIES (simulations)
# Attention: 
#     P1 corresponds to the attributable numbers with risk P1 and 
#        temperatures P2
#     P2 corresponds to the attributable numbers with risk P2 and 
#        temperatures P2
matsim <- matrix(NA,length(provincies_n_valides),4,dimnames=list(provincies_n_valides,
                                                                c("No_added_effect_P1","No_added_effect_P2",
                                                                  "Added_effect_P1", "Added_effect_P2")))

# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE INJURIES
arraysim <- array(NA,dim=c(length(provincies_n_valides),4,nsim),dimnames=list(provincies_n_valides,
                                                   c("No_added_effect_P1","No_added_effect_P2",
                                                     "Added_effect_P1", "Added_effect_P2")))


# Threshold according to the Plan (for each province) - we only consider 35 provinces, 
# where the Plan was activated in both periods
cen <- c(34,36,35,35,33,40,35,30.5,33,38,33,33,35,41,33,32,34,39,35,36,37,34,
         39,33,37,36,33,36.5,36,38,36,37,33,36,33,33,35,33,35,34,41,34,33,35,
         38,34,36,37,35,38)
cen <- cen[c(2:4,6:8,10:18,21:26,28:30,35,37:38,40:41,43,45:47,49:50)]

################################################################################

# RUN THE LOOP
for(i in 1:length(dlist_valides)){
  
    # PRINT
    cat(i,"")
    
    # EXTRACT THE DATA
    data_tot <- dlist_valides[[i]]
    
    
    # COMPUTE THE ATTRIBUTABLE MORTALITY
    # NB: THE REDUCED COEFFICIENTS ARE USED HERE

    ### Without added effect
      # Before and after
    matsim[i,"No_added_effect_P1"] <- attrdl_effectiveness_plan(data_tot,coef1=blupP1[[i]]$blup,
                                                     vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                                     vcov2=blupP2[[i]]$vcov,outcome="adeath",type="an",dir="back",tot=TRUE,
                                                     threshold=cen[i],range=NULL,sim=FALSE)[1]
    matsim[i,"No_added_effect_P2"] <- attrdl_effectiveness_plan(data_tot,coef1=blupP1[[i]]$blup,
                                                             vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                                             vcov2=blupP2[[i]]$vcov,outcome="adeath",type="an",dir="back",tot=TRUE,
                                                             threshold=cen[i],range=NULL,sim=FALSE)[2]
    
    ### With added effect
      # Before and after
    matsim[i,"Added_effect_P1"] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=coef_blupP1_cb[[i]],
                                                     vcov1=vcov_blupP1_cb[[i]],coef2=coef_blupP2_cb[[i]],
                                                     vcov2=vcov_blupP2_cb[[i]],coef_addeff_1=coef_blupP1_hw_plan[[i]],
                                                     coef_addeff_2=coef_blupP2_hw_plan[[i]],
								     vcov_addeff_1=vcov_blupP1_hw_plan[[i]],vcov_addeff_2=vcov_blupP2_hw_plan[[i]],
                                                     outcome="adeath",type="an",dir="back",tot=TRUE,
                                                     threshold=cen[i],range=NULL,sim=FALSE)[1]
    matsim[i,"Added_effect_P2"] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=coef_blupP1_cb[[i]],
                                                                       vcov1=vcov_blupP1_cb[[i]],coef2=coef_blupP2_cb[[i]],
                                                                       vcov2=vcov_blupP2_cb[[i]],coef_addeff_1=coef_blupP1_hw_plan[[i]],
                                                                       coef_addeff_2=coef_blupP2_hw_plan[[i]],
								                       vcov_addeff_1=vcov_blupP1_hw_plan[[i]],vcov_addeff_2=vcov_blupP2_hw_plan[[i]],
                                                                       outcome="adeath",type="an",dir="back",tot=TRUE,
                                                                       threshold=cen[i],range=NULL,sim=FALSE)[2]
    
        
    
    # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE INJURIES
    # USED TO DERIVE CONFIDENCE INTERVALS
    
    # Without added effect
      # Before and after
    arraysim[i,"No_added_effect_P1",] <- attrdl_effectiveness_plan(data_tot,coef1=blupP1[[i]]$blup,
                                                        vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                                        vcov2=blupP2[[i]]$vcov,outcome="adeath",type="an",dir="back",tot=TRUE,
                                                        threshold=cen[i],sim=T,nsim=nsim)[1:nsim]
    arraysim[i,"No_added_effect_P2",] <- attrdl_effectiveness_plan(data_tot,coef1=blupP1[[i]]$blup,
                                                                   vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                                                   vcov2=blupP2[[i]]$vcov,outcome="adeath",type="an",dir="back",tot=TRUE,
                                                                   threshold=cen[i],sim=T,nsim=nsim)[(nsim+1):(nsim*2)]
    
    # With added effect
      # Before and after
    arraysim[i,"Added_effect_P1",] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=coef_blupP1_cb[[i]],
                                                        vcov1=vcov_blupP1_cb[[i]],coef2=coef_blupP2_cb[[i]],
                                                        vcov2=vcov_blupP2_cb[[i]],coef_addeff_1=coef_blupP1_hw_plan[[i]],
                                                        coef_addeff_2=coef_blupP2_hw_plan[[i]],
								     vcov_addeff_1=vcov_blupP1_hw_plan[[i]],vcov_addeff_2=vcov_blupP2_hw_plan[[i]],
                                                        outcome="adeath",type="an",dir="back",tot=TRUE,
                                                        threshold=cen[i],sim=T,nsim=nsim)[1:nsim]
    arraysim[i,"Added_effect_P2",] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=coef_blupP1_cb[[i]],
                                                                             vcov1=vcov_blupP1_cb[[i]],coef2=coef_blupP2_cb[[i]],
                                                                             vcov2=vcov_blupP2_cb[[i]],coef_addeff_1=coef_blupP1_hw_plan[[i]],
                                                                             coef_addeff_2=coef_blupP2_hw_plan[[i]],
  								     				vcov_addeff_1=vcov_blupP1_hw_plan[[i]],vcov_addeff_2=vcov_blupP2_hw_plan[[i]],
                                                                             outcome="adeath",type="an",dir="back",tot=TRUE,
                                                                             threshold=cen[i],sim=T,nsim=nsim)[(nsim+1):(nsim*2)]
    
    
    # STORE THE TOTAL INJURIES (ACCOUNTING FOR MISSING)
    totdeath[i] <- sum(data_tot$adeath,na.rm=T)
    
    totdeath_1[i] <- sum(data_tot$adeath[which(data_tot$bef_aft==0)],na.rm=T) 
    totdeath_2[i] <- sum(data_tot$adeath[which(data_tot$bef_aft==1)],na.rm=T)
    
}

################################################################################
# ATTRIBUTABLE NUMBERS

#### OVERALL PERIOD ####
# CITY-SPECIFIC
ancity <- matsim
ancitylow <- apply(arraysim,c(1,2),quantile,0.025)
ancityhigh <- apply(arraysim,c(1,2),quantile,0.975)
rownames(ancity) <- rownames(ancitylow) <- rownames(ancityhigh) <- provincies_n_valides

# TOTAL
# NB: FIRST SUM THROUGH CITIES
antot <- colSums(matsim)
antotlow <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.025)
antothigh <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.975)


################################################################################
# TOTAL DEATHS

# BY COUNTRY. OVERALL
totdeathtot <- sum(totdeath_2)

# BEFORE PERIOD
totdeathtot_1 <- sum(tempDEATHS$adeath[which(tempDEATHS$bef_aft==0)])
# AFTER PERIOD
totdeathtot_2 <- sum(tempDEATHS$adeath[which(tempDEATHS$bef_aft==1)])

################################################################################
# ATTRIBUTABLE FRACTIONS


#### OVERALL PERIOD ####
# CITY-SPECIFIC
afcity <- ancity/totdeath_2*100
afcitylow <- ancitylow/totdeath_2*100
afcityhigh <- ancityhigh/totdeath_2*100

# TOTAL
aftot <- antot/totdeathtot_2*100
aftotlow <- antotlow/totdeathtot_2*100
aftothigh <- antothigh/totdeathtot_2*100

save(afcity,afcitylow,afcityhigh,file="af_cities_plan.RData")


