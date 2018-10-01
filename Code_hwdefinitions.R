
################################################################################
# "Impact of ambient temperatures on mortality in Spain (1993-2013)"
#   
#   CODE (INTERACTION TERM) AND HEAT WAVES INDICATOR
#
#   ISGlobal  
#   Januray 2018
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


#base_path <- "//fs.isglobal.lan/temperature/Mortality/02_Stata/03_data"
#base_path <- "/media/eloi/Linux Mint 17.1 Xfce 64-bit"
#setwd(paste0(base_path, "/doctorat/Mortality/02_Stata/03_data"))
#load(paste0("G:/doctorat/Mortality/02_Stata/03_data/tempDEATHS.Rdata"))
#load("//fs.isglobal.lan/temperature/Mortality/02_Stata/03_data/tempDEATHS.Rdata")
load("T:\\xavi\\etec\\erica\\mortality\\tempDEATHS.Rdata")


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
#(agafa el data frame i el converteix a llista de tants elements com provincies)
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


# FUNCTION TO CREATE AN HEAT WAVE INDICATOR FOR A TEMPERATURE SERIES
# 	BASED ON THE THRESHOLD AND THE DURATION, BY GROUPS (PROVINCE)
fun.hw.thr <- function(x,thr,dur,group=NULL) {
  as.numeric(apply(Lag(x>=thr,0:(dur-1),group=group),
                   1,sum,na.rm=T)>(dur-1))
}

# CREATE THE MATRICES TO STORE THE RESULTS
hw.N <- hw1.N <- hw2.N <- matrix(NA,length(dlist_total),12, dimnames=list(provincies_n_total,
                                                  paste("hw",rep(c(2,3,4),each=4),rep(c(90,925,95,975),2),sep=".")))

main.eff <- added.eff <- matrix(NA,length(provincies_n_total),24, 
                                dimnames=list(provincies_n_total,paste("hw",rep(c(2,3,4),each=8),
                                                           rep(c(90,925,95,975),each=2),c("est","sd"),sep=".")))
main.eff_before <- added.eff_before <- main.eff_after <- added.eff_after <- main.eff



### CODE FOR DIFFERENT HEAT WAVES DEFINITIONS

for(i in seq(length(dlist_total))) {
  
  # PRINT
  cat(i,"")
  
  # FIRST STAGE
  
  data <- dlist_total[[i]]
  percentiles <- quantile(data$tempmax_compl,c(75,90,92.5,95,97.5)/100,na.rm=T)
  range <- round(range(data$tempmax_compl,na.rm=T),0)
  
  # HW DEFINITIONS
  hw.def <- cbind(rep(percentiles[2:5],3),rep(c(2,3,4),c(4,4,4)))
   
  # RUN THE MODEL FOR EACH DEFINITION
  for(k in 1:nrow(hw.def)) {
    
    data1 <- subset(data, data$yyyy < 2003)
    data2 <- subset(data, data$yyyy > 2003)
    
    # CREATE HEATWAVE INDICATOR FOR THE SPECIFIC HW DEFINITION BEFORE
    hw <- fun.hw.thr(data$tempmax_compl,hw.def[k,1],hw.def[k,2],data$yyyy)
    hw.N[i,k] <- sum(hw)
    
    hw1 <- fun.hw.thr(data1$tempmax_compl,hw.def[k,1],hw.def[k,2],data1$yyyy)
    hw1.N[i,k] <- sum(hw1)
    hw2 <- fun.hw.thr(data2$tempmax_compl,hw.def[k,1],hw.def[k,2],data2$yyyy)
    hw2.N[i,k] <- sum(hw2)
    
    if (sum(hw2)>0 & sum(hw1)>0){
      ###### BEFORE
      # DEFINE THE CROSSBASIS
      argvar <- list(fun="bs",degree=2,knots=quantile(data1$tempmax_compl,varper/100,na.rm=T))
      
      # We define 2 internal knots placed at 75,90th percentiles
      cb1 <- crossbasis(data1$tempmax_compl, lag=lag, argvar=argvar, 
                       arglag=arglag)
      
      data1$t=1:dim(data1)[1]
      
      model1 <- glm(adeath ~ hw1 + cb1 + dow + hday + phday + total_influenza_h +
                      ns (doy,df=dfseas): factor(yyyy) + 
                      ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
                    data1, family=quasipoisson, na.action="na.exclude")
      
      
      ######## AFTER
      # DEFINE THE CROSSBASIS
      argvar <- list(fun="bs",degree=2,knots=quantile(data2$tempmax_compl,varper/100,na.rm=T))
      
      # We define 2 internal knots placed at 75,90th percentiles
      cb2 <- crossbasis(data2$tempmax_compl, lag=lag, argvar=argvar, 
                       arglag=arglag)
      
      data2$t=1:dim(data2)[1]
      
      model2 <- glm(adeath ~ hw2 + cb2 + dow + hday + phday + total_influenza_h +
                      ns (doy,df=dfseas): factor(yyyy) + 
                      ns(t, df=round(length(unique(yyyy))/dftrend/10)), 
                    data2, family=quasipoisson, na.action="na.exclude")
      
      }
      
      # SAVE MAIN EFFECT
      if(sum(hw2)>0 & sum(hw1)>0) {
        tmedian <- median(data$tempmax_compl[hw==1],na.rm=T)
        
        pred.before <- crosspred(cb1,model1,at=c((range[1]+1):(range[2]-1),tmedian),
                          cen=percentiles[1])
        pred.after <- crosspred(cb2,model2,at=c((range[1]+1):(range[2]-1),tmedian),
                          cen=percentiles[1])
        
        main.eff_before[i,c(k*2-1,k*2)] <- cbind(pred.before$allfit,
                                          pred.before$allse)[as.character(tmedian),]
        main.eff_after[i,c(k*2-1,k*2)] <- cbind(pred.after$allfit,
                                          pred.after$allse)[as.character(tmedian),]
       # SAVE ADDED EFFECT
        added.eff_before[i,c(k*2-1,k*2)] <- ci.lin(model1)["hw1",1:2]
        added.eff_after[i,c(k*2-1,k*2)] <- ci.lin(model2)["hw2",1:2]
      
      } else {main.eff_before[i,c(k*2-1,k*2)] <- c(NA,NA)
              main.eff_after[i,c(k*2-1,k*2)] <- c(NA,NA)
              }
      
    }

}  


########################################################################################################
# SECOND-STAGE ANALYSIS: UNIVARIATE META-ANALYSIS
########################################################################################################

label <- paste("hw",rep(c(2,3,4),each=4),rep(c(90,92.5,95,97.5),2),sep=".")
table1.before <- table1.after <- matrix(NA,12,7,dimnames=list(label,
                                      c("N comm","Est.main","95%CI.main","P-het.added","Est.added",
                                        "95%CI.added","P-het.added")))

for(i in 1:12) {	 
  
  # SET TO MISSING IF NO ESTIMATE FOR ADDED EFFECT
  added.eff_before[added.eff_before[,2*i]==0,c(2*i-1,2*i)] <- NA
  main.eff_before[is.na(added.eff_before[,2*i]),c(2*i-1,2*i)] <- NA
  added.eff_after[added.eff_after[,2*i]==0,c(2*i-1,2*i)] <- NA
  main.eff_after[is.na(added.eff_after[,2*i]),c(2*i-1,2*i)] <- NA
  
  
  # RUN THE META-ANALYSIS
  pool.main.before <- rma.uni(yi=main.eff_before[,2*i-1],sei=main.eff_before[,2*i])
  pool.added.before <- rma.uni(yi=added.eff_before[,2*i-1],sei=added.eff_before[,2*i])
  
  pool.main.after <- rma.uni(yi=main.eff_after[,2*i-1],sei=main.eff_after[,2*i])
  pool.added.after <- rma.uni(yi=added.eff_after[,2*i-1],sei=added.eff_after[,2*i])
  
  # FILL TABLE1
  table1.before[i,] <- c(sum(!is.na(added.eff_before[,2*i-1])),
                  round(exp(pool.main.before$b)*100-100,1),
                  paste(round(exp(pool.main.before$b-1.96*pool.main.before$se)*100-100,1),"to",
                        round(exp(pool.main.before$b+1.96*pool.main.before$se)*100-100,1)),
                  round(pool.main.before$QEp,3),
                  round(exp(pool.added.before$b)*100-100,1),
                  paste(round(exp(pool.added.before$b-1.96*pool.added.before$se)*100-100,1),"to",
                        round(exp(pool.added.before$b+1.96*pool.added.before$se)*100-100,1)),
                  round(pool.added.before$QEp,3))
  
  table1.after[i,] <- c(sum(!is.na(added.eff_after[,2*i-1])),
                  round(exp(pool.main.after$b)*100-100,1),
                  paste(round(exp(pool.main.after$b-1.96*pool.main.after$se)*100-100,1),"to",
                        round(exp(pool.main.after$b+1.96*pool.main.after$se)*100-100,1)),
                  round(pool.main.after$QEp,3),
                  round(exp(pool.added.after$b)*100-100,1),
                  paste(round(exp(pool.added.after$b-1.96*pool.added.after$se)*100-100,1),"to",
                        round(exp(pool.added.after$b+1.96*pool.added.after$se)*100-100,1)),
                  round(pool.added.after$QEp,3))
}



table1.before
table1.after



########################################################################################################
# ATTRIBUTABLE MEASURES (WITH ADDED EFFECTS)
########################################################################################################

#setwd("//fs.isglobal.lan/temperature/Mortality/02_Stata/01_do/R code/Temperatures_mortality-master/Functions - heatwaves/Erica")
#setwd("G:/doctorat/Mortality/02_Stata/01_do/R code/Temperatures_mortality-master/Functions - heatwaves/Erica")
setwd("T:\\Xavi\\ETEC\\Erica\\mortality\\separate models\\heatwaves\\erica")
setwd("H:/doctorat/Mortality/02_Stata/01_do/R code/Temperatures_mortality-master/separate models/heat waves/erica")
source("04_blups_per_attrm_hwdefinitions.R")

# per alguna rao no em funciona tal qual però si corro el codi de la funció si

blupsattr_plan <- blups_plan(provinces_total,provincies_n_total,dlist_total,varper,lag,arglag,"adeath",dfseas,dftrend)

save(blupsattr_plan,file="T:\\Xavi\\ETEC\\Erica\\mortality\\separate models\\heatwaves\\erica\\blupsattr_plan.RData")
save(blupsattr_plan,file="H:/doctorat/Mortality/02_Stata/01_do/R code/Temperatures_mortality-master/separate models/heat waves/erica/blupsattr_plan.RData")

# load(file="T:\\Xavi\\ETEC\\Erica\\mortality\\separate models\\heatwaves\\erica\\blupsattr_plan.RData")

coef_blupP1_hw_1 <- blupsattr_plan[[1]]; coef_blupP1_cb_hw_1 <- blupsattr_plan[[2]]
coef_blupP2_hw_1 <- blupsattr_plan[[3]]; coef_blupP2_cb_hw_1 <- blupsattr_plan[[4]]
vcov_blupP1_cb_hw_1 <- blupsattr_plan[[6]]; vcov_blupP2_cb_hw_1 <- blupsattr_plan[[8]]

coef_blupP1_hw_2 <- blupsattr_plan[[9]]; coef_blupP1_cb_hw_2 <- blupsattr_plan[[10]]
coef_blupP2_hw_2 <- blupsattr_plan[[11]]; coef_blupP2_cb_hw_2 <- blupsattr_plan[[12]]
vcov_blupP1_cb_hw_2 <- blupsattr_plan[[14]]; vcov_blupP2_cb_hw_2 <- blupsattr_plan[[16]]

coef_blupP1_hw_3 <- blupsattr_plan[[17]]; coef_blupP1_cb_hw_3 <- blupsattr_plan[[18]]
coef_blupP2_hw_3 <- blupsattr_plan[[19]]; coef_blupP2_cb_hw_3 <- blupsattr_plan[[20]]
vcov_blupP1_cb_hw_3 <- blupsattr_plan[[22]]; vcov_blupP2_cb_hw_3 <- blupsattr_plan[[24]]

coef_blupP1_hw_4 <- blupsattr_plan[[25]]; coef_blupP1_cb_hw_4 <- blupsattr_plan[[26]]
coef_blupP2_hw_4 <- blupsattr_plan[[27]]; coef_blupP2_cb_hw_4 <- blupsattr_plan[[28]]
vcov_blupP1_cb_hw_4 <- blupsattr_plan[[30]]; vcov_blupP2_cb_hw_4 <- blupsattr_plan[[32]]

coef_blupP1_hw_5 <- blupsattr_plan[[33]]; coef_blupP1_cb_hw_5 <- blupsattr_plan[[34]]
coef_blupP2_hw_5 <- blupsattr_plan[[35]]; coef_blupP2_cb_hw_5 <- blupsattr_plan[[36]]
vcov_blupP1_cb_hw_5 <- blupsattr_plan[[38]]; vcov_blupP2_cb_hw_5 <- blupsattr_plan[[40]]

coef_blupP1_hw_6 <- blupsattr_plan[[41]]; coef_blupP1_cb_hw_6 <- blupsattr_plan[[42]]
coef_blupP2_hw_6 <- blupsattr_plan[[43]]; coef_blupP2_cb_hw_6 <- blupsattr_plan[[44]]
vcov_blupP1_cb_hw_6 <- blupsattr_plan[[46]]; vcov_blupP2_cb_hw_6 <- blupsattr_plan[[48]]

coef_blupP1_hw_7 <- blupsattr_plan[[49]]; coef_blupP1_cb_hw_7 <- blupsattr_plan[[50]]
coef_blupP2_hw_7 <- blupsattr_plan[[51]]; coef_blupP2_cb_hw_7 <- blupsattr_plan[[52]]
vcov_blupP1_cb_hw_7 <- blupsattr_plan[[54]]; vcov_blupP2_cb_hw_7 <- blupsattr_plan[[56]]

coef_blupP1_hw_8 <- blupsattr_plan[[57]]; coef_blupP1_cb_hw_8 <- blupsattr_plan[[58]]
coef_blupP2_hw_8 <- blupsattr_plan[[59]]; coef_blupP2_cb_hw_8 <- blupsattr_plan[[60]]
vcov_blupP1_cb_hw_8 <- blupsattr_plan[[62]]; vcov_blupP2_cb_hw_8 <- blupsattr_plan[[64]]

coef_blupP1_hw_9 <- blupsattr_plan[[65]]; coef_blupP1_cb_hw_9 <- blupsattr_plan[[66]]
coef_blupP2_hw_9 <- blupsattr_plan[[67]]; coef_blupP2_cb_hw_9 <- blupsattr_plan[[68]]
vcov_blupP1_cb_hw_9 <- blupsattr_plan[[70]]; vcov_blupP2_cb_hw_9 <- blupsattr_plan[[72]]

coef_blupP1_hw_10 <- blupsattr_plan[[73]]; coef_blupP1_cb_hw_10 <- blupsattr_plan[[74]]
coef_blupP2_hw_10 <- blupsattr_plan[[75]]; coef_blupP2_cb_hw_10 <- blupsattr_plan[[76]]
vcov_blupP1_cb_hw_10 <- blupsattr_plan[[78]]; vcov_blupP2_cb_hw_10 <- blupsattr_plan[[80]]

coef_blupP1_hw_11 <- blupsattr_plan[[81]]; coef_blupP1_cb_hw_11 <- blupsattr_plan[[82]]
coef_blupP2_hw_11 <- blupsattr_plan[[83]]; coef_blupP2_cb_hw_11 <- blupsattr_plan[[84]]
vcov_blupP1_cb_hw_11 <- blupsattr_plan[[86]]; vcov_blupP2_cb_hw_11 <- blupsattr_plan[[88]]

coef_blupP1_hw_12 <- blupsattr_plan[[89]]; coef_blupP1_cb_hw_12 <- blupsattr_plan[[90]]
coef_blupP2_hw_12 <- blupsattr_plan[[91]]; coef_blupP2_cb_hw_12 <- blupsattr_plan[[92]]
vcov_blupP1_cb_hw_12 <- blupsattr_plan[[94]]; vcov_blupP2_cb_hw_12 <- blupsattr_plan[[96]]

# Added by XB:

vcov_blupP1_hw_1 <- blupsattr_plan[[5]]
vcov_blupP1_hw_2 <- blupsattr_plan[[13]]
vcov_blupP1_hw_3 <- blupsattr_plan[[21]]
vcov_blupP1_hw_4 <- blupsattr_plan[[29]]
vcov_blupP1_hw_5 <- blupsattr_plan[[37]]
vcov_blupP1_hw_6 <- blupsattr_plan[[45]]
vcov_blupP1_hw_7 <- blupsattr_plan[[53]]
vcov_blupP1_hw_8 <- blupsattr_plan[[61]]
vcov_blupP1_hw_9 <- blupsattr_plan[[69]]
vcov_blupP1_hw_10 <- blupsattr_plan[[77]]
vcov_blupP1_hw_11 <- blupsattr_plan[[85]]
vcov_blupP1_hw_12 <- blupsattr_plan[[93]]


vcov_blupP2_hw_1 <- blupsattr_plan[[7]]
vcov_blupP2_hw_2 <- blupsattr_plan[[15]]
vcov_blupP2_hw_3 <- blupsattr_plan[[23]]
vcov_blupP2_hw_4 <- blupsattr_plan[[31]]
vcov_blupP2_hw_5 <- blupsattr_plan[[39]]
vcov_blupP2_hw_6 <- blupsattr_plan[[47]]
vcov_blupP2_hw_7 <- blupsattr_plan[[55]]
vcov_blupP2_hw_8 <- blupsattr_plan[[63]]
vcov_blupP2_hw_9 <- blupsattr_plan[[71]]
vcov_blupP2_hw_10 <- blupsattr_plan[[79]]
vcov_blupP2_hw_11 <- blupsattr_plan[[87]]
vcov_blupP2_hw_12 <- blupsattr_plan[[95]]


# changed by XB:
#source("attrdl_effectiveness_added_effect_hwdef.R")
# XB: I made several changes to the function:
source("attrdl_effectiveness_added_effect_hwdef_XB.R")

# CREATE THE VECTORS TO STORE THE TOTAL MORTALITY (ACCOUNTING FOR MISSING)
totdeath <- rep(NA,length(provincies_n_total))
names(totdeath) <- provincies_n_total

totdeath_1 <- totdeath_2 <- totdeath

# CREATE THE MATRIX TO STORE THE ATTRIBUTABLE INJURIES (simulations)
# Attention: 
#     P1 corresponds to the attributable numbers with risk P1 and 
#        temperatures P2
#     P2 corresponds to the attributable numbers with risk P2 and 
#        temperatures P2

matsim <- lapply(1:12, function(x) matrix(NA, nrow=50, ncol=2, dimnames=list(provincies_n_total,
                                                                             c("Added_effect_P1", "Added_effect_P2"))))

# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE INJURIES
# arraysim1 <- array(NA,dim=c(length(provincies_n_total),2,nsim),dimnames=list(provincies_n_total,
#                                                                               c("Added_effect_P1", "Added_effect_P2")))
# arraysim2 <- arraysim3 <- arraysim4 <- arraysim5 <- arraysim6 <- arraysim7 <- arraysim8 <- arraysim9 <- 
#   arraysim10 <- arraysim11 <- arraysim12 <- arraysim1
arraysim <- lapply(1:12, function(x) array(NA, dim=c(length(provincies_n_total),2,nsim),dimnames=list(provincies_n_total,
                                                                              c("Added_effect_P1", "Added_effect_P2"))))


# RUN THE LOOP
for(i in 1:length(dlist_total)){
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data_tot <- dlist_total[[i]]
  
#  percentiles <- quantile(data$tempmax_compl,c(75,90,92.5,95,97.5)/100,na.rm=T)
#  range <- round(range(data$tempmax_compl,na.rm=T),0)
 
# changed by XB:
   percentiles <- quantile(data_tot$tempmax_compl,c(75,90,92.5,95,97.5)/100,na.rm=T)
   range <- round(range(data_tot$tempmax_compl,na.rm=T),0)
 
  # HW DEFINITIONS
#  hw.def <- cbind(rep(percentiles[2:5],2),rep(c(2,3,4),c(4,4,4)))

# Changed by XB:
   hw.def <- cbind(rep(percentiles[2:5],3),rep(c(2,3,4),c(4,4,4)))
 

  # RUN THE MODEL FOR EACH DEFINITION
  for(k in 1:nrow(hw.def)) {
  
    # COMPUTE THE ATTRIBUTABLE MORTALITY
    # NB: THE REDUCED COEFFICIENTS ARE USED HERE
    
     
    ### With added effect
    # Before and after

    # XB: in all calls to attrdl_effectiveness_plan_added_effect I made changes in the hw.def argument
    matsim[[k]][i,"Added_effect_P1"] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=get(paste0("coef_blupP1_cb_hw_",k))[[i]],
                                                     vcov1=get(paste0("vcov_blupP1_cb_hw_",k))[[i]],
                                                     coef2=get(paste0("coef_blupP2_cb_hw_",k))[[i]],
                                                     vcov2=get(paste0("vcov_blupP2_cb_hw_",k))[[i]],
                                                     coef_addeff_1=get(paste0("coef_blupP1_hw_",k))[[i]],
                                                     coef_addeff_2=get(paste0("coef_blupP2_hw_",k))[[i]],
								     vcov_addeff_1=get(paste0("vcov_blupP1_hw_",k))[[i]],
								     vcov_addeff_2=get(paste0("vcov_blupP2_hw_",k))[[i]],
                                                     outcome="adeath",type="an",dir="back",tot=TRUE,
                                                     threshold=hw.def[k,1], hw.def=hw.def[k,], range=NULL,sim=FALSE)[1]
    
    
    matsim[[k]][i,"Added_effect_P2"] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=get(paste0("coef_blupP1_cb_hw_",k))[[i]],
                                                                               vcov1=get(paste0("vcov_blupP1_cb_hw_",k))[[i]],
                                                                               coef2=get(paste0("coef_blupP2_cb_hw_",k))[[i]],
                                                                               vcov2=get(paste0("vcov_blupP2_cb_hw_",k))[[i]],
                                                                               coef_addeff_1=get(paste0("coef_blupP1_hw_",k))[[i]],
                                                                               coef_addeff_2=get(paste0("coef_blupP2_hw_",k))[[i]],
								     					 vcov_addeff_1=get(paste0("vcov_blupP1_hw_",k))[[i]],
								     					 vcov_addeff_2=get(paste0("vcov_blupP2_hw_",k))[[i]],
                                                                               outcome="adeath",type="an",dir="back",tot=TRUE,
                                                                               threshold=hw.def[k,1], hw.def=hw.def[k,], range=NULL,sim=FALSE)[2]
    
    
    
    # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE INJURIES
    # USED TO DERIVE CONFIDENCE INTERVALS
    
    # With added effect
    # Before and after
    arraysim[[k]][i,"Added_effect_P1",] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=get(paste0("coef_blupP1_cb_hw_",k))[[i]],
                                                                             vcov1=get(paste0("vcov_blupP1_cb_hw_",k))[[i]],
                                                                             coef2=get(paste0("coef_blupP2_cb_hw_",k))[[i]],
                                                                             vcov2=get(paste0("vcov_blupP2_cb_hw_",k))[[i]],
                                                                             coef_addeff_1=get(paste0("coef_blupP1_hw_",k))[[i]],
                                                                             coef_addeff_2=get(paste0("coef_blupP2_hw_",k))[[i]],
								     				     vcov_addeff_1=get(paste0("vcov_blupP1_hw_",k))[[i]],
								     				     vcov_addeff_2=get(paste0("vcov_blupP2_hw_",k))[[i]],
                                                                             outcome="adeath",type="an",dir="back",tot=TRUE,
                                                                             threshold=hw.def[k,1],hw.def=hw.def[k,],sim=T,nsim=nsim)[1:nsim]
    
    
    arraysim[[k]][i,"Added_effect_P2",] <- attrdl_effectiveness_plan_added_effect(data_tot,coef1=get(paste0("coef_blupP1_cb_hw_",k))[[i]],
                                                                             vcov1=get(paste0("vcov_blupP1_cb_hw_",k))[[i]],
                                                                             coef2=get(paste0("coef_blupP2_cb_hw_",k))[[i]],
                                                                             vcov2=get(paste0("vcov_blupP2_cb_hw_",k))[[i]],
                                                                             coef_addeff_1=get(paste0("coef_blupP1_hw_",k))[[i]],
                                                                             coef_addeff_2=get(paste0("coef_blupP2_hw_",k))[[i]],
								     				     vcov_addeff_1=get(paste0("vcov_blupP1_hw_",k))[[i]],
								     				     vcov_addeff_2=get(paste0("vcov_blupP2_hw_",k))[[i]],
                                                                             outcome="adeath",type="an",dir="back",tot=TRUE,
                                                                             threshold=hw.def[k,1],hw.def=hw.def[k,],sim=T,nsim=nsim)[(nsim+1):(nsim*2)]
    
    
    # STORE THE TOTAL INJURIES (ACCOUNTING FOR MISSING)
    totdeath[i] <- sum(data_tot$adeath,na.rm=T)
    
    totdeath_1[i] <- sum(data_tot$adeath[which(data_tot$bef_aft==0)],na.rm=T) 
    totdeath_2[i] <- sum(data_tot$adeath[which(data_tot$bef_aft==1)],na.rm=T)
  }
  
}


################################################################################
# ATTRIBUTABLE NUMBERS

#### OVERALL PERIOD ####
# CITY-SPECIFIC
ancity <- matsim
# changed by XB:
ancitylow=list()
for (i in 1:12) {
  ancitylow[[i]] <- apply(arraysim[[i]],c(1,2),quantile,0.025)
}
ancityhigh=list()
for (i in 1:12) {
  ancityhigh[[i]] <- apply(arraysim[[i]],c(1,2),quantile,0.975)
}

#rownames(ancity) <- rownames(ancitylow) <- rownames(ancityhigh) <- provincies_n_total

# TOTAL
# NB: FIRST SUM THROUGH CITIES

# changed by XB:
antot=list()
for (i in 1:12) {
  antot[[i]]=apply(ancity[[i]],2,sum)
}
#antot <- lapply(matsim, function(x) sum(x))
antotlow=list()
for (i in 1:12) {
  antotlow[[i]] <- apply(apply(arraysim[[i]],c(2,3),sum),1,quantile,0.025)
}
antothigh=list()
for (i in 1:12) {
 antothigh[[i]] <- apply(apply(arraysim[[i]],c(2,3),sum),1,quantile,0.975)
}

################################################################################
# TOTAL INJURIES

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
# changed by XB:

afcity=list()
for (i in 1:12) {
  afcity[[i]] <- ancity[[i]]/totdeath_2*100
}
afcitylow=list()
for (i in 1:12) {
  afcitylow[[i]] <- ancitylow[[i]]/totdeath_2*100
}

afcityhigh=list()
for (i in 1:12) {
  afcityhigh[[i]] <- ancityhigh[[i]]/totdeath_2*100
}

# TOTAL
aftot = list()
for (i in 1:12) {
  aftot[[i]] <- antot[[i]]/totdeathtot_2*100
}
aftotlow = list()
for (i in 1:12) {
  aftotlow[[i]] <- antotlow[[i]]/totdeathtot_2*100
}
aftothigh=list()
for (i in 1:12) {
  aftothigh[[i]] <- antothigh[[i]]/totdeathtot_2*100
}

mat.af=round(matrix(unlist(aftot),ncol=2,byrow=T),3)
mat.af.low=round(matrix(unlist(aftotlow),ncol=2,byrow=T),3)
mat.af.high=round(matrix(unlist(aftothigh),ncol=2,byrow=T),3)

res.before=matrix(paste0(mat.af[,1]," (",mat.af.low[,1],", ",mat.af.high[,1],")"),ncol=1)
res.after=matrix(paste0(mat.af[,2]," (",mat.af.low[,2],", ",mat.af.high[,2],")"),ncol=1)

print.data.frame(data.frame(res.before),quote=F,row.names=F)
print.data.frame(data.frame(res.after),quote=F,row.names=F)






