
################################################################################
# "Temporal changes in temperature-related mortality in Spain and effect of 
#        the implementation of a Heat Health Prevention Plan"
#   
#   ISGlobal  
#   June 2018
#   
#
################################################################################

# ATTRIBUTABLE MEASURES
# COMPUTE THE ATTRIBUTABLE INJURIES FOR EACH CITY, WITH EMPIRICAL CI
# ESTIMATED USING THE RE-CENTERED BASES
################################################################################

# LOAD THE FUNCTION FOR COMPUTING THE ATTRIBUTABLE RISK MEASURES
source("attrdl_effectiveness.R")

# open this file, result of running "blups per attrm.R"
load(file="blupP.RData")

library(mgcv)

# CREATE THE VECTORS TO STORE THE TOTAL MORTALITY (ACCOUNTING FOR MISSING)
totdeath <- rep(NA,length(provinces))
names(totdeath) <- provincies_n

totdeath_1 <- totdeath_2 <- totdeath

# CREATE THE MATRIX TO STORE THE ATTRIBUTABLE INJURIES (simulations)
# Attention: 
#     P1 corresponds to the attributable numbers with risk P1 and 
#        temperatures P2
#     P2 corresponds to the attributable numbers with risk P2 and 
#        temperatures P2
matsim <- matrix(NA,length(provinces),14,dimnames=list(provincies_n,
          c("glob_P1","cold_P1","heat_P1","extreme cold_P1","moderate cold_P1",
            "moderate heat_P1","extreme heat_P1","glob_P2","cold_P2","heat_P2",
            "extreme cold_P2","moderate cold_P2","moderate heat_P2",
            "extreme heat_P2")))

# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE INJURIES
arraysim <- array(NA,dim=c(length(provinces),14,nsim),dimnames=list(provincies_n,
            c("glob_P1","cold_P1","heat_P1","extreme cold_P1","moderate cold_P1",
              "moderate heat_P1","extreme heat_P1","glob_P2","cold_P2","heat_P2",
              "extreme cold_P2","moderate cold_P2","moderate heat_P2",
              "extreme heat_P2")))


## Calculate minimum temperature mortality for period 2
minperccityP2 <- mintempcityP2 <- rep(NA,length(dlist))
for(i in seq(length(dlist))) {
  data <- dlist[[i]]
  predvar <- quantile(data[data$bef_aft==1,]$tempmax_compl,1:99/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  
  # We define 3 knots at the 10th, 75th and 90th percentiles
  argvar <- list(x=predvar,fun="bs",
                 degree=2, knots=quantile(data[data$bef_aft==1,]$tempmax_compl, varper/100, na.rm=T),
                 Bound=range(data[data$bef_aft==1,]$tempmax_compl,na.rm=T))
  
  bvar <- do.call(onebasis,argvar)
  minperccityP2[i] <- (1:99)[which.min((bvar%*%blup2[[i]]$blup))]
  mintempcityP2[i] <- quantile(data$tempmax_compl,minperccityP2[i]/100,na.rm=T)
}

################################################################################

# RUN THE LOOP
for(i in 1:length(dlist)){
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data_tot <- dlist[[i]]
  
  # COMPUTE PERCENTILES 2.5 AND 97.5
  perc <- quantile(data_tot[data$bef_aft==1,]$tempmax_compl,c(0.025,0.975))
  
  # COMPUTE THE ATTRIBUTABLE MORTALITY
  # NB: THE REDUCED COEFFICIENTS ARE USED HERE
  
  #### COMPARING BEFORE AND AFTER PERIODS
  matsim[i,"glob_P1"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=NULL,sim=FALSE)[1]
  matsim[i,"glob_P2"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=NULL,sim=FALSE)[2]
  
  
  matsim[i,"cold_P1"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=c(-100,mintempcityP2[[i]]),
                            sim=FALSE)[1]
  matsim[i,"cold_P2"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=c(-100,mintempcityP2[[i]]),
                            sim=FALSE)[2]
  
  
  
  matsim[i,"heat_P1"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],100),
                            sim=FALSE)[1]
  matsim[i,"heat_P2"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],100),
                            sim=FALSE)[2]
  
  matsim[i,"extreme cold_P1"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                cen=mintempcityP2[[i]],range=c(-100,perc[1]),
                                sim=FALSE)[1]
  matsim[i,"extreme cold_P2"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                  vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                  vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                  cen=mintempcityP2[[i]],range=c(-100,perc[1]),
                                  sim=FALSE)[2]  

  matsim[i,"moderate cold_P1"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                  vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                  vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                  cen=mintempcityP2[[i]],range=c(perc[1],mintempcityP2[[i]]),
                                  sim=FALSE)[1]
  matsim[i,"moderate cold_P2"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                    vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                    vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                    cen=mintempcityP2[[i]],range=c(perc[1],mintempcityP2[[i]]),
                                    sim=FALSE)[2]  
    
  matsim[i,"moderate heat_P1"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                     vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                     vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                     cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],perc[2]),
                                     sim=FALSE)[1]
  matsim[i,"moderate heat_P2"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                     vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                     vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                     cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],perc[2]),
                                     sim=FALSE)[2]
    
  matsim[i,"extreme heat_P1"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                   vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                   vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                   cen=mintempcityP2[[i]],range=c(perc[2],100),
                                   sim=FALSE)[1]
  matsim[i,"extreme heat_P2"] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                   vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                   vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                   cen=mintempcityP2[[i]],range=c(perc[2],100),
                                   sim=FALSE)[2]
    

  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE INJURIES
  # USED TO DERIVE CONFIDENCE INTERVALS
  
  
  #### OVERALL PERIOD: 1993-2013 (EXCLUDING 2003)
  arraysim[i,"glob_P1",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                               vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                               vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                               cen=mintempcityP2[[i]],sim=T,nsim=nsim)[1:nsim]
  arraysim[i,"glob_P2",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                              vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                              vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                              cen=mintempcityP2[[i]],sim=T,nsim=nsim)[(nsim+1):(nsim*2)]
    
  arraysim[i,"cold_P1",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=c(-100,mintempcityP2[[i]]),
                            sim=T,nsim=nsim)[1:nsim]
  arraysim[i,"cold_P2",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                               vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                               vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                               cen=mintempcityP2[[i]],range=c(-100,mintempcityP2[[i]]),
                               sim=T,nsim=nsim)[(nsim+1):(nsim*2)]

  arraysim[i,"heat_P1",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                            vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                            vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                            cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],100),
                            sim=T,nsim=nsim)[1:nsim]
  arraysim[i,"heat_P2",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],100),
                                sim=T,nsim=nsim)[(nsim+1):(nsim*2)]
    
  arraysim[i,"extreme cold_P1",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                  vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                  vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                  cen=mintempcityP2[[i]],range=c(-100,perc[1]),
                                  sim=T,nsim=nsim)[1:nsim]
  arraysim[i,"extreme cold_P2",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                      vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                      vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                      cen=mintempcityP2[[i]],range=c(-100,perc[1]),
                                      sim=T,nsim=nsim)[(nsim+1):(nsim*2)]

  arraysim[i,"moderate cold_P1",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                       vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                       vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                       cen=mintempcityP2[[i]],range=c(perc[1],mintempcityP2[[i]]),
                                       sim=T,nsim=nsim)[1:nsim]
  arraysim[i,"moderate cold_P2",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                                          vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                                          vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                                          cen=mintempcityP2[[i]],range=c(perc[1],mintempcityP2[[i]]),
                                                          sim=T,nsim=nsim)[(nsim+1):(nsim*2)]

  arraysim[i,"moderate heat_P1",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                         vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                         vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                         cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],perc[2]),
                                         sim=T,nsim=nsim)[1:nsim]
  arraysim[i,"moderate heat_P2",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                         vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                         vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                         cen=mintempcityP2[[i]],range=c(mintempcityP2[[i]],perc[2]),
                                         sim=T,nsim=nsim)[(nsim+1):(nsim*2)]    
      
  arraysim[i,"extreme heat_P1",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                   vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                   vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                   cen=mintempcityP2[[i]],range=c(perc[2],100),
                                   sim=T,nsim=nsim)[1:nsim]
  arraysim[i,"extreme heat_P2",] <- attrdl_effectiveness(data_tot,coef1=blupP1[[i]]$blup,
                                       vcov1=blupP1[[i]]$vcov,coef2=blupP2[[i]]$blup,
                                       vcov2=blupP2[[i]]$vcov,type="an",dir="back",tot=TRUE,
                                       cen=mintempcityP2[[i]],range=c(perc[2],100),
                                       sim=T,nsim=nsim)[(nsim+1):(nsim*2)]
    

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
rownames(ancity) <- rownames(ancitylow) <- rownames(ancityhigh) <- provincies_n

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

#
## DIFFERENCES IN ATTR FRACTION (P2-P1)
af_differences <- matrix(NA,length(provinces),7,dimnames=list(provincies_n,
                            c("glob","cold","heat","extreme cold","moderate cold",
                            "moderate heat","extreme heat")))
af_differences[,1] <- afcity[,c(8)]-afcity[,c(1)]
af_differences[,2] <- afcity[,c(9)]-afcity[,c(2)]
af_differences[,3] <- afcity[,c(10)]-afcity[,c(3)]
af_differences[,4] <- afcity[,c(11)]-afcity[,c(4)]
af_differences[,5] <- afcity[,c(12)]-afcity[,c(5)]
af_differences[,6] <- afcity[,c(13)]-afcity[,c(6)] 
af_differences[,7] <- afcity[,c(14)]-afcity[,c(7)]


save(af_differences,ancity,ancitylow,ancityhigh,afcity,afcitylow,afcityhigh,
     aftot,aftothigh,aftotlow, afcity,afcitylow,afcityhigh,file="attrfraction_city_differences.csv")
load(file="attrfraction_city_differences.csv")


> row.names(af_differences)
 [1] "Alava"                  "Albacete"               "Alicante"              
 [4] "Almeria"                "Avila"                  "Badajoz"               
 [7] "Illes Balears"          "Barcelona"              "Burgos"                
[10] "Caceres"                "Cadiz"                  "Castellon"             
[13] "Ciudad Real"            "Cordoba"                "A Coruna"              
[16] "Cuenca"                 "Girona"                 "Granada"               
[19] "Guadalajara"            "Guipuzcoa"              "Huelva"                
[22] "Huesca"                 "Jaen"                   "Leon"                  
[25] "Lleida"                 "La Rioja"               "Lugo"                  
[28] "Madrid"                 "Malaga"                 "Murcia"                
[31] "Navarra"                "Ourense"                "Asturias"              
[34] "Palencia"               "Las Palmas"             "Pontevedra"            
[37] "Salamanca"              "Santa Cruz de Tenerife" "Cantabria"             
[40] "Segovia"                "Sevilla"                "Soria"                 
[43] "Tarragona"              "Teruel"                 "Toledo"                
[46] "Valencia"               "Valladolid"             "Vizcaya"               
[49] "Zamora"                 "Zaragoza"   

af_differences.df=data.frame(af_differences)
af_differences.df$ccaa=c(16,8,10,
	1,7,11,
	4,9,7,
	11,1,10,
	8,1,12,
	8,9,1,
	8,16,1,
	2,1,7,
	9,17,12,
	13,1,14,
	15,12,3,
	7,5,12,
	7,5,6,
	7,1,7,
	9,2,8,
	10,7,16,
	7,2)
af_differences.df$p1.1=NA
af_differences.df$p1.2=NA
af_differences.df$p1.3=NA
af_differences.df$p1.4=NA
af_differences.df$p1.5=NA
af_differences.df$p3.1=NA
af_differences.df$p3.2=NA
af_differences.df$p3.3=NA
af_differences.df$p3.4=NA
af_differences.df$p3.5=NA
af_differences.df$p4.1=NA
af_differences.df$p4.2=NA
af_differences.df$p4.3=NA
af_differences.df$p4.4=NA
af_differences.df$p4.5=NA
af_differences.df$p5.1=NA
af_differences.df$p5.2=NA
af_differences.df$p5.3=NA
af_differences.df$p5.4=NA
af_differences.df$p5.5=NA
af_differences.df$p6.1=NA
af_differences.df$p6.2=NA
af_differences.df$p6.3=NA
af_differences.df$p6.4=NA
af_differences.df$p6.5=NA
af_differences.df$p8.1=NA
af_differences.df$p8.2=NA
af_differences.df$p8.3=NA
af_differences.df$p8.4=NA
af_differences.df$p8.5=NA


af_differences.df[af_differences.df$ccaa==1,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,rep(0,4)),c(1,1,0,1,1),rep(1,5),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==1,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==1,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==2,9:38]=matrix(rep(c(c(0,0,0,0,1),rep(0,5),rep(0,5),rep(0,5),rep(0,5),rep(0,5)),dim(af_differences.df[af_differences.df$ccaa==2,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==2,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==3,9:38]=matrix(rep(c(c(0,1,0,1,1),c(0,1,0,0,0),c(1,0,0,0,0),rep(0,5),rep(0,5),c(1,1,1,0,0)),dim(af_differences.df[af_differences.df$ccaa==3,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==3,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==4,9:38]=matrix(rep(c(c(1,1,0,1,1),c(1,0,1,0,1),c(1,0,0,0,0),rep(0,5),rep(0,5),c(1,1,0,0,0)),dim(af_differences.df[af_differences.df$ccaa==4,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==4,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==5,9:38]=matrix(rep(c(c(0,0,0,0,1),c(1,0,0,0,0),c(1,0,0,0,0),c(1,0,0,0,0),rep(0,5),c(1,1,0,0,0)),dim(af_differences.df[af_differences.df$ccaa==5,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==5,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==6,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,0,0,0,0),c(1,1,0,1,1),rep(1,5),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==6,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==6,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==7,9:38]=matrix(rep(c(c(0,0,0,0,1),c(1,0,0,0,0),c(1,0,0,0,0),c(1,0,0,0,0),rep(0,5),c(1,1,0,0,0)),dim(af_differences.df[af_differences.df$ccaa==7,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==7,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==8,9:38]=matrix(rep(c(rep(1,5),c(1,1,1,0,0),c(1,0,0,0,0),c(1,1,0,1,1),c(1,0,1,1,1),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==8,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==8,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==9,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,0,0,0,0),c(1,1,0,1,1),rep(1,5),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==9,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==9,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==10,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,1,0,0,0),c(1,1,0,1,1),rep(1,5),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==10,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==10,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==11,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,0,0,0,0),c(1,0,0,1,0),c(0,1,0,0,0),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==11,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==11,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==12,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,0,0,0,0),c(1,1,0,1,1),rep(1,5),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==12,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==12,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==13,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,0,0,0,0),c(1,1,0,1,1),rep(1,5),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==13,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==13,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==14,9:38]=matrix(rep(c(rep(1,5),rep(0,5),c(1,0,0,0,0),rep(0,5),rep(0,5),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==14,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==14,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==15,9:38]=matrix(rep(c(c(1,1,0,0,1),c(0,0,0,0,1),c(1,0,0,0,0),c(1,0,0,0,0),rep(0,5),c(1,1,1,0,0)),dim(af_differences.df[af_differences.df$ccaa==15,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==15,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==16,9:38]=matrix(rep(c(c(0,0,0,0,1),c(0,0,0,0,0),c(1,0,0,0,0),c(1,0,0,0,0),rep(0,5),c(1,1,0,0,0)),dim(af_differences.df[af_differences.df$ccaa==16,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==16,])[1],byrow=T)
af_differences.df[af_differences.df$ccaa==17,9:38]=matrix(rep(c(rep(1,5),rep(1,5),c(1,0,0,0,0),c(1,1,0,1,1),c(0,0,1,1,1),rep(1,5)),dim(af_differences.df[af_differences.df$ccaa==17,])[1]),nrow=dim(af_differences.df[af_differences.df$ccaa==17,])[1],byrow=T)

# check
aggregate(af_differences.df[,9:38],by=list(ccaa=af_differences.df$ccaa),mean)

# We add 5 fot include those in 2.
apply(af_differences.df[,9:38],1,sum) + 5

se.p2.prov=(afcityhigh[,14]-afcity[,14])/1.96
se.p1.prov=(afcityhigh[,7]-afcity[,7])/1.96

# assuming independence, Var(x-y)=var(x)+var(y)
af_differences.df$se.diff=se.p2.prov+se.p1.prov

summary(lm(extreme.heat~ p1.1+p1.2+p1.3+p1.4+p1.5+
				p3.1+p3.2+p3.3+p3.4+p3.5+
				p4.1+p4.2+p4.3+p4.4+p4.5+
				p5.1+p5.2+p5.3+p5.4+p5.5+
				p6.1+p6.2+p6.3+p6.4+p6.5+
				p8.1+p8.2+p8.3+p8.4+p8.5,data=af_differences.df))
library(leaps)
dat=af_differences.df[,c(7,9:38)]
best.subset <- regsubsets(extreme.heat~., dat, nvmax=8)
best.subset.summary <- summary(best.subset)
best.subset.summary$outmat

best.subset.by.adjr2 <- which.max(best.subset.summary$adjr2)
best.subset.by.adjr2

coef(best.subset,5)

a=cor(dat)
a[1,]

af_differences.df$index=apply(af_differences.df[,9:38],1,sum) + 5

summary(lm(extreme.heat~index,data=af_differences.df))
summary(lm(heat~index,data=af_differences.df))
summary(lm(glob~index,data=af_differences.df))
summary(lm(cold~index,data=af_differences.df))


plot(af_differences.df$index,af_differences.df$extreme.heat,xlab="Index",ylab="")
abline(h=0,lty=2)
abline(lm(af_differences.df$extreme.heat~af_differences.df$index))


### Add differences by hw_plan
## obtained from code_plan_new.R
load(file="af_cities_plan.RData")

af.diff.plan=afcity[,4]-afcity[,3]
se.p2.plan=(afcityhigh[,4]-afcity[,4])/1.96
se.p1.plan=(afcityhigh[,3]-afcity[,3])/1.96

# assuming independence, Var(x-y)=var(x)+var(y)
se.diff.plan=se.p2.plan+se.p1.plan

df.diff.plan=data.frame(af.diff.plan=af.diff.plan,se.diff.plan=se.diff.plan)
af_differences.df$prov=row.names(af_differences.df)

af.diff.df.all=merge(af_differences.df,df.diff.plan,by="row.names",all.x=T)

cov.prov <- read.csv("t:\\xavi\\etec\\erica\\mortality\\separate models\\attrfraction_metaregression.csv",sep=";")
cov.prov <- read.csv("H:/doctorat/Mortality/02_Stata/04_results/Tables/attrfraction_metaregression.csv",sep=";")

cov.prov=cov.prov[,c(1,8:21)]
cov.prov$prov=as.character(cov.prov$provcode)
 

af.dif.cov=merge(af.diff.df.all,cov.prov,by="prov",all.x=T)

load("T:\\xavi\\etec\\erica\\mortality\\tempDEATHS.Rdata")
load("H:\\doctorat\\Mortality\\02_Stata\\03_data\\tempDEATHS.Rdata")
aa=aggregate(tempDEATHS$tempmax_compl,by=list(province=tempDEATHS$province),mean,na.rm=T)
names(aa)[2]="avg.tempmax"
aa$prov <- c("Alava", "Albacete", "Alicante", "Almeria", "Avila", 
  "Badajoz", "Illes Balears", "Barcelona", "Burgos", "Caceres", "Cadiz", 
  "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", 
  "Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen", "Leon", 
  "Lleida", "La Rioja", "Lugo", "Madrid", "Malaga", "Murcia", "Navarra", 
  "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
  "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", 
  "Tarragona", "Teruel", "Toledo", "Valencia", "Valladolid", "Vizcaya", 
  "Zamora", "Zaragoza")

af.dif.cov=merge(af.dif.cov,aa,by="prov",all.x=T)


# meta-regression
library(metafor)
af.dif.cov$var.diff=af.dif.cov$se.diff^2
af.dif.cov$var.diff.plan=af.dif.cov$se.diff.plan^2

res.extreme.heat <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index, data=af.dif.cov)
summary(res.extreme.heat)

size.extreme.heat=(1/(af.dif.cov$var.diff+res.extreme.heat$tau2))/10

plot(af.dif.cov$index,af.dif.cov$extreme.heat,xlab="Index",
	ylab="AF period 1 - AF Period 2",cex=size.extreme.heat,pch=19,las=1,col="gray50")
### add predicted values (and corresponding CI bounds)
preds <- predict(res.extreme.heat, newmods=c(6:31))
lines(6:31, preds$pred,lwd=2)
lines(6:31, preds$ci.lb, lty="dashed",lwd=2)
lines(6:31, preds$ci.ub, lty="dashed",lwd=2)
### dotted line at coef=0 (no difference between groups)
abline(h=0, lty="dotted")
#text(af.dif.cov$index+.05,af.dif.cov$extreme.heat+.03,
	labels=af.dif.cov$prov,cex=.5)



## add other covariates to meta-regression
############################################

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index+inhabprov_2011+prov_.unemploy_2001+prov_.aircond_2001+PromedioDeind_vulnerability+avg.tempmax, data=af.dif.cov)
summary(res.extreme.heat.cov)

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index, data=af.dif.cov)
summary(res.extreme.heat.cov)
# signif, zval -3.22

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ inhabprov_2011, data=af.dif.cov)
summary(res.extreme.heat.cov)

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ prov_.unemploy_2001, data=af.dif.cov)
summary(res.extreme.heat.cov)
# signif, zval -3.19

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ prov_.aircond_2001, data=af.dif.cov)
summary(res.extreme.heat.cov)
# signif, zval -3.51

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ PromedioDeind_vulnerability, data=af.dif.cov)
summary(res.extreme.heat.cov)
# signif, zval -3.62

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ avg.tempmax, data=af.dif.cov)
summary(res.extreme.heat.cov)
# signif, zval -4.04

# Index adjusted for others:
res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index + PromedioDeind_vulnerability, data=af.dif.cov)
summary(res.extreme.heat.cov)

# Index adjusted for others:
res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index + avg.tempmax, data=af.dif.cov)
summary(res.extreme.heat.cov)

# Index adjusted for others:
res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index + prov_.aircond_2001, data=af.dif.cov)
summary(res.extreme.heat.cov)


# MOre than 2 variables :
res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index + PromedioDeind_vulnerability+avg.tempmax+prov_.aircond_2001, data=af.dif.cov)
summary(res.extreme.heat.cov)

# final model?
res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ index + PromedioDeind_vulnerability+avg.tempmax, data=af.dif.cov)
summary(res.extreme.heat.cov)

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~  PromedioDeind_vulnerability+avg.tempmax, data=af.dif.cov)
summary(res.extreme.heat.cov)

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~  avg.tempmax, data=af.dif.cov)
summary(res.extreme.heat.cov)

cor(af.dif.cov[,c(42,59,61,58)],use="pairwise.complete.obs")

## Model for index:

summary(lm(index~avg.tempmax+PromedioDeind_vulnerability+prov_.aircond_2001,data=af.dif.cov,na.action=na.omit))
summary(lm(index~avg.tempmax+PromedioDeind_vulnerability,data=af.dif.cov,na.action=na.omit))

cor(af.dif.cov$index,af.dif.cov$avg.tempmax)
cor(af.dif.cov$index,af.dif.cov$PromedioDeind_vulnerability,use="complete.obs")



##########################
# Plots for Extreme heat   #
##########################

par(mfrow=c(1,2))

plot(af.dif.cov$index,af.dif.cov$extreme.heat,xlab="Index",
	ylab="AF period 1 - AF Period 2",cex=size.extreme.heat,pch=19,las=1,col="gray50")
### add predicted values (and corresponding CI bounds)
preds <- predict(res.extreme.heat, newmods=c(6:31))
lines(6:31, preds$pred,lwd=2)
lines(6:31, preds$ci.lb, lty="dashed",lwd=2)
lines(6:31, preds$ci.ub, lty="dashed",lwd=2)
### dotted line at coef=0 (no difference between groups)
abline(h=0, lty="dotted")
#text(8,-.8,"Cor = -0.45")
legend(10,-.6,"Cor = -0.45")
legend("topleft","A",text.font=2)

res.extreme.heat.cov <- rma(yi=extreme.heat, vi=var.diff, mods = ~ avg.tempmax, data=af.dif.cov)

plot(af.dif.cov$avg.tempmax,af.dif.cov$extreme.heat,xlab="Average Maximum Temperature",
	ylab="AF period 1 - AF Period 2",cex=size.extreme.heat,pch=19,las=1,col="gray50")
### add predicted values (and corresponding CI bounds)
preds <- predict(res.extreme.heat.cov, newmods=c(6:31))
lines(6:31, preds$pred,lwd=2)
lines(6:31, preds$ci.lb, lty="dashed",lwd=2)
lines(6:31, preds$ci.ub, lty="dashed",lwd=2)
### dotted line at coef=0 (no difference between groups)
abline(h=0, lty="dotted")
#text(18,-.8,"Cor = -0.56")
legend(18,-.6,"Cor = -0.56")
legend("topleft","B",text.font=2)

library(wCorr)
weightedCorr(af.dif.cov$extreme.heat, af.dif.cov$index,weights=(1/(af.dif.cov$var.diff+res.extreme.heat$tau2)),method="pearson")
weightedCorr(af.dif.cov$extreme.heat, af.dif.cov$avg.tempmax,weights=(1/(af.dif.cov$var.diff+res.extreme.heat$tau2)),method="pearson")




##############
### hw_plan
##############
res.hw.plan <- rma(yi=af.diff.plan, vi=var.diff.plan, mods = ~ index, data=af.dif.cov)
summary(res.hw.plan)

size.hw.plan=(1/(af.dif.cov$var.diff.plan+res.hw.plan$tau2))/200

plot(af.dif.cov$index,af.dif.cov$af.diff.plan,xlab="Index",
	ylab="",cex=size.hw.plan,pch=19,las=1,col="gray50")
### add predicted values (and corresponding CI bounds)
preds <- predict(res.hw.plan, newmods=c(6:31))
lines(6:31, preds$pred,lwd=2)
lines(6:31, preds$ci.lb, lty="dashed",lwd=2)
lines(6:31, preds$ci.ub, lty="dashed",lwd=2)
### dotted line at coef=0 (no difference between groups)
abline(h=0, lty="dotted")
#text(af.dif.cov$index+.05,af.dif.cov$af.diff.plan+.03,
	labels=af.dif.cov$prov,cex=.5)



#### hw_plan

res.hw.plan <- rma(yi=af.diff.plan, vi=var.diff.plan, mods = ~ index, data=af.dif.cov)
summary(res.hw.plan)
# signif, -2.2

res.hw.plan.cov <- rma(yi=af.diff.plan, vi=var.diff.plan, mods = ~prov_.aircond_2001 , data=af.dif.cov)
summary(res.hw.plan.cov)
# signif, -2.3

res.hw.plan.cov <- rma(yi=af.diff.plan, vi=var.diff.plan, mods = ~PromedioDeind_vulnerability , data=af.dif.cov)
summary(res.hw.plan.cov)

res.hw.plan.cov <- rma(yi=af.diff.plan, vi=var.diff.plan, mods = ~ avg.tempmax , data=af.dif.cov)
summary(res.hw.plan.cov)
# -1.84

## two variables

res.hw.plan.cov <- rma(yi=af.diff.plan, vi=var.diff.plan, mods = ~ index+avg.tempmax, data=af.dif.cov)
summary(res.hw.plan.cov)

res.hw.plan.cov <- rma(yi=af.diff.plan, vi=var.diff.plan, mods = ~ index+prov_.aircond_2001, data=af.dif.cov)
summary(res.hw.plan.cov)


plot(af_differences.df$index,af_differences.df$glob,xlab="Index",ylab="")
plot(af_differences.df$index,af_differences.df$heat,xlab="Index",ylab="")




