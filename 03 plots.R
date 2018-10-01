
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
# PLOTS
################################################################################

################################################################################


####################################################################
# 1. PLOT WITH TEMPERATURE PERCENTILE ####
indlab <- predper%in%c(0,1,10,50,90,99,100)

#pdf("figure2.pdf",width=4,height=4)
#layout(matrix(1:1,ncol=1))
#par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))

################################################################################
# POOLED: BEFORE AND AFTER PERIODS

# pdf("figure2.pdf",width=4,height=4)
# layout(1)
# par(mex=0.8,mgp=c(2.5,1,0),las=1)

period <- c("Period 1","Period 2")

plot(cp1,ylab="Percent change (%)",xlab="Temperature percentile",axes=F,
     ylim=c(0.8,2.2),lwd=2,col="chocolate1",ci.arg=list(density=20,col="chocolate1"))
lines(cp2,ci="area",lwd=2,col="olivedrab4",ci.arg=list(density=20,angle=-45,col="olivedrab4"))
axis(1,at=tmeancountry[indlab],labels=predper[indlab],cex.axis=0.9)
axis(2,cex.axis=0.9, at=c(0.8,1.00,1.2,1.4,1.6,1.8,2,2.2), 
     labels=c(-2,0,20,40,60,80,100,120))
#abline(v=c(tmeancountry[cenindcountry],tmeancountry[c("1.0%","99.0%")]),
#       lty=c(3,2,2))
abline(v=c(tmeancountry[c("1.0%","99.0%")]),
       lty=c(2,2))
legend("top",period,col=c("chocolate1","olivedrab4"),lwd=2,bg="white",cex=0.8,ncol=2)



#### OVERALL CUMULATIVE EXPOSURE-RESPONSE ASSOCIATIONS ####
## PERCENT CHANGE BY PROVINCES

lab <- c("Period 1","Period 2")

pdf("figureS1.pdf",width=9,height=9)
layout(matrix(seq(3*3),nrow=3,byrow=T))
par(mex=0.8,mgp=c(2.5,1,0),las=1)

for(i in seq(dlist)) {
  bvar <- onebasis(dlist[[i]]$tempmax_compl,fun="bs",degree=2,
                   knots=if(is.null(varper)) NULL else 
                     quantile(dlist[[i]]$tempmax_compl,varper/100,na.rm=T))
   #cen <- quantile(dlist[[i]]$tempmax_compl,cenpercountry/100,na.rm=T)
  cen1 <- quantile(dlist[[i]]$tempmax_compl,cenpercountry1/100,na.rm=T)
  cen2 <- quantile(dlist[[i]]$tempmax_compl,cenpercountry2/100,na.rm=T)
  
  perc <- quantile(dlist[[i]]$tempmax_compl, probs = seq(0, 1, by= 0.01),na.rm=T)
  
  pred1 <- crosspred(bvar,coef=blup1[[i]]$blup,vcov=blup1[[i]]$vcov,
                     model.link="log",cen=cen1)
  pred2 <- crosspred(bvar,coef=blup2[[i]]$blup,vcov=blup2[[i]]$vcov,
                     model.link="log",cen=cen2)
  
  # Maximum and minimum temperature for periods before and after
  max1 <- max(dlist[[i]]$tempmax_compl[dlist[[i]]$yyyy<2004])
  max2 <- max(dlist[[i]]$tempmax_compl[dlist[[i]]$yyyy>2003])
  min1 <- min(dlist[[i]]$tempmax_compl[dlist[[i]]$yyyy<2004])
  min2 <- min(dlist[[i]]$tempmax_compl[dlist[[i]]$yyyy>2003])
  
  # Fixing the maximum temperature to be represented in each curve
  if(max2>max1){
    pred1$predvar <- pred1$predvar[pred1$predvar<=max1]
    l <- length(pred1$predvar)
    pred1$allRRfit <- pred1$allRRfit[1:l]
    pred1$allRRlow <- pred1$allRRlow[1:l]
    pred1$allRRhigh <- pred1$allRRhigh[1:l]
    pred1$matRRfit <- pred1$matRRfit[1:l]
    pred1$matRRlow <- pred1$matRRlow[1:l]
    pred1$matRRhigh <- pred1$matRRhigh[1:l]
    pred1$matfit <- pred1$matfit[1:l]
    pred1$matse <- pred1$matse[1:l]
    pred1$allfit <- pred1$allfit[1:l]
    pred1$allse <- pred1$allse[1:l]
  }else{
    pred2$predvar <- pred2$predvar[pred2$predvar<=max2]
    l <- length(pred2$predvar)
    pred2$allRRfit <- pred2$allRRfit[1:l]
    pred2$allRRlow <- pred2$allRRlow[1:l]
    pred2$allRRhigh <- pred2$allRRhigh[1:l]
    pred2$matRRfit <- pred2$matRRfit[1:l]
    pred2$matRRlow <- pred2$matRRlow[1:l]
    pred2$matRRhigh <- pred2$matRRhigh[1:l]
    pred2$matfit <- pred2$matfit[1:l]
    pred2$matse <- pred2$matse[1:l]
    pred2$allfit <- pred2$allfit[1:l]
    pred2$allse <- pred2$allse[1:l]
  }
  
  plot(pred1,ylab="Percent change (%)",xlab="Temperature (ºC)",ylim=c(0.8,2.2),
       lwd=2,col="chocolate1",ci.arg=list(density=20,col="chocolate1"),yaxt='n',
       xlim = c(min(min1,min2),max(max1,max2)))
  axis(2,cex.axis=0.9, at=c(0.8,1.00,1.2,1.4,1.6,1.8,2,2.2), 
       labels=c(-2,0,20,40,60,80,100,120))
  lines(pred2,ci="area",lwd=2,col="olivedrab4",ci.arg=list(density=20,angle=-45,col="olivedrab4"))
  abline(v=perc[c("1%","99%")],lty=c(2,2))
  title(provincies_n[i])
  legend("top",lab,col=c("chocolate1","olivedrab4"),lwd=2,bg="white",cex=0.6,ncol=2,inset=0.05)
}

dev.off()



####################
##
## 4. MAP ATTRIBUTABLE FRACTION INJURIES - DIFFERENCES####
## (GADM data)
## Source of information: https://www.students.ncl.ac.uk/keith.newman/r/maps-in-r-using-gadm


library("sp")

setwd("H:/Doctorat/Dades/Maps regions/Gadm")
esp1 <- readRDS("ESP_adm2.rds")
esp1@data

# We introduce attributable fraction in each province (considering the same order than file esp1)
# for cold, heat and extreme heat

esp1@data$attrcold <- c(-6.56093434,-3.72560419,0.55486846,-4.1166136,-3.58492257,-0.70215719,3.35355276,
                        -4.33526916,-1.71342968,0.0498406,-0.2951487,-0.2951487,-1.19683438,-3.78436232,
                        -2.88797704,-0.59176809,-1.1958211,-0.10744298,-0.37793514,-3.94735051,-1.44736439,
                        -4.34860077,-1.40088507,-1.29820462,-2.50883155,-0.27938292,-4.55081325,-3.70190423,
                        -3.16823629,-0.48039605,-1.39401519,0,0,0.27659737,-4.02048013,-0.70192891,-1.42718345,
                        -0.45138984,-5.10659012,-3.78755356,-0.59229994,-3.79111945,-0.92413887,-1.85784147,
                        -0.8863206,3.04056067,-0.91340093,-0.91340093,-4.03804949,-4.03804949,-0.09308928,
                        -0.09308928,-0.98629731,0.49303025,0.39129435,2.72511244)


esp1@data$attrextremeheat <- c(-0.837077121,-0.567156438,-0.722973552,-0.224491874,-0.771560848,0.0545003,
                               -0.296582918,-0.405751468,-0.134677901,-0.018512831,0.143055198,0.143055198,
                               0.273607091,-0.516904137,-0.521742812,0.157080743,-0.232183756,-0.344623208,
                               0.006460784,0.156707813,-0.068455094,-0.025713113,-0.017190827,-0.012974686,
                               0.251527034,0.092947703,0.019771064,-0.124807281,-0.372590637,0.390734552,
                               -0.46272249,0,0,-0.227527437,0.488071013,-0.065497022,-0.148538733,0.081323171,
                               -0.227297942,-0.114877034,0.135475212,0.097462896,0.194545165,-0.305291073,
                               -0.101994548,0.234018109,-0.093383355,-0.093383355,0.325233864,0.325233864,
                               -0.011564877,-0.011564877,0.086369925,0.449350525,0.185842743,-0.143495722)

##################
## COLD
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrcold,c(-50,-1,0,1,50),
              labels=c('<-1', '(-1 - 0]', '(0-1]', '>1'))
#funcol <- colorRampPalette(c("#fef0d9", "#fdcc8a","#fc8d59","#e34a33", "#b30000"))
#funcol <- colorRampPalette(c("#feedde", "#fdbe85","#fd8d3c","#e6550d", "#a63603"))
funcol <- colorRampPalette(c("lightskyblue1","steelblue1", "steelblue4","lightblue4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
pdf("map_cold.pdf",width=6,height=5)

plot(esp1, col = col)
#text(coordinates(esp1)[-c(47,11,49,51),], labels = esp1[-c(47,11,49,51),]$NAME_2, cex=0.5)

legend(-12,39,paste0(c('<-1%', '(0% - -1%]', '[0% - 1%]', '>1%')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Difference in mortality attributable fraction \n in period 1 and 2 (cold)")

dev.off()


##################
## EXTREME HEAT
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrextremeheat,c(-50,-1,0,1,50),
              labels=c('<-1', '(-1 - 0]', '(0-1]', '>1'))
funcol <- colorRampPalette(c("khaki2","darkgoldenrod1", "chocolate","OrangeRed4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
pdf("map_extremeheat.pdf",width=6,height=5)
plot(esp1, col = col)

legend(-12,39,paste0(c('<-1%', '(0% - -1%]', '[0% - 1%]', '>1%')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Difference in mortality attributable fraction \n in period 1 and 2 (extreme heat)")

dev.off()



####################
##
## 5. MAP ATTRIBUTABLE FRACTION INJURIES - PERIOD 1####
## (GADM data)
## Source of information: https://www.students.ncl.ac.uk/keith.newman/r/maps-in-r-using-gadm


library("sp")

setwd("H:/Doctorat/Dades/Maps regions/Gadm")
esp1 <- readRDS("ESP_adm2.rds")
esp1@data

# We introduce attributable fraction in each province (considering the same order than file esp1)
# for cold, heat and extreme heat

esp1@data$attrcold <- c(7.5581425,7.4345704,4.4916013,6.9141991,5.6656479,1.4508683,4.4621408,6.4261121,
                        2.6422422,0.2617288,0.5389545,0.5389545,4.2142233,6.3475092,4.1621277,0.605068,
                        3.1088044,0.3924613,3.5353099,4.9096164,2.3214512,4.9460536,2.7692501,1.5626595,
                        3.5244719,2.632562,4.920208,5.8997933,5.2835209,0.9316378,4.6943723,0,0,2.9985866,
                        5.4546646,4.2850748,4.9648087,4.3928261,6.3112545,4.8881485,3.894882,6.2498022,5.8984811,
                        5.9491604,6.9164597,2.9812551,2.1611447,2.1611447,5.2581002,5.2581002,2.1463628,2.1463628,
                        2.4741082,3.445252,2.7175955,1.9551065)

esp1@data$attrextremeheat <- c(1.2871742,1.0770219,1.1126666,0.835088,1.4480063,0.6715872,1.0768447,1.059989,0.8164274,
                               0.9247279,0.6661909,0.6661909,0.1667613,1.0696716,1.299485,0.8917779,0.9020625,1.0808069,
                               0.678984,0.6135884,0.4791714,0.5433747,0.7785746,1.173891,0.7917428,0.5997346,0.8061137,
                               0.5416311,0.6397949,0.3133908,0.9420012,0,0,0.7534796,0.5009678,0.5062852,0.4395032,0.3055062,
                               0.8360581,1.0185481,0.2743854,0.4164556,0.1580746,1.0090594,0.7003616,0.4287722,0.8783432,
                               0.8783432,0.3527799,0.3527799,0.5970908,0.5970908,0.4817129,0.2730631,0.1035584,0.6014073)

##################
## COLD
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrcold,c(-50,2,4,6,10),
              labels=c('[0-2]', '(2-4]', '(4-6]', '>6'))
#funcol <- colorRampPalette(c("#fef0d9", "#fdcc8a","#fc8d59","#e34a33", "#b30000"))
#funcol <- colorRampPalette(c("#feedde", "#fdbe85","#fd8d3c","#e6550d", "#a63603"))
funcol <- colorRampPalette(c("lightskyblue1","steelblue1", "steelblue4","lightblue4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
pdf("map_cold.pdf",width=6,height=5)

plot(esp1, col = col)
#text(coordinates(esp1)[-c(47,11,49,51),], labels = esp1[-c(47,11,49,51),]$NAME_2, cex=0.5)

legend(-12,39,paste0(c('[0-2]', '(2-4]', '(4-6]', '>6')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Mortality attributable fraction \nPeriod 1 (cold)")

dev.off()



##################
## EXTREME HEAT
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrextremeheat,c(-50,0.4,0.8,1.2,10),
              labels=c('[0-0.4]', '(0.4-0.8]', '(0.8-1.2]', '>1.2'))
funcol <- colorRampPalette(c("khaki2","darkgoldenrod1", "chocolate","OrangeRed4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
pdf("map_extremeheat.pdf",width=6,height=5)
plot(esp1, col = col)

legend(-12,39,paste0(c('[0-0.4]', '(0.4-0.8]', '(0.8-1.2]', '>1.2')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Mortality attributable fraction \nPeriod 1 (extreme heat)")

dev.off()



####################
##
## 6. MAP ATTRIBUTABLE FRACTION INJURIES - PERIOD 2####
## (GADM data)
## Source of information: https://www.students.ncl.ac.uk/keith.newman/r/maps-in-r-using-gadm


library("sp")

setwd("H:/Doctorat/Dades/Maps regions/Gadm")
esp1 <- readRDS("ESP_adm2.rds")
esp1@data

# We introduce attributable fraction in each province (considering the same order than file esp1)
# for cold, heat and extreme heat

esp1@data$attrcold <- c(0.99720818,3.70896617,5.04646978,2.79758551,2.08072532,0.74871109,7.81569351,
                        2.09084296,0.92881252,0.3115694,0.24380582,0.24380582,3.01738897,2.5631469,
                        1.27415069,0.01329988,1.91298332,0.28501837,3.15737473,0.9622659,0.87408683,
                        0.59745281,1.36836504,0.26445493,1.01564035,2.35317904,0.36939473,2.19788911,
                        2.11528465,0.45124179,3.30035711,0,0,3.27518397,1.43418451,3.58314591,3.53762523,
                        3.94143631,1.20466439,1.10059497,3.30258205,2.45868272,4.97434222,4.09131897,
                        6.03013912,6.02181573,1.24774376,1.24774376,1.22005071,1.22005071,
                        2.05327356,2.05327356,1.48781092,3.93828225,3.10888986,4.68021898)

esp1@data$attrextremeheat <- c(0.4500971,0.5098654,0.3896931,0.6105961,0.6764455,0.7260875,0.7802618,0.6542375,
                               0.6817495,0.9062151,0.8092461,0.8092461,0.4403684,0.5527675,0.7777422,1.0488587,
                               0.6698787,0.7361837,0.6854448,0.7702962,0.4107163,0.5176616,0.7613838,1.1609163,
                               1.0432699,0.6926823,0.8258847,0.4168238,0.2672043,0.7041253,0.4792787,0,0,0.5259521,
                               0.9890388,0.4407882,0.2909644,0.3868294,0.6087602,0.903671,0.4098606,0.5139185,0.3526197,
                               0.7037683,0.5983671,0.6627903,0.7849599,0.7849599,0.6780137,0.6780137,0.5855259,0.5855259,
                               0.5680828,0.7224136,0.2894012,0.4579116)

##################
## COLD
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrcold,c(-50,2,4,6,10),
              labels=c('[0-2]', '(2-4]', '(4-6]', '>6'))
#funcol <- colorRampPalette(c("#fef0d9", "#fdcc8a","#fc8d59","#e34a33", "#b30000"))
#funcol <- colorRampPalette(c("#feedde", "#fdbe85","#fd8d3c","#e6550d", "#a63603"))
funcol <- colorRampPalette(c("lightskyblue1","steelblue1", "steelblue4","lightblue4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
pdf("map_cold.pdf",width=6,height=5)

plot(esp1, col = col)
#text(coordinates(esp1)[-c(47,11,49,51),], labels = esp1[-c(47,11,49,51),]$NAME_2, cex=0.5)

legend(-12,39,paste0(c('[0-2]', '(2-4]', '(4-6]', '>6')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Mortality attributable fraction \nPeriod 2 (cold)")

dev.off()



##################
## EXTREME HEAT
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrextremeheat,c(-50,0.4,0.8,1.2,10),
              labels=c('[0-0.4]', '(0.4-0.8]', '(0.8-1.2]', '>1.2'))
funcol <- colorRampPalette(c("khaki2","darkgoldenrod1", "chocolate","OrangeRed4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
pdf("map_extremeheat.pdf",width=6,height=5)
plot(esp1, col = col)

legend(-12,39,paste0(c('[0-0.4]', '(0.4-0.8]', '(0.8-1.2]', '>1.2')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Mortality attributable fraction \nPeriod 2 (extreme heat)")	

dev.off()



####################
##
## 7. MAP ATTRIBUTABLE FRACTION PLAN - PERIOD 1, 2 AND DIFFERENCES####
## (GADM data)
## Source of information: https://www.students.ncl.ac.uk/keith.newman/r/maps-in-r-using-gadm


library("sp")

setwd("H:/Doctorat/Dades/Maps regions/Gadm")
esp1 <- readRDS("ESP_adm2.rds")
esp1@data

# We introduce attributable fraction in each province (considering the same order than file esp1)
# for cold, heat and extreme heat

esp1@data$attrp1 <- c(1.66282125,0.8699196,0.69329923,0.09674808,1.21694473,0.42353253,2.6770967,
                      0.90856849,1.08456882,0,1.14642782,1.14642782,0,1.00824733,4.54338468,1.92089211,
                      0,0.95636013,0,0,0.03264295,0,0.22050719,1.11984533,0,0.04246696,0.03681686,0.92563672,
                      0.68849445,0.42661279,0.21028427,0,0,0.79413386,0,0.13243498,0.29329793,0.14311083,
                      0.83418198,1.0738213,0.07937538,0,0,0,0.03892552,1.31489518,2.08144197,2.08144197,0.17610938,
                      0.17610938,0,0,0,0,0,1.41491454)

esp1@data$attrp2 <- c(0.062177197,0.692847636,0.112820555,-0.001941201,0.505009391,-0.009826026,0.400189787,
                      -0.003883328,2.067932244,0,0.419351084,0.419351084,0,1.050039135,2.519196081,1.687552304,
                      0,0.607973697,0,0,0.00636277,0,0.171069102,0.675207992,0,0.025105283,0.028654764,1.313872529,
                      0.611678555,0.270566644,0.13723457,0,0,0.837806004,0,-0.048039324,0.039869455,0.13648286,
                      1.183733084,0.60399181,-0.011583042,0,0,0,-0.189555432,1.963345408,1.347397745,1.347397745,
                      0.041506744,0.041506744,0,0,0,0,0,0.094187859)

esp1@data$attrdiff <- c(-1.600644053,-0.177071964,-0.580478675,-0.098689281,-0.711935339,-0.433358556,-2.276906913,
                        -0.912451818,0.983363424,0,-0.727076736,-0.727076736,0,0.041791805,-2.024188599,-0.233339806,
                        0,-0.348386433,0,0,-0.02628018,0,-0.049438088,-0.444637338,0,-0.017361677,-0.008162096,
                        0.388235809,-0.076815895,-0.156046146,-0.0730497,0,0,0.043672144,0,-0.180474304,-0.253428475,
                        -0.00662797,0.349551104,-0.46982949,-0.090958422,0,0,0,-0.228480952,0.648450228,-0.734044225,
                        -0.734044225,-0.134602636,-0.134602636,0,0,0,0,0,-1.320726681)


##################
## PERIOD 1
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrp1,c(-50,0,1,2,10),
              labels=c('No activated','[0-1]', '(1-2]', '>2'))
funcol <- colorRampPalette(c("gray76","khaki2", "chocolate","OrangeRed4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
plot(esp1, col = col)

legend(-12,39,paste0(c('Plan no activated','[0-1]', '(1-2]', '>2')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Mortality attributable fraction \nPeriod 1")	



##################
## PERIOD 2
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrp2,c(-50,0,1,2,10),
              labels=c('No activated','[0-1]', '(1-2]', '>2'))
funcol <- colorRampPalette(c("gray76","khaki2", "chocolate","OrangeRed4"))
col <- funcol(length(levels(levels)))[levels]

# Plot map
plot(esp1, col = col)

legend(-12,39,paste0(c('Plan no activated','[0-1]', '(1-2]', '>2')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Mortality attributable fraction \nPeriod 2")	


##################
## DIFFERENCE
# Put levels to the attributable fraction
levels <- cut(esp1@data$attrdiff,c(-50,0,0.5,10),
              labels=c('<0', '[0-0.5]', '>0.5'))
funcol <- colorRampPalette(c("khaki2", "chocolate","OrangeRed4"))
col <- funcol(length(levels(levels)))[levels]
col[c(10,13,17,19,20,22,25,32,33,35,42,43,44,51,52,53,54,55)] <- "gray76"
  

# Plot map
plot(esp1, col = col)
legend(-12,39,paste0(c('Plan no activated','<0', '[0-0.5]', '>0.5')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=c("gray76","#EEE685","#D2691E","#8B2500"),bty="n",cex=0.6,
       title="Mortality attributable fraction \nDifference P2-P1")	


####################
##
## 8. MAP NUMBER OF DAYS THAT THE PLAN WAS ACTIVATED####
## (GADM data)
## Source of information: https://www.students.ncl.ac.uk/keith.newman/r/maps-in-r-using-gadm

library(maps)
map.text("state", regions=c("california", "texas"), labels=as.character(c(125, 250)))

library("sp")

setwd("H:/Doctorat/Dades/Maps regions/Gadm")
esp1 <- readRDS("ESP_adm2.rds")

esp1@data$nump1 <- c(116,26,58,10,52,36,132,73,145,0,48,48,0,73,232,51,0,128,0,0,5,0,5,49,0,10,5,65,10,5,15,
           0,0,28,0,16,21,46,40,57,10,0,0,0,25,20,33,33,20,20,0,0,0,0,0,46)
esp1@data$nump2 <- c(158,118,68,5,86,37,169,119,197,0,90,90,0,121,240,81,17,120,0,5,5,5,15,94,0,5,5,223,73,40,
           25,0,0,65,10,15,56,35,86,96,5,0,0,11,31,74,110,110,6,6,5,5,0,5,0,145)

# Period 1
levels <- cut(esp1@data$nump1,c(-50,10,50,300),
              labels=c('<10', '[10-50]', '>50'))
funcol <- colorRampPalette(c("chartreuse3", "darkorange","red2"))
col <- funcol(length(levels(levels)))[levels]

plot(esp1, col=col)
text(coordinates(esp1), labels = nump1, cex=0.3)
legend(-12,39,paste0(c('<10', '[10-50]', '>50')),
       xjust=0.5,yjust=0.5,pt.cex=1.5,
       pch=22,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.6,
       title="Number of days")	


# Period 2
levels <- cut(esp1@data$nump2,c(-50,10,50,300),
              labels=c('<10', '[10-50]', '>50'))
funcol <- colorRampPalette(c("chartreuse3", "darkorange","red2"))
col <- funcol(length(levels(levels)))[levels]

plot(esp1, col=col)
text(coordinates(esp1), labels = nump1, cex=0.3)
