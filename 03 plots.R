
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




