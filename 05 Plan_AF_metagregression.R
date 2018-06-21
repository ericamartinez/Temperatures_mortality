
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

attr_prov <- read.csv("attrfraction_metaregression.csv",sep=";")
plan_prov <- read.csv2("ccaa_prov.csv",sep=";")


attr_prov$provcode2 <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14",
                         "15","16","17","18","19","20","21","22","23","24","25","26","27","28",
                         "29","30","31","32","33","34","35","36","37","38","39","40",
                         "41","42","43","44","45","46","47","48","49","50")

plan_prov$provcode2 <- as.character(ifelse(plan_prov$provcode<10,paste0("0",plan_prov$provcode),plan_prov$provcode))
plan_prov$ccaacode2 <- as.character(ifelse(plan_prov$ccaacode<10,paste0("0",plan_prov$ccaacode),plan_prov$ccaacode))

# CREATE THE INDICATOR OF THE ACCOMPLISHMENT OF THE PLAN
indicator <- matrix(NA,nrow=17,ncol=2)
colnames(indicator) <- c("ccaa","indic")
indicator[,1] <- c("01","02","03","04","05","06","07","08","09","10","11","12",
                   "13","14","15","16","17")
indicator[,2] <- c(30,6,13,15,11,30,11,27,30,31,24,30,30,16,14,10,28)

rownames(indicator) <- indicator[,1]
all <- cbind(plan_prov, indicator[plan_prov$ccaacode2, ])

attr_prov_all <- merge(all,attr_prov,by=c("provcode2"))
attr_prov_all <- attr_prov_all[,c(1,5,7,3,8,10:29)]
attr_prov_all$indic <- as.numeric(attr_prov_all$indic)



###
# Differences in AF (P2 - P1)
attr_prov_all$diff_pc <- attr_prov_all$P2_extreme.heat - attr_prov_all$P1_extreme.heat

plot(attr_prov_all$P1_extreme.heat-attr_prov_all$P1_extreme.heat_low,attr_prov_all$P1_extreme.heat_high-attr_prov_all$P1_extreme.heat)
abline(a=0,b=1)
attr_prov_all$se.P1=(attr_prov_all$P1_extreme.heat_high-attr_prov_all$P1_extreme.heat_low)/(1.96*2)
attr_prov_all$se.P2=(attr_prov_all$P2_extreme.heat_high-attr_prov_all$P2_extreme.heat_low)/(1.96*2)

# assuming independence, Var(x-y)=var(x)+var(y)
attr_prov_all$se.diff=(attr_prov_all$se.P1+attr_prov_all$se.P2)


###
# META-ANALYSIS
library(metafor)

res <- rma(diff_pc, se.diff, mods = ~ prov_.unemploy_2011 + prov_.home1per_2011, data=attr_prov_all)
res


