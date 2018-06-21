
################################################################################
# "Temporal changes in temperature-related mortality in Spain and effect of 
#        the implementation of a Heat Health Prevention Plan"
#   
#   ISGlobal  
#   June 2018
#   
#
################################################################################

dd <- read.csv("results_plots.csv",sep = ";")


library(ggplot2)
library(gridExtra)


### ALL VARIABLES TOGETHER (SEX, AGE, CAUSE, URBAN/RURAL, VULNERABILITY INDEX)
dd_all <- dd[c(251:252,1:8,12:15,25:32,39:40,151:154,155:160,188:190,233:235),]

lab <- c("", "","Men","Women", "Men", "Women","16-64","65-74","75-84",">=85","16-64",
         "65-74","75-84",">=85","Diabetes","Diabetes","Mental diseases","Mental diseases",
         "Cardiovascular", "Cardiovascular","Respiratory", "Respiratory","External causes","External causes",
         "Rural","Urban","Rural","Urban","Low","Medium","High","Low","Medium","High",
         "Low","Medium","High","Low","Medium","High")
dd_all$category <- c("Overall","Overall","Sex", "Sex","Sex", "Sex","Age", "Age", "Age", "Age","Age", 
                   "Age", "Age", "Age","Causes of death","Causes of death","Causes of death",
                   "Causes of death","Causes of death","Causes of death","Causes of death",
                   "Causes of death","Causes of death","Causes of death","Environment","Environment",
                   "Environment","Environment","Vulnerability Index","Vulnerability Index","Vulnerability Index",
                   "Vulnerability Index","Vulnerability Index","Vulnerability Index","Air conditioning",
                   "Air conditioning","Air conditioning","Air conditioning","Air conditioning","Air conditioning")
dd_all$category <- factor(dd_all$category, levels = c("Overall","Sex", "Age","Causes of death","Environment",
                                                    "Vulnerability Index","Air conditioning"))
dd_all$lab <- c("", "","Men","Women", "Men", "Women","16-64","65-74","75-84",">=85","16-64",
                "65-74","75-84",">=85","Diabetes","Diabetes","Mental diseases","Mental diseases",
                "Cardiovascular", "Cardiovascular","Respiratory", "Respiratory","External causes","External causes",
                "Rural","Urban","Rural","Urban","Low","Medium","High","Low","Medium","High",
                "Low","Medium","High","Low","Medium","High")
dd_all$lab <- factor(dd_all$lab, levels = lab)
Periode_lab <- c("Period 1","Period 2")
dd_all$Periode_lab <- factor(dd_all$Periode, levels = Periode_lab)



# HEAT EFFECTS
h1 <- ggplot(dd_all, aes(x = lab, y = RR.99th, color=Periode_lab)) + ylab("Heat effects \n(Percent change and 95% CI)")
h1 <- h1 +  geom_point(size = 2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymax = IClow...99th, ymin = ICup...99th),  width=0.1,
                position=position_dodge(width=0.4))
h1 <- h1 + scale_color_manual("Periode", breaks=c("Period 1","Period 2"),values=c("chocolate1","olivedrab4"))
h1 <- h1 + geom_hline(aes(yintercept=1), linetype="dashed", colour="gray60")
h1 <- h1 + scale_y_log10(breaks=c(0.9,1,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80), 
                         labels=c(-5,0,10,20,30,40,50,60,70,80),limits = c(0.95,1.88))
h1 <- h1 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.title.y = element_text(size=11, face="bold"),
                 axis.title.x = element_blank(),
                 axis.text.y=element_text(size=10),
                 legend.title = element_blank(),
                 legend.position = "top")
h1 <- h1 + facet_grid(.~category, switch = NULL, scales = "free_x", space="free_x") +
  theme(strip.text.x = element_text(face="bold", size=9),legend.key = element_blank())
h1 <- h1 + theme(panel.background = element_rect(fill = "gray94"))

dodge <- 0.1
ann_line <- data.frame(lab=c(1,2),x=c(1-dodge,2-dodge),y=c(1.19,1.78),
                       xend=c(1-dodge,2-dodge),yend=c(1.84,1.84),
                       category=dd_all$category[c(15,17)])
h1 <- h1 + geom_segment(data=ann_line, aes(x=x,y=y,xend=xend,yend=yend), color="chocolate1",
                      arrow = arrow(length = unit(0.2, "cm")))

ann_line_r <- data.frame(lab=c(1),x=c(1-dodge),y=c(0.95),
                         xend=c(1-dodge),yend=c(1.19),
                         category=dd_all$category[c(15)])
h1 <- h1 + geom_segment(data=ann_line_r, aes(x=x,y=y,xend=xend,yend=yend), color="chocolate1",
                      arrow = arrow(ends="first",length = unit(0.2, "cm")))

ann_line_a <- data.frame(lab=c(2),x=c(2-dodge),y=c(1.58),
                         xend=c(2-dodge),yend=c(1.78),
                         category=dd_all$category[c(17)]) 
h1 <- h1 + geom_segment(data=ann_line_a, aes(x=x,y=y,xend=xend,yend=yend), color="chocolate1")
hline <- data.frame(tt=c(2-dodge), v=c(1.58),category=dd_all$category[c(17)])
h1 <- h1 + geom_point(data=hline, aes(tt,v), shape=95, size=3, color="chocolate1")
h1

# Annotations
annotation_text <- data.frame(lab=c(2.5), RR.99th=c(1.88), 
                              label=c("89 (71,109)"),
                              category = dd_all$category[c(18)]) # Cause of death
h1 <- h1 + geom_text(aes(x = lab, y = RR.99th, color=NULL), data = annotation_text, label=annotation_text$label, size=2.7)

annotation_point <- data.frame(lab=c(2+dodge), RR.99th=c(1.86),
                               category = dd_all$category[c(18)]) # Cause of death
h1 <- h1 + geom_point(aes(x = lab, y = RR.99th, color=NULL),data = annotation_point, shape= 17, size=2, colour="olivedrab4")
h1


# COLD EFFECTS
c1 <- ggplot(dd_all, aes(x = lab, y = RR.1st, color=Periode_lab)) + ylab("Cold effects \n(Percent change and 95% CI)")
c1 <- c1 +  geom_point(size = 2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymax = IClow...1st, ymin = ICup...1st),  width=0.1,
                position=position_dodge(width=0.4))
c1 <- c1 + scale_color_manual("Periode", breaks=c("Period 1","Period 2"),values=c("chocolate1","olivedrab4"))
c1 <- c1 + geom_hline(aes(yintercept=1), linetype="dashed", colour="gray60")
c1 <- c1 + scale_y_log10(breaks=c(0.9,1,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80), 
                         labels=c(-5,0,10,20,30,40,50,60,70,80),limits = c(0.95,1.86))
c1 <- c1 + theme(axis.text.x = element_text(angle=45, hjust=1, size=10),
                 axis.ticks.x=element_blank(),
                 axis.title.y = element_text(size=11, face="bold"),
                 axis.title.x = element_blank(),
                 axis.text.y=element_text(size=10),
                 legend.position = "none")
c1 <- c1 + facet_grid(.~category, switch = "x", scales = "free_x", space="free_x") +
  theme(strip.text.x = element_text(face="bold", size=9),legend.key = element_blank())
c1 <- c1 + theme(panel.background = element_rect(fill = "gray94"))
c1

plot_all <- grid.arrange(h1,c1)
ggsave(filename="Figure3.png", width = 44, height = 22, units = "cm", dpi=300, plot=plot_all)



