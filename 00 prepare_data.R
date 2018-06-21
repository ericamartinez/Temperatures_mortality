
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
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel)


load("tempDEATHS.Rdata")

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



