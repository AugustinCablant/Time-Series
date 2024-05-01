#INTRODUCTION
# Website link : "https://www.insee.fr/fr/statistiques/serie/010767732#Documentation"

#PACKAGES
install.packages("fUnitRoots") #Unit root test
require(tseries) #functions for time series
library(fUnitRoots)

#DATASET DOWNLOAD (we use our Github repository)
# getwd()
setwd("/Users/augustincablant/Documents/GitHub/Time-Series")
datafile <- "data.csv"
data <- read.csv(datafile, sep=";")

d <- length(data$Date)  #length of vector Date
dates <- as.yearmon(seq(from=2024+2/12, to=1990, by=-1/12))  #vector from february 2024 to january 1990
#indice
class(data$Indice)

#production
production.source <- zoo(data$Indice, order.by=dates)
T = length(data$Indice)
production <- production.source[1:(T-4)]
N <- length(production)
dates_rev <- rev(dates[5:T])