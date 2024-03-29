#INTRODUCTION
# Website link : "https://www.insee.fr/fr/statistiques/serie/010767732#Documentation"

#PACKAGES
install.packages("fUnitRoots") #Unit root test
require(tseries) #functions for time series
library(fUnitRoots)

datafile <- "data.csv"
data <- read.csv(datafile, sep=";")