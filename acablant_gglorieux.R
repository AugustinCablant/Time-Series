#INTRODUCTION
# Website link : "https://www.insee.fr/fr/statistiques/serie/010767832#Telechargement"

#PACKAGES
#install.packages("fUnitRoots") #Unit root test
#install.packages("zoo")
#install.packages("forecast")
require(tseries) #functions for time series
library(fUnitRoots)
require(forecast)
library(zoo)
library(astsa) # to properly represent confidence intervals of ARIMA

#DATASET DOWNLOAD 
setwd("/Users/augustincablant/Desktop/Time Series project") #set your path
datafile <- "data.csv"
data <- read.csv(datafile, sep=";")
data <- data[4:nrow(data), ]


class(data$Indice)
data$Indice <- as.numeric(data$Indice) 
dates <- as.yearmon(seq(from=2012+0/12, to=2024+2/12, by=1/12)) #from 2012 to 2024

#production
xm.source <- zoo(data$Indice, order.by=dates)
T <- length(data[[1]])
dates_r <- as.yearmon(seq(from=2012+0/12, to=2024+2/12, by=1/12)) # defines the time bounds of the truncated series
xm <- xm.source

# ---------------------------------------------------------------
#PART 1
#QUESTION 1


#QUESTION 2
#A) PRODUCTION STUDY 
#plot of the serie
plot(xm.source, type = "l", xlab = "Date", ylab = "Production", 
     main = "Figure 1 – Production index for the manufacture \n of essential oils from 2012 to 2024")

#Just to have a look 
#acf
acf(xm.source)  #not usefull

#pacf
pacf(xm.source)


#stationnary test
pp.test(xm.source)  

# Augmented Dickey Fuller test
# Some verifications
summary(lm(xm ~ dates_r))

# Function to perform Ljung-Box tests on the "series" series, with a lag of k, and fitdf = p+q (chosen by the user)
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value 
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

# Perform a series of ADF tests with 24 lags (which corresponds to 2 years) on the xm series. We take the constant + time trend type, consistent with our assumptions above.
series <- xm
kmax <- 24
adftype="ct"

adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2] 
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1
      cat("OK \n") 
    } else 
      cat("nope \n")
    k <- k+1
  }
  return(adf)
}

adf <- adfTest_valid(xm,24,adftype="ct")
adf


#Differenciation or order 1
dxm <- diff(xm,1)
dates_dxm <- dates_r[-1]
summary(lm(dxm ~ dates_dxm))

#ADF
adf <- adfTest_valid(dxm,24,"nc")
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients)) 
adf

#PP
pp.test(dxm)

#QUESTION 3
#Plot the differences between the two series : 
plot(cbind(xm,dxm), main = "Figure 2 - Serie Xm and serie \n xm differencied (order 1)")


# ---------------------------------------------------------------
#PART 2 ARMA(p,q)
#QUESTION 1

#autocorrelation function 
acf(dxm, main = "Figure 3A - Autocorrelation function \n of the differentiated series", 
    lag.max = 20)

pacf(dxm, main = "Figure 3B - Partial autocorrelation \n function of the differentiated series", 
     lag.max = 20)

# Question 2 : We create the ARIMA model and test the autocorrelations of the residuals. 

pmax=4; qmax=2
## Function to test the individual significances of coefficients
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef)) 
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2 
  return(rbind(coef,se,pval))
}
## Function to estimate an ARIMA model and check its fit (significance of coefficients) and validity (absence of residual autocorrelation)
modelchoice <- function(p,q,data=dxm, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") 
    return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA) )
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05 
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals, 24, length(estim$coef) - 1)[, 2] <= 0.05, na.rm = TRUE) == 0 
  checks <- c(arsignif, masignif, resnocorr)
  ok <- as.numeric(sum(checks, na.rm = TRUE) == (3 - sum(is.na(checks)))) 
  return(c("p" = p, "q" = q, "arsignif" = arsignif, "masignif" = masignif, "resnocorr" = resnocorr, "ok" = ok))
}
## Function to estimate and check all ARIMA(p,q) models with p<=pmax and q<=qmax
armamodelchoice <- function(pmax, qmax){
  pqs <- expand.grid(0:pmax, 0:qmax) 
  t(apply(matrix(1:dim(pqs)[1]), 1, function(row) {
    p <- pqs[row, 1]
    q <- pqs[row, 2] 
    cat(paste0("Computing ARMA(", p, ",", q, ") \n")) 
    modelchoice(p, q)
  })) 
}
armamodels <- armamodelchoice(pmax, qmax) # estimates all ARIMA models
selec <- armamodels[armamodels[, "ok"] == 1 & !is.na(armamodels[, "ok"]), ] # well-fitted and valid models
selec

# We want to choose the best model. To do so, we minimize the AIC and BIC criteria.
pqs <- apply(selec, 1, function(row) list("p" = as.numeric(row[1]), "q" = as.numeric(row[2]))) # creates a list of orders p and q for candidate models
names(pqs) <- paste0("arima(", selec[, 1], ",", selec[, 2], ")") # renames the elements of the list
models <- lapply(pqs, function(pq) arima(dxm, c(pq[["p"]], 0, pq[["q"]]))) # creates a list of estimated candidate models
vapply(models, FUN.VALUE = numeric(2), function(m) c("AIC" = AIC(m), "BIC" = BIC(m))) # computes the AIC and BIC of candidate models

#Call the ARIMA 
arima <- arima(dxm,c(0,0,1))
arima

# Part III : predictions 

#our modele 
modele=arima(xm, c(0,1,1))
modele

frequency(xm) <- 12
xm_ts <- as.ts(xm)
prevision <- sarima.for(xm_ts,n.ahead=2,p=0,d=1,q=1,P=0,D=0,Q=0,S=-1,tol = sqrt(.Machine$double.eps),
           no.constant = FALSE, plot = TRUE, plot.all = FALSE,
           xreg = NULL, newxreg = NULL, fixed = NULL)
