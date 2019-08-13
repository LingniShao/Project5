library(reshape2)  # data processing
library(ggplot2)  # data visualization
library(GGally)  # data visualization
library(imputeTS)  # impute missing values in univariate time series.
library(caret)  # package for creating cross validation folds
library(forecast)  # time series analysis



##------------------------------------------------------------------------------
# data loading and precossing


# set the working path
setwd("~/Documents/UoE/1_courses/dissertation/subproject05/code/data")

# define a method for changing string of specific forms to Date type.
setClass('gnDate')
setAs("character","gnDate", function(from) as.Date(from, format="%d/%m/%Y"))
setClass('flowDate')
setAs("character","flowDate", function(from) as.Date(from, format="%Y-%m-%d"))
setClass('levelDate')
setAs("character","levelDate", function(from) as.Date(from, format="%Y/%m/%d"))

# load data into workspace with formated date column
gnData <- read.csv('generation.csv', 
                   colClasses=c('gnDate', 'NULL', 'NULL', 'numeric', 'numeric', 
                                'numeric', 'numeric', 'numeric'))
flowDataD <- read.csv('don_parkhill_flow.csv', 
                      colClasses=c('flowDate', 'numeric', 'NULL'))
levelDataD <- read.csv('don_parkhill_level.csv', 
                       colClasses=c('levelDate', 'NULL', 'numeric', 'NULL'))
names(gnData) <- c('date','power','usl','dsl','tspeed','tflow')
names(levelDataD) <- c('date','level')


# get the daily generation data by suming the power and averaging the others
gnDataD <- aggregate(cbind(power,usl,dsl,tspeed,tflow)~date, gnData, mean)
# merge datasets by day
start <- max(gnDataD$date[1], flowDataD$date[1], levelDataD$date[1])
end <- min(gnDataD$date[nrow(gnDataD)], flowDataD$date[nrow(flowDataD)], 
           levelDataD$date[nrow(levelDataD)])
dates <- data.frame(date=seq(from=start, to=end, by='day'))
myDataD <- merge(merge(merge(dates,gnDataD, by='date', all.x=T), 
                       flowDataD, by='date', all.x=T), 
                 levelDataD, by='date', all.x=T)

# plot variates against time
mdata <- melt(data=myDataD,id='date')
ggplot(mdata, aes(x=date, y=value)) + 
  geom_line(aes(color=variable), size=0.4) +
  facet_grid(variable~., scales='free_y') +
  theme(legend.position = "none", axis.title=element_blank()) +
  scale_x_date(date_breaks='1 month', date_labels="%b")


##------------------------------------------------------------------------------
# exploratory data analysis
# delete the abnormal record
gnData <- gnData[-which.min(gnData$tspeed),]
# deal with noise

# get the daily generation data by suming the power and averaging the others
gnDataD <- aggregate(cbind(power,usl,dsl,tspeed,tflow)~date, gnData, mean)
# merge datasets by day
start <- max(gnDataD$date[1], flowDataD$date[1], levelDataD$date[1])
end <- min(gnDataD$date[nrow(gnDataD)], flowDataD$date[nrow(flowDataD)], 
           levelDataD$date[nrow(levelDataD)])
dates <- data.frame(date=seq(from=start, to=end, by='day'))
myDataD <- merge(merge(merge(dates,gnDataD, by='date', all.x=T), 
                       flowDataD, by='date', all.x=T), 
                 levelDataD, by='date', all.x=T)



# impute missing values and treat data with 0 power but flow > 6.76 as missing
myDataD_im <- myDataD
par(mfrow=c(1,1))
for (i in 2:8) {
  if (anyNA(myDataD[,i]) == T){
    myDataD_im[,i] <- na_kalman(myDataD[,i])
    plotNA.imputations(myDataD[,i],myDataD_im[,i], 
                       legend=F, ylab=colnames(myDataD)[i])
  }
}
summary(myDataD_im)
# plot variates against time
mdata <- melt(data=myDataD_im,id='date')
ggplot(mdata, aes(x=date, y=value)) + 
  geom_line(aes(color=variable), size=0.4) +
  facet_grid(variable~., scales='free_y') +
  theme(legend.position = "none", axis.title=element_blank()) +
  scale_x_date(date_breaks='1 month', date_labels="%b")


# scatter plots and correlation 
lowerFn <- function(data, mapping){
  p <- ggplot(data=data, mapping=mapping) +
    geom_point(color='aquamarine4', alpha=0.3, size=0.4) +
    geom_smooth(color='gray40', fill='gray40', method='auto', size=0.4)
  p
}
diagFn <- function(data,mapping){
  p <- ggplot(data,mapping=mapping) + 
    geom_histogram(aes(y=..density..),colour="gray40", fill="white") +
    geom_density(alpha=.3, color='gray40',fill="gold")
  p
}
ggpairs(myDataD_im, columns=c(7,8,3,4,5,6,2),
        lower=list(continuous=lowerFn),
        diag=list(continuous=diagFn),
        upper=list(continuous=wrap('cor', size=3)))

# transform data to achieve linearity
myDataD_im$logflow <- log(myDataD_im$flow)
myDataD_im$loglevel <- log(myDataD_im$level)

##------------------------------------------------------------------------------
# useful functions
getCV <- function(data, k, model='lm', segZ=NULL, npsi=NULL,...){
  # return cross valiadation rmse for linear model and lm/glm model
  # Args
  #   data: a dataframe of response(1st column) and covariates
  #   k: the fold for cross validation
  #   model: the model of interest
  #   segZ: segmented variable if model='seglm'
  #   npsi: the number of breaking points if model='seglm'
  y <- colnames(data)[1]
  rmses <- numeric(k)
  if (NCOL(data)==1){
    folds <- createFolds(data, k=k)
    if (model=='lm'){
      for (i in 1:k){
        tra <- data[-folds[[i]]]
        fit <- lm(tra~1)
        rmses[i] <- mean((data[folds[[i]]]-fit$coefficients)**2,na.rm=T)
      }
    } else if (model=='seglm'){
      print('unsupported calculation')
    }
  } else {
    folds <- createFolds(data[,1], k=k)
    if (model=='lm'){
      for (i in 1:k){
        fit <- lm(as.formula(paste(y, '~.')), data[-folds[[i]],])
        pred <- predict(fit, data[folds[[i]],])
        rmses[i] <- RMSE(pred, data[folds[[i]],1])
      }
    } else if (model=='seglm'){
      for (i in 1:k){
        fit <- segmented.lm(lm(as.formula(paste(y, '~.')), data[-folds[[i]],]),
                            seg.Z=as.formula(paste('~',paste(segZ,collapse='+'))),
                            npsi=npsi)
        pred <- predict(fit, cbind(data[folds[[i]],]))
        rmses[i] <- RMSE(pred, data[folds[[i]],1])
      }
    }
  }
  
  mean(rmses,na.rm=T)
}

farima <- function(y,h,order=NULL,model) {
  # Return a forecast object for models without regressors
  # Args:
  #   y: a univariate time series
  #   h: an integer means forecast lag
  #   model: time series model
  #   order: order of arima if method='arima
  if (model=='arima'){
    forecast(Arima(y, order=order), h=h)
  } else if (model=='auto.arima') {
    forecast(auto.arima(y), h=h) 
  }
}

farimaX <- function(y, h, xreg, model='auto.arima', order=NULL,...){
  # Return a forecast object for models with regressors
  # Args:
  #   y: a univariate time series
  #   h: an integer means forecast lag
  #   xreg: a matrix of regressors
  #   model: time series model
  #   order: order of arima if method='arima
  ncol <- NCOL(xreg)
  X <- matrix(xreg[1:length(y), ], ncol = ncol)
  if(NROW(xreg) < length(y) + h)
    stop("Not enough xreg data for forecasting")
  newX <- matrix(xreg[length(y) + (1:h), ], ncol = ncol)
  if (model=='auto.arima'){
    fit <- auto.arima(y, xreg=X)
    fcst <- forecast(fit, xreg = newX, h=h) 
  } else if (model=='arima'){
    fit <- Arima(y, order=order, xreg=X)
    fcst <- forecast(fit, xreg = newX, h=h) 
  } 
  fcst
}

ChzXreg <- function(data,model='lm',k=NULL, direction='backward', criteria='rmse',...){
  # Choose regressors by using backward/forward elimination based on rmse and return with AIC
  # Args: 
  #   data: a dataframe of reponse and regressors (the first column is response)
  #   model: the model of interest
  #   k: the parameter for getCV function if used
  N <- n <- ncol(data)-1
  can <- data.frame(NA)
  chZcol <- data.frame(NA)
  Rmses <- data.frame(name=rep(NA,N+1),rmse=rep(Inf,N+1), aic=rep(Inf,N+1))
  if (model=='lm'){
    if (direction=='backward'){
      y <- colnames(data)[1]
      Rmses$rmse[1] <- getCV(data, k=k, model='lm')
      Rmses$name[1] <- paste(unlist(colnames(data)[-1]), collapse=' ')
      Rmses$aic[1] <- AIC(lm(as.formula(paste(y,'~.')),data))
      Rmses$rmse[N+1] <- getCV(data[,1], k=k, model='lm')
      Rmses$aic[N+1] <- AIC(lm(as.formula(paste(y,'~1')),data))
      Rmses$name[N+1] <- 'No regressor'
      while (n>1){
        for (i in 1:n){
          Ddata <- data[,-i-1,drop=F]
          if (criteria=='rmse'){
            if (getCV(Ddata, k=k, model='lm') < Rmses$rmse[N+2-n]){
              can <- Ddata
              Rmses$name[N+2-n] <- paste(unlist(colnames(Ddata)[-1]), collapse=' ')
              Rmses$rmse[N+2-n] <- getCV(Ddata, k=k, model='lm')
              Rmses$aic[N+2-n] <- AIC(lm(as.formula(paste(y,'~.')), Ddata))
            }
          } else if (criteria=='aic'){
            if (AIC(lm(as.formula(paste(y,'~.')),Ddata)) < Rmses$aic[N+2-n]){
              can <- Ddata
              Rmses$name[N+2-n] <- paste(unlist(colnames(Ddata)[-1]), collapse=' ')
              Rmses$rmse[N+2-n] <- getCV(Ddata, k=k, model='lm')
              Rmses$aic[N+2-n] <- AIC(lm(as.formula(paste(y,'~.')), Ddata))
            }
          }
          
        }
        data <- can
        n <- ncol(data)-1
      }
    } else if (direction=='forward'){
      y <- colnames(data)[1]
      Rmses$rmse[N+1] <- getCV(data, k=k, model='lm')
      Rmses$name[N+1] <- paste(unlist(colnames(data)[-1]), collapse=' ')
      Rmses$aic[N+1] <- AIC(lm(as.formula(paste(y,'~.')),data))
      Rmses$rmse[1] <- getCV(data[,1], k=k, model='lm')
      Rmses$aic[1] <- AIC(lm(as.formula(paste(y,'~1')),data))
      Rmses$name[1] <- 'No regressor'
      while (n>1){
        for (i in 1:n){
          if (anyNA(chZcol)==T){
            Ddata <- data[,c(1,i+1),drop=F]
          } else {
            Ddata <- cbind(data[,c(1,i+1),drop=F],chZcol)
          }
          if (criteria=='rmse'){
            if (getCV(Ddata, k=k, model='lm') < Rmses$rmse[N+2-n]){
              can <- Ddata[,-1,drop=F]
              Rmses$name[N+2-n] <- paste(unlist(colnames(Ddata)[-1]), collapse=' ')
              Rmses$rmse[N+2-n] <- getCV(Ddata, k=k, model='lm')
              Rmses$aic[N+2-n] <- AIC(lm(as.formula(paste(y,'~.')), Ddata))
            }
          } else if (criteria=='aic'){
            if (AIC(lm(as.formula(paste(y,'~.')),Ddata))< Rmses$aic[N+2-n]){
              can <- Ddata[,-1,drop=F]
              Rmses$name[N+2-n] <- paste(unlist(colnames(Ddata)[-1]), collapse=' ')
              Rmses$rmse[N+2-n] <- getCV(Ddata, k=k, model='lm')
              Rmses$aic[N+2-n] <- AIC(lm(as.formula(paste(y,'~.')), Ddata))
            }
          }
        }
        chZcol <- can
        data <- data[,-which(names(data) %in% names(chZcol)),drop=F]
        n <- NCOL(data)-1
      }
    }
  } else if (model=='auto.arima'){
    if (direction=='backward'){
      y <- ts(data[,1])
      xreg <- data[,-1]
      Rmses$rmse[1] <- sqrt(mean(tsCV(y, farimaX, model=model, h=1,
                                      xreg=data.matrix(xreg))^2, na.rm=T))
      Rmses$aic[1] <- auto.arima(y,xreg=data.matrix(xreg))$aic
      Rmses$name[1] <- paste(unlist(colnames(xreg)), collapse=' ')
      Rmses$rmse[N+1] <- sqrt(mean(tsCV(y,farima,h=1, model=model)^2,
                                   na.rm=T))
      Rmses$aic[N+1] <- auto.arima(y)$aic
      Rmses$name[N+1] <- 'No regressor'
      while (n>1){
        for (i in 1:n){
          x <- xreg[,-i,drop=F]
          if (criteria=='rmse'){
            re <- sqrt(mean(tsCV(y, farimaX, model=model, h=1, 
                                 xreg=data.matrix(x))^2, na.rm=T))
            if ( re < Rmses$rmse[N+2-n]){
              can <- x
              Rmses$name[N+2-n] <- paste(unlist(colnames(x)), collapse=' ')
              Rmses$rmse[N+2-n] <- re
              Rmses$aic[N+2-n] <- auto.arima(y,xreg=data.matrix(x))$aic
            }
          } else if (criteria=='aic'){
            re <- auto.arima(y,xreg=data.matrix(x))$aic
            if ( re < Rmses$aic[N+2-n]){
              can <- x
              Rmses$name[N+2-n] <- paste(unlist(colnames(x)), collapse=' ')
              Rmses$rmse[N+2-n] <- sqrt(mean(tsCV(y, farimaX, model=model, h=1, 
                                                  xreg=data.matrix(x))^2, na.rm=T))
              Rmses$aic[N+2-n] <- re
            }
          }
        }
        xreg <- can
        n <- ncol(xreg)
      }
    } else if (direction=='forward'){
      y <- ts(data[,1])
      xreg <- data[,-1]
      Rmses$rmse[N+1] <- sqrt(mean(tsCV(y, farimaX, model=model, h=1,
                                      xreg=data.matrix(xreg))^2, na.rm=T))
      Rmses$aic[N+1] <- auto.arima(y,xreg=data.matrix(xreg))$aic
      Rmses$name[N+1] <- paste(unlist(colnames(xreg)), collapse=' ')
      Rmses$rmse[1] <- sqrt(mean(tsCV(y,farima,h=1, model=model)^2,
                                   na.rm=T))
      Rmses$aic[1] <- auto.arima(y)$aic
      Rmses$name[1] <- 'No regressor'
      while (n>1){
        for (i in 1:n){
          if (anyNA(chZcol)==T){
            x <- xreg[,i,drop=F]
          } else {
            x <- cbind(xreg[,i,drop=F],chZcol)
          }
          if (criteria=='rmse'){
            re <- sqrt(mean(tsCV(y, farimaX, model=model, h=1, 
                                 xreg=data.matrix(x))^2, na.rm=T))
            if ( re < Rmses$rmse[N+2-n]){
              can <- x
              Rmses$name[N+2-n] <- paste(unlist(colnames(x)), collapse=' ')
              Rmses$rmse[N+2-n] <- re
              Rmses$aic[N+2-n] <- auto.arima(y,xreg=data.matrix(x))$aic
            }
          } else if(criteria=='aic'){
            re <- auto.arima(y,xreg=data.matrix(x))$aic
            if ( re < Rmses$aic[N+2-n]){
              can <- x
              Rmses$name[N+2-n] <- paste(unlist(colnames(x)), collapse=' ')
              Rmses$rmse[N+2-n] <- sqrt(mean(tsCV(y, farimaX, model=model, h=1, 
                                                  xreg=data.matrix(x))^2, na.rm=T))
              Rmses$aic[N+2-n] <- re
            }
          }
        }
        chZcol <- can
        xreg <- xreg[,-which(names(xreg) %in% names(chZcol)),drop=F]
        n <- ncol(xreg)
      }
    }
  }
  Rmses
}

##______________________________________________________________________________
# fit model for power

# linear model
selectLmPowerF1 <- ChzXreg(data=myDataD_im[,c(2:6,9,10)],model='lm',k=10, 
                         direction='forward', criteria='rmse')
selectLmPowerF2 <- ChzXreg(data=myDataD_im[,c(2:6,9,10)],model='lm',k=10, 
                           direction='forward', criteria='aic')
selectLmPowerB1 <- ChzXreg(data=myDataD_im[,c(2:6,9,10)],model='lm',k=10, 
                          direction='backward',criteria='rmse')
selectLmPowerB2 <- ChzXreg(data=myDataD_im[,c(2:6,9,10)],model='lm',k=10, 
                           direction='backward',criteria='aic')
selectLmPowerF1
selectLmPowerF2
selectLmPowerB1
selectLmPowerB2
# fit selected lm model
fitlmPower <- lm(power~usl+dsl+tflow, myDataD_im)
summary(fitlmPower)
# residual plots of selected lm model
par(mfcol = c(1, 3))
plot(fitlmPower$residuals, ylab='Residuals')
plot(fitlmPower,which=c(1,2))
# calculate cross validation rmse
set.seed(1)
rmseLmPower <- getCV(myDataD_im[,c(2:4,6)], k=10, model='lm')


# fit time series model
selectArimaPowerF1 <- ChzXreg(data=myDataD_im[,c(2:8)],model='auto.arima',
                             direction='forward',criteria='rmse')
selectArimaPowerF2 <- ChzXreg(data=myDataD_im[,c(2:8)],model='auto.arima',
                              direction='forward',criteria='aic')
selectArimaPowerB1 <- ChzXreg(data=myDataD_im[,c(2:8)],model='auto.arima',
                             direction='backward',criteria='rmse')
selectArimaPowerB2 <- ChzXreg(data=myDataD_im[,c(2:8)],model='auto.arima',
                              direction='backward',criteria='aic')
selectArimaPowerF1
selectArimaPowerF2
selectArimaPowerB1
selectArimaPowerB2
# transform data type to ts
powerTS <- ts(data=myDataD_im$power)
# fit candidate ARIMA models
fitArimaPower1 <- auto.arima(powerTS, xreg=data.matrix(myDataD_im[,c(3,5,7,8)]))
fitArimaPower2 <- auto.arima(powerTS, xreg=data.matrix(myDataD_im[,c(3,5,8)]))
fitArimaPower3 <- auto.arima(powerTS, xreg=data.matrix(myDataD_im[,c(3,4,6)]))
fitArimaPower4 <- auto.arima(powerTS, xreg=data.matrix(myDataD_im[,c(3,4,5,8)]))
summary(fitArimaPower1)
summary(fitArimaPower2)
summary(fitArimaPower3)
summary(fitArimaPower4)
# calculate time series cross validation rmse and mae for candiate models
set.seed(1)
resArimaPower1 <- tsCV(powerTS, farimaX, model='arima',order=c(2,1,1), h=1,
                       xreg=data.matrix(myDataD_im[,c(3,5,7,8)]))
resArimaPower2 <- tsCV(powerTS, farimaX, model='arima',order=c(2,1,1), h=1,
                       xreg=data.matrix(myDataD_im[,c(3,5,8)]))
resArimaPower3 <- tsCV(powerTS, farimaX, model='arima',order=c(1,0,2), h=1,
                       xreg=data.matrix(myDataD_im[,c(3,4,6)]))
resArimaPower4 <- tsCV(powerTS, farimaX, model='arima',order=c(1,0,2), h=1,
                       xreg=data.matrix(myDataD_im[,c(3,4,5,8)]))

rmseArimaPower1 <- sqrt(mean(resArimaPower1^2, na.rm=T))
maeArimaPower1 <- mean(abs(resArimaPower1), na.rm=T)
rmseArimaPower2 <- sqrt(mean(resArimaPower2^2, na.rm=T))
maeArimaPower2 <- mean(abs(resArimaPower2), na.rm=T)
rmseArimaPower3 <- sqrt(mean(resArimaPower3^2, na.rm=T))
maeArimaPower3 <- mean(abs(resArimaPower3), na.rm=T)
rmseArimaPower4 <- sqrt(mean(resArimaPower4^2, na.rm=T))
maeArimaPower4 <- mean(abs(resArimaPower4), na.rm=T)
rmseArimaPower1; rmseArimaPower2; rmseArimaPower3; rmseArimaPower4
maeArimaPower1; maeArimaPower2; maeArimaPower3; maeArimaPower4
# model diagnostic
tsdiag(fitArimaPower3)

##------------------------------------------------------------------------------
# modeling tflow

# simple linear model
selectlmTflow <- ChzXreg(data=myDataD_im[,c(6,9,10)],model='lm',k=10,
                         direction='forward',criteria='aic')

selectlmTflow
# fit selected lm model
lfitlmTflow <- lm(tflow~loglevel, myDataD_im)
summary(fitlmTflow)
# residual plots of selected lm model
par(mfcol = c(1, 3))
plot(fitlmTflow$residuals, ylab='Residuals')
plot(fitlmTflow,which=c(1,2))



# arima model
selectArimaTflow <- ChzXreg(data=myDataD_im[,c(6,7,8)],model='auto.arima',
                             direction='forward', criteria='aic')

selectArimaTflow
# transform data type to ts
tflowTS <- ts(data=myDataD_im$tflow)
# fit candidate ARIMA models
fitArimaTflow1 <- auto.arima(tflowTS, xreg=data.matrix(myDataD_im[,8]))
fitArimaTflow2 <- auto.arima(tflowTS, xreg=data.matrix(myDataD_im[,c(7,8)]))
summary(fitArimaTflow1)
summary(fitArimaTflow2)
# calculate time series cross validation rmse and mae
set.seed(1)
resArimaTflow1 <- tsCV(tflowTS, farimaX, method='arima', order=c(0,1,2), h=1,
                       xreg=data.matrix(myDataD_im[,8]))
resArimaTflow2 <- tsCV(tflowTS, farimaX, method='arima', order=c(0,1,2), h=1,
                      xreg=data.matrix(myDataD_im[,c(7,8)]))
rmseArimaTflow1 <- sqrt(mean(resArimaTflow1^2, na.rm=T)) 
maeArimaTflow1 <- mean(abs(resArimaTflow1),na.rm=T)
rmseArimaTflow2 <- sqrt(mean(resArimaTflow2^2, na.rm=T)) 
maeArimaTflow2 <- mean(abs(resArimaTflow2),na.rm=T)
rmseArimaTflow1; rmseArimaTflow2
maeArimaTflow1; maeArimaTflow2


# model diagnostic
tsdiag(fitArimaTflow1)



##------------------------------------------------------------------------------
# modeling usl

# simple linear model
selectlmUsl <- ChzXreg(data=myDataD_im[,c(3,7,8)],model='lm',k=10,
                         direction='forward',criteria='aic')

selectlmUsl
# fit selected lm model
fitlmUsl<- lm(usl~level+flow, myDataD_im)
summary(fitlmUsl)
# residual plots of selected lm model
par(mfcol = c(1, 3))
plot(fitlmUsl$residuals, ylab='Residuals')
plot(fitlmUsl,which=c(1,2))



# select arima model
selectArimaUsl <- ChzXreg(data=myDataD_im[,c(3,7,8)],model='auto.arima',
                          direction='forward',criteria='aic')

selectArimaUsl
# transform data type to ts
uslTS <- ts(data=myDataD_im$usl)
# fit candidate ARIMA models
fitArimaUsl1 <- auto.arima(uslTS, xreg=data.matrix(myDataD_im[,8]))
fitArimaUsl2 <- auto.arima(uslTS, xreg=data.matrix(myDataD_im[,c(7,8)]))
summary(fitArimaUsl1)
summary(fitArimaUsl2)
# calculate time series cross validation rmse and mae for candidate models
set.seed(1)
resArimaUsl1 <- tsCV(uslTS, farimaX, method='arima', order=c(1,1,1), h=1,
                     xreg=data.matrix(myDataD_im[,8]))
resArimaUsl2 <- tsCV(uslTS, farimaX, method='arima', order=c(1,1,1), h=1,
                     xreg=data.matrix(myDataD_im[,c(7,8)]))
rmseArimaUsl1 <- sqrt(mean(resArimaUsl1^2, na.rm=T))
maeAArimaUsl1 <- mean(abs(resArimaUsl1), na.rm=T)
rmseArimaUsl2 <- sqrt(mean(resArimaUsl2^2, na.rm=T))
maeAArimaUsl2 <- mean(abs(resArimaUsl2), na.rm=T)
rmseArimaUsl1; rmseArimaUsl2
maeAArimaUsl1; maeAArimaUsl2


# model diagnostic
tsdiag(fitArimaUsl1)

##------------------------------------------------------------------------------
# modeling dsl

# simple linear model
selectLmDsl <- ChzXreg(data=myDataD_im[,c(4,7,8)],model='lm',k=10,
                         direction='forward',criteria='aic')


selectLmDsl
# fit selected lm model
fitlmDsl<- lm(dsl~level+flow, myDataD_im)
summary(fitlmDsl)
# residual plots of selected lm model
par(mfcol = c(1, 3))
plot(fitlmDsl$residuals, ylab='Residuals')
plot(fitlmDsl,which=c(1,2))



# select arima models
selectArimaDsl <- ChzXreg(data=myDataD_im[,c(4,7,8)],model='auto.arima',
                          direction='forward',criteria='aic')

selectArimaDsl
# transform data type to ts
dslTS <- ts(data=myDataD_im$dsl)
# fit candidate arima models
fitArimaDsl1 <- auto.arima(dslTS, xreg=data.matrix(myDataD_im[,8]))
fitArimaDsl2 <- auto.arima(dslTS, xreg=data.matrix(myDataD_im[,c(7,8)]))
summary(fitArimaDsl1)
summary(fitArimaDsl2)
# calculate time series cross validation rmse and mae for candiate models
set.seed(1)
resArimaDsl1 <- tsCV(dslTS, farimaX, method='arima', order=c(1,1,3), h=1,
                     xreg=data.matrix(myDataD_im[,8]))
resArimaDsl2 <- tsCV(dslTS, farimaX, method='arima', order=c(1,1,3), h=1,
                     xreg=data.matrix(myDataD_im[,c(7,8)]))
rmseArimaDsl1 <- sqrt(mean(resArimaDsl1^2, na.rm=T))
rmseArimaDsl2 <- sqrt(mean(resArimaDsl2^2, na.rm=T))
maeArimaDsl1 <- mean(abs(resArimaDsl1), na.rm=T)
maeArimaDsl2 <- mean(abs(resArimaDsl2), na.rm=T)
rmseArimaDsl1; rmseArimaDsl2
maeArimaDsl1; maeArimaDsl2


# model diagnostic
tsdiag(fitArimaDsl1)

##------------------------------------------------------------------------------
# modeling level data
# select data
Dates <- data.frame(date=seq.Date(as.Date('2017-03-31'), as.Date('2019-06-04'), by='day'))
my_level <- merge(Dates, levelDataD, by='date', all.x=T)
# summary data set
summary(my_level)
# impute missing values
my_level_im <- my_level
my_level_im$level<- na_kalman(my_level$level)
summary(my_level_im)
# change the class of level to ts 
levelTS <- ts(data=my_level_im$level)

# fit ARIMA model
fitArimaLevel <- auto.arima(levelTS,approximation = F)
summary(fitArimaLevel)
# model diagonistic
tsdiag(fitArimaLevel)
# calculate time series cross validation rmse and mae
resArimaLevel <- tsCV(levelTS,farima,h=1, order=c(3,1,2), model='arima')
rmseArimaLevel <- sqrt(mean(resArimaLevel^2, na.rm=T))
maeArimaLevel <- mean(abs(resArimaLevel), na.rm=T)
rmseArimaLevel
maeArimaLevel

##______________________________________________________________________________
# forecast

# forecast level
predLevel1 <- forecast(fitArimaLevel,h=27)
# forecast tflow, usl and dsl using observed level
predTflow0 <- forecast(fitArimaTflow1, h=27, xreg=levelDataD$level[792:818])
predUsl0 <- forecast(fitArimaUsl1, h=27, xreg=levelDataD$level[792:818])
predDsl0 <- forecast(fitArimaDsl1, h=27, xreg=levelDataD$level[792:818])
# forecast tflow, usl and dsl using forecast level
predTflow1 <- forecast(fitArimaTflow1, h=27, xreg=as.numeric(predLevel1$mean))
predUsl1 <- forecast(fitArimaUsl1, h=27, xreg=as.numeric(predLevel1$mean))
predDsl1 <- forecast(fitArimaDsl1, h=27, xreg=as.numeric(predLevel1$mean))
# forecast power using observed tflow, usl and dsl
predPower00 <- forecast(fitArimaPower3, h=27, 
                        xreg=data.matrix(gnDataD[326:352,c(3,4,6)]))
# forecast power using forecast tflow, usl and dsl, which are forecast using 
# observed level
predX0 <- data.frame(usl=as.numeric(predUsl0$mean),
                     dsl=as.numeric(predDsl0$mean),
                     tflow=as.numeric(predTflow0$mean))
predPower0 <- forecast(fitArimaPower3, h=27, xreg=data.matrix(predX0))
# forecast power using forecast tflow, usl and dsl, which are forecast using 
# forecast level
predX1 <- data.frame(usl=as.numeric(predUsl1$mean),
                    dsl=as.numeric(predDsl1$mean),
                    tflow=as.numeric(predTflow1$mean))
predPower1 <- forecast(fitArimaPower3, h=27, xreg=data.matrix(predX1))

# plot forecast results
par(mfrow=c(1,1))
autoplot(predLevel1,ylab='level',main=NULL) +
  geom_line(aes(797:825, levelDataD$level[792:820])) +
  geom_line(aes(1:796, as.numeric(fitArimaLevel$fitted)),color='orange')

autoplot(predTflow0,ylab='tflow',main=NULL) +
  geom_line(aes(328:354, gnDataD$tflow[326:352]))+
  geom_line(aes(1:327, as.numeric(fitArimaTflow1$fitted)),color='orange')
autoplot(predTflow1,ylab='tflow',main=NULL) +
  geom_line(aes(328:354, gnDataD$tflow[326:352]))+
  geom_line(aes(1:327, as.numeric(fitArimaTflow1$fitted)),color='orange')

autoplot(predUsl0,ylab='usl',main=NULL)+
  geom_line(aes(328:354, gnDataD$usl[326:352])) +
  geom_line(aes(1:327, as.numeric(fitArimaUsl1$fitted)),color='orange')
autoplot(predUsl1,ylab='usl',main=NULL)+
  geom_line(aes(328:354, gnDataD$usl[326:352])) +
  geom_line(aes(1:327, as.numeric(fitArimaUsl1$fitted)),color='orange')

autoplot(predDsl0,ylab='dsl',main=NULL)+
  geom_line(aes(328:354, gnDataD$dsl[326:352]))+
  geom_line(aes(1:327, as.numeric(fitArimaDsl1$fitted)),color='orange')
autoplot(predDsl1,ylab='dsl',main=NULL)+
  geom_line(aes(328:354, gnDataD$dsl[326:352]))+
  geom_line(aes(1:327, as.numeric(fitArimaDsl1$fitted)),color='orange')

autoplot(predPower0,ylab='power',main=NULL) +
  geom_line(aes(328:354, gnDataD$power[326:352])) +
  geom_line(aes(1:327, as.numeric(fitArimaPower3$fitted)),color='orange')
autoplot(predPower1,ylab='power',main=NULL) +
  geom_line(aes(328:354, gnDataD$power[326:352])) +
  geom_line(aes(1:327, as.numeric(fitArimaPower3$fitted)),color='orange')
autoplot(predPower00,ylab='power',main=NULL) +
  geom_line(aes(328:354, gnDataD$power[326:352])) +
  geom_line(aes(1:327, as.numeric(fitArimaPower3$fitted)),color='orange')


