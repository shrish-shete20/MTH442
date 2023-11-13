library(tseries)    # For adfuller 
library(stats)      # For acf, pacf
library(forecast)   # For auto.arima
library(ggplot2)    # For plotting
library(data.table) # For data manipulation
library(zoo)        # For time series manipulation
library(gridExtra)  # For subplot
library(lubridate)  # For date manipulation
library(DescTools) # for plotting the ACF and PACF

dat = read.csv("Data/MTH 442 data.csv")
df = data.frame(dat)
for(i in 1:1524){
  df[i+2, 1] = dat[3 + (i-1)*5,1]
  for(j in 2: 12){
    val = 0
    dv = 0
    for(k in 1: 5){
      if(dat[2 + (i-1)*5 + k,j] != ""){
        dv = dv +1;
        val = val + as.numeric(dat[2 + (i-1)*5 + k,j])
      }
    }
    df[i+2, j] = val/dv
  }
}

df = df[1:1522, ]
colnames(df)[1]=df[2,]$SDRs.per.currency.unit.for.the.period..January.01..1994...November.10..2023
colnames(df)[2]="Canadian Dollars(CAD)"
colnames(df)[3]=df[2,]$X.1
colnames(df)[4]=df[2,]$X.2
colnames(df)[5]=df[2,]$X.3
colnames(df)[6]=df[2,]$X.4
colnames(df)[7]=df[2,]$X.5
colnames(df)[8]=df[2,]$X.6
colnames(df)[9]=df[2,]$X.7
colnames(df)[10]=df[2,]$X.8
colnames(df)[11]=df[2,]$X.9
colnames(df)[12]=df[2,]$X.10
df=df[-1,]
df=df[-1,]



## Visualizations
date_vector<-df$Date
date_final<-as.Date(date_vector,format = "%d-%b-%Y")
new_df<-data.frame(dates=date_final,numbers=df$`Indian rupee   (INR)                     `)
labels <- c("1994", "2003", "2011","2020")
par(mfrow=c(2,2))
p1<-plot.ts(df$`Indian rupee   (INR)                     `,xaxt="n",ylab = "Indian Rupee(INR)")
axis(1,at=seq(0,length(df$`Indian rupee   (INR)                     `),500),labels=labels)
p2<-plot.ts(df$`U.S. dollar   (USD)                     `,xaxt="n",ylab = "U.S Dollars(USD)")
axis(1,at=seq(0,length(df$`Indian rupee   (INR)                     `),500),labels=labels)
p3<-plot.ts(df$`Japanese yen   (JPY)                     `,xaxt="n",ylab = "Japanese Yen (JPY)")
axis(1,at=seq(0,length(df$`Indian rupee   (INR)                     `),500),labels=labels)
p2<-plot.ts(df$`U.K. pound   (GBP)                     `,xaxt="n",ylab = "Pound(GBP)")
axis(1,at=seq(0,length(df$`Indian rupee   (INR)                     `),500),labels=labels)

### Convert all character columns to numeric
for (i in 2:12){
  df[,i] <- as.numeric(df[,i])
}
for(i in 2:12)
{
  df[,i]<-na.approx(df[,i])
}


### Plots of the trend, seasonal and random components
dats <- as.data.frame(cbind(df$Date, df$`Indian rupee   (INR)                     `))
colnames(dats) <- c("Date", "Indian_Rupee")

library(fpp)
time_se <- dats$Indian_Rupee
time_se <- as.numeric(time_se)
time_se <- na.omit(time_se)
plot(as.ts(time_se))

ts_tot <- ts(time_se, frequency = 52)
decomposition <- decompose(ts_tot, "additive")

plot(as.ts(decomposition$trend), ylab = "Trend")
plot(as.ts(decomposition$seasonal), ylab = "Seasonality")
plot(as.ts(decomposition$random), ylab = "Random")
plot(decomposition)

##### Tests for checking the stationarity

library(tseries)
adf.test(time_se)
kpss.test(time_se)
pp.test(time_se)
plot.ts(time_se)

##### PLotting the time series and ACF and PACF together

library(DescTools)
PlotACF(time_se, main = "INDIAN RUPEE")

####  Time Series for USD

time_series_USD <- df$`U.S. dollar   (USD)                     `
time_series_USD <- as.numeric(time_series_USD)
time_series_USD <- na.omit(time_series_USD)
plot(as.ts(time_series_USD), main = "USD Time Series")
adf.test(time_series_USD)
PlotACF(time_series_USD, main = "USD")    #### non stationary ARMA(2, 0) or ARMA(1, 1, 0)

######### However USD is clearly non stationary by the ADF test
##### Applying first order differencing and then applying the ADF test again

diff_1 <- diff(time_se, lag = 1)
adf.test(diff_1)
plot(as.ts(diff_1), main = "First Order Differencing", ylab = "USD Time Series")
#### By the plot, the first order differenced time series seems to be stationary


##### Plotting the time series, ACF and PACF
PlotACF(diff_1, main = "First Order Differenced USD")   #### stationary AR(1) process


### Validation actual
train<-df[1:1250,]
test<-df[1251:1520,]
model=auto.arima(train$`U.S. dollar   (USD)                     `)
fc<-forecast(model,h=length(test$`Indian rupee   (INR)                     `))
summary(fc)
plot(fc)


### Validation 
validation<-function(p,q,d){
  
  model_arima <- Arima(train$`U.S. dollar   (USD)                     `, order=c(p, d, q))
  
  fc_arima <- forecast(model_arima, h=length(test$`Indian rupee   (INR)                     `))
  
  plot(fc_arima)
  return(summary(fc_arima))
}


### AIC
p <- c(0, 1, 2, 3, 4)
d <- 1
q <- c(0, 1, 2, 3)

min_aic <- Inf
optimal_params <- NULL

for (i in p) {
  for (j in q) {
    tryCatch({
      fit <- arima(train$`U.S. dollar   (USD)                     `, order=c(i, d, j))
      current_aic <- AIC(fit)
      if (current_aic < min_aic) {
        min_aic <- current_aic
        optimal_params <- c(i, d, j)
      }
      print(paste("ARIMA(", i, ",", d, ",", j, ") - AIC:", current_aic, sep=""))
    }, error=function(e){})
  }
}

print(paste("Best params: (", optimal_params[1], ",", optimal_params[2], ",", optimal_params[3], "), min_aic: ", min_aic, sep=""))


### Residuals


predict_values<-c()

valid<-validation(optimal_params[1],optimal_params[3],optimal_params[2])

phi1<-valid$model$coef[1]
x_1250<-df$`U.S. dollar   (USD)                     `[1252]
x_tm<-x_1250
for(i in 1251:1520)
{
  predict_values<-c(predict_values,phi1*x_tm+rnorm(1))
  x_tm<-predict_values[i-1250]
}
residual<--predict_values+test$`U.S. dollar   (USD)                     `

adf.test(residual)
plot(residual,type = 'l')

qqnorm(residual)

