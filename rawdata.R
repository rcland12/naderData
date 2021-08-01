setwd("C:\\Users\\rl02898\\Documents\\Time Series\\Final Project")
library("randtests")
library("forecast")

HannanRisanen <- function(X, p, q){
  coef <- ar(X, method = "yule-walker")$ar
  n <- length(X)
  m <- length(coef)
  Xhat <- rep(0, n)
  for (i in (m+1):n) {
    Xhat[i] <- coef%*%rev(X[(i-m):(i-1)])
  }
  What <- X-Xhat
  Xmat <- matrix(NA, nrow = n-m-q, ncol = p)
  Wmat <- matrix(NA, nrow = n-m-q, ncol = q)
  for (i in 1:(n-m-q)) {
    Xmat[i,] <- X[(m+q+i-1):(m+q+i-p)]
    Wmat[i,] <- What[(m+q+i-1):(m+i)]
  }
  Z <- cbind(Xmat, Wmat)
  beta <- solve(t(Z)%*%Z, t(Z)%*%X[(m+q+1):n])
  return(list("phi" = beta[1:p], "theta" = beta[(p+1):(p+q)]))
}
getpsi <- function(AR, MA, m){
  n <- max(c(m+1,length(AR),length(MA)))
  phi <- c(AR,rep(0,n-length(AR)))
  theta <- c(MA,rep(0,n-length(MA)))
  psi <- rep(0,n+1)
  for(i in 1:(n+1)){
    if(i==1){
      psi[1] <- theta[1] + phi[1]
    }else{
      psi[i] <- theta[i] + sum(c(psi[(i-1):1],1)*phi[1:i])
    }
  }
  psi <- c(1,psi[1:(length(psi)-1)])
  return(psi[1:(m+1)])
}
myARMAacf <- function(AR, MA, m=50, ACVF = F){
  m1 <- max(100,2*m)
  psi <- getpsi(AR, MA, m1)
  acvf <- rep(0,m1)
  for(h in 1:m1){
    acvf[h] <- sum(psi[1:(m1-h+1)]*psi[h:m1])
  }
  if(ACVF){
    return(acvf[1:m])
  }
  return(acvf[1:m]/acvf[1])
}
innovations <- function(AR, MA, n){
  ACVF <- myARMAacf(AR, MA, m=max(50,2*n), ACVF = T)
  #gamma <- ACVF[2:(n+1)]
  Theta <- matrix(0, ncol=n, nrow=n)
  v <- vector(length=n+1)
  v[1] <- ACVF[1]
  Theta[1,1] <- ACVF[2]/v[1]
  v[2] <- ACVF[1]-Theta[1,1]^2*v[1]
  for(i in 2:n){
    for(k in 0:(i-1)){
      if(i-k+1==(n+1)){
        Theta[i,i-k] <- (ACVF[abs(i-k)+1]-sum(Theta[k,k:2]*Theta[i,i:(i-k)]*v[1:(k+1)]))/v[k+1]  
      }else{
        Theta[i,i-k] <- (ACVF[abs(i-k)+1]-sum(Theta[k,k:1]*Theta[i,i:(i-k+1)]*v[1:k]))/v[k+1]
      }
      
    }
    v[i+1] <- ACVF[1]-sum(Theta[i,i:1]^2*v[1:i])
  }
  list(Theta,v)
}
hStepPred <- function(x, H, AR, MA){
  n <- length(x)
  Theta <- innovations(AR, MA, n+H-1)[[1]]
  xhats <- rep(0,n); xhat_n_hs <- rep(0, H)
  for(i in 2:n){
    xhats[i] <- sum(Theta[i-1,1:(i-1)]*rev(x[1:(i-1)]-xhats[1:(i-1)]))
  }
  for(h in 1:H){
    xhat_n_hs[h] <- sum(Theta[n+h-1,h:(n+h-1)]*rev(x-xhats))
  }
  return(xhat_n_hs)
}

tornadoes <- c(
    7,  20,  21,  15,  61,  28,  23,  13,   3,   2,   4,   4,
    2,  10,   6,  26,  57,  76,  23,  27,   9,   2,  12,  10, 
   12,  27,  43,  37,  34,  34,  27,  16,   1,   0,   6,   3,
   14,  16,  40,  47,  94, 110,  32,  24,   5,   6,  12,  21,
    2,  17,  62, 113, 101, 107,  45,  49,  21,  14,   2,  17,
    3,   4,  42,  99, 147, 153,  49,  33,  15,  23,  20,   3,
    2,  47,  31,  85,  79,  65,  92,  42,  16,  29,   7,   9,
   17,   5,  38, 216, 227, 147,  55,  20,  17,  18,  59,  38,
   11,  20,  15,  76,  68, 128, 121,  46,  24,   9,  45,   1,
   16,  20,  43,  30, 226,  73,  63,  38,  58,  24,  11,   2,
    9,  28,  28,  70, 201, 125,  42,  48,  21,  18,  25,   1,
    1,  31, 124,  74, 137, 107,  77,  27,  53,  14,  36,  16,
   12,  25,  37,  41, 200, 171,  78,  51,  24,  11,   5,   2,
   15,   5,  49,  84,  71,  90,  62,  26,  33,  13,  15,   0,
   14,   2,  36, 157, 134, 137,  63,  79,  25,  22,  17,  18,
   21,  32,  34, 123, 273, 147,  85,  61,  64,  16,  34,   7,
    1,  28,  12,  80,  98, 126, 100,  58,  22,  29,  20,  11,
   39,   8,  42, 150, 116, 210,  90,  28, 139,  36,   8,  61,
    5,   7,  28, 102, 142, 136,  56,  66,  25,  14,  44,  32, 
    3,   5,   8,  68, 145, 137,  98,  70,  20,  26,   5,  23,
    9,  16,  25, 117,  88, 134,  81,  55,  54,  50,  10,  14,
   19,  83,  40,  75, 166, 199, 100,  50,  47,  38,  16,  56, 
   33,   7,  69,  96, 140, 114, 115,  59,  49,  34,  17,   8,
   33,  10,  80, 150, 250, 224,  80,  51,  69,  25,  81,  49,
   24,  23,  36, 267, 144, 194,  59, 107,  25,  45,  13,   8,
   52,  45,  84, 108, 188, 196,  79,  60,  34,  12,  39,  22,
   12,  36, 180, 113,  55, 169,  84,  38,  35,  11,   0,   1,
    5,  17,  64,  88, 228, 132,  99,  82,  65,  25,  24,  23,
   23,   7,  17, 107, 213, 148, 143,  65,  20,   7,   9,  30,
   16,   4,  53, 123, 112, 150, 132, 126,  69,  47,  21,   2,
    7,  11,  41, 137, 203, 217,  95,  73,  37,  43,   3,   2,
    2,  25,  33,  84, 187, 223,  98,  64,  26,  32,   7,   1,
   18,   3,  60, 150, 329, 196,  95,  34,  38,   9,  19,  96,
   13,  41,  71,  65, 249, 178,  99,  76,  19,  13,  49,  58,
    1,  27,  73, 176, 169, 242,  72,  47,  17,  49,  30,   4,
    2,   7,  38, 134, 182,  82,  51, 108,  40,  18,  19,   3,
    0,  30,  76,  84, 173, 134,  88,  67,  65,  26,  17,   5,
    6,  19,  38,  20, 126, 132, 163,  63,  19,   1,  55,  14,
   17,   4,  28,  58, 132,  63, 103,  61,  76,  19, 121,  20,
   14,  18,  43,  82, 231, 252,  59,  36,  31,  30,  57,   3,
   11,  57,  86, 108, 243, 329, 106,  60,  45,  35,  18,  35,
   29,  11, 157, 204, 335, 216,  64,  46,  26,  21,  20,   3,
   15,  29,  55,  53, 137, 399, 213, 115,  81,  34, 161,  20,
   17,  78,  48,  85, 177, 313, 242, 112,  65,  55,  19,   6,
   13,   9,  58, 205, 161, 234, 155, 120,  30,  51,  42,   4,
   36,   7,  49, 130, 392, 216, 162,  53,  19,  74,  79,  18,
   35,  14,  71, 177, 235, 128, 202,  72, 101,  68,  55,  15,
   70,  23, 102, 114, 225, 193, 188,  84,  32, 100,  25,  12,
   47,  72,  72, 182, 310, 376,  82,  61, 104,  86,  26,   6,
  216,  22,  56, 177, 310, 289, 102,  79,  56,  17,   7,  15,
   16,  56, 102, 135, 241, 135, 148,  52,  47,  64,  50,  26,
    5,  30,  34, 135, 240, 249, 120,  69,  84, 117, 110,  22,
    3,   2,  47, 117, 204,  97,  68,  86,  61,  58,  96,  99,
    0,  18,  43, 156, 542, 292, 167,  44,  32,  26,  53,   1,
    3,   9,  46, 125, 509, 268, 124, 179, 297,  79, 150,  26,
   33,  10,  62, 132, 122, 317, 138, 123, 133,  18, 149,  26,
   47,  12, 143, 244, 139, 139,  29,  31,  61,  74,  42,  40,
   22,  54, 181, 167, 251, 128,  69,  73,  51,  87,   7,  19,
   84, 147, 129, 189, 461, 294,  94, 101, 111,  21,  15,  46,
    6,  36, 115, 226, 202, 270, 120,  60,   8,  65,   3,  48,
   29,   1,  34, 139, 304, 325, 129,  55,  57, 108,  53,  32,
   16,  68,  75, 771, 322, 156, 102,  59,  51,  23,  46,  15,
   79,  55, 150, 205, 121, 112,  35,  40,  39,  36,   7,  53, 
   75,  39,  19,  80, 227, 126,  69,  46,  22,  61,  77,  14,
    4,  44,  20, 128, 130, 283,  90,  71,  42,  73,  23,  20,
   27,   2,  11, 171, 383, 184, 115,  47,  17,  40,  98,  83,
   18,  99,  84, 142, 218,  87, 108,  90,  40,  21,  50,  19,
  135,  69, 194, 215, 288, 142,  82, 116,  51,  72,  42,  12,
   15,  51,  55, 129, 169, 153,  91,  82, 108, 116,  88,  66,
   22,  27, 111, 276, 509, 179, 100,  78,  83,  62,  19,  56,
   88,  42,  79, 251, 126,  91,  96, 171,  37,  19,  24,  26,
   14,  10, 123
)

write.csv(tornadoes, file = "tornadoesMonths1950.csv", row.names = F)
write.csv(tornadoes[481:855], file = "tornadoesMonths1990.csv", row.names = F)
data <- ts(read.csv(file = "tornadoesMonths1990.csv"), start = 1990, frequency = 12)

# plot from 1950
data2 <- ts(read.csv(file = "tornadoesMonths1950.csv"), start = 1950, frequency = 12)
plot.ts(data2, main = "Monthly Confirmed Tornadoes in the US Since 1950", ylab = "Tornado Count", ylim = c(0, 850))
text(2002.5, 800, labels = "2011 Super Outbreak", cex = 0.5, font = 2)
arrows(2008.5, 800, x1 = 2010.5, y1 = 750, code = 0)
text(1960, 350, labels = "1965 Palm Sunday Outbreak", cex = 0.5, font = 2)
arrows(1962, 325, x1 = 1964.5, y1 = 250, code = 0)
text(1992, 650, labels = "Outbreak Sequence of May 2003 and 2004", cex = 0.5, font = 2)
arrows(1997, 625, x1 = 2002.5, y1 = 500, code = 0)
text(1970, 450, labels = "1974 Super Outbreak", cex = 0.5, font = 2)
arrows(1971, 425, x1 = 1973.95, y1 = 250, code = 0)

# plot from 1990
plot.ts(data, main = "Monthly Confirmed Tornadoes in the US Since 1990", ylim = c(0, 850))
text(2006, 800, labels = "2011 Super Outbreak", cex = 0.5, font = 2)
arrows(2009, 790, x1 = 2010.75, y1 = 750, code = 0)
text(2015, 600, labels = "2019 May Outbreak", cex = 0.5, font = 2)
arrows(2017, 575, x1 = 2018.9, y1 = 500, code = 0)
text(1995, 650, labels = "Outbreak Sequence of May 2003 and 2004", cex = 0.5, font = 2)
arrows(1997, 625, x1 = 2002.5, y1 = 500, code = 0)

# time series analysis
trend <- decompose(data)$trend
season <- decompose(data)$seasonal
noise <- decompose(data)$random
noise <- noise[!is.na(noise)]      
trend <- trend[!is.na(trend)]
plot(decompose(data))

# check for whiteness
Box.test(noise, type="Ljung-Box", lag=20)
turning.point.test(noise)
difference.sign.test(noise)
rank.test(noise)
qqnorm(noise)
qqline(noise)

# create a model
acf(noise, lag.max = 50)
acf(noise, lag.max = 50, type = "partial")

iterations <- 6    # maximum order we want to test for each p and q
AICCs <- matrix(NA, nrow = iterations + 1, ncol = iterations + 1)
bestAICC <- Inf
bestOrder <- c()
for (i in 0:iterations) { # loop over p
  for (j in 0:iterations) { # loop over q
    print(c(i, j))
    AICC <- arima(noise, order = c(i,0,j), include.mean = F, method = "ML")$aic
    if (AICC < bestAICC){
      bestOrder <- c(i,j) # (p,q)
      bestAICC <- AICC
      print(paste("The optimal number of AR and MA coefficients are", bestOrder[1], "and", bestOrder[2], sep = " "))
    }
    AICCs[i+1,j+1] <- AICC
  }
}
model <- HannanRisanen(noise, bestOrder[1], bestOrder[2])
phi_hat <- model$phi
theta_hat <- model$theta

# make predictions
pred2021 <- rep(0, 9)
for (i in 7:15) {
  n <- length(noise)
  xhat_n_hs <- hStepPred(noise, i, phi_hat, theta_hat)
  holt_fit <- holt(trend, h = i) 
  pred2021[i-6] <- holt_fit$mean[i]+season[i-3]+xhat_n_hs[i]
}
finalPredictions <- as.data.frame(pred2021)
colnames(finalPredictions) <- "Predictions"
rownames(finalPredictions) <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
finalPredictions

# plot the predictions
predTS <- ts(finalPredictions, start = c(2021, 4), frequency = 12)
plotTS <- window(data, start=c(2017,1), end = c(2021,3))
plot(NA, NA,
     xlab = "Year",
     ylab = "Tornado Count",
     xlim = c(2017,2022),
     ylim = c(0, max(as.matrix(plotTS), as.matrix(finalPredictions))),
     main = "Tornadoes in the US from 2017 to 2022")
points(plotTS, type = "l")
points(predTS, type = "p", pch = 25, col = "red")

# find the year with the cloest mean squared difference to 2021
year2021 <- c(data[373:375], finalPredictions[[1]])
diff <- rep(0, floor(length(data)/12))
for (i in 1:floor(length(data)/12)) {
  diff[i] <- sum((data[(i*12-11):(i*12)] - year2021)^2)
}
closestYears <- order(diff)

plot.ts(year2021,
        type = "l",
        xlab = "Month",
        ylab = "Tornado Count",
        main = "2021 vs. 2000 and 2007",
        ylim = c(0, 350),
        col = "black")
lines(data[(closestYears[1]*12-11):(closestYears[1]*12)],
      type = "l",
      col = "red")
lines(data[(closestYears[2]*12-11):(closestYears[2]*12)],
      type = "l",
      col = "blue")
legend(10, 240,
       legend=c("2021", "2000", "2007"),
       col=c("black", "red", "blue"),
       lty=1, cex=0.8, text.font=4,
       title="Year", bg='lightblue')

# plot the 2021 data
predTS2 <- ts(finalPredictions, start = c(2021, 4), frequency = 12)
plotTS2 <- window(data, start=c(2021,1), end = c(2021,3))
plot(NA, NA,
     xlab = "Year",
     ylab = "Tornado Count",
     xlim = c(2021,2022),
     ylim = c(0, max(as.matrix(plotTS2), as.matrix(finalPredictions))),
     main = "Tornado Forecasts for 2021")
points(plotTS2, type = "l")
points(predTS2, type = "p", pch = 25, col = "red")