boot_sd_ci <- function(x, confidence = 95, itr = 1000) { 
  
  # init iterations and sample 
  i <- itr 
  bt_avg <- NULL 
  
  # loop over iterations 
  while (i > 0) { 
    # sample randomly from x 
    spl <- sample(x, length(x), replace = TRUE) 
    
    # store mean 
    bt_avg <- c(bt_avg, sum(spl)/length(spl)) 
    
    # decrement i 
    i <- i-1 
  } 
  
  # mean over the bootstrapped samples 
  bt_est <- sum(bt_avg)/itr 
  
  # compute standard deviation 
  bt_sd <- sqrt((1/(length(x)-1)) * sum((bt_avg-bt_est)^2)) 
  
  # compute confidence interval 
  # sorting bt_avg numerically 
  st_avg <- sort(bt_avg) 
  
  # get the first value after the 2.5 first centiles 
  bt_95ci_inf <- st_avg[floor(0.5*(1-0.01*confidence)*itr)+1] 
  
  # get the last value before the 2.5 last centiles 
  bt_95ci_sup <- st_avg[floor((0.01*confidence+0.5*(1-0.01*confidence))*itr)-1] 
  
  res <- c(bt_sd, bt_95ci_inf, bt_95ci_sup) 
  return(res)  
  
} 

######## intra specific competition ######## 

# import data set 
# merged <- read.csv("mergedResults-intraSpecComp-predSpcf0-predOpnt0-pryConvRateRatio1-pryCtchProbRatio1-pryOfspRatio1-pryMaxConsRatio1.csv", header = T, stringsAsFactors = T, row.names = F) 
merged <- read.csv("../intraSpecComp/folder-intraSpecComp-prey2-ConvRate-over1/intraSpecComp-predSpcf0-predOpnt0-pryConvRateRatio2.5/statsintraSpecComp/mergedResults.csv", header = T, stringsAsFactors = T, row.names = NULL) 

# define tburn 
tburn = 400 

# nb of rep 
repNb <- max(merged$replicate) 

# split plot window in four 
dev.new() 
par(mfrow = c(2,2)) 

# no burning all data at once 
dataset <- merged 

# transform data to be able to fit Nt to deltaN 
deltaNt <- NULL 
Nt <- NULL 

for (k in 1:repNb) { 
  subrep <- subset(dataset, dataset$replicate == as.factor(k)) 
  
  Nt <- c(Nt, subrep$prey2PopulationSize[-length(subrep$prey2PopulationSize)]) 
  deltaNt <- c(deltaNt, subrep$prey2growthRate[-1]) 
} 

# fit model 
mod <- lm(formula = deltaNt/Nt ~ Nt) 
summary(mod) 

print(paste("LM : no burning all replicates at once, r2 = ", mod$coefficients[1], " a22 = ", -mod$coefficients[2]/mod$coefficients[1])) 
# plot  
plot(y = deltaNt/Nt, x = Nt, main = "no burning, all replicates at once") 
abline(a = mod$coefficients[1], b = mod$coefficients[2], col = "red")

# all data at once with burning 
subburn <- subset(merged, merged$timeStep >= tburn) 

dataset <- subburn 

# transform data to be able to fit Nt to deltaN 
deltaNt <- NULL 
Nt <- NULL 

for (k in 1:repNb) { 
  subrep <- subset(dataset, dataset$replicate == as.factor(k)) 
  
  Nt <- c(Nt, subrep$prey2PopulationSize[-length(subrep$prey2PopulationSize)]) 
  deltaNt <- c(deltaNt, subrep$prey2growthRate[-1]) 
} 

# fit model 
mod <- lm(formula = deltaNt/Nt ~ Nt) 
summary(mod) 

print(paste("LM: ", tburn, "burning all replicates at once, r2 = ", mod$coefficients[1], " a22 = ", -mod$coefficients[2]/mod$coefficients[1])) 
print(paste("1/K: ", tburn, "burning all replicates at once, ", "a22 = ", 1/mean(Nt)))  

alpha22 = -mod$coefficients[2]/mod$coefficients[1]

# plot  
plot(y = deltaNt/Nt, x = Nt, main = paste(tburn, "burning, all replicates at once")) 
abline(a = mod$coefficients[1], b = mod$coefficients[2], col = "red") 

# just 1 replicate 
subrep <- subset(merged, merged$replicate == as.factor(1)) 

dataset <- subrep 

Nt <- dataset$prey2PopulationSize[-length(dataset$prey2PopulationSize)] 
deltaNt <- dataset$prey2growthRate[-1]  

# fit model 
mod <- lm(formula = deltaNt/Nt ~ Nt) 
# summary(mod) 

print(paste("no burning, 1 replicate, r2 = ", mod$coefficients[1], " a22 = ", -mod$coefficients[2]/mod$coefficients[1])) 

# plot  
plot(y = deltaNt/Nt, x = Nt, main = "no burning, 1 replicate") 
abline(a = mod$coefficients[1], b = mod$coefficients[2], col = "red") 

# 1 replicate with burning 
dataset <- subset(subburn, subburn$replicate == as.factor(1)) 

Nt <- dataset$prey2PopulationSize[-length(dataset$prey2PopulationSize)] 
deltaNt <- dataset$prey2growthRate[-1]  

# fit model 
mod <- lm(formula = deltaNt/Nt ~ Nt) 
# summary(mod) 

print(paste(tburn, "burning 1 replicate, r2 = ", mod$coefficients[1], " a22 = ", -mod$coefficients[2]/mod$coefficients[1])) 
print(paste("1/K: ", tburn, "burning 1 replicate, ", "a22 = ", 1/mean(Nt))) 

# plot  
plot(y = deltaNt/Nt, x = Nt, main = paste(tburn, "burning, 1 replicate")) 
abline(a = mod$coefficients[1], b = mod$coefficients[2], col = "red") 

# dev.off() 

#### check which r2 is the most accurate #### 
par(mfrow=c(1,1)) 
plot(x = merged$timeStep, y = merged$prey2PopulationSize) 
plot(x = merged$timeStep, y = log(merged$prey2PopulationSize)) 

tInit = 100 

subinit <- subset(merged, merged$timeStep < tInit) 

plot(x = subinit$timeStep, y = subinit$prey2PopulationSize) 
plot(x = subinit$timeStep, y = log(subinit$prey2PopulationSize)) 

# fit model 
mod <- lm(formula = log(subinit$prey2PopulationSize) ~ subinit$timeStep) 
summary(mod) 

abline(a = log(subinit$prey2PopulationSize[1]), b = mod$coefficients[2], col = "red") 
abline(a = mod$coefficients[1], b = mod$coefficients[2], col = "blue") 

#### plot LV and LVc equations to check estimates #### 

plotFit <- function(dataN, dataT, tIntro, tMax, r, alpha, K, N0) { 
  
  # generate a time vector 
  time <- seq(0, tMax) 
  
  # generate N2 according to LV 
  
  # vector for N 
  N.LV = rep(NA, length(time)) 
  N.LV[1] = N0 
  dN = 0 
  
  for (i in 2:(length(time)-1))  
  { 
    
    # before predator intro 
    if (i < tIntro) { 
      dN = r * N.LV[i-1] * (1 - N.LV[i-1]/K) 
      N.LV[i] = N.LV[i-1] + dN 
    } 
    
    # plot(x = dataT, y = dataN, main = paste("r =", round(r, 3), "; alpha22 =", round(alpha22, 3)), xlab = "Time", ylab = "Population size") 
    # points(x = time, y = N.LV, type = "l", col = "blue") 
    
    # # according to LVc 
    # N.LVc = rep(NA, length(time)) 
    # N.LVc[1] = N0 
    # dN = 0 
    
    # for (i in 2:(length(time)-1)) { 
    else { 
      
      # dN = r * N.LVc[i-1] * (1 - alpha*N.LVc[i-1]) 
      # N.LVc[i] = N.LVc[i-1] + dN 
      dN = r * N.LV[i-1] * (1 - alpha*N.LV[i-1]) 
      N.LV[i] = N.LV[i-1] + dN 
      
    } 
  } 
  
  plot(x = dataT/10, y = dataN, main = paste("r =", round(r, 3), "; alpha22 =", round(alpha22, 4)), xlab = "Generations", ylab = "Population size") 
  points(x = time/10, y = N.LV, type = "l", col = "red") 
  
} 

# parameters 
tIntro = 200 
r2 = mod$coefficients[2] 
# r2 = 0.04 
# alpha22 = 0.0027 

subplat <- subset(merged, merged$timeStep > tInit & merged$timeStep <= tIntro) 
K2 = mean(subplat$prey2PopulationSize) 

plotFit(dataT = merged$timeStep, dataN = merged$prey2PopulationSize, tIntro = tIntro, tMax = max(merged$timeStep), r = r2, alpha = alpha22, K = K2, N0 = merged$prey2PopulationSize[1]) 

######## inter specific competition ######## 

# import data set 
bothPreys <- read.csv(file = "../localSA/folder-localSA-ConvRate-over1/localSA-predSpcf0-predOpnt0-pryConvRateRatio1-pryCtchProbRatio1-pryOfspRatio1-pryMaxConsRatio1/statslocalSA/mergedResults.csv") 


#### first method #### 

r1 = r2 
alpha11 = 0.002 
# alpha22 = 0.003 

# import stats 

repNb <- max(bothPreys$replicate) 
tburn = 200 

# no burning all data at once 
# dataset <- bothPreys 

# all data at once with burning 
subburn <- subset(bothPreys, bothPreys$timeStep >= tburn) 

dataset <- subburn 

# transform data to be able to fit Nt to deltaN 
t <- NULL 
deltaN1 <- NULL 
deltaN2 <- NULL 
deltaPr <- NULL
N1 <- NULL 
N2 <- NULL 
Pr <- NULL

for (k in 1:repNb) { 
  subrep <- subset(dataset, dataset$replicate == as.factor(k)) 
  
  t <- c(t, subrep$timeStep[-length(subrep$timeStep)]) 
  N1 <- c(N1, subrep$prey1PopulationSize[-length(subrep$prey1PopulationSize)]) 
  N2 <- c(N2, subrep$prey2PopulationSize[-length(subrep$prey2PopulationSize)]) 
  Pr <- c(Pr, subrep$predator1PopulationSize[-length(subrep$predator1PopulationSize)]) 
  deltaN1 <- c(deltaN1, subrep$prey1growthRate[-1]) 
  deltaN2 <- c(deltaN2, subrep$prey2growthRate[-1]) 
  deltaPr <- c(deltaPr, subrep$predator1growthRate[-1]) 
} 

alpha21est <- 1/N1 * (1 - alpha22*N2 - deltaN2 / (r2 * N2))  
alpha21 = mean(alpha21est) 
alpha21infCI = boot_sd_ci(alpha21est, itr = 10000)[2] 
alpha21supCI = boot_sd_ci(alpha21est, itr = 10000)[3] 

alpha12est <- 1/N2 * (1 - alpha11*N1 - deltaN1 / (r1 * N1))  
alpha12 = mean(alpha12est) 
alpha12infCI = boot_sd_ci(alpha12est)[2] 
alpha12supCI = boot_sd_ci(alpha12est)[3] 

shapiro.test(alpha21est) 
hist(alpha21est) 
hist(alpha21est, xlim = c(-0.005,0.005), breaks = 10) 
hist(alpha12est) 

plot(x = t, y = alpha21est, xlab = "Time", ylab = "Inter-specific competition coefficient", main = "alpha21 estimation") 
rect(xleft = t[1]-100, ybottom = alpha21-sd(alpha21est), xright = max(t)+100, ytop = alpha21+sd(alpha21est), col = rgb(1, 0.8, 0, 0.5), border = "coral") 
rect(xleft = t[1]-100, ybottom = alpha21infCI, xright = max(t)+100, ytop = alpha21supCI, col = rgb(1, 0.2, 0, 0.5), border = "coral") 
abline(h = 0, lty = 2, col = "blue") 
abline(h = alpha21, col = "red", lwd = 2) 

#### try and divide by alpha predators number ####

test12est <- alpha12est/Pr
test12 = mean(test12est) 
test12infCI = boot_sd_ci(test12est, itr = 1000)[2] 
test12supCI = boot_sd_ci(test12est, itr = 1000)[3] 

test21est <- alpha21est/Pr
test21 = mean(test21est) 
test21infCI = boot_sd_ci(test21est, itr = 1000)[2] 
test21supCI = boot_sd_ci(test21est, itr = 1000)[3] 

plot(x = t, y = test21est, xlab = "Time", ylab = "Inter-specific competition coefficient / Pred density", main = "alpha21 estimation") 
rect(xleft = t[1]-100, ybottom = test21-sd(test21est), xright = max(t)+100, ytop = test21+sd(test21est), col = rgb(1, 0.8, 0, 0.5), border = "coral") 
rect(xleft = t[1]-100, ybottom = test21infCI, xright = max(t)+100, ytop = test21supCI, col = rgb(1, 0.2, 0, 0.5), border = "coral") 
abline(h = 0, lty = 2, col = "blue") 
abline(h = test21, col = "red", lwd = 2) 

#### second method #### 

deltaN2solo <- deltaNt #[1:(length(deltaNt)/2)] 
deltaN2pres <- deltaN2 

alpha21est2nd <- (deltaN2solo - deltaN2pres) / (r2 * N1 * N2)  
alpha21 = mean(alpha21est2nd) 
alpha21infCI = boot_sd_ci(alpha21est2nd, itr = 10000)[2] 
alpha21supCI = boot_sd_ci(alpha21est2nd, itr = 10000)[3] 

shapiro.test(alpha21est2nd) 
hist(alpha21est2nd) # , xlim = c(-0.005,0.005), breaks = 10) 

plot(x = t, y = alpha21est2nd, xlab = "Time", ylab = "Inter-specific competition coefficient", main = "alpha21 estimation - 2nd method") 
rect(xleft = t[1]-100, ybottom = alpha21-sd(alpha21est2nd), xright = max(t)+100, ytop = alpha21+sd(alpha21est2nd), col = rgb(1, 0.8, 0, 0.5), border = "coral") 
rect(xleft = t[1]-100, ybottom = alpha21infCI, xright = max(t)+100, ytop = alpha21supCI, col = rgb(1, 0.2, 0, 0.5), border = "coral") 
abline(h = alpha21, col = "red", lwd = 2) 
abline(h = 0, lty = 2, col = "blue") 


#### plot LV and LVc equations to check estimates #### 

plotFit2preys <- function(dataN1, dataN2, dataP, dataT, tIntro, tMax, r1, r2, alpha11, alpha22, alpha21, alpha12, N10, N20) { 
  
  # generate a time vector 
  time <- seq(0, tMax) 
  
  # generate N2 according to LV 
  
  # vector for Ns 
  N1.LVc = rep(NA, length(time)) 
  N1.LVc[1] = N10 
  dN1 = 0 
  
  N2.LVc = rep(NA, length(time)) 
  N2.LVc[1] = N20 
  dN2 = 0 
  
  for (i in 2:length(time)) { 
    
    # before predator intro 
    if (i < tIntro) { 
      
      # compute delta Ns  
      dN1 = r1 * N1.LVc[i-1] * (1 - alpha11 * N1.LVc[i-1]) 
      dN2 = r2 * N2.LVc[i-1] * (1 - alpha11 * N2.LVc[i-1]) 
      
      # calculate Ns(t) 
      N1.LVc[i] = N1.LVc[i-1] + dN1 
      N2.LVc[i] = N2.LVc[i-1] + dN2 
      
    } else { 
      
      # compute delta Ns  
      dN1 = r1 * N1.LVc[i-1] * (1 - alpha11 * N1.LVc[i-1] - alpha12 * N2.LVc[i-1]) 
      dN2 = r2 * N2.LVc[i-1] * (1 - alpha22 * N2.LVc[i-1] - alpha21 * N1.LVc[i-1]) 
      
      # calculate Ns(t) 
      N1.LVc[i] = N1.LVc[i-1] + dN1 
      N2.LVc[i] = N2.LVc[i-1] + dN2 
    } 
  } 
  
  plot(x = 1, y = 1, main = "IBM vs LVc", xlab = "Generations", ylab = "Population size", type = "n",  
       xlim = c(dataT[1]/10, tMax/10), 
       ylim = c(min(c(min(dataN1), min(dataN2), min(dataP)))-10, max(c(max(dataN1), max(dataN2), max(dataP)))+10)) 
  points(x = dataT/10, y = dataN1, col = "red") 
  points(x = dataT/10, y = dataN2, col = "blue") 
  points(x = dataT/10, y = dataP, col = "orange") 
  points(x = time/10, y = N1.LVc, type = "l", col = "darkred") 
  points(x = time/10, y = N2.LVc, type = "l", col = "darkblue") 
  
  # # according to LVc 
  # N.LVc = rep(NA, length(time)) 
  # N.LVc[1] = N0 
  # dN = 0 
  #  
  # for (i in 2:(length(time)-1)) { 
  #   dN = r * N.LVc[i-1] * (1 - alpha*N.LVc[i-1]) 
  #   N.LVc[i] = N.LVc[i-1] + dN 
  # } 
  #  
  # points(x = time, y = N.LVc, type = "l", col = "red") 
  
} 

# parameters 
# r2 = mod$coefficients[2] 
# r2 = 0.033 
# alpha22 = 0.0027 
alpha11 = 0.002 
# alpha21 = 0.0005 
# alpha12 = 0.001 
# alpha21 = mean(alpha21est) 
# alpha12 = mean(alpha12est) 
N2init = bothPreys$prey2PopulationSize[1] 
N1init = bothPreys$prey1PopulationSize[1] 

plotFit2preys(dataT = bothPreys$timeStep, dataN1 = bothPreys$prey1PopulationSize, dataN2 = bothPreys$prey2PopulationSize, dataP = bothPreys$predator1PopulationSize, tIntro = 201, tMax = max(bothPreys$timeStep), r1 = r2, r2 = r2, alpha22 = alpha22, alpha11 = alpha11, alpha21 = alpha21, alpha12 = alpha12, N20 = N2init, N10 = N1init) 
