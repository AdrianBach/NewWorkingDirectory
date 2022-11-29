lvc <- function(tMax, r1, r2, alpha11, alpha22, alpha21, alpha12, N10, N20) {
  
  # generate a time vector
  time <- seq(0, tMax)
  
  # generate N2 according to LV
  
  # vector for Ns
  N1.solo = rep(NA, length(time))
  N1.solo[1] = N10
  dN1.solo = 0
  
  N1.LVc = rep(NA, length(time))
  N1.LVc[1] = N10
  dN1 = 0
  
  N2.solo = rep(NA, length(time))
  N2.solo[1] = N20
  dN2.solo = 0
  
  N2.LVc = rep(NA, length(time))
  N2.LVc[1] = N20
  dN2 = 0
  
  for (i in 2:length(time)) {
      
    # compute delta Ns 
    dN1.solo = r1 * N1.solo[i-1] * (1 - alpha11 * N1.solo[i-1])
    dN1 = r1 * N1.LVc[i-1] * (1 - alpha11 * N1.LVc[i-1] - alpha12 * N2.LVc[i-1])
    dN2.solo = r2 * N2.solo[i-1] * (1 - alpha22 * N2.solo[i-1])
    dN2 = r2 * N2.LVc[i-1] * (1 - alpha22 * N2.LVc[i-1] - alpha21 * N1.LVc[i-1])
    
    # calculate Ns(t)
    N1.solo[i] = N1.solo[i-1] + dN1.solo
    N1.LVc[i] = N1.LVc[i-1] + dN1
    N2.solo[i] = N2.solo[i-1] + dN2.solo
    N2.LVc[i] = N2.LVc[i-1] + dN2

  }
  
  plot(x = 1, y = 1, main = "LVc equation", xlab = "Generations", ylab = "Population size", type = "n", 
       xlim = c(min(time), max(time)),
       ylim = c(min(c(min(N1.solo), min(N2.solo)))-10, max(c(max(N1.solo), max(N2.solo)))+10))
  lines(x = time, y = N1.solo, lwd = 1.5, lty = 2, col = "darkred")
  lines(x = time, y = N2.solo, lwd = 1.5, lty = 2, col = "darkblue")
  lines(x = time, y = N1.LVc,  lwd = 1.5, lty = 1, col = "darkred")
  lines(x = time, y = N2.LVc,  lwd = 1.5, lty = 1, col = "darkblue")
  
}

# parameters
# r2 = mod$coefficients[2]
Alpha11 = 0.002
Alpha22 = 0.0027
Alpha21 = 0.001
Alpha12 = 0.001

lvc(tMax = 1000, r1 = 0.033, r2 = 0.033, alpha22 = Alpha22, alpha11 = Alpha11, alpha21 = Alpha21, alpha12 = Alpha12, N10 = 25, N20 = 25)

#### including predators ####
lvc.P <- function(tMax, r1, r2, rp, alpha11, alpha22, alpha21, alpha12, alphap1, alphap2, alpha1p, alpha2p, N10, N20, P0) {
  
  # generate a time vector
  time <- seq(0, tMax-1)
  
  # generate N2 according to LV
  
  # vector for Ns
  N1.solo = rep(NA, length(time))
  N1.solo[1] = N10
  dN1.solo = 0
  
  N1.LVc = rep(NA, length(time))
  N1.LVc[1] = N10
  dN1 = 0
  
  N2.solo = rep(NA, length(time))
  N2.solo[1] = N20
  dN2.solo = 0
  
  N2.LVc = rep(NA, length(time))
  N2.LVc[1] = N20
  dN2 = 0
  
  P.solo.N1 = rep(NA, length(time))
  P.solo.N1[1] = P0
  dP.solo.N1 = 0
  
  P.solo.N2 = rep(NA, length(time))
  P.solo.N2[1] = P0
  dP.solo.N2 = 0
 
  P.LVc = rep(NA, length(time))
  P.LVc[1] = P0
  dP = 0
  
  for (i in 2:length(time)) {
    
    # compute delta Ns 
    dN1.solo = r1 * N1.solo[i-1] * (1 - alpha11 * N1.solo[i-1] - alpha1p * P.solo.N1[i-1])
    dN1 = r1 * N1.LVc[i-1] * (1 - alpha11 * N1.LVc[i-1] - alpha12 * N2.LVc[i-1] - alpha1p * P.LVc[i-1])
    dN2.solo = r2 * N2.solo[i-1] * (1 - alpha22 * N2.solo[i-1] - alpha2p * P.solo.N2[i-1])
    dN2 = r2 * N2.LVc[i-1] * (1 - alpha22 * N2.LVc[i-1] - alpha21 * N1.LVc[i-1] - alpha2p * P.LVc[i-1])
    dP.solo.N1 = rp * P.solo.N1[i-1] * (1 - alphap1 * N1.solo[i-1])
    dP.solo.N2 = rp * P.solo.N2[i-1] * (1 - alphap2 * N2.solo[i-1])
    dP = rp * P.LVc[i-1] * (1 - alphap1 * N1.LVc[i-1] - alphap2 * N2.LVc[i-1])
    
    # calculate Ns(t)
    N1.solo[i] = N1.solo[i-1] + dN1.solo
    N1.LVc[i] = N1.LVc[i-1] + dN1
    N2.solo[i] = N2.solo[i-1] + dN2.solo
    N2.LVc[i] = N2.LVc[i-1] + dN2
    P.solo.N1[i] = P.solo.N1[i-1] + dP.solo.N1
    P.solo.N2[i] = P.solo.N2[i-1] + dP.solo.N2
    P.LVc[i] = P.LVc[i-1] + dP
    
    # round them
    N1.solo[i] = ifelse(N1.solo[i] > 0, round(N1.solo[i]) , 0)
    N1.LVc[i]  = ifelse(N1.LVc[i]  > 0, round(N1.LVc[i]), 0)
    N2.solo[i] = ifelse(N2.solo[i] > 0, round(N2.solo[i]), 0)
    N2.LVc[i]  = ifelse(N2.LVc[i]  > 0, round(N2.LVc[i]), 0)
    P.solo.N1[i] = ifelse(P.solo.N1[i] > 0, round(P.solo.N1[i]), 0)
    P.solo.N2[i] = ifelse(P.solo.N2[i] > 0, round(P.solo.N2[i]), 0)
    P.LVc[i] = ifelse(P.LVc[i] > 0, round(P.LVc[i]), 0)
    
  }
  
  plot(x = 1, y = 1, main = "LVc equation", xlab = "Generations", ylab = "Population size", type = "n", 
       xlim = c(min(time), max(time)),
       ylim = c(min(c(min(N1.solo), min(N2.solo), min(P.solo.N1), min(P.solo.N2)))-10, max(c(max(N1.solo), max(N2.solo), max(P.LVc)))+10))
  lines(x = time, lwd = 1.5, lty = 2, col = "darkred" , y = N1.solo)
  lines(x = time, lwd = 1.5, lty = 3, col = "darkblue", y = N2.solo)
  lines(x = time, lwd = 1.5, lty = 2, col = "orange" , y = P.solo.N1)
  lines(x = time, lwd = 1.5, lty = 3, col = "orange", y = P.solo.N2)
  lines(x = time, lwd = 1.5, lty = 1, col = "darkred" , y = N1.LVc )
  lines(x = time, lwd = 1.5, lty = 1, col = "darkblue", y = N2.LVc )
  lines(x = time, lwd = 1.5, lty = 1, col = "orange" , y = P.LVc )
  
}

# parameters
# r2 = mod$coefficients[2]
Alpha11 = 1/500
Alpha22 = 1/500
Alpha21 = 0.001
Alpha12 = 0.002
Alpha1p = 0.003
Alpha2p = 0.004
Alphap1 = -1.1
Alphap2 = -1.2

# lvc(tMax = 1000, 
    # r1 = 0.033, r2 = 0.033, rp = -1, 
    # alpha11 = Alpha11, alpha22 = Alpha22, alpha12 = Alpha12, alpha21 = Alpha21,
    # alphap1 = Alphap1, alphap2 = Alphap2, alpha1p = Alpha1p, alpha2p = Alpha2p,
    # N10 = 50, N20 = 50, P0 = 20)
lvc.P(tMax = 1000, 
      r1 = 0.033, r2 = 0.033, rp = -1, 
      alpha11 = 0.002, alpha22 = 0.002, alpha12 = 0, alpha21 = 0,
      alphap1 = 0.002, alphap2 = 0.003, alpha1p = 0.01, alpha2p = 0.02,
      N10 = 100, N20 = 100, P0 = 20)

