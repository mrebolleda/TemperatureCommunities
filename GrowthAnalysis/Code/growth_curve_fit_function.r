#######################################################
#### Growth curve fitting function ###################
######################################################

# Author: Djordje Bajic (D.Bajic@tudelft.nl)
# Last edited by: Maria Rebolleda-Gomez (mreboll1@uci.edu)

# In order to minimize assumptions imposed to the data, we are fitting
# a generalized additive model using Djordje Bajic's function. 

# The problem is the model becomes extremely sensitive to moments that 
# do not always reflect the overall trend. Those include:
  # the FIRST HOURS were noise is closer to the scale of the measured values
  # so now growth rate is only measured after a minimum time (2 hours). 
  # SMALL PEAKS within the overall trend, that we can smooth averaging 
  # the derivative 1 hour around the max. 


compute.gam <- function(x)
{
  
  id <- x %>% select(-t, -abs, -lOD) %>% unique()
  # First model the raw data, which is needed for MODEL THE RAW DATA, NEEDED FOR THE DERIVATIVE VAL. AT MAX GR
  gam0 <- gam(abs ~ s(t, bs = 'cs'), data=x)
  #gam0 <- gam(abs ~ s(t, bs = 'ad'), data=x)
  
  t <- unique(x$t)
  newd <- data.frame(t = t)
  pred0 <- predict(gam0, newd,
                   se.fit = TRUE) %>%
    as_tibble() %>%
    mutate(t = t)
  
  X0 <- predict(gam0, newd, type="lpmatrix")
  
  eps <- 1e-5 ## finite difference interval
  newd <- data.frame(t = t+eps)
  X1 <- predict(gam0, newd, type="lpmatrix")
  
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  
  df <- Xp%*%coef(gam0)              ## ith smooth derivative
  df.sd <- rowSums(Xp%*%gam0$Vp*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  #plot(t, df, type="l", ylim=range(c(df+2*df.sd,df-2*df.sd)))
  #lines(t,df+2*df.sd,lty=2)
  #lines(t,df-2*df.sd,lty=2)
  
  pred0 <- pred0 %>%
    mutate(deriv.fit = df[,1], deriv.sd = df.sd) %>%
    merge(id)
  
  # MODEL THE LOG OD
  #gam1 <- gam(lOD ~ s(t, bs = 'ad'), data=x)
  gam1 <- gam(lOD ~ s(t, bs = 'cs'), data=x)
  t <- unique(x$t)
  newd <- data.frame(t = t)
  pred <- predict(gam1, newd,
                  se.fit = TRUE) %>%
    as_tibble() %>%
    mutate(t = t)
  
  X0 <- predict(gam1, newd, type="lpmatrix")
  
  eps <- 1e-5 ## finite difference interval
  newd <- data.frame(t = t+eps)
  X1 <- predict(gam1, newd, type="lpmatrix")
  
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  
  
  df <- Xp%*%coef(gam1)              ## ith smooth derivative
  df.sd <- rowSums(Xp%*%gam1$Vp*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  #plot(t, df, type="l", ylim=range(c(df+2*df.sd,df-2*df.sd)))
  #lines(t,df+2*df.sd,lty=2)
  #lines(t,df-2*df.sd,lty=2)
  
  ## plot(gam1)
  
  # GET THE DATA
  pred <- pred %>%
    mutate(deriv.fit = df[,1], deriv.sd = df.sd,
           OD.frac = pred0$fit/max(pred0$fit),
           lOD = x$lOD) %>%
    merge(id)
  # one table has data plus fitted value, plus standard dev from fitted value
  # other table has summarized parameters - calculates growth rate, etc
  
  # GET GROWTH PARAMETERS
  ## compute lag
  t0 <- min(x$t)

  # Maximum growth rate is too sensitive to artifacts. Instead we can get the
  # median rate within an hour of the maximum.

  maxgr <- pred %>%
    filter(t >= 2) %>%
    filter(deriv.fit==max(deriv.fit))
  t.maxGr <- maxgr$t[1]
  
  OD.maxGr <- exp(maxgr$fit[1])
  
  p0 <- setDT(pred0)[t==t.maxGr]
  slope.maxGr <- p0$deriv.fit
  
  ## y = y_0 + slope * t, so:
  y_0 <- OD.maxGr - (slope.maxGr*t.maxGr)
  
  # I take as lag the intersection with y = min (OD), instead of y=0
  t_lag <- (min(pred0$fit)-y_0)/slope.maxGr # lag

  
  prm <- tibble(r = maxgr$deriv.fit[1],
                t.r = t.maxGr, # time at which it reaches max growth
                lag = t_lag, 
                maxOD.fit = max(pred0$fit)) %>%
    merge(id)
  
  return(list('data' = pred,
              'params' = prm))
}
