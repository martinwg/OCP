#!/usr/bin/Rscript

df <- read.csv('final_grid_points.csv')[,-c(1,2)]
head(df)

## If we want to model counts:
df$counts<-df$avg_FP*1000

## Regression model:
## h = sRKD, corr = estimated correlation off-diagonals, obs = N, p = dimension
## dist = distribution
reg<-lm(avg_FP ~ (h + corr + obs + p + dist)^3, data=df)
summary(reg)

## Binary Search Function
hBinarySearch <- function(l_pred_bound, u_pred_bound, pred_h, obs, p, dist = "Normal", s = 0.5, desired_FAP = 0.05, tol = 1e-03, niter = 1000, sim = 500)
{
  l <- l_pred_bound
  u <- u_pred_bound
  m <- pred_h
  #m <- (l_pred_bound + u_pred_bound)/2
  D <- avg_FP_calc (obs, p, m, s, dist, sim) 
  #D <- BinaryEntrp (xold, m) ## This is just the function we want to get to approximate root
  for (i in 1:niter)
  {
    if (D < desired_FAP)
    {
      u = m
      l = l
      m = (u+l)/2
    } 
    
    else
    {
      l = m
      u = u
      m = (u+l)/2
    }
    
    D <- avg_FP_calc (obs, p, m, s, dist, sim)  ## This is just the function we want to get to approximate root
    if (abs((D-desired_FAP)) < tol)
      return(list(h_value = m, FP_value = D)) 
  }
  #print (abs(D-desired_FAP))
  #stop ("Exceeded")
}

## Threshold Simulations
threshold_simulations <- function(obs, p, h1, s, dist, sim)
{
  FP <- c()
  for(i in 1:sim) 
  {
    set.seed(i)
    if (dist=="Normal") {
      gennorm_corr_S_M = function(obsin, p, s, obsout = 0, deltap = 0.023){
        require(mvtnorm)
        mu0 = (rep(0, p))
        p1 = sample(1:p, 1) ### Shifting Random Dimensions
        sigma <-matrix(s, p, p) + diag((1-s), p)
        mu1a = 20
        mu1a = rep(mu1a, p) 
        rsigns = sample(c(-1,1),p, replace = T)
        mu1b = sqrt(mu1a%*%mu1a/p1)
        mu1 = c(rep(mu1b, p1), rep(0,(p-p1))) ### Shifting Random Dimensions
        mu1 = mu1 * rsigns
        din = rmvnorm(obsin, mean=mu0, sigma = sigma)
        if (obsout > 0){
          dout = rmvnorm(obsout, mean=mu1, sigma = sigma)
        }
        else
        {
          dout = c()
        }
        
        return(list(datain = din, dataout = dout, data = rbind(din,dout), p1 = p1))
      }
      data <- gennorm_corr_S_M(obs, p, s)
      
    } else if (dist=="Lognormal") { 
      genlognorm_corr_S_L = function(obsin, p, s, obsout = 0,  deltap = 0.023){
        require(mvtnorm)
        mu0 = (rep(0, p))
        p1 = sample(1:p, 1) ### Shifting Random Dimensions
        sigma <-matrix(s, p, p) + diag((1-s), p)
        din = exp(rmvnorm(obsin, mean=mu0, sigma = sigma))
        if (obsout > 0){
          dout = exp(rmvnorm(obsout, mean=mu1, sigma = sigma))
        }
        else
        {
          dout = c()
        }
        return(list(datain = din, dataout = dout, data = rbind(din,dout), p1 = p1))
      }
      
      data <-  genlognorm_corr_S_L(obs, p, s)
      
    } else {
      gent_corr_S_L = function(obsin, p, s, obsout = 0, deltap = 0.023){
        require(mvtnorm)
        mu0 = (rep(0, p))
        p1 = sample(1:p, 1) ### Shifting Random Dimensions
        sigma <-matrix(s, p, p) + diag((1-s), p)
        #sigma <-matrix(0.2, p, p) + diag(0.8, p)
        #mu1a = qmvt(deltap, sigma  = sigma, df = 10)$quantile
        #mu1a = 20
        #mu1a = rep(mu1a, p) 
        #rsigns = sample(c(-1,1),p, replace = T)
        #mu1b = sqrt(mu1a%*%mu1a/p1)
        #mu1 = c(rep(mu1b, p1), rep(0,(p-p1))) ### Shifting Random Dimensions
        #mu1 = mu1 * rsigns
        din = rmvt(obsin, sigma = sigma, df = 10)
        if (obsout > 0){
          dout = rmvt(obsout, sigma = sigma, df = 10) + matrix(rep(mu1,obsout), ncol = p, byrow = T)
        }
        else
        {
          dout = c()
        }
        return(list(datain = din, dataout = dout, data = rbind(din,dout), p1 = p1))
      }
      
      data <- gent_corr_S_L(obs, p, s)
      
    }
    
    OCP = function(dataframe, pto = 2){
      require(kernlab)
      require(car)
      df_oc = data.frame(type=1,dataframe)
      N = dim(df_oc)[1]
      m = dim(df_oc)[1]
      p = dim(df_oc)[2]
      invsigma = 1/p
      d = 0
      o1 = 0
      bagdat = df_oc;
      while (m > pto)
      {
        o1 = o1+1;
        d = d+1 ;  
        if (class(bagdat[,-1]) =="numeric"){mu_OCP=mean(bagdat[,-1])} else {mu_OCP = colMeans(bagdat[,-1])}
        w = ksvm(type ~., data=bagdat, type="one-svc", kernel="rbfdot", kpar=list(sigma=invsigma),nu=0.0001,cross=10, scaled = F)
        sv = SVindex(w);
        bagdat = bagdat[-sv,]
        m = dim(bagdat)[1]  
      }
      if (class(dataframe) =="numeric"){dA = dataframe**2} else {dA = rowSums(dataframe*dataframe)}
      dB = sum(mu_OCP*mu_OCP)
      D = (dA+dB) - as.numeric(2*as.matrix(dataframe)%*%mu_OCP)
      D[D<0] = 0
      KS_OCP = 1 - exp(-D/(p*p))
      #MAD = 1.3238*median(abs(KS_OCP-median(KS_OCP)));
      #k<-powerTransform(KS_OCP, family="bcPower")
      #if (k$lambda<0){
      KS_OCP<-log(KS_OCP+0.001)
      #     }else{
      #      KS_OCP<-(KS_OCP+0.001)^k$lambda}
      MAD = median(abs(KS_OCP-median(KS_OCP)));
      RKD_OCP = (KS_OCP  - median(KS_OCP))/MAD;
      return(list(RKD_OCP= RKD_OCP, mu_OCP = mu_OCP))
    }
    
    #data <- gennorm_uncorr_S(obs, p)
    OCP_result <- OCP(data$data)
    FP[i] = sum(OCP_result$RKD_OCP > h1) / obs 
  }
  #proc.time() - ptm
  #iters<-do.call(rbind, iter)
  
  #ji <- ji + 1
  #avg_FP[ji] <- mean(FP)
  #avg_FP[ji] <- mean(iters)
  #}
  #proc.time() - ptm
  #avg_FP <- data.frame(avg_FP)
  #rownames(avg_FP) <- h1
  return(mean(FP))
}

avg_FP_calc <- function (obs, p, m, s, dist, sim)
{
  #library(doParallel)
  library(mvtnorm)
  #avg_FP <- rep(0, sim)
  #cl<-makeCluster(12)
  #registerDoParallel(cl)
  #ptm <- proc.time()
  #for (i in 1:sim){
  avg_FP<- threshold_simulations(obs, p, m, s, dist, sim) 
  #avg_FP[i] <- threshold_simulations(sim, obs=obs[i], p=p[i], h1=h[i], dist=dist[i],  s=corr[i])
   # print(i)
  #}
  #proc.time() - ptm
  #stopCluster(cl)
  return(avg_FP)
}


#Example
dist = "Normal"
obs = 500
p = 2000
corr = 0.25

## Inverse Prediction
library(investr)
FP_pred  <- invest(reg, y0=.05, interval=c("inversion"), level=0.99, x0.name="h", newdata = data.frame(dist= dist, obs=obs, p=p, corr=corr), lower=-10, upper=20)


## Example 
result <- hBinarySearch(FP_pred$lower, FP_pred$upper, FP_pred$estimate, obs, p,  dist, s = corr)
result$h_value
result$FP_value

write.csv(result, 'Normal_N500_p2000_c0_25.csv')