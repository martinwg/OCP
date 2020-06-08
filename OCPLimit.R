## These functions were used to generate the h limits for the OCP method in the paper.
## The function hBinarySearch refers to Alogrithim 2 in the paper.

#!/usr/bin/Rscript
library(doParallel)
library(mvnfast)
library(RhpcBLASctl)

df <- read.csv('final_grid_points.csv')[,-c(1,2)]
head(df)

generate_data <- function(obsin, obsout, p, s, dist, sim, h1) 
{
  FP <- c()
  cl=makeCluster(4)
  registerDoParallel(cl)
  iter<-foreach (i=1:sim) %dopar%
  {
    if (dist=="Normal") {
      gennorm = function(obsin, p, s, obsout = 0, deltap = 0.023){
        require(mvnfast)
        mu0 = (rep(0, p))
        p1 = sample(1:p, 1) ### Shifting Random Dimensions
        sigma <-matrix(s, p, p) + diag((1-s), p)
        din = rmvn(obsin, mu0, sigma,ncores=1)
        dout = c()
        return(list(datain = din, dataout = dout, data = rbind(din,dout)))
      }
      data <- gennorm(obsin, p, s, obsout)
      
    } else if (dist=="Lognormal") { 
      genlognorm = function(obsin, p, s, obsout = 0,  deltap = 0.023){
        require(mvnfast)
        mu0 = (rep(0, p))
        p1 = sample(1:p, 1) ### Shifting Random Dimensions
        sigma <-matrix(s, p, p) + diag((1-s), p)
        din = exp(rmvn(obsin, mu0, sigma,ncores=4))
        dout = c()
        return(list(datain = din, dataout = dout, data = rbind(din,dout)))
      }
      data <-  genlognorm(obsin, p, s, obsout)
      
    } else {
      gent = function(obsin, p, s, obsout = 0, deltap = 0.023){
        require(mvnfast)
        mu0 = (rep(0, p))
        p1 = sample(1:p, 1) ### Shifting Random Dimensions
        sigma <-matrix(s, p, p) + diag((1-s), p)
        din = mvnfast::rmvt(obsin, mu0,sigma, df = 10,ncores=4)
        dout = c()
        return(list(datain = din, dataout = dout, data = rbind(din,dout)))
      }
      data <- gent(obsin, p, s, obsout)
    }
    
    OCP = function(dataframe, peeled_to_obs = 2, h = NULL, plot = TRUE){
      require(kernlab)
      df_oc = data.frame(type=1,dataframe)
      N = dim(df_oc)[1]
      m = dim(df_oc)[1]
      p = dim(df_oc)[2]-1
      invsigma = 1/p
      d = 0
      o1 = 0
      bagdat = df_oc;
      while (m > peeled_to_obs)
      {
        o1 = o1+1;
        d = d+1 ;  
        if (class(bagdat[,-1]) =="numeric"){mu_OCP=mean(bagdat[,-1])} else {mu_OCP = colMeans(bagdat[,-1])}
        w = ksvm(type ~., data=bagdat, type="one-svc", kernel="rbfdot", kpar=list(sigma=invsigma), nu=0.0001, scaled = F)
        sv = SVindex(w);
        bagdat = bagdat[-sv,]
        m = dim(bagdat)[1]  
      }
      if (class(dataframe) =="numeric"){dA = dataframe**2} else {dA = rowSums(dataframe*dataframe)}
      dB = sum(mu_OCP*mu_OCP)
      D = (dA+dB) - as.numeric(2*as.matrix(dataframe)%*%mu_OCP)
      D[D<0] = 0
      KS_OCP = 1 - exp(-D/(p*p))
      MAD = median(abs(KS_OCP-median(KS_OCP)));
      RKD_OCP = (KS_OCP  - median(KS_OCP))/MAD;
      if(is.null(h)) {h <- quantile(RKD_OCP, 0.75) + 1.5 * IQR(RKD_OCP)}
      OCP_Flag = RKD_OCP > h
      if (plot == TRUE){
        plot (RKD_OCP, type = 'l', col = 'grey', lwd = 2,cex.main=1,  xlab="Observation Index", ylab="RKDs")
        par(new = TRUE)
        plot (RKD_OCP, col = ifelse(OCP_Flag > 0,'red', NA),  xaxt='n', ann=FALSE)
        par(new = TRUE)
        abline (h = h, lty=3, col = 'red', lwd = 1.8)
      }
      return(list(RKD_OCP= RKD_OCP, OCP_Flag = OCP_Flag))
    }
    
    ocp_res = OCP(data$data)
    FP[i] = sum(ocp_res$RKD_OCP > h1) / obsin
  }
  stopCluster(cl)
  iters<-do.call(rbind, iter)
  return(mean(iters))
}

## Binary Search Function
hBinarySearch <- function(l_pred_bound, u_pred_bound, pred_h, obs, p, dist, s, desired_FAP = 0.05, tol = 0.003, niter = 100, sim = 50)
{
  l <- l_pred_bound
  u <- u_pred_bound
  m <- pred_h
  D = generate_data(obs,0,p,s, dist, sim, m)
  test=abs((D-desired_FAP))
  while(test > tol){
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
    D = generate_data(obs,0,p,s, dist, sim, m)
    test=abs((D-desired_FAP))
    if(test<.005){
      sim=500
    }else if (test<.01){
      sim=100
    }
    cat("The value of h is: ", m, "\n")
    cat("Difference from desired FAP: ", test, "\n")
    cat("Actual FAP is: ",D, "\n")
    cat("\n")
    cat("\n")
  }
  return(list(h_value = m, FP_value = D))
}

## Example
dist = "Normal"
obs = 50
p = 200
corr = 0.1

## Inverse Prediction (if not converging STOP, provide new upper, lower and estimate values)
library(investr)
reg<-lm(avg_FP ~ (h + corr + obs + p + dist)^3, data=df)
summary(reg)
FP_pred  <- invest(reg, y0=.05, interval=c("inversion"), level=0.99, x0.name="h", newdata = data.frame(dist= dist, obs=obs, p=p, corr=corr), lower=-10, upper=20)

result <- hBinarySearch(FP_pred$lower, FP_pred$upper, FP_pred$estimate, obs, p,  dist, s = corr)
result$h_value
result$FP_value


