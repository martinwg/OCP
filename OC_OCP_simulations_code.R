library(MASS)
library(crossval)
library(doParallel)

#####################################################
###### Function to Generate the random data #########
#####################################################

generate_data <- function(obsin, obsout, p, s, dist, seed) 
{
  set.seed(seed)
  if (dist=="Normal") {
    gennorm = function(obsin, p, s, obsout = 0, deltap = 0.023){
      require(mvtnorm)
      mu0 = (rep(0, p))
      p1 = sample(1:p, 1) ### Shifting Random Dimensions
      sigma <-matrix(s, p, p) + diag((1-s), p)
      mu1a = qmvnorm(1-deltap, sigma  = sigma)$quantile
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
      return(list(datain = din, dataout = dout, data = rbind(din,dout), dim_shifted = p1))
    }
    data <- gennorm(obsin, p, s, obsout)
    
  } else if (dist=="Lognormal") { 
    genlognorm = function(obsin, p, s, obsout = 0,  deltap = 0.023){
      require(mvtnorm)
      mu0 = (rep(0, p))
      p1 = sample(1:p, 1) ### Shifting Random Dimensions
      sigma <-matrix(s, p, p) + diag((1-s), p)
      mu1a = qmvnorm(1-deltap, sigma  = sigma)$quantile
      mu1a = rep(mu1a, p) 
      rsigns = sample(c(-1,1),p, replace = T)
      mu1b = sqrt(mu1a%*%mu1a/p1)
      mu1 = c(rep(mu1b, p1), rep(0,(p-p1))) ### Shifting Random Dimensions
      mu1 = mu1 * rsigns
      din = exp(rmvnorm(obsin, mean=mu0, sigma = sigma))
      if (obsout > 0){
        dout = exp(rmvnorm(obsout, mean=mu1, sigma = sigma))
      }
      else
      {
        dout = c()
      }
      return(list(datain = din, dataout = dout, data = rbind(din,dout), dim_shifted = p1))
    }
    data <-  genlognorm(obsin, p, s, obsout)
    
  } else {
    gent = function(obsin, p, s, obsout = 0, deltap = 0.023){
      require(mvtnorm)
      mu0 = (rep(0, p))
      p1 = sample(1:p, 1) ### Shifting Random Dimensions
      sigma <-matrix(s, p, p) + diag((1-s), p)
      mu1a = qmvt(1-deltap, sigma  = sigma, df = 10)$quantile
      mu1a = rep(mu1a, p) 
      rsigns = sample(c(-1,1),p, replace = T)
      mu1b = sqrt(mu1a%*%mu1a/p1)
      mu1 = c(rep(mu1b, p1), rep(0,(p-p1))) ### Shifting Random Dimensions
      mu1 = mu1 * rsigns
      din = rmvt(obsin, sigma = sigma, df = 10)
      if (obsout > 0){
        dout = rmvt(obsout, sigma = sigma, df = 10) + matrix(rep(mu1,obsout), ncol = p, byrow = T)
      }
      else
      {
        dout = c()
      }
      return(list(datain = din, dataout = dout, data = rbind(din,dout), dim_shifted = p1))
    }
  data <- gent(obsin, p, s, obsout)
  }
  return(data)
}

#####################################################
######            OCP Function                  #####
#####################################################

OCP = function(dataframe, pto = 2){
  require(kernlab)
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
    w = ksvm(type ~., data=bagdat, type="one-svc", kernel="rbfdot", kpar=list(sigma=invsigma), nu=0.0001, cross=10, scaled = F)
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
  MAD = median(abs(KS_OCP-median(KS_OCP)))
  RKD_OCP = (KS_OCP  - median(KS_OCP))/MAD;
  return(list(RKD_OCP= RKD_OCP, mu_OCP = mu_OCP))
}


#####################################################
######      Performance Metrics Function        #####
#####################################################

perf_metrics <- function (RKD_OCP, cutoff, obsout)
{
  N = length(RKD_OCP)
  pred_v <- as.numeric(RKD_OCP > cutoff)
  actual_v <- c(rep(0, N - obsout), rep(1, obsout))
  diff = sum(abs(pred_v - actual_v))
  Percent_CC = (1 - (diff /  N)) * 100
  Detect_R = (sum(pred_v[(N-obsout+1):N]) / obsout) * 100
  return(list(Percent_CC = Percent_CC, Detect_R = Detect_R))
}


#####################################################
######                  Example                 #####
#####################################################

obsin = 100
obsout = 20
p = 200
s = 0.5
dist = "Lognormal"
seed = 12
cutoff = 4  ## This is the value obtained through simulation

dataset = generate_data(obsin, obsout, p, s, dist, seed)
ocp_res = OCP(dataset$data)
perf_metrics(ocp_res$RKD_OCP, cutoff, obsout)


#####################################################
######              Simulations                 #####
#####################################################

sims = 5
obs = 50
p = 100
s = 0.5
dist = "Normal"
cutoff = 4  ## This is the value obtained through simulation
percent_out = 0.05

obsout = round(percent_out * obs)
obsin = obs - obsout
metrics_combined = matrix(rep(0, 2*sims), ncol = 2)

for (i in 1:sims)
{
  dataset = generate_data(obsin, obsout, p, s, dist, i)
  ocp_res = OCP(dataset$data)
  met = perf_metrics(ocp_res$RKD_OCP, cutoff, obsout)
  metrics_combined[i,1] = met$Percent_CC
  metrics_combined[i,2] = met$Detect_R
}

colnames(metrics_combined) = c("Percent_CC", "Detect_R")
avgs = colMeans(metrics_combined)

