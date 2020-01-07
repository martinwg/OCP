################################################ OCP R Implementation ####################################################
#                                              One-Class Peeling Method                                                  #
##########################################################################################################################

## To determine a particular value of h for a dataset, use the binary search procedure in OCPLimits.R. 
## If h=NULL the robust limits are used by default (Q3+1.5(IQR)) where Q3 is found from the RKDs.
## Per default, data are scaled internally in the ksvm function to zero mean and unit variance.

## Example
## Inliers: Correlated Normal (mu = 0)  N = 100, p = 200, (corr = 0.25)  
## Outliers: 20 from Correlated Normal N = 100, p = 200 (corr = 0.25)  
# library(mvtnorm)
# p = 200
# corr = 0.25
# gamma = 0.023
# sigma <-matrix(corr, p, p) + diag((1-corr), p)
# nquantile = qmvnorm(1-gamma, sigma  = sigma)$quantile
# d1 = rmvnorm(100, mean=(rep(0, p)), sigma = sigma)
# d2 = rmvnorm(20, mean= rep(nquantile, p), sigma = sigma)
# datafinal <- rbind(d1,d2)
# OCP_result = OCP(datafinal, standardize = FALSE)

OCP = function(dataframe, peeled_to_obs = 2, h = NULL, standardize = TRUE, plot = TRUE){
  require(kernlab)
  if(standardize == TRUE){
    dataframe = data.frame(scale(dataframe, center = TRUE, scale = TRUE))
  }
  df_oc = data.frame(type=1,dataframe)
  N = dim(df_oc)[1]
  m = dim(df_oc)[1]
  p = dim(df_oc)[2]
  invsigma = 1/p
  d = 0
  o1 = 0
  bagdat = df_oc;
  while (m > peeled_to_obs)
  {
    o1 = o1+1;
    d = d+1 ;  
    if (class(bagdat[,-1]) =="numeric"){mu_OCP=mean(bagdat[,-1])} else {mu_OCP = colMeans(bagdat[,-1])}
    w = ksvm(type ~., data=bagdat, type="one-svc", kernel="rbfdot", kpar=list(sigma=invsigma), nu=0.0001, scaled = F)  ## The scaled
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
    plot (RKD_OCP, type = 'l', col = 'grey', lwd = 2,cex.main=1,  xlab="Observation Index", ylab="sRKDs")
    par(new = TRUE)
    plot (RKD_OCP, col = ifelse(OCP_Flag > 0,'red', NA),  xaxt='n', ann=FALSE)
    par(new = TRUE)
    abline (h = h, lty=3, col = 'red', lwd = 1.8)
  }
  return(list(RKD_OCP= RKD_OCP, OCP_Flag = OCP_Flag))
}
