setwd("/home/martinwg/OCP")
a = 12
write.csv(a, "a.csv")
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

#b = head(data)
#write.csv(b, "b.csv")

RMDP<-function(y,alpha=0.05,itertime=100)   ############RMDP procedure
{
  n<-nrow(y)
  p<-ncol(y)
  h<-round(n/2)+1
  init_h=2
  delta=alpha/2
  bestdet=0
  jvec<-array(0,dim=c(n,1))
  for(A in 1:itertime)
  {
    id=sample(n,init_h)
    ny=y[id,]
    mu_t=apply(ny,2,mean)
    var_t=apply(ny,2,var) 
    dist<-apply((t((t(y)-mu_t)/var_t^0.5))^2,1,sum)
    crit=10
    l=0
    while(crit!=0&l<=15)
    {	l=l+1
    ivec<-array(0,dim=c(n,1))
    dist_perm<-order(dist)
    ivec[dist_perm[1:h]]=1
    crit=sum(abs(ivec-jvec))
    jvec=ivec
    newy=y[dist_perm[1:h],]
    mu_t=apply(newy,2,mean)
    var_t=apply(newy,2,var) 
    dist<-apply((t((t(y)-mu_t)/var_t^0.5))^2,1,sum)
    }
    tempdet=prod(var_t)
    if(bestdet==0|tempdet<bestdet) 
    {bestdet=tempdet
    final_vec=jvec
    }
  }
  submcd<-seq(1,n)[final_vec!=0]
  ##Estimated Mean and Variance
  mu_t=apply(y[submcd,],2,mean)
  var_t=apply(y[submcd,],2,var) 
  ##Squared Disatances
  dist<-apply((t((t(y)-mu_t)/var_t^0.5))^2,1,sum)
  ##Squared Distances rescaled
  dist=dist*p/median(dist)
  ###Correlation matrix of the minimum volume subset
  ER=cor(y[submcd,])
  ##Squared correlation matrix of the minimum volume subset
  ER=ER%*%ER
  ##Calculating argument for inside of the Z{}function in eq (9)
  tr2_h=sum(diag(ER))
  tr2=tr2_h-p^2/h
  cpn_0=1+(tr2_h)/p^1.5
  ##Calculating weight from eq (9) w0 returns a TRUE or FALSE for each value in Sample
  w0=(dist-p)/sqrt(2*tr2*cpn_0)<qnorm(1-delta)
  nw=sum(w0)
  ##Data is subset based on TRUE w0 values
  sub=seq(1,n)[w0]	
  mu_t=apply(y[sub,],2,mean)
  var_t=apply(y[sub,],2,var) 
  ER=cor(y[sub,])
  ER=ER%*%ER
  tr2_h=sum(diag(ER))
  tr2=tr2_h-p^2/nw
  dist<-apply((t((t(y)-mu_t)/var_t^0.5))^2,1,sum)
  scale=1+1/sqrt(2*pi)*exp(-qnorm(1-delta)^2/2)/(1-delta)*sqrt(2*tr2)/p
  dist=dist/scale
  cpn_1=1+(tr2_h)/p^1.5
  w1=(dist-p)/sqrt(2*tr2*cpn_1)<qnorm(1-alpha)
  list(w1=w1)
  return(dis_rmdp = (dist-p)/sqrt(2*tr2*cpn_1))
}




## Binary Search Function
hBinarySearch <- function(l_pred_bound, u_pred_bound, pred_h, obs, p, dist = "Normal", s = 0.5, desired_FAP = 0.05, tol = 1e-04, niter = 100, sim = 500)
{
  l <- l_pred_bound
  u <- u_pred_bound
  m <- pred_h
  #m <- (l_pred_bound + u_pred_bound)/2
  D <- avg_FP_calc (obs, p, s, dist, m) 
  #D <- BinaryEntrp (xold, m) ## This is just the function we want to get to approximate root
  if (abs((D - desired_FAP)) < tol)
	return(list(h_value = m, FP_value = D))
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
    D <- avg_FP_calc (obs, p, s, dist, m)   ## This is just the function we want to get to approximate root
    #write.csv(D, "D.csv")
	if (abs((D-desired_FAP)) < tol)
      return(list(h_value = m, FP_value = D)) 
  }
  #print (abs(D-desired_FAP))
  #stop ("Exceeded")
}


avg_FP_calc <- function(obs, p, s, dist, h1)
{
  FP <- rep(0, 500)
  for (i in 1:500)
  {
    dataset = generate_data(obs, 0, p, s, dist, i)
    rmdp_res = RMDP(dataset$data)
    FP[i] = sum(rmdp_res > h1) / obs 
  }
  FP <- mean(FP)
}



## Example
dist = "Normal"
obs = 50
p = 50
corr = 0.25
lower = 0
upper = 10
est = 5

c = obs
#write.csv(c, "c.csv")

## Example 
result <- hBinarySearch(lower, upper, est, obs, p,  dist, s = corr)
result$h_value
result$FP_value
write.csv(result, "/home/martinwg/OCP/result.csv")
