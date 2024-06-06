llk <- function(X,mu,type,param)
{
  p = nrow(X)
  n = ncol(X)
  
  l = 0
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
      l = l - n*log(mean((X[j,]-mu[j,])^2))/2
    else if ( type[j] == 1 )
      l = l + sum(dbinom(X[j,],param[j],1/(1+exp(-mu[j,])),TRUE))
    else if ( type[j] == 2 )
      l = l + sum(dnbinom(X[j,],param[j],1/(1+exp(mu[j,])),log=TRUE))
    else if ( type[j] == 3 )
      l = l + sum(dpois(X[j,],exp(mu[j,]),TRUE))
  }
  return(l)
}


BIC_tuning <- function(dat,res,coff,type, param){
  n = dim(res$mu)[2]
  p = dim(res$mu)[1]
 # new w, new z and new m   
  #mu =  (abs(res$m)>0.05)*res$m +((abs(res$W)>0.05)*res$W)%*% ((abs(res$Z)>0.05)*res$Z)  
  BIC_res = -2*llk(dat,as.matrix(res$mu),type,param) + (sum(abs(res$W)>coff)+sum(abs(res$Z)>coff))*log(n*p)
  #BIC_res = -2*llk(dat,mu,type,param) + (sum(abs(res$W)>coff)+sum(abs(res$Z)>coff))*log(n*p)
  return(BIC_res)
}


