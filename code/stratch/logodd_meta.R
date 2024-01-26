
LogoddsMetaAnalysis <- function(logodds1,sigma1,logodds2,sigma2,logodds3,sigma3){
  sigma1.inv <- solve(sigma1)
  sigma2.inv <- solve(sigma2)
  sigma3.inv <- solve(sigma3)
  sigma.meta <- solve(sigma1.inv+sigma2.inv+sigma3.inv)
  logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2+
                                  sigma3.inv%*%logodds3)
  
  return(list(logodds.meta = logodds.meta,
              sigma.meta = sigma.meta))
  
}

