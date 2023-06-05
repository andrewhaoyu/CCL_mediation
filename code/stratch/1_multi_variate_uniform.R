library(MultiRNG)

cmat<-matrix(c(1,0.6691746,0.6494822,0.6691746,1,0.9127285,0.6494822,0.9127285,1), nrow=3, ncol=3)
mydata=draw.d.variate.uniform(no.row=3.53e07,d=3,cov.mat=cmat)
apply(mydata,2,mean)-rep(0.5,3)
cor(mydata)-cmat

sum(mydata[,3]<=1E-06)/3.53e07
library(ACAT)
test = ACAT(t(mydata))
head(test)
ACAT(mydata[1,])
result = cbind(mydata,test)
save(result,)