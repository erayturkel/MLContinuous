
generateData<-function(n,p,alpha,treatProb){
  X = matrix(rnorm(n * p), n, p)
  nu<-rnorm(n)
  t<-pmax(X[,1]*alpha+nu,0)
  for(i in 1:p){
    for(j in 1:p){
      X<-cbind(X,X[,i]*X[,j])  
    }
  }
  X<-cbind(X,t)
  return(X)
}

generateCoefficients<-function(treatStrength, sparsity,X){
  coef<-dim(X)[2]
  beta<-rnorm(coef)
  if(sparsity>0){
    indexzero<-sample(1:(coef-1))[1:floor(sparsity*(coef-1))]
    beta[indexzero]=0
  }
  beta[coef]<-treatStrength
  return(beta)
}

generateOutcome<-function(X,beta){
  n<-dim(X)[1]
  eps<-rnorm(n)
  Y= X%*%beta + eps
  return(Y)
}
n=5000
alpha=0.8
treatStrength=4
p=20
trtstr=3
sparsity=0.5
X<-generateData(n,p,alpha = alpha)
beta<-generateCoefficients(treatStrength,sparsity,X)
trtstr<-treatStrength
Y<-generateOutcome(X,beta)
t=X[,dim(X)[2]]
tmax<-as.numeric(quantile(t,0.9))
trueResponse<-seq(0,tmax,0.05)*trtstr
linmodel<-truncreg(t~ -1+X[,1:p])
coefsFromModel<-linmodel$coefficients
sigma<-coefsFromModel[p+1]
coefsFromModel<-coefsFromModel[1:p]
coefsFromModel<-as.vector(coefsFromModel)
sigma<-as.numeric(sigma)
meanVector<-X[,1:p]%*%coefsFromModel
propensities<-rep(0,length(t))
for(i in 1:length(t)){
  if(t[i]!=0){
    mu<-meanVector[i]  
    nom<-dnorm((t[i]-mu)/sigma)
    denom<- sigma
    propensities[i]<- nom/denom
  }else{
    mu<-meanVector[i]  
    propensities[i]<-pnorm(-mu/sigma)
  }
}
Estindex<-(propensities>0.025&propensities<0.975)
EstT<-t[Estindex]
EstX<-X[Estindex,]
EstY<-Y[Estindex]
EstP<-propensities[Estindex]
#HighOrderLinear
linmod<-lm(EstY~EstT+I(EstT^2)+EstP+I(EstP^2)+EstT*EstP)
EstResponse<-data.frame(t=0,est=0)
for(treat in seq(0,tmax,0.05)){
  if(treat!=0){
    EstP<-dnorm((rep(treat,length(Y))-meanVector)/sigma)/sigma
  }else{
    EstP<-pnorm(-meanVector/sigma)
  }
  Estindex<-(EstP>0.05&EstP<0.95)
  EstP<-EstP[Estindex]
  EstT<-rep(treat,length(EstP))
  EstTsq<-rep(treat^2,length(EstP))
  EstPsq<-EstP^2
  estResp<-mean(predict(linmod,newdata=data.frame(EstT=EstT,EstP=EstP)))
  EstResponse<-rbind(EstResponse,c(treat,estResp))
}
EstResponse<-EstResponse[-1,]
ResultErrors[rep,]$LINerror<-mean((trueResponse-EstResponse$est)^2)


##Gradient Boosting
library(gbm)
propensities<-dnorm(t-meanVector)
Estindex<-(propensities>0.05&propensities<0.95)
EstT<-t[Estindex]
EstX<-X[Estindex,]
EstY<-Y[Estindex]
EstP<-propensities[Estindex]
EstFrame<-data.frame(EstY=EstY,EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)
modelGBM<-gbm(EstY~.,data=EstFrame,distribution = "gaussian",n.trees = 5000,interaction.depth = 3,cv.folds = 5)
besttrees<-gbm.perf(modelGBM,method = "cv")
EstResponseGBM<-data.frame(t=0,est=0)
for(treat in seq(0,tmax,0.05)){
  if(treat!=0){
    EstP<-dnorm((rep(treat,length(Y))-meanVector)/sigma)/sigma
  }else{
    EstP<-pnorm(-meanVector/sigma)
  }
  Estindex<-(EstP>0.05&EstP<0.95)
  EstP<-EstP[Estindex]
  EstT<-rep(treat,length(EstP))
  EstFrame<-data.frame(EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)
  estResp<-mean(predict(modelGBM,newdata=EstFrame,num.trees=besttrees))
  EstResponseGBM<-rbind(EstResponseGBM,c(treat,estResp))
}
EstResponseGBM<-EstResponseGBM[-1,]
ResultErrors[rep,]$GBMerror<-mean((trueResponse-EstResponseGBM$est)^2)
#NEURAL NETS
propensities<-dnorm(t-meanVector)
Estindex<-(propensities>0.05&propensities<0.95)
EstT<-t[Estindex]
EstX<-X[Estindex,]
EstY<-Y[Estindex]
EstP<-propensities[Estindex]
EstFrame<-data.frame(EstY=EstY,EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)
modelNNET<-neuralnet(EstY~EstT+EstP+EstTsq+EstPsq+Intr, data= EstFrame ,rep=5,stepmax=100000,threshold=5000,act.fct="logistic",hidden = c(4),lifesign = "full",linear.output=TRUE)
EstResponseNNET<-data.frame(t=0,est=0)
for(treat in seq(0,tmax,0.05)){
  if(treat!=0){
    EstP<-dnorm((rep(treat,length(Y))-meanVector)/sigma)/sigma
  }else{
    EstP<-pnorm(-meanVector/sigma)
  }
  Estindex<-(EstP>0.05&EstP<0.95)
  EstP<-EstP[Estindex]
  EstY<-Y[Estindex]
  EstT<-rep(treat,length(EstP))
  EstFrame<-data.frame(EstY=EstY,EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)
  estResp<-mean(compute(modelNNET,covariate=EstFrame[,2:6])$net.result)
  EstResponseNNET<-rbind(EstResponseNNET,c(treat,estResp))
}
EstResponseNNET<-EstResponseNNET[-1,]
ResultErrors[rep,]$NNETerror<-mean((trueResponse-EstResponseNNET$est)^2)


ggplot()+geom_point(aes(x=seq(0,tmax,0.05),y=seq(0,tmax,0.05)*trtstr,color="True Response"))+geom_point(aes(x=EstResponseGBM$t,y=EstResponseGBM$est,color="GBM"))+geom_point(aes(x=EstResponse$t,y=EstResponse$est,color="Linear Regression"))+geom_point(aes(x=EstResponseNNET$t,y=EstResponseNNET$est,color="NN"))




#NAIVE PREDICTION

EstT<-t
EstX<-X[,1:p]
EstY<-Y
EstFrame<-data.frame(EstY=Y,EstT=t,EstX=EstX)
modelNNET<-neuralnet(paste("EstY","~",paste(colnames(EstFrame)[2:22],collapse="+")), data= EstFrame ,rep=5,stepmax=100000,threshold=5000,act.fct="logistic",hidden = c(4),lifesign = "full",linear.output=TRUE)
EstResponseNaive<-data.frame(t=0,est=0)
for(treat in seq(0,tmax,0.05)){
  EstY<-Y
  EstT<-rep(treat,length(EstY))
  EstFrame<-data.frame(EstY=Y,EstT=EstT,EstX=EstX)
  estResp<-mean(compute(modelNNET,covariate=EstFrame[,2:22])$net.result)
  EstResponseNaive<-rbind(EstResponseNaive,c(treat,estResp))
}
EstResponseNaive<-EstResponseNaive[-1,]


ggplot()+geom_point(aes(x=seq(0,tmax,0.05),y=seq(0,tmax,0.05)*trtstr,color="True Response"))+geom_point(aes(x=EstResponseGBM$t,y=EstResponseGBM$est,color="GBM"))+geom_point(aes(x=EstResponse$t,y=EstResponse$est,color="Linear Regression"))+geom_point(aes(x=EstResponseNNET$t,y=EstResponseNNET$est,color="NN"))+geom_point(aes(x=EstResponseNaive$t,y=EstResponseNaive$est,color="Naive Predictive NN"))



