
generateData<-function(n,p,alpha){
  X = matrix(rnorm(n * p), n, p)
  nu<-rnorm(n)
  t<-pmax(X[,1]*alpha+nu,0)
  for(i in 1:p){
    for(j in 1:p){
      X<-cbind(X,X[,i]*X[,j])  
    }
  }
  X<-cbind(X,t)
  X<-cbind(X,t^2)
  X<-cbind(X,t^3)
  return(X)
}

generateCoefficients<-function(treatStrength, sparsity,X){
  coef<-dim(X)[2]
  beta<-rnorm(coef)
  if(sparsity>0){
    indexzero<-sample(1:(coef-3))[1:floor(sparsity*(coef-3))]
    beta[indexzero]=0
  }
  beta[coef-2]<-treatStrength
  beta[coef-1]<-1/2*treatStrength
  beta[coef]<-1/8*treatStrength
  return(beta)
}

generateOutcome<-function(X,beta){
  n<-dim(X)[1]
  eps<-rnorm(n)
  Y= X%*%beta + eps
  return(Y)
}

simulateExc<-function(n,p,sparsity,treatStrength,alpha,simrep){
  ResultErrors<-data.frame(LINerror=rep(0,simrep),GBMerror=rep(0,simrep),NNETerror=rep(0,simrep))
  
  for(rep in 1:simrep){
    X<-generateData(n,p,alpha = alpha)
    beta<-generateCoefficients(treatStrength,sparsity,X)
    Y<-generateOutcome(X,beta)
    t=X[,dim(X)[2]-2]
    tmax<-as.numeric(quantile(t,0.9))
    trueResponse<-seq(0,tmax,0.05)*trtstr+seq(0,tmax,0.05)^2*trtstr*1/2+seq(0,tmax,0.05)^3*trtstr*(1/8)
    linmodel<-lm(t~ -1+X[,1:p])
    coefsFromModel<-linmodel$coefficients
    coefsFromModel<-as.vector(coefsFromModel)
    meanVector<-X[,1:p]%*%coefsFromModel
    propensities<-dnorm(t-meanVector)
    Estindex<-(propensities>0.025&propensities<0.975)
    EstT<-t[Estindex]
    EstX<-X[Estindex,]
    EstY<-Y[Estindex]
    EstP<-propensities[Estindex]
    #HighOrderLinear
    linmod<-lm(EstY~EstT+I(EstT^2)+EstP+I(EstP^2)+EstT*EstP)
    EstResponse<-data.frame(t=0,est=0)
    for(treat in seq(0,tmax,0.05)){
      EstP<-dnorm(rep(treat,length(Y))-meanVector)
      Estindex<-(EstP>0.05&EstP<0.95)
      EstP<-EstP[Estindex]
      EstT<-rep(treat,length(EstP))
      EstTsq<-rep(treat^2,length(EstP))
      EstPsq<-EstP^2
      estResp<-mean(predict(linmod,newdata=data.frame(EstT=EstT,EstP=EstP)))
      EstResponse<-rbind(EstResponse,c(treat,estResp))
    }
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
      EstP<-dnorm(rep(treat,length(Y))-meanVector)
      Estindex<-(EstP>0.05&EstP<0.95)
      EstP<-EstP[Estindex]
      EstT<-rep(treat,length(EstP))
      EstFrame<-data.frame(EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)
      estResp<-mean(predict(modelGBM,newdata=EstFrame,num.trees=besttrees))
      EstResponseGBM<-rbind(EstResponseGBM,c(treat,estResp))
    }
    ResultErrors[rep,]$GBMerror<-mean((trueResponse-EstResponseGBM$est)^2)
    #NEURAL NETS
    propensities<-dnorm(t-meanVector)
    Estindex<-(propensities>0.05&propensities<0.95)
    EstT<-t[Estindex]
    EstX<-X[Estindex,]
    EstY<-Y[Estindex]
    EstP<-propensities[Estindex]
    EstFrame<-data.frame(EstY=EstY,EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)
    modelNNET<-neuralnet(EstY~EstT+EstP+EstTsq+EstPsq+Intr, data= EstFrame ,rep=3,stepmax=500000,threshold=30,act.fct="tanh",hidden = c(4,2),lifesign = "full",linear.output=TRUE)
    EstResponseNNET<-data.frame(t=0,est=0)
    for(treat in seq(0,tmax,0.05)){
      EstP<-dnorm(rep(treat,length(Y))-meanVector)
      Estindex<-(EstP>0.05&EstP<0.95)
      EstP<-EstP[Estindex]
      EstY<-Y[Estindex]
      EstT<-rep(treat,length(EstP))
      EstFrame<-data.frame(EstY=EstY,EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)
      estResp<-mean(compute(modelNNET,covariate=EstFrame[,2:6])$net.result)
      EstResponseNNET<-rbind(EstResponseNNET,c(treat,estResp))
    }
    ResultErrors[rep,]$NNETerror<-mean((trueResponse-EstResponseNNET$est)^2)
  }
  return(ResultErrors)
  
}