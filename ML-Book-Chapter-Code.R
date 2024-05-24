library(e1071)
library(caTools)
library(survival)
library(dplyr)

smsurv <-function(Time,Status,X,beta,w,model){  # function to estimate baseline hazard using cox ph  
  death_point <- sort(unique(subset(Time, Status==1)))
  if(model=='ph') coxexp <- exp((beta)%*%t(X[,-1]))  
  lambda <- numeric()
  event <- numeric()
  for(i in 1: length(death_point)){
    event[i] <- sum(Status*as.numeric(Time==death_point[i]))
    if(model=='ph')  temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
    #if(model=='aft')  temp <- sum(as.numeric(Time>=death_point[i])*w)
    temp1 <- event[i]
    lambda[i] <- temp1/temp
  }
  HHazard <- numeric()
  for(i in 1:length(Time)){
    HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
    if(Time[i]>max(death_point))HHazard[i] <- Inf
    if(Time[i]<min(death_point))HHazard[i] <- 0
  }
  survival <- exp(-HHazard)
  list(survival=survival)
}

em.svm.RC <-function(Time,Status,X,Z,offsetvar,uncureprob,beta,emmax,eps,data) # EM function considering proportional hazards structure
{     
  w <- Status	# conditional probability of the ith individual remains uncured at the mth iteration. We use Status as initial value
  n <- length(Status)
  s <- smsurv(Time,Status,X,beta,w,model="ph")$survival # Estimate baseline hazards using COX PH
  
  convergence<- 1000;i <-1
  while (convergence > eps & i < emmax){ 
    #print(i)
    survival<-drop(s^(exp((beta)%*%t(X[,-1])))) # Estimate susceptible survival probabilities
    ## E step 
    w <- Status+(1-Status)*(uncureprob*survival)/((1-uncureprob)+uncureprob*survival) # conditional probability of the ith individual remains uncured at the mth iteration
    ## M step
    multipleuncureprob=matrix(1:5*n, nrow=n,ncol=5) #multiple imputation 
    for (j in 1:n){multipleuncureprob[j,]<-rbinom(5,size = 1,prob=w[j])}
    uncureprob1<-c(1,1)
    uncureprob2<-c(1,1)
    uncureprob3<-c(1,1)
    uncureprob4<-c(1,1)
    uncureprob5<-c(1,1)
    for (j in 1:n){uncureprob1[j]=multipleuncureprob[j,1]}
    for (j in 1:n){uncureprob2[j]=multipleuncureprob[j,2]}
    for (j in 1:n){uncureprob3[j]=multipleuncureprob[j,3]}
    for (j in 1:n){uncureprob4[j]=multipleuncureprob[j,4]}
    for (j in 1:n){uncureprob5[j]=multipleuncureprob[j,5]}
    
    for (j in 1:n){uncureprob1[j]=uncureprob1[j]*2-1}
    for (j in 1:n){uncureprob2[j]=uncureprob2[j]*2-1}
    for (j in 1:n){uncureprob3[j]=uncureprob3[j]*2-1}
    for (j in 1:n){uncureprob4[j]=uncureprob4[j]*2-1}
    for (j in 1:n){uncureprob5[j]=uncureprob5[j]*2-1}
    
    uncureprob1<-as.factor(uncureprob1)
    uncureprob2<-as.factor(uncureprob2)
    uncureprob3<-as.factor(uncureprob3)
    uncureprob4<-as.factor(uncureprob4)
    uncureprob5<-as.factor(uncureprob5)
    update_cureb<-c(1,1)
    # Estimate of incidence part using SVM
    obj<-tune(svm,uncureprob1~Z[,-1],data=data,kernel="radial", ranges=list(gamma = 2^(-1:1),cost=2^(2:4)),tunecontrol=tune.control(sampling = "fix"))
    bg<-obj$best.parameters[1]  
    bc<-obj$best.parameters[2]
    
    mod1<-svm(Z[,-1],uncureprob1,method = "C-classification",kernel="radial",gamma=bg[[1]],cost=bc[[1]], probability=TRUE)
    pred1<-predict(mod1,Z[,-1],probability = TRUE)
    proba1<-attr(pred1, "probabilities")
    update_cureb1<-c(1,1)
    for (z in 1:n){update_cureb1[z]<-proba1[z,colnames(proba1)==1]}
    uncureprob1<-as.numeric(as.character(uncureprob1))
    
    mod2<-svm(Z[,-1],uncureprob2,method = "C-classification",kernel="radial",gamma=bg[[1]],cost=bc[[1]], probability=TRUE)
    pred2<-predict(mod2,Z[,-1],probability = TRUE)
    proba2<-attr(pred2, "probabilities")
    update_cureb2<-c(1,1)
    for (z in 1:n){update_cureb2[z]<-proba2[z,colnames(proba2)==1]}
    uncureprob2<-as.numeric(as.character(uncureprob2))
    
    mod3<-svm(Z[,-1],uncureprob3,method = "C-classification",kernel="radial",gamma=bg[[1]],cost=bc[[1]], probability=TRUE)
    pred3<-predict(mod3,Z[,-1],probability = TRUE)
    proba3<-attr(pred3, "probabilities")
    update_cureb3<-c(1,1)
    for (z in 1:n){update_cureb3[z]<-proba3[z,colnames(proba3)==1]}
    uncureprob3<-as.numeric(as.character(uncureprob3))
    
    mod4<-svm(Z[,-1],uncureprob4,method = "C-classification",kernel="radial",gamma=bg[[1]],cost=bc[[1]], probability=TRUE)
    pred4<-predict(mod4,Z[,-1],probability = TRUE)
    proba4<-attr(pred4, "probabilities")
    update_cureb4<-c(1,1)
    for (z in 1:n){update_cureb4[z]<-proba4[z,colnames(proba4)==1]}
    uncureprob4<-as.numeric(as.character(uncureprob4))
    
    mod5<-svm(Z[,-1],uncureprob5,method = "C-classification",kernel="radial",gamma=bg[[1]],cost=bc[[1]], probability=TRUE)
    pred5<-predict(mod5,Z[,-1],probability = TRUE)
    proba5<-attr(pred5, "probabilities")
    update_cureb5<-c(1,1)
    for (z in 1:n){update_cureb5[z]<-proba5[z,colnames(proba5)==1]}
    uncureprob5<-as.numeric(as.character(uncureprob5))
    
    for (z in 1:n){update_cureb[z]<-(update_cureb1[z]+update_cureb2[z]+update_cureb3[z]+update_cureb4[z]+update_cureb5[z])/5}
    
    # Estimate of latency part using non-parametric coxph based on the Breslow method
    update_beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)),subset=w!=0, method="breslow")$coef
    update_s <-smsurv(Time,Status,X,beta,w,model="ph")$survival
    #convergence<-sum((beta-update_beta)^2)+sum((s-update_s)^2)+sum((uncureprob-update_cureb)^2)
    
    convergence<-sum(c(  mean(update_cureb)-mean(uncureprob),update_beta-beta,mean(update_s)-mean(s)      )^2   )
    
    uncureprob <- update_cureb
    beta <- update_beta 
    s<-update_s
    
    i <- i+1
    #print(i)
  }
  
  S1 = drop(s^(exp((beta)%*%t(X[,-1])))) # survival function of susceptible group
  Sp = (1-uncureprob)+(uncureprob*S1)  # survival function of overall population
  
  em.svm <- list(latencyfit=beta,Uncureprob=uncureprob,S0=s,S1=S1,Sp=Sp,tau=convergence, Mod1=mod1, Mod2=mod2, Mod3=mod3, Mod4=mod4, Mod5=mod5)
}

smcure <-function(formula,cureform,offset=NULL,data,na.action=na.omit,Var1=T,emmax=500,eps=1e-3,nboot=100){
  call <- match.call()
  cat("Program is running..be patient...")
  ## prepare data
  data <- na.action(data)
  n <- dim(data)[1]
  mf <- model.frame(formula,data)
  cvars <- all.vars(cureform)
  Z <- as.matrix(cbind(rep(1,n),data[,cvars]))
  colnames(Z) <- c("(Intercept)",cvars)
  if(!is.null(offset)) {
    offsetvar <- all.vars(offset)
    offsetvar<-data[,offsetvar]
  }else offsetvar <- NULL
  Y <- model.extract(mf,"response")
  X <- model.matrix(attr(mf,"terms"), mf)
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  Time <- Y[,1]
  Status <- Y[,2]
  bnm <- colnames(Z)
  nb <- ncol(Z)
  betanm <- colnames(X)[-1]
  nbeta <- ncol(X)-1
  ## initial values
  w <- Status
  nw<-c(1,1)
  for(i in 1: n) {nw[i]= w[i]*2-1}
  nw <- as.factor(nw)
  obj<-tune(svm,nw~Z[,-1], data=data,kernel="radial",ranges=list(gamma=2^(-1:1),cost=2^(2:4)),tunecontrol=tune.control(sampling = "fix"))
  bg<-obj$best.parameters[1]  
  bc<-obj$best.parameters[2]
  mod <- svm(Z[,-1],nw,method = "C-classification",kernel="radial",gamma=bg[[1]],cost=bc[[1]], probability=TRUE)
  pred <- predict(mod,Z[,-1], probability = TRUE)
  proba<-attr(pred, "probabilities")
  uncureprob<-c(1,1)
  for (z in 1:n){uncureprob[z]<-proba[z,colnames(proba)==1]}
  nw<-as.numeric(as.character(nw))
  
  beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)), subset=w!=0, method="breslow")$coef
  w <- as.factor(w)
  
  ## do EM algo
  emfit <- em.svm.RC(Time,Status,X,Z,offsetvar,uncureprob,beta,emmax,eps,data)#Time,Status,X,Z,offsetvar,uncureprob,beta,emmax,eps,data
  beta.est = emfit$latencyfit
  UN<-emfit$Uncureprob
  MOD1<-emfit$Mod1
  MOD2<-emfit$Mod2
  MOD3<-emfit$Mod3
  MOD4<-emfit$Mod4
  MOD5<-emfit$Mod5
  S1 = emfit$S1
  Sp = emfit$Sp
  S0 = emfit$S0
  
  if(Var1){ # Estimate standard deviation of latency parameters using bootstrapping
    uncure_boot<-matrix(rep(0,nboot*n), nrow=nboot)
    latency_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
    iter <- matrix(rep(0,nboot),ncol=1)
    tempdata <- cbind(Time,Status,X,Z)
    data1<-subset(tempdata,Status==1);data0<-subset(tempdata,Status==0)
    n1<-nrow(data1);n0<-nrow(data0)  
    i<-1
    while (i<=nboot){
      print(i)
      id1<-sample(1:n1,n1,replace=TRUE);id0<-sample(1:n0,n0,replace=TRUE)
      bootdata<-rbind(data1[id1,],data0[id0,])
      bootZ <- bootdata[,bnm]
      bootX <- as.matrix(cbind(rep(1,n),bootdata[,betanm]))
      
      bootfit <- em.svm.RC(bootdata[,1],bootdata[,2],bootX,bootZ,offsetvar,UN,beta.est,emmax,eps,data)
      latency_boot[i,] <- bootfit$latencyfit
      uncure_boot[i,] = bootfit$Uncureprob
      if (bootfit$tau<eps){
        i<-i+1}
    }#end of while
    
    latency_var <- apply(latency_boot, 2, var) 
    uncure_var = apply(uncure_boot,2,var)
    latency_sd <- sqrt(latency_var)
    uncure_sd = sqrt(uncure_var)
    lower_uncure = UN - (1.96*uncure_sd)
    upper_uncure = UN + (1.96*uncure_sd)
  }
  fit<-list()
  class(fit) <- c("smcure")
  
  fit$latency <- beta.est
  if(Var1){
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$uncure_sd = uncure_sd
    fit$latency_zvalue <- fit$latency/latency_sd
    fit$latency_pvalue <- (1-pnorm(abs(fit$latency_zvalue)))*2
    fit$lower_uncure = lower_uncure
    fit$upper_uncure = upper_uncure
  }
  cat(" done.\n")
  fit$call <- call
  fit$bnm <- bnm
  fit$betanm <- betanm
  #fit$s <- s
  fit$Time <- Time
  fit$UN<- UN
  fit$MOD1<- MOD1
  fit$Sp = Sp
  fit$S1 = S1
  fit$S0 = S0
  fit
}

