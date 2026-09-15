
rm(list=ls())

## The functions below cp.JG and RenewGLM.JG are copied and pasted from the R package
# RenewGLM and corrected so that they run without error for batch updates with batch sample
# size value 1

# For any R libraries please install them first as needed

cp.JG <- function(X, y, w) {
  X <- as.matrix(X)
  w <- as.numeric(w)
  crossprod(sqrt(w) * X)
}

invlink <- function(X, beta, type){
    #  if (length(beta)==1){
    #    eta <- X * beta
    #  }else{
    eta <- drop(X %*% beta)
    #  }
    if(type=="gaussian"){ out <- eta }
    if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta)) }
    if(type=="poisson") { out <- exp(eta) }
    out
  }

invlinkdiv <- function(X, beta, type){
    #  if (length(beta)==1){
    #    eta <- X * beta
    #  }else{
    eta <- drop(X %*% beta)
    #  }
    if(type=="gaussian"){ out <- 1 }
    if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta))^2 }
    if(type=="poisson") { out <- exp(eta) }
    out
  }

RenewGLM.JG <-function (X, y, type, betahat, infomats, intercept, s, phi) 
{
  tol <- 1e-06
  max_iter <- 100
  p <- dim(as.matrix(X))[2]
  betahat_old <- betahat
  W <- invlinkdiv(X, betahat_old, type = type)
  H <- cp.JG(X, y, W)
  #H <- matrix(H, nrow = ncol(infomats), ncol = ncol(infomats))
  U = chol(infomats + H)
  L = t(U)
  for (r in 1:max_iter) {
    g_0 = t(crossprod((y - invlink(X, betahat, type)), X))
    g_1 = -t(crossprod((betahat - betahat_old), infomats))
    g <- g_0 + g_1
    d_beta <- backsolve(U, forwardsolve(L, g))
    df_beta <- crossprod(g, d_beta)
    if (abs(df_beta) < tol) {
      break
    }
    else {
      betahat <- betahat + drop(d_beta)
    }
  }
  if (type == "gaussian") {
    s <- s + length(y)
    phi <- (s - length(y) - p)/(s - p) * phi + 1/(s - p) * 
      t(betahat_old) %*% infomats %*% (betahat_old - betahat) + 
      t(y) %*% (y - X %*% betahat)/(s - p)
  }
  W <- invlinkdiv(X, betahat, type = type)
  H_new <- cp.JG(X, y, W)
  #H_new <- matrix(H_new, nrow = ncol(infomats), ncol = ncol(infomats))
  infomats <- infomats + H_new
  if (type == "gaussian") {
    return(list(betahat, infomats, s, phi))
  }
  else {
    return(list(betahat, infomats))
  }
}

library(MASS)

library(devtools)

#install_github("luolsph/RenewGLM_pkg",force=TRUE)
#install_github("luolsph/RenewGLM_pkg")
library("RenewGLM")

# loads training data with first 4 rows corresponding to prior and sample size 1800
# not including prior rows
load(file="data.train.RData") 

n <- c(304, rep(1, 9))      # Short demonstration run with 10 batches
# n <- c(304, rep(1, 1500)) # Comment out to use 1501 batches for full analysis in the paper
N <- sum(n) # Cumulative training sample size
B <- length(n) # Total number of batches

# loads test data with sample size 1000
load(file="data.pred.RData")

names(data.train) # First column is intercept
dim(data.train)
#[1] 1804   20


names(data.pred) # First column is intercept
dim(data.pred)
#[1] 1000   20

source("allmodels.R")

p <- ncol(data.train)-1 # number of predictors including intercept
# 19

mat.all <- allmodels(p-1)
dim(mat.all)#262144    18

n.mc <- 1000 # Monte Carlo sample size 
n.top <- 5000 # Number of top models to keep track of
nsims <- 1 # One replicate

#Begin:  Online 2 ###
infomats.online2 <- vector("list",2^(p-1)) # stores infomats of diff dims for renewglm
beta.online2 <- vector("list",2^(p-1)) # stores betahat of diff dims for renewglm
incprob.online2 <- vector("list",nsims) # stores incprobs by sims and batches
prob.online2 <- vector("list",nsims) # stores model probs by sims and batches
betabma.online2 <- vector("list",nsims) # stores betabma by sims and batches
betasamp.online2 <- vector("list",nsims) # stores beta Monte Carlo 
online.time <- vector("list",nsims)


        for (i in 1:nsims)    
      {
        print(i)
        for (m in 1:(2^(p-1)))
        {
          dim <- sum(mat.all[m,])+1  
          beta.online2[[m]] <- rep(0,dim)
          infomats.online2[[m]] <- matrix(0,nrow=dim,ncol=dim)
        }
        online.time[[i]] <- rep(NA,B)
        betabma.online2[[i]] <- matrix(0,B,p)
        betasamp.online2[[i]] <- array(0,dim=c(B,n.mc,p))
        incprob.online2[[i]] <- matrix(0,B,p)
        prob.online2[[i]] <- matrix(0,B,2^(p-1))
        loglik.online2 <- rep(0,2^(p-1))
        sumy <- 0
        endind <- 0
        for (b in 1:1)
         online.time[[i]][b] <- system.time(
        {  
          print(c(i,b))
          betahat.m <- matrix(0,2^(p-1),p) # stores betahat for all models for current batch to make BMA estimate calculation easy
          bic.online2 <- rep(0,2^(p-1))
          #
          startind <- endind + 1 
          endind <- sum(n[1:b]) 
          selrows <- c(startind:endind)
          alldata <- data.train[selrows,]
          y <- matrix(alldata$is_correct,nrow=n[b],ncol=1)
          sumy <-  sumy + sum(y)
          sampleprop <- sumy/endind #sumy/(b*n)
          beta.online2[[1]] <- log(sampleprop/(1-sampleprop))
          infomats.online2[[1]] <- endind*sampleprop*(1-sampleprop)
          betahat.m[1,1] <- beta.online2[[1]]
          loglik.online2[1] <- loglik.online2[1] + sum(dbinom(x=y, size=1, prob=sampleprop, log=TRUE)) # null model log lik
          bic.online2[1] <- (-2)*loglik.online2[1] + log(endind)*1 #log(b*n)*1 # null model BIC so # of params=1
          #  
          for (m in 2:(2^(p-1)))
          {
            betahat  <- beta.online2[[m]]
            infomats <- infomats.online2[[m]]    
            which.vars <- which(c(1,mat.all[m,])!=0) # variable numbers in design matrix without intercept column
            n.pars<- length(which.vars) # dim including intercept
            X <-as.matrix(alldata[,1:p][,which.vars])
            summary<-RenewGLM_out(X,y,"binomial",betahat,infomats,intercept=FALSE)
            beta.online2[[m]] <- summary[[1]]
            infomats.online2[[m]] <- summary[[2]]
            online2.term <- 0.5*t(beta.online2[[m]]-betahat)%*%infomats%*%(beta.online2[[m]]-betahat)
            probs <- exp(X%*%beta.online2[[m]])/(1+exp(X%*%beta.online2[[m]]))
            loglik.online2[m] <- loglik.online2[m] + sum(dbinom(x=y, size=1, prob=probs, log=TRUE)) - online2.term
            #
            betahat.m[m,][which.vars] <- beta.online2[[m]]
            bic.online2[m] <- (-2)*loglik.online2[m] + log(endind)*n.pars 
          }
          #
          bic.online2 <- bic.online2-min(bic.online2)
          postprobs <- exp(-bic.online2/2)
          postprobs <- postprobs/sum(postprobs)
          #
          incprobs <- t(postprobs)%*%mat.all
          incprobs <- c(1,incprobs)
          incprob.online2[[i]][b,] <- incprobs
          prob.online2[[i]][b,] <- postprobs 
          betabma <- t(postprobs)%*%betahat.m
          betabma.online2[[i]][b,] <- betabma
          samp.mc <- sample(1:(2^(p-1)),prob = postprobs, replace = TRUE, size=n.mc)
          for (mc in 1:n.mc)
          {
            which.mc <- which(mat.all[samp.mc[mc],]==1)
            betasamp.online2[[i]][b,mc,c(1,which.mc+1)] <- mvrnorm(n=1,mu=beta.online2[[samp.mc[mc]]],Sigma=solve(infomats.online2[[samp.mc[mc]]]))
          }
          #
          top.which <- order(bic.online2,decreasing=FALSE)[1:n.top]
       })[3]
        for (b in 2:B)
          online.time[[i]][b] <- system.time(
        {  
          print(c(i,b))
          betahat.m <- matrix(0,2^(p-1),p) # stores betahat for all models for current batch to make BMA estimate calculation easy
          bic.online2 <- rep(NA,2^(p-1))
          #
          startind <- endind + 1 
          endind <- sum(n[1:b]) 
          selrows <- c(startind:endind)
          alldata <- data.train[selrows,]
          y <- matrix(alldata$is_correct,nrow=n[b],ncol=1)
          if (1%in%top.which)
              {
          sumy <-  sumy + sum(y)
          sampleprop <- sumy/endind #sumy/(b*n)
          beta.online2[[1]] <- log(sampleprop/(1-sampleprop))
          infomats.online2[[1]] <- endind*sampleprop*(1-sampleprop)
          betahat.m[1,1] <- beta.online2[[1]]
          loglik.online2[1] <- loglik.online2[1] + sum(dbinom(x=y, size=1, prob=sampleprop, log=TRUE)) # null model log lik
          bic.online2[1] <- (-2)*loglik.online2[1] + log(endind)*1 #log(b*n)*1 # null model BIC so # of params=1
          }
          #  
             for (m in top.which[top.which!=1])  
          {
            betahat  <- beta.online2[[m]]
            infomats <- infomats.online2[[m]]    
            which.vars <- which(c(1,mat.all[m,])!=0) # variable numbers in design matrix without intercept column
            n.pars<- length(which.vars) # dim including intercept
            X <-matrix(as.matrix(alldata[,1:p,drop=FALSE][,which.vars,drop=FALSE]),nrow=n[b],ncol=n.pars)
            summary<-RenewGLM.JG(X,y,"binomial",betahat,infomats,intercept=FALSE)
            beta.online2[[m]] <- summary[[1]]
            infomats.online2[[m]] <- summary[[2]]
            online2.term <- 0.5*t(beta.online2[[m]]-betahat)%*%infomats%*%(beta.online2[[m]]-betahat)
            probs <- exp(X%*%beta.online2[[m]])/(1+exp(X%*%beta.online2[[m]]))
            loglik.online2[m] <- loglik.online2[m] + sum(dbinom(x=y, size=1, prob=probs, log=TRUE)) - online2.term
            #
            betahat.m[m,][which.vars] <- beta.online2[[m]]
            bic.online2[m] <- (-2)*loglik.online2[m] + log(endind)*n.pars #log(b*n)*n.pars
          }
          #
          bic.online2 <- bic.online2-min(na.omit(bic.online2))
          postprobs <- exp(-bic.online2/2)
          postprobs[is.na(postprobs)] <- 0
          postprobs <- postprobs/sum(postprobs)
          #
          incprobs <- t(postprobs)%*%mat.all
          incprobs <- c(1,incprobs)
          incprob.online2[[i]][b,] <- incprobs
          prob.online2[[i]][b,] <- postprobs 
          betabma <- t(postprobs)%*%betahat.m
          betabma.online2[[i]][b,] <- betabma
          #
          samp.mc <- sample(1:(2^(p-1)),prob = postprobs, replace = TRUE, size=n.mc)
          for (mc in 1:n.mc)
          {
            which.mc <- which(mat.all[samp.mc[mc],]==1)
            betasamp.online2[[i]][b,mc,c(1,which.mc+1)] <- mvrnorm(n=1,mu=beta.online2[[samp.mc[mc]]],Sigma=solve(infomats.online2[[samp.mc[mc]]]))
          }
          })[3]
            } 
        

# Save output as RData files
  save(prob.online2,file="prob.online2.RData") 
  save(betabma.online2,file="betabma.online2.RData")
  save(incprob.online2,file="incprob.online2.RData")
  save(betasamp.online2,file="betasamp.online2.RData")
  save(online.time,file="online.time.RData")
 
 
  # Remove R objects from R workspace to free up RAM   
  rm(prob.online2)
  rm(betabma.online2)
  rm(incprob.online2)
  rm(betasamp.online2)
  rm(online.time)

 

  ## Begin: Offline ##
  infomats.off <- vector("list",2^(p-1)) # stores infomats of diff dims for renewglm
  betamle.off <- vector("list",2^(p-1)) # stores betahat of diff dims for renewglm
  betabma.off <- vector("list",length=nsims)
  incprob.off <- vector("list",length=nsims)
  prob.off <- vector("list",length=nsims)
  postprobs.off <- matrix(NA,2^(p-1),B)
  betasamp.off <- vector("list",nsims) # stores beta Monte Carlo 
  offline.time <- vector("list",nsims)
  

      for (i in 1:nsims)
      {
        offline.time[[i]] <- rep(NA,B)
        betabma.off[[i]] <- matrix(NA,B,p)
        betasamp.off[[i]] <- array(0,dim=c(B,n.mc,p))
        incprob.off[[i]] <- matrix(NA,B,p)
        prob.off[[i]] <- matrix(NA,B,2^(p-1))
        for (m in 1:(2^(p-1)))
        {
          dim <- sum(mat.all[m,])+1  
          betamle.off[[m]] <- rep(0,dim)
          infomats.off[[m]] <- matrix(0,nrow=dim,ncol=dim)
        }
        
          for (b in 1:1)
            offline.time[[i]][b] <- system.time(
        {
          print(c(i,b))       
          alldata <- data.train[1:sum(n[1:b]),]
          betahat.m <- matrix(0,2^(p-1),p)
          bic.v <- rep(0,2^(p-1))
          sampleprop <- sum(alldata$is_correct)/nrow(alldata)
          betahat.m[1,1] <- log(sampleprop/(1-sampleprop))
          betamle.off[[1]] <- betahat.m[1,1]
          infomats.off[[1]] <- sum(n[1:b])*sampleprop*(1-sampleprop)
          loglik.null <- sum(dbinom(x=alldata$is_correct, size=1, prob=sampleprop, log=TRUE))
          bic.v[1] <- (-2)*loglik.null + log(sum(n[1:b]))*1 #log(b*n)*1
          for (m in 2:(2^(p-1)))
          {
            which.vars <- which(c(1,mat.all[m,])!=0)
            n.pars<- length(which.vars)
            y<-alldata$is_correct
            X<-as.matrix(alldata[,which.vars])
            summary<-RenewGLM_out(X,y,"binomial",rep(0,n.pars),matrix(0,n.pars,n.pars),intercept=FALSE)
            betamle.off[[m]] <- summary[[1]]
            infomats.off[[m]] <- summary[[2]]
            probs <- exp(X%*%betamle.off[[m]])/(1+exp(X%*%betamle.off[[m]]))
            loglik <- sum(dbinom(x=y, size=1, prob=probs, log=TRUE)) 
            betahat.m[m,][which.vars] <- betamle.off[[m]]
            bic.v[m] <- (-2)*loglik + log(sum(n[1:b]))*n.pars #log(b*n)*n.pars
          }
          bic.v <- bic.v-min(bic.v)
          postprobs.v <- exp(-bic.v/2)
          postprobs.v <- postprobs.v/sum(postprobs.v)
          postprobs.off[,b] <- postprobs.v
          #
          incprobs <- t(postprobs.v)%*%mat.all
          incprobs <- c(1,incprobs)
          incprob.off[[i]][b,] <- incprobs
          prob.off[[i]][b,] <- postprobs.v 
          betabma <- t(postprobs.v)%*%betahat.m
          betabma.off[[i]][b,] <- betabma
          samp.mc <- sample(1:(2^(p-1)),prob = postprobs.v, replace = TRUE, size=n.mc)
          for (mc in 1:n.mc)
          {
            which.mc <- which(mat.all[samp.mc[mc],]==1)
            betasamp.off[[i]][b,mc,c(1,which.mc+1)] <- mvrnorm(n=1,mu=betamle.off[[samp.mc[mc]]],Sigma=solve(infomats.off[[samp.mc[mc]]]))
          }
          top.which.off <- order(bic.v,decreasing=FALSE)[1:n.top]
          })[3]
        for (b in 2:B)
          offline.time[[i]][b] <- system.time(
            {
          print(c(i,b))       
          alldata <- data.train[1:sum(n[1:b]),]
          betahat.m <- matrix(0,2^(p-1),p)
          bic.v <- rep(NA,2^(p-1))
          if (1%in%top.which.off)
          {
          sampleprop <- sum(alldata$is_correct)/nrow(alldata)
          betahat.m[1,1] <- log(sampleprop/(1-sampleprop))
          betamle.off[[1]] <- betahat.m[1,1]
          infomats.off[[1]] <- sum(n[1:b])*sampleprop*(1-sampleprop)
          loglik.null <- sum(dbinom(x=alldata$is_correct, size=1, prob=sampleprop, log=TRUE))
          bic.v[1] <- (-2)*loglik.null + log(sum(n[1:b]))*1 
          }
          for (m in top.which.off[top.which.off!=1])  
          {
            which.vars <- which(c(1,mat.all[m,])!=0)
            n.pars<- length(which.vars)
            y<-alldata$is_correct
            X<-as.matrix(alldata[,which.vars])
            summary<-RenewGLM_out(X,y,"binomial",rep(0,n.pars),matrix(0,n.pars,n.pars),intercept=FALSE)
            betamle.off[[m]] <- summary[[1]]
            infomats.off[[m]] <- summary[[2]]
            probs <- exp(X%*%betamle.off[[m]])/(1+exp(X%*%betamle.off[[m]]))
            loglik <- sum(dbinom(x=y, size=1, prob=probs, log=TRUE)) 
            betahat.m[m,][which.vars] <- betamle.off[[m]]
            bic.v[m] <- (-2)*loglik + log(sum(n[1:b]))*n.pars 
          }
          bic.v <- bic.v-min(na.omit(bic.v))
          postprobs.v <- exp(-bic.v/2)
          postprobs.v[is.na(postprobs.v)]=0
          postprobs.v <- postprobs.v/sum(postprobs.v)
          postprobs.off[,b] <- postprobs.v
          #
          incprobs <- t(postprobs.v)%*%mat.all
          incprobs <- c(1,incprobs)
          incprob.off[[i]][b,] <- incprobs
          prob.off[[i]][b,] <- postprobs.v 
          betabma <- t(postprobs.v)%*%betahat.m
          betabma.off[[i]][b,] <- betabma
          samp.mc <- sample(1:(2^(p-1)),prob = postprobs.v, replace = TRUE, size=n.mc)
          for (mc in 1:n.mc)
          {
            which.mc <- which(mat.all[samp.mc[mc],]==1)
            betasamp.off[[i]][b,mc,c(1,which.mc+1)] <- mvrnorm(n=1,mu=betamle.off[[samp.mc[mc]]],Sigma=solve(infomats.off[[samp.mc[mc]]]))
          }
            } 
          )[3]
      } 
  
  
  # Save output as RData files 
  save(offline.time,file="offline.time.RData")
  save(prob.off,file="prob.off.RData")  
  save(betabma.off,file="betabma.off.RData")
  save(incprob.off,file="incprob.off.RData")
  save(betasamp.off,file="betasamp.off.RData")
  
  # Remove R objects from R workspace to free up RAM  
  rm(offline.time)
  rm(prob.off)
  rm(betabma.off)
  rm(incprob.off)
  rm(betasamp.off)
  

  # Load saved output for analysis and plots
 
    load(file="prob.off.RData")  
    load(file="betabma.off.RData")
    load(file="incprob.off.RData")
    load(file="betasamp.off.RData")
    load(file="offline.time.RData")
 
    load(file="prob.online2.RData")  
    load(file="betabma.online2.RData")
    load(file="incprob.online2.RData")
    load(file="betasamp.online2.RData")
    load(file="online.time.RData")
   
    pdf("runningtime.pdf",height=5,width=5)
    plot(rep(2:B,2),c(online.time[[1]][-1],offline.time[[1]][-1]),type="n",xlab="Batches",ylab="Running time (seconds)",xlim=c(2,B))
    points(2:B,online.time[[1]][-1],col="blue",pch=1,cex=0.5)
    points(2:B,offline.time[[1]][-1],col="red",pch=4,cex=0.5)
    # legend position may need readjusting
    #legend("topleft",legend=c("Offline","Online"),col=c("red","blue"),pch=c(4,1),cex=1.2)
    dev.off()  
   
    
    n.xpred <- nrow(data.pred)
    probs.ci.online2 <- array(0,dim=c(B,n.xpred,2))
    probs.mean <- matrix(NA,B,n.xpred)
    
    # Online 
    
    for (b in 1:B)
    {
      if (b%%100==0) {print(b)}
      xbeta <- as.matrix(data.pred[,1:p])%*%t(betasamp.online2[[1]][b,,])
      xbeta[xbeta>=709] <- 709
      probs <- exp(xbeta)/(1+exp(xbeta))
      probs.ci.online2[b,,] <- t(apply(probs,1,function(x){quantile(x,prob=c(0.05,0.95))}))
      probs.mean[b,] <- apply(probs,1,mean)
    }
    
   
    probs.ci.online2.width <- matrix(0,n.xpred,B)
    
    for (b in 1:B)
    {
      probs.ci.online2.width[,b] <- probs.ci.online2[b,,2] - probs.ci.online2[b,,1]
    }
    
    probs.ci.online2.width <- data.frame(probs.ci.online2.width)
   
    #install.packages("pROC")
    library(pROC)
    
    auc.all.2 <- rep(NA,B)
    
    for (b in 1:B)
    {
      if (b%%100==0) {print(b)}
      auc.all.2[b] <- auc(roc(response=data.pred$is_correct,predictor=probs.mean[b,],levels=c(0,1),direction="<",smooth=TRUE,auc=TRUE))
    }
    
    
   
      
# Offline 

    n.xpred <- nrow(data.pred)
    probs.ci.off <- array(0,dim=c(B,n.xpred,2))
    probs.mean.off <- matrix(NA,B,n.xpred)
    
    
    
    for (b in 1:B)
    {
      if (b%%100==0) {print(b)}
      xbeta <- as.matrix(data.pred[,1:p])%*%t(betasamp.off[[1]][b,,])
      xbeta[xbeta>=709] <- 709
      probs.off <- exp(xbeta)/(1+exp(xbeta))
      probs.ci.off[b,,] <- t(apply(probs.off,1,function(x){quantile(x,prob=c(0.05,0.95))}))
      probs.mean.off[b,] <- apply(probs.off,1,mean)
    }
    
   
    probs.ci.off.width <- matrix(0,n.xpred,B)
    
    for (b in 1:B)
    {
      probs.ci.off.width[,b] <- probs.ci.off[b,,2] - probs.ci.off[b,,1]
    }
    
    probs.ci.off.width <- data.frame(probs.ci.off.width)
    
   
    #library(pROC)
    auc.all.off.2 <- rep(NA,B)
    
    for (b in 1:B)
    {
      if (b%%100==0) {print(b)}
      auc.all.off.2[b] <- auc(roc(response=data.pred$is_correct,predictor=probs.mean.off[b,],levels=c(0,1),direction="<",smooth=TRUE,auc=TRUE))
    }
    
      pdf("auc.pdf")    
      plot(rep(1:B,2),c(auc.all.off.2,auc.all.2),type="n",xlab="Batches",ylab="AUC")
      points(1:B,auc.all.off.2,pch=3,col="red",cex=0.4)
      points(1:B,auc.all.2,pch=1,col="blue",cex=0.4)
      #legend(x=100,y=0.770,legend=c("Offline","Online"),col=c("red","blue"),pch=c(3,1),cex=1.2)
      dev.off()   


    #install.packages("plotrix")  
    library(plotrix)

offset <- 20

# Note for this plot need to run for 1501 batches first
pdf("probint.pdf",height=8,width=6.5)
par(mfrow=c(4,1))
for (i in c(96, 140, 180, 337)) 
    {
 subset <- seq(1,B,by=250)     
 x.on <- c(1:B)[subset]+offset
 y.on <- probs.mean[,i][subset]
 lower.on <- probs.ci.online2[,i,1][subset]
 upper.on <-  probs.ci.online2[,i,2][subset]

 subset <- seq(1,B,by=250)     
 x.off <- c(1:B)[subset]-offset
 y.off <- probs.mean.off[,i][subset]
 lower.off <- probs.ci.off[,i,1][subset]
 upper.off <-  probs.ci.off[,i,2][subset]

 yrange <- range(c(lower.on, upper.on, lower.off, upper.off))

 plotCI(x=x.off, y=y.off, li = lower.off, ui = upper.off, col = "red", slty = 1, pch = 4, ylab = "90% credible intervals", xlab = "Batches", ylim = c(0,1), xlim = range(c(x.on, x.off, 1800)), 
        main=paste("Probability that Q" ,i,sep="", " in test sample is correctly answered (true y = ",data.pred$is_correct[i],")"),cex=0.4)
 plotCI(x=x.on, y=y.on, li = lower.on, ui = upper.on, slty = 2, pch = 1, col="blue", add = TRUE,cex=0.4)
 legend("bottomright", legend = c("Offline","Online"),  col = c("red","blue"),  pch = c(4,1),  lty = c(1, 2))
 }
dev.off()

subset <- seq(1,B,by=3)     

pdf("ciwidth.pdf",height=5,width=5)
plot((1:B)[subset],apply(probs.ci.online2.width,2,mean)[subset],
     ylim = range(c(apply(probs.ci.online2.width,2,mean)[subset],apply(probs.ci.off.width,2,mean)[subset])),
     col="blue",pch=1,cex=0.5,cex.lab=1,
     xlab="Batches",ylab="Width of 90% credible intervals")
points((1:B)[subset],apply(probs.ci.off.width,2,mean)[subset],col="red",pch=4,cex=0.5)
legend("topright",legend=c("Offline","Online"),col=c("red","blue"),pch=c(4,1),cex=1.2)
dev.off()
#

# Run the section below only after running the full 1501 batch analysis above. This section will
# return an error for the short 10 batch run as it requires output from batch 1501.
if (FALSE)
{
## Top models analysis begin
a.off <- cbind(mat.all[order(prob.off[[1]][B,],decreasing=TRUE)[1:5],],prob.off[[1]][B,][order(prob.off[[1]][B,],decreasing=TRUE)[1:5]])
colnames(a.off) <- c(names(data.train[,-1])[-p],"post.prob")
round(cbind(t(a.off),c(incprob.off[[1]][B,-1],NA)),2)
# Commented out results are for 1501 batches 
#                                [,1] [,2] [,3] [,4] [,5] [,6]
#problem_number                  0.00 0.00 0.00 0.00 0.00 0.03
#exercise_problem_repeat_session 0.00 0.00 0.00 0.00 0.00 0.03
#level_prev                      1.00 0.00 0.00 1.00 0.00 0.47
#points_prev                     0.00 0.00 0.00 0.00 1.00 0.05
#badges_cnt_prev                 0.00 0.00 0.00 0.00 0.00 0.04
#user_grade                      0.00 0.00 0.00 0.00 0.00 0.04
#has_teacher_cnt                 0.00 0.00 0.00 0.00 0.00 0.03
#is_self_coachTRUE               0.00 0.00 0.00 0.00 0.00 0.03
#has_student_cnt                 0.00 0.00 0.00 0.00 0.00 0.02
#belongs_to_class_cnt            0.00 0.00 0.00 0.00 0.00 0.03
#has_class_cnt                   0.00 0.00 0.00 0.00 0.00 0.05
#learning_stagejunior            0.00 0.00 1.00 1.00 0.00 0.26
#correct_rate                    1.00 1.00 1.00 1.00 1.00 1.00
#difficulty_lvl                  0.00 0.00 0.00 0.00 0.00 0.03
#login_dates_cnt                 0.00 0.00 0.00 0.00 0.00 0.04
#problem_correct_rate            1.00 1.00 1.00 1.00 1.00 1.00
#level2_id_count                 0.00 0.00 0.00 0.00 0.00 0.03
#level3_id_count                 0.00 0.00 0.00 0.00 0.00 0.03
#post.prob                       0.24 0.21 0.11 0.04 0.02   NA

# This code needs running for 1501 batches
gamma.1501.off <- betasamp.off[[1]][1501,,]
gamma.1501.off[gamma.1501.off!=0] <- 1
gamma.1501.off.string <- apply(gamma.1501.off, 1, paste0, collapse = "")
      top5.off.string <- apply(cbind(1,mat.all[order(prob.off[[1]][B,],decreasing=TRUE)[1:5],]), 1, paste0, collapse = "")
top5.off.string
#Commented out results are for 1501 batches 
#[1] "1001000000000100100" "1000000000000100100" "1000000000001100100" "1001000000001100100"
#[5] "1000100000000100100"
  
#which(top5.off.string[1]==gamma.1501.off.string)


      apply(betasamp.off[[1]][1501,which(top5.off.string[1]==gamma.1501.off.string),c(1,4,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})

      apply(betasamp.off[[1]][1501,which(top5.off.string[2]==gamma.1501.off.string),c(1,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})

      apply(betasamp.off[[1]][1501,which(top5.off.string[3]==gamma.1501.off.string),c(1,13,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})
      
      apply(betasamp.off[[1]][1501,which(top5.off.string[4]==gamma.1501.off.string),c(1,4,13,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})

       apply(betasamp.off[[1]][1501,which(top5.off.string[5]==gamma.1501.off.string),c(1,5,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))}) 

#       Commented out results are for 1501 batches 
prob.off[[1]][B,][order(prob.off[[1]][B,],decreasing=TRUE)[1:5]]
#[1] 0.24107976 0.21286666 0.11027076 0.04395149 0.01944506    

      
      
a.on <- cbind(mat.all[order(prob.online2[[1]][B,],decreasing=TRUE)[1:5],],prob.online2[[1]][B,][order(prob.online2[[1]][B,],decreasing=TRUE)[1:5]])
colnames(a.on) <- c(names(data.train[,-1])[-p],"post.prob")
round(cbind(t(a.on),c(incprob.online2[[1]][B,-1],NA)),2)
#Commented out results are for 1501 batches 
#                               [,1] [,2] [,3] [,4] [,5] [,6]
#
#problem_number                   0.0 0.00  0.0 0.00 0.00 0.03
#exercise_problem_repeat_session  0.0 0.00  0.0 0.00 0.00 0.03
#level_prev                       1.0 0.00  0.0 1.00 0.00 0.45
#points_prev                      0.0 0.00  0.0 0.00 0.00 0.06
#badges_cnt_prev                  0.0 0.00  0.0 0.00 1.00 0.08
#user_grade                       0.0 0.00  0.0 0.00 0.00 0.09
#has_teacher_cnt                  0.0 0.00  0.0 0.00 0.00 0.04
#is_self_coachTRUE                0.0 0.00  0.0 0.00 0.00 0.03
#has_student_cnt                  0.0 0.00  0.0 0.00 0.00 0.03
#belongs_to_class_cnt             0.0 0.00  0.0 0.00 0.00 0.03
#has_class_cnt                    0.0 0.00  0.0 0.00 0.00 0.05
#learning_stagejunior             0.0 0.00  1.0 1.00 0.00 0.30
#correct_rate                     1.0 1.00  1.0 1.00 1.00 1.00
#difficulty_lvl                   0.0 0.00  0.0 0.00 0.00 0.03
#login_dates_cnt                  0.0 0.00  0.0 0.00 0.00 0.04
#problem_correct_rate             1.0 1.00  1.0 1.00 1.00 1.00
#level2_id_count                  0.0 0.00  0.0 0.00 0.00 0.04
#level3_id_count                  0.0 0.00  0.0 0.00 0.00 0.04
#post.prob                        0.2 0.18  0.1 0.04 0.02   NA
gamma.1501.on <- betasamp.online2[[1]][1501,,]
gamma.1501.on[gamma.1501.on!=0] <- 1
gamma.1501.on.string <- apply(gamma.1501.on, 1, paste0, collapse = "")
      top5.on.string <- apply(cbind(1,mat.all[order(prob.online2[[1]][B,],decreasing=TRUE)[1:5],]), 1, paste0, collapse = "")
top5.on.string
#Commented out results are for 1501 batches 
#[1] "1001000000000100100" "1000000000000100100" "1000000000001100100" "1001000000001100100"
#[5] "1000010000000100100"
  
#which(top5.on.string[1]==gamma.1501.on.string)


      apply(betasamp.online2[[1]][1501,which(top5.on.string[1]==gamma.1501.on.string),c(1,4,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})

      apply(betasamp.online2[[1]][1501,which(top5.on.string[2]==gamma.1501.on.string),c(1,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})

      apply(betasamp.online2[[1]][1501,which(top5.on.string[3]==gamma.1501.on.string),c(1,13,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})
      
      apply(betasamp.online2[[1]][1501,which(top5.on.string[4]==gamma.1501.on.string),c(1,4,13,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))})

       apply(betasamp.online2[[1]][1501,which(top5.on.string[5]==gamma.1501.on.string),c(1,6,14,17)],2,function(x){quantile(x,prob=c(0.05,0.5,0.95))}) 

#       Commented out results are for 1501 batches 
prob.online2[[1]][B,][order(prob.online2[[1]][B,],decreasing=TRUE)[1:5]]
#[1] 0.20299681 0.17746490 0.10139686 0.04082954 0.02329259
}
 
