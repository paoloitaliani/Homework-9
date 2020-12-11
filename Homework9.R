########EXERCISE 1

library(leaps)
library(glmnet)
library(pls)
library(knitr)

covid.weekly=read.table("covidweeklygrowths.txt",header = T)
covid.weekly <- covid.weekly[,5:170]
P <- 165
n <- 175
covid.weekly$inc.var <- covid.weekly$longestinc <- covid.weekly$inc.longest <-
  covid.weekly$longestdec <- covid.weekly$dec.longest <-
  covid.weekly$nphases <- covid.weekly$period1 <- covid.weekly$period2 <-
  covid.weekly$period3 <- covid.weekly$period4 <- covid.weekly$period5 <-
  covid.weekly$nmax10 <- rep(NA,n)
daysdata <- as.matrix(covid.weekly[,2:(P+1)])


q90 <- quantile(as.vector(daysdata),probs=0.9) # For finding highest 10%
increases <- matrix(NA,nrow=n,ncol=P)
for(i in 1:n){
  covid.weekly$inc.var[i] <- var(daysdata[i,]) # Variance for country i
  covid.weekly$nmax10[i] <- sum(daysdata[i,]>q90) # Number of values in highest 10%
  for(j in 2:P) # Do values increase from day j to j+1?
    increases[i,j] <- as.integer(daysdata[i,j]>daysdata[i,j-1])
  runinfo <- rle(increases[i,2:P])
  # rle: "Compute the lengths and values of runs of equal values"
  runloc <- c(1,cumsum(runinfo$lengths)) # locations of run starts
  covid.weekly$nphases[i] <- sum(runinfo$values==1) # Number of runs of increases
  covid.weekly$longestinc[i] <- max(runinfo$lengths[runinfo$values==1]) # longest run
  covid.weekly$longestdec[i] <- max(runinfo$lengths[runinfo$values==0])
  inclengths <- declengths <- runinfo$lengths
  inclengths[runinfo$values==0] <- 0
  declengths[runinfo$values==1] <- 0
  incloc <- 1+which.max(inclengths) # incloc element indicating increase run end
  decloc <- 1+which.max(declengths) # same for decrease run
  covid.weekly$inc.longest[i] <-
    mean(daysdata[i,runloc[incloc-1]:runloc[incloc]]) # mean value of increase run
  covid.weekly$dec.longest[i] <-
    mean(daysdata[i,runloc[decloc-1]:runloc[decloc]]) # same for decrease run
  covid.weekly$period1[i] <- sum(daysdata[i,1:33]) # Summing over five periods
  covid.weekly$period2[i] <- sum(daysdata[i,34:66])
  covid.weekly$period3[i] <- sum(daysdata[i,67:99])
  covid.weekly$period4[i] <- sum(daysdata[i,100:132])
  covid.weekly$period5[i] <- sum(daysdata[i,133:165])
}

cvdata<-covid.weekly
cvdata=cvdata[,c(1,167:178)]
varnames <- names(cvdata)[-1]
response <- names(cvdata)[1]

cvdata<-cvdata[sample(nrow(cvdata)),]
indexes <- cut(seq(1,nrow(cvdata)),breaks=10,labels=FALSE)
fold=list()
for(i in 1:10){
  fold[[i]] <- which(indexes==i,arr.ind=TRUE)
  
}

nmodels <- 12
P=12
response <- 1
names(cvdata)[response] <- "y"
cvdata$y=log(cvdata$y)
full.model=as.formula(paste("y","~", paste(varnames,collapse=" + ")))
n=175
yhat <- besti_sub <- numeric(0)
fmbesti_sub <- mbesti_sub <- bestmodel <- modforcvi <- list()
yhati <- sqrlossi <- matrix(NA,nrow=nmodels,ncol=n)
foldnumber=10
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  yhati <- matrix(NA,nrow=nmodels,ncol=n)
  modforcv <- list()
  for (j in (1:foldnumber)[-i]){
    cat("Inner fold ",j,"\n")
    bestsubi <- regsubsets(full.model,data=cvdata[-c(fold[[i]],fold[[j]]),],
                           nvmax=nmodels,method="forward") # leaving folds i and j out
    sbesti_sub <- summary(bestsubi)
    for (k in 1:nmodels){
      modforcv[[k]] <- as.formula(paste("y~", paste(varnames[sbesti_sub$which[k,2:(P+1)]],
                                                    collapse=" + "))) # extract best model
      fmi <- lm(modforcv[[k]], data=cvdata[-c(fold[[i]],fold[[j]]),]) # fit it
      yhati[k,fold[[j]]] <- predict(fmi,cvdata[fold[[j]],])
      # predict fold j for finding best k
    } # end for k (models)
  } # end for j (inner loop) outer loop still running
  for (k in 1:nmodels)
    sqrlossi[k,i] <- sqrt(mean((yhati[k,]-cvdata$y)^2,na.rm=TRUE))
  besti_sub[i] <- which.min(sqrlossi[,i]) # Best model chosen without fold i
  bestmodel[[i]] <- regsubsets(full.model,data=cvdata[-fold[[i]],],
                               nvmax=besti_sub[i],method="forward") # run forward selection on data without fold i
  sbesti_sub <- summary(bestmodel[[i]])
  modforcvi[[i]] <- as.formula(paste("y~", paste(varnames[sbesti_sub$which[besti_sub[[i]],2:(P+1)]],
                                                 collapse=" + ")))
  # Extract best model as found in inner loop
  fmbesti_sub[[i]] <- lm(modforcvi[[i]], data=cvdata[-fold[[i]],]) # Fit this
  yhat[fold[[i]]] <- predict(fmbesti_sub[[i]],cvdata[fold[[i]],])
  # Predict fold i data with best model selected without fold i.
}
sqrlossbest4 <- sqrt(mean((yhat-cvdata$y)^2))
besti_sub

#########LASSO
# Use existing folds for 10-fold for outer loop, to make result
# comparable with forward:
foldidorig <- rep(NA,n)
for(i in 1:foldnumber){
  foldidorig[fold[[i]]] <- i # cv.glmnet allows fold indicating vector as input.
}

# The outer loop needs to be the same, but not the inner loop!
response <- 1
yhat <- bestp <- bestlambda<-besti_LASSO <- bestcvm <- numeric(0)
lmodel <- list()
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  lmodel[[i]] <- glmnet(x=as.matrix(cvdata[varnames])[-fold[[i]],],y=cvdata$y[-fold[[i]]])
  # Need this in oder to unify lambdas for all fitted models.
  lambdai <- lmodel[[i]]$lambda
  foldidx <- foldidorig # Preparation of removing fold i
  foldidx[foldidx>=i] <- foldidx[foldidx>=i]-1 # Update fold numbers if larger than i
  foldidi <- foldidx[-fold[[i]]] # remove fold i from foldid vector
  foldilasso <- cv.glmnet(x=as.matrix(cvdata[varnames])[-fold[[i]],],
                          y=cvdata$y[-fold[[i]]],foldid=foldidi,lambda=lambdai)
  besti_LASSO[i] <- which.min(foldilasso$cvm)
  bestp[i] <- foldilasso$nzero[besti_LASSO[i]] # Number of nonzero variables
  bestlambda[i] <- foldilasso$lambda[besti_LASSO[i]] # lambda of best model
  bestcvm[i] <- foldilasso$cvm[besti_LASSO[i]] # cross validation loss of best model
  yhat[fold[[i]]] <- predict(foldilasso$glmnet.fit,
                             as.matrix(cvdata[varnames][fold[[i]],]),s=bestlambda[i])
}

sqrlosslasso4 <- sqrt(mean((yhat-cvdata$y)^2))
bestp
##############PCR
response <- 1
yhat<- besti_PCR<-numeric(0)
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  pcrLBP <- pcr(full.model, ncomp=12,data=cvdata[-fold[[i]],],scale=TRUE,validation="LOO")
  besti_PCR[i] <- which.min(as.data.frame(RMSEP(pcrLBP)$val)[1,])-1
  
  yhat[fold[[i]]] <-predict(pcrLBP,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
}

sqrlossPCR4 <- sqrt(mean((yhat-cvdata$y)^2))
besti_PCR


#############PLS
yhat<- besti_PLS<-numeric(0)
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  pcrLBP <- plsr(full.model,ncomp=12,data=cvdata[-fold[[i]],],scale=TRUE,validation="LOO")
  besti_PLS[i] <- which.min(as.data.frame(RMSEP(pcrLBP)$val)[1,])-1
  
  yhat[fold[[i]]] <-predict(pcrLBP,cvdata[fold[[i]],] , ncomp = besti_PLS[i])
}

sqrlossPLS <- sqrt(mean((yhat-cvdata$y)^2))
besti_PLS







#########mean
yhat=numeric(0)
for (i in 1:foldnumber){
  model=lm(y~1,data=cvdata[-fold[[i]],],)
  yhat[fold[[i]]] <-predict(model,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
  
}
sqrlossmean4<- sqrt(mean((yhat-cvdata$y)^2))

#########à fullmodel
yhat=numeric(0)
for (i in 1:foldnumber){
  model=lm(y~.,data=cvdata[-fold[[i]],],)
  yhat[fold[[i]]] <-predict(model,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
  
}
sqrlossfull.model<- sqrt(mean((yhat-cvdata$y)^2))



##############Exercise 4
crime <- read.csv("london-borough-profilesf.csv")


i=1
c=c()
pb <- txtProgressBar(min = 0, max = 85, style = 3)
while(i < 85){
  if (class(crime[,i])=="factor"){
    c=append(c,i)
    i=i+1
  }
  else{
    i=i+1
  }
  setTxtProgressBar(pb, i)
}

crime_new=crime[,-c]
crime_new=crime_new[,-c(31,32)]
names(crime_new)
crime_new=crime_new[,-c(66,65,60,21,22,27,31,33,47,50,52)]
names(crime_new)
cor=cor(crime_new)
cor=as.data.frame(cor)

varlist=list()
corlist=list()
for(i in 1:71){
  varlist[[i]]=row.names(cor[cor[i,]>0.9,])
  corlist[[i]]=cor[cor[i,]>0.95,][,i]
}
