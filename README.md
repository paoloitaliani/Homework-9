

# Exercise 1


```r
setwd("~/Desktop")
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

set.seed(1234)
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

response <- 1
yhat<- besti_PCR<-numeric(0)
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  pcrLBP <- pcr(full.model, ncomp=12,data=cvdata[-fold[[i]],],scale=TRUE,validation="LOO")
  besti_PCR[i] <- which.min(as.data.frame(RMSEP(pcrLBP)$val)[1,])-1
  
  yhat[fold[[i]]] <-predict(pcrLBP,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
}

sqrlossPCR4 <- sqrt(mean((yhat-cvdata$y)^2))

yhat<- besti_PLS<-numeric(0)
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  pcrLBP <- plsr(full.model,ncomp=12,data=cvdata[-fold[[i]],],scale=TRUE,validation="LOO")
  besti_PLS[i] <- which.min(as.data.frame(RMSEP(pcrLBP)$val)[1,])-1
  
  yhat[fold[[i]]] <-predict(pcrLBP,cvdata[fold[[i]],] , ncomp = besti_PLS[i])
}

sqrlossPLS <- sqrt(mean((yhat-cvdata$y)^2))


######mean
yhat=numeric(0)
for (i in 1:foldnumber){
  model=lm(y~1,data=cvdata[-fold[[i]],],)
  yhat[fold[[i]]] <-predict(model,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
  
}
sqrlossmean4<- sqrt(mean((yhat-cvdata$y)^2))

######### fullmodel
yhat=numeric(0)
for (i in 1:foldnumber){
  model=lm(y~.,data=cvdata[-fold[[i]],],)
  yhat[fold[[i]]] <-predict(model,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
  
}
sqrlossfull.model<- sqrt(mean((yhat-cvdata$y)^2))
```


```r
trial <-matrix(c(sqrlossbest4,sqrlosslasso4,sqrlossPCR4,sqrlossPLS,sqrlossfull.model, sqrlossmean4) ,ncol=6,nrow=1)
colnames(trial) <-c("regsubsets","LASSO","PCR","PLS","fullmodel", "mean")
rownames(trial)="prediction error"
kable(trial)
```



|                 | regsubsets|   LASSO|      PCR|      PLS| fullmodel|    mean|
|:----------------|----------:|-------:|--------:|--------:|---------:|-------:|
|prediction error |   1.407029| 1.36724| 1.397812| 1.399506|  1.397812| 1.68628|

According to the square loss function the best model is built using Lasso, nevertheless the results of the four methods are similar and by running again the analysis using different folds, other methods give a minimum for the square loss function. All methods propose models that are better than the mean, but not always better than the full model. 


```r
trial <-matrix(c(besti_sub,bestp,besti_PCR,besti_PLS) ,ncol=10,nrow=4)
rownames(trial) <-c("regsubsets","LASSO","PCR","PLS" )
kable(trial,caption = "number of variables")
```




|           |   |   |   |   |   |   |   |   |   |   |
|:----------|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|regsubsets | 12|  8|  8| 12| 12| 12| 12| 12| 10| 12|
|LASSO      |  9| 11| 10| 10| 12| 12| 12| 12|  9| 12|
|PCR        | 12| 10| 12| 12| 12| 12| 12|  8| 10| 12|
|PLS        |  8| 10| 12| 12| 12| 12| 12| 12| 11|  9|

The number of variables included in the model are not always different from the case of the full model, this explains why the square loss function of the full model is slightly different from the one of the proposed models.
In this case with such small number of variables and so little improvement a dimension reduction procedure perhaps it's not advisable.

# Exercise 2

$$
\begin{aligned}
 \\[0,3in]
 \beta_c=(C'C)^{-1}C'Y=(XA'AX')^{-1}XA'Y=(XX')^{-1}XA'Y
\end{aligned}
$$

# Exercise 3

$$
\begin{aligned}
 IF(X,P,T)&=\lim_{x\to\infty}\frac{T((1-\epsilon)E_p[(X)^2]+\epsilon\delta_x)-T(P)}{} \\
 &=\lim_{x\to\infty}\frac{(1-\epsilon)-1}{\epsilon}=-1
\end{aligned}
$$




# Exercise 4


```r
crime <- read.csv("london-borough-profilesf.csv")


i=1
c=c()
pb <- txtProgressBar(min = 0, max = 85, style = 3)
```

`

```r
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
```


```r
crime_new=crime[,-c]
crime_new=crime_new[,-c(31,32)]
crime_new=crime_new[,-c(66,65,60,21,22,27,31,33,47,50,52)]

crime_new=crime_new[,-c(6,7,8,3,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,
                        26,27,28,29,30,36,37,38,39,40,41,43,44,45,45,46,
                        48,49,57,58,59,60)]
crime_new$GLA=(crime_new$GLA.Population.Estimate.2015+crime_new$GLA.Household.Estimate.2015)/2
crime_new$NoAnxiety.score.2011.14..out.of.10.=10-crime_new$Anxiety.score.2011.14..out.of.10.
crime_new$Psychological_status=(crime_new$Life.satisfaction.score.2011.14..out.of.10.+
                                  crime_new$Worthwhileness.score.2011.14..out.of.10.+
                                  crime_new$NoAnxiety.score.2011.14..out.of.10.+crime_new$Happiness.score.2011.14..out.of.10.)/4
crime_new$life.expectancy=(crime_new$Female.life.expectancy...2011.13.+
                             crime_new$Male.life.expectancy...2011.13.)/2

crime_new=crime_new[,-c(22,1:2,14:15,16:19)]

names(crime_new)
```

```
##  [1] "Population.density..per.hectare..2015"             
##  [2] "Average.Age..2015"                                 
##  [3] "Net.internal.migration..2014."                     
##  [4] "Unemployment.rate..2014."                          
##  [5] "Two.year.business.survival.rates..started.in.2011."
##  [6] "Crime.rates.per.thousand.population.2014.15"       
##  [7] "Fires.per.thousand.population..2014."              
##  [8] "Ambulance.incidents.per.hundred.population..2014." 
##  [9] "Median.House.Price..2014"                          
## [10] "Total.carbon.emissions..2013."                     
## [11] "Rates.of.Children.Looked.After..2014."             
## [12] "Mortality.rate.from.causes.considered.preventable" 
## [13] "GLA"                                               
## [14] "Psychological_status"                              
## [15] "life.expectancy"
```

This are the variables 15 that I suggest to include in the model. A lot of variables have been excluded because they are irrelavant or don't add much information to the variables considered. 12 variables are the original ones and 3 are the result of combining the original variables to obtain an indicator that can summarize them:

-GLA is given by  (GLA.Population.Estimate.2015 + GLA.Household.Estimate.2015)/2
-Psychological_status is given by  (Worthwhileness.score + Happiness.score + (10- Anxiety.score) + Life.satisfaction.score)/4.
-life.expectancy is given by (Male.life.expectancy+Female.life.expectancy)/2
