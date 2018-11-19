# Change appropriately the folder path where simdatn2.RData 
# and the expectedrates.RT.dat are located
#setwd()

# Install the needed packages
reqPcks <- c("relsurv", "survival","geepack")

for(p in reqPcks){
  if(!require(p, character.only=TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)}
}

packageVersion("relsurv")
packageVersion("survival")

# Load data
load("data_pseudo_tutorial2.RData")

str(simdatn2)
head(simdatn2)
str(expectedrates.RT)
N <- nrow(simdatn2)

# Convert (1) age at diagnosis, (2) time from 01-01-1960 until diag date in 
# days.
#-----------------------------------------------------------------------------
# (1) Age at diagnosis in days 
simdatn2$agediagdays <- round(simdatn2$age*365.241)

# (2) Transform continuous year of diagnosis the date of diagnosis 
# and subtract it from 01-01-1960.
simdatn2$diag_year<- floor(simdatn2$year)
simdatn2$diag_moCont <- as.numeric(substr(simdatn2$year,5,
                                          nchar(simdatn2$year)))
simdatn2$diagDate <- as.Date(paste(simdatn2$diag_year,
                                   "-01-01",sep=""))+
  round(simdatn2$diag_moCont*365.241)

simdatn2$diagdays1960 <- as.numeric(simdatn2$diagDate-
                                      as.Date("01/01/1960",
                                              format="%d/%m/%Y"))

# (3) Survival time in days
simdatn2$timesurvD <- floor(simdatn2$timesurv*365.241)
times_c <- quantile(simdatn2$timesurv, probs=seq(0.15,1,0.15))
times_c <- round(times_c,1)
print(times_c)

#------------------------------------------------------------------
# L E A V E  O N E  O U T  E S T I M A T O R
#------------------------------------------------------------------

# Thetas based on the whole sample
fit_all <- cmp.rel(Surv(timesurvD,vstat)~1+ratetable(AGE.RT=agediagdays,
                                                     SEX.RT=sex,
                                                     YEAR.RT=diagdays1960),
                   ratetable=expectedrates.RT,data=simdatn2,tau=3652.41,
                   conf.int=0.95)

results_relsurv<- list(summary(fit_all, times = times_c)$est, 
                       cbind(fit_all$causeSpec$area, fit_all$population$area))


ls <- list()

#--------------------------------------------------------------------------------
# Tip: you can speed the process by replacing the for loop for example
# with foreach from package doSNOW
    install.packages("doSNOW")
    library(doSNOW)
    clustersNum <- 8 #number of clusters selected
    cl<-makeCluster(8, type = "SOCK")
    registerDoSNOW(cl)
    
    ls_foreach<- foreach(y = 1:nrow(simdatn2),.combine = 'rbind') %dopar% {
      library(relsurv)

      #Fit the model

      fit <- cmp.rel(Surv(timesurvD,vstat)~ratetable(AGE.RT=agediagdays,SEX.RT=sex,
                                                     YEAR.RT=diagdays1960),
                     ratetable=expectedrates.RT,data=simdatn2[-y,],tau=3652.41,
                     conf.int=0.95)
      #This will be stored
      list(summary(fit, times = times_c)$est,
                  cbind(fit$causeSpec$area, fit$population$area))
      
    }
    stopCluster(cl)
    
    #Rearrange results as to run smoothly with the rest of the code
    ls<- lapply(1:N, function(x)list(ls_foreach[[x]], ls_foreach[[N+x]]))
#--------------------------------------------------------------------------------

# 
# for (y in 1:nrow(simdatn2)){
#   fit <- cmp.rel(Surv(timesurvD,vstat)~ratetable(AGE.RT=agediagdays,SEX.RT=sex,
#                                                  YEAR.RT=diagdays1960),
#                  ratetable=expectedrates.RT,data=simdatn2[-y,],tau=3652.41,
#                  conf.int=0.95)
#   ls[[y]] <- list(summary(fit, times = times_c)$est, 
#                   cbind(fit$causeSpec$area, fit$population$area)) # to be stored
# }

# Separate estimates based on the indicator and the cause 

CPr.1 <- t(sapply(1:N, function (x) ls[[x]][[1]][1,])) 
#dataframe dimensions: (N x times_c)
CPr.2 <- t(sapply(1:N, function (x) ls[[x]][[1]][2,]))

lyl.1 <- sapply(1:N, function (x) ls[[x]][[2]][,1]) 
#dataframe dimensions: (N x 1)
lyl.2 <- sapply(1:N, function (x) ls[[x]][[2]][,2])


# Final step: leave - one - out estimator

pseudo_CPr.1<-data.frame(matrix(1,nrow(simdatn2),length(times_c)))
pseudo_CPr.2<-data.frame(matrix(1,nrow(simdatn2),length(times_c)))
colnames(pseudo_CPr.1)<- colnames(pseudo_CPr.2)<- 
  paste(times_c,"y",sep="")

pseudo_lyl.1<-data.frame(matrix(1,nrow(simdatn2),1))
pseudo_lyl.2<-data.frame(matrix(1,nrow(simdatn2),1))
colnames(pseudo_lyl.1)<- colnames(pseudo_lyl.2)<- 
  paste(max(times_c),"y",sep="")

for(y in 1:length(times_c)){
  for (x in 1:N){
    pseudo_CPr.1[x,y]<- N*results_relsurv[[1]][1,y]-(N-1)*CPr.1[x,y]
    pseudo_CPr.2[x,y]<- N*results_relsurv[[1]][2,y]-(N-1)*CPr.2[x,y]
  }
}

for (x in 1:N){
  pseudo_lyl.1[x,]<- N*results_relsurv[[2]][,1]-(N-1)*lyl.1[x]
  pseudo_lyl.2[x,]<- N*results_relsurv[[2]][,2]-(N-1)*lyl.2[x]
  
}

# Crude Probability of Death
#---------------------------
linkfunc <- "identity"
covastr <- "independence"

for (h in 1:2){
  pseudo <- get(paste("pseudo_CPr",h,sep="."))
  
  b <- NULL
  for (it in 1:length(times_c)) {
    b <- rbind(b, cbind(simdatn2,
                        pseudo = pseudo[,it],
                        tpseudo = times_c[it],
                        id = 1:nrow(simdatn2)))
  }
  b <- b[order(b$id), ]
  assign(paste("b_cpd", h, sep="."),b)
  
  
  #Put sex always 2nd and interaction before any variables!
  
  pseudo_fit <- geese(pseudo ~ as.factor(tpseudo) +
                        agecr+sex+yearcr, data = b, id = id,
                      jack = TRUE, scale.fix = TRUE,
                      family = gaussian, mean.link = linkfunc,
                      corstr = covastr)
  assign(paste("pseudo_fit_cpd", h, sep="."),pseudo_fit)
}

print(paste("Cause: Cancer"))
print(summary(pseudo_fit_cpd.1)$mean)
print(paste("Cause: Other causes"))
print(summary(pseudo_fit_cpd.2)$mean)

# Life years lost
#---------------------------
linkfunc <- "identity"
covastr <- "independence"
times_lyl<- max(times_c)

for (h in 1:2){
  pseudo <- get(paste("pseudo_lyl",h,sep="."))
  
  b <- NULL
  for (it in 1:length(times_lyl)) {
    b <- rbind(b, cbind(simdatn2,
                        pseudo = pseudo[,it],
                        tpseudo = times_lyl[it],
                        id = 1:nrow(simdatn2)))
  }
  b <- b[order(b$id), ]
  assign(paste("b_lyl", h, sep="."),b)
  
  pseudo_fit <- geese(pseudo ~ agecr+sex+yearcr, data = b, id = id,
                      jack = TRUE, scale.fix = TRUE,
                      family = gaussian, mean.link = linkfunc,
                      corstr = covastr)
  assign(paste("pseudo_fit_lyl", h, sep="."),pseudo_fit)
}

print(paste("Cause: Cancer"))
print(summary(pseudo_fit_lyl.1)$mean)
print(paste("Cause: Other causes"))
print(summary(pseudo_fit_lyl.2)$mean)


pseudo_pred <- function(mod, tpseudo, var1, var2, var3, data, p, times){
  
  xis  <-data.frame(data[,c(tpseudo,var1, var2,var3)])
  
  Ntim <- length(times)
  N <- nrow(data)/Ntim
  xiss <-data.frame(matrix(1,nrow(data),Ntim))
  xis<-cbind(xis,xiss)
  
  betas<- mod$beta
  cov_vars <- mod$vbeta.ajs
  
  pred<-data.frame(matrix(1,N,Ntim))
  var<-NULL
  
  y=p
  
  if (Ntim!=1){
    for (i in 2:Ntim){
      u<-y+1+i
      xis[,u]<- ifelse(xis[,tpseudo]==times[i],1,0)
    }
    
    xis_transf <- as.matrix(cbind(xis[,c((ncol(xis)-
                                            Ntim+1):ncol(xis))],
                                  xis[,c(2:(ncol(xis)-Ntim))]))
    #Define separately the first point
    xis_transf.1 <- subset(xis_transf,xis_transf[,1]==1 &
                             rowSums(xis_transf[,2:Ntim])==0)
    weights.1 <- matrix(rep(1/N,nrow(xis_transf.1)), nrow =1)
    new_xis_trans.1 <- weights.1%*%xis_transf.1
    var[1]<- as.numeric(new_xis_trans.1%*%cov_vars%*%t(new_xis_trans.1))
    pred[,1] <- as.numeric(betas%*%t(xis_transf.1))
    
    for (l in 2:Ntim){
      xis_transf.l <- subset(xis_transf, xis_transf[,l]==1)
      weights.l <- matrix(rep(1/N,nrow(xis_transf.l)), nrow =1)
      new_xis_trans.l <- weights.l%*%xis_transf.l
      var.l<- as.numeric(new_xis_trans.l%*%cov_vars%*%t(new_xis_trans.l))
      pred.l <- as.numeric(betas%*%t(xis_transf.l))
      
      
      pred[,l] <- pred.l
      var[l] <- var.l
    }
    
  } else if (Ntim==1){
    xis_transf <- as.matrix(cbind(xis[,ncol(xis)],
                                  xis[,c(2:(ncol(xis)-1))]))
    
    weights.1 <- matrix(rep(1/N,nrow(xis_transf)), nrow =1)
    new_xis_trans.1 <- weights.1%*%xis_transf
    var[1]<- as.numeric(new_xis_trans.1%*%cov_vars%*%t(new_xis_trans.1))
    pred[,1] <- as.numeric(betas%*%t(new_xis_trans.1))
  }   
  
  if (pseudo_fit$model$mean.link=="cloglog"){
    est <- colMeans(1-exp(-exp(pred)))
    lower <- 1-exp(-exp(log(-log(1-est))-qnorm(0.975)*sqrt(var)))
    upper <- 1-exp(-exp(log(-log(1-est))+qnorm(0.975)*sqrt(var)))
    var1 <- var*((1-est)*log(1-est))^2
  } else if (pseudo_fit$model$mean.link=="log"){
    est <- colMeans(exp(pred))
    upper <- exp(log(est)+qnorm(0.975)*sqrt(var))
    lower<-  exp(log(est)-qnorm(0.975)*sqrt(var))
    var1 <- var*(est^2 )
  } else if (pseudo_fit$model$mean.link=="identity"){
    est <- colMeans(pred)
    lower <- est-qnorm(0.975)*sqrt(var)
    upper <- est+qnorm(0.975)*sqrt(var)
    var1 <- var
  }
  return(as.data.frame(t(data.frame(rbind(est,lower,upper,var1)))))
}


# N o n - P a r a m e t r i c
NP_fit <- cmp.rel(Surv(timesurvD,vstat)~ratetable(AGE.RT=agediagdays,
                                                  SEX.RT=sex,
                                                  YEAR.RT=diagdays1960),
                  ratetable=expectedrates.RT,data=simdatn2,
                  tau=max(times_c)*365.241,conf.int=0.95)


# C r u d e  p r o b a b i l i t i e s 

for (h in 1:2){  
  pseudo_fit <- get(paste("pseudo_fit_cpd",h,sep="."))
  b <- get(paste("b_cpd",h,sep="."))
  
  pseudo_data<-pseudo_pred(mod=pseudo_fit, tpseudo="tpseudo",
                           var1="agecr", var2="sex", var3="yearcr",
                           data=b,p=3, times=times_c)
  rownames(pseudo_data)<- paste(paste(paste(paste("CPr",h,sep="."), 
                                            "at",sep=" "), times_c, sep=" "), 
                                "y",sep="")
  assign(paste("pseudo_data_cpd", h, sep="."),pseudo_data)
}


# L i f e  Y e a r s  L o s t

for (h in 1:2){  
  pseudo_fit <- get(paste("pseudo_fit_lyl",h,sep="."))
  b <- get(paste("b_lyl",h,sep="."))
  
  pseudo_data<-pseudo_pred(mod=pseudo_fit, tpseudo="tpseudo",
                           var1="agecr", var2="sex", var3="yearcr",
                           data=b,p=3, times=max(times_c))
  assign(paste("pseudo_data_lyl", h, sep="."),pseudo_data)
}


# Plot to compare non-parametric and model estimates for CPr
plot(fit, conf.int = TRUE, ylab="Crude probabilities of death", ylim=c(0,0.5))
points(times_c*365.241,pseudo_data_cpd.1$est, col="red", pch=19)
arrows(times_c*365.241,pseudo_data_cpd.1$lower,
       times_c*365.241,pseudo_data_cpd.1$upper,
       length=0.05, angle=90, code=3, col="red")
points(times_c*365.241,pseudo_data_cpd.2$est, col="red", pch=19)
arrows(times_c*365.241,pseudo_data_cpd.2$lower,
       times_c*365.241,pseudo_data_cpd.2$upper,
       length=0.05, angle=90, code=3, col="red")