#-.-.-.-.-.-
# STEP 1
#-.-.-.-.-.-

cDataDesignOptim <- function (n=NULL, cens.admin, ydiagmin, ydiagmax){
  
  id 			<- c(1:n)
  sex 		<- ifelse(runif(n) < 0.3, 2, 1)
  IsexH 		<- ifelse(sex==1,1,0)
  year 		<- runif(n, min = ydiagmin, max = ydiagmax)
  age   		<- c(runif(n*0.3, min = 30, max = 65),
               runif(n*0.3, min = 65, max = 75),
               runif(n*0.40, min = 75, max = 85))
  
  tab2              <- data.frame(id, sex, IsexH,year, age)
  tab2$agecr        <- (tab2$age-70)/10
  tab2$yearcr       <- tab2$year-2002
  tab2 <- tab2[order(tab2$id),]
  
  return(list(tab2 = tab2, n=n, cens.admin = cens.admin, ydiagmin=ydiagmin,
              ydiagmax=ydiagmax))
}


Simdat.exa<- cDataDesignOptim(n = 10, cens.admin = 2014, ydiagmin=2000, 
                              ydiagmax=2003)$tab2


###############################################################
# Main function for simulating survival times
sim_rel <- function(sdh,betaagec=NULL, betasex=NULL, betayearcr=NULL, kappa=NULL,
                    rho=NULL, alpha=NULL, adm.cens=NULL, drop.out=NULL, lifetable.DF=NULL,
                    lifetable.RT= NULL, Simdata=Simdat.exa){
  
  reqPcks <- c("doSNOW", "lubridate","plyr","survival", "statmod")
  for(p in reqPcks){
    if(!require(p, character.only=TRUE)) {
      install.packages(p)
    }
  }
  
  library(doSNOW)
  cl<-makeCluster(10, type = "SOCK")
  registerDoSNOW(cl)
  
  Simdata$lambda1 <- Simdata$lambda2 <- Simdata$cause <- Simdata$finaltime <- 
    Simdata$timesurv<- Simdata$finalcause <- Simdata$cens<- 999
  
  
  SimdataALL.df <- foreach(iloop = 1:nrow(Simdata),.combine = "rbind") %dopar% {
    
    library(lubridate)
    library(plyr)
    library(survival)
    library(statmod)
    
    Xi= c(Simdata$agecr[iloop],Simdata$IsexH[iloop], Simdata$yearcr[iloop])
    
    #-.-.-.-.-.-
    # STEP 2
    #-.-.-.-.-.-
    
    lambda2_ft <- Vectorize(function(t){
      lambda2 <- subset(lifetable.DF, select="exprate", 
                        subset=(AGERT==trunc(Simdata$age[iloop]+t) & 
                                  YEARRT==trunc(Simdata$year[iloop]+t) &
                                  SEXRT==Simdata$sex[iloop])
      )
      return(lambda2)
    })
    
    cumlam2_ft <- Vectorize(function(t){
      ddiag<-format(date_decimal(Simdata$year[iloop]), "%d-%m-%Y")
      yearout <- format(date_decimal(Simdata$year[iloop]+t), "%d-%m-%Y")
      agediag<-Simdata$age[iloop]
      sex<-Simdata$sex[iloop]
      fuptimeday1<-as.Date(yearout, format="%d-%m-%Y")-as.Date(ddiag, 
                                                               format="%d-%m-%Y")
      fuptimeday <- ifelse(fuptimeday1==0,1,fuptimeday1)
      #number of days
      ddiagday<-as.numeric(as.Date(ddiag, format="%d-%m-%Y")-as.Date("01-01-1960", 
                                                                     format="%d-%m-%Y"))
      
      agediagday<-Simdata$age[iloop]*365.25
      
      #should be changed according to the attributes of the object 'lifetable'
      cumlam2_x<- -log(survexp(~ ratetable(AGE.RT = agediagday, SEX.RT = sex, 
                                           YEAR.RT = ddiagday),
                               ratetable = lifetable.RT,
                               times = fuptimeday)$surv)
      return(cumlam2_x)
    })
    
    #-.-.-.-.-.-
    # STEP 3
    #-.-.-.-.-.-
    
    gamma_0<- Vectorize(function(t) {kappa*(rho^kappa)*(t^(kappa-1))/(1+
                                                                        ((rho*t)^kappa)/alpha)})
    gamma_ft <- Vectorize(function(t) {exp(c(betaagec, betasex,
                                             betayearcr)%*%Xi)*gamma_0(t)})
    cumgam_ft<- Vectorize(function(lo,up) {integrate(gamma_ft, lower = lo,
                                                     upper = up)[1]$value})
    
    GLw1 <- gauss.quad(n=1, kind="legendre")
    GLw10 <- gauss.quad(n=10, kind="legendre")
    GLw <- gauss.quad(n=30, kind="legendre")
    Rescale <- function(gl,a,b){
      gl$mynod <- gl$nodes*(b-a)/2+(a+b)/2
      gl$myw <- gl$weights*(b-a)/2
      return(gl)
    }
    
    lambda1_ft <- Vectorize(function(t){
      Numer <- gamma_ft(t)*exp(-cumgam_ft(0,t)+cumlam2_ft(t))
      fDenom <- Vectorize(function(z) {gamma_ft(z)*exp(-cumgam_ft(0,z)+cumlam2_ft(z))})
      if (t<4) {gg <- Rescale(GLw10,0,t)} else {gg <- Rescale(GLw,0,t)}
      if (t<0.1){gg <- Rescale(GLw1,0,t)}
      aa <- fDenom(gg$mynod)
      Denom <- 1-sum(gg$myw*aa)
      return(Numer/Denom)
    })
    
    cumlam1_ft <- function(t){
      if (t<4) {gg <- Rescale(GLw10,0,t)} else {gg <- Rescale(GLw,0,t)}
      if (t<0.1){gg <- Rescale(GLw1,0,t)}
      bb=lambda1_ft(gg$mynod)
      cum.GL=sum(gg$myw*bb)
      if (t<0.003) cum.GL<-0.00001
      return(cum.GL)
    }
  } 
  
  #-.-.-.-.-.-
  # STEP 4
  #-.-.-.-.-.-
  
  lambda_tot <- function(t) {lambda1_ft(t)+lambda2_ft(t)}
  
  cumlam_tot <- function(t) {cumlam1_ft(t) + cumlam2_ft(t)}
  
  u <- runif(1)
  temp1 <- Vectorize(function(t,u){if (t==0) {myf=-u} else 
    myf=(1-exp(-cumlam_tot(t))- u) 
  return(myf)})
  
  res <- try(uniroot(temp1, interval=c(0.001,adm.cens), u=u, tol=0.001), silent=T)
  if (class(res)=="try-error"){stime=adm.cens} else {stime=res$root}
  
  
  lambda1<- lambda1_ft(stime)
  lambda2<- lambda2_ft(stime)
  
  Cause <- sample(1:2, 1, prob = c(lambda1, lambda2))
  Simdata$lambda1[iloop] <- lambda1
  Simdata$lambda2[iloop] <- lambda2
  
  
  Simdata$cause[iloop] <- Cause
  Simdata$finaltime[iloop] <- stime
  
  
  #-.-.-.-.-.-
  # STEP 5
  #-.-.-.-.-.-
  
  #Censoring: administrative+drop outs
  Simdata$cause[iloop] <- ifelse(stime==adm.cens,0,Cause)
  temp.ut2 <- runif(nrow(Simdata))
  Simdata$cens<- temp.ut2/drop.out
  
  
  #-.-.-.-.-.-
  # STEP 6
  #-.-.-.-.-.-
  Simdata$timesurv <- sapply(1:nrow(Simdata), function (x) min(Simdata$finaltime[x], 
                                                               Simdata$cens[x]))
  Simdata$finalcause <- sapply(1:nrow(Simdata), function (x) 
    ifelse(Simdata$finaltime[x]==Simdata$timesurv[x], Simdata$cause[x], 0) )
  
  return(Simdata[iloop,])
  
}

stopCluster(cl)
return(SimdataALL.df)
}
