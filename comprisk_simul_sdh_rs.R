#-.-.-.-.-.--.-.-.-.-.-.--.-.-.-.-.-.--.-.-.--.-.-.--.
# Main function used for simulating survival times using
# a subdistribution hazard (gamma_0) for the cause of interest
# and population life tables for the other causes.
# The function computes also the cause-specific hazard for the cause
# of interest in order to determine the occurence of any event
#-.-.-.-.-.--.-.-.-.-.-.--.-.-.-.-.-.--.-.-.--.-.-.--.

#ncl          :   Number of clusters (used in doSNOW)
#betaagec     :   BETA for age
#betasex      :   BETA for sex
#betayearcr   :   BETA for year
#kappa,rho, alpha : parameters for subdistribution hazard
#adm.cens     : administrative censoring
#drop.out     : lost to follow-up
#lifetable.DF : population life tables in data.frame format (names: AGERT, SEXRT, YEARRT)
#lifetable.RT : population life tables in rate.table format (names: AGE.RT, SEX.RT, YEAR.RT)
#Simdata      : population data for which we want to generate the survival times


#Note 1: betaagec, betasex, betayearcr correspond to the beta coef. used in
#        the subdistribution hazard model for the coause of interest. The corresponding 
#        covariates are identified in the data as (agecr, IsexH, yearcr) and they are 
#        defined as [age-mean(age)]/10, sex(0,1), [year-mean(year)]/10.

#Note 2: the function as it is accounts for life tables which are stratified by  
#        age, sex, year. The variables are identified must be able to identify 
#        in the data also as age (continuous), sex (binary), year(continuous)).
#        The user can account for different life tables by
#        adapting STEP 1 from the code below. 


# required packages
reqPcks <- c("doSNOW", "lubridate","plyr","survival", "statmod")
for(p in reqPcks){
  if(!require(p, character.only=TRUE)) {
    install.packages(p)
  }
}

# main function
sim_rel <- function(ncl,betaagec, betasex, betayearcr, kappa,
                    rho, alpha, adm.cens, drop.out, lifetable.DF,
                    lifetable.RT, Simdata){
  
  #we use doSNOW to parallelise the procedure
  library(doSNOW)
  cl<-makeCluster(ncl, type = "SOCK")
  registerDoSNOW(cl)
  
  Simdata$lambda1 <- Simdata$lambda2 <- Simdata$cause <- Simdata$finaltime <- 
    Simdata$timesurv<- Simdata$finalcause <- Simdata$cens<- 999
  
  
  SimdataALL.df <- foreach(iloop = 1:nrow(Simdata),.combine = "rbind") %dopar% {
    
    library(lubridate)
    library(plyr)
    library(survival)
    library(statmod)
    
    Xi= c(Simdata$agecr[iloop],Simdata$IsexH[iloop], Simdata$yearcr[iloop])

    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    # STEP 1 : Population mortality (attributed to other causes)
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
    
    # Note: age, year, sex are continuous
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
    
    #-.-.-.-.-.--.-.-.--.-.-.--.-.-.-
    # STEP 3: Subdistribution hazard
    #-.-.-.-.-.--.-.-.--.-.-.--.-.-.-
    
    gamma_0<- Vectorize(function(t) {kappa*(rho^kappa)*(t^(kappa-1))/(1+
                                                                        ((rho*t)^kappa)/alpha)})
    gamma_ft <- Vectorize(function(t) {exp(c(betaagec, betasex,
                                             betayearcr)%*%Xi)*gamma_0(t)})
    cumgam_ft<- Vectorize(function(lo,up) {integrate(gamma_ft, lower = lo,
                                                     upper = up)[1]$value})
    
    #-.-.-.-.-.--.-.-.--.-.-.--.-.-.-
    # STEP 4: Cancer-specific hazard
    #-.-.-.-.-.--.-.-.--.-.-.--.-.-.-
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
  
  
    #-.-.-.-.-.--.-.-.--.-.-.--
    # STEP 5: Determine event
    #-.-.-.-.-.--.-.-.--.-.-.--
    
    lambda_tot <- function(t) {lambda1_ft(t)+lambda2_ft(t)}
    
    cumlam_tot <- function(t) {cumlam1_ft(t) + cumlam2_ft(t)}
    
    u <- runif(1)
    temp1 <- Vectorize(function(t,u){if (t==0) {myf=-u} else 
      myf=(1-exp(-cumlam_tot(t))- u) 
    return(myf)})
    
    res <- try(uniroot(temp1, interval=c(0.001,adm.cens), u=u, tol=0.001), silent=T)
    if (class(res)=="try-error"){stime=adm.cens} else {stime=res$root}
    
    
    lambda1<- lambda1_ft(stime)
    lambda2<- as.numeric(lambda2_ft(stime))
    
    Cause <- sample(1:2, 1, prob = c(lambda1, lambda2))
    Simdata$lambda1[iloop] <- lambda1
    Simdata$lambda2[iloop] <- lambda2
    
    
    Simdata$cause[iloop] <- Cause
    Simdata$finaltime[iloop] <- stime
    
  
  #-.-.-.-.-.--.-.-.-.--.-.-
  # STEP 6: Apply Censoring
  #-.-.-.-.-.--.-.-.-.--.-.-
  
  #Censoring: administrative+drop outs
  Simdata$cause[iloop] <- ifelse(stime==adm.cens,0,Cause)
  temp.ut2 <- runif(nrow(Simdata))
  Simdata$cens<- temp.ut2/drop.out
  
  
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
  # STEP 6: Final survival time
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
  
  Simdata$timesurv <- sapply(1:nrow(Simdata), function (x) min(Simdata$finaltime[x], 
                                                               Simdata$cens[x]))
  Simdata$finalcause <- sapply(1:nrow(Simdata), function (x) 
    ifelse(Simdata$finaltime[x]==Simdata$timesurv[x], Simdata$cause[x], 0) )
  
  return(Simdata[iloop,])
  
  }
  return(SimdataALL.df)
}

#------------------------------------------------------------------------------
#----------------------------- T E S T ----------------------------------------
#------------------------------------------------------------------------------
               #FOR THE TEST USE THE examp_data_for_simulation.R.RData
                               
                          
Simulasion_data <- sim_rel(ncl = 20,
                           betaagec = 0.22,betasex = 0.26, betayearcr = 0.009,
                           kappa=2,rho=1.6, alpha=0.05,
                           adm.cens=10, drop.out=0.035,
                           lifetable.DF=exprates.DF.UKdep9,lifetable.RT= exprates.RT.UKdep9, Simdata = Simdat.exa)

