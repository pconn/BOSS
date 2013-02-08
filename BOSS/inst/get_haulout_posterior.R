library(BOSSPow)
data(glmmLDTS.ribbon.fit)
data(glmmLDTS.bearded.fit)
data(glmmLDTS.spotted.fit)

set.seed(123456)
load("c:/users/paul.conn/git/BOSS/BOSS/data/Effort_points.Rdat")

### 1: single value for each species
#draw random sample of date-times from on effort flights
effort_data=effort_data[-which(effort_data[["flightid"]]=="AeroFl10"),]  
DT=as.POSIXlt(effort_data[["gps_dt"]])
effort_data=effort_data[which(DT$mon==3 & DT$mday>20 & DT$mday<29),] #mon goes 0-11
DT.samp=sample(effort_data[["gps_dt"]],10000,replace=TRUE)
DT.samp=as.POSIXlt(DT.samp)

Day <- DT.samp$yday/365 
Hour<-DT.samp$hour

#ribbon seals
RibbonFit <- rep(0,10000)
RibbonH <- rep(0,10000)
RibbonCV <- rep(0,10000)
for(i in 1:10000) {
    hour <- Hour[i]
    day <- Day[i]

    DateFit <- glmmLDTS.ribbon.fit$fixed.effect[
      glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
    glmmLDTS.ribbon.fit$fixed.effect[
      glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
      glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day +
    glmmLDTS.ribbon.fit$fixed.effect[
      glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
      glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day^2 +
    glmmLDTS.ribbon.fit$fixed.effect[
      glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
      glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day^3

		L <- (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "intercept")*1 +
		(glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
      glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day +
		(glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
      glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day^2 +
		(glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
      glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day^3

		ell <- sum(glmmLDTS.ribbon.fit$fixed.effect$estimate*L)
		h <- 1/(1+exp(-ell))
		Cell <- t(L) %*% glmmLDTS.ribbon.fit$covb %*% L
		varh <- (exp(ell)/(1+exp(ell))^2)^2*Cell
		
    RibbonFit[i] <- DateFit
		RibbonH[i] <- h
		RibbonCV[i] <- sqrt(varh)/h
}

h.ribbon=mean(RibbonH)
cv.ribbon=mean(RibbonCV)


#spotted seals
SpottedFit <- rep(0,10000)
SpottedH <- rep(0,10000)
SpottedCV <- rep(0,10000)
for(i in 1:10000) {
  hour <- Hour[i]
  day <- Day[i]
  
  DateFit <- glmmLDTS.spotted.fit$fixed.effect[
    glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
    glmmLDTS.spotted.fit$fixed.effect[
      glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
        glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day +
    glmmLDTS.spotted.fit$fixed.effect[
      glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
        glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day^2 +
    glmmLDTS.spotted.fit$fixed.effect[
      glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
        glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day^3
  
  L <- (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "intercept")*1 +
    (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
       glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day +
    (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
       glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day^2 +
    (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
       glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day^3
  
  ell <- sum(glmmLDTS.spotted.fit$fixed.effect$estimate*L)
  h <- 1/(1+exp(-ell))
  Cell <- t(L) %*% glmmLDTS.spotted.fit$covb %*% L
  varh <- (exp(ell)/(1+exp(ell))^2)^2*Cell
  
  SpottedFit[i] <- DateFit
  SpottedH[i] <- h
  SpottedCV[i] <- sqrt(varh)/h
}

h.spotted=mean(SpottedH)
cv.spotted=mean(SpottedCV)

#bearded seals
BeardedFit <- rep(0,10000)
BeardedH <- rep(0,10000)
BeardedCV <- rep(0,10000)
for(i in 1:10000) {
  hour <- Hour[i]
  day <- Day[i]
  
  DateFit <- glmmLDTS.bearded.fit$fixed.effect[
    glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
    glmmLDTS.bearded.fit$fixed.effect[
      glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
        glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day +
    glmmLDTS.bearded.fit$fixed.effect[
      glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
        glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day^2 +
    glmmLDTS.bearded.fit$fixed.effect[
      glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
        glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""),
      "estimate"]*day^3
  
  L <- (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "intercept")*1 +
    (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
       glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day +
    (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
       glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day^2 +
    (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
       glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(hour,",", sep = ""))*day^3
  
  ell <- sum(glmmLDTS.bearded.fit$fixed.effect$estimate*L)
  h <- 1/(1+exp(-ell))
  Cell <- t(L) %*% glmmLDTS.bearded.fit$covb %*% L
  varh <- (exp(ell)/(1+exp(ell))^2)^2*Cell
  
  BeardedFit[i] <- DateFit
  BeardedH[i] <- h
  BeardedCV[i] <- sqrt(varh)/h
}

h.bearded=mean(BeardedH)
cv.bearded=mean(BeardedCV)

### 2: lookup table for each species

startDate <- "2006-04-21 01:00:00"
startDay <- (as.POSIXlt(as.POSIXct(startDate))$yday+1)/365
endDate <- "2006-04-28 01:00:00"
endDay <- (as.POSIXlt(as.POSIXct(endDate))$yday+1)/365
DateRange <- startDay + (0:7)/7*(endDay - startDay)
HourRange <-  as.character(0:23)
fracYear2POSIX(2006,DateRange)


# ------------------------------------------------------------------------------
#                   RIBBON SEALS
# ------------------------------------------------------------------------------



Fit.ribbon <- matrix(NA, nrow = length(DateRange), ncol = 24)
L=rep(NA,length(glmmLDTS.ribbon.fit$fixed.effect[,"effect"]))
for(i in 1:length(DateRange)) {
  for (j in 1:24) {
    
    Hour <- HourRange[j]
    Day <- DateRange[i]
    
    Fit.ribbon[i,j] <- glmmLDTS.ribbon.fit$fixed.effect[
      glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
      glmmLDTS.ribbon.fit$fixed.effect[
        glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
          glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day +
      glmmLDTS.ribbon.fit$fixed.effect[
        glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
          glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^2 +
      glmmLDTS.ribbon.fit$fixed.effect[
        glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
          glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^3
    
    L <- rbind(L, (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "intercept")*1 +
      (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
         glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day +
      (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
         glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^2 +
      (glmmLDTS.ribbon.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
         glmmLDTS.ribbon.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^3)    
  }
}
L=L[-1,]
Ell <- L%*%glmmLDTS.ribbon.fit$fixed.effect$estimate
Cell <- L %*% glmmLDTS.ribbon.fit$covb %*% t(L)
#generate predictions on the logit scale and back transform
H.ribbon=array(0,dim=c(dim(Fit.ribbon),1000))
for(i in 1:1000){
  if(i%%10==0)cat(paste('iter ',i,'\n'))
  H.ribbon[,,i]=matrix(rmvnorm(1,mean=Ell,sigma=Cell,method="svd"),nrow=length(DateRange),ncol=24,byrow=TRUE)
}
H.ribbon=1/(1+exp(-H.ribbon)) #back transform

# ------------------------------------------------------------------------------
#                   Bearded SEALS
# ------------------------------------------------------------------------------



Fit.bearded <- matrix(NA, nrow = length(DateRange), ncol = 24)
L=rep(NA,length(glmmLDTS.bearded.fit$fixed.effect[,"effect"]))
for(i in 1:length(DateRange)) {
  for (j in 1:24) {
    
    Hour <- HourRange[j]
    Day <- DateRange[i]
    
    Fit.bearded[i,j] <- glmmLDTS.bearded.fit$fixed.effect[
      glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
      glmmLDTS.bearded.fit$fixed.effect[
        glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
          glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day +
      glmmLDTS.bearded.fit$fixed.effect[
        glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
          glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^2 +
      glmmLDTS.bearded.fit$fixed.effect[
        glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
          glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^3
    
    L <- rbind(L, (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "intercept")*1 +
                 (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
                    glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day +
                 (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
                    glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^2 +
                 (glmmLDTS.bearded.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
                    glmmLDTS.bearded.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^3)    
  }
}
L=L[-1,]
Ell <- L%*%glmmLDTS.bearded.fit$fixed.effect$estimate
Cell <- L %*% glmmLDTS.bearded.fit$covb %*% t(L)
#generate predictions on the logit scale and back transform
H.bearded=array(0,dim=c(dim(Fit.bearded),1000))
for(i in 1:1000){
  if(i%%10==0)cat(paste('iter ',i,'\n'))
  H.bearded[,,i]=matrix(rmvnorm(1,mean=Ell,sigma=Cell,method="svd"),nrow=length(DateRange),ncol=24,byrow=TRUE)
}
H.bearded=1/(1+exp(-H.bearded)) #back transform

# ------------------------------------------------------------------------------
#                   Spotted SEALS
# ------------------------------------------------------------------------------



Fit.spotted <- matrix(NA, nrow = length(DateRange), ncol = 24)
L=rep(NA,length(glmmLDTS.spotted.fit$fixed.effect[,"effect"]))
for(i in 1:length(DateRange)) {
  for (j in 1:24) {
    
    Hour <- HourRange[j]
    Day <- DateRange[i]
    
    Fit.spotted[i,j] <- glmmLDTS.spotted.fit$fixed.effect[
      glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "intercept", "estimate"]  +
      glmmLDTS.spotted.fit$fixed.effect[
        glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
          glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day +
      glmmLDTS.spotted.fit$fixed.effect[
        glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
          glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^2 +
      glmmLDTS.spotted.fit$fixed.effect[
        glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
          glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""),
        "estimate"]*Day^3
    
    L <- rbind(L, (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "intercept")*1 +
                 (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear" &
                    glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day +
                 (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear2" &
                    glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^2 +
                 (glmmLDTS.spotted.fit$fixed.effect[,"effect"] == "hourCat:DayofYear3" &
                    glmmLDTS.spotted.fit$fixed.effect[,"levels"] == paste(Hour,",", sep = ""))*Day^3)    
  }
}
L=L[-1,]
Ell <- L%*%glmmLDTS.spotted.fit$fixed.effect$estimate
Cell <- L %*% glmmLDTS.spotted.fit$covb %*% t(L)
#generate predictions on the logit scale and back transform
H.spotted=array(0,dim=c(dim(Fit.spotted),1000))
for(i in 1:1000){
  if(i%%10==0)cat(paste('iter ',i,'\n'))
  H.spotted[,,i]=matrix(rmvnorm(1,mean=Ell,sigma=Cell,method="svd"),nrow=length(DateRange),ncol=24,byrow=TRUE)
}
H.spotted=1/(1+exp(-H.spotted)) #back transform

Haulout.samples=list(spotted=H.spotted,bearded=H.bearded,ribbon=H.ribbon)
#output lookup table
save(Haulout.samples,file="c:/users/paul.conn/git/BOSS/BOSS/data/Haulout_samples.Rdat")
