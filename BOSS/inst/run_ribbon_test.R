# script to compare ribbon seal estimation with that from a straight log-Gaussian Cox process (in spatial prediction)
#
# Paul B. Conn
# 

library(sp)
library(rgeos)
library(doBy)

load('c:/users/paul.conn/git/spatpred/Ribbon_data.Rda')

#also include quadratic term for ice_conc in dataset
ice_conc2=(Cur.grid@data[["ice_conc"]])^2
Cur.grid@data=cbind(Cur.grid@data,ice_conc2)
#attach Ribbon seal counts to data frame
#Cur.grid@data=cbind(Cur.grid@data,Count)
Which.not.sampled=c(1:S)
Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]

#add in some zero counts in cells with no ice
Which.no.ice=which(Cur.grid@data[,"ice_conc"]<0.001)
Which.no.ice=Which.no.ice[which((Which.no.ice%in%Which.not.sampled)==TRUE)]
Mapping=c(Mapping,Which.no.ice)
Area.trans=c(Area.photo,rep(0.1,length(Which.no.ice)))
Count=Ribbon.count[Mapping]

Dat=data.frame(Transect=rep(1,sum(Count)),Photo=rep(1,sum(Count)),Obs=rep(1,sum(Count)),Grp=rep(1,sum(Count)))
irow=1
for(iobs in 1:length(Count)){
  if(Count[iobs]>0){
    for(iind in 1:Count[iobs]){
      Dat[irow,]=c(iobs,1,1,1)
      irow=irow+1
    }
  }
}

n.transects=length(Mapping)

Hab.cov=Cur.grid@data


#8) Reformat Psi array using Brett's misID data (species and certainty are both in different order)
library(hierarchicalDS)
source('c:/users/paul.conn/git/BOSS/BOSS/R/util_funcs.R')
source('c:/users/paul.conn/git/BOSS/BOSS/R/hierarchical_boss.R')
source('c:/users/paul.conn/git/BOSS/BOSS/R/mcmc_boss.R')

n.species=1
photo.adjust=1  #change to photo=0 with this probability (for examining power with a reduced proportion photographs)
Which.photo=which(Dat[,"Photo"]==1)
New.photo=rbinom(length(Which.photo),1,photo.adjust)
Dat[Which.photo,"Photo"]=New.photo

set.seed(22222) 
ZIP=FALSE 
misID=FALSE
spat.ind=FALSE #do not make spatially independent; i.e. estimate spatial autocorrelation!  


#generate thinning priors using haulout and det prob data
Thin=array(1,dim=c(1,8,24,1000))
DayHour=matrix(1,n.transects,2)
colnames(DayHour)=c("day","hour")
Psi=array(1,dim=c(1,1,1000))


#2) declare inputs and call hierarchical model; 
#Obs.cov=matrix(0,n.transects,1)
n.obs.cov=0 
Prop.photo=rep(1,n.transects)  #proportion of surveyed area in each transect that is photographed (used in post. loss calcs) 
Obs.cov=NULL
Hab.cov=Cur.grid@data #covariate on abundance intensity same as used to generate data
Hab.pois.formula=list("vector",1)
for(i in 1:1)Hab.pois.formula[[i]]=~1
Hab.bern.formula=c(~1,~1,~1,~1,~1)  #formula for Bernoulli part of ZIP model
Cov.prior.parms=array(0,dim=c(n.species,2,1))
Cov.prior.parms[,1,1]=0.1  
Cov.prior.parms[,2,1]=0.1
#Cov.prior.parms[,1,1]=c(0.1,0.2,0.3,0.4,0)
Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
Cov.prior.pdf=matrix(0,n.species,1)  
Cov.prior.pdf[,1]="pois1"  #model group size as a zero truncated poisson
Cov.prior.n=matrix(2,n.species,1)
fix.tau.nu=FALSE
srr=TRUE
srr.tol=0.8
grps=TRUE
post.loss=FALSE
#Control=list(iter=600000,burnin=100000,thin=250,MH.nu=matrix(.2,n.species,S),adapt=400)
Control=list(iter=20000,burnin=10000,thin=1,MH.nu=matrix(.2,n.species,S),adapt=1000)
#hab.pois=matrix(0,n.species,2) #covariates are intercept, index
#hab.pois[,1]=log(Sim$G.tot/S+10) #start 'near' true value
#hab.bern=matrix(0,n.species,1)
#hab.bern[,1]=0.5
#hab.pois[1,]=c(log(39),1) 
#hab.pois[2,]=c(log(9),0)
#hab.pois[3,]=c(log(20),-1)
#hab.pois[4,]=c(log(10),0.5)
#hab.pois[5,]=c(log(15),0)
#hab.bern=matrix(.5,1,1)
Inits=list(tau.nu=rep(100,n.species))
#Inits=list(hab.pois=hab.pois,hab.bern=hab.bern,tau.nu=rep(100,n.species)) #provide some initial values to ensure MCMC doesn't start out at weird place
Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.tau=0.01) #(1,.01) prior makes it closer to a uniform distribution near the origin
adapt=TRUE
set.seed(12345)   #chain1
Out=hierarchical_boss(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,DayHour=DayHour,Thin=Thin,Prop.photo=Prop.photo,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.pois.formula=Hab.pois.formula,Hab.bern.formula=Hab.bern.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,ZIP=ZIP,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,Inits=Inits,grps=grps,n.species=n.species,Control=Control,adapt=adapt,Prior.pars=Prior.pars,Psi=Psi,post.loss=post.loss)

save(Out,file='USBering_ribbon_test.Rdat')

#png('ice_eff2.png')
#plot_covar(DM=Out$DM.hab.pois,MCMC=Out$MCMC,Vars=c("ice_conc"),n.species=n.species,n.points=20,Sp.names=c("Bearded","Ribbon","Ringed","Spotted","Other"),const.tau=100,bern=FALSE)
#dev.off()

Tmp.grid=Cur.grid
Tmp<-Tmp.grid
iter.start=1
New.dat=matrix(0,S,n.species)
for(isp in 1:n.species){
  New.dat[,isp]=apply(Out$Post$N[isp,iter.start:(dim(Out$Post$N)[2]),],2,'mean')
}
Tmp@data=cbind(Tmp.grid@data,New.dat)
library(ggplot2)
library(plyr)
library(grid)
Tmp@data$id=rownames(Tmp@data)
tmp1<-fortify(Tmp,region='id')
tmp2<-join(tmp1,Tmp@data,by="id")
new.colnames=colnames(tmp2)
new.colnames[1:2]=c("Easting","Northing")
colnames(tmp2)=new.colnames
#qplot(long, lat, data = tmp2, fill = bearded, geom = "raster")
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
#tmp2[which(tmp2[,"New.dat"]>400),"New.dat"]=400
p1=ggplot(tmp2)+aes(Easting,Northing,fill=New.dat)+geom_raster()+tmp.theme
p1

plot(Out$MCMC[,1])
#   pdf('dist_contour.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=dist_contour)+geom_raster()
#   dev.off()
#   pdf('dist_edge.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=dist_edge)+geom_raster()
#   dev.off()
