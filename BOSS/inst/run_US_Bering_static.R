# script to estimate abundance of seals in the U.S. Bering, April 20-27
#
# Paul B. Conn
# 

  library(sp)
  library(rgeos)
  library(doBy)
  
  #First, pull in data
  load('c:/users/paul.conn/git/BOSS/BOSS/data/Transect_Data.Rdat') #read in effort data (held in "flight_segs" SpLinesDF)
  load('c:/users/paul.conn/git/BOSS/BOSS/data/AlaskaBeringData.Rdat') #read in "Data" holding grid info
  load('c:/users/paul.conn/git/BOSS/BOSS/data/Hotspot_Data.Rdat')  #read in hotspot data (held in 'hotspots' SpPointsDF)
  load('c:/users/paul.conn/git/BOSS/BOSS/data/p13.Rdata')  #read in confusion array
  load('c:/users/paul.conn/git/BOSS/BOSS/data/Haulout_samples.Rdat')  #read in confusion array
  load('c:/users/paul.conn/git/BOSS/BOSS/data/Effort_points.Rdat')  #read in Sp Points DF "effort_data" for on effort points (use for covariates)
  
  #do some cleaning up of hotspot records
  if(sum(is.na(hotspots[["hotspot_found"]]))>0)hotspots=hotspots[-which(is.na(hotspots[["hotspot_found"]])),] #remove records with a missing hotspot_found entry
  if(sum(is.na(hotspots[["hotspot_type"]]))>0)hotspots=hotspots[-which(is.na(hotspots[["hotspot_type"]])),]
  #temporarily get rid of records for which hotspot_type=seal but species=NA
  Tmp.seal=hotspots[["hotspot_type"]]
  Tmp.seal[is.na(Tmp.seal)]="NA"
  if(sum(is.na(hotspots[["species"]]) & Tmp.seal=="seal")>0)hotspots=hotspots[-which(is.na(hotspots[["species"]]) & Tmp.seal=="seal"),]
  #delete records for which species="o" 
  if(length(which(hotspots[["species"]]=="o"))>0)hotspots=hotspots[-which(hotspots[["species"]]=="o"),]
  
  #assemble a list of unique flight names and first/last on effort dates for each
  Flt.ids=unique(flight_segs[["flightid"]])
  Flt.table=data.frame(id=Flt.ids,min.time=rep(flight_segs[["min_dt"]][1],length(Flt.ids)),max.time=rep(flight_segs[["max_dt"]][1],length(Flt.ids)))
  for(i in 1:nrow(Flt.table)){
    Cur.which=which(flight_segs[["flightid"]]==Flt.table[i,"id"])
    Flt.table[i,"min.time"]=min(flight_segs[["min_dt"]][Cur.which])
    Flt.table[i,"max.time"]=max(flight_segs[["max_dt"]][Cur.which])
  }
  #select flights which began between April 20 and April 27
  min.ti=as.POSIXct("2012-04-20",tz="GMT")
  max.ti=as.POSIXct("2012-04-28",tz="GMT")
  Flt.anal=as.character(Flt.table[which(Flt.table[,"min.time"]>min.ti & Flt.table[,"min.time"]<max.ti),1])
  Flt.anal=Flt.anal[-which(Flt.anal=="AeroFl10")]  #didn't include this in our list of 10 flights for some reason
  Flights=flight_segs[flight_segs[["flightid"]]%in%Flt.anal,]
  #Flights=Flights[-which(Flights[["distance"]]==0),]
  plot(Flights)
  #output table of duration, distance by flight
  Fl.table=summaryBy(duration+distance~flightid,data=Flights@data,FUN=sum)
  Fl.table[,"distance.sum"]=Fl.table[,"distance.sum"]/1000
  write.csv(Fl.table,file="FLIR_effort_summary_power.csv")
  
  ##create a single grid object for this analysis, averaging time varying covariates over the dates considered
  Cur.grid=Data$Grid[[1]]
  S=nrow(Cur.grid)  #number of cells
  #first object corresponds to April 4 so start with 17th SpDF
  Mean.covs=matrix(0,S,3)
  for(i in 17:24)Mean.covs=Mean.covs+Data$Grid[[i]]@data[,5:7]    
  Mean.covs=Mean.covs/8
  Cur.grid@data[,5:7]=Mean.covs
  rownames(Cur.grid@data)=c(1:S)
  Cur.grid<-spChFIDs(Cur.grid,as.character(c(1:S)))  #change ID names to go from 1:S
  
  #turn land cover covariate into Poisson intensity modifier
  Area.hab=1-Cur.grid@data[,"land_cover"]
  
  #overlay flight tracks over grid cells, calculating cumulative distance flown in each
  int<-gIntersects(Flights,Cur.grid,byid=TRUE)
  vec <- vector(mode="list", length=nrow(Flights))
  for (i in seq(along=vec)) vec[[i]] <- try(gIntersection(Flights[i,],Cur.grid[int[,i],], byid=TRUE))
  out <- do.call("rbind", vec)
  rn <- row.names(out)
  nrn <- do.call("rbind", strsplit(rn, " "))
  Length.df <- data.frame(Fl=nrn[,1], poly=as.numeric(as.character(nrn[,2])), len=gLength(out,byid=TRUE))
  #calculate area surveyed for each of these flight_segment * grid cell combos
  Row.index=rep(0,nrow(Length.df))
  for(i in 1:length(Row.index))Row.index[i]=which(Flights[["seg_id"]]==as.character(Length.df[i,"Fl"]))
  Diam=Flights[["swath_diam"]][Row.index]
  DT=Flights[["min_dt"]][Row.index]+(Flights[["max_dt"]][Row.index]-Flights[["min_dt"]][Row.index])/2
  Length.df["day"]=as.POSIXlt(DT)$yday
  Length.df["day"]=Length.df["day"]-min(Length.df["day"])+1
  Length.df["hour"]=as.POSIXlt(DT)$hour+1
  
  
  #intersect on effort spatial points with grid
  effort_data=effort_data[effort_data[["flightid"]]%in%Flt.anal,]
  int<-gIntersects(effort_data,Cur.grid,byid=TRUE)
  effort_data=effort_data[apply(int,2,'sum')==1,]
  which.cell=function(x)which(x==1)
  Cell.id=unlist(apply(int,2,which.cell))
  Date.lt=as.POSIXlt(effort_data[["gps_dt"]])
  Day=Date.lt$mday-20 #make april 21st day 1
  Hour=Date.lt$hour+1
  Swath=effort_data[["diameter"]]
  
  #overwrite swath diameter using on effort spatial points if they overlap in a cell (flight segment values are used if all points data are corrupted for a given cell)
  for(i in 1:nrow(Length.df)){
    Which.pts=which(Cell.id==Length.df[i,"poly"])
    if(length(Which.pts)>0){
      Diam[i]=mean(Swath[which(Cell.id==Length.df$poly[i])])
      
    }
  }
  #total area surveyed in each flight * grid cell combo
  Length.df["area"]=Length.df[,"len"]*0.001*(Diam*.001)/(625*Area.hab[Length.df[,"poly"]])  #proportional area surveyed for each cell
  
  Mapping=as.numeric(as.character(unique(Length.df[,"poly"])))

  #sum total proportion of the habitable portion of each strata that is surveyed (for each surveyed cell)
  n.transects=length(Mapping)
  Area.trans=rep(0,n.transects)
  for(i in 1:n.transects)Area.trans[i]=sum(Length.df[Length.df["poly"]==Mapping[i],"area"])
  
  #determine a day and hour for each surveyed cell to match with 'Thin' (for time & date specific haulout correction)
  DayHour=data.frame(day=rep(0,n.transects),hour=rep(0,n.transects))
  for(i in 1:length(Mapping)){
    Which.pts=which(Cell.id==Mapping[i])
    if(length(Which.pts)>0){
      DayHour[i,"day"]=round(mean(Day[Which.pts]))
      Cur.hour=Hour[Which.pts]-1
      Hour.trans=2*pi*Cur.hour/24
      x.bar=mean(sin(Hour.trans))
      y.bar=mean(cos(Hour.trans))
      r=sqrt(x.bar^2+y.bar^2)
      if(x.bar>=0)alpha.bar=acos(y.bar/r)
      else alpha.bar=2*pi-acos(y.bar/r)
      hour.bar=round(12*alpha.bar/pi)  
      DayHour[i,"hour"]=hour.bar+1  #average using circular statistics so that the mean of hour 1 and hour 23 is 0 not 12
    }
    else{ #for cells that don't have an on effort point because of data corruption
      DayHour[i,"day"]=round(mean(Length.df[which(Length.df[,"poly"]==Mapping[i]),"day"]))
      Cur.hour=Length.df[which(Length.df[,"poly"]==Mapping[i]),"hour"]
      Hour.trans=2*pi*Cur.hour/24
      x.bar=mean(sin(Hour.trans))
      y.bar=mean(cos(Hour.trans))
      r=sqrt(x.bar^2+y.bar^2)
      if(x.bar>0)alpha.bar=acos(y.bar/r)
      else alpha.bar=2*pi-acos(y.bar/r)
      hour.bar=round(12*alpha.bar/pi)  
      DayHour[i,"hour"]=hour.bar+1  #average using circular statistics so that the mean of hour 1 and hour 23 is 0 not 12
    }
  } 
  if(sum(DayHour[,"hour"]==25)>0)DayHour[DayHour[,"hour"]==25,"hour"]=1
  #Thin.pl=(DayHour[,"hour"]-1)*8+DayHour[,"day"] #maybe easiest to just one index
  
  Adj=Data$Adj
  rm(Data)
    
  #Format hotspot data to be compatible with R BOSS package for analysis
  
  #1) attach a transect to each hotspot
  Tmp=gIntersects(hotspots,Cur.grid,byid=TRUE)
  #remove several hotspots that weren't on the grid
  I.intersects=apply(Tmp,2,'sum')
  hotspots=hotspots[-which(I.intersects==0),]
  Tmp=Tmp[,-which(I.intersects==0)]
  fun<-function(x)which(x==1)
  Transects=apply(Tmp,2,'fun')
  #associate points that don't occur in surveyed cells with closest surveyed cell
  I.no.survey=1-(Transects %in% Mapping)
  if(sum(I.no.survey)>0){
    Which.no.survey=which(I.no.survey==TRUE)
    for(i in 1:sum(I.no.survey)){
      Cur.dist=gDistance(hotspots[Which.no.survey[i],],Cur.grid[Mapping,],byid=TRUE)
      Transects[Which.no.survey[i]]=Mapping[which(Cur.dist==min(Cur.dist))]
    }
  }
  #associate each point with a "transect" [vector index for the Mapping object]
  fun<-function(x,Mapping)which(Mapping==x)
  Transects=unlist(sapply(Transects,'fun',Mapping=Mapping))
  
  #2) assign a numeric value for 'hotspot_found'
  Photo=ifelse(hotspots[["hotspot_found"]]=="yes",1,0)
  
  #3) assign a numeric species value (with max integer being unknown)
  Sp=hotspots[["species"]]
  Sp[Sp=="bd"]=1
  Sp[Sp=="rn"]=2
  Sp[Sp=="rd"]=3
  Sp[Sp=="sd"]=4
  Sp[which(hotspots[["hotspot_type"]] %in% c("other_animal","dirty_ice_anomaly","dirty_ice_anom","unknown","o"))]=5
  Sp[Sp=="unk"]=6
  Sp[which(hotspots[["hotspot_type"]]=="seal_evidence")]=6
  Sp=as.numeric(Sp)
  
  #4) format observer uncertainty into certain, likely, guess
  Cert=rep(1,length(Sp))
  Cert[hotspots[["species_conf"]]=="likely"]=2
  Cert[hotspots[["species_conf"]]=="guess"]=3
  #change any 'likely' or 'guess' of 'unknown' species to 'certain' (no unknown unknowns!)
  Cert[which(Sp==6)]=1
  #change all anomalies to have certainty 1
  Cert[which(Sp==5)]=1
  
  #5) get group numbers right
  Grp=rep(1,length(Sp))
  Grp[which(hotspots[["hotspot_type"]]=="seal")]=as.numeric(hotspots[["numseals"]])[which(hotspots[["hotspot_type"]]=="seal")]
  
  #6) assemble final data.frame
  Dat=data.frame(Transect=Transects,Photo=Photo,Obs=3*(Sp-1)+Cert,Grp=Grp)
  
  #Do some cleaning up
  rm(Cert,Diam,DT,effort_data,flight_segs,Flights,Flt.anal,Fl.table,Flt.ids,Flt.table,hotspots,Hour,hour.bar,Hour.trans,int,Photo,r,rn,nrn,Which.pts,Swath,Day,Cell.id,Tmp)
  
  library(hierarchicalDS)
	source('c:/users/paul.conn/git/BOSS/BOSS/R/util_funcs.R')
	source('c:/users/paul.conn/git/BOSS/BOSS/R/hierarchical_boss.R')
	source('c:/users/paul.conn/git/BOSS/BOSS/R/mcmc_boss.R')
	
	set.seed(22222) 
	ZIP=FALSE 
  misID=TRUE
	spat.ind=TRUE #do not make spatially independent; i.e. estimate spatial autocorrelation!  
  
  #Reformat Psi array using Brett's misID data (species and certainty are both in different order)
  n.species=5
  Psi=array(0,dim=c(n.species,(n.species+1)*3,dim(p13)[3]))
  Psi[1,1:3,]=p13[3,9:7,]
  Psi[1,4:6,]=p13[3,6:4,]
  Psi[1,7:9,]=p13[3,12:10,]
  Psi[1,10:12,]=p13[3,3:1,]
  Psi[1,16,]=p13[3,13,]
  Psi[2,1:3,]=p13[2,9:7,]
  Psi[2,4:6,]=p13[2,6:4,]
  Psi[2,7:9,]=p13[2,12:10,]
  Psi[2,10:12,]=p13[2,3:1,]
  Psi[2,16,]=p13[2,13,]
  Psi[3,1:3,]=p13[4,9:7,]
  Psi[3,4:6,]=p13[4,6:4,]
  Psi[3,7:9,]=p13[4,12:10,]
  Psi[3,10:12,]=p13[4,3:1,]
  Psi[3,16,]=p13[4,13,]
  Psi[4,1:3,]=p13[1,9:7,]
  Psi[4,4:6,]=p13[1,6:4,]
  Psi[4,7:9,]=p13[1,12:10,]
  Psi[4,10:12,]=p13[1,3:1,]
  Psi[4,16,]=p13[1,13,]
  Psi[5,13,]=1
  rm(p13)
  
  #generate thinning priors using haulout and det prob data
  Thin=array(1,dim=c(n.species,8,24,1000))
  P=rbeta(1000,67,5)  #conjugate beta(1,1) for binomial detection data (66/70 successes)
  for(isp in 1:4){
    for(iday in 1:8){
      for(ihr in 1:24){
        Thin[isp,iday,ihr,]=P
        if(isp==1)Thin[isp,iday,ihr,]=Thin[isp,iday,ihr,]*Haulout.samples$bearded[iday,ihr,]
        if(isp==2)Thin[isp,iday,ihr,]=Thin[isp,iday,ihr,]*Haulout.samples$ribbon[iday,ihr,]
        if(isp==4)Thin[isp,iday,ihr,]=Thin[isp,iday,ihr,]*Haulout.samples$spotted[iday,ihr,]
      }
    }
  }
  rm(Haulout.samples)  
  
  #2) declare inputs and call hierarchical model; 
  #Obs.cov=matrix(0,n.transects,1)
	n.obs.cov=0 #1 observer covariate, "Seat" is included in the dataset even though it isn't modeled
  Prop.photo=rep(NA,n.transects)  #proportion of surveyed area in each transect that is photographed (used in post. loss calcs) NOT IMPLEMENTED YET
	Hab.cov=Cur.grid@data #covariate on abundance intensity same as used to generate data
  Hab.cov2=Hab.cov^2
  paste_fun=function(x)paste(x,2,sep='')
  colnames(Hab.cov2)=sapply(colnames(Hab.cov),paste_fun)
  Hab.cov=cbind(Hab.cov,Hab.cov2)
  Hab.pois.formula=list("vector",n.species)
  for(i in 1:n.species)Hab.pois.formula[[i]]=~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_contour+dist_edge #the density of species two isn't affected by the habitat covariate but we'll estimate this effect anyway
  Hab.bern.formula=c(~1,~1,~1,~1,~1)  #formula for Bernoulli part of ZIP model
	Cov.prior.parms=array(0,dim=c(n.species,2,1))
	Cov.prior.parms[,1,1]=0.1  
	Cov.prior.parms[,2,1]=0.1
  #Cov.prior.parms[,1,1]=c(0.1,0.2,0.3,0.4,0)
	Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
	Cov.prior.pdf=matrix(0,n.species,1)  
	Cov.prior.pdf[,1]="pois1"  #model group size as a zero truncated poisson
	Cov.prior.n=matrix(2,n.species,1)
  fix.tau.nu=TRUE
	srr=TRUE
	srr.tol=0.95
	grps=TRUE
  post.loss=FALSE 
  Control=list(iter=110000,burnin=10000,thin=20,MH.nu=matrix(.2,n.species,S),adapt=300)
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
	set.seed(8327329)   #chain1
	Out=hierarchical_boss(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,DayHour=DayHour,Thin=Thin,Prop.photo=Prop.photo,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.pois.formula=Hab.pois.formula,Hab.bern.formula=Hab.bern.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,ZIP=ZIP,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,Inits=Inits,grps=grps,n.species=n.species,Control=Control,adapt=adapt,Prior.pars=Prior.pars,Psi=Psi,post.loss=post.loss)

  #Post hoc correction for density at sea concentration = 0; note: will need to adjust N in MCMC object too
  Out$Post$N[,,which(Hab.cov[,"ice_conc"]<.0001)]=0
  
  ###3) plot and summarize results; note that chain would need to be run a lot longer to summarize the posterior very well!!!
  #plot(Out$MCMC)
	#summary_N(Out)
  #post_loss(Out)
	Tmp<-Cur.grid
  iter.start=1
  New.dat=matrix(0,S,n.species)
  for(isp in 1:n.species){
    New.dat[,isp]=apply(Out$Post$N[isp,iter.start:(dim(Out$Post$N)[2]),],2,'median')
  }
  colnames(New.dat)=c("bearded","ribbon","ringed","spotted","other")
  Tmp@data=cbind(Cur.grid@data,New.dat)
  library(ggplot2)
  Tmp@data$id=rownames(Tmp@data)
  tmp1<-fortify(Tmp,region='id')
  tmp2<-join(tmp1,Tmp@data,by="id")
  #qplot(long, lat, data = tmp2, fill = bearded, geom = "raster")
  ggplot(tmp2)+aes(long,lat,fill=bearded)+geom_raster()
  ggplot(tmp2)+aes(long,lat,fill=ribbon)+geom_raster()
  ggplot(tmp2)+aes(long,lat,fill=ringed)+geom_raster()
  ggplot(tmp2)+aes(long,lat,fill=spotted)+geom_raster()
  ggplot(tmp2)+aes(long,lat,fill=other)+geom_raster()
#   pdf('dist_mainland.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=dist_mainland)+geom_raster()
#   dev.off()
#   pdf('dist_land.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=dist_land)+geom_raster()
#   dev.off()
#   pdf('dist_shelf.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=dist_shelf)+geom_raster()
#   dev.off()
#   pdf('ice_conc.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=ice_conc)+geom_raster()
#   dev.off()
#   pdf('dist_contour.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=dist_contour)+geom_raster()
#   dev.off()
#   pdf('dist_edge.pdf',height=2.9, width=5)    
#   ggplot(tmp2)+aes(long,lat,fill=dist_edge)+geom_raster()
#   dev.off()

