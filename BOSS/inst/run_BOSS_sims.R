#' script to simulate distance sampling data (with misID) and run a simple hierarchical example analysis using said dataset
#' #@return see help for hierarchical_DS.R
#' #@export
#' #@keywords simulation
#' #@author Paul B. Conn
#' 
  #1) First, simulate data
	S=900 #Number of grid cells; this needs to be a square number (square grid assumed)
  library(hierarchicalDS)
	source('c:/users/paul.conn/git/BOSS/BOSS/R/util_funcs.R')
  source('c:/users/paul.conn/git/BOSS/BOSS/inst/simulate_data.R')
	source('c:/users/paul.conn/git/BOSS/BOSS/R/hierarchical_boss.R')
	source('c:/users/paul.conn/git/BOSS/BOSS/R/mcmc_boss.R')
	source('c:/users/paul.conn/git/hierarchicalDS/hierarchicalDS/R/spat_funcs.R')
	
	n.transects=200 #one transect per cell
  n.species=5
	set.seed(11111) 
	ZIP=FALSE 
  misID=TRUE
	spat.ind=FALSE #do not make spatially independent; i.e. estimate spatial autocorrelation!
  tau.pois=10000
	if(spat.ind==FALSE)tau.pois=20
  tau.bern=10000
  #if(spat.ind==FALSE)tau.bern=20
  prop.photo=0.8
	Sim=simulate_BOSS(S=S,prop.photo=prop.photo,n.sampled=n.transects,misID=misID,ZIP=ZIP,tau.pois=tau.pois,tau.bern=tau.bern)
  #plot 'true' abundance by species
	Plot.dat=data.frame(cbind(rep(c(1:sqrt(S)),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),t(Sim$G.true)))
  colnames(Plot.dat)=c("Northing","Easting","Species1","Species2","Species3","Species4","Species5")  
	library(ggplot2)
	library(plyr)
	library(grid)
	pdf(file="sim_maps.pdf")
	pushViewport(viewport(layout=grid.layout(3,2)))
	tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
	p1=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species1)+geom_raster()+tmp.theme
	print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
	p2=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species2)+geom_raster()+tmp.theme
	print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
	p3=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species3)+geom_raster()+tmp.theme
	print(p3,vp=viewport(layout.pos.row=2,layout.pos.col=1))
	p4=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species4)+geom_raster()+tmp.theme
	print(p4,vp=viewport(layout.pos.row=2,layout.pos.col=2))
	p5=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species5)+geom_raster()+tmp.theme
	print(p5,vp=viewport(layout.pos.row=3,layout.pos.col=1))
	dev.off()

  Dat=Sim$Dat
  Psi=array(Sim$Psi,dim=c(nrow(Sim$Psi),ncol(Sim$Psi),20))
  for(i in 2:20){
    Psi[,,i]=Psi[,,1]*exp(rnorm(prod(dim(Psi)[1:2]),0,0.02)) #small amount of lognormal variability around true relationship for testing
    Psi[,,i]=Psi[,,i]/apply(Psi[,,i],1,'sum')
  }
  tabulate(Dat[,"Species"])
  #2) declare inputs and call hierarchical model; 
  #Obs.cov=matrix(0,n.transects,1)
	n.obs.cov=0 #1 observer covariate, "Seat" is included in the dataset even though it isn't modeled
	Adj=rect_adj_RW2(sqrt(S),sqrt(S))
	Mapping=Sim$Sampled.cells
	Area.trans=rep(0.1,n.transects)
	Area.hab=rep(1,S)
  Prop.photo=rep(prop.photo,n.transects)  #proportion of surveyed area in each transect that is photographed (used in post. loss calcs)
	DayHour=data.frame(day=rep(1,n.transects),hour=rep(1,n.transects))
	Thin=array(1,dim=c(n.species,2,2,1000))	#need >1 day, hour for matrix subscripting to work right in mcmc_boss
  Hab.cov=data.frame(Sim$X[,2:ncol(Sim$X)]) #covariate on abundance intensity same as used to generate data
  colnames(Hab.cov)=c("Easting","Northing","Matern.cov")
	Hab.pois.formula=c(~Easting+Northing+Matern.cov,~Easting+Northing+Matern.cov,~Easting+Northing+Matern.cov,~Easting+Northing+Matern.cov,~Easting+Northing+Matern.cov) #the density of species two isn't affected by the habitat covariate but we'll estimate this effect anyway
  Hab.bern.formula=c(~1,~1,~1,~1,~1)  #formula for Bernoulli part of ZIP model
	Cov.prior.parms=array(0,dim=c(n.species,2,1))
	Cov.prior.parms[,1,1]=2  #we'll put priors a little off from their true values; expected group sizes are 4 and 2 for each species
	Cov.prior.parms[,2,1]=1
  #Cov.prior.parms[,1,1]=c(0.1,0.2,0.3,0.4,0)
	Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
	Cov.prior.pdf=matrix(0,n.species,1)  
	Cov.prior.pdf[,1]="pois1"  #model group size as a zero truncated poisson
	Cov.prior.n=matrix(2,n.species,1)
  fix.tau.nu=TRUE
	srr=TRUE
	srr.tol=50
	grps=TRUE
  post.loss=FALSE 
  Control=list(iter=600000,burnin=100000,thin=100,MH.nu=matrix(.2,n.species,S),adapt=400)
  hab.pois=matrix(0,n.species,4) #covariates are intercept, index
	hab.pois[,1]=log(Sim$G.tot/S+10) #start 'near' true value
  #hab.bern=matrix(0,n.species,1)
  #hab.bern[,1]=0.5
	#hab.pois[1,]=c(log(39),1) 
	#hab.pois[2,]=c(log(9),0)
	#hab.pois[3,]=c(log(20),-1)
	#hab.pois[4,]=c(log(10),0.5)
	#hab.pois[5,]=c(log(15),0)
	hab.bern=matrix(c(1,.5,0,.4,0),5,1)
	
	Inits=list(hab.pois=hab.pois,hab.bern=hab.bern,tau.nu=rep(500,n.species)) #provide some initial values to ensure MCMC doesn't start out at weird place
	Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.tau=0.01) #(1,.01) prior makes it closer to a uniform distribution near the origin
	adapt=TRUE	
	set.seed(8327329)   #chain1
	Out=hierarchical_boss(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Prop.photo=Prop.photo,DayHour=DayHour,Thin=Thin,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.species=n.species,n.obs.cov=n.obs.cov,Hab.pois.formula=Hab.pois.formula,Hab.bern.formula=Hab.bern.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,ZIP=ZIP,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,Inits=Inits,grps=grps,Control=Control,adapt=adapt,Prior.pars=Prior.pars,Psi=Psi,post.loss=post.loss)
  save.image("sim_Euring2.rdat")
  ###3) plot and summarize results; note that chain would need to be run a lot longer to summarize the posterior very well!!!
  #plot(Out$MCMC)
	#summary_N(Out)
  #post_loss(Out)
  
	N.tot=matrix(0,n.species,dim(Out$Post$N)[2])
	for(iiter in 1:dim(Out$Post$N)[2])N.tot[,iiter]=apply(Out$Post$N[,iiter,],1,'sum')
	par(mfrow=c(2,2))
	hist(N.tot[1,],main='Species 1',xlab='',freq=FALSE,ylab='Posterior density',cex.lab=1.3,cex.axis=1.3)
  abline(v=Sim$G.tot[1]*1.1,col='red',lwd=2)
	hist(N.tot[2,],main='Species 2',xlab='',freq=FALSE,ylab='Posterior density',cex.lab=1.3,cex.axis=1.3)
	abline(v=Sim$G.tot[2]*1.2,col='red',lwd=2)
  hist(N.tot[3,],main='Species 3',xlab='Abundance',freq=FALSE,ylab='Posterior density',cex.lab=1.3,cex.axis=1.3)
	abline(v=Sim$G.tot[3]*1.3,col='red',lwd=2)
  hist(N.tot[4,],main='Species 4',xlab='Abundance',freq=FALSE,ylab='Posterior density',cex.lab=1.3,cex.axis=1.3)
	abline(v=Sim$G.tot[4]*1.4,col='red',lwd=2)
	
  
  iter.start=1
	New.dat=matrix(0,S,n.species)
	for(isp in 1:n.species){
	  New.dat[,isp]=apply(Out$Post$N[isp,iter.start:(dim(Out$Post$N)[2]),],2,'mean')
	}
	New.dat=cbind(New.dat,c(1:S),rep(c(1:sqrt(S)),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)))
  I.sample=rep(0,S)
  I.sample[Mapping]=1
  New.dat=cbind(New.dat,I.sample)
	colnames(New.dat)=c("Species1","Species2","Species3","Species4","Species5","id","Northing","Easting","Sampled")
  tmp2=as.data.frame(New.dat)
	Sampled.df=tmp2[tmp2[,"Sampled"]==1,]
	tmp2[,"Sampled"]=factor(tmp2[,"Sampled"])
  
  library(ggplot2)
	library(plyr)
	library(grid)
  
  #png(file="sim_sampled.png")
	theme_xy=theme(axis.ticks = element_blank(), axis.text = element_blank(),panel.grid.major = element_blank(),
	                panel.grid.minor = element_blank(),panel.background = element_blank(),legend.key.size=unit(0.5,'cm'))
	theme_y=theme(axis.ticks = element_blank(), axis.text = element_blank(),panel.grid.major = element_blank(),
	               panel.grid.minor = element_blank(),panel.background = element_blank(),legend.key.size=unit(0.5,'cm'),axis.title.x=element_blank())
	theme_x=theme(axis.ticks = element_blank(), axis.text = element_blank(),panel.grid.major = element_blank(),
	               panel.grid.minor = element_blank(),panel.background = element_blank(),legend.key.size=unit(0.5,'cm'),axis.title.y=element_blank())
	theme_no=theme(axis.ticks = element_blank(), axis.text = element_blank(),panel.grid.major = element_blank(),
	               panel.grid.minor = element_blank(),panel.background = element_blank(),legend.key.size=unit(0.5,'cm'),axis.title.x=element_blank(),axis.title.y=element_blank())
  p0=ggplot(tmp2)+aes(Easting,Northing,fill=Sampled)+geom_raster()+scale_fill_manual(breaks=c("0","1"),labels=c("N","Y"),values=c("grey","red"))+tmp.theme
  #p0
  #dev.off()
  
  png(file="sim_estimated2.png")
	pushViewport(viewport(layout=grid.layout(5,2)))
	p1=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species1)+geom_raster()+scale_fill_gradient2(limits=c(0,90),midpoint=45,low="white",mid="blue",high="black")+theme_xy
	print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
	p2=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species2)+geom_raster()+scale_fill_gradient2(limits=c(0,30),midpoint=15,low="white",mid="blue",high="black")+theme_xy
  print(p2,vp=viewport(layout.pos.row=2,layout.pos.col=1))
	p3=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species3)+geom_raster()+scale_fill_gradient2(limits=c(0,150),midpoint=75,low="white",mid="blue",high="black")+theme_xy
	print(p3,vp=viewport(layout.pos.row=3,layout.pos.col=1))
	p4=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species4)+geom_raster()+scale_fill_gradient2(limits=c(0,180),midpoint=90,low="white",mid="blue",high="black")+theme_xy
	print(p4,vp=viewport(layout.pos.row=4,layout.pos.col=1))
	p5=ggplot(Plot.dat)+aes(Easting,Northing,fill=Species5)+geom_raster()+scale_fill_gradient2(limits=c(0,70),midpoint=35,low="white",mid="blue",high="black")+theme_xy
	print(p5,vp=viewport(layout.pos.row=5,layout.pos.col=1))

  #note the following geom_point can simply be changed to geom_rect for outlines but looks pretty clumped
	p6=ggplot(tmp2)+aes(Easting,Northing,fill=Species1)+geom_raster()+scale_fill_gradient2(limits=c(0,90),midpoint=45,low="white",mid="blue",high="black")+theme_xy
	p6=p6+geom_point(data=Sampled.df,size=1,fill=NA,colour="red",aes(xmin=Easting-0.5,xmax=Easting+0.5,ymin=Northing-0.5,ymax=Northing+0.5))	
	print(p6,vp=viewport(layout.pos.row=1,layout.pos.col=2))
	p7=ggplot(tmp2)+aes(Easting,Northing,fill=Species2)+geom_raster()+scale_fill_gradient2(limits=c(0,30),midpoint=15,low="white",mid="blue",high="black")+theme_xy
	p7=p7+geom_point(data=Sampled.df,size=1,fill=NA,colour="red",aes(xmin=Easting-0.5,xmax=Easting+0.5,ymin=Northing-0.5,ymax=Northing+0.5))
	print(p7,vp=viewport(layout.pos.row=2,layout.pos.col=2))
	p8=ggplot(tmp2)+aes(Easting,Northing,fill=Species3)+geom_raster()+scale_fill_gradient2(limits=c(0,150),midpoint=75,low="white",mid="blue",high="black")+theme_xy
	p8=p8+geom_point(data=Sampled.df,size=1,fill=NA,colour="red",aes(xmin=Easting-0.5,xmax=Easting+0.5,ymin=Northing-0.5,ymax=Northing+0.5))
	print(p8,vp=viewport(layout.pos.row=3,layout.pos.col=2))
	p9=ggplot(tmp2)+aes(Easting,Northing,fill=Species4)+geom_raster()+scale_fill_gradient2(limits=c(0,180),midpoint=90,low="white",mid="blue",high="black")+theme_xy
	p9=p9+geom_point(data=Sampled.df,size=1,fill=NA,colour="red",aes(xmin=Easting-0.5,xmax=Easting+0.5,ymin=Northing-0.5,ymax=Northing+0.5))
	print(p9,vp=viewport(layout.pos.row=4,layout.pos.col=2))
	p10=ggplot(tmp2)+aes(Easting,Northing,fill=Species5)+geom_raster()+scale_fill_gradient2(limits=c(0,70),midpoint=35,low="white",mid="blue",high="black")+theme_xy
	p10=p10+geom_point(data=Sampled.df,size=1,fill=NA,colour="red",aes(xmin=Easting-0.5,xmax=Easting+0.5,ymin=Northing-0.5,ymax=Northing+0.5))
	print(p10,vp=viewport(layout.pos.row=5,layout.pos.col=2))
	dev.off()
	

