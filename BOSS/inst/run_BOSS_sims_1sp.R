#' script to simulate distance sampling data (with misID) and run a simple hierarchical example analysis using said dataset
#' #@return see help for hierarchical_DS.R
#' #@export
#' #@keywords simulation
#' #@author Paul B. Conn
#' 
  #1) First, simulate data
	S=100 #Number of grid cells; this needs to be a square number (square grid assumed)
  library(hierarchicalDS)
	source('c:/users/paul.conn/git/BOSS/BOSS/R/util_funcs.R')
  source('c:/users/paul.conn/git/BOSS/BOSS/inst/simulate_data_1sp.R')
	source('c:/users/paul.conn/git/BOSS/BOSS/R/hierarchical_boss.R')
	source('c:/users/paul.conn/git/BOSS/BOSS/R/mcmc_boss.R')
	
	n.transects=100 #one transect per cell
	set.seed(22222) 
	ZIP=TRUE 
  misID=FALSE
	spat.ind=TRUE #do not make spatially independent; i.e. estimate spatial autocorrelation!
  tau.pois=10000
	if(spat.ind==FALSE)tau.pois=15
  tau.bern=10000
  if(spat.ind==FALSE)tau.bern=20
  prop.photo=0.5
	Sim=simulate_BOSS(S=S,prop.photo=prop.photo,n.sampled=n.transects,misID=misID,ZIP=ZIP,tau.pois=tau.pois,tau.bern=tau.bern)
	Dat=Sim$Dat
  
  Psi=NULL
  #2) declare inputs and call hierarchical model; 
  #Obs.cov=matrix(0,n.transects,1)
	n.obs.cov=0 #1 observer covariate, "Seat" is included in the dataset even though it isn't modeled
	Adj=square_adj(sqrt(S))
	Mapping=Sim$Sampled.cells
	Area.trans=rep(0.1,n.transects)
	Area.hab=rep(1,S)
  Prop.photo=rep(prop.photo,n.transects)  #proportion of surveyed area in each transect that is photographed (used in post. loss calcs)
	Hab.cov=data.frame(rep(log(c(1:sqrt(S)/sqrt(S))),each=sqrt(S))) #covariate on abundance intensity same as used to generate data
  colnames(Hab.cov)=c("Cov1")
	Hab.pois.formula=c(~Cov1,~Cov1,~Cov1,~Cov1,~Cov1) #the density of species two isn't affected by the habitat covariate but we'll estimate this effect anyway
  Hab.bern.formula=c(~1,~1,~1,~1,~1)  #formula for Bernoulli part of ZIP model
	n.species=1
	Cov.prior.parms=array(0,dim=c(n.species,2,1))
	Cov.prior.parms[,1,1]=2  #we'll put priors a little off from their true values; expected group sizes are 4 and 2 for each species
	Cov.prior.parms[,2,1]=1
  #Cov.prior.parms[,1,1]=c(0.1,0.2,0.3,0.4,0)
	Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
	Cov.prior.pdf=matrix(0,n.species,1)  
	Cov.prior.pdf[,1]="pois1"  #model group size as a zero truncated poisson
	Cov.prior.n=matrix(2,n.species,1)
  fix.tau.nu=FALSE
	srr=TRUE
	srr.tol=0.5
	grps=TRUE
  post.loss=FALSE 
  Control=list(iter=420000,burnin=20000,thin=100,MH.nu=matrix(.2,n.species,S),adapt=400)
  hab.pois=matrix(0,n.species,2) #covariates are intercept, index
	hab.pois[,1]=log(Sim$G.tot/S+10) #start 'near' true value
  #hab.bern=matrix(0,n.species,1)
  #hab.bern[,1]=0.5
	#hab.pois[1,]=c(log(39),1) 
	#hab.pois[2,]=c(log(9),0)
	#hab.pois[3,]=c(log(20),-1)
	#hab.pois[4,]=c(log(10),0.5)
	#hab.pois[5,]=c(log(15),0)
	hab.bern=matrix(.5,1,1)
	
	Inits=list(hab.pois=hab.pois,hab.bern=hab.bern,tau.nu=rep(100,n.species)) #provide some initial values to ensure MCMC doesn't start out at weird place
	Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.tau=0.01) #(1,.01) prior makes it closer to a uniform distribution near the origin
	adapt=TRUE	
	set.seed(8327329)   #chain1
	Out=hierarchical_boss(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Prop.photo=Prop.photo,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.pois.formula=Hab.pois.formula,Hab.bern.formula=Hab.bern.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,ZIP=ZIP,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,Inits=Inits,grps=grps,Control=Control,adapt=adapt,Prior.pars=Prior.pars,Psi=Psi,post.loss=post.loss)

  ###3) plot and summarize results; note that chain would need to be run a lot longer to summarize the posterior very well!!!
  #plot(Out$MCMC)
	#summary_N(Out)
  #post_loss(Out)
	

