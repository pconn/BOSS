#' Primary function for hierarchical, areal analysis of distance sampling data.  This function
#' pre-processes data and calls other functions to perform the analysis, and is the only function
#' the user needs to call themselves. 
#'
#' @param Dat 	A data frame with the following columns:
#' 		(1)numeric transect ID (1 to # of transects)
#' 		(2)photo obtained (0/1) ? 
#' 		(3) Observation type (integer - the max integer being 'unknown' if applicable) [NOTE: modeled as factor, but need to be input as integers to account for unknown species observations]
#' 		(4-x)(Observer covariates); (things like survey conditions or observer skill that affect misID; things that don't change during a transect.  Note these also need to be provided in Obs.cov)
#' 		(x+1-??)(Group size and other individual covariates thought to influence detection; if group size is one of them, it's assumed to be column x+1); could also include observer certainty here
#' 		Note that column names can be used to tag covariates, and that object types (e.g. numeric, factor) will be preserved in analysis (note all covariates must be categorical right now!)
#' @param Adj   Adjacency matrix for habitat cells (diagonal matrix implies spatial independence)
#' @param Area.hab   A vector giving the area of each geographical strata (default is equal area)
#' @param Mapping  A vector giving the habitat cell id number for each transect
#' @param Area.trans	A vector giving the effective area covered by each transect as fraction of total area in the strata it is located
#' @param DayHour A (n.transect X 2) matrix providing row and column entries into the Thin array. Each row corresponds to an entry in Mapping
#' @param Thin An (n.species X n.days X n.hours X n.iter) array providing n.iter posterior samples of the thinning parameters
#' @param Prop.photo A vector giving the proportion of of the area sampled in each transect that is photographed
#' @param n.species An integer giving the true number of species
#' @param n.obs.cov	Number of observer covariates (e.g., visibility, etc.)
#' @param Hab.cov	A data.frame object giving covariates thought to influence abundance intensity at strata level; column names index individual covariates
#' @param Obs.cov  A (# of transects X # of observer covariates) size matrix giving observer covariate values for each transect flown
#' @param Hab.pois.formula	A formula vector giving the specific model for Poisson abundance intensity at the strata level (e.g., ~Vegetation+Latitude) for each species
#' @param Hab.bern.formula  If ZIP=TRUE, a formula vector giving the specific model for the zero component for abundance intensity at the strata level (e.g., ~Vegetation+Latitude) for each species
#' @param Cov.prior.pdf	If individual covariates are provided, this character matrix gives the form of the prior pdfs for each species & covariate
#'		  current possibilities are "poisson", "pois1","poisson_ln","pois1_ln",uniform.disc","multinom","uniform.cont", or "normal".
#'		  "pois1" is 1+x where x~poisson; "poisson_ln" and "pois1_ln" are lognormal poisson models that incorporate overdispersion. 
#' @param Cov.prior.parms	A (s X k X n) array where s is the number of species, n is the number of individual covariates (other than distance), and
#' 		k is the maximum number of parameters considered for a single covariate (NAs can be used to fill this matrix
#'      out for covariate priors that have <k parameters).  If Cov.prior.fixed=1 for a given entry, the prior parameters supplied
#'      in each column apply to the prior pdf itself, and are treated as fixed.  If Cov.prior.fixed=0, the model will attempt
#'  	to estimate the posterior distribution of model parameters, given hyperpriors.  In this case, it is actually the hyperpriors
#'      that are being specified.  For "poisson", and "pois1", it is assumed that lambda~gamma(alpha,beta), so alpha
#' 		and beta must be supplied.  For "poisson_ln", and "pois1_ln", the model is lambda_i=exp(-sigma*Z_i+theta), so it is priors
#' 		for theta and sigma that are specified (in that order).  Theta is assumed to have a normal(mu,s^2) distribution,
#' 		and sigma is assumed to have a uniform(0,a) distribution; thus, priors are specified for these models as (mu,s, and a).
#' 		For the multinomial pdf, prior parameters of the dirichlet distribution must be specified if Cov.prior.fixed=1.
#' @param Cov.prior.fixed  An indicator matrix specifying which (if any) individual covariate distributions should be fixed during estimation
#' @param Cov.prior.n  An (# species X # indiv. covariates) matrix giving the number of parameters in each covariate pdf
#' @param ZIP  If TRUE, estimate ZIP model for abundance that includes a Bernoulli model for zeros and a Poisson + 1 model for positive values (default is FALSE)
#' @param spat.ind	If TRUE, assumes spatial independence (no spatial random effects on abundance intensity) default is FALSE
#' @param fix.tau.nu  If TRUE, fixes tau.nu during estimation (the value to fix it to can be provided in "Inits")
#' @param srr  If TRUE, uses spatially retricted regression, where smoothing occurs on residuals and all spatial effects are orthogonal to the linear predictors (by default, analysis is limited to the highest 50 eigenvalues of the decomposition of the residual projection matrix to reduce computing time)
#' @param srr.tol Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation (default is 0.5)
#' @param Psi An array holding posterior samples of the confusion matrix (dim = #species,#obs,#mcmc.iter)
#' @param Thin An array holding posterior samples of thinning probabilities (dim = # species, # days surveyed, # hours in day, #iterations)
#' @param Thin.pointer A matrix indicating which element of Thin each surveyed cell belongs to
#' @param grps 	If FALSE, detections are assumed to all be of individual animals
#' @param Control	A list object including the following objects:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'	"adapt": if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal ranges prior to final MCMC run; 
#'	"MH.nu": MH tuning parameters for Nu parameters (dimension = # species X # of unique strata sampled)
#' @param Inits	An (optional) list object providing initial values for model parameters, with the following objects:
#'  "hab.pois": Initial values for habitat linear predictor parameters for poisson model;
#'  "hab.bern": If ZIP=TRUE, initial values for habitat linear predictor parameters for bernoulli zero model;
#'	"Nu": Gives log(lambda) for each spatial strata;
#'	"Eta.pois": If spat.ind==FALSE, spatial random effects for Poisson abundance model; one for each cell and for each species
#'  "Eta.bern": If spat.ind==FALSE & ZIP=TRUE, spatial random effects for Bernoulli abundance model; one for each cell and for each species
#'	"tau.eta.pois": If spat.ind==FALSE, precision for spatial ICAR model(s) for the Poisson component
#'  "tau.eta.bern": If spat.ind==FALSE & ZIP=TRUE, precision for spatial ICAR model(s) for the Bernoulli component
#'	"tau.nu": Precision for Nu (overdispersion relative to the Poisson distribution)
#'  One need not specify an initial value for all parameter types (if less are specified, the others are generated randomly)
#' @param adapt	If adapt==TRUE, run an additional Control$adapt number of MCMC iterations to optimize MCMC proposal distributions prior to primary MCMC
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following objects
#'	"a.eta": alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'  "b.eta": beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'	"a.nu": alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))
#'	"b.nu": beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu)) 
#'	"beta.tau": prior precision for regression coefficients (assumed Normal(0,(beta.tau*X'X)^(-1))
#' @param post.loss If TRUE, calculates observed values and posterior predictions for detection data to use with posterior predictive loss functions
#' @return returns a list with the following objecs: 
#' 	MCMC: A list object containing posterior samples;
#'  Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis-Hastings;
#'  Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#' @export
#' @import Matrix
#' @keywords areal model, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn \email{paul.conn@@noaa.gov} 
#' @examples print("example analysis included in the script run_BOSS_sims.R")
hierarchical_boss<-function(Dat,Adj,Area.hab=1,Mapping,Area.trans,DayHour,Thin,Prop.photo=Prop.photo,Hab.cov,Obs.cov,Hab.pois.formula,Hab.bern.formula=NULL,Cov.prior.pdf,Cov.prior.parms,Cov.prior.fixed,Cov.prior.n,n.species=1,n.obs.cov=0,ZIP=FALSE,spat.ind=FALSE,fix.tau.nu=FALSE,srr=TRUE,srr.tol=0.5,Psi,Inits=NULL,grps=FALSE,M,Control,adapt=TRUE,Prior.pars,post.loss=TRUE){
	require(mvtnorm)
	require(Matrix)
	require(truncnorm)
	require(mc2d)
	require(MCMCpack)
	DEBUG=FALSE
	
	Adj=as.matrix(Adj)  #just in case the adjacency matrix = 1 (for 1 transect)
	S=length(Area.hab)
	n.transects=length(Area.trans)
	n.ind.cov=ncol(Dat)-(3+n.obs.cov) #number of individual covariates 
  Obs.NArm=Dat[,3][-which(is.na(Dat[,3]))]
  n.obs.types=length(unique(Obs.NArm))
	
	#By no means exhaustive checking to make sure input values are internally consistent
	if(nrow(Control$MH.nu)!=n.species)cat("\n ERROR: Control$MH.nu does not have # rows = number of species \n")
  #More later...	
  
	#if(length(unique(Dat[,3]))>1)Dat[,3]=as.factor(Dat[,3])  #convert species to factors if not already
	cur.colnames=colnames(Dat)
  cur.colnames[1:3]=c("Transect","Photo","Species")
	if(grps==TRUE)cur.colnames[4+n.obs.cov]="Group"
	if(length(colnames(Dat))!=length(cur.colnames))cat("\n ERROR: mismatch between dimension of Dat and expected columns: check to make sure n.obs.cov, etc. correct")
	colnames(Dat)=cur.colnames

  #convert character entries to factors
  for(i in 1:ncol(Dat))if(is.character(Dat[,i])&length(unique(Dat[,i]))>1)Dat[,i]=as.factor(Dat[,i])
  
  #for factors, determine levels, save labels, and convert to numeric
	factor.ind=sapply(Dat[1,],is.factor)
  which.factors=which(factor.ind==1)
  n.factors=sum(factor.ind)
  
  Factor.labels=vector("list",n.factors)
  if(n.factors>0){
    for(i in 1:n.factors){
      Factor.labels[[i]]=levels(Dat[,which.factors[i]])
    }
  }
  
	Dat.num=Dat
  #if(sum(Dat[,"Obs"]==0)>0)Dat.num[Dat[,"Obs"]==0,"Species"]=0  #if a missing obs, set species=0
  Levels=NULL
  if(n.factors>0){
    for(icol in which.factors){
		  Dat.num[,icol]=as.numeric((Dat[,icol]))
	  }
    Levels=vector("list",n.factors)
	  for(i in 1:n.factors){
	    Levels[[i]]=sort(unique(Dat.num[,which.factors[i]]))
	  }
	  names(Levels)=colnames(Dat[,which.factors])
  }
  #update observer covariate values to reflect new factor values going from 1,2,...
  if(n.obs.cov>0){
	  for(icov in 1:n.obs.cov){
	    if((icov+3)%in%which.factors){
        Obs.cov[,icov]=as.factor(as.numeric(as.factor(Obs.cov[,icov])))
	    }
	  }	
  }
	  
	#add an additional column for "True species" and fill
	True.sp=as.numeric(as.character(Dat[,"Species"]))
  n.missed=sum(is.na(True.sp))
  Which.photo=which(is.na(True.sp)==FALSE)
  if(n.missed>0){
	  True.sp.miss=which(is.na(True.sp)==1)
    if(n.species==1)True.sp[True.sp.miss[1]]=1
	  else True.sp[True.sp.miss[1]]=sample(c(1:(n.species-1)),1)
    if(n.missed>1){
      for(i in 2:length(True.sp.miss)){
        if(n.species==1)True.sp[True.sp.miss[i]]=1
        else True.sp[True.sp.miss[i]]=sample(c(1:n.species),1)
      }
    }
	  #fill group sizes for unphotographed animals 
	  for(imiss in 1:length(True.sp.miss))Dat.num[True.sp.miss[imiss],"Group"]=switch_sample(1,Cov.prior.pdf[True.sp[True.sp.miss[imiss]],1],Cov.prior.parms[True.sp[True.sp.miss[imiss]],,1])
  }
  if(misID){
    for(iind in 1:length(Which.photo)){
      True.sp[Which.photo[iind]]=sample(c(1:n.species),1,prob=Psi[,True.sp[Which.photo[iind]],1])
    }
  }
	  
	if(DEBUG==TRUE)True.sp=Sim$True.species #for debugging
	Dat.num=cbind(Dat.num[1:3],True.sp,Dat.num[4:ncol(Dat.num)])
	
	G.transect=matrix(0,n.species,n.transects)  #number of groups by transect; each row gives results for separate species
	N.transect=G.transect #total abundance by transect
  for(isp in 1:n.species){
    for(itrans in 1:n.transects){
      Cur.ind=which(Dat.num[,"Transect"]==itrans & Dat.num[,"True.sp"]==isp)
      if(length(Cur.ind)>0){
        G.transect[isp,itrans]=length(Cur.ind)
        N.transect[isp,itrans]=sum(Dat.num[Cur.ind,"Group"])
      }
    }
  }
  Dat=Dat.num
  cur.colnames=colnames(Dat)
  cur.colnames[3]="Obs"
  cur.colnames[4]="Species"
  colnames(Dat)=cur.colnames

	N.hab.pois.par=rep(0,n.species)
	DM.hab.pois=vector('list',n.species)
	if(1==1){
		if(is.null(Hab.cov)|Hab.pois.formula[[1]]==~1){
			DM.hab.pois[[1]]=as.matrix(rep(1,S),ncol=1)
			colnames(DM.hab.pois[[1]])="Intercept"
		}
		else DM.hab.pois[[1]]=model.matrix(Hab.pois.formula[[1]],data=Hab.cov)
	}
	N.hab.pois.par[1]=ncol(DM.hab.pois[[1]])
	if(n.species>1){
		for(i in 2:n.species){  #create design matrices for each species. e.g., name for first species will be DM.hab.pois1
			if(is.null(Hab.cov)|Hab.pois.formula[[i]]==~1){
				DM.hab.pois[[i]]=as.matrix(rep(1,S),ncol=1)
				colnames(DM.hab.pois[[i]])="Intercept"
			}
			else DM.hab.pois[[i]]=model.matrix(Hab.pois.formula[[i]],data=Hab.cov)
			N.hab.pois.par[i]=ncol(DM.hab.pois[[i]])
		}
	}
  DM.hab.bern=NULL
  N.hab.bern.par=NULL
	if(ZIP){
	  DM.hab.bern=vector('list',n.species)
	  N.hab.bern.par=rep(0,n.species)
	  if(1==1){
	    if(is.null(Hab.cov)|Hab.bern.formula[[1]]==~1){
	      DM.hab.bern[[1]]=as.matrix(rep(1,S),ncol=1)
	      colnames(DM.hab.bern[[1]])="Intercept"
	    }
	    else DM.hab.bern[[1]]=model.matrix(Hab.bern.formula[[1]],data=Hab.cov)
	  }
	  N.hab.bern.par[1]=ncol(DM.hab.bern[[1]])
	  if(n.species>1){
	    for(i in 2:n.species){  #create design matrices for each species. e.g., name for first species will be DM.hab.bern1
	      if(is.null(Hab.cov)|Hab.bern.formula[[i]]==~1){
	        DM.hab.bern[[i]]=as.matrix(rep(1,S),ncol=1)
	        colnames(DM.hab.bern[[i]])="Intercept"
	      }
	      else DM.hab.bern[[i]]=model.matrix(Hab.bern.formula[[i]],data=Hab.cov)
	      N.hab.bern.par[i]=ncol(DM.hab.bern[[i]])
	    }
	  }
	}	  

	Par=generate_inits_BOSS(DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,N.hab.pois.par=N.hab.pois.par,N.hab.bern.par=N.hab.bern.par,G.transect=G.transect,Area.trans=Area.trans,Area.hab=Area.hab,Mapping=Mapping,spat.ind=spat.ind,grp.mean=Cov.prior.parms[,1,1])	
  if(misID)Par$Psi=Psi[,,sample(dim(Psi)[3],1)]
  if(is.null(Inits)==FALSE){  #replace random inits with user provided inits for all parameters specified
		I.init=names(Inits)
		for(ipar in 1:length(I.init)){
			eval(parse(text=paste("Par$",names(Inits)[ipar],"=Inits$",names(Inits[ipar]))))
		}
	}
	#start Nu out at a compatible level
  for(isp in 1:n.species){
    if(ncol(Par$hab.pois)<ncol(DM.hab.pois[[isp]]))cat('Error: Abundance intensity model has parameter/formula mismatch')
    Par$Nu[isp,]=DM.hab.pois[[isp]]%*%Par$hab.pois[isp,1:N.hab.pois.par[isp]]
  }
  if(length(Par$tau.nu)!=n.species)cat('Error: length of initial value vector for tau.nu should be equal to # of species')
  
	#get initial individual covariate parameter values
	Par$Cov.par=Cov.prior.parms 
	for(i in 1:n.ind.cov){	
		for(j in 1:n.species){
			if(Cov.prior.fixed[j,i]==1)Par$Cov.par[j,,i]=Cov.prior.parms[j,,i]
			else{
				temp=switch_sample_prior(Cov.prior.pdf[j,i],Cov.prior.parms[j,,i])
				Par$Cov.par[j,1:length(temp),i]=temp
			}
		}
	}
	

	n.hab.cov=ifelse(is.null(Hab.cov)==1 | length(Hab.cov)==1,0,ncol(Hab.cov))
		
	i.Covered=c(1:S)%in%Mapping
	Covered.area=rep(0,S)
	for(i in 1:S){
		if(i.Covered[i]==1){
			Covered.area[i]=sum(Area.trans[which(Mapping==i)])
		}
	}

	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q=Matrix(Q)	
  
	Meta=list(n.transects=n.transects,n.species=n.species,S=S,spat.ind=spat.ind,Area.hab=Area.hab,Area.trans=Area.trans,
      DayHour=DayHour,Thin=Thin,Prop.photo=Prop.photo,Adj=Adj,Mapping=Mapping,Covered.area=Covered.area,
			factor.ind=factor.ind,Levels=Levels,misID=misID,
			G.transect=G.transect,N.transect=N.transect,grps=grps,n.ind.cov=n.ind.cov,
			Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,ZIP=ZIP,fix.tau.nu=fix.tau.nu,
			srr=srr,srr.tol=srr.tol,N.hab.pois.par=N.hab.pois.par,N.hab.bern.par=N.hab.bern.par,post.loss=post.loss)

	if(adapt==TRUE){
		cat('\n Beginning adapt phase \n')
		Out=mcmc_boss(Par=Par,Dat=Dat,Psi=Psi,cur.iter=Control$adapt,adapt=1,Control=Control,DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_boss(Par=Par,Dat=Dat,Psi=Psi,cur.iter=Control$iter,adapt=0,Control=Out$Control,DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	else{
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_boss(Par=Par,Dat=Dat,Psi=Psi,cur.iter=Control$iter,adapt=0,Control=Control,DM.hab.pois=DM.hab.pois,DM.hab.bern=DM.hab.bern,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	Out	
}
