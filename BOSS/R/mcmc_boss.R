#' Function for MCMC analysis of BOSS surveys 
#' 
#' @param Par 	A list comprised of the following parameters:
#' 		"hab.pois": a vector giving the current iteration's linear model parameters for Poisson abundance intensity; each row gives parameters for a particular species
#'   	"hab.bern": a vector giving the current iteration's linear model parameters for Bernoulli part of ZIP model for abundance (if Meta$ZIP=TRUE)
#' 		"Nu": a vector giving the log of the abundance intensity for each strata;
#'    "Eta.pois": If Meta$spat.ind==FALSE, spatial random effects for Poisson abundance model; one for each cell and for each species
#'    "Eta.bern": If Meta$spat.ind==FALSE & Meta$ZIP=TRUE, spatial random effects for Bernoulli abundance model; one for each cell and for each species
#'	  "tau.eta.pois": If Meta$spat.ind==FALSE, precision for spatial ICAR model(s) for the Poisson component
#'    "tau.eta.bern": If Meta$spat.ind==FALSE & Meta$ZIP=TRUE, precision for spatial ICAR model(s) for the Bernoulli component
#'	  "tau.nu": Precision for Nu (overdispersion relative to the Poisson distribution)
#' 		"G": a vector giving the number of groups of animals in each strata; 
#' 		"N": a vector giving the number of animals in each strata
#' 		"Psi": a matrix holding observation assignment probabilities (true species on rows and observation type on columns)
#' 		"Cov.par": an (n.species X n X n.ind.cov)  array holding parameters of individual covariate distributions.
#' @param Dat   A matrix with the first row giving the transect, the second a binary indicator for whether the observation 
#'      had an accompanying photograph, the third gives observation type (if photographed), the fourth holds latent species values,
#'      the fifth to fifth+n.obs.cov give observer/survey condition covariates.  The final columns aftter this hold individual level
#'      covariates (starting with group size)
#' @param Psi An array holding posterior predictions of classification probabilities (from another analysis); dimension (n.species,n.obs.types,# posterior samples)
#' @param cur.iter   Number of iterations to run
#' @param adapt	If adapt==TRUE, run MCMC in adapt mode, optimizing MCMC proposal distributions prior to primary MCMC
#' @param Control	A list object including the following objects:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'	"adapt": if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal 
#' 				ranges prior to final MCMC run; 
#'	"MH.nu": MH tuning parameter for Nu parameters (Langevin-Hastings multivariate update);
#' @param DM.pois.hab	A design matrix for the Poisson model for abundance intensity (log scale)
#' @param DM.bern.hab If Meta$ZIP=TRUE, a design matrix for the Bernoulli zero model (probit scale)
#' @param Q			An inverse precision matrix for the spatial ICAR process
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following objects
#'	"a.eta": alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'  "b.eta": beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'	"a.nu": alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))
#'	"b.nu": beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu)) 
#'	"beta.tau": Prior precision for regression coefficients 
#' @param Meta	A list object giving a number of other features of the dataset, including:
#' 	"n.transects"	Number of transects
#'  "n.species"     Number of species
#' 	"S"				Number of strata cells
#'  "spat.ind"		Indicator for spatial dependence
#'  "Area.hab"		Vector giving relative area covered by each strata
#'  "Area.trans"	Vector giving fraction of area of relevant strata covered by each transect
#'  "Prop.photo"  Vector giving proportion of area in each transect that is photographed
#' 	"Adj"			Adjacency matrix giving connectivity of spatial grid cells
#'  "Mapping" 		Vector mapping each transect into a parent strata
#'  "Covered.area"	Vector giving the fraction of each strata covered by transects
#'  "factor.ind"	Indicator vector specifying whether data columns are factors (1) or continuous (0)
#'  "Levels"		a list object, whose elements are comprised of detection model names; each element gives total # of levels in the combined dataset
#'  "G.transect"	vector holding current number of groups of animals present in area covered by each transect		
#'  "N.transect"    vector holding current number of animals present in covered area by each transect
#'  "grps"			indicator for whether observations are for groups rather than individuals
#'  "n.ind.cov" 	Number of individual covariates (distance is not included in this total, but group size is)
#'  "Cov.prior.pdf" character vector giving the probability density function associated with each individual covariate (type ? hierarchical_DS for more info)
#'  "Cov.prior.parms"	An (n.species X n X n.ind.cov) array providing "pseudo-prior" parameters for individual covarate distributions (only the first row used if a signle parameter distribution)
#'  "Cov.prior.fixed" indicator vector for whether parameters of each covariate distribution should be fixed within estimation routine
#'  "Cov.prior.n" 	(#species X #covariates) Matrix giving number of parameters in each covariate pdf 
#'  "ZIP"  If TRUE, fit a ZIP model to abundance
#'  "fix.tau.nu"	Indicator for whether tau.nu should be fixed (1) or estimated(0)
#'  "srr"			Indicator for whether a spatially restricted regression model should be employed (1) or not (0)
#'  "srr.tol"		Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation
#'  "N.hab.pois.par"	    A vector specifying the number of parameters needed for each species' Poisson abundance model
#'  "N.hab.bern.par"    If fitting a ZIP model, this vector specifying the number of parameters needed for each species' Bernoulli zero model
#'  "post.loss"  If TRUE, observed and predicted detections are compiled for posterior predictive loss 
#' @return returns a list with the following objects: 
#' 	"MCMC": An 'mcmc' object (see 'coda' R package) containing posterior samples;
#'  "Accept": A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
#'  "Control": A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used) 
#'  "Obs.N":  Records latent abundance in each transect; dimension is (n.species X # samples X # transects)
#'  "Pred.N": Posterior predictive distribution for abundance in each transect; obtained by sampling a Poisson distribution given current parameter values (with possible zero inflation)
#'  "Post": Holds posterior samples for strata specific group sizes ("Post$G") and abundance ("Post$N")
#'  "Obs.det":  if Meta$post.loss=TRUE, a matrix holding observed detection types for posterior predictive loss calculations dim = c(n.transects,n.obs.types) 
#'  "Pred.det": if Meta$post.loss=TRUE, an array holding predicted detection types for posterior predictive loss calculations dim = c(n.mcmc.iter,n.transects,n.obs.types)
#' @export
#' @import Matrix
#' @keywords areal, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn

mcmc_boss<-function(Par,Dat,Psi,cur.iter,adapt,Control,DM.hab.pois,DM.hab.bern=NULL,Q,Prior.pars,Meta){	
	#require(mvtnorm)
	#require(Matrix)
	#require(truncnorm)
 	Lam.index=c(1:Meta$S)
	grp.pl=NULL
	if(Meta$grps==TRUE)grp.pl=which(colnames(Dat)=="Group")
  
	
	#initialize G.obs (number of groups observed per transect)
  n.obs=nrow(Dat)
	g.tot.obs=colSums(Meta$G.transect) #total number of observations of animals seen at least once
	n.obs.cov=ncol(Dat)-Meta$n.ind.cov-4
  i.na=(sum(is.na(Dat[,"Species"]))>0)
  if(i.na==0)Obs.NArm=Dat[,"Species"]
  if(i.na==1)Obs.NArm=Dat[,"Species"][-which(is.na(Dat[,"Species"]))]
 	n.obs.types=length(unique(Obs.NArm))
 	
	Z=matrix(1,Meta$n.species,Meta$S)  #need to define for non-ZIP models
  #in case of ZIP model, initialize Z, Z.tilde
  if(Meta$ZIP){
    #start all zeros as arising from Bernoulli component
    Z[Par$G==0]=0    
    Z.tilde=matrix(0,Meta$n.species,Meta$S)
    G.gt0=(Par$G>0)
    G.eq0=(Par$G==0)
    for(isp in 1:Meta$n.species){
      n.gt0=sum(G.gt0[isp,])
      n.eq0=sum(G.eq0[isp,])
      ExpZ=DM.hab.bern[[isp]]%*%Par$hab.bern[isp,]
      if(n.gt0>0)Z.tilde[isp,which(G.gt0[isp,]==1)]=rtruncnorm(n.gt0,a=0,b=Inf,ExpZ[G.gt0[isp,]],1)
      if(n.eq0>0)Z.tilde[isp,which(G.eq0[isp,]==1)]=rtruncnorm(n.eq0,a=-Inf,b=0,ExpZ[G.eq0[isp,]],1)
    }
  }

	#initialize lambda
	Lambda=matrix(0,Meta$n.species,Meta$S)
	Lambda.trans=matrix(0,Meta$n.species,Meta$n.transects)
	for(isp in 1:Meta$n.species){
		Lambda[isp,]=exp(Par$Nu[isp,])*Meta$Area.hab
		Lambda.trans[isp,]=Lambda[isp,Meta$Mapping]*Meta$Area.trans
	}
	grp.lam=rep(0,Meta$n.species)
	
	
	#initialize statistics/matrices needed for MCMC updates
	XpXinv.pois=vector('list',Meta$n.species)
	XpXinvXp.pois=XpXinv.pois
	for(isp in 1:Meta$n.species){
		XpXinv.pois[[isp]]=solve(crossprod(DM.hab.pois[[isp]]))
		XpXinvXp.pois[[isp]]=XpXinv.pois[[isp]]%*%t(DM.hab.pois[[isp]])
	}
  
  if(Meta$ZIP){
    XpXinv.bern=vector('list',Meta$n.species)
    XpXinvXp.bern=XpXinv.bern
    for(isp in 1:Meta$n.species){
      XpXinv.bern[[isp]]=solve(crossprod(DM.hab.bern[[isp]]))
      XpXinvXp.bern[[isp]]=XpXinv.bern[[isp]]%*%t(DM.hab.bern[[isp]])
    }
  }

	if(Meta$srr){
		L.t.pois=XpXinv.pois
		L.pois=L.t.pois
		Qt.pois=L.t.pois
		cross.L.pois=L.t.pois
		Theta.pois=L.t.pois
		N.theta.pois=rep(0,Meta$n.species)
		for(isp in 1:Meta$n.species){
			P.c=diag(Meta$S)-DM.hab.pois[[isp]]%*%solve(crossprod(DM.hab.pois[[isp]]),t(DM.hab.pois[[isp]]))
			Omega=(P.c%*%Meta$Adj%*%P.c)*(Meta$S/sum(Meta$Adj))
			Eigen=eigen(Omega)
			if(max(Eigen$values)<Meta$srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
			Ind=which(Eigen$values>Meta$srr.tol)
			L.t.pois[[isp]]=Eigen$vectors[,Ind]
			cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
			L.pois[[isp]]=t(L.t.pois[[isp]])
			Qt.pois[[isp]]=L.pois[[isp]]%*%Q%*%L.t.pois[[isp]]
			cross.L.pois[[isp]]=L.pois[[isp]]%*%L.t.pois[[isp]]
			N.theta.pois[isp]=nrow(Qt.pois[[isp]])
			Theta.pois[[isp]]=rnorm(N.theta.pois[isp],0,sqrt(1/Par$tau.eta.pois[isp]))
		}
    if(Meta$ZIP){
      L.t.bern=XpXinv.bern
      L.bern=L.t.bern
      Qt.bern=L.t.bern
      cross.L.bern=L.t.bern
      Theta.bern=L.t.bern
      N.theta.bern=rep(0,Meta$n.species)
      for(isp in 1:Meta$n.species){
        P.c=diag(Meta$S)-DM.hab.bern[[isp]]%*%solve(crossprod(DM.hab.bern[[isp]]),t(DM.hab.bern[[isp]]))
        Omega=(P.c%*%Meta$Adj%*%P.c)*(Meta$S/sum(Meta$Adj))
        Eigen=eigen(Omega)
        if(max(Eigen$values)<Meta$srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
        Ind=which(Eigen$values>Meta$srr.tol)
        L.t.bern[[isp]]=Eigen$vectors[,Ind]
        cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
        L.bern[[isp]]=t(L.t.bern[[isp]])
        Qt.bern[[isp]]=L.bern[[isp]]%*%Q%*%L.t.bern[[isp]]
        cross.L.bern[[isp]]=L.bern[[isp]]%*%L.t.bern[[isp]]
        N.theta.bern[isp]=nrow(Qt.bern[[isp]])
        Theta.bern[[isp]]=rnorm(N.theta.bern[isp],0,sqrt(1/Par$tau.eta.bern[isp]))
      }
    }
	}
  
  
	Sampled=unique(Meta$Mapping)
	n.unique=length(Sampled)
	Sampled.area.by.strata=rep(0,n.unique)
	for(i in 1:Meta$n.transects)Sampled.area.by.strata[Sampled==Meta$Mapping[i]]=Sampled.area.by.strata[which(Sampled==Meta$Mapping[i])]+Meta$Area.trans[i]
	
	#initialize MCMC, Acceptance rate matrices
	mcmc.length=(Control$iter-Control$burnin)/Control$thin
  MCMC=list(Psi=array(0,dim=c(dim(Psi)[1:2],mcmc.length)),N.tot=matrix(0,Meta$n.species,mcmc.length),N=array(0,dim=c(Meta$n.species,mcmc.length,Meta$S)),G=array(0,dim=c(Meta$n.species,mcmc.length,Meta$S)),Hab.pois=array(0,dim=c(Meta$n.species,mcmc.length,ncol(Par$hab.pois))),tau.eta.pois=matrix(0,Meta$n.species,mcmc.length),tau.nu=matrix(0,Meta$n.species,mcmc.length),Cov.par=array(0,dim=c(Meta$n.species,mcmc.length,length(Par$Cov.par[1,,]))))
  if(Meta$ZIP){
    MCMC$Hab.bern=array(0,dim=c(Meta$n.species,mcmc.length,ncol(Par$hab.bern)))
    MCMC$tau.eta.bern=matrix(0,Meta$n.species,mcmc.length)
  }
	Accept=list(N=matrix(0,Meta$n.species,Meta$n.transects),Nu=matrix(0,Meta$n.species,n.unique),Psi=rep(0,Meta$n.transects))
	Pred.N=array(0,dim=c(Meta$n.species,mcmc.length,Meta$n.transects))
	Obs.N=Pred.N

	Obs.det=NA
	Pred.det=NA
	Which.photo=which(Dat[,"Photo"]==1)
	Which.no.photo=which(Dat[,"Photo"]==0)
	n.photo=length(Which.photo)
	n.misID.updates=round(n.photo*0.5)
	n.no.photo=length(Which.no.photo) 
  Index.misID=c(1:n.misID.updates)
  Index.photo=c(1:n.photo)
  
  #establish framework for individual covariate contributions to species updates (all possible combos of individual covariates)
  if(Meta$n.ind.cov>0){
    if(Meta$n.ind.cov==1)Cov.pointer=matrix(unique(Dat[Which.photo,4+n.obs.cov+1]),length(unique(Dat[Which.photo,4+n.obs.cov+1])),1)
    else{
      Unique.tags=apply(Dat[Which.photo,(4+n.obs.cov+1):(4+n.obs.cov+Meta$n.ind.cov)],1,'paste',collapse='')
      Unique.tags=which(!duplicated(Unique.tags))
      Cov.pointer=Dat[Which.photo,(4+n.obs.cov+1):(4+n.obs.cov+Meta$n.ind.cov)][Unique.tags,]
    }
    colnames(Cov.pointer)=colnames(Dat)[(4+n.obs.cov+1):ncol(Dat)]   
    Cov.logL=matrix(0,nrow(Cov.pointer),Meta$n.species)  #this will hold log likelihoods for each species
    #fill log likelihoods
    for(isp in 1:Meta$n.species){
      for(icov in 1:Meta$n.ind.cov){
        Cov.logL[,isp]=Cov.logL[,isp]+sapply(Cov.pointer[,icov],'switch_pdf',pdf=Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,,icov],RE=0)       
      }
    }
    Cov.pl=apply(matrix(Dat[Which.photo,(4+n.obs.cov+1):(4+n.obs.cov+Meta$n.ind.cov)],ncol=Meta$n.ind.cov),1,'get_place',Pointer=Cov.pointer)      
  }
     

	if(Meta$post.loss){ #calculate observed counts of different detection types, initialize prediction arrays
    Obs.det=matrix(0,Meta$n.transects,n.obs.types+1) #col=n.obs.types+1 is for hotspots w/o accompanying photographs
    Pred.det=array(0,dim=c(mcmc.length,Meta$n.transects,n.obs.types+1))
    for(iobs in 1:n.photo)Obs.det[Dat[Which.photo[iobs],"Transect"],Dat[Which.photo[iobs],"Obs"]]=Obs.det[Dat[Which.photo[iobs],"Transect"],Dat[Which.photo[iobs],"Obs"]]+Dat[Which.photo[iobs],"Group"]
    Obs.det[,n.obs.types+1]=tabulate(Dat[Which.no.photo,"Transect"],nbins=n.transects)     
  }
	
	#initialize random effect matrices for individual covariates if required
	if(sum(1-Meta$Cov.prior.fixed)>0)RE.cov=matrix(0,n.obs,Meta$n.ind.cov)
  
  #make copy of Dat, with transect as a factor
  Dat2=Dat #establish levels so xtabs works when there are zeros 
  Dat2[,"Transect"]=factor(Dat[,"Transect"],levels=c(1:Meta$n.transects))  
  Dat2[,"Species"]=factor(Dat2[,"Species"],levels=c(1:Meta$n.species))  
	
	PROFILE=FALSE  #outputs time associated with updating different groups of parameters
	DEBUG=FALSE
	if(DEBUG){
		Par$misID[[1]]=2
		Par$misID[[2]]=1
		Par$misID[[3]]=3
	}
	st <- Sys.time()
	##################################################
	############   Begin MCMC Algorithm ##############
	##################################################
	for(iiter in 1:cur.iter){
		#cat(paste('\n ', iiter))
		for(isp in 1:Meta$n.species){		
		  #update Z,Z.tilde if ZIP model specified
		  if(Meta$ZIP){
		    Z[isp,which(Par$G[isp,]>0)]=1
		    Which.0=which(Par$G[isp,]==0)
		    Cur.p=pnorm(DM.hab.bern[[isp]]%*%Par$hab.bern[isp,]+Par$Eta.bern[isp,]) #bernoulli success prob
		    #Cur.pois0=Cur.p*exp(-Meta$Area.hab*exp(Mu))
		    Cur.pois0=Cur.p*exp(-Meta$Area.hab*exp(Par$Nu[isp,]))
		    Cur.p=Cur.pois0/((1-Cur.p)+Cur.pois0)
		    if(length(Which.0)>0)Z[isp,which(Par$G[isp,]==0)]=rbern(length(Which.0),Cur.p[Which.0])
		    Z.tilde[isp,]=rtruncnorm(S,a=ifelse(Z[isp,]==1,0,-Inf),b=ifelse(Z[isp,]==1,Inf,0),DM.hab.bern[[isp]]%*%Par$hab.bern[isp,]+Par$Eta.bern[isp,])       
		  }
		  
		  #update nu parameters (log lambda)
	 	  #1) for sampled cells for which z-tilde>0 (if ZIP)
			Hab.pois=Par$hab.pois[isp,1:Meta$N.hab.pois.par[isp]]
	 	  Eta.pois=Par$Eta.pois[isp,]
		  Mu=DM.hab.pois[[isp]]%*%Hab.pois+Eta.pois
		  G.sampled=rep(0,n.unique) #total number of groups currently in each sampled strata
		  for(i in 1:Meta$n.transects)G.sampled[Sampled==Meta$Mapping[i]]=G.sampled[Sampled==Meta$Mapping[i]]+Meta$G.transect[isp,i]
      for(i in 1:n.unique){
        if(Z[isp,Sampled[i]]==1){
		      prop=Par$Nu[isp,Sampled[i]]+runif(1,-Control$MH.nu[isp,i],Control$MH.nu[isp,i])				
		      old.post=dnorm(Par$Nu[isp,Sampled[i]],Mu[Sampled[i]],1/sqrt(Par$tau.nu[isp]),log=TRUE)+dpois(G.sampled[i],Sampled.area.by.strata[i]*exp(Par$Nu[isp,Sampled[i]]),log=TRUE)
		      new.post=dnorm(prop,Mu[Sampled[i]],1/sqrt(Par$tau.nu[isp]),log=TRUE)+dpois(G.sampled[i],Sampled.area.by.strata[i]*exp(prop),log=TRUE)
		      if(runif(1)<exp(new.post-old.post)){
		        Par$Nu[isp,Sampled[i]]=prop
		        Accept$Nu[isp,i]=Accept$Nu[isp,i]+1
		      }
        }
		  }
		
			#2) simulate nu for areas not sampled
			Par$Nu[isp,-Sampled]=rnorm(Meta$S-n.unique,Mu[-Sampled],1/sqrt(Par$tau.nu[isp]))
		  #3) if ZIP, sample nu for areas where Z=0
      if(Meta$ZIP){
        which.Z.eq0=which(Z[isp,]==0)   
        sampled.Z0=which.Z.eq0[which.Z.eq0%in%Sampled]
        if(length(sampled.Z0)>0)Par$Nu[isp,sampled.Z0]=rnorm(length(sampled.Z0),Mu[sampled.Z0],1/sqrt(Par$tau.nu[isp]))
      }
		  if(PROFILE==TRUE){
		    cat(paste("Nu: ", (Sys.time()-st),'\n'))
		    st=Sys.time()
		  } 
			########## update abundance parameters at the strata scale   ################
      #if(iiter==1035){
      #  cat("\n Par$G \n")
      #  print(Par$G)
      #}

      #update spatial random effects
			if(Meta$spat.ind==FALSE){
				if(Meta$srr==FALSE){
					#update eta parameters (spatial random effects)
					V.eta.inv <- Par$tau.nu[isp]*diag(Meta$S) + Par$tau.eta.pois[isp]*Q
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*(Par$Nu[isp,]-DM.hab.pois[[isp]]%*%Hab.pois))		
					Par$Eta.pois[isp,]<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
					Par$Eta.pois[isp,]=Par$Eta.pois[isp,]-mean(Par$Eta.pois[isp,])  #centering
					#update tau_eta  (precision of spatial process)
					Par$tau.eta.pois[isp] <- rgamma(1, (Meta$S-1)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta.pois[isp,], Q %*% Par$Eta.pois[isp,])*0.5) + Prior.pars$b.eta)
          if(Meta$ZIP){
            Hab.bern=Par$hab.bern[isp,]
            V.eta.inv <- diag(Meta$S) + Par$tau.eta.bern[isp]*Q
            M.eta <- solve(V.eta.inv, Z.tilde[isp,]-DM.hab.bern[[isp]]%*%Hab.bern)		
            Par$Eta.bern[isp,]<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
            Par$Eta.bern[isp,]=Par$Eta.bern[isp,]-mean(Par$Eta.bern[isp,])  #centering
            #update tau_eta  (precision of spatial process)
            Par$tau.eta.bern[isp] <- rgamma(1, (Meta$S-1)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta.bern[isp,], Q %*% Par$Eta.bern[isp,])*0.5) + Prior.pars$b.eta)            
					}
				}
				else{
					#Update Theta
					Dat.minus.Exp=Par$Nu[isp,]-DM.hab.pois[[isp]]%*%Hab.pois
					V.eta.inv <- cross.L.pois[[isp]]*Par$tau.nu[isp] + Par$tau.eta.pois[isp]*Qt.pois[[isp]]
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*L.pois[[isp]]%*%Dat.minus.Exp)
					Theta.pois <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(N.theta.pois[isp],0,1))
					Par$Eta.pois[isp,]=as.numeric(L.t.pois[[isp]]%*%Theta.pois)					
					#update tau.eta
					Par$tau.eta.pois[isp] <- rgamma(1, N.theta.pois[isp]*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Theta.pois, Qt.pois[[isp]] %*% Theta.pois)*0.5) + Prior.pars$b.eta)
          if(Meta$ZIP){
            Hab.bern=Par$hab.bern[isp,]
            Dat.minus.Exp=Z.tilde[isp,]-DM.hab.bern[[isp]]%*%Hab.bern
            V.eta.inv <- cross.L.bern[[isp]] + Par$tau.eta.bern[isp]*Qt.bern[[isp]]
            M.eta <- solve(V.eta.inv, L.bern[[isp]]%*%Dat.minus.Exp)
            Theta.bern <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(N.theta.bern[isp],0,1))
            Par$Eta.bern[isp,]=as.numeric(L.t.bern[[isp]]%*%Theta.bern)	
            #update tau.eta
            Par$tau.eta.bern[isp] <- rgamma(1, N.theta.bern[isp]*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Theta.bern, Qt.bern[[isp]] %*% Theta.bern)*0.5) + Prior.pars$b.eta)
          }
				}
			}
			#update tau_nu	 (precision for Poisson overdispersion)
			Mu=DM.hab.pois[[isp]]%*%Hab.pois+Par$Eta.pois[isp,]
			if(Meta$fix.tau.nu==FALSE){
        Cur.ind=c(Sampled,which(Z[isp,]==1))
        Cur.ind=Cur.ind[duplicated(Cur.ind)]
				Diff=Par$Nu[isp,Cur.ind]-Mu[Cur.ind]
				Par$tau.nu[isp] <- rgamma(1,length(Cur.ind)/2 + Prior.pars$a.nu, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.nu)
			}
			if(PROFILE==TRUE){
				cat(paste("Tau nu: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}

			#translate to lambda scale 
      Lambda[isp,]=exp(Par$Nu[isp,])*Meta$Area.hab  
		  if(Meta$ZIP)Lambda[isp,]=Lambda[isp,]*Z[isp,] 
			Lambda.trans[isp,]=Lambda[isp,Meta$Mapping]*Meta$Area.trans

			if(DEBUG==FALSE){
				#update Betas for habitat relationships
			  Par$hab.pois[isp,1:Meta$N.hab.pois.par[isp]]=rmvnorm(1,XpXinvXp.pois[[isp]]%*%(Par$Nu[isp,]-Par$Eta.pois[isp,]),XpXinv.pois[[isp]]/(Par$tau.nu[isp]+Prior.pars$beta.tau))
        if(Meta$ZIP)Par$hab.bern[isp,1:Meta$N.hab.bern.par[isp]]=rmvnorm(1,XpXinvXp.bern[[isp]]%*%(Z.tilde[isp,]-Par$Eta.bern[isp,]),XpXinv.bern[[isp]]/(1+Prior.pars$beta.tau))
			}
		
			if(PROFILE==TRUE){
				cat(paste("Habitat vars, etc.: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
		}
    
    #if(iiter==1568){
      #cat(paste("iiter ",iiter,'\n'))
      #cat("G.transect \n")
      #cat(G.transect)
    #}

    if(length(Which.no.photo)>0){
      ########## update species for observations without photos
      Dat[Which.no.photo,"Species"]=sapply(Dat[Which.no.photo,"Transect"],sample_nophoto_sp,Lam=Lambda.trans,n.sp=Meta$n.species)
      Dat2[,"Species"]=Dat[,"Species"]
      Meta$G.transect=xtabs(~Species+Transect,data=Dat2)  #recalculate number of groups per species/transect combo
      ########## update ind. covariates for observations without photos  #############
      if(Meta$n.ind.cov>0){
        for(icov in 1:Meta$n.ind.cov){
          if(Meta$Cov.prior.pdf[icov]=='poisson_ln' | Meta$Cov.prior.pdf[icov]=='pois1_ln')cur.RE=RE.cov[,icov]
          else cur.RE=0
          for(isp in 1:Meta$n.species){
            I.cur=(Dat[,"Photo"]==0)*(Dat[,"Species"]==isp)
            if(length(cur.RE>0))cur.RE=cur.RE[I.cur==1]
            Dat[I.cur==1,4+n.obs.cov+icov]=switch_sample(n=sum(I.cur),pdf=Meta$Cov.prior.pdf[icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
          }
        }
      }	
    }
    
	
#     if(iiter==1568){
#      # cat("\n Lambda \n")
#       print(Lambda)
#       cat("\n G.transect \n")
#       print(G.transect)
#       cat("\n Par$G \n")
#       print(Par$G)
#       cat("\n Lambda.trans \n")
#       print(Lambda.trans)
		#     }
		##### update true species for observed animals ######
		Which.sampled=sample(n.photo,n.misID.updates) #replace needs to be false here or there's issues with using Old.sp
		Old.sp=Dat[Which.photo[Which.sampled],"Species"]
		Prop.sp=sample(c(1:Meta$n.species),n.misID.updates,replace=TRUE)
		Cur.obs=Dat[Which.photo[Which.sampled],"Obs"]
		#covariate contributions
		MH=sapply(Index.misID,"get_mat_entries",Mat=Cov.logL,Row=Cov.pl[Which.sampled],Col=Prop.sp)-sapply(Index.misID,"get_mat_entries",Mat=Cov.logL,Row=Cov.pl[Which.sampled],Col=Old.sp)
		#confusion matrix contributions
		if(Meta$misID)MH=MH+log(sapply(Index.misID,"get_mat_entries",Mat=Par$Psi,Row=Prop.sp,Col=Dat[Which.photo[Which.sampled],"Obs"]))-log(sapply(Index.misID,"get_mat_entries",Mat=Par$Psi,Row=Dat[Which.photo[Which.sampled],"Species"],Col=Dat[Which.photo[Which.sampled],"Obs"]))
		#Poisson abundance model contributions (new state)
		MH=MH+log(sapply(Index.misID,"get_mat_entries",Mat=Lambda.trans,Row=Prop.sp,Col=Dat[Which.photo[Which.sampled],"Transect"]))
		Rsamp=runif(n.misID.updates)
		#       if(iiter==1568){
		#       cat("\n MH \n")
		#       print(MH)
		#       cat("\n Dat \n")
		#       print(Dat[Which.photo[Which.sampled],])
		#       cat("\n Prop.sp \n")
		#       print(Prop.sp)
		#       cat("\n Mapping \n")
		#       print(Mapping)
		#       }
		for(isamp in 1:n.misID.updates){
		  mh=MH[isamp]-log(Lambda.trans[Dat[Which.photo[Which.sampled[isamp]],"Species"],Dat[Which.photo[Which.sampled[isamp]],"Transect"]])
		  if(Rsamp[isamp]<exp(mh)){
		    Meta$G.transect[Old.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]=Meta$G.transect[Old.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]-1
		    Meta$G.transect[Prop.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]=Meta$G.transect[Prop.sp[isamp],Dat[Which.photo[Which.sampled[isamp]],"Transect"]]+1
		    Dat[Which.photo[Which.sampled[isamp]],"Species"]=Prop.sp[isamp]
		  }
		}
		
		if(PROFILE==TRUE){
		  cat(paste("True species: ", (Sys.time()-st),'\n'))
		  st=Sys.time()
		}
    #if(iiter==1568)cat("\n made it 1 \n")
    
		########## recalculate abundance totals
    Dat2[,"Group"]=Dat[,"Group"]
    Dat2[,"Species"]=factor(Dat[,"Species"],levels=c(1:Meta$n.species))
    Meta$N.transect=xtabs(Group~Species+Transect,data=Dat2)

		#if(iiter==1568){
      #cat("\n Meta$N.transect \n")
      #print(Meta$N.transect)
      #cat("\n made it 2 \n")
		#}
    for(isp in 1:Meta$n.species){
      #if(iiter==1568)cat(paste("\n isp ",isp,'\n'))
      ########## update group abundance at strata level using most recent updates to covered area abundance
      Par$G[isp,]=rpois(Meta$S,Lambda[isp,]*(1-Meta$Covered.area))
      grp.lam[isp]=ifelse(Meta$Cov.prior.pdf[isp,1] %in% c("pois1_ln","poisson_ln"),exp(Par$Cov.par[isp,1,1]+(Par$Cov.par[isp,2,1])^2/2),Par$Cov.par[isp,1,1])
      Par$N[isp,]=Par$G[isp,]+rpois(Meta$S,grp.lam[isp]*Par$G[isp,]) #add the Par$G since number in group > 0 
      for(ipl in 1:length(Meta$Mapping)){
        #if(iiter==1568)cat(paste("\n ipl ",ipl,'\n'))
        Par$G[isp,Meta$Mapping[ipl]]=Par$G[isp,Meta$Mapping[ipl]]+Meta$G.transect[isp,ipl]
        Par$N[isp,Meta$Mapping[ipl]]=Par$N[isp,Meta$Mapping[ipl]]+Meta$N.transect[isp,ipl]
      }
    }		
    
		#if(iiter==1568)cat("\n made it 3 \n")
		
    ##### update misID parameters in independence chain
    if(Meta$misID){
      Psi.prop=Psi[,,sample(dim(Psi)[3],1)]
		  mh=sum(log(sapply(Index.photo,"get_mat_entries",Mat=Psi.prop,Row=Dat[Which.photo,"Species"],Col=Dat[Which.photo,"Obs"]))-log(sapply(Index.photo,"get_mat_entries",Mat=Par$Psi,Row=Dat[Which.photo,"Species"],Col=Dat[Which.photo,"Obs"])))
      if(runif(1)<exp(mh)){
        Par$Psi=Psi.prop
        Accept$Psi=Accept$Psi+1
      }
    }
		
		#####update parameters of individual covariate distributions (if fixed=0)
    Cov.logL=Cov.logL*0
		for(isp in 1:Meta$n.species){
      if(sum(Meta$G.transect[isp,])>0){
        for(icov in 1:Meta$n.ind.cov){
          Cur.cov=Dat[Dat[,"Species"]==isp & Dat[,"Photo"]==1,4+n.obs.cov+icov] 
          if(Meta$Cov.prior.fixed[isp,icov]==0){
            if(Meta$Cov.prior.pdf[isp,icov]=="normal")cat("\n Warning: hyper-priors not yet implemented for normal dist. \n")
            if(Meta$Cov.prior.pdf[isp,icov]=="poisson"){
              Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="pois1"){
              Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)-length(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="poisson_ln" | Meta$Cov.prior.pdf[isp,icov]=="pois1_ln"){
              cat('RE on covariates not currently implemented for BOSS analysis')
            }
            if(Meta$Cov.prior.pdf[isp,icov]=="multinom"){
              Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov]=rdirichlet(1,Meta$Cov.prior.parms[isp,1:Meta$Cov.prior.n[isp,icov],icov]+tabulate(factor(Cur.cov)))
            }
          }
        }
      }
      #update log likelihood contributions for combinations of different covariate values
      for(icov in 1:Meta$n.ind.cov){
        Cov.logL[,isp]=Cov.logL[,isp]+sapply(Cov.pointer[,icov],'switch_pdf',pdf=Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,,icov],RE=0)       
      }    
		}
	
		if(PROFILE==TRUE){
			cat(paste("Ind cov pars: ", (Sys.time()-st),'\n'))
			st=Sys.time()
		}
		#if(iiter==1568)cat("\n made it 4 \n")
		
		
		if(adapt==TRUE){
			if(iiter%%100==0){
				for(ipar in 1:Meta$n.species){
					for(i in 1:n.unique)
					if(Accept$Nu[ipar,i]<30)Control$MH.nu[ipar,i]=Control$MH.nu[ipar,i]*.95
					if(Accept$Nu[ipar,i]>40)Control$MH.nu[ipar,i]=Control$MH.nu[ipar,i]*1.053
				}
				Accept$Nu=Accept$Nu*0
				#Accept$MisID=lapply(Accept$MisID,function(x)x*0)
			}
		}
		
		
		#store results if applicable
		if(iiter>Control$burnin & iiter%%Control$thin==0){
      if(Meta$misID==TRUE)MCMC$Psi[,,(iiter-Control$burnin)/Control$thin]=Par$Psi
      for(isp in 1:Meta$n.species){
				MCMC$G[isp,(iiter-Control$burnin)/Control$thin,]=Par$G[isp,]
				MCMC$N[isp,(iiter-Control$burnin)/Control$thin,]=Par$N[isp,]
				MCMC$N.tot[isp,(iiter-Control$burnin)/Control$thin]=sum(Par$N[isp,])
				MCMC$Hab.pois[isp,(iiter-Control$burnin)/Control$thin,]=Par$hab.pois[isp,]
				MCMC$tau.eta.pois[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.eta.pois[isp]
        if(Meta$ZIP){
				  MCMC$Hab.bern[isp,(iiter-Control$burnin)/Control$thin,]=Par$hab.bern[isp,]
				  MCMC$tau.eta.bern[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.eta.bern[isp]
        }
				MCMC$tau.nu[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.nu[isp]
				MCMC$Cov.par[isp,(iiter-Control$burnin)/Control$thin,]=Par$Cov.par[isp,,]
				Obs.N[isp,(iiter-Control$burnin)/Control$thin,]=Meta$N.transect[isp,]
				Temp.G=Meta$Area.hab[Meta$Mapping]*Meta$Area.trans*exp(rnorm(Meta$n.transects,(DM.hab.pois[[isp]]%*%Par$hab.pois[isp,1:Meta$N.hab.pois.par[isp]]+Par$Eta.pois[isp,])[Meta$Mapping],sqrt(1/Par$tau.nu[isp])))
        if(Meta$ZIP)Temp.G=Temp.G*rbern(Meta$n.transects,pnorm((DM.hab.bern[[isp]]%*%Par$hab.bern[isp,1:Meta$N.hab.bern.par[isp]]+Par$Eta.bern[isp,])[Meta$Mapping]))
        Pred.N[isp,(iiter-Control$burnin)/Control$thin,]=Temp.G+rpois(Meta$n.transects,grp.lam[isp]*Temp.G)	
			}
      
			#if(Meta$post.loss){ #calculate observed counts of different detection types, initialize prediction arrays
			#  Obs.det=matrix(0,Meta$n.transects,n.obs.types+1) #col=n.obs.types+1 is for hotspots w/o accompanying photographs
			#  Pred.det=array(0,dim=c(mcmc.length,Meta$n.transects,n.obs.types+1))
			#  for(iobs in 1:n.photo)Obs.det[Dat[Which.photo[iobs],"Transect"],Dat[Which.photo[iobs],"Obs"]]=Obs.det[Dat[Which.photo[iobs],"Transect"],Dat[Which.photo[iobs],"Obs"]]+Dat[Which.photo[iobs],"Group"]
			#  Obs.det[,n.obs.types+1]=tabulate(Dat[Which.no.photo,"Transect"],nbins=n.transects)     
			#}
			
       #posterior predictions of data given nu & misclassification parameters
			if(Meta$post.loss){
        Cur.G=matrix(rpois(Meta$n.species*Meta$n.transects,exp(Par$Nu[,Sampled])),Meta$n.species,Meta$n.transects)		
        if(Meta$ZIP)Cur.G=Z[,Sampled]*Cur.G         
        for(itrans in 1:Meta$n.transects){
          for(isp in 1:Meta$n.species){
            if(Cur.G[isp,itrans]>0){
              Cur.dat=matrix(0,Cur.G[isp,itrans],ncol(Dat))
              Cur.dat[,4]=isp
              #fill observer covariates
              if(n.obs.cov>0){
                cat('Warning: Observer/survey level covariates not yet entirely implemented for BOSS surveys!!! \n') 
              }
              #sample from individual covariate distributions 
              if(Meta$n.ind.cov>0){
                for(icov in 1:Meta$n.ind.cov){
                  #if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=rep(rnorm(Cur.G[isp,itrans],0,1),each=Meta$n.Observers[itrans])
                  #else 
                  cur.RE=0
                  rsamp=switch_sample(n=Cur.G[isp,itrans],pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
                  Cur.dat[,4+n.obs.cov+icov]=rsamp
                }
              }
              #could be sped up by using apply
              for(iind in 1:Cur.G[isp,itrans]){
                if(n.species>1)cur.obs=sample(n.obs.types,1,prob=Par$Psi[isp,])
                else cur.obs=1
                if(runif(1)<Prop.photo[itrans])Pred.det[(iiter-Control$burnin)/Control$thin,itrans,cur.obs]=Pred.det[(iiter-Control$burnin)/Control$thin,itrans,cur.obs]+1    
                else Pred.det[(iiter-Control$burnin)/Control$thin,itrans,n.obs.types+1]=Pred.det[(iiter-Control$burnin)/Control$thin,itrans,n.obs.types+1]+1    
              }           
            }
          }
			  }
			}	
		}
		
		if(iiter==100){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/100
			ttc <- round((cur.iter-100)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}	
		if((iiter%%1000)==1)cat(paste('iteration ', iiter,' of ',cur.iter,' completed \n'))
	}
	cat(paste('\n total elapsed time: ',difftime(Sys.time(),st,units="mins"),' minutes \n'))

	Post=list(N=MCMC$N,G=MCMC$G)
	#convert Out$MCMC into mcmc object for use with coda, S3 methods
	Hab.pois.names=vector("list",Meta$n.species)
	for(isp in 1:Meta$n.species)Hab.pois.names[[isp]]=colnames(DM.hab.pois[[isp]])
  if(Meta$ZIP){
    Hab.bern.names=vector("list",Meta$n.species)
    for(isp in 1:Meta$n.species)Hab.bern.names[[isp]]=colnames(DM.hab.bern[[isp]])
  }
	Cov.names=vector("list",Meta$n.ind.cov)
	Cov.par.n=0
	if(Meta$n.ind.cov>0){
		for(icov in 1:Meta$n.ind.cov){
			Par.name=switch(Meta$Cov.prior.pdf[icov],pois1_ln=c("mean.minus.1","sd"),poisson_ln=c("mean","sd"),multinom=paste("prop.cell.",c(1:(Meta$Cov.prior.n[isp,icov]-1)),sep=''),normal="mean",pois1="mean.minus.1",poisson="mean")
			Cov.names[[icov]]=paste(Meta$stacked.names[Meta$dist.pl+icov],".",Par.name,sep='')
			Cov.par.n=Cov.par.n+length(Cov.names[[icov]])
		}
	}
  if(1==1){
    if(Meta$ZIP){
      N.hab.bern.par=Meta$N.hab.bern.par
      Hab.bern.names=Hab.bern.names
    }
    else{
      N.hab.bern.par=NA
      Hab.bern.names=NA
    }
  }
  if(n.species>1)Post$Psi=MCMC$Psi
	MCMC=convert.BOSS.to.mcmc(MCMC=MCMC,N.hab.pois.par=Meta$N.hab.pois.par,N.hab.bern.par=N.hab.bern.par,Cov.par.n=Cov.par.n,Hab.pois.names=Hab.pois.names,Hab.bern.names=Hab.bern.names,Cov.names=Cov.names,fix.tau.nu=Meta$fix.tau.nu,spat.ind=Meta$spat.ind)
	Out=list(Post=Post,MCMC=MCMC,Accept=Accept,Control=Control,Obs.N=Obs.N,Pred.N=Pred.N,Obs.det=Obs.det,Pred.det=Pred.det)
	Out
}
