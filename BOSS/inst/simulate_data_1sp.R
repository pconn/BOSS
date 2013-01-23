#' function to simulate distance sampling data from simple model with increasing abundance
#' intensity, no assumed spatial structure, and point independence
#' @param S number of spatial strata (a single transect is placed in each strata and assumed to cover the whole strata)
#' @param n.sampled  number of sampled cells
#' @return a distance sampling dataset
#' @export
#' @keywords distance sampling, simulation
#' @author Paul B. Conn
simulate_BOSS<-function(S,prop.photo=1,n.sampled,misID=1,ZIP=TRUE,tau.pois=15,tau.bern=20){
  #note: currently hardwired for 4 species and an anomaly category
	require(mvtnorm)
  library(hierarchicalDS)
	#set.seed(122412)
  Sampled.cells=sample(S,n.sampled)
  Area=0.1 #proportion of each cell surveyed
  n.species=1
	
  if((sqrt(S)%%1)!=0)cat("Error: S must be a square #")
	Adj=square_adj(sqrt(S))
	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q.pois=Matrix(tau.pois*Q)
  Q.bern=Matrix(tau.bern*Q)
	#simulate icar process
  Eta.pois=matrix(0,S,1)
  Eta.bern=Eta.pois
	SP.pois=array(0,dim=c(1,sqrt(S),sqrt(S)))
	SP.bern=SP.pois
  for(i in 1:1){
	  Eta.pois[,i]=rrw(Q.pois)
    Eta.bern[,i]=rrw(Q.bern)
    SP.pois[i,,]=matrix(Eta.pois[,i],sqrt(S),sqrt(S))
    SP.bern[i,,]=matrix(Eta.bern[,i],sqrt(S),sqrt(S))
  }
  if(ZIP==FALSE){
    SP.bern=0*SP.bern
	  Eta.bern=0*Eta.bern
  }

	#process parameters
	Lambda.grp=0.1
	X.site=cbind(rep(1,S),rep(log(c(1:sqrt(S)/sqrt(S))),each=sqrt(S))) #covariate on abundance intensity, sp 1
	Beta.pois=matrix(0,1,2)
	Beta.pois[1,]=c(log(100),0)
  X.bern=matrix(1,S,1)
  Beta.bern=matrix(0,1,1)
  if(ZIP==FALSE)Beta.bern=matrix(10,1,1)
	#plot(exp(X.site1%*%Beta.site1))
	#points(exp(X.site2%*%Beta.site2))
			
	#N=rpois(S,0.5*exp(X.site%*%Beta.site))
  N=matrix(0,1,S)
  for(i in 1:1)N[i,]=rbern(S,pnorm(X.bern%*%Beta.bern[i,]+Eta.bern[,i]))*(rpois(S,exp(X.site%*%Beta.pois[i,]+Eta.pois[,i])))
  N.sampled=matrix(0,5,n.sampled)
  for(i in 1:1)N.sampled[i,]=rbinom(n.sampled,N[i,Sampled.cells],Area)
  True.species=rep(1,sum(N.sampled[1,]))
 
  G.tot=apply(N,1,'sum')
	cat(paste("\n True G.tot= ",G.tot,'\n'))
		
	Dat=matrix(NA,sum(N.sampled),4)
  
	pl=1
	for(i in 1:5){
    for(j in 1:n.sampled){
      if(N.sampled[i,j]>0){
        for(k in 1:N.sampled[i,j]){
          Dat[pl,1]=j
          Dat[pl,2]=rbern(1,prop.photo)
          if(Dat[pl,2]==1){
            Dat[pl,3]=i
		        Dat[pl,4]=rpois(1,Lambda.grp[i])+1
          }
          pl=pl+1
        }
      }
    }
	}
          
          
	Dat=as.data.frame(Dat)
	colnames(Dat)=c("Transect","Photo","Species","Group")
	
  stacked.names=colnames(Dat)
	factor.ind=sapply(Dat[1,],is.factor)
	which.factors=which(factor.ind==1)
	n.factors=sum(factor.ind)
	
  Levels=NA
  if(n.factors>0){
    Factor.labels=vector("list",n.factors)
    for(i in 1:n.factors){
      Factor.labels[[i]]=levels(Dat[,which.factors[i]])
    }
    
    Dat.num=Dat
    #if(sum(Dat[,"Obs"]==0)>0)Dat.num[Dat[,"Obs"]==0,"Species"]=0  #if a missing obs, set species=0
    for(icol in which.factors){
      Dat.num[,icol]=as.numeric((Dat[,icol]))
    }
    Levels=vector("list",n.factors)
    for(i in 1:n.factors){
      Levels[[i]]=sort(unique(Dat.num[,which.factors[i]]))
    }
    names(Levels)=colnames(Dat[,which.factors])
  }
  if(misID==FALSE)Psi=NULL	
	Out=list(Dat=Dat,G.tot=G.tot,G.true=N,True.species=True.species,Sampled.cells=Sampled.cells,Psi=Psi)
}

