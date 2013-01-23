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
  n.species=5
	
  if((sqrt(S)%%1)!=0)cat("Error: S must be a square #")
	Adj=square_adj(sqrt(S))
	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q.pois=Matrix(tau.pois*Q)
  Q.bern=Matrix(tau.bern*Q)
	#simulate icar process
  Eta.pois=matrix(0,S,5)
  Eta.bern=Eta.pois
	SP.pois=array(0,dim=c(5,sqrt(S),sqrt(S)))
	SP.bern=SP.pois
  for(i in 1:5){
	  Eta.pois[,i]=rrw(Q.pois)
    Eta.bern[,i]=rrw(Q.bern)
    SP.pois[i,,]=matrix(Eta.pois[,i],sqrt(S),sqrt(S))
    SP.bern[i,,]=matrix(Eta.bern[,i],sqrt(S),sqrt(S))
  }
  if(ZIP==FALSE)SP.bern=0*SP.bern
	
	#process parameters
	Lambda.grp=c(.1,.2,.3,.4,0)
	X.site=cbind(rep(1,S),rep(log(c(1:sqrt(S)/sqrt(S))),each=sqrt(S))) #covariate on abundance intensity, sp 1
	Beta.pois=matrix(0,5,2)
  Beta.pois[1,]=c(log(39),1) 
	Beta.pois[2,]=c(log(9),0)
	Beta.pois[3,]=c(log(20),-1)
	Beta.pois[4,]=c(log(10),0.5)
	Beta.pois[5,]=c(log(15),0)
  X.bern=matrix(1,S,1)
  Beta.bern=matrix(c(1,.5,0,.4,0),5,1)
  if(ZIP==FALSE)Beta.bern=matrix(10,5,1)
	#plot(exp(X.site1%*%Beta.site1))
	#points(exp(X.site2%*%Beta.site2))
			
	#N=rpois(S,0.5*exp(X.site%*%Beta.site))
  N=matrix(0,5,S)
  for(i in 1:5)N[i,]=rbern(S,pnorm(X.bern%*%Beta.bern[i,]+Eta.bern[,i]))*(rpois(S,exp(X.site%*%Beta.pois[i,]+Eta.pois[,i])))
  N.sampled=matrix(0,5,n.sampled)
  for(i in 1:5)N.sampled[i,]=rbinom(n.sampled,N[i,Sampled.cells],Area)
  True.species=rep(1,sum(N.sampled[1,]))
  for(i in 2:5)True.species=c(True.species,rep(i,sum(N.sampled[i,])))
 
  G.tot=apply(N,1,'sum')
	cat(paste("\n True G.tot= ",G.tot,'\n'))
	
  Cur.mat=matrix(0,n.species,3*(n.species+1))
  Alpha=matrix(0,n.species,n.species+1)
	C=array(0,dim=c(n.species,n.species+1,3))
  for(isp in 1:(n.species-1)){
    for(iobs in 1:(n.species-1)){
      if(isp==iobs){
        Alpha[isp,iobs]=0.6
        C[isp,iobs,]=c(0.5,0.35,0.15)
      }
      else{
        Alpha[isp,iobs]=0.08
        C[isp,iobs,]=c(0,0.3,0.7)
      } 
    }
  }
  C[n.species,n.species,]=c(1,0,0)
  C[1:(n.species-1),n.species+1,1]=1
  Alpha[n.species,n.species]=1
	Alpha[,n.species+1]=1-rowSums(Alpha)
	
  Psi=get_misID_mat(Cur.mat,Alpha,C,n.species,n.conf=3)
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
  	
	if(misID==TRUE){
		# Now, put in partial observation process
    for(i in 1:nrow(Dat)){
      if(Dat[i,2]==1)Dat[i,3]=sample(c(1:(3*(n.species+1))),1,prob=Psi[Dat[i,3],])
    }
	}
	Out=list(Dat=Dat,G.tot=G.tot,G.true=N,True.species=True.species,Sampled.cells=Sampled.cells,Psi=Psi)
}

