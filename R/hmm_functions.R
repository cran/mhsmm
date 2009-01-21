
print.hmmspec <- function(x, ...){
  cat("Hidden Markov Model specification:\n")
  cat(sprintf("J (number of states): \n%i \n", x$J))
  cat("init:\n")
  print(x$init)
  cat ("transition:\n")
  print(x$transition)
  cat("emission:\n")
  print(x$emission)
  return(invisible(x))
}

summary.hmm <- function (object, ...) 
{
    cat("init: \n", round(object$model$init, 3), "\n")

    cat ("\ntransition:\n")
    print(round(object$model$trans, 3)) ## FIXME: This is **Very** fragile
    cat("\nemission:\n")
    print(object$model$emission)
    return(invisible(object))
}


simulate.hmmspec <- function(object, nsim, seed=NULL, ...) {

  if(!is.null(seed)) set.seed(seed)
  if(length(nsim)==1) {
   s1 = sim.mc(object$init,object$transition, nsim)
   x = sapply(s1, object$r, model=object)
   ret = list(s = s1, x = x, N = nsim)
   class(ret) <- "hsmm.data"
   ret
  }
  else  .sim.mhmm(object,nsim)
}

.sim.mhmm <- function(model,N) {
  s1 = sim.mc(model$init,model$transition,N)#the hidden states
  x = sapply(s1,model$r,model) #simulate the observed state sequence
  ret = list(s=s1,x=x,N=N)
  class(ret) <- "mhsmm.data"
  ret                                            
}


hmmspec <- function(init, trans, emission,r=rnorm.hsmm,f=dnorm.hsmm) {
 ans <- list(J=length(init), init = init, transition = trans, emission = emission,r=r,f=f)
 class(ans) <- "hmmspec"
 return(ans)
}

print.hmm <- function(x, ...) {
 cat("hmm object contains the following slots:\n")
 print(names(x))
 return(invisible(x))
}

hmmfit <- function(x,start.val,f=dnorm.hsmm,mstep=mstep.norm,tol=1e-08,maxit=1000) 
{
  model = start.val
  K = nrow(model$trans)
  if(class(x)=="numeric" | class(x)=="integer") {
  	warning('x is a primitive vector.  Assuming single sequence.')
  	N = NN = NROW(x)
  }
  else{
	  N = NN = x$N
  	x = x$x
  }
  
  if(K<2) stop("K must be larger than one.")	
  if(any(dim(model$trans)!=K)) stop("dimensions of a incorrect")
  if(length(model$init)!=K) stop("dimensions of st incorrect")
  if(NROW(x)!=sum(N)) stop("dimensions incorrect")
  if(length(N)==1) NN=1
  else NN = length(N)

  loglik=numeric(maxit)
  loglik=NA
  loglik[1]=-Inf  
  gam = double(K*sum(N))
  for(i in 1:maxit) {  
    p = sapply(1:K,fn <- function(state) f(x,state,model))

    #estep duh	
    test = .C("mo_estep_hmm",a=as.double(t(model$transition)),pi=as.double(t(model$init)),p=as.double(t(p)),N=as.integer(N),nsequences=as.integer(NN),
      K=as.integer(K),gam=gam,ll=double(1),PACKAGE='mhsmm')
    #mstep
    loglik[i]=test$ll                                   
  if(i>1)    if(abs(loglik[i]-loglik[i-1])<tol) break("Converged")
#    if((loglik[i]-loglik[i-1])<(-tol)) stop(paste("loglikelihood has decreased on iteration",i))
    gam = matrix(test$gam,ncol=K)
    model$emission = mstep(x,gam)
    model$transition=matrix(test$a,nrow=K,byrow=TRUE)
    model$init=test$pi
  }  
  ret = list(model=model,K=K,f=f,mstep=mstep,gam=gam,loglik=loglik[!is.na(loglik)],N=N,p=gam,yhat=apply(gam,1,which.max))
  class(ret) <- "hmm"
  return(ret)	
}


predict.hmm <- function(object,x,method="viterbi",...) {
  nseq=1  
  if(class(x)=="numeric" | class(x)=="integer") {
  	warning('x is a primitive vector.  Assuming single sequence.')
  	N = NROW(x)
  	NN = c(0,N)
  	if(N<1) stop("N less than one")
  }
  else{
  	N = x$N
  	NN = cumsum(c(0,x$N))
  	x0 = x$x
  }
#   state = integer(obj
  if(method=="viterbi") {
    K = object$K
    p = sapply(1:K,fn <- function(state) object$f(x$x,state,object$model))
    p[p==0]= 1e-200
    tmp = object$model$trans
    tmp[!tmp>0] =  1e-200
    logtrans = as.double(log(t(tmp)))
    tmp = object$model$init
    tmp[!tmp>0] =  1e-20
    logpi = as.double(log(t(tmp)))
    state = integer(sum(x$N))
    for(i in 1:length(x$N)) {    
         state[(NN[i]+1):NN[i+1]] = .C("viterbi_hmm",a=logtrans,pi=logpi,p=as.double(log(t(p[(NN[i]+1):NN[i+1],]))),N=as.integer(x$N[i]),NN=as.integer(nseq),K=as.integer(object$K),
                q=as.integer(rep(-1,x$N[i])),PACKAGE='mhsmm')$q+1
    }
    ans <- list(s=state,x=x$x,N=x$N)
  }
  else if(method=="smoothed") {
    tmp <- hmmfit(x,object$model,object$f,object$mstep,maxit=1)
    ans <- list(s=tmp$yhat,x=x$x,N=x$N,p=tmp$p)
  }
  else stop("Unavailable prediction method")
  class(ans) <- "hsmm.data"
  ans
#  void viterbi(double *a,double *start,double *p,int *T,int *nsequences,int *nstates,int *q) {
}


