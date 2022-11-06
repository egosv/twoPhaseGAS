### These new functions do NOT consider the phi (dispersion parameter) as a quantity of inference
fnc_wgts <- function(formula,family,dat0,theta,q,disp,G_,Z_,uniqZ){
  N <- wg0 <- id_ <- NULL # suggested by https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  resp<-as.character(formula[[2]]); terms_q<-c(G_,Z_)
  pG_<-aggregate(as.formula(paste0("q~",paste(G_,collapse="+"))),q,FUN=sum); names(pG_)[ncol(pG_)]<-"q_g"
  dat0$ord<-1:nrow(dat0)
  dat0bis<-merge(as.data.table(dat0),as.data.table(q[,c(terms_q,"q")]),by=terms_q,all.x=T,sort = F)
  dat0bis<-merge(dat0bis,as.data.table(pG_),by=G_,all.x=T,sort = F)

  dat0bis$ind<-ifelse(dat0bis[[Z_]]%in%uniqZ,1,0)
  dat0bis$q<-dat0bis$ind*dat0bis$q + (1-dat0bis$ind)*dat0bis$q_g

  X0 = model.matrix(formula,dat0bis)
  eta = as.vector(if (NCOL(X0) == 1L){ X0 * theta }else X0 %*% theta)
  mu <- family$linkinv(eta)

  if(family$family=="gaussian"){
    dat0bis$wg0 <- dnorm(x=dat0bis[[resp]],mean=mu,sd=sqrt(disp))*dat0bis[["q"]]
  }else if(family$family=="binomial"){
    dat0bis$wg0 <- dbinom(x=dat0bis[[resp]],size=1,prob=mu)*dat0bis[["q"]]
  }else if(family$family=="poisson"){
    dat0bis$wg0 <- dpois(x=dat0bis[[resp]],lambda=mu)*dat0bis[["q"]]
  }else if(family$family=="Gamma"){
    dat0bis$wg0 <- dgamma(x=dat0bis[[resp]],shape=1/disp,scale=mu)*dat0bis[["q"]]
  }else stop("family not one of: gaussian, binomial, Gamma or poisson.")

  dat0bis[, N:=sum(wg0), by = id_]

  # dat0bis1 <- dat0bis

  # if(family$family=="gaussian"){
  #   dat0bis$lwg0 <- dnorm(x=dat0bis[[resp]],mean=mu,sd=sqrt(disp),log=TRUE) + log(dat0bis[["q"]])
  # }else if(family$family=="binomial"){
  #   dat0bis$lwg0 <- dbinom(x=dat0bis[[resp]],size=1,prob=mu,log=TRUE) + log(dat0bis[["q"]])
  # }else if(family$family=="poisson"){
  #   dat0bis$lwg0 <- dpois(x=dat0bis[[resp]],lambda=mu,log=TRUE) + log(dat0bis[["q"]])
  # }else if(family$family=="Gamma"){
  #   dat0bis$lwg0<-dgamma(x=dat0bis[[resp]],shape=1/disp,scale=mu,log=TRUE) + log(dat0bis[["q"]])
  # }else stop("family not one of: gaussian, binomial, Gamma or poisson.")
  #
  # dat0bis$wg0 <- exp(dat0bis$lwg0)
  #
  # dat0bis[, N:=sum(wg0), by = id_] ## fast way to add a column (N) with the sum of the weigths by id_
  #
  # #    ## Procedure to normalize probabilities under small values taken from http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
  # if( any(dat0bis$N==0) ){
  #   zerindx <- which(dat0bis$N==0)
  #   dfwg0 <- dat0bis[zerindx,c("id_","ord","lwg0"), with=FALSE]
  #   dfwg0[, max:=max(lwg0),by=id_]
  #   if( !all(dfwg0$ord==dat0bis$ord[zerindx]) ) stop("Something went wrong in the probability normalization process.")
  #   dat0bis$wg0[zerindx] <- exp(dfwg0$lwg0-dfwg0$max)
  #   dat0bis[zerindx, N:=sum(wg0), by = id_]
  # }
  #
  # dat0bis1$wg <- dat0bis1$wg0/dat0bis1$N
  # dat0bis1 <- dat0bis1[order(dat0bis1$ord),]

  dat0bis$wg <- dat0bis$wg0/dat0bis$N
  dat0bis <- dat0bis[order(dat0bis$ord),]

  return(dat0bis$wg)
  #return(cbind(dat0bis$wg,dat0bis1$wg,dat0bis$N,dat0bis1$N,dat0bis$wg0,dat0bis1$wg0))
}

.loglik<-function(theta,q,disp,formula,Y_,G_,Z_,dat,family){
  datbis=merge(dat,q,by=c(G_,Z_))
  In<-t(fac2sparse(datbis$id_))#class.ind2.Mat(datbis$id_)
  mu<-family$linkinv(model.matrix(formula,datbis)%*%theta)

  if(family$family=="gaussian"){
    sig<-sqrt(disp)
    ll <- sum( log( colSums( dnorm(datbis[,Y_], mu, sig)*datbis$q*In ) ) )
  }else if(family$family=="poisson"){
    ll <- sum( log( colSums( dpois(datbis[,Y_], lambda=mu)*datbis$q*In ) ) )
  }else if(family$family=="Gamma"){
    ll <- sum( log( colSums( dgamma(datbis[,Y_],shape=1/disp, scale=mu)*datbis$q*In ) ) )
  }else if(family$family=="binomial"){
    ll <- sum( log( colSums( dbinom(datbis[,Y_],size=1,prob=mu)*datbis$q*In ) ) )
  }else stop("family not one of: gaussian, binomial, Gamma, inverse.gaussian or poisson.")

  return(ll)
}

score_glm<-function(theta,disp,q,formula,Y_,data,family){
  q_j<-q[,"q"]; nq<-nrow(q)
  y <- data[,Y_]
  x <- as.matrix(model.matrix(formula,data))

  nobs <- NROW(y); nvars <- NCOL(x)

  variance <- family$variance
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta

  eta <- as.vector(if (NCOL(x) == 1L){ x * theta }else x %*% theta)
  mu <- linkinv(eta)
  varmu <- variance(mu)
  mu.eta.val <- mu.eta(eta)
  if (any(is.na(varmu))) #anyNA(varmu)
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  if (any(is.na(mu.eta.val)))
    stop("NAs in d(mu)/d(eta)")

  w <- (mu.eta.val^2)/varmu
  S <- x[, , drop = FALSE] * w * (y - mu)/(disp*mu.eta.val)

  #return(list(S=S,X=x[, , drop = FALSE],Wx2=w^2*varmu/(disp*mu.eta.val^2)))
  return(list(S=S,X=x[, , drop = FALSE],varmu=varmu,Y=y,mu=mu,w=w,mu.eta.val=mu.eta.val))
}

FIM.glm <- function(theta,qG,disp,formula,Y_,dat,family,glmfit){
  nq<-nrow(qG); nbetas<-length(theta)
  phi<-c(theta,log(qG[,"q"])) #send the log(q_j) here to deal with the numerical issue on the hessian computation
  n_phi<-length(phi)

  Qd1=sapply(dat$k_,function(k){1/qG[qG[,"k_"]==k,"q"]})
  Qd2=sapply(dat$k_,function(k){-1/qG[qG[,"k_"]==k,"q"]^2})
  Qmat = as.matrix(t(fac2sparse(dat$k_)))#class.ind2(dat$k_)
  Qmat.d1=Qmat*Qd1
  Psi=dat$wg_

  ## second derivatives (i.e. hessian of log P(Y_i|g_j) + log q_j or jacobian of the gradient)
  I1_0<-crossprod(qr.R(glmfit$qr))/disp # upper left matrix with the fisher information matrix for betas
  I2_0 <-diag(colSums(Psi * Qmat * Qd2)) #second derivatives for an i,jth observation

  I1<-rbind(cbind(I1_0, matrix(0, nrow=nrow(I1_0), ncol=ncol(I2_0))),
            cbind(matrix(0, nrow=nrow(I2_0), ncol=ncol(I1_0)), -I2_0))

  ## first derivatives (i.e. gradient of log P(Y_i|g_j,z_j) + log q_j)
  Sc=score_glm(theta,disp,qG,formula,Y_,dat,family)
  U<-colSums(Psi*cbind(Sc$S,Qmat.d1)) ##cbind(ptheta,pq_j)

  ## calculate second matrix, i.e. expected value of cross score functions
  X=Sc$X
  Wx2<- (Sc$w * (Sc$Y - Sc$mu)/(disp*Sc$mu.eta.val))^2
  Mx2<- crossprod( X,  Psi*Wx2*X )

  Wxq2<- Sc$w * (Sc$Y - Sc$mu)/(disp*Sc$mu.eta.val)
  Mxq2<-crossprod( X,  Psi*Wxq2*Qmat.d1)

  Mq2<- crossprod( Qmat.d1,  Psi*Qmat.d1 )

  obsI2<-rbind(cbind(Mx2, Mxq2),
               cbind(t(Mxq2),Mq2))

  ## Third matrix, i.e. outer product of expected score. Note this uses the W's from the second matrix and weight them accordingly
  IIt <- crossprod(fac2sparse(dat$id_))#tcrossprod(class.ind2.Mat(dat$id_)) #IIt <- tcrossprod(class.ind2.Mat(dat$id_)) #IIt<-tcrossprod(class.ind2(dat$id_)) ## much slower, specially for larga matrices

  #Mx3<- crossprod( sign(Wx2)*sqrt(abs(Wx2))*X*Psi, IIt %*% (sqrt(abs(Wx2))*X*Psi) )
  Mx3 <- as.matrix(crossprod( Wxq2*X*Psi, IIt %*% (Wxq2*X*Psi) ))
  #Mxq3<-crossprod( sign(Wxq2)*sqrt(abs(Wxq2))*X*Psi, IIt %*% (sqrt(abs(Wxq2))*Qmat.d1*Psi) )
  Mxq3<- as.matrix(crossprod( Wxq2*X*Psi, IIt %*% (Qmat.d1*Psi) ))
  Mq3<- as.matrix(crossprod( Qmat.d1*Psi, IIt %*% (Qmat.d1*Psi) ))

  obsI3<-as.matrix(rbind(cbind(Mx3, Mxq3),
                         cbind(t(Mxq3),Mq3)))

  obsFIM=I1-obsI2+obsI3
  return(list(FIM=obsFIM,I1=I1,I2=obsI2,I3=obsI3,U=U))
}

expFIM.glm<-function(theta,qG,disp,formula,Y_,dat,family,glmfit){ ## for this case Psi is the vector of weights (not the matrix)
  nq<-nrow(qG); nbetas<-length(theta)
  phi<-c(theta,log(qG[,"q"]))
  n_phi<-length(phi) # n<-nrow(Psi);

  Qd1=sapply(dat$k_,function(k){1/qG[qG[,"k_"]==k,"q"]})
  Qd2=sapply(dat$k_,function(k){-1/qG[qG[,"k_"]==k,"q"]^2})
  Qmat=as.matrix(t(fac2sparse(dat$k_)))#class.ind2(dat$k_)
  Qmat.d1=Qmat*Qd1
  Qmat.d2=Qmat*Qd2
  Psi=dat$wg_

  ## second derivatives (i.e. hessian of log P(Y_i|g_j) + log q_j or jacobian of the gradient)
  I1_0<-crossprod(qr.R(glmfit$qr))/disp # upper left matrix with the (expected) Fisher information matrix for betas
  I2_0<-colSums(Psi * Qmat.d2)
  I2_0 <-diag(I2_0,NROW(I2_0))
  #I2_0 <-diag(colSums(Psi * Qmat.d2)) ##old implementation (note that this brings a problem when the length  inside diag is equal to 1, i.e. one estimate for qG )

  expI1<-rbind(cbind(I1_0, matrix(0, nrow=nrow(I1_0), ncol=ncol(I2_0))),
               cbind(matrix(0, nrow=nrow(I2_0), ncol=ncol(I1_0)), -I2_0))

  ## first derivatives (i.e. gradient of log P(Y_i|g_j,z_j) + log q_j)
  Sc<-score_glm(theta,disp,qG,formula,Y_,dat,family)
  lst_indx<-nq+nbetas+1
  I2_0_<-cbind(Sc$S,Qmat.d1,Psi)
  U<-colSums(I2_0_[,lst_indx]*I2_0_[,-lst_indx])

  ## calculate second matrix, i.e. expected value of cross score functions
  X=Sc$X
  Wx2<- Sc$w^2*Sc$varmu/(disp*Sc$mu.eta.val^2)
  Mx2<- crossprod( X,  Psi*Wx2*X )

  Wxq2<- 0
  Mxq2<-crossprod( X,  Psi*Wxq2*Qmat.d1)

  Mq2<- crossprod( Qmat.d1,  Psi*Qmat.d1 )

  expI2<-rbind(cbind(Mx2, Mxq2),
               cbind(t(Mxq2),Mq2))

  ## Third matrix, i.e. outer product of expected score. Note this uses the W's from the second matrix and weight them accordingly
  IIt <- crossprod(fac2sparse(dat$id_)) #tcrossprod(class.ind2.Mat(dat$id_)) #IIt<-tcrossprod(class.ind2(dat$id_)) ## much slower, specially for larga matrices

  Mx3<- as.matrix(crossprod( sign(Wx2)*sqrt(abs(Wx2))*X*Psi, IIt %*% (sqrt(abs(Wx2))*X*Psi) ))
  Mxq3<-as.matrix(crossprod( sign(Wx2)*sqrt(abs(Wxq2))*X*Psi, IIt %*% (sqrt(abs(Wxq2))*Qmat.d1*Psi) ))
  Mq3<- as.matrix(crossprod( Qmat.d1*Psi, IIt %*% (Qmat.d1*Psi) ))

  expI3<-as.matrix(rbind(cbind(Mx3, Mxq3),
                         cbind(t(Mxq3),Mq3)))

  expFIM=expI1-expI2+expI3

  return(list(FIM=expFIM,I1=expI1,I2=expI2,I3=expI3,U=U))
}

Testing_EM_joint0 <- function(theta,q,Betas_ids=NULL,G_,FIM.obj){
  if(is.null(Betas_ids)) Betas_ids<- which(names(theta)%in%G_)
  nbetas<-length(theta)
  phi<-c(theta,q[,"q"])
  n_phi<-length(phi)

  Q1<-FIM.obj$FIM
  #Vcov matrix
  D<-rbind(diag(n_phi-1),0); F1 <-  t(D) %*% Q1 %*% D
  Omega1 <- tryCatch(solve( F1 ), error=function(e){ ginv(F1) })

  ###  Score test (following Lin's paper)
  IndU<-rep(F,n_phi); IndU[Betas_ids]<-T
  U<-FIM.obj$U
  U1<-U[IndU]

  Ind<-rep(F,n_phi-1); Ind[Betas_ids]<-T
  invF1_noInd<-tryCatch(solve(F1[!Ind,!Ind]), error=function(e){ ginv(F1[!Ind,!Ind]) })
  V1<-F1[Ind,Ind] - F1[Ind,!Ind] %*% invF1_noInd %*% F1[!Ind,Ind]

  invOmega1_Ind<- tryCatch(solve(Omega1[Ind,Ind]), error=function(e){ ginv(Omega1[Ind,Ind]) })
  W<-as.numeric(t(phi[IndU]) %*% invOmega1_Ind %*% phi[IndU])
  invV1<- tryCatch(solve( V1 ), error=function(e){ ginv(V1) })
  S<-as.numeric(t(U1) %*% invV1 %*% U1)

  return(list(Var=diag(Omega1)[1:nbetas],W=W,S=S,df=length(Betas_ids)))
}

Testing_EM_joint<-function(...){
  tryCatch(Testing_EM_joint0(...),error = function(e){list(Var=NA,W=NA,S=NA,df=NA)} )
}

# svdsubsel from http://bwlewis.github.io/GLM/svdss.html
#
# Input m*p matrix A, m >= p.
# Number of output columns k<=p.
# Returns an index subset of columns of A that *estimates* the k most linearly
# independent columns of A.

svdsubsel <- function(A,k=ncol(A))
{
  S <- svd(scale(A,center=FALSE,scale=TRUE), k)
  n <- which(svd(A)$d < 2*.Machine$double.eps)[1]
  if(!is.na(n) && k>=n)
  {
    k <- n - 1
    warning("k was reduced to match the rank of A")
  }
  Q <- qr( t(S$v[,1:k]) ,LAPACK=TRUE)
  sort(Q$pivot[1:k],decreasing=FALSE)
}

### taken from https://bwlewis.github.io/GLM/ and solves some inestability when zero weights are present and other edge cases in glm estimation (note this implementation does not handle offsets yet)
irls_svdnewton = function(A, b, family=binomial, maxit=25, tol=1e-08, weights=rep(1,nrow(A)), rank_deficiency=c("select columns","minimum norm","error")) {
  rank_deficiency = match.arg(rank_deficiency)
  m = nrow(A)
  n = ncol(A)
  select = 1:n
  zw = weights==0
  if(any(zw)) A[zw,]=0
  S = svd(A)
  tiny_singular_values = S$d/S$d[1] < tol
  k = sum(tiny_singular_values)
  if(k>0)
  {
    if(rank_deficiency=="select columns")
    {
      warning("Numerically rank-deficient model matrix")
      # NB This is a different selection method than R's default glm.fit uses.
      # See https://bwlewis.github.io/GLM and https://bwlewis/github.io/GLM/svdss.html
      select = svdsubsel(A,n-k)
      S = svd(A[,select,drop=FALSE])
    } else if(rank_deficiency=="error")
    {
      stop("Near rank-deficient model matrix")
    }
  }
  t = rep(0,m)
  s = rep(0,length(select))
  good = weights > 0
  for(j in 1:maxit)
  {
    g       = family$linkinv(t[good])
    varg    = family$variance(g)
    if(any(is.na(varg))) stop("NAs in variance of the inverse link function")
    if(any(varg==0)) stop("Zero value in variance of the inverse link function")
    gprime  = family$mu.eta(t[good])
    if(any(is.na(gprime))) stop("NAs in the inverse link function derivative")
    z       = rep(0,m)
    W       = rep(0,m)
    z[good] = t[good] + (b[good] - g) / gprime ### from glm.fit (may help with offsets z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good])
    W[good] = weights[good] * as.vector(gprime^2 / varg)

    good   = W > .Machine$double.eps*2
    if(sum(good)<m) warning("Tiny weights encountered")
    s_old   = s
    C   = chol(crossprod(S$u[good,,drop=FALSE], W[good]*S$u[good,,drop=FALSE]))
    s   = forwardsolve(t(C), crossprod(S$u[good,,drop=FALSE],W[good]*z[good]))
    s   = backsolve(C,s)
    t   = rep(0,m)
    t[good] = S$u[good,,drop=FALSE] %*% s
    if(sqrt(crossprod(s - s_old)) < tol) break
  }
  if( j>=maxit ) warning("Maximum number of irls iterations reached.")
  x = rep(NA, n)
  if(rank_deficiency=="minimum norm") S$d[tiny_singular_values] = Inf
  x[select] = S$v %*% ((1/S$d) * crossprod(S$u[good,],t[good]))
  names(x)=colnames(A)
  eta = as.vector(if (NCOL(A) == 1L){ A * x }else A %*% x)
  mu = family$linkinv(eta)
  residuals <- (b - mu)/family$mu.eta(eta)
  qr=qr(A*as.numeric(sqrt(W)))
  list(coefficients=x,iterations=j,residuals=residuals,qr=qr)
}
