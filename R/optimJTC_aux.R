
check_alloc<-function(alloc_exp,Obs_samp){#Check for allocation higher than the observed elements in strata (if that exists keep all subjects in that strata and adjust to keep the sample size fixed)
  while( any(alloc_exp>Obs_samp) ){
    alloc_init<-alloc_exp
    indx1= alloc_exp>Obs_samp
    to_allocate<-sum( alloc_exp[indx1]-Obs_samp[indx1] )
    alloc_exp[indx1]<-Obs_samp[indx1]
    indx2<-!alloc_init>=Obs_samp & alloc_init>0
    alloc_exp[indx2]<- alloc_init[indx2] + floor(to_allocate*alloc_init[indx2]/sum(alloc_init[indx2] ))
  }
  alloc_exp
}


BSS <- function(samp.fracs,n,data,stratf)  {


  ### samp.fracs: vector of sampling fractions with length equal to the number of strata determined by stratf (and in the same order). Must add up to 1.
  ### n: phase 2 sample size
  ### data: phase 1 data, MUST contain stratification factors
  ### stratf: right hand side of a formula with stratification factors (Z and Yst in this case)


  tot_df=data.frame(xtabs(stratf,data=data))
  tot_df$n_strat=round(n*samp.fracs)


  # xtabs(n_stat~Z+S,data=tot_df)
  # xtabs(Freq~Z+S,data=tot_df)


  tot_df$n_strat<-check_alloc(tot_df$n_strat,tot_df$Freq)

  if( any(tot_df$n_strat>tot_df$Freq) ) stop("Sampling fraction greater than 1 for some stratum in the allocation")

  if( sum(round(n*samp.fracs))-sum(tot_df$n_strat)>2 )  {


    to_allocate<-sum(round(n*samp.fracs))-sum(tot_df$n_strat)
    alloc_init<-tot_df$n_strat
    indx <- tot_df$n_strat/tot_df$Freq<1 & !is.na(tot_df$n_strat/tot_df$Freq)
    q<-alloc_init[indx]/sum(alloc_init[indx])
    tot_df$n_strat[indx] <- alloc_init[indx] + ceiling(to_allocate*ifelse(is.nan(q),1/length(alloc_init[indx]),q))


  }


  tot_df$n_strat<-check_alloc(tot_df$n_strat,tot_df$Freq)

  ### Sampling begins
  #p2 <- tot_df$n_strat/tot_df$Freq
  data$R_ <- rep(0, NROW(data))
  stratf_<-all.vars(stratf)
  data$ord_ <-1:NROW(data)
  tot_df$k_<-1:NROW(tot_df)

  df=merge(data,tot_df[,c(stratf_,"k_")],by=stratf_,all.x=TRUE,sort=FALSE)

  for(k in 1:NROW(tot_df)){ #k<-3
    n_k=tot_df[k,"n_strat"]
    s <- sample(which( df$k_==k ), n_k)
    df$R_[s] <- 1
  }

  list(R=df[order(df$ord_),"R_"],eff.samp.fracs=tot_df$n_strat/n)


}


.optMsure<-function(type){
  Dopt<-function(vcov,k){
    if( is.null(k) ){
      return(det(vcov))
    }else{
      return(det(vcov[k,k]))
    }
  }
  Aopt<-function(vcov,k){
    if( is.null(k) ){
      return(sum(diag(vcov)))
    }else{
      return(sum(diag(vcov[k,k])))
    }
  }
  Param<-function(vcov,k){
    vcov[k,k]
  }
  switch(type,'D-opt'=Dopt,'A-opt'=Aopt,"Par-spec"=Param)
}


FIMf <- function(formula,family,dat,beta,p_gz,disp,G_,Z_){
  X <- model.matrix(formula, dat)
  eta = as.vector(if (NCOL(X) == 1L){ X * beta }else X %*% beta)
  mu = family$linkinv(eta)
  ## other quantities of interest
  varmu <- family$variance(mu)
  mu.eta.val <- family$mu.eta(eta)
  if (any(is.na(varmu))) #anyNA(varmu)
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  if (any(is.na(mu.eta.val)))
    stop("NAs in d(mu)/d(eta)")
  w <- (mu.eta.val^2)/varmu
  qr <- qr(X*as.numeric(sqrt(dat$wg_ * w)))

  #Create objects needed in the testing part
  # qG <- aggregate(as.formula(paste0("q~",paste(G_,collapse="+"))),p_gz,FUN=sum); names(qG)[which(names(qG)=="q")]<-"qG"
  qG <- p_gz; names(qG)[which(names(qG)=="q")]<-"qG"

  qG$k_ <- 1:nrow(qG)
  dat$ord_ <- 1:nrow(dat)

  # dat <- as.data.frame(merge(as.data.table(dat),as.data.table(qG[,c(G_,"qG","k_")]),by=G_, sort = F)) # merge(dat,qG[,c(G_,"qG","k_")],by=G_) ## this one may be sped up
  dat <- as.data.frame(merge(as.data.table(dat),as.data.table(qG[,c(G_,Z_,"qG","k_")]),by=c(G_,Z_), sort = F)) # merge(dat,qG[,c(G_,"qG","k_")],by=G_) ## this one may be sped up

  dat <- dat[order(dat$ord_),]

  Qd1=1/dat$qG
  Qd2=-1/dat$qG^2

  # Qmat <- crossprod(sapply(dat$k_,function(x)1*x==qG$k_), diag(nrow(qG))) ## Alternative way
  # Qmat.d1 <- Qmat*Qd1
  Qmat <- as.matrix(t(fac2sparse(dat$k_))) #class.ind2(dat$k_) ## old way
  Qmat.d1 = Qmat*Qd1
  Qmat.d2 = Qmat*Qd2

  # Qd1 <- 1/qG$qG
  # Qd2 <- -1/qG$qG^2

  # colsmatch <- c(G_,Z_)[c(G_,Z_) %in% names(dat)]
  # sapply(dat[,colsmatch],function(z){qG[,colsmatch] %in% z})

  # Ind <- sapply(dat[,Z_],function(z){qG[,Z_] %in% z})
  #
  # Qmat.d1 <- tcrossprod(model.matrix(~1,dat),Qd1)*t(Ind)
  # colnames(Qmat.d1) <- 1:nrow(qG)
  # Qmat.d2 <- tcrossprod(model.matrix(~1,dat),Qd2)*t(Ind)
  # colnames(Qmat.d2) <- 1:nrow(qG)

  Omega=dat$wg_

  ### Expected information
  {
    ## first matrix, or second derivatives (i.e. hessian of log P(Y_i|g_j) + log q_j or jacobian of the gradient)
    I1_0 <- crossprod(qr.R(qr))/disp # upper left matrix with the (expected) Fisher information matrix for betas
    I2_0 <- colSums(Omega * Qmat.d2)
    I2_0 <- diag(I2_0,NROW(I2_0))

    expI1 <- rbind(cbind(I1_0, matrix(0, nrow=nrow(I1_0), ncol=ncol(I2_0))),
                 cbind(matrix(0, nrow=nrow(I2_0), ncol=ncol(I1_0)), -I2_0))


    ## calculate second matrix, i.e. expected value of cross score functions
    X = X[,,drop=FALSE]
     # w * (Y - mu)/(disp*mu.eta.val)
    Wx2 <- w^2*varmu/(disp*mu.eta.val^2)
    Mx2<- crossprod( X,  Omega*Wx2*X )

    Wxq2 <- 0
    Mxq2 <- crossprod( X,  Omega*Wxq2*Qmat.d1)

    Mq2 <- crossprod( Qmat.d1,  Omega*Qmat.d1 )

    expI2 <- rbind(cbind(Mx2, Mxq2),
                 cbind(t(Mxq2),Mq2))


    ## Third matrix, i.e. outer product of expected score. Note this uses the W's from the second matrix and weight them accordingly
    IIt <- crossprod(fac2sparse(dat$id_)) #tcrossprod(class.ind2.Mat(dat$id_)) #IIt<-tcrossprod(class.ind2(dat$id_)) ## much slower, specially for larga matrices

    Mx3<- as.matrix(crossprod( sign(Wx2)*sqrt(abs(Wx2))*X*Omega, IIt %*% (sqrt(abs(Wx2))*X*Omega) ))
    Mxq3<-as.matrix(crossprod( sign(Wx2)*sqrt(abs(Wxq2))*X*Omega, IIt %*% (sqrt(abs(Wxq2))*Qmat.d1*Omega) ))
    Mq3<- as.matrix(crossprod( Qmat.d1*Omega, IIt %*% (Qmat.d1*Omega) ))

    expI3<-as.matrix(rbind(cbind(Mx3, Mxq3),
                           cbind(t(Mxq3),Mq3)))

    expFIM=expI1-expI2+expI3
    return(list(FIM=expFIM,expI1=expI1,expI2=expI2,expI3=expI3))
  }
}

# another version for obsIMf1 that calculates the observed information by individual faster that
# calling obsIMf1 each time (needed for optimTP.GA), note that by.id=FALSE should render the
# exact same results as obsIMf1
obsIMf1 <- function(formula,family,dat,beta,p_gz,disp,G_,Z_,by.id=FALSE){
  Y <- dat[,as.character(formula[[2]])]
  X <- model.matrix(formula, dat)
  eta = as.vector(if (NCOL(X) == 1L){ X * beta }else X %*% beta)
  mu = family$linkinv(eta)
  ## other quantities of interest
  varmu <- family$variance(mu)
  mu.eta.val <- family$mu.eta(eta)

  if (any(is.na(varmu))) #anyNA(varmu)
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  if (any(is.na(mu.eta.val)))
    stop("NAs in d(mu)/d(eta)")
  Omega <- dat$wg_
  w <- (mu.eta.val^2)/varmu

  # qG <- aggregate(as.formula(paste0("q~",paste(G_,collapse="+"))),p_gz,FUN=sum); names(qG)[which(names(qG)=="q")]<-"qG"
  qG <- p_gz; names(qG)[which(names(qG)=="q")]<-"qG"

  qG$k_ <- 1:nrow(qG)

  # dat$ord_ <- 1:nrow(dat)
  # # dat <- as.data.frame(merge(as.data.table(dat),as.data.table(qG[,c(G_,"qG","k_")]),by=G_, sort = F)) # merge(dat,qG[,c(G_,"qG","k_")],by=G_) ## this one may be sped up
  # dat <- as.data.frame(merge(as.data.table(dat),as.data.table(qG[,c(G_,Z_,"qG","k_")]),by=c(G_,Z_), sort = F)) # merge(dat,qG[,c(G_,"qG","k_")],by=G_) ## this one may be sped up
  # dat <- dat[order(dat$ord_),]
  # ## Need the 3 lines below to always have the same number of columns in Qmat irrespective of the subset
  # Qd1 <- 1/dat$qG
  # Qmat0 <- crossprod(sapply(dat$k_,function(x)1*x==qG$k_), diag(nrow(qG)))
  # Qmat.d1bis <- Qmat0*Qd1

  Qmat <- tcrossprod(matrix(rep(1,nrow(dat)),ncol=1),1/qG$qG)
  matchrows <- if( all(c(G_,Z_)%in%names(dat)) & all(c(G_,Z_)%in%names(qG)) ){
    match(interaction(dat[,c(G_,Z_)]), interaction(qG[,c(G_,Z_)]))
  }else if( Z_ %in% names(qG) & Z_ %in% names(dat) ){
    match(dat[,c(Z_)],qG[,c(Z_)])
  }else if( G_ %in% names(qG) & G_ %in% names(dat)){
    match(dat[,c(G_)],qG[,c(G_)])
  }else stop("No matching columns b/t p_gz and dat")

  Qmat.d1 <- Qmat*t(sapply(matchrows,function(s){1*s==1:nrow(qG)}))
  # Qmat.d1b <- crossprod(sapply(dat[,Z_],function(z){1*(p_gz[,Z_]==z)/p_gz[,"q"]}), diag(nrow(p_gz)))
  colnames(Qmat.d1) <- 1:nrow(qG)

  ## Original way assuming all the values of qG are in dat$k_
  # Qmat <- as.matrix(t(fac2sparse(dat$k_)))   #class.ind2(dat$k_)
  # Qmat.d1 <- Qmat*Qd1

  X <- X[, , drop = FALSE]

  Wxq2 <- w * (Y - mu)/(disp*mu.eta.val)

  U <- cbind(Wxq2*X,Qmat.d1)

  Ind_id <- fac2sparse(dat$id_)
  EU <- Ind_id %*% (Omega*U)

  if( !by.id ){
    ## Second matrix, i.e. expected value of cross score functions

    obsI2 <- crossprod( sqrt(Omega)*U,  sqrt(Omega)*U )

    # Mx2 <- crossprod( Wxq2*X,  Omega*Wxq2*X )
    # Mxq2 <- crossprod( Wxq2*X,  Omega*Qmat.d1)
    #
    # Mq2 <- crossprod( Qmat.d1,  Omega*Qmat.d1 )
    #
    # obsI2b <- rbind(cbind(Mx2, Mxq2),
    #                cbind(t(Mxq2),Mq2))

    ## Third matrix, i.e. outer product of expected scores.
    ## Expected score
    obsI3 <- as.matrix(crossprod( EU,  EU ))

    # IIt <- crossprod(Ind_id) #tcrossprod(class.ind2.Mat(dat$id_)) #IIt<-tcrossprod(class.ind2(dat$id_)) ## much slower, specially for large matrices
    #
    # signWq2 <- sign(Wxq2)
    # sqrtWq2 <- sqrt(abs(Wxq2))
    #
    # Mx3 <- as.matrix(crossprod( Wxq2*X*Omega, IIt %*% (Wxq2*X*Omega) ))
    # Mxq3 <- as.matrix(crossprod( signWq2*sqrtWq2*X*Omega, IIt %*% (sqrtWq2*Qmat.d1*Omega) ))
    # Mq3 <- as.matrix(crossprod( Qmat.d1*Omega, IIt %*% (Qmat.d1*Omega) ))
    #
    # obsI3b <- as.matrix(rbind(cbind(Mx3, Mxq3),
    #                        cbind(t(Mxq3),Mq3)))
    ### all.equal(obsI3,obsI3b)
    # obsI3 <- as.matrix(crossprod( U*Omega, IIt %*% U*Omega) ) ## not symmetric
  }else{
    ### this function implement row-wise outer product
    rwouter <- function(mat){t(apply( mat, 1, function(x) x%o%x))}
    ## calculate second matrix, i.e. expected value of cross score functions
    obsI2 <- as.matrix( Ind_id %*% rwouter( U*sqrt(Omega) ) )
    ## Third matrix, i.e. outer product of expected score
    obsI3 <- rwouter( EU )
  }

  return(list(I2=obsI2,I3=obsI3))
}

obsIMf2 <- function(formula,family,dat,beta,p_gz,disp,G_,Z_,by.id=FALSE){
  #### use enrich package to do this for now
  linkf <- enrich(make.link(family$link), with = "all")
  family <- enrich(family, with = "all")

  Y <- dat[,as.character(formula[[2]])]
  X <- model.matrix(formula, dat)
  eta = as.vector(if (NCOL(X) == 1L){ X * beta }else X %*% beta)
  mu = family$linkinv(eta)
  ## other quantities of interest
  varmu <- family$variance(mu)
  mu.eta.val <- family$mu.eta(eta)
  d2mus <- linkf$d2mu.deta(eta)
  d1variances <- family$d1variance(mu)

  if (any(is.na(varmu))) #anyNA(varmu)
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  if (any(is.na(mu.eta.val)))
    stop("NAs in d(mu)/d(eta)")
  Omega <- dat$wg_
  w <- (mu.eta.val^2)/varmu
  qrX <- qr(X*as.numeric(sqrt(Omega * w)))
  w1 <- Omega * (d2mus / varmu - mu.eta.val^2 * d1variances / varmu^2) * (Y - mu)

  # qG <- aggregate(as.formula(paste0("q~",paste(G_,collapse="+"))),p_gz,FUN=sum); names(qG)[which(names(qG)=="q")]<-"qG"
  qG <- p_gz; names(qG)[which(names(qG)=="q")]<-"qG"

  ###
  # Qd1 <- 1/qG$qG
  # Qmat.d1 <- tcrossprod(matrix(rep(1,nrow(dat)),ncol=1),Qd1) #wrong since doesn not discriminate by Z, all values are present

  # Qmat <- as.matrix(t(fac2sparse(dat$k_)))#class.ind2(dat$k_)
  # Qmat <- crossprod(sapply(dat$k_,function(x)1*x==qG$k_), diag(nrow(qG)))
  # Qmat.d1 <- Qmat*Qd1
  # Qmat.d1 <- tcrossprod(model.matrix(~1,dat),Qd1)
  # Qmat <- tcrossprod(matrix(rep(1,nrow(dat)),ncol=1),1/qG$qG)
  # matchrows <- if( all(c(G_,Z_)%in%names(dat)) & all(c(G_,Z_)%in%names(qG)) ){
  #   match(interaction(dat[,c(G_,Z_)]), interaction(qG[,c(G_,Z_)]))
  # }else if( Z_ %in% names(qG) & Z_ %in% names(dat) ){
  #   match(dat[,c(Z_)],qG[,c(Z_)])
  # }else if( G_ %in% names(qG) & G_ %in% names(dat)){
  #   match(dat[,c(G_)],qG[,c(G_)])
  # }else stop("No matching columns b/t p_gz and dat")
  # Qmat.d1 <- Qmat*t(sapply(matchrows,function(s){1*s==1:nrow(qG)}))

  Qmat.d1 <- crossprod(sapply(dat[,Z_],function(z){1*(p_gz[,Z_]==z)/p_gz[,"q"]}), diag(nrow(p_gz)))

  colnames(Qmat.d1) <- 1:nrow(qG)
  if( !by.id ){
    ## calculate second matrix, i.e. expected value of cross score functions
    X <- X[, , drop = FALSE]

    Wxq2 <- w * (Y - mu)/(disp*mu.eta.val)

    # U <- cbind(Wxq2*X,Qmat.d1)
    #
    # obsI2 <- crossprod( U,  Omega*U )

    Mx2 <- crossprod( Wxq2*X,  Omega*Wxq2*X )
    Mxq2 <- crossprod( X,  Omega*Wxq2*Qmat.d1)

    Mq2 <- crossprod( Qmat.d1,  Omega*Qmat.d1 )
    Mq2 <- diag(diag(Mq2),nrow=nrow(Mq2))

    obsI2 <- rbind(cbind(Mx2, Mxq2),
                   cbind(t(Mxq2),Mq2))

    ## Third matrix, i.e. outer product of expected score.
    IIt <- crossprod(fac2sparse(dat$id_)) #tcrossprod(class.ind2.Mat(dat$id_)) #IIt<-tcrossprod(class.ind2(dat$id_)) ## much slower, specially for large matrices
    signWxq2 <- sign(Wxq2)
    sqrtWxq2 <- sqrt(abs(Wxq2))

    Mx3 <- as.matrix(crossprod( Wxq2*X*Omega, IIt %*% (Wxq2*X*Omega) ))
    Mxq3 <- as.matrix(crossprod( signWxq2*sqrtWxq2*X*Omega, IIt %*% (sqrtWxq2*Qmat.d1*Omega) ))
    Mq3 <- as.matrix(crossprod( Qmat.d1*Omega, IIt %*% (Qmat.d1*Omega) ))

    obsI3 <- as.matrix(rbind(cbind(Mx3, Mxq3),
                             cbind(t(Mxq3),Mq3)))

    # obsI3 <- as.matrix(crossprod( U*Omega, IIt %*% U*Omega) ) ## not symmetric

    return(list(I2=obsI2,I3=obsI3))
  }else{

    X <- X[, , drop = FALSE]

    ### these functions implement row-wise outer product
    rwouter <- function(mat){t(apply( mat, 1, function(x) x%o%x))}
    rwouterxy <- function(mat1,mat2){nx=NCOL(mat1);ny=NCOL(mat2);t(apply( cbind(mat1,mat2), 1, function(x) x[1:nx]%o%x[(nx+1):(nx+ny)]))}

    ## calculate second matrix, i.e. expected value of cross score functions
    Wxq2 <- w * (Y - mu)/(disp*mu.eta.val)
    X_ <- cbind(Wxq2*X*sqrt(Omega),Qmat.d1*sqrt(Omega))
    obsI2 <- do.call(rbind,by(rwouter( X_),list(id_=dat$id_),colSums))

    ## Third matrix, i.e. outer product of expected score. Note this uses the W's from the second matrix and weight them accordingly
    IIt <- crossprod(fac2sparse(dat$id_))

    X_ <- cbind(Wxq2*X*Omega,Qmat.d1*Omega)

    obsI3 <- do.call(rbind,by(rwouterxy(  X_, as.matrix(IIt %*%  X_)),list(id_=dat$id_),colSums))

    return(list(I2=obsI2,I3=obsI3))
  }
}

Dmu_fnc <- function(formula,family,dat,beta,disp){
  #### use enrich package to do this for now
  Y <- dat[,as.character(formula[[2]])]
  X <- model.matrix(formula, dat)
  eta = as.vector(if (NCOL(X) == 1L){ X * beta }else X %*% beta)
  mu = family$linkinv(eta)
  ## other quantities of interest
  varmu <- family$variance(mu)
  mu.eta.val <- family$mu.eta(eta)

  if (any(is.na(varmu))) #anyNA(varmu)
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  if (any(is.na(mu.eta.val)))
    stop("NAs in d(mu)/d(eta)")

  w <- (mu.eta.val^2)/varmu

  Dmu <- w * (Y - mu)/(disp*mu.eta.val)

  return(Dmu)
}

weigthedloglik <-  function(theta,nbetas,q,disp,formula,Y_,G_,Z_,dat,family)  {
  betas <- theta[1:nbetas]
  q$qG_ <- theta[(nbetas+1):length(theta)]
  datbis <- merge(dat,q,by=c(G_))
  In <- t(fac2sparse(datbis$id_))#class.ind2.Mat(datbis$id_)
  mu <- family$linkinv(model.matrix(formula,datbis)%*%betas)

  if(family$family=="gaussian"){
    sig<-sqrt(disp)
    ll <- sum(  datbis$wg_*log( colSums( dnorm(datbis[,Y_], mu, sig)*datbis$qG_*In ) ) )

  }else if(family$family=="binomial"){
    ll <- sum( datbis$wg_*log( colSums( dbinom(datbis[,Y_],size=1,prob=mu)*datbis$qG_*In ) ) )
  }else stop("family not one of: gaussian, binomial, Gamma, inverse.gaussian or poisson.")

  return(ll)
}

# dat0=dat; q0=p_gz
# dat0=dat_; q0=p_gz
omega1_g <- function(formula,family,dat0,beta,q0,disp,G_,Z_){

  N <- wg0 <- id_ <- lwg0 <- NULL # as suggested by https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/

  resp<-as.character(formula[[2]]); terms_q<-c(G_,Z_)
  dat0$ord<-1:nrow(dat0)
  q_ <- names(q0)[!names(q0) %in% terms_q]
  if( q_ %in% names(dat0) ){
    namsold<-names(dat0)[which(names(dat0)==q_)]
    names(dat0)[which(names(dat0)==q_)]<-paste0(q_,".old")
  }
  dat0bis<-merge(as.data.table(dat0),as.data.table(q0[,c(terms_q,"q")]),by=terms_q,all.x=T,sort = F)

  X0 = model.matrix(formula,dat0bis)
  eta = as.vector(if (NCOL(X0) == 1L){ X0 * beta }else X0 %*% beta)
  mu <-family$linkinv(eta)

  if(family$family=="gaussian"){
    dat0bis$lwg0 <- dnorm(x=dat0bis[[resp]],mean=mu,sd=sqrt(disp), log=TRUE) + log(dat0bis[["q"]])
  }else if(family$family=="binomial"){
    dat0bis$lwg0 <- dbinom(x=dat0bis[[resp]],size=1,prob=mu, log=TRUE) + log(dat0bis[["q"]])
  }else if(family$family=="poisson"){
    dat0bis$lwg0 <- dpois(x=dat0bis[[resp]],lambda=mu, log=TRUE) + log(dat0bis[["q"]])
  }else if(family$family=="Gamma"){
    dat0bis$lwg0<-dgamma(x=dat0bis[[resp]],shape=1/disp,scale=mu, log=TRUE) + log(dat0bis[["q"]])
  }else stop("family not one of: gaussian, binomial, Gamma or poisson.")

  dat0bis$wg0 <- exp(dat0bis$lwg0)

  dat0bis[, N:=sum(wg0,na.rm = TRUE), by = id_] ## fast way to add a column (N) with the sum of the weigths by id_

  ##    ## Procedure to normalize probabilities under small values taken from http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
  if( any(dat0bis$N==0) ){
    zerindx <- which(dat0bis$N==0)
    dfwg0 <- dat0bis[zerindx,c("id_","ord","lwg0"), with=FALSE]
    dfwg0[, max:=max(lwg0),by=id_]
    if( !all(dfwg0$ord==dat0bis$ord[zerindx]) ) stop("Something went wrong in the probability normalization process.")
    dat0bis$wg0[zerindx] <- exp(dfwg0$lwg0-dfwg0$max)
    dat0bis[zerindx, N:=sum(wg0), by = id_]
  }

  dat0bis$wg <- dat0bis$wg0/dat0bis$N
  dat0bis <- dat0bis[order(dat0bis$ord),]

  return(list(wg=dat0bis$wg,pYZ=dat0bis$N))
  #return(dat0bis$wg)
}

# function fitnessTP
fitnessTP  <-  function(obsIM.R,Rj,optimMeasure,K.idx=NULL)  {

  FIM <- obsIM.R$IM1 + matrix(colSums(obsIM.R$IM2_id*Rj),ncol=obsIM.R$npars) - matrix(colSums(obsIM.R$IM3_id*Rj),ncol=obsIM.R$npars)

  D <- rbind(diag(obsIM.R$npars-1),0)
  F1 <-  crossprod(D, FIM %*% D)

  ### Function that determines the optimality measure (D-, A- or parameter-specfic for now)
  ### if optimMeasure=="Par-spec" then K.idx MUST be not null
  if( optimMeasure=="Par-spec" & is.null(K.idx) ) stop("For a parameter-specific criterion K.idx must be provided.")

  ## asign these functions to the execution environment
  optMsuref=.optMsure(optimMeasure)

  if( obsIM.R$Null ){
    Ind <- obsIM.R$Ind

    invF1_noInd <- ginv(F1[!Ind,!Ind])
    V1 <- F1[Ind,Ind] - F1[Ind,!Ind] %*% invF1_noInd %*% F1[!Ind,Ind]
    invV1 <- solve( V1 )

    return( optMsuref(invV1,K.idx) )

  }else {
    vcov = solve(F1)

    return(optMsuref(vcov,K.idx))
  }

}

# function obsIM.R
obsIM.R  <-  function(formula,miscov,auxvar,family,data,beta,p_gz,disp=NULL)  {
  ### Indicators (0,1) for phase 2 belonging, same length as nrow(data)
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  terms_dat <- unique( c(all.vars(formula), attr(terms(miscov), "term.labels") , attr(terms(auxvar), "term.labels")) )
  Y_<-as.character(formula[[2]])
  G_<-all.vars(miscov)
  Z_<-all.vars(auxvar)

  ## remove enties with zero or negative probabilities
  p_gz <- p_gz[p_gz$q>0,]

  ### check if names in G_ and Z_ are in the joint ditribution data.frame p_gz
  g.vals<-unique(p_gz[,G_]);
  N <- NROW(data)

  if( is.null(disp) ){
    if( family$family %in% c("gaussian","Gamma") ){
      disp = var(data[,Y_])
    }else disp=1
  }

  dat_ <- cbind(id_=rep(seq(1,N,by=1), each=NROW(g.vals)),data[rep(seq_len(N), each=NROW(g.vals)), ],matrix(rep(t(g.vals),N), ncol=NCOL(g.vals), byrow=TRUE))
  names(dat_)[(NCOL(dat_)-NCOL(g.vals)+1):ncol(dat_)] <- G_
  dat_ <- droplevels(dat_)

  wgs <- omega1_g(formula,family,dat_,beta,p_gz,disp,G_,Z_)

  dat_$wg_ <- wgs$wg
  dat_ <- dat_[!is.na(dat_$wg_),]

  IM1 <- obsIMf1(formula,family,dat_,beta,p_gz,disp,G_,Z_)$I3/N

  IM23_id <- obsIMf1(formula,family,dat_,beta,p_gz,disp,G_,Z_,by.id = TRUE)

  IM2_id <- IM23_id$I2/N

  IM3_id <- IM23_id$I3/N

  npars <- NCOL(IM1)

  Ind <- rep(FALSE,npars-1)
  Betas_ids <- which(grepl(paste0(paste0("^",G_),collapse="|"),colnames(model.matrix(formula,dat_[1,]))))
  Ind[Betas_ids] <- TRUE

  Null <- all(beta[Ind]==0)

  return( list(IM1=IM1,IM2_id=IM2_id,IM3_id=IM3_id,npars=npars, Ind=Ind, Null=Null) )
}


## p_Z <- data.frame(xtabs(~Z,dat_sim)/nrow(dat_sim)); LD.r <- 0.75; maf_G <- 0.2
### This function calculates p_gz given a distribution for Z (using, say xtabs), a MAF for G and a LD value (R).

# This function calculates p_gz given a distribution for Z (using, say xtabs), a MAF for G and a LD value (R).
# INTERNAL FUNCTION
p_gz_func <- function(p_Z, maf_G, LD.r, G_, Z_)  {
  ### I'm gonna assume here that G has 3 level genotypes (GG, Gg and gg), I will look further for generalizations later. Also, this only takes into account a single G. Further extensions for multiple G's are warranted.
  nZ <- NROW(p_Z)
  P_d <- maf_G; ## MAF of G
  P_a <- if(nZ==1){
    p_Z
  } else ifelse(nZ==3, p_Z[p_Z[,1]==2,"Freq"] + p_Z[p_Z[,1]==1,"Freq"]/2, p_Z[p_Z[,1]==1,"Freq"]) ## MAF of Z
  P_D <- 1 - P_d; P_A <- 1 - P_a ;

  Rsq <- LD.r^2
  LD <- LD.r * sqrt(P_D*P_d*P_A*P_a)
  # Dpr <- LD /( (LD<0)*min(P_D*P_A,P_d*P_a) + (LD>=0)*min(P_d*P_A,P_D*P_a) )
  #
  # if( Dpr > 1 | Dpr < -1 ) stop("This combination of MAF for G and genotype for Z seems incompatible. D' > 1 or D'< -1.")
  Dpr <- LD /( (LD<0)*max(-P_D*P_A,-P_d*P_a) + (LD>=0)*min(P_D*P_a,P_d*P_A) )
  if( Dpr > 1 | Dpr < 0 ) stop("This combination of MAF for G and genotype for Z seems incompatible. D' > 1 or D'< 0.")

  ## Haplotype frequencies (two SNPs)
  h.freqs <- rep(0, 4)
  h.freqs[1] <- LD + (1-P_d)*(1-P_a) ## DA
  h.freqs[2] <- 1 - P_d - h.freqs[1] ## Da
  h.freqs[3] <- 1 - P_a - h.freqs[1] ## dA
  h.freqs[4] <- LD + P_d*P_a ##  da
  names(h.freqs) <- c()

  ## Genotype frequencies (two SNPs)
  cross <- h.freqs %*% t(h.freqs)
  genot <- matrix(rep(0,3*3), nrow=3, ncol=3)
  genot[1,1] <- cross[1,1]
  genot[1,2] <- cross[1,2] + cross[2,1]
  genot[1,3] <- cross[2,2]
  genot[2,1] <- cross[1,3] + cross[3,1]
  genot[2,2] <- cross[1,4] + cross[2,3] + cross[3,2] + cross[4,1]
  genot[2,3] <- cross[2,4] + cross[4,2]
  genot[3,1] <- cross[3,3]
  genot[3,2] <- cross[3,4] + cross[4,3]
  genot[3,3] <- cross[4,4]
  rownames(genot) <- 0:2
  colnames(genot) <- 0:2
  crosstable <- data.frame(genot,check.names = FALSE)
  p_gz <- data.frame(rows = rownames(crosstable), stack(crosstable))
  names(p_gz)[names(p_gz)=="rows"] <- G_
  names(p_gz)[names(p_gz)=="ind"] <- Z_
  names(p_gz)[names(p_gz)=="values"] <- "q"
  p_gz[,G_] <- as.numeric(as.character(p_gz[,G_]))
  p_gz[,Z_] <- as.numeric(as.character(p_gz[,Z_]))
  return(p_gz[,c(G_,Z_,"q")])
}


### taken from gtools to avoid dependency
perms <- function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE){
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) !=
      0)
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) !=
      0)
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n)
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n)
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed)
    sub <- function(n, r, v) {
      if (r == 1)
        matrix(v, n, 1)
      else if (n == 1)
        matrix(v, 1, r)
      else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner),
                                                  ncol = ncol(inner), nrow = nrow(inner) * n,
                                                  byrow = TRUE))
      }
    }
  else sub <- function(n, r, v) {
    if (r == 1)
      matrix(v, n, 1)
    else if (n == 1)
      matrix(v, 1, r)
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n -
                                                        1, r - 1, v[-i])))
      X
    }
  }
  sub(n, r, v[1:n])
}


# This function calculates multilocus genotype frequencies (p_gz) in which each individual locus is biallelic given allele frequencies and linkage disequilibrium among loci. Adapted from the original by Liam J. Revell 2018
## Osvaldo Espin-Garcia, 2021
# @param nloci number of locus involved or a vector with SNP names.
# @param p  (Major) allele frequencies of the dominant (in this case, merely uppercase) allele at each locus.
# @param LD Linkage disequilibrium info across pairs or loci. Flat LD matrix with as many rows as unique pairs of nloci. It must contain 3 columns: locus1, locus2 and r value (Pearson correlation).
# INTERNAL FUNCTION
multilocus.Pgz <- function(nloci=2,p=NULL,LD=NULL){
  if(is.null(p)) p <- rep(0.5, nloci)
  if(is.null(LD)) LD <- cbind(t(combn(1:nloci,2)), 0)
  if(length(p)!=nloci) nloci<-length(p)
  allelesfreqs <- cbind(p,1-p)
  COMBN <- perms(n=2,r=nloci,set=T,repeats.allowed=T)
  ALLELESp <- LETTERS[1:nloci]
  ALLELESq <- letters[1:nloci]
  ALLELESnams <- vector()
  FREQ <- rep(1,nrow(COMBN))
  for( i in 1:nrow(COMBN) ){
    for( j in 1:nloci ){
      FREQ[i] <- FREQ[i]*allelesfreqs[j,COMBN[i,j]]

      atype <- if(COMBN[i,j]==1) ALLELESp[j]
      else if(COMBN[i,j]==2) ALLELESq[j]

      ALLELESnams[i] <- if(j==1) atype else paste(ALLELESnams[i],atype,sep="")
    }
  }
  ### Add/Subtract LD
  FREQ2 <- FREQ
  pairs <- t(combn(1:nloci,2))
  if( nrow(pairs)!=nrow(LD) ){
    stop("The number of rows in LD matrix should correspond to the number of possible combinations: nrows(LD):",nrow(LD),"!=", nrow(pairs), " (num. possible combinations).")
  }
  for( k in 1:nrow(pairs) ){
    LD.pair <- LD[pairs[k,1]==LD[,1] & pairs[k,2]==LD[,2], 3]
    kthfreqs <- allelesfreqs[t(pairs[k,]),]
    LD.pair <- LD.pair * sqrt(prod(kthfreqs))

    Dpr <- LD.pair / ( (LD.pair<0)*max(-prod(kthfreqs[,1]),-prod(kthfreqs[,2])) + (LD.pair>=0)*min(prod(diag(kthfreqs)),prod(1-diag(kthfreqs))) )
    if( Dpr > 1 | Dpr < 0 ) message("One combination of MAFs and LD seems incompatible. i.e. D' > 1 or D'< 0. The frequencies will be adjusted/normalized in the results accordingly. Zero frequencies can be expected. SNP Pair: ", paste0(pairs[k,],collapse=", "))

    ind <- COMBN[,pairs[k,1]]==COMBN[,pairs[k,2]]
    FREQ2 <- FREQ2 + ind*LD.pair - (1-ind)*LD.pair
  }

  ### Generate Genotypes
  cross <- FREQ2 %*% t(FREQ2)
  colnames(cross) <- rownames(cross) <- ALLELESnams
  crossnames <- matrix(do.call(paste0, expand.grid(ALLELESnams,ALLELESnams)), ncol=length(ALLELESnams))
  crossnames_sorted <- apply(crossnames,c(1,2), function(x) sapply(lapply(strsplit(x, NULL), sort), paste, collapse=""))

  GENOTYPE <- vector()
  GENO <- perms(n=3,r=nloci,set=T,repeats.allowed=T)
  qval <- rep(0, nrow(GENO))
  for( i in 1:nrow(GENO) ){
    for( j in 1:nloci ){

      gtype <- if(GENO[i,j]==1) paste(rep(ALLELESp[j],2),collapse="")
      else if(GENO[i,j]==2) paste(c(ALLELESp[j],ALLELESq[j]),collapse="")
      else if(GENO[i,j]==3) paste(rep(ALLELESq[j],2),collapse="")

      GENOTYPE[i]<-if(j==1) gtype else paste(GENOTYPE[i],gtype,sep="")
    }
    ### match strings
    indMat <- which(crossnames_sorted %in% sapply(lapply(strsplit(GENOTYPE[i], NULL), sort), paste, collapse=""))
    # cat(indMat,"\n")
    qval[i] <- qval[i] + sum(cross[indMat])
  }

  ### Ensure values are between 0 and 1
  if( any(qval<0 | qval>1) ) message("At least one computed frequency was found to be <0 or >1. The frequencies will be adjusted/normalized accordingly. Zero frequencies can be expected.")
  qval <- pmin(pmax(qval,0),1)
  qval <- qval/sum(qval)

  res <- cbind(GENO-1, q=qval)
  res <- data.frame(res)
  cols <- if( nloci==2 ){
    c("Z","G")
  }else c("Z",paste0("G",1:(nloci-1)))
  colnames(res)[1:nloci] <- cols
  return(res)
}






