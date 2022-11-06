#' function optimTP.LM
#'
#'
#' @param formula formula1name
#' @param miscov xxxxxxxxcode
#' @param auxvar auxvar1
#' @param strata strata1
#' @param family family1
#' @param n n1
#' @param data data1
#' @param beta beta1
#' @param p_gz p_gz1
#' @param disp disp1
#' @param optimMeasure optimMeasure1
#' @param K.idx K.idx1
#' @param min.nk min.nk1
#' @param logical.sub logical.sub1
#'
#' @details details at here
#'
#' @examples
#'
#' print(1)
#'
#' @export
#' @importFrom nloptr cobyla
#' @importFrom dfoptim hjkb
#' @importFrom utils combn stack tail
#' @importFrom enrichwith enrich
optimTP.LM  <-  function(formula,miscov,auxvar,strata,family,n,data,beta,p_gz,disp=NULL,optimMeasure,K.idx=NULL,min.nk=NULL,logical.sub=NULL){

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  ### TODO: check if all unique values in data[,Z_] are in p_gz

  terms_dat <- unique( c(all.vars(formula), attr(terms(miscov), "term.labels") , attr(terms(auxvar), "term.labels"), attr(terms(strata), "term.labels") , attr(terms(auxvar), "term.labels")) )

  Y_ <- as.character(formula[[2]])
  G_ <- all.vars(miscov)
  Z_ <- all.vars(auxvar)
  S_ <- all.vars(strata)

  if( is.null(disp) ){
    if( family$family %in% c("gaussian","Gamma") ){
      disp = var(data[,Y_])
    }else disp=1
  }

  ## remove enties with zero or negative probabilities
  p_gz <- p_gz[p_gz$q>0,]

  N=NROW(data)
  rho=n/N
  if( rho<0.1 )warning("The phase 2 sample size (n2) is smaller than recommended. Please proceed with caution.")

  stratadf <- data.frame(xtabs(strata,data=data))
  #Xs <- model.matrix(strata, stratadf)
  stratadf[S_] <- lapply(stratadf[S_], function(x)as.numeric(as.character(x))) ## may need it for data.table
  g.vals <- unique(p_gz[,G_])
  z.vals <- unique(p_gz[,Z_])

  if( !is.null(logical.sub) ){
    data <- data[logical.sub,]
    N <- NROW(data)
  }

  dat_ <- cbind(id_=rep(seq(1,N,by=1), each=NROW(g.vals)),data[rep(seq_len(N), each=NROW(g.vals)), ],matrix(rep(t(g.vals),N), ncol=NCOL(g.vals), byrow=TRUE))
  names(dat_)[(ncol(dat_)-NCOL(g.vals)+1):ncol(dat_)] <- G_
  dat_ <- droplevels(dat_)
  dat_$ord_ <- 1:nrow(dat_)

  wgs <- omega1_g(formula,family,dat_,beta,p_gz,disp,G_,Z_) # P(G|Y,Z) and P(Y,Z)

  dat_$wg_ <- wgs$wg
  dat_ <- dat_[!is.na(dat_$wg_),]

  M1 <- obsIMf1(formula,family,dat_,beta,p_gz,disp,G_,Z_)$I3
  M1 <- M1/N

  IM3_subs <- do.call(rbind,by(dat_, INDICES = dat_[,S_],function(df.by){
    Nk <- N#length(unique(df.by$id_))
    IM <- obsIMf1(formula,family,df.by,beta,p_gz,disp,G_,Z_)
    data.frame(df.by[1,S_,drop=FALSE],t(as.vector(IM$I3)/Nk))
  }))


  IM2_subs <- do.call(rbind,by(dat_, INDICES = dat_[,S_],function(df.by){
    Nk <- N#length(unique(df.by$id_))
    IM <- obsIMf1(formula,family,df.by,beta,p_gz,disp,G_,Z_)
    data.frame(df.by[1,S_,drop=FALSE],t(as.vector(IM$I2)/Nk))
  }))

  lenS_ <- length(S_)
  npars <- NCOL(M1)

  ### Indicator for the beta parameters of interest
  Ind <- rep(FALSE,npars-1)
  Betas_ids <- which(grepl(paste0(paste0("^",G_),collapse="|"),colnames(model.matrix(formula,dat_[1,]))))
  Ind[Betas_ids] <- TRUE

  if( optimMeasure=="Par-spec" & is.null(K.idx) ) stop("For a parameter-specific criterion K.idx must be provided.")
  ## asign these functions to the execution environment
  optMsuref=.optMsure(optimMeasure);  environment(optMsuref)=environment()

  ### parameters for the constained optimization
  # if( !all.equal(pS_[S_],stratadf[,S_],attributes=FALSE,check.attributes = FALSE) ) stop("Strata groups in pS and stratadf don't match!")
  # ais <- (stratadf$Freq)#/sum(stratadf$Freq)
  # ais <- N*pS_$pS_#/sum(stratadf$Freq)
  npis <- nrow(stratadf)-1
  rho <- n
  upperLM <- as.numeric(stratadf$Freq) #rep(1,npis+1)
  # M1inv <- MASS::ginv(M1)
  # In <- diag(npars)
  D <- rbind(diag(npars-1),0)
  # Dinv <- MASS::ginv(D)
  # Dtinv  <- MASS::ginv((t(D)))
  ### result on how to calcuate the inverse of A+B
  # https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
  # (A+B)^{-1} = A^{-1}(I-(I+BA^{-1})^{-1}BA^{-1})


  ### The LM approach is going to get optimal values for nk/P(R=1,Z,K) as opposed to pr(R=1|Y,Z). In additio the  objective function is determined according to whether the design regression parameters are under the null or alternative

  if( all(beta[Ind]==0) ){
    ObjFun <- function(pis){ # pis = (upperLM*n/N)[-1]
      pis1 <- c(0,pis)
      pis1[1] <- (rho-sum(pis1))
      stratadf$prReq1YZ <- (pis1/N)/stratadf$Freq

      IM2b_subs <- as.data.frame(merge(as.data.table(IM2_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))
      IM3b_subs <- as.data.frame(merge(as.data.table(IM3_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))

      FIM <- M1 + matrix(colSums(IM2b_subs[,(lenS_+1):(npars^2+lenS_)]*IM2b_subs$prReq1YZ),ncol=npars) - matrix(colSums(IM3b_subs[,(lenS_+1):(npars^2+lenS_)]*IM3b_subs$prReq1YZ),ncol=npars)

      ## to deal with the constraint on the p_g's
      F1 <-  crossprod(D, FIM %*% D)

      ###  Score test variance
      invF1_noInd <- ginv(F1[!Ind,!Ind])
      V1 <- F1[Ind,Ind] - ( F1[Ind,!Ind] %*% invF1_noInd %*% F1[!Ind,Ind] )
      invV1 <- solve( V1 )

      return( optMsuref(invV1,K.idx) )
    }
  }else {
    ObjFun <- function(pis){ # pis = (upperLM*n/N)[-1]
      pis1 <- c(0,pis)
      pis1[1] <- (rho-sum(pis1))
      stratadf$prReq1YZ <- (pis1/N)/stratadf$Freq

      IM2b_subs <- as.data.frame(merge(as.data.table(IM2_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))
      IM3b_subs <- as.data.frame(merge(as.data.table(IM3_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))

      FIM <- M1 + matrix(colSums(IM2b_subs[,(lenS_+1):(npars^2+lenS_)]*IM2b_subs$prReq1YZ),ncol=npars) - matrix(colSums(IM3b_subs[,(lenS_+1):(npars^2+lenS_)]*IM3b_subs$prReq1YZ),ncol=npars)

      ## to deal with the constraint on the p_g's
      F1 <-  crossprod(D, FIM %*% D)

      # B <- matrix(colSums(IM2b_subs[,(lenS_+1):(npars^2+lenS_)]*IM2b_subs$prReq1YZ),ncol=npars) - matrix(colSums(IM3b_subs[,(lenS_+1):(npars^2+lenS_)]*IM3b_subs$prReq1YZ),ncol=npars)
      #
      #
      # FIMinv <- M1inv %*% (In - solve(In+B%*%M1inv)) %*% (B%*%M1inv)
      #
      # vcov = Dinv %*% FIMinv %*% Dtinv

      vcov = solve(F1)

      return( optMsuref(vcov,K.idx) )

    }
  }


  if( is.null(min.nk) ){
    minLB.val <- 0 #ifelse(rho<0.05,rho*0.01,0.05)
    minLB <- rep(minLB.val,npis+1)
  }else{
    if( !length(min.nk) %in% c(1,npis+1) ) stop("The number of elements of min.nk is either 1 or ",npis+1)
    if( any(min.nk<0) )stop("The minimum stratum size must greater or equal than zero")
    if( any(min.nk>stratadf$Freq) )stop("The minimum stratum size must not be greater than the strata sample sizes determined by strata")
    minLB <- min.nk
  }
  initX0 <- (upperLM*n/N)[-1] #rep(rho,npis)
  initX1 <- (upperLM*n/N)  #rep(rho,npis+1)

  ObjFun_in <- function(pis){
    ui <- rbind(rep(-1,npis),rep(1,npis))
    ci <- c(-rho,0)
    return(as.numeric(ui%*%pis-ci))
  }

  ObjFun_eq <- function(pis1){
    return(sum(pis1)-rho)
  }

  ### Mind the warning: "For consistency with the rest of the package the inequality sign may be switched from >= to <= in a future nloptr version."

  ## To disable this warning and avoid printing do the following. However, new versions of the package may return an error when  this is enforced.
  options(nloptr.show.inequality.warning = FALSE)

  ### Note that the function can be changed when the change is implemented using
  # hin <- function(x) (-1)*f2(x, ...)  # NLOPT expects hin <= 0

  ### use this one for base, if an error then use the rest
  sol <- cobyla(initX0, fn = ObjFun, lower = minLB[-1], upper = upperLM[-1], hin = ObjFun_in, nl.info = FALSE, control = list(xtol_rel = 1e-6, maxeval = 10000))

  # return the option to the original state
  options(nloptr.show.inequality.warning = TRUE)

  if( !(sol$convergence>0 & sol$convergence<5) ){
    if( sol$convergence==5 ) message("The maximum number of evaluations in nloptr:::cobyla() has been reached. Trying dfoptim:::hjkb() now. You could also increase the number of evaluations via maxeval (although it is fairly large already.)")
    if( sol$convergence<0 ) message(paste0("Convergence in nloptr:::cobyla() was not achieved; error code ",sol$convergence,". Trying dfoptim:::hjkb() now."))

    ObjFunb <- function(pisA, k){ # pisA <- (upperLM*n/N)
      nstar <- sum(pisA)
      return( ObjFun(pisA[-1]) + k*(nstar-rho)^2 )
    }

    sol1 <- hjkb(initX1, fn=ObjFunb, lower=minLB, upper=upperLM, k=0.01, control=list(maxfeval=10000))
    sol1 <- tryCatch(hjkb(sol1$par, fn=ObjFunb, lower=minLB, upper=upperLM, k=100, control=list(maxfeval=10000)), error=function(e){ NULL })
    if( is.null(sol1) ){
      solB <- c(0,sol$par)
      solB[1] <- (rho-sum(solB))

      sol1 <- hjkb(solB, fn=ObjFunb, lower=minLB, upper=upperLM, k=100, control=list(maxfeval=10000))
    }

    if( sol1$convergence!=0 ){
      warning(paste0("Convergence in dfoptim:::hjkb() was not achieved. Error code ",sol1$convergence,". Returning last par value."))
    }

    stratadf$convergence <- c(sol1$convergence,rep(NA,nrow(stratadf)-1))
    pisopt <- sol1$par

  } else {
    stratadf$convergence <- c(sol$convergence,rep(NA,nrow(stratadf)-1))
    pisopt <- c(0,sol$par)
    pisopt[1] <- (rho-sum(pisopt))
  }

  stratadf$prR_cond_optim <- pisopt/sum(pisopt)

  return(stratadf)
}


#' function optimTP.GA
#'
#' @param ncores ncores1
#' @param formula the formula
#' @param miscov miscov1
#' @param auxvar auxvar1
#' @param family family1
#' @param n n1
#' @param data the data gggg
#' @param beta beta1
#' @param p_gz p_gz1
#' @param disp disp1
#' @param ga.popsize ga.popsize1
#' @param ga.propelit ga.propelit1
#' @param ga.proptourney ga.proptourney1
#' @param ga.ngen ga.ngen1
#' @param ga.mutrate ga.mutrate1
#' @param ga.initpop ga.initpop1
#' @param optimMeasure optimMeasure1
#' @param K.idx K.idx1
#' @param seed seed1
#' @param verbose verbose1
#'
#' @details details here
#'
#' @examples
#'
#' print(1)
#'
#' @export
optimTP.GA  <-  function(ncores,formula,miscov,auxvar,family,n,data,beta,p_gz,disp=NULL,ga.popsize,ga.propelit,ga.proptourney,ga.ngen,ga.mutrate,ga.initpop=NULL,optimMeasure,K.idx=NULL,seed=1,verbose=0){
  ptm0 <- Sys.time()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  terms_dat<-unique( c(all.vars(formula), attr(terms(miscov), "term.labels") , attr(terms(auxvar), "term.labels")) )
  Y_<-as.character(formula[[2]]); G_<-all.vars(miscov); Z_<-all.vars(auxvar)

  ## remove entries with zero or negative probabilities
  p_gz <- p_gz[p_gz$q>0,]

  ### check if names in G_ and Z_ are in the joint ditribution data.frame p_gz
  g.vals<-unique(p_gz[,G_]);
  N <- NROW(data)

  if( is.null(disp) ){
    if( family$family %in% c("gaussian","Gamma") ){
      disp = var(data[,Y_])
    }else disp=1
  }
  ### Function that determines the optimality measure (D-, A- or parameter-specfic for now)
  ### if optimMeasure=="Par-spec" then K.idx MUST be not null
  if( optimMeasure=="Par-spec" & is.null(K.idx) ) stop("For a parameter-specific criterion K.idx must be provided.")

  ## assign these functions to the execution environment
  omega1_g=omega1_g; optMsuref=.optMsure(optimMeasure)#; class.ind2=class.ind2; class.ind2.Mat=class.ind2.Mat

  # envup=environment(omega1_g)=environment(class.ind2)=environment(class.ind2.Mat)=environment(optMsuref)=environment()
  envup=environment(omega1_g)=environment(optMsuref)=environment()

  uniqZ<-unique(p_gz[,Z_]); uniqZ<-uniqZ[order(uniqZ)]

  ### calculate FIMs for all samples (hopefully, a faster way than before)
  dat_ <- cbind(id_=rep(seq(1,N,by=1), each=NROW(g.vals)),data[rep(seq_len(N), each=NROW(g.vals)), ],matrix(rep(t(g.vals),N), ncol=NCOL(g.vals), byrow=TRUE))
  names(dat_)[(ncol(dat_)-NCOL(g.vals)+1):ncol(dat_)] <- G_
  dat_ <- droplevels(dat_)

  wgs <- omega1_g(formula,family,dat_,beta,p_gz,disp,G_,Z_)

  dat_$wg_ <- wgs$wg
  dat_ <- dat_[!is.na(dat_$wg_),]

  IM1 <- obsIMf1(formula,family,dat_,beta,p_gz,disp,G_,Z_)$I3/N

  obsIM23_id <- obsIMf1(formula,family,dat_,beta,p_gz,disp,G_,Z_,by.id = TRUE)

  npars <- NCOL(IM1)

  IM2_id = obsIM23_id$I2/N
  IM3_id = obsIM23_id$I3/N

  Ind <- rep(FALSE,npars-1)
  Betas_ids <- which(grepl(paste0(paste0("^",G_),collapse="|"),colnames(model.matrix(formula,dat_[1,]))))
  Ind[Betas_ids] <- TRUE
  D <- rbind(diag(npars-1),0)

  if( all(beta[Ind]==0) ){
    Assess_fitness <- function(R){

      FIM <- IM1 + matrix(colSums(IM2_id*R),ncol=npars) - matrix(colSums(IM3_id*R),ncol=npars)
      F1 <-  t(D) %*% FIM %*% D

      invF1_noInd <- ginv(F1[!Ind,!Ind])
      V1 <- F1[Ind,Ind] - F1[Ind,!Ind] %*% invF1_noInd %*% F1[!Ind,Ind]
      invV1 <- solve( V1 )

      return( optMsuref(invV1,K.idx) )
    }


  } else {
    Assess_fitness <- function(R){

      FIM <- IM1 + matrix(colSums(IM2_id*R),ncol=npars) - matrix(colSums(IM3_id*R),ncol=npars)
      F1 <-  t(D) %*% FIM %*% D

      vcov = solve(F1)
      return(optMsuref(vcov,K.idx))

    }
  }

  #### R<-rep(0,N); R[sample(N,n)] <- 1
  R0s <- rep(0,N)
  ObjFun <- function(v){
    R <- R0s; R[v]<-1
    Assess_fitness(R)
  }

  time.elapsed0 <- as.numeric(difftime(Sys.time(), ptm0, units="secs"))

  ptm <- Sys.time()

  out <- kofnGA.mc(verbose=verbose, ncores, n = N, k = n, OF = ObjFun, popsize = ga.popsize,
                   keepbest = ceiling(ga.popsize*ga.propelit),
                   tourneysize = max(ceiling(ga.popsize*ga.proptourney), 2),
                   ngen = ga.ngen, mutprob = ga.mutrate, initpop = ga.initpop,
                   varsExport=c("Assess_fitness","optMsuref","K.idx","IM1","IM2_id","IM3_id","D","R0s"),envup=envup)

  time.elapsed <- as.numeric(difftime(Sys.time(), ptm, units="secs"))
  out$time.elapsed <- c(preproc=time.elapsed0,GA=time.elapsed)

  return(out)
}





