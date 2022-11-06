
#' Performs inference on two-phase studies data via semiparametric maximum likelihood.
#'
#'
#' @param formula  regression formula, note that if it does not contain the missing-by-design variable, miscov, it will return results under the null hypothesis. Hypothesis testing corresponds to the score statistic. Otherwise, estimates and hypothesis testing ocurr under the alternative hypothesis leading to Wald statistics. All the elements in formula except miscov must be present in data0 and data1.
#' @param miscov right hand side formula with the missing-by-design covariate(s), i.e. the potential causal locus (loci). Must be present in data1 but absent in data0.
#' @param auxvar  right hand side formula with the auxiliary variable(s), i.e. the GWAS SNP from phase 1. Must be present in data0 and data1.
#' @param family member of the exponential family (see \code{\link{family}}, 'quasi' models not available). Default gaussian().
#' @param data0 a dataframe with the complement of the phase 2 data. Must contain the unique elements in formula and auxvar but NOT miscov.
#' @param data1 a dataframe with the phase 2 data. Must contain the unique elements in formula, auxvar, and miscov.
#' @param start.values a named list with initial values for the regression parameters and joint distribution between miscov and auxvar (only one can be specified). Defaults to NULL
#' @param verbose verbose output? logical, defaults to FALSE.
#'
#' @details these are some additional details
#'
#' @return  A list of objects
#'
#'
#' @examples
#'
#' data = DataGeneration_TPD()
#' set.seed(1)
#' R = rep(0, nrow(data)); R[sample(nrow(data),500)] <- 1 # random phase 2 subsample of 500.
#' data0 = data[R==0,c('Y','Z')]
#' data1 = data[R==1,c('Y','Z','G1')]
#' res_Ho = twoPhaseSPML(formula =  Y ~ Z,
#' miscov = ~ G1,
#' auxvar = ~ Z,
#' data0 = data0, data1 = data1)
#' res_Ha = twoPhaseSPML(formula =  Y ~ Z + G1,
#' miscov = ~ G1,
#' auxvar = ~ Z,
#' data0 = data0, data1 = data1)
#' @export
#' @import stats
#' @import data.table
#' @importFrom MASS ginv
#' @importFrom Matrix fac2sparse t crossprod colSums
twoPhaseSPML <- function(formula, miscov, auxvar, family=gaussian, data0, data1, start.values=NULL, verbose=FALSE)
{
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
  terms_data<-if( all(G_ %in% all.vars(formula)) ){ c(all.vars(formula)[-which(all.vars(formula)%in%G_)],Z_) }else c(all.vars(formula),Z_)
  n1<-nrow(data1); n0<-nrow(data0); n<-n1+n0
  nG_ = length(G_)
  uniqZ <- unique(data1[,Z_]); uniqZ<-uniqZ[order(uniqZ)]


  uniqG <- unique(data1[,G_])
  dat0 <- cbind(id_=rep(seq(n1+1,n,by=1), each=NROW(uniqG)),
                data0[rep(seq_len(n0), each=NROW(uniqG)), ],
                matrix(rep(t(uniqG),n0), ncol=NCOL(uniqG), byrow=TRUE))
  names(dat0)[(ncol(dat0)-NCOL(uniqG)+1):ncol(dat0)] <- G_
  qform <- as.formula(paste0("wg_ ~",paste(G_,collapse="+"),"+",Z_)) # may need to modify this if multiple Z's


  dat <- rbind(cbind(id_=1:n1,data1[,terms_dat]),dat0[,c("id_",terms_dat)]) # pseudo-data

  #initialize the non parametric terms, i.e. p_g,z, the goint distribution between miscov and auxvar
  uniqGZ<-unique(data1[,c(G_,Z_)])
  q0<-cbind(uniqGZ,q=1/nrow(uniqGZ))
  Zuniq<-unique(c(data1[,Z_],data0[,Z_]))
  q<-cbind(matrix(rep(t(uniqG),length(Zuniq)), ncol=NCOL(uniqG), byrow=TRUE),Z=rep(Zuniq,each=NROW(uniqG))); colnames(q)<-c(G_,Z_)
  q<-merge(q,q0,all.x=T,sort=F)
  q[is.na(q$q),"q"]<-0
  nq<-nrow(q)
  if(!is.null(start.values$q)){
    q1<-merge(q,start.values$q,by=c(G_,Z_),all.x=T,sort=F)
    q1<-q1[,c(G_,Z_,"q.y")]; names(q1)[3]<-"q"
    q<-q1
  }

  #initialize the theta parameters, i.e. regression parameters (with variance term in the gaussian case)
  X <- model.matrix(formula,dat) #We'll need this model matrix for the fitting

  if( is.null(start.values$betas) ){
    theta.o <- glm(formula, data = dat, family = family)$coef
    } else theta.o <- start.values$betas

  if(family$family %in% c("gaussian","Gamma")){
    disp=var(c(data1[,Y_],data0[,Y_]))
  }else disp=1

  #some needed values
  nparms<-length(theta.o)
  ng<-NROW(uniqG)
  n1ones <- rep(1,n1)

  #initialize the parameters for the EM stopping rules (may consider introducing a list with control options for flexibility.)
  maxiter=1000
  tol=1e-05
  iter<-0
  theta<-theta.o; theta.o<-theta+1;
  q.o=q; q.o[,"q"]=q.o[,"q"]-1;
  if(verbose) cat("Initial iteration:",iter,"theta =",theta,"pGZ =",as.numeric(q[,"q"]),"\n")

  while( (iter<=maxiter & (max(abs(theta-theta.o))>tol | max(abs(q[,"q"]-q.o[,"q"]))>tol)) ){

    theta.o<-theta; q.o<-q; iter<-iter+1

    #E-step
    dat$wg_ <- c(n1ones, fnc_wgts(formula,family,dat0,theta,q,disp,G_,Z_,uniqZ))

    #M-step
    glmfit <- suppressWarnings( irls_svdnewton(X, dat[,Y_], family=family,
                                               weights=dat$wg_, rank_deficiency="minimum norm") )
    theta <- glmfit$coef

    if( family$family %in% c("gaussian","Gamma") ){
      disp = sum(dat$wg_*glmfit$residuals^2)/n ##sigma^2
    }
    q.it <- data.frame(xtabs(qform,data=dat)/n)

    ordq <- match(apply(q[,c(G_,Z_)],1,paste,collapse="-"),
                  apply(q.it[,c(G_,Z_)],1,paste,collapse="-"))

    q[,"q"] <- q.it[ordq,"Freq"]
    if(verbose) cat("Iteration:",iter,"theta =",theta,"pGZ =",as.numeric(q[,"q"]),"\n")
  }

  if(verbose) cat("Total iterations:",iter,"theta_end=",theta,"pGZ_end=",as.numeric(q[,"q"]),"\n")

  #Create objects needed in the testing part
  qG <- aggregate(as.formula(paste0("q~",paste(G_,collapse="+"))),q,FUN=sum)
  qG$k_ <- 1:nrow(qG)
  dat$ord_ <- 1:nrow(dat)
  dat <- merge(dat,qG[,c(G_,"k_")],by=G_)
  dat <- dat[order(dat$ord_),]

  if( Ho<-!all(G_ %in% names(theta)) ){ ## under the null (no G in formula)
    indx_thetas_H0 <- which(!G_ %in% names(theta))
    nG_ <- length(indx_thetas_H0)
    val <- c(theta, rep(0,nG_))
    id <- c(seq_along(theta), 1+seq(0.1,0.9,length.out=nG_))
    theta1 <- val[order(id)]
    names(theta1)[which(names(theta1)=="")] <- G_[indx_thetas_H0]
    formula1 <- update(formula,paste0(".~",paste(G_[indx_thetas_H0],collapse="+"),"+ ."))
    X1 <- model.matrix(formula1,dat)
    eta1 <- X1%*%theta1
    W <- family$mu.eta(eta1)^2/family$variance(family$linkinv(eta1))
    glmfit1 <- list(); glmfit1$qr <- qr(X1*sqrt(as.numeric(dat$wg_*W)))

    FIM = FIM.glm(theta1,qG,disp,formula1,Y_,dat,family,glmfit1)
    expFIM = expFIM.glm(theta1,qG,disp,formula1,Y_,dat,family,glmfit1)

    res_SobsIM <- Testing_EM_joint(theta1,qG,NULL,G_[indx_thetas_H0],FIM)
    res_SexpIM <- Testing_EM_joint(theta1,qG,NULL,G_[indx_thetas_H0],expFIM)

  }else{ ## under the aternative (G in formula)
    FIM = FIM.glm(theta,qG,disp,formula,Y_,dat,family,glmfit)
    expFIM = expFIM.glm(theta,qG,disp,formula,Y_,dat,family,glmfit)

    res_VobsIM <- Testing_EM_joint(theta,qG,NULL,G_,FIM)
    res_VexpIM <- Testing_EM_joint(theta,qG,NULL,G_,expFIM)
  }

  ll <- .loglik(theta,q,disp,formula,Y_,G_,Z_,dat,family) ## log-likelihood

  if(Ho){
    return(list(theta=theta,qGZ=q,dispersion=disp,Sobs=res_SobsIM$S,Sexp=res_SexpIM$S,df=res_SexpIM$df,qG=qG,iter=iter,ll=ll))
  }else{
    return(list(theta=theta,qGZ=q,dispersion=disp,var_theta=res_VobsIM$Var,Wobs=res_VobsIM$W,Wexp=res_VexpIM$W,df=res_VexpIM$df,qG=qG,iter=iter,ll=ll))
  }
}

