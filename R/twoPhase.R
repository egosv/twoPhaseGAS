#' Internal function to optimize over a range of maf, LD and betas.
#'
#'
#' @param design_formula Formula for the regression model, default is Y ~ G + Z. Z is the GWAS SNP, G is the sequence variant. Y is outcome.
#' @param beta  Vector of betas (length of 3, corresponding to intercept, effect size for G, effect size for Z).
#' @param maf_G Maf for G (numeric, ranging from 0 to 1).
#' @param LD Correlation r betweek G and Z (numeric, ranging from -1 to 1, r value).
#' @param data Data frame with Y and Z variables.
#' @param n2 Phase 2 sample size
#' @param family Distrubution of the outcome (default: gaussian()  ). Familes available for glm can be used here. See help(stats::family) for examples.
#' @param useGeneticAlgorithm   If TRUE, use genetic algorithm in addition to Lagrange multiplier approach (slower). Default: FALSE
#'
#'
#'
#' @return  sth
#'
#'
#' @examples
#'
#' data = twoPhase()
#'
#' @export
twoPhase <- function(
  beta = c( 1,1,1 ),
  maf_G = 0.1,
  LD = 0.3 ,
  data = NA,
  n2 = NA,
  design_formula  =  Y ~  G  +  Z,
  family = gaussian(),
  useGeneticAlgorithm = FALSE
){

  cat("Optimizing for a range of maf and LD\n")

  # maf_G = c(0.2,0.4)
  #LD = c(0.1,0.3,0.5)


  if( class(beta) != "list") beta = list( beta )


  values = expand.grid(nbeta = 1:length(beta), maf_G = maf_G, LD = LD)


  L = nrow(values)


  res = list()
  nres = 1


  for(l in 1:L ) {

      val_beta = beta[[values$nbeta[l]]]
      val_maf_G = values$maf_G[l]
      val_LD = values$LD[l]

      cat("  - running beta=", paste(val_beta,sep=",",collapse=",") , " maf_G=", val_maf_G, " LD=", val_LD,"\n", sep="")

      result <- try( twoPhaseDesign( beta = val_beta , maf_G = val_maf_G, LD = val_LD, data = data, n2 = n2,  design_formula = design_formula,
                                   family = family, design = ifelse(useGeneticAlgorithm,"optGA","optLM") ) )


      res[[nres]] = result
      nres = nres+1


      names(res)[l] = sprintf("beta=%s maf_G=%.3f LD=%.3f", paste(val_beta,collapse=","), val_maf_G , val_LD )

  }


  # res = res[ sapply(res,class) != "try-error" ] ## better do this filtering outside the function to identify which combinations broke

  return( res )
}




#' Compute sample allocations for a two phase study design.
#'
#'
#' @param design_formula Formula for the regression model, default is Y ~ G + Z, where Z is the GWAS SNP, G is the sequence variant. Y is outcome. Rename the variables in your data.frame to match Y and Z, G is the seq-SNP not present in the data.frame.
#' @param beta  Vector of betas (length of 3, corresponding to intercept, effect size for G, effect size for Z).
#' @param maf_G minor allele frequency for G (numeric, ranging from 0 to 1). Numeric or vector of possible values.
#' @param LD Correlation r betweek G and Z (numeric, ranging from -1 to 1, r value). Numeric or vector of possible values.
#' @param data Data frame with Y and Z variables.
#' @param n2 Phase 2 sample size.
#' @param family Distrubution of the outcome. Default: gaussian(). Familes available for glm can be used here. See help(stats::family) for examples.
#' @param S stratification of the outcome Y. Optional. Should be a numeric vector with strata categories (e.g. 1 1 2 2 3 3 ). If present, its length must be equal to the number of rows in data. Default: NULL. Needed when Y does not render itself into strata, e.g. Gaussian, Poisson, Gamma.
#' @param perc_Y vector of percentiles in increasing order for which the outcome Y will be stratified. Default: c(1/5,4/5). Only used when S is NA. Note that setting up S is strongly suggested.
#' @param p_gz data frame with the joint distribution between G and Z. See examples for the right format. Default: NULL. If present, values of maf_G and LD are disregarded in the analysis.
#' @param design string for the design to use for phase 2 sample selection. One of residual-dependent sampling ("RDS"), optimal as defined by Tao, Zheng and Lin (2019) ("TZL"), optimal via Lagrange multipliers ("LM"), optimal via genetic algorithm ("GA"), probability proportional to size ("PPS"), balanced ("BAL") or combined ("COM") allocations. Default: "RDS". See details for a more explanations.
#' @param ndraws integer that determines the number of draws to examine when design is one of "pps", "bal" or "comb" and the design parameter combinations is greater than 1. Default: 10
#' @param optimCriterion string denoting the optimality criterion used during the optimization. One of "Par-spec" (default), "A-opt" or "D-opt". For parameter-specific, A-optimality or D-optimality, respectively.
#' @param overallMethod string denoting the method to select the overall design when multiple design parameters are given. One of "med-max" (default) or "cumm" for median-maximum and cummulative frequencies, respectively.
#'
#' Note that in this version strata are always defined in terms of Z and S, i.e. a joint design, future implementations may relax this by allowing for only S or only Z (marginal designs, outcome- or covariate-dependent, respectively)
#' @return  sth
#'
#'
#' @export
twoPhaseDesign <- function(

      beta,
      maf_G,
      LD,
      data,
      n2,

      design_formula  =  Y ~  G  +  Z,
      family = gaussian(),
      S,
      perc_Y =  c(1/5,4/5),
      p_gz,

      #useGeneticAlgorithm = FALSE,
      design = c("RDS","LM","GA","PPS","BAL","COM","TZL"),
      ndraws = 10,
      optimCriterion = c("Par-spec","A-opt","D-opt"),
      overallMethod = c("med-max","cumm")

)
{

  #### Perform some checkings

  if( missing(beta) | missing(maf_G) | missing(LD) | missing(data) | missing(n2)) stop("none of beta, maf_G, LD, data or n2 can be missing.")

  if( class(data) != "data.frame" ) stop("data should be a data.frame")

  terms1 = all.vars(design_formula)
  terms2 = setdiff(terms1,c("G","Z","Y") )

  if( ! ( "G" %in% terms1 ) ) { stop("G should be included in the design_formula") }
  if( ! ( "Z" %in% terms1 ) ) { stop("Z should be included in the design_formula") }
  if( ! ( "Y" %in% terms1 ) ) { stop("Y should be included in the design_formula") }


  if( sum(terms2 %in% colnames(data) ) != length(terms2)  ) {

    s = which( ! (terms2 %in% colnames(data) ) ) [1]

    stop(sprintf("The variable %s is missing from the data.frame data", terms2[s]) )

  }


  data_tmp  = data[1:2,]; if(  !( "G" %in% colnames(data_tmp) ) )  data_tmp$G = 0.1;
  beta_names = colnames( model.matrix(design_formula,data_tmp) )

  if ( class(beta) == "list" ) {

    tmp = sapply(beta,length)
    beta_test = beta[[1]]

    if( min(tmp) != max(tmp) ) stop("Length of betas is not equal for all elements of list\n");

  } else  beta_test = beta

  if( length(beta_names) != length(beta_test) ) {

    stop( sprintf("Number of terms in design formula and length of beta do not match: %d terms in beta, %d terms in formula\n",
                  length(beta_test), length(beta_names) ) )

  }


  N = nrow(data) ## phase 1 sample size

  Z = data$Z

  maf_Z = mean(Z, na.rm=T)/2

  Y = data$Y

  if ( missing(S) ){

    if( !is.numeric(perc_Y) || !all(perc_Y<1,perc_Y>0) || is.unsorted(perc_Y) ) stop("perc_Y must all be numeric values between 0 and 1 and in increasing order." )
    S = as.numeric( cut(Y, breaks=quantile(Y, probs=c(0,perc_Y,1), na.rm=TRUE),
            include.lowest=TRUE) )

  } else {
    S_uniq = unique(S)
    S_uniq_sort = S_uniq[order(S_uniq)]
    if( !is.numeric(S_uniq_sort) | any(diff(S_uniq_sort)!=1) | S_uniq_sort[1]!=1 ) stop("S must contain only consecutive integers starting from 1")
  }

  data$S = S

  design <- match.arg(design)

  # These parameters are used in the functions
  strataformula <- ~ Z+S ## formula that determines the stratification factors, i.e. how the strata is defined (other formulas may be used, for example strataformula = ~S or strataformmula= ~Z). These must be discrete numeric variables in the phase 1 data

  ## optimality criterion (so far only 3 implemented)
  optMethod <- match.arg(optimCriterion)
  overallMethod <- match.arg(overallMethod)

  if( optMethod=="Par-spec" ){
    Kind <- which(beta_names=="G")
  } else Kind <- NULL

  if( !missing(p_gz) ){
    message("Since 'p_gz' is specified, values of maf_G and LD are disregarded from the analysis.")
    maf_G <- NULL
    LD <- NULL
  }

  # Optimize over a series of parameters

  if( length(maf_G) > 1 || length(LD) > 1 || class( beta ) == "list" ) {

    if( design %in% c("optLM","optGA") ){
      result = twoPhase( beta = beta ,
                                 maf_G = maf_G,
                                 LD = LD,
                                 data = data,
                                 n2 = n2,
                                 design_formula = design_formula,
                                 family = family,
                                 useGeneticAlgorithm = design=="optGA" )

      indx_noerror <- sapply(result,class) != "try-error"
      result = result[ indx_noerror ]
      nams = names(result)
      L = sum(indx_noerror)
      result_tmp = matrix( NA, nrow = N, ncol = L )
      for( l in 1:L ) {
        result_tmp[,l] = if( design=="optGA" ){
                      result[[l]]$R.opt.GA
                    }else result[[l]]$R.opt.LM
      }
      result = result_tmp
      colnames(result) = nams
      rm(result_tmp)

      if( class(beta) != "list" ) beta = list( beta )

      values = expand.grid(nbeta = 1:length(beta), maf_G = maf_G, LD = LD)

    } else{

      result = twoPhaseHeuristic(
        design = design,
        ndraws = ndraws,
        data = data,
        n2 = n2,
        family = family,
        S = data$S
      )

      L = ndraws


      if( class(beta) != "list" ) beta = list( beta )

      values = expand.grid(nbeta = 1:length(beta), maf_G = maf_G, LD = LD)

      indx_noerror = sapply(1:nrow(values),function(i){
        maf_G = values[i,"maf_G"]
        LD = values[i,"LD"]
        res = try(p_gz_func(maf_Z, maf_G, LD, "G", "Z"))
        if( class(res) == "try-error") cat("  - skipping this MAFs-LD combination: maf_Z=", maf_Z, " maf_G=", maf_G, " LD=", LD,"\n", sep="")

        class(res) != "try-error"
      })
    }

    values = values[ indx_noerror, ,drop=FALSE]
    L2 = nrow(values)

    # To select a single sample across all design values, we propose two methods (but there may be more)

    ### method 1 calculate the fitness and perform med-max
    if ( overallMethod == "med-max" ){

      FitnessVals <- matrix( NA, nrow = L, ncol = L2 )

      for( k in 1:L2 ){ # k <- 2
        val_beta_b = beta[[values$nbeta[k]]]
        val_maf_G_b = values$maf_G[k]
        val_LD_b = values$LD[k]
        p_gz0_b = p_gz_func(maf_Z, val_maf_G_b, val_LD_b, "G", "Z")

        obsIM.Rk <- obsIM.R(formula = design_formula,
                            miscov = ~G,
                            auxvar = ~Z,
                            family= family,
                            data = data,
                            beta = val_beta_b,
                            p_gz = p_gz0_b,
                            disp=NULL)

        for( l in 1:L ) {
          Rl = result[,l]
          FitnessVals[l,k] <- fitnessTP(obsIM.R = obsIM.Rk,
                                        Rj = Rl,
                                        optimMeasure = optMethod, K.idx = Kind)
        }
      }

      x <- apply(FitnessVals,1,median)
      x1 <- which(x==min(x))
      if( length(x1)>1 ) x1 <- x1[1]
      indx <- x1

      R.overall <- result[,indx]

    ### method 2, group the indicators to inform a new sampling
    } else if ( overallMethod == "cumm" ){

      x <- NA
      Freqs <- rowSums(result)
      indx <- sample(1:N,n2,prob=Freqs)
      R.overall <- rep(0, N)
      R.overall[indx] <- 1

    } else stop("Method for overall design not identified. Select one of med-max or cumm.")

    ### TODO: implement the sequential approach discussed with Radu on July 8, 2019

    data$R = R.overall
    R.name = switch(design,
                    optLM="R.opt.LM",
                    optGA="R.opt.GA",
                    pps="R.pps",
                    bal="R.bal",
                    comb="R.comb")
    names(data)[which(names(data)=="R")] = R.name

    return(list(overall_design=data,indiv_design=result,medianFitness=x))
  }


  # example data
  if(0) {

    maf_G = 0.1
    LD = 0.3
    data = data.frame(  Y = rnorm(1000),  Z = rbinom(1000, 2, 0.2) )
    beta = c( 1,1,1 )
    n2 = 200
    family = gaussian
    design_formula =  Y ~  G  +  Z
  }

  #** Single draw from heuristic design ----
  if( design %in% c("pps","bal","comb")){

    data.ret = twoPhaseHeuristic(design = design,
                                 ndraws = 1,
                                 data = data,
                                 n2 = n2,
                                 family = family,
                                 S = data$S)

    data.ret = cbind(data,R=data.ret)
    R.name = switch(design,pps="R.pps",bal="R.bal",comb="R.comb")
    names(data.ret)[which(names(data.ret)=="R")] = R.name
    return ( data.ret )
  }
  # TODO: Check data for validity and data types.


  ## determine "design" quantities to select the phase 2 data

  p_Z <- data.frame( xtabs(~Z,data)/nrow(data) )  ## calculates genotype distribution of Z from data

  LD.r <- LD ## specify LD between GWAS SNP, Z, and potentially "causal" locus

  ### if the joint distribution between g and z is provided, do not calculate it
  if( !missing(p_gz) ){

    p_gz0 = p_gz
  ### else calculate p_gz, i.e. the joint distribution of G and Z using maf_G, LD.r and available data from Z using function p_gz_func() in optimJTC_v1.R. This procedure sometimes fails, if that's the case we may cheat a bit and use the complete data.
  } else p_gz0 <-  p_gz_func(p_Z, maf_G, LD.r, "G", "Z")


  #** Phase 2 sample selection via optimal designs ----
  ## R1-R4 essentially contain a vector of size N with zeros (unselected) and ones (selected)

  ### check number of strata
  strat_df = data.frame(xtabs(strataformula,data=data))
  n_strat = nrow(strat_df)
  S_unique = unique(data$S)

  sel.prob.bal <- rep(1/n_strat,n_strat) ## selection probabilities for a balanced design

  ### To generalize the concept of a combined design for more than 3 groups for Y, we'll always set to selection probability for the median value(s) to 0

  ### Find the strata that will have zero probability of being selected (the median)
  if(  (length(S_unique)%%2) == 1  ){ ### calculate the median depending on whether the number of strata for S is even or odd
    med_stra <- median(unique(data$S))
  }else {
    med_stra <- as.numeric(quantile(S_unique,0.5,type=1))
    med_stra <- c(med_stra,med_stra+1)
  }
  indx_com <- which(!strat_df$S%in%med_stra)
  n_strat_com <- length(indx_com)


  sel.prob.com <- rep(0,n_strat) ## selection probabilities for a combined design
  sel.prob.com[indx_com] <- 1/n_strat_com


  # Basic stratified sampling for the "intuitive allocations"
  # R1 <- BSS(samp.fracs=sel.prob.bal,n2,data,strataformula)$R ## indicators for balanced design
  # R2 <- BSS(samp.fracs=sel.prob.com,n2,data,strataformula)$R ## indicators for combined design




  if(1)

    opt.prop <- optimTP.LM(  formula = design_formula,
                               miscov = ~ G,
                               auxvar = ~ Z,

                               strata = strataformula,

                               family = family,

                               n = n2,

                               data = data,

                               beta = beta,

                               p_gz = p_gz0,


                               disp = NULL,

                               optimMeasure  =  optMethod,

                               K.idx = Kind,
                               min.nk = NULL
                            )




  R3 <- BSS(samp.fracs = opt.prop$prR_cond_optim, n2,  data, strataformula) ## indicators for LM design



  ### TODO::::::::::::::  Add fitness of the LM results for testing



  if( design  == "optGA" )  {

  ### TODO:::::::::::::: Add more control to the GA parameters, e.g. whether initial sample is provided or not. If provided we should give means to let the user initialize thm

    ## For the genetic algorithm (GA), instead of providing a random initialization, I give an informed initialization based on combined, balanced and LM designs (a third each).


    pop1 <- sapply(1:20,function(x)  {


      sa <- which(BSS(samp.fracs=opt.prop$prR_cond_optim, n2, data, strataformula)==1)
      len <- length(sa)
      if( len==n2 ){ return(sa)
      }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
      }else return( c(sa,sample((1:N)[-sa],n2-len)) )


    })



    ## combined


    pop2 <- sapply(1:20,function(x)  {

      sa <- which(BSS(samp.fracs=sel.prob.com,n2,data,strataformula)$R==1)
      len <- length(sa)
      if( len==n2 ){ return(sa)
      }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
      }else return( c(sa,sample((1:N)[-sa],n2-len)) )

    })



    ## balanced


    pop3 <- sapply(1:20,function(x)  {

      sa <- which(BSS(samp.fracs=sel.prob.bal,n2,data,strataformula)$R==1)
      len <- length(sa)
      if( len==n2 ){ return(sa)
      }else if( len>n2 ){ return( sa[-sample(len, len-n2)] )
      }else return( c(sa,sample((1:N)[-sa],n2-len)) )


    } )


    ## Genetic Algorithm approach. This function may take a bit to run (in the order of a few minutes).
    GA.sol <- optimTP.GA(  ncores=1,
                            formula=design_formula,
                            miscov=~G,
                            auxvar=~Z,
                            family=family,
                            n=n2,
                            data,
                            beta=beta,
                            p_gz=p_gz0,
                            disp=NULL,

                ### TODO:::: ideally these parameters could be specified by the user
                            ga.popsize = 60,
                            ga.propelit = 0.8,
                            ga.proptourney = 0.8,
                            ga.ngen = 300,
                            ga.mutrate = 0.001,
                            ga.initpop = t(cbind(pop1,pop2,pop3)),

                            optimMeasure = optMethod,
                            K.idx = Kind,
                            seed = 1

              )


    R4 <- rep(0, N ); R4[GA.sol$bestsol] <- 1  ## indicators for GA design

    data$R.opt.GA = R4

  }



  if( design  == "optLM" ) {

     data$R.opt.LM = R3
    # print(opt.prop)
  }


  return ( data )


}

#' Select samples for phase 2 under heuristic designs.
#'
#' @param design Heuristic design to use for phase 2 sample selection. One of probability proportional to size ("pps"), balanced ("bal") or combined ("comb") allocations. Default: "pps".
#' @param ndraws Number of draws of the heuristic design to generate Default: 1.
#' @param data Data frame with Y and Z variables.
#' @param n2 Phase 2 sample size.
#' @param family Distrubution of the outcome. Default: gaussian(). Familes available for glm can be used here. See help(stats::family) for examples.
#' @param S stratification of the outcome Y. Optional. Should be a numeric vector with strata categories (e.g. 1 1 2 2 3 3 ). If present, its length must be equal to the number of rows in data. Default: NULL. Needed when Y does not render itself into strata, e.g. Gaussian, Poisson, Gamma.
#' @param perc_Y vector of percentiles in increasing order for which the outcome Y will be stratified. Default: c(1/5,4/5). Only used when S is NA. Note that setting up S is strongly suggested.
#'
#' Note that in this version strata are always defined in terms of Z and S, i.e. a joint design, future implementations may relax this by allowing for only S or only Z (marginal designs, outcome- or covariate-dependent, respectively)
#' @return  sth
#'
#'
#' @export
twoPhaseHeuristic <- function(
  design = c("pps","bal","comb"),
  ndraws = 1,
  data = NA,
  n2 = NA,

  family = gaussian(),
  S = NULL,
  perc_Y =  c(1/5,4/5)
)
{

  #### Perform some checkings

  if( class(data) != "data.frame" ) stop("data should be a data.frame")

  terms1 = all.vars(as.formula("Y~Z"))
  terms2 = setdiff(terms1,c("Z","Y") )

  if( ! ( "Z" %in% terms1 ) ) { stop("Z should be included in the design_formula") }
  if( ! ( "Y" %in% terms1 ) ) { stop("Y should be included in the design_formula") }


  if( sum(terms2 %in% colnames(data) ) != length(terms2)  ) {

    s = which( ! (terms2 %in% colnames(data) ) ) [1]

    stop(sprintf("The variable %s is missing from the data.frame data", terms2[s]) )

  }

  N = nrow(data) ## phase 1 sample size

  Y = data$Y

  if ( is.null(S) ){

    if( !all(perc_Y<1,perc_Y>0) || is.unsorted(perc_Y) ) stop("perc_Y must all be values between 0 and 1 and in increasing order." )
    S = as.numeric( cut(Y, breaks=quantile(Y, probs=c(0,perc_Y,1), na.rm=TRUE),
                        include.lowest=TRUE) )

  } else {
    S_uniq = unique(S)
    S_uniq_sort = S_uniq[order(S_uniq)]
    if( !is.numeric(S_uniq_sort) | any(diff(S_uniq_sort)!=1) | S_uniq_sort[1]!=1 ) stop("S must contain only consecutive integers starting from 1")
  }

  data$S = S

  design <- match.arg(design)

  ### TODO: check that S and Z are discrete


  # These parameters are used in the functions
  strataformula <- ~ Z+S ## formula that determines the stratification factors, i.e. how the strata is defined (other formulas may be used, for example strataformula = ~S or strataformmula= ~Z). These must be discrete numeric variables in the phase 1 data

  ### check number of strata
  strat_df = data.frame(xtabs(strataformula,data=data))
  n_strat = nrow(strat_df)
  S_unique = unique(data$S)

  ### draw the phase 2 sample indicators accroding to the selected design
  if( design == "pps" ){
    sel.prob.pps = rep(n2/N,n_strat) ## selection probabilities for a design with selection probabilities proportional to size

    Rs = sapply(1:ndraws,function(i)BSS(samp.fracs=sel.prob.pps,n2,data,strataformula) )

  } else if ( design == "bal" ){
    sel.prob.bal <- rep(1/n_strat,n_strat) ## selection probabilities for a balanced design

    Rs = sapply(1:ndraws, function(x)BSS(samp.fracs=sel.prob.bal,n2,data,strataformula)$R)

  } else if (design == "comb" ){

    ### To generalize the concept of a combined design for more than 3 groups for Y, we'll always set to selection probability for the median value(s) to 0

    ### Find the strata that will have zero probability of being selected (the median)
    if(  (length(S_unique)%%2) == 1  ){ ### calculate the median depending on whether the number of strata for S is even or odd
      med_stra <- median(unique(data$S))
    }else {
      med_stra <- as.numeric(quantile(S_unique,0.5,type=1))
      med_stra <- c(med_stra,med_stra+1)
    }
    indx_com <- which(!strat_df$S%in%med_stra)
    n_strat_com <- length(indx_com)

    sel.prob.com <- rep(0,n_strat) ## selection probabilities for a combined design
    sel.prob.com[indx_com] <- 1/n_strat_com

    Rs = sapply(1:ndraws, function(x)BSS(samp.fracs=sel.prob.com,n2,data,strataformula)$R)
  }

  ### Rs above should correspond to a matrix of N rows and ndraws columns with zeros and ones
  return ( Rs[,,drop=FALSE] )
}

