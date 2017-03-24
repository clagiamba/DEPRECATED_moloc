#library(mvtnorm)

##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##' 
##' @title Estimate trait variance, internal function
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##' 
##' @author Chris Wallace
sdY.est <- function(vbeta, maf, n) {
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  if(coef(m)[["oneover"]] < 0)
    stop("Trying to estimate trait variance from betas, and getting negative estimate.  Something is wrong.  You can 'fix' this by supplying an estimate of trait standard deviation yourself, as sdY=<value> in the dataset list.")
  return(sqrt(coef(m)[["oneover"]]))
}

#' variance of MLE of beta for quantitative trait, assuming var(y)=0
#'
#' Internal function
#' @title Var.data
#' @param f minor allele freq
#' @param N sample number
#' @return variance of MLE beta
#' @author Claudia Giambartolomei
#' @keywords internal
Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

#' Internal function, stats_p
#'
#' Get beta, SE from p
#' @title Internal function, stats_p
#' @param p p value
#' @param f MAF
#' @param type "quant" or "cc"
#' @param N sample size
#' @param s proportion of samples that are cases, ignored if type=="quant"
#' @param suffix suffix to append to column names of returned data.frame
#' @return data.frame containing beta, SE
#' @author Claudia Giambartolomei
#' @keywords internal
stats_p <- function(p, n, maf) {
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    var_mle <- 1/(2*maf*(1-maf) * ( n + z^2))
    # var_mle <- 1 / (2 * maf * (1 - maf) * n)
    SE <- sqrt(var_mle)
    BETA = z * SE
    # vars = 2 * maf * ( 1 - maf) * n * SE^2 * (n -1) + 2 * maf * ( 1- maf) * n * BETA^2 
    # vars = sqrt(median(vars/(n-1)))
    # Neff_est <- vars^2/(2*maf*(1-maf)*SE^2) - (BETA^2/SE^2) +1
    df <- cbind.data.frame(BETA, SE)
    return(df)
}
 
#' Internal function, stats_p_cc
#'
#' Internal function
#' @title stats_p_cc
#' @inheritParams Var.data
#' @param s proportion of samples that are cases
#' @param f minor allele freq
#' @param N sample number
#' @return variance of MLE beta
#' @examples
#' stats = stats_p_cc(listData[[1]]$PVAL, listData[[1]]$N, listData[[1]]$F, listData[[1]]$Ncases/listData[[1]]$N)
#' @author Claudia Giambartolomei
#' @keywords internal
stats_p_cc <- function(p, n, maf, s) {
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    var_mle <- (s*(1-s))/(2*maf*(1-maf) * ( n + z^2))
    SE <- sqrt(var_mle)
    BETA = z * SE
    df <- cbind.data.frame(BETA, SE)
    return(df)
}

#' Internal function, stats_p_cc
#'
#' Internal function
#' @title stats_p_cc
#' @inheritParams stats_p
#' @inheritParams stats_p_cc
#' @param s proportion of samples that are cases
#' @param f minor allele freq
#' @param N sample number
#' @return variance of MLE beta
#' @keywords internal
get_stats <- function(x) {
      message("Var and Beta from pvalues")
      if (!all(c("MAF", "N") %in% names(x))) stop("Must provide freq and N")
      if (all(c("N", "Ncases") %in% names(x))) {
         message("This is a case-control")
         x <- cbind.data.frame(x, stats_p_cc(x$PVAL, x$N, x$MAF, x$Ncases/x$N))
      } else {
         message("This is not a case-control")
         x <- cbind.data.frame(x, stats_p(x$PVAL, x$N, x$MAF))
      }
}

#' Internal function, approx.bf.p
#'
#' Calculate approximate Bayes Factors
#' @title Internal function, approx.bf.p
#' @param p p value
#' @param f MAF
#' @param type "quant" or "cc"
#' @param N sample size
#' @param s proportion of samples that are cases, ignored if type=="quant"
#' @param suffix suffix to append to column names of returned data.frame
#' @return data.frame containing lABF and intermediate calculations
#' @author Claudia Giambartolomei
#' @keywords internal
approx.bf.p <- function(p,f,type, N, s, suffix=NULL) {
  if(type=="quant") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  r <- sd.prior^2 / (sd.prior^2 + V)
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}

#' Internal function, logsum
#'
#' This function calculates the log of the sum of the exponentiated
#' logs taking out the max, i.e. insuring that the sum is not Inf
#' @title logsum
#' @param x numeric vector
#' @return max(x) + log(sum(exp(x - max(x))))
#' @author Claudia Giambartolomei, Vincent Plagnol
#' @keywords internal
logsum <- function(x) {
  my.max <- max(x)                              #take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

#' Internal function, logdiff
#'
#' This function calculates the log of the difference of the exponentiated
#' logs taking out the max, i.e. insuring that the difference is not negative
#' @title logdiff
#' @param x numeric
#' @param y numeric
#' @return max(x) + log(exp(x - max(x,y)) - exp(y-max(x,y)))
#' @author Chris Wallace
#' @keywords internal
logdiff <- function(x,y) {
  my.max <- max(x,y)                              #take out the maximum value in log form
  if (x>y) {
    my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  }
  if (x<y) {
    my.res <- my.max + log(exp(y-my.max) - exp(x - my.max ))
  }
  # If both x and y are zero, return zero?
  if (sum(x,y)==0) my.res =0
  return(my.res)
}

#' Internal function, is_pos_def
#'
#' This function tests if a matrix is positive definite
#' @title is_pos_def
#' @param A matrix
#' @author Claudia Giambartolomei, Jimmy Liu
#' @keywords internal
is_pos_def <- function(A) {
    cholStatus <- try(u <- chol(A), silent = TRUE)
    cholError <- ifelse(class(cholStatus) == "try-error", FALSE, TRUE)
    return(cholError)
    }


#' Internal function, get_psd
#'
#' This function attempts to make a matrix positive definite
#' @title get_psd
#' @param A matrix
#' @author Claudia Giambartolomei, Jimmy Liu
#' @keywords internal
get_psd <- function(A) {
    d = diag(A)
    # detA = np.linalg.det(A)
    while (!is_pos_def(A)) {
        # while detA <= 0:
        A = 0.999 * A
        for (i in 1:length(A[1,])){
            A[i, i] = d[i]
            # detA = np.linalg.det(A)
            }
        }
    return(A)
    }


#' Internal function, get_combos
#'
#' Find all possible combinations for the traits
#' @title get_combos
#' @param x names of each data frame
#' @return character vector with all combinations of names
#' @author Claudia Giambartolomei
#' @keywords internal
get_combos <- function(x) {
  x2 <- unlist(lapply(seq_along(x), function(m) 
    sapply(combn(x, m, simplify=FALSE), paste0, collapse='')))
  x3 <- expand.grid(rep(list(x2), length(x)))
  x4 <- sapply(apply(x3, 1, unique), function(x) paste0(sort(x), collapse=','))
  unique(grep('.*([^,]).*\\1', x4, val=TRUE, invert=TRUE))
}


#' Internal function, make_sigma
#'
#' This function creates Wakefield's ABF for each SNP, taking account of the correlation between traits
#' @title make_sigma
#' @param cor_mat matrix, correlation between the traits for each SNP
#'     If assuming no overlap, cor_mat is identity matrix.
#' @return matrix of ABF adjusted for correlations
#' @author Jimmy Liu
#' @keywords internal
make_sigma <- function(cor_mat, v, w){
    sigma = diag(v+w)
    for(i in 1:(length(v)-1)) {
        for(j in (i+1):length(v)) {
            c = cor_mat[i, j]
            sigma[i, j] = c * sqrt((v[i] + w[i]) * (v[j] + w[j]))
            sigma[j, i] = c * sqrt((v[i] + w[i]) * (v[j] + w[j]))
        }
    }
    if (!is_pos_def(sigma)) {
            message("Matrix is not positive definite")
            sigma <- get_psd(sigma)
        }

    return(sigma)
}


#' Adjusted Bayes factors for each SNP and each single association and combinations of
#' sharing/not sharing of causal variants across the datasets in the input data
#'
#' @title adjust_bfs
#' @param listData A list of data frames to be analyzed. Each data frame
#'     must contain columns named SNP (the SNP name consistent across all the data frames), 
#'     BETA, 
#'     SE
#' @param overlap Logical, do the individuals in the datasets overlap?
#' @param prior_var One or more numbers for the prior variance of the ABF.
#' @param compute_sdY Estimate standard deviation of Y. If false assume sdY is 1
#' @return An array containing the log adjusted Bayes Factors for each SNP and
#'         for each SNP and each configuration combination
#' @examples
#' ABF <- adjust_bfs(listData, overlap=FALSE, prior_var=c(0.01, 0.1, 0.5))
#' to view all configuration combinations for a SNP called "bp38878088"
#' ABF["bp38878088",,]
#' 
#' @keywords internal
#' @author Jimmy Liu, Claudia Giambartolomei
adjust_bfs <- function(listData, overlap=FALSE, prior_var=0.15, compute_sdY = FALSE, from_p=FALSE){
  # keep only common SNPs in all data: 
  listData <- lapply(listData, function(x) x[x$SNP %in% Reduce(intersect, Map("[[", listData, "SNP")), ])
  if (nrow(listData[[1]])==0) stop("There are no common SNPs in the datasets: check that SNP names are consistent")       
  listData <- lapply(listData, function(df){ df[order(df$SNP),]})
  ##### This is not right
  if (from_p) {
    listData <- lapply(listData, function(x) {
         y = get_stats(x)
         return(y)
         })    # if (!all(c("BETA", "SE") %in% names(x)))
  }       
  #######
  n_files <- length(listData)
  d <- letters[1:n_files]
  configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
  configs <- do.call(paste0, configs_cases)

    # adjusted bfs
    # correlation matrix
    # if assuming no overlap, cor_mat is identity matrix
    cor_mat <- diag(1, n_files)
    if (overlap) { 
        # for each combination of traits, i.e. 
        #loop through unique combination
        for(i in 1:(n_files -1)) {
            for(j in (i+1):n_files) {
                beta1 = listData[[i]]$BETA
                beta2 = listData[[j]]$BETA
                cortest = cor.test(beta1,beta2) # pvalue extracted from cor.test is not exact...
                if (cortest$p.value >= 0.01) { 
                    cor_mat[i, j] = 0.0
                    cor_mat[j, i] = 0.0
                } else {
                    cor_mat[i, j] = cortest$estimate
                    cor_mat[j, i] = cortest$estimate
                } 
        }
       }
    }
    # configs - different adjusted ABF depending on configs
    # get adjusted ABF for each config combo
    snps <- as.character(listData[[1]]$SNP)
    ABF <- array(,dim=c(length(snps),1,length(configs)))
    dimnames(ABF) <- list(snps,NULL,configs) 
    # grid of priors?
    # grid_priors <- matrix(0, nrow(configs_cases), ncol(configs_cases))
    if (compute_sdY) {
        sdY = mapply(sdY.est, lapply(lapply(listData, "[[", "SE"), '^',2), lapply(listData, "[[", "MAF"),  lapply(listData, "[[", "N"), SIMPLIFY = FALSE)
        varY = lapply(sdY, '^',2)
        names(varY) = d
        varY_configs <- apply(configs_cases, 1, function(x) prod(unlist(varY[unlist(x)])) )
    }

    names(listData) <- d
    var = lapply(listData, "[[", "SE")
    var = lapply(var, '^',2)
    var = do.call(cbind, var)

    betas = lapply(listData, "[[", "BETA")
    betas = do.call(cbind, betas)

    # data_sets <- apply(configs_cases,1, function(x) do.call(cbind.data.frame, listData[unlist(x[x!=''])]))
    means <- rep(0, n_files)
    w_0 <- rep(0, n_files)
    sigma_h0 <- lapply(split(var, seq(NROW(var))), FUN=make_sigma, cor_mat=cor_mat, w=w_0)

    pdf_null <- c()
    for (j in 1:length(sigma_h0)) {
         x <- dmvnorm(betas[j,], means, sigma_h0[[j]])
         pdf_null <- c(pdf_null, x)
    }

    for (i in 1:nrow(configs_cases)) {
               config <- configs[i]
               adj_i_average <- data.frame(matrix(ncol=0,nrow=length(snps)),stringsAsFactors=FALSE)
               for (w_i in prior_var) {
                  w <- as.numeric(ifelse(nchar(as.matrix(configs_cases[i,])), w_i, 0))
                  if (compute_sdY) {
                      w <- w * varY_configs[i]
                  }

                  sigma_h1 <- lapply(split(var, seq(NROW(var))), FUN=make_sigma, cor_mat=cor_mat, w=w)
                  pdf_alt <- c()
                  for (j in 1:length(sigma_h0)) {
                      x <- dmvnorm(betas[j,], means, sigma_h1[[j]])
                      pdf_alt <- c(pdf_alt, x)
                  }
                  adj_bf <- pdf_alt / pdf_null

                  adj_bf <- mapply(function(x, y)
                          if (y ==0 | y<1e-300) res <- 1.0 else res <- x / y, 
                          pdf_alt, pdf_null)

                  adj_i_average<- cbind(adj_i_average, adj_bf)
               }  
               adj_bf <- apply(adj_i_average, 1, mean)

               if (any(log(adj_bf)== "Inf") | any(log(adj_bf)== "-Inf") ) stop("Division impossible for adj_bf")

          ABF[,1,config] <- log(adj_bf)
          }
    return(ABF)
    }

#' Likelihood frame and posterior probability for combinations of sharing/not sharing
#' of causal variants across the datasets in the input data
#'
#' @title config_coloc
#' @param ABF is an array containing the single and coloc logBFs
#' @param n_files is one number, i.e. the number of traits to be analyzed
#' @param priors are numbers, the prior for one variant to be associated with 1 trait and with each additional trait
#' @return A data frame containing the likelihoods and posteriors for each configuration combination
#' @examples
#' lkl <- config_coloc(ABF, n_files=3, priors=c(1e-04, 1e-06, 1e-07))
#' 
#' @keywords internal
#' @author Jimmy Liu, Claudia Giambartolomei
config_coloc <- function(ABF, n_files, priors){
    # n_files = max(nchar(ABF))
    # have as many priors as traits
    d <- letters[1:n_files]
    configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
    configs <- do.call(paste0, configs_cases)
    final_configs <- get_combos(d) # includes the non-coloc configs

    # store sum of bfs for each configurarion, prior, and prior*sum(bfs_config)
    config_lkl <- array(,dim=c(length(final_configs),3,1), dimnames=list(final_configs,NULL,NULL))
    nsnps <- nrow(ABF)
    
    # likelihoods for colocalized configurations
    for (config in final_configs) {
        if (config %in% configs) {

        prior <- priors[nchar(config)]
        config_lkl[config,1,] <- prior
        #single_bfs[config,,] <- prior * sum(trait_bfs)
        config_lkl[config,2,] <- logsum(ABF[,,config])
        }
        }
    # likelihoods for configurations that can be derived from colocalized configurations
     for (config in final_configs) {
        if (!config %in% configs) {
         # print(config)
         # lH3.abf <-  logdiff(lH1.abf + lH2.abf, lH4.abf)
            # left side
            left_trait_bfs <- 0
            prior_num <- 1.0
            composite_bfs <- unlist(strsplit(config, split=","))
            for (i in composite_bfs) { # gather info on relevant single_bfs
               # print(i)
               left_trait_bfs <- left_trait_bfs + config_lkl[i,2,]
               prior_num <- prior_num * config_lkl[i,1,]
               }
            # relevant coloc configuration
            coloc_config <- gsub(",", "", config)
            # must sort otherwise it doesn't find the match
            coloc_config <- as.character(lapply(lapply(strsplit(coloc_config,NULL),sort),paste,collapse=""))
            # right_trait_bfs <- ( prior_num / single_bfs[coloc_config,1,] ) * single_bfs[coloc_config,3,]
            # right_trait_bfs <- prior_num * single_bfs[coloc_config,3,]
            right_trait_bfs <- config_lkl[coloc_config,2,]
            # print(left_trait_bfs)
            # print(right_trait_bfs)
            config_bf <- logdiff(left_trait_bfs, right_trait_bfs)
            config_lkl[config,2,] <- config_bf
            config_lkl[config,1,] <- prior_num
        }
        }
      config_lkl[,3,] <- log(config_lkl[,1,]) + config_lkl[,2,] # prior * Sum(BFs)

    # add config where nothing is associated
    config_ppas <-  as.data.frame(config_lkl)
    names(config_ppas)=c("prior", "sumbf", "loglkl")
    # add zero
    config_ppas <- rbind(config_ppas, zero=c(0.999, 0, 0))

    all.abf <- config_ppas$loglkl
    my.denom.log.abf <- logsum(all.abf)
    config_ppas$PPA <- exp(all.abf - my.denom.log.abf)

    # print(config_ppas)
    return(list(config_ppas, nsnps))
}


#' Posterior probability that each SNP is THE causal variant for a shared signal
#'
#' @title snp_ppa
#' @param ABF is an array containing the single and coloc logBFs
#' @param n_files is one number, i.e. the number of traits to be analyzed
#' @param priors is the prior for one variant to be associated with 1 trait and with each additional trait
#' @return A data frame containing the likelihoods and posteriors for each configuration combination
#' @examples
#' snp <- snp_ppa(ABF, n_files=3, config_ppas=lkl[[1]])
#' 
#' @keywords internal
#' @author Jimmy Liu, Claudia Giambartolomei
snp_ppa <- function(ABF, n_files, config_ppas){
    d <- letters[1:n_files]
    configs_cases <- do.call(expand.grid, lapply(d, function(x) c("", x)))[-1,]
    configs <- do.call(paste0, configs_cases)
    final_configs <- get_combos(d) # includes the non-coloc configs
    SNP.PP.H4.df <- as.data.frame(sapply(as.data.frame(ABF), function(x) exp(x-logsum(x)) ))
    names(SNP.PP.H4.df) = dimnames(ABF)[[3]]
    rownames(SNP.PP.H4.df) = dimnames(ABF)[[1]]

    #best.t <- exp(rowSums(as.data.frame(ABF)[,grep(i, colnames(as.data.frame(ABF)))]) - my.denom.log.abf)

    # printout
    # print out colocalization posteriors for each trait
    coloc_ppas = numeric()
    best.snp.coloc = character()
    # for (i in d) {
    for (i in configs) {
        final_config_ls <- strsplit(final_configs, ",")
        trait_coloc_list <- Filter(function(x) any(nchar(grep(i, x, value = TRUE))>1), final_config_ls)
        trait_coloc <- sapply(trait_coloc_list,function(x) paste(x, collapse=","))
        trait_coloc <- c(trait_coloc, grep("a[a-z]c", configs, perl=T, value=T))
        trait_coloc <- trait_coloc[!duplicated(trait_coloc)]
        # print(trait_coloc)
        hh4 <- sum(config_ppas[rownames(config_ppas) %in% trait_coloc,"PPA"])

        # what is the denominator here, other coloc models or all the coloc and non coloc?
        if (sum(names(SNP.PP.H4.df) %in% trait_coloc)>1) {
            SNP.coloc <- rowSums(SNP.PP.H4.df[,names(SNP.PP.H4.df) %in% trait_coloc]) 
            } else {
            SNP.coloc <- SNP.PP.H4.df[,names(SNP.PP.H4.df) %in% trait_coloc]
            }
        names(SNP.coloc) <- rownames(SNP.PP.H4.df)
        s = names(which.max(SNP.coloc))
        message("Best SNP per trait ", i, ": ", s)
        best.snp.coloc <- c(best.snp.coloc, s)
        #if (names(which.max(SNP.pp4)) != causal.1) browser("***causal not recognized")
        coloc_ppas <- c(coloc_ppas, hh4)
        message("Probability that trait ", i, " colocalizes with at least one other trait = ", signif(hh4, digits = 2))
     }
     #SNP.pp4.any <- rowSums(ABF[,,dimnames(ABF)[[3]] %in% any_coloc])/denom
     #message("Best SNP for any trait colocalizing: ", names(which.max(SNP.pp4.any)))
     
     # names(coloc_ppas) <- d
     best_snp = cbind.data.frame(coloc_ppas, best.snp.coloc)
     rownames(best_snp) <- configs
     return(best_snp)
}


#' Bayesian multiple trait colocalization analysis using list of data.frames
#'
#' Runs \code{\link{adjust_bfs}}, \code{\link{config_coloc}}, \code{\link{snp_ppa}}
#' 
#' @title moloc_test
#' @param listData A list of data frames to be analyzed. Each data frame
#'     must contain columns named SNP (the SNP name consistent across all the data frames), 
#'     BETA, 
#'     SE
#' @param overlap Logical, do the individuals in the datasets overlap?
#' @param prior_var One or more numbers for the prior variance of the ABF.
#' @param priors are numbers, the prior for one variant to be associated with 1 trait and with each additional trait
#' @param compute_sdY Estimate standard deviation of Y. If false assume sdY is 1
#' 
#' @param ... parameters passed to \code{\link{adjust_bfs}}, \code{\link{config_coloc}}, \code{\link{snp_ppa}}
#' @return A list of three elements: 
#'         First is a data frame with 4 variables: priors ('prior'), likelihoods ('sumbf') and 
#'           posteriors ('loglkl' and 'PPA') for each configuration; 
#'         Second is a number, the number of SNPs in common in the region;
#'         Third is a data frame with 2 variables: 
#'           SNP with the best posterior for each scenario ('coloc_ppas'), 
#'           SNP name with the best posterior for each scenario ('best.snp.coloc').
#' @export
#' @author Claudia Giambartolomei
#' @examples
#' moloc <- moloc_test(listData) # uses default priors
#' 
#' @author Claudia Giambartolomei
moloc_test <- function(listData, overlap=FALSE, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07), compute_sdY=FALSE, from_p=FALSE) {
    n_files <- length(listData)
    if(missing(priors)) {
      priors <- 10^-(seq(from=4, to=4+n_files-1, by=1))
      message("Use default priors: ", paste(priors, collapse=","))
    } else {
      if (length(priors)!=n_files) stop("Priors need to be the same length as the number of traits")
      priors <- as.numeric(priors)
    }

    ABF <- adjust_bfs(listData, overlap, prior_var, from_p)
    lkl <- config_coloc(ABF, n_files, priors)
    snp <- snp_ppa(ABF, n_files=n_files, config_ppas=lkl[[1]])
    nsnps <- lkl[[2]]
    lkl <- lkl[[1]]
    res <- list(lkl, nsnps, snp)
    names(res) <- c("priors_lkl_ppa", "nsnps", "best_snp")
    return(res)
}
