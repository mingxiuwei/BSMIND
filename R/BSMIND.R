#'
#' BSMIND Model Main Function
#'
#' BSMIND: Bayesian integration of spatial multiomics for consensus domain identification
#'
#' The function BSMIND implements a multiomics domain identification method for spatial multiomics data.
#'     It normalizes raw count matriceso f multiple omics sharing identical
#'     spatial coordinates without taking logarithm throughout the pipeline. The method can
#'     simultaneously quantify and output the contribution weight of each omics layer to spatial region identification.
#'


#'
#' Internal auxiliary function
#' @noRd
trunc_rbeta <- function(num, s1, s2, w1, w2, by = 10^(-7)) {
  x <- seq(w1, w2, by = by)
  x <- x[-1]
  lf <- (s1 - 1) * log(x) + (s2 - 1) * log(1 - x)
  max_lf <- max(lf)
  ss <- sum(exp(lf - max_lf))
  prob <- exp(lf - max_lf) / ss
  elements <- sample(x,
                     size = num,
                     prob = prob,
                     replace = TRUE)
  return(elements)
}

#' @noRd
ComputePottsDist <- function(tmpBeta, c_vec, nei_idxs, numOfData) {
  sum_nei_m = 0
  for (i in 1:numOfData) {
    tmpNei = nei_idxs[[i]]
    sum_nei_m = sum_nei_m + sum(c_vec[tmpNei] == c_vec[i])
  }
  return(sum_nei_m / 2 * tmpBeta)
}

#'
#' Internal auxiliary function
#' @noRd
FindNeighbors_R <- function(i, X_loc, Y_loc, platform) {
  x_tmp <- X_loc[i]
  y_tmp <- Y_loc[i]
  numOfData <- length(X_loc)
  
  if (platform == "ST") {
    candidates_idx <- which((X_loc == x_tmp & abs(Y_loc - y_tmp) == 1) |
                              (Y_loc == y_tmp &
                                 abs(X_loc - x_tmp) == 1))
  } else if (platform == "Visium") {
    if (y_tmp %% 2 == 0) {
      candidates_idx <- which(
        (X_loc == x_tmp & Y_loc == y_tmp + 1) |
          (X_loc == x_tmp & Y_loc == y_tmp - 1) |
          (X_loc == x_tmp - 1 & Y_loc == y_tmp) |
          (X_loc == x_tmp + 1 & Y_loc == y_tmp) |
          (X_loc == x_tmp + 1 & Y_loc == y_tmp + 1) |
          (X_loc == x_tmp - 1 & Y_loc == y_tmp - 1)
      )
    } else {
      candidates_idx <- which(
        (X_loc == x_tmp & Y_loc == y_tmp + 1) |
          (X_loc == x_tmp & Y_loc == y_tmp - 1) |
          (X_loc == x_tmp - 1 & Y_loc == y_tmp) |
          (X_loc == x_tmp + 1 & Y_loc == y_tmp) |
          (X_loc == x_tmp - 1 & Y_loc == y_tmp + 1) |
          (X_loc == x_tmp + 1 & Y_loc == y_tmp - 1)
      )
    }
  } else {
    stop("platform must be 'ST' or 'Visium'")
  }
  return(as.integer(candidates_idx))
}

#'
#' Internal auxiliary function
#' @noRd
GetAllNeighbors <- function(X_loc, Y_loc, platform) {
  numOfData <- length(X_loc)
  nei_idxs <- vector("list", numOfData)
  
  for (i in 1:numOfData) {
    nei_idxs[[i]] <- FindNeighbors_R(i, X_loc, Y_loc, platform)
  }
  return(nei_idxs)
}


#' BSMIND function
#' @importFrom stats kmeans rnorm runif dnorm
#' @importFrom MCMCpack rinvgamma rdirichlet ddirichlet
#' @importFrom truncnorm rtruncnorm dtruncnorm
#' @param data_1 First n*p1 data matrix (spots x genes), gene expression modality.
#' @param data_2 Second n*p2 data matrix (spots x features), secondary feature modality.
#' @param coord Data frame of spot coordinates, 2 columns: x, y coordinates.
#' @param platform Spatial sequencing platform; "ST" (square grid) or "Visium" (hex grid).
#' @param K_fixed Fixed number of regions for all latent layers. Default is 3.
#' @param a_mu Mean of normal prior for gene expression cluster means mu_gk. Default is 0.
#' @param b_mu Standard deviation of normal prior for mu_gk. Default is 1.
#' @param IGkappa Shape parameter of inverse-gamma prior for gene expression sigma_g. Default is 2.
#' @param IGtau Scale parameter of inverse-gamma prior for sigma_g. Default is 10.
#' @param a_mu2 Mean of normal prior for second modality cluster means mu_pk. Default is 2.
#' @param b_mu2 Standard deviation of normal prior for mu_pk. Default is 1.
#' @param IGkappa2 Shape parameter of inverse-gamma prior for second modality sigma_p. Default is 2.
#' @param IGtau2 Scale parameter of inverse-gamma prior for sigma_p. Default is 10.
#' @param a_beta Prior mean for spatial Potts interaction parameter beta. Default is 0.7.
#' @param tau_beta Prior standard deviation for beta. Default is 1.
#' @param dpAlpha Dirichlet concentration parameter for cluster proportion prior pi_vec. Default is 1.
#' @param minPsi1 Lower bound of uniform prior for gene expression contribution weight psi1. Default is 0.8.
#' @param maxPsi1 Upper bound of uniform prior for gene expression contribution weight psi1. Default is 0.95.
#' @param minPsi2 Lower bound of uniform prior for secondary feature modality contribution weight psi2. Default is 0.8.
#' @param maxPsi2 Upper bound of uniform prior for secondary feature modality contribution weight psi2. Default is 0.95.
#' @param tau0 Scale parameter for Dirichlet proposal distribution of pi_vec. Default is 2.
#' @param tau1 Proposal standard deviation for Potts parameter beta. Default is 0.05.
#' @param numOfMCMC Total MCMC iterations including burn-in. Default is 5000.
#' @param burnIn Number of burn-in iterations (discarded for posterior inference). Default is 2500.
#' @param Is_print Logical; if TRUE, print iteration progress. Default is TRUE.
#' @param print_gap Iteration interval for progress printing. Default is 10.
#' @param Is_random_seed Logical; if TRUE, set fixed random seed for reproducibility. Default is FALSE.
#' @param random_seed Random seed number used if Is_random_seed=TRUE. Default is 2234.
#'
#' @return List containing full and post-burn-in MCMC posterior samples:
#' \item{mu_gk_mcmc}{Posterior samples of gene expression cluster means}
#' \item{sigma_g_mcmc}{Posterior samples of gene expression standard deviations}
#' \item{mu_pk_mcmc}{Posterior samples of second modality cluster means}
#' \item{sigma_p_mcmc}{Posterior samples of second modality standard deviations}
#' \item{clIds_mcmc}{Posterior samples of integrative cluster labels}
#' \item{L1_Ids_mcmc}{Posterior samples of gene-specific cluster labels}
#' \item{L2_Ids_mcmc}{Posterior samples of second modality-specific cluster labels}
#' \item{pottsBeta_mcmc}{Posterior samples of spatial interaction parameter beta}
#' \item{psi1_mcmc}{Posterior samples of modality 1 weight psi1}
#' \item{psi2_mcmc}{Posterior samples of modality 2 weight psi2}
#' \item{pi_mcmc}{Posterior samples of cluster proportion vector}
#' \item{all_mu_gk_mcmc}{List, all posterior samples including burn-in of gene expression cluster means}
#' \item{all_sigma_g_mcmc}{List, all posterior samples including burn-in of gene expression standard deviations}
#' \item{all_mu_pk_mcmc}{List, all posterior samples including burn-in of second modality cluster means}
#' \item{all_sigma_p_mcmc}{List, all posterior samples including burn-in of second modality standard deviations}
#' \item{all_clIds_mcmc}{List, all posterior samples including burn-in of integrative cluster labels}
#' \item{all_L1_Ids_mcmc}{List, all posterior samples including burn-in of gene-specific cluster labels}
#' \item{all_L2_Ids_mcmc}{List, all posterior samples including burn-in of second modality-specific cluster labels}
#' \item{all_pottsBeta_mcmc}{Vector, all posterior samples including burn-in of spatial interaction parameter beta}
#' \item{all_psi1_mcmc}{Vector, all posterior samples including burn-in of modality 1 weight psi1}
#' \item{all_psi2_mcmc}{Vector, all posterior samples including burn-in of modality 2 weight psi2}
#' \item{all_pi_mcmc}{List, all posterior samples including burn-in of cluster proportion vector}
#' \item{exeTime}{Total model execution time}
#' @export
BSMIND <- function(data_1,
                   data_2,
                   coord,
                   platform = c("ST", "Visium"),
                   K_fixed = 3,
                   a_mu = 0,
                   b_mu = 1,
                   IGkappa = 2,
                   IGtau = 10,
                   a_mu2 = 2,
                   b_mu2 = 1,
                   IGkappa2 = 2,
                   IGtau2 = 10,
                   a_beta = 0.7,
                   tau_beta = 1,
                   dpAlpha = 1,
                   minPsi1 = 0.8,
                   maxPsi1 = 0.95,
                   minPsi2 = 0.8,
                   maxPsi2 = 0.95,
                   tau0 = 2,
                   tau1 = 0.05,
                   numOfMCMC = 1000,
                   burnIn = 500,
                   Is_print = TRUE,
                   print_gap = 10,
                   Is_random_seed = F,
                   random_seed = 2234) {
  # two types of data
  gene_data = as.matrix(data_1)   # spots x genes
  data_2 = as.matrix(data_2)
  
  # location coordinates
  X_loc = as.integer(coord[, 1])
  Y_loc = as.integer(coord[, 2])
  
  # data shapes
  numOfData = nrow(gene_data)
  G = ncol(gene_data)
  P = ncol(data_2)
  
  if (Is_random_seed)
    set.seed(random_seed)
  
  
  # warm start using KMeans
  kmres_L1 <- kmeans(x = gene_data, centers = K_fixed)
  kmres_L2 <- kmeans(x = data_2, centers = K_fixed)
  
  # label switching
  init_label1 = kmres_L1[["cluster"]]
  uniq = unique(init_label1)
  for (i in 1:length(uniq)) {
    init_label1[kmres_L1[["cluster"]] == uniq[i]] = i
  }
  
  init_label2 = kmres_L2[["cluster"]]
  uniq = unique(init_label2)
  for (i in 1:length(uniq)) {
    init_label2[kmres_L2[["cluster"]] == uniq[i]] = i
  }
  L1_Ids = init_label1
  L2_Ids = init_label2
  clIds  = L1_Ids
  
  K1_max = max(L1_Ids)
  K2_max = max(L2_Ids)
  
  # initialize
  mu_gk <- matrix(
    data = rnorm(G * K1_max, mean = a_mu, sd = b_mu),
    nrow = G,
    ncol = K1_max
  )
  sigma_g = sqrt(rinvgamma(G, shape = IGkappa, scale = IGtau))
  mu_pk <- matrix(
    data = rnorm(P * K2_max, mean = a_mu2, sd = b_mu2),
    nrow = P,
    ncol = K2_max
  )
  sigma_p = sqrt(rinvgamma(P, shape = IGkappa2, scale = IGtau2))
  pottsBeta = rtruncnorm(1, a = 0, mean = a_beta, sd = tau_beta)
  psi1 = runif(1, minPsi1, maxPsi1)
  psi2 = runif(1, minPsi2, maxPsi2)
  pi_vec <- rdirichlet(1, rep(dpAlpha, K_fixed))
  
  
  # store MCMC samples (numOfMCMC - burnIn)
  mu_gk_mcmc <- vector(mode = "list", length = numOfMCMC)
  sigma_g_mcmc <- vector(mode = "list", length = numOfMCMC)
  mu_pk_mcmc <- vector(mode = "list", length = numOfMCMC)
  sigma_p_mcmc <- vector(mode = "list", length = numOfMCMC)
  clIds_mcmc <- vector(mode = "list", length = numOfMCMC)
  L1_Ids_mcmc <- vector(mode = "list", length = numOfMCMC)
  L2_Ids_mcmc <- vector(mode = "list", length = numOfMCMC)
  pottsBeta_mcmc = array(0, dim = numOfMCMC)
  psi1_mcmc = array(0, dim = numOfMCMC)
  psi2_mcmc = array(0, dim = numOfMCMC)
  pi_mcmc <- vector(mode = "list", length = numOfMCMC)
  
  nei_idxs <- GetAllNeighbors(X_loc = X_loc,
                              Y_loc = Y_loc,
                              platform = platform)
  
  
  
  cat(paste0("=== MCMC Iterations ===\n"))
  ssTime = Sys.time()
  
  for (mcmc in 1:numOfMCMC) {
    #------------------------------------------------------------------------------
    ### ==== TBC
    # tmplen = length(dpXi)
    # dpXi_proposal = numeric(tmplen)
    # dpXi_proposal[1:M0] = rtruncnorm(M0, a = 0, b = 1, mean = dpXi[1:M0], sd = tau0)
    pi_proposal = rdirichlet(1, tau0 * pi_vec)
    beta_proposal = rtruncnorm(1, a = 0, mean = pottsBeta, sd = tau1)
    # if (tmplen > M0) dpXi_proposal[(M0 + 1):tmplen] = rbeta(tmplen - M0, 1, dpAlpha)
    
    if (!beta_proposal < 0) {
      # Step 2: Propose new C* from P(C | beta*, {pi_k*})
      aux_c = clIds
      for (i in 1:numOfData) {
        tmpNei = nei_idxs[[i]]
        aux_p = numeric(K_fixed)
        for (m in 1:K_fixed) {
          aux_p[m] = pi_proposal[m] * exp(beta_proposal * sum(aux_c[tmpNei] == m))
        }
        aux_c[i] = sample(1:K_fixed, 1, prob = aux_p)
      }
      
      
      # Step 3: accept proposal with probability min(1,r)
      # log ratio
      # prior
      logfrac1 = log(ddirichlet(pi_proposal, rep(dpAlpha, K_fixed))) +
        log(dtruncnorm(
          beta_proposal,
          a = 0,
          mean = a_beta,
          sd = tau_beta
        ))
      logfrac2 = log(ddirichlet(pi_vec, rep(dpAlpha, K_fixed))) +
        log(dtruncnorm(
          pottsBeta,
          a = 0,
          mean = a_beta,
          sd = tau_beta
        ))
      
      # C
      logfrac1 = logfrac1 +
        sum(log(pi_proposal[clIds])) + sum(log(pi_vec[aux_c])) +
        ComputePottsDist(beta_proposal, clIds, nei_idxs, numOfData) + ComputePottsDist(pottsBeta, aux_c, nei_idxs, numOfData)
      logfrac2 = logfrac2 +
        sum(log(pi_vec[clIds])) + sum(log(pi_proposal[aux_c])) +
        ComputePottsDist(pottsBeta, clIds, nei_idxs, numOfData) + ComputePottsDist(beta_proposal, aux_c, nei_idxs, numOfData)
      
      # proposal
      logfrac1 = logfrac1 +
        log(ddirichlet(pi_vec, tau0 * pi_proposal)) +
        log(dtruncnorm(
          pottsBeta,
          a = 0,
          mean = beta_proposal,
          sd = tau1
        ))
      logfrac2 = logfrac2 +
        log(ddirichlet(pi_proposal, tau0 * pi_vec)) +
        log(dtruncnorm(
          beta_proposal,
          a = 0,
          mean = pottsBeta,
          sd = tau1
        ))
      
      ratio = logfrac1 - logfrac2
      ratio = exp(ratio)
      prob = min(1, ratio)
      
      tmpU = runif(1)
      if (!is.na(prob) && tmpU < prob) {
        pottsBeta = beta_proposal
        pi_vec = pi_proposal
      } else {
        # print("reject")
      }
    }
    
    
    # update psi1 and psi2
    num_eq1 = sum(L1_Ids == clIds)
    num_eq2 = sum(L2_Ids == clIds)
    psi1 = trunc_rbeta(1, 1 + num_eq1, 1 + numOfData - num_eq1, minPsi1, maxPsi1)
    psi2 = trunc_rbeta(1, 1 + num_eq2, 1 + numOfData - num_eq2, minPsi2, maxPsi2)
    
    ### ==== TBC
    
    
    
    ### update C
    for (i in 1:numOfData) {
      unnorm_prob <- numeric(K_fixed)
      
      for (k in 1:K_fixed) {
        nei_j <- nei_idxs[[i]]
        
        if (length(nei_j) > 0) {
          sum_nei_k <- sum(clIds[nei_j] == k)
        } else {
          sum_nei_k <- 0
        }
        
        log_spatial_potts <- pottsBeta * sum_nei_k
        
        if (L1_Ids[i] == k) {
          Prob_L1 <- psi1
        } else {
          Prob_L1 <- (1 - psi1) / (K_fixed - 1)
        }
        
        if (L2_Ids[i] == k) {
          Prob_L2 <- psi2
        } else {
          Prob_L2 <- (1 - psi2) / (K_fixed - 1)
        }
        unnorm_prob[k] <- log(pi_vec[k]) + log_spatial_potts + log(Prob_L1) + log(Prob_L2)
      }
      unnorm_prob = unnorm_prob - max(unnorm_prob)
      norm_prob <- exp(unnorm_prob) / sum(exp(unnorm_prob))
      clIds[i] <- sample(x = 1:K_fixed,
                         size = 1,
                         prob = norm_prob)
    }
    
    
    ### update L1
    for (i in 1:numOfData) {
      unnorm_prob_L1 <- numeric(K_fixed)
      
      for (k in 1:K_fixed) {
        if (k == clIds[i]) {
          Prob_L1 <- psi1
        } else {
          Prob_L1 <- (1 - psi1) / (K_fixed - 1)
        }
        
        log_likelihood_sum <- sum(dnorm(
          x = gene_data[i, ],
          mean = mu_gk[, k],
          sd = sigma_g,
          log = TRUE
        ))
        
        unnorm_prob_L1[k] <- log(Prob_L1) + log_likelihood_sum
        
      }
      unnorm_prob_L1 = unnorm_prob_L1 - max(unnorm_prob_L1)
      norm_prob_L1 = exp(unnorm_prob_L1) / sum(exp(unnorm_prob_L1))
      L1_Ids[i] = sample(1:K_fixed,
                         size = 1,
                         replace = TRUE,
                         prob = norm_prob_L1)
    }
    
    
    ### update L2
    for (i in 1:numOfData) {
      unnorm_prob_L2 <- numeric(K_fixed)
      
      
      for (k in 1:K_fixed) {
        if (k == clIds[i]) {
          Prob_L2 <- psi2
        } else {
          Prob_L2 <- (1 - psi2) / (K_fixed - 1)
        }
        
        log_likelihood_sum <- sum(dnorm(
          data_2[i, ],
          mean = mu_pk[, k],
          sd = sigma_p,
          log = TRUE
        ))
        
        unnorm_prob_L2[k] <- log(Prob_L2) + log_likelihood_sum
        
      }
      
      unnorm_prob_L2 = unnorm_prob_L2 - max(unnorm_prob_L2)
      norm_prob_L2 = exp(unnorm_prob_L2) / sum(exp(unnorm_prob_L2))
      L2_Ids[i] = sample(1:K_fixed,
                         size = 1,
                         replace = TRUE,
                         prob = norm_prob_L2)
      
    }
    
    
    ### update mu_gk
    for (g in 1:G) {
      for (k in 1:K1_max) {
        spot_ind_1 = which(L1_Ids == k)
        Nk = length(spot_ind_1)
        
        if (Nk != 0) {
          tmpSum = if (Nk > 1) {
            colSums(gene_data[spot_ind_1, ])
          } else{
            gene_data[spot_ind_1, ]
          }
          
          b_mu2_tilde = 1 / (1 / (b_mu)^2 + Nk / (sigma_g[g])^2)
          a_mu_tilde =  b_mu2_tilde * (a_mu / (b_mu)^2 + tmpSum[g]   / (sigma_g[g])^2)
          mu_gk[g, k] = rnorm(1, mean = a_mu_tilde, sd = sqrt(b_mu2_tilde))
          
        } else{
          b_mu2_tilde = 1 / (1 / b_mu^2 + Nk / (sigma_g[g])^2)
          mu_gk[g, k] = rnorm(1, mean = a_mu, sd = b_mu)
        }
        
      }
    }
    
    
    ###  update sigma_g
    for (g in 1:G) {
      tmpMu = mu_gk[g, L1_Ids]
      tmpSumSqu = sum((gene_data[, g] - tmpMu)^2) / 2
      sigma_g[g] = sqrt(rinvgamma(
        1,
        shape = numOfData / 2 + IGkappa,
        scale = IGtau + tmpSumSqu
      )) # sigma_g is std!!!
    }
    
    
    ### mu_pk
    for (p in 1:P) {
      for (k in 1:K2_max) {
        spot_ind_2 = which(L2_Ids == k)  # L2_Ids??
        Nk = length(spot_ind_2)
        
        if (Nk != 0) {
          tmpSum = if (Nk > 1) {
            colSums(data_2[spot_ind_2, ])
          } else{
            data_2[spot_ind_2, ]
          }
          
          b_mu2_tilde_p = 1 / (1 / (b_mu2)^2 + Nk / (sigma_p[p])^2)
          a_mu_tilde_p =  b_mu2_tilde_p * (a_mu2 / (b_mu2)^2 + tmpSum[p] / (sigma_p[p])^2)
          mu_pk[p, k] = rnorm(1, mean = a_mu_tilde_p, sd = sqrt(b_mu2_tilde_p))
          
        } else{
          b_mu2_tilde_p = 1 / (1 / (b_mu2)^2 + Nk / (sigma_p[p])^2)
          mu_pk[p, k] = rnorm(1, mean = a_mu2, sd = b_mu2)
        }
        
      }
    }
    
    
    ###  update sigma_p
    for (p in 1:P) {
      tmpMu = mu_pk[p, L2_Ids]
      tmpSumSqu = sum((data_2[, p] - tmpMu)^2) / 2
      sigma_p[p] = sqrt(rinvgamma(
        1,
        shape = numOfData / 2 + IGkappa2,
        scale = IGtau2 + tmpSumSqu
      )) # sigma_g is std!!!
    }
    
    
    ### Results
    clIds_mcmc[[mcmc]] = clIds
    L1_Ids_mcmc[[mcmc]] = L1_Ids
    L2_Ids_mcmc[[mcmc]] = L2_Ids
    
    mu_gk_mcmc[[mcmc]] <- mu_gk
    sigma_g_mcmc[[mcmc]] = sigma_g
    mu_pk_mcmc[[mcmc]] <- mu_pk
    sigma_p_mcmc[[mcmc]] = sigma_p
    
    pottsBeta_mcmc[mcmc] = pottsBeta
    psi1_mcmc[mcmc] = psi1
    psi2_mcmc[mcmc] = psi2
    pi_mcmc[[mcmc]] = pi_vec
    
    
    ### Output
    if (mcmc <= burnIn) {
      if (mcmc == 0)
        cat(" Burn-in:")
      else if (mcmc / print_gap == floor(mcmc / print_gap))
        cat(paste0("+++", mcmc))
    } else {
      if (mcmc == burnIn + 1)
        cat("\n MCMC sampling:")
      else if (mcmc > burnIn &
               mcmc / print_gap == floor(mcmc / print_gap))
        cat(paste0("...", mcmc))
    }
    
  }
  
  
  eeTime = Sys.time()
  exeTime = eeTime - ssTime
  cat(paste0("\n=== End Train ===\n"))
  
  
  ### return results list
  res_list = vector("list")
  
  # All MCMC iteration results
  res_list[["all_mu_gk_mcmc"]]      = mu_gk_mcmc
  res_list[["all_sigma_g_mcmc"]]    = sigma_g_mcmc
  res_list[["all_mu_pk_mcmc"]]      = mu_pk_mcmc
  res_list[["all_sigma_p_mcmc"]]    = sigma_p_mcmc
  res_list[["all_clIds_mcmc"]]      = clIds_mcmc
  res_list[["all_L1_Ids_mcmc"]]     = L1_Ids_mcmc
  res_list[["all_L2_Ids_mcmc"]]     = L2_Ids_mcmc
  res_list[["all_pottsBeta_mcmc"]]  = pottsBeta_mcmc
  res_list[["all_psi1_mcmc"]]       = psi1_mcmc
  res_list[["all_psi2_mcmc"]]       = psi2_mcmc
  res_list[["all_pi_mcmc"]]         = pi_mcmc
  
  # 2. Results after burn-in only
  res_list[["mu_gk_mcmc"]]      = mu_gk_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["sigma_g_mcmc"]]    = sigma_g_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["mu_pk_mcmc"]]      = mu_pk_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["sigma_p_mcmc"]]    = sigma_p_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["clIds_mcmc"]]      = clIds_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["L1_Ids_mcmc"]]     = L1_Ids_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["L2_Ids_mcmc"]]     = L2_Ids_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["pottsBeta_mcmc"]]  = pottsBeta_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["psi1_mcmc"]]       = psi1_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["psi2_mcmc"]]       = psi2_mcmc[(burnIn + 1):numOfMCMC]
  res_list[["pi_mcmc"]]         = pi_mcmc[(burnIn + 1):numOfMCMC]
  
  res_list[["exeTime"]]         = exeTime
  return(res_list)
}
