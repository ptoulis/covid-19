#' Analysis of COVID-19 serology studies in the US (Partial Identification approach).
#' https://bfi.uchicago.edu/wp-content/uploads/BFI_WP_202054.pdf
#' Panos Toulis (panos.toulis@chicagobooth.edu)
#' Version 5.0 -- May, 2020.
#' 
rm(list=ls())
require(plyr)
require(dplyr)

#' To try main test.
#' -----------------------
#' m = create_Model(1)  # 1=Santa Clara study. See "create_Model" for other options.
#' theta0 = c(1/100,90/100, 12)
#' include_theta_Exact(theta0, m, vis=TRUE)   # Test H0: 1% FPR, 90% TPR, 12 true infected.
#'
#' RETURNS: list(incl, strength, reduction_pct)
#'    where incl = whether we include theta0 or not in confidence set.
#'          strength = the larger than 0.05 the stronger the acceptance.
#'          reduction_pct = % of computational reduction by pruning the sample space (see Section 5.1)
#'
#' To try likelihood ratio test       
#' --------------------------------
#' LR_test(theta0, m)
#' 
#' Procedure 1 and 3 from Chen et al (2018) -- mainly for Santa Clara study.
#' -----------------------------------------
#' Proc1_chen(m, num_mcmc=20000)
#' Proc3_chen(m)



#' Creates an object that contains all information for analysis:
#' @param use_Sharp FALSE if to use classical (conservative construction), TRUE to use sharper construction.
#' @return List of (CalibrationStudy, MainStudy, Theta)
#' CalibrationStudy has validation data on test outcomes for which the patient condition is known.
#' MainStudy has the data from the main study where the underlying condition is not known.
#' Theta is the parameter space to be explored -- comprised of (FPR, TPR, prevalence)
create_Model = function(use_data_no, use_Sharp=TRUE, SC_version=1) {
  #' Calibration study. Remains fixed.
  stopifnot(use_data_no %in% c(1, 2, 3, 4, 7))
  #' Codes for different datasets.
  # 1 = Santa Clara
  # 2 =  LA county
  # 3 = Santa Clara + LA County
  # 4 = NYC
  # 7 = Santa Clara + LA County + NYC
  #' Calibration study.
  #' Based on v2 by the Santa Clara authors.
  Calibration = list(Neg=list(total=401, obs=2),
                     Pos=list(total=37+75+85, obs=25+75+78))
  #
  #' Based on version 2 by the Santa Clara authors.
  #' https://www.medrxiv.org/content/10.1101/2020.04.14.20062463v2
  #' See also https://statmodeling.stat.columbia.edu/2020/04/30/updated-santa-clara-study-of-coronavirus-infection/
  if(SC_version==2) {
    Calibration = list(Neg=list(total=3324, obs=16),
                       Pos=list(total=157, obs=130))
  }
  print(sprintf("> [INFO] Using Santa Clara validation version %d", SC_version))
  print(str(Calibration))
  #' Main study.
  #' 
  N_values =   c(3330, 846, 3330 + 846, 3000, NA, NA, 3330 + 846 + 3000)
  n_obs_values = c(50, 35,   50 + 35,   420, NA,  NA,  50 + 35 + 420)
  #                     # max prevalence 4%, 7%, 5%,...
  max_infect = as.integer(c(4, 7, 5, 20, 0, 0, 9) * N_values/100) 
  # c(100, 70, 120, 650, NA, NA, 800)
  #' Names of studies. SC+LA = Santa Clara and LA county combined.
  names = c("SC", "LA", "SC+LA", "NY", NA, NA, "SC+LA+NY")
  full_names = c("Santa Clara", "LA County", "Santa Clara + LA County", "New York State", NA, NA, 
                 "SC + LA + NY")
  N = N_values[use_data_no]
  n_obs = n_obs_values[use_data_no] # observed positives in study
  if(is.na(N) | is.na(n_obs)) stop("Invalid USE_DATA.")
  # Main study statistics.
  MainStudy = list(total=N, obs=n_obs)
  
  #' Parameter space.
  Theta = expand.grid(FPR=seq(0, 0.05, length.out=101), 
                      TPR=seq(0.6, 1, length.out=61), 
                      num_infect=seq(0, max_infect[use_data_no]))

  # Return MODEL object.
  return(list(CalibrationStudy=Calibration,
              MainStudy=MainStudy,
              Theta=Theta,
              name=names[use_data_no], 
              full_name=full_names[use_data_no],
              sc_version=SC_version, 
              use_Sharp=use_Sharp))
}


#' Main function. This performs one single test for specific values of the model parameters.
#' Implements the tests in Equation (8) and Equation (10) of paper.
#'
#' theta = (FPR, TPR, num_infect)
#' @FPR False positive rate: Pr(test=1 | state=0)
#' @TPR True positive rate: Pr(test=1 | state=1)
#' @num_infect Number of true infected individuals in main study.
#' @model Object model that contains information about calibration/main studies, and Theta (param space).
include_theta_Exact = function(theta, model, 
                               alpha_level=0.05, 
                               visualize=FALSE, dens_cutoff=100) {
  FPR = theta[1]; TPR = theta[2]; num_infect = theta[3]
  # p = FPR, q= TPR.
  stopifnot(FPR <= 1 & FPR >= 0 & TPR <= 1 & TPR >= 0)
  # p = 0.03; q = 1; num_infect=335
  # Calibration study, total trials
  Nc_neg = model$CalibrationStudy$Neg$total # Nc- = true negatives in validation study (true negatives)
  Nc_pos = model$CalibrationStudy$Pos$total # Nc+ = true positives.
  
  Dc_neg = dbinom(seq(0, Nc_neg), size=Nc_neg, prob=FPR, log=T)  # P(Sc- = 0, 1, ...Nc-)
  Dc_pos = dbinom(seq(0, Nc_pos), size=Nc_pos, prob=TPR, log=T)  # P(Sc+ = 0, 1..., Nc+)
  
  Nm = model$MainStudy$total  # total trials in main study.
  d1t_all = dbinom(seq(0, Nm), size=num_infect, prob=TPR, log=T)  # P(Sm+ = 0, 1, ...N), true positives
  d1f_all = dbinom(seq(0, Nm), size=Nm - num_infect, prob=FPR, log=T) # P(Sm- = 0, 1...N), false positives.
  
  log_dens_CalibrationStudy = function(sc_neg, sc_pos) (Dc_neg[sc_neg + 1] + Dc_pos[sc_pos + 1])
  
  log_dens_MainStudy = function(s1) {
    # Note that S1 = S1+  +  S1-. Need convolution density.
    x = seq(0, s1)  # S1+ = x, S1- = s1-x
    d1t = d1t_all[x + 1]
    d1f = d1f_all[s1 - x + 1] 
    log(sum(exp(d1t + d1f)))
  }
  
  # 1. observed density
  S1_obs = model$MainStudy$obs
  S0n_obs = model$CalibrationStudy$Neg$obs
  S0p_obs = model$CalibrationStudy$Pos$obs
  f_obs = exp(log_dens_MainStudy(S1_obs) + log_dens_CalibrationStudy(S0n_obs, S0p_obs))
  
  # Calculate full density
  incl = NA; stren = NA; reduction = NA
  if(!model$use_Sharp) {
    d1 = sapply(seq(0, Nm), log_dens_MainStudy)  # log P(S1=0, 1, ...N)
    d0n = Dc_neg
    d0p = Dc_pos
    
    cutoff = -dens_cutoff  # runs things faster.
    d1 = d1[d1 > cutoff]
    d0n = d0n[d0n > cutoff]  # take only support, f>0
    d0p = d0p[d0p > cutoff]  #
    
    dens = expand.grid(d0n=d0n, d0p=d0p, d1=d1)
    
    dens$f = exp(dens$d0n + dens$d0p + dens$d1)  # Because of test independence.
    # 2. Calculate nu. (critical number.)
    nu =  nrow(subset(dens, f <= f_obs & f > 0))  # all those have P > 0 and P <= s.
    # 3. nu function for test decision.
    incl = ((nu * f_obs) > alpha_level)
    stren = nu*f_obs
  } else {
    #' 
    d1 = sapply(seq(0, Nm), log_dens_MainStudy)  # log P(S1=0, 1, ...N)
    d0n = Dc_neg
    d0p = Dc_pos
    
    # 
    d1_reduct = d1[d1 > log(f_obs)]
    d0n_reduct = d0n[d0n > log(f_obs)]  
    d0p_reduct = d0p[d0p > log(f_obs)]
    
    # How much is the search space reduced?
    C1 = length(d1_reduct) / length(d1)
    C2 = length(d0n_reduct) / length(d0n)
    C3 = length(d0p_reduct) / length(d0p)
    # print(paste("Reduction numbers ", C1, " ", C2, " ", C3))
    # print(d0p)
    
    reduction = C1*C2*C3
    
    dens = expand.grid(d0n=d0n_reduct, d0p=d0p_reduct, d1=d1_reduct)
    dens$f = exp(dens$d0n + dens$d0p + dens$d1)  # Because of test independence.
    #dens00 = subset(dens, f > 0 & f <= f_obs)
    #print(sum(dens00$f))#
    dens = subset(dens, f > f_obs)

    cutoff = max(0, 1 - sum(dens$f))
    incl = (cutoff > alpha_level)
    stren = cutoff
  }
  ###
  retObj = list(incl=incl, strength=stren, reduction_pct=100 - 100*reduction)
  
  
  if(visualize) {
    print("> Visualization ON.")
    if(model$use_Sharp) {
      print("> Using sharp construction.")
      print(sprintf("> Compare I(f<= fobs) * f = %f with %f", 
                    cutoff, alpha_level))
    } else {
      print("> Using classic construction.")
      print(sprintf("> Compare nu*f_obs = %f with %f", stren, alpha_level))
    }
    # Create density  again 
    N = model$MainStudy$total; 
    M = model$CalibrationStudy$Neg$total; 
    Mpos = model$CalibrationStudy$Pos$obs
    m = seq(0, M)
    n = seq(0, N)
    d1 = sapply(n, log_dens_MainStudy)
    d0 = sapply(m, function(i) log_dens_CalibrationStudy(i, Mpos))
    
    #d0 = d0n_all
    dens = expand.grid(d0=d0, d1=d1)
    arg = expand.grid(s0=m, s1=n)
    dens = cbind(arg, dens)
    dens$f = exp(dens$d0 + dens$d1)
    sub = subset(dens, f >= .0001*max(f))  # plot only around the mode.
    nu2 = nrow(subset(dens, f <= f_obs & f > 0))
    # nrow(subset(dens, f <= f_obs & f > 0))
    print(sprintf("Density at observed point, f_obs = %.10f", f_obs))
    print(sprintf("Total %d points of support with <= f_obs", nu2))
    
    m_obs = model$CalibrationStudy$Neg$obs
    n_obs = model$MainStudy$obs
    # Revise to make plot better.
    m = seq(min(c(sub$s0,  m_obs)), max(c(sub$s0, m_obs)))
    n = seq(min(c(sub$s1,  n_obs)), max(c(sub$s1, n_obs)))
    arg = expand.grid(s0=m, s1=n)
    d1 = sapply(n, function(i) log_dens_MainStudy(i))
    d0 = sapply(m, function(i) log_dens_CalibrationStudy(i, Mpos))
    viz_dens = expand.grid(d0=d0, d1=d1)
    viz_dens = cbind(arg, viz_dens)
    viz_dens$f = exp(viz_dens$d0 + viz_dens$d1)
    
    # viz_dens = (s0, s1, d0, d1, f)
    
    #out = expand.grid(s0=m, s1=n
    z = matrix(viz_dens$f, nrow=length(m), ncol=length(n))
    require(plot3D)
    par(mfrow=c(1, 1))
    
    main_vals =sprintf("(%d,%d,%d)", model$CalibratioStudy$Neg$obs, 
                       model$CalibrationStudy$Pos$obs, model$MainStudy$obs)
    smain = expression(paste("(", s[c][obs]^"-",", ", s[c][obs]^"+", ", ", s[m][obs]^"", ") = (2, 178, 50)"))
    if(model$sc_version==2) {
      smain = expression(paste("(", s[c][obs]^"-",", ", s[c][obs]^"+", ", ", s[m][obs]^"", ") = (16, 130, 50)"))
    }
    sx = expression(paste("FP in calibration study ", S[c]))
    persp3D(m, n, z, border = "black", axes=TRUE, ticktype="detailed",
            xlab="false positives in calibration study (Sc-)", 
            ylab="positives in main study (Sm)", theta=-50, phi=25, 
            expand=0.5, zlab="density", cex.main=1.5,
            cex.lab=1.3, main=smain)
    
    # text3D(10, 0, 0, rot=c(1, 1, 1, 0.5), labels = sx, add = TRUE, adj=1)
    
    # sprintf("%s study. *Data: (s0n, s1)=(%d, %d)", model$full_name, m_obs, n_obs))
    #points3D(m, m, m)
    m_obs = m_obs; n_obs = n_obs
    d1 = subset(dens, s0==m_obs & s1==n_obs)$f
    # print(d1)
    stopifnot(length(d1) > 0)
    scatter3D(m_obs, n_obs, 1.2*d1, col="red", add=T, pch="*", cex=5)
    
  }
  
  return(retObj)
}



Sample_Data = function(model, theta_0, nsampl=1000) {
  FPR = theta_0[1]; TPR = theta_0[2]; num_infect = theta_0[3]
  stopifnot(FPR <= 1 & FPR >= 0 & TPR <= 1 & TPR >= 0)
  
  Nc_neg = model$CalibrationStudy$Neg$total # Nc- = true negatives in validation study (true negatives)
  Nc_pos = model$CalibrationStudy$Pos$total # Nc+ = true positives.
  Nm = model$MainStudy$total
  
  neg_obs = rbinom(nsampl, size=Nc_neg, prob=FPR)
  pos_obs = rbinom(nsampl, size=Nc_pos, prob=TPR)
  obs = rbinom(nsampl, size=num_infect, prob=TPR) + rbinom(nsampl, size=Nm - num_infect, prob=FPR)
  
  all_samples = list()
  for(i in 1:nsampl) {
    m = model
    m$CalibrationStudy$Neg$obs = neg_obs[i]
    m$CalibrationStudy$Pos$obs = pos_obs[i]
    m$MainStudy$obs = obs[i]
    all_samples[[i]] = m
  }
  return(all_samples)
}


#' Log likelighood for given model and parameters.
#' @param model It is the data model. See covid19_serology.R
#' @param theta_0 Parameter vector (FPR, TPR, #infected)
#' 
loglik = function(model, theta_0) {
  FPR = theta_0[1]; TPR = theta_0[2]; num_infect = theta_0[3]
  GOOD = (FPR <= 1 & FPR >= 0 & TPR <= 1 & TPR >= 0)
  if(!GOOD) {
    print(theta_0)
    stop("")
  }
  Nc_neg = model$CalibrationStudy$Neg$total # Nc- = true negatives in validation study (true negatives)
  obs_neg = model$CalibrationStudy$Neg$obs # observed positives in Neg-study
  
  Nc_pos = model$CalibrationStudy$Pos$total # Nc+ = true positives.
  obs_pos = model$CalibrationStudy$Pos$obs # observed positives in Pos-study
  
  log_dens_CalibrationStudy = dbinom(obs_neg, size=Nc_neg, prob=FPR, log=T) + 
    dbinom(obs_pos, size=Nc_pos, prob=TPR, log=T)
  
  
  Nm = model$MainStudy$total  # total trials in main study.
  d1t_all = dbinom(seq(0, Nm), size=num_infect, prob=TPR, log=T)  # P(Sm+ = 0, 1, ...N), true positives
  d1f_all = dbinom(seq(0, Nm), size=Nm - num_infect, prob=FPR, log=T) # P(Sm- = 0, 1...N), false positives.
  
  obs = model$MainStudy$obs # obs positives in main study.
  x = seq(0, obs)  # S1+ = x, S1- = s1-x
  d1t = d1t_all[x + 1]
  d1f = d1f_all[obs - x + 1] 
  log_dens_MainStudy = log(sum(exp(d1t + d1f)))
  
  f_obs = log_dens_MainStudy + log_dens_CalibrationStudy
  # f_obs = log_dens_MainStudy
  # warning("only main study")
  return(f_obs)
}


# Max-likelihood routines
expit = function(x) {
  if(x > 50) return(1)
  if(x < -50) return(0)
  exp(x) / (1+exp(x))
}

logit = function(z) {
  stopifnot(z >= 0 & z <= 1)
  e = 1e-6
  z = max(e, min(z, 1-e))
  return(log(z / (1-z)))
}

#' Transfers from natural parameterization to theta-parameterization.
theta_param = function(natural_param, Nm) {
  b0 = expit(natural_param[1])
  b1 = expit(natural_param[2])
  b2 = as.integer(expit(natural_param[3]) * Nm)
  
  c(b0, b1, b2)
}

#' Maps from theta-parameterization to natural one.
natural_param = function(theta_param, Nm) {
  fpr = theta_param[1]; tpr = theta_param[2]; num_infect = theta_param[3]
  q = (num_infect/Nm)
  
  par0 = c(logit(fpr), logit(tpr), logit(q))
  return(par0)
}

#' Maximum likelihood for this particular model.
#' @param model Data model to use (see covid19_serology.R for model definitions.)
max_loglik = function(model) {
  
  Nm = model$MainStudy$total
  
  obj = function(par) {
    # b0 = expit(par[1])
    # b1 = expit(par[2])
    # b2 = as.integer(expit(par[3]) * model$MainStudy$total)
    theta = theta_param(par, Nm)
    if(any(is.na(theta))) {
      print(theta)
      model$Theta = NULL
      print(model)
      print(par)
      stop("ERROR")
    }
    ret = -loglik(model, theta)
    
    if(is.infinite(ret)) { return(1e20)}
    return(ret)
  }
  
  fpr_hat = model$CalibrationStudy$Neg$obs/  model$CalibrationStudy$Neg$total
  tpr_hat = model$CalibrationStudy$Pos$obs/  model$CalibrationStudy$Pos$total
  num_infect_hat = model$MainStudy$obs
  
  theta_hat = c(fpr_hat, tpr_hat, num_infect_hat)
  par0 = natural_param(theta_hat, Nm)
  # print(par0)
  out = optim(par=par0, fn = obj, method="BFGS", control=list(maxit=100))
  # print("done")
  theta_best = theta_param(out$par, Nm)
  return(list(theta=theta_best, val=-out$val))
}

#' Likelihood ratio test.
#' @param model It is the data model. See covid19_serology.R
#' @param theta_0 Parameter vector (FPR, TPR, #infected)
#' 
LR_test = function(theta_0, model, nsampl=1000, vis=F) {
  # th = apply(model$Theta, 1, function(row) all(row==theta_0))
  # stopifnot(sum(th)==1)
  # Lmax = max(apply(model$Theta[-which(th), ], 1, function(row) Likelihood(model, row)))
  best = max_loglik(model)
  
  Lmax = best$val
  L0 = loglik(model, theta_0)
  
  if(is.infinite(L0)) {
    return(list(incl=FALSE, pvalue=0))
  }
  Tobs = exp(L0  - Lmax)
  
  s = Sample_Data(model, theta_0, nsampl = nsampl)
  tvals = sapply(1:length(s), function(i) {
    m = s[[i]]
    Lmax = max_loglik(m)$val
    exp(loglik(m, theta_0) - Lmax)
  })
  
  if(vis) {
    print(L0)
    print(tvals)
    hist(tvals, breaks=30, col="blue")
    abline(v=Tobs, col="red", lwd=2)
  }
  m = mean(Tobs >= tvals)
  pval = min(m, 1-m)
  return(list(incl=(pval >= 0.025), pvalue=pval))
}



Proc1_chen = function(model = create_Model(1, TRUE), num_mcmc=200000, sigma_mcmc = 0.5) {
  #' Implements Procedure 1 from https://arxiv.org/pdf/1605.00499.pdf
  #' Chen et. al. (2018)
  #' 
  Nm = model$MainStudy$total
  # where to start from?
  # out = max_loglik(model)
  # NEW.
  # num_infect0 = 40
  # PQ = expand.grid(fpr=seq(0, 0.05, length.out=200), tpr=seq(0.7, 1, length.out=100))  
  # z = apply(PQ, 1, function(row) {
  #   th = c(as.numeric(row), num_infect0)
  #   loglik(model, th)
  # })
  theta0 = sample_MCMC_theta()
  print("> Theta-hat initialized.")
  # chain = matrix(theta_hat, nrow=3, ncol=1)
  chain = matrix(theta0, nrow=3, ncol=1)
  # par_hat = natural_param(theta_hat, Nm)
  
  num_acc = 0
  t0 = proc.time()[3]
  print(" > Running MCMC")
  for(j in 1:num_mcmc) {
    theta_old = chain[, ncol(chain)]
    # par_old = natural_param(theta_old, Nm)
    # par_new = rnorm(3, mean=par_old, sd=sigma_mcmc)
    # theta_new = theta_param(par_new, Nm)
    theta_new = sample_MCMC_theta()
    acc = min(1, exp(loglik(model, theta_new) - loglik(model, theta_old)))
    if(runif(1) < acc) {
      # print("accept")
      chain = cbind(chain, theta_new)
      num_acc = num_acc + 1
    }
    t1 = proc.time()[3]
    if(t1 - t0 > 10) {
      print(sprintf("j=%d/%d -- Acceptance ratio = %.2f%%", j, num_mcmc, 100*num_acc/j))
      t0 = t1
      plot(as.mcmc(t(chain))) # plot the chain
    }
  }
  print(sprintf("> MCMC Done. Acceptance ratio = %.2f%%", 100*num_acc/num_mcmc))
  
  n = ncol(chain)
  chain = chain[, tail(1:n, as.integer(0.8*n))]  # cut first 20%
  chain = t(chain)
  colnames(chain) = c("FPR", "TPR", "num_infect")
  rownames(chain) = NULL
  require(coda)
  
  plot(as.mcmc(chain)) # plot the chain
  
  100 *quantile(chain[,3] / model$MainStudy$total, probs=c(0.025, 0.975))
  ## Build CI
  lvals = sapply(1:nrow(chain), function(j) {
    theta = as.numeric(chain[j,])
    loglik(model, theta)
  })
  
  Q = as.numeric(quantile(lvals, probs = 0.95))
  t0 = proc.time()[3]
  incl = rep(FALSE, nrow(model$Theta))
  for(i in 1:nrow(model$Theta)) {
    
    theta = as.numeric(model$Theta[i,])
    # print(theta)
    ll0 = loglik(model, theta)
    # print(paste("ll=", ll0, " incl? ", ll0 > Q))
    incl[i] = as.logical(ll0  >= Q)
    t1 = proc.time()[3]
    if(t1 - t0 > 10) {
      t0=t1
      print(paste(" Done ", 100*i/nrow(model$Theta), " %."))
      # print(model$Theta[which(incl),])
    }
    
  }
  
  Results = model$Theta[which(incl),]
  save(Results, file="Prevalence_SC_MCMC.rdata")
  save(chain, file="Prevalence_SC_MCMC_chain.rdata")
  plot_Results_ribbon2(Results, model)
}

Proc3_chen = function(model) {
  # grid for (FPR, TPR)
  PQ = expand.grid(fpr=seq(0, 0.05, length.out=250), tpr=seq(0.7, 1, length.out=100))  
  
  # 1. Find some good \theta^
  num_infect0 = 40
  z = apply(PQ, 1, function(row) {
    th = c(as.numeric(row), num_infect0)
    loglik(model, th)
  })
  
  theta_hat = c(as.numeric(PQ[which.max(z),]), num_infect0)
  Lhat = loglik(model, theta_hat)
  print(paste("Lhat = ", Lhat))
  # 2. grid for prevalence + Q function.
  # n = model$CalibrationStudy$Neg$total + model$CalibrationStudy$Pos$total + model$MainStudy$total
  Q = function(theta) {
    2*(Lhat - loglik(model, theta))
  }
  
  # Main iteration
  n_infect = seq(0, 60)
  incl = rep(0, length(n_infect))
  t0 = proc.time()[3]
  for(i in 1:length(n_infect)) {
    # preval = round(n_infect/model$MainStudy$total)
    num_infect = n_infect[i] # as.integer(preval * model$MainStudy$total)
    y = min(apply(PQ, 1, function(row) {
      th = c(as.numeric(row), num_infect)
      Q(th)
    }))
    
    # 3. Decision.
    incl[i] = (y <= qchisq(0.95, df=1))
    t1 = proc.time()[3]
    if(t1 - t0 > 5) {
      ci = n_infect[which(incl==1)] / model$MainStudy$total
      
      print(sum(incl))
      print(paste("i = ", i, "/", length(n_infect), " num_infect=", num_infect, " ", paste(min(ci), ", ", max(ci))))
      t0 = t1
    }
  }
}

sample_MCMC_theta = function() {
  fpr = rbeta(1, shape1 = 6, shape2 = 1000)
  tpr = rbeta(1, shape1 = 85, shape2 = 10)
  all = 0:80
  num_infect = sample(all, size=1, replace=T)
  #prob=sapply(all, function(i) min(i, 80-i)))
  num_infect = sample(all, size=1, replace=T, prob=(80-all)^2)
  
  return(c(fpr, tpr, num_infect))
}

