packages <- c("dplyr", "netmeta", "parallel", "foreach", "doParallel", "doRNG")
lapply(packages, library, character.only = TRUE)

lor2prob <- function(p1, lor){
  p2 <- p1/(p1 + exp(lor)*(1-p1))
  return(p2)
}

wide2long <- function(MTCdata){
  N <- max(MTCdata$Number.of.arms)
  study <- rep(MTCdata$Study.number,N)
  t <- NULL
  n <- NULL
  r <- NULL
  for(i in c(2,1,3)){
    r <- c(r, eval(parse(text = paste0("MTCdata$Number.of.Event.in.arm.",i, sep = ""))))
    n <- c(n, eval(parse(text = paste0("MTCdata$Total.number.in.arm.",i, sep = ""))))
    t <- c(t, eval(parse(text = paste0("MTCdata$Arm.",i, sep = ""))))
  }
  res <- data.frame(id = study, t = t, r = r, n = n)
  res <- res %>% dplyr::filter(!is.na(n)) %>% arrange(id)
  res
}




new_study <- function(p_neg, p_pos, p_new1, sigma, pig_alloc = c(50,50,50)){
  r_vec <- rbinom(3, size = pig_alloc, prob = c(p_neg, p_pos, p_new1))
  n_vec <- pig_alloc
  new_s <- data.frame(id = rep(100,3), t = c("No active control", "Enrofloxacin", "New Treatment"),
                      r = r_vec, n = n_vec)
  return(new_s)
}

###non-inferiority and superiority with previous network
bio_equal <- function(p_neg, p_pos, p_new1, sigma, re_sigma = 0, pig_alloc = c(50, 50, 50), data_prev, n_iter){
  res_TE <- numeric(n_iter)
  res_seTE <- numeric(n_iter)
  power <- numeric(n_iter)
  sup_power_enro <- numeric(n_iter)
  sup_power_new <- numeric(n_iter)
    for(i in 1:n_iter){
      new_s <- new_study(p_neg = p_neg, p_pos = p_pos, p_new1 = p_new1, sigma = sigma, pig_alloc = pig_alloc)
      data_final <- rbind(data_prev, new_s)
      BRD_new <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = data_final, allstudies = T, sm = "OR")
      nma_res <- netmeta::netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new,sm="OR",comb.fixed =F,comb.random = T, tol.multiarm.se = 0.01)
      test_stat <- abs(nma_res$TE.fixed[8,5])/nma_res$seTE.fixed[8,5]
      power[i] <- nma_res$TE.fixed[8,5] + qnorm(0.95) * nma_res$seTE.fixed[8,5] < 0.2
      sup_power_enro[i] <- abs(nma_res$TE.fixed[5, 9])/nma_res$seTE.fixed[5,9] > qnorm(0.975)
      sup_power_new[i] <- abs(nma_res$TE.fixed[8, 9])/nma_res$seTE.fixed[8,9] > qnorm(0.975)
    }
  res <- matrix(c(sup_power_enro, sup_power_new, power), ncol = 3)
  return(res)
}

###non-inferiority and superiority without previous network
bioeq_single_s <- function(pig_alloc = c(50, 50, 50), p_nac, p_enro, p_new1 ,n_iter){
  power <- numeric(n_iter)
  sup_power_enro <- numeric(n_iter)
  sup_power_new <- numeric(n_iter)
  for(i in 1:n_iter){
    BRD_s <- data.frame(id = rep(100,3), t = c( "Enrofloxacin", "New Treatment", "No active control"), r = rbinom(3, pig_alloc, c(p_enro, p_new1, p_nac)), n = pig_alloc)
    BRD_new_s <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_s, allstudies = T, sm = "OR")
    nma_s <- netmeta::netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_s,sm="OR",comb.fixed =F,comb.random = T, tol.multiarm.se = 0.01)
    power[i] <- nma_s$TE.fixed[2,1] + qnorm(0.95) * nma_s$seTE.fixed[2,1] < 0.2
    sup_power_enro[i] <- abs(nma_s$TE.fixed[1, 3])/nma_s$seTE.fixed[1,3] > qnorm(0.975)
    sup_power_new[i] <- abs(nma_s$TE.fixed[2, 3])/nma_s$seTE.fixed[2,3] > qnorm(0.975)
  }
  res <- matrix(c(sup_power_enro, sup_power_new, power), ncol = 3)
  return(res)
}

