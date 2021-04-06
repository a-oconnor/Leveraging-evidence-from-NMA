packages <- c("MASS","NlcOptim","netmeta","dplyr")
lapply(packages, library, character.only = TRUE)

SolveSampleSize_Withprev <- function(p1, p2, p3, sigma,power_level){
  beta1 = log(p1/(1-p1))
  beta2 = log(p2/(1-p2)) - log(p1/(1-p1))
  beta3 = log(p3/(1-p3)) - log(p1/(1-p1))
  sigma_prev <- sigma
  n0=c(100,100,100)
  
  power_withprev <- function(n){
    var <- 1/(p2 * n[2]) + 1/(p3 * n[3]) - 1/(p2^2 * n[2]^2 *(sigma_prev^2 + 1/(p1 * n[1]) + 1/(p2 * n[2])))
    se <- sqrt(var)
    z <- (beta3-beta2)/se
    power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
    return(power)
  }
  
  confun_withprev <- function(n){
    f = power_level-power_withprev(n)
    f = rbind(f,-n[1])
    f = rbind(f,-n[2])
    f = rbind(f, -n[3])
    return(list(ceq=NULL,c=f))
  }
  
  objfun=function(n){
    n[1]+n[2]+n[3]
  }
  
  solution_temp <- solnl(n0,objfun=objfun,confun=confun_withprev)$par
  solution_temp_int <- round(solution_temp,0)
  # get the integer solution around this
  dat_para <- expand.grid(n1=c(max(solution_temp_int[1]-2,1):(solution_temp_int[1]+2)),n2=c(max(solution_temp_int[2]-2,1):(solution_temp_int[2]+2)),
                          n3=c(max(solution_temp_int[3]-2,1):(solution_temp_int[3]+2)))
  
  for(c in 1:nrow(dat_para)){
    dat_para[c,4] <- power_withprev(c(dat_para$n1[c],dat_para$n2[c],dat_para$n3[c]))
  }
  
  dat_para <- dat_para[dat_para$V4>=power_level,]
  dat_para$n <- dat_para$n1+dat_para$n2 + dat_para$n3
  nmin <- min(dat_para$n)
  dat_para <- dat_para[dat_para$n==nmin,]
  dat_para <- dat_para[order(dat_para$V4,decreasing = T),]
  return(as.numeric(dat_para[1,1:3]))
}




SolveSampleSize_Withprev_equal <- function(p1, p2, p3, sigma,power_level){
  # cal the total sample size when we added one more condition: sample sizes are equal for both treatment groups
  beta1 = log(p1/(1-p1))
  beta2 = log(p2/(1-p2)) - log(p1/(1-p1))
  beta3 = log(p3/(1-p3)) - log(p1/(1-p1))
  sigma_prev <- sigma
  n0=c(100, 100, 100)
  
  power_withprev <- function(n){
    size <- n/3
    var <- 1/(p2 * size) + 1/(p3 * size) - 1/(p2^2  * size^2 *(sigma_prev^2 + 1/(p1 * size) + 1/(p2 * size)))
    se <- sqrt(var)
    z <- (beta3-beta2)/se
    power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
    return(power)
  }
  
  confun_withprev <- function(n){
    f = power_level-power_withprev(n)
    f = rbind(f,-n)
    return(list(ceq=NULL,c=f))
  }
  
  objfun=function(n){
    return(n)
  }
  
  solution_temp <- solnl(n0,objfun=objfun,confun=confun_withprev)$par
  solution_temp_int <- round(solution_temp[1],0)
  
  if(solution_temp_int %% 3 ==0){
    res <- solution_temp_int
  }else if(solution_temp_int %% 3 == 1){
    res <- solution_temp_int + 2
  } else {
    res <- solution_temp_int + 1
  } 
  
  return(res)
}





SolveSampleSize_Single_equal <- function(p1, p2, p3,power_level){
  beta1 = log(p1/(1-p1))
  beta2 = log(p2/(1-p2)) - log(p1/(1-p1))
  beta3 = log(p3/(1-p3)) - log(p1/(1-p1))
  n0=c(100,100,100)
  
  power_single <- function(n){
    size <- n/3
    var <- 1/(p2 * size) + 1/(p3 * size)
    se <- sqrt(var)
    z <- (beta3-beta2)/se
    power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
    return(power)
  }
  
  confun_single <- function(n){
    f = power_level-power_single(n)
    f = rbind(f,-n)
    return(list(ceq=NULL,c=f))
  }
  
  objfun=function(n){
    return(n)
  }
  
  solution_temp <- solnl(n0,objfun=objfun,confun=confun_single)$par
  solution_temp_int <- round(solution_temp[1],0)
  
  
  if(solution_temp_int %% 3 ==0){
    res <- solution_temp_int
  }else if(solution_temp_int %% 3 == 1){
    res <- solution_temp_int + 2
  } else {
    res <- solution_temp_int + 1
  } 
  return(res)
}


# rearrange the data to the long type: one arm one study per row
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

# log(k1/k2)=lor
lor2prob <- function(p1, lor){
  p2 <- p1/(p1 + exp(lor)*(1-p1))
  return(p2)
}
