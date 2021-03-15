library(DEoptimR)

BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)

sigma2_prev <- nma_old$seTE.fixed[8,5]^2

var_fn <- function(n, beta1, beta2, beta3, sigma2_prev, s){
  mu_1 <- exp(beta1)/(1+exp(beta1))^2
  mu_2 <- exp(beta1 + beta2)/(1+exp(beta1 + beta2))^2
  mu_3 <- exp(beta1 + beta3)/(1+exp(beta1 + beta3))^2
  1/(mu_2 * n[2]) + 1/(mu_3 * n[3]) - 1/(mu_2^2 * n[2]^2 *(sigma2_prev + 1/(mu_1 * n[1]) + 1/(mu_2 * n[2])))
}

# var_fn_neg <- function(n, beta1, beta2, beta3, sigma2_prev, s){
#   mu_1 <- exp(beta1)/(1+exp(beta1))^2
#   mu_2 <- exp(beta1 + beta2)/(1+exp(beta1 + beta2))^2
#   mu_3 <- exp(beta1 + beta3)/(1+exp(beta1 + beta3))^2
#   #((1/(mu_1 *n[1]) + 1/(mu_2 * n[2])) * (1/(mu_2 * n[2]) + 1/(mu_3 * n[3])) - 1/(mu_2 * n[2])^2 +
#   #(1/(mu_1 *n[1]) -  1/(mu_3 * n[3])) * sigma2_prev)/(sigma2_prev + 1/(mu_1 * n[1]) + 1/(mu_2 * n[2]))
#   1/(mu_1 * n[1]) + 1/(mu_3 * n[3]) - 1/(mu_1^2 * n[1]^2 *(sigma2_prev + 1/(mu_1 * n[1]) + 1/(mu_2 * n[2])))
# }

con_fn <- function(n, beta1, beta2, beta3, sigma2_prev, s){
  n1 <- n[1]
  n2 <- n[2]
  n3 <- n[3]
  n1 + n2 + n3 - s
}

p_nac <- BRD_long %>% filter(t == "No active control") %>% summarise(p = sum(r)/sum(n)) %>% pull
lor_nac_2_enro <- nma_old$TE.fixed[8, 5]
m <- p_nac/(1-p_nac) * exp(-lor_nac_2_enro)
p_enro <- m/(1+m)

# optimal pig allocation, s is the number of pigs
s = 3600
tmp <- JDEoptim(lower = c(0,0,0), upper = c(s,s,s), fn = var_fn, constr = con_fn, meq = 1,trace = T, triter = 50, beta1 = log(p_nac/(1-p_nac)), beta2 = log(p_enro/(1-p_enro)) - log(p_nac/(1-p_nac)), beta3 = log(p_enro/(1-p_enro)) - log(p_nac/(1-p_nac)), sigma2_prev = sigma2_prev, s = s)

# sqrt(var_fn(c(23, 1, 26), beta1 = log(p_nac/(1-p_nac)), beta2 = log(p_enro/(1-p_enro)) - log(p_nac/(1-p_nac)), beta3 = log(p_enro/(1-p_enro)) - log(p_nac/(1-p_nac)), sigma2_prev = sigma2_prev, s = s))
# 
# sqrt(var_fn(c(0, 25, 25), beta1 = log(p_nac/(1-p_nac)), beta2 = log(p_enro/(1-p_enro)) - log(p_nac/(1-p_nac)), beta3 = log(p_enro/(1-p_enro)) - log(p_nac/(1-p_nac)), sigma2_prev = sigma2_prev, s = s))
