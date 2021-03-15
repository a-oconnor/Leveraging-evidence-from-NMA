source(file = "./functions.R")

BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)

p_nac <- BRD_long %>% filter(t == "No active control") %>% summarise(p = sum(r)/sum(n)) %>% pull
lor_nac_2_enro <- nma_old$TE.fixed[8, 5]
m <- p_nac/(1-p_nac) * exp(-lor_nac_2_enro)
p_enro <- m/(1+m)

cl <- makeCluster(16)
registerDoParallel(cl)
set.seed(20200921)
tmp <- foreach(i = 1:200) %dorng%{
  bio_equal(p_neg = p_nac, p_pos = p_enro,  p_new1 = p_enro, sigma = 0, re_sigma = 0, data_prev = BRD_long,n_iter = 50, pig_alloc = c(800, 800, 800))
}
stopCluster(cl)

BE <- do.call("rbind", tmp)

BE <- as.data.frame(BE)
names(BE) <- c("superiority enro to nac", "superiority new to nac", "non inferiority")

write.csv(BE, file = "./res/alloc_2400_800_800_800.csv", row.names = F)



