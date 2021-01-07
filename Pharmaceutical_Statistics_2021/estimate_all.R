library(rstan)
library(dplyr)

S <- 3
K <- 15
L <- 2
Hyp1 <- 0.5
Hyp2 <- 1.5

d <- read.csv('input/data.csv', stringsAsFactors=FALSE)
d_ptxd <- read.csv('input/data_ptxd.csv', stringsAsFactors=FALSE)

## prepare data.frame for data matrix
d <- d %>%
  tidyr::unite(col=AECOS_GRAD, c(AECOS, GRAD), sep='_', remove=FALSE) %>%
  filter(!AECOS %in% c('ALOPECIA'))
subj_str <- as.character(sort(unique(d$SUBJID)))
conv_subj <- setNames(seq_along(subj_str), subj_str)
ae_str <- sort(unique(d$AECOS_GRAD))
conv_ae <- setNames(seq_along(ae_str), ae_str)
E <- length(conv_ae)  # determination of E
d <- d %>%
  mutate(SUBJID=conv_subj[as.character(SUBJID)], AECOS_GRAD_N=conv_ae[AECOS_GRAD])  # renumbering

## preparation for severities
d_ae <- data.frame(num=seq_along(conv_ae), ae=gsub('_[1-4]', '', names(conv_ae)))
d_ae_count <- d_ae %>%
  group_by(ae) %>% summarize(start=min(num), end=max(num)) %>%
  arrange(start) %>% ungroup()
GX <- d_ae_count %>% filter(start != end) %>%
  select(start, end) %>% as.matrix()
E_GX <- nrow(GX)
GX0 <- d_ae_count %>% filter(start == end) %>% .$start %>% c(GX[,1]) %>% sort()
E_GX0 <- length(GX0)

## prepare data.frame for the number of patients
d_ptxd <- d_ptxd %>% filter(SUBJID %in% subj_str) %>%
  mutate(SUBJID=conv_subj[as.character(SUBJID)])
subjs <- lapply(1:10, function(trtid) {
  d_ptxd %>%
    filter(TRTID == trtid) %>% arrange(SUBJID) %>% .$SUBJID %>% conv_subj[.]
})

## calculate the numbers of patients
subjs[[1]] <- subjs[[2]]
subjs[[5]] <- subjs[[6]]
subjs[[8]] <- subjs[[9]]
N <- sapply(subjs, length)
N_A1     <- N[1]
N_A1_AC  <- N[2]
N_A1_T   <- N[3]
N_A2_FAC <- N[4]
N_A3     <- N[5]
N_A3_A   <- N[6]
N_A3_CMF <- N[7]
N_A4     <- N[8]
N_A4_AC  <- N[9]
N_A4_CMF <- N[10]
Nx_A1_Re <- which(subjs[[1]] %in% subjs[[3]])
Nx_A3_Re <- which(subjs[[5]] %in% subjs[[7]])
Nx_A4_Re <- which(subjs[[8]] %in% subjs[[10]])

## construct AE data matrix
Y_all <- reshape2::acast(d, SUBJID ~ AECOS_GRAD_N ~ TRTID)
Y_A1_S1  <- t(Y_all[subjs[[1]],,1])
Y_A1_AC  <- t(Y_all[subjs[[2]],,2])
Y_A1_T   <- t(Y_all[subjs[[3]],,3])
Y_A2_FAC <- t(Y_all[subjs[[4]],,4])
Y_A3_S3  <- t(Y_all[subjs[[5]],,5])
Y_A3_A   <- t(Y_all[subjs[[6]],,6])
Y_A3_CMF <- t(Y_all[subjs[[7]],,7])
Y_A4_S3  <- t(Y_all[subjs[[8]],,8])
Y_A4_AC  <- t(Y_all[subjs[[9]],,9])
Y_A4_CMF <- t(Y_all[subjs[[10]],,10])

## calculate the numbers of cycles
Y_Cyc <- reshape2::acast(d_ptxd, SUBJID ~ TRTID, value.var='cycles')
Cyc_A1_AC  <- Y_Cyc[subjs[[2]],1]
Cyc_A1_T   <- Y_Cyc[subjs[[3]],2]
Cyc_A2_FAC <- Y_Cyc[subjs[[4]],3]
Cyc_A3_A   <- Y_Cyc[subjs[[6]],4]
Cyc_A3_CMF <- Y_Cyc[subjs[[7]],5]
Cyc_A4_AC  <- Y_Cyc[subjs[[9]],6]
Cyc_A4_CMF <- Y_Cyc[subjs[[10]],7]

## estimation
stanmodel <- stan_model(file='model/model_all.stan')
fit <- optimizing(stanmodel, seed=7, iter=100000, verbose=TRUE, as_vector=FALSE)
save(fit, file='output/result.RData')
