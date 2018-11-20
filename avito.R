set.seed(25)

library(tidyverse)
library(lubridate)
library(Rcpp)
library(BH)
library(bbmle)

U <- read_tsv('UserInfo.tsv')
U <- U %>% mutate(fingerprint = interaction(UserDeviceID, UserAgentFamilyID))
U_sample <- U %>%
  sample_n(5000)

S <- read_tsv('SearchInfo.tsv')


############################
# fit model parameters to data excluding test sample
loc_props <- S %>% 
  filter(!(UserID %in% U_sample$UserID)) %>% 
  group_by(LocationID) %>% 
  summarise(n=n()) %>% ungroup
loc_prior <- loc_props$n + 1550
loc_prior <- loc_prior / sum(loc_prior)

fit_nb <- S %>% 
  filter(!(UserID %in% U_sample$UserID)) %>% 
  group_by(UserID) %>% 
  summarise(n=n()-1) %>% 
  .$n %>% sample(., 10000) %>% MASS::fitdistr(densfun = 'negative binomial')

exp_tables2 <- function(alpha, theta, nn) {
  et1_log <- lgamma(theta + nn + alpha) + lgamma(theta + 1) -
    log(alpha) - lgamma(theta + nn) - lgamma(theta + alpha)
  et <- exp(et1_log) - theta/alpha
  return(et)
}
fp_per_user <- U %>% filter(!(UserID %in% U_sample$UserID)) %>% .$fingerprint
N_seq <- round(10**(seq(log10(10), log10(100000), length.out = 100)))
ngroups_sim <- sapply(N_seq, function(x) length(unique(sample(fp_per_user, x, replace=FALSE))))
fit_crp <- nls(ngroups_sim ~ exp_tables2(alpha, -alpha + theta, N_seq), 
               algorithm = 'port', control = nls.control(maxiter = 100, warnOnly = TRUE),
               start = list(alpha = 0.5, theta = .5),
               lower = list(alpha = 0, theta = 0), upper = list(alpha = 1, theta = Inf))

############################

S_sample <- S %>%
  filter(UserID %in% U_sample$UserID)

rm(S, U); gc()

C <- read_tsv('Category.tsv')
L <- read_tsv('Location.tsv')

S_sample <- S_sample %>%
  left_join(U_sample, by='UserID') %>%
  left_join(C, by='CategoryID') %>%
  left_join(L, by='LocationID')

colnames(S_sample)[15] <- 'CatLevel'
colnames(S_sample)[18] <- 'LocLevel'

rm(C,L,U_sample); gc()

############################################################
############################################################

S_sample <- S_sample %>%
  arrange(fingerprint, SearchDate) %>%
  mutate(group = cumsum(ifelse(fingerprint != lag(fingerprint, default=TRUE), 1, 0)))
group <- as.integer(S_sample$group)
group_sizes <- S_sample %>%
  group_by(group) %>%
  summarise(n=n()) %>%
  .$n
group_start_inds <- as.integer(cumsum(c(0, group_sizes[1:(length(group_sizes)-1)])))
N_group <- S_sample %>% group_by(fingerprint) %>%
  summarise(n = n_distinct(UserID)) %>% .$n

S_sample <- S_sample %>% 
  mutate(mins = difftime(SearchDate, min(SearchDate), units="mins"), 
         rel_time = as.numeric(mins) / 
           (as.numeric(difftime(max(SearchDate), min(SearchDate), units="mins")))) %>% 
  arrange(fingerprint, rel_time)
## add jitter to avoid ties in SearchDate
rel_second = 1.0 / as.numeric(difftime(max(S_sample$SearchDate), min(S_sample$SearchDate), unit="secs"))
S_sample <- S_sample %>% mutate(rel_time = rel_time + runif(n(), 0, rel_second))
S_sample <- S_sample %>% arrange(fingerprint, rel_time)

time_ <- S_sample$rel_time

S_sample <- S_sample %>% 
  inner_join(
    loc_props %>% select(LocationID) %>% mutate(loc = 1:n())
  )
loc <- S_sample$loc
################################

sourceCpp('mcmcem.cpp', verbose=TRUE, rebuild=TRUE)

params <- list(nb_shape1 = coef(fit_nb)[[1]], nb_shape2 = coef(fit_nb)[[1]] / coef(fit_nb)[[2]],
               crp_discount = coef(fit_crp)[[1]], crp_A = -coef(fit_crp)[[1]] + coef(fit_crp)[[2]],
               dircat_discount = 0.32, rho = NULL)

fit_model1 <- mcmcem(250, 1e6, FALSE, FALSE, params,
                     group, group_sizes, group_start_inds, time_, loc, loc_prior, c(0L), c(0L), c(0L))
fit_model2 <- mcmcem(250, 1e6, FALSE, TRUE, params,
                     group, group_sizes, group_start_inds, time_, loc, loc_prior, c(0L), c(0L), c(0L))

fit_model3 <- mcmcem(250, 1e6, TRUE, FALSE, params,
                     group, group_sizes, group_start_inds, time_, loc, loc_prior, c(0L), c(0L), c(0L))


################################
# estimate performance metrics
sample_pair <- function(group_start_inds, group_sizes, group_weighted = FALSE) {
  if (group_weighted) {
    g <- sample.int(length(group_sizes), 1, prob = replace(group_sizes, which(group_sizes==1), 0))
    pair <- group_start_inds[g] + sample.int(group_sizes[g], 2)
  } else {
    done <- FALSE
    while (!done) {
      id1 <- sample.int(sum(group_sizes), 1)
      g <- max(which(group_start_inds < id1))
      if (group_sizes[g] >= 2)
        done = TRUE
    }
    pair <- group_start_inds[g] + sample.int(group_sizes[g], 2)
  }
  return(sort(pair))
}

check_pair <- function(x) {
  if (S_sample$UserID[x[1]] == S_sample$UserID[x[2]]) {
    return(1)
  } else {
    return(0)
  }
}

log_loss <- function(pred, obs) {
  return(-mean(obs * log(pmax(pmin(pred, 1-1e-6), 1e-6)) + (1-obs) * log(pmax(pmin(1-pred, 1-1e-6), 1e-6))))
}

set.seed(25)

scoreSample <- t(sapply(1:5000, function(i) sample_pair(group_start_inds, group_sizes)))
scoreSample <- rbind(scoreSample,
                     t(sapply(1:5000, function(i) sample_pair(group_start_inds, group_sizes, TRUE))))

obs <- apply(scoreSample, 1, check_pair)
# benchmarks
brier_score_naive <- mean((1.0 - obs[1:1000])^2)
brier_score_naive_weighted <- mean((1.0 - obs[1001:2000])^2)
logloss_naive <- log_loss(1.0, obs[1:1000])
logloss_naive_weighted <- log_loss(1.0, obs[1001:2000])

out_model1 <- mcmcem(1, 50*1e6, FALSE, FALSE,
                     params, group, group_sizes, group_start_inds, time_, loc, loc_prior, 
                     fit_model1$L1, scoreSample[,1], scoreSample[,2])
params2 <- params
params2$rho <- fit_model2$rho_vec[length(fit_model2$rho_vec)]
out_model2 <- mcmcem(1, 50*1e6, FALSE, TRUE, params2,
                     group, group_sizes, group_start_inds, time_, loc, loc_prior, 
                     fit_model2$L1, scoreSample[,1], scoreSample[,2])
out_model3 <- mcmcem(1, 50*1e6, TRUE, FALSE, params,
                     group, group_sizes, group_start_inds, time_, loc, loc_prior, 
                     fit_model3$L1, scoreSample[,1], scoreSample[,2])

brierscore_model1 <- mean((out_model1$ScoreSample[1:1000] - obs[1:1000])^2)
logloss_model1 <- log_loss(out_model1$ScoreSample[1:1000], obs[1:1000])
brierscore_model1_weighted <- mean((out_model1$ScoreSample[1001:2000] - obs[1001:2000])^2)
logloss_model1_weighted <- log_loss(out_model1$ScoreSample[1001:2000], obs[1001:2000])
brierscore_model2 <- mean((out_model2$ScoreSample[1:1000] - obs[1:1000])^2)
logloss_model2 <- log_loss(out_model2$ScoreSample[1:1000], obs[1:1000])
brierscore_model2_weighted <- mean((out_model2$ScoreSample[1001:2000] - obs[1001:2000])^2)
logloss_model2_weighted <- log_loss(out_model2$ScoreSample[1001:2000], obs[1001:2000])
brierscore_model3 <- mean((out_model3$ScoreSample[1:1000] - obs[1:1000])^2)
logloss_model3 <- log_loss(out_model3$ScoreSample[1:1000], obs[1:1000])
brierscore_model3_weighted <- mean((out_model3$ScoreSample[1001:2000] - obs[1001:2000])^2)
logloss_model3_weighted <- log_loss(out_model3$ScoreSample[1001:2000], obs[1001:2000])


########################
# make some plots
pdf("calibration-plots.pdf", width = 7.7, height = 3.0)
par(mfrow = c(1, 3), mar = c(1, 1, 2, 1), oma = c(4, 4, 0, 0), asp=1)
# par(pty="s")
plot(c(0,1), c(0,1), type = "n", xlim = c(0,1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i", xlab = NA, ylab = NA, axes = FALSE)
axis(1, labels = TRUE)
axis(2, labels = TRUE)
title("Model 1")
box()
abline(0, 1, lty = 2)
lines(lowess(out_model1$ScoreSample, obs, iter=0))
# par(pty="s")
plot(c(0,1), c(0,1), type = "n", xlim = c(0,1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i", xlab = NA, ylab = NA, axes = FALSE, title = "model1")
axis(1, labels = TRUE)
axis(2, labels = FALSE)
title("Model 2")
box()
abline(0, 1, lty = 2)
lines(lowess(out_model2$ScoreSample, obs, iter=0))
# par(pty="s")
plot(c(0,1), c(0,1), type = "n", xlim = c(0,1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i", xlab = NA, ylab = NA, axes = FALSE)
axis(1, labels = TRUE)
axis(2, labels = FALSE)
title("Model 3")
box()
abline(0, 1, lty = 2)
lines(lowess(out_model3$ScoreSample, obs, iter=0))
mtext(text= "Estimated correspondence probability", side=1, line=2, outer = TRUE)
mtext(text= "Actual probability",side=2,line=2,outer=TRUE)
dev.off()



