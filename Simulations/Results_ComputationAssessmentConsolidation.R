
library(GMCB)

#### GMCB-MH ####

mh_dp2002 <- readRDS("gmcb_mh_n100p5q5_DP2002_3A_XJ2019_1_ComputationAssessment.rds")
mh_mod <- readRDS("gmcb_mh_n100p5q5_moderatesparseAR_ComputationAssessment.rds")
mh_small <- readRDS("gmcb_mh_n100p5q5_smallsignalsCS_ComputationAssessment.rds")

reps <- 2000

# extract the ESS values
mh_dp2002_ess <- sapply(1:reps, function(x) mh_dp2002$ests[[x]]$postmean_ess)
mh_mod_ess <- sapply(1:reps, function(x) mh_mod$ests[[x]]$postmean_ess)
mh_small_ess <- sapply(1:reps, function(x) mh_small$ests[[x]]$postmean_ess)
mh_dp2002_qs_ess <- sapply(1:reps, function(x) mh_dp2002$ests[[x]]$steinquad_ess)
mh_mod_qs_ess <- sapply(1:reps, function(x) mh_mod$ests[[x]]$steinquad_ess)
mh_small_qs_ess <- sapply(1:reps, function(x) mh_small$ests[[x]]$steinquad_ess)

# extract time
mh_dp2002_time <- lapply(1:reps, function(x) mh_dp2002$ests[[x]]$time)
mh_dp2002_time_durationunitsplit <- lapply(mh_dp2002_time, function(x) strsplit(x, " "))
mh_dp2002_time_durationsum <- lapply(mh_dp2002_time_durationunitsplit,
                                     function(x) as.numeric(x$`Pre-sampling computations`[1]) +
                                       as.numeric(x$Sampling[1]))
summary(as.factor(sapply(mh_dp2002_time_durationunitsplit, 
                         function(x) c(x$`Pre-sampling computations`[2],
                                       x$Sampling[2])))) 
# all secs

mh_mod_time <- lapply(1:reps, function(x) mh_mod$ests[[x]]$time)
mh_mod_time_durationunitsplit <- lapply(mh_mod_time, function(x) strsplit(x, " "))
mh_mod_time_durationsum <- lapply(mh_mod_time_durationunitsplit,
                                     function(x) as.numeric(x$`Pre-sampling computations`[1]) +
                                       as.numeric(x$Sampling[1]))
summary(as.factor(sapply(mh_mod_time_durationunitsplit, 
                         function(x) c(x$`Pre-sampling computations`[2],
                                       x$Sampling[2])))) 
# all secs

mh_small_time <- lapply(1:reps, function(x) mh_small$ests[[x]]$time)
mh_small_time_durationunitsplit <- lapply(mh_small_time, function(x) strsplit(x, " "))
mh_small_time_durationsum <- lapply(mh_small_time_durationunitsplit,
                                     function(x) as.numeric(x$`Pre-sampling computations`[1]) +
                                       as.numeric(x$Sampling[1]))
summary(as.factor(sapply(mh_small_time_durationunitsplit, 
                         function(x) c(x$`Pre-sampling computations`[2],
                                       x$Sampling[2])))) 
# all secs

#### GMCB-SMN ####

smn_dp2002 <- readRDS("gmcb_smn_n100p5q5_DP2002_3A_XJ2019_1_ComputationAssessment.rds")
smn_mod <- readRDS("gmcb_smn_n100p5q5_moderatesparseAR_ComputationAssessment.rds")
smn_small <- readRDS("gmcb_smn_n100p5q5_smallsignalsCS_ComputationAssessment.rds")

# ess
smn_dp2002_f_ess <- sapply(1:reps, function(x) smn_dp2002$ests[[x]]$postmean_ess)
smn_dp2002_qs_ess <- sapply(1:reps, function(x) smn_dp2002$ests[[x]]$steinquad_ess)
smn_mod_f_ess <- sapply(1:reps, function(x) smn_mod$ests[[x]]$postmean_ess)
smn_mod_qs_ess <- sapply(1:reps, function(x) smn_mod$ests[[x]]$steinquad_ess)
smn_small_f_ess <- sapply(1:reps, function(x) smn_small$ests[[x]]$postmean_ess)
smn_small_qs_ess <- sapply(1:reps, function(x) smn_small$ests[[x]]$steinquad_ess)

# extract time
smn_dp2002_time <- lapply(1:reps, function(x) smn_dp2002$ests[[x]]$time)
smn_dp2002_time_durationunitsplit <- lapply(smn_dp2002_time, function(x) strsplit(x, " "))
smn_dp2002_time_durationsum <- lapply(smn_dp2002_time_durationunitsplit,
                                     function(x) as.numeric(x$`Pre-sampling computations`[1]) +
                                       as.numeric(x$Sampling[1]))
summary(as.factor(sapply(smn_dp2002_time_durationunitsplit, 
                         function(x) c(x$`Pre-sampling computations`[2],
                                       x$Sampling[2])))) 
# all secs

smn_mod_time <- lapply(1:reps, function(x) smn_mod$ests[[x]]$time)
smn_mod_time_durationunitsplit <- lapply(smn_mod_time, function(x) strsplit(x, " "))
smn_mod_time_durationsum <- lapply(smn_mod_time_durationunitsplit,
                                  function(x) as.numeric(x$`Pre-sampling computations`[1]) +
                                    as.numeric(x$Sampling[1]))
summary(as.factor(sapply(smn_mod_time_durationunitsplit, 
                         function(x) c(x$`Pre-sampling computations`[2],
                                       x$Sampling[2])))) 
# all secs

smn_small_time <- lapply(1:reps, function(x) smn_small$ests[[x]]$time)
smn_small_time_durationunitsplit <- lapply(smn_small_time, function(x) strsplit(x, " "))
smn_small_time_durationsum <- lapply(smn_small_time_durationunitsplit,
                                    function(x) as.numeric(x$`Pre-sampling computations`[1]) +
                                      as.numeric(x$Sampling[1]))
summary(as.factor(sapply(smn_small_time_durationunitsplit, 
                         function(x) c(x$`Pre-sampling computations`[2],
                                       x$Sampling[2])))) 
# all secs

# checking whether any ess estimate exceeds number of iterations
length(which(smn_dp2002_f_ess > 1e5)) # 0
length(which(smn_dp2002_qs_ess > 1e5)) # 0
length(which(smn_mod_f_ess > 1e5)) # 0
length(which(smn_mod_qs_ess > 1e5)) # 0
length(which(smn_small_f_ess > 1e5)) # 0
length(which(smn_small_qs_ess > 1e5)) # 0

#### HSGHS ####

# times in seconds
hsghs_dp2002_timing <- read.csv("hsghs_n100p5q5_DP2002_3A_XJ2019_1_timinginsecs.csv",
                         header = FALSE)
hsghs_mod_timing <- read.csv("hsghs_n100p5q5_moderatesparsesignalsAR_timinginsecs.csv",
                             header = FALSE)
hsghs_small_timing <- read.csv("hsghs_n100p5q5_smallsignalsCS_timinginsecs.csv",
                               header = FALSE)

hsghs_dp2002 <- readRDS("ComputationEffortResults_DP2002_3A_XJ2019_1.rds")
hsghs_mod <- readRDS("ComputationEffortResults_moderatesparsesignalsAR.rds")
hsghs_small <- readRDS("ComputationEffortResults_smallsignals.rds")

# extract the ESS values
hsghs_dp2002_ess <- sapply(1:reps, function(x) hsghs_dp2002[[x]]$postmean_ess)
hsghs_mod_ess <- sapply(1:reps, function(x) hsghs_mod[[x]]$postmean_ess)
hsghs_small_ess <- sapply(1:reps, function(x) hsghs_small[[x]]$postmean_ess)
hsghs_dp2002_qs_ess <- sapply(1:reps, function(x) hsghs_dp2002[[x]]$steinquad_ess)
hsghs_mod_qs_ess <- sapply(1:reps, function(x) hsghs_mod[[x]]$steinquad_ess)
hsghs_small_qs_ess <- sapply(1:reps, function(x) hsghs_small[[x]]$steinquad_ess)

#### Collecting results in the same dataframe ####

# collecting into dataframe
nscenarios <- 3
results <- data.frame(Algorithm = rep(rep(rep(c("GMCB-MH", "GMCB-SMN", "HSGHS"), 
                                              rep(reps, 3)), nscenarios), 2),
                      Estimator = c(rep("F", reps*3*nscenarios), rep("QS", reps*3*nscenarios)),
                      scenario = rep(rep(1:nscenarios, rep(reps*3, nscenarios)), 2),
                      ess = c(mh_small_ess, smn_small_f_ess, hsghs_small_ess,
                              mh_mod_ess, smn_mod_f_ess, hsghs_mod_ess,
                              mh_dp2002_ess, smn_dp2002_f_ess, hsghs_dp2002_ess,
                              mh_small_qs_ess, smn_small_qs_ess, hsghs_small_qs_ess,
                              mh_mod_qs_ess, smn_mod_qs_ess, hsghs_mod_qs_ess,
                              mh_dp2002_qs_ess, smn_dp2002_qs_ess, hsghs_dp2002_qs_ess))
results$Algorithm <- as.factor(results$Algorithm)
results$scenario <- as.factor(results$scenario)

timeresults <- data.frame(Algorithm = c(rep("GMCB-MH", reps*nscenarios),
                                        rep("GMCB-SMN", reps*nscenarios),
                                        rep("HSGHS", reps*nscenarios)),
                          scenario = rep(rep(1:nscenarios, rep(reps, nscenarios)), 3),
                          time = c(unlist(mh_small_time_durationsum), 
                                   unlist(mh_mod_time_durationsum), 
                                   unlist(mh_dp2002_time_durationsum), 
                                   unlist(smn_small_time_durationsum),
                                   unlist(smn_mod_time_durationsum),
                                   unlist(smn_dp2002_time_durationsum),
                                   unlist(hsghs_small_timing),
                                   unlist(hsghs_mod_timing),
                                   unlist(hsghs_dp2002_timing)))
timeresults$scenario <- as.factor(timeresults$scenario)
timeresults$Algorithm <- as.factor(timeresults$Algorithm)

#### Generating Values for Tables ####

# multivariate ESS
(ess_mean <- aggregate(results$ess, list(Scenario = results$scenario, Algorithm = results$Algorithm,
                            Estimator = results$Estimator), mean))

(ess_se <- aggregate(results$ess, list(Scenario = results$scenario, Algorithm = results$Algorithm,
                            Estimator = results$Estimator), function(x) sd(x)/2000))

max(ess_se$x)

# computation time - in seconds
(time_mean <- aggregate(timeresults$time, list(Scenario = timeresults$scenario,
                            Algorithm = timeresults$Algorithm), mean))

(time_se <- aggregate(timeresults$time, list(Scenario = timeresults$scenario,
                                 Algorithm = timeresults$Algorithm), function(x) sd(x)/2000))

max(time_se$x) 

#### Plotting ####

library(ggplot2)
library(ggpubr)

##### Total Computation time #####

# only GMCB
timeplot_gmcb <- ggplot(timeresults[timeresults$Algorithm != "HSGHS",]) + 
  geom_boxplot(aes(x = scenario, y = time, linetype = Algorithm)) + 
  ylab("") + xlab("Scenario") +
  theme(legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=14), legend.title = element_text(size = 16),
        axis.title.x = element_text(size=16), axis.text = element_text(size = 14))
timeplot_gmcb

#### ACF plots ####

p <- 5
q <- 5

##### small signals, CS #####

smallcs_mh <- readRDS("acf_gmcb_mh_n100p5q5_smallsignalsCS.rds")
smallcs_smn <- readRDS("acf_gmcb_smn_n100p5q5_smallsignalsCS.rds")

set.seed(1)
(b.index <- sample(1:(p*q), 1)) # 25 
matrix(1:25, nrow = 5) # 25 corresponds to b[5,5]
(d.index <- sample(1:(q*(q-1)/2), 1)) # 4 
ref.t <- matrix(0, nrow = q, ncol = q)
ref.t[upper.tri(ref.t)] <- 1:(q*(q-1)/2)
t(ref.t) # 4 corresponds to delta[3,1]

par(mfrow = c(1, 2)) 
acf(smallcs_mh$mcmc$b[,b.index], lag.max = 50, main = bquote("B"[55]))
acf(smallcs_mh$mcmc$delta[,d.index], lag.max = 50, main = bquote(delta[4][","][1]))

acf(smallcs_smn$mcmc$b[,b.index], lag.max = 50, main = bquote("B"[55]))
acf(smallcs_smn$mcmc$delta[,d.index], lag.max = 50, main = bquote(~delta[4][","][1]))

##### moderate signals, AR #####

modar_mh <- readRDS("acf_gmcb_mh_n100p5q5_moderatesparsesignalsAR.rds")
modar_smn <- readRDS("acf_gmcb_smn_n100p5q5_moderatesparsesignalsAR.rds")

set.seed(2)
(b.index <- sample(1:(p*q), 1)) # 21 
matrix(1:25, nrow = 5) # 21 corresponds to b[5,1]
(d.index <- sample(1:(q*(q-1)/2), 1)) # 6 
ref.t <- matrix(0, nrow = q, ncol = q)
ref.t[upper.tri(ref.t)] <- 1:(q*(q-1)/2)
t(ref.t) # 6 corresponds to delta[4, 3]

par(mfrow = c(1, 2)) 
acf(modar_mh$mcmc$b[,b.index], lag.max = 50, main = bquote("B"[51]))
acf(modar_mh$mcmc$delta[,d.index], lag.max = 50, main = bquote(delta[4][","][3]))

acf(modar_smn$mcmc$b[,b.index], lag.max = 50, main = bquote("B"[51]))
acf(modar_smn$mcmc$delta[,d.index], lag.max = 50, main = bquote(delta[4][","][3]))

##### DP2002 3A, XJ2019 1 #####

dp2002_mh <- readRDS("acf_gmcb_mh_n100p5q5_DP2002_3A_XJ2019_1.rds")
dp2002_smn <- readRDS("acf_gmcb_smn_n100p5q5_DP2002_3A_XJ2019_1.rds")

set.seed(3)
(b.index <- sample(1:(p*q), 1)) # 5 
matrix(1:25, nrow = 5) # 5 corresponds to b[1,5]
(d.index <- sample(1:(q*(q-1)/2), 1)) # 10 
ref.t <- matrix(0, nrow = q, ncol = q)
ref.t[upper.tri(ref.t)] <- 1:(q*(q-1)/2)
t(ref.t) # 10 corresponds to delta[5, 4]

par(mfrow = c(1, 2)) 
acf(dp2002_mh$mcmc$b[,b.index], lag.max = 50, main = bquote("B"[15]))
acf(dp2002_mh$mcmc$delta[,d.index], lag.max = 50, main = bquote(delta[5][","][4]))

acf(dp2002_smn$mcmc$b[,b.index], lag.max = 50, main = bquote("B"[15]))
acf(dp2002_smn$mcmc$delta[,d.index], lag.max = 50, main = bquote(delta[5][","][4]))
