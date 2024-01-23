
#### GMCB timing results ####

reps <- 50

gmcb_mh_time <- vector("list", reps)
gmcb_smn_time <- vector("list", reps)
for (i in 1:reps) {
  cat("i = ", i, "\n")
  filename_mh <- paste("GMCBMH_run", i, ".rds", sep = "")
  filename_smn <- paste("GMCBSMN_run", i, ".rds", sep = "")
  run_mh <- readRDS(filename_mh)
  run_smn <- readRDS(filename_smn)
  
  gmcb_mh_time[[i]] <- run_mh$timing
  gmcb_smn_time[[i]] <- run_smn$timing
}

gmcb_mh_time_durationunitsplit <- lapply(gmcb_mh_time, function(x) strsplit(x, " "))
gmcb_smn_time_durationunitsplit <- lapply(gmcb_smn_time, function(x) strsplit(x, " "))

gmcb_mh_time_durationunitsplit <- lapply(gmcb_mh_time, function(x) strsplit(x, " "))
gmcb_smn_time_durationunitsplit <- lapply(gmcb_smn_time, function(x) strsplit(x, " "))

sapply(gmcb_mh_time_durationunitsplit, 
       function(x) c(x$`Pre-sampling computations`[2], x$Sampling[2]))
sapply(gmcb_smn_time_durationunitsplit, 
       function(x) c(x$`Pre-sampling computations`[2], x$Sampling[2]))
# all pre-sampling computations completed in seconds
# sampling completed in hours for both

# set in hours for ease of comparison with HSGHS
gmcb_mh_time_durationsum <- sapply(gmcb_mh_time_durationunitsplit,
                                     function(x) as.numeric(x$`Pre-sampling computations`[1])/3600 +
                                       as.numeric(x$Sampling[1]))
gmcb_smn_time_durationsum <- sapply(gmcb_smn_time_durationunitsplit,
                                     function(x) as.numeric(x$`Pre-sampling computations`[1])/3600 +
                                       as.numeric(x$Sampling[1]))

#### HSGHS timing results ####

hsghstime <- read.csv("hsghs_n100p120q50_timinginhours.csv", 
					  header = FALSE)

library(ggplot2)
timeresults <- data.frame(Algorithm = c(rep("GMCB-MH", reps), rep("GMCB-SMN", reps), rep("HS-GHS", reps)),
						  time = c(gmcb_mh_time_durationsum, gmcb_smn_time_durationsum, unlist(hsghstime)))
timeresults$Algorithm <- as.factor(timeresults$Algorithm)
ggplot(timeresults) + 
  geom_boxplot(aes(x = Algorithm, y = time)) + 
  ylab("Total Computation Time (hours)") + xlab("Algorithm")

# computation time mean
(time_mean <- data.frame(GMCB_MH = mean(gmcb_mh_time_durationsum),
           GMCB_SMN = mean(gmcb_smn_time_durationsum),
           HSGHS = mean(unlist(hsghstime))))

# standard errors
(time_se <- data.frame(GMCB_MH = sd(gmcb_mh_time_durationsum)/reps,
           GMCB_SMN = sd(gmcb_smn_time_durationsum)/reps,
           HSGHS = sd(unlist(hsghstime))/reps))

max(time_se)
