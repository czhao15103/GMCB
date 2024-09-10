
#### convert data to .csv for HSGHS ####

data <- readRDS("Data_n40p30q50.rds")
library(mvtnorm)

write.csv(data$x.stand, file = "n40p30q50_x.csv", row.names = FALSE)
write.csv(data$b, file = "n40p40q50_b.csv", row.names = FALSE)
write.csv(data$true.sigma, file = "n40p30q50_sigma.csv", row.names = FALSE)
write.csv(data$true.omega, file = "n40p30q50_omega.csv", row.names = FALSE)

for (i in 1:length(data$y.centered)) {
	write.csv(data$y.centered[[i]], file = paste("n40p30q50_y", i, ".csv", sep = ""), row.names = FALSE)
}

