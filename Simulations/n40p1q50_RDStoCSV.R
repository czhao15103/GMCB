
#### Creating the .csv files for HSGHS ####

data <- readRDS("Data_n40q50_meancovariance_longrange.rds")

for (i in 1:length(data$y)) {
	write.csv(data$y[[i]], file = paste("n40p1q50_longrange_y", i, ".csv", sep = ""), row.names = FALSE)
}

#### If we run GHS on this same data ####

data <- readRDS("Data_n40q50_meancovariance_longrange.rds")

for (i in 1:length(data$y)) {
	ghsy <- sweep(data$y[[i]], 2, colMeans(data$y[[i]]))
	write.csv(data$y[[i]], file = paste("n40p1q50_meancovariance_longrange_ghs_y", i, ".csv", sep = ""), row.names = FALSE)
}

