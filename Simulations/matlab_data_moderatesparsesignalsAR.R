
#### Turning Simulation Datasets into usable files for Matlab ####

# This is the moderate sparse B, AR(1) Sigma setting

data <- readRDS("Data_n100p5q5_moderatesparsesignalsAR.rds")

nsets <- length(data$y)

# x and b
write.csv(data$x.stand, file = "moderatesparsesignalsAR_x.csv", row.names = FALSE)
write.csv(data$b, file = "moderatesparsesignalsAR_b.csv", row.names = FALSE)

for (i in 1:nsets) {
	y.name <- paste("moderatesparsesignalsAR_y", i, ".csv", sep = "")

	# save to csv
	write.csv(data$y.centered[[i]], file = y.name, row.names = FALSE)
}



