
#### Turning Simulation Datasets into usable files for Matlab ####

# DP 2002 3A XJ 2019 1

data <- readRDS("Data_n100p5q5_DP2002_3A_XJ2019_1.rds")

nsets <- length(data$y)

# x and b
write.csv(data$x.stand, file = "n100p5q5_DP2002_3A_XJ2019_1_x.csv", row.names = FALSE)
write.csv(data$b, file = "n100p5q5_DP2002_3A_XJ2019_1_b.csv", row.names = FALSE)
write.csv(data$y.cov, file = "n100p5q5_DP2002_3A_XJ2019_1_cov.csv", row.names = FALSE)
write.csv(data$y.prec, file = "n100p5q5_DP2002_3A_XJ2019_1_prec.csv", row.names = FALSE)

for (i in 1:nsets) {
	y.name <- paste("n100p5q5_DP2002_3A_XJ2019_1_y", i, ".csv", sep = "")

	# save to csv
	write.csv(data$y.centered[[i]], file = y.name, row.names = FALSE)
}



