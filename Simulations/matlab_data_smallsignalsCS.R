
#### Turning Simulation Datasets into usable files for Matlab ####

# This is the dense B, compound symmetric Sigma setting

# functions and packages
data <- readRDS("Data_n100p5q5_smallsignalsCS.rds")

nsets <- length(data$y)

# x and b
write.csv(data$x.stand, file = "smallsignalsCS_x.csv", row.names = FALSE)
write.csv(data$b, file = "smallsignalsCS_b.csv", row.names = FALSE)

for (i in 1:nsets) {
	y.name <- paste("smallsignalsCS_y", i, ".csv", sep = "")

	# save to csv
	write.csv(data$y.centered[[i]], file = y.name, row.names = FALSE)
}



