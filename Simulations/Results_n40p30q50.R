
#### n = 40, p = 30, q = 50 Results ####

mssl <- readRDS("mSSL_n40p30q50_losses.rds")
smn <- readRDS("n40p30q50_SMNlosses.rds")
hsghs <- readRDS("HSGHS_n40p30q50_losses.rds")

rbind(smn = colMeans(smn), 	  hsghs = colMeans(hsghs))

colMeans(mssl)


rbind(sapply(1:ncol(smn), function(x) sd(smn[,x])/nrow(smn)),
	  sapply(1:ncol(hsghs), function(x) sd(hsghs[,x])/nrow(hsghs)))

sapply(1:ncol(mssl), function(x) sd(mssl[,x])/nrow(mssl))

allsds <- c(sapply(1:4, function(x) sd(smn[,x])/nrow(smn)),
	  sapply(1:4, function(x) sd(hsghs[,x])/nrow(hsghs)),
	  sapply(1:4, function(x) sd(mssl[,x])/nrow(mssl)))
max(allsds)
