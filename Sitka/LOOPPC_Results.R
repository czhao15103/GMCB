
#### LOO PPC results ####

setwd("C:/Users/clz15/Desktop/ResearchCode/RealData/SitkaDataset")

mh0112001 <- readRDS("LOOPPC_MH_0112001.rds")
mh114005 <- readRDS("LOOPPC_MH_114005.rds")
smn0112001 <- readRDS("LOOPPC_SMN_0112001.rds")
smn114005 <- readRDS("LOOPPC_SMN_114005.rds")

results <- data.frame(mh0112001 = mh0112001,
                      mh114005 = mh114005,
                      smn0112001 = smn0112001,
                      smn114005 = smn114005)

colMeans(results) # -0.5692786 -0.5693516 -0.5746820 -0.5731782 
apply(results, 2, function(x) sd(x)/nrow(results)) # 0.01722958 0.01722822 0.01712938 0.01717993 
