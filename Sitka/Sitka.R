
#### 1988 dataset ####

library(MASS)
data("Sitka")
colnames(Sitka) 
summary(as.factor(Sitka$tree)) # 5 measurements for all 79 trees

unique(Sitka$Time) # 152 174 201 227 258 - 5 unique time points

###### Rearranging for analysis ######

# response matrix will be 79 x 5
# covariate matrix will be 79 x 2, with intercept and indicator for ozone treatment

times <- unique(Sitka$Time)

n <- length(unique(Sitka$tree))
p <- 2 
q <- length(times)

y <- matrix(nrow = n, ncol = q)
x <- matrix(nrow = n, ncol = p)
x[,1] <- 1
for (i in 1:n) {
  for (j in 1:q) {
    y[i,j] <- Sitka$size[Sitka$tree == i & Sitka$Time == times[j]]
  }
  x[i,2] <- ifelse(Sitka$treat[Sitka$tree == i][1] == "ozone", 1, 0)
}

###### assessing normality ######

# marginal normality, without taking into consideration covariates
par(mfrow = c(2,3))
for (i in 1:q) {
  qqnorm(y[,i], main = bquote(Y[.(i)]))
  qqline(y[,i])
}

y.means <- colMeans(y)
y.sds <- apply(y, 2, sd)
y.centered <- sweep(y, 2, y.means)
y.stand <- sweep(y.centered, 2, y.sds, FUN = "/")

# conditional normality of autoregressions
par(mfrow = c(2,3))
l <- lm(y.stand[,1] ~ x)
qqnorm(residuals(l), main = bquote(Residuals~of~autoregression~of~Y[1]))
qqline(residuals(l))
for (i in 2:q) { 
  l <- lm(y.stand[,i] ~ y.stand[,1:(i-1)] + x)
  qqnorm(residuals(l), main = bquote(Residuals~of~autoregression~of~Y[.(i)])) 
  qqline(residuals(l))
}

###### save files ######

saveRDS(list(y = y, y.stand = y.stand, y.centered = y.centered, x = x),
        file = "Sitka1988_arranged.rds")

###### growth curves ######

ozonemeanresponse <- colMeans(y[x[,2] == 1,,drop = FALSE])
controlmeanresponse <- colMeans(y[x[,2] == 0,,drop = FALSE])

# plotting the growth patterns by growth condition
par(mfrow = c(1,2))
ozone <- which(x[,2] == 1)
control <- which(x[,2] == 0)
plot(times, y[ozone[1],], type = "l", main = "Ozone", xlab = "Days since January 1, 1988",
     ylab = "Log Size", ylim = c(2,8))
for (i in 2:length(ozone)) {
  lines(times, y[ozone[i],])
}
plot(times, y[control[1],], type = "l", main = "Control", xlab = "Days since January 1, 1988",
     ylab = "Log Size", ylim = c(2,8))
for (i in 2:length(control)) {
  lines(times, y[control[i],])
}
