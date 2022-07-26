Syy <- cov(multivardata[, 5:6])
Sxy <- cov(cbind(multivardata[, 1], multivardata[, 2], multivardata[, 5:6])) 
Syx <- cov(cbind(multivardata[, 5:6], multivardata[, 1], multivardata[, 2]))
Sxx <- cov(cbind(multivardata[, 1], multivardata[, 2]))

rsquare <- solve(Syy) %*% Syx 
rsquare

x <- cbind(multivardata[, 1], multivardata[, 2])
y <- multivardata[, 5:7]

bhat1 <- solve(t(x) %*% x)
bhat2 <-  t(x) %*% y
