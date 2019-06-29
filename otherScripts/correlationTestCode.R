range = quantile(ut2017$ELEVATION, c(0, seq(0.25, 0.75, .1), 1))
range = quantile(ut2017$ELEVATION, seq(0.25, 0.75, .1), 1)
range[1] <- range[1] - 1
range[length(range)] <- range[length(range)] + 1
results <- results2 <- results3 <- matrix(0, ncol = 2, nrow = length(range))
for(i in 1:length(range)){
  results[i, 1] <- range[i]
  results[i, 2] <- check(ut2017, range[i], method = "kendall")

  results2[i, 1] <- range[i]
  results2[i, 2] <- check(ut1992, range[i], method = "kendall")

  results3[i, 1] <- range[i]
  results3[i, 2] <- check(id2015, range[i], method = "kendall")
}

print(results[which.min(results[,2]),1])
print(results2[which.min(results2[,2]),1])

par(mfrow = c(1, 3))
plot(results[,1], results[,2], type = "l")
plot(results2[,1], results2[,2], type = "l")
plot(results3[,1], results3[,2], type = "l")
abline(v = 1219.2, col = "red")
par(mfrow = c(1, 1))
