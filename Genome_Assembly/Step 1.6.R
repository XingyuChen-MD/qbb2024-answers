genome_coverage_30x <- scan("genome_coverage_30x.txt")

coverage_hist_30x <- hist(genome_coverage_30x, breaks = 50, plot = FALSE)

hist(genome_coverage_30x, breaks = 50, probability = TRUE, 
     main = "Coverage Distribution for 30x Coverage",
     xlab = "Coverage", ylab = "Density", col = "gray", border = "black")

lambda_30x <- 30
x_vals_30x <- 0:max(coverage_hist_30x$breaks)
poisson_probs_30x <- dpois(x_vals_30x, lambda_30x)
lines(x_vals_30x, poisson_probs_30x, type = "h", col = "red", lwd = 2)


mean_val_30x <- 30
std_dev_30x <- sqrt(30)
normal_probs_30x <- dnorm(x_vals_30x, mean = mean_val_30x, sd = std_dev_30x)
lines(x_vals_30x, normal_probs_30x, col = "blue", lwd = 2)

legend("topright", legend = c("Histogram", "Poisson (λ=30)", "Normal (mean=30, sd=5.47)"), 
       col = c("gray", "red", "blue"), lwd = 2)
png("ex1_30x_cov.png")
dev.off()

zero_coverage_count_30x <- sum(genome_coverage_30x == 0)
total_bases_30x <- length(genome_coverage_30x)
percentage_unsequenced_30x <- (zero_coverage_count_30x / total_bases_30x) * 100
print(percentage_unsequenced_30x)

