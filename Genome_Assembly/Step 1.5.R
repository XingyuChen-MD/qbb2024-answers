


genome_coverage_10x <- scan("genome_coverage_10x.txt")

coverage_hist_10x <- hist(genome_coverage_10x, breaks = 50, plot = FALSE)

hist(genome_coverage_10x, breaks = 50, probability = TRUE, 
     main = "Coverage Distribution for 10x Coverage",
     xlab = "Coverage", ylab = "Density", col = "gray", border = "black")

lambda_10x <- 10
x_vals_10x <- 0:max(coverage_hist_10x$breaks)
poisson_probs_10x <- dpois(x_vals_10x, lambda_10x)
lines(x_vals_10x, poisson_probs_10x, type = "h", col = "red", lwd = 2)

mean_val_10x <- 10
std_dev_10x <- sqrt(10)
normal_probs_10x <- dnorm(x_vals_10x, mean = mean_val_10x, sd = std_dev_10x)
lines(x_vals_10x, normal_probs_10x, col = "blue", lwd = 2)

legend("topright", legend = c("Histogram", "Poisson (λ=10)", "Normal (mean=10, sd=sqrt(10))"), 
       col = c("gray", "red", "blue"), lwd = 2)
png("ex1_10x_cov.png")
# [Insert the plotting code here]
dev.off()
