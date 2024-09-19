

genome_coverage <- scan("genome_coverage.txt")

coverage_hist <- hist(genome_coverage, breaks = 50, plot = FALSE)

hist(genome_coverage, breaks = 50, probability = TRUE, 
     main = "Coverage Distribution for 3x Coverage",
     xlab = "Coverage", ylab = "Density", col = "gray", border = "black")

lambda <- 3
x_vals <- 0:max(coverage_hist$breaks)
poisson_probs <- dpois(x_vals, lambda)
lines(x_vals, poisson_probs, type = "h", col = "red", lwd = 2)

# Overlay Normal distribution with mean = 3 and std. dev. = sqrt(3)
mean_val <- 3
std_dev <- sqrt(3)
normal_probs <- dnorm(x_vals, mean = mean_val, sd = std_dev)
lines(x_vals, normal_probs, col = "blue", lwd = 2)

# Add legend
legend("topright", legend = c("Histogram", "Poisson (λ=3)", "Normal (mean=3, sd=sqrt(3))"), 
       col = c("gray", "red", "blue"), lwd = 2)

png("ex1_3x_cov.png")
dev.off()

zero_coverage_count <- sum(genome_coverage == 0)
total_bases <- length(genome_coverage)
percentage_unsequenced <- (zero_coverage_count / total_bases) * 100
print(percentage_unsequenced)
