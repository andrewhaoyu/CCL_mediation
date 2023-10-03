# Generate some example data
set.seed(123)

n1 <- 100
n2 <- 150

data1 <- data.frame(x = rnorm(n1), y = rnorm(n1))
data2 <- data.frame(x = rnorm(n2), y = rnorm(n2))

# Fit the regression models
fit1 <- lm(y ~ x, data = data1)
fit2 <- lm(y ~ x, data = data2)

beta1 <- coef(fit1)["x"]
beta2 <- coef(fit2)["x"]


library(boot)

# Bootstrap function
bootstrap_meta <- function(data, indices) {
  sample1 <- data1[indices[1:n1], ]
  sample2 <- data2[indices[(n1+1):(n1+n2)], ]
  
  fit1 <- lm(y ~ x, data = sample1)
  fit2 <- lm(y ~ x, data = sample2)
  
  beta1 <- coef(fit1)["x"]
  beta2 <- coef(fit2)["x"]
  
  w1 <- 1 / summary(fit1)$coefficients["x", "Std. Error"]^2
  w2 <- 1 / summary(fit2)$coefficients["x", "Std. Error"]^2
  
  combined_beta <- (beta1 * w1 + beta2 * w2) / (w1 + w2)
  
  return(combined_beta)
}

# Combine the datasets for bootstrapping
combined_data <- rbind(data1, data2)

# Run bootstrap
set.seed(456)
results <- boot(data = combined_data, statistic = bootstrap_meta, R = 1000)


# Confidence interval
conf_int <- quantile(results$t, c(0.025, 0.975))

# Two-tailed p-value for combined effect size being different from 0
p_val <- mean(abs(results$t) >= abs(mean(results$t)))

list(confidence_interval = conf_int, p_value = p_val)
