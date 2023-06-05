
result = read.csv("/Users/zhangh24/Downloads/LDpred_CT_15sets.csv")
result = result[c(3,4,7,8,11,12,15),]
n = result$n_eff
y = result$ratio_var
n_data = nrow(result)
#assume n_effect more than 50000, the ratio gets 1
n[n_data+(1:6)] = c(50000,60000,70000,80000,90000,100000)
y[n_data+c(1:6)] = 1.00 + abs(rnorm(6,0,1E-06))

library(minpack.lm)
# Define the models
exp_decay <- function(n, a, b) {
  a * exp(-b * n) + 1
}

power_law <- function(n, a, b) {
  a * n^(-b) + 1
}


logistic <- function(n, a, b, c) {
  a / (1 + exp(-b * (n - c))) + (1 - a)
}

gompertz <- function(n, a, b, c) {
  a * exp(-b * exp(-c * n)) + 1
}


weibull <- function(n, a, b, c) {
  a * exp(-b * n^c) + 1
}
# Fit the models
exp_decay_fit <- nlsLM(y ~ exp_decay(n, a, b), start = list(a = 3.263, b = 0.001))
power_law_fit <- nlsLM(y ~ power_law(n, a, b), start = list(a = 3.263, b = 0.5),)
logistic_fit <- nlsLM(y ~ logistic(n, a, b, c), start = list(a = 3.263, b = 0.001, c = 1500))
gompertz_fit <- nlsLM(y ~ gompertz(n, a, b, c), start = list(a = 4.263, b = 0.02, c = 0.0001))
weibull_fit <- nlsLM(y ~ weibull(n, a, b, c), start = list(a = 4.263, b = 0.02, c = 0.05))

# Calculate R-squared values for each model
r2 <- function(fit) 1 - sum(residuals(fit)^2) / sum((y - mean(y))^2)

r_squared <- c(
  exp_decay = r2(exp_decay_fit),
  power_law = r2(power_law_fit),
  logistic = r2(logistic_fit),

  gompertz = r2(gompertz_fit),
  weibull = r2(weibull_fit)
)

# Determine the best-fitting model
best_model_name <- names(which.max(r_squared))
best_model_fit <- get(paste0(best_model_name, "_fit"))
best_model_coef <- coef(best_model_fit)
# Print best model name and R-squared
cat("Best Model:", best_model_name, "\n")
cat("R-squared:", r_squared[best_model_name], "\n")

# Compute fitted values
fitted_values <- fitted(best_model_fit)

# Determine the best-fitting model
best_model <- best_model_fit  

# Compute fitted values
fitted_values <- fitted(best_model)

# Create a plot of the original data points
plot(n[1:n_data], y[1:n_data], main = "Data with Best Fitting Model",
     xlab = "n", ylab = "y",
     pch = 19, col = "blue")

# Add the best-fitting curve to the plot
lines(n, fitted_values, col = "red", lwd = 2)

# Add a legend
legend("topleft", legend = c("Data", "Best Fit"),
       col = c("blue", "red"), lty = 1, pch = c(19, NA), lwd = c(NA, 2))





n_values <- seq(75, 100000, length.out = 1000)
y_weibull <- weibull(n_values, a = 3.51102, b = 0.04858, c= 0.42930)

# Extract the best fit parameters for the logistic function
logistic_params <- coef(logistic_fit)
y_logistic <- logistic(n_values, logistic_params["a"], logistic_params["b"], logistic_params["c"])

# Plot the results
plot(n_values, y_weibull, type = "l", col = "red", ylim = c(0, 4.5), xlab = "Effective Sample Size", ylab = "Ratio of PRS Phenotypic Variance", main = "Comparison of Fitted Functions")
lines(n_values, y_logistic, col = "blue")
legend("topright", legend = c("Weibull", "Logistic"), col = c("red", "blue"), lty = 1, bty = "n")

