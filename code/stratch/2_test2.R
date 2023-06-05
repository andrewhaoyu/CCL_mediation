#goal: test the polytomous model to get results

n_rep = 1000
n = 15000
prob = c(1/3,1/3,1/3)
set.seed(123)  # For reproducibility

# Sample size
beta_1 = rep(0,n_rep)
beta_2  = rep(0,n_rep)
for(k in 1:n_rep){
  if(k%%10){print(k)}
  # Generate the control group (0), case type 1 (1), and case type 2 (2)
  group <- factor(sample(0:2, n, replace = TRUE, prob = c(0.4, 0.3, 0.3)))
  
  # Generate predictors (e.g., x1 and x2)
  x1 <- rnorm(n)
  idx <- which(group==0|group==1)
  y_new = group[idx]
  x_new = x1[idx]
  model = glm(y_new ~ x_new, family = "binomial")
  beta_1[k] = coef(model)[2]
  idx <- which(group==0|group==2)
  y_new = group[idx]
  x_new = x1[idx]
  model = glm(y_new ~ x_new, family = "binomial")
  
  beta_2[k] = coef(model)[2]
  
}
cor(beta_1,beta_2)
n1 = n*0.3
n2 = n*0.3
n0 = n*0.4
sqrt(n1*n2/((n0+n2)*(n0+n1)))
