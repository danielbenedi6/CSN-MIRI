library(ggplot2)
library(stats)

read_dataset <- function(language) {
  dataset = read.table(paste("./data/",language,"_dependency_tree_metrics.txt",sep=""), header=FALSE)
  colnames(dataset) = c("vertices","degree_2nd_moment", "mean_length")
  dataset = dataset[order(dataset$vertices),]
  dataset
}

generate_table_entry <- function(language, dataset) {
  sprintf(r"(%s & %d & %.3f & %.3f & %.3f & %.3f \\)",
          language,
          length(dataset),
          mean(dataset$vertices),
          sd(dataset$vertices),
          mean(dataset$mean_length),
          sd(dataset$mean_length)
          )
}

print_row_type2 <- function(language, model0, model1, model2, model3, model4, model1p, model2p, model3p, model4p){
  sprintf(r"( %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\)",
          language,
          model0,  # Model 0
          model1,  # Model 1
          model2,  # Model 2
          model3,  # Model 3
          model4,  # Model 4
          model1p, # Model 1+
          model2p, # Model 2+
          model3p, # Model 3+
          model4p  # Model 4+
  )
}

model1 <- function(n,b) {
  (n/2)**b
}
log_model1 <- function(n,b) {
  b * log(n/2)
}
model2 <- function(n,a,b) {
  a*n**b
}
log_model2 <- function(n,a,b) {
  log(a) + b*log(n)
}
model3 <- function(n,a,c) {
  a*exp(c*n)
}
log_model3 <- function(n,a,c) {
  log(a) + c*n
}
model4 <- function(n,a) {
  a*log(n)
}
log_model4 <- function(n,a) {
  log(a) + log(log(n))
}
model1_plus <- function(n.b,d) {
  (n/2)**b +d
}
log_model1_plus <- function(n,b,d) {
  log((n/2)**b + d)
}
model2_plus <- function(n,a,b,d) {
  a*n**b+d
}
log_model2_plus <- function(n,a,b,d) {
  log(a*n**b+d)
}
model3 <- function(n,a,c,d) {
  a*exp(c*n)+d
}
log_model3 <- function(n,a,c,d) {
  log(a*exp(c*n)+d)
}
model4 <- function(n,a) {
  a*log(n) + d
}
log_model4 <- function(n,a) {
  log(a*log(n) + d)
}


languages = c("Arabic", "Basque", "Catalan",
              "Chinese", "Czech", "English",
              "Greek", "Hungarian", "Italian",
              "Turkish")
table_summary <- r"(\begin{table}[!htb]
\centering
\resizebox{\columnwidth}{!}{
\begin{tabular}{llllll}{
Language & N & $\mu_n$ & $\sigma_n$ & $\mu_x$ & $\sigma_x$ \\ \hline
)"
table_res_se <- r"(\begin{table}[!htb]
\centering
\resizebox{\columnwidth}{!}{
\begin{tabular}{llllllllll}{
Language & 0 & 1 & 2 & 3 & 4 & 1+ & 2+ & 3+ & 4+ \\ \hline
)"
table_aic <- r"(\begin{table}[!htb]
\centering
\resizebox{\columnwidth}{!}{
\begin{tabular}{llllllllll}{
Language & 0 & 1 & 2 & 3 & 4 & 1+ & 2+ & 3+ & 4+ \\ \hline
)"
table_aic_diff <- r"(\begin{table}[!htb]
\centering
\resizebox{\columnwidth}{!}{
\begin{tabular}{llllllllll}{
Language & 0 & 1 & 2 & 3 & 4 & 1+ & 2+ & 3+ & 4+ \\ \hline
)"
table_params <- r"(\begin{table}[!htb]
\centering
\resizebox{\columnwidth}{!}{
\begin{tabular}{lllllllllllllllll}{
         & \multicolumn{16}{l}{Model} \\ \cline{2-17} 
         & \multicolumn{1}{l|}{1} & \multicolumn{2}{l|}{2}     & \multicolumn{2}{l|}{3}     & \multicolumn{1}{l|}{4} & \multicolumn{2}{l|}{1+}    & \multicolumn{3}{l|}{2+}        & \multicolumn{3}{l|}{3+}        & \multicolumn{2}{l}{4+} \\
Language & \multicolumn{1}{l|}{b} & a & \multicolumn{1}{l|}{b} & a & \multicolumn{1}{l|}{c} & \multicolumn{1}{l|}{a} & b & \multicolumn{1}{l|}{d} & a & b & \multicolumn{1}{l|}{d} & a & c & \multicolumn{1}{l|}{d} & a & d \\ \hline
)"

for(language in languages){
  data = read_dataset(language)
  table_summary <- paste(table_summary,generate_table_entry(language,data), sep="\n")
  
  means = aggregate(data, list(data$vertices), mean)
  
  plt <- ggplot() + geom_point(data=dataset, aes(x=vertices,y=mean_length))
  plt <- plt + geom_line(data=means, aes(x=vertices,y=mean_length), color="green")
  plt <- plt + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
  plt <- plt + labs(title=paste("Models for mean length dependency for ", language))
  plt <- plt + labs(colour="Model")
  plt <- plt + xlab("Number of vertices")
  plt <- plt + ylab("Mean length")
  
  table_res_se <- paste(table_res_se, print_row_type2(language,
                                                      0., # Model 0
                                                      0., # Model 1
                                                      0., # Model 2
                                                      0., # Model 3
                                                      0., # Model 4
                                                      0., # Model 1+
                                                      0., # Model 2+
                                                      0., # Model 3+
                                                      0.  # Model 4+
  ), sep="\n")
  table_aic <- paste(table_aic, print_row_type2(language,
                                                      0., # Model 0
                                                      0., # Model 1
                                                      0., # Model 2
                                                      0., # Model 3
                                                      0., # Model 4
                                                      0., # Model 1+
                                                      0., # Model 2+
                                                      0., # Model 3+
                                                      0.  # Model 4+
  ), sep="\n")
  table_aic_diff <- paste(table_aic_diff, print_row_type2(language,
                                                      0., # Model 0
                                                      0., # Model 1
                                                      0., # Model 2
                                                      0., # Model 3
                                                      0., # Model 4
                                                      0., # Model 1+
                                                      0., # Model 2+
                                                      0., # Model 3+
                                                      0.  # Model 4+
  ), sep="\n")
  table_params <- paste(table_params, sprintf(r"( %s & \multicolumn{1}{l|}{%.3f} & %.3f & \multicolumn{1}{l|}{%.3f} & %.3f & \multicolumn{1}{l|}{%.3f} & \multicolumn{1}{l|}{%.3f} & %.3f & \multicolumn{1}{l|}{%.3f} & %.3f & %.3f & \multicolumn{1}{l|}{%.3f} & %.3f & %.3f & \multicolumn{1}{l|}{%.3f} & %.3f & %.3f \\)",
                                              language,
                                              param_b_model1,
                                              param_a_model2, param_b_model2,
                                              param_a_model3, param_c_model3,
                                              param_a_model4,
                                              param_b_model1p, param_d_model1p,
                                              param_a_model2p, param_b_model2p, param_d_model2p,
                                              param_a_model3p, param_c_model3p, param_d_model3p,
                                              param_a_model4p, param_d_model4p
  ), sep="\n")
}
table_summary <- paste(table_summary, r"(\end{tabular}}
\caption{Summary of the properties of the degree sequences \label{tab:summary}}
\end{table})",sep="\n")
table_res_se <- paste(table_res_se, r"(\end{tabular}}
\caption{Residual standard error for each model \label{tab:res_se}}
\end{table})",sep="\n")
table_aic <- paste(table_aic, r"(\end{tabular}}
\caption{Akaike information criterion (AIC) of each model \label{tab:aic}}
\end{table})",sep="\n")
table_aic_diff <- paste(table_aic_diff, r"(\end{tabular}}
\caption{AIC differences \label{tab:aic_diff}}
\end{table})",sep="\n")
table_params <- paste(table_params, r"(\end{tabular}}
\caption{Parameters of each fitted model \label{tab:params}}
\end{table})",sep="\n")


# Data Analysis

# 0. loads info of catalan dependency trees, sorting it increasingly by # of vertices
Catalan = read.table("./data/Catalan_dependency_tree_metrics.txt", header = FALSE)
colnames(Catalan) = c("vertices","degree_2nd_moment", "mean_length")
Catalan = Catalan[order(Catalan$vertices), ]

# 1. preliminary plot
plot(Catalan$vertices, Catalan$mean_length,
     xlab = "vertices", ylab = "mean dependency length")
# same but taking logs on both axes --> 
## suggest a power-law dependency between mean length and number of vertices,
## in spite of the high dispersion. 
plot(log(Catalan$vertices), log(Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean dependency length)")
# even clearer, it can be seen by averaging the mean lengths for a given # of 
# vertices
mean_Catalan = aggregate(Catalan, list(Catalan$vertices), mean)
plot(mean_Catalan$vertices, mean_Catalan$mean_length,
     xlab = "vertices", ylab = "mean mean dependency length")
# in log scale
plot(log(mean_Catalan$vertices), log(mean_Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean mean dependency length)")

# NOTE: intuition about how far the real scaling of the mean dependency length
# {d} is from the random linear arrangement can be seen comparing both plots
# -> random linear arrangement coming from expected mean length = (n+1)/3 --> NULL MODEL
## real mean dependency length -> green
## expected mean ...           -> red
## (plot in double log scale)
## NOTE: plot for {d} versus n suggests an almost power-law dependency.
plot(log(Catalan$vertices), log(Catalan$mean_length),
       xlab = "vertices", ylab = "mean dependency length")
lines(log(mean_Catalan$vertices),log(mean_Catalan$mean_length), col = "green")
lines(log(mean_Catalan$vertices),log((mean_Catalan$vertices+1)/3), col = "red")

# 2. Ensemble of models..
# Non-linear regression -> fit of the models

# Model 2: an^b (a power-law model)
a_initial = 4
b_initial = 4
nonlinear_model = nls(mean_length ~ a*vertices^b,data=Catalan,
                      start = list(a = a_initial, b = b_initial), trace = TRUE)

# ... before, arbitrary initial values 
# --> obtention of good initial values for the params
linear_model = lm(log(mean_length)~log(vertices), Catalan)
# retrieving the coefs of the linear regression model done to obtain the initial
# values for the params -> a + bx 
a_initial = exp(coef(linear_model)[1])
b_initial = coef(linear_model)[2]

# Running again the non-linear regression with "better" initial values for the 
# params --> faster convergence
nonlinear_model = nls(mean_length~a*vertices^b,data=Catalan,
                      start = list(a = a_initial, b = b_initial), trace = TRUE)

# Errors of the non-linear regression model -> to measure the quality of the fit
deviance(nonlinear_model)
AIC(nonlinear_model)

# residual standard error 
# - measures the standard deviation of the residuals in a regression model
# - the smaller the residual standard error, the better a regression model fits 
#   a dataset
# - sqrt(RSS / df) -> sum of squared residuals (RSS) / degrees of freedom 
#   (= # of observations - # of model parameters)

sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model))

# Params giving the best fit for the model
coef(nonlinear_model)
coef(nonlinear_model)["a"]
coef(nonlinear_model)["b"]

# The RSS, s and AIC for a non-parametric model (such as the null model)
# {d} = n/3 + 1/3
RSS <- sum((Catalan$mean_length-(Catalan$vertices+1)/3)^2)
n <- length(Catalan$vertices)
p <- 0
s <- sqrt(RSS/(n - p))
AIC <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)


# Plot the empirical data and the curve for the best fit:
## empirical data plot
plot(log(Catalan$vertices), log(Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean dependency length)")
## best fit plot
lines(log(Catalan$vertices), log(fitted(nonlinear_model)), col = "green")


################################################################################
# Fitting the models for each of the languages
languages = c("Arabic", "Basque", "Catalan",
              "Chinese", "Czech", "English",
              "Greek", "Hungarian", "Italian",
              "Turkish")

for(language in languages) {
  cat("---------------------------------------------------------------------\n")
  cat("Language: ", language, "\n")
  language = "Arabic"
  data = read_dataset(language)
  # Homocesdasticity test -> checking the assumption of the homogenity of variance
  # ...
  # Aggregation of the data in the case assumption does not hold (it does not for
  # any language)
  data = aggregate(data, list(data$vertices), mean)
  
  models <- list()
  s_list <- list()
  AIC_list <- list()
  # Fit of the models
  # Initial values for the params -> 1 for all the params: a, b and d
  # a_initial = 1
  # b_initial = 1
  # d_initial = 1
  cat("################################## Model 1  ##################################")
  # f(n) = (n/2)^b
  b_initial = 1
  model_1 = nls(mean_length~(vertices/2)^b,data=data,
                        start = list(b = b_initial), trace = TRUE)
  model_1
  models[["model_1"]] <- model_1
  # ---------------------------------------------------------------------------#
  cat("... with optimized initial values: \n")
  # Opt: better init value for b param -> faster convergence
  linear_model = lm(log(mean_length)~log(vertices/2), data=data)
  b_init = coef(linear_model)[2]
  nls(mean_length~(vertices/2)^b,data=data,
                        start = list(b = b_init), trace = TRUE)
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  coef(model_1)
  coef(model_1)["b"]
  # Errors of the non-linear regression model
  s <- sqrt(deviance(model_1)/df.residual(model_1))
  s_list <- append(s_list, s)
  cat("Residual Standard Error (s): \n")
  print(s)
  AIC_value <- AIC(model_1)
  AIC_list <- append(AIC_list, AIC_value)
  cat("AIC: \n")
  print(AIC_value)
  
  cat("################################## Model 2  ##################################")
  # f(n) = an^b - a power law model
  a_initial = 1
  b_initial = 1
  model_2 = nls(mean_length~a*vertices^b,data=data,
                        start = list(a = a_initial, b = b_initial), trace = TRUE)
  models[["model_2"]] <- model_2
  # ---------------------------------------------------------------------------#
  cat("... with optimized initial values: \n")
  # Opt: better init values -> faster convergence
  linear_model = lm(log(mean_length)~log(vertices), data=data)
  a_init =  exp(coef(linear_model)[1])
  b_init = coef(linear_model)[2]
  nls(mean_length~a*vertices^b,data=data,
      start = list(a = a_init, b = b_init), trace = TRUE)
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  coef(model_2)
  coef(model_2)["a"]
  coef(model_2)["b"]
  # Errors of the non-linear regression model
  s <- sqrt(deviance(model_2)/df.residual(model_2))
  s_list <- append(s_list, s)
  cat("Residual Standard Error (s): \n")
  print(s)
  AIC_value <- AIC(model_2)
  AIC_list <- append(AIC_list, AIC_value)
  cat("AIC: \n")
  print(AIC_value)

  cat("################################## Model 3  ##################################")
  # f(n) = ae^(cn) - an exponential model
  a_initial = 1
  b_initial = 0
  model_3 = nls(mean_length~a*exp(b*vertices),data=data,
                        start = list(a = a_initial, b = b_initial), trace = TRUE)
  models[["model_3"]] <- model_3
  # ---------------------------------------------------------------------------#
  cat("... with optimized initial values: \n")
  # Opt: better init values -> faster convergence
  linear_model = lm(log(mean_length)~vertices, data=data)
  linear_model
  b_init = coef(linear_model)[2] # slope 
  a_init = exp(coef(linear_model)[1]) # e^intercept = a
  nls(mean_length~a*exp(b*vertices),data=data,
      start = list(a = a_init, b = b_init), trace = TRUE)
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  coef(model_3)
  coef(model_3)["a"]
  coef(model_3)["b"] # b is the same as c in this case
  # Errors of the non-linear regression model
  s <- sqrt(deviance(model_3)/df.residual(model_3))
  s_list <- append(s_list, s)
  cat("Residual Standard Error (s): \n")
  print(s)
  AIC_value <- AIC(model_3)
  AIC_list <- append(AIC_list, AIC_value)
  cat("AIC: \n")
  print(AIC_value)
  
  cat("################################## Model 1+  ##################################")
  # f(n) = (n/2)^b + d
  b_initial = 1
  d_initial = 1
  model_1s = nls(mean_length~(vertices/2)^b + d,data=data,
                        start = list(b = b_initial, d = d_initial), trace = TRUE)
  models[["model_1s"]] <- model_1s
  # ---------------------------------------------------------------------------#
  cat("... with optimized initial values: \n")
  # Opt: better init values -> faster convergence
  # take d_init as 1/2 of the min of the mean_length observations: 
  # good reference: 
  # https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
  d_init = min(data$mean_length) / 2
  linear_model = lm(log(mean_length - d_init)~log(vertices/2), data=data)
  linear_model
  b_init = coef(linear_model)[2]
  nls(mean_length~(vertices/2)^b + d,data=data,
      start = list(b = b_init, d = d_init), trace = TRUE)
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  coef(model_1s)
  coef(model_1s)["b"]
  coef(model_1s)["d"] 
  # Errors of the non-linear regression model
  s <- sqrt(deviance(model_1s)/df.residual(model_1s))
  s_list <- append(s_list, s)
  cat("Residual Standard Error (s): \n")
  print(s)
  AIC_value <- AIC(model_1s)
  AIC_list <- append(AIC_list, AIC_value)
  cat("AIC: \n")
  print(AIC_value)
  
  cat("################################## Model 2+  ##################################")
  # f(n) = ae^(cn) - an exponential model
  a_initial = 1
  b_initial = 1
  d_initial = 1
  model_2s = nls(mean_length~a*vertices^b + d,data=data,
                        start = list(a = a_initial, b = b_initial, d = d_initial),
                        trace = TRUE)
  models[["model_2s"]] <- model_2s
  # ---------------------------------------------------------------------------#
  # NOTE: SLOWER THAN WITH INITIAL = 1 VALUE OF PARAMETERS!!
  cat("... with optimized initial values: \n")
  # Opt: better init values -> faster convergence
  # take d_init as 1/2 of the min of the mean_length observations: 
  d_init = min(data$mean_length) / 2
  linear_model = lm(log(mean_length - d_init)~log(vertices), data=data)
  linear_model
  a_init = exp(coef(linear_model)[1]) # a = e^intercept
  b_init = coef(linear_model)[2] # b = slope
  nls(mean_length~a*vertices^b + d,data=data,
      start = list(a = a_init, b = b_init, d = d_init), trace = TRUE)
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  coef(model_2s)
  coef(model_2s)["a"]
  coef(model_2s)["b"]
  coef(model_2s)["d"]
  # Errors of the non-linear regression model
  s <- sqrt(deviance(model_2s)/df.residual(model_2s))
  s_list <- append(s_list, s)
  cat("Residual Standard Error (s): \n")
  print(s)
  AIC_value <- AIC(model_2s)
  AIC_list <- append(AIC_list, AIC_value)
  cat("AIC: \n")
  print(AIC_value)
  
  cat("################################## Null model ################################")
  # The RSS, s and AIC for a non-parametric model (such as the null model)
  # {d} = n/3 + 1/3
  RSS <- sum((data$mean_length-(data$vertices+1)/3)^2)
  n <- length(data$vertices)
  p <- 0
  s <- sqrt(RSS/(n - p))
  s
  AIC <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)
  AIC
  ########
  # Plot the empirical data and the curve for the best fit:
  # Best model - the one with lower error values
  modelNames <- names(models)
  
  ### wrt to s_value
  cat("Model with best s error: ", modelNames[which.min(s_list)], "\n")
  ## empirical data plot
  plot(log(data$vertices), log(data$mean_length),
       xlab = "log(vertices)", ylab = "log(mean dependency length)")
  ## best fit plot
  lines(log(data$vertices), log(fitted(models[[which.min(s_list)]])), col = "green")
  
  
  ###  wrt to AIC value
  cat("Model with best AIC: ", modelNames[which.min(AIC_list)], "\n")
  ## empirical data plot
  plot(log(data$vertices), log(data$mean_length),
       xlab = "log(vertices)", ylab = "log(mean dependency length)")
  ## best fit plot
  lines(log(data$vertices), log(fitted(models[[which.min(AIC_list)]])), col = "green")
  
}
