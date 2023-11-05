library(stats)
library(lmtest)

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
  cat("---------------------------------------------------------------------\n")
  cat("Language: ", language, "\n")
  dataset = read_dataset(language)
  data <- dataset
  
  # First condition 4-6-n <= <k²> <= n - 1
  data = subset(subset(data, 4-6/vertices <= degree_2nd_moment),degree_2nd_moment <= vertices - 1)
  # Second condition n/(8(n-1))<k²> +1/2 <= <d> <= n - 1
  data = subset(subset(data, vertices/(8*(vertices-1)) * degree_2nd_moment + 1/2 <= mean_length),mean_length <= vertices - 1)
  
  table_summary <- paste(table_summary,generate_table_entry(language,data), sep="\n")
  
  if(bptest(lm(mean_length ~ vertices, data=data))$p.value < 0.05) {
    cat("    Heteroscedasticity\n")
    data = aggregate(data, list(data$vertices), mean)
  }  else {
    cat("    Homoscedasticity\n")
  }
  models <- list()
  s_list <- list()
  AIC_list <- list()
  # Fit of the models
  cat("################################## Model 1  ##################################\n")
  # f(n) = (n/2)^b
  linear_model = lm(log(mean_length)~log(vertices/2), data=data)
  b_init = coef(linear_model)[2]
  model_1 = nls(mean_length~(vertices/2)^b,data=data,
                start = list(b = b_init), trace = FALSE)
  models[["model_1"]] <- model_1
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_b_model1 = coef(model_1)["b"]
  # Errors of the non-linear regression model
  s_model1 <- sqrt(deviance(model_1)/df.residual(model_1))
  s_list <- append(s_list, s_model1)
  cat("Residual Standard Error (s): \n")
  print(s_model1)
  aic_model1 <- AIC(model_1)
  AIC_list <- append(AIC_list, aic_model1)
  cat("AIC: \n")
  print(aic_model1)
  
  cat("################################## Model 2  ##################################\n")
  # f(n) = an^b - a power law model
  linear_model = lm(log(mean_length)~log(vertices), data=data)
  a_init =  exp(coef(linear_model)[1])
  b_init = coef(linear_model)[2]
  model_2 = nls(mean_length~a*vertices^b,data=data,
                start = list(a = a_init, b = b_init), trace = FALSE)
  models[["model_2"]] <- model_2
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_a_model2 = coef(model_2)["a"]
  param_b_model2 = coef(model_2)["b"]
  # Errors of the non-linear regression model
  s_model2 <- sqrt(deviance(model_2)/df.residual(model_2))
  s_list <- append(s_list, s_model2)
  cat("Residual Standard Error (s): \n")
  print(s_model2)
  aic_model2 <- AIC(model_2)
  AIC_list <- append(AIC_list, aic_model2)
  cat("AIC: \n")
  print(aic_model2)
  
  cat("################################## Model 3  ##################################\n")
  # f(n) = ae^(cn) - an exponential model
  linear_model = lm(log(mean_length)~vertices, data=data)
  c_init = coef(linear_model)[2] # slope 
  a_init = exp(coef(linear_model)[1]) # e^intercept = a
  model_3 = nls(mean_length~a*exp(c*vertices),data=data,
                start = list(a = a_init, c = c_init), trace = FALSE)
  models[["model_3"]] <- model_3
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_a_model3 <- coef(model_3)["a"]
  param_c_model3 <- coef(model_3)["c"]
  # Errors of the non-linear regression model
  s_model3 <- sqrt(deviance(model_3)/df.residual(model_3))
  s_list <- append(s_list, s_model3)
  cat("Residual Standard Error (s): \n")
  print(s_model3)
  aic_model3 <- AIC(model_3)
  AIC_list <- append(AIC_list, aic_model3)
  cat("AIC: \n")
  print(aic_model3)
  
  cat("################################## Model 4  ##################################\n")
  # f(n) = a log(n) - a logarithmic model
  linear_model = lm(log(mean_length)~vertices, data=data)
  a_init = exp(coef(linear_model)[1]) # e^intercept = a
  model_4 = nls(mean_length~a*log(vertices),data=data,
                start = list(a = a_init), trace = FALSE)
  models[["model_4"]] <- model_4
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_a_model4 <- coef(model_4)["a"]
  # Errors of the non-linear regression model
  s_model4 <- sqrt(deviance(model_4)/df.residual(model_4))
  s_list <- append(s_list, s_model4)
  cat("Residual Standard Error (s): \n")
  print(s_model4)
  aic_model4 <- AIC(model_4)
  AIC_list <- append(AIC_list, aic_model4)
  cat("AIC: \n")
  print(aic_model4)
  
  cat("################################## Model 1+  ##################################\n")
  # f(n) = (n/2)^b + d
  # Take d_init as something smaller than the min of the mean_length observations: 
  # good reference: 
  # https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
  d_init = min(data$mean_length) - 0.001
  linear_model = lm(log(mean_length - d_init)~log(vertices/2), data=data)
  b_init = coef(linear_model)[2]
  model_1p <- nls(mean_length~(vertices/2)^b + d,data=data,
                  start = list(b = b_init, d = d_init), trace = FALSE)
  models[["model_1p"]] <- model_3
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_b_model1p <- coef(model_1p)["b"]
  param_d_model1p <- coef(model_1p)["d"]
  # Errors of the non-linear regression model
  s_model1p <- sqrt(deviance(model_1p)/df.residual(model_1p))
  s_list <- append(s_list, s_model1p)
  cat("Residual Standard Error (s): \n")
  print(s_model1p)
  aic_model1p <- AIC(model_1p)
  AIC_list <- append(AIC_list, aic_model1p)
  cat("AIC: \n")
  print(aic_model1p)
  
  
  cat("################################## Model 2+  ##################################\n")
  # f(n) = an^b + d
  d_init = min(data$mean_length) - 0.001
  linear_model = lm(log(mean_length - d_init)~log(vertices), data=data)
  a_init = exp(coef(linear_model)[1]) # a = e^intercept
  b_init = coef(linear_model)[2] # b = slope
  model_2p = nls(mean_length~a*vertices^b + d,data=data,
                 start = list(a = a_init, b = b_init, d = d_init), trace = FALSE)
  models[["model_2p"]] <- model_2p
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_a_model2p = coef(model_2)["a"]
  param_b_model2p = coef(model_2)["b"]
  param_d_model2p = coef(model_2)["d"]
  # Errors of the non-linear regression model
  s_model2p <- sqrt(deviance(model_2p)/df.residual(model_2p))
  s_list <- append(s_list, s_model2p)
  cat("Residual Standard Error (s): \n")
  print(s_model2p)
  aic_model2p <- AIC(model_2p)
  AIC_list <- append(AIC_list, aic_model2p)
  cat("AIC: \n")
  print(aic_model2p)
  
  cat("################################## Model 3+  ##################################\n")
  # f(n) = ae^(cn) + d - an exponential model
  d_init = min(data$mean_length) - 0.001
  linear_model = lm(log(mean_length - d_init)~vertices, data=data)
  c_init = coef(linear_model)[2] # slope 
  a_init = exp(coef(linear_model)[1]) # e^intercept = a
  model_3p = nls(mean_length~a*exp(c*vertices)+d,data=data,
                start = list(a = a_init, c = c_init, d = d_init), trace = FALSE,control=nls.control(maxiter=100))
  models[["model_3p"]] <- model_3p
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_a_model3p <- coef(model_3p)["a"]
  param_c_model3p <- coef(model_3p)["c"]
  param_d_model3p <- coef(model_3p)["d"]
  # Errors of the non-linear regression model
  s_model3p <- sqrt(deviance(model_3p)/df.residual(model_3p))
  s_list <- append(s_list, s_model3p)
  cat("Residual Standard Error (s): \n")
  print(s_model3p)
  aic_model3p <- AIC(model_3p)
  AIC_list <- append(AIC_list, aic_model3p)
  cat("AIC: \n")
  print(aic_model3p)
  
  cat("################################## Model 4+  ##################################\n")
  # f(n) = a log(n) - a logarithmic model
  d_init = min(data$mean_length) - 0.001
  linear_model = lm(log(mean_length - d_init)~vertices, data=data)
  d_init = coef(linear_model)[2] # slope 
  a_init = exp(coef(linear_model)[1]) # e^intercept = a
  model_4p = nls(mean_length~a*log(vertices)+d,data=data,
                start = list(a = a_init, d = d_init), trace = FALSE)
  models[["model_4p"]] <- model_4p
  # ---------------------------------------------------------------------------# 
  # Params giving the best fit for the model
  param_a_model4p <- coef(model_4p)["a"]
  param_d_model4p <- coef(model_4p)["d"]
  # Errors of the non-linear regression model
  s_model4p <- sqrt(deviance(model_4p)/df.residual(model_4p))
  s_list <- append(s_list, s_model4p)
  cat("Residual Standard Error (s): \n")
  print(s_model4p)
  aic_model4p <- AIC(model_4p)
  AIC_list <- append(AIC_list, aic_model4p)
  cat("AIC: \n")
  print(aic_model4p)
  
  cat("################################## Null model ################################\n")
  # The RSS, s and AIC for a non-parametric model (such as the null model)
  # {d} = n/3 + 1/3
  model0 <- mean_length ~ ((vertices+1)/3)^2
  RSS <- sum((data$mean_length-(data$vertices+1)/3)^2)
  n <- length(data$vertices)
  p <- 0
  s_model0 <- sqrt(RSS/(n - p))
  s_list <- append(s_list, s_model0)
  cat("Residual Standard Error (s): \n")
  print(s_model0)
  aic_model0 <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)
  AIC_list <- append(AIC_list, aic_model0)
  cat("AIC: \n")
  print(aic_model0)
  
  ########
  # Plot the empirical data and the curve for the best fit:
  # Best model - the one with lower error values
  modelNames <- names(models)
  
  ### wrt to s_value
  cat("Model with best s error: ", modelNames[which.min(s_list)], "\n")
  ## empirical data plot
  pdf(paste("./plot/",language,"_best_s.pdf", sep=""))
  plot(data$vertices, data$mean_length, log="xy",
       xlab = "Vertices", ylab = "Mean dependency length", 
       main = paste("Best s for", language, modelNames[which.min(s_list)], sep=" "))
  ## best fit plot
  lines(data$vertices, fitted(models[[which.min(s_list)]]), col = "green")
  dev.off()
  
  
  ###  wrt to AIC value
  cat("Model with best AIC: ", modelNames[which.min(AIC_list)], "\n")
  ## empirical data plot
  pdf(paste("./plot/",language,"_best_aic.pdf", sep=""))
  plot(data$vertices, (data$mean_length), log="xy",
       xlab = "Vertices", ylab = "Mean dependency length", 
       main = paste("Best AIC for", language, modelNames[which.min(AIC_list)], sep=" "))
  ## best fit plot
  lines(log(data$vertices), log(fitted(models[[which.min(AIC_list)]])), col = "green")
  dev.off()
  
  #############################################
  #############################################
  #############################################
  
  pdf(paste("./plot/",language,"_all_models.pdf", sep=""))
  #plot(dataset$vertices, dataset$mean_length, ylab="Mean Length", xlab="Vertices", log="xy")
  plot(data$vertices, data$mean_length, ylab="Mean Length", xlab="Vertices", log="xy")
  #lines(data$vertices, data$mean_length, type="l", lty=2, col="green")
  lines(data$vertices, (data$vertices + 1)^2/9, type="l", lty=1, col="grey")
  lines(data$vertices, fitted(model_1), type="l", lty=1, col=1)
  lines(data$vertices, fitted(model_2), type="l", lty=1, col=2)
  lines(data$vertices, fitted(model_3), type="l", lty=1, col=3)
  lines(data$vertices, fitted(model_4), type="l", lty=1, col=4)
  lines(data$vertices, fitted(model_1p), type="l", lty=1, col=5)
  lines(data$vertices, fitted(model_2p), type="l", lty=1, col=6)
  lines(data$vertices, fitted(model_3p), type="l", lty=1, col=7)
  lines(data$vertices, fitted(model_4p), type="l", lty=1, col=8)
  legend("topleft",legend = c("Null Hypothesis",
                              "Model 1",
                              "Model 2",
                              "Model 3",
                              "Model 4",
                              "Model 1+",
                              "Model 2+",
                              "Model 3+",
                              "Model 4+"
                              ),
         col=append(c("grey"),1:8), lty=1
        )
  title(main=paste("Models for mean length dependency for ", language))
  dev.off()
  
  table_res_se <- paste(table_res_se, print_row_type2(language,
                                                      s_model0, # Model 0
                                                      s_model1, # Model 1
                                                      s_model2, # Model 2
                                                      s_model3, # Model 3
                                                      s_model4, # Model 4
                                                      s_model1p, # Model 1+
                                                      s_model2p, # Model 2+
                                                      s_model3p, # Model 3+
                                                      s_model4p  # Model 4+
  ), sep="\n")
  table_aic <- paste(table_aic, print_row_type2(language,
                                                      aic_model0, # Model 0
                                                      aic_model1, # Model 1
                                                      aic_model2, # Model 2
                                                      aic_model3, # Model 3
                                                      aic_model4, # Model 4
                                                      aic_model1p, # Model 1+
                                                      aic_model2p, # Model 2+
                                                      aic_model3p, # Model 3+
                                                      aic_model4p  # Model 4+
  ), sep="\n")
  best_AIC <- AIC_list[[which.min(AIC_list)]]
  table_aic_diff <- paste(table_aic_diff, print_row_type2(language,
                                                          best_AIC - aic_model0, # Model 0
                                                          best_AIC - aic_model1, # Model 1
                                                          best_AIC - aic_model2, # Model 2
                                                          best_AIC - aic_model3, # Model 3
                                                          best_AIC - aic_model4, # Model 4
                                                          best_AIC - aic_model1p, # Model 1+
                                                          best_AIC - aic_model2p, # Model 2+
                                                          best_AIC - aic_model3p, # Model 3+
                                                          best_AIC - aic_model4p  # Model 4+
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

write(table_summary, file="./table/summary.tex")
write(table_res_se, file="./table/res_se.tex")
write(table_aic, file="./table/aic.tex")
write(table_aic_diff, file="./table/aic_diff.tex")
write(table_params, file="./table/params.tex")

# residual standard error 
# - measures the standard deviation of the residuals in a regression model
# - the smaller the residual standard error, the better a regression model fits 
#   a dataset
# - sqrt(RSS / df) -> sum of squared residuals (RSS) / degrees of freedom 
#   (= # of observations - # of model parameters)

