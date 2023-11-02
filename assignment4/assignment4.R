library(ggplot2)
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
# 
for(language in languages){
  cat(language)
  data = read_dataset(language)
  
  # First condition 4-6-n <= <k²> <= n - 1
  data = subset(subset(data, 4-6/vertices <= degree_2nd_moment),degree_2nd_moment <= vertices - 1)
  # Second condition n/(8(n-1))<k²> +1/2 <= <d> <= n - 1
  data = subset(subset(data, vertices/(8*(vertices-1)) * degree_2nd_moment + 1/2 <= mean_length),mean_length <= vertices - 1)
  
  table_summary <- paste(table_summary,generate_table_entry(language,data), sep="\n")
  
  if(bptest(lm(mean_length ~ vertices, data=data))$p.value < 0.05) {
    cat(" Heteroscedasticity\n")
    data = aggregate(data, list(data$vertices), mean)
  }  else {
    cat(" Homoscedasticity\n")
  }
  
  #DO FIT FOR EACH MODEL
  
  model_fit_1 = nls(mean_length ~ (vertices/2)^b, start = list(b = 1),data=data)
  param_b_model1 = coef(model_fit_1)[1]
  deviance_model1 = deviance(model_fit_1)
  aic_model1 = AIC(model_fit_1)
  s_model1 = sqrt(deviance_model1/df.residual(model_fit_1))
  
  model_fit_2 = nls(mean_length ~ a*vertices^b, start = list(a = 1, b = 1),data=data)
  param_a_model2 = coef(model_fit_2)[1]
  param_b_model2 = coef(model_fit_2)[2]
  deviance_model2 = deviance(model_fit_2)
  aic_model2 = AIC(model_fit_2)
  s_model2 = sqrt(deviance_model2/df.residual(model_fit_2))
  
  plt <- ggplot() + geom_point(data=dataset, aes(x=vertices,y=mean_length))
  plt <- plt + geom_line(data=means, aes(x=vertices,y=mean_length), color="green")
  plt <- plt + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
  plt <- plt + labs(title=paste("Models for mean length dependency for ", language))
  plt <- plt + labs(colour="Model")
  plt <- plt + xlab("Number of vertices")
  plt <- plt + ylab("Mean length")
  
  next;
  
  table_res_se <- paste(table_res_se, print_row_type2(language,
                                                      0., # Model 0
                                                      s_model1, # Model 1
                                                      s_model2, # Model 2
                                                      0., # Model 3
                                                      0., # Model 4
                                                      0., # Model 1+
                                                      0., # Model 2+
                                                      0., # Model 3+
                                                      0.  # Model 4+
  ), sep="\n")
  table_aic <- paste(table_aic, print_row_type2(language,
                                                      0., # Model 0
                                                      aic_model1, # Model 1
                                                      aic_model2, # Model 2
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