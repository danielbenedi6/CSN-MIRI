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

# 2. Ensemble of models





