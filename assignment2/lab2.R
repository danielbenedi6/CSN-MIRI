#lab 2

install.packages("stats4")
install.packages("VGAM")
install.packages("latex2exp")

require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function
require("latex2exp") #to use latex in legend
require("plyr") # for the laply function

#Write the probability distribution functions
displazed_poisson <- function( k , lambda ){
  lambda**k * exp(-lambda) / factorial(k) / ( 1 - exp(-lambda))
}

displazed_geometric <- function( k , q ){
  (1-q)**k*q
}

zeta_2 <- function( k ){
  k**(-2)/zeta(2)
}

zeta_g <- function( k , gamma ){
  k**(-gamma)/zeta(gamma)
}

right_truncated_zeta <- function( k , gamma , km ){
  k**(-gamma)/sum( (1:floor(km)) ** (-gamma) )
}

altmann <- function( k , n , gamma , delta ){
  c = 1 / sum( (1:n) **(-gamma)*exp(-(1:n)*delta ) )
  c*k**(-gamma)*exp(-delta*k)
}

# Fit each probability distribution for each language

languages = c("Arabic", "Basque", "Catalan",
              "Chinese", "Czech", "English",
              "Greek", "Hungarian", "Italian",
              "Turkish")

table <- r"(\begin{table}[!htb]
\centering
\begin{tabular}{llllllll}
\multicolumn{1}{r}{\textbf{}} & \multicolumn{5}{c}{\textbf{Model}} \\ \cline{2-8} 
\multicolumn{1}{c}{\textbf{}} &
  \multicolumn{1}{c}{D. Poisson} &
  \multicolumn{1}{c}{D. Geom} &
  \multicolumn{1}{c}{Zeta} &
  \multicolumn{2}{c}{Truncated Z} &
  \multicolumn{2}{c}{Altmann} \\
\multicolumn{1}{c}{\textbf{Language}} &
  \multicolumn{1}{c}{$\lambda$} &
  \multicolumn{1}{c}{$q$} &
  \multicolumn{1}{c}{$\gamma_1$} &
  \multicolumn{1}{c}{$\gamma_2$} &
  \multicolumn{1}{c}{$k_{max}$}  &
  \multicolumn{1}{c}{$\gamma_3$} &
  \multicolumn{1}{c}{$\delta$} \\ \hline)"

table_aic <- r"(\begin{table}[!htb]
\centering
\resizebox{\columnwidth}{!}{
\begin{tabular}{lllllll}
\multicolumn{1}{r}{\textbf{}} & \multicolumn{5}{c|}{\textbf{Model}}    &  \\ \cline{2-6}
\multicolumn{1}{c}{\textbf{Language}} &
  \multicolumn{1}{c}{D. Poisson} &
  \multicolumn{1}{c}{D. Geom} &
  \multicolumn{1}{c}{Zeta} &
  \multicolumn{1}{c}{Truncated Z} &
  \multicolumn{1}{l|}{Altmann} &
  $AIC_{best}$ \\ \hline)"

for(language in languages) {
  degree_sequence = read.table(paste("./data/",language,"_in-degree_sequence.txt",sep=""), header = FALSE)
  degree_spectrum = table(degree_sequence)
  
  ds = data.frame( degree_spectrum )
  for( j in 1:length(degree_spectrum) ){
    #Yes, this is weird, but how R manages factors is weirder
    ds$V2[j] = as.numeric(as.character(ds$V1[j])) 
  }
  ds = ds[c("V2","Freq")]
  ds$p = ds$Freq/sum( degree_spectrum )
  colnames(ds) = c("k","f","p")
  
  #Get the graph information
  M = sum( degree_sequence )
  N = sum( degree_spectrum )
  kmax = max(ds$k)
  p = degree_spectrum / N
  C = 0
  for( k in degree_sequence$V1 ){
    if( k > 1 ) { C = C + sum( log(1:k)) }
  }
  Mp = sum( log(degree_sequence) )

  
  # Write the likelihoods
  minus_log_likelihood_displazed_poisson <- function( lambda ){
    -M*log(lambda) + N*(lambda+log(1-exp(-lambda))) + C
  }
  minus_log_likelihood_displazed_geometric <- function( q ){
    (N-M)*log( 1 - q ) - N*log(q) 
  }
  minus_log_likelihood_zeta_2 <- function( ){
    2*Mp + N*log(pi**2/6)
  }
  minus_log_likelihood_zeta <- function(gamma) {
    gamma*Mp + N*log( zeta(gamma) )
  }
  minus_log_likelihood_right_truncated_zeta <- function( gamma , km ){
    gamma * Mp + N * log( sum( (1:km) ** (-gamma)) )
  }
  minus_log_likelihood_altmann <- function( gamma , delta ){
    d = sum( (1:N) **(-gamma)*exp(-(1:N)*delta ) )
    N*log(d) + gamma*Mp + delta*M
  }
  
  
  # Define function to compute the AIC of a model
  get_AIC <- function(m2logL,K,N) {
    if(N - K > 1000){
      m2logL + 2*K  # AIC without correction when sample is big enough
    } else {
      m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
    }
  }
  
  #Get all the mle parameters
  mle_displazed_poisson <- mle( minus_log_likelihood_displazed_poisson , 
                                start = list( lambda = 1) , method = "L-BFGS-B" , 
                                lower = c(0.00000001) )
  mle_dp_lambda <- attributes(summary(mle_displazed_poisson))$coef[1]
  mle_dp_m2logL <- attributes(summary(mle_displazed_poisson))$m2logL
  mle_dp_aic <- get_AIC(mle_dp_m2logL,1,N)
  
  mle_displazed_geometric <- mle( minus_log_likelihood_displazed_geometric , 
                                  start = list( q = 0.5) , method = "L-BFGS-B" , 
                                  lower = c(0.00000001) , upper = c(0.99999999))
  mle_dg_q <- attributes(summary(mle_displazed_geometric))$coef[1]
  mle_dg_m2logL <- attributes(summary(mle_displazed_geometric))$m2logL
  mle_dg_aic <- get_AIC(mle_dg_m2logL,1,N)
  
  mle_displazed_zeta_g <- mle( minus_log_likelihood_zeta ,
                               start = list( gamma = 2 ) , method = "L-BFGS-B" , 
                               lower = c(1.00000001) )
  mle_dz_gamma <- attributes(summary(mle_displazed_zeta_g))$coef[1]
  mle_dz_m2logL <- attributes(summary(mle_displazed_zeta_g))$m2logL
  mle_dz_aic <- get_AIC(mle_dz_m2logL,1,N)
  
  #mle_drtz_km <- ds$k[ds$f == 1][1]

  mle_drtz_m2logL_prev <-  1e7
  mle_drtz_m2logL <-  1e6
  mle_drtz_km  <- kmax
  while(abs(mle_drtz_m2logL-mle_drtz_m2logL_prev)>1e-6){
    mle_drtz_km<-mle_drtz_km+round(0.05*kmax)
    likelihood <- function(gamma) {
      minus_log_likelihood_right_truncated_zeta(gamma, mle_drtz_km)
    }
    mle_displazed_right_truncated_zeta <- mle( likelihood ,
                                               start = list( gamma = 2 ) , method = "L-BFGS-B" , 
                                               lower = c( 0.00000001 ) , upper = c( 100 ) )
    mle_drtz_m2logL <- attributes(summary(mle_displazed_right_truncated_zeta))$m2logL
    mle_drtz_gamma <- attributes(summary(mle_displazed_right_truncated_zeta))$coef[1]
    mle_drtz_m2logL_prev <- mle_drtz_m2logL
  }
  mle_drtz_aic <- get_AIC(mle_drtz_m2logL,2,N)
  
  # TODO: Find proper values to bounds and start delta and gamma
  mle_altmann <- mle( minus_log_likelihood_altmann ,
                      start = list(gamma = 2 , delta = 2 ), 
                      method = "L-BFGS-B" ,
                      lower = c( 0.00000001 , 0.0000001) , 
                      upper = c( 100 , 100) )
  mle_altm_gamma <- attributes(summary(mle_altmann))$coef[1]
  mle_altm_delta <- attributes(summary(mle_altmann))$coef[2]
  mle_altm_m2logL <- attributes(summary(mle_altmann))$m2logL
  mle_altm_aic <- get_AIC(mle_altm_m2logL,2,N)
  
  #plot all the lines
  pdf(paste("./plot/",language,"_fit.pdf",sep=""))
  matplot(ds$k,ds$p, main = language,
          xlab = "degree", ylab = "number of vertices" , log = "xy" , pch = 19 , cex.lab=1.6,  cex.main=1.8, cex.names=1.4 )
  
  lines( 1:kmax , displazed_poisson( 1:kmax , mle_dp_lambda) , col = 2 )
  lines( 1:kmax , displazed_geometric( 1:kmax , mle_dg_q) , col = 3)
  lines( 1:kmax , zeta_2( 1:kmax) , col = 4 )
  lines( 1:kmax , zeta_g( 1:kmax , mle_dz_gamma) , col = 5)
  lines( 1:kmax , right_truncated_zeta( 1:kmax , mle_drtz_gamma, mle_drtz_km) , col = 6)
  lines( 1:kmax , altmann(1:kmax , N , mle_altm_gamma , mle_altm_delta ) , col = 7)
  legend( 10**(0.45*(log10(max(ds$k))+log10(min(ds$k)))) , max(ds$p) , 
          legend = c(TeX(sprintf(r'(D. Poisson($\lambda = %0.3f$))', mle_dp_lambda)), 
                     TeX(sprintf(r'(D. Geom($\q = %0.3f$))', mle_dg_q)), 
                     "Zeta(2)", 
                     TeX(sprintf(r'(Zeta($\gamma = %0.3f$))', mle_dz_gamma)),
                     TeX(sprintf(r'(Trunc Zeta($\gamma = %0.3f, k_{max} = %d$))', mle_drtz_gamma, mle_drtz_km)),
                     TeX(sprintf(r'(Altmann($\gamma = %0.3f, \delta = %0.3f$))', mle_altm_gamma, mle_altm_delta))
          ), 
          col =2:7 , cex = 1 , title = "Fit lines" , lty = rep(1,4) 
        )  
  dev.off()
  
  # Write resume table with the fittted parameters
  table <- paste(table,sprintf(r"(\multicolumn{1}{c}{%s} & %.3f & %.3f & %.3f & %.3f & %d  & %.3f & %.3f  \\)",language,mle_dp_lambda,mle_dg_q,mle_dz_gamma,mle_drtz_gamma,mle_drtz_km,mle_altm_gamma,mle_altm_delta),sep="\n")
  
  # Obtain the best AIC and compute the AIC difference
  AIC_best <- min(mle_dp_aic, mle_dg_aic, mle_dz_aic, mle_drtz_aic, mle_altm_aic)
  table_aic <- paste(table_aic,sprintf(r"(\multicolumn{1}{c}{%s} & %.3f & %.3f & %.3f & %.3f & \multicolumn{1}{l|}{%.3f} & %.3f \\)",language,
                                    mle_dp_aic - AIC_best,
                                    mle_dg_aic - AIC_best,
                                    mle_dz_aic - AIC_best,
                                    mle_drtz_aic - AIC_best,
                                    mle_altm_aic -  AIC_best,
                                    AIC_best),
                      sep="\n")
  
}

table <- paste(table,r"(\end{tabular}
\caption{Values of fitted parameters per distribution and language} \label{tab:fitted}
\end{table})",sep="\n")
table_aic <- paste(table_aic,r"(\end{tabular}}
\caption{AIC difference of the models per language} \label{tab:aic}
\end{table})",sep="\n")

cat(table)
cat(table_aic)



