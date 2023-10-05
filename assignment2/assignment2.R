languages = c("Arabic", "Basque", "Catalan",
              "Chinese", "Czech", "English",
              "Greek", "Hungarian", "Italian",
              "Turkish")

## PART 1

write_summary <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  #cat(paste(language,length(degree_sequence$V1),max(degree_sequence$V1),sum(degree_sequence$V1)/length(degree_sequence$V1),length(degree_sequence$V1)/sum(degree_sequence$V1), sep=" & "), "\\\\ \n")
  cat(sprintf(r"(%s & %d & %d & %.3f & %.6f \\ )",
              language,
              length(degree_sequence$V1),
              max(degree_sequence$V1),
              sum(degree_sequence$V1)/length(degree_sequence$V1),
              length(degree_sequence$V1)/sum(degree_sequence$V1)), sep="\n")
}

write_dataset_table <-function(languages, suffix){
  cat("\\begin{table}\n\\centering\n\\begin{tabular}{c c c c c}\n")
  cat("Language & $N$ & Maximum degree & $M/N$ & $N/M$ \\\\ \n")
  cat("\\hline  \\\\ \n")
  for(language in languages) {
    write_summary(language, paste("./data/",language,suffix,sep=""))
  }
  cat("\\end{tabular}\n\\caption{Summary of dataset}\\label{tab:data}\n\\end{table}\n")
}

write_dataset_table(languages, "_in-degree_sequence.txt")

## PART 2

plot_spectrum <- function(languages, suffix) {
  for(language in languages) {
    degree_sequence = read.table(paste("./data/",language,suffix,sep=""), header = FALSE)
    degree_spectrum = table(degree_sequence)
    
    ds = data.frame( degree_spectrum )
    for( j in 1:length(degree_spectrum) ){
      #Yes, this is weird, but how R manages factors is weirder
      ds$V2[j] = as.numeric(as.character(ds$V1[j])) 
    }
    ds = ds[c("V2","Freq")]
    ds$p = ds$Freq/sum( degree_spectrum )
    colnames(ds) = c("k","f","p")
    
    pdf(paste("./plot/",language,"_spectrum.pdf", sep=""))
    barplot(degree_spectrum, main = language, xlab = "Degree", ylab = "Number of Vertices", cex.lab=1.6,  cex.main=1.8, cex.names=1.4)
    dev.off()
    pdf(paste("./plot/",language,"_spectrum_loglog.pdf", sep = ""))
    matplot(ds$k,ds$p, main = language, xlab = "Degree", ylab = "Number of Vertices", log="xy" , pch = 19, cex.lab=1.6, cex.main=1.8, cex.names=1.4)
    dev.off()
  }
}


plot_spectrum(languages, "_in-degree_sequence.txt")

## PART 3
require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function
# require("plyr")

displaced_poisson_likelihood <- function(M, Mprime, N, C, kmax, lambda) {
  M*log(lambda) - N*(lambda+log(1-exp(-lambda)))-C
}

displaced_geometric_likelihood <- function(M, Mprime, N, C, kmax, q) {
  (M-N)*log(1-q) + N*log(q)
}

zeta_likelihood <- function(M, Mprime, N, C, kmax, gamma) {
  - gamma * M - N*log(zeta(gamma)) 
}

right_truncated_zeta_likelihood <- function(M, Mprime, N, C, kmax, gamma) {
  -gamma * M - N*log(sum((1:kmax)^(-gamma)))
}

likelihoods = c(displaced_poisson_likelihood, displacedgeometriclikelihood, zeta_likelihood, right_truncated_zeta_likelihood)


for(language in languages) {
  degree_sequence = read.table(paste("./data/",language,suffix,sep=""), header = FALSE)
  
  for(likelihood in likelihoods) {
    func <- function(param) {
      -likelihood(
        sum(degree_sequence),
        sum(log(degree_sequence)),
        nrow(degree_sequence),
        sum(laply(degree_sequence$V1, function(k){sum(log(1:k))})),
        max(degree_sequence),
        param
      )
    }
    mle_likelihood <- mle(func,
                          start = list(param = 2),
                          method = "L-BFGS-B",
                          lower=c(1.000000001))
    param <- attributes(mle_likelihood)$coef[1]
    
  } 
}