library(igraph)

p <- 10^(seq(-4,0,0.2))

N = 100
NODES = 1000

#for (neigh in c(2,3,4,5) ){
for (neigh in c(2,3,4,5,6,8,10,20) ){
  L <- rep(0,length(p))
  C <- rep(0,length(p))
  
  for (i in 1:length(p)) {
    for (j in 1:N) {
      # Note that the smaller the population and the larger the neighbourhood, then the clustering coefficient will be higher
      g <- sample_smallworld(1,NODES,neigh,p[i])
      L[i] = L[i] + transitivity(g, "global")
      C[i] = C[i] + mean_distance(g)
    }
  
    L[i] = L[i]/N
    C[i] = C[i]/N
  }

  L <- L/L[1]
  C <- C/C[1]
  
  
  pdf(paste("./plots/WS_neigh_",neigh,".pdf", sep=""))
  plot(p, L, ylim = c(0,1), ylab='coeff', xlab='p', pch=0, log='x', , main= paste("Transitivity and Mean distance for WS with ",NODES," nodes and neigh = ",neigh,sep=""))
  points(p, C, ylim = c(0,1), ylab='coeff', xlab='p', pch=16)
  legend("bottomleft", inset=.02,c('C(p)/C(0)','L(p)/L(0)'), pch=c(0,16))
  dev.off()
}

#
#
#
#
#

num_nodes = round(10^(seq(0.1,5,0.2)))
apl_gnp <- function(n) mean_distance(sample_gnp(n, (1+0.0001)*log(n)/n))
mean_apl <- function(n) sum(rep(apl_gnp(n), N))/N

ASP <- Map(mean_apl,num_nodes)

ASP <- rep(0,length(num_nodes))
for ( i in 1:length(num_nodes) ){
  for (j in 1:N) {
    g <- sample_gnp(num_nodes[i], (1+0.001)*log(num_nodes[i])/num_nodes[i])
    ASP[i] = ASP[i] + average.path.length(g)
  }
  
  ASP[i] <- ASP[i]/N
  print(i)
}

plot(num_nodes, ASP)
