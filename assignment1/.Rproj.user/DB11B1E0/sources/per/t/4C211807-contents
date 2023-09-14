library(igraph)

p <- 10^(seq(-4,0,0.01))

N = 100

L <- rep(0,length(p))
C <- rep(0,length(p))

for (i in 1:length(p)) {
  for (j in 1:N) {
    g <- sample_smallworld(1,100,4,p)
    #g <- watts.strogatz.game(1,100,4,p[i])
    L[i] = L[i] + transitivity(g, "undirected")/N
    C[i] = C[i] + average.path.length(g)/N
  }
  
  #L[i] = L[i]/N
  #C[i] = C[i]/N
}

L <- L/L[1]
C <- C/C[1]

plot(p, L, ylim = c(0,1), ylab='coeff', xlab='p', pch=0, log='x')
points(p, C, ylim = c(0,1), ylab='coeff', xlab='p', pch=16)
legend("bottomleft", inset=.02,c('C(p)/C(0)','L(p)/L(0)'), pch=c(0,16))


#
#
#
#
#

num_nodes = x <- 10^(seq(0,6,0.2))
ASP <- c()
for( n in num_nodes ){
  ASP <- append(ASP, 0)
  for(i in 1:N) {
    g <- sample_gnp(n, (1+0.0001)*log(n)/n)
    ASP[length(ASP)] = ASP[length(ASP)] + average.path.length(g)
  }
  
  ASP[length(ASP)] = ASP[length(ASP)]/N
}

plot(num_nodes, ASP)
