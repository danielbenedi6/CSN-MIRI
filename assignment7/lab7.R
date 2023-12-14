library(igraph)
library(xtable)

simulate <- function(G, gamma, beta, p0) {
  curr <- rep(0., vcount(G))
  curr[sample(1:vcount(G), p0)] <- 1.
  prev <- rep(0, vcount(G))
  
  evolution <- c(sum(curr)/vcount(G))
  t <- 0
  
  prevNorm <- Inf
  while(t < 10 || abs(prevNorm - sqrt(sum((curr-prev)**2))) > 1e-3) {
    prevNorm <- sqrt(sum((curr-prev)**2))
    prev <- curr
    curr <- rep(0., vcount(G))
    
    for(node in 1:vcount(G)) {
      if(prev[node] == 1) {
        curr[node] = runif(1) > gamma
      } else {
        K <- 0
        for(neighbor in neighbors(G, node)) {
          if(prev[neighbor] == 1) {
            K <- K + 1
          }
        }
        curr[node] = runif(1) >  (1-beta)**K
      }
    }
    
    t <- t + 1
    evolution <- append(evolution, sum(curr)/vcount(G))
  }
  
  evolution
}

repetition <- function(num_rep, G, gamma, beta, p0) {
  sum_evols <- c(0,0)
  for(rep in 1:num_rep) {
    evol <- simulate(G, gamma, beta, p0)
    
    if(length(sum_evols) < length(evol)) {
      sum_evols <- append(sum_evols, rep(sum_evols[length(sum_evols)], length(evol)-length(sum_evols)))    
    } else if(length(sum_evols) > length(evol)) {
      evol <- append(evol, rep(evol[length(evol)], length(sum_evols)-length(evol)))  
    }
    
    sum_evols <- sum_evols + evol
  }
  
  sum_evols/num_rep
}

plot_erdos_renyi_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,1,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- erdos.renyi.game(N, p)
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs)
  title("Boxplot of the  Eigenalues of the Erdös-Renyi Graphs")
  data
}

plot_erdos_renyi <- function(N,step,num_rep, gamma) {
  plot.new()
  axis(1)
  axis(2)
  c <- 0
  for(p in seq(step,1,step)) {
    c <- c + 1
    G <- erdos.renyi.game(N, p)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta = gamma / max_eigenvalue + 0.05
    res <- repetition(num_rep, G, gamma, beta, 100)
    lines(seq(0,1,1/(length(res)-1)),res, col=c)
  }
  title("SIS Simulation on different Erdös-Renyi Graph varying probability")
  legend("bottomright",legend=seq(step,1,step),col=1:c, lty=1)
}

N <-1000
step <- 0.2
num_rep <- 10
gamma <- 0.3

# This experiment  allows us to see how SIS converges  depending  on  the
# probabiblity of Erdös-Renyi graphs
plot_erdos_renyi(N,step,num_rep,gamma)

# This experiment  allows us to see how much variance there is in the eigenvaleus
# for different Erdös-Renyi graphs.
plot_erdos_renyi(N,step,num_rep,gamma)


