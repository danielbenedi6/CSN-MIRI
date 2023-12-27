library(igraph)
library(clustAnalytics)
library(xtable)
library(ggplot2)

simulate <- function(G, gamma, beta, p0) {
  curr <- rep(0., vcount(G))
  curr[sample(1:vcount(G), p0)] <- 1.
  prev <- rep(0, vcount(G))
  
  evolution <- c(sum(curr)/vcount(G))
  t <- 0
  
  prevNorm <- Inf
  while(t < 10 || abs(prevNorm - sqrt(sum((curr-prev)**2))) > 1e-1) {
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


#####ERDOS-RENYI FUNCTIONS######
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
  title("Boxplot of the  Eigenalues of the Erdös-Renyi Graph")
  data
}

plot_erdos_renyi <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 1, step)) {
    G <- erdos.renyi.game(N, p)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p, x = seq(0, 1, 1 / (length(res) - 1)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes') +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Erdös-Renyi Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Spread Probability')
}



####BARABASI-ALBERT FUNCTIONS######
BarabasiAlbert <- function(N,pb) {
  B <- matrix(c(1, 0.2, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 0.2, 1), ncol=4)
  barabasi_albert_blocks(m=8, p=c(pb, pb, pb, pb), B=B, t_max=N, type="Hajek", sample_with_replacement = FALSE)
}


plot_barabasi_albert_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,1,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- BarabasiAlbert(N, p)
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs)
  title("Boxplot of the  Eigenalues of the Barabasi-Albert Graph")
  data
}



plot_barabasi_albert <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 1, step)) {
    G <- BarabasiAlbert(N, p)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p, x = seq(0, 1, 1 / (length(res) - 1)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes') +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Barabasi-Albert Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Spread Probability')
}

####SMALL-WORLD (WATTS-STROGATZ) FUNCTIONS####
plot_small_world_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,1,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- sample_smallworld(dim=1, size=N, nei=5, p=p)
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs)
  title("Boxplot of the  Eigenalues of the Small-World Graph")
  data
}

plot_small_world <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 1, step)) {
    G <- sample_smallworld(dim=1, size=N, nei=5, p=p)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p, x = seq(0, 1, 1 / (length(res) - 1)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes') +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Small-World Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Spread Probability')
}

####TREE FUNCTIONS####
plot_tree_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,1,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- make_tree(N, children = 2, mode = "undirected")
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs)
  title("Boxplot of the  Eigenalues of the Tree Graph")
  data
}



plot_tree <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 1, step)) {
    G <- make_tree(N, children = 2, mode = "undirected")
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p, x = seq(0, 1, 1 / (length(res) - 1)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'SIS Simulation Result') +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Tree Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Spread Probability')
}


####STAR FUNCTIONS####
plot_star_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,1,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- make_star(N, mode = "undirected", center = 1)
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs)
  title("Boxplot of the  Eigenalues of the Star Graph")
  data
}


plot_star <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 1, step)) {
    G <- make_star(N, mode = "undirected", center = 1)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p, x = seq(0, 1, 1 / (length(res) - 1)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes') +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Star Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Spread Probability')
}


####COMPLETE GRAPHS FUNCTIONS (VERY SIMILAR TO ERDOS-RENYI)####
plot_complete_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,1,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- make_full_graph(N, directed = FALSE, loops = FALSE)
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs)
  title("Boxplot of the  Eigenalues of the Complete Graph")
  data
}


plot_complete <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  for (p in seq(step, 1, step)) {
    G <- make_full_graph(N, directed = FALSE, loops = FALSE)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p, x = seq(0, 1, 1 / (length(res) - 1)), y = res))
  }
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes') +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Complete Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Spread Probability')
}


####LATTICE GRAPHS FUNCTIONS (Ring is lattice with dimension=1)####
#make_lattice(dimvector = NULL, length = NULL, dim = NULL, nei = 1, directed = FALSE, mutual = FALSE, circular = FALSE) fOR MORE THAN ONE DIMENSION USE THIS INSTEAD OF RING

plot_lattice_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,1,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- make_ring(n=N, directed = FALSE)
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs)
  title("Boxplot of the  Eigenalues of the Lattice Graph")
  data
}


plot_lattice <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  for (p in seq(step, 1, step)) {
    G <- make_ring(n=N, directed = FALSE,)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p, x = seq(0, 1, 1 / (length(res) - 1)), y = res))
  }
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes') +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Lattice Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Spread Probability')
}


##filepath to store plots
filepath<-'/home/didac/Documents/CSN/LAB/CSN-MIRI/assignment7/img'

N <-1000
step <- 0.2
num_rep <- 10
#gamma <- 0.3
listgamma <- c(0.15, 0.3, 0.45, 0.6, 0.75, 0.9)

#######################
##### EXERCISE 1 ######
#######################

#####GRAPHS SELECTION (SEARCH FOR DIFFERENT BEHAVIOURS IN RESULTING PLOTS)#####
#plot_erdos_renyi(N,step,num_rep,gamma)
#plot_barabasi_albert(N,step,num_rep,gamma)
#plot_small_world(N,step,num_rep,gamma)
#plot_tree(N,step,num_rep,gamma)
#plot_star(N,step,num_rep,gamma)
#plot_complete(N,step,num_rep,gamma)    Similar to erdos-renyi
#plot_lattice(N,step,num_rep,gamma)     Similar to tree


for (gamma in listgamma){
  # This experiment  allows us to see how SIS converges  depending  on  the
  # probabiblity of each graph
  ErdoRenyiConvergence <- plot_erdos_renyi(N,step,num_rep,gamma)
  filename<-paste0('ErdosRenyiConvergence_', gamma, '.png')
  ggsave(filename, path=filepath, device='png', plot = ErdoRenyiConvergence, width = 6, height = 4, units = 'in')
  BarabasiAlbertConvergence <- plot_barabasi_albert(N,step,num_rep,gamma)
  filename<-paste0('BarabasiAlbertConvergence_', gamma, '.png')
  ggsave(filename, path=filepath, device='png', plot = BarabasiAlbertConvergence, width = 6, height = 4, units = 'in')
  SmallWorldConvergence <- plot_small_world(N,step,num_rep,gamma)
  filename<-paste0('SmallWorldConvergence_', gamma, '.png')
  ggsave(filename, path=filepath, device='png', plot = SmallWorldConvergence, width = 6, height = 4, units = 'in')
  TreeConvergence <- plot_tree(N,step,num_rep,gamma)
  filename<-paste0('TreeConvergence_', gamma, '.png')
  ggsave(filename, path=filepath, device='png', plot = TreeConvergence, width = 6, height = 4, units = 'in')
  StarConvergence <- plot_star(N,step,num_rep,gamma)
  filename<-paste0('StarConvergence_', gamma, '.png')
  ggsave(filename, path=filepath, device='png', plot = StarConvergence, width = 6, height = 4, units = 'in')
  
  # This experiment  allows us to see how much variance there is in the eigenvalues
  # for each graph.
  #
  filename<-paste0(filepath, '/ErdosRenyiEigenvalues_', gamma, '.png')
  png(filename)
  ErdoRenyiEigenvalues <- plot_erdos_renyi_eigen(N, step, num_rep)
  dev.off()
  filename<-paste0(filepath, '/BarabasiAlbertEigenvalues_', gamma, '.png')
  png(filename)
  BarabasiAlbertEigenvalues <- plot_barabasi_albert_eigen(N, step, num_rep)
  dev.off()
  filename<-paste0(filepath, '/SmallWorldEigenvalues_', gamma, '.png')
  png(filename)
  SmallWorldEigenvalues <- plot_small_world_eigen(N, step, num_rep)
  dev.off()
  filename<-paste0(filepath, '/TreeEigenvalues_', gamma, '.png')
  png(filename)
  TreeEigenvalues <- plot_tree_eigen(N, step, num_rep)
  dev.off()
  filename<-paste0(filepath, '/StarEigenvalues_', gamma, '.png')
  png(filename)
  StarEigenvalues <- plot_star_eigen(N, step, num_rep)
  dev.off()
}

#######################
##### EXERCISE 2 ######
#######################

#....
