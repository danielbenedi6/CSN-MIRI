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
    data <- rbind(data, data.frame(prob = p,x = seq(1/length(res), 1, 1 /length(res)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Erdös-Renyi Graph varying graph density') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Expected Density')
}



####BARABASI-ALBERT FUNCTIONS######
BarabasiAlbert <- function(N,num_clusters) {
  B <- diag(0.8, num_clusters) + matrix(0.2,num_clusters,num_clusters)
  barabasi_albert_blocks(m=8, p=rep(1/num_clusters,num_clusters), B=B, t_max=N)
}


plot_barabasi_albert_eigen <- function(N,step,num_rep)  {
  probs  <- rep(seq(step,5*step,step), num_rep)
  eigenvals <- rep(0,length(probs))
  c <- 0
  for(p in probs) {
    c <- c + 1
    G <- BarabasiAlbert(N, p)
    eigenvals[c] <-max(eigen(as_adj(G))$values)
  }
  data <- data.frame(probs, eigenvals)
  boxplot(eigenvals ~ probs, xlab="Number of clusters")
  title("Boxplot of the  Eigenalues of the Barabasi-Albert Graph")
  data
}

plot_barabasi_albert <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 5*step, step)) {
    G <- BarabasiAlbert(N, p)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p,x = seq(1/length(res), 1, 1 /length(res)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Barabasi-Albert Graph varying number of clusters') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Number of clusters')
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
    data <- rbind(data, data.frame(prob = p,x = seq(1/length(res), 1, 1 /length(res)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Small-World Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Rewiring Probability')
}

####TREE FUNCTIONS####
## No sense to do this experiment because tree is deterministic
# plot_tree_eigen <- function(N,step,num_rep)  {
#   probs  <- rep(seq(step,1,step), num_rep)
#   eigenvals <- rep(0,length(probs))
#   c <- 0
#   for(p in probs) {
#     c <- c + 1
#     G <- make_tree(N, children = 2, mode = "undirected")
#     eigenvals[c] <-max(eigen(as_adj(G))$values)
#   }
#   
#   data <- data.frame(probs, eigenvals)
#   boxplot(eigenvals ~ probs)
#   title("Boxplot of the  Eigenalues of the Tree Graph")
#   data
# }

plot_tree <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 5*step, step)) {
    G <- make_tree(N, children = p, mode = "undirected")
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p,x = seq(1/length(res), 1, 1 /length(res)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'SIS Simulation Result', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Tree Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Number of children')
}


####STAR FUNCTIONS####
## No sense to do this experiment because star is deterministic
# plot_star_eigen <- function(N,step,num_rep)  {
#   probs  <- rep(seq(step,1,step), num_rep)
#   eigenvals <- rep(0,length(probs))
#   c <- 0
#   for(p in probs) {
#     c <- c + 1
#     G <- make_star(N, mode = "undirected", center = 1)
#     eigenvals[c] <-max(eigen(as_adj(G))$values)
#   }
#   
#   data <- data.frame(probs, eigenvals)
#   boxplot(eigenvals ~ probs)
#   title("Boxplot of the  Eigenalues of the Star Graph")
#   data
# }


plot_star <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  G <- make_star(N, mode = "undirected", center = 1)
  max_eigenvalue <- max(eigen(as_adj(G))$values)
  beta <- gamma / max_eigenvalue + epsilon
  res <- repetition(num_rep, G, gamma, beta, 100)
  data <- rbind(data, data.frame(x = seq(1/length(res), 1, 1/length(res) ), y = res))
  
  ggplot(data, aes(x = x, y = y)) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Star Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line()
}


####COMPLETE GRAPHS FUNCTIONS (VERY SIMILAR TO ERDOS-RENYI)####
## No sense to do this experiment because complete graphs are deterministic
# plot_complete_eigen <- function(N,step,num_rep)  {
#   probs  <- rep(seq(step,1,step), num_rep)
#   eigenvals <- rep(0,length(probs))
#   c <- 0
#   for(p in probs) {
#     c <- c + 1
#     G <- make_full_graph(N, directed = FALSE, loops = FALSE)
#     eigenvals[c] <-max(eigen(as_adj(G))$values)
#   }
#   data <- data.frame(probs, eigenvals)
#   boxplot(eigenvals ~ probs)
#   title("Boxplot of the  Eigenalues of the Complete Graph")
#   data
# }


plot_complete <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  G <- make_full_graph(N, directed = FALSE, loops = FALSE)
  max_eigenvalue <- max(eigen(as_adj(G))$values)
  beta <- gamma / max_eigenvalue + epsilon
  res <- repetition(num_rep, G, gamma, beta, 100)
  data <- rbind(data, data.frame(x = seq(1/length(res), 1, 1/length(res) ), y = res))
    
  ggplot(data, aes(x = x, y = y)) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Complete Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line()
}


####LATTICE GRAPHS FUNCTIONS (Ring is lattice with dimension=1)####
## No sense to do this experiment because lattice is deterministic
# plot_lattice_eigen <- function(N,step,num_rep)  {
#   probs  <- rep(seq(step,1,step), num_rep)
#   eigenvals <- rep(0,length(probs))
#   c <- 0
#   for(p in probs) {
#     c <- c + 1
#     G <- make_ring(n=N, directed = FALSE)
#     eigenvals[c] <-max(eigen(as_adj(G))$values)
#   }
#   data <- data.frame(probs, eigenvals)
#   boxplot(eigenvals ~ probs)
#   title("Boxplot of the  Eigenalues of the Lattice Graph")
#   data
# }


plot_lattice <- function(N, step, num_rep, gamma, epsilon) {
  data <- data.frame()
  
  for (p in seq(step, 5*step, step)) {
    G <- make_lattice(length=ceiling(N**(1/step)) , dim=step, directed = FALSE,)
    max_eigenvalue <- max(eigen(as_adj(G))$values)
    beta <- gamma / max_eigenvalue + epsilon
    res <- repetition(num_rep, G, gamma, beta, 100)
    data <- rbind(data, data.frame(prob = p,x = seq(1/length(res), 1, 1 /length(res)), y = res))
  }
  
  ggplot(data, aes(x = x, y = y, color = factor(prob))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation on the Lattice Graph varying probability') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Number of dimensions')
}


##filepath to store plots
filepath<-'./img'

N <-1000
step <- 0.2
num_rep <- 10
gamma <- 0.3
listgamma <- c(0.15, 0.3, 0.45, 0.6, 0.75, 0.9)
listp0 <- c(0.1, 0.25, 0.5, 0.75, 0.9)

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

for (p0 in listp0) { # We will experiment varying the percentage of initial inffected
  print(sprintf("Simulating with %.1f%% infected initialy...", p0*100))
  infected <- floor(p0*N)
  gamma <- 0.3 # Recovery probability
  beta <- 0.4  # Spread probability
  
  data <- data.frame()
  eigens <- data.frame()
  
  G <- erdos.renyi.game(N, 0.4)
  max_eigen <- max(eigen(as_adj(G))$values)
  if(beta > gamma/max_eigen) {
    print("    Erdös-Renyi is above threshold")
  } else {
    print("    Erdös-Renyi is below threshold")
  }
  res_er <- repetition(num_rep, G, gamma, beta, infected)
  eigens <- rbind(eigens, data.frame(graph="Erdös-Renyi", max_eigenvalue=c(max_eigen), threshold=c(gamma/max_eigen)))
  
  num_clusters <- 4
  B <- diag(0.9,num_clusters) + matrix(0.1,num_clusters,num_clusters)
  G <- barabasi_albert_blocks(m=8, p=rep(1/num_clusters, num_clusters), B=B, t_max=N)
  max_eigen <- max(eigen(as_adj(G))$values)
  if(beta > gamma/max_eigen) {
    print("    Barbási-Albert is above threshold")
  } else {
    print("    Barabási-Albert is below threshold")
  }
  res_ba <- repetition(num_rep, G, gamma, beta, infected)
  eigens <- rbind(eigens, data.frame(graph="Barabási-Albert", max_eigenvalue=c(max_eigen), threshold=c(gamma/max_eigen)))
  
  G <- sample_smallworld(dim=1, size=N, nei=5, p=0.25)
  max_eigen <- max(eigen(as_adj(G))$values)
  if(beta > gamma/max_eigen) {
    print("    Watts-Strogatz is above threshold")
  } else {
    print("    Watts-Strogatz is below threshold")
  }
  res_ws <- repetition(num_rep, G, gamma, beta, infected)
  eigens <- rbind(eigens, data.frame(graph="Watts-Strogatz", max_eigenvalue=c(max_eigen), threshold=c(gamma/max_eigen)))
  
  G <- make_tree(N, children = 4, mode = "undirected")
  max_eigen <- max(eigen(as_adj(G))$values)
  if(beta > gamma/max_eigen) {
    print("    Tree is above threshold")
  } else {
    print("    Tree is below threshold")
  }
  res_tree <- repetition(num_rep, G, gamma, beta, infected)
  eigens <- rbind(eigens, data.frame(graph="Tree", max_eigenvalue=c(max_eigen), threshold=c(gamma/max_eigen)))
  
  G <- make_star(N, mode = "undirected", center = 1)
  max_eigen <- max(eigen(as_adj(G))$values)
  if(beta > gamma/max_eigen) {
    print("    Star is above threshold")
  } else {
    print("    Star is below threshold")
  }
  res_star <- repetition(num_rep, G, gamma, beta, infected)
  eigens <- rbind(eigens, data.frame(graph="Star", max_eigenvalue=c(max_eigen), threshold=c(gamma/max_eigen)))
  
  print("    Saving results...")
  
  print(xtable(eigens, type="latex", auto=TRUE),  include.rownames=FALSE, file=sprintf("./tables/eigens_%.2f.tex",p0))
  
  max_len <- max(length(res_er), length(res_ba), length(res_ws), length(res_tree), length(res_star))
  res_er <- append(res_er, rep(res_er[length(res_er)], max_len-length(res_er)))
  res_ba <- append(res_ba, rep(res_ba[length(res_ba)], max_len-length(res_ba)))
  res_ws <- append(res_ws, rep(res_ws[length(res_ws)], max_len-length(res_ws)))
  res_tree <- append(res_tree, rep(res_tree[length(res_tree)], max_len-length(res_tree)))
  res_star <- append(res_star, rep(res_star[length(res_star)], max_len-length(res_star)))
  
  data <- rbind(data, data.frame(graph="Erdös-Renyi", x = seq(1 / max_len, 1, 1 / max_len), y = res_er))
  data <- rbind(data, data.frame(graph="Barabási-Albert", x = seq(1 / max_len, 1, 1 / max_len), y = res_ba))
  data <- rbind(data, data.frame(graph="Watts-Strogatz", x = seq(1 / max_len, 1, 1 / max_len), y = res_ws))
  data <- rbind(data, data.frame(graph="Tree", x = seq(1 / max_len, 1, 1 / max_len), y = res_tree))
  data <- rbind(data, data.frame(graph="Star", x = seq(1 / max_len, 1, 1 / max_len), y = res_star))
  
  ggplot(data, aes(x = x, y = y, color = factor(graph))) +
    scale_x_log10(name = 'Normalized Time') +
    scale_y_continuous(name = 'Proportion of infected nodes', limits=c(0,1)) +
    #theme_minimal() +
    labs(title = 'SIS Simulation') +
    theme(legend.position = 'bottom') + geom_line() +
    scale_color_discrete(name = 'Network type')
  ggsave(sprintf("./img/AllNetworks_%.2f.pdf",p0))
}

# This experiment  allows us to see how much variance there is in the eigenvalues
# for each graph.
#
filename<-paste0(filepath, '/ErdosRenyiEigenvalues.pdf')
pdf(filename)
ErdoRenyiEigenvalues <- plot_erdos_renyi_eigen(N, step, num_rep)
dev.off()
filename<-paste0(filepath, '/BarabasiAlbertEigenvalues.pdf')
pdf(filename)
BarabasiAlbertEigenvalues <- plot_barabasi_albert_eigen(N, 2, num_rep)
dev.off()
filename<-paste0(filepath, '/SmallWorldEigenvalues.pdf')
pdf(filename)
SmallWorldEigenvalues <- plot_small_world_eigen(N, step, num_rep)
dev.off()
# filename<-paste0(filepath, '/TreeEigenvalues.pdf')
# pdf(filename)
# TreeEigenvalues <- plot_tree_eigen(N, step, num_rep)
# dev.off()
# filename<-paste0(filepath, '/StarEigenvalues.pdf')
# pdf(filename)
# StarEigenvalues <- plot_star_eigen(N, step, num_rep)
# dev.off()

#######################
##### EXERCISE 2 ######
#######################

for (gamma in listgamma){
  print(sprintf("Simulating with gamma=%.2f ...", gamma))
  
  # This experiment  allows us to see how SIS converges  depending  on  the
  # probabiblity of each graph
  print("    Erdös-Renyi Pandemic")
  ErdoRenyiConvergence <- plot_erdos_renyi(N,step,num_rep,gamma,0.05)
  filename<-paste0('ErdosRenyiConvergence_Above_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = ErdoRenyiConvergence, width = 6, height = 4, units = 'in')
  print("    Barabási-Albert Pandemic")
  BarabasiAlbertConvergence <- plot_barabasi_albert(N,2,num_rep,gamma,0.05)
  filename<-paste0('BarabasiAlbertConvergence_Above_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = BarabasiAlbertConvergence, width = 6, height = 4, units = 'in')
  print("    Watts-Strogatz Pandemic")
  SmallWorldConvergence <- plot_small_world(N,step,num_rep,gamma,0.05)
  filename<-paste0('SmallWorldConvergence_Above_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = SmallWorldConvergence, width = 6, height = 4, units = 'in')
  print("    Tree Pandemic")
  TreeConvergence <- plot_tree(N,2,num_rep,gamma,0.05)
  filename<-paste0('TreeConvergence_Above_', gamma, '.pdf')
  print("    Star Pandemic")
  ggsave(filename, path=filepath, device='pdf', plot = TreeConvergence, width = 6, height = 4, units = 'in')
  StarConvergence <- plot_star(N,step,num_rep,gamma,0.05)
  filename<-paste0('StarConvergence_Above_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = StarConvergence, width = 6, height = 4, units = 'in')
  
  
  # This experiment  allows us to see how SIS converges  depending  on  the
  # probabiblity of each graph
  print("    Erdös-Renyi Controlled")
  ErdoRenyiConvergence <- plot_erdos_renyi(N,step,num_rep,gamma,-0.05)
  filename<-paste0('ErdosRenyiConvergence_Below_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = ErdoRenyiConvergence, width = 6, height = 4, units = 'in')
  print("    Barabási-Albert Controlled")
  BarabasiAlbertConvergence <- plot_barabasi_albert(N,2,num_rep,gamma,-0.05)
  filename<-paste0('BarabasiAlbertConvergence_Below_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = BarabasiAlbertConvergence, width = 6, height = 4, units = 'in')
  print("    Watts-Strogatz Controlled")
  SmallWorldConvergence <- plot_small_world(N,step,num_rep,gamma,-0.05)
  filename<-paste0('SmallWorldConvergence_Below_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = SmallWorldConvergence, width = 6, height = 4, units = 'in')
  print("    Tree Controlled")
  TreeConvergence <- plot_tree(N,2,num_rep,gamma,-0.05)
  filename<-paste0('TreeConvergence_Below_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = TreeConvergence, width = 6, height = 4, units = 'in')
  print("    Star Controlled")
  StarConvergence <- plot_star(N,step,num_rep,gamma,-0.05)
  filename<-paste0('StarConvergence_Below_', gamma, '.pdf')
  ggsave(filename, path=filepath, device='pdf', plot = StarConvergence, width = 6, height = 4, units = 'in')
}
