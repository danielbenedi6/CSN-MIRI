# Lab1: Reproduce the graphs that were introduced in the first lecture.
library(igraph)

#(a) plotting the clustering coefficient and the average shortest - path as a
# function of the parameter p of the WS model

p <- 10^(seq(-4,0,0.2))

N = 100
NODES = 1000
for (neigh in c(2,3,4,5,6,8,10,20) ){
  L <- rep(0,length(p))
  C <- rep(0,length(p))

  for (i in 1:length(p)) {
    for (j in 1:N) {
      # Note that the smaller the population and the larger the neighborhood,
      # then the clustering coefficient will be higher
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
  plot(p, L, ylim = c(0,1), ylab='coeff', xlab='p', pch=0, log='x', ,
       main= paste("Transitivity and Mean distance for WS with ",NODES,
                   " nodes and neigh = ",neigh,sep=""))
  points(p, C, ylim = c(0,1), ylab='coeff', xlab='p', pch=16)
  legend("bottomleft", inset=.02,c('C(p)/C(0)','L(p)/L(0)'), pch=c(0,16))
  dev.off()
}

#b) Plotting the average shortest-path length as a function of the network size 
# of the ER model
average_path_length_ER <- function(n, p) {
  g <- erdos.renyi.game(n, p)
  # The example code for the first lab lead to errors for me when trying to plot
  # the graph, therefore I am using another code version for the ER model.
  # g <- sample_gnm(n,p)
  # If graph is not connected, calculate for the largest connected component
  if (!is.connected(g)) {
    g <- induced.subgraph(g, which.max(clusters(g)$csize))
  } # TODO: Write this down in the report
  return(mean_distance(g, directed = FALSE))
}

# Generating a sequence of steps that is similar to the one shown in the lecture
max_nodes <- 100000
# The amount of steps shown in the according task graph is 21
n_steps <- 21 
power <- 2
linear_seq <- seq(0, 1, length.out = n_steps)
non_linear_seq <- linear_seq^power
network_sizes <- round(1 + (max_nodes - 1) * non_linear_seq)

# fixed probability
#p <- 0.1  
# Static p
#avg_lengths <- sapply(network_sizes, function(n) average_path_length_ER(n, p))
#Dynamic p
avg_lengths <- sapply(network_sizes, function(n) average_path_length_ER(n, (2*log(n)/n) ))

plot(network_sizes, avg_lengths, type="b", xlab="Number of Nodes", ylab="Average Shortest Path Length",
     main=paste("ER model with p = 2*log(n)/n"))
