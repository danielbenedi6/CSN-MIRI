library(igraph)
library(clustAnalytics)
library(stringr)
library(xtable)

load_network <- function(network, ground_truth) {
  G <- read_graph(network, format="edgelist", n=0, directed=FALSE)
  G <- simplify(G)
  
  if(!missing(ground_truth)) {
    # Read ground truth file as list of lists
    GT <- sapply(str_split(readLines(ground_truth),"\t"),as.numeric)
    # Note:
    # We consider that communities are no overlapping
    G$Faction <- rep(-1, vcount(G))
    
    # This is too slow, anything better?
    for(cluster in 1:length(GT)) {
      for(vertex in GT[[cluster]]) {
        G$Faction[vertex] <- cluster
      }
    }
  }
  
  return(G)
}

ensure_connectedness <- function(network) {
  if(!is_connected(network)) {
    comps <- components(network)
    lcc <- which.max(comps$csize)
    network <- induced_subgraph(network, V(network)[comps$membership == lcc])
  }
  
  network
}

jaccard_sim <- function(left, right) {
  M = matrix(nrow=length(left), ncol=length(right))

  for(cluster1 in 1:length(left) ) {
    for(cluster2 in  1:length(right)  ) {
      c1 = unlist(groups(left)[cluster1])
      c2 = unlist(groups(right)[cluster2])
      i = length(intersect(c1,c2))
      u = length(union(c1,c2))
      M[cluster1,cluster2] <- i/u
    } 
  }
  
  M <- as.data.frame(M)
  colnames(M) <-  1:length(right)
  M
}

best_indices <- function(JC) {
  return(apply(JC,1,which.max))
}

cluster_weights <- function(clustering) {
  sizes(clustering) / sum(sizes(clustering))
}

match_clusters <- function(JC, name1, name2) {
  best_values <- apply(JC,1,max)
  df <- as.data.frame(rbind(best_values), stringsAsFactors = FALSE)
  best_id <- best_indices(JC)
  names <- sapply(1:nrow(JC), function(i) {paste(name1,".", i, ", ", name2,".", best_id[i], sep="")})
  colnames(df) <- names
  return(df)
}

Wmean <- function(MC, weights) {
  weighted.mean(MC[1,], weights)
}

evaluate <- function(graph, methods, name) {
  D <- data.frame(matrix(ncol = 0, nrow = 1))
  
  # Note: Barabasi-Albert does not store in Faction
  #       it should be adapted the graph before calling
  #       this method.
  if(is.null(V(graph)$Faction)) {
    # TODO: Change to proper selection method
    BASELINE <- unlist(methods)[[1]](graph)
    label2 <- D$method[1]
  } else {
    V(graph)$label <- V(graph)$Faction
    BASELINE <- make_clusters(graph, membership = V(graph)$Faction)
    label2 <- "Ground Truth"
  }
  
  pdf(paste("./plot/",name,"_baseline.pdf", sep=""))
  plot(BASELINE, graph)
  dev.off()
  
  for(method_name in names(methods)) {
    clustered <- methods[[method_name]](graph)
    global_jaccard <- Wmean(match_clusters(jaccard_sim(clustered, BASELINE), method_name, label2), cluster_weights(clustered))
    D[method_name] = c(global_jaccard)
    
    pdf(paste("./plot/",name,"_",method_name,".pdf", sep=""))
    plot(clustered, graph)
    dev.off()
  }
  return(D)
}

# Case 1

data(karate, package="igraphdata")
karate <- upgrade_graph(karate)

E_karate <- evaluate(
  karate,
  list(
    "Louvain" = cluster_louvain,
    "Label Propagation" = cluster_label_prop,
    "Walktrap" = cluster_walktrap,
    "Edge Betweenness" = cluster_edge_betweenness,
    "Fast Greedy" = cluster_fast_greedy,
    "Spin-Glass" = cluster_spinglass
  ),
  "karate"
)

# Case 2

num_vertices = 100
added_edges = 4 # num_vertices = added_edges * num_vertices
num_clusters = 2

# All clusters have the same probability, we could change this
cluster_probabilities = rep(1/num_clusters, num_clusters)

# Affinity blocks
B <- matrix(c(1, 0.1, 0.1, 1), ncol=num_clusters)

BA <- barabasi_albert_blocks(
       m=added_edges,
       p=cluster_probabilities,
       B=B,
       t_max=num_vertices,
       type="Hajek",
       sample_with_replacement = FALSE
     )
V(BA)$Faction <- V(BA)$label

E_BA <- evaluate(
  BA,
  list(
    "Louvain" = cluster_louvain,
    "Label Propagation" = cluster_label_prop,
    "Walktrap" = cluster_walktrap,
    "Edge Betweenness" = cluster_edge_betweenness,
    "Fast Greedy" = cluster_fast_greedy,
    "Spin-Glass" = cluster_spinglass
  ),
  "BA"
)

# Case 3

data(enron, package="igraphdata")
enron <- upgrade_graph(enron)
enron <- as.undirected(enron, mode="collapse")
enron <- simplify(enron)

# ENRON is not connected
# There are only two nodes disconnected.
# We opted to simply remove them
enron <- ensure_connectedness(enron)

E_enron <- evaluate(
  enron,
  list(
    "Louvain" = cluster_louvain,
    "Label Propagation" = cluster_label_prop,
    "Walktrap" = cluster_walktrap,
    "Edge Betweenness" = cluster_edge_betweenness,
    "Fast Greedy" = cluster_fast_greedy,
    "Spin-Glass" = cluster_spinglass
  ),
  "enron"
)

# Case 4

#DBLP <- load_network("./data/com-dblp.ungraph.txt","./data/com-dblp.top5000.cmty.txt")
DBLP <- load_network("./data/com-dblp.ungraph.txt")

# DBLP is not connected
# All the disconnected components have size 1
# We opted to simply remove them
DBLP <- ensure_connectedness(DBLP)

# Random subgraph to easy execution
DBLP <- induced_subgraph(DBLP, sample(1:vcount(DBLP), floor(0.1*vcount(DBLP))))
DBLP <- ensure_connectedness(DBLP)

E_DBLP<- evaluate(
  DBLP,
  list(
    "Louvain" = cluster_louvain,
    "Label Propagation" = cluster_label_prop,
    "Walktrap" = cluster_walktrap,
    "Edge Betweenness" = cluster_edge_betweenness,
    "Fast Greedy" = cluster_fast_greedy,
    "Spin-Glass" = cluster_spinglass
  ),
  "dblp"
)

# Table

evals <- rbind(E_karate,E_BA, E_enron, E_DBLP)
rownames(evals) <- c("karate", "Barabasi-Albert", "ENRON", "DBLP")

print(xtable(evals, type="latex", auto=TRUE), file="./table/evaluations.tex")
