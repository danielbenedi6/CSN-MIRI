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

evaluate <- function(graph, methods, name, baseline_method = "Louvain") {
  D <- data.frame(matrix(ncol = 0, nrow = 1))
  
  # Note: Barabasi-Albert does not store in Faction
  #       it should be adapted the graph before calling
  #       this method.
  if(is.null(V(graph)$Faction)) {
    # TODO: Change to proper selection method
    BASELINE <- methods[[baseline_method]](graph)
    label2 <- D$method[1]
  } else {
    V(graph)$label <- V(graph)$Faction
    BASELINE <- make_clusters(graph, membership = setNames(V(graph)$Faction, names(V(graph))))
    label2 <- "Ground Truth"
  }
  
  pdf(paste("./plot/",name,"_baseline.pdf", sep=""))
  plot(BASELINE, graph)
  dev.off()
  
  for(method_name in names(methods)) {
    print(paste("computing ",method_name))
    
    clustered <- methods[[method_name]](graph)
    global_jaccard <- Wmean(match_clusters(jaccard_sim(clustered, BASELINE), method_name, label2), cluster_weights(clustered))
    D[method_name] = c(global_jaccard)
    print(paste("computing pdf for ",method_name))
    
    pdf(paste("./plot/",name,"_",method_name,".pdf", sep=""))
    plot(clustered, graph)
    dev.off()
  }
  return(D)
}

set.seed(0)

summary <- data.frame(
              "Network"=character(),
              "Vertices"=integer(),
              "Edges"=integer(),
              "Mean Degree"=double(),
              "Density"=double()
          )
clustering_algorithms = list(
  "Louvain" = cluster_louvain,
  "Label Propagation" = cluster_label_prop,
  "Walktrap" = cluster_walktrap,
  "Edge Betweenness" = cluster_edge_betweenness,
  "Fast Greedy" = cluster_fast_greedy,
  "Spin-Glass" = cluster_spinglass
)

# Case 1
print("Calculating karate...")
data(karate, package="igraphdata")
karate <- upgrade_graph(karate)

summary[nrow(summary) + 1,] = c(
    "Karate",
    length(V(karate)),
    length(E(karate)),
    round(mean(degree(karate)),digits=4),
    round(edge_density(karate),digits=4)
)

E_karate <- evaluate(
  karate,
  clustering_algorithms,
  "karate"
)

# Case 2
print("Calculating barabasi albert...")

num_vertices = 200
added_edges = 4 # num_vertices = added_edges * num_vertices
num_clusters = 4

# All clusters have the same probability, we could change this
cluster_probabilities = rep(1/num_clusters, num_clusters)

# Affinity blocks
B <- matrix(c(1, 0.1, 0.1, 1,
              1, 0.1, 0.1, 1,
              1, 0.1, 0.1, 1,
              1, 0.1, 0.1, 1), ncol=num_clusters)

BA <- barabasi_albert_blocks(
       m=added_edges,
       p=cluster_probabilities,
       B=B,
       t_max=num_vertices,
       type="Hajek",
       sample_with_replacement = FALSE
     )
V(BA)$Faction <- V(BA)$label

summary[nrow(summary) + 1,] = c(
  "Barabasi-Albert",
  length(V(BA)),
  length(E(BA)),
  round(mean(degree(BA)),digits=4),
  round(edge_density(BA),digits=4)
)

E_BA <- evaluate(
  BA,
  clustering_algorithms,
  "BA"
)

# Case 3
print("Calculating enron...")

data(enron, package="igraphdata")
enron <- upgrade_graph(enron)
enron <- as.undirected(enron, mode="collapse")
enron <- simplify(enron)

# ENRON is not connected
# There are only two nodes disconnected.
# We opted to simply remove them
enron <- ensure_connectedness(enron)

summary[nrow(summary) + 1,] = c(
  "ENRON",
  length(V(enron)),
  length(E(enron)),
  round(mean(degree(enron)),digits=4),
  round(edge_density(enron),digits=4)
)
E_enron <- evaluate(
  enron,
  clustering_algorithms,
  "enron",
  "Label Propagation"
)

# Case 4
print("Calculating DBLP...")

#DBLP <- load_network("./data/com-dblp.ungraph.txt","./data/com-dblp.top5000.cmty.txt")
DBLP <- load_network("./data/com-dblp.ungraph.txt")
print("ensure_connectedness...")
# DBLP is not connected
# All the disconnected components have size 1
# We opted to simply remove them
DBLP <- ensure_connectedness(DBLP)
print("induced subgraph...")

# Random subgraph to easy execution
DBLP = induced_subgraph(
                  DBLP,
                  random_walk(
                     DBLP,
                     start=1,
                     steps=floor(0.01*vcount(DBLP))
                   )
                 )

#DBLP <- induced_subgraph(DBLP, sample(1:vcount(DBLP), floor(0.1*vcount(DBLP))))
print("ensure_connectedness again...")

DBLP <- ensure_connectedness(DBLP)
print("compute summary...")

summary[nrow(summary) + 1,] = c(
  "DBLP",
  length(V(DBLP)),
  length(E(DBLP)),
  round(mean(degree(DBLP)),digits=4),
  round(edge_density(DBLP),digits=4)
)
print("evaluating cluster algorithms")
E_DBLP<- evaluate(
  DBLP,
  clustering_algorithms,
  "dblp"
)

# Table

evals <- rbind(E_karate,E_BA, E_enron, E_DBLP)
rownames(evals) <- c("karate", "Barabasi-Albert", "ENRON", "DBLP")

print(xtable(summary, type="latex", auto=TRUE),  include.rownames=FALSE, file="./table/summary.tex")
print(xtable(evals, type="latex", auto=TRUE), file="./table/evaluations.tex")

karate_significance = evaluate_significance(karate, clustering_algorithms, gt_clustering = V(karate)$Faction)
print(xtable(karate_significance, type="latex", auto=TRUE), file="./table/karate_significance.tex")

ba_significance = evaluate_significance(BA, clustering_algorithms, gt_clustering = V(BA)$Faction)
print(xtable(ba_significance, type="latex", auto=TRUE), file="./table/barabasi_albert_significance.tex")

enron_significance = evaluate_significance(enron, clustering_algorithms)
print(xtable(enron_significance, type="latex", auto=TRUE), file="./table/enron_significance.tex")

DBLP_significance = evaluate_significance(DBLP, clustering_algorithms)
print(xtable(DBLP_significance, type="latex", auto=TRUE), file="./table/DBLP_significance.tex")

