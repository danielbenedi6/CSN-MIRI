library(igraph)
library(clustAnalytics)
library(stringr)

load_network <- function(network, ground_truth) {
  G <- read_graph(network, format="edgelist", n=0, directed=FALSE)
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

evaluate <- function(graph, methods) {
  D <- data.frame(matrix(ncol = 0, nrow = 1))
  
  # Note: Barabasi-Albert does not store in Faction
  #       it should be adapted the graph before calling
  #       this method.
  if(is.null(V(graph)$Faction)) {
    # TODO: Change to proper selection method
    BASELINE <- unlist(methods)[[1]](graph)
    label2 <- D$method[1]
  } else {
    BASELINE <- make_clusters(graph, membership = setNames(V(graph)$Faction, names(V(graph))))
    label2 <- "Ground Truth"
  }
  
  for(method_name in names(methods)) {
    clustered <- methods[[method_name]](graph)
    global_jaccard <- Wmean(match_clusters(jaccard_sim(clustered, BASELINE), method_name, "baseline"), cluster_weights(clustered))
    D[method_name] = c(global_jaccard)
  }
  return(D)
}

E <- evaluate(
        karate,
        list(
          "Louvain" = cluster_louvain,
          "Label Propagation" = cluster_label_prop,
          "Walktrap" = cluster_walktrap,
          "Edge Betweenness" = cluster_edge_betweenness,
          "Fast Greedy" = cluster_fast_greedy,
          "Spin-Glass" = cluster_spinglass
        )
)

# Wmean(match_clusters(jaccard_sim(fc,wc), "FC", "WC"),cluster_weights(fc) )
# 
# data(karate, package="igraphdata")
# karate <- upgrade_graph(karate)
# 
# 
# wc <- walktrap.community(karate)
# modularity(wc)
# unname(membership(wc))
# plot(wc,karate)
# 
# 
# fc <- fastgreedy.community(karate)
# modularity(fc)
# dendPlot(fc)
# plot(fc, karate)


# 
# M <- as_adjacency_matrix(as.undirected(karate, mode="each"))
# 
# evaluate_significance(karate, 
#                     alg_list=list(
#                        Louvain=cluster_louvain,
#                        "label prop"=cluster_label_prop,
#                        walktrap=cluster_walktrap,
#                        fastgreedy=fastgreedy.community
#                      ),
#                      gt_clustering = V(karate)$Faction
# )

# plot(karate, vertex.color=V(karate)$Faction)
# 
# B <- matrix(c(1, 0.2, 0.2, 1), ncol=2)
# G <- barabasi_albert_blocks(
#         m=4,
#         p=c(0.5,0.5),
#         B=B,
#         t_max=100,
#         type="Hajek",
#         sample_with_replacement = FALSE
#       )
# plot(G, vertex.color=(V(G)$label),vertex.label=NA,vertex.size=10)

