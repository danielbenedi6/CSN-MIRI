library(igraph)
library(clustAnalytics)

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
  
  D <- data.frame(
    method = names(methods),
    result = sapply(methods, function(f) f(graph), simplify = FALSE)
  )
   
  # louvain <- cluster_louvain(graph)
  # label_prop <- cluster_label_prop(graph)
  # walktrap <- cluster_walktrap(graph)
  # edge_betweenness <- cluster_edge_betweenness(graph)
  # fast_greedy <- cluster_fast_greedy(graph)
  # spinglass <- cluster_spinglass(graph)
  
  if(is.nul(V(G)$Faction)) {
    # TODO: Change to proper selection method
    BASELINE <- D$result[1]
    label2 <- D$method[1]
  } else {
    BASELINE <- V(G)$Faction
    label2 <- "Ground Truth"
  }
  
  for(method1 in methods) {
     # Compare method1 with BASELINE
  }
  
  
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

Wmean(match_clusters(jaccard_sim(fc,wc), "FC", "WC"),cluster_weights(fc) )

data(karate, package="igraphdata")
karate <- upgrade_graph(karate)


wc <- walktrap.community(karate)
modularity(wc)
unname(membership(wc))
plot(wc,karate)


fc <- fastgreedy.community(karate)
modularity(fc)
dendPlot(fc)
plot(fc, karate)


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

