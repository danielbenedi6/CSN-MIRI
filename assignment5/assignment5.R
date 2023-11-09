library(igraph)
library(clustAnalytics)

jaccard_sim <- function(left, right) {
  M = matrix(nrow=length(left), ncol=length(right))
  
  for(cluster1 in unique(unname(membership(left))) ) {
    for(cluster2 in unique(unname(membership(right))) ) {
      intersection = length(intersect(left[cluster1], right[cluster2]))
      union = length(left[cluster1]) + length(right[cluster2]) / intersection
      M[cluster1,cluster2] <- intersection/union
    } 
  }
  
  M <- as.data.frame(M)
  colnames(M) <-  1:length(right)
  M
}

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


M <- as_adjacency_matrix(as.undirected(karate, mode="each"))

evaluate_significance(karate, 
                      alg_list=list(
                        Louvain=cluster_louvain,
                        "label prop"=cluster_label_prop,
                        walktrap=cluster_walktrap,
                        fastgreedy=fastgreedy.community
                      ),
                      gt_clustering = V(karate)$Faction
)

plot(karate, vertex.color=V(karate)$Faction)

B <- matrix(c(1, 0.2, 0.2, 1), ncol=2)
G <- barabasi_albert_blocks(
        m=4,
        p=c(0.5,0.5),
        B=B,
        t_max=100,
        type="Hajek",
        sample_with_replacement = FALSE
      )
plot(G, vertex.color=(V(G)$label),vertex.label=NA,vertex.size=10)

