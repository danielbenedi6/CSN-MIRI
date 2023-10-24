library(ggplot2)
library(stats)

read_dataset <- function(language) {
  dataset = read.table(paste("./data/",language,"_dependency_tree_metrics.txt",sep=""), header=FALSE)
  colnames(dataset) = c("vertices","degree_2nd_moment", "mean_length")
  dataset = dataset[order(dataset$vertices),]
  dataset
}

languages = c("Arabic", "Basque", "Catalan",
              "Chinese", "Czech", "English",
              "Greek", "Hungarian", "Italian",
              "Turkish")
for(language in languages){
  data = read_dataset(language)
  means = aggregate(data, list(data$vertices), mean)
  
  plt <- ggplot() + geom_point(data=dataset, aes(x=vertices,y=mean_length))
  plt <- plt + geom_line(data=means, aes(x=vertices,y=mean_length), color="green")
  plt <- plt + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
  
  # lines(log(mean_Catalan$vertices),log((mean_Catalan$vertices+1)/3), col = "red")
}