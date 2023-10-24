library(ggplot2)
library(stats)

read_dataset <- function(language) {
  dataset = read.table(paste("./data/",language,"_dependency_tree_metrics.txt",sep=""), header=FALSE)
  colnames(dataset) = c("vertices","degree_2nd_moment", "mean_length")
  dataset = dataset[order(dataset$vertices),]
  dataset
}

generate_table_entry <- function(language, dataset) {
  sprintf(r"(%s & %d & %.3f & %.3f & %.3f & %.3f \\)",
          language,
          length(dataset),
          mean(dataset$vertices),
          sd(dataset$vertices),
          mean(dataset$mean_length),
          sd(dataset$mean_length)
          )
}

model1 <- function(n,b) {
  (n/2)**b
}
log_model1 <- function(n,b) {
  b * log(n/2)
}
model2 <- function(n,a,b) {
  a*n**b
}
log_model2 <- function(n,a,b) {
  log(a) + b*log(n)
}
model3 <- function(n,a,c) {
  a*exp(c*n)
}
log_model3 <- function(n,a,c) {
  log(a) + c*n
}
model4 <- function(n,a) {
  a*log(n)
}
log_model4 <- function(n,a) {
  log(a) + log(log(n))
}
model1_plus <- function(n.b,d) {
  (n/2)**b +d
}
log_model1_plus <- function(n,b,d) {
  log((n/2)**b + d)
}
model2_plus <- function(n,a,b,d) {
  a*n**b+d
}
log_model2_plus <- function(n,a,b,d) {
  log(a*n**b+d)
}
model3 <- function(n,a,c,d) {
  a*exp(c*n)+d
}
log_model3 <- function(n,a,c,d) {
  log(a*exp(c*n)+d)
}
model4 <- function(n,a) {
  a*log(n) + d
}
log_model4 <- function(n,a) {
  log(a*log(n) + d)
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
  plt <- plt + labs(title=paste("Models for mean length dependency for ", language))
  plt <- plt + labs(colour="Model")
  plt <- plt + xlab("Number of vertices")
  plt <- plt + ylab("Mean length")
  
  # lines(log(mean_Catalan$vertices),log((mean_Catalan$vertices+1)/3), col = "red")
}