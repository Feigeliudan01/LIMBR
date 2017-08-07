library(igraph)

demo(package="igraph")

nodes <- read.csv("pep_cor_nodes.csv", header=T, as.is=T)

edges <- read.csv("pep_cor_edges.csv", header=T, row.names=1)
edges <- 4*abs(as.matrix(edges))



g <- graph_from_incidence_matrix(edges,weighted=T)

V(g)$x <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
V(g)$y <- c(5, 4, 3, 2, 1, 5, 4, 3, 2, 1)
V(g)$color <- c('lightblue', 'lightblue', 'lightblue', 'lightblue', 'lightblue', 'red', 'red', 'red', 'red', 'red')

E(g)$width <- E(g)$weight

plot(g,vertex.size=30)

