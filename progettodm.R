library(igraph)
library(MASS)
library(Matrix)

#Lettura dataset

nodes<- read.delim("nodes.txt",
                    header = TRUE,
                    sep = "\t",
                    colClasses = "character")

expression <- read.delim("expression_data.txt",
                          header = TRUE,
                          sep = "\t",
                          colClasses = "character")

clinical <- read.delim("clinical_data.txt",
                       header = TRUE,
                       sep = "\t",
                       colClasses = "character")

edges<- read.delim("edges.txt",
                    header = TRUE,
                    sep = "\t",
                    colClasses = "character")


#Preparo la matrice

edges$Type <- NULL
edges <- edges[!(edges$Subtype!="ACTIVATION" & edges$Subtype!="INHIBITION"),]
as.character(edges$Subtype)
edges$Subtype[edges$Subtype=="ACTIVATION"] <- "1"
edges$Subtype[edges$Subtype=="INHIBITION"] <- "-1"
edges$Subtype = as.numeric(as.character(edges$Subtype))
colnames(edges)[3] <- "weight"


#ed = data.frame(Genes=unique(as.character(unlist(edges[, c("X.Start", "End")])))); 
#nodesNew <- merge(ed, nodes, by = c(nodes[,1]), all=TRUE) MERGE DI DUE DATASET


graph <- graph_from_data_frame(edges, directed=TRUE)
#simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
#         edge.attr.comb = igraph_opt("edge.attr.comb"))

precision_matrix <- as_adjacency_matrix(graph, 
                                        names = TRUE,
                                        attr = "weight", 
                                        sparse=T)

precision_matrix <- precision_matrix[2100:2150, 2100:2150]

write.matrix(precision_matrix,
             file="matrix.txt",
             sep = "\t")

matrix <- read.delim("matrix.txt",
                     header = TRUE,
                     sep = "\t")


#CALCOLO DELLA PRECISION MATRIX

n <- nrow(precision_matrix) # Numero di vertici, inoltre numero righe = numero colonne nella precision_matrix


for(k in 1:n){
  for(i in 1:n){
    for(j in 1:n){
      if(precision_matrix[i,j]==0){
        precision_matrix[i,j] = (precision_matrix[i,k]*precision_matrix[k,j])
      }
    }
  }
}                   # ??(n^3) temporale ??(n^2) spaziale

# Creo il grafo di correlazione

correlation_graph <- graph_from_adjacency_matrix(precision_matrix, 
                                                 mode = "directed",
                                                 weighted = TRUE,
                                                 )

#Stampa del grafo di correlazione

plot.igraph(correlation_graph, edge.label=round(E(correlation_graph)$weight, 3), vertex.size=10)

plot( correlation_graph, layout = layout.reingold.tilford,
      edge.width = 1,
      edge.arrow.width = 0.3,
      vertex.size = 5,
      edge.arrow.size = 0.5,
      vertex.size2 = 3,
      vertex.label.cex = 1,
      asp = 0.35,
      margin = -0.1,
      edge.label=round(E(correlation_graph)$weight, 3),
      vertex.label=NA)
