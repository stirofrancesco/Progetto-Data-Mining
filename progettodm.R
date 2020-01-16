library(igraph)
library(DNetFinder)
library(stringi)

#Creazione del grafo di correlazione


#Lettura dataset di espressione
expression <- read.delim("expression_data.txt",
                          header = TRUE,
                          sep = "\t",
                          colClasses = "character")

#Caricamento nodi METAPATHWAY
nodes<- read.delim("nodes.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = "character")


#Caricamento archi METAPATHWAY
edges<- read.delim("edges.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = "character")


#Tengo solamente Attivazione e Inibizione e assegno i pesi negli archi della METAPATHWAY.
{
edges$Type <- NULL
edges <- edges[!(edges$Subtype!="ACTIVATION" & edges$Subtype!="INHIBITION"),]
edges$Subtype[edges$Subtype=="ACTIVATION"] <- "1"
edges$Subtype[edges$Subtype=="INHIBITION"] <- "-1" 
edges$Subtype = as.numeric(as.character(edges$Subtype))
colnames(edges)[3] <- "weight"
}

#Considero solo i geni presenti in tutti i dataset.
{
  
  v1 <- rownames(expression)
  v2 <- unique(c(edges$X.Start,edges$End))
  v3 <- nodes$X.Id
  
  common <- c()
  for(i in 1:length(v1)){
    if(v1[i] %in% v2 && v1[i] %in% v3) common <- c(common,v1[i])
  }
  
  expression <- expression[common,]
  edges <- edges[(edges$X.Start %in% common) & (edges$End %in% common),]
  nodes <- nodes[nodes$X.Id %in% common,]
  
  v1 <- rownames(expression)
  v2<- unique(c(edges$X.Start,edges$End))
  v3 <- nodes$X.Id
  
  common <- c()
  for(i in 1:length(v1)){
    if(v1[i] %in% v2 && v1[i] %in% v3) common <- c(common,v1[i])
  }
  
  expression <- expression[common,]
  nodes <- nodes[nodes$X.Id %in% common,]
  
  v1 <- rownames(expression)
  v3 <- nodes$X.Id
}

#Creazione grafo METAPATHWAY
metapathway <- graph_from_data_frame(edges,
                                     directed=TRUE)

#Creazione matrice di adiacenza METAPATHWAY
adj_metapathway <- as.matrix(as_adjacency_matrix(metapathway))


# CREAZIONE DEL GRAFO DI CORRELAZIONE

#Sampling delle matrici
{
  numrow <- 300
  numcol <- 300
}

#Carico la matrice di espressione dal dataset
expression_matrix <- data.matrix(expression)

# Calcolo della trasposta della matrice di espressione
expression_matrix <- t(expression_matrix)

#Ordino i geni sulla base della metapathway
expression_matrix <- expression_matrix[,rownames(adj_metapathway)]

#Riduco numero di righe e colonne
{
  expression_matrix <- expression_matrix[1:numrow,1:numcol]
  adj_metapathway <- adj_metapathway[1:numrow,1:numcol]
}

#Calcolo i coefficenti di regressione (scegliere tra lassoGGM() e lassoNPN() )
coefGGM <- lassoNPN(expression_matrix)

#Setto a 0 le colonne in cui lasso non ha restituito un valore.
coefGGM[is.na(coefGGM)] <- 0

#Creo la matrice di adiacenza, inizialmente con tutti 0.
precision_matrix <- matrix(0,nrow = numcol, ncol = numcol)

#Introduco nuovamente gli ID dei geni
{
rownames(precision_matrix) <- colnames(expression_matrix)
colnames(precision_matrix) <- colnames(expression_matrix)
}

#Riempimento della matrice di precisione.
for(i in 1:numcol){
  for(j in 1:numcol){
    if(i<j){
      if(coefGGM[i,j]>0) precision_matrix[i,j] = 1
      else if(coefGGM[i,j]<0) precision_matrix[i,j] = -1
    }
    else if(i>j){
      if(coefGGM[i-1,j]>0) precision_matrix[i,j] = 1
      else if(coefGGM[i-1,j]<0) precision_matrix[i,j] = -1
    }  
  }
}

#Creazione del grafo di correlazione
correlation_graph <- graph_from_adjacency_matrix(precision_matrix, 
                                                  mode = "directed",
                                                  weighted = TRUE
                                                )

#Stampa del grafo ci correlazione
plot(correlation_graph)



#Confronto con la METAPATHWAY


#Calcolo la matrice differenza tra la matrice di adiacenza della METAPATHWAY e tra la matrice di adiacenza del grafo di correlazione
diffmatrix <- adj_metapathway - precision_matrix

#Setto a 0 i self-loop
for(i in 1:length(rownames(diffmatrix))){
  diffmatrix[i,i] <- 0
}

#Creazione delle reti per la divisione dei pesi
{
  network1 <- matrix(0, nrow = length(which(diffmatrix==1)), ncol = 2 )
  network2 <- matrix(0,nrow = length(which(diffmatrix==2)), ncol = 2 )
  network_neg1 <- matrix(0,nrow = length(which(diffmatrix==-1)), ncol = 2 )
  network_neg2 <- matrix(0,nrow = length(which(diffmatrix==-2)), ncol = 2 )
  colnames(network1) <- c("Start","End")
  colnames(network2) <- c("Start","End")
  colnames(network_neg1) <- c("Start","End")
  colnames(network_neg2) <- c("Start","End")
}

#Riempimento di ciascuna matrice per ogni relativo peso
{
  index <- c(1,1,1,1)
  for(i in 1:length(rownames(diffmatrix))){
    for (j in 1:length(colnames(diffmatrix))){
      if(diffmatrix[i,j]==1){
        network1[index[1],] <- c(rownames(diffmatrix)[i],colnames(diffmatrix)[j])
        index[1] <- index[1]+1
      }
      else if(diffmatrix[i,j]==2){
        network2[index[2],] <- c(rownames(diffmatrix)[i],colnames(diffmatrix)[j])
        index[2] <- index[2]+1
      }
      else if(diffmatrix[i,j]==-1){
        network_neg1[index[3],] <- c(rownames(diffmatrix)[i],colnames(diffmatrix)[j])
        index[3] <- index[3]+1
      }
      else if(diffmatrix[i,j]==-2){
        network1[index[4],] <- c(rownames(diffmatrix)[i],colnames(diffmatrix)[j])
        index[4] <- index[4]+1
      }
    }  
  }
}


#Conversione delle matrici a DataFrame
{
  network1 <- as.data.frame(network1)
  network2 <- as.data.frame(network2)
  network_neg1 <- as.data.frame(network_neg1)
  network_neg2 <- as.data.frame(network_neg2)
}

#Grafo rete pesi 1
graph_net1 <- graph_from_data_frame(network1, 
                                    directed = TRUE, 
                                    vertices = nodes[nodes$X.Id %in% unique(c(as.vector.factor(network1$Start),as.vector.factor(network1$End))),])
#Grafo rete pesi 2
graph_net2 <- graph_from_data_frame(network2, 
                                    directed = TRUE, 
                                    vertices = nodes[nodes$X.Id %in% unique(c(as.vector.factor(network2$Start),as.vector.factor(network2$End))),])
#Grafo rete pesi -1
graph_neg_net1 <- graph_from_data_frame(network_neg1, 
                                    directed = TRUE, 
                                    vertices = nodes[nodes$X.Id %in% unique(c(as.vector.factor(network_neg1$Start),as.vector.factor(network_neg1$End))),])
#Grafo rete pesi -2
graph_neg_net2 <- graph_from_data_frame(network_neg2, 
                                    directed = TRUE, 
                                    vertices = nodes[nodes$X.Id %in% unique(c(as.vector.factor(network_neg2$Start),as.vector.factor(network_neg2$End))),])


#Funzione per la normalizzazione
rescale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}

#Visualizzazione Grafo 1
{
  deg <- degree(graph_net1, mode = "all")
  V(graph_net1)$size <- rescale(deg,min(deg),max(deg),2,35)
  E(graph_net1)$arrow.size <- .2

  tkplot(graph_net1, vertex.label = ifelse(degree(graph_net1) >= 13, sapply(strsplit(c(V(graph_net1)$Name),","), `[`, 1), NA) )
}

#Visualizzazione Grafo 2
{
  deg <- degree(graph_net2, mode = "all")
  V(graph_net2)$size <- rescale(deg,min(deg),max(deg),2,35)
  E(graph_net2)$arrow.size <- .2
  
  tkplot(graph_net2, vertex.label = ifelse(degree(graph_net2) >= 13, sapply(strsplit(c(V(graph_net2)$Name),","), `[`, 1), NA) )
}

plot(graph_net2)

#Visualizzazione Grafo -1
{
  deg <- degree(graph_neg_net1, mode = "all")
  V(graph_neg_net1)$size <- rescale(deg,min(deg),max(deg),2,35)
  E(graph_neg_net1)$arrow.size <- .2
  
  tkplot(graph_neg_net1, vertex.label = ifelse(degree(graph_neg_net1) >= 13, sapply(strsplit(c(V(graph_neg_net1)$Name),","), `[`, 1), NA) )
}

#Visualizzazione Grafo -2
{
  deg <- degree(graph_neg_net2, mode = "all")
  V(graph_neg_net2)$size <- rescale(deg,min(deg),max(deg),2,35)
  E(graph_neg_net2)$arrow.size <- .2
  
  tkplot(graph_neg_net2, vertex.label = ifelse(degree(graph_neg_net2) >= 13, sapply(strsplit(c(V(graph_neg_net2)$Name),","), `[`, 1), NA) )
}

