library(igraph)
library(matrix)
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

#Considero solo i geni presenti sia nella METAPATHWAY sia nel dataset di espressione.
{
  
  v1 <- rownames(expression)
  v2<- unique(c(edges$X.Start,edges$End))
  
  common <- c()
  for(i in 1:length(v1)){
    if(v1[i] %in% v2) common <- c(common,v1[i])
  }
  
  expression <- expression[common,]
  edges <- edges[(edges$X.Start %in% common) & (edges$End %in% common),]
  
  v1 <- rownames(expression)
  v2<- unique(c(edges$X.Start,edges$End))
  
  common <- c()
  for(i in 1:length(v1)){
    if(v1[i] %in% v2) common <- c(common,v1[i])
  }
  
  expression <- expression[common,]
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

#Introduco nuovamente i nomi dei geni
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



#Confronto con la METAPATHWAY (Non finito)



#Calcolo la matrice differenza tra la matrice di adiacenza della METAPATHWAY e tra la matrice di adiacenza del grafo di correlazione
diffmatrix <- adj_metapathway - precision_matrix
