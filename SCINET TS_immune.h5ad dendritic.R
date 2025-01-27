# https://github.com/shmohammadi86/SCINET/blob/master/demo/Run_SCINET_PBMC.Rmd

library(anndata)
library(SCINET)
library(Matrix)


# Read TS data

immune <- read_h5ad("/Users/caeden/Downloads/TS_immune.h5ad")
immune

# Meta data (optional)
immune$obs
immune$var

unique(immune$obs["cell_ontology_class"])

for (i in t(unique(immune$obs["cell_ontology_class"]))) {
  print(paste0("Working on... ", i))
  if("dendritic cell" == i){
    
    
    # Write all your processing code below, in this for loop
    
    # Take a subset of the matrix
    immune_subset = immune[immune$obs["cell_ontology_class"] == i]
    
    immune_subset_mat = t(as.matrix(immune_subset$X)) 
    immune_subset_mat
  }
}

immune_subset$obs
rownames(immune_subset_mat) <- colnames(immune_subset_mat) <- NULL
immune_subset_mat

celltype_labels = immune_subset$obs$cell_ontology_class
celltype = "nk cell"
celltype.cols = which(celltype_labels == celltype)

activity.scores = SCINET::compute_gene_activities(A = immune_subset_mat, samples = celltype.cols, thread_no = 1)
rownames(activity.scores) = as.array(immune_subset$var$gene_symbol)

activity.scores = SCINET::compute_gene_activities(A = immune_subset_mat, samples = celltype.cols, thread_no = 1)
rownames(activity.scores) = as.array(immune_subset$var$gene_symbol)

G = readRDS('/Users/caeden/Downloads/GAAD_net.RDS') # Link: https://github.com/shmohammadi86/SCINET/blob/master/demo/GAAD_net.RDS
paired.datasets = pair.datasets(G, activity.scores)
EL = get.edgelist(G, names = FALSE)
G.adj = as(get.adjacency(paired.datasets$net), 'dgTMatrix')
edge.idx = (G.adj@i)*nrow(G.adj) + (G.adj@j+1)

immune_cells = construct_cell_networks(net = G.adj, gene_activities = paired.datasets$activity.scores[, 1:3], thread_no = 1)
immune_cells

Mu.ref = Reduce("+", x = immune_cells) / length(immune_cells)
Mu.ref

Immune.withRef.graph = graph_from_adjacency_matrix(Mu.ref, mode = "undirected", weighted = TRUE)
Immune.withRef.graph

topo.spec.withRef = topo.spec(Immune.withRef.graph, sample_no = 100)
topo.spec.withRef

print(paired.datasets$genes[order(topo.spec.withRef, decreasing = TRUE)[1:20]]) #taking top 20 genes, change to top 50

