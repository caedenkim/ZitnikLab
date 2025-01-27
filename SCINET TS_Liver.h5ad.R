library(anndata)
library(SCINET)
library(Matrix)


#########################
#
# Read TS data
#
#########################

liver <- read_h5ad("/Users/caeden/Downloads/TS_Liver.h5ad")
liver

# Meta data (optional)
liver$obs # for each cell
liver$var # for each gene

# Process raw count matrix
liver_mat = t(as.matrix(liver$X))
rownames(liver_mat) <- colnames(liver_mat) <- NULL # Remove all row and column names
liver_mat


#########################
#
# Pick a cell type
#
#########################

celltype_labels = liver$obs$cell_ontology_class

celltype = "monocyte"
celltype.cols = which(celltype_labels == celltype)

#########################
#
# Run SCINET
#
#########################

# Calculate activity scores

activity.scores = SCINET::compute_gene_activities(A = liver_mat, samples = celltype.cols, thread_no = 1)
rownames(activity.scores) = as.array(liver$var$gene_symbol)

# Load reference network

G = readRDS('/Users/caeden/Downloads/GAAD_net.RDS') # Link: https://github.com/shmohammadi86/SCINET/blob/master/demo/GAAD_net.RDS
paired.datasets = pair.datasets(G, activity.scores)
EL = get.edgelist(G, names = FALSE)
G.adj = as(get.adjacency(paired.datasets$net), 'dgTMatrix')
edge.idx = (G.adj@i)*nrow(G.adj) + (G.adj@j+1)

liver_cells = construct_cell_networks(net = G.adj, gene_activities = paired.datasets$activity.scores[, 1:10], thread_no = 1)
liver_cells

Mu.ref = Reduce("+", x = liver_cells) / length(liver_cells)
Mu.ref

Liver.withRef.graph = graph_from_adjacency_matrix(Mu.ref, mode = "undirected", weighted = TRUE)
Liver.withRef.graph

topo.spec.withRef = topo.spec(Liver.withRef.graph, sample_no = 100)
topo.spec.withRef

print(paired.datasets$genes[order(topo.spec.withRef, decreasing = TRUE)[1:100]]) #taking top 20 genes, change to top 50
