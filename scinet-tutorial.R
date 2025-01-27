# https://github.com/shmohammadi86/SCINET/blob/master/demo/Run_SCINET_PBMC.Rmd

library(anndata)
library(SCINET)
library(Matrix)


#########################
#
# Read tutorial data
#
#########################

sce = readRDS('/Users/michelle/caeden_kim/reduced_sce_PBMC_annotated.RDS') # Link: https://github.com/shmohammadi86/SCINET/blob/master/demo/reduced_sce_PBMC_annotated.RDS
ACTIONet.out = readRDS('/Users/michelle/caeden_kim/ACTIONet_out_PBMC.RDS') # Link: https://github.com/shmohammadi86/SCINET/blob/master/demo/ACTIONet_out_PBMC.RDS

Bcell.cols = which(sce$Labels == "B cell")
marker.genes = c('CD19', 'CD22', 'CD4', 'CD8A', 'CD14')
marker.rows = match(marker.genes, rownames(sce))

W = ACTIONet.out$signature.profile[, ACTIONet.out$core.out$core.archs]
H = ACTIONet.out$reconstruct.out$H_stacked[ACTIONet.out$core.out$core.archs, ]
A =  W %*% H


#########################
#
# SCINET
#
#########################

# Calculate activity scores

activity.scores = SCINET::compute_gene_activities(A = A, samples = Bcell.cols, thread_no = 1)

# Format activity scores
X = t(activity.scores[marker.rows, ])
colnames(X) = marker.genes
df = reshape2::melt(X)
colnames(df) = c('id', 'gene', 'activity')

require(ggplot2)
ggplot(df, aes(x=activity, fill=gene)) +
  geom_density()

rownames(activity.scores) = as.array(rownames(sce))

# Load reference network

G = readRDS('/Users/michelle/caeden_kim/GAAD_net.RDS') # Link: https://github.com/shmohammadi86/SCINET/blob/master/demo/GAAD_net.RDS

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

print(paired.datasets$genes[order(topo.spec.withRef, decreasing = TRUE)[1:20]])

