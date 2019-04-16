library(igraph)
library(MCODER)

### draw an whole network (directed), Example 1
adj_mat_directed <- adjacency.matrix(raw_table2,v1='outdegree_node',v2='indegree_node',direct=T,multiple.edge=F)
# vertex value table (label,K-core,node density, node score)
clusters <- list()
clusters[[1]] <- rownames(adj_mat_directed)
clusters_edges <- list()
clusters_edges[[1]] <- adj_mat_directed
##prepare for essential plotting input
graph_par_vertex <- graph.par.vertex(x=clusters)
graph_par_edge <- graph.par.edge(x=clusters_edges)
# transform to network
IG <- network.format(x=clusters_edges,vertices.par=graph_par_vertex,edges.par=graph_par_edge,direct=T)
# plotting
plot(plot.network(IG[1])



### draw clusters with multiple edges, example 2
adj_mat_directed <- adjacency.matrix(raw_table2,v1='outdegree_node',v2='indegree_node',direct=F,multiple.edge=T)
vertex_table <- calc.vertex.value(x=adj_mat_directed,loop=F)

# finding clusters
clusters <- find.cluster(vertex_table,adj_mat_directed,node.score.cutoff=0.2)

# Fluff
clusters <- post.fluff(x=clusters,adj_mat_directed,vertex_table)
clusters <- post.fluff.bind(x=clusters)

# haircut
clusters <- post.haircut(clusters,adj_mat_directed)

# edge matrix order by cluster number
clusters_edges <- cluster.adj.matrix(x=clusters,adj_mat_directed)

#prepare for essential plotting input
graph_par_vertex <- graph.par.vertex(x=clusters)
graph_par_edge <- graph.par.edge(x=clusters_edges) #option here for multiple edge, direction

# network format
IG <- network.format(x=clusters_edges,vertices.par=graph_par_vertex,edges.par=graph_par_edge,direct=T)

# plotting
for(i in 1:length(IG)) { plot.network(IG[i]) }

