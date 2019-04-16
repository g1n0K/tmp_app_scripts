# Gino Sungjin, Kwon 
# gino.kwon@gmail.com
# Yonsei Genomic center, Biomedical research institute, Yonsei university medical center

#### MCODE library
# import required packages
if(requireNamespace("sna",quietly = TRUE)) {
	kcores <- sna::kcores
} else {
	"Required package: SNA"
}
if(requireNamespace("igraph",quietly = TRUE)) {
	library(igraph)
} else {
	"Required package: igraph"
}

# shape adjust
shapes <- function() {
mycircle <- function(coords, v=NULL, params) {
	vertex.color <- params("vertex", "color")
	if (length(vertex.color) != 1 && !is.null(v)) { vertex.color <- vertex.color[v] }
	vertex.size  <- 1/250 * params("vertex", "size")
	if (length(vertex.size) != 1 && !is.null(v)) { vertex.size <- vertex.size[v] }
	vertex.frame.color <- params("vertex", "frame.color")
	if (length(vertex.frame.color) != 1 && !is.null(v)) { vertex.frame.color <- vertex.frame.color[v] }
	vertex.frame.width <- params("vertex", "frame.width")
	if (length(vertex.frame.width) != 1 && !is.null(v)) { vertex.frame.width <- vertex.frame.width[v] }
	mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,vertex.size, vertex.frame.width,FUN=function(x, y, bg, fg, size, lwd) { symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,circles=size, add=TRUE, inches=FALSE) })
}
mysquare <- function(coords, v=NULL, params) {
	vertex.color <- params("vertex", "color")
	if (length(vertex.color) != 1 && !is.null(v)) { vertex.color <- vertex.color[v] }
	vertex.size  <- 1/150 * params("vertex", "size")
	if (length(vertex.size) != 1 && !is.null(v)) { vertex.size <- vertex.size[v] }
	vertex.frame.color <- params("vertex", "frame.color")
	if (length(vertex.frame.color) != 1 && !is.null(v)) { vertex.frame.color <- vertex.frame.color[v] }
	vertex.frame.width <- params("vertex", "frame.width")
	if (length(vertex.frame.width) != 1 && !is.null(v)) { vertex.frame.width <- vertex.frame.width[v] }
	mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color, vertex.size, vertex.frame.width,FUN=function(x, y, bg, fg, size, lwd) { symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,squares=size, add=TRUE, inches=FALSE) })
}
add.vertex.shape("circle", clip=igraph.shape.noclip, plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=2))
add.vertex.shape("square", clip=igraph.shape.noclip, plot=mysquare, parameters=list(vertex.frame.color=1, vertex.frame.width=2))
}


### build vertex(node) matrix
# 'input_tb' is vertex information table
# first and second columns are candidates, otherwise set names of target column.
# this matrix is upper triangle, undirected graph. plus, ignore self loop and duplicated edges
adjacency.matrix <- function(input_tb,v1=NULL,v2=NULL,sampling=NULL,remove.loop=T,directed=F,multiple.edge=F) {
	samples <- sampling
	# get vertex column location
	if( (length(v1) >= 1) & (length(v2) >= 1) ) {
		loc1 <- which(colnames(input_tb) == v1)
		loc2 <- which(colnames(input_tb) == v2)
	} else {
		loc1 <- 1
		loc2 <- 2
	}
	# sampling from raw data, char
	if( length(samples) >= 1 ) {	input_tb <- input_tb[which( (input_tb[,loc1] %in% samples) & (input_tb[,loc2] %in% samples) ),]	}
	
	# remove isolated vertex(null) and remove duplicated nodes
	# Self-loop removing
	if(remove.loop == T) {
		del_list <- which( unlist( lapply( 1:nrow(input_tb), function(i) if( ( is.na(input_tb[i,loc1]) ) | ( is.na(input_tb[i,loc2]) ) | (input_tb[i,loc1] == input_tb[i,loc2]) ) {'DEL'} else {'N'} ) ) == 'DEL' )
		if( length(del_list) >= 1 ) { input_tb <- input_tb[-del_list,] }
	}
	# vertex_labels(list)
	vertex_labels <- sort(unique(append(input_tb[,loc1],input_tb[,loc2])))
	# make empty matrix for vertex adjacency matrix
	adjacency_matrix <- matrix(0,nrow=length(vertex_labels),ncol=length(vertex_labels))
	# make flag on adjacency matrix, 1 = edge
	id1 <- unlist( lapply(1:nrow(input_tb), function(x) which(input_tb[x,loc1] == vertex_labels) ) )
	id2 <- unlist( lapply(1:nrow(input_tb), function(x) which(input_tb[x,loc2] == vertex_labels) ) )
	# directed or indirected
	for( i in 1:length(id1) ) {
		if(directed == T) {
			adjacency_matrix[id1[i],id2[i]] <- adjacency_matrix[id1[i],id2[i]] + 1
		} else {
			tmp_id1 <- sort(c(id1[i],id2[i]))[1]
			tmp_id2 <- sort(c(id1[i],id2[i]))[2]
			if( tmp_id1 == tmp_id2 ) {
				adjacency_matrix[tmp_id1,tmp_id2] <- adjacency_matrix[tmp_id1,tmp_id2] + 1
			} else {
			adjacency_matrix[tmp_id1,tmp_id2] <- adjacency_matrix[tmp_id1,tmp_id2] + 1
			adjacency_matrix[tmp_id2,tmp_id1] <- adjacency_matrix[tmp_id2,tmp_id1] + 1
			}
		}
	}
	# remove multiple edges
	if( multiple.edge == F ) {
		adjacency_matrix[which(adjacency_matrix != 0)] <- 1
	}
	rownames(adjacency_matrix) <- vertex_labels
	colnames(adjacency_matrix) <- vertex_labels
	return(adjacency_matrix)
}

### staring with adj matrix -> nodes inf tables
# input - adjc square table with labels, minimum degree cut-off
# calc density(weight), score, local k-cores
calc.node.scores <- function(x,adj.mat,minimum.edges=2,loop=F) {
	seed <- x
	flag <- adj.mat
	flag[which(flag > 0)] <- 1
	if(minimum.edges < 2) {minimum.edges <- 2} # minimum edge must be at least 2
	min_edges_num <- minimum.edges
	if( (degree_cutoff != min_edges_num) & (degree_cutoff >= 3) ) { min_edges_num <- degree_cutoff } else if( (degree_cutoff != min_edges_num) & (degree_cutoff < 2) ) { stop('Minimum Degree-cutoff is 2') }
	edges <- flag
	seed <- as.character(seed)
	id <- which(rownames(edges) == seed)
	# get neighbours nodes
	neighbours <- append(id,which(edges[id,] != 0))
	#tar <- append(id,which(edges[id,] != 0))
	
	# minimum edge cut_off
	if( (sum(edges[id,]) <= min_edges_num-1) && (sum(edges[id,]) >= 1) ) {
		Den <- 0
		Sco <- 0
		K <- 1
	} else if(sum(edges[id,]) == 0) {  
		Den <- 0
		Sco <- 0
		K <- 0
	} else {
		tmp_mat <- flag[neighbours,neighbours]
		# calc nighbour k-cores, then cut out. ks[1] is seed.
		ks <- kcores(tmp_mat,mode = "graph")
		tmp_mat <- tmp_mat[which(ks >= ks[1]),which(ks >= ks[1])]
		# result1 is K-core value by the SNA package
		K <- as.numeric(ks[1])
		# calc edge sum,node sum Act numdes, Poss nodes, density, score
		E <- sum( as.matrix(tmp_mat) ) / 2
		V <- nrow(tmp_mat)
		Den <- (2 * E) / (V*(V-1)) #exclude loop
		if( loop==T ) { Den <- (2 * E) / (V*(V+1)) } #include loop
		Sco <- Den * K
	}
	
	rst <- list(Den,Sco,K)

	names(rst) <- c('density','score','k')
	return(rst)
}

## node table, calculating vertex values: k, density, score
# flag = matrix
calc.vertex.value <- function(x,loop=F,minimum.edges=2) {
	flag <- x
	flag[which(flag > 0)] <- 1
	# input needs undirect edge matrix with row and col labels
	if(is.matrix(flag) == FALSE) {stop('Require matrix of edge adjacency matrix')}
	# vertexs is verte value table which includes k-core,density,score
	vertexs <- as.data.frame(rownames(flag))
	colnames(vertexs) <- 'Label'
	vertexs$K_core <- as.numeric( kcores(flag,mode = "graph") ) # k-cores undirect, MCODE is undirected graph algorithm
	
	# function calculate node score
	node_values <- lapply(1:nrow(vertexs), function(x) calc.node.scores(vertexs$Label[x],adj.mat=flag,minimum.edges=minimum.edges,loop=loop) )
	vertexs$density <- sapply(1:nrow(vertexs), function(i) node_values[[i]]$density )
	vertexs$score <- sapply(1:nrow(vertexs), function(i) node_values[[i]]$score )
	vertexs$local_k <- sapply(1:nrow(vertexs), function(i) node_values[[i]]$k )
	return(vertexs)
}

## finding clusters
find.cluster <- function(vertex.score,adj.mat,node.score.cutoff=0.2,min.k.core=2,max.depth=100) {
	vertex <- vertex.score
	flag <- adj.mat
	flag[which(flag > 0)] <- 1
	cut_off <- node.score.cutoff
	min_core <- min.k.core
	max_depth <- max.depth
	ids <- vertex$score[which(vertex$local_k != 0)]
	names(ids) <- which(vertex$local_k != 0)

	ids <- sort(ids,decreasing = T)
	tmp_df <- as.data.frame(ids)
	colnames(tmp_df) <- 'score'
	# make table for cluster target
	tmp_df$id <- names(ids)
	#tmp_df$symbol <- vertex$NODE[as.numeric(names(ids))]
	tmp_df$label <- vertex$Label[as.numeric(names(ids))]
	# -99 as NULL
	tmp_df$cluster <- -99
	# fr the highest to lowest
	z = 1
	while( length( which(tmp_df$cluster < 0) ) != 0 )
	{
		# get a seed node. which has not cluster number
		i <- which( (tmp_df$score == max( tmp_df$score[which(tmp_df$cluster < 0 ) ] ) ) & (tmp_df$cluster < 0) )[1]
		
		seed <- as.numeric( tmp_df$id[i]) # vertex number
		score_cut <- (1-cut_off) * tmp_df$score[i] # vertex score cut-off
		
		### clustering step
		# neighbour list
		tar <- append(as.numeric( tmp_df$id[i] ) ,which( flag[as.numeric( tmp_df$id[i] ),] != 0) )
		loc <- which(tmp_df$id %in% tar )
		# rm which already did
		loc <- loc[ which(tmp_df$cluster[loc] < 0) ]
		# score cut off
		pass <- loc[which( tmp_df$score[loc] >= score_cut)]
		if( length(pass) >= 1 ) 
		{
			# give cluster num
			tmp_df$cluster[pass] <- z 
			# loop until max depth
			for( j in 1:(max_depth-1))
			{
				tar <- unique (unlist( lapply(1:length(pass), function(x) which(flag[as.numeric(tmp_df$id[pass[x]]),] != 0 ) ) ) )
				tar <- tar[which(!tar %in% as.numeric( tmp_df$id[pass] ) ) ]
				loc <- which(tmp_df$id %in% tar )
				loc <- loc[ which( tmp_df$cluster[loc] < 0 ) ]
				pass <- loc[which( tmp_df$score[loc] >= score_cut)]
				if( length(pass) >= 1 ) { tmp_df$cluster[pass] <- z }
				if( length(pass) < 1 ) { break }
			}
		}
		z <- z+1
	}
	# list up cluster list
	cluster_ls <- unique(tmp_df$cluster)
	tmp_len <- unlist(lapply(1:length(cluster_ls), function(x) length(which(cluster_ls[x] == tmp_df$cluster))))
	# minimum core. number of node > core + 1. 
	cluster_ls <- cluster_ls[which( tmp_len > min_core )]
	
	rst <- lapply(1:length(cluster_ls), function(x) sort( as.character( tmp_df$label[ which( cluster_ls[x] == tmp_df$cluster) ] ) ) )
	i=1
	del.ls <- rep(0,length(cluster_ls))
	while( i <= length(cluster_ls))
	{
		ids <- as.numeric( tmp_df$id[ which( tmp_df$cluster %in% cluster_ls[i] ) ] )
		mat <- flag[ids,ids]
		k <- kcores(mat,mode = "graph")
		if( max(k) < min_core )	{ del.ls[i] <- 1 }
		i <- i+1
	}
	#rst <- rst[-which(del.ls == 1)]
	rst <- rst[which(del.ls != 1)]
	rst <- rst[lapply(rst,length)>0]
	# get cluster list
	return(rst)
}
# EOF 


### Post PROC, HAIRCUT
post.haircut <- function(x,adj.mat) {
	list_clust <- x
	vertex <- adj.mat
	vertex[which(vertex > 0)] <- 1
	diag(vertex) <- 0
	len <- length(list_clust)
	for( z in 1:len )
	{
		components <- list_clust[[z]]
		ids <- which(rownames(vertex) %in% components)
		mat <- vertex[ids,ids]
		k <- kcores(mat,mode = "graph")
		# siginly connected node equals K=1
		del.ls <- names(k)[which( k == 1 )]
		if( length(del.ls) >= 1 )
		{
			components <- components[which(!components %in% del.ls)]
			list_clust[[z]] <- components
		}
	}
	return(list_clust)
}
# EOF 


### Post PROC, FLUFF
post.fluff <- function(x,adj.mat,vertex.score.table,density.cutoff=0.6) {
	cut_off <- density.cutoff
	flag <- adj.mat
	flag[which(flag > 0)] <- 1
	diag(flag) <- 0
	list_clust <- x
	len <- length(list_clust)
	node_vals <- vertex.score.table
	node_vals$cluster <- 0
	for( z in 1:len )
	{
		node_vals$cluster[ which( node_vals$Label %in% list_clust[[z]] ) ] <- z
	}
	for( z in 1:len )
	{
		components <- list_clust[[z]]
		ids_ver <- which( rownames(flag) %in% components )
		tmp_mat <- as.matrix(flag[ids_ver,])
		# get neighbours of seed cluster except other clusters and itself
		neighbours <- unique( unlist( lapply(1:nrow(tmp_mat), function(x) colnames(tmp_mat)[ which( tmp_mat[x,] != 0 ) ] ) ) )
		neighbours <- neighbours[ which(!neighbours %in% components) ]
		neighbours <- neighbours[ which( neighbours %in% as.character( node_vals$Label[ which( node_vals$cluster == 0 ) ] ) ) ]

		# get table of possible fluff nodes
		DF <- node_vals[ which(node_vals$Label %in% neighbours), ]
		# calc local Density for fluff
		lc_dens <- unlist( lapply(1:nrow(DF), function(x) calc.local.density(DF$Label[x],flag) ) )

		# find Fluff target (bigger then Density cut off)
		add_vertex_by_fluff <- as.character( DF$Label[ which(lc_dens > cut_off) ] )
		if( length(add_vertex_by_fluff) >= 1 )
		{
			components <- sort( append( components,add_vertex_by_fluff ) )
			list_clust[[z]] <- components
		}
	}
	return(list_clust)
}
# EOF 

### after FLUFF, some case has common vertexes. but, Cytoscape App MCODE igonore this case and print as each different groups
post.fluff.bind <- function(x) {
	nodes <- as.data.frame(unlist(x),stringsAsFactors=F)
	colnames(nodes) <- 'nodes'
	tmp <- sapply(1:length(x), function(i) length(x[[i]]) )
	nodes$num <- unlist( lapply(1:length(tmp) , function(i) rep(i,tmp[i]) ) )
	nodes$dup <- duplicated(nodes[,1])
	if( length(which(nodes$dup == 'TRUE') ) >= 1 ) {
		dup_id <- which(nodes$dup == 'TRUE')
		dup_num <- nodes$num[dup_id]
		#nodes$dup[ which( nodes$nodes %in% nodes$nodes[which(nodes$dup == 'TRUE')] ) ] <- 'TRUE'
		for(i in 1:length(dup_id)) {
			tmp <- which(nodes$nodes == nodes$nodes[dup_id[i]] )
			tmp <- tmp[which(tmp != dup_id[i] ) ]
			group_id <- nodes$num[tmp]
			nodes$num[ which(nodes$num == dup_num[i] ) ] <- group_id
		}
		nodes <- nodes[which(nodes$dup != 'TRUE'),]
		rst <- lapply(1:length(unique(nodes$num)), function(x) nodes$nodes[nodes$num == unique(nodes$num)[x]] )
	}	else {
		rst <- x
	}
	return(rst)
}


### local density
calc.local.density <- function(x,vertex) {
	seed <- x
	seed_id <- which( rownames(vertex) == seed )
	ids <- as.numeric( which(vertex[ seed_id, ] != 0) )
	tmp_df <- as.matrix( vertex[append(seed_id,ids) , append(seed_id,ids)] )
	Den <- (2 * ( sum( as.matrix(tmp_df) ) / 2) ) / (nrow(tmp_df)*(nrow(tmp_df)-1)) 
	return(Den)
}
# EOF 

### make individual adjacency matrix by clusters
# att edge matrix to cluster list
cluster.adj.matrix <- function(x,adj.mat) {
	cls <- x
	ver <- adj.mat
	rst <- list()
	for(i in 1:length(cls))
	{
		ids <- sapply(1:length(cls[[i]]), function(x) which( rownames(ver) == cls[[i]][x] ) )
		mat <- ver[ids,ids]
		rst[[i]] <- mat
	}
	return(rst)
}
# EOF 

# plot input table for nodes
graph.par.vertex <- function(x) {
	if(length(x) < 1) {stop('Need clusters list')}
	tmp <- as.data.frame(unlist(x),stringsAsFactors=F)
	colnames(tmp) <- 'x'
	tmp$number <- NA
	for( i in 1:length(x) ) {
		ids <- which( tmp$x %in% unlist(x[[i]]) )
		tmp$number[ids] <- i
	}
	tmp$colour <- 'white' # colour of vertex frame; inside or vertex. defualt = 'white'
	tmp$size <- 27 # vertex frame size, defualt = 27
	tmp$font <- 1.0 #vertex label font size, defualt = 1
	tmp$shape <- 'circle' # vertex shape, defualt = 'circle'
	tmp$lab_col <- 'black' # label colour, defualt = 'black'
	tmp$border <- 'red' # vertex border colour, defualt = 'red'
	return(tmp)
}

# edges information data frame build , input is 'clusters.edges'
graph.par.edge <- function(x,directed=F,multiple.edge=F) {
	for( i in 1:length(x)) {
		# get cluster's edge information
		tmp <- x[[i]]
		labs <- rownames(tmp)
		
		if(directed==F) {
		tmp[which(lower.tri(tmp))] <- 0
		ids <- which(tmp != 0)
		rvs <- ids %% nrow(tmp)
		cvs <- as.integer(ids/nrow(tmp)) + 1
		if(multiple.edge==T) {
			len <- tmp[ids]
			rvs <- unlist(lapply(1:length(rvs), function(k) rep(rvs[k],len[k]) ) )
			cvs <- unlist(lapply(1:length(cvs), function(k) rep(cvs[k],len[k]) ) )
		}
		tmp_rst <- cbind( as.data.frame(labs[rvs],stringsAsFactors=F),as.data.frame(labs[cvs],stringsAsFactors=F) )
		tmp_rst$number <- i
		colnames(tmp_rst)[1:2] <- c('x','y')
		} else {
		rvs <- NA
		cvs <- NA
			for(j in 1:nrow(tmp)) {
				vs <- which(tmp[j,] != 0)
				if(length(vs) >= 1) {
					x_node <- rep(labs[j],length(vs))
					y_node <- labs[vs]
					if(multiple.edge==T) {
						len <- tmp[j,vs]
						x_node <- unlist(lapply(1:length(x_node), function(k) rep(x_node[k],len[k]) ) )
						y_node <- unlist(lapply(1:length(y_node), function(k) rep(y_node[k],len[k]) ) )
					}
					rvs <- append(rvs,x_node)
					cvs <- append(cvs,y_node)
				}
			}
		rvs <- rvs[-1]
		cvs <- cvs[-1]
		tmp_rst <- cbind( as.data.frame(rvs,stringsAsFactors=F),as.data.frame(cvs,stringsAsFactors=F) )
		tmp_rst$number <- i
		colnames(tmp_rst)[1:2] <- c('x','y')
		}
		if( i == 1 ) {rst <- tmp_rst} else {rst <- rbind(rst,tmp_rst)}
	}
	
	# set edges graphical options, defualt = 2
	rst$colour = 'blue' #edge colour, defualt = 'black'
	rst$width = 1 #edge width, defualt = 2

	return(rst)
}

### Add annotation function
add.annotation <- function(x,vertices.par,edges.par,annotation.table) {
	if( "FALSE" %in% c(exists("vertices.par"),exists("x"),exists("edges.par"),exists("annotation.table")) ) { stop("cluster.edges, vertex.table, edges.tables and annotation tables are required") }
	for(i in 1:length(x) ) { if( length(which( rownames(x[[i]]) %in% annotation.table[,1]) ) >= 1 ) {	
			add <- which(rownames(x[[i]]) %in% annotation.table[,1])
			add_desc <- annotation.table[which(annotation.table[,1] %in% rownames(x[[i]])[add] ),2]
			tmp_mat <- matrix(0,nrow=nrow(x[[i]])+length(add_desc),ncol=nrow(x[[i]])+length(add_desc) )
			colnames(tmp_mat) <- append(rownames(x[[i]]),add_desc)
			rownames(tmp_mat) <- append(rownames(x[[i]]),add_desc)
			for(j in 1:length(add_desc)) {
				tmp_mat[add[j],(nrow(x[[i]])+j) ] <- 1
				tmp_mat[(nrow(x[[i]])+j),add[j] ] <- 1
				vertices.par[nrow(vertices.par)+1,] <- NA
				vertices.par[nrow(vertices.par),1] <- add_desc[j]
				vertices.par[nrow(vertices.par),2] <- i
				vertices.par[nrow(vertices.par),3] <- 'white'
				vertices.par[nrow(vertices.par),4] <- 18
				vertices.par[nrow(vertices.par),5] <- 0.6
				vertices.par[nrow(vertices.par),6] <- 'rectangle'
				vertices.par[nrow(vertices.par),7] <- 'black'
				vertices.par[nrow(vertices.par),8] <- 'black'
				edges.par[nrow(edges.par)+1,] <- NA
				edges.par[nrow(edges.par),3] <- i
				edges.par[nrow(edges.par),2] <- rownames(x[[i]])[add[j]]
				edges.par[nrow(edges.par),1] <- add_desc[j]
				edges.par[nrow(edges.par),4] <- 'black'
				edges.par[nrow(edges.par),5] <- 1		
			}
			for(j in 1:nrow(x[[i]]) ) {
				for(k in j:nrow(x[[i]]) ) {
					tmp_mat[j,k] <- x[[i]][j,k]
					tmp_mat[k,j] <- x[[i]][k,j]
				}
			}
			# duplicated node
			if( length( which(duplicated(rownames(tmp_mat))) ) >= 1 ) {
				dup_id <- which(duplicated(rownames(tmp_mat)))
				dup_lab <- unique(rownames(tmp_mat)[dup_id])
				for(j in 1:length(dup_lab)) {
					tmp_id <- which(rownames(tmp_mat) == dup_lab[j])
					for(k in 1:length(tmp_id)) {
						if(k == 1) {
							base_array <- as.numeric(tmp_mat[,tmp_id[1]])
						} else {
							base_array <- base_array + as.numeric(tmp_mat[,tmp_id[k]])
						}
					}
					tmp_mat[tmp_id[1],] <- base_array
					tmp_mat[,tmp_id[1]] <- base_array
					tmp_mat <- tmp_mat[-(tmp_id[-1]),-(tmp_id[-1])]
				}
			}
			x[[i]] <- tmp_mat
		}
	}
	# dupplicated node
	if( length( which( duplicated( paste(vertices.par$x,vertices.par$number,sep='___'))) ) >= 1 ) {
		del_ls <- which( duplicated( paste(vertices.par$x,vertices.par$number,sep='___')))
		vertices.par <- vertices.par[-del_ls,]
	}
	return(list(x,vertices.par,edges.par))
}

# data frame to network format
network.format <- function(x,vertices.par=NULL,edges.par=NULL,direct=F) {
	# 'vertices' is vertices parameter table
	# 'edges' is edges parameter table
	
	x.components <- unlist(lapply(1:length(x) , function(i) rownames(x[[i]]) ) )
	len <- length(x.components)
	x.list <- lapply(1:length(x), function(i) rownames(x[[i]]) )
	
	# clusters.edges to 'igraph'
	igraphs <- lapply(1:length(x), function(i) igraph.form( vertices.par[which(vertices.par$number==i),],edges.par[which(edges.par$number==i),],direct) )
	return(igraphs)
}

# belonging to function 'network.format'
# edge matrix to igraph form, directed graph
igraph.form <- function(vertices.par,edges.par,direct=F) {
	rst <- graph(unlist(lapply(1:nrow(edges.par), function(k) append(as.character(edges.par[k,1]),as.character(edges.par[k,2])) )),directed=direct)
	E(rst)$colour <- edges.par$colour
	E(rst)$width <- edges.par$width
	E(rst)$weight <- 1
	labs <- V(rst)$name
	
	V(rst)$colour <- unlist(lapply(1:length(labs), function(k)  vertices.par$colour[which(vertices.par[,1] == labs[k])] ))
	V(rst)$font_size <- unlist(lapply(1:length(labs), function(k)  vertices.par$font[which(vertices.par[,1] == labs[k])] ))
	V(rst)$shape <- unlist(lapply(1:length(labs), function(k)  vertices.par$shape[which(vertices.par[,1] == labs[k])] ))
	V(rst)$size <- unlist(lapply(1:length(labs), function(k)  vertices.par$size[which(vertices.par[,1] == labs[k])] ))
	V(rst)$lab_col <- unlist(lapply(1:length(labs), function(k)  vertices.par$lab_col[which(vertices.par[,1] == labs[k])] ))
	V(rst)$frame_col <- unlist(lapply(1:length(labs), function(k)  vertices.par$border[which(vertices.par[,1] == labs[k])] ))
	return(rst)
}


### drawing network
# draw single network with full-frame
plot.network <- function(input,network.layout='layout.fruchterman.reingold',arrow.size=0.5) {
	for(z in 1:length(input)) {
		x <- input[[z]]
		plot(x,layout=get(network.layout), vertex.shape=V(x)$shape, vertex.label = V(x)$name, vertex.size=V(x)$size,vertex.color=V(x)$colour,vertex.frame.color=V(x)$frame_col,vertex.label.color=V(x)$lab_col,edge.width=E(x)$width,edge.color=E(x)$colour,vertex.label.cex=V(x)$font_size,rescale=TRUE,edge.arrow.size=arrow.size)
	}
}













