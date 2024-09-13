library(igraph)
library(dplyr)
library(qgraph)

### Credit: Shushrruth Sai Srinivasan ### 

# Figure 5d: EXN TF Network
path_to_tf <- "/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/TF_Peak_Gene_networks_for_Sai_2M/"

tf_exh = read.csv(paste(path_to_tf,'EXC.csv',sep=''))

exc_genes <- c('SMARCC1','EGR2','TAL2')
tf_exh <- tf_exh[tf_exh$TF %in% c('SMARCC1','EGR2','TAL2'),]
x <- tf_exh[tf_exh$Correlation > 0.65,]
y <- tf_exh[tf_exh$Correlation > 0.7,]

ptsd_gwas <- unique(y[y$PTSD_GWAS=='TRUE',]$geneName)
up_genes <- unique(y[y$DEG=='UP',]$geneName)
down_genes <- unique(y[y$DEG=='DOWN',]$geneName)

df <- x[c('TF','geneName')]

df_exc <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df_exc) <- c('TF', 'geneName')

for (tf in unique(df$TF)){
    df_exc <- rbind(df_exc,unique(df[df$TF==tf,]))
}

gr <- graph.data.frame(d = df_exc, directed = FALSE)

V(gr)$size[V(gr)$name %in% df$geneName] = 1
V(gr)$color[V(gr)$name %in% df$geneName] = 'white'
V(gr)$frame.color[V(gr)$name %in% df$geneName] = 'gray68'
                           
V(gr)$size[V(gr)$name %in% up_genes] = 4
V(gr)$color[V(gr)$name %in% up_genes] = "#EF3B2C"

V(gr)$size[V(gr)$name %in% down_genes] = 4
V(gr)$color[V(gr)$name %in% down_genes] = "#4292C6"

V(gr)$size[V(gr)$name %in% c('TAL2','SMARCC1')] = 15
V(gr)$size[V(gr)$name %in% 'EGR2'] = 20
V(gr)$color[V(gr)$name %in% exc_genes] = '#b22222'

V(gr)$size[V(gr)$name %in% ptsd_gwas] = 4
V(gr)$color[V(gr)$name %in% ptsd_gwas] = "gold"

for (deg in unique(up_genes)){
    V(gr)$label[V(gr)$name == deg] = deg
    V(gr)$label.cex[V(gr)$name == deg] = 0.6
    V(gr)$label.color[V(gr)$name == deg] = 'black'
    V(gr)$frame.color[V(gr)$name == deg] = 'gray28'
}

for (deg in unique(down_genes)){
    V(gr)$label[V(gr)$name == deg] = deg
    V(gr)$label.cex[V(gr)$name == deg] = 0.6
    V(gr)$label.color[V(gr)$name == deg] = 'black'
    V(gr)$frame.color[V(gr)$name == deg] = 'gray28'
}

for (tf_gene in unique(c(exc_genes))){
    V(gr)$label[V(gr)$name == tf_gene] = tf_gene
    V(gr)$label.cex[V(gr)$name == tf_gene] = 2
    V(gr)$label.color[V(gr)$name == tf_gene] = 'black'
    V(gr)$frame.color[V(gr)$name == tf_gene] = 'gray28'
}

for (tf_gene in unique(c(ptsd_gwas))){
    V(gr)$label[V(gr)$name == tf_gene] = tf_gene
    V(gr)$label.cex[V(gr)$name == tf_gene] = 0.6
    V(gr)$label.color[V(gr)$name == tf_gene] = 'black'
    V(gr)$frame.color[V(gr)$name == tf_gene] = 'gray28'
}

E(gr)$color = 'gray78' #   mediumpurple1   pink
E(gr)$width = 0.0001
E(gr)$arrow.size = 0.0001

e <- get.edgelist(gr,names=FALSE)
layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(gr))
plot(gr, layout = layout, vertex.label.family='Helvetica', vertex.frame.width=0.5)

layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(gr))
pdf('/home/ah2428/EXC_TF.pdf', width=10, height=10)
plot(gr, layout = layout, vertex.label.family='Helvetica', vertex.frame.width=0.5,asp=0)
dev.off()


# Figure 5g: IN TF Network
path_to_tf <- "/gpfs/gibbs/pi/gerstein/jz435/ShareZhangLab/PTSD/ATAC/data/TF_Peak_Gene_networks_for_Sai_2M/"

tf_inh = read.csv(paste(path_to_tf,'INH.csv',sep=''))
x <- tf_inh[tf_inh$TF=='TFAP4',]
deg <- unique(x[x$DEG %in% c('UP','DOWN'),]$geneName)

inh_genes <- c('WT1','TFAP4','ZNF238')
tf_inh <- tf_inh[tf_inh$TF %in% c('WT1','TFAP4','ZNF238'),]
x <- tf_inh[tf_inh$Correlation > 0.75,]

ptsd_gwas <- unique(x[x$PTSD_GWAS=='TRUE',]$geneName)
up_genes <- unique(x[x$DEG=='UP',]$geneName)
down_genes <- unique(x[x$DEG=='DOWN',]$geneName)

df <- x[c('TF','geneName')]

df_inh <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df_inh) <- c('TF', 'geneName')

for (tf in unique(df$TF)){
    df_inh <- rbind(df_inh,unique(df[df$TF==tf,]))
}

gr <- graph.data.frame(d = df_inh, directed = FALSE)

V(gr)$size[V(gr)$name %in% df$geneName] = 1
V(gr)$color	[V(gr)$name %in% df$geneName] = 'white'
V(gr)$frame.color[V(gr)$name %in% df$geneName] = 'gray68'

V(gr)$size[V(gr)$name %in% up_genes] = 4
V(gr)$color[V(gr)$name %in% up_genes] = "#EF3B2C"

V(gr)$size[V(gr)$name %in% down_genes] = 4
V(gr)$color[V(gr)$name %in% down_genes] = "#4292C6"

V(gr)$size[V(gr)$name %in% c('WT1','ZNF238')] = 15
V(gr)$size[V(gr)$name %in% 'TFAP4'] = 20
V(gr)$color[V(gr)$name %in% inh_genes] = '#2E8B57'

V(gr)$size[V(gr)$name %in% ptsd_gwas] = 4
V(gr)$color[V(gr)$name %in% ptsd_gwas] = "gold"

V(gr)$size[V(gr)$name %in% ptsd_gwas] = 4
V(gr)$color[V(gr)$name %in% ptsd_gwas] = "gold"

for (deg in unique(up_genes)){
    V(gr)$label[V(gr)$name == deg] = deg
    V(gr)$label.cex[V(gr)$name == deg] = 0.6
    V(gr)$label.color[V(gr)$name == deg] = 'black'
    V(gr)$frame.color[V(gr)$name == deg] = 'gray28'
}

for (deg in unique(down_genes)){
    V(gr)$label[V(gr)$name == deg] = deg
    V(gr)$label.cex[V(gr)$name == deg] = 0.6
    V(gr)$label.color[V(gr)$name == deg] = 'black'
    V(gr)$frame.color[V(gr)$name == deg] = 'gray28'
}

for (tf_gene in unique(c(inh_genes))){
    V(gr)$label[V(gr)$name == tf_gene] = tf_gene
    V(gr)$label.cex[V(gr)$name == tf_gene] = 2
    V(gr)$label.color[V(gr)$name == tf_gene] = 'black'
    V(gr)$frame.color[V(gr)$name == tf_gene] = 'gray28'
    
}

for (tf_gene in unique(c(ptsd_gwas))){
    V(gr)$label[V(gr)$name == tf_gene] = tf_gene
    V(gr)$label.cex[V(gr)$name == tf_gene] = 0.6
    V(gr)$label.color[V(gr)$name == tf_gene] = 'black'
    V(gr)$frame.color[V(gr)$name == tf_gene] = 'gray28'
    
}

E(gr)$color = 'gray78' #   mediumpurple1   pink
E(gr)$width = 0.0001
E(gr)$arrow.size = 0.0001

e <- get.edgelist(gr,names=FALSE)
layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(gr))
plot(gr, layout = layout, vertex.label.family='Helvetica', vertex.frame.width=0.5)

layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(gr))
pdf('/home/ah2428/INH_TF_2.pdf', width=10, height=10)
plot(gr, layout = layout, vertex.label.family='Helvetica', vertex.frame.width=0.5,asp=0)
dev.off()
