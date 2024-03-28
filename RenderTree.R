# load libraries
library(ggtree)
library(treeio)
library(ggplot2)

# read parameters
args <- commandArgs(trailingOnly = TRUE)

tree_file = args[1] # './RSV_A.nwk'
out_grp = args[2] #'RSVA/Human/USA/78G-104-01/1978|KU316155'
anno_file <- args[3] #'RSV_A.csv'
fig_name <- args[4] #'RSV_A.png'

# read tree and anno data
tree <- read.newick(tree_file)
strain_data <- read.csv(anno_file, header = FALSE)
colnames(strain_data) <- c('label','Type')
# reroot tree
root_tree <- root(tree, outgroup = out_grp)
# plot tree
p <- ggtree(root_tree, size=2) %<+% strain_data 
node_data = p$data
max_x <- max(node_data$x)
p <- p + 
  geom_tiplab(size=6, color="black",geom = 'label') + 
  geom_tippoint(aes(color = Type), size = 4) + 
  geom_treescale(linesize=1.5, color='black', x = max_x*0.8) + 
  theme(legend.position = 'bottom', legend.background = element_rect(), legend.key = element_blank(), legend.key.size = unit(0.8, 'cm'), legend.text = element_text(size = 15), title = element_text(size = 15)) + 
  xlim(0, max_x*1.2)
# save tree to PNG
ggsave(plot = p, filename = fig_name, width = 13, height = 15)
