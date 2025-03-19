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
color_file <- args[5] # color.csv

print(args)

# read tree and anno data
tree <- read.newick(tree_file)
strain_data <- read.csv(anno_file, header = FALSE)
num_strain = dim(strain_data)[1]
colnames(strain_data) <- c('label','Clade','Type')
strain_data$Type <- factor(strain_data$Type, levels = c('Reference','Query'))
# read color
col <- read.csv(color_file, header = FALSE)
color_set <- col$V2
names(color_set) <- col$V1
# reroot tree
root_tree <- root(tree, outgroup = out_grp)
# plot tree
p <- ggtree(root_tree, size=2) %<+% strain_data 
node_data = p$data
max_x <- max(node_data$x)
p <- p + 
  geom_tiplab(size=6, color="black",geom = 'label') + 
  geom_tippoint(aes(color = Clade, shape = Type, size = Type)) + 
  geom_treescale(linesize=1.5, color='black', x = max_x*0.8) + 
  theme(legend.position = 'right', legend.background = element_rect(), legend.key = element_blank(), legend.key.size = unit(0.8, 'cm'), legend.text = element_text(size = 15), title = element_text(size = 15)) + 
  xlim(0, max_x*1.8) + 
  scale_color_manual(values = color_set) + scale_shape_manual(values = c(20,17)) + scale_size_manual(values = c(7,9))
# save tree to PNG
if(num_strain > 30)
{
  size_factor = sqrt( num_strain / 30 )
  ggsave(plot = p, filename = fig_name, width = 15*size_factor, height = 15*size_factor)
} else {
  ggsave(plot = p, filename = fig_name, width = 15*1.2, height = 15*1.2)
}

