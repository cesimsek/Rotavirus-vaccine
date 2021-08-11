library(ggtree)
library(treeio)
library(Biostrings)
library(ape)
library(phytools)


setwd("~/Documents/PhD/RVAB/RVA_study/RV_Phylogenies/vp7")

treefile <- "vp7_gtree.nwk"
tre <- read.tree(treefile)


#get the tip number
#MRCA(tre, tip=c('RVA/Human-wt/BEL/F09825/2017/G3P8'))


######### Make a genotype heatmap next to the tree

genotype <- read.table('nodes_genotypes_vp7.txt',sep="\t", stringsAsFactor=F )


cols <- c('#E64B35FF',"#59260B" , '#67D89A', '#228B22', '#FF9900',  
          '#1E5687', '#9CDEF6', '#d7bdec', '#DC52BF', '#0DA3BA', '#6549DA', 
          "#ab0903",'#ffbaba', '#086972','#00A087FF')


# Visualize the tree + genotypes, arrange the colors beforehand
p <- ggtree(tre) + geom_treescale(y=-2) + geom_tiplab(align = T, size = 4)
p <- rotate(p, 127)

gheatmap(p, genotype, offset = 0.4, width=0.4, colnames_position = "top", font.size = 3) +
scale_fill_manual(values=cols) 



#check the tip labels
ggtree(tre) + geom_text(aes(label=node), hjust=-.3)






# you can also read raxml with bootstrap values - check the nodes with bootstrap values above 70

rax <- read.tree(file = 'vp7_new.raxml.support')
ggtree(rax) + geom_nodelab(aes(label=label, subset=as.numeric(label) > 70))

#reroot tree to G6
rax_rooted <- rax
rax_rooted@phylo <- root(rax_rooted@phylo, outgroup = "RVA/Roe-deer-wt/SLO/D38-14/2014/G6P15")
rax_rooted <- reroot(as.phylo(rax_rooted), outgroup)


p1 <- ggtree(rax) + geom_treescale(y=-2) + geom_nodelab(aes(label=label,subset=as.numeric(label) > 70), size=3) +
  geom_tiplab(align = T, size = 4)

gheatmap(p1, genotype, offset = 0.4, width=0.4, colnames_position = "top", font.size = 3) +
  scale_fill_manual(values=cols) 




