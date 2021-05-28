library(ggtree)
library(treeio)
library(Biostrings)
library(ape)
library(phytools)


setwd("~/Documents/PhD/RVAB/RVA_study/RV_Phylogenies/vp7")

treefile <- "vp7_gtree.nwk"
tre <- read.tree(treefile)


#get the tip number - didnt work
#MRCA(tre, tip=c('RVA/Human-wt/BEL/F09825/2017/G3P8'))

######## OK: VISUALIZE A TREE
ggtree(tre) + geom_treescale() + geom_tiplab(size = 2)



#########  OK: Make a genotype heatmap next to the tree

genotype <- read.table('nodes_genotypes_vp7.txt',sep="\t", stringsAsFactor=F )


cols <- c('#E64B35FF',"#59260B" , '#67D89A', '#228B22', '#FF9900',  
          '#1E5687', '#9CDEF6', '#d7bdec', '#DC52BF', '#0DA3BA', '#6549DA', 
          "#ab0903",'#ffbaba', '#086972','#00A087FF')


# SEE THE TREE WITH GENOTYPES, ARRANGE THE COLORS!
p <- ggtree(tre) + geom_treescale(y=-2) + geom_tiplab(align = T, size = 4)
p <- rotate(p, 127)

gheatmap(p, genotype, offset = 0.4, width=0.4, colnames_position = "top", font.size = 3) +
scale_fill_manual(values=cols) 

#save at 17*20?


#check the tip labels
ggtree(tre) + geom_text(aes(label=node), hjust=-.3)





###### not tried yet but can be done later!

#treeio
read.codeml()



#read raxml with bootstrap values - only get the bootstrap values above 70

rax <- read.tree(file = 'vp7_new_japan.raxml.support')
ggtree(rax) + geom_nodelab(aes(label=label, subset=as.numeric(label) > 70))

#reroot tree to G6
rax_rooted <- rax
rax_rooted@phylo <- root(rax_rooted@phylo, outgroup = "RVA/Roe-deer-wt/SLO/D38-14/2014/G6P15")
rax_rooted <- reroot(as.phylo(rax_rooted), outgroup)


outgroup <- findMRCA(as.phylo(rax_rooted), c("RVA/Roe-deer-wt/SLO/D38-14/2014/G6P15"))


p1 <- ggtree(rax) + geom_treescale(y=-2) + geom_nodelab(aes(label=label,subset=as.numeric(label) > 70), size=3) +
  geom_tiplab(align = T, size = 4)

gheatmap(p1, genotype, offset = 0.4, width=0.4, colnames_position = "top", font.size = 3) +
  scale_fill_manual(values=cols) 

#ggtree(rax) + geom_text(aes(label=node), hjust=-.3)


