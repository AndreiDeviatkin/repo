library(ggplot2)
library (scales)
library(igraph)
library(treeio)
library(ape)
library(phytools)

#Specify path to a folder with files
setwd("D:\\R\\data\\subtypes")

#read output of computeNeiGojobori.py
df=read.csv("tick_nt_0.98_al_ORF_codon.fas_dnds.txt")
#calculate raw dn/ds ratio
df$ratio_dnds=(df[,4]/df[,3])
#format table, generate a three column table 
dnds=as.data.frame(cbind(as.character(df$name1), as.character(df$name2), as.numeric(df$ratio_dnds)))
colnames(dnds)=c("id1", "id2", "dnds_ratio")

#generate matrix for pairwise dn/ds ratios
dnds_for_tree=data.frame(cbind(dnds$id1, dnds$id2, dnds$dnds_ratio))
colnames(dnds_for_tree)=c("id1", "id2", "dnds_ratio")
dnds_for_tree$dnds_ratio=as.numeric(dnds_for_tree$dnds_ratio)
dnds_for_tree <- graph.data.frame(dnds_for_tree, directed=FALSE)
dnds_for_tree_tmp=get.adjacency(dnds_for_tree, attr="dnds_ratio", sparse=FALSE)
tmp2=as.numeric(dnds_for_tree_tmp)
dnds_for_tree=matrix(tmp2, ncol=nrow(dnds_for_tree_tmp), nrow=nrow(dnds_for_tree_tmp))
colnames(dnds_for_tree)=colnames(dnds_for_tree_tmp)
rownames(dnds_for_tree)=rownames(dnds_for_tree_tmp)
final= as.matrix(dnds_for_tree)

#construct phylogenetic tree for pairwise dn/ds ratios

tree <- nj(final)
tree <- ladderize(tree)
tree=midpoint.root(tree)

#save  phylogenetic tree in Newick format
write.tree(tree, file = "dnds.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)