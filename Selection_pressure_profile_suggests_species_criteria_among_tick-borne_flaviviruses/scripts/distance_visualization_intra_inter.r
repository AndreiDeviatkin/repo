library(ape)
library(otuSummary)
library(phangorn)
library(ggplot2)
library(scales)


####

#Specify path to a folder with files
setwd("D:\\R\\data\\subtypes")


#read FASTA file for nucleic acid sequences
dna <- read.dna( "tick_nt_0.98_al_ORF.fas", format="fasta")
dna_sl1=dna[1:nrow(dna), seq(from = 1, to = length(dna[1,]), by=1)]
#calculate pairwise distance, save the result as matrix
dist_sl1 = dist.gene(dna_sl1, method = "percentage",  pairwise.deletion = TRUE)
#convert  matrix to three column table (id1, id2, nucleotide distance)
nucl_dist1 <- matrixConvert(dist_sl1, colname = c("sp1", "sp2", "nucl_distance")) #convert matrix to a list of values

#read FASTA file for nucleic acid sequences
protein=read.aa("tick_aa_0.98_al_ORF.fas", format = "fasta")
#calculate pairwise distance, save the result as matrix
protein2=as.matrix(dist.hamming(protein))
#convert  matrix to three column table (id1, id2, amino acid distance)
prot_dist1 <- matrixConvert(protein2, colname = c("sp1", "sp2", "prot_distance")) 

#merge two three column tables to one four column table  (id1, id2, nucleotide distance, amino acid distance)  
dist1=cbind( nucl_dist1, prot_dist1$prot_distance)

#change format of the table
final=as.data.frame(cbind(as.character(dist1$sp1),as.character(dist1$sp2),as.numeric(dist1$nucl_distance), as.numeric(prot_dist1$prot_distance)))
colnames(final)=c("id1", "id2", "nucl_distance", "prot_distance") 

#Use regular expression to select a taxon name from FASTA header. This code parse
# "genus_mammalflavi_species_TBEV_subtype_TBEVSib_DQ486861" to "TBEV" and save this
# in a new column entitled "sp1"
final$sp1=gsub('_subtype_.*','', gsub('.*_species_','', final$id1)) 
#Use regular expression to select a taxon name from FASTA header. This code parse
# "genus_mammalflavi_species_GGV_subtype_GGV2_MN830233" to "GGV" and save this
# in a new column entitled "sp2"
final$sp2=gsub('_subtype_.*','', gsub('.*_species_','', final$id2))  

#This regular expression were used to select the title of the subtype
#final$subtype1=gsub('_.*','', gsub('.*_subtype_','', final$id1)) 
#final$subtype2=gsub('_.*','', gsub('.*_subtype_','', final$id2)) 
#intersubtype_intraspecies=as.data.frame(subset (intraspecies, intraspecies$subtype1!=intraspecies$subtype2))
#intrasubtype_intraspecies=as.data.frame(subset (intraspecies, intraspecies$subtype1==intraspecies$subtype2))


#Select rows from the table with identical values in columns "sp1" and "sp2"
intraspecies = as.data.frame(subset (final, final$sp1==final$sp2))
#Select rows from the table with not identical values in columns "sp1" and "sp2"
interspecies = as.data.frame(subset (final, final$sp1!=final$sp2))


#generate a plot demonstrating the distribution of pairwise distances in amino acid and nucleotide sequences
dna_protein_inter_ORF=ggplot()+
  geom_point(data=interspecies, aes(as.numeric(interspecies[,3])*100,as.numeric(interspecies[,4])*100),
             colour = "darkred", size=7, alpha=0.03, shape=4)+
  geom_smooth(method = "lm", se = FALSE, data=interspecies, aes(as.numeric(interspecies[,3])*100,as.numeric(interspecies[,4])*100),
              colour = "red",  linewidth=1.5, linetype = 1)+
  geom_point(data=intraspecies, aes(as.numeric(intraspecies[,3])*100,as.numeric(intraspecies[,4])*100),
             colour = "darkblue", size=7, alpha=0.03, shape=4)+
  geom_smooth(method = "lm", se = FALSE, data=intraspecies, aes(as.numeric(intraspecies[,3])*100,as.numeric(intraspecies[,4])*100),
              colour = "blue",  linewidth=1.5, linetype = 1)+
  xlab('nucleotide distance, %')+ylab('protein distance, %')+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 53)) + 
  scale_y_continuous(expand = c(0, 0),limits = c(0, 53), breaks= pretty_breaks())+
  theme(
    axis.line = element_line(colour = 'black', linewidth = 2),
    axis.title.x = element_text( size = 30, angle = 0, hjust = 0.5, vjust = -0.5, face = "plain"),
    axis.title.y = element_text( size = 30, angle = 90, hjust = 0.5, vjust = 0.5, face = "plain"),
    axis.text.x = element_text(color = "black", size = 25,  hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(color = "black", size = 25,  hjust = .5, vjust = .5, face = "plain"),
    axis.ticks.length=unit(.25, "cm"),
    plot.background=element_rect(fill = "white"),
    panel.background = element_rect(fill = "white")
    
  )

#visualize a plot
dna_protein_inter_ORF
