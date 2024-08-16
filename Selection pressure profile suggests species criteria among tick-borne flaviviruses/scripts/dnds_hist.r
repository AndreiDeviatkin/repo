library(ggplot2)
library (scales)

#Specify path to a folder with files
setwd("D:\\R\\data\\subtypes")

#read output of computeNeiGojobori.py
df=read.csv("tick_nt_0.98_al_ORF_codon.fas_dnds.txt")
#calculate raw dn/ds ratio
df$ratio_dnds=(df[,8]/df[,7])
#format table, generate a three column table 
dnds=as.data.frame(cbind(as.character(df$name1), as.character(df$name2), as.numeric(df$ratio_dnds)))
colnames(dnds)=c("id1", "id2", "dnds_ratio")

#Use regular expression to select a taxon name from FASTA header. This code parse
# "genus_mammalflavi_species_GGV_subtype_GGV2_MN830233" to "GGV" and save this
# in a new column entitled "id1_species"
dnds$id1_species=gsub('_subtype_.*','', gsub('.*_species_','', dnds$id1))
#Use regular expression to select a taxon name from FASTA header. This code parse
# "genus_mammalflavi_species_TBEV_subtype_TBEVSib_DQ486861" to "TBEV" and save this
# in a new column entitled "id2_species"
dnds$id2_species=gsub('_subtype_.*','', gsub('.*_species_','', dnds$id2)) 

#Select rows from the table with identical values in columns "id1_species" and "id2_species"
intraspecies = as.data.frame(subset (dnds, dnds$id1_species==dnds$id2_species))
#Select rows from the table with not identical values in columns "id1_species" and "id2_species"
interspecies = as.data.frame(subset (dnds, dnds$id1_species!=dnds$id2_species))

intra2=data.frame(as.numeric(intraspecies$dnds_ratio))
intra2$status <- 'intraspecies'
inter2=data.frame(as.numeric(interspecies$dnds_ratio))
inter2$status <- 'interspecies'
colnames(inter2)=c("ratio", "status")
colnames(intra2)=c("ratio", "status")
combo <- rbind(inter2, intra2)


#generate a plot demonstrating the distribution of pairwise dn/ds distances
ggplot(combo, aes(ratio, fill = status)) + geom_histogram(alpha = 0.8, binwidth = 0.01)+  
  scale_x_continuous(expand = c(0, 0),limits = c(0, max(combo$ratio)*1.05), breaks= pretty_breaks())+
  scale_y_continuous(expand = c(0, 0), 
  breaks= pretty_breaks())+
  theme(
    axis.line = element_line(colour = 'black', linewidth = 1),
    axis.title.x = element_text( size = 35, angle = 0, hjust = .5, vjust = -0.5, face = "plain"),
    axis.title.y = element_text( size = 35, angle = 90, hjust = .5, vjust = 0.5, face = "plain"),
    axis.text.x = element_text(color = "black", size = 30,  hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(color = "black", size = 30,  hjust = .5, vjust = .5, face = "plain"),
    axis.ticks.length=unit(.25, "cm"),
    plot.background=element_rect(fill = "white"),
    panel.background = element_rect(fill = "white")
    )+  geom_vline(xintercept = 0.37, linetype="dashed",
             color = "black", linewidth=1.5)+
  theme(legend.position = "none")
