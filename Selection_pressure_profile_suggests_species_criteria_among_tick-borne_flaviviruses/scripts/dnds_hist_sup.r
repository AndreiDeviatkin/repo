library(ggplot2)
library (scales)
library(ggpubr)
library(dplyr)


#Specify path to a folder with Yang-Nielsen dN/dS tables
setwd("../data/")

# list of genome regions
proteins <- c("ORF", "E", "NS3")
# list of vectors
vectors = c('tick','mosquito')


for (vector in vectors){
    plots_list = list()
    for (x in proteins) {
      if (vector=='tick'){
        title_vector = "TBFV"
      }
      if (vector=='mosquito'){
        title_vector = "MBFV"
      }

      # read table with pairwise Yang-Nielsen dN/dS values. These tables are the result of parse_yn.py script 
      # which processes the output of yn0.exe (PAML)
      fname = paste0(vector,"_nt_0.98_al_", x, "_codonYN_ed.txt")
      df=read.csv(fname)
    
      coldS = 10
      coldN = 8
      df[,coldN] = as.numeric(df[,coldN])
      df[,coldS] = as.numeric(df[,coldS])
      df$ratio_dnds=(df[,coldN]/df[,coldS])

      
      dnds=as.data.frame(cbind(as.character(df$name1), as.character(df$name2), as.numeric(df$ratio_dnds)))
      colnames(dnds)=c("id1", "id2", "dnds_ratio")
      dnds$dnds_ratio = as.numeric(df$ratio_dnds)

      #Use regular expression to select a taxon name from FASTA header. This code parse
      # "genus_mammalflavi_species_GGV_subtype_GGV2_MN830233" to "GGV" and save this
      # in a new column entitled "id1_species"
      dnds$id1_species=gsub('_subtype_.*','', gsub('.*_species_','', dnds$id1)) 
      dnds$id2_species=gsub('_subtype_.*','', gsub('.*_species_','', dnds$id2)) 

      # add new column that indicates whether comparison is inter or intraspecies
      dnds = dnds %>% dplyr::mutate(status = ifelse(id1_species == id2_species, "intraspecies", "interspecies"))

      #generate a histogram for the distribution of pairwise Yang-Nielsen dN/dS
      hist_intersp = ggplot(dnds, aes(dnds_ratio, fill = status)) + geom_histogram(alpha = 0.8, binwidth = 0.01)+  
        scale_x_continuous(expand = c(0, 0), breaks= pretty_breaks())+
        scale_y_continuous(expand = c(0, 0), 
                           breaks= pretty_breaks())+
        ggtitle(x)+
        xlab("Yang-Nielsen dN/dS ratio") +
        theme(
          plot.title = element_text(size = 30),
          axis.line = element_line(colour = 'black', linewidth = 1),
          axis.title.x = element_text( size = 20, angle = 0, hjust = .5, vjust = -0.5, face = "plain"),
          axis.title.y = element_text( size = 20, angle = 90, hjust = .5, vjust = 0.5, face = "plain"),
          axis.text.x = element_text(color = "black", size = 20,  hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "black", size = 20,  hjust = .5, vjust = .5, face = "plain"),
          axis.ticks.length=unit(.25, "cm"),
          plot.background=element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")
        )+ 
        theme(legend.position = "none")
      
      plots_list = append(plots_list, list(hist_intersp))
 
    }
  
  plot_all = ggarrange(plotlist=plots_list, ncol = 3, nrow=1)
  ggsave(filename=paste0(vector,"_YN.jpg"),plot_all ,width=15, height=5, dpi=600)
}
