# ml conda_R/4.3
# R
###########
# library(Seurat)
# library(SummarizedExperiment)
# library(SpatialExperiment)
# library(Banksy)
library(scuttle)
# library(scater)
library(cowplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

###########
dir_out = "/dcs04/hansen/data/ywang/xenium_joint_runs/out/"
setwd(dir_out)

###########
list_pval_allGenes_perClu = readRDS("list_pval_allGenes_perClu.rds")
list_fc_allGenes_perClu = readRDS("list_fc_allGenes_perClu.rds")
list_DEG_perClu = readRDS("list_DEG_perClu.rds")
# genes_deg_tmp = c("LYPD6B", "CDH13","MC3R","CDH4","RXRG")
genes_deg_tmp = c("LYPD6B", "CDH13","CDH4")


###########
clus_tmp = "33"

df = data.frame(logfc = list_fc_allGenes_perClu[[clus_tmp]],
                pval = list_pval_allGenes_perClu[[clus_tmp]],
                gene = names(list_fc_allGenes_perClu[[clus_tmp]]))
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "not significant"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$logfc > 0.6 & df$pval < 0.05] <- "higher in male"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$logfc < -0.6 & df$pval < 0.05] <- "higher in female"
# head(df[order(df$padj) & df$diffexpressed == 'DOWN', ])
df$diffexpressed = factor(df$diffexpressed, levels = c("higher in male","not significant","higher in female"))
df$delabel <- ifelse(df$gene %in% genes_deg_tmp, df$gene, NA)

pdf(paste0("representative_plots/vocano_clu_",clus_tmp,".pdf"),width=5,height=5)
myvolcanoplot <- ggplot(data = df, aes(x = logfc, y = -log10(pval), col = diffexpressed, label = delabel))  +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed')  +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed')  +
  geom_point(size = 1) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable<br />
                     labels = c("higher in male", "not significant", "higher in female")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)<br />
  # coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits<br />
  # coord_cartesian(ylim = c(0, 4), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits<br />
  labs(color = 'Differential expression', #legend_title,<br />
       x = expression("log fold change"), y = expression("-log p-value"))  +
  # scale_x_continuous(breaks = seq(-10, 10, 2)) + 
  ggtitle(paste0("cell cluster ", clus_tmp))  + # to customise the breaks in the x axis<br />
  # ggtitle('Thf-like cells in severe COVID vs healthy patients') + # Plot title<br />
  geom_text_repel(max.overlaps = Inf, min.segment.length=0, seed = 22, box.padding = 0.5,show.legend  = FALSE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))#+
  # geom_text(show.legend  = FALSE)

print(myvolcanoplot)
dev.off()
