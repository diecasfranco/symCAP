library(ggplot2)
library(reshape2)
library(forcats)
library(tidyverse)
library(ggnewscale)



genes <- read.table("file_for_genes_plot.txt", header = TRUE, sep = "\t")
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t")
#colnames(metadata)
gene.molten <- melt(genes, measure.vars = 2:47)
gene.metadata <- merge(gene.molten, metadata)

colors <- c("Issidae" = "#6eb66c23", "Achilidae" = "#a379272c", "Delphacidae" = "#db19281a", "Caliscelidae" = "#ff99373b", 
            "Dictyopharidae" = "#006eff20", "Cixiidae" = "#00808022", "Meenoplidae" = "#f0072c44", "Derbidae" = "#ff80e720", 
            "Fulgoridae" ="#a379272c", "Tropiduchidae"="#2b00002c", "Lophopidae"="#c8b7b745", "Tettigometridae"="#ff2a581d", "white"="white")

sample_order <- c("VFASICLA","VFUGYOPS","VFSTEMIN","VFKELRIB","VFKELVIT","VFHYALUT","VFPENROR","VFOLIH","VFCIXIUS","VFCIXNER","VFMEEALB","VFTRYOCC","VFAKOQUE","VFCIXPIL","VFMALBOS","VFDER3","VFPYRCLA","VFPYRLAN","VFPYRVIR","VFPENVAR","VFRANSCY","VFPARPLA","VFRANEDI","VFCALKRU","VFDICMUL","VFDICEUR","VFDICPAN5","VFRR05105","VFTETSUL","VFTETHEX","VFTETIMP","VFTETMAC","VFPARIOC","VFOMMLON","VFOMMDIS","VFCALBON1","VFISSCOL1","VFISSKAZ","VFLIBAN","VFSCODIS","VFAGAHAV","VFZOPTEN","VFMYCCUN")

gene.metadata$samples <- factor(gene.molten$samples, levels = sample_order)

ggplot(data=gene.metadata, aes(x = variable, y = fct_rev(samples), fill= factor (value), color= factor(value))) + 
  geom_point(size = 6, shape = 21) + scale_fill_manual(values = c("white","darkgrey","black"), guide = FALSE) + 
  scale_color_manual(values = c("black","black","black"), guide = FALSE) +
  scale_fill_manual(values = colors, guide = FALSE, na.value = "white") +
  scale_x_discrete(position = "top") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,size=10)) + theme(axis.text.y=element_text(size=14)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "", y = "", fill = "")
  



ggplot(data=gene.metadata, aes(x = variable, y = fct_rev(samples), fill=factor (value))) +
  geom_tile(color = "white", size = 0.5) + geom_point(aes(color= factor(value),size = 6, shape =factor(value))) +
  scale_shape_manual(values = c("white","darkgrey","black")) +
  theme_bw() +
  scale_x_discrete(position = "top") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "", y = "", fill = "")


gene.metadata$family <- factor(gene.metadata$family, levels = names(colors))

ggplot(data=gene.metadata, aes(x = variable, y = fct_rev(samples))) + 
  geom_point(aes(fill = family), size = 7, shape = 21) + scale_fill_manual(values = c(colors)) 
+ scale_x_discrete(position = "top") + theme(axis.text.x=element_text(angle=90,hjust=1,size=16)) 
+ theme(axis.title.x = element_text(size=20))



ggplot(data = gene.molten, aes(x = variable, y = fct_rev(samples))) +
  geom_point(aes(fill = factor(ifelse(value == 0.5, "GRAY", ifelse(value > 0, "TRUE", "FALSE")))), size = 7, shape = 21) +
  scale_fill_manual(values = c("white", "black", "gray"), drop = FALSE) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16)) +
  theme(axis.title.x = element_text(size = 20))

########

library(ggplot2)
library(reshape2)

genes <- read.table("file_for_genes_plot.txt", header = TRUE, sep = "\t")
gene.molten <- melt(genes, measure.vars = 2:47)

ggplot(data = gene.molten, aes(x = variable, y = fct_rev(samples))) +
  geom_point(aes(shape = ifelse(value == 0.5, "half", "full"), 
                 fill = factor(ifelse(value > 0, "TRUE", "FALSE"))), 
             size = 7) +
  scale_shape_manual(values = c(full = 21, half = 22)) +
  scale_fill_manual(values = c("white", "black"), drop = FALSE) +
  scale_x_discrete(position = "bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16)) +
  theme(axis.title.x = element_text(size = 20))




# set the order of the samples
sample_order <- c("VFASICLA","VFUGYOPS","VFSTEMIN","VFKELRIB","VFKELVIT","VFHYALUT","VFPENROR","VFOLIH","VFCIXIUS","VFCIXNER","VFMEEALB","VFTRYOCC","VFAKOQUE","VFCIXPIL","VFMALBOS","VFDER3","VFPYRCLA","VFPYRLAN","VFPYRVIR","VFPENVAR","VFRANSCY","VFPARPLA","VFRANEDI","VFCALKRU","VFDICMUL","VFDICEUR","VFDICPAN5","VFRR05105","VFTETSUL","VFTETHEX","VFTETIMP","VFTETMAC","VFPARIOC","VFOMMLON","VFOMMDIS","VFCALBON1","VFISSCOL1","VFISSKAZ","VFLIBAN","VFSCODIS","VFAGAHAV","VFZOPTEN","VFMYCCUN")
gene.metadata$samples <- factor(gene.metadata$samples, levels = sample_order)

# create a new column for fill color based on the value
gene.metadata$fill_color <- ifelse(gene.metadata$value == "TRUE", gene.metadata$family, 
                                   ifelse(gene.metadata$value == "FALSE", "white",
                                          ifelse(gene.metadata$value == "HALF", paste0(gene.metadata$family, "80"), 
                                                 "gray")))



# plot the data
ggplot(data = gene.metadata, aes(x = variable, y = fct_rev(samples))) +
  geom_point(size = 6, shape = 21, aes(fill = fill_color), color = "black") +
  scale_fill_manual(values =colors, guide = FALSE)

ggplot(data = gene.metadata, aes(x = variable, y = fct_rev(samples))) +
  geom_point(size = 6, shape = 21, aes(fill = fill_color), color = "black") +
  scale_fill_manual(values = c("FALSE" = "white", "HALF" = "gray",colors), guide = FALSE) +
  scale_x_discrete(position = "top") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,size=10)) + theme(axis.text.y=element_text(size=14)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),panel.background = element_blank()) +
  labs(x = "", y = "", fill = "")




### Plot diamonds if HALF
ggplot(data = gene.metadata, aes(x = variable, y = fct_rev(samples))) +
  geom_point(size = 6, shape = 21, aes(fill = fill_color), color = "black", 
             data = subset(gene.metadata, value %in% c(TRUE, FALSE))) +
  geom_point(size = 4, shape = 23, aes(fill = fill_color, color = "black"), 
             data = subset(gene.metadata, value == "HALF")) +
  scale_fill_manual(values = colors, guide = FALSE)+
  scale_x_discrete(position = "top") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,size=10)) + theme(axis.text.y=element_text(size=13)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),panel.background = element_blank()) +
  labs(x = "", y = "", fill = "") + theme(legend.position="none")
