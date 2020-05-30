##### Script for analysing bulk ATAC-seq of reprogramming intermediates
# Note: all files can be downloaded from http://hrpi.ddnetbio.com/
# Note: all files are assumed to be in the data folder

### Clear workspace and load libraries
rm(list=ls())
library(dplyr)
library(edgeR)
library(irlba)
library(sjstats)
library(Mfuzz)
# library(limma)
# library(dplyr)
# library(made4)
# library(ggplot2)
# library(ggrepel)
# library(stringr)
# library(reshape2)
# library(plyr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Mfuzz")

# source("https://bioconductor.org/biocLite.R")
# biocLite("VariantAnnotation")
# biocLite("systemPipeR")

# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Glimma")
# library(Glimma)
# library(ggfortify)
# library(GenomicRanges)
# source("http://bioconductor.org/biocLite.R")
# biocLite("annotatr")
# library(annotatr)
# install.packages("xlsx")
# library("xlsx")

# Define colour palettes

colorMedia = c("black","darkorange","blue","darkolivegreen2","red2","pink2")
names(colorMedia) = c("D0, 3, 7","Primed","t2iLGoY","5iLAF","NHSM","RSeT")
shapeTime = c(0,1,2,15,16,17,18)
names(shapeTime) = c("D0","D3","D7","D13","D21","P3","P10")

# Exploratory Analysis - PCA
## Import read counts of consensus peak set.

peak.counts <-
  read.table(file = "data/consensus_peak_set_counts.txt",
             header = T,
             stringsAsFactors = F)

### Reorder samples 
order <- colnames(peak.counts)[7:ncol(peak.counts)]
order <- order[c(1, 2, 11:14, 3, 5, 7, 9, 19, 21, 15, 17, 4, 6, 8, 10, 20, 22, 16, 18)]
order <- factor(order, levels = order)

### Peak information
peak.info <- colnames(peak.counts)[1:6]

### Peak counts table reordered

peak.counts <- peak.counts[, c(peak.info, levels(order))]

### Change Geneid column name to PeakID

colnames(peak.counts)[1] <- "PeakID"
rownames(peak.counts) <- peak.counts$PeakID

## Within and between library normalization.

### Quantile normalised log2(RPKM + 1)

peaks.fpkm <-
  rpkm(peak.counts[, 7:ncol(peak.counts)],
       log = F, gene.length = peak.counts$Length)

###Keep peaks with more than 5 FPKMs in at least 2 samples for exploratory analysis.

keep <- rowSums(peaks.fpkm > 5) >= 2

peaks.fpkm <- peaks.fpkm[keep, ]
peaks.fpkm <- log2(peaks.fpkm + 1)
peaks.fpkm <- normalizeQuantiles(peaks.fpkm)

### PCA 
pca.plot.annotation <- data.table(sample_id = colnames(peaks.fpkm))

pca.plot.annotation$sample_labels <- c(
  paste(rep("FM_", 3), rep("D", 3), rep(c(0, 3, 7), each = 2), rep(c("_R1", "_R2")), sep = ""),
  paste(rep("PM_", 2), rep("D", 2), rep(c(13, 21), each = 2), rep(c("_R1", "_R2")), sep = ""),
  paste(rep("PM_", 2), rep("P", 2), rep(c(3, 10), each = 2), rep(c("_R1", "_R2")), sep = ""),
  paste(rep("NM_", 2), rep("D", 2), rep(c(13, 21), each = 2), rep(c("_R1", "_R2")), sep = ""),
  paste(rep("NM_", 2), rep("P", 2), rep(c(3, 10), each = 2), rep(c("_R1", "_R2")), sep = "")
)

pca.plot.annotation$media = "D0, 3, 7"
pca.plot.annotation[grep("PM", sample_labels)]$media = "Primed"
pca.plot.annotation[grep("NM", sample_labels)]$media = "t2iLGoY"
pca.plot.annotation$media = factor(pca.plot.annotation$media, 
                           levels = c("D0, 3, 7", "Primed", "t2iLGoY"))
pca.plot.annotation$timept = tstrsplit(pca.plot.annotation$sample_labels, "_")[[2]]
pca.plot.annotation$timept = factor(pca.plot.annotation$timept, 
                            levels = unique(pca.plot.annotation$timept))

pca <- prcomp_irlba(t(peaks.fpkm))

pca.plot.annotation$PC1 = pca$x[,1]
pca.plot.annotation$PC2 = pca$x[,2]
pca.plot.annotation$PC3 = pca$x[,3]

ggplot(pca.plot.annotation, aes(PC1, PC2, color=media, label=timept)) +
  geom_point(size = 3) + geom_text_repel(size = 5) + 
  scale_color_manual(values = colorMedia[1:3]) +
  scale_shape_manual(values = c(16,15)) + theme_classic(base_size = 24)

ggplot(pca.plot.annotation, aes(PC1, PC3, color=media, label=timept)) +
  geom_point(size = 3) + geom_text_repel(size = 5) + 
  scale_color_manual(values = colorMedia[1:3]) +
  scale_shape_manual(values = c(16,15)) + theme_classic(base_size = 24)


## Fuzzy clustering.
### Aggregate sample counts means and RPKM.
# Object to use 

peak.counts

conditions.labels <-
  c(
    "D0",
    "D3",
    "D7",
    "D13.*Primed",
    "D21.*Primed",
    "P3.*Primed",
    "P10.*Primed",
    "D13.*Smith",
    "D21.*Smith",
    "P3.*Smith",
    "P10.*Smith"
  )



aggregated.peak.counts <-
  sapply(conditions.labels, function(x)
    rowSums(peak.counts[, grep(x,
                               colnames(peak.counts),
                               value = TRUE)]))

colnames(aggregated.peak.counts) <-
  colnames(aggregated.peak.counts) %>% gsub("\\.\\*", "_", .)

aggregated.peak.fpkm <-
  rpkm(aggregated.peak.counts,
       log = F,
       gene.length = peak.counts$Length)


aggregated.keep <- rowSums(aggregated.peak.fpkm > 10) >= 1

aggregated.peak.fpkm <- aggregated.peak.fpkm[aggregated.keep, ]
aggregated.peak.fpkm <- log2(aggregated.peak.fpkm + 1)
aggregated.peak.fpkm <- normalizeQuantiles(aggregated.peak.fpkm)

### Discard low coefficient of variation peaks.

cv.20 <- apply(aggregated.peak.fpkm, 1, function(x) cv(x) * 100) > 20
cv.20 %>% table()

aggregated.peak.fpkm.high.cv <- aggregated.peak.fpkm[cv.20, ]

### Fuzzy clustering

exprSet <- ExpressionSet(assayData = aggregated.peak.fpkm.high.cv)

s.exprSet <- standardise(exprSet)

m <- mestimate(s.exprSet)

set.seed(1234)
cl.s.exprSet <- mfuzz(s.exprSet, c = 8, m = m)

#### Filter high affinity peaks and annotate

core.exprSet <- acore(s.exprSet,
                      cl = cl.s.exprSet,
                      min.acore = 0.8)

core.exprSet.cluster.counts <- lapply(core.exprSet, nrow) %>%
  unlist() %>% as.data.frame() %>% 
  mutate(Cluster = c(1:8))

colnames(core.exprSet.cluster.counts)[1] <- "members"


core.exprSet.peaks.ids <- sapply(c(1:8), function(x) as.vector(core.exprSet[[x]]$NAME)
)


colnames(cl.s.exprSet.high.cv.peak.ids)[1] <- "Cluster"
s.exprSet.high.cv.ann <- as.data.frame(s.exprSet.high.cv@assayData$exprs) %>% 
  dplyr::mutate(., PeakID = rownames(.))

cl.s.exprSet.high.cv.ann <- dplyr::inner_join(s.exprSet.high.cv.ann,
                                              cl.s.exprSet.high.cv.peak.ids,
                                              by = "PeakID")




high.cluster.affinity.peak.ids <- as.character(unlist(high.cluster.affinity.peak.ids))

cl.s.exprSet.high.cv.high.affinity.peaks.ann <- dplyr::filter(cl.s.exprSet.high.cv.ann,
                                                              PeakID %in% high.cluster.affinity.peak.ids)

cl.s.exprSet.high.cv.high.affinity.peaks.ann.homer <- dplyr::inner_join(cl.s.exprSet.high.cv.high.affinity.peaks.ann,
                                                                        hg19.atac.seq.intersect.biol.rep.peak.details, by = "PeakID")

comparelists(core.exprSet.high.cv[[2]]$NAME,
             filter(cl.s.exprSet.high.cv.high.affinity.peaks.ann, Cluster == 2)$PeakID)


lapply(c(1:8), 
       function(x) write.table(dplyr::filter(cl.s.exprSet.high.cv.high.affinity.peaks.ann.homer, Cluster == x), 
                               file = paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/fuzzy_clustering/combined_sets/", "high_affinity_peaks_cluster_",x,".tsv",sep=""), 
                               append=F,
                               quote=F,
                               sep="\t", 
                               col.names=T,
                               row.names=F, 
                               eol="\n")
)
### Stable cluster
cl.s.exprSet.high.cv.high.affinity.peaks.ann.homer <- dplyr::inner_join(cl.s.exprSet.high.cv.high.affinity.peaks.ann,
                                                                        hg19.atac.seq.intersect.biol.rep.peak.details, by = "PeakID")

```
#### Stable cluster processing.
```{r stable cluster processing}
hg19.atac.seq.intersect.biol.rep.sum.peak.counts.rpkm.log2.qnorm.low.cv

exprSet.low.cv <- ExpressionSet(assayData = hg19.atac.seq.intersect.biol.rep.sum.peak.counts.rpkm.log2.qnorm.low.cv)

s.exprSet.low.cv <- standardise(exprSet.low.cv)
s.exprSet.low.cv.ann <- as.data.frame(s.exprSet.low.cv@assayData$exprs)
s.exprSet.low.cv.ann$PeakID <- rownames(s.exprSet.low.cv.ann)
s.exprSet.low.cv.ann.melt <- melt(s.exprSet.low.cv.ann, id = "PeakID")
### ADD MUTATE IFELSE FOR MEDIA AND STAGE POST MELT


s.exprSet.low.cv.ann.melt <-
  s.exprSet.low.cv.ann.melt %>% mutate(stage = ifelse(variable == "D0", "D0", ifelse(
    variable == "D3", "D3", ifelse(
      variable == "D7",
      "D7",
      ifelse(
        variable == "D13_Primed" |
          variable == "D13_Smith",
        "D13",
        ifelse(
          variable == "D21_Primed" |
            variable == "D21_Smith",
          "D21",
          ifelse(variable == "P3_Primed" |
                   variable == "P3_Smith", "P3", "P10")
        )
      )
    )
  )))

s.exprSet.low.cv.ann.melt <- s.exprSet.low.cv.ann.melt %>%
  mutate(media = case_when(
    grepl("D0|D3|D7", variable) ~ "FM",
    grepl("Smith", variable) ~ "NM",
    grepl("Primed", variable) ~ "PM"
  ))

s.exprSet.low.cv.ann.melt$stage <-
  factor(s.exprSet.low.cv.ann.melt$stage,
         levels = c("D0", "D3", "D7", "D13", "D21", "P3", "P10"))
s.exprSet.low.cv.ann %>% dim()

ggplot(data = s.exprSet.low.cv.ann.melt %>% filter(media %in% c("FM", "NM")), 
       aes(stage, value, group = PeakID, colour = media)) +
  geom_line(alpha = 1/500) +
  stat_summary(aes(group = 1), fun.y = mean,   geom = "line",size = 0.75, color = "black") +
  stat_summary(aes(group = 1), fun.y = mean,   geom = "point",size = 2) +
  geom_line(data = s.exprSet.low.cv.ann.melt %>% filter(media %in% c("FM", "NM")),
            alpha = 1/500) +
  stat_summary(data = s.exprSet.low.cv.ann.melt %>% filter(media %in% c("FM", "PM")),
               aes(group = 1), fun.y = mean,   geom = "line", color = "black", size = 0.75) +
  stat_summary(data = s.exprSet.low.cv.ann.melt %>% filter(media %in% c("FM", "PM")),
               aes(group = 1), fun.y = mean,   geom = "point",size = 2) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#2302FE", "PM" = "#FF8D03"), name = "Culture Media") +
  ylab("Gene expression (z-score)") +
  xlab("Reprogramming stages") +
  ylim(c(-3, 3)) +
  # ggtitle("Fuzzy clustering")  +
  geom_text(data = NULL, aes(x = 1, y = 3, label=paste("n = ", "26361")),colour="black", inherit.aes=FALSE, parse=FALSE, size = 3, fontface = "bold") +
  gplot.theme

ggsave("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/fuzzy_clustering/fuzzy_clusters_combined_sets.pdf", height = 20, width = 25, units = "cm", dpi = 150, scale = 1)

```
##### Export annotated stable clustering for Homer processing.
```{r}
s.exprSet.low.cv.ann$Cluster <- 9
s.exprSet.low.cv.ann.homer <- dplyr::inner_join(s.exprSet.low.cv.ann,
                                                hg19.atac.seq.intersect.biol.rep.peak.details,
                                                by = "PeakID")
lapply(c(9), 
       function(x) write.table(dplyr::filter(s.exprSet.low.cv.ann.homer, Cluster == x), 
                               file = paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/fuzzy_clustering/combined_sets/", "high_affinity_peaks_cluster_",x,".tsv",sep=""), 
                               append=F,
                               quote=F,
                               sep="\t", 
                               col.names=T,
                               row.names=F, 
                               eol="\n")
)
```

### Combine high CV and stable cluster HOMER objects
```{r high cv plus stable}
head(s.exprSet.low.cv.ann.homer) %>% dim()
head(cl.s.exprSet.high.cv.high.affinity.peaks.ann.homer)
cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.homer <-
  dplyr::bind_rows(cl.s.exprSet.high.cv.high.affinity.peaks.ann.homer,
                   s.exprSet.low.cv.ann.homer)

```

## Integrate annotated peaks with cluster information.

```{r integrate annotated peaks}
homer.annotated.clusters <- lapply(Sys.glob("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/fuzzy_clustering/combined_sets/cluster_analysis/annotated_cluster_peaks/*_annotated_table.tsv"), function(x) read.delim(x, sep = "\t", header = T, stringsAsFactors = F))

homer.annotated.clusters.names <- Sys.glob("/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/fuzzy_clustering/combined_sets/cluster_analysis/annotated_cluster_peaks/*_annotated_table.tsv") %>% unlist() %>% gsub("/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/fuzzy_clustering/combined_sets/cluster_analysis/annotated_cluster_peaks/", "", .) %>% gsub("_homer_annotated_table.tsv", "", .)
names(homer.annotated.clusters) <- homer.annotated.clusters.names
names(homer.annotated.clusters)

### Combined annotated dynamic clusters with stable clusters
cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann <- dplyr::bind_rows(cl.s.exprSet.high.cv.high.affinity.peaks.ann,
                                                                        s.exprSet.low.cv.ann)

homer.annotated.clusters.details <- lapply(names(homer.annotated.clusters), function(x) inner_join(homer.annotated.clusters[[x]], filter(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann[, c(12:13, 1:11)], Cluster == gsub(".*\\_", "", x)), by = "PeakID"))

names(homer.annotated.clusters.details) <- names(homer.annotated.clusters)


homer.annotated.gene.by.atac.seq.cluster <- ldply(homer.annotated.clusters.details, rbind) %>% filter(!Nearest.Ensembl == "" ) %>% select(c("Nearest.Ensembl", "Cluster"))
colnames(homer.annotated.gene.by.atac.seq.cluster)[1] <- "Ensembl_id"

homer.annotated.clusters.details.table <- ldply(homer.annotated.clusters.details, rbind) %>% filter(!Nearest.Ensembl == "" ) %>% .[, -1]

homer.annotated.gene.by.atac.seq.changing.cluster <- homer.annotated.gene.by.atac.seq.cluster %>% filter(Cluster != 9)
homer.annotated.changing.clusters.details.table <- homer.annotated.clusters.details.table %>% filter(Cluster != 9)

```
## Subset annotated peaks closes to the nearest TSS
### All clusters.
```{r subsetting annotated clusters to TSS}
homer.annotated.clusters.details.table.nearest.TSS <- as.data.frame(homer.annotated.clusters.details.table) %>%
  dplyr::group_by(Nearest.Ensembl) %>% 
  dplyr::slice(which.min(abs(Distance.to.TSS)))

dplyr::filter(homer.annotated.clusters.details.table, Nearest.Ensembl == "ENSG00000146674")
dplyr::filter(homer.annotated.clusters.details.table.nearest.TSS, Nearest.Ensembl == "ENSG00000146674")

homer.annotated.clusters.details.table.nearest.TSS.cluster.ids <- homer.annotated.clusters.details.table.nearest.TSS %>% dplyr::select("Nearest.Ensembl", "Cluster")
colnames(homer.annotated.clusters.details.table.nearest.TSS.cluster.ids)[1] <- "Ensembl_id"



dplyr::filter(homer.annotated.clusters.details.table, Nearest.Ensembl %in% grep(homer.annotated.clusters.details.table.gene.ids[1], homer.annotated.clusters.details.table$Nearest.Ensembl, value = T))

dplyr::filter(homer.annotated.clusters.details.table.nearest.TSS, Nearest.Ensembl == "ENSG00000109182")
```
#### Distance to TSS plots
```{r}
distance.to.tss.table <- homer.annotated.clusters.details.table.nearest.TSS[, c("PeakID", "Distance.to.TSS", "Cluster")]

distance.to.tss.table <- distance.to.tss.table %>% mutate(log10_abs_dist_TSS = log10(abs(Distance.to.TSS) + 1))
distance.to.tss.table.cluster.counts <- distance.to.tss.table %>% group_by(Cluster) %>% dplyr::count(Cluster)

ggplot(data = distance.to.tss.table, aes(x = log10_abs_dist_TSS)) +
  geom_histogram(bins = 50) +
  scale_x_discrete(limits = c(0, 1, 2, 3, 4, 5, 6), 
                   labels = c( "TSS", "10b", "100b", "1Kb", "10Kb", "100kb", "1Mb")) +
  facet_wrap(~ Cluster, ncol = 1, nrow = 9, scales = "free") +
  labs(x = "Absolute distance to TSS (log10 TSS + 1", y = "Peaks (no.)") + 
  # geom_vline(xintercept = 3) +
  geom_text(
    data = distance.to.tss.table.cluster.counts,
    aes(
      x = 0,
      y = Inf,
      vjust = 1.25,
      label = paste("n = ", n)
    ),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  gplot.theme
```
### Changing clusters.

```{r subsetting changing annotated clusters to TSS}
homer.annotated.changing.clusters.details.table.nearest.TSS <- as.data.frame(homer.annotated.changing.clusters.details.table) %>%
  dplyr::group_by(Nearest.Ensembl) %>% 
  dplyr::slice(which.min(abs(Distance.to.TSS)))

dplyr::filter(homer.annotated.changing.clusters.details.table, Nearest.Ensembl == "ENSG00000146674")
dplyr::filter(homer.annotated.changing.clusters.details.table.nearest.TSS, Nearest.Ensembl == "ENSG00000146674")

homer.annotated.changing.clusters.details.table.nearest.TSS.cluster.ids <- homer.annotated.changing.clusters.details.table.nearest.TSS %>% dplyr::select("Nearest.Ensembl", "Cluster")
colnames(homer.annotated.changing.clusters.details.table.nearest.TSS.cluster.ids)[1] <- "Ensembl_id"

dplyr::filter(homer.annotated.changing.clusters.details.table, Nearest.Ensembl %in% grep(homer.annotated.clusters.details.table.gene.ids[1], homer.annotated.changing.clusters.details.table$Nearest.Ensembl, value = T))

dplyr::filter(homer.annotated.changing.clusters.details.table.nearest.TSS, Nearest.Ensembl == "ENSG00000109182")
```
#### Distance to TSS plots
```{r}
distance.to.tss.changing.clusters.table <- homer.annotated.changing.clusters.details.table.nearest.TSS[, c("PeakID", "Distance.to.TSS", "Cluster")]

distance.to.tss.changing.clusters.table <- distance.to.tss.changing.clusters.table %>% mutate(log10_abs_dist_TSS = log10(abs(Distance.to.TSS) + 1))
distance.to.tss.changing.clusters.table.counts <- distance.to.tss.changing.clusters.table %>% group_by(Cluster) %>% dplyr::count(Cluster)

distance.to.tss.changing.clusters.table$Cluster <- factor(distance.to.tss.changing.clusters.table$Cluster, levels = c(1, 2, 5, 4, 7, 6, 3, 8))


distance.to.tss.changing.clusters.labels <- c("1" = "C 1 (n = 2784)", "2" = "C 2 (n = 2017)", "3" = "C 7 (n = 1219)", "4" = "C 4 (n = 602)", "5" = "C 3 (n = 1333)", "6" = "C 6 (n = 1963)", "7" = "C 5 (n = 2133)", "8" = "C 8 (n = 2323)")

distance.to.tss.changing.clusters.plot <-
  ggplot(data = distance.to.tss.changing.clusters.table, aes(x = log10_abs_dist_TSS)) +
  geom_histogram(bins = 50) +
  scale_x_discrete(
    limits = c(0, 1, 2, 3, 4, 5, 6),
    labels = c("TSS", "10b", "100b", "1Kb", "10Kb", "100kb", "1Mb")
  ) +
  facet_wrap(
    ~ Cluster,
    ncol = 1,
    nrow = 8,
    # scales = "free",
    labeller = labeller(Cluster = distance.to.tss.changing.clusters.labels)
  ) +
  labs(x = "Absolute distance to TSS\n(log10 TSS + 1)", y = "Peaks (no.)") +
  # # geom_vline(xintercept = 3) +
  # geom_text(
  # data = distance.to.tss.table.cluster.counts,
  # aes(
  # x = 0,
  # y = Inf,
  # vjust = 1.25,
  # label = paste("n = ", n)
  # ),
  # colour = "black",
  # inherit.aes = FALSE,
# parse = FALSE,
# size = 3,
# fontface = "bold"
# ) +
gplot.theme + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7))

ggsave(plot = distance.to.tss.changing.clusters.plot,
       "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Supp_Figure_5D_ATAC-seq_distance_to_tss_reordered_resized.pdf",
       height = 18,
       width = 7.5,
       units = "cm",
       dpi = 200,
       scale = 1,
       device = "pdf"
)
```



# Integration of ATAC-seq with bulk RNA-seq
## RPKM from raw counts
```{r RPKM from raw counts}
bulk.RNA.seq.raw.counts <-
  read.table(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_RNA-seq_Human_Reprogramming/analysis/data/bulkRNA_rawCounts.tab",
             header = T,
             stringsAsFactors = F)

keep.min.bulk.RNA.seq.raw.counts <-
  apply(bulk.RNA.seq.raw.counts[, 4:ncol(bulk.RNA.seq.raw.counts)], 1 , max) >= 10
table(keep.min.bulk.RNA.seq.raw.counts)

keep.cpm.bulk.RNA.seq.raw.counts <- rowSums(cpm(bulk.RNA.seq.raw.counts[, 4:ncol(bulk.RNA.seq.raw.counts)]) > 2.0) >= 2           

table(keep.cpm.bulk.RNA.seq.raw.counts)      

keep.bulk.RNA.seq.raw.counts <- keep.min.bulk.RNA.seq.raw.counts & keep.cpm.bulk.RNA.seq.raw.counts

table(keep.bulk.RNA.seq.raw.counts)

bulk.RNA.seq.filtered.low.count.genes <- bulk.RNA.seq.raw.counts[keep.bulk.RNA.seq.raw.counts, ]

colnames(bulk.RNA.seq.filtered.low.count.genes)
rownames(bulk.RNA.seq.filtered.low.count.genes) <- bulk.RNA.seq.filtered.low.count.genes$geneID
dim(bulk.RNA.seq.filtered.low.count.genes)

y.bulk.RNA.seq.filtered.low.count.genes <-
  DGEList(counts = bulk.RNA.seq.filtered.low.count.genes[, 4:ncol(bulk.RNA.seq.filtered.low.count.genes)],
          genes = bulk.RNA.seq.filtered.low.count.genes[, 1:3])

y.bulk.RNA.seq.filtered.low.count.genes <- calcNormFactors(y.bulk.RNA.seq.filtered.low.count.genes)

log2.rpkm.bulk.RNA.seq.filtered.low.count.genes <-
  rpkm(
    y.bulk.RNA.seq.filtered.low.count.genes,
    normalized.lib.sizes = T,
    log = T,
    prior.count = 1
  )

```
### Mean log RPKM
```{r mean log2 RPKM}
bulk.RNA.seq.conditions.labels <- c("Fibroblast_D00", "Fibroblast_D03", "Fibroblast_D07", "Primed_D13", "Primed_D21", "Primed_P03", "Primed_P10", "t2iLGoY_D13", "t2iLGoY_D21", "t2iLGoY_P03", "t2iLGoY_P10")



mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes <-
  sapply(bulk.RNA.seq.conditions.labels, function(x)
    rowMeans(log2.rpkm.bulk.RNA.seq.filtered.low.count.genes[, grep(x,
                                                                    colnames(log2.rpkm.bulk.RNA.seq.filtered.low.count.genes),
                                                                    value = TRUE)])) %>% as.data.frame()


mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes$Ensembl_id <- rownames(log2.rpkm.bulk.RNA.seq.filtered.low.count.genes)

plotMDS(mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes[, -12],
        top = 15000,
        gene.selection = "common")

hist(mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes[, 1])

write.table(mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes,
            file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_RNA-seq_Human_Reprogramming/analysis/data/mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.tab",
            sep = "\t",
            row.names = F,
            quote = F)

```
### Z-score of all genes (not discarding low CV) all clusters.
```{r}
exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes <- ExpressionSet(assayData = as.matrix(mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes[, -12]))

exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes@assayData$exprs %>% head()


s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes <-
  standardise(exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes)

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes <-
  as.data.frame(s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes@assayData$exprs) %>%
  dplyr::mutate(., Ensembl_id = rownames(.))


s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster <-
  dplyr::inner_join(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes,
    homer.annotated.gene.by.atac.seq.cluster,
    by = "Ensembl_id"
  )

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster <-
  dplyr::inner_join(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes,
    homer.annotated.gene.by.atac.seq.changing.cluster,
    by = "Ensembl_id"
  )
```
#### All clusters. Select peaks closest to TSS.
```{r closest to tss all genes}

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss <-  dplyr::inner_join(s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes,
                                                                                                          homer.annotated.clusters.details.table.nearest.TSS.cluster.ids, by = "Ensembl_id")


```
#### Changing clusters. Select peaks closest to TSS.
```{r closest to tss all genes}

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss <-  dplyr::inner_join(s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes,
                                                                                                                   homer.annotated.changing.clusters.details.table.nearest.TSS.cluster.ids, by = "Ensembl_id")

intersect(filter(as.data.frame(homer.annotated.changing.clusters.details.table.nearest.TSS.cluster.ids), Cluster == "1") %>% .$Ensembl_id,
          s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes$Ensembl_id)

```

### Z-score of genes with more than 10 %CV - RPKM based.
Standardise and annotate based on ATAC-seq cluster association. RPKM
```{r standardise rna-seq association dicard low cv}
## Discard low CV genes
RNA.seq.bulk.cv.index <- apply(mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes[, -12], 1, function(x) cv(x) * 100) < 10

mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv <- mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes[!RNA.seq.bulk.cv.index, ]

exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv <- ExpressionSet(assayData = as.matrix(mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv[, -12]))

exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv@assayData$exprs %>% head()


s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv <-
  standardise(exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv)

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv <-
  as.data.frame(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv@assayData$exprs
  ) %>%
  dplyr::mutate(., Ensembl_id = rownames(.))


s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv.atac.seq.cluster <-
  dplyr::inner_join(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv,
    homer.annotated.gene.by.atac.seq.cluster,
    by = "Ensembl_id"
  )


```
#### Closest to TSS peaks - RPKM based analysis
Genes with a CV of more than 10%.
```{r subsetted closest to TSS peaks}

colnames(homer.annotated.clusters.details.table.nearest.TSS.cluster.ids)[1] <- "Ensembl_id"


s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv.atac.seq.cluster.tss <- dplyr::inner_join(s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.high.cv,
                                                                                                                 homer.annotated.clusters.details.table.nearest.TSS.cluster.ids, by = "Ensembl_id")

```

# Plots
## GGPlots theme
```{r}
gplot.theme <- theme(panel.border=element_blank(), 
                     panel.background=element_blank(),
                     axis.line=element_line(),
                     axis.line.y=element_line("black", size=1),
                     axis.line.x=element_line("black", size=1),
                     axis.text=(element_text(colour="black", size = 6, face="bold")),
                     axis.title=(element_text(colour="black", size = 8, face="bold")),
                     axis.ticks = (element_line("black", size=0.5)),
                     legend.title = element_text(colour="black", size = 9, face="bold"),
                     legend.text = element_text(colour="black", size= 8, face="bold"),
                     legend.key.size=unit(0.1,"cm"),
                     # legend.margin = unit(0, "cm"),
                     legend.position="none",
                     plot.margin=unit(c(0,2,0,0),"mm"),
                     strip.text.x = element_text(size = 6, colour = "black", face = "bold", margin = margin(0,0,0,0, "cm"))
)

```

## MDS plots
```{r mds plots}
plotMDS(hg19.atac.seq.intersect.biol.rep.rpkm.log2.qnorm,
        gene.selection = "common",
        top = nrow(hg19.atac.seq.intersect.biol.rep.rpkm.log2.qnorm),
        dim.plot = c(1, 2)
)
MDS.plot <- plotMDS(
  # log2.cpm,
  # log2.cpm[, -c(15:18)],
  # log2.cpm[, fib.prime],
  # log2.cpm[, fib.smith],
  # log2.cpm[, f.32],
  # log2.cpm[, f.55],
  # hg19.atac.seq.intersect.biol.rep.cpm,
  # hg19.atac.seq.intersect.biol.rep.strict.cpm,
  hg19.atac.seq.intersect.biol.rep.rpkm.log2.qnorm[, levels(samples.order)],
  gene.selection = "common",
  dim.plot = c(1, 2),
  top = nrow(hg19.atac.seq.intersect.biol.rep.rpkm.log2.qnorm)
)
variance.contribution <- cmdscale(as.dist(MDS.plot$distance.matrix), eig = T, k = 21)

variance.contribution$eig[2]/sum(variance.contribution$eig)

glMDSPlot(hg19.atac.seq.intersect.biol.rep.rpkm.log2.qnorm,
          gene.selection = "common",
          top = nrow(hg19.atac.seq.intersect.biol.rep.rpkm.log2.qnorm
          ))

MDS.plot.coordinates <- data.frame(x = MDS.plot$x,
                                   y = MDS.plot$y)
MDS.plot.coordinates$sample.name <- rownames(MDS.plot.coordinates)
MDS.plot.coordinates$stage <- MDS.plot.coordinates$sample.name %>%  gsub("ATAC[_:.]", "", .) %>% gsub("_rep2", "", .) %>% gsub("\\.", "\\_", .) %>% gsub("\\_15.*", "", .) %>% gsub("Smith\\_R", "SmithR", .) %>% gsub("\\_KSR", "", .) %>% gsub("_HDFa", "", .) %>% gsub("_32F|_55F", "", .) %>% gsub("_SmithR|_Primed", "", .)



MDS.plot.coordinates$sample.name <- factor(MDS.plot.coordinates$sample.name, levels = MDS.plot.coordinates$sample.name)
MDS.plot.coordinates$colour <- as.character(MDS.plot.coordinates$sample.name)
MDS.plot.coordinates$colour[grep("D0|D3|D7",MDS.plot.coordinates$colour)] <- "black"
MDS.plot.coordinates$colour[grep("Primed",MDS.plot.coordinates$colour)] <- "#FF8C00"
MDS.plot.coordinates$colour[grep("Smith",MDS.plot.coordinates$colour)] <- "#0000FF"

MDS.plot.coordinates$media <- as.character(MDS.plot.coordinates$sample.name)
MDS.plot.coordinates$media[grep("D0|D3|D7",MDS.plot.coordinates$media)] <- "FM"
MDS.plot.coordinates$media[grep("Primed",MDS.plot.coordinates$media)] <- "PM"
MDS.plot.coordinates$media[grep("Smith",MDS.plot.coordinates$media)] <- "NM"


MDS.ggplot <- ggplot(data=NULL, aes(x = x, y = y)) +
  # geom_line(data=bezier.points.2d.f[1:48,], aes(x,y), color=factor(bezier.points.2d.f.color[1:48]), lwd=1, alpha=0.95) + 
  geom_point(data = MDS.plot.coordinates, aes(colour = media), pch = 19, size = 2 , alpha= 1) +
  scale_color_manual(breaks = MDS.plot.coordinates$media, values = c("FM" = "black", "PM" = "#FF8C00", "NM" = "#0000FF"), name = "Culture Media") +
  # xlim(-1,1) +
  # ylim(-0.75,1) +
  xlab("PC 1 (44.1%)") +
  # ylab("PC 3 (10.7%)") +
  ylab("PC 2 (17.8%)") +
  # ggtitle("MDS plot of D0-D7 Fs, Naive and Prime\n(counts of biol. reps. intersection−rescued peaks)") +
  geom_text_repel(data = MDS.plot.coordinates, aes(label= MDS.plot.coordinates$stage), 
                  fontface = "plain",
                  # fontface = "bold",
                  # size = 2,
                  size = 2.5,
                  # segment.color = "darkblue", 
                  # color="darkblue", 
                  segment.size = 0.25,# C 1 and C 3
                  point.padding = unit(0.5, 'lines')#, C 1 and C 3
                  # force = 1.5#,
                  # force = 1#,
                  # arrow = arrow(length = unit(1, 'mm'))
  ) +
  labs(title = NULL) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  theme(panel.border=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.line.y=element_line("black", size=0.35),
        axis.line.x=element_line("black", size=0.35),
        axis.text=(element_text(colour="black", size = 8, face="bold")),
        axis.title=(element_text(colour="black", size = 8, face="bold")),
        axis.ticks = (element_line("black", size=0.25)),
        legend.title = element_text(colour="black", size = 9, face="bold"),
        legend.text = element_text(colour="black", size= 8, face="bold"),
        legend.key.size=unit(0.1,"cm"),
        legend.margin = unit(0, "cm"),
        legend.position="none",
        plot.margin=unit(c(0,2,0,0),"mm")
  )

print(MDS.ggplot)
ggsave("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Rebuttal/20191111/MDS_plot_ATAC-seq_c_1_c2.pdf", height = 8, width = 11, units = "cm", dpi = 200, scale = 1)
```

## Fuzzy clustering plots
### Changing clusters
```{r fuzz clustering plots}
cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt <- melt(cl.s.exprSet.high.cv.high.affinity.peaks.ann, id = c("PeakID", "Cluster"))
### ADD MUTATE IFELSE FOR MEDIA AND STAGE POST MELT
cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt.prime <- dplyr::filter(cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt, variable %in% 
                                                                           c("D0", "D3", "D7", "D13_Primed", "D21_Primed", "P3_Primed", "P10_Primed"))
cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt.prime %>% filter(PeakID %in% paste("consensus_peak_", 561:570, sep = ""))
cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt.naive <- dplyr::filter(cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt, variable %in% 
                                                                           c("D0", "D3", "D7", "D13_Smith", "D21_Smith", "P3_Smith", "P10_Smith"))

cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt <- cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt %>% mutate(stage = ifelse(variable == "D0", "D0",ifelse(variable == "D3", "D3", ifelse(variable == "D7", "D7", ifelse(variable == "D13_Primed" | variable == "D13_Smith", "D13", ifelse(variable == "D21_Primed" | variable == "D21_Smith", "D21", ifelse(variable == "P3_Primed" | variable == "P3_Smith", "P3", "P10")))))))

cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt <- cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt %>% 
  mutate(media = case_when(grepl("D0|D3|D7", variable) ~ "FM",
                           grepl("Smith", variable) ~ "NM",
                           grepl("Primed", variable) ~ "PM"))

cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt$stage <- factor(cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt$stage, 
                                                                  levels = c("D0", "D3", "D7", "D13", "D21", "P3", "P10"))


cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt$Cluster <- factor(cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt$Cluster, levels = c(1, 2, 5, 4, 7, 6, 3, 8))
cluster.labels.atac.seq <- c("1" = "C 1 (n = 12024)", "2" = "C 2 (n = 7779)", "3" = "C 7 (n = 4885)", "4" = "C 4 (n = 3334)", "5" = "C 3 (n = 5077)", "6" = "C 6 (n = 10129)", "7" = "C 5 (n = 9117)", "8" = "C 8 (n = 7739)")

all.changing.clusters.ggplot <-
  ggplot(data = cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt %>% filter(media %in% c("FM", "PM")),
         aes(stage, value, group = PeakID, colour = media)) +
  geom_line(alpha = 1 / 50) +
  stat_summary(
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    size = 0.75,
    color = "black"
  ) +
  stat_summary(aes(group = 1),
               fun.y = mean,
               geom = "point",
               size = 2) +
  geom_line(data = cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt %>% filter(media %in% c("FM", "NM")),
            alpha = 1 / 50) +
  stat_summary(
    data = cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt %>% filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    color = "black",
    size = 0.75
  ) +
  stat_summary(
    data = cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt %>% filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "point",
    size = 2
  ) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#2302FE", "PM" = "#FF8C00"),
                     name = "Culture Media") +
  facet_wrap(~ Cluster, nrow = 1, ncol = 8, labeller = labeller(Cluster = cluster.labels.atac.seq)) +
  # facet_wrap( ~ Cluster, nrow = 8, ncol = 1, labeller = paste("C", 1:8)) +
  ylab("Chromatin accessibility\n(Z-scaling)") +
  xlab("Reprogramming stages") +
  # ggtitle("Fuzzy clustering")  +
  # geom_text(
  # data = gene.counts.core.exprSet.high.cv,
  # aes(
  # x = 2,
  # y = 3,
  # label = paste("n = ", members)
  # ),
  # colour = "black",
  # inherit.aes = FALSE,
  # parse = FALSE,
# size = 3,
# fontface = "bold"
# ) +
gplot.theme + theme(legend.position = "none")

print(all.changing.clusters.ggplot)
ggsave(plot = all.changing.clusters.ggplot,
       "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_4_ATAC-seq_fuzzy_clustering_resized_horizontal.pdf",
       height = 5,
       width = 28,
       units = "cm",
       dpi = 200,
       scale = 1,
       device = "pdf"
)
```
### All clusters.
```{r changing clusters + stable}

cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt
s.exprSet.low.cv.ann.melt$Cluster <- 9
s.exprSet.low.cv.ann.melt <- s.exprSet.low.cv.ann.melt[, c(1, 6, 2:5)]

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt <-
  bind_rows(cl.s.exprSet.high.cv.high.affinity.peaks.ann.melt,
            s.exprSet.low.cv.ann.melt)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt.counts <- cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt %>% dplyr::count(Cluster) %>% mutate(n = n/levels(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt$variable) %>% length())

ggplot(data = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt %>% filter(media %in% c("FM", "PM")),
       aes(stage, value, group = PeakID, colour = media)) +
  geom_line(alpha = 1 / 500) +
  stat_summary(
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    size = 0.75,
    color = "black"
  ) +
  stat_summary(aes(group = 1),
               fun.y = mean,
               geom = "point",
               size = 2) +
  geom_line(data = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt %>% filter(media %in% c("FM", "NM")),
            alpha = 1 / 500) +
  stat_summary(
    data = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt %>% filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    color = "black",
    size = 0.75
  ) +
  stat_summary(
    data = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt %>% filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "point",
    size = 2
  ) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#2302FE", "PM" = "#FF8D03"),
                     name = "Culture Media") +
  facet_wrap( ~ Cluster, nrow = 8, ncol = 1) +
  ylab("Gene expression (z-score)") +
  xlab("Reprogramming stages") +
  # ggtitle("Fuzzy clustering")  +
  geom_text(
    data = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.melt.counts,
    aes(
      x = 2,
      y = 3,
      label = paste("n = ", members)
    ),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  gplot.theme

```

## Bulk RNA-seq expression based on ATAC-seq clustering
### RPKM
#### All Clusters 
##### All Peaks.
```{r RPKM all peaks}
s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster %>% dim()

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt <-
  melt(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster,
    id = c("Ensembl_id", "Cluster")
  )

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt <-
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt %>%
  mutate(media = ifelse(grepl("Fibroblast", variable), "FM",
                        ifelse(grepl("Primed", variable), "PM", "NM")))  %>%
  mutate(stage = gsub(".*_", "", variable))

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt$stage <-
  factor(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt$stage,
    levels = c("D00", "D03", "D07", "D13", "D21", "P03", "P10"),
    labels = c("D0", "D3", "D7", "D13", "D21", "P3", "P10")
  )

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.counts <-
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster %>%
  group_by(Cluster) %>%
  dplyr::count()

ggplot(data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt %>%
         filter(media %in% c("FM", "PM")),
       aes(stage, value, group = Ensembl_id, colour = media)) +
  # geom_line(alpha = 1 / 250) +
  stat_summary(aes(group = 1),
               fun.y = mean,
               geom = "line",
               size = 1) +
  # geom_line(data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt %>%
  #             filter(media %in% c("FM", "NM")),
  # alpha = 1 / 250) +
  stat_summary(
    data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.melt %>%
      filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    size = 1
  ) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#0000FF", "PM" = "#FF8C00"),
                     name = "Culture Media") +
  facet_wrap( ~ Cluster, nrow = 9, ncol = 1) +
  ylab("Gene expression (z-score)") +
  xlab("Reprogramming stages") +
  # ggtitle("Fuzzy clustering")  +
  geom_text(
    data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.counts,
    aes(
      x = 2,
      y = 0.6,
      label = paste("n = ", n)
    ),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  gplot.theme

```

##### Peaks closest to TSS - All detected genes.

```{r all detected genes closest to TSS peaks}

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt <-
  melt(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss,
    id = c("Ensembl_id", "Cluster")
  )

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt <-
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt %>%
  mutate(media = ifelse(grepl("Fibroblast", variable), "FM",
                        ifelse(grepl("Primed", variable), "PM", "NM")))  %>%
  mutate(stage = gsub(".*_", "", variable))

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt$stage <-
  factor(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt$stage,
    levels = c("D00", "D03", "D07", "D13", "D21", "P03", "P10"),
    labels = c("D0", "D3", "D7", "D13", "D21", "P3", "P10")
  )

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.counts <-
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss %>%
  group_by(Cluster) %>%
  dplyr::count()

ggplot(data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt %>%
         filter(media %in% c("FM", "PM")),
       aes(stage, value, group = Ensembl_id, colour = media)) +
  # geom_line(alpha = 1 / 250) +
  stat_summary(aes(group = 1),
               fun.y = mean,
               geom = "line",
               size = 1) +
  # geom_line(data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt %>%
  #             filter(media %in% c("FM", "NM")),
  # alpha = 1 / 250) +
  stat_summary(
    data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.melt %>%
      filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    size = 1
  ) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#0000FF", "PM" = "#FF8C00"),
                     name = "Culture Media") +
  facet_wrap( ~ Cluster, nrow = 9, ncol = 1) +
  ylab("Gene expression (z-score)") +
  xlab("Reprogramming stages") +
  # ggtitle("Fuzzy clustering")  +
  geom_text(
    data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.cluster.tss.counts,
    aes(
      x = 2,
      y = 0.6,
      label = paste("n = ", n)
    ),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 2.5,
    fontface = "bold"
  ) +
  gplot.theme

```

##### Peaks closest to TSS - High CV detected genes. 

#### Changing clusters - 
##### Peaks closest to TSS - All detected genes. 

```{r all detected genes closest to TSS peaks}

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss

gene.ids <- bulk.RNA.seq.raw.counts[, 1:2]
colnames(gene.ids)[1] <- "Ensembl_id"

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.annotated <-
  dplyr::inner_join(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss,
    gene.ids,
    by = "Ensembl_id"
  )
for (i in c(1:8)){
  write.table(s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.annotated %>% dplyr::filter(Cluster == i) %>% dplyr::select(c("Ensembl_id", "Cluster", "geneName")),
              file = paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/fuzzy_clustering/combined_sets/cluster_analysis/bulk_RNA-seq_genes_over_ATAC-seq_changing_clusters_peaks_TSS/",
                           "cluster_",
                           i,
                           "_gene_ids.csv",sep = ""),
              row.names = F,
              sep = "\t",
              quote = F)
}


s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt <-
  melt(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss,
    id = c("Ensembl_id", "Cluster")
  )

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt <-
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt %>%
  mutate(media = ifelse(grepl("Fibroblast", variable), "FM",
                        ifelse(grepl("Primed", variable), "PM", "NM")))  %>%
  mutate(stage = gsub(".*_", "", variable))

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt$stage <-
  factor(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt$stage,
    levels = c("D00", "D03", "D07", "D13", "D21", "P03", "P10"),
    labels = c("D0", "D3", "D7", "D13", "D21", "P3", "P10")
  )

s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.counts <-
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss %>%
  group_by(Cluster) %>%
  dplyr::count()



s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt$Cluster <- factor(s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt$Cluster, levels = c(1, 2, 5, 4, 7, 6, 3, 8))
cluster.labels.rna.seq <- c("1" = "C 1 (n = 1823)", "2" = "C 2 (n = 1347)", "3" = "C 7 (n = 790)", "4" = "C 4 (n = 339)", "5" = "C 3 (n = 892)", "6" = "C 6 (n = 1027)", "7" = "C 5 (n = 1239)", "8" = "C 8 (n = 1547)")
atac.seq.clustes.tss.peaks.bulk.RNA.seq.all.genes.ggplot <-
  ggplot(data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt %>%
           dplyr::filter(media %in% c("FM", "PM")),
         aes(stage, value, group = Ensembl_id, colour = media)) +
  stat_summary(aes(group = 1),
               fun.y = mean,
               geom = "line",
               size = 1) +
  stat_summary(
    data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt %>%
      dplyr::filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    size = 1
  ) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#0000FF", "PM" = "#FF8C00"),
                     name = "Culture Media") +
  facet_wrap(~ Cluster, nrow = 1, ncol = 8) +
  # facet_wrap(~ Cluster, nrow = 1, ncol = 8, labeller = labeller(Cluster = cluster.labels.rna.seq)) +
  ylab("Gene expression (z-score)") +
  xlab("Reprogramming stages") +
  # ggtitle("Fuzzy clustering")  +
  # geom_text(
  # data = s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.counts,
  # aes(
  # x = 2,
  # y = -0.4,
  # label = paste("n = ", n)
  # ),
  # colour = "black",
  # inherit.aes = FALSE,
  # parse = FALSE,
# size = 2.5,
# fontface = "bold"
# ) +
gplot.theme
print(atac.seq.clustes.tss.peaks.bulk.RNA.seq.all.genes.ggplot)
ggsave(
  plot = atac.seq.clustes.tss.peaks.bulk.RNA.seq.all.genes.ggplot,
  filename = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_4_bulk_RNA-seq_over_ATAC-seq_changing_clusters_peaks_tss_resized_reordered_horizontal_no_labels.pdf",
  height = 5,
  width = 28,
  units = "cm",
  dpi = 200,
  scale = 1,
  device = "pdf"
)

extended.figure.6.c.source.data <- select(
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt,
  -variable
)

extended.figure.6.c.source.data <- extended.figure.6.c.source.data %>% mutate(
  Cluster_Label = case_when(
    Cluster == "1" ~ "C 1",
    Cluster == "2" ~ "C 2",
    Cluster == "3" ~ "C 7",
    Cluster == "4" ~ "C 4",
    Cluster == "5" ~ "C 3",
    Cluster == "6" ~ "C 6",
    Cluster == "7" ~ "C 5",
    TRUE ~ "C 8"
  )
)
extended.figure.6.c.source.data$Cluster_Label  <-
  factor(x = extended.figure.6.c.source.data$Cluster_Label,
         levels = c("C 1", "C 2", "C 3", "C 4", "C 5", "C 6", "C 7", "C 8"))

# Plot test

ggplot(data = extended.figure.6.c.source.data %>%
         dplyr::filter(media %in% c("FM", "PM")),
       aes(stage, value, group = Ensembl_id, colour = media)) +
  stat_summary(aes(group = 1),
               fun.y = mean,
               geom = "line",
               size = 1) +
  stat_summary(
    data = extended.figure.6.c.source.data %>%
      dplyr::filter(media %in% c("FM", "NM")),
    aes(group = 1),
    fun.y = mean,
    geom = "line",
    size = 1
  ) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#0000FF", "PM" = "#FF8C00"),
                     name = "Culture Media") +
  facet_wrap(~ Cluster_Label, nrow = 1, ncol = 8) +
  ylab("Gene expression (z-score)") +
  xlab("Reprogramming stages") +
  gplot.theme

write_excel_csv(x = select(extended.figure.6.c.source.data,-Cluster),
                path = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/Source_Data/extended.figure.6.c.source.data.csv")
```

##### ATAC-seq peaks dynamics closest to TSS selected genes
```{r atac-seq peaks dynamics closest to tss seleceted genes}

homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene <- homer.annotated.changing.clusters.details.table.nearest.TSS %>% ungroup %>% dplyr::select("Gene.Name", "D0", "D3", "D7", "D13_Primed", "D21_Primed", "P3_Primed", "P10_Primed", "D13_Smith", "D21_Smith", "P3_Smith", "P10_Smith")

homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt <- melt(homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene, id = c("Gene.Name"))
### ADD MUTATE IFELSE FOR MEDIA AND STAGE POST MELT

homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt <- homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt %>% mutate(stage = ifelse(variable == "D0", "D0",ifelse(variable == "D3", "D3", ifelse(variable == "D7", "D7", ifelse(variable == "D13_Primed" | variable == "D13_Smith", "D13", ifelse(variable == "D21_Primed" | variable == "D21_Smith", "D21", ifelse(variable == "P3_Primed" | variable == "P3_Smith", "P3", "P10")))))))

homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt <- homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt %>% 
  mutate(media = case_when(grepl("D0|D3|D7", variable) ~ "FM",
                           grepl("Smith", variable) ~ "NM",
                           grepl("Primed", variable) ~ "PM"))

homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt$stage <-
  factor(
    homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt$stage,
    levels = c("D0", "D3", "D7", "D13", "D21", "P3", "P10")
  )

homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt.candidates <- dplyr::filter(homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt, Gene.Name %in% levels(
  z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl$geneName))

# homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt.candidates$Gene.Name <- factor(homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt.candidates$Gene.Name, levels = levels(
#   z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl$geneName))

homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt.candidates$Gene.Name <- factor(homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt.candidates$Gene.Name, levels = c("SNAI1", "TWIST2", "ZEB2", "BHLHE40", "F11R", "ZIC5", "OCLN", "NLRP7"))


rna.seq.over.atac.seq.dynamics <- ggplot(
  data = homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt.candidates %>%
    filter(
      media %in% c("FM", "PM")),
  aes(stage, value, group = Gene.Name, colour = media)
) +
  geom_line(size = 1.5) +
  geom_line(
    data = homer.annotated.changing.clusters.details.table.nearest.TSS.by.gene.melt.candidates %>%
      filter(
        media %in% c("FM", "NM")
      ),
    aes(stage, value, group = Gene.Name, colour = media), size = 1.5
  ) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#0000FF", "PM" = "#FF8C00"),
                     name = "Culture Media") +
  facet_wrap( ~ Gene.Name, nrow = 8, ncol = 1) +
  ylab("Chromatin Accessiblity (z-score)") +
  xlab("Reprogramming stages") +
  gplot.theme
ggsave(plot = rna.seq.over.atac.seq.dynamics,
       "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_4_bulk_RNA-seq_over_ATAC-seq_changing_clusters_line_plots_closest_to_tss_dynamics_reordered.pdf",
       height = 15,
       width = 4.5,
       units = "cm",
       dpi = 300,
       scale = 1
)
```
##### Selected gene expression plots - Cluster basd
```{r}
figure.4.bulk.RNA.seq.candidates.canonical.ensembl <- queryMany(figure.4.bulk.RNA.seq.candidates.canonical$hg19.kgXref.geneSymbol, fields = "ensembl.gene", scope = "symbol", species = "human")

z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl <-
  filter(
    s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt,
    Ensembl_id %in% c(
      figure.4.bulk.RNA.seq.candidates.canonical.ensembl$ensembl %>% unlist()
    )
  )
z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl <- inner_join(z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl, gene.ids, by = "Ensembl_id")
z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl$Cluster <- factor(z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl$Cluster, levels = c(1:8))

z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl$geneName <- factor(z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl$geneName,
                                                                              levels = c("SNAI1", "TWIST2", "ZEB2", "BHLHE40", "F11R", "ZIC5", "OCLN", "NLRP7"))
selected.genes.line.plots <- ggplot(data = z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl %>% 
                                      filter(media %in% c("FM", "PM")), 
                                    aes(stage, value, group = geneName, colour = media)) +
  geom_line(size = 1.5) +
  geom_line(data = z.score.figure.4.bulk.RNA.seq.candidates.canonical.ensembl %>% 
              filter(media %in% c("FM", "NM")), aes(stage, value, group = geneName, colour = media), size = 1.5) +
  scale_color_manual(values = c("FM" = "black", "NM" = "#0000FF", "PM" = "#FF8C00"), name = "Culture Media") +
  facet_wrap(~geneName, nrow = 8, ncol = 1) +
  ylab("Gene expression (z-score)") +
  xlab("Reprogramming stages") +
  gplot.theme
ggsave(plot = selected.genes.line.plots,
       "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_4_bulk_RNA-seq_over_ATAC-seq_changing_clusters_line_plots_reordered.pdf",
       height = 15,
       width = 4.5,
       units = "cm",
       dpi = 300,
       scale = 1
)
```

## ATAC-seq coverage plots (Wiggle plots)
### Ensembl annotation
```{r atacseq wiggle plots}
listMarts()
useMart("ensembl")
ensembl.mart <- useMart("ensembl", host = "grch37.ensembl.org")
listDatasets(ensembl.mart)
ensembl.grch37 <- useDataset("hsapiens_gene_ensembl", mart = ensembl.mart)
ensembl.grch37

attributes.grch37 <- listAttributes(ensembl.grch37)

selected.attributes = c("ensembl_transcript_id", "ensembl_gene_id", 
                        "external_gene_name", "strand", 
                        "gene_biotype", "transcript_biotype")
transcript.metadata.grch37 = getBM(attributes = selected.attributes, mart = ensembl.grch37)

transcript.metadata.grch37 <- dplyr::rename(transcript.metadata.grch37,
                                            transcript_id = ensembl_transcript_id,
                                            gene_id = ensembl_gene_id,
                                            gene_name = external_gene_name)

saveRDS(transcript.metadata.grch37, file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/objects/transcript_metadata_grch37.rds")

txdb.grch37 <- makeTxDbFromBiomart(#biomart = "ensembl",
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl", 
  host = "grch37.ensembl.org")

exons.grch37 <- exonsBy(txdb.grch37, by = "tx", use.names = T)
cdss.grch37 <- cdsBy(txdb.grch37, by = "tx", use.names = T)


plotTranscripts(exons = exons.grch37["ENST00000334384"], transcript_annotations = transcript.metadata.grch37)
plotTranscripts(exons = exons.grch37["ENST00000259915"], cdss = cdss.grch37["ENST00000259915"], 
                transcript_annotations = transcript.metadata.grch37, 
                rescale_introns = F, flanking_length = c(2500, 2500))

```
### UCSC annotation
```{r}
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb.hg19.ucsc <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons.hg19.ucsc <- exonsBy(txdb.hg19.ucsc, by = "tx", use.names = T)
cdss.hg19.ucsc <- cdsBy(txdb.hg19.ucsc, by = "tx", use.names = T)
introns.hg19.ucsc <- intronsByTranscript(txdb.hg19.ucsc, use.names = T)
transcripts.hg19.ucsc <- transcriptsBy(txdb.hg19.ucsc, by = "gene")

pou5f1.coordinates.hg19.ucsc <- transcripts.hg19.ucsc["5460"]
pou5f1.canonical.introns.hg19.ucsc <- introns.hg19.ucsc["uc003nsv.3"]
pou5f1.canonical.exons.hg19.ucsc <- exons.hg19.ucsc["uc003nsv.3"]
pou5f1.canonical.cdss.hg19.ucsc <- cdss.hg19.ucsc["uc003nsv.3"]

save(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/objects/pou5f1.coordinates.hg19.ucsc.Rdata", pou5f1.coordinates.hg19.ucsc)

save(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/objects/pou5f1.canonical.introns.hg19.ucsc.Rdata", pou5f1.canonical.introns.hg19.ucsc)

save(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/objects/pou5f1.canonical.exons.hg19.ucsc.Rdata", pou5f1.canonical.exons.hg19.ucsc)

save(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/objects/pou5f1.canonical.cdss.hg19.ucsc.Rdata", pou5f1.canonical.cdss.hg19.ucsc)


plotTranscripts(exons = exons.hg19.ucsc["uc003eef.3"], cdss = cdss.hg19.ucsc["uc003eef.3"])

canonical.selected.transcripts <- read.table(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/figure_1_canonical_transcripts_ensembl_its.txt",
                                             sep = "\t",
                                             header = T, 
                                             stringsAsFactors = F)

ncol(canonical.selected.transcripts)
canonical.selected.transcripts.metadata <- dplyr::data_frame(transcript_id = canonical.selected.transcripts$hg19.knownToEnsembl.name,
                                                             gene_name = canonical.selected.transcripts$hg19.kgXref.geneSymbol,
                                                             strand = 1)
canonical.selected.transcripts.metadata <- filter(canonical.selected.transcripts.metadata, !transcript_id == "n/a")

canonical.selected.transcripts.metadata[c(3, 6:8), "strand"] <- 0

canonical.selected.transcripts.2 <- read.table(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/figure_5_transcripts_ucsc_ids.txt",
                                               sep = "\t",
                                               header = T, 
                                               stringsAsFactors = F)

ncol(canonical.selected.transcripts)
canonical.selected.transcripts.2.metadata <- dplyr::data_frame(transcript_id = canonical.selected.transcripts.2$hg19.knownCanonical.transcript,
                                                               gene_name = canonical.selected.transcripts.2$hg19.kgXref.geneSymbol,
                                                               strand = 1)
canonical.selected.transcripts.2.metadata <- filter(canonical.selected.transcripts.2.metadata, !transcript_id == "n/a")

canonical.selected.transcripts.2.metadata[c(2), "strand"] <- 0

figure.4.bulk.RNA.seq.candidates <- read.table(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/figure_4_all_transcripts_ucsc_ids.txt",
                                               sep = "\t",
                                               header = T, 
                                               stringsAsFactors = F)



figure.4.bulk.RNA.seq.candidates.canonical <- figure.4.bulk.RNA.seq.candidates %>% filter(hg19.knownCanonical.transcript != "n/a")
figure.4.bulk.RNA.seq.candidates.canonical.metadata <-
  dplyr::data_frame(
    transcript_id = figure.4.bulk.RNA.seq.candidates.canonical$hg19.knownToEnsembl.name,
    gene_name = figure.4.bulk.RNA.seq.candidates.canonical$hg19.kgXref.geneSymbol,
    strand = ifelse(
      figure.4.bulk.RNA.seq.candidates.canonical$hg19.knownGene.strand == "+",
      1,
      0
    )
  )

figure.4.bulk.RNA.seq.candidates.canonical <- figure.4.bulk.RNA.seq.candidates %>% filter(hg19.knownCanonical.transcript != "n/a")


par(cex = 2)
plotCoverage(exons = exons.hg19.ucsc["uc003nsv.3"], cdss = cdss.hg19.ucsc["uc003nsv.3"],
             track_data = track.data, transcript_annotations = canonical.selected.transcripts.metadata,
             heights = c(0.75, 0.1), plot_fraction = 0.5,
             fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
             rescale_introns = F, 
             flanking_length = c(5000, 5000)#,
             # region_coords = c(31128965, 31142579)
)
```

### Create sample data and coverage plots - UCSC annotation
```{r}
sample.conditions <- samples.labels %>% gsub("_32F|_55F", "", .) %>% gsub("_HDFa", "", .) %>% gsub("SmithR", "Naive", .)
sample.data <- dplyr::data_frame(sample_id = samples.labels,
                                 condition = factor(sample.conditions, levels = unique(sample.conditions)),
                                 scaling_factor = 1)



sample.path <- as.data.frame(list.files(path = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/UCSC_hg19/peaks/macs2_bw", 
                                        pattern = "pileup",
                                        full.names = T,
                                        include.dirs = T))

colnames(sample.path)[1] <- "bigWig"
sample.path$bigWig <- as.character(sample.path$bigWig)

sample.path$bigW %>% gsub(".*/", "", .) %>% gsub("_treat_.*", "", .) %>% gsub("-", "_", .) %>% gsub("ATAC[_:.]", "", .) %>% gsub("_rep2", "", .) %>% gsub("\\.", "\\_", .) %>% gsub("\\_15.*", "", .) %>% gsub("Smith\\_R", "SmithR", .) %>% gsub("\\_KSR", "", .)

sample.path$sample_id <- as.character(sample.path$bigWig %>% gsub(".*/", "", .) %>% gsub("_treat_.*", "", .) %>% gsub("-", "_", .) %>% gsub("ATAC[_:.]", "", .) %>% gsub("_rep2", "", .) %>% gsub("\\.", "\\_", .) %>% gsub("\\_15.*", "", .) %>% gsub("Smith\\_R", "SmithR", .) %>% gsub("\\_KSR", "", .))

sample.data <- inner_join(sample.data, sample.path, by = "sample_id")

track.data  <- dplyr::mutate(sample.data, track_id = condition, colour_group = condition)
```
### Coverage Plots
```{r coverage plots}
for (i in canonical.selected.transcripts.metadata$transcript_id) {
  pdf(paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/", i, ".pdf", sep = ""), paper = "special",
      height = 10, width = 4)
  
  p <- plotCoverage(exons = exons.hg19.ucsc[i], cdss = cdss.hg19.ucsc[i],
                    track_data = track.data, transcript_annotations = canonical.selected.transcripts.metadata,
                    heights = c(0.75, 0.1), plot_fraction = 0.5,
                    fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                    rescale_introns = F, flanking_length = c(5000, 7500))
  print(p)
  dev.off()
}

for (i in figure.4.bulk.RNA.seq.candidates.canonical.metadata$transcript_id) {
  pdf(paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/", 
            "figure_4_bulk_RNA-seq_candidates_ATAC-seq_tracks_", i, ".pdf", sep = ""), paper = "special",
      height = 10, width = 4)
  
  p <- plotCoverage(exons = exons.hg19.ucsc[i], cdss = cdss.hg19.ucsc[i],
                    track_data = track.data, transcript_annotations = figure.4.bulk.RNA.seq.candidates.canonical.metadata,
                    heights = c(0.75, 0.1), plot_fraction = 0.5,
                    fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                    rescale_introns = F, flanking_length = c(5000, 7500))
  print(p)
  dev.off()
}
introns.coverage.plot <- plotCoverage(exons = introns.hg19.ucsc["uc003nsv.3"], cdss = introns.hg19.ucsc["uc003nsv.3"][1],
                                      track_data = track.data, transcript_annotations = canonical.selected.transcripts.metadata,
                                      heights = c(0.75, 0.1), plot_fraction = 0.5,
                                      fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                                      rescale_introns = F, connect_exons = F, flanking_length = c(500, 500)
) #,

print(introns.coverage.plot)

## Canonical transcripts ATAC-seq 2
### NANOG: chr12:7,936,214-7,949,170 (to show wider region to see the peak clearly). I have used 6000 up 1000 down TTS

pdf(paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_5_ATAC-seq_coverage_plots/atac-seq_coverage_plot_","nanog" , ".pdf", sep = ""), paper = "special",
    height = 10, width = 4)

nanog.atac.seq.coverage.p <- plotCoverage(exons = exons.hg19.ucsc["uc009zfy.1"], cdss = cdss.hg19.ucsc["uc009zfy.1"],
                                          track_data = track.data, transcript_annotations = canonical.selected.transcripts.2.metadata,
                                          heights = c(0.75, 0.1), plot_fraction = 0.5,
                                          fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                                          rescale_introns = F, flanking_length = c(6000, 1000))
print(nanog.atac.seq.coverage.p)
dev.off()
# KLF17 has two TSS?,  chr1:44,581,761-44,601,780  or chr1:44,509,786-44,599,226

pdf(paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_5_ATAC-seq_coverage_plots/atac-seq_coverage_plot_","KLF17" , ".pdf", sep = ""), paper = "special",
    height = 10, width = 4)

klf17.atac.seq.coverage.p <- plotCoverage(exons = exons.hg19.ucsc["uc001clp.3"], cdss = cdss.hg19.ucsc["uc001clp.3"],
                                          track_data = track.data, transcript_annotations = canonical.selected.transcripts.2.metadata,
                                          heights = c(0.75, 0.1), plot_fraction = 0.5,
                                          fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                                          rescale_introns = F, flanking_length = c(3000, 1000))
print(klf17.atac.seq.coverage.p)
dev.off()

# PRMD14
pdf(paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_5_ATAC-seq_coverage_plots/atac-seq_coverage_plot_","PRMD14" , ".pdf", sep = ""), paper = "special",
    height = 10, width = 4)

prmd14.atac.seq.coverage.p <- plotCoverage(exons = exons.hg19.ucsc["uc003xym.3"], cdss = cdss.hg19.ucsc["uc003xym.3"],
                                           track_data = track.data, transcript_annotations = canonical.selected.transcripts.2.metadata,
                                           heights = c(0.75, 0.1), plot_fraction = 0.5,
                                           fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                                           rescale_introns = F, flanking_length = c(2500, 2500))
print(prmd14.atac.seq.coverage.p)
dev.off()

# ZNF729: chr19:22,459,849-22,502,709

pdf(paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_5_ATAC-seq_coverage_plots/atac-seq_coverage_plot_","ZNF729" , ".pdf", sep = ""), paper = "special",
    height = 10, width = 4)

znf729.atac.seq.coverage.p <- plotCoverage(exons = exons.hg19.ucsc["uc021urs.1"], cdss = cdss.hg19.ucsc["uc021urs.1"],
                                           track_data = track.data, transcript_annotations = canonical.selected.transcripts.2.metadata,
                                           heights = c(0.75, 0.1), plot_fraction = 0.5,
                                           fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                                           rescale_introns = F, flanking_length = c(9500, 2750))
print(znf729.atac.seq.coverage.p)
dev.off()

# SUSD2: chr22:24,575,513-24,586,967

pdf(paste("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_5_ATAC-seq_coverage_plots/atac-seq_coverage_plot_","SUSD2" , ".pdf", sep = ""), paper = "special",
    height = 10, width = 4)

susd2.atac.seq.coverage.p <- plotCoverage(exons = exons.hg19.ucsc["uc002zzn.1"], cdss = cdss.hg19.ucsc["uc002zzn.1"],
                                          track_data = track.data, transcript_annotations = canonical.selected.transcripts.2.metadata,
                                          heights = c(0.75, 0.1), plot_fraction = 0.5,
                                          fill_palette = c(rep("black", 3), rep("#FF8C00", 4), rep("#0000FF", 4)),
                                          rescale_introns = F, flanking_length = c(2000, 2000))
print(susd2.atac.seq.coverage.p)
dev.off()


```




## ATAC-seq peaks annotation by stage
```{r annotation by stage}
stage.specific.peaks.annotated <- lapply(Sys.glob("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/UCSC_hg19/peaks/replicate_intersection/strong_drive/reduced_peaks/*_annotation_stats_table.txt"), function(x) read.delim(x, sep = "\t", header = T, stringsAsFactors = F))

stage.specific.peaks.annotated.names <- Sys.glob("~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/UCSC_hg19/peaks/replicate_intersection/strong_drive/reduced_peaks/*_annotation_stats_table.txt") %>% unlist() %>% gsub("/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/UCSC_hg19/peaks/replicate_intersection/strong_drive/reduced_peaks/", "", .) %>% 
  gsub("_intersect_strong_merged_reduced_annotation_stats_table.txt", "", .) %>% gsub("ATAC-", "", .) %>% gsub("-Smith", "_NM", .) %>% 
  gsub("-Primed", "_PM", .)
stage.specific.peaks.annotated.names[c(1, 6, 7)] <- c("D0_FM", "D3_FM", "D7_FM")
names(stage.specific.peaks.annotated) <- stage.specific.peaks.annotated.names

stage.specific.peaks.annotated.subset <- lapply(stage.specific.peaks.annotated, "[", c(1:2))
stage.specific.peaks.annotated.peak.no <- join_all(stage.specific.peaks.annotated.subset, by = "Annotation")

colnames(stage.specific.peaks.annotated.peak.no)[2:ncol(stage.specific.peaks.annotated.peak.no)] <- names(stage.specific.peaks.annotated.subset)

rownames(stage.specific.peaks.annotated.peak.no) <- stage.specific.peaks.annotated.peak.no$Annotation
stage.specific.peaks.annotated.peak.no.freq <- as.data.frame(sapply(stage.specific.peaks.annotated.peak.no[2:ncol(stage.specific.peaks.annotated.peak.no)], function(x) x/sum(x)))

stage.specific.peaks.annotated.peak.no.freq$annotation <- stage.specific.peaks.annotated.peak.no$Annotation

stage.specific.peaks.annotated.peak.no.freq.melted <- melt(stage.specific.peaks.annotated.peak.no.freq, by = "annotation")

stage.specific.peaks.annotated.peak.no.freq.melted$variable <- factor(stage.specific.peaks.annotated.peak.no.freq.melted$variable,
                                                                      levels = c("D0_FM", "D3_FM", "D7_FM", "D13_PM", "D21_PM", "P3_PM", "P10_PM", "D13_NM", "D21_NM", "P3_NM", "P10_NM"))

ggplot(data = stage.specific.peaks.annotated.peak.no.freq.melted, aes( x = variable, y = value, fill = annotation)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_flip() + scale_x_discrete(limits = rev(levels(stage.specific.peaks.annotated.peak.no.freq.melted$variable))) +
  scale_fill_brewer(palette = "Paired") +
  xlab("Reprogramming Stages") +
  ylab("Proportion of stage specific ATAC-seq peaks class") +
  labs(fill = "Class") +
  gplot.theme
```
## Create a Granges object and use annotatr
### All clusters
```{r granges object from clusters}

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges <-
  makeGRangesFromDataFrame(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.homer,
                           keep.extra.columns = T, gen)

builtin_annotations() %>% grep("hg19", . , value = T)
hg19.enhancers <- c("hg19_enhancers_fantom")
hg19.gral <- c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
               'hg19_genes_intronexonboundaries')
hg19.basic.genes <- c('hg19_basicgenes', 'hg19_genes_intergenic')

hg19.enhancers.annotations <- build_annotations(genome = "hg19", annotations = hg19.enhancers)
hg19.basic <- build_annotations(genome = "hg19", annotations = hg19.basic.genes)
hg19.gral <- build_annotations(genome = "hg19", annotations = hg19.gral)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated <-
  annotate_regions(
    regions = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges,
    annotations = hg19.enhancers.annotations,
    ignore.strand = T,
    quiet = FALSE
  )

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic <-
  annotate_regions(
    regions = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges,
    annotations = hg19.basic,
    ignore.strand = T,
    quiet = FALSE
  )

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.gral <-
  annotate_regions(
    regions = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges,
    annotations = hg19.gral,
    ignore.strand = T,
    quiet = FALSE
  )
print(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated)
print(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.gral)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.df <- data.frame(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.gral.df <- data.frame(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.gral)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.df <- 
  data.frame(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic)


cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.summary <-
  summarize_annotations(annotated_regions = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated,
                        quiet = TRUE)
print(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.summary)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.plot <- cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.df %>% 
  dplyr::select(PeakID, Cluster, annot.type) %>% 
  dplyr::distinct(PeakID, annot.type, .keep_all = T) %>% 
  filter(!annot.type == "hg19_genes_promoters")

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.plot.freqs <-
  cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.plot %>%
  group_by(Cluster, annot.type) %>%
  summarise(n = n()) %>%
  mutate(total = sum(n), rel.freq = n / total)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.plot.freqs$Cluster <- factor(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.plot.freqs$Cluster)

ggplot(data = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.plot.freqs, aes(x = Cluster, y = rel.freq, fill = annot.type)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_flip() + scale_x_discrete(limits = rev(levels(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.basic.plot.freqs$Cluster))) +
  scale_fill_brewer(palette = "Paired") +
  xlab("ATAC-seq clusters") +
  ylab("Proportion of Cluster specific ATAC-seq peaks class") +
  labs(fill = "Class") +
  gplot.theme
```
### Changing Clusters
```{r basic annotation changing clusters}

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges <- 
  makeGRangesFromDataFrame(cl.s.exprSet.high.cv.high.affinity.peaks.ann.homer,
                           keep.extra.columns = T)

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic <-
  annotate_regions(
    regions = cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges,
    annotations = hg19.basic,
    ignore.strand = T,
    quiet = FALSE
  )

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.df <- 
  data.frame(cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic)

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot <- cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.df %>% 
  dplyr::select(PeakID, Cluster, annot.type) %>% 
  dplyr::distinct(PeakID, annot.type, .keep_all = T)# %>% 
# filter(!annot.type == "hg19_genes_promoters")

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs <-
  cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot %>%
  group_by(Cluster, annot.type) %>%
  summarise(n = n()) %>%
  mutate(total = sum(n), rel.freq = n / total)

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs$Cluster <-
  factor(cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs$Cluster, levels = c(1, 2, 5, 4, 7, 6, 3, 8))

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs.labels <- c("1" = "C 1", "2" = "C 2", "3" = "C 7", "4" = "C 4", "5" = "C 3", "6" = "C 6", "7" = "C 5", "8" = "C 8")

high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot <- ggplot(data = cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs, aes(x = Cluster, y = rel.freq, fill = annot.type)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_flip() + scale_x_discrete(limits = rev(levels(cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs$Cluster)),
                                  labels = cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs.labels) +
  scale_fill_brewer(palette = "Paired") +
  xlab("ATAC-seq clusters") +
  ylab("Proportion of Cluster specific ATAC-seq peaks class") +
  labs(fill = "Class") +
  gplot.theme + theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 9), axis.title = element_text(size = 10), axis.text = element_text(size = 9))

ggsave(
  plot = high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot,
  file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/plots/Figure_PDFs/Figure_5E.pdf",
  height = 8,
  width = 12.5,
  units = "cm",
  dpi = 200,
  scale = 1,
  device = "pdf"
)

## Add labels for source data

extended.figure.6.a.source.data <- cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.basic.plot.freqs                           

extended.figure.6.a.source.data <-
  extended.figure.6.a.source.data %>% mutate(
    Cluster_Label = case_when(
      Cluster == "1" ~ "C 1",
      Cluster == "2" ~ "C 2",
      Cluster == "3" ~ "C 7",
      Cluster == "4" ~ "C 4",
      Cluster == "5" ~ "C 3",
      Cluster == "6" ~ "C 6",
      Cluster == "7" ~ "C 5",
      TRUE ~ "C 8"
    )
  )
extended.figure.6.a.source.data$Cluster_Label  <-
  factor(x = extended.figure.6.a.source.data$Cluster_Label,
         levels = c("C 1", "C 2", "C 3", "C 4", "C 5", "C 6", "C 7", "C 8"))

ggplot(data = extended.figure.6.a.source.data, aes(x = Cluster_Label, y = rel.freq, fill = annot.type)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_flip() + scale_x_discrete(limits = rev(levels(extended.figure.6.source.data$Cluster_Label))) +
  scale_fill_brewer(palette = "Paired") +
  xlab("ATAC-seq clusters") +
  ylab("Proportion of Cluster specific ATAC-seq peaks class") +
  labs(fill = "Class") +
  gplot.theme + theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 9), axis.title = element_text(size = 10), axis.text = element_text(size = 9))

write_excel_csv(
  x = extended.figure.6.a.source.data %>% ungroup() %>% select(Cluster_Label, rel.freq, annot.type),
  path = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/Source_Data/extended.figure.6.a.source.data.csv"
)
filter(extended.figure.6.a.source.data, Cluster == "1") %>% dplyr::select(n) %>% .$n %>% sum()
```


#### Annotations from John (FANTOM5 enhancers)
```{r annotations from John}
hg19.enhancers.john <-
  read.table(
    file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/0F5enhancer_FibEscIps.tab",
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )

hg19.enhancers.john <- hg19.enhancers.john %>% separate(geneID, c("chr", "start", "end", sep = ":")) %>% select(-":")
hg19.enhancers.john <- hg19.enhancers.john %>% separate(col = "roadmap", into = c("tx_id", "id"), sep = "\\|")

hg19.enhancers.john.granges <-
  makeGRangesFromDataFrame(hg19.enhancers.john,
                           keep.extra.columns = T)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john <-
  annotate_regions(
    regions = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges,
    annotations = hg19.enhancers.john.granges,
    ignore.strand = T,
    quiet = FALSE
  )

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.john <- 
  annotate_regions(
    regions = cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges,
    annotations = hg19.enhancers.john.granges,
    ignore.strand = T,
    quiet = FALSE
  )

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john.df <- data.frame(cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john)

cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.john.df <- as.data.frame(cl.s.exprSet.high.cv.high.affinity.peaks.ann.granges.annotated.john)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john.df$annot.type %>% table()

GenomicRanges::intersect(
  cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges,
  hg19.enhancers.john.granges,
  ignore.strand = T
)

hg19.enhancers.john.granges[queryHits(findOverlaps(hg19.enhancers.john.granges, cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges, type="any", ignore.strand = T)),] 

hg19.enhancers.john.granges[overlapsAny(hg19.enhancers.john.granges,
                                        cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges, type = "any", ignore.strand = T)]

subsetByOverlaps()

GenomicRanges::width(hg19.enhancers.john.granges)

cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john.df.freqs <-
  cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john.df %>%
  group_by(Cluster, annot.type) %>%
  summarise(n = n()) %>%
  mutate(total = sum(n), rel.freq = n / total)

ggplot(data = cl.s.exprSet.high.cv.high.affinity.stable.peaks.ann.granges.annotated.john.df.freqs, aes(x = Cluster, y = rel.freq, fill = annot.type)) +
  geom_bar(position = "stack", stat = "identity") +
  coord_flip() + scale_x_discrete(limits = rev(levels(test$Cluster))) +
  scale_fill_brewer(palette = "Paired") +
  xlab("ATAC-seq clusters") +
  ylab("Proportion of stage specific ATAC-seq peaks class") +
  labs(fill = "Class") +
  gplot.theme


```
# ATAC-seq hg19 biological TSS

```{r import peak counts}
hg19.atac.seq.tss.biol.reps <- read.table(file = "/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/UCSC_hg19/peaks/counts_over_tss/HR_ATAC-seq_hg19_counts_over_tss_10_100_table.txt",
                                          header = T,
                                          # # skip = 2,
                                          stringsAsFactors = F
                                          # colClasses = c(rep("character", 5), rep("numeric", 10))
)

hg19.atac.seq.tss.biol.reps %>% head()
```

### CPM
```{r RPKM}
rownames(hg19.atac.seq.tss.biol.reps) <- hg19.atac.seq.tss.biol.reps$Geneid


hg19.atac.seq.tss.biol.reps.cpm <- cpm(hg19.atac.seq.tss.biol.reps[, 7:ncol(hg19.atac.seq.tss.biol.reps)], log = F)

which((rowSums(hg19.atac.seq.tss.biol.reps.cpm > 10) >= 2) == TRUE) %>% length()

keep.hg19.atac.seq.tss.biol.reps.cpm<- which((rowSums(hg19.atac.seq.tss.biol.reps.cpm > 10) >= 2) == TRUE)

hg19.atac.seq.tss.biol.reps.cpm <- hg19.atac.seq.tss.biol.reps.cpm[keep.hg19.atac.seq.tss.biol.reps.cpm, ]

hg19.atac.seq.tss.biol.reps.cpm.log2 <- log2(hg19.atac.seq.tss.biol.reps.cpm + 1)
boxplot(hg19.atac.seq.tss.biol.reps.cpm.log2)

hg19.atac.seq.tss.biol.reps.cpm.log2.qnorm <- normalizeQuantiles(hg19.atac.seq.tss.biol.reps.cpm.log2)
boxplot(hg19.atac.seq.tss.biol.reps.cpm.log2.qnorm)


plotMDS(hg19.atac.seq.tss.biol.reps.cpm.log2.qnorm,
        top = nrow(hg19.atac.seq.tss.biol.reps.cpm.log2.qnorm),
        gene.selection = "common")
save(file = "~/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/objects/hg19.atac.seq.tss.biol.reps.cpm.Rdata", hg19.atac.seq.tss.biol.reps.cpm)

```

# Supplementary Tables
## Suplementary Table 9 (Fuzzy clustering of ATAC-seq signals)
```{r supplementary table 9}
supplementary.table.9 <- cl.s.exprSet.high.cv.high.affinity.peaks.ann.homer %>% 
  dplyr::mutate(Cluster_label = ifelse(Cluster == "1",                                                                                            "C 1",
                                       ifelse(
                                         Cluster == "2",
                                         "C 2",
                                         ifelse(Cluster == "3", "C 7",
                                                ifelse(
                                                  Cluster == "4", " C 4",
                                                  ifelse(Cluster == "5", " C 3",
                                                         ifelse(
                                                           Cluster == "6", "C 6",
                                                           ifelse(Cluster == "7", " C 5", "C 8")
                                                         ))
                                                ))
                                       )))

colnames(supplementary.table.9)[1:11] <- 
  c(paste("FM", "_", "D", c(0, 3, 7), sep = ""),
    paste("PM", "_", "D", c(13, 21), sep = ""),
    paste("PM", "_", "P", c(3, 10), sep = ""),
    paste("NM", "_", "D", c(13, 21), sep = ""),
    paste("NM", "_", "P", c(3, 10), sep = ""))

supplementary.table.9 %>% select(PeakID, Cluster_label, Chr, Start, End, Length, everything())

write.xlsx2(
  x = supplementary.table.9 %>% dplyr::filter(Cluster == 1) %>% dplyr::select(-c("Cluster", "Strand")) %>% select(PeakID, Cluster_label, Chr, Start, End, Length, everything()),
  file = "/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/Supplementary_Tables/Supplementary_Table_9/Supplementary_Table_9.xlsx",
  sheetName = paste("Cluster_", 1),
  row.names = F,
  append = FALSE
)

for (i in c(2, 5, 4, 7, 6, 3, 8)) {
  write.xlsx2(
    x = supplementary.table.9 %>% dplyr::filter(Cluster == i) %>% dplyr::select(-c("Cluster", "Strand")) %>% select(PeakID, Cluster_label, Chr, Start, End, Length, everything()),
    file = "/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/Supplementary_Tables/Supplementary_Table_9/Supplementary_Table_9.xlsx",
    sheetName = paste("Cluster_", i),
    row.names = F,
    append = TRUE
  )
}
```

## Supplementary Table 10 (Standardized gene expression (averaged) for genes closest to ATAC-seq cluster peaks
)
```{r supplementary table 10}
### Reordering code As reminder
# factor(s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.melt$Cluster, levels = c(1, 2, 5, 4, 7, 6, 3, 8))
# cluster.labels.rna.seq <- c("1" = "C 1 (n = 1823)", "2" = "C 2 (n = 1347)", "3" = "C 7 (n = 790)", "4" = "C 4 (n = 339)", "5" = "C 3 (n = 892)", "6" = "C 6 (n = 1027)", "7" = "C 5 (n = 1239)", "8" = "C 8 (n = 1547)")
supplementary.table.10 <-
  s.exprSet.mean.log2.rpkm.bulk.RNA.seq.filtered.low.count.genes.atac.seq.changing.cluster.tss.annotated %>% dplyr::mutate(Cluster_label = ifelse(Cluster == "1",                                                                                            "C 1",
                                                                                                                                                  ifelse(
                                                                                                                                                    Cluster == "2",
                                                                                                                                                    "C 2",
                                                                                                                                                    ifelse(Cluster == "3", "C 7",
                                                                                                                                                           ifelse(
                                                                                                                                                             Cluster == "4", " C 4",
                                                                                                                                                             ifelse(Cluster == "5", " C 3",
                                                                                                                                                                    ifelse(
                                                                                                                                                                      Cluster == "6", "C 6",
                                                                                                                                                                      ifelse(Cluster == "7", " C 5", "C 8")
                                                                                                                                                                    ))
                                                                                                                                                           ))
                                                                                                                                                  )))

supplementary.table.10 %>% colnames()

colnames(supplementary.table.10)[1:11] <- 
  c(paste("FM", "_", "D", c(0, 3, 7), sep = ""),
    paste("PM", "_", "D", c(13, 21), sep = ""),
    paste("PM", "_", "P", c(3, 10), sep = ""),
    paste("NM", "_", "D", c(13, 21), sep = ""),
    paste("NM", "_", "P", c(3, 10), sep = ""))


colnames(supplementary.table.10)[14] <- "Symbol"

# xlsx pacakge tutorial here http://www.sthda.com/english/wiki/r-xlsx-package-a-quick-start-guide-to-manipulate-excel-files-in-r#write-data-to-an-excel-file

write.xlsx2(
  x = supplementary.table.10 %>% dplyr::filter(Cluster == 1) %>% dplyr::select(-c("Cluster")) %>% select(Ensembl_id, Symbol, Cluster_label, everything()),
  file = "/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/Supplementary_Tables/Supplementary_Table_10/Supplementary_Table_10.xlsx",
  sheetName = paste("Cluster_", 1),
  row.names = F,
  append = FALSE,
  fon
)

for (i in c(2, 5, 4, 7, 6, 3, 8)) {
  write.xlsx2(
    x = supplementary.table.10 %>% dplyr::filter(Cluster == i) %>% dplyr::select(-c("Cluster")) %>% select(Ensembl_id, Symbol, Cluster_label, everything()),
    file = "/home/fernandr/projects/Polo_Group/Ethan_Liu/Ethan_Liu_ATAC-seq_Human_Reprogramming_human_samples/analysis/R/data/Supplementary_Tables/Supplementary_Table_10/Supplementary_Table_10.xlsx",
    sheetName = paste("Cluster_", i),
    row.names = F,
    append = TRUE
  )
}




```

