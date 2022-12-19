#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/rprojectnb/lasvchal/Jacquelyn/papers/alpha-spring2021")
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper variables
cols.voc <- c(Alpha="#e41a1c", Beta="#377eb8", Gamma="#4daf4a", 
              "Non-VoC"="#525252", Failed="#969696")
cols.int <- c(Absent="white", Subconsensus="#e41a1c", 
              Consensus="grey60", Fixed="black")
cols <- circlize::colorRamp2(c(0, 1, 2, 3), 
                             c("#de2d26", "#000000", "#bdbdbd", "#ffffff"))

## inputs ----------------------------------------------------------------------
# sample info
meta <- read.csv("sampleinfo.csv",
                 na.strings="",
                 colClasses=c(Collection.date="Date", 
                              Collection.week="Date"))

# difference matrix
difs <- read.csv("clustering/difs.csv", row.names=1) %>%
        as.matrix()

# clusters 
clst <- read.csv("clustering/clusters.csv")

## figure 1 --------------------------------------------------------------------
# rolling daily average
seq(min(meta$Collection.date)+3,
    max(meta$Collection.date)-3, 
    by=1) %>%
  lapply(function(i) {
    x <- meta %>%
         filter(Collection.date >= i-3,
                Collection.date <= i+3) %>%
         group_by(Collection.date) %>%
         summarise(Cases=n(),
                   .groups="drop")
    x <- mean(x$Cases)
    data.frame(Collection.date=i,
               Cases=x)
  }) %>%
  do.call(rbind, .) %>%
  ggplot(aes(Collection.date, Cases)) +
  geom_line(size=2) +
  labs(x=element_blank(),
       y="Case count",
       title="BU rolling 7-day average cases")
ggsave("analysis/fig1a.png", units="cm", width=11, height=8)

# weekly VoC breakdown
# first we need to annotate the VoC column
meta$Label <- meta$VoC
# fill in non-voc
meta[is.na(meta$VoC), "Label"] <- "Non-VoC"
# sequencing failed when lineage is NA
meta[is.na(meta$Lineage), "Label"] <- "Failed"
# calculate weekly cases and plot
meta %>%
  group_by(Collection.week, Label) %>%
  summarise(Cases=n(),
            .groups="drop") %>%
  mutate(Label=factor(Label, levels=rev(names(cols.voc)))) %>%
  ungroup() %>%
  ggplot(aes(Collection.week, Cases, fill=Label)) +
  geom_col(col="black") +
  scale_fill_manual(values=cols.voc) +
  labs(x=element_blank(),
       fill="VoC",
       title="BU VoC frequency by week") +
  theme(legend.position="bottom")
ggsave("analysis/fig1b.png", units="cm", width=11, height=10)

## supp. figure 1 --------------------------------------------------------------
# rolling 7-day VoC cases 
seq(min(meta$Collection.date)+3,
    max(meta$Collection.date)-3, 
    by=1) %>%
  lapply(function(i) {
    # get 7-day spanning cases
    meta %>%
      filter(Collection.date >= i-3,
            Collection.date <= i+3) %>%
      # group by date of collection and VoC label, then count
      group_by(Collection.date, Label) %>%
      summarise(Cases=n(),
              .groups="drop") %>%
      # regroup by label only and find mean
      group_by(Label) %>%
      summarise(Cases=mean(Cases),
               .groups="drop") %>%
      # add back midpoint date
      mutate(Collection.date=i)
  }) %>%
  do.call(rbind, .) %>%
  # add in missing values as zeros
  right_join(expand.grid(Collection.date=unique(.$Collection.date),
                         Label=unique(.$Label)),
             by=c("Label", "Collection.date")) %>%
  replace_na(list(Cases=0)) %>%
  # set label order via factor
  mutate(Label=factor(Label, levels=c("Failed", "Non-VoC", 
                                      "Alpha", "Beta", "Gamma"))) %>%
  ggplot(aes(Collection.date, Cases, col=Label)) +
  geom_line(size=2) +
  facet_wrap(~Label, ncol=1) +
  scale_color_manual(values=cols.voc) +
  labs(x=element_blank(),
       y="Case count",
       title="BU 7-day rolling average cases by VoC") +
  theme(legend.position="none")
ggsave("analysis/sup1a.png", units="cm", width=12, height=20)

## figure 2 --------------------------------------------------------------------
# genotype set sizes
x <- clst %>%
     group_by(Cluster.dif0) %>%
     summarise(Size=n(),
               .groups="drop") %>%
     group_by(Size) %>%
     summarise(Clusters=n(),
               .groups="drop")
# as histogram
x %>%
  ggplot(aes(Size, Clusters)) +
  geom_col(fill="black") +
  ggbreak::scale_y_break(c(50, 680)) +
  scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50, 680, 690, 700), 
                     limits=c(0, 702)) +
  labs(x="Cluster size",
         y="Number of unique genotypes")
ggsave("analysis/fig2a.png", units="cm", width=10, height=6)
# as tree map
x %>%
  ggplot(aes(area=Clusters, fill=as.factor(Size))) +
  treemapify::geom_treemap(col="black") +
  scale_fill_brewer(palette="Paired") +
  labs(fill="Cases",
       title="Unique genotype cluster sizes")
ggsave("analysis/fig2b.png", units="cm", width=8, height=6)
rm(x)

# all cases heatmap with VoC colorbar
# subset metadata to only sequenced samples
meta <- filter(meta, Barcode %in% rownames(difs))
difs <- difs[meta$Barcode, meta$Barcode]
# make VoC annotation
ann <- ComplexHeatmap::rowAnnotation(VoC=meta$VoC,
                                     col=list(VoC=cols.voc[1:3]), 
                                     border=TRUE)
png("analysis/fig2c.png", 
    units="cm", res=300, width=10, height=8)
ComplexHeatmap::Heatmap(difs, name="# difs", 
                        col=cols,
                        row_dend_width=unit(3, "cm"),
                        left_annotation=ann, 
                        border=TRUE,
                        show_row_names=FALSE,
                        show_row_dend=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
rm(ann)

# introduction and circulation "ladder"
clst %>%
  left_join(meta, by="Barcode") %>%
  group_by(Cluster.dif0, Collection.date) %>%
  summarise(Cases=n(), 
           .groups="drop") %>%
  arrange(Collection.date) %>%
  mutate(Cluster.dif0=factor(Cluster.dif0, levels=unique(Cluster.dif0)),
        Cluster.dif0=as.integer(Cluster.dif0)) %>%
  mutate(Cases=factor(Cases, levels=sort(unique(Cases)))) %>%
  ggplot(aes(Collection.date, Cluster.dif0, group=Cluster.dif0)) +
  geom_line() +
  geom_point(aes(size=Cases, fill=Cases), pch=21) +
  scale_fill_manual(values=c("15"="red"), na.value="grey60") +
  scale_size_manual(values=c("1"=0.5, "2"=1.2, "3"=2.2, 
                             "4"=3.2, "5"=4.2, "15"=6)) +
  labs(x="Collection date",
       y="Genotypes",
       title="Introductions and duration at BU")
ggsave("analysis/fig2d.png", units="cm", width=15, height=22)

## figure 3 --------------------------------------------------------------------
# Greek life cluster is 373 identical, 296 dif1, and 271 dif2
# working with dif1 for now!

# subset cluster info, metadata, and difference matrix
grk.meta <- clst %>%
            filter(Cluster.dif1==296) %>%
            left_join(meta, by="Barcode") %>%
            arrange(Collection.date, Cluster.dif0) %>%
            rownames_to_column("Sample.order") %>%
            mutate(Sample.order=as.integer(Sample.order),
                   Genotype=factor(Cluster.dif0,
                                   labels=LETTERS[1:6])) 
grk.difs <- difs[grk.meta$Barcode,
                 grk.meta$Barcode]

# figure 3A: cluster heatmap
png("analysis/fig3a.png", 
    units="cm", res=300, width=9, height=8)
ComplexHeatmap::Heatmap(grk.difs, name="# difs", 
                        col=cols,
                        row_dend_width=unit(3, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_row_dend=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()

# figure 3B: accumulation
grk.meta %>%
  ggplot(aes(Collection.date, Sample.order, 
             fill=Genotype, 
             group=Cluster.dif1)) +
  geom_line() +
  geom_point(size=2, pch=21) +
  scale_fill_brewer(palette="Paired") +
  labs(x=element_blank(),
       y="Cummulative cases")
ggsave("analysis/fig3b.png", units="cm", width=10, height=10)

# load all consensus SNVs
snvs <- grk.meta$Barcode %>%
        lapply(function(i) {
          paste0("data/", i, ".csv") %>%
            read.csv() %>% 
            mutate(Barcode=i) %>%
            select(Barcode, NT.ID, AA.ID, Frequency)
        }) %>%
        do.call(rbind, .)
# remove consensus SNVs found in all samples
# and annotate as fixed, consensus, subconsensus, or not present
x <- snvs %>%
     filter(Frequency > 0.5) %>%
     group_by(NT.ID) %>%
     summarise(Cases=n(),
               .groups="drop") %>%
     filter(Cases < dim(grk.meta)[1]) %>%
     select(NT.ID) %>%
     unlist()
snvs <- snvs %>%
        filter(NT.ID %in% x) %>%
        mutate(Intensity=cut(Frequency, 
                             breaks=c(0, 0.05, 0.5, 0.95, 1), 
                             labels=c("Absent", "Subconsensus", 
                                      "Consensus", "Fixed")))
rm(x)

# figure 3C: 488-G-A frequency plot
snvs %>%
  filter(NT.ID=="488-G-A") %>% 
  ggplot(aes(Intensity, Frequency, fill=Intensity)) +
  geom_jitter(pch=21, width=0.2, height=0, size=3, alpha=0.8) +
  scale_fill_manual(values=cols.int) +
  ylim(0, 1) +
  labs(x="NSP1 D75N",
       y="Cummulative cases") +
  theme(legend.position="none")
ggsave("analysis/fig3c.png", units="cm", width=9, height=10)

# figure 3D: 488-G-A grows out of index case
snvs %>%
  filter(NT.ID=="488-G-A") %>%
  left_join(grk.meta, by="Barcode") %>% 
  ggplot(aes(Collection.date, Sample.order, fill=Intensity, 
             group=Cluster.dif1)) +
  geom_line() +
  geom_point(size=2, pch=21) +
  scale_fill_manual(values=cols.int) +
  labs(x=element_blank(),
       y="Cummulative cases",
       fill="NSP1 D75N")
ggsave("analysis/fig3d.png", units="cm", width=11.5, height=10)
rm(grk.meta, grk.difs, snvs)

## figure 4 --------------------------------------------------------------------
# Greek life cluster is 238 identical, 202 dif1, and 191 dif2
# working with dif1 from I238 for now!
x <- difs[, "bu0033"]
x <- names(x[x <= 1])

# subset cluster info, metadata, and difference matrix
grk.meta <- clst %>%
            filter(Barcode %in% x) %>%
            left_join(meta, by="Barcode") %>%
            arrange(Collection.date, Cluster.dif0) %>%
            rownames_to_column("Sample.order") %>%
            mutate(Sample.order=as.integer(Sample.order),
                   Genotype=factor(Cluster.dif0,
                                   labels=LETTERS[1:10])) 
grk.difs <- difs[grk.meta$Barcode,
                 grk.meta$Barcode]
rm(x)

# figure 4A: cluster heatmap
png("analysis/fig4a.png", 
    units="cm", res=300, width=9, height=8)
ComplexHeatmap::Heatmap(grk.difs, name="# difs", 
                        col=cols,
                        row_dend_width=unit(3, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_row_dend=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()

# figure 4B: accumulation
grk.meta %>%
  ggplot(aes(Collection.date, Sample.order, 
             fill=Genotype, 
             group=Cluster.dif1)) +
  geom_line() +
  geom_point(size=2, pch=21) +
  scale_fill_brewer(palette="Paired") +
  labs(x=element_blank(),
       y="Cummulative cases")
ggsave("analysis/fig4b.png", units="cm", width=10, height=10)

# load all consensus SNVs
snvs <- grk.meta$Barcode %>%
        lapply(function(i) {
          paste0("data/", i, ".csv") %>%
            read.csv() %>% 
            mutate(Barcode=i) %>%
            select(Barcode, NT.ID, AA.ID, Frequency)
        }) %>%
        do.call(rbind, .)
# remove consensus SNVs found in all samples
# and annotate as fixed, consensus, subconsensus, or not present
x <- snvs %>%
     filter(Frequency > 0.5) %>%
     group_by(NT.ID) %>%
     summarise(Cases=n(),
               .groups="drop") %>%
     filter(Cases < dim(grk.meta)[1]) %>%
     select(NT.ID) %>%
     unlist()
snvs <- snvs %>%
        filter(NT.ID %in% x) %>%
        mutate(Intensity=cut(Frequency, 
                             breaks=c(0, 0.05, 0.5, 0.95, 1), 
                             labels=c("Absent", "Subconsensus", 
                                      "Consensus", "Fixed")))
rm(x)

# 23194-T-C: all but 1 at consensus
# 5301-C-T: all but 1 at consensus
# 28378-G-T: all but 2 at consensus
# figure 4C: frequency plot
snvs %>%
  filter(NT.ID=="28378-G-T") %>% 
  right_join(grk.meta, by="Barcode") %>%
  replace_na(list(Frequency=0,
                  Intensity="Absent")) %>%
  ggplot(aes(Intensity, Frequency, fill=Intensity)) +
  geom_jitter(pch=21, width=0.2, height=0, size=3, alpha=0.8) +
  scale_fill_manual(values=cols.int) +
  ylim(0, 1) +
  labs(x="INSERT AA ID",
       y="Cummulative cases") +
  theme(legend.position="none")
ggsave("analysis/fig4c.png", units="cm", width=9, height=10)

# figure 4D: grows out of index case
snvs %>%
  filter(NT.ID=="28378-G-T") %>% 
  right_join(grk.meta, by="Barcode") %>%
  replace_na(list(Frequency=0,
                  Intensity="Absent")) %>%
  ggplot(aes(Collection.date, Sample.order, fill=Intensity, 
             group=Cluster.dif1)) +
  geom_line() +
  geom_point(size=2, pch=21) +
  scale_fill_manual(values=cols.int) +
  labs(x=element_blank(),
       y="Cummulative cases",
       fill="INSERT AA ID")
ggsave("analysis/fig4d.png", units="cm", width=11.5, height=10)
rm(grk.meta, grk.difs, snvs)
