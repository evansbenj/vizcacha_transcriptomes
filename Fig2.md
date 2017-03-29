# Here's the script for Fig2 (scatterplot.R):
```R
setwd('/projects/submitted/Tympanoctomys/2016_WGS_kmer_analysis/AO245_reparc_29')
library (ggplot2)
pdf("contig_length_and_coverage.pdf",w=8, h=4, version="1.4", bg="transparent")
dat = read.table("high_abundance_kmer_contigs_coverage_and_length.txt", header=TRUE)
dat$color<-"red"
dat$labelz<-""
dat$labelz[which(dat$species %in% "tympa")[1]]<-"T. barrerae"
dat$color[which(dat$species == "octomys")]<-"blue"
dat$labelz[which(dat$species %in% "octomys")[1]]<-"O. mimax"
# octomys mtdna (length above 200, so length here is 28 bp too small)
dat$color[which(dat$length == "16302")]<-"black"
dat$color[which(dat$length == "208" & dat$coverage == "3.980769")]<-"black"
# tympa mtdna
dat$color[which(dat$length == "10273")]<-"black"
dat$color[which(dat$length == "4580" & dat$coverage == "3.623363")]<-"black"
dat$color[which(dat$length == "740" & dat$coverage == "3.687838")]<-"black"
dat$color[which(dat$length == "265" & dat$coverage == "3.060377")]<-"black"
dat$color[which(dat$length == "178" & dat$coverage == "3.556180")]<-"black"

# add 28 to get the length in base pairs
dat$length <-dat$length + 28

# reorder so the mtDNA prints last
dat <- dat[order(-as.numeric(factor(dat$color))),]




d<-ggplot(dat, aes(x=length, y=coverage, colour=color, size = color, fill=factor(color))) +
  # make it clean
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  # label axis 
  labs(x=expression("Length in base pairs"), y=expression("Coverage")) +
  # remove the legend
  theme(legend.position="none") +
  # add points
  geom_point(alpha = 1, shape=21) +
  # fill the dots differently
  scale_fill_manual(values=c("gray46","royalblue1", "red1")) + 
  # make the outline of the mtDNA dots black
  scale_color_manual(values=c("black","royalblue1", "red1"))+
  # make the mtdna larger
  scale_size_manual(values=c(2,2,2)) +
  # plot in seperate facets
  facet_grid(. ~ species) +
  # add text
  geom_text(size=5, aes(x=8000,y=225,fontface="italic", label=labelz)) +
  # remove the strips on top
  theme(strip.background = element_blank(), strip.text.x = element_blank()) 
d
dev.off()


```
