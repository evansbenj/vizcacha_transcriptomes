# Here's the script for Fig1

```R
setwd('/projects/submitted/Tympanoctomys/2016_WGS_kmer_analysis')

library(ggplot2)

# make a pdf
pdf("35mer_barplot_plot.pdf",w=6, h=3, version="1.4", bg="transparent")

# load the data
tymp<-read.table("AO245_jelly_count_afterquake_35mers.histo", header=FALSE)
tymp$species <- 'tymp'
oct<-read.table("AO248_jelly_count_afterquake_35mers.histo", header=FALSE)
oct$species <- 'oct'

# make a subset of these values up to 100
tymp.sub <- subset(tymp, V1 > 1 & V1 < 102)
oct.sub <- subset(oct, V1 > 1 & V1 < 102)
# Now replace the value of V2[101] with the sum of this column starting at
# that value
tymp.sub$V2[100]<-sum(as.numeric(tymp$V2[101:10001]),na.rm=TRUE)
oct.sub$V2[100]<-sum(as.numeric(oct$V2[101:10001]),na.rm=TRUE)
# combine the data
dat <- rbind(tymp.sub, oct.sub)

# make labels
dat$labelz<-""
dat$labelz[which(dat$species %in% "tymp")[1]]<-"T. barrerae"
dat$labelz[which(dat$species %in% "oct")[1]]<-"O. mimax"

# plot
ggplot(dat, aes(x=V1, y=V2, fill=species)) + 
  scale_fill_manual(values=c("royalblue1", "red")) +
  geom_bar(stat="identity",position="dodge", width=2) +
  # remove the legend
  guides(fill=FALSE) +
  # x label
  xlab("Count") +
  # y label
  ylab("Occurrence") +
  # plot them seperately (in facets)
  facet_grid(. ~ species) +
  # get rid of extra stuff
  theme_classic() +
  # add text
  geom_text(size=5, aes(x=50,y=2.0e+08,fontface="italic", label=labelz)) +
  # remove the strips on top
  theme(strip.background = element_blank(), strip.text.x = element_blank()) 

dev.off()
# DONE!


```
