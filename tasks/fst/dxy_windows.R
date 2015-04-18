# 
# look at dxy compared to Fst
#

setwd("c:/work/slavici-clanek/")

# granges must be loaded first beacuse it masks dplyr functions otherwise
# with some useless variants
library(GenomicRanges)

library(dplyr)
library(ggplot2)
library(gtools)

sortchrom <- function(df) df %>% mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev))

# load the data and fix chromosome order
d <- read.delim("data/variant-table.tsv") %>% sortchrom
   
# load the dxy data
dxy <- read.delim('data/dxy.tsv', col.names=c("chrom", "pos", "ndiff", "ncomp", "dxy")) %>%
  filter(!is.na(ndiff))

# join on chrom, pos
da <- d %>% inner_join(dxy)

# pick 10 biggest chroms
bigchroms <- da %>% 
  filter(chrom != "chrUn") %>%
  group_by(chrom) %>%
  summarize(chrom_len=max(pos)) %>%
  arrange(dplyr::desc(chrom_len)) %>%
  head(n=10) %>%
  .$chrom

# plot the variance ~ ncomparisons
dxy %>% ggplot(aes(ncomp, dxy)) + geom_point(alpha=0.1)
dxy %>% ggplot(aes(ncomp)) + geom_histogram()

# dot plot along chromosomes
# pick only reasonable values of ncomp based on previous plot
da %>% 
  filter(chrom %in% bigchroms, ncomp > 50) %>%
  ggplot(aes(pos, dxy)) + 
    geom_point(colour="#cccccc") + 
    facet_wrap(~chrom, ncol=1)

# use windows.R::find_variants to assign variants to sliding windows
daf_dxy <- da %>% filter(ncomp > 30)
ovr_dxy <- daf_dxy %>% find_variants

# calculate smoothed dxy
# improved since windows.R
smoothed_values <- function(hits, values) {
  hits %>%
    mutate(fst=values[subjectHits]) %>%
    group_by(chrom, zf_pos) %>%
    summarize(smooth=mean(fst, na.rm=T), nvars=n())
}

tdxy <- smoothed_values(ovr_dxy, daf_dxy$dxy) %>% dplyr::rename(dxy_smooth=smooth)
tdxy %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, dxy_smooth)) + 
    geom_line(colour="green") +
    facet_wrap(~chrom, ncol=1)
  
# get fst smoothers
dfst <- read.delim('data/lp2-var-filtered.weir.fst', col.names=c("chrom", "pos", "fst"))

daf_fst <- d %>% left_join(dfst) %>% filter(qual > 10)
ovr_fst <- daf_fst %>% find_variants
tfst <- smoothed_values(ovr_fst, daf_fst$fst) %>% dplyr::rename(fst_smooth=smooth)

# plot smooth fst and dxy
# looks i'm loosing the resolutino of fst by filtering the vars, go the other way around..
tdxy %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=dxy_smooth), colour="green") +
  geom_line(aes(y=fst_smooth), data=tfst %>% filter(chrom %in% bigchroms), colour="blue") +
  facet_wrap(~chrom, ncol=1)

tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, dxy_smooth)) + 
  geom_line(colour="green") +
  ylim(c(0, .0007)) +
  facet_wrap(~chrom, ncol=1) +
  ggtitle("per site dxy scaled by containing contig length")

ggsave('results/dxy_scaled.pdf', width=20, height=16)

# weight dxy values by size of contig thei're in
# -> few big spikes probably due to missing data, let's check it
tdxy_scaled <- smoothed_values(ovr_dxy, daf_dxy$dxy / daf_dxy$contig_size) %>% dplyr::rename(dxy_smooth=smooth)
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=dxy_smooth), colour="green") +
  geom_line(aes(y=nvars), colour="#bbbbbb") +
  geom_line(aes(y=fst_smooth), data=tfst %>% filter(chrom %in% bigchroms), colour="blue") +  
  facet_wrap(~chrom, ncol=1) +
  scale_y_log10() + 
  ggtitle("response of fst(blue) and dxy (green) to variant site density (gray)")

ggsave('results/metrics_var_density.pdf', width=20, height=16)

# try to do rescaled plots

# rescale numeric vector into (0, 1) interval
# clip everything outside the range 
rescale <- function(vec, lims=range(vec), clip=c(0, 1)) {
  # find the coeficients of transforming linear equation
  # that maps the lims range to (0, 1)
  slope <- (1 - 0) / (lims[2] - lims[1])
  intercept <- - slope * lims[1]
  
  xformed <- slope * vec + intercept
  
  # do the clipping
  xformed[xformed < 0] <- clip[1]
  xformed[xformed > 1] <- clip[2]
  
  xformed
}

# use brewer pallette for the overlapping lines
# to get same 'visual intensity'
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=rescale(dxy_smooth, lims=c(0, .0007)), colour="Dxy")) +
  geom_line(aes(y=rescale(fst_smooth), colour="Fst"), data=tfst %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(nvars), colour="nvars")) +
  facet_wrap(~chrom, ncol=1) +
  ggtitle("response fst and dxy to variant site density, rescaled to 1") +
  ylab("rescaled values") +
  xlab("zebra finch chromosome position") +
  scale_colour_brewer(type="qual", palette="Set1")

ggsave('results/metrics_var_density-scaled.pdf', width=20, height=16)


# check if the fst data is ok, find good ranges in histograms
# 
ggplot(tboot, aes(fst_boot)) + geom_histogram()
library(tidyr)
tt <- gather(tboot, type, fst, smooth, fst_boot)
tt %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, fst, colour=type)) + 
  geom_line() + 
  geom_point(y=0.2, colour="blue", shape=15, size=2, data=tboot %>% filter(smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1)

ggplot(tt, aes(fst, fill=type)) + geom_histogram()

tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>% 
  ggplot(aes(nvars)) + geom_histogram()

# add bootstrap for fst
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=rescale(dxy_smooth, lims=c(0, .0007)), colour="Dxy")) +
  geom_line(aes(y=rescale(smooth, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(fst_boot, lims=c(0, .25)), colour="Fst_boot"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(nvars), colour="nvars")) +
  geom_point(y=1, colour="blue", shape=15, size=2, data=tboot %>% filter(smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  ggtitle("response fst and dxy to variant site density, rescaled to 1") +
  ylab("rescaled values") +
  xlab("zebra finch chromosome position") +
  scale_colour_brewer(type="qual", palette="Set1")

ggsave('results/metrics_var_density-scaled.pdf', width=20, height=16)

# check if shifted tracks is more legible
# add means per valua and chromosome
tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=rescale(dxy_smooth, lims=c(0, .0007)) + 1, colour="Dxy")) +
  geom_line(aes(y=rescale(smooth, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(fst_boot, lims=c(0, .25)), colour="Fst_boot"), data=tboot %>% filter(chrom %in% bigchroms)) +  
  geom_line(aes(y=rescale(nvars, lims=c(0, 700)) + 2, colour="nvars")) +
  geom_point(y=0, colour="blue", shape=15, size=2, data=tboot %>% filter(smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, .0007)) + 1, colour="Dxy"), data=tdxy_scaled %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=mean(dxy_smooth))) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, .25)), colour="Fst"), data=tboot %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=mean(smooth))) +
  geom_hline(aes(yintercept=rescale(cmean, lims=c(0, 700)) + 2, colour="nvars"), data=tdxy_scaled %>% filter(chrom %in% bigchroms) %>% group_by(chrom) %>% summarize(cmean=median(nvars))) +
  ggtitle("response fst and dxy to variant site density, rescaled to 1") +
  ylab("rescaled values") +
  xlab("zebra finch chromosome position") +
  scale_colour_brewer(type="qual", palette="Set1") +
  theme(legend.title=element_blank())

ggsave('results/metrics_var_density-scaled-shift.pdf', width=20, height=16)


tdxy_scaled %>% 
  filter(chrom %in% bigchroms) %>% 
  group_by(chrom) %>% 
  summarize(cmean=mean(dxy_smooth) %>% rescale(c(0, .0007))) %>%
  ggplot(aes(chrom, cmean)) + geom_point() + geom_hline(aes(yintercept=0.3))

# TODO: center tracks around the mean?

# questions
# - how to deal with uneven sample coverage in contigs and across contigs
# - 