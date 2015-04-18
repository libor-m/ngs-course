# 
# do some basic checks on the data and 
# calculate running window statistics for variants
#
setwd("c:/work/slavici-clanek/")

library(dplyr)
library(ggplot2)
library(gtools)

# load the data and fix chromosome order

d <- read.delim("data/variant-table.tsv") %>% 
  mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev)) 


# start with some checks 
# does quality depend on read depth?
d %>% 
  filter(qual < 999) %>%
  ggplot(aes(raw_depth, qual)) + geom_point(alpha=0.05) + scale_x_log10()

# qual ~ read depth plot
d %>% 
  filter(qual < 999) %>%
  mutate(dpsum=ref_f + ref_r + alt_r + alt_f) %>%
  ggplot(aes(dpsum, qual)) + geom_point(alpha=0.1) + scale_x_log10()

# histogram of var qualities
d %>% 
  filter(qual < 999) %>%
  ggplot(aes(qual)) + geom_histogram()

# 'gene' sizes in zebra finch
d %>% 
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%
  ggplot(aes(zf_len)) + geom_histogram()

# gene sizes per chromosome
d %>%
  group_by(chrom, contig_name) %>%
  summarize(zf_max=max(zf_max), zf_min=min(zf_min)) %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%  
  ggplot(aes(chrom, zf_len)) + geom_boxplot() + coord_flip()

# max chrom sizes
d %>% 
  group_by(chrom) %>%
  summarize(chrom_len=max(zf_max)) %>%
  ggplot(aes(chrom, chrom_len)) + geom_bar(stat="identity") + coord_flip()

# mvz annotations
d %>%
  ggplot(aes(mvz_seqpart)) + geom_histogram()


# merge in fst data
dfst <- read.delim('data/lp2-var-filtered.weir.fst')
da <- d %>% left_join(dfst, c("chrom" = "CHROM", "pos" = "POS"))

# save for further reuse
colnames(da)[19] <- "fst"
# get rid of positions without fst
daf <- da %>% filter(!is.na(fst))

save(daf, file='data/daf.RData')
write.table(daf, file='data/daf.tsv', row.names=F, quote=F, sep="\t")

load('data/daf.RData')

# check if the fst outside the mapped exons differ
ggplot(daf, aes(fst, fill=is.na(zf_pos))) + geom_density(alpha=0.6, colour=NA)
# no difference

# check check distribution of fst in badly mapped contigs
ggplot(daf, aes(fst, fill=zf_max - zf_min < 5e5)) + geom_density(alpha=0.6, colour=NA)
ggplot(daf, aes(zf_max - zf_min < 5e5, fst)) + geom_boxplot() + coord_flip()
# maybe a bit lower fst values for the badly mapped exons

# plot few chromosomes
chroms <- c("chr1", "chrZ")
da %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST), chrom %in% chroms) %>%
  ggplot(aes(pos, WEIR_AND_COCKERHAM_FST)) + geom_point(colour="#bbbbbb") + facet_wrap(~chrom, ncol=1)

# do the means over 1 Mb window in zebra finch coordinates
# x is 'chrom', 'zf_pos'
# dd is the original dataset
window_average <- function(x, dd, window_half) {
  dd %>%
    filter(chrom == x[[1]]) %>%  # pick only variants on the same chromosome
    mutate(var_dist = abs(zf_pos - (x[[2]])) ) %>% # calculate distance to all variants from the current variant
    filter(var_dist < window_half) %>%
    summarize(mean_fst = mean(fst, na.rm=T)) %>% 
    .$mean_fst
}

# pick only zf positioned and real fst rows
daf <- da %>% filter(!is.na(zf_pos), !is.na(WEIR_AND_COCKERHAM_FST))
poslist <- split(daf[,c("chrom", "zf_pos")], seq_along(daf[,1]))

fst_smooth <- sapply(lpoints, window_average, daf, win_size/2)
daf$fst_smooth <-fst_smooth

# looks the calculation, apart from being hell slow,
# ended with many NAs in the results
chroms <- c("chr1", "chrZ")
daf %>%
  filter(chrom %in% chroms) %>%
  ggplot(aes(zf_pos, WEIR_AND_COCKERHAM_FST)) + 
    geom_point(colour="#bbbbbb") + 
    geom_line(aes(y=fst_smooth, group=chrom)) +
    facet_wrap(~chrom, ncol=1)

daf %>%
  filter(chrom %in% chroms) %>%
  ggplot(aes(zf_pos, fst_smooth)) + 
  geom_line() +
  facet_wrap(~chrom, ncol=1)

# move the data to python do do the window calculations
write.table(daf, file="data/vars-filtered-fst.tsv", sep="\t", row.names=F, quote=F)

daf <- read.delim('data/daf-fst.tsv') %>% mutate(chrom=chrom %>% factor(levels=chrom %>% levels %>% mixedsort %>% rev)) 
ggplot(daf %>% filter(chrom == "chrZ"), aes(zf_pos, fst_1mb)) + geom_line()

# a faster approach - an equal sampling along the chromosomes
# say 1k points along the longest chromosome
win_size <- 1e6
# required points
req_points <- 1e3
# ensure that the whole chrom is covered
stride <- as.integer(max(daf$zf_max) %>% { .  / max(. / win_size, req_points) })

# create equally spaced poins along all chromosomes, spaced by `stride`
dfpoints <- daf %>% 
  filter(chrom != "chrUn") %>%         # take only the known chromosomes
  mutate(chrom=as.character(chrom)) %>% # do not copy around the factor levels
  group_by(chrom) %>%                   # calculate chromosome sizes
  summarize(zf_size=max(zf_max)) %>%    # ..
  filter(zf_size > win_size) %>%        # pick chromosomes bigger than requested window
  split(seq_len(dim(.)[1])) %>%         # make data suitable for lapply
  lapply(function(x) 
    data.frame(chrom=x[[1]], zf_pos=seq(from=win_size/2, to=x[[2]], by=stride), stringsAsFactors=F)) %>%
  rbind_all                         # bind the dataframes together

# add the calculated values to the 'query' points
fst_smooth <- dfpoints %>% split(seq_len(dim(.)[1])) %>% sapply(window_average, daf, win_size/2)
dfst_smooth <- dfpoints %>% mutate(fst_smooth = fst_smooth)
ggplot(dfst_smooth, aes(zf_pos, fst_smooth)) + geom_line() + facet_wrap(~chrom)

bigchroms <- c("chr1", "chr1A", "chr2", "chr3", "chr4", "chrZ")
dfst_smooth %>% filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos, fst_smooth)) + geom_line() + facet_wrap(~chrom, ncol=1)

# qual - fst plot
daf %>% filter(chrom %in% bigchroms, qual < 900) %>%
  ggplot(aes(qual, fst)) + geom_point(alpha=0.1) + facet_wrap(~chrom)

# point + line plot
# does not make sens in ggplot, because there is no dual axes in ggplot, and the
# scales are different

# fst smooth on quality filtered data
dafq <- daf %>% filter(qual > 10)
fst_smooth <- dfpoints %>% split(seq_len(dim(.)[1])) %>% sapply(window_average, dafq, win_size/2)
dfst_smooth <- dfpoints %>% mutate(fst_smooth = fst_smooth)

# randomize the fst values in each chromosome
dfr <- dafq %>% 
  select(chrom, zf_pos, fst) %>%
  group_by(chrom) %>%
  mutate(fst=sample(fst))

# 
ptm <- proc.time()
rfst_smooth <- dfpoints %>% split(seq_len(dim(.)[1])) %>% sapply(window_average, dfr, win_size/2)
proc.time() - ptm
dfst_smooth <- dfpoints %>% mutate(fst_smooth = fst_smooth, fst_random = as.numeric(rfst_smooth))

dfst_smooth %>% filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
    geom_line(aes(y=fst_smooth), colour="blue") + 
    geom_line(aes(y=fst_random), colour="yellow") +
    facet_wrap(~chrom, ncol=1)


# attempts to speed things up (naive approach takes 180 s on single core)
# - range query
library(GenomicRanges)

# this has to be done only once, subsequent queries
# with bootstrapped fst can use the result
# vars is a df with chrom, zf_pos, zf_max, zf_min
find_variants <- function(vars, win_size=1e6, req_points=1e3) {
  # ensure that the whole biggest chrom is covered with windows
  stride <- as.integer(max(vars$zf_max) %>% { .  / max(. / win_size, req_points) })

  # create equally spaced poins along all chromosomes, spaced by `stride`
  dfpoints <- vars %>% 
    filter(chrom != "chrUn") %>%         # take only the known chromosomes
    mutate(chrom=as.character(chrom)) %>% # do not copy around the factor levels
    group_by(chrom) %>%                   # calculate chromosome sizes
    summarize(zf_size=max(zf_max)) %>%    # ..
    filter(zf_size > win_size) %>%        # pick chromosomes bigger than requested window
    split(seq_len(dim(.)[1])) %>%         # make data suitable for lapply
    lapply(function(x) 
      data.frame(chrom=x[[1]], zf_pos=seq(from=win_size/2, to=x[[2]], by=stride), stringsAsFactors=F)) %>%
    rbind_all                         # bind the dataframes together

  # place variants without direct mapping in zebra finch into wide range of all exons in given contig
  ir <- IRanges(start=ifelse(is.na(vars$zf_pos), vars$zf_min, vars$zf_pos), 
              end=ifelse(is.na(vars$zf_pos), vars$zf_max, vars$zf_pos + 1))

  # create interval forest for faster querying
  gr <- GRanges(seqnames=vars$chrom, ranges=ir)
  gf <- GIntervalTree(gr)

  # create whole set of query ranges
  irq <- IRanges(start=dfpoints$zf_pos - win_size / 2, end=dfpoints$zf_pos + win_size / 2)
  gq <- GRanges(seqnames=dfpoints$chrom, ranges=irq)

  # query the intervalforest
  # annotate the results with chromosome and midpoint of query window
  rv <- findOverlaps(gq, gf) %>% 
    as.data.frame %>%
    mutate(chrom=as.factor(seqnames(gf)[subjectHits]), 
           zf_pos=gq %>% ranges %>% start %>% {.[queryHits] + win_size / 2}  )
  
  rv
}

smoothed_values <- function(hits, values) {
  hits %>%
    mutate(fst=values[subjectHits]) %>%
    group_by(chrom, zf_pos) %>%
    summarize(fst_smooth=mean(fst, na.rm=T))
}

# choose top n bigest chroms
bigchroms <- daf %>% 
  filter(chrom != "chrUn") %>%
  group_by(chrom) %>%
  summarize(chrom_len=max(zf_max)) %>%
  arrange(dplyr::desc(chrom_len)) %>%
  head(n=10) %>%
  .$chrom

bigchroms <- c("chr1", "chr1A", "chr2", "chr3", "chr4", "chrZ")
fst_plot <- function(x, chroms) 
  x %>% 
  filter(chrom %in% chroms) %>%
  ggplot(aes(zf_pos, fst_smooth)) +
  geom_line() +
  facet_wrap(~chrom, ncol=1)

# fst randomized by sampling the same chromosome
rand_fst <- function(d)
  d %>% 
  select(chrom, zf_pos, fst) %>%
  group_by(chrom) %>%
  mutate(fst=sample(fst)) %>%
  .[,"fst"] %>%
  .[[1]]

ptm <- proc.time()
ovr <- find_variants(daf)
proc.time() - ptm
# 11 seconds

ovr %>% View

ptm <- proc.time()
t <- smoothed_values(ovr, daf$fst)
proc.time() - ptm
# 0.5 seconds

ptm <- proc.time()
t <- smoothed_values(ovr, rand_fst(daf))
proc.time() - ptm
# <1 second
# that is 25k in few hours

t %>% fst_plot(bigchroms)

reps <- 1:25000
lt <- sapply(reps, function(x) smoothed_values(ovr, rand_fst(daf))$fst_smooth)
save(lt, file='data/bootstraps.RData')
t$fst_boot <- apply(lt, 1, max)
write.table(t %>% select(-smooth), "data/tfst0.tsv", sep="\t", quote=F, row.names=F)

# recalculate real smooth values without badly mapped contigs
daff <- daf %>% filter(zf_max - zf_min < 5e5)
ovr <- find_variants(daff)
t2 <- smoothed_values(ovr, daff$fst)

tm <- left_join(t %>% select(-fst_smooth), t2 %>% select(chrom, zf_pos, fst_smooth))

# almost final plot, showing max bootstrapped value
# the real smoothed value and the detected 'islands'
# t %>% 
tm %>%
  filter(chrom %in% bigchroms) %>%
  ggplot(aes(zf_pos)) + 
  geom_line(aes(y=fst_boot), colour="yellow") +
  geom_line(aes(y=fst_smooth), colour="blue") + 
  geom_point(aes(y=zf_pos), y=0.2, colour="blue", size=2, data=t %>% filter(fst_smooth > fst_boot, chrom %in% bigchroms)) +
  facet_wrap(~chrom, ncol=1) +
  ylim(c(-.1, .2)) +
  ggtitle("nightingale speciation islands as mapped to zebra finch chromosomes, 25k bootstrap")

# filter out too wide contig targets
daf %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len < 5e4) %>%
  ggplot(aes(zf_len)) + geom_histogram()


# 1426 contigs are filtered out as 'too long' at 5e4
# 1005 at 1e5
# checking the zebra finch annotation, only 6 genes is longer than 500k (ensGenes)
daf %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len > 1e6) %>%
  .[,"contig_name"] %>%
  as.character %>% factor %>%
  levels %>% length

# this filters out 16k additional variants
daf %>%
  mutate(zf_len=zf_max - zf_min) %>%
  filter(zf_len > 5e4) %>%
  summarize(nvars=n())

# TODO - final filtering:
# - badly mapped contigs
# - low quality variants
# TODO - compare convolution with 1M box kernel with the fst_boot
# i guess the bootstrap value is negatively correlated with the var density
# TODO - fst according to mvz seqpart

# - subset data - query only chromX data for chromX positions
# - parallelize from the outside

