# 
# replicate some FastQC plots with ggplot2
# 
####

library(dplyr)
library(ggplot2)

setwd("c:/work/ngs-course/ngs-course-repo/tasks//fastq")

# sequence quality lines, not very informative
# first check one, then more..
d <- read.delim("quals.tsv", col.names=c("seq", "pos", "end_pos", "qual"), header=F)
d <- read.delim("quals2.tsv", col.names=c("seq", "pos", "end_pos", "base", "qual"), header=F)

sel <- levels(d$seq)[1:10]
ggplot(d %>% filter(seq %in% sel), aes(pos, qual, colour=seq, group=seq)) + geom_line()

# fastqc uses bins with varying size: 
# 1-9 by one, up to 75 by 5, up to 300 by 50, rest by 100
# the real bin sizes are a bit weird, use some nice approximation

breaks <- c(0:9, seq(14, 50, by=5), seq(59, 100, by=10), seq(100, 300, by=50), seq(400, 1000, by=100))

# create nice labels for the intervals
labs <- data.frame(l=breaks[1:length(breaks)-1], r=breaks[2:length(breaks)]) %>%
  mutate(diff=r-l, lab=ifelse(diff > 1, paste0(l+1, "-", r), as.character(r)))

# use bins and labels to discretize sequence positions
dm <- d %>% mutate(bin=cut(pos, breaks, labels=labs$lab))

# data for quality zones
quals <- data.frame(ymin=c(0, 20, 28), ymax=c(20, 28, 40), colour=c("red", "orange", "green"))

# check the quals only
ggplot(quals, aes(ymin=ymin, ymax=ymax, fill=colour)) + 
  geom_rect(alpha=0.3, xmin=-Inf, xmax=Inf) + 
  scale_fill_identity() + 
  scale_x_discrete()

# plot binned quality
ggplot(dm, aes(bin, qual)) +
  geom_boxplot(outlier.colour=NA) + 
  ylim(c(0, 45))

# replicate the whole plot with 
# - quality zones
# - colours
# - smoother
ggplot(quals) + 
  geom_rect(xmin=-Inf, xmax=Inf, aes(ymin=ymin, ymax=ymax, fill=colour), alpha=0.3) + 
  scale_fill_identity() +
  geom_boxplot(data=dm, aes(bin, qual), outlier.colour=NA, fill="yellow") +
  geom_smooth(data=dm, aes(bin, qual, group=1), colour="blue")

# reiterate, use dm as the base data, which is conceptually more correct;)
ggplot(dm) + 
  geom_rect(xmin=-Inf, xmax=Inf, data=quals, aes(ymin=ymin, ymax=ymax, fill=colour), alpha=0.3) + 
  scale_fill_identity() +
  geom_boxplot(aes(bin, qual), outlier.colour=NA, fill="yellow") +
  geom_smooth(aes(bin, qual, group=1), colour="blue")

# this is at the limit of reasonable data size, and it's only 
# 1000 sequences - how it is possible to do it for many more?
# - the statistics for each bin are calculated on the fly, 
#   as the sequences are read from the file
# - streaming mean is easy to explain - all you need is sum and number of values

#
# replicate per_base_sequence_content.png 
#
# first count distinct pairs of bin and base
t <- dm %>% mutate(baseu=toupper(base)) %>% select(baseu, bin) %>% table %>% data.frame

# plot it
ggplot(t, aes(bin, Freq, colour=baseu, group=baseu)) + geom_line()

# we're almost there, only need to normalize the numbers in each column
# that is - divide by sum 
# first check, if the work has already been done by someone else..
ggplot(t, aes(bin, Freq, fill=baseu, group=baseu)) + geom_bar(stat="identity")

# changing the position to fill automatically calculates relative abundances
ggplot(t, aes(bin, Freq, fill=baseu, group=baseu)) + geom_bar(stat="identity", position="fill")

# this actually looks better than the line chart from fastq (i'm getting lost in the lines)
# finally put the N to the last position
levs <- rev(c("A", "C", "G", "T", "N"))
t %>% 
  mutate(baseuo=factor(baseu, levels=levs)) %>% 
  ggplot(aes(bin, Freq, fill=baseuo, order=factor(baseuo, levels=rev(levs)) )) + geom_bar(stat="identity", position="fill")

# and what if we want the line chart..?
# this looks a bit like magic, but is caused by R's very own "vector recycling"
# ie sum(Freq) which is of length() 1 is recycled, until it has the length() of the other operand
# http://cran.r-project.org/web/packages/dplyr/vignettes/window-functions.html
tn <- t %>% group_by(bin) %>% mutate(Freqn=Freq / sum(Freq))
tn %>%
  mutate(baseuo=factor(baseu, levels=levs)) %>% 
  ggplot(aes(bin, Freqn, colour=baseuo, group=baseuo)) + geom_line(size=1.3)

# what is better about this plot compared to the barchart, what is worse? 
# (line: better - can see absolute proportion of base, worse - a bit messy at both ends)
# (bar: better - can easily see the disproportion between bins, worse - difficult to find the most abundant)