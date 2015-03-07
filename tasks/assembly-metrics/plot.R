#
# contig snakes
#

library(dplyr)
library(ggplot2)

# read the data
d <- read.delim("genomes/all.tsv", col.names=c("file", "contig", "len"))

dm <- d %>% group_by(file) %>% arrange(len) %>% mutate(nr=row_number(), cumlen=cumsum(len))

# check if everything is ok
dm %>% group_by(file) %>% summarize(m=max(cumlen))

# plot the snake v1
ggplot(dm, aes(nr, cumlen, colour=file)) + geom_line(size=2)

# calculate Nx
# order by descending contig length
# calculate relative portion of the whole assembly covered by given contig
# discretize portions to 100 bins
# pick the shortest contig from each bin as N(x)
dnx <- d %>% 
  group_by(file) %>%
  arrange(desc(len)) %>%
  mutate(cumlen=cumsum(len), cumrel=cumlen/max(cumlen), nx=cut(100*cumrel, 0:100, labels=1:100)) %>%
  group_by(file, nx) %>% 
  summarize(mlen=min(len))

ggplot(dnx, aes(as.numeric(nx), mlen, colour=file)) + 
  geom_line(size=2) +
  geom_vline(xintercept=50, colour="gray") + 
  geom_vline(xintercept=90, colour="gray") +
  xlab("N(x)") + ylab("base pairs") + ggtitle("N(x) values for different assemblies")

# compare the information provided with a historam
ggplot(d, aes(len)) + geom_histogram() + facet_wrap(~file) + scale_y_log10()
