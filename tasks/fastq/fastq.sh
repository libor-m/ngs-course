# visualize base qualities of first 1k reads

# http://maq.sourceforge.net/fastq.shtml
# better http://en.wikipedia.org/wiki/FASTQ_format
# this is the 'translation table'
awk 'BEGIN{for(i=33;i<127;i++){printf("%c\t%d\n", i, i - 33)}}' 

# create translation script - this won't work, because of special chars in sed
awk 'BEGIN{for(i=33;i<127;i++){printf("s/\t%c/\t%d/\n", i, i - 33)}}' > q2num.sed

# pick only qualities of reads > 50 bases
<00-raw/G59B7NP01.fastq sed '1~4s/^@//;2~4d;3~4d' | paste - - | awk 'length($2) > 50' |
  # pick first 1000
  head -1000 |
  # melt into long format
  awk 'BEGIN{OFS="\t"; 
         for(i=33;i<127;i++) quals[sprintf("%c", i)] = i - 33;
       }
       { 
         l = length($2)
         for(i=1;i<=l;i++) { 
           print $1, i, l - i, quals[substr($2, i, 1)];}
       }'\
> quals.tsv

R <<EOF
library(dplyr)
library(ggplot2)

# sequence quality lines, not very informative
# first check one, then more..
d <- read.delim("quals.tsv", col.names=c("seq", "pos", "end_pos", "qual"), header=F)
sel <- levels(d$seq)[1:10]
ggplot(d %>% filter(seq %in% sel), aes(pos, qual, colour=seq, group=seq)) + geom_line()

# 
EOF

# length of fasta sequences to get histogram
<file.fa \
  awk '$0 ~ /^>/{printf("\n%s ", $1)} $0 !~ /^>/ {printf("%s",$0)}'| 
  awk 'NR > 1{print $1 "\t"  length($2)}'|head

# quick summary from scrimer minor revision
# - to be found in brno3


http://www.smashingmagazine.com/2014/03/28/design-principles-visual-perception-and-the-principles-of-gestalt/
http://www.vanseodesign.com/web-design/design-unity/
http://gastonsanchez.com/blog/archive/