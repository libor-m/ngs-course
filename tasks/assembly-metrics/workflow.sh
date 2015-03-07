# https://wiki.nbic.nl/index.php/Assembly_metrics_explained
# generate a plot of scaffold size vs contig number (sorted by size)
# or x:N(x)

# calculate lengths from fasta file
IN=22-newbler/454Isotigs.fna.filtered

<$IN mawk '
  /^>/ { 
      if(NR > 1) print name "\t" len;
      name = $1;
      len = 0;
  } 
  $0 !~ /^>/ {
      len += length($0)
  }' | 
  cut -c2- \
> ${IN%%/*}.lens

# merge lengths into single file to visualize (pick only contigs > 500)
awk '$2 > 500 {print FILENAME "\t" $1 "\t" $2}' *.lens | 
  cut -c4- | 
  sed 's/\.lens//'\
> all.tsv

# genomes do not have prefix in their name
INS=genomes/*.lens
awk '{print FILENAME "\t" $1 "\t" $2}' $INS | 
  sed 's/\.lens//'\
> all.tsv
