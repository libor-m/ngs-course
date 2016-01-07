# combine variant and annotation info into one convenient table
# magic numbers (Fst..) will be joined in R
python dump_vars.py data/lp2.sorted.gff3 data/lp2-var-filtered.vcf.gz > data/variant-table.tsv

python dump_vars.py data/lp2.sorted.gff3.gz data/lp2-var-filtered.vcf.gz | pv -l > data/variant-table.tsv

# find the longest gene models in zebra finch ensemble annotation
<taeGut1/annot/ensGenes.sorted.bed.gz zcat|mawk '{print $3 - $2;}'|sort -rn|head -20
# ensGenes max is <1M, 6 models is >500k
# refSeq has one 1036509, the second is 495217

ANNOTS=~/brno3/genomes/taeGut1/annot
ANN="$ANNOTS/ensGenes.bed.gz $ANNOTS/refSeqGenes.bed.gz"
zcat $ANN | 
  bedtools sort | 
  bedtools intersect -wa -sorted -a - -b <( <80-islands/islands.bed tr -d "\r" | bedtools sort ) \
> 80-islands/genes-in-islands.bed

<80-islands/genes-in-islands.bed cut -f4 | grep ENS > 80-islands/genes-in-islands.ensg
<80-islands/genes-in-islands.bed cut -f4 | grep -v ENS > 80-islands/genes-in-islands.refseq


# pick song related genes
# the paper states that blue, orange and dark green modules are 'song modules'
<song/neuron_10977_mmc3.txt awk '$6 == "blue" || $6 = "orange" || $6 == "darkgreen"' > song/song-modules.tsv
# 689 probes

# go search it in the genome, because many of the genes were not annotated back then..
module add blat-suite-34

# extract the probe sequences
<81-song/song-modules.tsv awk '{print ">seq" NR; print $3;}' > 81-song/song-modules.fa
GENOMES=~/brno3/genomes
GENOME=$GENOMES/taeGut1/taeGut1.2bit

# search! 
blat -fastMap -dots=10 $GENOME 81-song/song-modules.fa 81-song/probe-blat.psl

# check the number of hits for fastMap
<81-song/probe-blat.psl tail -n +6 | cut -f10|sort|uniq|wc -l
# 476 - not good for 689 probes

# search again, consider -fine and -noHead
blat -dots=10 $GENOME 81-song/song-modules.fa 81-song/probe-blat-full.psl
# 634 matches looks better

# dxy calculation
python dxy.py data/lp2-var-filtered.vcf.gz populations > data/dxy.tsv