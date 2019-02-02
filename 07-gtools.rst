Genomic tools session
=====================

A lot of command line tools available for genomics, e.g.:

**Read alignment data:**
 * samtools (https://samtools.github.io)

**Variant data:**
 * vcftools (https://vcftools.github.io/index.html)
 * bcftools (https://samtools.github.io/bcftools/)

**Annotation data (genome arithmetics):**
 * bedtools (https://bedtools.readthedocs.io/en/latest/)
 * bedops (https://github.com/bedops/bedops)

**Sequence/Alignment/Tree data:**
 * newick-utilities (https://github.com/tjunier/newick_utils/wiki)
 * BuddySuite (https://github.com/biologyguy/BuddySuite)

Exercise
--------

Get a population differentiation calculated as Fst between *M. m. musculus*
and *M. m. domesticus* within a given sliding window and find candidate
genes within highly differentiated regions:

	1. use ``vcftools`` to filter data and calculate Fst for individual SNPs
	2. use ``bedtools makewindows`` to create sliding windows of three sizes:

		a) 100 kb + 10 kb step
		b) 500 kb + 50 kb step
		c) 1 Mb + 100 kb step

	3. calculate average Fst for each window
	4. use R-Studio and ggplot2 to plot Fst values across the genome
	5. use R or ``tabtk`` to obtain the 99th percentile and use it to obtain a set of candidate genomic regions
	6. use ``bedtools intersect`` to get a list of candidate genes

Extract genotype data for European mouse individuals and filter out
variants having more than one missing genotype and minor allele frequency 0.2
(we have already started - you should have prepared VCF file with European samples
and filtered out variants with missing genomes and low minor allele frequency).

.. code-block:: bash

    mkdir -p ~/projects/fst
    
    cd ~/projects/fst

    IN=/data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz
    SAMPLES=/data-shared/mus_mda/00-popdata/euro_samps.txt

	vcftools --gzvcf $IN \
	--keep $SAMPLES \
	--recode --stdout |
	vcftools --vcf - \
	--max-missing 1 \
	--maf 0.2 \
	--recode \
	--stdout \
	> popdata_mda_euro.vcf

Calculate Fst values for variants between *M. m. musculus*
and *M. m. domesticus* populations (populations specified in
``musculus_samps.txt`` and ``domesticus_samps.txt``):

.. code-block:: bash

    MUS=/data-shared/mus_mda/00-popdata/musculus_samps.txt
    DOM=/data-shared/mus_mda/00-popdata/domesticus_samps.txt
    IN=popdata_mda_euro.vcf 

	vcftools --vcf $IN \
	--weir-fst-pop $MUS \
	--weir-fst-pop $DOM \
	--stdout |
	tail -n +2 |
	awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$2-1,$2,$1":"$2,$3}' \
	> popdata_mda_euro_fst.bed

Make three sets of sliding windows (100 kb, 500 kb, 1 Mb)
and concatenate them into a single file:

.. code-block:: bash

    average_fst() {
        
        bedtools makewindows \
	       -g $1 \
	       -w $2 \
	       -s $3 |
        awk -v win=$2 -F $'\t' 'BEGIN{OFS=FS}{ print $0,win }' |
        bedtools intersect \
	       -a - \
	       -b $4 \
           -wa -wb |
        sort -k4,4 -k1,1 -k2,2n |
        groupBy -i - \
	       -g 4,1,2,3 \
	       -c 9 \
	       -o mean
           
    }
    
    ## Average Fst
    
    IN=popdata_mda_euro_fst.bed
    
    grep -E '^2|^11' /data-shared/mus_mda/02-windows/genome.fa.fai > genome-fst.fa.fai
    
    GENOME=genome-fst.fa.fai
    
    # 1 Mb sliding windows with 100 kb step
    
    WIN=1000000
    STEP=100000
    
    average_fst $GENOME $WIN $STEP $IN > fst_1000kb.bed
    
    # 500 kb sliding windows with 50 kb step
    
    WIN=500000
    STEP=50000
    
    average_fst $GENOME $WIN $STEP $IN > fst_500kb.bed
    
    # 100 kb sliding windows with 10 kb step
    
    WIN=100000
    STEP=10000

    average_fst $GENOME $WIN $STEP $IN > fst_100kb.bed
    
    cat fst*.bed > windows_mean_fst.tsv

Visualize the average Fst values within the sliding windows of the three sizes
between the two house mouse subspecies in `R-Studio <http://localhost:8787>`_.
Plot the distribution of the Fst values for the three window sizes and
also plot the average Fst values along the chromosomes.

.. note:: R ggplot2 commands to plot population differentiation

	.. code-block:: bash

		library(tidyverse)

		setwd("~/projects/fst")

		## Read Fst file and rename names in header
		read_tsv('windows_mean_fst.tsv', col_names=F) -> fst

		names(fst) <- c("win_size", "chrom", "start", "end", "avg_fst" )

		# Reorder levels for window size
		fst %>%
		  mutate(win_size = factor(win_size, levels=c("100kb", "500kb", "1000kb"))) ->
		  fst

		# Plot density distribution for average Fst values across windows
		ggplot(fst, aes(avg_fst)) +
			geom_density(fill=I("blue")) +
			facet_wrap(~win_size)

	.. image:: _static/fst_dist.png
			:align: center

	.. code-block:: bash

		## Plot Fst values along physical position
		ggplot(fst, aes(y=avg_fst, x=start, colour=win_size)) +
			geom_line() +
			facet_wrap(~chrom, nrow=2) +
			scale_colour_manual(name="Window size", values=c("green", "blue","red"))

		## Retrieve 99% quantiles
		fst %>%
			group_by(win_size) %>%
			summarize(p=quantile(avg_fst,probs=0.99)) -> fst_quantiles

		## Add 99% quantiles for 500kb window
		ggplot(fst, aes(y=avg_fst, x=start, colour=win_size)) +
			geom_line() +
			facet_wrap(~chrom, nrow=2) +
			geom_hline(yintercept=as.numeric(fst_quantiles[2,2]), colour="black") +
			scale_colour_manual(name="Window size", values=c("green", "blue","red"))

	.. image:: _static/fst_on_chroms.png
			:align: center

Find the 99th percentile of genome-wide distribution of Fst values
in order to guess possible outlier genome regions. 99th percentile
can be obtained running R as command line or by using ``tabtk``.
The output would be a list of windows having Fst higher
than or equal to 99% of the data.

.. code-block:: bash

	## Calculate the 99 % quantile for average Fst for 500 kb windows
	Q=$( grep '500kb' windows_mean_fst.tsv | tabtk num -c5 -q0.99 )

	## Use of variables in AWK: -v q=value

	grep 500kb windows_mean_fst.tsv |
	  awk -v q=$Q -F $'\t' 'BEGIN{OFS=FS}$5>=q{print $2,$3,$4}' |
	  sortBed |
	  bedtools merge -i stdin \
		> signif_500kb.bed

Use the mouse gene annotation file to retrieve genes within
the windows of high Fst (i.e. putative reproductive isolation loci).

.. code-block:: bash

	bedtools intersect \
		-a signif_500kb.bed \
		-b /data-shared/bed_examples/Ensembl.NCBIM37.67.bed -wa -wb | \
		cut -f4-7 | \
		tr ";" "\t" | \
		column -t | less
