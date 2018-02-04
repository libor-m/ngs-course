Genomic tools session
=====================

Genome feature arithmetics & summary
------------------------------------

**Explore bedtools & bedops functionality**

- http://bedtools.readthedocs.io/
- https://bedops.readthedocs.io/

1. Count the number of merged open chromatin regions overlapping with genes or are within 1000 bp window on each side of a gene

In the first exercise we will work with open chromatin regions
based on DNaseI hypersensitive sites in file ``encode.bed`` obtained
from ENCODE database. We would like to parse and count those open
chromatin regions which overlap with known genes retrieved from Ensembl
database.

.. code-block:: bash

	# Explore the Ensembl.NCBIM37.67.bed file
	less /data-shared/bed_examples/Ensembl.NCBIM37.67.bed

	# Explore the encode.bed file
	less /data-shared/bed_examples/encode.bed

	# Count the number of open chromatin regions overlapping with genes:
	bedtools intersect \
	-a <( sortBed -i encode.bed ) \
	-b <( sortBed -i /data-shared/bed_examples/Ensembl.NCBIM37.67.bed ) |
	wc -l

	## Count the number of open chromatin regions overlapping with genes and within 1000 bp window on each side
	bedtools window -w 1000 \
	-a <( sortBed -i encode.bed ) \
	-b <( sortBed -i /data-shared/bed_examples/Ensembl.NCBIM37.67.bed ) |
	wc -l

2. Make two sets of sliding windows across mouse genome (1 Mb, 5 Mb)
with the step size 0.2 by the size of the window and obtain gene density
within these sliding windows. To speed up the process we focus only on chromosome X.

.. code-block:: bash

	# Explore fasta index file
	less /data-shared/bed_examples/genome.fa.fai

	# Make 1Mb sliding windows (step 200kb)
	bedtools makewindows \
	-g <( grep '^X' /data-shared/bed_examples/genome.fa.fai ) \
	-w 1000000 \
	-s 200000 \
	-i winnum \
	> windows_1mb.bed

	# Make 5Mb sliding windows (step 1Mb)
	bedtools makewindows \
	-g <( grep '^X' /data-shared/bed_examples/genome.fa.fai ) \
	-w 5000000 \
	-s 1000000 \
	-i winnum \
	> windows_5mb.bed

	# Obtain densities of genes within individual windows
	bedtools coverage \
	-a windows_1mb.bed \
	-b <( sortBed -i /data-shared/bed_examples/Ensembl.NCBIM37.67.bed ) \
	> gdens_windows_1mb.tab

	bedtools coverage \
	-a windows_5mb.bed \
	-b <( sortBed -i /data-shared/bed_examples/Ensembl.NCBIM37.67.bed ) \
	> gdens_windows_5mb.tab

The gene density can be visualized in R-Studio.

VCFtools
--------

**Explore vcftools functionality**

- http://vcftools.sourceforge.net

Prepare working directory ``projects/fst``:

.. code-block:: bash

	cd
	mkdir projects/fst
	cd projects/fst

Obtaining the basic file statistics (number of variants & number of samples):

.. code-block:: bash

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz

Viewing and printing out the content of the VCF file:

.. code-block:: bash

	# To print out the content of the VCF file

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz \
	--recode \
	--out new_vcf

	# To view the content directly

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz \
	--recode \
	--stdout | less -S

Basic data filtering - use of appropriate flags:

.. code-block:: bash

	--keep ind.txt # Keep these individuals
	--remove ind.txt # Remove these individuals
	--snps snps.txt # Keep these SNPs
	--snps snps.txt –-exclude # Remove these SNPs

To select a subset of samples:

.. code-block:: bash

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz \
	--keep /data-shared/mus_mda/00-popdata/euro_samps.txt \
	--recode \
	--stdout |
	less -S

Select subset of samples and SNPs based on physical position in genome:

.. code-block:: bash

	# Flags you can use:
	--chr 11 # Keep just this chromosome
	--not-chr 11 # Remove this chromosome
	--not-chr 11 –not-chr 2 # Remove these two chromosomes
	--from-bp 20000000 # Keep SNPs from this position
	--to-bp 22000000 # Keep SNPs to this position
	--bed keep.bed # Keep only SNPs overlapping with locations listed in a file
	--exclude-bed remove.bed # The opposite of the previous

.. code-block:: bash

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz \
	--chr 11 \
	--from-bp 22000000 \
	--to-bp 23000000 \
	--keep /data-shared/mus_mda/00-popdata/euro_samps.txt \
	--recode \
	--stdout |
	less -S

Select subset of samples and then select SNPs with no missing data
and with minor allele frequency (MAF) no less than 0.2:

.. code-block:: bash

	# Flags you can use:
	--maf 0.2 # Keep just variants with Minor Allele Freq higher than 0.2
	--hwe 0.05 # Keep just variants which do not deviate from HW equilibrium (p-value = 0.05)
	--max-missing (0-1) # Remove SNPs with given proportion of missing data (0 = allowed completely missing, 1 = no missing data allowed)
	--minQ 20 # Minimal quality allowed (Phred score)

.. code-block:: bash

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz \
	--keep /data-shared/mus_mda/00-popdata/euro_samps.txt \
	--recode \
	--stdout |
	vcftools \
	--vcf - \
	--max-missing 1 \
	--maf 0.2 \
	--recode \
	--stdout |
	less -S

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz \
	--keep /data-shared/mus_mda/00-popdata/euro_samps.txt \
	--recode \
	--stdout |
	vcftools --vcf - \
	--max-missing 1 \
	--maf 0.2 \
	--recode \
	--stdout \
	> popdata_mda_euro.vcf

Use the newly created ``popdata_mda_euro.vcf`` representing variants
only for a subset of individuals and variants to calculate Fst index.
In order for vcftools to calculate Fst index the populations
have to be specified in the output - each one with a separate file
(``--weir-fst-pop pop1.txt`` and ``--weir-fst-pop pop2.txt``).

.. code-block:: bash

	# Flags you can use:
	--site-pi # Calculates per-site nucleotide diversity (π)
	--window-pi 1000000 --window-pi-step 250000 # Calculates per-site nucleotide diversity for windows of 1Mb with 250Kb step
	--weir-fst-pop pop1.txt --weir-fst-pop pop2.txt # Calculates Weir & Cockerham's Fst
	--fst-window-size 1000000 –-fst-window-step 250000 # Calculates Fst for windows of 1Mb with 250Kb step

.. code-block:: bash

	vcftools --vcf popdata_mda_euro.vcf \
	--weir-fst-pop /data-shared/mus_mda/00-popdata/musculus_samps.txt \
	--weir-fst-pop /data-shared/mus_mda/00-popdata/domesticus_samps.txt \
	--stdout |
	less -S

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

	cd ~/projects/fst

	vcftools --gzvcf /data-shared/mus_mda/00-popdata/popdata_mda.vcf.gz \
	--keep /data-shared/mus_mda/00-popdata/euro_samps.txt \
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

	vcftools --vcf popdata_mda_euro.vcf \
	--weir-fst-pop /data-shared/mus_mda/00-popdata/musculus_samps.txt   \
	--weir-fst-pop /data-shared/mus_mda/00-popdata/domesticus_samps.txt \
	--stdout |
	tail -n +2 |
	awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$2-1,$2,$1":"$2,$3}' \
	> popdata_mda_euro_fst.bed

Make the three sets of sliding windows (100 kb, 500 kb, 1 Mb)
and concatenate them into a single file:

.. code-block:: bash

	## Create windows of 1 Mb with 100 kb step
	bedtools makewindows \
	-g <(grep -E '^2|^11' /data-shared/mus_mda/02-windows/genome.fa.fai) \
	-w 1000000 \
	-s 100000 |
	awk -F $'\t' 'BEGIN{OFS=FS}{print $0,"1000kb"}' \
	> windows_1000kb.bed

	## Create windows of 500 kb with 500 kb step
	bedtools makewindows \
	-g <(grep -E '^2|^11' /data-shared/mus_mda/02-windows/genome.fa.fai) \
	-w 500000 \
	-s 50000 |
	awk -F $'\t' 'BEGIN{OFS=FS}{print $0,"500kb"}' \
	> windows_500kb.bed

	## Create windows of 100 kb with 10 kb step
	bedtools makewindows \
	-g <(grep -E '^2|^11' /data-shared/mus_mda/02-windows/genome.fa.fai) \
	-w 100000 \
	-s 10000 | \
	awk -F $'\t' 'BEGIN{OFS=FS}{print $0,"100kb"}' \
	> windows_100kb.bed

	## Concatenate windows of all sizes
	cat windows_*.bed > windows.bed

Calculate average Fst within the sliding windows:

.. code-block:: bash

	## Input files for bedtools groupby need to be sorted

	# Join Fst values and the 'windows.bed' file
	bedtools intersect \
	  -a <( sortBed -i windows.bed ) \
	  -b <( sortBed -i popdata_mda_euro_fst.bed ) -wa -wb \
	> windows_fst.tsv

	# Run bedtools groupby command to obtain average values of Fst
	sort -k4,4 -k1,1 -k2,2n windows_fst.tsv |
	groupBy -i - \
	-g 4,1,2,3 \
	-c 9 \
	-o mean > windows_mean_fst.tsv

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

	## Use of variables in AWK: -v q=value

	grep 500kb windows_mean_fst.tsv |
	  awk -v q=0.9166656 -F $'\t' 'BEGIN{OFS=FS}$5>=q{print $2,$3,$4}' |
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
