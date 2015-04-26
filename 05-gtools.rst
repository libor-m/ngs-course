Genomic tools session
=====================

**Explore bedtools & bedops functionality**

- http://bedtools.readthedocs.org/en/
- http://bedops.readthedocs.org/en/

.. code-block:: bash
	
	## Prepare files (features.bed, genes.bed, my.genome)
	
	cd
	
	mkdir data/bed
	cp /data/bed/* data/bed/.
	cd data/bed/
	
	## Get parts of features that overlap
	
	bedops --intersect genes.bed features.bed
	bedtools intersect -a genes.bed -b features.bed
	
	## Merge entire features
	
	bedops --merge genes.bed features.bed
	
	cat *.bed | sortBed > features2.bed
	bedtools merge -i features2.bed
	
	## Get complement features
	
	bedops --complement genes.bed features.bed
	
	cat *.bed | sortBed > features2.bed
	bedtools complement -i features2.bed -g my.genome
	
	## Report A which overlaps B
	
	bedops --element-of 1 genes.bed features.bed
	bedtools intersect -u -a genes.bed -b features.bed
	
	## Report B which overlpas A
	
	bedops --element-of 1 features.bed genes.bed
	bedtools intersect -u -a features.bed -b genes.bed
	
	## Report A,B which overlap each other
	
	bedtools intersect -wa -wb -a genes.bed -b features.bed
	
	## What is the base coverage of features within genes?
	
	bedmap --echo --count --bases-uniq genes.bed features.bed
	coverageBed -b genes.bed -a features.bed
	
**Explore vcftools functionality**

- http://vcftools.sourceforge.net

.. code-block:: bash

	## Prepare data files
	
	cd
	mkdir data/diff

	cp /data/mus_mda/00-popdata/*.txt data/diff/.
	cp /data/mus_mda/00-popdata/popdata_mda.vcf.gz data/diff/.

	cd data/diff/
	
	## vcf file statistics - i.e. number of samples, number of SNPs

	vcftools --gzvcf popdata_mda.vcf.gz

	## Open compressed (.gz) vcf file and view it in less
	
	vcftools --gzvcf popdata_mda.vcf.gz --recode --stdout | less -S
	
	## Open compressed (.gz) vcf file and save it as a new file
	
	vcftools --gzvcf popdata_mda.vcf.gz --recode --out new_vcf
	
	## Select subset of samples

	vcftools --gzvcf popdata_mda.vcf.gz --keep euro_samps.txt --recode --stdout | less -S

	## Select subset of samples and SNPs based on physical position in genome

	vcftools --gzvcf popdata_mda.vcf.gz --chr 11 --from-bp 22000000 --to-bp 23000000 --keep euro_samps.txt --recode --stdout | less -S

	## Select subset of samples and then select SNPs with no missing data and with minor allele frequency (MAF) no less than 0.2

	vcftools --gzvcf popdata_mda.vcf.gz --keep euro_samps.txt --recode --stdout | vcftools --vcf - --max-missing 1 --maf 0.2 --recode --stdout | less -S

	vcftools --gzvcf popdata_mda.vcf.gz --keep euro_samps.txt --recode --stdout | vcftools --vcf - --max-missing 1 --maf 0.2 --recode --stdout > popdata_mda_euro.vcf

	## Calculate Fst
	
	vcftools --vcf popdata_mda_euro.vcf --weir-fst-pop musculus_samps.txt --weir-fst-pop domesticus_samps.txt --stdout | less -S
	
**Exercise: Population differentiation**

.. code-block:: bash

	vcftools --gzvcf popdata_mda.vcf.gz --keep euro_samps.txt --recode --stdout | vcftools --vcf - --max-missing 1 --maf 0.2 --recode --stdout > popdata_mda_euro.vcf

.. code-block:: bash

	vcftools --vcf popdata_mda_euro.vcf --weir-fst-pop musculus_samps.txt  --weir-fst-pop domesticus_samps.txt --stdout | tail -n +2 | awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$2-$1,$2,$1":"$2,$3}' > popdata_mda_euro_fst.bed

.. code-block:: bash

	cp /data/mus_mda/02-windows/genome.fa.fai .

	## Create windows of 1 Mb with 100 kb step
	bedtools makewindows -g <(grep '^2\|^11' genome.fa.fai) -w 1000000 -s 100000 -i winnum | awk '{print $0":1000kb"}' > windows_1000kb.bed

	## Create windows of 500 kb with 500 kb step
	bedtools makewindows -g <(grep '^2\|^11' genome.fa.fai) -w 500000 -s 50000 -i winnum | awk '{print $0":500kb"}' > windows_500kb.bed
	
	## Create windows of 100 kb with 10 kb step		
	bedtools makewindows -g <(grep '^2\|^11' genome.fa.fai) -w 100000 -s 10000 -i winnum | awk '{print $0":100kb"}' > windows_100kb.bed
	
.. code-block:: bash

	## Concatenate windows of all sizes
	cat windows_*.bed > windows.bed

.. code-block:: bash

	## Input files for bedops need to be sorted
	sort-bed windows.bed > windows_sorted.bed
	sort-bed popdata_mda_euro_fst.bed > popdata_mda_euro_fst_sorted.bed

	bedmap --echo --mean --count windows_sorted.bed popdata_mda_euro_fst_sorted.bed | grep -v NA | tr "|:" "\t" > windows2snps_fst.bed

.. note:: R ggplot2 commands to plot population differentiation

	Get to the Rstudio by typing `localhost:8787` in your web browser.

	.. code-block:: bash

		library(ggplot2)

		setwd("~/data/diff")

		fst <- read.table("windows2snps_fst.bed", header=F,sep="\t")

		names(fst) <- c("chrom", "start", "end", "win_id", "win_size", "fst", "cnt_snps")

		fst$win_size <- factor(fst$win_size, levels=c("100kb", "500kb", "1000kb"))

		qplot(fst, data=fst, geom="density",fill=I("blue")) + facet_wrap(~win_size)
	
	.. code-block:: bash	
	
		ggplot(fst, aes(y=fst, x=start, colour=win_size)) + 
			geom_line() + 
			facet_wrap(~chrom, nrow=2) + 
			scale_colour_manual(name="Window size", values=c("green", "blue","red"))

		q <- quantile(subset(fst,win_size=="500kb",select="fst")[,1],prob=0.99)[[1]]

		ggplot(fst, aes(y=fst, x=start, colour=win_size)) + 
			geom_line() + 
			facet_wrap(~chrom, nrow=2) + 
			geom_hline(yintercept=q,colout="black") +
			scale_colour_manual(name="Window size", values=c("green", "blue","red"))
		
.. code-block:: bash

	## Use of variables: var=value
	## $() can be used to assign output of command as a variable
	## do not use ` (backticks) please, they're depracated and confusing..:)
	q500=$( grep 500kb windows2snps_fst.bed | cut -f 6 | Rscript -e 'quantile(as.numeric(readLines("stdin")),probs=0.99)[[1]]' | cut -d " " -f 2 )
	
	## Call variable
	echo $q500
	
	grep 500kb windows2snps_fst.bed | awk -v a=$q500 -F $'\t' 'BEGIN{OFS=FS}{if($6 >= a){print $1,$2,$3}}' | bedtools merge -i stdin > signif_500kb.bed
	
.. code-block:: bash
	
	cp /data/mus_mda/05-fst2genes/Mus_musculus.NCBIM37.67.gtf .
	
	bedtools intersect -a signif_500kb.bed -b Mus_musculus.NCBIM37.67.gtf -wa -wb | grep protein_coding | cut -f 1,2,3,4,13 | cut -d ' ' -f 1,3,9 | tr -d '";' | sort -u > fst2genes.tab
