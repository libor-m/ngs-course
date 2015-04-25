Genomic tools session
=====================

**Explore vcftools's functionality**

.. code-block:: bash

	vcftools –-gzvcf popdata_mda.vcf.gz –-recode –-stdout | less -S

.. code-block:: bash

	vcftools –-gzvcf popdata_mda.vcf.gz –-recode –-out new_vcf

.. code-block:: bash

	vcftools –-gzvcf popdata_mda.vcf.gz –-recode –-stdout –-keep euro_samples.txt | less -S

.. code-block:: bash

	vcftools –-gzvcf popdata_mda.vcf.gz –-recode –-stdout –-keep euro_samples.txt –-chr 11 –-from-bp 22000000 –-to-bp 23000000 | less -S

.. code-block:: bash

	vcftools –-gzvcf popdata_mda.vcf.gz –-recode –-stdout –-keep euro_samples.txt | vcftools –-vcf - --recode –-stdout –-max-missing 1 –maf 0.2 | less -S

.. code-block:: bash

	vcftools –-gzvcf popdata_mda.vcf.gz –-recode –-stdout –-keep euro_samples.txt | vcftools –-vcf - --recode –-stdout –-max-missing 1 –maf 0.2 > popdata_mda_euro.vcf

.. code-block:: bash

	vcftools –-vcf popdata_mda_euro.vcf --stdout –-weir-fst-pop musculus_samps.txt –-weir-fst-pop domesticus_samps.txt | less -S

**Exercise: Population differentiation**

.. code-block:: bash

	vcftools –-gzvcf popdata_mda.vcf.gz –-recode –-stdout –-keep euro_samples.txt | vcftools –-vcf - --recode –-stdout –-max-missing 1 –maf 0.2 > popdata_mda_euro.vcf

.. code-block:: bash

	vcftools –-vcf popdata_mda_euro.vcf --stdout –-weir-fst-pop musculus_samps.txt –-weir-fst-popdomesticus_samps.txt | tail -n +2 | awk -F $'\t' 'BEGIN{OFS=FS}{ print $1,$2-1,$2,$1":"$2,$3}' > popdata_mda_euro_fst.bed

.. code-block:: bash

	bedtools makewindows -g <(grep '^2\|^11' genome.fa.fai) -w 1000000 -s 100000 -i winnum | awk '{ print $0":1000kb" }' > windows_1000kb.bed

.. code-block:: bash

	cat windows_*.bed > windows.bed

.. code-block:: bash

	## Input files for bedops need to be sorted
	sort-bed windows.bed > windows_sorted.bed
	sort-bed popdata_mda_euro_fst.bed > popdata_mda_euro_fst_sorted.bed

	bedmap --echo --mean –-count windows_sorted.bed popdata_mda_euro_fst_sorted.bed | grep -v NA | tr "|:" "\t" > windows2snps_fst.bed

.. note:: R ggplot2 commands to plot population differentiation

	.. code-block:: bash

		library(ggplot2)

		setwd("~/Data/projects/unix_workshop_data")

		fst <- read.table("windows2snps_fst.bed", header=F,sep="\t")

		names(fst) <- c("chrom", "start", "end", "win_id", "win_size", "fst", "cnt_snps")

		fst$win_size <- factor(fst$win_size, levels=c("100kb", "500kb", "1000kb"))

		qplot(fst, data=fst, geom="density",fill=I("blue")) + facet_wrap(~win_size)
	
	.. code-block:: bash	
	
		ggplot(fst, aes(y=fst, x=start, colour=win_size)) + 	geom_line() + 
			facet_wrap(~chrom, nrow=2) + 
			scale_colour_manual(name="Window size", values=c("green", "blue","red"))

		q <- quantile(subset(fst,win_size=="500kb",select="fst")[,1],prob=0.99)[[1]]

		ggplot(fst, aes(y=fst, x=start, colour=win_size)) + 	geom_line() + 
			facet_wrap(~chrom, nrow=2) + 	geom_hline(yintercept=q,colout="black") +
			scale_colour_manual(name="Window size", values=c("green", "blue","red"))
		
.. code-block:: bash

	## Use of variables: var=value
	## `` can be used to assign output of command as a variable
	q500=`grep 500kb windows2snps_fst.bed | cut -f 6 | Rscript -e 'quantile(as.numeric(readLines("stdin")),p=c(0.99))[[1]]' | cut -d " " -f 2`

	## Call variable
	echo $q500

	grep 500kb windows2snps_fst.bed | awk -v a=$q500 -F $'\t' 'BEGIN{OFS=FS}{ if($6 >= a){print $1,$2,$3} }' | bedtools merge -i stdin > signif_500kb.bed

.. code-block:: bash

	bedtools intersect –a signif.bed –b Mus_musculus.NCBIM37.67.gtf -wa -wb | grep protein_coding | cut -f 1,2,3,4,13 | cut -d ' ' -f 1,3,9 | tr -d '"";' | sort | uniq > fst2genes.tab


