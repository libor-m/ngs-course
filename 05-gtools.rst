Genomic tools session
=====================

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

.. code-block:: bash

	## Use of variables: var=value
	## `` can be used to assign output of command
	q500=`grep 500kb windows2snps_fst.bed | cut -f 6 | Rscript -e 'quantile(as.numeric(readLines("stdin")),p=c(0.99))[[1]]' | cut -d " " -f 2`

	## Call variable
	echo $q500

	grep 500kb windows2snps_fst.bed | awk -v a=$q500 -F $'\t' 'BEGIN{OFS=FS}{ if($6 >= a){print $1,$2,$3} }' | bedtools merge -i stdin > signif_500kb.bed

.. code-block:: bash

	bedtools intersect –a signif.bed –b Mus_musculus.NCBIM37.67.gtf -wa -wb | grep protein_coding | cut -f 1,2,3,4,13 | cut -d ' ' -f 1,3,9 | tr -d '"";' | sort | uniq > fst2genes.tab


