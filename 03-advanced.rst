Advanced UNIX session
=====================

List of Tasks:
--------------

1. How many records in the GTF file
2. Explore the 'group' column (column 9) in the GTF file
3. Get list of chromosomes (column 1)
4. Get list of features (column 3)
5. Get the number of genes mapping onto chromosomes in total
6. Get the number of protein coding genes mapping onto chromosomes
7. Get the number of protein coding genes on chromosome X and Y
8. Get the number of transcripts of protein coding genes mapping onto chromosomes
9. Get the gene with the highest number of transcripts
10. Get the gene with the highest number of exons
11. What is the total size (in Mb) of coding sequences
12. Get the longest gene

.. note:: During the afternoon session we are going to use these commands:

	.. code-block:: bash
	
		cut
		sort
		uniq
		grep
		tr
		sed
		awk

**1. How many records in the GTF file

	.. code-block:: bash
	
	cat Mus_musculus.NCBIM37.67.gtf | wc -l

**2. Explore the 'group' column (column 9) in the GTF file

	.. code-block:: bash
	
	cut -f 9 Mus_musculus.NCBIM37.67.gtf | less -S
	
**3. Get list of chromosomes (column 1)

	.. code-block:: bash
	
	cut -f 1 Mus_musculus.NCBIM37.67.gtf | sort | uniq
	
**4. Get list of features (column 3)

	.. code-block:: bash
	
	cut -f 3 Mus_musculus.NCBIM37.67.gtf | sort | uniq
	
**5. Get the number of genes mapping onto chromosomes in total

	.. code-block:: bash
	
	grep -v ^NT Mus_musculus.NCBIM37.67.gtf | cut -f 9 | cut -d ";" -f 1 | sort | uniq | wc -l
	
**6. Get the number of protein coding genes mapping onto chromosomes

	.. code-block:: bash
	
	grep -v ^NT Mus_musculus.NCBIM37.67.gtf | grep protein_coding | cut -f 9 | cut -d ";" -f 1 | sort | uniq | wc -l
	
**7. Get the number of protein coding genes on chromosome X and Y

	.. code-block:: bash
	
	grep ^[XY] Mus_musculus.NCBIM37.67.gtf | grep protein_coding | cut -f 1,9 | cut -d ';' -f 1 | sort | uniq | cut -f 1 | sort | uniq -c
	
**8. Get the number of transcripts of protein coding genes mapping onto chromosomes

	.. code-block:: bash
	
	grep -v ^NT Mus_musculus.NCBIM37.67.gtf | grep protein_coding | cut -f 9 | cut -d ";" -f 2 | sort | uniq | wc -l
	
**9. Get the gene with the highest number of transcripts

	.. code-block:: bash
	
	grep -v ^NT Mus_musculus.NCBIM37.67.gtf | grep protein_coding | cut -f 9 | cut -d " " -f 3,5,9 | tr -d '";' | sort -k1,1 | uniq | cut -d ' ' -f 1,3 | uniq -c | sed 's/^ *//' | tr ' ' "\t" | sort -nr -k1,1 | head
	
**10. Get the gene with the highest number of exons

	.. code-block:: bash
	
	grep -v ^NT Mus_musculus.NCBIM37.67.gtf | grep protein_coding | grep exon | cut -f 9 | cut -d " " -f 3,5,9 | tr -d '";' | sort | uniq -c | sed 's/^ *//g' | tr " " "\t" | sort -rn -k1,1 | head
	
**11. What is the total size (in Mb) of coding sequences

	.. code-block:: bash
	
	grep CDS Mus_musculus.NCBIM37.67.gtf | awk -F $'\t' 'BEGIN{OFS=FS;t=0}{s=$5-$4+1;t+=s}END{print t/1000000" Mb"}'
	
**12. Get the longest gene

	.. code-block:: bash
	
	grep protein_coding Mus_musculus.NCBIM37.67.gtf | grep exon | cut -f 1,4,5,9 | cut -d " " -f 1,3 | tr -d '";' | sort -k4,4 -k2,2n > exons.bed < exons.bed
	
	awk -F $'\t' 'BEGIN{ OFS=FS }{if(NR==1){ gene=$4; chrom=$1; gene_start=$2; gene_end=$3 }else{ if(gene==$4){if(gene_end<=$3){gene_end=$3}}else{ print gene,chrom,gene_start,gene_end,gene_end-gene_start; gene=$4;chrom=$1;gene_start=$2;gene_end=$3; }}}END{print gene,chrom,gene_start,gene_end,gene_end-gene_start }' | sort -rn -k5,5 | head