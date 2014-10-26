Important NGS formats
=====================

A selection of the most commonly used formats in NSG data processing pipelines.

FASTQ - Short reads
-------------------

Sequencing instruments produce not only base calls, but usualy can assign some quality score to each called base.
Fastq contains multiple sequences, and each sequence is paired with quality scores for each base. The quality scores
are encoded in text form.

  - http://maq.sourceforge.net/fastq.shtml

SAM - Reads mapped to reference
-------------------------------
SAM stands for Sequence Alignment/Mapping format. It includes parts of the original reads, that were mapped
to a reference genome, together with the position where they belong to. There is an effective binary encoded
counterpart called **BAM**.

  - http://samtools.github.io/hts-specs/SAMv1.pdf

BED and GFF - Annotations
-------------------------
Annotations are regions in given reference genome with some optional additional information. BED is very simple
and thus easy to work with for small tasks, GFF (General Feature Format) is a comprehensive format allowing feature
nesting, arbitrary data fields for each feature and so on.

  - http://genome.ucsc.edu/FAQ/FAQformat.html#format1
  - http://www.ensembl.org/info/website/upload/gff.html

VCF - Variants in individuals
-----------------------------
VCF stands for Variant Call Format. Given a reference and a set of sequenced individuals, 
VCF is a format to store the differences in these individuals, compared to the reference, efficiently.
There is also a binary counterpart **BCF**.

  - http://samtools.github.io/hts-specs/VCFv4.2.pdf