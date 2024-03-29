Session 2: Genomics data
========================

During this course we are going to use various genomics data stored in standard formats
that can be easily read and processed with Unix shell commands. This session provides
overview of these common data formats. Specialized tools that are optimized for work
with these data will be mentioned here as well.

Complex tasks are more easy to be conducted with specialized genomics tools rather
than trying to achieve that solely by using built-in Unix command line tools.
Not mentioning here very specific bioinformatic tasks like read mapping,
assembly and others that simply need a software developed for such purposes.

The advantage is that majority of the specialized genomics software can be run
in the command line. The most efficient work is based on combination of both
built-in Unix shell commands and specialized genomics tools. Built-in Unix shell
commands can be used to direct flow of the data in the bioinformatics pipeline,
carry basic processing of the data, and run genomics tools with specified parameters.

Common genomics tools
---------------------

Below is the list of the most widely used genomics tools that can be used for genomics data
processing based on type of task that they carry:

Read alignment data
^^^^^^^^^^^^^^^^^^^
 - samtools (https://samtools.github.io)

Variant data
^^^^^^^^^^^^
 - vcftools (https://vcftools.github.io/index.html)
 - bcftools (https://samtools.github.io/bcftools/)

Annotation data (genome arithmetics)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 - bedtools (https://bedtools.readthedocs.io/en/latest/)
 - bedops (https://github.com/bedops/bedops)

Sequence/Alignment/Tree data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 - newick-utilities (https://github.com/tjunier/newick_utils/wiki)
 - BuddySuite (https://github.com/biologyguy/BuddySuite)


Standard genomic data formats
-----------------------------


High throughput specs hub: https://samtools.github.io/hts-specs/

.. image:: _static/ngs-flows.png

FASTQ - Short reads
^^^^^^^^^^^^^^^^^^^
Sequencing instruments produce not only base calls, but usually can assign
some quality score to each called base. Fastq contains multiple sequences, and
each sequence is paired with quality scores for each base. The quality scores
are encoded in text form.

  - http://maq.sourceforge.net/fastq.shtml

SAM - Reads mapped to reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
SAM stands for Sequence Alignment/Mapping format. It includes parts of the
original reads, that were mapped to a reference genome, together with the
position where they belong to. There is an effective binary encoded
counterpart called **BAM** and even more compressed **CRAM**.

  - http://samtools.github.io/hts-specs/SAMv1.pdf
  - https://samtools.github.io/hts-specs/CRAMv3.pdf

BED and GFF - Annotations
^^^^^^^^^^^^^^^^^^^^^^^^^
Annotations are regions in given reference genome with some optional
additional information. BED is very simple and thus easy to work with for
small tasks, GFF (General Feature Format) is a comprehensive format allowing
feature nesting, arbitrary data fields for each feature and so on.

  - http://genome.ucsc.edu/FAQ/FAQformat.html#format1
  - http://www.ensembl.org/info/website/upload/gff.html

VCF - Variants in individuals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
VCF stands for Variant Call Format. Given a reference and a set of sequenced
individuals, VCF is a format to store the differences in these individuals,
compared to the reference, efficiently. There is also a binary counterpart
BCF.

  - http://samtools.github.io/hts-specs/VCFv4.2.pdf
