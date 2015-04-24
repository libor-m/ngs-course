Quality session
===============

Many steps of genomic data processing have some associated quality value for
their results. Here we will briefly check the first and last of those. But
there is no simple way to set your quality thresholds. You have to recognize
completely bad data. But after that there is a continuum. Sometimes you just
need an idea of the underlying biology. Find some variants for further
screening. Sometimes you're trying to pinpoint particular variant causing a
disease.


Read quality
^^^^^^^^^^^^
Each read that comes out of the (now common) sequencing machines like Illumina
or Ion Torrent has a quality score assigned with each of the bases. This is not 
true for the upcoming NanoPore or almost forgotten SOLiD machines, that are reading
more bases a time.

Phred encoding
--------------
The quality in Fastq files is encoded in `Phred quality score <http://en.wikipedia.org/wiki/Phred_quality_score>`_,
a number on a logarithmic scale, similar to decibels. 

  +---------------+-----------------------+
  | Phred quality | Probability of error  |
  +---------------+-----------------------+
  |            20 | 1 in 100              |
  +---------------+-----------------------+
  |            40 | 1 in 10,000           |
  +---------------+-----------------------+
  |            60 | 1 in 1,000,000        |
  +---------------+-----------------------+

Each Phred number is in turn encoded as a single character, so there is
straightforward mapping between the bases and the quality scores. The 
most common mapping is ``Phred+33``::

  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
  | |                       |    |         |
  0.2......................26...31........41

FastQC
------
FastQC is a nice tool that you run with your data and get nice graphical 
reports on the quality. It is not completely installed in your images,
so you can try to run a tool that was just unpacked (this is how this 
tool is distributed, there is no install script - a daily bread of a 
bioinformatician;). You have to use the full path to the tool to run it::

.. code-block:: bash

   # make a project directory for the qualities
   cd
   mkdir data/quality
   cd data

   ~/sw/FastQC/fastqc -o quality --noextract fastq/*.fastq

Now transfer the ``.html`` files from the virtual machine to yours.
Open the files on your machine. You should see a report with plots
like this:

.. image:: _static/fqc_per_base_quality.png 

.. image:: _static/fqc_per_base_sequence_content.png 

Parsing Fastq and decoding Phred
--------------------------------
To understand better what is in the FastQC plots, we will try to do the same
using UNIX and ggplot. You should be able to understand the following
pipeline, at least by taking it apart with the help of head or less. A brief
description:

- ``sed`` replaces the leading '@' with an empty string in 
  the first of every four lines and deletes the third of every four lines 
  (the useless '+' line)
- ``paste`` merges every three lines
- ``awk`` selects only reads longer than 50 bases
- ``head`` takes first 1,000 sequences
- ``awk`` creates a Phred decoding table, then uses it to decode the values,
  outputs one row for each base (see 'tidy data')

.. code-block:: bash

    IN=fastq/G59B7NP01.fastq

    <$IN sed '1~4s/^@//;3~4d' |
      paste - - - |              
      awk 'length($2) > 50' |
      head -1000 |
      awk 'BEGIN{OFS="\t"; 
             for(i=33;i<127;i++) quals[sprintf("%c", i)] = i - 33;
           }
           { 
             l = length($2)
             for(i=1;i<=l;i++) { 
               print $1, i, l - i, substr($2, i, 1), quals[substr($3, i, 1)];}
           }'\
    > quality/quals.tsv

Quality by position
-------------------
The first of the FastQC plots shows a summary of base qualities
according to position in the read.

Variant quality
^^^^^^^^^^^^^^^

