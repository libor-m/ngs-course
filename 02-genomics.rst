Unix - Advanced I
=================

This session focuses on plain text file data extraction/modification
using build-in Unix tools.

Pattern search & regular expressions
------------------------------------

``grep`` useful tool to search pattern using regular expressions.

.. block-code:: bash

  # ^A         # match A at the beginning of line
  # A$         # match A at the end of line
  # [0-9]      # match numerical characters
  # [A-Z]      # match alphabetical characters
  # [ATGC]     # match A or T or G or C
  # .          # match any character
  # A*         # match A letter 0 or more times
  # A\{2\}     # match A letter exactly 2 times
  # A\{1,\}    # match A letter 1 or more times
  # A+         # match A letter 1 or more times (extended regular expressions)
  # A\{1,3\}   # match A letter at least 1 times but no more than 3 times
  # AATT\|TTAA # match AATT or TTAA
  # \s         # match whitespace (also TAB)

*Use mouse annotation file (GTF)*

.. block-code:: bash

  cd ~
  sudo cp /data/mus_mda/05-fst2genes/Mus_musculus.NCBIM37.67.gtf.gz ~/data/.
  gunzip data/Mus_musculus.NCBIM37.67.gtf.gz
  less -S data/Mus_musculus.NCBIM37.67.gtf

1. Count the number of records on the chromosome X

.. block-code:: bash

  < data/Mus_musculus.NCBIM37.67.gtf grep '^X' | wc -l

2. Count the number of records on chromosome X and Y

.. block-code:: bash

  < data/Mus_musculus.NCBIM37.67.gtf grep '^[XY]'
  < data/Mus_musculus.NCBIM37.67.gtf grep '^X\|^Y' | wc -l
  < data/Mus_musculus.NCBIM37.67.gtf grep -E '^X' -E '^Y' | wc -l

3. Count the number of 'CDS' on the chromosome X

.. block-code:: bash

  < data/Mus_musculus.NCBIM37.67.gtf grep 'CDS' | grep '^X' | wc -l

*Use nightingale variant call file (VCF)*

.. block-code:: bash

  cd ~
  sudo cp /data/mus_mda/05-fst2genes/luscinia_vars_flags.vcf.gz ~/data/.
  gunzip data/luscinia_vars_flags.vcf.gz
  less -S data/luscinia_vars_flags.vcf

1. Count the number variants in the file

.. block-code:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | wc -l

2. Count the number of variants passing/failing the quality threshold

.. block-code:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'PASS' | wc -l
  < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'FAIL' | wc -l

3. Count the number of variants on the chromosome Z passing the quality threshold

.. block-code:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'PASS' | grep '^chrZ\s' | wc -l

4. Count the number of records on large autosomes which passed quality threshold

.. block-code:: bash

 < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'PASS' | grep '^chr[1-9]\{1,2\}\s' | wc -l


Cutting out, sorting and replacing text
---------------------------------------

We are going to use these commands: ``cut``, ``sort``, ``uniq``, ``tr``, ``sed``.

*Use nightingale variant call file (VCF)*

1. Which chromosome has the highest and the least number of variants?

.. block-code:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | cut -f 1 | sort | uniq -c | sed 's/^ \{1,\}//' | tr " " "\t" | sort -k1,1nr

2. What is the number of samples in the VCF file?

.. block-code:: bash

  < data/luscinia_vars_flags.vcf grep -v '^##' | head -n1 | cut --complement -f 1-9 | tr "\t" "\n" | wc -l

Joining multiple file + subshell
--------------------------------

``paste``, ``join``

*Use nightingale FASTQ file*

1. Join all nightingale FASTQ files and create a TAB separated file with one line per read

  < cat *.fastq | paste - - - - | cut -f 1-3 | less

2. Make a TAB-separated file having four columns:
    1. chromosome name
    2. number of variants in total for given chromosome
    3. number of variants which pass
    4. number of variants which fails

.. block-code:: bash

  # Command 1
  < data/luscinia_vars_flags.vcf grep -v '^#' | cut -f 1 | sort | uniq -c | sed 's/^ \{1,\}//' | tr " " "\t" > count_vars_chrom.txt

  # Command 2
  < data/luscinia_vars_flags.vcf grep -v '^#' | cut -f 1,7 | sort -r | \
  uniq -c | sed 's/^ \{1,\}//' | tr " " "\t" | paste - - | cut --complement -f 2,3,6 > count_vars_pass_fail.txt

  # Command 3
  join -1 2 -2 3 count_vars_chrom.txt count_vars_pass_fail.txt | wc -l

  # How many lines did you retrieved?

  # You have to sort the data before sending to ``join`` - subshell
  join -1 2 -2 3 <( sort -k2,2 count_vars_chrom.txt ) <( sort -k3,3 count_vars_pass_fail.txt ) | tr " " "\t" > count_all.txt

All three commands together using subshell:

.. block-code:: bash

  join -1 2 -2 3 <( < lp2-var-filtered-rand2.vcf grep -v '^#' | cut -f 1 | sort | uniq -c | \
  sed 's/^ \{1,\}//' | tr " " "\t" | sort -k2,2 ) \
  <( < lp2-var-filtered-rand2.vcf grep -v '^#' | cut -f 1,7 | sort -r | uniq -c | \
  sed 's/^ \{1,\}//' | tr " " "\t" | paste - - | cut --complement -f 2,3,6 | \
  sort -k3,3  ) | tr " " "\t" > count_all.txt


Exercise
--------

How many bases were sequenced?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``wc`` can count characters (think bases) as well. But to get a reasonable number,
we have to get rid of the other lines that are not bases.

One way to do it is to pick only lines comprising of letters A, C, G, T and N.
There is a ubiquitous mini-language called `regular expressions` that can be used
to define text patterns. `A line comprising only of few possible letters` is
a text pattern. ``grep`` is the basic tool for using regular expressions::

  cat *.fastq | grep '^[ACGTN]*$' | less -S

Check if the output looks as expected. This is a very common way to work - build a part of
the pipeline, check the output with ``less`` or ``head`` and fix it or add more commands.

Now a short explanation of the ``^[ACGTN]*$`` pattern (``grep`` works one line a time):

- ``^`` marks beginning of the line - otherwise ``grep`` would search anywhere in the line
- the square brackets (``[]``) are a `character class`, meaning one character of the list, ``[Gg]rep``
  matches ``Grep`` and ``grep``
- the ``*`` is a count suffix for the square brackets, saying there should be zero or more of such characters
- ``$`` marks end of the line - that means the whole line has to match the pattern

To count the bases read, we extend our pipeline::

  cat *.fastq | grep '^[ACGTN]*$' | wc -c

The thing is that this count is not correct. ``wc -c`` counts every character,
and the end of each line is marked by a special character written as ``\n`` (n
for newline). To get rid of this character, we can use another tool, ``tr``
(transliterate). ``tr`` can substitute one letter with another  (imagine you
need to lowercase all your data, or mask lowercase bases in your Fasta file).
Additionally ``tr -d`` (delete) can remove characters::

  cat *.fastq | grep '^[ACGTN]*$' | tr -d "\n" | wc -c

.. note::  If you like regular expressions, you can hone your skills at https://regex.alf.nu/.
