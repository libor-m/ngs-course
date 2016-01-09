Unix - Advanced I
=================

This session focuses on plain text file data extraction/modification
using built-in Unix tools.

Pattern search & regular expressions
------------------------------------

``grep`` is a useful tool to search for patterns using a mini-language called **regular expressions**.

.. code-block:: bash

  ^A         # match A at the beginning of line
  A$         # match A at the end of line
  [0-9]      # match numerical characters
  [A-Z]      # match alphabetical characters
  [ATGC]     # match A or T or G or C
  .          # match any character
  A*         # match A letter 0 or more times
  A\{2\}     # match A letter exactly 2 times
  A\{1,\}    # match A letter 1 or more times
  A+         # match A letter 1 or more times (extended regular expressions)
  A\{1,3\}   # match A letter at least 1 times but no more than 3 times
  AATT\|TTAA # match AATT or TTAA
  \s         # match whitespace (also TAB)

*Use mouse annotation file (GTF)*

.. code-block:: bash

  cd ~
  cp /data/mus_mda/05-fst2genes/Mus_musculus.NCBIM37.67.gtf.gz ~/data
  gunzip data/Mus_musculus.NCBIM37.67.gtf.gz
  less -S data/Mus_musculus.NCBIM37.67.gtf

.. note::

  You can check file permissions by typing ``ll`` instead of ``ls``.
  ``rwx`` stand for *Read*, *Write*, *eXecute*, and are repeated three times,
  for *User*, *Group*, and *Others*. The two names you see next to the 
  permissions are file's owner user and group. 

  You can change the permissions - if you have the permission to do so - 
  by e.g. ``chmod go+w`` - "add write permission to group and others".

1. Count the number of records on the chromosome X

.. code-block:: bash

  < data/Mus_musculus.NCBIM37.67.gtf grep '^X' | wc -l

2. Count the number of records on chromosome X and Y

.. code-block:: bash

  < data/Mus_musculus.NCBIM37.67.gtf grep '^[XY]' | wc -l
  < data/Mus_musculus.NCBIM37.67.gtf grep '^X\|^Y' | wc -l

3. Count the number of 'CDS' on the chromosome X

.. code-block:: bash

  < data/Mus_musculus.NCBIM37.67.gtf grep 'CDS' | grep '^X' | wc -l

*Use nightingale variant call file (VCF)*

.. code-block:: bash

  cd ~
  cp /data/mus_mda/05-fst2genes/luscinia_vars_flags.vcf.gz ~/data
  gunzip data/luscinia_vars_flags.vcf.gz
  less -S data/luscinia_vars_flags.vcf

1. Count the number variants in the file

.. code-block:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | wc -l

2. Count the number of variants passing/failing the quality threshold

.. code-block:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'PASS' | wc -l
  < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'FAIL' | wc -l

3. Count the number of variants on the chromosome Z passing the quality threshold

.. code-block:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'PASS' | grep '^chrZ\s' | wc -l

4. Count the number of records on large autosomes which passed quality threshold

.. code-block:: bash

 < data/luscinia_vars_flags.vcf grep -v '^#' | grep 'PASS' | grep '^chr[1-9]\{1,2\}\s' | wc -l


Cutting out, sorting and replacing text
---------------------------------------

We are going to use these commands: ``cut``, ``sort``, ``uniq``, ``tr``, ``sed``.

.. note::

  ``sed`` (text Stream EDitor) can do a lot of things, however,
  pattern replacement is the best thing to use it for. The 'sed language'
  consists of single character commands, and is no fun to code and even less 
  fun to read (what does ``sed 'h;G;s/\n//'`` do?;). Use ``awk`` for more 
  complex processing.

  General syntax:

  .. code-block:: bash

    sed 's/pattern/replacement/'

    # Replace one or more A or C or G or T by N
    sed 's/^[AGCT]\{1,\}/N/'

    # The same thing using extended regular expressions:
    sed -r 's/^[AGCT]+/N/'


*Use nightingale variant call file (VCF)*

1. Which chromosome has the highest and the least number of variants?

.. code-block:: bash

  < data/luscinia_vars_flags.vcf grep -v '^#' | cut -f 1 | \
  sort | uniq -c | sed 's/^ \{1,\}//' | tr " " "\t" | sort -k1,1nr

  # sed -r (extended regular expressions)
  < data/luscinia_vars_flags.vcf grep -v '^#' | cut -f 1 | \
  sort | uniq -c | sed -r 's/^ +//' | tr " " "\t" | sort -k1,1nr

2. What is the number of samples in the VCF file?

.. code-block:: bash

  < data/luscinia_vars_flags.vcf grep -v '^##' | head -n1 | \
  cut --complement -f 1-9 | tr "\t" "\n" | wc -l

.. note::

  Difference between ``sed`` and ``tr``:

  ``tr`` (from TRansliterate) replaces (or deletes) individual characters:
  Ideal for removing line ends (``tr -d "\n"``) or replacing some 
  separator to TAB (``tr ";" "\t"``).

  ``sed`` replaces (or deletes) complex patterns.

Joining multiple files + subshell
---------------------------------

Use ``paste``, ``join`` commands.

.. note::

  Shell substitution is a nice way to pass a pipeline in a place where a file
  is expected, be it input or output file (Just use the appropriate sign).
  Multiple pipelines can be used in a single command:

  .. code-block:: bash

    cat <( cut -f 1 file.txt | sort -n ) <( cut -f 1 file2.txt | sort -n ) | less

*Use nightingale FASTQ file*

1. Join all nightingale FASTQ files and create a TAB separated file with one line per read

.. code-block:: bash

  # repeating input in paste causes it to take more lines from the same source
  < cat *.fastq | paste - - - - | cut -f 1-3 | less

2. Make a TAB-separated file having four columns:

    1. chromosome name
    2. number of variants in total for given chromosome
    3. number of variants which pass
    4. number of variants which fails

.. code-block:: bash

  # Command 1
  < data/luscinia_vars_flags.vcf grep -v '^#' | cut -f 1 | \
  sort | uniq -c | sed 's/^ \{1,\}//' | tr " " "\t" > count_vars_chrom.txt

  # Command 2
  < data/luscinia_vars_flags.vcf grep -v '^#' | cut -f 1,7 | sort -r | \
  uniq -c | sed 's/^ \{1,\}//' | tr " " "\t" | paste - - | \
  cut --complement -f 2,3,6 > count_vars_pass_fail.txt

  # Command 3
  join -1 2 -2 3 count_vars_chrom.txt count_vars_pass_fail.txt | wc -l

  # How many lines did you retrieved?

  # You have to sort the data before sending to ``join`` - subshell
  join -1 2 -2 3 <( sort -k2,2 count_vars_chrom.txt ) \
  <( sort -k3,3 count_vars_pass_fail.txt ) | tr " " "\t" > count_all.txt

All three commands together using subshell:

.. code-block:: bash

  join -1 2 -2 3 <( < lp2-var-filtered-rand2.vcf grep -v '^#' | cut -f 1 | \
  sort | uniq -c | sed 's/^ \{1,\}//' | tr " " "\t" | sort -k2,2 ) \
  <( < lp2-var-filtered-rand2.vcf grep -v '^#' | cut -f 1,7 | sort -r | uniq -c | \
  sed 's/^ \{1,\}//' | tr " " "\t" | paste - - | cut --complement -f 2,3,6 | \
  sort -k3,3  ) | tr " " "\t" > count_all.txt

  # and indented a bit more nicely
  IN=lp2-var-filtered-rand2.vcf
  join -1 2 -2 3 \
      <( <$IN  grep -v '^#' |
        cut -f 1 | 
        sort | 
        uniq -c | 
        sed 's/^ \{1,\}//' | 
        tr " " "\t" | 
        sort -k2,2 ) \
      <( <$IN grep -v '^#' | 
        cut -f 1,7 | 
        sort -r | 
        uniq -c | 
        sed 's/^ \{1,\}//' | 
        tr " " "\t" | 
        paste - - | 
        cut --complement -f 2,3,6 | 
        sort -k3,3  ) |
    tr " " "\t" \
  > count_all.txt

Exercise
--------

How many bases were sequenced?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``wc`` can count characters (think bases) as well. But to get a reasonable number,
we have to get rid of the other lines that are not bases.

One way to do it is to pick only lines comprising of letters A, C, G, T and N.
There is a ubiquitous mini-language called `regular expressions` that can be used
to define text patterns. `A line comprising only of few possible letters` is
a text pattern. ``grep`` is the basic tool for using regular expressions:

.. code-block:: bash

  cat *.fastq | grep '^[ACGTN]*$' | less -S

Check if the output looks as expected. This is a very common way to work - build a part of
the pipeline, check the output with ``less`` or ``head`` and fix it or add more commands.

Now a short explanation of the ``^[ACGTN]*$`` pattern (``grep`` works one line a time):

- ``^`` marks beginning of the line - otherwise ``grep`` would search anywhere in the line
- the square brackets (``[]``) are a `character class`, meaning one character of the list, ``[Gg]rep``
  matches ``Grep`` and ``grep``
- the ``*`` is a count suffix for the square brackets, saying there should be zero or more of such characters
- ``$`` marks end of the line - that means the whole line has to match the pattern

To count the bases read, we extend our pipeline:

.. code-block:: bash

  cat *.fastq | grep '^[ACGTN]*$' | wc -c

The thing is that this count is not correct. ``wc -c`` counts every character,
and the end of each line is marked by a special character written as ``\n`` (n
for newline). To get rid of this character, we can use another tool, ``tr``
(transliterate). ``tr`` can substitute one letter with another  (imagine you
need to lowercase all your data, or mask lowercase bases in your Fasta file).
Additionally ``tr -d`` (delete) can remove characters:

.. code-block:: bash

  cat *.fastq | grep '^[ACGTN]*$' | tr -d "\n" | wc -c

.. note::  If you like regular expressions, you can hone your skills at http://regex.alf.nu/.
