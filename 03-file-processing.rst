Session 3: Plain text file processing in Unix
=============================================

This session focuses on plain text file data extraction/modification
using built-in Unix tools.

Pattern search & regular expressions
------------------------------------

``grep`` command can be used to efficiently match and also retrieve pattern in text files.

.. code-block:: bash

  grep pattern file.txt # Returns lines matching a pattern

  grep -v pattern file.txt # Returns lines not matching a pattern

  grep -E regex file.txt # Returns lines not matching a regex

  grep -c pattern file.txt # Returns number of lines matching a pattern

  grep -B pattern file.txt # Returns number of lines before a line matching a pattern

  grep -o pattern file.txt # Returns only matching part of lines

  man grep # For other options

But what if we want match a variable pattern, e.g. differing length or content?
**Regular expressions** are exactly the tool that we need.

.. code-block:: bash

  ^A         # match A at the beginning of line
  A$         # match A at the end of line
  [0-9]      # match numerical characters
  [A-Z]      # match alphabetical characters
  [ATGC]     # match A or T or G or C
  .          # match any character
  A*         # match A letter 0 or more times
  A{2}     # match A letter exactly 2 times
  A{1,}    # match A letter 1 or more times
  A+         # match A letter 1 or more times (extended regular expressions)
  A{1,3}   # match A letter at least 1 times but no more than 3 times
  AATT|TTAA # match AATT or TTAA
  \s         # match whitespace (also TAB)


``grep -E`` is a useful tool to search for patterns using a mini-language called
**regular expressions**. You can use the ``egrep`` shorthand, which means the same.

Word and line count
-------------------

``wc`` command states for *word count* provides a quick summary statistics on the plain text file
content. Bytes, words and lines can be counted in defined files.

.. code-block:: bash

  wc -c file.txt # Number of bytes in a file
  wc -w file.txt # Number of words in a file
  wc -l file.txt # Number of lines in a file

To calculate the number of bytes, words and lines in every file listed:

.. code-block:: bash

  wc -c *.txt
  wc -w *.txt
  wc -l *.txt

1. Count the number of variants in the file

.. code-block:: bash

  < /data-shared/vcf_examples/luscinia_vars_flags.vcf.gz zcat |
  grep -v '^#' |
  wc -l

2. Count the number of variants passing/failing the quality threshold

.. code-block:: bash

  < /data-shared/vcf_examples/luscinia_vars_flags.vcf.gz zcat |
  grep -v '^#' |
  grep 'PASS' |
  wc -l

  < /data-shared/vcf_examples/luscinia_vars_flags.vcf.gz zcat |
  grep -v '^#' |
  grep 'FAIL' |
  wc -l

3. Count the number of variants on the chromosome Z passing the quality threshold

.. code-block:: bash

  < /data-shared/vcf_examples/luscinia_vars_flags.vcf.gz zcat |
  grep -v '^#' |
  grep 'PASS' |
  grep '^chrZ\s' |
  wc -l


Retrieve & count unique records
-------------------------------

Often times we face a problem of how many unique records we have
in a file or how many there are instances of every unique item.

Unix provides efficient way to cut (``cut``) desired columns and retrieve unique
records for selected column (``sort -u``). Additionaly, we can count the instances (``uniq -c``).

Typical examle of use in bioinformatics is to count the number of genes
or SNPs per contig or chromosome.

To select specified columns ``cut`` command can be used. By default, ``cut``
use whitespace as separator. When there is need to distinguish between
standard whitespece and ``TAB`` (i.e. ``TAB``-separated files) then ``-d $'\t'``
delimiter has to be set explicitly. When all columns except a specific column(s)
are supposed to be selected ``--complement`` flag can be used.

.. code-block:: bash

  cut -f1-3 file.txt
  cut -d $'\t' -f1-3 file.txt
  cut --complement -f4 file.txt # Select all columns except column 4

Content of the file can be sorted based on the specified column (``-k1,1``)
or range of columns (``-k1,3``). When the data needs to be sorted numerically (``-n``)
or in reverse order (``-r``) appropriate flags need to be added. Similarly to ``cut``
command ``sort`` recognizes as separator any whitespace. When ``TAB`` is used as a separator,
to enforce distinction from the possible whitespaces used in the file,
``-d $'\t'`` flag has to be used explicitly.

.. code-block:: bash

  sort -k1,1 file.txt # Sort based on first column
  sort -k1,1 -k2,2nr file.txt # Sort first based on first column, then on second column numerically in reversed order
  sort -k1,3 file.txt # Sort based on range of columns

``sort`` command can also be used to retrieve the unique records using flag ``-u``.
When counts of instances for every unique items are supposed to be provided ``uniq -c``
command should be used in combination with ``sort`` as records before sent to ``uniq``
have to be sorted.

.. code-block:: bash

  sort -u file.txt # Retrieve unique records in the file
  < file.txt sort | uniq -c # Count number of instances for every unique item


String extraction and replacement
---------------------------------

Another common task in bioinformatics is a string extraction and/or replacement.
Often times we need to extract specific piece of information from a complex data.

Typical example is the extraction of a specific value according to a TAG in gff3
or VCF file. As positions of TAGS can differ from line to line, using ``cut``
command is simply not possible. Matching using ``sed`` based on a TAG is
the only possibility. ``regex`` can be used to match appropriate pattern using ``sed``.

Another typical task is a replacing of delimiters. ``tr`` command is very
well suited for this task. ``-d`` flag can be used to remove specific characters
from the file. The whole classes can be replaced which can be for instance
used to change uppercase to lowercase or vice versa.

For extraction of repeating strings ``grep -o`` is the most efficient tool. In bioinformatics 
it can be easily used to match and retrieve microsatellites from the sequence for instance.

.. note::

  Difference between ``sed`` and ``tr``:

  ``tr`` (from TRansliterate) replaces (or deletes) **individual characters**:
  Ideal for removing line ends (``tr -d "\n"``) or replacing some
  separator to TAB (``tr ";" "\t"``).

  ``sed`` replaces (or deletes) **complex patterns**.

Typical usage of ``tr`` is as follows:

.. code-block:: bash

  tr "\t" ";" file.txt # To replace TAB separator to semicolon
  tr -d "\n" file.txt # Remove new line characters (``\n``)
  tr "[A-Z]" "[a-z]" # Replace uppercase to lowercase

To match a specific string in the file ``sed`` can use ``regex`` similar
to ``grep`` command. However, to use full ``regex`` functionality
and simplify the regex syntax, **extended regular expression** flag
``r`` (``-E`` for Mac OSX) has to be used.

Comparison of use of standard ``sed`` and use of ``sed`` with extended
regular expressions ``sed -r``:

.. code-block:: bash

  sed 's/pattern/replacement/'

  # Replace one or more A or C or G or T by N

  # Standard sed
  sed 's/^[^ACGT]\{1,\}/N/'

  # The same thing using extended regular expressions:
  sed -r 's/^[^ACGT]{1,}/N/'
  sed -r 's/^[^ACGT]+/N/'

``sed`` can be used also for string extraction. Matched string designated
to be extracted has to be marked in rounded brackets ``(string)``
and passed to the output with following notation: ``\#`` where # character
states for the position starting with 1 in the matched string (i.e. there can be
multiple extractions from one matched string).

.. code-block:: bash

  # Returns TTTGGG
  echo 'AAATTTCCCGGG' | sed -r 's/A+(T+)C+(G+)/\1\2/'

.. note::

  ``sed -r`` (text Stream EDitor) can do a lot of things, however,
  pattern replacement and extraction is the best thing to use it for.
  The 'sed language' consists of single character commands, and it is no fun
  to code and even less fun to read (what does ``sed 'h;G;s/\n//'`` do?;).
  Use ``awk`` for more complex processing (*see next session*).

``grep -o`` extracts only matching parts of the string. This command can be used
to exctract repeating patterns (i.e. very usefull for extraction of microsatellite 
sequences).

.. code-block:: bash

  # Match AT di-nucleotide twice or more times
  grep -o -E "(AT){2,}"

  # Match GTC tri-nuleotide twice or more times
  grep -o -E "(GTC){2,}"

  # Match any repeating pattern
  grep -o -E "([ATGC]{1,})\1+"


*Exercies: Use nightingale variant call file (VCF)*

1. Which chromosome has the highest and the least number of variants?

.. code-block:: bash

  < data-shared/luscinia_vars_flags.vcf grep -v '^#' |
  cut -f 1 |
  sort |
  uniq -c |
  sed -r 's/^ +//' |
  tr " " "\t" |
  sort -k1,1nr

2. What is the number of samples in the VCF file?

.. code-block:: bash

  < data-shared/luscinia_vars_flags.vcf grep -v '^##' |
  head -n1 |
  cut --complement -f 1-9 |
  tr "\t" "\n" |
  wc -l

Figure out alternative solution for exercise 2.

Join & paste data
-----------------

The final part of this session is joining and pasting data. Here, we seek to merge multiple
files into one. ``join`` command corresponds to standard ``JOIN`` command known from ``SQL`` language.
It joins two files based on specific key column. It assumes that both files contain a column representing
keys the are in both files. **Both files must be sorted by the key before any joining task**.
``join`` command has the same functionality as a standard ``JOIN`` meaning that supports ``INNER JOIN``,
``LEFT JOIN``, ``RIGHT JOIN`` and ``FULL OUTER JOIN`` (`Join types <http://www.sql-join.com/sql-join-types>`_).

By default the column considered to be **key** is the first column in both input files. As already mentioned
the key column needs to be sorted in a same way in both files.

.. code-block:: bash

  sort -k1,1 file1.txt > file1.tmp
  sort -k1,1 file2.txt > file2.tmp
  join file1.tmp file2.tmp > joined-file.txt

If **key** column is at different position it needs to be specified on the input
using ``-1`` and ``-2`` flags:

.. code-block:: bash

  sort -k2,2 file1.txt > file1.tmp # key column on the 2nd position
  sort -k3,3 file2.txt > file2.tmp # key column on the 3rd position
  join -12 -23 file1.tmp file2.tmp > joined-file.txt

To specify that the ``join`` is supposed to print **unpaired** lines corresponding to **left, right
and full outer join**, specification of the file to print unpaired lines from has to be done using
``-a`` flag. Also, ``-e`` flag sets value for missing values

.. code-block:: bash

  # Left join
  join -a1 -e NA file1.tmp file2.tmp > joined-file.txt

  # Right join
  join -a2 -e NA file1.tmp file2.tmp > joined-file.txt

  # Full outer join
  join -a1 -a2 -e NA file1.tmp file2.tmp > joined-file.txt

Another command that can be used to merge two or more files is ``paste``. ``paste`` as opposed to
``join`` simply align files by column (corresponding to ``cbind`` in ``R``). No **key** column
is needed as it assumes **one to one correspondence** between the two files.

.. code-block:: bash

  paste file1.txt file2.txt > file-merged.txt

``paste`` command can also be used for smart file transpositions. ``paste`` by default
expects input multiple files per one line. However, when only one file provided with other
possible file possitions filed with ``-`` the command ``paste`` takes the further columns
from next lines of the only file provided. This feature enables to transpose multiple lines
into one line.

Example:

.. code-block:: bash

  file.txt

  item-line1
  item-line2
  item-line3
  item-line4

  < filte.txt paste - -

  item-line1  item-line2
  item-line3  item-line4

**This feature can be used in bioinformatics to convert ``.fastq`` files into ``.tab``
type separated files with one read per line.** We will use this functionality in upcoming
session.


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
