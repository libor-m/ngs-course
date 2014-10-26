Reference manual to UNIX introduction
=====================================

Basic orientation in UNIX
-------------------------

**Multiple windows (screen)**

  You're all used to work with multiple windows (in MS Windows;). You can have them in (remote) UNIX as well. The main benefit, however, is that you can log off and your programs keep running.

  To go into a screen mode type:

  .. code-block:: bash
  
    screen

  Once in a screen mode you can control it after you press the master key: ``ctrl+a key``. To create a new window within the screen mode, press ``ctrl+a c`` (create). To flip among your windows press ``ctrl+a space`` (you flip windows often, it's the biggest key available). To detach screen (i.e. keep your programs running and go home), press ``ctrl+a d`` (detach).

  Once outside the screen mode, to restore screen type:
  
  .. code-block:: bash
  
    screen -r

  To list running screens, type:
  
  .. code-block:: bash
  
    screen -ls

**Controlling processes (htop/top)**

  ``htop`` or ``top`` serve to see actual memory burden for each running process.


**Getting help (man)**

  Type man and name of the command:

  .. code-block:: bash
  
    man screen

.. _moving_around:

Moving around & manipulation with files and directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  Basic commands to move around and manipulate files/directories.

  .. code-block:: bash
  
    pwd    # prints current directory path
    cd     # changes current directory path
    ls     # lists current directory contents
    ll     # lists detailed contents of current directory
    mkdir  # creates a directory
    rm     # removes a file
    rm -r  # removes a directory
    cp     # copies a file/directory
    mv     # moves a file/directory
    locate # tries to find a file by name

  Usage:

  *cd*

  To change into a specific subdirectory, and make it our present working directory:

  .. code-block:: bash

    cd go/into/specific/subdirectory

  To change to parent directory:
  
  .. code-block:: bash
  
    cd ..
  
  To change to home directory:
  
  .. code-block:: bash
  
    cd
  
  To go up one level to the parent directory then down into the directory2:
  
  .. code-block:: bash
  
    cd ../directory2
  
  To go up two levels:
  
  .. code-block:: bash
  
    cd ../../

  *ls*

  To list all (including hidden) files and directories (``-a``) in current in given folder along with human readable (``-h``) size of files (``-s``), type:
  
  .. code-block:: bash
  
    ls -ash

  *mv*

  To move a file data.fastq from current working directory to directory ``/home/directory/fastq_files``, type:
  
  .. code-block:: bash
  
    mv data.fastq /home/directory/fastq_files/data.fastq

  *cp*

  To copy a file data.fastq from current working directory to directory ``/home/directory/fastq_files``, type:
  
  .. code-block:: bash
  
    cp data.fastq /home/directory/fastq_files/data.fastq

  *locate*

  This command enables to find any string on system. It helps find location of given files.

  So to locate file data.fastq type:
  
  .. code-block:: bash
  
    locate data.fastq

  This commands uses database of files and directories which updates just once a day. When you look for recent files you may not find them. So to rearch for these files one has to update database before:
  
  .. code-block:: bash
  
    sudo updatedb

**Symbolic links**

  Symbolic links refer by their names to some files or directories in different location. It is useful when one wants to work with some general files accessible to more users but at the same time to have them in local directory. Also, it is usefull when one works at multiple projects and uses same files (especially large ones). Instead of copying them into each project directory one can use simply symbolic links.

  Symbolic link can be create by:
  
  .. code-block:: bash
  
    ln -s /data/genomes/luscinia/genome.fa genome/genome.fasta

  This command creates symbolic link on file in general location (``/data/genomes/luscinia/genome.fa``) and the link is created in subdirectory to the current working directory (``genome/genome.fasta``).



Exploring and basic manipulation with data
------------------------------------------

  *less*

  Program to view (but not to change) and navigate throughout the contents of text files. As it reads only part of a file on the screen (i.e. does not have to read entire file before starting), it has fast load times for large files.

  To view text file while disabling line wrap and add line numbers add options ``-S`` and ``-N``, respectively.

  .. code-block:: bash
  
    less -SN data.fasta

  To navigate within the text file while viewing use:
  
  
    +-----------+-----------------+
    |  Key      | Command         |
    +===========+=================+
    | Space bar | Next page       |
    +-----------+-----------------+
    | b         | Previous page   |
    +-----------+-----------------+
    | Enter key | Next line       |
    +-----------+-----------------+
    | /<string> | Look for string |
    +-----------+-----------------+
    | <n>G      | Go to line <n>  |
    +-----------+-----------------+
    | h         | Help            |
    +-----------+-----------------+
    | q         |  Quit           |
    +-----------+-----------------+

  *cat*

  Utility which outputs the contents of a specific file and can be used to concatenate and list files.

  .. code-block:: bash
  
    cat seq1_a.fasta seq1_b.fasta > seq1.fasta

  *head*

  By default, this utility prints first 10 lines. The number of first n lines can be specified by ``-n`` option.

  To print first 50 lines type:
  
  .. code-block:: bash
  
    head -n 50 data.txt

  *tail*

  By default, this utility prints last 10 lines. The number of last n lines can be specified by ``-n`` option as in case of head.

  To print last 20 lines type:
  
  .. code-block:: bash
  
    tail -n 20 data.txt

  To skip first few lines in the file (e.g. to remove header line of the file):
  
  .. code-block:: bash
  
    tail -n +2 data.txt

  *grep*

  This utility enables you to search text file(s) for lines matching text patterns. To match given pattern it uses either specific string or regular expressions. Regular expressions enable for certain amount of character variability in searched strings (similar to globbing).

  To obtain one file with list of sequence IDs in multiple fasta files type:
  
  .. code-block:: bash
  
    grep '>' *.fasta > seq_ids.txt


  To print all but #-starting lines from the vcf file use option ``-v`` (invert-match):
  
  .. code-block:: bash
  
    grep -v ^# snps.vcf > snps.tab

  The ^ mark specifies beginning of line (i.e. it skips all # which are not at the beginning of line).
  
  *wc*

  This utility generates set of statistics on either standard input or list of text files. It provides these statistics:
  * line count (``-l``)
  * word count (``-w``)
  * character count (``-m``)
  * byte count (``-c``)
  * length of the longest line (``-L``)

  If specific word provided it returns count of this word in a given file.

  To obtain number of files in a given directory type:
  
  .. code-block:: bash
  
    ls | wc -l

  The ``|`` symbol is explained in further section.
  
  *cut*

  Cut out specific columns (fields/bytes) out of a file. By default, fields are separated by TAB. Otherwise, change delimiter using ``-d`` option. To select specific fields out of a file use ``-f`` option (position of selected fields/columns separated by commas). If needed to complement selected fields (i.e. keep all but selected fields) use ``--complement`` option.

  Out of large matrix select all but first column and row representing IDs of rows and columns, respectively:
  
  .. code-block:: bash
  
    < matrix1.txt tail -n +2 | cut --complement -f 1 > matrix2.txt

  *sort*

  This utility sorts a file based on whole lines or selected columns. To sort numerically use ``-n`` option. Range of columns used as sorting criterion is specified by ``-k`` option.

  Extract list of SNPs with their IDs and coordinates in genome from vcf file and sort them based on chromosome and physical position:
  
  .. code-block:: bash
  
    < snps.vcf grep ^# | cut -f 1-4 | sort -n -k2,2 -k3,3 > snps.tab

  *uniq*

  This utility takes sorted lists and provides unique records and also counts of non-unique records (``-c``). To have more numerous records on top of output use ``-r`` option for ``sort`` command.

  Find out count of SNPs on each chromosome:
  
  .. code-block:: bash
  
    < snps.vcf grep ^# | cut -f 2 | sort | uniq -c > chromosomes.tab

  *tr*

  Replaces or removes specific sets of characters within files.

  To replace a characters a and b in the entire file for characters c and d  type:
  
  .. code-block:: bash
  
    tr 'ab' 'cd' < file1.txt > file2.txt
    
  Multiple consecutive occurrences of specific character can be replaced by single character using ``-s`` option. To remove empty lines type:
  
  .. code-block:: bash
  
    tr -s '\n' < file1.txt > file2.txt

  To replace lower case to upper case in fasta sequence type:
  
  .. code-block:: bash
  
    tr "[:lower:]" "[:upper:]" < file1.txt > file2.txt


Building commands
-----------------

**Globbing**

  Refers to manipulating (searching/listing/etc.) files based on pattern matching using specific characters.
  
  Example:
  
  .. code-block:: bash
  
    :~$ ls
    a.bed b.bed seq1_a.fasta seq1_b.fasta seq2_a.fasta seq2_b.fasta
    :~$ ls *.fasta
    seq1_a.fasta seq1_b.fasta seq2_a.fasta seq2_b.fasta


  Character ``*`` in previous example replaces any number of any characters and it indicates to ``ls`` command to list any file ending with ".fasta".

  However, if we look for fastq instead, we got no result:
  
  .. code-block:: bash
  
    :~$ ls *.fastq
    :~$


  Character ``?`` in following example replaces just right the one character (a/b) and it indicates to ls functions to list files containing `seq2_` at the beginning, any single character in the middle (a/b) and ending with ".fasta"

  .. code-block:: bash

    :~$ ls
    a.bed b.bed seq1_a.fasta seq1_b.fasta seq2_a.fasta seq2_b.fasta
    :~$ ls seq2_?.fasta
    seq2_a.fasta seq2_b.fasta
    
  .. code-block:: bash
    
    :~$ ls
    a.bed b.bed seq1_a.fasta seq1_b.fasta seq2_a.fasta seq2_b.fasta
    :~$ ls seq2_[ab].fasta
    seq2_a.fasta seq2_b.fasta

  One can specifically list altering characters (a,b) using brackets ``[]``. One may also be more general and list all files having any alphabetical character ``[a-z]`` or any numerical character ``[0-9]``:

  .. code-block:: bash

    :~$ ls
    a.bed b.bed seq1_a.fasta seq1_b.fasta seq2_a.fasta seq2_b.fasta
    :~$ ls seq[0-9]_[a-z].fasta
    seq1_a.fasta seq1_b.fasta seq2_a.fasta seq2_b.fasta
    

**TAB completition**

  Using key TAB one can finish unique file names or paths without having to fully type them. (try and see)

  From this perspective it is important to think about names for directories in advance as it can spare you a lot time in future. For instance, when processing data with multiple steps one can use numbers at beginnings of names:

  * 00-beginning
  * 01-first-processing
  * 02-second-processsing
  * ...

**Variables**

  UNIX environment enables to use shell variables. To set primer sequence ``'GATACGCTACGTGC'`` to variable ``PRIMER1`` in a command line and print it on screen using ``echo``, type:
  
  .. code-block:: bash
  
    :~$ PRIMER1=GATACGCTACGTGC
    :~$ echo $PRIMER1
    GATACGCTACGTGC

.. note:: It is good habit in UNIX to use capitalized names for variables: ``PRIMER1`` not ``primer1``.

**Pipes**

  UNIX environment enables to chain commands using pipe symbol ``|``. Standard output of the first command serves as standard input of the second one, and so on.

  .. code-block:: bash
  
    ls | head -n 20

**Subshell**

  Subshell enables to run two commands in parallel and merge the outputs. It can be helpful in dealing with data files headers. Use of subshell enables to remove header, run the set of operations on the data, and later insert the header back to file. The basic syntax is:

  .. code-block:: bash

    (command1 file1.txt && command2 file1.txt) > file2.txt

  To sort data file based on two columns without including header type:
  
  .. code-block:: bash
  
    (head -n 1 file1.txt && tail -n +2 file1.txt | sort -n -k1,1 -k2,2) > file2.txt

  Another use of subshell to produce input on the fly (saving useless intermediate files):

  .. code-block:: bash

    paste <(< file1.txt tr ' ' '\t') <(<file2.txt tr '' '\t') > file3.txt


Advanced text manipulation (sed)
--------------------------------

``sed`` editor enables to manipulate easily text. It provides the same functionality as ``tr`` command does but it enables even more advanced manipulations.

More complex data manipulation (awk)
------------------------------------

``awk`` enables to manipulate text data in a very complex way. In fact, it is a simple programming language with functionality similar to regular programming languages. As such it enables enormous variability in ways of how to process text data.

It can be used to write a short script and which can be chained along with UNIX commands in one pipeline. awk script goes line by line and runs given operation on each line separately. The basic structure of the script is divided into three parts and any of these three parts may or may not be included in the script (according to the intention of user). The first part ``'BEGIN{}'`` conducts operation before going through the input file, the middle part ``'{}'`` goes throughout the input file and conducts operations on each line separately. The last part ``'END{}'`` conducts operation after going through the input file.

The basic syntax:

  .. code-block:: bash
  
    < data.txt awk 'BEGIN{<before data processing>}{<data processing>}END{<after data processing>}' > output.txt

**Built-in variables**

  awk has several built-in variables which can be used to track and process data without having to program specific feature.

  The basic four built-in variables:
  
  * ``FS`` - input field separator
  * ``OFS`` - output field separator
  * ``NR`` - record (line) number
  * ``NF`` - number of fields in record (in line)

More advanced built-in variables are: ``RS``, ``ORS``, ``FILENAME``, ``FNR``

Use of built-in variables:

awk recognizes individual columns based on white symbol (i.e. space). When used different delimiter (e.g. TAB) it has to be specified at the beginning of a script using ``-F`` option. This option passes given character into the script which has to be further passed into built-in variable within the BEGIN part:

  .. code-block:: bash
  
    < data.txt awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$2}' > output.txt

  This command takes file data.txt, extract first two TAB delimited columns of the input file and print them TAB delimited into the output file output.txt. When we look more closely on the syntax we see that the TAB delimiter was set using ``-F`` option. This option corresponds to the ``FS`` built-in variable. As we want TAB delimited columns in the output file we pass ``FS`` to ``OFS`` (i.e. ouput field separator) in the ``BEGIN`` section. Further, in the middle section we print out first two columns which can be extracted by numbers with ``$`` symbol (``$1``, ``$2``). The numbers correspond to position of the column in the input file. We could, of course, use for this operation the ``tr`` command which is even simpler. However, the awk enables to conduct any other operation on given data.

  .. note:: To work with the whole line one can approach it by ``$0``.


The ``NR`` built-in variable can be used to capture each second line in a file type:

  .. code-block:: bash
  
    < data.txt awk '{ if(NR % 2 == 0){ print $0 }}' > output.txt

  The ``%`` symbol represents modulo operator which returns the remainder of division. The ``if()`` condition is used to decide on whether the modulo is 0 or not.

  Here is a bit more complex example of how to use ``awk``. We write a command which retrieves coordinates of introns from coordinates of exons.

  Example of input file:
  
  .. code-block:: bash
  
    GeneID            Chromosome   Exon_Start   Exon_End
    ENSG00000139618   chr13        32315474     32315667
    ENSG00000139618   chr13        32316422     32316527
    ENSG00000139618   chr13        32319077     32319325
    ENSG00000139618   chr13        32325076     32325184
    ...               ...          ...          ...

  The command is going to be as follows:
  
  When we look at the command step by step we first remove header and sort data based on GeneID and Exon_Start columns:
  
  .. code-block:: bash
  
    < exons.txt tail -n +2 | sort -k1,1 -k3,3n | ...

  Further, we write a short script using awk to obtain coordinates of introns:
  
  .. code-block:: bash
  
    ... | awk -F $'\t' 'BEGIN{OFS=FS}{ 
             if(NR==1){ 
               x=$1; end1=$4+1;
             }else{ 
               if(x==$1) {
                   print $1,$2,end1,$3-1; end1=$4+1; 
               }else{ 
                   x=$1; end1=$4+1;
               }
             }
           }' > introns.txt

  In the ``BEGIN{}`` part we set TAB as output field separator. Further, using ``NR==1`` test we set GeneID for first line into ``x`` variable and intron start into end1 variable. Otherwise we do nothing. For others records ``NR > 1`` condition ``x==$1`` test if we are still within the same gene. If so we print exon end from previous line (``end1``) as intron start and exon start of current line we use as intron end. Next, we set new intron start (i.e. exon end from current line) into end1. If we have already moved into new one ``x<>$1``) we repeat procedure for the first line and print nothing waiting for next line.


  