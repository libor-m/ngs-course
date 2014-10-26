Practical part (Unix introduction)
==================================

Basic orientation
-----------------

1. Use of multiple windows
2. Use some 'move around' commands to explore directory structure
3. Prepare data in your home directory

.. note:: Before we start, please, install this software:
  
  .. code-block:: bash
  
    sudo apt-get install htop


**1. Use multiple windows**

  First, type ``screen`` in your terminal:

  .. code-block:: bash

    screen

  You created first window. To create another one press ``ctrl+a c``. To switch between them press ``ctrl+a space``.
 
  To see how useful it is you can try to run ``htop`` in one of these windows

  .. code-block:: bash

    htop

  Htop displays CPU and memory utilization of the (virtual) computer. Continue your work in another one (``ctrl+a space``).
  You can switch back to htop window to monitor progress of some lengthy operation.

**2. Use some 'move around' commands to explore directory structure**

  Your bash session has a `working directory`. All paths you use are supposed to start in your 
  working directory, unless you start them with ``/``. 
  Try the following command in the order they are provided, and figure out what they do.
  Then use your knowledge to explore the directory structure of the virtual machine.

  .. code-block:: bash

    pwd
    ls
    ls /
    cd /
    pwd
    ls
    cd
    pwd

  In case you don't know ask or go to :ref:`moving_around`.

**3. Prepare data in your home directory**

  All data we are going to use today are located in directory ``/data`` and as that they are potentially available to all users. However, we want to have them in our own data directory. Without a need to copy them all to our directory we can have them there using symbolic link (``ln -s``). This keeps data in their original location but creates reference to them.

  In case you did not know where the data are but you knew their names or suffixes you could use ``locate`` command. We know they have fastq, vcf and gff3 suffices:

  .. code-block:: bash

    locate fastq
    locate vcf
    locate gff3

  Once we know their actual position we can create symbolic links:

  .. code-block:: bash

    mkdir data # create directory data
    cd data # go to your new data directory
    ln -s /data/00-reads 00-reads
    ln -s /data/01-genome 01-genome
    ln -s /data/02-variants 02-variants

  We created symbolic links to three places: ``00-reads``, ``01-genome`` and ``02-variants``. We can check them by typing:

  .. code-block:: bash

    ls -l


Installing software
-------------------
The easiest way to install software is via a package manager as you've seen in the beginning (``apt-get`` for all Debian variants).
When the required software is not in the repositories, or one needs the latest version, it's necessary to take the more diffucult path.
The canonical UNIX way is::

  wget -O - ..url.. | tar xvz   # download the 'tarball' from internet
  cd ..unpacked directory..     # set working directory to the project directory
  ./configure                   # check your system and choose the way to build it
  make && sudo make install     # convert source code to machine code and if successful, copy the results to your system

In our example, some steps are ommited. We'll install ``bedtools`` program from a github repository. 
User installed software can be found in ``~/sw`` directory. To install a new software go to this directory:

  .. code-block:: bash

    cd ~/sw

  When the software source code is in a single file (`tarball`), ``wget`` command is the best option to get the file. The latest versions are usually not packaged, and many of the tools can be found at GitHub. To get stuff from GitHub, ``git clone`` command is usually the easiest option.
 
  .. code-block:: bash

    git clone https://github.com/arq5x/bedtools2

  This creates a `clone` of the online repository in directory ``bedtools2``.

  .. code-block:: bash

    cd bedtools2

  To compile (convert from text source form to machine executable form) software on UNIX use the ``make`` command:

  .. code-block:: bash

    make

  It should take a while, so you can flip to your `htop` window with ``ctrl-a space`` and watch the CPU spin;)

  When ``bedtools`` is compiled you have to copy bedtools binaries to ``/usr/local/bin`` directory for UNIX system to find the program when calling from any place in the system.

  .. warning:: Before you use command below to copy binaries make sure you are really in directory you want to be!
 
  .. code-block:: bash
    
    cd bedtools2/bin
    ls # Check that you are really in directory you want to be!
    sudo cp * /usr/local/bin

  We used two commands: ``sudo`` and ``cp``. The sudo command tells the system that we want to make changes in system directories and as such we are asked for password. This step prevents us from harming system. The ``cp`` command copies all bedtools binaries from local bin directory to the system binary repository.
  
.. note:: We used the ``*`` symbol which tells the system that all files in the current directory should be selected. We explain this later.


FASTQ
-----

Explore lengths of short reads in FASTQ files:

1. Explore FASTQ files
2. Find out how many reads are there in each file and in total
3. Calculate summary statistics of read lengths
4. Find primers in FASTQ files

**1. Explore FASTQ files**

  To view contents of FASTQ file go to ``data/00-reads`` directory see contents of the directory using ``ls`` command and view file using ``less``:

  .. code-block:: bash

    cd data
    ls
    less -SN 00-reads/GS60IET02.RL1.fastq

  .. note:: You don't have to type the whole file name. Try to use TAB completition!
  
  Try use and unuse ``-S`` and ``-N`` options and see what's the difference.

  You can also use ``head`` command to view first lines:

  .. code-block:: bash

    head 00-reads/GS60IET02.RL1.fastq        # the default is to show 10 lines
    head -n 20 00-reads/GS60IET02.RL1.fastq  # to show first 20 lines, use -n 20

  or ``tail`` command to view last 20 lines:

  .. code-block:: bash

    tail -n 20 00-reads/GS60IET02.RL1.fastq

**2. How many reads are there?**
  
  We found out that FASTQ files has certain structure and that each read takes 
  four lines (ID, sequence, quality, +; see NGS formats). To obtain number of 
  sequences in each FASTQ file we build a pipeline by combination of 
  ``grep`` and ``wc`` commands along with UNIX feature called *globbing*.

  First, let's try to see what's globbing:

  .. code-block:: bash

    echo 00-reads/*.fastq

  When bash encounters a special character like ``*`` or ``?``, it tries to match 
  filename patterns in the directory structure, where ``?`` match for any single 
  character and ``*`` for 0 or more any characters, respectively. It can, however, 
  match more complex patterns.

  So let's get back to counts of reads...

  We can first try to get counts of *lines* in each file simply by typing:

  .. code-block:: bash

    wc -l 00-reads/*.fastq

  However, to obtain counts of *reads* in each file we have to select just ID lines using ``grep`` command:

  .. code-block:: bash

    grep "^@[0-9A-Z]*$" 00-reads/*.fastq | wc -l


  Command ``grep`` enables to search file for specific character or string of characters.
  Here, we used so-called regular expressions to specify the pattern ``grep`` is supposed 
  to search for.   Regular expressions is a very concise and 'magical' way to describe text
  ptterns. Let's go through our expression piece by piece.

  - ``^`` marks beginning of line - otherwise grep would search anywhere in the line
  - the square brackets (``[]``) represent a character of given class (0 to 9 or A to Z)
  - the ``*`` is a count suffix for the square brackets, saying there should be zero or more of such characters
  - ``$`` marks end of line - that means the whole line has to match the pattern

  If you like regular expressions, you can hone your skills at https://regex.alf.nu/.

**3. Calculate summary statistics of read lengths**

  In this particular task we will need first two lines (ID, sequence) of FASTQ files for each sequence.
  The simplest way to do that is to use UNIX built-in programmatic interface called ``awk``. This program
  enables to efficiently handle the data of various complexity. We build a bit more complex pipeline which 
  is going to combine awk tool along with other commands (``tr``, ``tail``, ``tabtk``).

  This is how the whole pipeline looks like:

  .. code-block:: bash

    awk '{ if( (NR+3) % 4 == 0 || (NR+2) % 4 == 0 ){print $0} }' 00-reads/*.fastq  | tr '\n@' '\t\n' | tail -n +2 | awk -F $'\t' 'BEGIN{OFS=FS}{ print $1,length($2)}' | tabtk num -c 2

  Now we can go step by step through the proces of building it (this is how we did it, there's no other magic):

  In the first step we are going to send all FASTQ files to command written in ``awk``. This command is supposed to return just ID and sequence for each read (i.e first and second line).

  .. code-block:: bash

    awk '{ if( (NR+3) % 4 == 0 || (NR+2) % 4 == 0 ){print $0} }' 00-reads/*.fastq | less -S

  ``NR`` is an ``awk`` built-in variable set to the number of current line (see reference for others built-in variables).

  Now, we created a file with read IDs and sequences where these two alter by line. To create a file with IDs and sequences on the same line we take advantage of the structure of the file. We use ``tr`` command which replaces and deletes characters in file. 

  .. code-block:: bash

    awk '{ if( (NR+3) % 4 == 0 || (NR+2) % 4 == 0 ){print $0} }' 00-reads/*.fastq | tr '\n@' '\t\n' | tail -n +2 | head

  First we replace symbol for newlines (``\n``) with symbol with TAB (``\t``). This concatenates all lines into one, each one separated by TAB. Second, we want to have record for each read (i.e. ID, sequence) in one line. Thus, we introduce newline symbol (``\n``) instead of @ symbol. Lastly, as we find out that first line is empty, we remove it by invoking tail command. This command with -n +2 option takes all lines throughout the file starting at line two.

  Now, we have TAB delimited file with two columns. The first one is for read ID and second one is the read sequence. However, we are interested in  the length of sequence. So we use awk again to calculate the length for each read:

  .. code-block:: bash

    awk '{ if( (NR+3) % 4 == 0 || (NR+2) % 4 == 0 ){print $0} }' 00-reads/*.fastq | tr '\n@' '\t\n' | tail -n +2 | awk -F $'\t' 'BEGIN{OFS=FS}{ print $1,length($2)}' | head


  TThe syntax of this command is simple. First, we need to set TAB as separator because by default awk considers white space as separator. To set TAB as input and output field separator we use two other built-in variables (``FS``, ``OFS``). The input field separator (``FS``) is set by ``-F`` option. The output field  separator is set in the ``BEGIN{}`` part by passing value of ``FS`` to ``OFS``. Next, in the middle section we print for each line (i.e. each read) the first column (read ID) and length of sequence. The length of sequence is obtained using awk built-in function ``length()``. The ``$'\t'`` is a way how to pass TAB character - because if you just press it on the keyboard, it invokes bash autocompletition and does not type the character.

  Lastly, we calculate read length summary statistics using program ``tabtk`` we installed at the beginning:

  .. code-block:: bash

    awk '{ if( (NR+3) % 4 == 0 || (NR+2) % 4 == 0 ){print $0} }' 00-reads/*.fastq | tr '\n@' '\t\n' | tail -n +2 | awk -F $'\t' 'BEGIN{OFS=FS}{ print $1,length($2)}' | tabtk num -c 2

**4. Find primers in FASTQ files**

  Reads in FASTQ files contain adaptors that were used for reverse transcription of the mRNA. 
  Try to identify them and visualize them using basic UNIX commands.

  First, we store the primer sequences into shell variables which we use later. This steps help us to get families to what it is and how to work with shell variables in UNIX environment.

  To set primer sequences into PRIMER# variable type:

  .. code-block:: bash

    PRIMER1="AAGCAGTGGTATCAACGCAGAGTACGCGGG"
    PRIMER2="AAGCAGTGGTATCAACGCAGAGT"

  To interpret a string as shell variable name, prefix it with ``$``:

  .. code-block:: bash

    echo $PRIMER1
    # AAGCAGTGGTATCAACGCAGAGT

  The ``echo`` command printed contents of the variable. However, the variable can be used in any other command in UNIX. We use them in searching for primers in FASTQ files:

  .. code-block:: bash

    grep --color=always $PRIMER1 00-reads/*.fastq | less -RS

  Here, the ``grep``'s coloured output was sent to ``less`` which kept the colors of the matched primers. To colour matches add ``--color=always`` in ``grep`` command and ``-R`` option in ``less``.

GFF, VCF, BED
-------------

Find SNPs and INDELs identified using reads which overlap with 5' untranslated regions.

1. Explore GFF file
2. Create BED file for 5' untranslated regions
3. Explore VCF files
4. Create BED file for SNPs and INDELs
5. Join the two BED files using BEDTools

**1. Explore GFF file (less)**

**2. Create BED file for 5' untranslated regions**

  The whole command looks like this:

  .. code-block:: bash

    grep 5utr 01-genome/luscinia_small.gff3 | tr '; ' '\t' | sed 's/Name=//' | awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$4-1,$5,$10}' > 01-genome/utrs.bed

  Let's go step by step:

  First, we need to filter out records corresponding to 5' UTRs in the GFF file. For this task we can use ``grep`` function and ``less`` to see the results:

  .. code-block:: bash

    grep 5utr 01-genome/luscinia_small.gff3 | less -S

  Having just 5' UTR records we need to remove and resort some columns. We use combination of ``sed``, ``tr`` and ``awk`` commands:

  .. code-block:: bash

    grep 5utr 01-genome/luscinia_small.gff3 | tr '; ' '\t' | sed 's/Name=//' | less -S

  First, we use ``tr`` command to extract gene ID. We replace semicolon and white space by TAB separator. These replacements cause the INFO column to split into three. Subsequently we delete ``'Name='`` part in the gene ID column using ``sed`` command.

  .. code-block:: bash

    grep 5utr 01-genome/luscinia_small.gff3 | tr '; ' '\t' | sed 's/Name=//' | awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$4-1,$5,$10}' | less -S

  Further, as BEDTools assume zero based coordinate system, we use ``awk`` to subtract one from all start coordinates.

  We can print the whole output into utrs.bed file:

  .. code-block:: bash

    grep 5utr 01-genome/luscinia_small.gff3 | tr '; ' '\t' | sed 's/Name=//' | awk -F $'\t' 'BEGIN{OFS=FS}{print $1,$4-1,$5,$10}' > 01-genome/utrs.bed

**3. Explore VCF file (less)**

**4. Create BED file out of VCF file for SNPs and INDELs**

  .. code-block:: bash

    grep -hv ^# 02-variants/*.vcf | awk -F $'\t' 'BEGIN{OFS=FS}{ if(length($4)==1){ print $1,($2-1),($2+length($4)-1),"SNP"}else{ print $1,($2-1),($2+length($4)-1),"INDEL"} }' > 02-variants/variants.bed

  First, we use inverted grep command (``-v`` option) to remove INFO lines (beginning with ``#`` symbol). Also, as we grep from multiple files (i.e. ``*`` globbing) we use option ``-h`` to suppress file names in the output. Try run grep with and without ``-h`` option:

  .. code-block:: bash

    grep -hv ^# 02-variants/*.vcf | head

  Second, we want to distinguish between SNPs and INDELs and create BED file. The difference is in length of REF column in VCF files. SNPs have always only single character, whereas INDELs have always at least two. So we can use easy ``if()`` condition in ``awk`` based on length of REF column. Also, as in the VCF file is only first position of the variant, when creating BED file one has to calculate the second coordinate. So the start position of a SNP is one minus the actual position, whereas the end position is the actual position:

  .. code-block:: bash

    grep -hv ^# 02-variants/*.vcf | awk -F $'\t' 'BEGIN{OFS=FS}{ if(length($4)==1){ print $1,($2-1),$2,"SNP"}else{ print $1,($2-1),($2+length($4)-1),"INDEL"} }' | head

  Finally, the output can be printed into variants.bed


**5. Join the two BED files using BEDTools**

  Finally, we are interested in how many of SNPs and INDELs are located in 5' UTRs. For this task we use BEDTools that represent a suite of tools to do easily so-called "genome arithmetic".

  Full pipeline:

  .. code-block:: bash

    bedtools intersect -a 01-genome/utrs.bed -b 02-variants/variants.bed -wa -wb | cut -f 4,8 |  sort -k2,2 | bedtools groupby -g 2 -c 1 -o count


  First, we use BEDTools tool ``intersect`` to find an overlap between SNPs, INDELs and 5' UTRs.

  .. code-block:: bash

    bedtools intersect -a 01-genome/utrs.bed -b 02-variants/variants.bed -wa -wb | head

  Here, the ``-a`` and ``-b`` options state for file a and file b. Also, it is necessary to specify which of the two files (or both of them) to print in the output (``-wa``, ``-wb``).

  As you may notice, the output contains eight columnts (i.e. four for each file). For us, however, what is important is only information on gene ID and type of variant (SNPs or INDELs). So we cut out only these two columns using ``cut`` command:

  .. code-block:: bash

    bedtools intersect -a 01-genome/utrs.bed -b 02-variants/variants.bed -wa -wb | cut -f 4,8 | head

  The ``-f`` option in the ``cut`` command states for specification of columns which are supposed to be cut out.

  Now, we want to obtain counts of SNPs and INDELs overlapping with 5' UTRs. We use another BEDTools tool - ``groupby``. This tool enables to group data based on column of choice and to do some summary statistics on another another one. Before grouping, however, we need to sort the data according to the column which we use as a grouping column:

  .. code-block:: bash

    bedtools intersect -a 01-genome/utrs.bed -b 02-variants/variants.bed -wa -wb | cut -f 4,8 | sort -k2,2 | bedtools groupby -g 2 -c 1 -o count

  To sort based on certain column one has to use ``-k`` option along with specification of range (in columns) of sorting. If we want to sort based on one column - as in the case above - we specify range using column position. Here, we sort based on second column so we specify range as ``-k2,2``. The BEDTools tool groupby has several options. ``-g`` option specifies column based on which we group, ``-c`` option specifies column to which we apply summary statistics and ``-o`` option specifies type of summary statistics (see manual at http://bedtools.readthedocs.org).


  
  