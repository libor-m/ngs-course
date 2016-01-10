Genome assembly
===============

We'll be assembling a genome of *E. coli* from paired Illumina reads. We're
using ``Velvet``, which is a bit obsolete, but nice for teaching purposes.
Runs quite fast and does not consume that much memory. In reality we'd use
something like ``Megahit`` for prokaryotic organisms or ``Spades`` for the
rest today.

You've already got a project directory for the assembly: ``~/projects/assembly``.
You can find linked ``00-reads`` there with the fastq files, and two assembly 
results should the process would fail for anyone.

.. code-block:: bash

  cd ~/projects/assembly
  
  # look what's around
  ll

Velvet is used in two phases, the first phase prepares the reads, the second
phase  does the assembly itself. Open the `Velvet manual
<https://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf>`_. When using anything more
complex than a simple notepad you will actually save your time by reading the
manual. Surprised?;)

Also - run ``screen`` at this point, because you want to continue working, 
while Velvet will be blocking one of the consoles.

.. code-block:: bash

  # load and prepare the reads
  velveth 03-velvet-k21 21 -fastq -short -separate 00-reads/MiSeq_Ecoli_MG1655_50x_R1.fastq 00-reads/MiSeq_Ecoli_MG1655_50x_R2.fastq

  # do the assembly - you can flip through the manual in the meantime..
  velvetg 03-velvet-k21

Running the whole assembly process is actually that simple. What is not simple
is deciding, whether your assembly is correct and whether it is the best one
you can get with your data. There is actually a lot of trial and error involved
if you're decided to get the best one. People usually combine several assemblers,
test several settings for each assembler and then combine the best runs from each 
of the assemblers with another assembler..;) 

The best criteria for evaluating your assembly are usually external - N50 is
nice, but does not tell you much about chimeric contigs for example. So
overlaps with another  assembly of some close species, the number of genes
that can be found using protein sequences from a close species are good metrics.

.. code-block:: bash

  # check the expected (assembly) coverage
  ~/sw/velvet_1.2.10/contrib/estimate-exp_cov/velvet-estimate-exp_cov.pl 03-velvet-k21/stats.txt | less

On the other hand when it's bad, any metrics will do - the reported N50 of 94
basepairs means there is something wrong. Let's try to use the paired reads.

.. code-block:: bash

  velveth 04-velvet-k31-paired 31 -fastq -shortPaired -separate 00-reads/MiSeq_Ecoli_MG1655_50x_R1.fastq 00-reads/MiSeq_Ecoli_MG1655_50x_R2.fastq
  
  velvetg 04-velvet-k31-paired -exp_cov auto -ins_length 150

  # check observerd insert length
  ~/sw/velvet_1.2.10/contrib/observed-insert-length.pl/observed-insert-length.pl 03-velvet-k21
