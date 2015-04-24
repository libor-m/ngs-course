Graphics session
================
.. pull-quote:: A picture is worth a thousand words. 

Especially when your data is big. We'll try to show you one of the easiest
ways to get nice pictures from  your UNIX. We'll be using R, but we're not
trying to teach you R. R Project is huge, and mostly a huge mess. We're cherry
picking just the best bits;)

Mouse variants
^^^^^^^^^^^^^^
R is best for working with 'tables'. That means data, where each line 
contains the same amount of 'fields', delimited by some special character
like ``;`` or ``<tab>``. The first row can contain column names. VCF is 
almost a nice tabular file. The delimiter is ``<tab>``, but there is some mess
in the beginning of the file::

  </data/mus_mda/00-popdata/popdata_mda_euro.vcf less -S


Prepare the input file
----------------------
We want to get rid of the comment lines starting with ``##``, and keep the 
line starting with ``#`` as column names (getting rid of the ``#`` itself):

.. code-block:: bash

   # create a new 'project' directory in data
   mkdir ~/data/plotvcf

   # we'll be reusing the same long file name, store it into a variable
   IN=/data/mus_mda/00-popdata/popdata_mda_euro.vcf

   # get rid of the '##' lines (quotes have to be there, otherwise
   # '#' means a comment in bash)
   <$IN grep -v '##' | less -S

   # good, now trim the first #
   <$IN grep -v '##' | tail -c +2 | less -S

   # all looks ok, store it (tsv for tab separated values)
   <$IN grep -v '##' | tail -c +2 > ~/data/plotvcf/popdata_mda_euro.tsv

Now we will switch to R Studio. You can just click here: `Open RStudio <http://localhost:8787>`_.

In R Studio choose ``File > New file > R Script``. R has a working directory as well.
You can change it with ``setwd``. Type this code into the newly created file::

  setwd('~/data/plotvcf')

With the cursor still in the ``setwd`` line, press ``ctrl+enter``. This copies the command
to the console and executes it. Now press ``ctrl+s``, and save your script as ``plots.R``.
It is a better practice to write all your commands in the script window, and execute with 
``ctrl+enter``. You can comment them easily, you'll find them faster...

Load and check the input
------------------------
Tabular data is loaded by ``read.table`` and it's shorthands. On a new line, type
``read.table`` and press ``F1``. Help should pop up. We'll be using the ``read.delim`` 
shorthand, that is preset for loading ``<tab>`` separated data with US decimal separator::

  d <- read.delim('popdata_mda_euro.tsv')

A new entry should show up in the 'Environment' tab. Click the arrow and explore. Click the 
'd' letter itself.

You can see that ``CHROM`` was encoded as a number only and it was loaded as
``integer``. But in fact it is a factor, not a number (remember e.g.
chromosome X). Fix this in the ``read.delim`` command, loading the data again
and overwriting `d`. The plotting would not work otherwise::

  d <- read.delim('popdata_mda_euro.tsv', colClasses=c("CHROM"="factor"))

First plot
----------
We will use the ``ggplot2`` library. The 'grammatical' structure of the
command says what to plot, and how to represent the values. There are some
sensible defaults - e.g. ``geom_bar`` of a factor sums the observations::

  library(ggplot2)
  ggplot(d, aes(CHROM)) + geom_bar()

This shows the number of variants in each chromosome. You can see here, that
we've included only a subset of the data, comprising chromosomes 2 and 11.

Summarize the data
------------------
We're interested in variant density along the chromosomes. We can simply
break the chromosome into equal sized chunks, and count variants in each of them
as a measure of density.

There is a function ``round_any`` in the package ``plyr``, which given
precision rounds the numbers. We will use it to round the variant position to
1x10^6 (million base pairs), and then use this rounded position as the block
identifier. Because the same positions repeat on each chromosome, we need to
calculate it once per each chromosome. This is guaranteed by ``group_by``.
``mutate`` just adds a column to the data.

You're already used to pipes from the previous exercises. While it's not
common in R, it is possible to build your commands in a similar way thanks to
the ``magrittr`` package. The name of the package is an homage to the Belgian
surrealist René Magritte and his most popular painting.

.. image:: _static/magritte.jpg
   :align: center
   :alt: Ceci n'est pas une pipe. This is not a pipe.

Although the magrittr ``%>%`` operator is not a pipe, it behaves like one. You
can chain your commands like when building a bash pipeline:

.. code-block:: r

   dc <- d %>% group_by(CHROM) %>% mutate(POS_block=round_any(POS, 1e6))

   # the above command is equivalent to 
   dc <- mutate(group_by(d, CHROM), POS_block=round_any(POS, 1e6))


  