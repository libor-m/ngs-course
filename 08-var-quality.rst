Session 7: Variant quality
==========================

In this part you will be working on your own. You're already familiar with the
VCF format and some reformatting and plotting tools. There is a file with
variants from several nightingale individuals::

  /data-shared/vcf_examples/luscinia_vars.vcf.gz

Your task now is:

- set up a new project directory
- pick only data for chromosomes ``chr1`` and ``chrZ``
- extract the sequencing depth ``DP`` from the ``INFO`` column
- extract variant type by checking if the ``INFO`` column contains ``INDEL`` string
- merge these two columns together with the first six columns of the VCF
- create a new git repository and commit your script
- push the repository to GitHub

And a bit of guidance here:

- create a new project directory in your ``projects``
- get rid of the comments (they start with ``#``, that is ``^#`` regular expression)
- filter lines based on chromosomes (``grep -e 'chr1\s' -e 'chrZ\s'``)
- extact the first 6 columns (``cut -f1-6``)
- extract ``DP`` column (``egrep -o 'DP=[^;]*' | sed 's/DP=//'``)
- check each line for ``INDEL`` (``awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}'``)
- merge the data (columns) before loading to R (``paste``)

.. pull-quote:: Good luck! (We will help you;)

Extra R task
------------

- load the data into R
- explore graphically
  - barchart of variant types
  - boxplot of qualities for INDELs and SNPs (use ``scale_y_log10()`` if you don't like the outliers)
  - histogram of qualities for INDELs and SNPs (use ``scale_x_log10()``, ``facet_wrap()``) - what is the problem?

Hints:
- add column names while loading the data with ``read_tsv(..., col_names=c(...))``
