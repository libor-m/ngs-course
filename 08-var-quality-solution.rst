.. _varq_solution:

Solution to the variant quality task
====================================
This is Libor's solution, YMMV, as you can sometimes see in some 
tech documents (Your mileage may vary).

.. code-block:: bash
    
    mkdir -p ~/projects/qual-exercise/data
    cd ~/projects/qual-exercise
    </data-shared/vcf_examples/luscinia_vars.vcf.gz zcat | grep -v '^#' > data/no-headers.vcf
    
    IN=data/no-headers.vcf
    <$IN cut -f 1-6 > data/cols1-6.tsv
    <$IN egrep -o 'DP=[^;]*' | sed 's/DP=//' > data/col-dp.tsv
    <$IN grep -v '^#' | awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' > data/col-type.tsv

    # check if all the files are of the same length
    wc -l data/*.tsv
    paste data/cols1-6.tsv data/col-dp.tsv data/col-type.tsv > data/cols-all.tsv

The data is ready, switch to R to visualize.

.. code-block:: r

    library(tidyverse)

    setwd("~/projects/qual-exercise")
    read_tsv("cols-all.tsv", 
             col.names=c("chrom", "pos", "dot", "ref", "alt", "qual", "DP", "TYPE")) ->
             d
    
    # few plots to try
    d %>%
      ggplot(aes(TYPE)) + 
      geom_bar()

    d %>%
      ggplot(aes(TYPE, DP)) + 
      geom_boxplot() +
      scale_y_log10()

    d %>%
      ggplot(aes(DP)) + 
      geom_histogram() +
      scale_x_log10()

    d %>%
      filter(qual < 900) %>%
      ggplot(aes(qual)) + 
      geom_histogram() +
      scale_x_log10() + 
      facet_wrap(~TYPE, scales="free_y")
