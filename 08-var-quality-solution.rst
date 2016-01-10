.. _varq_solution:

Solution to the variant quality task
====================================
This is Libor's solution, YMMV, as you can sometimes see in some 
tech documents (Your mileage may vary).

.. code-block:: bash
    
    cd ~/data
    mkdir varq
    cd varq
    cat /data/slavici/02-variants/*.vcf | grep -v '^#' > pure-data.vcf
    
    IN=/data/vcf_examples/luscinia_vars.vcf.gz
    <$IN zcat | cut -f 1-6 > cols1-6.tsv
    <$IN zcat | egrep -o 'DP=[^;]*' | sed 's/DP=//' > col-dp.tsv
    # <$IN zcat | egrep -o 'TYPE=[^;]*' | sed 's/TYPE=//' > col-type.tsv
    <$IN zcat | grep -v '^#' | mawk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' > col-type.tsv

    # check if all the files are of the same length
    wc -l *.tsv
    paste cols1-6.tsv col-dp.tsv col-type.tsv > cols-all.tsv
    

The data is ready, switch to R to visualize.

.. code-block:: r

    library(ggplot2)
    library(dplyr)

    setwd("~/data/varq/")
    d <- read.delim("cols-all.tsv", 
                    col.names=c("chrom", "pos", "dot", "ref", "alt", "qual", "DP", "TYPE"))
    
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
