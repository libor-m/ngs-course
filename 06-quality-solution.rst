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
    
    <pure-data.vcf cut -f 1-6 > cols1-6.tsv
    <pure-data.vcf egrep -o 'DP=[^;]*' | sed 's/DP=//' > col-dp.tsv
    <pure-data.vcf egrep -o 'TYPE=[^;]*' | sed 's/TYPE=//' > col-type.tsv
    
    # check if all the files are of the same length
    wc -l *.tsv
    #  62313 col-dp.tsv
    #  62313 col-type.tsv
    #  62313 cols1-6.tsv
    paste cols1-6.tsv col-dp.tsv col-type.tsv > cols-all.tsv
    

The data is ready, switch to R to visualize.

.. code-block:: r

    library(ggplot2)
    library(dplyr)

    setwd("~/data/varq/")
    d <- read.delim("cols-all.tsv", 
                    col.names=c("chrom", "pos", "dot", "ref", "alt", "qual", "DP", "TYPE"))
    
    # alpha=0.1 makes the points transparent
    # and thus helps with overplotting (too much points in the same place)
    d %>%
      filter(!grepl(",", TYPE)) %>%
      ggplot(aes(DP, qual)) + 
      geom_point(alpha=0.1) + 
      facet_wrap(~TYPE) +
      scale_x_log10() +
      scale_y_log10()
      
