Other speakers' materials
=========================

Petri Kemppainen - LDna
-----------------------

To install the LDna package in the Virtual machine, start the machine, log in with PuTTY 
and run following commands (takes ~10 minutes to build):

.. code:: bash

  sudo su   # password is 'user'
  apt-get install libcurl3-dev
  R --no-save -q <<< 'install.packages(c("devtools"), repos="http://cran.rstudio.com/")'
  R --no-save -q <<< $'options(repos=structure(c(CRAN="http://cran.rstudio.com/")))\ndevtools::install_github("petrikemppainen/LDna")'

If internet connection does not work in your virtual machine, ask for a new image (Libor or Petri).

If you're interested in installing `LDna` into your very own R/RStudio, you can check the instructions here: https://github.com/petrikemppainen/LDna.

Jean Francois Martin
--------------------
Nothing here yet.