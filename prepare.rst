
Course materials preparation
============================

VirtualBox image
----------------

Download Debian net install image - use i386 so there is as few problems with virtualization as possible.
Not all machines can virtualize x64.

https://www.debian.org/CD/netinst/

Create new VirtualBox machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Linux/Debian (32 bit)
- 1 GB RAM - this can be changed at the users machine, if enough RAM is available
- 12 GB HDD as system drive (need space for basic system, gcc, rstudio and some data)
- name the machine 'node'
- users: root:debian, user:user
- setup port forwarding
  - 22 to 22 (ssh)
  - 8787 to 8787 (rstudio server)

Log in as root:

.. code-block:: bash

  apt-get install sudo
  usermod -a -G sudo user

Login as user:

.. code-block:: bash

  # colrize prompt - uncomment force_color_prompt=yes
  # add ll alias - uncomment alias ll='ls -l'
  # fast sort and uniq
  # export LC_ALL=C 
  # maximal width of man
  # export MANWIDTH=120
  nano ~/.bashrc
  . ~/.bashrc

  # everyone likes git and screen
  sudo apt-get install git screen
  
  # add important stuff to python
  sudo apt-get install python-dev python-pip python-virtualenv

  # make a vbox snapshot here 'usable system'

  # install CloudBioLinux into virtual environment (not to pollute whole system)
  mkdir sw
  virtualenv py-cbl
  . py-cbl/bin/activate

  # cloudbiolinux installation (https://github.com/chapmanb/cloudbiolinux)
  git clone git://github.com/chapmanb/cloudbiolinux.git
  cd cloudbiolinux
  python setup.py build
  python setup.py install

  # fix problem with distribution (new debian wheezy not yet supported?)
  sudo bash -c "echo DISTRIB_CODENAME=wheezy >> /etc/os-release"

  # choose a minimal flavor for installing
  fab -f fabfile.py -H localhost -p user install_biolinux:flavor=ngs_pipeline_minimal 

  # to ease problem debugging
  sudo apt-get install strace
  sudo updatedb

  # fix /usr/local/lib missing in search path
  sudo bash -c "echo /usr/local/lib >> /etc/ld.so.conf.d/local.conf"
  # rebuild the cache
  sudo ldconfig

R Studio server:

.. code-block:: bash

  # install rstudio server
  # https://www.rstudio.com/ide/download/server.html
  sudo apt-get install gdebi-core
  wget http://download2.rstudio.org/rstudio-server-0.98.507-amd64.deb
  # get old openssl
  wget http://ftp.de.debian.org/debian/pool/main/o/openssl/libssl0.9.8_0.9.8o-4squeeze14_amd64.deb
  dpkg -i libssl0.9.8_0.9.8o-4squeeze14_amd64.deb

  # update R to some decent version
  # http://cran.r-project.org/bin/linux/debian/README.html
  sudo bash -c "echo 'deb http://mirrors.nic.cz/R/bin/linux/debian wheezy-cran3/' >> /etc/apt/sources.list"
  sudo apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480
  sudo apt-get update
  sudo apt-get install r-base
  sudo R
  >  update.packages(.libPaths(), checkBuilt=TRUE, ask=F)

  # add some packages by hand
  curl http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2|tar xvj
  cd parallel-20140422/
  ./configure
  make && sudo make install


Prepare data
------------

Create a subset of nightingale data on other machine:


.. code-block::bash

  # see read counts for chromosomes
  samtools view 41-map-smalt/alldup.bam | mawk '{cnt[$3]++;} END{for(c in cnt) print c, cnt[c];}' | sort --key=2rn,2
  # extract readnames that mapped to chromosome 1 or chromosome Z
  mkdir -p kurz/00-reads
  samtools view 41-map-smalt/alldup.bam | mawk '($3 == "chr1" || $3 == "chrZ"){print $1;}' | sort > kurz/readnames
  parallel "fgrep -A 3 -f kurz/readnames {} | grep -v '^--$' > kurz/00-reads/{/}" ::: 10-mid-split/*.fastq

  # reduce the genome as well
  # http://edwards.sdsu.edu/labsite/index.php/robert/381-perl-one-liner-to-extract-sequences-by-their-identifer-from-a-fasta-file
  perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(chr1 chrZ)}print if $c' 51-liftover-all/lp2.fasta > kurz/20-genome/luscinia_small.fasta

Transfer them to VirtualBox:

.. code-block:: bash

  sudo mkdir /data
  sudo chown user:user /data

Create documentation
--------------------
This was not done in the virtual machine, but belongs to the course preparation...

.. code-block:: bash

  mkdir ngs-course-2014
  cd ngs-course-2014
  
  # use default answers to all the questions
  sphinx-quickstart

  # track the progress with git
  git init
  git commit -a -m "empty docs and slide"

Spare parts
-----------
If instaling from remote machine:
Use fabricant to install cloudbiolinux:
need the 127.0.0.1 otherwise it does not use ssh

.. code-block:: bash

  fab -f fabfile.py -H 127.0.0.1 --port=2222 -u user -p user install_biolinux:flavor=ngs_pipeline_minimal

  # full install - does not work
  fab -f fabfile.py -H localhost install_biolinux
