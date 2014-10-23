Course materials preparation
============================ 
This section contains the steps that we did to produce the materials that course participants
got ready-made. That is the **linux machine image**, **online documentation** and the **slide deck**.

Online documentation
--------------------
Login to https://github.com. Create a new project called `ngs-course-nhrady`, with a default readme file.


Clone the project to local machine and initialize `sphinx` docs. Choose ``SSH`` clone link in GitHub.

.. code-block:: bash

  git clone git@github.com:libor-m/ngs-course-nhrady.git

  cd ngs-course-nhrady
  
  # use default answers to all the questions
  # enter project name and version 1.0
  sphinx-quickstart


Now track all files created by `sphinx-quickstart` in current directory with `git` and publish 
to GitHub.

.. code-block:: bash
  
  git add .
  git commit -m 'empty sphinx project'

  # ignore _build directory in git
  echo _build >> .gitignore
  git add .gitignore
  git commit -m 'ignore _build directory'
  
  # publish the first docs
  # setting up argument less git pull with '-u'
  git push -u origin master
  
To get live view of the documents, login to https://readthedocs.org. Your `GitHub` account can be paired with 
`Read the Docs` account in `Edit Profile/Social Accounts`, then you can simply 'import' new projects 
from your GitHub with one click. Import the new project and wait for it to build. After the build
the docs can be found at http://ngs-course-nhrady.readthedocs.org (or click the ``View`` button).
  
Now write the docs, commit and push. Rinse and repeat. Try to keep the commits small, just one change a time.

.. code-block:: bash
  
  git add _whatever_new_files_
  git commit -m '_your meaningful description of what you did here_'
  git push

References that may come handy:

- `Thomas Cokelaer's cheat sheet <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_

VirtualBox image
----------------
Create new VirtualBox machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Linux/Debian (32 bit)
- 1 GB RAM - this can be changed at the users machine, if enough RAM is available
- 12 GB HDD as system drive (need space for basic system, gcc, rstudio and some data)
- setup port forwarding
  - 2222 to 22 (ssh, avoiding possible collisions on linux machines with sshd running)
  - 8787 to 8787 (rstudio server)

Install Debian
^^^^^^^^^^^^^^
Download Debian net install image - use i386 so there is as few problems with virtualization as possible.
Not all machines can virtualize x64.

https://www.debian.org/CD/netinst/

Connect the iso to IDE in the virtual machine. Start the machine. Choose ``Install``.

Mostly the default settings will do.
- English language (it will cause less problems)
- Pacific time zone (it is connected with language, no easy free choice;)
- hostname ``node``, domain ``vbox``
- users: root:debian, user:user
- simple partitioning (all in one partition, no LVM)
- Czech mirror to get fast installer file downloads
- pick only SSH server and Standard system utilities

Log in as root:

.. code-block:: bash

  apt-get install sudo
  usermod -a -G sudo user

Login as user (can be done by ``su user`` in root shell):

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
  sudo apt-get install git screen pv
  
  # add important stuff to python
  sudo apt-get install python-dev python-pip python-virtualenv

This is what it takes to create a basic usable system in VirtualBox.
We can shut it down now with ``sudo shutdown -h now`` and take a snapshot of the machine.
If any installation goes haywire from now on, it's easy to revert to this basic system.

Install additional software
^^^^^^^^^^^^^^^^^^^^^^^^^^^

R is best used in RStudio - server version can be used in web browser.

.. code-block:: bash

  mkdir sw
  cd sw

  # install latest R
  # http://cran.r-project.org/bin/linux/debian/README.html
  sudo bash -c "echo 'deb http://mirrors.nic.cz/R/bin/linux/debian wheezy-cran3/' >> /etc/apt/sources.list"
  sudo apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480
  sudo apt-get update
  sudo apt-get install r-base
  sudo R
  > update.packages(.libPaths(), checkBuilt=TRUE, ask=F)
  > install.packages(c("ggplot2", "dplyr", "reshape2", "GGally", "stringr", "vegan", "svd", "tsne"))

  # RStudio with prerequisities
  wget http://ftp.us.debian.org/debian/pool/main/o/openssl/libssl0.9.8_0.9.8o-4squeeze14_i386.deb
  sudo dpkg -i libssl0.9.8_0.9.8o-4squeeze14_i386.deb
  sudo apt-get install gdebi-core
  wget http://download2.rstudio.org/rstudio-server-0.98.1081-i386.deb
  sudo gdebi rstudio-server-0.98.1081-i386.deb

There are packages that are not in the standard repos, or the versions in the repos is very obsolete.
It's worth it to install such packages by hand, when there is not much dependencies.

.. code-block:: bash

  # pipe viewer
  wget -O - http://www.ivarch.com/programs/sources/pv-1.5.7.tar.bz2 | tar xvj
  cd pv-1.5.7/
  ./configure
  make
  sudo make install
  cd ..

  # parallel
  wget -O - http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2|tar xvj
  cd parallel-20141022/
  ./configure
  make
  sudo make install

  # tabtk
  git clone https://github.com/lh3/tabtk.git
  cd tabtk/
  # no configure in the directory
  make
  # no installation procedure defined in makefile
  # just copy the executable to a suitable location
  sudo cp tabtk /usr/local/bin


Sample datasets
^^^^^^^^^^^^^^^
Use data from my nightingale project, subset the data for two selected chromosomes.

.. code-block:: bash

  # see read counts for chromosomes
  samtools view 41-map-smalt/alldup.bam | mawk '{cnt[$3]++;} END{for(c in cnt) print c, cnt[c];}' | sort --key=2rn,2
  # extract readnames that mapped to chromosome 1 or chromosome Z
  mkdir -p kurz/00-reads
  samtools view 41-map-smalt/alldup.bam | mawk '($3 == "chr1" || $3 == "chrZ"){print $1;}' | sort > kurz/readnames
  parallel "fgrep -A 3 -f kurz/readnames {} | grep -v '^--$' > kurz/00-reads/{/}" ::: 10-mid-split/*.fastq

  # reduce the genome as well
  # http://edwards.sdsu.edu/labsite/index.php/robert/381-perl-one-liner-to-extract-sequences-by-their-identifer-from-a-fasta-file
  perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(chr1 chrZ)}print if $c' 51-liftover-all/lp2.fasta > kurz/20-genome/luscinia_small.fasta

  # subset the vcf file with grep
  # [the command got lost;]

Prepare the ``/data`` folder.

.. code-block:: bash

  sudo mkdir /data
  sudo chmod user:user /data

Transfer the files to the VirtualBox image, /data directory using WinSCP.

Slide deck
----------
The slide deck was created using Adobe InDesign.
