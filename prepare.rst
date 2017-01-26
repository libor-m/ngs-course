Course materials preparation
============================
This section contains the steps that we did to produce the materials that course participants
got ready-made. That is the **linux machine image**, **online documentation** and the **slide deck**.

Online documentation
--------------------
Login to https://github.com. Create a new project called `ngs-course`, with a default readme file.


Clone the project to local machine and initialize `sphinx` docs. Choose ``SSH`` clone link in GitHub.

.. code-block:: bash

  git clone git@github.com:libor-m/ngs-course.git

  cd ngs-course

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
the docs can be found at http://ngs-course.readthedocs.org (or click the ``View`` button).

Now write the docs, commit and push. Rinse and repeat. Try to keep the commits small, just one change a time.

.. code-block:: bash

  git add _whatever_new_files_
  git commit -m '_your meaningful description of what you did here_'
  git push

References that may come handy:

- `Thomas Cokelaer's cheat sheet <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_

Use http://goo.gl to shorten a link to www.seznam.cz, to get a tracking counter
url for the network connectivity test.

Adding another instance
^^^^^^^^^^^^^^^^^^^^^^^
Check out the version which will serve as starting material, create and publish new branch.

.. code-block:: bash

  git pull
  git checkout praha-january-2016
  git checkout -b praha-january-2017
  git push -u origin praha-january-2017:praha-january-2017

Log in to `Read the Docs`, go to `Admin > Versions
<https://readthedocs.org/dashboard/ngs-course/versions/>`_,
make the new version 'Active', set as the default version.

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
  - 5690 to 5690 (rstudio + shiny)

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
  # # wget impersonating normal browser
  # # good for being tracked with goo.gl for example
  # alias wgets='H="--header"; wget $H="Accept-Language: en-us,en;q=0.5" $H="Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8" $H="Connection: keep-alive" -U "Mozilla/5.0 (Windows NT 5.1; rv:10.0.2) Gecko/20100101 Firefox/10.0.2" --referer=/ '
  nano ~/.bashrc
  . ~/.bashrc

  # set timezone so the time is displayed correctly
  echo "TZ='Europe/Prague'; export TZ" >> ~/.profile

  # some screen settings
  cat > ~/.screenrc << 'EOF'
  hardstatus alwayslastline
  hardstatus string '%{= kG}[%{G}%H%? %1`%?%{g}][%= %{= kw}%-w%{+b yk} %n*%t%?(%u)%? %{-}%+w %=%{g}][%{B}%d.%m. %{W}%c%{g}]'

  defscrollback 20000

  startup_message off
  EOF

  # MOTD
  sudo su
  cat > /etc/motd <<EOF

    _ __   __ _ ___        ___ ___  _   _ _ __ ___  ___
   | '_ \ / _` / __|_____ / __/ _ \| | | | '__/ __|/ _ \
   | | | | (_| \__ \_____| (_| (_) | |_| | |  \__ \  __/
   |_| |_|\__, |___/      \___\___/ \__,_|_|  |___/\___|
          |___/
  EOF
  exit

  # everyone likes git and screen
  sudo apt-get install git screen pv curl wget jq

  # build tools
  sudo apt-get install build-essential pkg-config autoconf

  # add important stuff to python
  sudo apt-get install python-dev python-pip python-virtualenv

  # java because of fastqc
  sudo apt-get install openjdk-7-jre-headless

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
  sudo bash -c "echo 'deb http://mirrors.nic.cz/R/bin/linux/debian jessie-cran3/' >> /etc/apt/sources.list"
  sudo apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480
  sudo apt-get update
  sudo apt-get install r-base

  sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev
  sudo R
  > update.packages(.libPaths(), checkBuilt=TRUE, ask=F)
  > install.packages(c("tidyverse", "shiny", "reshape2", "vegan"))
  > exit()

  # RStudio with prerequisities
  sudo apt-get install gdebi-core
  wget https://download2.rstudio.org/rstudio-server-1.0.136-i386.deb
  sudo gdebi rstudio-server-1.0.136-i386.deb

There are packages that are not in the standard repos, or the versions in the repos is very obsolete.
It's worth it to install such packages by hand, when there is not much dependencies.

.. code-block:: bash

  # pipe viewer
  cd ~/sw
  wget -O - http://www.ivarch.com/programs/sources/pv-1.6.0.tar.bz2 | tar xj
  cd pv-*
  ./configure
  make && sudo make install

  # parallel
  cd ~/sw
  wget -O - http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2|tar xj
  cd parallel-*/
  ./configure
  make && sudo make install

  # tabtk
  cd ~/sw
  git clone https://github.com/lh3/tabtk.git
  cd tabtk/
  # no configure in the directory
  make
  # no installation procedure defined in makefile
  # just copy the executable to a suitable location
  sudo cp tabtk /usr/local/bin

  # fastqc
  cd ~/sw
  wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
  unzip fastqc_*.zip
  rm fastqc_*.zip
  chmod +x FastQC/fastqc

  # vcftools
  cd ~/sw
  wget -O - https://github.com/vcftools/vcftools/tarball/master | tar xz
  cd vcftools*
  ./autogen.sh
  ./configure
  make && sudo make install

  # samtools
  cd ~/sw
  wget -O - https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj
  cd samtools*
  ./configure
  make && sudo make install

  # bcftools
  cd ~/sw
  wget -O - https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 | tar xj
  cd bcftools*
  ./configure
  make && sudo make install

  # htslib (tabix)
  cd ~/sw
  wget -O - https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 | tar xj
  cd htslib*
  ./configure
  make && sudo make install

  # bwa
  cd ~/sw
  wget -O - https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 | tar xj
  cd bwa*
  # add -msse2 to CFLAGS
  nano Makefile
  sudo cp bwa /usr/local/bin
  # copy the man
  sudo bash -c "<bwa.1 gzip > /usr/share/man/man1/bwa.1.gz"

  # velvet
  cd ~/sw
  wget -O - https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz | tar xz
  cd velvet*
  # comment out the -m64 line, we're on x86
  nano Makefile
  make
  sudo cp velveth velvetg /usr/bin/local

  # bedtools
  cd ~/sw
  wget -O - https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz | tar xz
  cd bedtools2/
  make && sudo make install

TODO - future proofing of the installs with getting the latest - but release -
quality code with something like this (does not work with tags yet)::

  gh-get-release() { echo $1 | cut -d/ -f4,5 | xargs -I{} curl -s https://api.github.com/repos/{}/releases/latest | jq -r .tarball_url | xargs -I{} curl -Ls {} | tar xz ;}

Check what are the largest packages::

  dpkg-query -Wf '${Installed-Size}\t${Package}\n' | sort -n

Sample datasets
^^^^^^^^^^^^^^^
Pull some data from SRA. This is prety difficult;)

.. code-block:: bash

  sraget () { wget -O $1.sra http://sra-download.ncbi.nlm.nih.gov/srapub/$1 ;}


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

Transfer the files to the VirtualBox image, /data directory using WinSCP.

Do the quality checks:

.. code-block:: bash

  cd /data/slavici
  ~/sw/FastQC/fastqc -o 04-read-qc --noextract 00-reads/*

  # update the file database
  sudo updatedb

Packing the image
^^^^^^^^^^^^^^^^^

Now shut down the VM and click in VirtualBox main window ``File > Export appliance``. Upload the file to a file sharing
service, and use the `goo.gl` url shortener to track the downloads.

Slide deck
----------
Libor's slide deck was created using Adobe InDesign (you can get the CS2 version almost legally for free).
Vasek's slide deck was created with Microsoft Powerpoint. Images are shamelessly taken from the internet,
with the 'fair use for teaching' policy ;)
