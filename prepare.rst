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
make the new version 'Active', set as the default version in `Admin > Advanced
<https://readthedocs.org/dashboard/ngs-course/advanced/>`_.

If the version (branch) is not visible yet, do a force build of some previous
version to get a fresh checkout.

Check if webhooks are set up both in `ReadTheDocs > Project > Admin > Integratinos`
and in `GitHub > Settings > Webhooks`.

Slack
^^^^^
http://ngs-course.slack.com/

Add a new channel every year. Add the channel to
`defaults <https://ngs-course.slack.com/admin/settings#default_channels>`_
in slack Admin.

Update invite link in `index.rst` (30 day validity).

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

  apt install sudo
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
  cat > /etc/motd <<"EOF"

    _ __   __ _ ___        ___ ___  _   _ _ __ ___  ___
   | '_ \ / _` / __|_____ / __/ _ \| | | | '__/ __|/ _ \
   | | | | (_| \__ \_____| (_| (_) | |_| | |  \__ \  __/
   |_| |_|\__, |___/      \___\___/ \__,_|_|  |___/\___|
          |___/

  EOF
  exit

  # everyone likes git and screen
  sudo apt install git screen pv curl wget jq locate

  # build tools
  sudo apt install build-essential pkg-config autoconf

  # add important stuff to python
  sudo apt install python-dev python-pip python-virtualenv

  # java because of fastqc
  sudo apt install openjdk-8-jre-headless

This is what it takes to create a basic usable system in VirtualBox. We can shut
it down now with ``sudo shutdown -h now`` and take a snapshot of the machine. If
any installation goes haywire from now on, it's easy to revert to this basic
system.

Install additional software
^^^^^^^^^^^^^^^^^^^^^^^^^^^

R is best used in RStudio - server version can be used in web browser.

.. code-block:: bash

  mkdir sw
  cd sw

  # install latest R
  # https://cran.r-project.org/bin/linux/debian/
  sudo bash -c "echo 'deb http://mirrors.nic.cz/R/bin/linux/debian buster-cran35/' >> /etc/apt/sources.list"
  sudo apt install dirmngr
  sudo apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
  sudo apt update
  sudo apt install r-base

  sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev
  sudo R
  > update.packages(.libPaths(), checkBuilt=TRUE, ask=F)
  > install.packages(c("tidyverse", "shiny", "reshape2", "vegan"))
  > quit(save="no")

  # RStudio with prerequisities
  sudo apt install gdebi-core

  # 1.1.463 is the latest 32 bit version, no more updates...
  # https://support.rstudio.com/hc/en-us/articles/206569407-Older-Versions-of-RStudio
  wget https://download2.rstudio.org/rstudio-server-1.1.463-i386.deb

  # https://rstudio.com/products/rstudio/download-server/debian-ubuntu/
  # 64 bit
  wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.2.5019-amd64.deb

  # occasionally it's necessary to install older libssl
  # see https://unix.stackexchange.com/a/394462
  # go to https://packages.debian.org/jessie/i386/libssl1.0.0/download
  # copy .deb the link there, do gdebi .deb
  sudo gdebi rstudio-server-*.deb
  # and fix upstart config
  # https://support.rstudio.com/hc/en-us/community/posts/200780986-Errors-during-startup-asio-netdb-error-1-Host-not-found-authoritative-
  # remove 2 from [2345]
  sudo nano /usr/lib/rstudio-server/extras/upstart/rstudio-server.conf
  rm rstudio-server-*.deb

Open http://localhost:8787 and reconfigure layout and colors.

There are packages that are not in the standard repos, or the versions in the
repos is very obsolete. It's worth it to install such packages by hand, when
there is not much dependencies.

.. code-block:: bash

  mkdir ~/sw

  # install a tar with the most common method
  inst-tar() {
    cd ~/sw
    wget -O - "$1" | tar xj
    # extract possible dir name from the tar path
    cd $( echo "$1" | egrep -o '/[^-/]+-' |  sed 's/^.//;s/$/*/' )
    ./configure
    make && sudo make install
  }

  # pipe viewer
  inst-tar http://www.ivarch.com/programs/sources/pv-1.6.6.tar.bz2

  # parallel
  inst-tar http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2

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
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
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
  inst-tar https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2

  # bcftools
  inst-tar https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2

  # htslib (tabix)
  inst-tar https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2

  # bwa
  cd ~/sw
  wget -O - https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar xj
  cd bwa*
  # 32 bit: add -msse2 to CFLAGS
  # nano Makefile
  make
  sudo cp bwa /usr/local/bin
  # copy the man
  sudo bash -c "<bwa.1 gzip > /usr/share/man/man1/bwa.1.gz"

  # velvet
  cd ~/sw
  wget -O - https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz | tar xz
  cd velvet*
  # 32 bit: comment out the -m64 line, we're on x86
  # nano Makefile
  make
  sudo cp velveth velvetg /usr/local/bin

  # bedtools
  cd ~/sw
  wget -O - https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz | tar xz
  cd bedtools2/
  make && sudo make install

  # htop if network fails
  wget http://ftp.cz.debian.org/debian/pool/main/h/htop/htop_2.2.0-2_i386.deb
  wget http://ftp.cz.debian.org/debian/pool/main/h/htop/htop_2.2.0-2_amd64.deb
  # then gdebi htop* at the lesson

  # clean up
  rm -rf bcftools-*/ bedtools2/ bwa-*/ htslib-*/ parallel-*/ pv-*/ samtools-*/ tabtk/ vcftools-vcftools-*/

TODO - future proofing of the installs with getting the latest - but release -
quality code with something like this (does not work with tags yet)::

  gh-get-release() { echo $1 | cut -d/ -f4,5 | xargs -I{} curl -s https://api.github.com/repos/{}/releases/latest | jq -r .tarball_url | xargs -I{} curl -Ls {} | tar xz ;}

Check what are the largest packages::

  dpkg-query -Wf '${Installed-Size}\t${Package}\n' | sort -n

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

Transfer the data to `user` directory (`root` cannot log in remotely):

.. code-block:: bash

  # on host machine
  cd somewhere.../data-pack
  scp -P 2222 -r data-shared user@localhost:~
  scp -P 2222 -r home/user/projects user@localhost:~

  # hyperv non-localhost
  VM=192.168.62.71
  scp -r data-shared "user@$VM:~"
  scp -r home/user/projects "user@$VM:~"

Back on the guest machine.

.. code-block:: bash

  # make the shared data 'shared'
  sudo mv ~/data-shared /

  # change permissons back to 'read only' for user
  sudo chown -R root:root /data-shared

Cleanup
^^^^^^^

.. code-block:: bash

  # update the file database
  sudo updatedb

  # remove history not to confuse users
  sudo su
  history -cw

  # ctrl-d
  history -cw

Packing the image
^^^^^^^^^^^^^^^^^
Now shut down the VM and click in VirtualBox main window ``File > Export
appliance``. Upload the file to a file sharing service, and use the `goo.gl` url
shortener to track the downloads.

Slide deck
----------
Libor's slide deck was created using Adobe InDesign (you can get the CS2 version
almost legally for free). Vasek's slide deck was created with Microsoft
Powerpoint. Images are shamelessly taken from the internet, with the 'fair use
for teaching' policy ;)
