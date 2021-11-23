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

Cloud image
-----------
Create new machine
^^^^^^^^^^^^^^^^^^
We expect ~16 participants. To make things simple we'll host them all on a single instance.

Follow the Meta Cloud `quick start <https://cloud.gitlab-pages.ics.muni.cz/documentation/quick-start/>`_.
Briefly:

- add ssh keys
- add SSH and ICMP security rules (more rules later)
- `Compute > Instance > Launch instance`, fill this in the wizard dialog

    - Debian (64 bit)
    - flavor `hpc.16core-32ram`
    - 32 GB RAM - little less than 2 GB per user
    - 16 vCPUs - keep 2 of the allowed 18 for the testing instance
    - 160 GB HDD as system drive (need space for basic system, gcc, rstudio and produced data * N participants)

- more rules in security group

  - HTTP to set up let's encrypt cert
  - 443 for secured RStudio
  - 60k-61k for mosh
  - 5690 rstudio + shiny



Debian conifg
^^^^^^^^^^^^^
SSH to the machine - read the IP in the OpenStack interface and log in with `debian`
user name.

.. code-block:: bash

  ssh debian@${INSTANCE_IP}

  # start as super user
  sudo su

  # Prague time zone
  dpkg-reconfigure tzdata

  # find fastest mirror
  apt install netselect-apt

  # patch it in sources.list
  vi /etc/sources.list

  # upgrade all
  apt update
  apt upgrade

  # keep the sources list over reboot
  # +apt_preserve_sources_list: true
  vi /etc/cloud/cloud.cfg

  # install the basic tools for more configuration work
  apt install vim screen mosh git

  # log in as debian
  su debian

  # create an ssh key
  ssh-keygen -t ed25519

  # checkout dotfiles
  git clone git@github.com:libor-m/dotfiles.git

  # link vim config
  ln -s dotfiles/vim/.vimrc .

  # back to root shell
  exit

  # link vim config for root
  cd
  ln -s ~debian/dotfiles/vim/.vimrc .

Now it should be easy to work as `debian` user, with vim configured even for sudo.

Tiny fixes to make work as `debian` pleasurable.

.. code-block:: bash

  # colrize prompt - uncomment force_color_prompt=yes
  # add ll alias - uncomment alias ll='ls -l'
  # export MANWIDTH=120
  vi ~/.bashrc
  . ~/.bashrc

Set up the user skeleton, so the newly created users will be set up as needed.
Fancy login message will sure help;)

.. code-block:: bash

  sudo su

  # colrize prompt - uncomment force_color_prompt=yes
  # add ll alias - uncomment alias ll='ls -l'
  # fast sort and uniq
  # export LC_ALL=C
  # maximal width of man
  # export MANWIDTH=120
  # # wget impersonating normal browser
  # # good for being tracked with goo.gl for example
  # alias wgets='H="--header"; wget $H="Accept-Language: en-us,en;q=0.5" $H="Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8" $H="Connection: keep-alive" -U "Mozilla/5.0 (Windows NT 5.1; rv:10.0.2) Gecko/20100101 Firefox/10.0.2" --referer=/ '
  vi /etc/skel/.bashrc

  # some screen settings
  cat > /etc/skel/.screenrc << 'EOF'
  hardstatus alwayslastline
  hardstatus string '%{= kG}[%{G}%H%? %1`%?%{g}][%= %{= kw}%-w%{+b yk} %n*%t%?(%u)%? %{-}%+w %=%{g}][%{B}%d.%m. %{W}%c%{g}]'

  defscrollback 20000

  startup_message off
  EOF

  # basic RStudio ide config
  # obtained by configuring one instance for liborm and then copying the
  # resulting file
  mkdir -p /etc/skel/.config/rstudio
  cat > /etc/skel/.config/rstudio/rstudio-prefs.json <<'EOF'
  {
      "save_workspace": "never",
      "font_size_points": 11,
      "editor_theme": "Solarized Dark",
      "panes": {
          "quadrants": [
              "TabSet1",
              "TabSet2",
              "Source",
              "Console"
          ],
          "tabSet1": [
              "Environment",
              "History",
              "Files",
              "Connections",
              "Build",
              "VCS",
              "Tutorial",
              "Presentation"
          ],
          "tabSet2": [
              "Plots",
              "Packages",
              "Help",
              "Viewer"
          ],
          "console_left_on_top": false,
          "console_right_on_top": false
      },
      "posix_terminal_shell": "bash"
  }
  EOF

  # MOTD
  cat > /etc/motd <<"EOF"

    _ __   __ _ ___        ___ ___  _   _ _ __ ___  ___
   | '_ \ / _` / __|_____ / __/ _ \| | | | '__/ __|/ _ \
   | | | | (_| \__ \_____| (_| (_) | |_| | |  \__ \  __/
   |_| |_|\__, |___/      \___\___/ \__,_|_|  |___/\___|
          |___/

  EOF
  exit


Install some basic software

.. code-block:: bash

  sudo apt install pv curl wget jq locate

  # build tools
  sudo apt install build-essential pkg-config autoconf

  # add important stuff to python
  sudo apt install python-dev python-pip python-virtualenv

  # java because of fastqc
  # sudo apt install openjdk-8-jre-headless

  # let's try default jre
  sudo apt install default-jre-headless

Set up a dynamic DNS to get some nice login name.

.. code-block:: bash

  cd
  ln -s dotfiles/duckdns

  cat duckdns/duck.cron
  # add the printed line to crontab
  crontab -e

This is what it takes to create a basic usable system in VirtualBox. We can shut
it down now with ``sudo shutdown -h now`` and take a snapshot of the machine. If
any installation goes haywire from now on, it's easy to revert to this basic
system.

Install R and RStudio
^^^^^^^^^^^^^^^^^^^^^

R is best used in RStudio - server version can be used in web browser.

.. code-block:: bash

  mkdir ~/sw
  cd ~/sw

  # install latest R
  # https://cran.r-project.org/bin/linux/debian/
  sudo bash -c "echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' > /etc/apt/sources.list.d/cran.list"
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
  wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.3.1093-amd64.deb
  sudo gdebi rstudio-server-*.deb

  # and fix upstart config
  # https://support.rstudio.com/hc/en-us/community/posts/200780986-Errors-during-startup-asio-netdb-error-1-Host-not-found-authoritative-
  # remove 2 from [2345]
  sudo nano /usr/lib/rstudio-server/extras/upstart/rstudio-server.conf

  # install nginx as a front end
  # snapd is needed for certbot ;(
  sudo apt install nginx snapd

  # test if http is accessible from local browser

  # simple nginx proxy config for rstudio
  sudo su
  cat > /etc/nginx/sites-enabled/ngs-course.duckdns.org <<'EOF'
    map $http_upgrade $connection_upgrade {
    default upgrade;
    ''      close;
    }

    server {
    location / {
        proxy_pass http://localhost:8787;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection $connection_upgrade;
        proxy_read_timeout 20d;
    }

    server_name ngs-course.duckdns.org;

    listen 80;

    }
  EOF

  # remove the default site
  rm /etc/nginx/sites-enabled/default

  # test and reload
  nginx -t
  nginx -s reload

  # test if RStudio login page is visible at http
  # .. we'll use the non-sudo account to access rstudio later

  # secure with certbot
  # (snap paths are somehow broken..and restarting the whole system is soo windows98)
  /snap/bin/certbot --nginx

Install additional software
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are packages that are not in the standard repos, or the versions in the
repos is very obsolete. It's worth it to install such packages by hand, when
there is not much dependencies.

.. code-block:: bash

  mkdir -p ~/sw

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
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
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
  inst-tar https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2

  # bcftools
  inst-tar https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2

  # htslib (tabix)
  inst-tar https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2

  # bwa
  cd ~/sw
  wget -O - https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar xj
  cd bwa*
  make
  sudo cp bwa /usr/local/bin
  # copy the man
  sudo bash -c "<bwa.1 gzip > /usr/share/man/man1/bwa.1.gz"

  # velvet
  cd ~/sw
  wget -O - https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz | tar xz
  cd velvet*
  make
  sudo cp velveth velvetg /usr/local/bin

  # bedtools
  cd ~/sw
  wget -O - https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz | tar xz
  cd bedtools2/
  make && sudo make install

  # clean up
  rm -rf bcftools-*/ bedtools2/ bwa-*/ htslib-*/ parallel-*/ pv-*/ samtools-*/ tabtk/ vcftools-vcftools-*/

TODO - future proofing of the installs with getting the latest - but release -
quality code with something like this (does not work with tags yet)::

  gh-get-release() { echo $1 | cut -d/ -f4,5 | xargs -I{} curl -s https://api.github.com/repos/{}/releases/latest | jq -r .tarball_url | xargs -I{} curl -Ls {} | tar xz ;}

Check what are the largest packages::

  dpkg-query -Wf '${Installed-Size}\t${Package}\n' | sort -n

Create the user accounts
^^^^^^^^^^^^^^^^^^^^^^^^
For a multi-user machine, we need the low-privileged accounts and at least a quota
to prevent DoS by overfilling the disk.

Name the accounts `user01` to `user22`:

.. code-block:: bash

  sudo su
  cd

  # aptitude search '?provides(wordlist)'
  apt install wamerican

  # generate some funny passwords
  </usr/share/dict/words egrep "^[a-z]{5,8}$" |
    sort -R |
    paste -d' ' - - - |
    head -22 |
    nl -w2 -n'rz' |
    sed 's/^/user/' \
  > users.tsv

  # use `adduser` as debian alternative
  # --gecos '' --disabled-password to get unattended run
  adduser --gecos '' --disabled-password liborm
  adduser --gecos '' --disabled-password janouse1
  usermod -a -G sudo liborm
  usermod -a -G sudo janouse1

  # normal users
  <users.tsv cut -f1 | xargs -n1 adduser --gecos '' --disabled-password

  # use chpasswd to update the passwords
  <users.tsv tr "\t" ":" | chpasswd

  # add quotas
  # https://www.digitalocean.com/community/tutorials/how-to-set-filesystem-quotas-on-debian-10
  apt install quota
  # add ,usrquota to / mount
  vi /etc/fstab
  mount -o remount /
  quotacheck -ugm /
  quotaon -v /
  <users.tsv cut -f1 | xargs -I{} setquota -u {} 8G 10G 0 0 /

  # copy-paste users.tsv to shared google sheet
  # delete on disk
  rm users.tsv

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

  VM=ngs-course.duckdns.org
  scp -r data-shared "debian@${VM}:~"
  scp -r home/user/projects "debian@${VM}:~"

On the remote machine:

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

Update the machine
^^^^^^^^^^^^^^^^^^
When Debian + RStudio are reasonably updatable, we can keep the previous image.
Hostname is derived from instance name via `cloud-init`, so renaming the instance in
OpenStack should do the trick. Still `/etc/hosts` need to be edited to make `sudo` happy.

.. code-block:: bash

  # as root
  sudo su

  # general update
  # (add new CRAN key)
  apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'

  apt update
  apt upgrade

  # update certificates
  snap refresh
  certbot certonly --nginx
  systemctl restart nginx

  # update R packages
  R
  > update.packages(lib.loc=.libPaths()[1], ask=F, checkBuilt=T, Ncpus=15)

  # update rstudio as normal user
  cd ~/sw
  wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2021.09.1-372-amd64.deb
  sudo rstudio-server active-sessions
  sudo rstudio-server offline
  sudo gdebi rstudio-server-2021.09.1-372-amd64.deb
  sudo rstudio-server online


Slide deck
----------
Libor's slide deck was created using Adobe InDesign (you can get the CS2 version
almost legally for free). Vasek's slide deck was created with Microsoft
Powerpoint. Images are shamelessly taken from the internet, with the 'fair use
for teaching' policy ;)
