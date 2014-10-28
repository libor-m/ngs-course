Doing useful stuff
==================

Moving files around
-------------------
In the ``/data`` directory you've got sample data with some precomputed results
for the case some computations fail. You don't want to overwrite those,
so you will create a 'clean' directory with only the input data.

Using links you can access to the same data from different locations:

.. warning:: Linux uses case sensitive filesystems - File is not file.

.. code-block:: bash

  # create a sandbox directory
  cd /data
  mkdir slavici_sandbox
  cd slavici_sandbox

  # link the data from the original directory
  ln -s ../slavici/00-reads
  ln -s ../slavici/01-genome

  # readgroups is small, we can copy it
  cp ../slavici/20-smalt/readgroups.txt

  # check if there is everythig we need
  ll
  ll 00-reads

Installing software
-------------------
There is a canonical software install procedure in UNIX. It can be summarized as

.. code-block:: bash

  wget -O - http://some.site.com/package.tar.gz | tar xvz
  cd package
  less R<tab> # usually README or README.txt, type capital R!
  ./configure
  make
  sudo make install

Looks easy .. ? Some packages do not have ``configure`` file, you just skip the 
configure step then. And some - usualy biological - packages are just weird.
Then you have to look for information in ``README.*`` or ``INSTALL.txt``.

Let's try to install two packages - ``Pipe viewer`` and ``vcflib``.
Pipe viewer is a nice tool you can use to watch the progress of your operations. 
It is distributed in standard ``.tar.gz`` form. Vcflib lives at GitHub - this is where
a lot of current open source software resides nowadays.

Pipe viewer: go to google, enter ``pipe viewer``. Click the ``ivarch.com`` link.
Look for downloads. Right click the ``pv-1.5.2.tar.gz`` link, select ``copy address`` 
or something similar (depends on your browser). Go to PuTTY, type ``wget -O -<space>``
and right-click your mouse. Then type `` | tar xvz <enter>``

.. code-block:: bash

    # go to a directory with software
    cd ~/sw

    # this is a spoiler, you should create the first line yourself
    # and do not copy it here
    wget -O - http://www.ivarch.com/programs/sources/pv-1.5.2.tar.gz | tar xvz
    cd pv<tab>
    ./configure
    make
    sudo make install

    # test pv
    </dev/zero pv > /dev/null


vcflib: go to google, type ``vcflib``, choose the GitHub link. Find ``clone url``.
Click the clipboard button. Go to PuTTY, type ``git clone<space>`` and right-click your mouse
in PuTTY window. Press <enter>.

.. code-block:: bash

    cd ~/sw
    git clone --recursive https://github.com/ekg/vcflib.git
    cd vcflib
    make

vcflib does not have ``make install``. We need to copy the binaries to ``$PATH``
manually.

.. code-block:: bash

  sudo cp bin/* /usr/local/bin

