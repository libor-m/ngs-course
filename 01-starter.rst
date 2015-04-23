Starter session
===============

This session will give you all the basics that you need 
to smoothly move around when using a UNIX system (in the text mode!).

Basic orientation
^^^^^^^^^^^^^^^^^

Use multiple windows (and be safe when the network fails)
---------------------------------------------------------

First, type ``screen`` in your terminal::

  screen

Screen creates the first window for you. To create another one press 
``ctrl+a c``. To switch between the windows press ``ctrl+a space``.

.. note:: 

   Keyboard shortcuts notation: ``ctrl+a space`` means press ``ctrl`` key and ``a`` key
   simultaneously and ``space`` key after you release both of the previous keys.

Check what the computer is doing
--------------------------------

Run ``htop`` in one of your screen windows::

  htop

Htop displays CPU and memory utilization of the (virtual) computer. Continue your 
work in another window (``ctrl+a space``). You can switch back to the htop window to
monitor progress of some lengthy operation.

Move around the directory structure
-----------------------------------

Unlike 'drives' in MS Windows, UNIX has a single directory tree 
that starts in ``/`` (called root). Everything can be reached from the root.
The next important directory is ``~`` (called user's home directory). It is 
a shortcut for ``/home/user``.

Your bash session has a `working directory`. All filenames and paths you 
type refer to your working directory, unless you start them with ``/``. 

Try the following commands in the order they are provided, and figure out what they do.
Then use your knowledge to explore the directory structure of the virtual machine.

.. code-block:: bash

    pwd
    ls
    ls /
    cd /
    pwd
    ls
    cd
    pwd


A neat trick to go back where you've been before the last `cd` command::

  cd -

More in :ref:`moving_around`.

Prepare data in your home directory
-----------------------------------

All the data we are going to use are located in ``/data``. However, we want to have it 
to be conveniently located in our own directory (``~``, ``/home/user``). Without the need
to copy all of it, we can use a symbolic link (``ln -s``). This keeps the data in their original
location but creates a reference.

In case you do not know where your files are but you do know some part of the name,
use the ``locate`` command. We know the names contained fastq, vcf and gff3 suffices:

.. code-block:: bash

    locate fastq
    locate vcf
    locate gff3

Once we know their actual position we can create symbolic links:

.. code-block:: bash

    # create directory data
    # <- this marks a comment - anything after first # is ignored
    mkdir data 
    
    # go to your new data directory
    cd data 
    
    # create a link to the nightingale reads
    # and name it 'fastq'
    ln -s /data/slavici/00-reads fastq


You created a `symbolic link` named ``fastq`` with ``/data/slavici/00-reads`` as a `target`.
Check it by typing::

  ls -l

.. note:: 

   You should use bash `autocomplete` feature, when typing paths. It is easier, faster
   and less prone to error. Type a part of the path, like ``/da`` and press the ``tab``
   key. When nothing appears, press ``tab`` once more. There is either no possible completion
   or more possibilities, that will be displayed on the second press.

It is possible to create a bad link. There is no validation on the target::

  ln -s /nothing_here bad-link

  # the bad link has a different color in the output
  ls -l

  # get rid of the bad link
  rm bad-link


Installing software
^^^^^^^^^^^^^^^^^^^
The easiest way to install software is via a package manager (eg. ``apt-get`` for all Debian
variants). When the required software is not in the repositories, or one needs the latest
version, it's necessary to take the more difficult path. The canonical UNIX way is::

  wget -O - ..url.. | tar xvz   # download and unpack the 'tarball' from internet
  cd ..unpacked directory..     # set working directory to the project directory
  ./configure                   # check your system and choose the way to build it
  make && sudo make install     # convert source code to machine code and if successful, copy the results to your system

First we'll get the latest pipe viewer. Google ``pipe viewer``, choose the ivarch.com site. 
Check the current version number on the site. Now check the version in your image::

  pv --version

.. note::

   It is a good habit to include ``--version`` option for a command. You need to check 
   the version of given tool in your system when you're trying to use some new features.

The version found at the site should be higher then the one in your image. A good reason for 
update;) Copy the link for the ``.tar.bz2`` file on the site.

.. code-block:: bash

   # go to the directory where software installations live
   cd ~/sw

   wget -O - ..paste the link here .. | tar xvj

   # the complete command from above for those who are cheating
   wget -O - http://www.ivarch.com/programs/sources/pv-1.6.0.tar.bz2 | tar xvj

   # do not copy this, try the autocompletion
   # cd pv<tab> <tab> <6> <tab> <enter>

   ls
   # you can see green configure script in the listing

   # to run something in current directory, the path has
   # to be given
   ./configure
   make

   # to make changes system wide, super user 'powers' have to be used
   sudo make install


.. figure:: _static/sandwich.png

   Normal users cannot change (and break) the (UNIX) system. There is one special 
   user in each system called ``root``, who has the rights to make system wide changes.
   You can either directly log in as root, or use ``sudo`` (super user do) to execute
   one command as ``root``.


In our example, some steps are ommited. We'll install ``bedtools`` program from a github repository. 
User installed software can be found in ``~/sw`` directory. To install a new software go to this directory:

.. code-block:: bash

    

  When the software source code is in a single file (`tarball`), ``wget`` command is the best option to get the file.
  The latest versions are usually not packaged, and many of the tools can be found at GitHub. To get stuff from GitHub,
  ``git clone`` command is usually the easiest option.
 
  .. code-block:: bash

    git clone https://github.com/arq5x/bedtools2

  For those without internet access in their virtual machines - you need to download the content to your
  normal computer and then transfer it to the virtual machine.

    - download https://github.com/arq5x/bedtools2/archive/master.zip
    - transfer it to the virtual machine with WinSCP (Windows) or scp (Mac or Linux)
    - unpack the file with ``unzip``
    - rename the folder with ``mv`` to bedtools2

  This creates a `clone` of the online repository in directory ``bedtools2``.

  .. code-block:: bash

    cd bedtools2

  To compile (convert from text source form to machine executable form) software on UNIX use the ``make`` command:

  .. code-block:: bash

    make

  It should take a while, so you can flip to your `htop` window with ``ctrl-a space`` and watch the CPU spin;)

  When ``bedtools`` is compiled you have to copy bedtools binaries to ``/usr/local/bin`` directory for UNIX system to
  find the program when calling from any place in the system.

  .. warning:: Before you use command below to copy binaries make sure you are really in directory you want to be!
 
  .. code-block:: bash
    
    cd bedtools2/bin
    ls # Check that you are really in directory you want to be!
    sudo cp * /usr/local/bin

  We used two commands: ``sudo`` and ``cp``. The sudo command tells the system that we want to make changes in system
  directories and as such we are asked for password. This step prevents us from harming system. The ``cp`` command
  copies all bedtools binaries from local bin directory to the system binary repository.
  
.. note:: We used the ``*`` symbol which tells the system that all files in the current directory should be selected. We explain this later.


.. note:: 

   To paste text into PuTTY just click right mouse button anywhere in the window.
   To copy text to clipboard, just select it. No keyboard shortcuts are necessary.
