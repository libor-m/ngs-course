Session 1: Unix - Introduction
==============================

This session will give you all the basics that you need
to smoothly move around when using a UNIX system (in the text mode!).

Check your keyboard
-------------------
Before we do any serious typing, make sure you know where are the important keys.
I'd suggest using English keyboard, if you don't want to constantly press right alt
and five random letters before you find the one you need.
You will definitely need those keys:

.. code-block:: bash

  [] - squared brackets
  {} - curly brackets
  <> - angle brackets (smaller-than, bigger-than sign)
  () - parentheses
  ~ - tilde
  / - slash
  \ - back slash
  | - pipe
  ^ - caret
  $ - dollar sign
  : - colon
  ; - semicolon
  . - dot
  , - comma
  # - hash
  _ - underscore
  - - dash
  * - asterisk
  ! - exclamation mark
  ? - question mark
  & - ampersand
  @ - at sign
  '' - quotation mark single
  "" - quotation mark double

Basic orientation - directory structure
---------------------------------------

Unlike 'drives' in MS Windows, UNIX has a single directory tree
that starts in ``/`` (called root directory). Everything can be reached from the root directory.
The next important directory is ``~`` (called user's home directory). It is
a shortcut for ``/home/user`` here, ``/home/..your login name..`` in general.

Your bash session has a `working directory` that can be changed with ``cd`` (change directory)
and printed with ``pwd`` (print working directory). All filenames and paths you
type refer to your working directory (relative paths), unless you start them with ``/`` (absolute paths).

Try the following commands in the order they are provided, and figure out what they do.
Then use your knowledge to explore the directory structure of the virtual machine.

Figure out what these commands do:

.. code-block:: bash

    pwd
    ls
    ls /
    ls ..
    ls ~
    cd
    cd /
    cd ..
    cd ~

A neat trick to go back where you've been before the last ``cd`` command:

.. code-block:: bash

  cd -

More in :ref:`moving_around`.

Helpful commands (dir content and its size, disc usage)
-------------------------------------------------------

.. code-block:: bash

  ls -shaR
  du -sh
  df -h

Moving/copying files/directories
--------------------------------

.. code-block:: bash

  touch # make a file
  mkdir -p # make a directory (``-p`` makes missing directory above)
  rm -r # remove a file/directory
  mv # move a file/directory
  cp -r # copy a file/directory

Prepare data directory in your HOME directory
and copy FASTQ data from common repository:

.. code-block:: bash

  cd ~
  mkdir -p data/fastq
  sudo cp -r /data/fastq/fastq.tar.gz data/fastq/.
  ls data/fastq

Uncompressing files
-------------------

.. code-block:: bash

  tar -xzvf data/fastq/fastq.tar.gz

Viewing plain text file content
-------------------------------

.. code-block:: bash

  less -SN
  tail -n 5
  head -n 5
  cat
  nano

Exercise (What does these commands do?):

.. code-block:: bash

  cd ~
  head -n 8 data/fastq/HRTMUOC01.RL12.00.fastq
  tail -n 8 data/fastq/HRTMUOC01.RL12.00.fastq

Pipes
-----

Using the ``|`` (pipe) character you instruct the shell to take the output of the first command
and use it as an input for the second command. You can also use ``less`` as a part of the
pipeline::

  head -4 G59B7NP01.fastq | less -S

The complement to ``head`` is ``tail``. It displays last lines of the input.
It can be readily combined with ``head`` to show the second sequence in the file.

.. code-block:: bash

    head -8 G59B7NP01.fastq | tail -4 | less -S

    # or the third sequence data ;)
    head -12 G59B7NP01.fastq | tail -4 | less -S

**Exercise (How many reads are there?)**::

We found out that FASTQ files have a particular structure (four lines per read).
To find the total number of reads in our data, we will use another tool, ``wc``
(stands for `word count`, not for a toilet at the end of the pipeline;). ``wc``
counts words, lines and characters.

Our data is in three separate files. To merge them on the fly we'll use another tool,
``cat`` (for conCATenate). ``cat`` takes a list of file names and outputs a continuous
stream of the data that was in the files (there is no way to tell where one file ends
from the stream).

.. code-block:: bash

    ls

# now double click on each file name in the listing,
# and click right mouse button to paste (insert space in between)
cat G59B7NP01.fastq GS60IET02.RL1.fastq GS60IET02.RL2.fastq | wc -l

The number that appeared is four times the number of sequences (each sequence takes
four lines). And there is even a built-in calculator in bash::

  echo $(( 788640 / 4 ))
  expr 788640 / 4

Imagine you've got 40 FASTQ files instead of 3. You don't want to copy and paste all
the names! There is a feature that comes to rescue. It's called `globbing`. It allows
you to specify more filenames at once by defining some common pattern. All your
read files have ``.fastq`` extension::

echo *.fastq

``echo`` is no magic, it outputs whatever you give it (try ``echo ahoj``). The magic
is done by bash - whenever it sees an asterisk (``*``), it tries to expand it by
matching to the files and directories. ``*.fastq`` means *a file named by any number of
characters followed by '.fastq'*.

Globbing
--------

.. code-block:: bash

  ls *.fastq
  ls *.fast?

Producing list in Unix
----------------------

.. code-block:: bash

  touch dir-{1..12}
  touch dir-0{1..9} dir-{10..12}
  touch {2013..2015}-0{1..9} {2013..2015}-{10..12}
  touch {2013..2014}-0{1..9}-{a..c} {2013..2014}-{10..12}-{a..c}

**Exercise**::

Program runs 20 runs of simulations for three datasets (hm, ss, mm) using
three different sets of values: small (sm), medium sized (md) and large (lg).
There are three groups of output files, which should go into subdirectory A, B and C.
Make a directory for each dataset-set of parameters-run-subdirectory.
Count the number of directories.

Producing lists of subdirectories

.. code-block:: bash

  mkdir –p {2013..2015}/{A..C}
  mkdir –p {2013..2015}/0{1..9}/{A..C} {2013..2015}/{10..12}/{A..C}

Variables & Loops
-----------------

.. code-block:: bash

  list=`ls HRTMUOC01.RL12.0{1..9}.fastq`

  for i in $list
  do
    echo $i
  done

  for i in $list
  do
    head -n1 $i
  done

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

.. note::

  To find the name of the command that does what you need (``grep``), use google::

    linux search for string

  Once you know the name of the command that does what you need, all the
  details are easily accessible using ``man``. To get all possible help about
  finding text patterns in files do:

  .. code-block:: bash

    man grep

.. note::

  Some useful keyboard shortcuts::

    ctrl+c  - kills current running program (except for bash, nano, vim, ...)
            - clears the command line in bash

    ctrl+d  - means end of input (if you run e.g. bc interactively)
            - end of input means logout in bash

    ctrl+r  - starts history search in bash, just type a part of a long command
              and it will come back (ctrl+c to the rescue;)

    ctrl+k  - clears the command line from cursor to the end,
              you will need this while exploring long pipelines...


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

.. note::

   To paste text into PuTTY just click right mouse button anywhere in the window.
   To copy text to clipboard, just select it. No keyboard shortcuts are necessary.

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

It is possible to create a bad link. There is no validation on the target:

.. code-block:: bash

  ln -s /nothing_here bad-link

  # the bad link has a different color in the output
  ls -l

  # get rid of the bad link
  rm bad-link

Installing software
-------------------
The easiest way to install software is via a package manager (eg. ``apt-get`` for all Debian
variants). When the required software is not in the repositories, or one needs the latest
version, it's necessary to take the more difficult path. The canonical UNIX way is::

  wget -O - ..url.. | tar xvz   # download and unpack the 'tarball' from internet
  cd ..unpacked directory..     # set working directory to the project directory
  ./configure                   # check your system and choose the way to build it
  make && sudo make install     # convert source code to machine code and if successful, copy the results to your system

htop
^^^^

Pipe viewer
^^^^^^^^^^^
First we'll get the latest pipe viewer. Pipe viewer can show you how
much of the data was already processed in your `pipeline`. Google ``pipe viewer``,
choose the ivarch.com site. Check the current version number on the site.
Now check the version in your image::

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


.. note::

   Normal users cannot change (and break) the (UNIX) system. There is one special
   user in each system called ``root``, who has the rights to make system wide changes.
   You can either directly log in as root, or use ``sudo`` (super user do) to execute
   one command as ``root``.


   .. image:: _static/sandwich.png
      :align: center

Bedtools
^^^^^^^^
Another common place where you find a lot of software is `GitHub`. We'll install
``bedtools`` from a GitHub repository:

.. code-block:: bash

  cd ~/sw

  # get the latest bedtools
  git clone https://github.com/arq5x/bedtools2

This creates a `clone` of the online repository in ``bedtools2`` directory.

.. code-block:: bash

   cd bedtools2
   make

The compilation should take a while, so you can flip to your `htop` window with
``ctrl-a space`` and watch the CPU spin;)
