Unix - Basics
=============

This session will give you all the basics that you need
to smoothly move around when using a UNIX system (in the text mode!).

Basic orientation - begining
----------------------------

Check your keyboard
^^^^^^^^^^^^^^^^^^^

Before we do any serious typing, make sure you know where are the important keys.
I'd suggest using English keyboard, if you don't want to constantly press right alt
and five random letters before you find the one you need.
You will definitely need those keys::

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

Directory structure
^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  ls -shaR # list content of a directory
  du -sh # disc usage (by directory)
  df -h # disc free
  ls | wc -l # what does this command do?
  locate # find a file/program

Moving or copying files and directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  touch # make a file
  mkdir -p # make a directory (``-p`` makes missing directory above)
  rm -r # remove a file/directory
  mv # move a file/directory
  cp -r # copy a file/directory
  ln -s # make a symbolic link

Prepare data directory in your HOME directory
and copy FASTQ data from common repository:

.. code-block:: bash

  cd ~
  mkdir -p data/fastq
  sudo cp -r /data/fastq/fastq.tar.gz data/fastq/.
  cd data/fastq
  ls

.. note::

   Normal users cannot change (and break) the (UNIX) system. There is one special
   user in each system called ``root``, who has the rights to make system wide changes.
   You can either directly log in as root, or use ``sudo`` (super user do) to execute
   one command as ``root``.

   .. image:: _static/sandwich.png
      :align: center

Uncompressing files
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  # Compressed tarball archives
  tar -xzvf fastq.tar.gz

  # gzipped files
  gunzip file.txt.gz

Viewing plain text file content
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  less -SN
  tail -n 5
  head -n 5
  cat
  nano

Try these commands:

.. code-block:: bash

  cd ~
  head -n 8 HRTMUOC01.RL12.00.fastq
  tail -n 8 HRTMUOC01.RL12.00.fastq

Pipes
^^^^^

Using the ``|`` (pipe) character you instruct the shell to take the output of the first command
and use it as an input for the second command.

The complement to ``head`` is ``tail``. It displays last lines of the input.
It can be readily combined with ``head`` to show the second sequence in the file.

.. code-block:: bash

    head -8 HRTMUOC01.RL12.00.fastq | tail -4 | less -S

    # or the third sequence data ;)
    < HRTMUOC01.RL12.00.fastq head -8 | tail -4 | less -S

**Exercise (How many reads are there?)**:

We found out that FASTQ files have a particular structure (four lines per read).
To find the total number of reads in our data, we will use another tool, ``wc``
(stands for `word count`, not for a toilet at the end of the pipeline;). ``wc``
counts words, lines and characters.

Our data is in three separate files. To merge them on the fly we'll use another tool,
``cat`` (for conCATenate). ``cat`` takes a list of file names and outputs a continuous
stream of the data that was in the files (there is no way to tell where one file ends
from the stream).

# now double click on each file name in the listing,
# and click right mouse button to paste (insert space in between)

.. code-block:: bash

  cat HRTMUOC01.RL12.00.fastq | wc -l

The number that appeared is four times the number of sequences (each sequence takes
four lines). And there is even a built-in calculator in bash:

.. code-block:: bash

  echo $(( 788640 / 4 ))
  expr XXXX / 4

Globbing
^^^^^^^^

Imagine you've got 40 FASTQ files instead of 3. You don't want to copy and paste all
the names! There is a feature that comes to rescue. It's called `globbing`. It allows
you to specify more filenames at once by defining some common pattern. All your
read files have ``.fastq`` extension. ``*.fastq`` means *a file named by any number of
characters followed by '.fastq'*.

.. code-block:: bash

  cat HRTMUOC01.RL12.*.fastq | wc -l
  expr XXXX / 4

  cat HRTMUOC01.RL12.0?.fastq | wc -l
  expr XXXX / 4

Producing lists
^^^^^^^^^^^^^^^

What do these commands do?

.. code-block:: bash

  touch file-0{1..9}.txt file-{10..20}.txt
  touch 0{1..9}-{a..f}.txt {10..12}-{a..f}.txt
  touch 0{1..9}-{jan,feb,mar}.txt {10..12}-{jan,feb,mar}.txt

**Exercise**:

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
^^^^^^^^^^^^^^^^^

.. code-block:: bash

  CPU=4
  echo $CPU

  FILE=data/fastq/HRTMUOC01.RL12.00.fastq
  echo $FILE

  FILES=`ls ~/data/fastq/*.fastq`
  echo $FILES

.. code-block:: bash

  list=`ls ~/data/fastq/HRTMUOC01.RL12.0{1..9}.fastq`

  for i in $list
  do
    echo $i
  done

  for i in $list
  do
    head -n1 $i | wc -c
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

Installing software
-------------------

The easiest way to install software is via a package manager (eg. ``apt-get`` for all Debian
variants). When the required software is not in the repositories, or one needs the latest
version, it's necessary to take the more difficult path. The canonical UNIX way is::

  wget -O - ..url.. | tar xvz   # download and unpack the 'tarball' from internet
  cd ..unpacked directory..     # set working directory to the project directory
  ./configure                   # check your system and choose the way to build it
  make                          # convert source code to machine code (compile it)
  sudo make install             # copy the results to your system

htop
^^^^

Installing software from common repository:

.. code-block:: bash

  sudo apt-get install htop

Bedtools
^^^^^^^^

Install software which is not in the common repository. You just need to find
a source code and compile it:

.. code-block:: bash

  wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
  tar -zxvf bedtools-2.25.0.tar.gz
  cd bedtools2
  make


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

Exercise
--------

.. note::

  1. What is the output of this command ``cd ~/ && ls | wc -l``?

    a) The total count of files in subdirectories in home directory
    b) The count of lines in files in home directory
    c) The count of files/directories in home directory
    d) The count of files/directories in current directory

  2. How many directories this command ``mkdir {1999-2001}-{1st,2nd,3rd,4th}-{1..5}`` makes?

    a) 56
    b) 60
    c) 64
    d) 72

  3. When files created using this command ``touch file0{1..9}.txt file{10..30}.txt``, how many files matched by ``ls file?.txt`` and ``ls file*0.txt``

    a) 30 and 0
    b) 0 and 30
    c) 30 and 4
    d) 0 and 3

  4. Which file would match this pattern ``ls *0?0.*``?

    a) file36500.tab
    b) file456030
    c) 5460230.txt
    d) 456000.tab

  5. Where do we get with this command ``cd ~/ && cd ../..``?

    a) two levels below home directory
    b) one level above home directory
    c) to root directory
    d) two levels above root directory

  6. What number does this command ``< file.txt head -n10 | tail -n+9 | wc -l`` print? (Assume the file.txt is not empty)

    a) 0
    b) 1
    c) 2
    d) 3
