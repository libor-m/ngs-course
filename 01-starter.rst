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
that starts in ``/`` (called root directory). Everything can be reached from the root directory.
The next important directory is ``~`` (called user's home directory). It is 
a shortcut for ``/home/user`` here, ``/home/..your login name..`` in general.

Your bash session has a `working directory` that can be changed with ``cd`` (change directory) 
and printed with ``pwd`` (print working directory). All filenames and paths you 
type refer to your working directory (relative paths), unless you start them with ``/`` (absolute paths). 

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
^^^^^^^^^^^^^^^^^^^
The easiest way to install software is via a package manager (eg. ``apt-get`` for all Debian
variants). When the required software is not in the repositories, or one needs the latest
version, it's necessary to take the more difficult path. The canonical UNIX way is::

  wget -O - ..url.. | tar xvz   # download and unpack the 'tarball' from internet
  cd ..unpacked directory..     # set working directory to the project directory
  ./configure                   # check your system and choose the way to build it
  make && sudo make install     # convert source code to machine code and if successful, copy the results to your system

Pipe viewer
-----------
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
--------
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


Show me the data!
^^^^^^^^^^^^^^^^^
Until now we were working with files and directories. But the real data is
inside the files. 

Explore FASTQ files
-------------------

The ``less`` tool is used to list through contents of a text file.  We will check some 
of the FASTQ files linked in our ``~/data`` directory.

.. code-block:: bash

   # cd by itself means cd ~ (that is cd /home/user here)
   # this will get you to your home directory, wherever you are
   cd

   # a file can be referenced in various ways
   # option 1: absolute path (<q> to quit the viewer)
   less /home/user/data/fastq/G59B7NP01.fastq

   # option 2: relative path from working directory
   less data/fastq/G59B7NP01.fastq

   # option 3: move 'closer' to the file
   cd data/fastq
   less G59B7NP01.fastq

.. note:: Reminder: you don't have to type the whole file name. Try to use TAB auto-completion!

The data you see looks like mess. One of the reasons is there are long lines, that
get wrapped so you see all the letters. But then you don't see the file structure.
Add the ``-S`` option, and see the four different line types in the FASTQ file::

  less -S G59B7NP01.fastq

The lines are:

  1. sequence name
  2. dna letters
  3. ``+`` sign
  4. encoded quality scores

The options can be given either one by one - which is more legible, or combined. Another interesting
option is ``-N``, showing the line numbers::

  less -S -N G59B7NP01.fastq

  # this is the same as above
  less -SN G59B7NP01.fastq

.. note:: If you forgot to type ``-S`` at the prompt, you can type ``-S`` also while in ``less``. Try it!

UNIX Pipes
----------
For a quick glance over the contents of the file, you can also use the ``head`` command::

  head G59B7NP01.fastq

The problem with the wrapped lines comes back again. ``head`` is not meant to be a file viewer,
so it does not have any text wrapping options. Instead you can combine two tools. ``cut`` allows you 
to choose only a part of each line.

  .. code-block:: bash

    # show up to 50 characters from each 
    # of the first 10 lines in the file
    head G59B7NP01.fastq | cut -c 50

    # we can get only first four lines
    head -4 G59B7NP01.fastq | cut -c 50

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

How many reads are there?
-------------------------
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
  # and click right mouse button to paste (isert space in between)
  cat G59B7NP01.fastq GS60IET02.RL1.fastq GS60IET02.RL2.fastq | wc -l

The number that appeared is four times the number of sequences (each sequence takes 
four lines). And there is even a built-in calculator in bash::

  echo $(( 788640 / 4 ))

Imagine you've got 40 FASTQ files instead of 3. You don't want to copy and paste all
the names! There is a feature that comes to rescue. It's called `globbing`. It allows 
you to specify more filenames at once by defining some common pattern. All your 
read files have ``.fastq`` extension::

  echo *.fastq

``echo`` is no magic, it outputs whatever you give it (try ``echo ahoj``). The magic
is done by bash - whenever it sees an asterisk (``*``), it tries to expand it by 
matching to the files and directories. ``*.fastq`` means *a file named by any number of 
characters followed by '.fastq'*.

Globbing works even across directories, try::

  cd ..
  echo fastq/*.fastq

Now we can use it in our read counting pipeline to make it shorter and more versatile::

  cd fastq
  cat *.fastq | wc -l

How many bases were sequenced?
------------------------------
``wc`` can count characters (think bases) as well. But to get a reasonable number,
we have to get rid of the other lines that are not bases.

One way to do it is to pick only lines comprising of letters A, C, G, T and N.
There is a ubiquitous mini-language called `regular expressions` that can be used
to define text patterns. `A line comprising only of few possible letters` is 
a text pattern. ``grep`` is the basic tool for using regular expressions::

  cat *.fastq | grep '^[ACGTN]*$' | less -S

Check if the output looks as expected. This is a very common way to work - build a part of 
the pipeline, check the output with ``less`` or ``head`` and fix it or add more commands.

Now a short explanation of the ``^[ACGTN]*$`` pattern (``grep`` works one line a time):

- ``^`` marks beginning of the line - otherwise ``grep`` would search anywhere in the line
- the square brackets (``[]``) are a `character class`, meaning one character of the list, ``[Gg]rep`` 
  matches ``Grep`` and ``grep``
- the ``*`` is a count suffix for the square brackets, saying there should be zero or more of such characters
- ``$`` marks end of the line - that means the whole line has to match the pattern

To count the bases read, we extend our pipeline::

  cat *.fastq | grep '^[ACGTN]*$' | wc -c

The thing is that this count is not correct. ``wc -c`` counts every character, and the end of the line
is marked by a special character written as ``\n`` (n for newline). To get rid of this character,
we can use another tool, ``tr`` (transliterate). ``tr`` can substitute one letter with another 
(imagine you need to lowercase all your data, or mask lowercase bases in your Fasta file). Additionally
``tr -d`` (delete) can remove characters::

  cat *.fastq | grep '^[ACGTN]*$' | tr -d "\n" | wc -c

.. note::  If you like regular expressions, you can hone your skills at https://regex.alf.nu/.
