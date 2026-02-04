A few more tips
===============

This is a collection of tips, that may help to overcome the initial barrier of
working with a 'foreign' system. There is a lot of ways to achieve the
solution, those presented here are not the only correct ones, but some that
proved beneficial to the authors.

Easiest ways to get Unix
------------------------

To get the most basic Unix tools, you can download an install  `Git for
Windows <https://git-scm.com/download/win>`_. It comes with a nice terminal
emulator, and installs to your right-click menu as 'Git Bash here' - which
runs terminal in the folder that you clicked. Git itself is meant for managing
versions of directories, but it cannot live without the Unix environment, so
someone did the hard work and packaged it all nicely together.

If you need more complete Unix environment, there are currently several options.
If you have a recent version of Windows 10 (yes, there are different versions
of Window 10), you can enable 'Windows Subsystem for Linux (WSL)' and then install
Ubuntu or Debian from the Windows Store. It's a marvel of engineering to connect
two operating systems, you're getting a 'real' Linux.


In older Windows,  you can use `Cygwin <https://www.cygwin.com/>`_. It is quite
complete, but it can't replace native Unix, as you'll find sooner or later.

Another easy way of getting a Unix environment in Windows is to install a basic Linux
into a virtual machine. Our previous courses used this method. It's much more
convenient that the dual boot configurations, and the risk of completely
breaking your computer is lower. You can be using Unix while having all your
familiar stuff at hand. The only downside is that you have to transfer all the
data as if the image was a remote machine.

It's much more convenient to use a normal terminal like PuTTY to connect to
the  machine rather than typing the commands into the virtual screen of
VirtualBox - It's usually lacking clipboard support, you cannot change the size
of the window, etc.

Mac OS X and Linux are Unix based, you just have to know how to start your
terminal program (``konsole``, ``xterm``, ``Terminal`` ...).

Essentials
----------
Always use ``screen`` for any remote work. Not using screen will cause your
jobs being interrupted when the network link fails (given you're working remotely),
and it will make you keep your home computer running even if your calculation is running
on a remote server.

Track system resources usage with ``htop``. System that is running low on memory won't
perform fast. System with many cores where only one core ('CPU') is used should be used for
more tasks - or can finish your task much faster, if used correctly.

Data organization
-----------------
Make a new directory for each project. Put all your data into subdirectories. Use
symbolic links to reference huge data that are reused by more projects in your current
project directory.
Prefix your directory names with 2 digit numbers, if your projects have more than few
subdirectories. Increase the number as the data inside is more and more 'processed'.
Keep the code in the top directory. It is easy to distinguish data references just by
having ``[0-9]{2}-`` prefix.

Example of genomic pipeline data directory follows:

.. code::

    00-raw --> /data/slavici/all-reads
    01-fastqc
    10-trim-adaptors
    13-fastqc
    20-assembly-newbler
    30-map-reads-gmap
    31-map-reads-bwa
    50-variants-samtools

Take care to note all the code used to produce all the intermediate data files.
This has two benefits:
1) your results will be really **reproducible**
2) it will **save you much work** when doing the same again, or trying different settings

If you feel geeky, use ``git`` to track your code files. It will save you from having 20 versions
of one script - and you being completely lost a year later, when trying to figure out which one
was the one that was actually working.

Building command lines
----------------------
Build the pipelines command by command, keeping ``| less -S`` (or ``| head`` if you don't expect lines
of the output to be longer than your terminal width) at the end. Every time you check if the
output is what you expect, and only after that add the next command. If there is a ``sort`` in
your pipeline, you have to put ``head`` in front of the ``sort``, because otherwise sort has to process
all the data before it gives out any output.

I (Libor) do prefer the 'input first' syntax (``<file command | command2 | command3
>out``) which improves legibility, resembles the real world pipeline (garden
hose, input tap -> garden hose -> garden sprinkler) more, and when changing
the input file names when reusing the pipeline, the names are easier to find.

Wrap your long pipelines on ``|`` - copy and paste to bash still works, because bash knows there
has to be something after ``|`` at the end of the line. Only the last line has to be escaped with ``\``,
otherwise all your output would go to the screen instead of a file.

.. code-block:: bash

  <infile sort -k3,3 |
    uniq -c -s64 |
    sort -k1rn,1 \
  >out

You can get a nice progress bar if you use ``pv`` (pipe viewer) instead of ``cat`` at the beginning
of the pipeline. But again, if there is a ``sort`` in your pipeline, it has to consume all the data
before it starts to work.

Use variables instead of hard-coded file names / arguments, especially when the name is used more times
in the process, or the argument is supposed to be tuned:

.. code-block:: bash

  FILE=/data/00-reads/GS60IET02.RL1.fastq
  THRESHOLD=300

  # count sequences in file
  <$FILE awk '(NR % 4 == 2)' | wc -l
  # 42308

  # count sequences longer that
  <$FILE awk '(NR % 4 == 2 && length($0) > $THRESHOLD)' | wc -l
  # 14190


Parallelization
---------------

Many tasks, especially in Big Data and NGS, are 'data parallel' - that means
you can split the data in pieces, compute the results on each piece separately
and then combine the results to get the complete result. This makes very easy
to exploit the full power of modern multi core machines, speeding up your
processing e.g. 10 times. ``GNU parallel`` is a nice tool that helps to
parallelize bash pipelines, check the manual for some examples: ``man
parallel_tutorial``.
