Unix - Advanced II
==================
You will learn to:

- keep your code in a file and data in a different place
- write handy scripts in ``awk``
- tidy your code by using functions and script files
- scale to multiple files and speed up your processing

Keep your stuff ordered
-----------------------
Let's pretend we're starting to work on something serious, a new project::

  cd ~/projects

  # a new project dir, separate dir for data in one shot
  mkdir -p unix-advanced/data

  # set your work dir to the 'current' project dir
  cd unix-advanced

  # definitely run screen to be safe
  screen -S advanced-2

  # prepare source data
  cp /data-shared/bed_examples/Ensembl.NCBIM37.67.bed data/

  # add a new window with ctrl+a c

  # note and save all the (working) code
  nano workflow.sh

Now you'll be typing your code in the nano window and pasting it to the
other window where the shell is running. You'll be writing code on your own
now, so you need to keep track of what you created. I'll occasionally post
solutions to `Slack <https://ngs-course.slack.com/>`_.

awk (pronounced [auk])
----------------------

``awk`` is most often used instead of ``cut``, when the fields are separated
by spaces and padded to a fixed width ``awk`` can ignore the whitespace -
and where ``cut`` also falls short, ``awk`` can reorder the columns::

.. code-block:: bash

  # test the smart <tab> auto-complete!
  INPUT=data/Ensembl.NCBIM37.67.bed

  # $INPUT contains some genome annotations
  # look aronund the file a bit
  # - there are chromosome ids in the first column
  # - we want to count the annotations per each chromosome

  <$INPUT cut -f # complete the rest!!

But what if you wanted to have a table which mentions chromosomes first and
then the counts? Enter ``awk``. Every line (called ``record``) is split
into ``fields``, which are assigned to variables ``$1`` for first field,
``$2`` for second etc. The whole line is in ``$0`` if needed. ``awk`` expects
the program as first argument::

.. code-block:: bash

  # note the single quotes, they're important because of $
  your code here | awk '{print $2, $1}'

Additionally you can add conditions when the code is executed::

  .. | awk '($2 < 10) {print $2, $1}'

Or even leave out the code body, which is the same as one ``{print $0}``
statement - that is print the matching lines::

  .. | awk '($2 < 10)'

There are some other variables pre-filled for each line, like
record number ``NR`` (starting at 1) and number of fields ``NF``.

.. code-block:: bash

  # NF comes handy when checking if it's okay to
  # process a file with (say) cut
  <$INPUT awk '{print NF}' | uniq


Let's play with some fastq files. Extract first five files to ``data``::

  INPUT=/data-shared/fastq/fastq.tar.gz
  <$INPUT tar tz | head -5 | xargs tar xvf $INPUT -C data

Look at the data with ``less`` - these are reads from 454, with varying read lengths.
Let's check the lengths::

  <data/HRTMUOC01.RL12.01.fastq paste - - - - | awk '{print $1, length($2)}' | head

We could do a length histogram easily now... But let's filter on the length::

  <data/HRTMUOC01.RL12.01.fastq paste - - - - | # can you figure out?

  # and we'd like to have a valid fastq file on the output
  # - what if we replaced all the \t with \n (hint: tr)

Functions in the Shell
----------------------
Create a command ``uniqt`` that will behave as ``uniq -c``, but there
will be no padding (spaces) in front of the numbers, and numbers will
be separated by <tab>, so e. g. ``cut`` will work.

Do not use the same name as the original command, otherwise you'll create
an endless loop.

.. code-block:: bash

    uniqt() { uniq -c | sed -r 's/^ *([0-9]+) /\1\t/'  ;}

You can see that the basics of the syntax are ``your-name() { normal commands ;}``.
What about creating a function called ``fastq-min-length``, with one argument
(use ``$1`` in the body of the function) giving the minimal length::

  fastq-min-length() { paste - - - - | ... ;}

We'll go through the 'quoting hell' and some methods to solve it here briefly.
Awk uses ``$1`` for something else than the shell, we need to protect it with
single quotes, but we still need to get through shell's ``$1`` somehow...
Awk's ``-v`` argument helps in this case.


Shell Scripts
-------------
Another way to organize your code is to put it into a separate file
called a 'script file'. It begins with a ``shebang`` line, telling the computer
which language is the script in. Bash shebang is ``#! /bin/bash``.
Take care to give a descriptive name to your script::

    nano fastq-filter-length.sh

.. note::

   Let's pop-open the matryoshka. What is terminal, what is a shell, what is
   Bash?

   The program which takes care of collecting your keystrokes and
   rendering the colored characters which come from the server is called a
   terminal. Famous terminals are ``mintty`` (that's what you're using in
   Windows now), ``Konsole``, ``Terminal App``... The next doll inside is
   ``ssh``. It takes care of encrypted communication with the remote server.
   An interesting alternative for geeks is ``mosh`` (google it yourself;). Now
   you need a program to talk to on the other side - that is the shell,
   running on the remote side. We're in ``bash``, sometimes you can meet the
   simpler cousin ``sh``, and the kool kids are doing ``zsh``. To recap, Bash
   is to shell what Firefox is to browser.

Then collect your code from before and paste it below the shebang.

.. code-block:: bash

    #!/bin/bash

    # your code comes here
    # to stay with the 'tool concept' output the results to stdout

We need to mark the file as executable:

.. code-block:: bash

    chmod +x filter_fastq.sh

    # check with ls, filter_fastq.sh should be green now
    # and using ll you should see the 'x' (eXecutable) permission
    ls
    ll


Multi- processing
-----------------
Multi-file processing is best done with ``find`` and ``xargs``. That's basic
UNIX. If you install ``parallel``, it substitutes ``xargs`` and does much
better job, having 'nicer' syntax, and makes multi-file multi-core processing
a breeze.

Let's check the basic concepts - ``find`` converts directory structure to
'data' (stdout), ``xargs`` converts stdin to command line(s).

.. code-block:: bash

  # Investigate!

  find data -type f

  find data -type f | xargs echo

  find data -type f | xargs -I{} echo File: {} found!

``parallel`` runs one instance of the command per each CPU in your machine.
Regrettably your **virtual** machine has only one CPU, so this won't help
much. But modern machines do have  four and more CPUs, and then it really
helps.

Do control the number of jobs (``-j``) only when sharing the machine with
someone, or when you're sure that your task is IO bound. Otherwise
``parallel`` does a good job choosing the number of tasks to run for you.

.. note::

  Parallelizing things **IS** difficult. There's no discussion about that.
  There are some rules of thumb, which can help - but if you want to squeeze
  out the maximum performance from your machine, it's still a lot of
  'try-monitor performance-try again' cycles.

  In general, you need a work unit which takes much longer to calculate than
  it takes to load the data from the hard drive (compare times of
  ``pv data > /dev/null`` to ``pv data | your-task > /dev/null``) ... TODO

.. code-block:: bash

  # some example here

There is a lot of magic to be done with ``{.}, {/}, {#}`` placeholders,
check ``man parallel``. If your data is a single file, but the processing
of one line is not dependent on the other lines, ``split`` will help.
