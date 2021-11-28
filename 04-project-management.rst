Session 4: Unix - "Real" Project
================================

Here you will learn several concepts needed for running real projects:

- use one screen session per project (if you're juggling more projects)
- keep your code in a text file
- keep your data files separate from your code
- make your intent more legible by wrapping code in functions and script files
- make your research shareable and reproducible via ``git``
- write handy scripts in ``awk``
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
  # screen sessions can be named, which helps you
  screen -S advanced

  # prepare source data
  cp /data-shared/bed_examples/Ensembl.NCBIM37.67.bed data/

  # add a new window with ctrl+a c

  # note and save all the (working) code
  nano workflow.sh

.. note::
  The nano bottom bar has a lot of weird symbols. ``^O`` means type ``ctrl-o``,
  pressing both of the keys at once. ``M-U`` means ``alt-u``, again pressing both
  simultaneously. For some reason, ``Write Out`` means **Save** in usual parlance.

Now you'll be typing your code in the nano window and pasting it to the
other window where the shell is running. You'll be writing code on your own
now, so you need to keep track of what you created. I'll occasionally post
solutions to `Slack <https://ngs-course.slack.com/>`_.

awk (pronounced [auk])
----------------------

``awk`` is most often used instead of ``cut``, when the fields are separated
by spaces and padded to a fixed width ``awk`` can ignore the whitespace -
and where ``cut`` also falls short, ``awk`` can reorder the columns:

.. topic:: Hands on!

  .. code-block:: bash

    # test the smart <tab> auto-complete!
    INPUT=data/Ensembl.NCBIM37.67.bed

    # $INPUT contains some genome annotations
    # look around the file a bit
    # - there are chromosome ids in the first column
    # - we want to count the annotations per each chromosome

    <$INPUT cut -f # complete the rest!!

But what if you wanted to have a table which mentions chromosomes first and
then the counts? Enter ``awk``. Every line (called ``record``) is split
into ``fields``, which are assigned to variables ``$1`` for first field,
``$2`` for second etc. The whole line is in ``$0`` if needed. ``awk`` expects
the program as first argument:

.. topic:: Hands on!

  .. code-block:: bash

    # note the single quotes, they're important because of $
    your_code_here | awk '{print $2, $1}'

Additionally you can add conditions when the code is executed:

.. topic:: Hands on!

  .. code-block:: bash

    your_code_here | awk '($2 < 10) {print $2, $1}'

Or even leave out the code body, which is the same as one ``{print $0}``
statement - that is print the matching lines:

.. topic:: Hands on!

  .. code-block:: bash

    your_code_here | awk '($2 < 10)'

There are some other variables pre-filled for each line, like
record number ``NR`` (starting at 1) and number of fields ``NF``.

.. code-block:: bash

  # NF comes handy when checking if it's okay to
  # process a file with (say) cut
  <$INPUT awk '{print NF}' | uniq

Let's play with some fastq files. Extract first five files to ``data``:

.. code-block:: bash

  INPUT=/data-shared/fastq/fastq.tar.gz
  <$INPUT tar tz | head -5 | xargs tar xvf $INPUT -C data

Look at the data with ``less`` - these are reads from 454, with varying read lengths.
Let's check the lengths:

.. code-block:: bash

  <data/HRTMUOC01.RL12.01.fastq paste - - - - | awk '{print $1, length($2)}' | head

We could do a length histogram easily now... But let's filter on the length:

.. topic:: Hands on!

  .. code-block:: bash

    <data/HRTMUOC01.RL12.01.fastq paste - - - - | # can you figure out?

    # and we'd like to have a valid fastq file on the output
    # - what if we replaced all the \t with \n (hint: tr)

Functions in the Shell
----------------------

This creates a command called ``uniqt`` that will behave as ``uniq -c``, but
there will be no padding (spaces) in front of the numbers, and numbers will be
separated by <tab>, so you can use it with ``cut`` will work.

.. code-block:: bash

  uniqt() { uniq -c | sed -r 's/^ *([0-9]+) /\1\t/' ;}

Now test it::

  <data/Ensembl.NCBIM37.67.bed cut -f1 | sort | uniqt | head

You can see that the basics of the syntax are ``your-name() { command pipeline ;}``.
If you want to pass some arguments into the function, use ``$1``, ``$2`` etc.::

  test-function() { echo First argument: $1 ;}
  test-function my-argument

Now create a function called ``fastq-min-length``, with one argument
(use ``$1`` in the body of the function) giving the minimal length:

.. topic:: Hands on!

  .. code-block:: bash

    fastq-min-length() { paste - - - - | your_code_here ;}

    # which will be used like this:
    <data/HRTMUOC01.RL12.01.fastq fastq-min-length 90 > data/filtered.fastq

We'll go through the 'quoting hell' and some methods to solve it here briefly.
Awk uses ``$1`` for something else than the shell, we need to protect it with
single quotes, but we still need to get through shell's ``$1`` somehow...
Awk's ``-v`` argument helps in this case - use it like ``awk -v min_len=$1
'(length($2) > min_len)'``.

.. note::

   Let's pop-open the matryoshka. What is terminal, what is a shell, what is
   Bash?

   The program which takes care of collecting your keystrokes and rendering
   the colored characters which come from the server is called a **terminal**.
   Famous terminals are ``mintty`` (that's what you're using in Windows now),
   ``Konsole``, ``Terminal App``... The next doll inside is ``ssh``. It takes
   care of encrypted communication with the remote server. An interesting
   alternative for geeks is ``mosh`` (google it yourself;). Now you need a
   program to talk to on the remote side - that is the **shell**. We're using
   ``bash`` now, sometimes you can meet the simpler cousin ``sh``, and the kool
   kids are doing ``zsh``. To recap, Bash is to shell what Firefox is to
   browser.

Shell Scripts
-------------
Another way to organize your code is to put it into a separate file
called a 'script file'. It begins with a ``shebang`` line, telling the computer
which language is the script in. Bash shebang is ``#! /bin/bash``.
Take care to give a descriptive name to your script::

    nano fastq-filter-length.sh

Copy and paste the following code block into the nano editor, save it with ``ctrl+o``
and switch to another bash window in screen.

.. topic:: Hands on!

  .. code-block:: bash

    #!/bin/bash

    # your_code_here

    echo Replace me with real code!
    echo Arguments: $1 $2

    # to stay with the 'tool concept'
    # expect input on stdin and output the results to stdout

We need to mark the file as executable and test it:

.. code-block:: bash

    chmod +x fastq-filter-length.sh

    # check with ls, filter_fastq.sh should be green now
    # and using ll you should see the 'x' (eXecutable) permission
    ls
    ll

    # and run it (the ./ is important!)
    ./fastq-filter-length.sh

.. note::

  You can check file permissions by typing ``ll`` instead of ``ls``.
  ``rwx`` stand for *Read*, *Write*, *eXecute*, and are repeated three times,
  for *User*, *Group*, and *Others*. The two names you see next to the
  permissions are file's owner user and group.

  You can change the permissions - if you have the permission to do so -
  by e.g. ``chmod go+w`` - "add write permission to group and others".

Now collect your code from above (contents of your function, not the whole
function) and paste it below the shebang. Don't forget to remove the debug ``echo``
parts - otherwise your script will spoil it's output with some useless chatter.

.. code-block:: bash

    # when the final code is there, you need to give it input (and maybe save the output):
    <data/HRTMUOC01.RL12.01.fastq ./fastq-filter-length.sh 90 > data/filtered.fastq


Code management and sharing via GIT
-----------------------------------
Once you start writing your scripts, you'll soon find yourself handling files
named like ``script.sh``, ``script_bak.sh``, ``script_previous_previous.sh``,
``script_last_working.sh`` .. etc. There is one cure for them all: ``git``.
It was originally created by Linus Torvalds, the author of the Linux kernel,
for managing the source code of Linux, it slowly gained popularity in other
communities.

What ``git`` does is managing versions of a directory tree. The managed subtree
is called a **repository**. Each saved version is called a **commit**.
You usually create a commit when you got a working version of your code. By
allowing you to go back to any committed version git effectively removes the
need for all ``_previous_working.sh`` copies of your code.

.. note::
  Please note the difference between ``git`` and Google Docs - Google Docs keeps
  track of all versions of a particular file. ``git`` keeps track of manually
  selected versions (snapshots) of a whole directory. That makes sense when
  ``script1.sh`` calls ``script2.sh``, and they have to match.

You will be using ``git`` to hand in the final exam, so please take care to set
it up correctly. Once on every new machine you need to tell git who you are,
because the commits are 'signed' by the author.

.. code-block:: bash

  git config --global user.email "you@example.com"
  git config --global user.name "Your Name"

To be able to share your code, you need an account on a 'social coding' site
like `GitHub <https://github.com>`_. If you don't have a GitHub account, please
get one. You don't have to do this right now, but if you believe could use some
help, now it's the best time.

To upload your code, you need to add a key to GitHub
`here <https://github.com/settings/keys>`_. This is how you generate the needed
key:

.. code-block:: bash

  # generate the private and public keys
  ssh-keygen -t ed25519

  # show the public key to copy and paste to github
  cat ~/.ssh/id_ed25519.pub

And now after a tedious setup let's reap the benefits. We'll store the current
version of your scripts in ``unix-advanced`` project, ignoring the data.

.. code-block:: bash

  # make sure we're in ~/projects/unix-advanced
  pwd

  # tell git we want to track this directory
  git init

  # tell git that we don't want to track and store the huge data
  # (git is not good at storing big data)
  echo 'data*' >> .gitignore

  # check what git sees in our brand new repo
  git status

  # add a single file
  git add workflow.sh

  # make a commit
  # add some descriptive message
  git commit -m 'solution to ngs-course exercises"

This is all you need to track your versions locally. If you want to publish your
creations, make a backup for yourself, or move your code to some shared machine
which will do bigger computations, you need to **push** it somewhere. If you
already have the GitHub account, you can create a repo on the website, add it as
a **remote** to your local repo and `push`. You don't have to be shy, GitHub allows
you to create a private repo.

.. code-block:: bash

  # use the commands suggested by GitHub to add a remote
  # ...

  # then push
  git push

When using git, you can gradually learn about more concepts and commands, as
you find the need for them. To give you a head start:

 - ``git pull`` updates your local repo if the remote is newer
 - by pulling other's changes over yours, you'll soon encounter **merge**
 - ``git stash`` can be used to "hide" local changes during ``pull``
   to avoid a commit and following merge, ``git stash pop`` brings them back
 - ``git checkout -b new-name`` and ``git branch some-name`` allow you to
   keep more simultaneous versions in one repo and switch between them

Multi-file, multi-core processing
---------------------------------
Multi-file processing is best done with ``find`` and ``xargs``. That's basic
Unix. If you install ``parallel``, it substitutes ``xargs`` and does much
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
  '*try - monitor performance - try again*' cycles.

  To get good performance it is important to know what happens during data processing:
  First the data is loaded from hard drive to memory, then from memory to the CPU,
  the CPU does the calculation, then the results have to get to the memory and saved
  to the hard drive again. Different workloads take different amounts of time in each step.

  .. image:: _static/data-processing.png
     :align: center

  In general, you need a work unit which takes much longer to calculate than
  it takes to load the data from the hard drive (compare times of ``pv data >
  /dev/null`` to ``pv data | your-task > /dev/null``), usually  a good work
  unit takes on the order of minutes. When disk access seems to be  the
  limiting factor, you can try to compress the data with some fast compressor
  like ``lz4``. **Do not** parallelize disk intensive tasks, it will make
  things only  slower! If you still want to use ``parallel``'s syntax, use
  ``parallel -j1`` to use only single core.

The most powerful thing about parallel is it's substitution strings like
``{.}``, ``{/}``, ``{#}`` - check ``man parallel``.

.. code-block:: bash

  parallel echo Ahoj ::: A B C

  parallel --dry-run echo Ahoj ::: A B C

  parallel echo File: {} found! ::: data/*.fastq

  parallel echo File: {/} found! ::: data/*.fastq

  parallel echo File: {/.} found! ::: data/*.fastq

.. note::

  If your data is a single file, but the processing of one line is not
  dependent on the other lines, you can use the ``split`` command to create
  several files each with defined number of lines from the original file.
