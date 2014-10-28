UNIX primer
===========
You need to connect to the system, type in commands, keep stuff running
while you're disconnected and then pick up your results.

Connecting to the system
------------------------
To connect from MS Windows we'll use PuTTY. PuTTY is a simple window that 
sends everything you type to the remote computer and displays anything the 
remote computer sends back. 

.. warning:: Clipboard works differently in PuTTY! When you select text
    PuTTY assumes you want to copy it - so it is automatically copied to 
    clipboard. To paste text you need to press right mouse button.

.. warning:: Trying to use windows shortcuts - especially Ctrl-C kills 
    your current running program.

Run putty and enter following information::

  Host: localhost
  Port: 2222

The shell
---------
What you see now is the shell. Shell is a program for entering commands.
Your shell is *bash*. Bash is to shell what MS Word is to text editor.
You can choose your shell if you need, but most people use bash.

Multiple windows
----------------
You're all used to work with multiple windows (in MS Windows;). You can 
have them in (remote) Linux as well.

.. code-block:: bash

  screen

.. note:: The additional benefit is that you can log off, and your programs 
  keep running.

``Screen`` is controled after you press the master key - ``ctrl-a``. The next key you 
press is a command to ``screen``.

To create a new window, press ``ctrl-a c`` (create). To flip among your windows press ``ctrl-a space``
(you flip windows often, it's the biggest key available).
To detach ``screen`` - "keep your programs running and go home" - press ``ctrl-a d`` (detach).

Coming back to work you need to connect to your ``screen`` (-r is for restore).

.. code-block:: bash

  screen -r

Moving around
-------------
You need to type these commands in bash - keep your eye on the ``prompt``
- the beginning of the line where you type. Different programs present different 
prompts.

.. code-block:: bash

	pwd    # prints current directory path
	cd     # changes current directory path
	ls     # lists current directory contents
	ll     # lists detailed contents of current directory
	mkdir  # creates a directory
	rm     # removes a file
	rm -r  # removes a directory
	cp     # copies a file/directory
	mv     # moves a file/directory
	locate # tries to find a file by name
	ln -s  # create symbolic link

Getting help
------------
Call me or Vaclav to get any help ;)

Once you know the name of the command that does what 
you need, all the details are easily accessible using ``man``.
To get all possible help about finding text do:

.. code-block:: bash

  man grep

To find the name of the command that does what you need, use google::

  linux search for string

Viewing files
-------------
``less`` is the command::

  less /data/slavici/00-reads/GSVZDOM02.fastq

Toggle line wrapping by typing ``-S<enter>``.

Search for sequence ``ACGT`` by typing ``/ACGT<enter>``. Press ``n`` (next) to 

Exit less by typing ``q``.


Chaining commands
-----------------
You know how to display whole file (``less``). What if you want
to display just specific information from the file?

Change the directory, so we don't have to type so much (press ``<tab>`` often
to spare some typing in bash)::

  cd /data/slavici

View only first 100 lines::

  <00-reads/GSVZDOM02.fastq head -100 | less

View only sequence names (they all start with @)::

  <00-reads/GSVZDOM02.fastq grep ^@ | less

We can see that the simple assumption was not correct - 
not only sequence names start with @. Let's display every fourth 
line - names are only on fourth lines::

  <00-reads/GSVZDOM02.fastq awk '(NR % 4 == 1)' | less

Writing to file instead of looking at it is easy::

  <00-reads/GSVZDOM02.fastq awk '(NR % 4 == 1)' > test-file

  # check if the data is there ;)
  less test-file

  # get rid of the file
  rm test-file

Chaining is not limited to two commands. I need first 1000 sequence names without
the @::

  <00-reads/GSVZDOM02.fastq awk '(NR % 4 == 1)' | cut -c2- | head -1000 > second-test
  less second-test

Getting out of it all
---------------------
.. code-block:: bash
	
	exit   # quits current session
