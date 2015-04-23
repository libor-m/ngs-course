Starter session
===============

This session will give you all the basics that you need 
to smoothly move around when using a UNIX system (in the text mode!).

Use multiple windows (and be safe when the network fails)
---------------------------------------------------------

First, type ``screen`` in your terminal::

  screen

Screen creates the first window for you. To create another one press 
``ctrl+a c``. To switch between the windows press ``ctrl+a space``.

Check what the computer is doing
--------------------------------

Run ``htop`` in one of your screen windows::

  htop

Htop displays CPU and memory utilization of the (virtual) computer. Continue your 
work in another window (``ctrl+a space``). You can switch back to the htop window to
monitor progress of some lengthy operation.

.. note:: 

   Keyboard shortcuts notation: `ctrl+a space` means press `ctrl` key and `a` letter
   simultaneously and `space` key after you release both of the previous keys.

Move around the directory structure
-----------------------------------

Unlike 'drives' in MS Windows, UNIX has a single directory tree 
that starts in `/` (called root). Everything can be reached from the root.
The next important directory is `~` (called user's home directory). It is 
a shortcut for `/home/user`.

Your bash session has a `working directory`. All filenames and paths you 
type refer to your working directory, unless you start them with `/`. 

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

Installing software
===================

blah

.. note:: 

   To paste text into PuTTY just click right mouse button anywhere in the window.
   To copy text to clipboard, just select it. No keyboard shortcuts are necessary.
