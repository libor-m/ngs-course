Connecting to the virtual machine
=================================
.. note:: 
  You need to start the virtual machine first!

.. _ssh_connect:

Connect to control the machine
------------------------------
To control the machine, you need to connect to the ssh service. 
This is also referred to as 'logging in'.

In Windows this is done with PuTTY.

- start PuTTY
- fill Host Name: ``localhost``
- fill Port: ``2222``
- click Open or press <Enter>

.. image:: _static/putty-config.png

In the black wnidow that appears, type your credentials:

- login as: ``user``
- user@localhost’s password: ``user``

.. image:: _static/putty.png

In Mac OS X or Linux, you start your terminal program first ('Terminal', 'Konsole', 'xterm').
In the terminal window your shell is running (probably 'bash'). Here you use ssh to connect 
to the virtual machine::

  ssh -p 2222 user@localhost

Testing the Internet connection
-------------------------------
When you're logged in, check your internet connection from the virtual machine. Your main
computer has to be connected to the internet, of course. Copy the following command, and 
paste it to the command prompt (click right mouse button in PuTTY window).

.. code-block:: bash

  wget -q -O - http://goo.gl/n8XK2Y | head -1
  # <!DOCTYPE html>

If the ``<!DOCTYPE html>`` does not appear, something is probably wrong with the connection.

Connect to copy files
---------------------
In Windows, WinSCP is used to copy files to Linux machines. You use the same information
as for PuTTY to log in.

.. image:: _static/winscp.png

In Mac OS X or Linux, the most simple command to copy a file into 
a home directory of ``user`` on a running virtual machine is::

  scp -P 2222 myfile user@localhost:~

Connect to RStudio
-----------------
This is the easiest one, just click this link: `Open RStudio <http://localhost:8787>`_.
Login with the same credentials (user, user).

