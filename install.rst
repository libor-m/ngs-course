Installation instructions
=========================

We will be using a virtual computer with Debian Linux, all the software tools
and all the data necessary for the exercises.

.. note:: 
  You need to install the image even if your main system is Linux / Mac OS X!

Installation steps (it should take about 10 minutes):

- Install VirtualBox (https://www.virtualbox.org/wiki/Downloads). It works on Linux and Mac too.
- Download the virtual machine image from this link: http://goo.gl/FwFk8Z You'll get a single
  file with ``.ova`` extension on your hard drive.
- You can either double click the ``.ova`` file, or run VirtualBox, and choose ``File > Import Appliance``.
  Follow the instructions after the import is started.

After successful installation you should see something like this (only the machine list will contain just one machine).
Check whether you can start the virtual machine: click ``Start`` in the main VirtualBox window:

.. image:: _static/vbox-main.png

After a while you should see something like this:

.. image:: _static/vbox.png

You don't need to type anything into that window, just checking that it looks like the screen shot is enough.

In case of any problems contact the tutors, we'll try to resolve all problems before the course. As the last resort
we can try to resolve remaining problems in the beginning session.

How to access the machine
-------------------------
Because it is much more comfortable to use a native terminal application than the small VM screen,
you will connect to the machine depending on what system you are using.

Windows
-------
Install PuTTY and WinSCP. PuTTY will be used to control the virtual computer. WinSCP will be used to transfer
files between your computer and the virtual computer.

- PuTTY (http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html - look for putty.exe) 
- WinSCP (http://winscp.net/eng/download.php - look for Installation package). 

Mac OS X and Linux
------------------
Ssh is used to control the virtual computer. It should be installed in your computer.

Files can be transferred with ``scp``, ``rsync`` or ``lftp`` (recommended) 
from the command line. `Scp` and `rsync` could be already installed in your system, 
if you want to use `lftp`, you'll probably have to install it yourself.

Mac users that prefer grapical clients can use something like `CyberDuck`. See
http://apple.stackexchange.com/questions/25661/whats-a-good-graphical-sftp-utility-for-os-x .

Machine configuration details
-----------------------------

- Administrative user: `root`, password: `debian`
- Normal user: `user`, password: `user`
- ssh on port 2222
- RStudio on port 8787

Time to log in!
---------------
Try to log in following the instructions in :ref:`ssh_connect`. 
