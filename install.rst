Installation instructions
=========================

We will be all connecting to a single remote server hosted by
`MetaCentrum <https://www.metacentrum.cz/>`_.

.. note::
    You will connect to the machine differently, depending on what
    operating system you are using.

Windows
^^^^^^^
Install `Git for Windows <https://git-scm.com/download/win>`_. We'll use it to
control the remote computer.

When installing Git Bash, go with the default settings, just
make sure you check ``Git Bash Here``. Keep in mind that the images are from an
older version of the installer.

.. image:: _static/git-for-win-1.png

.. image:: _static/git-for-win-2.png

To set up your terminal run the ``Git Bash`` from Start menu,
run this and exit the terminal (``exit``)::

  curl -sL https://owncloud.cesnet.cz/index.php/s/1B1NnrI4lqQtY9Q/download > ~/.minttyrc

Install `WinSCP <http://winscp.net/eng/download.php>`_ (look for
``Installation package``). WinSCP will be used to transfer files between your
computer and the remote computer.

macOS and Linux
^^^^^^^^^^^^^^^
``ssh`` is used to control the virtual computer. It should be already installed in your computer.

Files can be transferred with ``scp``, ``rsync`` or ``lftp`` (recommended)
from the command line. `scp` and `rsync` could be already installed in your system,
if you want to try `lftp`, you'll probably have to install it yourself.

Mac users that prefer grapical clients can use something like `CyberDuck`. Check
the `answers here
<http://apple.stackexchange.com/questions/25661/whats-a-good-graphical-sftp-utility-for-os-x>`_.

Time to log in!
---------------
Try to log in following the instructions in :ref:`ssh_connect`.
