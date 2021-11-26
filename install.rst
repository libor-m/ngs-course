Installation instructions
=========================

We will be all connecting to a single remote server hosted by
`MetaCentrum <https://www.metacentrum.cz/>`_.

.. note::
    You will connect to the machine differently, depending on what
    operating system you are using.

Windows
^^^^^^^
We need to install two packages:

  - `Git for Windows <https://git-scm.com/download/win>`_
  - `WinSCP <http://winscp.net/eng/download.php>`_

We'll use **Git Bash** from the Git for Windows package to control the remote
computer.

When installing Git for Windows, keep the default settings, just
make sure you check ``Git Bash Here``. Keep in mind that the images are from an
older version of the installer.

.. image:: _static/git-for-win-1.png

.. image:: _static/git-for-win-2.png

To set up your terminal run the ``Git Bash`` from Start menu,
run this and exit the terminal (``exit``)::

  curl -sL https://owncloud.cesnet.cz/index.php/s/1B1NnrI4lqQtY9Q/download > ~/.minttyrc

**WinSCP** will be used to transfer files between your computer and the remote
computer. Look for ``Installation package`` on the
`WinSCP download<http://winscp.net/eng/download.php>`_ page.

.. note::
    Microsoft recently developed an app called `Windows Terminal`, which seems
    superior in functionality to `mintty`, but it is still more difficult to set
    up. That's why we stick with `Git Bash`.

macOS and Linux
^^^^^^^^^^^^^^^
``ssh`` is used to control the virtual computer. It should be already installed in your computer.

Files can be transferred with ``scp``, ``rsync`` or ``lftp`` (recommended)
from the command line. `scp` and `rsync` could be already installed in your system,
if you want to try `lftp`, you'll probably have to install it yourself.

Mac users that prefer grapical clients can use something like `CyberDuck`. Check
the `answers here
<http://apple.stackexchange.com/questions/25661/whats-a-good-graphical-sftp-utility-for-os-x>`_.


Chromebook
^^^^^^^^^^
On Chromebook, you need to install the
`Secure Shell App <https://chrome.google.com/webstore/detail/secure-shell-app/pnhechapfaindjhompbnflcldabbghjo?hl=en>`_.

Another option is to install the `full blown Linux <https://chromeos.dev/en/linux>`_, if your device allows it.

Time to log in!
^^^^^^^^^^^^^^^
Try to log in following the instructions in :ref:`ssh_connect`.
