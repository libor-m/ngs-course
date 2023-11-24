Session 7: Git and GitHub
=========================

Once you start writing your scripts, you'll soon find yourself handling files
named like ``script.sh``, ``script_bak.sh``, ``script_previous_previous.sh``,
``script_last_working.sh`` .. etc. There is one cure for them all: ``git``.
It was originally created by Linus Torvalds, the author of the Linux kernel,
for managing the source code of Linux, it slowly gained popularity in other
communities.

What ``git`` does is managing versions of a directory tree. The managed subtree
is called a **repository**. Each saved version is called a **commit**. You
usually create a commit when the code you're working on behaves as expected 🙂.
By allowing you to go back to any committed version, ``git`` effectively removes the
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

  git config --global user.email "ferda@mraveniste.cz"
  git config --global user.name "Ferda Mravenec"

To be able to share your code, you need an account on a 'social coding' site
like `GitHub <https://github.com>`_. If you don't have a GitHub account, please
create one.

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
