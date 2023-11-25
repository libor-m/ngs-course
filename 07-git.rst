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

Configure git
-------------
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

Create a local repository
-------------------------
And now after the tedious setup let's enjoy the benefits. You have to do the
setup only once per machine, and then you use it for all your projects.

Let's store the current version of your scripts in ``unix-advanced`` project.
Git repo is not a good place for data, we'll ignore the `data` directory.

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
  git commit -m "solution to ngs-course exercises"

This is all you need to track your versions *locally*. If you want to share your
creations with others, transfer it to another machine or just make a backup for
yourself, you need to **push** it somewhere.

Upload repository to GitHub
---------------------------
In your GitHub account, create a repository and use the suggested commands to add it as
a **remote** to your local repo and `push`. Select *Private* if you don't want your code
to be publicly accessible.

.. code-block:: bash

  # use the commands suggested by GitHub to add a remote
  git remote add origin git@github.com:{user}/{repo}.git
  git push -u origin master

When using git, you can gradually learn about more concepts and commands, as
you find the need for them. To give you a head start:

git commands: Basics
--------------------

.. glossary::

  working copy
    the directory you're currently working in

  ``git status``
    check the working copy for changes versus last commit

  ``git add {file}``
    add file to the next commit

  ``git commit``
    create a commit from all added files (opens editor for commit message)

  ``git commit -am '{what i did}'``
    shortcut for adding all changed files (previously in repo) and committing them

  ``git push``
    upload the commits to a remote repository

  ``git pull``
    download new commits from a remote repository, merge them into your working copy

  ``git stash``
    can be used to "hide" local changes during ``pull``

  ``git stash pop``
    brings the "hidden" changes back

  ``git log``
    show the previous commits

  ``git checkout -- {filename}``
    overwrite the file with the version from the last commit

git commands: Branches
----------------------

.. glossary::

  ``git checkout -b {new-name}``
    create a new branch from the current one and switch to it

  ``git push -u origin {branch-name}``
    upload commits on current branch to a remote repository

  ``git checkout {branch-name}``
    switch to another branch

git commands: Merge
-------------------

.. glossary::

  ``git pull``
    if the remote branch has new commits and your local brach has some other commits
    ``pull`` will do a merge for you

  ``git checkout --theirs``
    in conflict, choose the version from the remote branch

  ``git checkout --ours``
    in conflict, choose the version from the local branch

  ``git merge --continue`` or ``git commit``
    continue the merge after resolving conflicts

  ``git merge --abort``
    abort the merge
