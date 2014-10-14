Course materials preparation
============================ 
This section contains the steps that we did to produce the materials that course participants
got ready-made. That is the **linux machine image**, **online documentation** and the **slide deck**.

Online documentation
--------------------
Login to https://github.com. Create a new project called `ngs-course-nhrady`, with a default readme file.


Clone the project to local machine and initialize `sphinx` docs. Choose ``SSH`` clone link in GitHub.

.. code-block:: bash

  git clone git@github.com:libor-m/ngs-course-nhrady.git

  cd ngs-course-nhrady
  
  # use default answers to all the questions
  # enter project name and version 1.0
  sphinx-quickstart


Now track all files created by `sphinx-quickstart` in current directory with `git` and publish 
to GitHub.

.. code-block:: bash
  
  git add .
  git commit -m 'empty sphinx project'

  # ignore _build directory in git
  echo _build >> .gitignore
  git add .gitignore
  git commit -m 'ignore _build directory'
  
  # publish the first docs
  # setting up argument less git pull with '-u'
  git push -u origin master
  
To get live view of the documents, login to https://readthedocs.org. Your `GitHub` account can be paired with 
`Read the Docs` account in `Edit Profile/Social Accounts`, then you can simply 'import' new projects 
from your GitHub with one click. Import the new project and wait for it to build. After the build
the docs can be found at http://ngs-course-nhrady.readthedocs.org (or click the ``View`` button).
  
Now write the docs, commit and push. Rinse and repeat. Try to keep the commits small, just one change a time.

.. code-block:: bash
  
  git add _whatever_new_files_
  git commit -m '_your meaningful description of what you did here_'
  git push

References that may come handy:

- `Thomas Cokelaer's cheat sheet <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
