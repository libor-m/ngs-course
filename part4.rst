Graphical tools
===============

IGV - Integrative Genomics Viewer
---------------------------------
Use ``WinScp`` to copy your data from the virtual machine. Run IGV.
Load the data to IGV and look around.

Galaxy
------
The same data that you see in IGV can be visualized online by uploading
to ``Galaxy`` server. I uploaded the data beforehand, so we don't upload 10x 500 MB
at once.

To see the loaded data, go to::

  https://usegalaxy.org/u/liborm/v/ngs-course-2014

To get the data to Galaxy interface, go to::

  https://usegalaxy.org/u/liborm/h/ngs-course-2014

To see Galaxy, go to::

  http://usegalaxy.org

MetaCentrum
===========
The best thing is that now you know almost everything to use MetaCentrum.
It is the same as using PuTTY to acces a local virtual machine.

Key differences
---------------
* you have to register
* you don't use localhost:2222, but something like skirit.ics.muni.cz:22
* you cannot use ``sudo``
* you need to allocate computers using ``qsub`` command
* your data is somewhere else than in ``/data``
* you can have 100 cores instead of 1