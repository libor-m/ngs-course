Best practice
=============

This is a collection of tips, that may help to overcome the initial barrier of working with a 'foreign' system.
There is a lot of ways to achieve the soulution, those presented here are not the only correct ones, but some
that proved beneficial to the authors.

.. note:: Sorry, this is just a document skeleton now.

Easiest ways to get UNIX
------------------------
- VirtualBox with PuTTy for Windows
- Mac OS X and Linux are UNIX based

Essentials
----------
Always use screen session for any serious work. You can track progress of your 
operations with ``pv`` if you're parsing single file, or at least control the 
resource usage in ``htop``. If no CPU is spinning an no IO is happening, maybe 
something is wrong..

Data organization
-----------------
- waterfall sysem
- ``workflow.sh`` in the top directory of waterfall, all data in subdirectories

Building command lines
----------------------
- write down your precious pipelines carefully - you'll reuse them often!
- input first syntax (``<file command | comm2 | comm3 > out) - improves legibility
- wrap your long pipelines on ``|``

.. code-block:: bash

  <infile sort -k3,3 |
    uniq -c -s64 |
    sort -k1rn,1 \
  >out
  
- check pieces of pipeline often with head
- enjoy progress indicator with pv
