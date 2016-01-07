Additional exercises
====================
These are tasks that do not fit any particular session, but we still
consider them interesting enough to share them with you.

Counting heads
^^^^^^^^^^^^^^
This is a nice example where bash can be used to solve a combinatorial 
problem by enumerating all the possibilities. And it is a long pipeline,
so you have a lot of code to examine;)

.. pull-quote:: 

   Eight people are sitting around a circular table. Each has a coin.
   They all flip their coins. What is the probability that no two adjacent
   people will get heads?

The basic idea is that there is not much possibilities (only 2 to the power of 8,
that is 256). We can just enumerate all the combinations and check if there is 
two adjacent heads.

This is the final solution, take your time to take it apart to see what each 
piece does.

.. code-block:: bash

    (echo "obase=2;"; printf "%d\n" {0..255}) | bc |  # generate numbers 0-255 in binary
      sed 's/^/0000000/' | egrep -o '.{8}$' |         # pad the output to 8 characters
      sed 'h;G;s/\n//' |                              # double the pattern on each line to simulate ring
      grep -v 11 |                                    # ignore all patterns where two heads (1) are next to each other
      wc -l                                           # count the rest


To get a more understandable code, we can split it to functional parts. Then 
we can just play and try different implementations of the parts:

.. code-block:: bash

   generate () { (echo "obase=2;"; printf "%d\n" {0..255}) | bc ;}
   pad () { sed 's/^/0000000/' | egrep -o '.{8}$' ;}
   ring () { sed p | paste - - | tr -d "\t" ;}
   
   generate | pad | ring | grep -v 11 | wc -l

These are alternative solutions - you can paste them one by one, 
and check if the pipe is still working.

.. code-block:: bash

    generate () { (echo "obase=2;"; echo {0..255} | tr " " "\n") | bc ;}
    generate () { (echo "obase=2;"; echo {0..255} | xargs -n1) | bc ;}
    generate () { (echo "obase=2;"; printf "%d\n" {0..255}) | bc ;}

    pad () { awk '{print substr(sprintf("0000000%s", $0), length);}' ;}
    pad () { sed 's/^/0000000/' | rev | cut -c-8 | rev ;}
    pad () { sed 's/^/0000000/' | egrep -o '.{8}$' ;}

    ring () { sed p | paste - - | tr -d "\t" ;}
    ring () { sed 'h;G;s/\n//' ;}

    generate | pad | ring | grep -v 11 | wc -l

The question was asking for the probability, thats one more division::

  echo "scale=3;$( generate | pad | ring | grep -v 11 | wc -l ) / 256" | bc

Solutions by participants
-------------------------
One way to get a shorter (but much slower) solution is to ignore the binary
conversion altogether, just use a huge list of decimal numbers and filter out
anything that does not look like binary. Few variants follow:

.. code-block:: bash

  seq -w 0 11111111 | grep ^[01]*$ | awk '!/11/ && !/^1.*1$/' | wc -l  
  seq -w 0 11111111 | grep ^[01]*$ | grep -v -e 11 -e ^1.*1$ | wc -l
  seq -w 0 11111111 | awk '/^[01]*$/ && !/11/ && !/^1.*1$/' | wc -l
  seq -w 0 11111111 | awk '!/[2-9]/ && !/11/ && !/^1.*1$/' | wc -l
  seq -w 0 11111111 | grep -v -e [2-9] -e 11 -e ^1.*1$ | wc -l

I believe there are still a few ways to make it shorter;)

Dimensionality reduction
^^^^^^^^^^^^^^^^^^^^^^^^
Methods like PCA and MDS (sometimes called PCoA to increase the confusion)
are usually regarded as black box by many. Here we try to present a simple example
that should help with getting a better idea on what are these magic boxes 
actually doing.

Load and visualize your data set
--------------------------------
Let's link the project directory to your own data directory::

  ln -s /data/banana/ ~/data

Now you can go to R and load the data:

.. code-block:: r

  setwd('~/data/banana')
  d <- read.csv("webapp/data/rotated.csv")

Plot the data to look what we've got:

.. code-block:: r

  library(ggplot2)
  ggplot(d, aes(x, y)) + geom_point() + coord_equal()


Correct the distortion
----------------------
Maybe you can already recognize what's in your data. But it appears to be a
bit .. rotated. Here is a code for 3d rotation of points, copy, paste and run it
in your R session:

.. code-block:: r

    # create a 3d rotation matrix
    # https://www.math.duke.edu/education/ccp/materials/linalg/rotation/rotm3.html
    rotX <- function(t) matrix(c(cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1), nrow=3)
    rotY <- function(t) matrix(c(1, 0, 0, 0, cos(t), sin(t), 0, -sin(t), cos(t)), nrow=3)
    rotZ <- function(t) matrix(c(cos(t), 0, -sin(t), 0, 1, 0, sin(t), 0, cos(t)), nrow=3)
    rot3d <- function(tx, ty, tz) rotX(tx) %*% rotY(ty) %*% rotZ(tz)

    # rotate a data frame with points in rows
    rot3d_df <- function(df, tx, ty, tz) {
      rmx <- rot3d(tx, ty, tz)
      res <- data.frame(t(apply(df, 1, function(x) rmx %*% as.numeric(x))))
      colnames(res) <- colnames(df)
      res
    }

Now try to rotate the object a bit, so we can see it better. Try to find good values 
for the rotation yourself (numbers are in radians, 0..2*PI makes sense):

.. code-block:: r

    dr <- rot3d_df(d, .9, .1, 2)
    ggplot(dr, aes(x, y)) + geom_point() + coord_equal()

Enter PCA. It actually finds the best rotation for you. Even in a way that the 
first axis has the most variability (longest side of the object), the second axis
has the maximum of the remaining variability etc.

.. code-block:: r

  pc <- prcomp(as.matrix(dr))
  ggplot(data.frame(pc$x), aes(PC1, PC2)) + geom_point() + coord_equal()
  ggplot(data.frame(pc$x), aes(PC1, PC3)) + geom_point() + coord_equal()
  ggplot(data.frame(pc$x), aes(PC2, PC3)) + geom_point() + coord_equal()

MDS
---
Metric MDS (multidimensional scaling) with `euclidean` distance equals to PCA. We will 
use the non-metric variant here, which tries to keep only the order of pairwise 
distances, not the distances themselves. You prefer MDS when you want to use a different
distance than `euclidean` - we're using `manhattan` (`taxicab`) distance here:

.. code-block:: r

    library(MASS)
    dmx <- dist(dr, "manhattan")
    mds <- isoMDS(dmx)
    ggplot(data.frame(mds$points), aes(X1, X2)) + geom_point() + coord_equal()

Shiny
-----
And now there is something you definitely wanted, while you were trying to find 
the good values for rotation of your object::

  setwd('/data/banana/webapp/')

Now ``File > Open``, and open ``server.R``. There should be a green ``Run App`` 
button at the top right of the editor window. Click that button!