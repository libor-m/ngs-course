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
