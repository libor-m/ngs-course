Unix - Advanced II
==================

Scripting session: 

- scripting in one line using ``awk``
- writing functions and scripts in shell
- speeding up your processing by running in parallel


Scripting in one line (awk)
---------------------------

*Use nightingale FASTQ files:*

1. Extract IDs of a FASTQ file and count the number of reads

.. code-block:: bash

	< data/fastq/HRTMUOC01.RL12.00.fastq awk '{ if( (NR + 3) % 4 == 0 ){ print $0 } }' | wc -l

2. Make a file with read ID and read lengths in one line

.. code-block:: bash

	< data/fastq/HRTMUOC01.RL12.00.fastq \
	awk 'BEGIN{
		OFS="\t"
	}{
		if(( NR + 3 ) % 4 == 0 ){
			id = $0
		}else{
			if( (NR + 3) % 4 == 1 ){
				print id,length($0)
			}
		}
	}' | less

3. Get average read length

.. code-block:: bash

	< data/fastq/HRTMUOC01.RL12.00.fastq \
	awk 'BEGIN{
		OFS="\t"; l=0; n=0
	}{
		if( ( NR + 3 ) % 4 == 1 ){
			l = l + length($0);
			n = n + 1;
		}
	}END{
		print "Average read length:", l/n
	}'

4. Filter out short sequences (set the minimum size allowed)

.. code-block:: bash

	< data/fastq/HRTMUOC01.RL12.00.fastq \
	awk -v l=80 '{
	  if( (NR + 3) % 4 == 0 ){
	    id=$0;
	  }else if( (NR + 3) % 4 == 1 ){
	    seq=$0;
	  }else if( (NR + 3) % 4 == 2 ){
	    q=$0;
	  }else{
	    if( length(seq) >= l ){
	      print id"\n"seq"\n"q"\n+";
	    }
	  }
	}' | less

Functions in Shell
------------------
Create a command ``uniqt`` that will behave as ``uniq``, but there
will be no padding (spaces) in front of the numbers, and numbers will 
be separated by <tab>, so eg. ``cut`` will work.

Do not use the same name as the original command, otherwise you'll create
an endless loop.

.. code-block:: bash

	uniqt(){ uniq -c | sed -r 's/^ *([0-9]+) /\1\t/'  ;}

Shell Scripts
-------------

.. code-block:: bash

	nano script.sh

Make a script ``filter_fastq.sh`` which reads a FASTQ file, filters out short
sequences and saves to a file named ``$INPUT-filtered``:

.. code-block:: bash

	#!/bin/sh

	FILE=$1
	LENGTH=$2
	OUT=$1-filtered

	< $FILE awk -v l=$LENGTH '{
		if( (NR + 3) % 4 == 0 ){
			id=$0;
		}else if( (NR + 3) % 4 == 1 ){
			seq=$0;
		}else if( (NR + 3) % 4 == 2 ){
			q=$0;
		}else{
			if( length(seq) >= l ){
				print id"\n"seq"\n"q"\n+";
			}
		}
	}' > $OUT

	echo File `basename $FILE` done

To run the script:

.. code-block:: bash

    chmod +x filter_fastq.sh
    # check with ls, filter_fastq.sh should be green now
    # and using ll you should see the 'x' (eXecutable) permission
	./filter_fastq.sh data/fastq/HRTMUOC01.RL12.00.fastq 80

	# or, without a need for the shebang line (#!) in the file
	# and without +x permission
	bash filter_fastq.sh data/fastq/HRTMUOC01.RL12.00.fastq 80

Parallel
--------

Runs one instance of the command per each CPU in your machine. Regretably your
**virtual** machine has only one CPU, so this won't help much. But modern
machines do have  four and more CPUs, and then it really helps.

Do control the number of jobs (``-j``) only when sharing the machine with
someone, or when you're sure that yout task is IO bound. Otherwise
``parallel`` does a good job choosing the number of tasks to run.

.. code-block:: bash

  parallel 'bash script.sh {} > {}.out' ::: {1..10}

Run the ``filter_fastq.sh`` in parallel:

.. code-block:: bash

  parallel 'bash filter_fastq.sh {} 80' ::: data/fastq/*.fastq

There is a lot of magic to be done with ``{.}, {/}, {#}`` placeholders,
check ``man parallel``. If your data is a single file, but the processing 
of one line is not dependent on the other lines, ``split`` will help.