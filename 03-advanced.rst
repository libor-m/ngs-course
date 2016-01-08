Unix - Advanced II
==================

Scripting session: scripting in one line using ``awk``,
writing functions and scripts in shell, and running procedures in parallel.


Scripting in one line (awk)
---------------------------

Use nightingales FASTQ files:

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

	LENGTH=80

	< data/fastq/HRTMUOC01.RL12.00.fastq \
	awk -v l=$LENGTH '{
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

.. code-block:: bash

	uniq(){ uniq -c | sed -r 's/^ *([0-9]+) /\1\t/'; }

Shell Scripts
-------------

.. code-block:: bash

	nano script.sh

.. code-block::bash

	#!/bin/sh

	FILE=$1
	LENGTH=$2
	OUT=$1-filtered

	< data/fastq/HRTMUOC01.RL12.00.fastq \
	awk -v l=$LENGTH '{
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

Parallel
--------
