Scrimer
=======
Scrimer is a pipeline for designing primers from transcriptome.
Most of the steps are general to NGS data analysis. 

Go to::

  http://scrimer.rtfd.org

In the documentation you will find commands to perform individual steps
of NGS analysis. 

The data in your virtual machine are filtered to a size where the 
steps won't take too long - so you should be able to run them all 
during the session.

Looking at quality checks
-------------------------

Steps to take
-------------
You should try to understand what the 
commands do (you don't have to understand *how* they do it).

Go to ``/data/slavici_sandbox`` directory in your virtual machine.

Set the number of CPUs that can be used (there's only one in the VM)::

  CPUS=1

The following headings can be found in the **Scrimer manual**. Follow instructions
there, with additional info given here.

Remove cDNA synthesis adaptors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set the adaptor sequences, that were used 
during the library preparation. They will be removed from the sequences.

.. code-block:: bash

    # primers used to synthetize cDNA
    # (sequences were found in .pdf report from the company that did the normalization)
    PRIMER1=AAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTTNN
    PRIMER2=AAGCAGTGGTATCAACGCAGAGTACGCGGG
    PRIMER3=AAGCAGTGGTATCAACGCAGAGT

Run the QC on the raw data and on the trimmed data - so you can see the difference.
QC reports produced by ``fastqc`` are html pages. You need to get them to your machine,
because your Linux does not have any graphical display (only text).

Use WinScp to copy them to your machine, after you generate them.

Map reads to reference assembly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Follow all the instructions on the page.

Detect and choose variants - Call variants with FreeBayes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Skip the line following ``# the rest, ignore order``.
