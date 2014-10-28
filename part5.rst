Visualizing your data
=====================
We'll use RStudio that is installed in the virtual machine. To start RStudio
go to (in your browser)::

  http://localhost:8787

Look around the program - if you ever user R withou RStudio, the difference is big!

Get the data::

  # download it here to your machine
  https://owncloud.cesnet.cz/public.php?service=files&t=aab865a16555adc995b50e33b148318a

  # use WinScp to copy it to virtual machine

  # unpack the data
  unzip data_viz.zip

Then use RStudio to navigate to the folder where you unpacked the data. 
Open ``multires_profiles.R``. 

Use package manager to install ``gtools`` library. (It's easy.)

Run the scripts, and see GC profiles of various genomes.

You can suggest other types of plots on this kind of data - we can try to create them.

Extra UNIX excersise
--------------------
If you want to plot the same profiles for ``/data/slavici``, there is a script 
``base_counts.py`` that can sum it for you. But I wrote it in a hurry for a certain
type of data - one gzipped chromosome per file. You can try to convert ``luscinia_small.fasta``
to this format and run the script.
