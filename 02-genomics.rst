Unix - Advanced I
=================

This session focuses on plain text file data extraction/modification
using build-in Unix tools.



The commands for this session are in the talk's slide deck.

How many bases were sequenced? => GREP EXERCISE FOR SATURDAY
------------------------------
``wc`` can count characters (think bases) as well. But to get a reasonable number,
we have to get rid of the other lines that are not bases.

One way to do it is to pick only lines comprising of letters A, C, G, T and N.
There is a ubiquitous mini-language called `regular expressions` that can be used
to define text patterns. `A line comprising only of few possible letters` is
a text pattern. ``grep`` is the basic tool for using regular expressions::

  cat *.fastq | grep '^[ACGTN]*$' | less -S

Check if the output looks as expected. This is a very common way to work - build a part of
the pipeline, check the output with ``less`` or ``head`` and fix it or add more commands.

Now a short explanation of the ``^[ACGTN]*$`` pattern (``grep`` works one line a time):

- ``^`` marks beginning of the line - otherwise ``grep`` would search anywhere in the line
- the square brackets (``[]``) are a `character class`, meaning one character of the list, ``[Gg]rep``
  matches ``Grep`` and ``grep``
- the ``*`` is a count suffix for the square brackets, saying there should be zero or more of such characters
- ``$`` marks end of the line - that means the whole line has to match the pattern

To count the bases read, we extend our pipeline::

  cat *.fastq | grep '^[ACGTN]*$' | wc -c

The thing is that this count is not correct. ``wc -c`` counts every character,
and the end of each line is marked by a special character written as ``\n`` (n
for newline). To get rid of this character, we can use another tool, ``tr``
(transliterate). ``tr`` can substitute one letter with another  (imagine you
need to lowercase all your data, or mask lowercase bases in your Fasta file).
Additionally ``tr -d`` (delete) can remove characters::

  cat *.fastq | grep '^[ACGTN]*$' | tr -d "\n" | wc -c

.. note::  If you like regular expressions, you can hone your skills at https://regex.alf.nu/.
