# fiat
FASTA In A Terminal

## Introduction

*fiat.py* is a text-mode viewer for DNA sequences read from FASTA files.
It offers a powerful search feature that allows for mismatches and ambiguous
bases, and the ability to define regions in the sequence, characterized by
a color and a name.

Fiat can handle multiple sequences at the same time, and allows switching between
sequences at any time. Search results and regions are specific to each sequence.

## Startup

The fiat.py command line is:

```bash
$ fiat.py fastafiles...
```

where fastafiles represents one or more files in FASTA format. Files can contain
multiple sequences. Each sequence contained in the input files will be loaded
into fiat. The program will then display the contents of the first loaded sequence,
starting at position 1.

## Display

The fiat.py display is divided into four regions:

* The top three rows contain the ruler, showing base positions from 1 to 60, followed
by an empty line.

* The rest of the screen, except for the last two rows, contains the sequence. The format
of each line is: position of the first base in the line, 60 bases, and optionally the name
of all hits and all regions present in the line.

* The last-but-one row is the status row. It displays the name of the currently displayed
sequence, its length, the range of bases currently visible, the number of define regions
and the number of hits (see Regions and Searching below).

* The last row is the /message/ row. It is used to display messages or to accept user input.

See ![demo screenshot](demo/screenshot.png).
