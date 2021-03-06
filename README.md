# fiat
FASTA In A Terminal

## Introduction

*fiat.py* is a text-mode viewer for DNA sequences read from FASTA/Genbank files.
It offers a powerful search feature that allows for mismatches and ambiguous
bases, and the ability to define regions in the sequence, characterized by
a color and a name.

Fiat can handle multiple sequences at the same time, and allows switching between
sequences at any time. Search results and regions are specific to each sequence.

## Startup

The fiat.py command line is:

```bash
$ fiat.py files...
```

where files represents one or more files in FASTA or Genbank format. Files in FASTA
format can contain multiple sequences. Each sequence contained in the input files
will be loaded into fiat. The program will then display the contents of the first loaded sequence, starting at position 1.

## Display

The fiat.py display is divided into four regions:

* The top three rows contain the ruler, showing base positions from 1 to 60, followed
by an empty line.

* The rest of the screen, except for the last two rows, contains the sequence. The format
of each line is: position of the first base in the line, 60 bases, and optionally the name
of all hits and all regions present in the line.

* The last-but-one row (in reverse) is the status row. It displays the name of the currently displayed
sequence, its length, the range of bases currently visible, the number of define regions
and the number of hits (see Regions and Searching below).

* The last row is the /message/ row. It is used to display messages or to accept user input.

The following image shows an example of what the display looks like:

![demo screenshot](demo/screenshot.png)

## Usage

When displaying a sequence, use the `?` key to display all available commands. Use the `q` command
to quit the program.

## Basic commands

You can use the arrow keys, Home/End, PgUp/PgDn, Enter and Spacebar to navigate through the sequence.
Left arrow scrolls back by 10 lines, right arrow scrolls forward by 10 lines.

The `g` key (for 'go to') asks for a position in the sequence, and scrolls the display so that the
line containing that position is at the top of the screen.

The `s` key prints some statistics on the sequence in the message row. Currently these consists of
the counts and percentages for each base, and the GC%.

The `o` key is used to set options. Currently these are:
* The maximum number of mismatches allowed in a search (with the `m` key);
* The search mode (with the `s` key) - this can be any combination of the four letters f, r, c, d,
  meaning forward, reverse, complement, and reverse-complement respectively. For example, if the mode
  is `fr` the program will search for sequences in both forward and reverse orientations. See the Searching
  section for more details. 

## Regions

Regions are arbitrary subsequences characterized by a start and end position,
a name, and a color. Regions can be defined manually (using the `a` key) or
as a result of a search. Use `<` and `>` to focus the previous / next region
respectively. When a region is focused, the message line starts with `* Region`
and shows the number of the current region and its coordinates.

If you use the `a` key to enter a region manually, the program will prompt for the
region coordinates, its name, and its color (red by default). Coordinates can be entered
in one of two formats: `start-end`, or `start+length`. For example, `1000-1500` and `1000+500`
indicate the same region.

Key | Command
----|--------
<, >         | jump to previous / next region
r            | rename the currently focused region
c            | change color of the currently focused region
d            | delete the currently focused region

The following five colors can be used for regions, identified by their initial:
(r)ed, (g)reen, (b)lue, (m)agenta, (c)yan.

## Searching

The search command finds all occurrences of a specified pattern (hits). Hits are
displayed in black on yellow in the sequence, and their names are displayed in
yellow. The name of a search hit is automatically generated and cannot be changed. 

The pattern can be specified using the full IUPAC nucleotide alphabet, shown here:

```
A      K = G/T       B = C/G/T
C      Y = C/T       D = A/G/T
G      S = C/G       H = A/C/T
T      W = A/T       V = A/C/G
       R = A/G
       M = A/C       N = A/C/G/T
```

For example, the pattern CAST will match CACT or CAGT.

Matching is also affected by the maximum number of mismatches allowed (that can be
set using the `o` command followed by `m`) and by the search mode (set with `o`
followed by `s`), which can be any combination of the following:

  f = forward, r = reverse, c = complement, d = reverse complement

For example, if the search mode is fd, GATTA will match both GATTA and TAATC.

Use `,` and `.` to focus the previous / next hit respectively. When a hit is focused, 
the message line starts with `* Hit` and shows the number of the current hit, its 
coordinates, and its orientation.

Key | Function
----|---------
.    | jump to the previous hit
,    | jump to the next hit
a    | save this hit as a region (prompts for name and color)
A    | save all hits as regions (prompts for color)

To clear the list of hits, enter `-` as the search string.
