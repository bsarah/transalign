= Transalign - calculate transitive alignments 

This program calculates _transitive alignments_, it takes as its input
a set of alignmnets from a set of _query_ sequences to an
_intermediate database_, a set of aligmnets from the intermediate
database to a _target database_, and outputs a set of alignments from
the queries to the targets.

Typically, the target database is well described but distantly related
to the query sequences, and the intermediate database is large, but
less well described. Using the intermediate to provide a large
evolutionary _context_ allows `transalign` to detect distant
relationships with a higher sensitivity than direct alignments, and
without needing to construct explicit stochastic models of the
targets.

== Running transalign

First, you need BLAST results in XML format, and an installed
`transalign` executable (see the next section for this).  This will
require a bit of disk space.  We're going to use UniRef 50 as the
intermediate database, and SwissProt (sprot.fa) as the final target.
The process will look something like this:

    formatdb -i uniref50.fa -pT
	formatdb -i sprot.fa -pT
	blastx -i input.fa -d uniref50.fa -e1e-4 -a 8 -m 7 -o inp_vs_u50.xml
	blastp -i uniref50.fa -d sprot.fa -e1e-4 -a 8 -m 7 -o u50_vs_sp.xml

(Options are: `-e` limits the evalue of hits to avoid generating an excess of false
positives, `-a` specifies the number of parallel threads, and `-m`
specifies the output format, in this case XML.)
You should now have the necessary input data, and you can run

    transalign inp_vs_u50.xml u50_vs_sp.xml > inp_vs_sp.txt

== Transalign options

You can display the brief, built-in help by running 
`transalign --help`.  This gives the following output:

transalign v0.1, ©2012 Ketil Malde

    transalign [OPTIONS] [BLASTXML FILES]
    Calculate sequence alignments transitively.

    Common flags:
      -l --long             long output detailing matching positions
      -c --cache            generate alignment cache for initial alignment
      -n --limit=INT        max number of alignments to consider
      -o --outfile=ITEM     output to file
      -e --extract=ITEM     explicit list of alignments to extract
         --cite             output citation information
      -b --blastfilter=NUM  exclude intermediate alignment with per-column score
                            less than this value (not using this option disables
                            the filter)
      -d --debug            show debug output
      -? --help             Display help message
      -V --version          Print version information
         --numeric-version  Print just the version number
      -v --verbose          Loud verbosity
      -q --quiet            Quiet verbosity

Long output produces a large table matching query positions with
target positions, while the default is to output a table similar to
BLAST tabular output.

Sometimes BLAST will generate a large number of alignments (for
instance will very repetitive proteins generate many alternative
pairwise matches) which causes `transalign` to consume a substantial
amount of memory.  You can limit the number of considered alignments
using `-n`.

Using the `-e` option, `transalign` outputs only alignments for the
requested query sequences.  This is most useful when the alignment
caches are already generated.

To speed up operation and avoid doing the same work over, `transalign`
builds alignment caches (for BLAST results `foo.xml`, it will create a
directory `foo.xml.d` containing the cache).  You can also construct
this cache separately, using `blastextract foo.xml`.  There is also a
program, `showcache`, to inspect a cached alignment.  The default is
to build a cache for the second step, but not the first, `-c` builds
caches for both input files.

== Good and bad practices

As the BLAST XML output can sometimes be large, transalign will
parse these files once to generate a cache for them.  This will
generate a large number of files in a single directory, so you need a
file system that handles this.  Some do, but NFS is _not_ one of
them.

== Examining the output

The output is in a table format, somewhat similar to BLAST's (-m
8). The columns are: Query, Target, Score, Alignment length, Average
score, Query start, Query end, Target start, Target end.

== Downloading and installing transalign

The program is written in Haskell, and distributed as source code.
This means you need a working Haskell compiler and environment, and
optionally, the necessary tools to download the source from its
`darcs` repository.

See the generic
[installation instructions](http://biohaskell.org/Installation) for
details on the various ways to acquire and install `transalign`.



