## Sampler: dataset reduction utility for phylogenetic analyses

Sampler takes a multi-FASTA (not an alignment!) or an EMBOSS-formatted
distance matrix and samples a requested number of sequences from it,
trying to be as comprehensive as possible. If a FASTA file is supplied
and the distance matrix is not, the distances are calculated by
Scoredist metric.

Returns a reduced distance matrix and, if possible, a reduced FASTA file.



### Installation and running:

    git clone https://github.com/SynedraAcus/sampler
    cd sampler

To reduce a dataset, call `sampler.py`. Either `-f`, `-d` or both should
always be supplied. Common uses:

Sample half of the sequences from aminoacid FASTA
Builds distance matrix internally

    ./sampler.py -f sequences.fasta -p -n 0.5

Sample exactly ten sequences from the same file

    ./sampler.py -f sequences.fasta -p -n 10

Sample half of the sequences, using the supplied matrix
FASTA file is used only for creating output files and can be
omitted

    ./sampler.py -f sequences.fasta -d matrix.txt -p -n 0.5

Get only IDs of the selected sequences

    ./sampler.py -f sequences.fasta -d matrix.txt -p -i -n 0.5

If FASTA is not supplied, FASTA output cannot be produced and -i
is added implicitly

    ./sampler.py -d matrix.txt -n 0.5


#### Dependencies:

* Python 3

* Biopython

* numpy (tested on 1.8.2 but should probably work on older versions)

* EMBOSS 6.6.0.0

Make sure EMBOSS `needle` is available in your $PATH. If it is not
possible (but you do have EMBOSS installed), add EMBOSS/bin directory
path to the NEEDLE_CALL variable in reductor.py

Sampler uses `needle` for pairwise alignment, a crucial step in the
distance matrix construction. As such, without needle it can only work
with supplied distance matrices.

### Options:

`-f` Multi-FASTA file

`-d` PHYLIP-formatted square distance matrix in lower triangular format

`-n` Number (int) or proportion (float) of sequences to sample

`-p` FASTA contains aminoacid sequences. Defaults to False

`-i` Print only IDs of selected sequences. Defaults to False

`-h` Print option list and quit

### Using reductor library

The sampler script is just a wrapper around the `reductor.py` library,
which can be used directly in your Python3 projects. Two classes can be
imported from this library: `DistanceMatrix`, which is a numpy-based
distance matrix implementation supporting the sampling procedure, and a
`DistanceMatrixFactory` to create or read matrices. The description of
the API is available from classes' and methods' docstrings.

### LICENSING
This software is distributed under the terms of CC-BY 4.0 license. The
complete text of the license is available at:

https://creativecommons.org/licenses/by/4.0/legalcode

When used in a scientific publication, please cite this paper:

A.A. Morozov, Y.P. Galachyants 2017 "Distant joining: a sequence
sampling method for the complex evolutionary histories".
