#! /usr/bin/env python3
__author__ = 'morozov'

import argparse
from Bio import SeqIO
from reductor import *
import sys

arg_parser = argparse.ArgumentParser(description='Reduce sequence dataset using Distant Joining approach')
arg_parser.add_argument('-f', type=str, help='Fasta file with the sequence set to be reduced')
arg_parser.add_argument('-d', type=str, help='Distance matrix file in PHYLIP lower-triangular format')
arg_parser.add_argument('-p', help='Using protein sequence set. Defaults to False',
                        action='store_true')
arg_parser.add_argument('-n', type=str, help='Number (int) or percentage (float) of sequences to be sampled\nDefaults to 0.2',
                        default='0.2')
arg_parser.add_argument('-i', help='Write IDs of reduced set and quit without writing distance matrix or FASTA\nDefaults to False',
                        action='store_true')
args = arg_parser.parse_args()

#  READ SEQUENCES, IF ANY
#  EITHER READ OR CALCULATE MATRIX
m = MatrixFactory()
if args.d:
    # If matrix was supplied
    distmat = m.read(open(args.d))
    length = len(distmat.indices.keys())  #  Yes, I could've given matrix __len__
    if args.f:
        # Check that the same sequences are in matrix and FASTA file
        fasta_ids = [x.id for x in SeqIO.parse(open(args.f), format='fasta')]
        if not sorted(fasta_ids) == sorted(distmat.indices.keys()):
            raise ValueError('Different sequence collections in matrix and FASTA')
elif args.f:
    # Build matrix
    if args.p:
        distmat = m.create_aminoacid_matrix(args.f)
    else:
        distmat = m.create_nucleotide_matrix(args.f)
    length = len(distmat.indices)
else:
    sys.stderr.write('At least one of -f or -d should be used.\n')
    quit()


#  How many sequences do we need?
if '.' in args.n:
    final_number = int(float(args.n)*length)
    if final_number > length:
        raise ValueError('Cannot sample more sequences than there are!')
else:
    final_number = int(args.n)

#  Actually subsampling
sample_ids = distmat.dj(final_count=final_number)
if not args.i:
    sample_matrix = distmat.submatrix(sample_ids)

#  OUTPUT

#  Generating handle base to which extensions for distance matrix and shit would be added
if args.d:
    handle_base = args.d
else:
    handle_base = args.f

#  Simplified ID output
if args.i or not args.f:
    with open('{0}.ids'.format(handle_base), mode='w') as id_handle:
        for x in sample_ids:
            print(x, file=id_handle)
    quit()

#  It's not args.i, so we need to print the rest
#  Defining matrix handle
if not args.d:
    matrix_handle = open('{0}.dist.reduced'.format(args.f), mode='w')
else:
    matrix_handle = open('{0}.reduced'.format(args.d), mode='w')

# Writing matrix
sample_matrix.write(matrix_handle)

# Writing unreduced matrix, if it was not in the input data
if not args.d:
    full_matrix_handle = open('{0}.dist'.format(args.f), mode='w')
    distmat.write(full_matrix_handle)

# Writing reduced fasta, if any
if args.f:
    with open(args.f) as old_fasta_handle:
        sample_seqs = []
        for seq in SeqIO.parse(old_fasta_handle, format='fasta'):
            if seq.id in sample_ids:
                sample_seqs.append(seq)
    with open('{0}.reduced'.format(args.f), mode='w') as new_fasta_handle:
        SeqIO.write(sample_seqs, handle=new_fasta_handle, format='fasta')
