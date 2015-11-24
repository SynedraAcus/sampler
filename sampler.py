#! /usr/bin/env python3
__author__ = 'morozov'

import argparse
from Bio import SeqIO
from reductor import *
import sys

arg_parser = argparse.ArgumentParser(description='Reduce sequence dataset using Distant Joining approach')
arg_parser.add_argument('-f', type=str, help='Fasta file with sequence set to be reduced')
arg_parser.add_argument('-d', type=str, help='Distance matrix file in PHYLIP lower-triangular format')
arg_parser.add_argument('-p', type=bool, help='Set to true if using protein sequence set',
                        default=True)
arg_parser.add_argument('-n', type=str, help='Number (int) or percentage (float) of sequences to be sampled\nDefault is 0.2',\
                        default='0.2')
arg_parser.add_argument('-i', type=bool, help='Write IDs of reduced set and quit without writing distance matrix or FASTA',
                        default=False)
args = arg_parser.parse_args()
if not(args.f or args.d):
    sys.stderr.write('Either -f or -d option should be used')


#  READ SEQUENCES IF ANY
#  EITHER READ OR CALCULATE MATRIX
if args.d:
    distmat = DistanceMatrix(handle=open(args.d))
    length = len(distmat.ids)  #  Yes, I could've given matrix __len__
if args.f and args.d:
    if not sorted(seqs.sequences.keys()) == sorted(distmat.ids):
        raise ValueError('Different sequence collections in matrix and FASTA')
if args.f and args.p and not args.d:
    with open(args.f) as fasta_fh:
        distmat = make_aminoacid_matrix(fasta_fh)
        length = len(distmat.ids)
if args.f and not args.p and not args.d:
    raise NotImplementedError('DNA distance calculation is not yet supported')


#  How many sequences do we need?
if '.' in args.n:
    final_number = int(float(args.n)*length)
else:
    final_number = int(args.n)

#  Actually subsampling, creating distance matrix if needed
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
if args.i:
    with open('{0}.ids'.format(handle_base), mode='w') as id_handle:
        for x in sample_ids:
            print(x, file=id_handle)
    quit()

#  Assume it's not args.i, so we need to print the rest
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