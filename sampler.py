#! /usr/bin/env python3
__author__ = 'morozov'

import argparse
from Bio import SeqIO
from reductor import *
import sys

arg_parser = argparse.ArgumentParser(description='Reduce sequence dataset using Distant Joining approach')
arg_parser.add_argument('-f', type=str, help='Fasta file with sequence set to be reduced')
arg_parser.add_argument('-d', type=str, help='Distance matrix file in PHYLIP lower-triangular format')
arg_parser.add_argument('-n', type=str, help='Number (int) or percentage (float) of sequences to be sampled\nDefault is 0.2',\
                        default='0.2')
args = arg_parser.parse_args()
if not(args.f or args.d):
    sys.stderr.write('Either -f or -d option should be used')


#  READ SEQUENCES IF ANY
#  EITHER READ OR CALCULATE MATRIX

if args.d:
    distmat = DistanceMatrix()
    distmat.from_handle(open(args.d))
    length = len(distmat.ids)  #  Yes, I could've given matrix __len__
if args.f and args.d:
    seqs = SequenceSet(handle=open(args.f), find_distances=False)
    seqs.matrix = distmat
    length = len(seqs)
if args.f and not args.d:
    seqs = SequenceSet(handle=open(args.f), find_distances=True)
    length = len(seqs)

#  Write matrix to STDERR. Add proper command later

# for seqid in seqs.matrix.ids:
#     sys.stderr.write("{0}\t{1}\n".format(seqid,
#                                        '\t'.join((str(seqs.matrix[seqid, x]) for x in seqs.matrix.ids))))

#  HOW MANY SEQUENCES DO WE NEED?
if '.' in args.n:
    final_number = int(float(args.n)*length)
else:
    final_number = int(args.n)

#  Actual subsampling, creating distance matrix

if args.f:
    final = seqs.matrix.dj(final_number)
else:
    final = distmat.dj(final_number)

#  Output

#  Defining matrix handle
if not args.d:
    matrix_handle = open('{0}.dist.reduced'.format(args.f), mode='w')
else:
    matrix_handle = open('{0}.reduced'.format(args.d), mode='w')

# Writing matrix
final.to_handle(matrix_handle)

# Writing fasta, if any
if args.f:
    fasta_handle = open('{0}.reduced'.format(args.f), mode='w')
    for seqid in final.ids:
        SeqIO.write(seqs[seqid], fasta_handle, 'fasta')