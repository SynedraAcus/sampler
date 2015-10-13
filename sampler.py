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

#  Sequence/distmat input
if args.d:
    raise NotImplementedError('Distance input is not implemented yet')
if args.f:
    seqs = SequenceSet(handle=open(args.f))
    length = len(seqs)

#  Reading number of sequences needed

if '.' in args.n:
    final_number = int(float(args.n)*length)
else:
    final_number = int(args.n)

#  Actual subsampling

final = seqs.matrix.dj(final_number)

#  Output

if args.f:
    for seqid in final.ids:
        SeqIO.write(seqs[seqid], sys.stdout, 'fasta')
else:
    print('\n'.join(final.ids))