#! /usr/bin/env python3
__author__ = 'morozov'

import argparse
from Bio import SeqIO
from reductor import *
import sys

arg_parser = argparse.ArgumentParser(description='Reduce sequence dataset using Distant Joining approach')
arg_parser.add_argument('-f', type=str, help='Fasta file with sequence set to be reduced')
arg_parser.add_argument('-d', type=str, help='Distance matrix file in PHYLIP square format')
arg_parser.add_argument('-n', type=str, help='Number (int) or percentage (float) of sequences to be sampled\nDefault is 0.2',\
                        default='0.2')
args = arg_parser.parse_args()
if not(args.f or args.d):
    sys.stderr.write('Either -f or -d option should be used')

if args.f:
    seqs = SequenceSet(handle=open(args.f))
#  Insert matrix/seq parsing here so we can calculate uniform number

if '.' in args.n:
    final_number = float(args.n)
else:
    final_number = int(args.n)