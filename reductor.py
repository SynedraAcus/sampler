__author__ = 'morozov'

from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import subprocess
from math import log
import sys
########
NWALIGN = '/home/morozov/tools/NWalign/NWalign' # Location of NWalign software. REMOVE HARDCODE!
########

#  Converting BLOSUM 62 to square form
matrix = MatrixInfo.blosum62
keys = list(matrix.keys())
for a in keys:
    matrix[(a[1], a[0])] = matrix[a]

def scoredist(seq1, seq2, matrix=matrix, correction=1.337):
    """
    Return a Scoredist distance [Sonnhammer, Hollich 2005] for an alignment of two protein sequences.
    Note: only matched positions affect score: gaps are discarded.
    :param seq1: SeqRecord
    :param seq2: SeqRecord
    :param matrix: dict
    :param correction: float
    :return: float
    """
    EXPECT = -0.5209
    sub = subprocess.check_output([NWALIGN, str(seq1.seq), str(seq2.seq), '3'],
                                  stderr=subprocess.PIPE,
                                  universal_newlines=True)
    lines = sub.split('\n')
    j = 0
    seqs = ['', '']
    for a in lines[7:]:
        if ':' in a:
            j = 1
        elif '1' in a:
            break # If numeric lines started then all sequence was printed already
        else:
            seqs[j]+=a.rstrip()
    #print(seqs)
    del lines, sub #Close this stuff
    score = 0
    score_a = 0
    score_b = 0
    length = 0
    for col in range(len(seqs[0])):
        if seqs[1][col]=='-' or seqs[0][col]=='-':
            continue  # Discard gapped columns
        score_a += matrix[(seqs[0][col], seqs[0][col])]
        score_b += matrix[(seqs[1][col], seqs[1][col])]
        score += matrix[(seqs[1][col], seqs[0][col])]
        length += 1.0
    sigma_n = score - length * EXPECT
    sigma_un = (score_a + score_b)/2 - length * EXPECT
    dist = -1*log(sigma_n/sigma_un)*100*correction
    return dist

class SequenceSet(object):
    '''
    Container class for sequences.
    '''
    def __init__(self, handle=None, sequences=None):
        self.sequences = {}
        self.matrix = DistanceMatrix()
        if handle is not None:
            self._load_fasta(handle)

    def _load_fasta(self, handle):
        for record in SeqIO.parse(handle, format='fasta'):
            self.sequences.update({record.name: record})
            dist_row = [scoredist(record, self[x]) for x in self.matrix.ids]
            self.matrix.add_row(record.id, dist_row)

    def __getitem__(self, item):
        '''
        Return a sequence with the given ID
        '''
        return self.sequences[item]

    def append(self, item):
        '''
        Append a sequence to Sequence set and recalculate distance
        '''
        pass

    def remove(self, item):
        pass

class DistanceMatrix(object):
    '''
    Distance matrix class. Contains a list of IDs and a distance matrix for the same IDs
    '''
    def __init__(self, handle=None):
        self.matrix = {}
        self.ids = []
        if not handle is None:
            self._from_handle(handle)

    def _from_handle(self, handle):
        '''
        Read distance matrix from filehandle
        :param handle: filehandle to read from
        :return:
        '''
        line_list = []
        for line in handle:
            line_list.append(line.split('\t'))
        self.ids = [x[0] for x in line_list]
        for row in line_list:
            del(row[0])
            #  IDs went to the attribute, no use for them here
        for seq1_number in range(len(line_list)):
            for seq2_number in range(len(line_list)):
            #  Order of the lines & items in line is the same as that of IDs
            #  Creating both of items to save __getitem__/__setitem__ hassle. Maybe add that later
                self.matrix[(self.ids[seq1_number], self.ids[seq2_number])] = line_list[seq1_number][seq2_number]
                self.matrix[(self.irangeds[seq2_number], self.ids[seq1_number])] = line_list[seq1_number][seq2_number]

    def dj(self, final_count):
        pass

    def __getitem__(self, item):
        '''
        Get a distance between two sequences
        :param item: a tuple of sequence names
        :return:
        '''
        return self.matrix[item]

    def add_row(self, new_id, dist_list):
        '''
        Add a sequence to the matrix
        :param new_id: sequence id
        :param dist_list: list of distances to all the other sequences
        :return:
        '''
        for pos in range(len(self.ids)):
            self.matrix[(new_id, self.ids[pos])] = dist_list[pos]
            self.matrix[(self.ids[pos], new_id)] = dist_list[pos]
        self.ids.append(new_id)

    def keys(self):
        '''
        Return the list of sequence pairs in this matrix object
        :return: A list of tuples
        '''
        return self.matrix.keys()