__author__ = 'morozov'

from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import subprocess
from math import log
import copy
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
    def __init__(self, handle=None, find_distances = True):
        self.sequences = {}
        self.find_distances = find_distances
        self.matrix = DistanceMatrix()
        if handle is not None:
            self._load_fasta(handle)

    def _load_fasta(self, handle):
        '''
        Load sequences from multi-FASTA and populate DistanceMatrix
        :param handle: FASTA filehandle
        '''
        for record in SeqIO.parse(handle, format='fasta'):
            self.append(record)

    def __getitem__(self, item):
        '''
        Return a sequence with the given ID
        :param item: sequence ID
        '''
        return self.sequences[item]

    def append(self, item):
        '''
        Append a sequence to Sequence set and recalculate distance, if needed
        :param item: SeqRecord object
        '''
        self.sequences.update({item.name: item})
        if self.find_distances:
            dist_row = [scoredist(item, self[x]) for x in self.matrix.ids]
            self.matrix.add_row(item.id, dist_row)

    def __len__(self):
        return len(self.sequences)


class DistanceMatrix(object):
    '''
    Distance matrix class. Contains a list of IDs and a distance matrix for the same IDs
    '''
    def __init__(self, handle=None):
        self.matrix = {}
        self.ids = []
        if not handle is None:
            self.from_handle(handle)

    def from_handle(self, handle):
        '''
        Read distance matrix from filehandle
        :param handle: filehandle to read from
        :return:
        '''
        line_list = []
        for line in handle:
            l = line.split('\t')
            if len(l)>1:
                line_list.append(l)
        self.ids = [x[0] for x in line_list]
        for x in range(len(line_list)):
            del(line_list[x][0])
            line_list[x] = [float(j) for j in line_list[x]]
            #  IDs went to the attribute, no use for them here
        for seq1_number in range(len(line_list)):
            for seq2_number in range(len(line_list)):
            #  Order of the lines & items in line is the same as that of IDs
            #  Creating both of items to save __getitem__/__setitem__ hassle. Maybe add that later
                self.matrix[(self.ids[seq1_number], self.ids[seq2_number])] = line_list[seq1_number][seq2_number]
                self.matrix[(self.ids[seq2_number], self.ids[seq1_number])] = line_list[seq1_number][seq2_number]

    def to_handle(self, handle):
        '''
        Write distance matrix to filehandle
        :param handle: Filehandle to write matrix to
        :return:
        '''
        handle.write('{0}\n'.format(len(self.ids)))
        for seqid in self.ids:
            line = seqid+'\t'+'\t'.join((str(self.matrix[(seqid, x)]) for x in self.ids))+'\n'
            handle.write(line)

    def dj(self, final_count):
        '''
        Subsample from this distance matrix
        :param final_count: Number of sequences that should be in a final sample
        :return: DistanceMatrix object
        '''
        if final_count > len(self.matrix):
            raise ValueError('Cannot sample more sequences than there are in a dataset!')
        subset = DistanceMatrix()
        non_sampled = copy.deepcopy(self.ids)
        max_sum = 0
        sum_dists = {}
        leader = None
        for x in self.ids:
            sum_dists[x] = sum(self[(x, y)] for y in self.ids if x != y )
            if sum_dists[x] > max_sum:
                max_sum = sum_dists[x]
                leader = x
        subset.add_row(leader, [])
        non_sampled.remove(leader)
        #  We just added the item that is furthest from others as the first item in subsample
        #  From now on sum_dists is a sum of distances from every element IN SUBSET, starting from first one
        sum_dists = {x: self[(x, leader)] for x in non_sampled}
        for j in range(final_count - 1):
            max_sum = 0
            for candidate in non_sampled:
                #  Looking for furthest one again
                if sum_dists[candidate]>max_sum:
                    max_sum = sum_dists[candidate]
                    leader = candidate
            #  He's in leader now
            subset.add_row(leader, [self[(leader, x)] for x in subset.ids])
            non_sampled.remove(leader)
            for x in non_sampled:
                sum_dists[x]+=self[(x, leader)]
        return subset


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
        # sys.stderr.write(len(self.ids))
        #  print('{0} {1}'.format(new_id, ' '.join((str(x) for x in dist_list))))
        for pos in range(len(self.ids)):
            self.matrix[(new_id, self.ids[pos])] = dist_list[pos]
            self.matrix[(self.ids[pos], new_id)] = dist_list[pos]
        self.matrix[(new_id, new_id)] = 0
        self.ids.append(new_id)

    def keys(self):
        '''
        Return the list of sequence pairs in this matrix object
        :return: A list of tuples
        '''
        return self.matrix.keys()