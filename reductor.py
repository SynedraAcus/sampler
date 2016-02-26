__author__ = 'morozov'

from Bio import SeqIO, Alphabet, AlignIO
from Bio.SubsMat import MatrixInfo
from io import StringIO
import subprocess
from math import log
import copy
import re
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

def make_aminoacid_matrix(fasta_handle):
    """
    Create a distance matrix for all sequences within FASTA filehandle
    Assumes sequences to be all aminoacid
    :param fasta_handle: filehandle
    :return: DistanceMatrix
    """
    ret = DistanceMatrix()
    seq_list = []
    for sequence in SeqIO.parse(fasta_handle, format='fasta'):
        seq_list.append(sequence)
    while len(seq_list)>1:
        sequence1 = seq_list.pop()
        for sequence2 in seq_list:
            ret[(sequence1.id, sequence2.id)] = scoredist(sequence1, sequence2)
    # adding zeros
    ret[(seq_list[0].id, seq_list[0].id)] = 0.0
    for x in ret.ids:
        ret[(x, x)] = 0.0
    return ret

class MatrixFactory(object):
    '''
    Distance matrix factory. Take        print(aln)
s a FASTA handle and AA/nucleotide boolean as an input, produces DistanceMatrix
    '''
    #  EMBOSS call lines defined as class variables

    #  needle call line. format(fasta_file, asequenceID, bsequenceID)
    NEEDLE_CALL = 'needle -asequence {0}:{1} -bsequence {0}:{2} -gapopen 10.0 -gapextend 0.5 -aformat3 fasta -outfile stdout'
    def __init__(self):
        pass

    def create_aminoacid_matrix(self, fasta_file=None):
        """
        Create a distance matrix for all sequences within FASTA filehandle
        Assumes sequences to be all aminoacid
        :param fasta_file: FASTA filename
        :return: DistanceMatrix
        """
        ret = DistanceMatrix()
        seq_list = []
        for sequence in SeqIO.parse(open(fasta_file), format='fasta'):
            seq_list.append(sequence)
        while len(seq_list)>1:
            sequence1 = seq_list.pop()
            for sequence2 in seq_list:
                ret[(sequence1.id, sequence2.id)] = self._scoredist(fasta_file, sequence1.id, sequence2.id)
        # adding zeros
        ret[(seq_list[0].id, seq_list[0].id)] = 0.0
        for x in ret.ids:
            ret[(x, x)] = 0.0
        return ret
        pass

    def _scoredist(self, fasta, id1, id2, matrix=matrix, correction=1.337):
        EXPECT = -0.5209
        byte_aln = subprocess.check_output(self.NEEDLE_CALL.format(fasta, id1, id2),
                                           shell=True,
                                           stderr=subprocess.DEVNULL)
        stream = StringIO(byte_aln.decode())
        aln = AlignIO.read(stream, format='fasta')
        score = 0
        score_a = 0
        score_b = 0
        length = 0
        for col in range(len(aln[0])):
            if '-' in aln[:, col]:
                continue  # Ignore gapped columns
            score_a += matrix[(aln[0][col], aln[0][col])]
            score_b += matrix[(aln[1][col], aln[1][col])]
            score += matrix[(aln[1][col], aln[0][col])]
            length += 1.0
        sigma_n = score - length * EXPECT
        sigma_un = (score_a + score_b)/2 - length * EXPECT
        dist = -1*log(sigma_n/sigma_un)*100*correction
        return dist

    def create_nucleotide_matrix(self, fasta_file=None):
        pass


class DistanceMatrix(object):
    '''
    Distance matrix class. Contains a list of IDs and a distance matrix for the same IDs
    '''
    def __init__(self, handle=None):
        self.matrix = {}
        self.ids = []
        if not handle is None:
            self.read(handle)

    #  I/O & matrix creation methods

    def read(self, filehandle):
        """
        Read distance matrix from filehandle. Expects EMBOSS format (upper-right).
        Works using some voodoo
        :param filehandle:
        :return:
        """
        id_dict = {}
        num_matrix = {}  # Tmp matrix w/numeric IDs
        for line in filehandle:
            if '\t' not in line:
                #  Skip lines that are not tab-separated, headers and such
                continue
            line=line.rstrip()
            # l_arr = re.split('\s+', line)
            l_arr = line.split('\t')
            if l_arr[1].lstrip() == '1':
                #  Skip numbers line, it's useless
                continue
            ids = l_arr.pop(-1)
            (al_id, num_id) = re.split('\s+', ids)
            num_id = int(num_id)
            id_dict[num_id] = al_id
            for num in range(num_id, len(l_arr)-1):
                num_matrix[(num, num_id)] = float(l_arr[num])
                num_matrix[(num_id, num)] = float(l_arr[num])
        #  Changing num_matrix IDs to proper names
        #  And writing to object
        self.ids = list(id_dict.values())
        for num_id1 in id_dict.keys():
            for num_id2 in id_dict.keys():
                if num_id1 != num_id2:
                    self[(id_dict[num_id1], id_dict[num_id2])] = num_matrix[(num_id1, num_id2)]
                else:
                    self[(id_dict[num_id2], id_dict[num_id2])] = 0.0

    def write(self, fh):
        """
        Write distance matrix to filehandle using EMBOSS format (upper-right)
        :param filehandle:
        :return:
        """
        print('Distance matrix\n---------------\n\nReduced by sampler\n', file=fh)
        print('\t'+'\t'.join(str(x) for x in range(1, len(self.ids)+1)), file=fh)
        num_code = {self.ids[x-1]: x for x in range(len(self.ids)+1)}
        for id in self.ids:
            fh.write('\t'*num_code[id])
            for x in range(num_code[id]-1, len(self.ids)):
                fh.write('{0:.2f}'.format(self[(id, self.ids[x])])+'\t')
            fh.write('\t{0} {1}\n '.format(id, num_code[id]))

    def dj(self, final_count):
        """
        Return a reduced ID list using distant joining algorithm
        Reduces the initial matrix ID list to final_count elements
        :param final_count: int
        :return: list of strings
        """
        if final_count > len(self.ids):
            raise ValueError ('Cannot reduce matrix to more elements than it has')
        reduced_list = []
        not_sampled = copy.deepcopy(self.ids)
        # Finding the initial object that is furthest from all in not_sampled
        max_dist = 0.0
        leader = None
        for candidate in not_sampled:
            dist = sum((self[(candidate, x)] for x in not_sampled))
            if dist >= max_dist:
                leader = candidate
                max_dist = dist
        not_sampled.remove(leader)
        reduced_list.append(leader)
        sum_dist = {x: self[(leader, x)] for x in self.ids}
        #  Iterative addition of the elements
        while len(reduced_list) < final_count:
            leader = max(sum_dist.items(), key=lambda a:a[1])
            reduced_list.append(leader[0])
            del sum_dist[leader[0]]
            not_sampled.remove(leader[0])
            for i in not_sampled:
                sum_dist[i] += self[(leader[0], i)]
        return reduced_list

    def submatrix(self, ids):
        """
        Create a submatrix that contains only the items in the list
        :param ids: list of strings
        :return:
        """
        ret = DistanceMatrix()
        id_list = copy.deepcopy(ids)
        while len(id_list)>1:
            id = id_list.pop()
            ret[(id, id)] = 0.0
            for x in id_list:
                ret[(id, x)] = self[(id, x)]
        #  Add zero distance from last element to itself
        ret[(id_list[0], id_list[0])] = 0.0
        return ret

    #  Dict-like behaviour

    def __getitem__(self, item):
        '''
        Get a distance between two sequences
        :param item: a tuple of sequence names
        :return:
        '''
        assert type(item) is tuple
        assert len(item) == 2
        if item in self.keys():
            return self.matrix[item]
        elif (item[1], item[0]) in self.keys():
            return self.matrix[(item[1], item[0])]
        else:
            raise KeyError('Invalid matrix key: {0}'.format(item))

    def __setitem__(self, key, value):
        '''
        Add an item to the matrix
        :param key: a 2-item tuple
        :param item: float
        '''
        assert type(key) is tuple
        assert len(key) == 2
        assert type(value) is float
        # Check if ID list of matrix is incomplete and correct it, if so
        if key[0] not in self.ids:
            self.ids.append(key[0])
        if key[1] not in self.ids:
            self.ids.append(key[1])
        #  Check if there is a place w/oppposite ID order
        #  Not doing so doubles use of space
        if (key[1], key[0]) in self.keys():
            self.matrix[(key[1], key[0])] = value
            return
        self.matrix[key] = value


    def keys(self):
        '''
        Return the list of sequence pairs in this matrix object
        :return: A list of tuples
        '''
        return self.matrix.keys()