__author__ = 'morozov'

from Bio import SeqIO, Alphabet, AlignIO
from Bio.SubsMat import MatrixInfo
import numpy as np
from io import StringIO
import subprocess
from math import log
import copy
import re

#  Converting BLOSUM 62 to square form
matrix = MatrixInfo.blosum62
keys = list(matrix.keys())
for a in keys:
    matrix[(a[1], a[0])] = matrix[a]


class MatrixFactory(object):
    """
    Distance matrix factory.
    Supports reading distance matrices from EMBOSS-formatted files and
    building them from sequence collections using EMBOSS needle alignment
    and either Scoredist (aminoacid) or kimura2param (nucleotide) distances.
    """

    #  needle call line. format(fasta_file, asequenceID, bsequenceID)
    NEEDLE_CALL = 'needle -asequence {0}:{1} -bsequence {0}:{2} -gapopen 10.0 -gapextend 0.5 -aformat3 fasta -outfile stdout'
    def __init__(self):
        pass

    def read(self, filehandle, format='emboss'):
        """
        Read distance matrix from filehandle in either emboss or clustalo format
        Works using some voodoo
        :param filehandle: Filehandle for the matrix file
        :param format: str, one of 'emboss' or 'clustalo'
        :return: DistanceMatrix
        """
        if format == 'emboss':
            return self.read_emboss(filehandle)
        elif format == 'clustalo':
            return self.read_clustalo(filehandle)
        else:
            raise ValueError('Invalid format')
    
    def read_clustalo(self, filehandle):
        """
        Read distance matrix from clustalo file
        :param filehandle:
        :return:
        """
        r = DistanceMatrix()
        count = int(filehandle.readline().rstrip('\n'))
        r.init_nd(count)
        counter = 0
        for line in filehandle:
            arr = line.rstrip().split()
            r.indices[arr[0]] = counter
            for index, value in enumerate(arr[1:]):
                r.array[(index, counter)] = float(value)
                if index == counter:
                    break
                    # Loading only lower-left part of the matrix to reduce
                    # memory. DistanceMatrix.__getitem__ would deal with NaN's
            counter += 1
        return r
        
    def read_emboss(self, filehandle):
        ret = DistanceMatrix()
        # id_dict = {}
        num_matrix = {}  # Tmp matrix w/numeric IDs
        for line in filehandle:
            if '\t' not in line:
                #  Skip lines that are not tab-separated, headers and such
                continue
            line = line.rstrip()
            l_arr = line.split('\t')
            if l_arr[1].lstrip() == '1':
                #  Counts start with 1, not zero. This is not important, as we
                # can use any numeric IDS so long as we write matrix.indices correctly
                array_size = l_arr[-1]
                ret.init_nd(int(array_size)+1)
                # ID line should be used only once, to define the array size.
                # The rest of procedure is a dist parser
                continue
            #  Get item ID from the last element and add it to distance matrix
            # 'indices' property
            ids = l_arr.pop(-1)
            (al_id, num_id) = re.split('\s+', ids)
            num_id = int(num_id)
            ret.indices[al_id] = num_id
            # id_dict[num_id] = al_id
            for num in range(num_id, len(l_arr)-1):
                ret.array[num, num_id] = float(l_arr[num])
                num_matrix[(num, num_id)] = float(l_arr[num])
                num_matrix[(num_id, num)] = float(l_arr[num])
        return ret

    def _build_alignment(self, fasta_file, id1, id2):
        '''
        Internal function that calls EMBOSS needle on the default params to build
        a pairwise alignment
        :param fasta_file: FASTA filename
        :param id1: ID of sequence 1
        :param id2: ID of sequence 2
        :return: Bio.Align instance
        '''
        screened_id1 = re.sub(r'([\|\[\]\\\/ ])', lambda a: '\{}'.format(a.groups(1)[0]), id1)
        screened_id2 = re.sub(r'([\|\[\]\\\/ ])', lambda a: '\{}'.format(a.groups(1)[0]), id2)
        call = self.NEEDLE_CALL.format(fasta_file, screened_id1, screened_id2)
        # print(screened_call)
        # quit()
        byte_aln = subprocess.check_output(call,
                                           shell=True,
                                           stderr=subprocess.DEVNULL)
        stream = StringIO(byte_aln.decode())
        aln = AlignIO.read(stream, format='fasta')
        return aln

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
        ret.init_nd(len(seq_list))
        ret.indices = {seq_list[x].id: x for x in range(len(seq_list))}
        while len(seq_list)>1:
            sequence1 = seq_list.pop()
            for sequence2 in seq_list:
                ret[(sequence1.description, sequence2.description)] =\
                    self._scoredist(fasta_file, sequence1.description, sequence2.description)
        # adding zeros
        ret[(seq_list[0].id, seq_list[0].id)] = 0.0
        for x in ret.indices.keys():
            ret[(x, x)] = 0.0
        return ret
        pass

    def _scoredist(self, fasta, id1, id2, matrix=matrix, correction=1.337):
        # 1.337 is not a 1337 joke, it's an actual correction value for blosum62
        EXPECT = -0.5209
        aln = self._build_alignment(fasta, id1, id2)
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

    #  This differs from create_aminoacid_matrix only in the function used.
    #  These two methods should be refactored into a single one
    def create_nucleotide_matrix(self, fasta_file=None):
        '''
        Return a Distance Matrix built for all the DNA sequences in fasta_file
        :param fasta_file: FASTA filename
        :return: DistanceMatrix
        '''
        ret = DistanceMatrix()
        seq_list = []
        for sequence in SeqIO.parse(open(fasta_file), format='fasta'):
            seq_list.append(sequence)
        ret.init_nd(len(seq_list))
        ret.indices = {seq_list[x].id: x for x in range(len(seq_list))}
        while len(seq_list)>1:
            sequence1 = seq_list.pop()
            for sequence2 in seq_list:
                ret[(sequence1.id, sequence2.id)] = self._kimura(fasta_file, sequence1.description, sequence2.description)
        # adding zeros
        ret[(seq_list[0].id, seq_list[0].id)] = 0.0
        for x in ret.ids:
            ret[(x, x)] = 0.0
        return ret
        pass

    def _kimura(self, fasta, id1, id2):
        aln = self._build_alignment(fasta, id1, id2)
        p = 0
        q = 0
        transition_sets=[set('A', 'T'), set('C', 'G')]
        for col in range(len(aln[0])):
            if set((aln[0][col], aln[1][col])) in transition_sets:
                p += 1
            else:
                q += 1
        w1 = 1 - 2*p - q
        w2 = 1 - 2*q
        d = log(w1)/2.0 - log(w2)/4.0
        return d


class DistanceMatrix(object):
    '''
    Distance matrix class.
    Contains a list of IDs and a distance matrix for the same IDs. Can produce
    a submatrix when given an ID list, or calculate ID lists of a desired size
    using DJ algorithm
    '''
    def __init__(self, handle=None):
        #  ID-to-index dictionary
        self.indices={}
        #  ndarray with the matrix itself. As it cannot be defined before we know how many sequences there are,
        #  the attribute is declared with None. Actually defining array and loading data is done by factory
        self.array = None
        if not handle is None:
            self.read(handle)

    def init_nd(self, size):
        self.array = np.empty(shape=(size, size), dtype=np.float32)
        self.array.fill(np.nan)
    #  I/O & matrix creation methods

    def write(self, fh):
        """
        Write distance matrix to filehandle using EMBOSS format (upper-right)
        :param filehandle:
        :return:
        """
        print('Distance matrix\n---------------\n\nReduced by sampler\n', file=fh)
        #  The following is a relic from the time when there was a self.ids list attribute, and probably should be rewritten
        ids = list(self.indices.keys())
        print('\t'+'\t'.join(str(x) for x in range(1, len(ids)+1)), file=fh)
        num_code = {ids[x-1]: x for x in range(len(ids)+1)}
        for id in ids:
            fh.write('\t'*num_code[id])
            for x in range(num_code[id]-1, len(ids)):
                fh.write('{0:.2f}'.format(self[(id, ids[x])])+'\t')
            fh.write('\t{0} {1}\n '.format(id, num_code[id]))

    def dj(self, final_count):
        """
        Return a reduced ID list using distant joining algorithm
        Reduces the initial matrix ID list to final_count elements
        :param final_count: int
        :return: list of strings
        """
        if final_count > len(self.indices.keys()):
            raise ValueError('Cannot reduce matrix to more elements than it has')
        reduced_list = []
        not_sampled = list(self.indices.keys())
        # Finding the initial object that is furthest from all in not_sampled
        max_dist = 0.0
        leader = not_sampled[0]
        for candidate in not_sampled[1:]:
            dist = sum((self[(candidate, x)] for x in not_sampled))
            if dist >= max_dist:
                leader = candidate
                max_dist = dist
        not_sampled.remove(leader)
        reduced_list.append(leader)
        minima = {x: self[(leader, x)] for x in self.indices.keys()}
        #  Iterative addition of the elements
        # For each sequence in not_sampled, take the minimum distance to ones
        # in set. EG if there are three seqs in reduced_list, and for x the
        # distances are [4,4,3], its minimum is 3.
        while len(reduced_list) < final_count:
            leader = max(minima.items(), key=lambda a: a[1])
            reduced_list.append(leader[0])
            del minima[leader[0]]
            not_sampled.remove(leader[0])
            for i in not_sampled:
                d = self[(leader[0], i)]
                if d < minima[i]:
                    minima[i] = d
        return reduced_list

    def submatrix(self, ids):
        """
        Create a submatrix that contains only the items in the list
        :param ids: list of strings
        :return:
        """
        ret = DistanceMatrix()
        ret.init_nd(len(ids))
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
        index = [self.indices[x] for x in item]
        if np.isnan(self.array[index[0], index[1]]):
            return self.array[index[1], index[0]]
        else:
            return self.array[index[0], index[1]]

    def __setitem__(self, key, value):
        '''
        Add an item to the matrix
        :param key: a 2-item tuple
        :param item: float
        '''
        assert type(key) is tuple
        assert len(key) == 2
        # assert type(value) is float
        if key[0] not in self.indices.keys():
            self.indices.update({key[0]: len(self.indices.keys())})
        if key[1] not in self.indices.keys():
            self.indices.update({key[1]: len(self.indices.keys())})
        self.array[self.indices[key[0]], self.indices[key[1]]] = value
