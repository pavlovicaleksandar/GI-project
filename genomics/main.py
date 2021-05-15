import argparse
import numpy
import argparse


# BWT & FM Index

def rotations(t):
    """ Return list of rotations of input string t """
    tt = t * 2
    return [tt[i:i + len(t)] for i in range(0, len(t))]


def bwm(t):
    """ Return lexicographically sorted list of t’s rotations """
    return sorted(rotations(t))


def bwt_via_bwm(t):
    """ Given T, returns BWT(T) by way of the BWM """
    return ''.join(map(lambda x: x[-1], bwm(t)))


def make_occurrences_matrix(bw):
    """ Given BWT string bw, returns a map of lists. Keys are
    characters and lists are cumulative # of occurrences up to and
    including the row. """
    tots = {}
    occ_matrix = {}
    for c in bw:
        if c not in tots:
            tots[c] = 0
            occ_matrix[c] = []
    for c in bw:
        tots[c] += 1
        for c in tots.keys():
            occ_matrix[c].append(tots[c])
    return occ_matrix, tots


def first_col(tots):
    """ Return map from character to the range of rows prefixed by the
    character. """
    first = {}
    totc = 0
    for c, count in sorted(tots.items()):
        first[c] = (totc, totc + count - 1)
        totc += count
    return first


def make_suffix_array(s):
    """ Given T return suffix array SA(T).  We use Python's sorted
        function here for simplicity, but we can do better. """
    satups = sorted([(s[i:], i) for i in range(len(s))])
    # Extract and return just the offsets
    return list(map(lambda x: x[1], satups))


def bwt_via_sa(t):
    """ Given T, returns BWT(T) by way of the suffix array. """
    bw = []
    for si in make_suffix_array(t):
        if si == 0:
            bw.append('$')
        else:
            bw.append(t[si - 1])
    return ''.join(bw)


def calculate_start_end_range(c, occ, query):
    try:
        reversed_query = query[::-1]
        start, end = c[reversed_query[0]]
        for ix in range(1, len(reversed_query)):
            if start > end:
                raise ValueError
            start = c[reversed_query[ix]][0] + occ[reversed_query[ix]][start - 1]
            end = c[reversed_query[ix]][0] + occ[reversed_query[ix]][end] - 1
        return start, end
    except Exception as e:
        raise e


def find_all_query_positions_in_word_via_suffix_arr(start, end, suff_arr):
    query_positions_in_word = []
    while start <= end:
        query_positions_in_word.append(suff_arr[start])
        start += 1
    return sorted(query_positions_in_word)


# Global alignment - Needle - Wunsch + scoring table

def scoring_matrix_inplace(this_c, that_c):
    """ Generate a scoring matrix such that matching has the highest value, transitions are penalised,
    transversions are penalised even more, and gaps are penalised the most. Reminder: A-G, C-T transitions,
    other combinations ase transversions. Scoring points are held globally. """
    # todo: a better solution for scoring points probably exists
    if this_c == that_c:
        return scoring_points['M']
    if this_c == '_' or that_c == '_':
        return scoring_points['G']
    maxb, minb = max(this_c, that_c), min(this_c, that_c)
    if (minb == 'A' and maxb == 'G') or (minb == 'C' and maxb == 'T'):
        return scoring_points['Ti']
    else:
        return scoring_points['Tv']


def global_alignment(this, that, scoring_matrix):
    """ Firstly, fill in V(i, 0) and V(0, j) with gap values, then the whole distance matrix. """
    distance_matrix = numpy.zeros((len(this) + 1, len(that) + 1), dtype=int)

    for i in range(1, len(this) + 1):
        distance_matrix[i, 0] = distance_matrix[i - 1, 0] + scoring_matrix(this[i - 1], '_')
    for j in range(1, len(that) + 1):
        distance_matrix[0, j] = distance_matrix[0, j - 1] + scoring_matrix('_', that[j - 1])

    for i in range(1, len(this) + 1):
        for j in range(1, len(that) + 1):
            distance_matrix[i, j] = max(distance_matrix[i - 1, j] + scoring_matrix(this[i - 1], '_'),
                                        distance_matrix[i, j - 1] + scoring_matrix('_', that[j - 1]),
                                        distance_matrix[i - 1, j - 1] + scoring_matrix(this[i - 1], that[j - 1]))

    # function returns table and global alignment score
    # alignment score is in cell (n,m) of the matrix
    return distance_matrix, distance_matrix[len(this), len(that)]


def traceback(this, that, distance_matrix, scoring_matrix):
    # initializing starting position cell(n,m)
    i = len(this)
    j = len(that)

    # initializing strings we use to represent alignments in x, y, edit transcript and global alignment
    ax, ay, am, tr = '', '', '', ''

    # exit condition is when we reach cell (0,0)
    while i > 0 or j > 0:

        # calculating diagonal, horizontal and vertical scores for current cell
        d, v, h = -100, -100, -100

        if i > 0 and j > 0:
            delta = 1 if this[i - 1] == that[j - 1] else 0
            d = distance_matrix[i - 1, j - 1] + scoring_matrix(this[i - 1], that[j - 1])  # diagonal movement
        if i > 0:
            v = distance_matrix[i - 1, j] + scoring_matrix(this[i - 1], '_')  # vertical movement
        if j > 0:
            h = distance_matrix[i, j - 1] + scoring_matrix('_', that[j - 1])  # horizontal movement

        # backtracking to next (previous) cell
        if d >= v and d >= h:
            ax += this[i - 1]
            ay += that[j - 1]
            if delta == 1:
                tr += 'M'
                am += '|'
            else:
                tr += 'R'
                am += ' '
            i -= 1
            j -= 1
        elif v >= h:
            ax += this[i - 1]
            ay += '_'
            tr += 'D'
            am += ' '
            i -= 1
        else:
            ay += that[j - 1]
            ax += '_'
            tr += 'I'
            am += ' '
            j -= 1

    alignment = '\n'.join([ax[::-1], am[::-1], ay[::-1]])
    return alignment, tr[::-1]


def define_arg_parser():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--fasta', type=str, required=True, help="Path to the fasta file")
    arg_parser.add_argument('--fastq', type=str, required=True, help="Path to the fastq file")
    arg_parser.add_argument('--margin', default=2)
    arg_parser.add_argument('--match', type=int, nargs="+", default=[0, 1, 2])
    arg_parser.add_argument('--mismatch', type=int, nargs="+", default=[-3, -2])
    # arg_parser.add_argument('--Ti', type=int, nargs="+", default=-1)
    # arg_parser.add_argument('--Tv', type=int, nargs="+", default=-3)
    arg_parser.add_argument('--gap', type=int, nargs="+", default=[-7, -5])
    arg_parser.add_argument('--seed-length', type=int, default=10)
    return arg_parser


def import_or_generate_scoring_points():
    # todo : create this in define_arg_parser()
    parser = argparse.ArgumentParser(description='Process scoring points.')
    parser.add_argument('integers', type=int, nargs='*',
                        help='4 scoring points', default=[1, -1, -2, -7])
    args = parser.parse_args()
    scoring_points_input = sorted(args.integers)
    return {'M': scoring_points_input[3], 'Ti': scoring_points_input[2], 'Tv': scoring_points_input[1],
            'G': scoring_points_input[0]}


# args = define_arg_parser().parse_args()

scoring_points = import_or_generate_scoring_points()
word = 'banana$'
query = 'ana'
last_column = bwt_via_bwm(word)
occ_matrix, tots = make_occurrences_matrix(last_column)
c = first_col(tots)
start, end = calculate_start_end_range(c, occ_matrix, query)
bwt_via_sa(word)
suff_arr = make_suffix_array(word)
query_pos = find_all_query_positions_in_word_via_suffix_arr(start, end, suff_arr)

x = 'TACGTCAGC'
y = 'TATGTCATGC'
distances, alignment_score = global_alignment(x, y, scoring_matrix_inplace)
alignment, transcript = traceback(x, y, distances, scoring_matrix_inplace)
