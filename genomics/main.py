import argparse
import numpy


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

def scoring_matrix_inplace(a, b):
    """ Generate a scoring matrix such that matching has the highest value, transitions are penalised,
    transversions are penalised even more, and gaps are penalised the most. Reminder: A-G, C-T transitions,
    other combinations ase transversions. Scoring points are held globally. """
    # todo: a better solution for scoring points probably exists
    if a == b:
        return scoring_points['M']
    if a == '_' or b == '_':
        return scoring_points['G']
    maxb, minb = max(a, b), min(a, b)
    if (minb == 'A' and maxb == 'G') or (minb == 'C' and maxb == 'T'):
        return scoring_points['Ti']
    else:
        return scoring_points['Tv']


def global_alignment(x, y, s):
    """ Firstly, fill in V(i, 0) and V(0, j) with gap values, then the whole distance matrix. """
    distance_matrix = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)

    for i in range(1, len(x) + 1):
        distance_matrix[i, 0] = distance_matrix[i - 1, 0] + s(x[i - 1], '_')
    for j in range(1, len(y) + 1):
        distance_matrix[0, j] = distance_matrix[0, j - 1] + s('_', y[j - 1])

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            distance_matrix[i, j] = max(distance_matrix[i - 1, j] + s(x[i - 1], '_'),
                                        distance_matrix[i, j - 1] + s('_', y[j - 1]),
                                        distance_matrix[i - 1, j - 1] + s(x[i - 1], y[j - 1]))

    # function returns table and global alignment score
    # alignment score is in cell (n,m) of the matrix
    return distance_matrix, distance_matrix[len(x), len(y)]


def import_or_generate_scoring_points():
    parser = argparse.ArgumentParser(description='Process scoring points.')
    parser.add_argument('integers', type=int, nargs='*',
                        help='4 scoring points', default=[1, -1, -2, -7])
    args = parser.parse_args()
    scoring_points_input = sorted(args.integers)
    return {'M': scoring_points_input[3], 'Ti': scoring_points_input[2], 'Tv': scoring_points_input[1],
                      'G': scoring_points_input[0]}


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

test_x = 'TACGTCAGC'
test_y = 'TATGTCATGC'
distances, alginment = global_alignment(test_x, test_y, scoring_matrix_inplace)
print(distances)