import logging

logger = logging.getLogger(__name__)


def rotations(t):
    """ Return list of rotations of input string t """
    b = len(t)
    all_rotations = []
    for i in range(b):
        c = t[i:]+t[:i]
        all_rotations.append(c)
    return all_rotations


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
                return -1, -1
            start = c[reversed_query[ix]][0] + occ[reversed_query[ix]][start - 1]
            end = c[reversed_query[ix]][0] + occ[reversed_query[ix]][end] - 1
        return start, end
    except Exception as e:
        return -1, -1


def find_all_query_positions_in_word_via_suffix_arr(start, end, suff_arr):
    query_positions_in_word = []
    while start <= end:
        query_positions_in_word.append(suff_arr[start])
        start += 1
    return sorted(query_positions_in_word)
