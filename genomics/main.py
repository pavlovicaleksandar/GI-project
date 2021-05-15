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


word = 'banana$'
query = 'ana'
last_column = bwt_via_bwm(word)
occ_matrix, tots = make_occurrences_matrix(last_column)
c = first_col(tots)
start, end = calculate_start_end_range(c, occ_matrix, query)
