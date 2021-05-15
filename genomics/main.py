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


def rank_bwt(bw):
    """ Given BWT string bw, return parallel list of B-ranks.  Also
        returns tots: map from character to # times it appears. """
    tots = dict()
    ranks = []
    for c in bw:
        if c not in tots: tots[c] = 0
        ranks.append(tots[c])
        tots[c] += 1
    return ranks, tots


def first_col(tots):
    """ Return map from character to the range of rows prefixed by the
    character. """
    first = {}
    totc = 0
    for c, count in sorted(tots.items()):
        first[c] = (totc, totc + count)
        totc += count
    return first



word = 'banana$'
last_column = bwt_via_bwm(word)
ranks, tots = rank_bwt(last_column)
last_column_with_ranks = zip(last_column, ranks)
pass