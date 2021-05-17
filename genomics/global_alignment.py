import numpy
import logging

logger = logging.getLogger(__name__)


def generate_scoring_points(match, mismatch, gap):
    return {
        'M': match,
        'Ti': mismatch,
        'Tv': mismatch,
        'G': gap
    }


def scoring_matrix_inplace(this_c, that_c, scoring_points):
    """ Generate a scoring matrix such that matching has the highest value, transitions are penalised,
    transversions are penalised even more, and gaps are penalised the most. Reminder: A-G, C-T transitions,
    other combinations ase transversions. """
    if this_c == that_c:
        return scoring_points['M']
    if this_c == '_' or that_c == '_':
        return scoring_points['G']
    maxb, minb = max(this_c, that_c), min(this_c, that_c)
    if (minb == 'A' and maxb == 'G') or (minb == 'C' and maxb == 'T'):
        return scoring_points['Ti']
    else:
        return scoring_points['Tv']


def global_alignment(this, that, scoring_points):
    """ Firstly, fill in V(i, 0) and V(0, j) with gap values, then the whole distance matrix. """
    distance_matrix = numpy.zeros((len(this) + 1, len(that) + 1), dtype=int)

    for i in range(1, len(this) + 1):
        distance_matrix[i, 0] = distance_matrix[i - 1, 0] + scoring_matrix_inplace(this[i - 1], '_', scoring_points)
    for j in range(1, len(that) + 1):
        distance_matrix[0, j] = distance_matrix[0, j - 1] + scoring_matrix_inplace('_', that[j - 1], scoring_points)

    for i in range(1, len(this) + 1):
        for j in range(1, len(that) + 1):
            distance_matrix[i, j] = max(
                distance_matrix[i - 1, j] + scoring_matrix_inplace(this[i - 1], '_', scoring_points),
                distance_matrix[i, j - 1] + scoring_matrix_inplace('_', that[j - 1], scoring_points),
                distance_matrix[i - 1, j - 1] + scoring_matrix_inplace(this[i - 1], that[j - 1], scoring_points)
            )

    # function returns table and global alignment score
    # alignment score is in cell (n,m) of the matrix
    return distance_matrix, distance_matrix[len(this), len(that)]


def traceback(this, that, distance_matrix, scoring_points):
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
            d = distance_matrix[i - 1, j - 1] + scoring_matrix_inplace(this[i - 1], that[j - 1], scoring_points)  # diagonal movement
        if i > 0:
            v = distance_matrix[i - 1, j] + scoring_matrix_inplace(this[i - 1], '_', scoring_points)  # vertical movement
        if j > 0:
            h = distance_matrix[i, j - 1] + scoring_matrix_inplace('_', that[j - 1], scoring_points)  # horizontal movement

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
