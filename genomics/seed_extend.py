from genomics.burrows import *
from genomics.global_alignment import *
import logging

logger = logging.getLogger(__name__)


def reverse_and_complement(this):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(map(lambda x: complements[x], this[::-1]))


def get_seed(seed_length, read):
    return read[0: seed_length]


# make coarse then do fine refine
def seed_and_extend(scoring_points, references, reads, margin=3, seed_length=10, reference_index=0):

    # init
    reference = references[reference_index]+'$'

    logger.info('Starting bwt via bwm')
    bwt_reference = bwt_via_bwm(reference)
    logger.info('Finished bwt via bwm')

    occ_matrix, tots = make_occurrences_matrix(bwt_reference)
    c = first_col(tots)
    suff_arr = make_suffix_array(reference)
    result = []
    for read in reads:
        # for each read extract substrins - seeds and reverse complement
        rc_read = reverse_and_complement(read)

        seed = get_seed(seed_length, read)
        rc_seed = get_seed(seed_length, rc_read)

        result.extend(extract(c, margin, occ_matrix, read, reference, seed, seed_length, scoring_points, suff_arr))
        result.extend(extract(c, margin, occ_matrix, rc_read, reference, rc_seed, seed_length, scoring_points, suff_arr))

    return sorted(result, key=lambda tup: tup[2], reverse=True)


def extract(c, margin, occ_matrix, read, reference, seed, seed_length, scoring_points, suff_arr):
    start, end = calculate_start_end_range(c, occ_matrix, seed)
    if (start, end) == (-1, -1):
        return []
    seed_positions = find_all_query_positions_in_word_via_suffix_arr(start, end, suff_arr)
    # ref     "ctagtcgt agctagctgatcg"
    # read    "gctgt cgtc"
    # s       "gctgt"
    result = []
    for start_of_seed in seed_positions:
        # start/end of reference
        start_pos = start_of_seed + seed_length
        end_pos = start_pos + len(read) - seed_length + margin + 1
        end_pos = len(reference) if end_pos > len(reference) else end_pos

        distances, alignment_score = global_alignment(read[seed_length:], reference[start_pos:end_pos],
                                                      scoring_points)
        alignment, transcript = traceback(read[seed_length:], reference[start_pos:end_pos], distances,
                                          scoring_points)
        result.append((start, end, alignment_score, transcript))
    return result
