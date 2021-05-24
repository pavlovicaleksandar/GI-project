from genomics.burrows import *
from genomics.global_alignment import *


def reverse_and_complement(this):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(map(lambda x: complements[x], this[::-1]))


def get_seed(seed_length, read):
    return read[0: seed_length]


# make coarse then do fine refine
def seed_and_extend(reference, reads, occ_matrix, c, suff_arr, scoring_points, margin=3, seed_length=10):
    result = []
    for read in reads:
        # for each read extract substrings - seeds and reverse complement
        rc_read = reverse_and_complement(read)

        seed = get_seed(seed_length, read)
        rc_seed = get_seed(seed_length, rc_read)

        # Added direction to be able to have more info in the result
        fwd_result = extract(c, margin, occ_matrix, read, reference, seed,
                             seed_length, scoring_points, suff_arr, dir="fwd")
        rc_result = extract(c, margin, occ_matrix, rc_read, reference, rc_seed,
                            seed_length, scoring_points, suff_arr, dir="rc")

        alignment_score_index_in_result_tuple = 3
        if fwd_result[alignment_score_index_in_result_tuple] >= rc_result[alignment_score_index_in_result_tuple] \
                and fwd_result[alignment_score_index_in_result_tuple] != -1000:
            result.append(fwd_result)
        elif rc_result[alignment_score_index_in_result_tuple] != -1000:
            # skip if both are -1000
            result.append(rc_result)

    return sorted(result, key=lambda tup: tup[alignment_score_index_in_result_tuple], reverse=True)


def extract(c, margin, occ_matrix, read, reference, seed, seed_length, scoring_points, suff_arr, dir):
    # Empirically discovered...
    best_alignment_score = -1000

    start, end = calculate_start_end_range(c, occ_matrix, seed)
    result = (start, end, dir, best_alignment_score, '')
    if (start, end) == (-1, -1):
        # return default result to avoid extracting data from empty tuple 
        return result
    
    # ref     "ctgctgt agctagctgatcg"
    # read    "gctgt cgtc"
    # s       "gctgt"

    seed_positions = find_all_query_positions_in_word_via_suffix_arr(start, end, suff_arr)    
    for start_of_seed in seed_positions:
        # start/end of reference
        start_pos = start_of_seed + seed_length
        end_pos = start_pos + len(read) - seed_length + margin + 1
        end_pos = len(reference) if end_pos > len(reference) else end_pos

        distances, alignment_score = global_alignment(read[seed_length:], reference[start_pos:end_pos],
                                                      scoring_points)
        alignment, transcript = traceback(read[seed_length:], reference[start_pos:end_pos], distances,
                                          scoring_points)

        if best_alignment_score < alignment_score:
            result = (start, end, dir, alignment_score, transcript)

    return result
