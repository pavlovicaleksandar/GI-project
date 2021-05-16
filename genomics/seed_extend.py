from genomics.burrows import *
from genomics.global_alignment import *
import csv
import logging
import time


logger = logging.getLogger(__name__)


def reverse_and_complement(this):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(map(lambda x: complements[x], this[::-1]))


def get_seed(seed_length, read):
    return read[0: seed_length]


def generate_scoring_points(match, mismatch, gap):
    return {
        'M': match,
        'Ti': mismatch,
        'Tv': mismatch,
        'G': gap
    }


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

    matches = [0, 1, 2]
    mismatches = [-3, -2]
    gaps = [-7, -5]
    for match in matches:
        for mismatch in mismatches:
            for gap in gaps:
                scoring_points = generate_scoring_points(match, mismatch, gap)
                file_out = handle_reads(reads, seed_length, c, margin, occ_matrix, reference, scoring_points, suff_arr)
                write_results_to_csv_file2(f'results-{time.time()}-match-{match}-mismatch{mismatch}-gap{gap}.csv', file_out)


def write_results_to_csv_file2(file_name, results):
    with open(file_name, mode='w') as result_file:
        writer = csv.writer(result_file)
        writer.writerow(['start', 'end', 'alignment-score', 'transcription'])
        for result in results:
            writer.writerow(list(result))


def handle_reads(reads, seed_length, c, margin, occ_matrix, reference, scoring_points, suff_arr):
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
