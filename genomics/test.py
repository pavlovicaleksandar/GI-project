from genomics.main import *
import pytest


# BWT & FM Index
@pytest.mark.parametrize("example, result", [
    ('', []),
    ('$', ['$']),
    ('banana$', ['banana$', 'anana$b', 'nana$ba', 'ana$ban', 'na$bana', 'a$banan', '$banana'])
])
def test_rotations(example, result):
    assert rotations(example) == result


@pytest.mark.parametrize("example, result", [
    ('', []),
    ('$', ['$']),
    ('banana$', ['$banana', 'a$banan', 'ana$ban', 'anana$b', 'banana$', 'na$bana', 'nana$ba'])
])
def test_bwm(example, result):
    assert bwm(example) == result


def test_make_occurrences_matrix():
    test_occ_matrix, test_tots = make_occurrences_matrix(bwt_via_bwm('banana$'))
    assert test_occ_matrix == {
        'a': [1, 1, 1, 1, 1, 2, 3],
        'n': [0, 1, 2, 2, 2, 2, 2],
        'b': [0, 0, 0, 1, 1, 1, 1],
        '$': [0, 0, 0, 0, 1, 1, 1]
    }
    assert test_tots == {'a': 3, 'n': 2, 'b': 1, '$': 1}


@pytest.mark.parametrize("example, result", [
    ({'a': 3, 'n': 2, 'b': 1, '$': 1}, {'$': (0, 0), 'a': (1, 3), 'b': (4, 4), 'n': (5, 6)})
])
def test_first_col(example, result):
    assert first_col(example) == result

@pytest.mark.parametrize("example_word, example_query, result_start, result_end", [
    ('banana$', 'ana', 2, 3),
    ('bananaban$', 'ban', 5, 6)
])
def test_calculate_start_end_range(example_word, example_query, result_start, result_end):
    test_occ_matrix, test_tots = make_occurrences_matrix(bwt_via_bwm(example_word))
    test_c = first_col(test_tots)
    test_start, test_end = calculate_start_end_range(test_c, test_occ_matrix, example_query)
    assert (test_start, test_end) == (result_start, result_end)

# Global Alignment
#
# assert scoring_matrix_inplace('C', '_') == -7
# assert scoring_matrix_inplace('C', 'G') == -2
# assert scoring_matrix_inplace('C', 'A') == -2
# assert scoring_matrix_inplace('C', 'T') == -1
# assert scoring_matrix_inplace('C', 'C') == 1
# test_x = 'TACGTCAGC'
# test_y = 'TATGTCATGC'
# assert scoring_matrix_inplace(test_x[6], test_y[6]) == 1
#
# # todo: probably has a better solution. Py knowledge critically low :(
# # scoring_points = {'M': 2, 'Ti': -3, 'Tv': -6, 'G': -9}
# # assert scoring_matrix_inplace('C', '_') == -9
# # assert scoring_matrix_inplace('C', 'G') == -6
# # assert scoring_matrix_inplace('C', 'A') == -6
# # assert scoring_matrix_inplace('C', 'T') == -3
# # assert scoring_matrix_inplace('C', 'C') == 2
# # assert scoring_matrix_inplace(test_x[6], test_y[6]) == 2
#
# expected_distance = numpy.array([[0, -7, -14, -21, -28, -35, -42, -49, -56, -63, -70],
#                                  [-7, 1, -6, -13, -20, -27, -34, -41, -48, -55, -62],
#                                  [-14, -6, 2, -5, -12, -19, -26, -33, -40, -47, -54],
#                                  [-21, -13, -5, 1, -6, -13, -18, -25, -32, -39, -46],
#                                  [-28, -20, -12, -6, 2, -5, -12, -19, -26, -31, -38],
#                                  [-35, -27, -19, -11, -5, 3, -4, -11, -18, -25, -32],
#                                  [-42, -34, -26, -18, -12, -4, 4, -3, -10, -17, -24],
#                                  [-49, -41, -33, -25, -19, -11, -3, 5, -2, -9, -16],
#                                  [-56, -48, -40, -32, -24, -18, -10, -2, 3, -1, -8],
#                                  [-63, -55, -47, -39, -31, -25, -17, -9, -3, 1, 0]])
# expected_alignment_score = 0
# distances, alignment_score = global_alignment(test_x, test_y, scoring_matrix_inplace)
# assert (expected_distance == distances).all()
# assert expected_alignment_score == alignment_score
#
# alignment, transcript = traceback(x, y, distances, scoring_matrix_inplace)
# expected_alignment = """TACGTCA_GC
# || |||| ||
# TATGTCATGC"""
# assert "MMRMMMMIMM" == transcript
# assert expected_alignment == alignment
#
# expected_fasta_len = 1
# expected_fastq_len = 20000
# test_references, test_reads = import_fasta_fastq()
# assert expected_fasta_len == len(test_references)
# assert expected_fastq_len == len(test_reads)
