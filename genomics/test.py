from genomics.main import *

# BWT & FM Index

assert rotations('') == []
assert rotations('$') == ['$']
assert rotations('banana$') == ['banana$', 'anana$b', 'nana$ba', 'ana$ban', 'na$bana', 'a$banan', '$banana']

assert bwm('') == []
assert bwm('$') == ['$']
assert bwm('banana$') == ['$banana', 'a$banan', 'ana$ban', 'anana$b', 'banana$', 'na$bana', 'nana$ba']

assert bwt_via_bwm('banana$') == 'annb$aa'
test_occ_matrix, test_tots = make_occurrences_matrix(bwt_via_bwm('banana$'))
assert test_occ_matrix == {'a': [1, 1, 1, 1, 1, 2, 3], 'n': [0, 1, 2, 2, 2, 2, 2], 'b': [0, 0, 0, 1, 1, 1, 1],
                           '$': [0, 0, 0, 0, 1, 1, 1]}
assert test_tots == {'a': 3, 'n': 2, 'b': 1, '$': 1}
assert first_col(test_tots) == {'$': (0, 0), 'a': (1, 3), 'b': (4, 4), 'n': (5, 6)}

test_word = 'banana$'
test_query = 'ana'
test_occ_matrix, test_tots = make_occurrences_matrix(bwt_via_bwm(test_word))
test_c = first_col(test_tots)
test_start, test_end = calculate_start_end_range(test_c, test_occ_matrix, test_query)
assert test_start == 2
assert test_end == 3

test_word = 'bananaban$'
test_query = 'ban'
test_occ_matrix, test_tots = make_occurrences_matrix(bwt_via_bwm(test_word))
test_c = first_col(test_tots)
test_start, test_end = calculate_start_end_range(test_c, test_occ_matrix, test_query)
assert test_start == 5
assert test_end == 6

# Global Alignment

assert scoring_matrix_inplace('C', '_') == -7
assert scoring_matrix_inplace('C', 'G') == -2
assert scoring_matrix_inplace('C', 'A') == -2
assert scoring_matrix_inplace('C', 'T') == -1
assert scoring_matrix_inplace('C', 'C') == 1
test_x = 'TACGTCAGC'
test_y = 'TATGTCATGC'
assert scoring_matrix_inplace(test_x[6], test_y[6]) == 1

expected_distance = numpy.array([[0, -7, -14, -21, -28, -35, -42, -49, -56, -63, -70],
              [-7, 1, -6, -13, -20, -27, -34, -41, -48, -55, -62],
              [-14, -6, 2, -5, -12, -19, -26, -33, -40, -47, -54],
              [-21, -13, -5, 1, -6, -13, -18, -25, -32, -39, -46],
              [-28, -20, -12, -6, 2, -5, -12, -19, -26, -31, -38],
              [-35, -27, -19, -11, -5, 3, -4, -11, -18, -25, -32],
              [-42, -34, -26, -18, -12, -4, 4, -3, -10, -17, -24],
              [-49, -41, -33, -25, -19, -11, -3, 5, -2, -9, -16],
              [-56, -48, -40, -32, -24, -18, -10, -2, 3, -1, -8],
              [-63, -55, -47, -39, -31, -25, -17, -9, -3, 1, 0]])
expected_alignment_score = 0
distances, alignment_score = global_alignment(test_x, test_y, scoring_matrix_inplace)
assert (expected_distance == distances).all()
assert alignment_score == expected_alignment_score
