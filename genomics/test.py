from genomics.main import *

assert rotations('') == []
assert rotations('$') == ['$']
assert rotations('banana$') == ['banana$', 'anana$b', 'nana$ba', 'ana$ban', 'na$bana','a$banan', '$banana']

assert bwm('') == []
assert bwm('$') == ['$']
assert bwm('banana$') == ['$banana', 'a$banan', 'ana$ban', 'anana$b', 'banana$', 'na$bana', 'nana$ba']

assert bwt_via_bwm('banana$') == 'annb$aa'
test_occ_matrix, test_tots = make_occurrences_matrix(bwt_via_bwm('banana$'))
assert test_occ_matrix == {'a': [1, 1, 1, 1, 1, 2, 3], 'n': [0, 1, 2, 2, 2, 2, 2], 'b': [0, 0, 0, 1, 1, 1, 1], '$': [0, 0, 0, 0, 1, 1, 1]}
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
