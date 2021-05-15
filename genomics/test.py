from genomics.main import *

assert rotations('') == []
assert rotations('$') == ['$']
assert rotations('banana$') == ['banana$', 'anana$b', 'nana$ba', 'ana$ban', 'na$bana','a$banan', '$banana']

assert bwm('') == []
assert bwm('$') == ['$']
assert bwm('banana$') == ['$banana', 'a$banan', 'ana$ban', 'anana$b', 'banana$', 'na$bana', 'nana$ba']

assert bwt_via_bwm('banana$') == 'annb$aa'
occ_matrix, tots = make_occurrences_matrix(bwt_via_bwm('banana$'))
assert occ_matrix == {'a': [1, 1, 1, 1, 1, 2, 3], 'n': [0, 1, 2, 2, 2, 2, 2], 'b': [0, 0, 0, 1, 1, 1, 1], '$': [0, 0, 0, 0, 1, 1, 1]}
assert tots == {'a': 3, 'n': 2, 'b': 1, '$': 1}
