from genomics.main import *

assert rotations('') == []
assert rotations('$') == ['$']
assert rotations('banana$') == ['banana$', 'anana$b', 'nana$ba', 'ana$ban', 'na$bana','a$banan', '$banana']

assert bwm('') == []
assert bwm('$') == ['$']
assert bwm('banana$') == ['$banana', 'a$banan', 'ana$ban', 'anana$b', 'banana$', 'na$bana', 'nana$ba']

assert bwt_via_bwm('banana$') == 'annb$aa'
