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


scoring_points = {'M': 1, 'Ti': -1, 'Tv': -3, 'G': -7}
x = 'TACGTCAGC'
y = 'TATGTCATGC'
expected_distance = numpy.array([[0, -7, -14, -21, -28, -35, -42, -49, -56, -63, -70],
                                 [-7, 1, -6, -13, -20, -27, -34, -41, -48, -55, -62],
                                 [-14, -6, 2, -5, -12, -19, -26, -33, -40, -47, -54],
                                 [-21, -13, -5, 1, -6, -13, -18, -25, -32, -39, -46],
                                 [-28, -20, -12, -6, 2, -5, -12, -19, -26, -31, -38],
                                 [-35, -27, -19, -11, -5, 3, -4, -11, -18, -25, -32],
                                 [-42, -34, -26, -18, -12, -4, 4, -3, -10, -17, -24],
                                 [-49, -41, -33, -25, -19, -11, -3, 5, -2, -9, -16],
                                 [-56, -48, -40, -32, -24, -18, -10,  -2,   2,  -1,  -8],
                                 [-63, -55, -47, -39, -31, -25, -17,  -9,  -3,  -1,   0]])


@pytest.mark.parametrize("example_this, example_that, result", [
    ('C', '_', -7),
    ('C', 'G', -3),
    ('C', 'A', -3),
    ('C', 'T', -1),
    ('C', 'C', 1),

])
def test_scoring_matrix_inplace(example_this, example_that, result):
    assert scoring_matrix_inplace(example_this, example_that, scoring_points) == result


def test_global_alignment():
    expected_alignment_score = 0
    distances, alignment_score = global_alignment(x, y, scoring_points)
    assert (expected_distance == distances).all()
    assert expected_alignment_score == alignment_score


def test_traceback():
    alignment, transcript = traceback(x, y, expected_distance, scoring_points)
    expected_alignment =\
        '''TACGTCA_GC
|| |||| ||
TATGTCATGC'''  # really unfortunate formatting, consider removing since i think that this check is pointless
    assert "MMRMMMMIMM" == transcript
    assert expected_alignment == alignment


def test_fasta_fastq_import():
    test_references, test_reads = import_fasta_fastq()
    assert 1 == len(test_references)
    assert 172 == len(test_reads)


@pytest.mark.parametrize("example_this, example_that, result", [
    ('ACTGTAAC', 'GTTACAGT', True),
    ('AA', 'TT', True),
    ('', '', True),
    ('', 'A', False),
    ('A', 'C', False),
    ('A', 'G', False),
    ('ACTGTAAC', 'TGACATTG', False),

])
def test_reverse_and_complement(example_this, example_that, result):
    assert (reverse_and_complement(example_this) == example_that) == result
