from genomics.burrows import *
from genomics.global_alignment import *
from genomics.seed_extend import *
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import numpy

# CMD Line Args


def define_arg_parser():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--fasta', type=str, required=True, help="Path to the fasta file")
    arg_parser.add_argument('--fastq', type=str, required=True, help="Path to the fastq file")
    arg_parser.add_argument('--margin', default=2)
    arg_parser.add_argument('--match', type=int, nargs="+", default=[0, 1, 2])
    arg_parser.add_argument('--mismatch', type=int, nargs="+", default=[-3, -2])
    # arg_parser.add_argument('--Ti', type=int, nargs="+", default=-1)
    # arg_parser.add_argument('--Tv', type=int, nargs="+", default=-3)
    arg_parser.add_argument('--gap', type=int, nargs="+", default=[-7, -5])
    arg_parser.add_argument('--seed-length', type=int, default=10)
    return arg_parser


def import_or_generate_scoring_points():
    # todo : create this in define_arg_parser()
    parser = argparse.ArgumentParser(description='Process scoring points.')
    parser.add_argument('integers', type=int, nargs='*',
                        help='4 scoring points', default=[1, -1, -2, -7])
    args = parser.parse_args()
    scoring_points_input = sorted(args.integers)
    return {
        'M': scoring_points_input[3],
        'Ti': scoring_points_input[2],
        'Tv': scoring_points_input[1],
        'G': scoring_points_input[0]
    }


# Read Fasta Fastq


def import_file(path, file_type):
    return list(map(lambda r: str(r.seq), SeqIO.parse(path, file_type)))


def import_fasta_fastq(fasta_path='../data/example_human_reference.fasta', fastq_path='../data/example_human_Illumina.pe_1.fastq'):
    return import_file(fasta_path, file_type='fasta'), import_file(fastq_path, file_type='fastq')





# args = define_arg_parser().parse_args()
scoring_points = import_or_generate_scoring_points()
word = 'banana$'
query = 'ana'
last_column = bwt_via_bwm(word)
occ_matrix, tots = make_occurrences_matrix(last_column)
c = first_col(tots)
start, end = calculate_start_end_range(c, occ_matrix, query)
bwt_via_sa(word)
suff_arr = make_suffix_array(word)
query_pos = find_all_query_positions_in_word_via_suffix_arr(start, end, suff_arr)

x = 'TACGTCAGC'
y = 'TATGTCATGC'
distances, alignment_score = global_alignment(x, y, scoring_points)
alignment, transcript = traceback(x, y, distances, scoring_points)

references, reads = import_fasta_fastq()

print(seed_and_extend(scoring_points, references, reads))
