import pysam
from matplotlib import pyplot as plt

# constants
reversed_complement_flag = 0x10
segment_unmapped_flag = 0x04


def complement_direction_analysis(sam_file_path):
    sam_file = pysam.AlignmentFile(sam_file_path, 'r')

    num_reversed = 0
    num_normal = 0
    total_number_of_reads = 0
    for read in sam_file:
        from_reversed = read.flag & reversed_complement_flag > 0
        unmapped = read.flag & segment_unmapped_flag > 0

        if (from_reversed and not unmapped):
            num_reversed += 1
        if not from_reversed and not unmapped:
            num_normal += 1

        total_number_of_reads += 1

    reversed_percentage = (num_reversed / total_number_of_reads) * 100
    normal_percentage = (num_normal / total_number_of_reads) * 100

    return (normal_percentage, reversed_percentage)


def plot_alignments_sorted_by_alignment_score(sam_file_path):
    sam_file = pysam.AlignmentFile(sam_file_path, 'r')
    alignment_scores = map(lambda item: item.tags[2][1], sam_file)
    alignment_scores = sorted(list(alignment_scores), reverse=True)
    plt.plot(alignment_scores)
    plt.show()


sam_file = "_1_reads.sam"
print(complement_direction_analysis(sam_file))

plot_alignments_sorted_by_alignment_score(sam_file)