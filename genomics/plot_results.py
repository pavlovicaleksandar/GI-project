import matplotlib.pyplot as plt
# todo - please check requirements.txt


def import_bwa_mem(file_path):
    # todo - import sam file
    result = []
    return result


def parse_sam_file(file):
    # todo - take as and reads
    result = []
    return result


def import_single_giska_output(file_path):
    # todo - import one file
    result = []
    return result


def import_all_giska_outputs():
    # todo - for each giska output
    result = []
    return result


def draw_plot(title, legend, read_alignment_scores):
    # x - reads, y - alignment score
    # title
    # legend - margin, seed, etv
    plt.figure()
    # plt.legend(legend) todo
    plt.suptitle(title)
    plt.xlabel("Reads")
    plt.ylabel("Alignment Score")
    i = range(len(read_alignment_scores))
    plt.plot(i, read_alignment_scores)
    plt.show()


draw_plot("proba", "margin: 2", [1, 2, 3, 4, 5, 9])
