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


def draw_single_plot(graph_num, pars, read):
    # x - reads, y - alignment score
    plt.figure()
    plt.suptitle(f'Graph number: {graph_num}')
    plt.xlabel("Reads")
    plt.ylabel("Alignment Score")
    plt.figtext(0.5, .8,
                pars,
                horizontalalignment="center",
                wrap=True, fontsize=10,
                bbox={'facecolor': 'yellow',
                      'alpha': 0.3, 'pad': 5})
    plt.plot(range(len(read)), read)
    plt.show()


def draw_all_plot(reads):
    # x - reads, y - alignment score
    # label - match, mismatch, gap, margin, seed, etc
    plt.figure()
    plt.suptitle("All in one")
    plt.xlabel("Reads")
    plt.ylabel("Alignment Score")
    for r in range(len(reads)):
        plt.plot(range(len(reads[r])), reads[r], label=r)
    plt.legend()
    plt.show()


# dummy test data
reads = [[7, 2, 5, 3, 5, 9], [1, 2, 3, 4, 5, 9], [11, 5, 7, 1, 0, 3]]
draw_single_plot(1, "match: 2; mismatch: -2; gap: -5", reads[0])
draw_all_plot(reads)
