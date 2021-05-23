import matplotlib.pyplot as plt
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def import_single_giska_output(match, mismatch, gap, margin, seed_length):
    version = f'match_{match}_mismatch_{mismatch}_gap_{gap}_margin_{margin}_seed-length_{seed_length}'
    file_name = f'results_{version}.csv'
    with open(file_name) as f:
        txt = pd.read_csv(f, skiprows=1)
    logger.info(f'Successfully read {file_name}')
    return txt["alignment-score"]


def import_draw_all_giska_outputs():
    margin = 2
    seed_length = 10
    all_reads = []
    i = 0
    for match in [0, 1, 2]:
        for mismatch in [-3, -2]:
            for gap in [-7, -5]:
                single_read = import_single_giska_output(match=match, mismatch=mismatch, gap=gap, margin=margin, seed_length=seed_length)
                i += 1
                draw_single_plot(i, f'match: {match}, mismatch: {mismatch}, gap: {gap}', single_read)
                all_reads.append(single_read)
    draw_all_plot(all_reads)
    return all_reads


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
# reads = [[7, 2, 5, 3, 5, 9], [1, 2, 3, 4, 5, 9], [11, 5, 7, 1, 0, 3]]
# draw_single_plot(1, "match: 2; mismatch: -2; gap: -5", reads[0])
# draw_all_plot(reads)
read = import_single_giska_output(match=1, mismatch=-2, gap=-7, margin=2, seed_length=10)
draw_single_plot(1, f'match: {1}, mismatch: {-2}, gap: {-7}', read)
