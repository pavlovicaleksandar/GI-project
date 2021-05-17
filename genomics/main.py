from genomics.seed_extend import *
from genomics.burrows import *
from Bio import SeqIO
import logging
import click
import time
import csv
import sys
import ast

logger = logging.getLogger(__name__)


def import_file(path, file_type):
    return list(map(lambda r: str(r.seq), SeqIO.parse(path, file_type)))


def import_fasta_fastq(fasta_path='../data/example_human_reference.fasta', fastq_path='../data/example_human_Illumina.pe_1.fastq'):
    return import_file(fasta_path, file_type='fasta'), import_file(fastq_path, file_type='fastq')


def write_results_to_csv_file(file_name, results):
    with open(file_name, mode='w') as result_file:
        writer = csv.writer(result_file)
        writer.writerow(['start', 'end', 'alignment-score', 'transcription'])
        for result in results:
            writer.writerow(list(result))


def write_data_to_file(data, file_name):
    with open(file_name, 'w') as om:
        # TODO - not happy with this, feels like I am duplicating the data. Discuss
        om.write(str(data))


def read_from_file_to_dict(file_name):
    with open(file_name) as f:
        txt = f.read()
    # TODO - again, I feel like I'm loading the same thing twice, but I know no better. Discuss
    generated_dict = ast.literal_eval(txt)
    logger.info(f'Successfully read {file_name}')
    return generated_dict


def enable_logging(level=logging.INFO):
    logging.basicConfig(filename='example.log', level=logging.INFO)
    _logger = logging.getLogger()
    _logger.setLevel(level=level)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level=level)
    _logger.addHandler(handler)


@click.group()
def cli():
    pass


@cli.command()
@click.option('--fasta_path', type=str, required=True)
def giska(fasta_path):
    try:
        enable_logging()
        logger.info(f'Starting giseed giska for {fasta_path}')

        reference = import_file(path=fasta_path, file_type='fasta')[0]
        logger.info('Finished importing fasta file')

        occurrences_matrix, totals = get_occurrences_matrix_and_totals(reference)
        logger.info(f'Finished generating occurrences matrix and totals for {fasta_path} fasta file')

        write_data_to_file(occurrences_matrix, 'occurrences_matrix.txt')
        logger.info('Finished writing occurrences matrix to txt file')

        c = first_col(totals)
        logger.info("Generated c array from totals")

        write_data_to_file(c, 'c.txt')
        logger.info("Finished writing to c.txt")

        suffix_array = make_suffix_array(reference)
        logger.info(f"Created suffix array from reference {fasta_path}")

        write_data_to_file(suffix_array, 'suffix_array.txt')
        logger.info("Finished writing to suffix_array.txt")

    except Exception as exc:
        logger.error(f'Following exception occurred {exc}')


@cli.command()
@click.option('--fasta_path', type=str, required=True)
@click.option('--fastq_path', type=str, required=True)
@click.option('--occurrences_matrix_path', type=str, required=True)
@click.option('--c_path', type=str, required=True)
@click.option('--suffix_array_path', type=str, required=True)
@click.option('--margin', type=int, default=2)
@click.option('--match', type=int, required=True)
@click.option('--mismatch', type=int, required=True)
@click.option('--gap', type=int, required=True)
@click.option('--seed_length', type=int, default=10)
def ekstendovic(fasta_path, fastq_path, occurrences_matrix_path, c_path,
                suffix_array_path, margin, match, mismatch, gap, seed_length):
    try:
        enable_logging()

        scoring_points = generate_scoring_points(match, mismatch, gap)
        logger.info(f'Generated scoring points for match: {match}, mismatch: {mismatch}, gap: {gap}')

        references, reads = import_fasta_fastq(fasta_path, fastq_path)
        logger.info(f'Imported fasta: {fastq_path} and fastq: {fastq_path}')

        occurrences_matrix = read_from_file_to_dict(occurrences_matrix_path)
        c = read_from_file_to_dict(c_path)
        suffix_array = read_from_file_to_dict(suffix_array_path)

        results = seed_and_extend(references[0], reads, occurrences_matrix, c,
                                  suffix_array, scoring_points, margin, seed_length)
        logger.info('Finished with seed and extend.')

        write_results_to_csv_file(f'results-{time.time()}.csv', results)
        logger.info('Finished writing results to csv file')

    except Exception as exc:
        logger.error(f'Following exception occurred {exc}')


if __name__ == '__main__':
    cli()
