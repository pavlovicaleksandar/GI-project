from genomics.seed_extend import *
from Bio import SeqIO
import logging
import click
import time
import csv
import sys

logger = logging.getLogger(__name__)


def import_or_generate_scoring_points(match, mismatch, gap):
    return {
        'M': match,
        'Ti': mismatch,
        'Tv': mismatch,
        'G': gap
    }


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


def enable_logging(level=logging.INFO):
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
@click.option('--fastq_path', type=str, required=True)
@click.option('--margin', type=int, default=2)
@click.option('--match', type=int, required=True)
@click.option('--mismatch', type=int, required=True)
@click.option('--gap', type=int, required=True)
@click.option('--seed_length', type=int, default=10)
def run(fasta_path, fastq_path, margin, match, mismatch, gap, seed_length):
    try:
        enable_logging()
        logger.info(f'Starting giseed with params fasta_path={fastq_path} '
                    f'fastq_path={fastq_path} margin={margin} match={match} '
                    f'mismatch={mismatch} gap={gap} seed_length={seed_length}')

        scoring_points = import_or_generate_scoring_points(match, mismatch, gap)

        logger.info('Importing fasta and fastq files')
        references, reads = import_fasta_fastq(fasta_path, fastq_path)
        logger.info('Finished importing fasta and fastq files')

        logger.info('Starting seed and extend')
        results = seed_and_extend(scoring_points, references, reads, margin, seed_length)
        logger.info('Finished seed and extend')

        logger.info('Writing results to csv file')
        # write_results_to_csv_file(f'results-{time.time()}.csv', results)
        logger.info('Finished writing results to csv file')
    except Exception as exc:
        logger.error('Following exception occurred ', exc)


if __name__ == '__main__':
    cli()
