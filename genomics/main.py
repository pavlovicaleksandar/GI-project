from genomics.seed_extend import *
from Bio import SeqIO
import click
import time
import csv


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
    scoring_points = import_or_generate_scoring_points(match, mismatch, gap)
    references, reads = import_fasta_fastq(fasta_path, fastq_path)
    results = seed_and_extend(scoring_points, references, reads, margin, seed_length)
    write_results_to_csv_file(f'results-{time.time()}.csv', results)


if __name__ == '__main__':
    cli()
