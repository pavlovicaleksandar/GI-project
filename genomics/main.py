from genomics.seed_extend import *
from Bio import SeqIO
import click


def import_or_generate_scoring_points(match, missmatch, gap):
    return {
        'M': match,
        'Ti': missmatch,
        'Tv': missmatch,
        'G': gap
    }


def import_file(path, file_type):
    return list(map(lambda r: str(r.seq), SeqIO.parse(path, file_type)))


def import_fasta_fastq(fasta_path='../data/example_human_reference.fasta', fastq_path='../data/example_human_Illumina.pe_1.fastq'):
    return import_file(fasta_path, file_type='fasta'), import_file(fastq_path, file_type='fastq')


@click.group()
def cli():
    pass


@cli.command()
@click.option('--fasta_path', type=str, required=True)
@click.option('--fastq_path', type=str, required=True)
@click.option('--margin', type=int, default=2)
@click.option('--match', type=int, required=True)
@click.option('--missmatch', type=int, required=True)
@click.option('--gap', type=int, required=True)
@click.option('--seed_length', type=int, default=10)
def run(fasta_path, fastq_path, margin, match, missmatch, gap, seed_length):
    scoring_points = import_or_generate_scoring_points(match, missmatch, gap)
    references, reads = import_fasta_fastq(fasta_path, fastq_path)
    seed_and_extend(scoring_points, references, reads, margin, seed_length)


if __name__ == '__main__':
    cli()
