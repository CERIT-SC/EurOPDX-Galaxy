#! /usr/bin/env python3
"""
For generating filtered tsv output for each sample.
Version: 1.1.0
"""

import argparse as ap
import logging
import os
import sys

import pandas

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def parse_arguments():
    """Read and parse commandline arguments"""
    parent_parser = ap.ArgumentParser(add_help=False)
    parent_parser.add_argument('--sample_id', nargs=1, required=True)
    parent_parser.add_argument('--input_file', nargs=1, required=True)
    parent_parser.add_argument('--output_file', nargs=1, required=True)

    parser = ap.ArgumentParser(prog='create_tsv_output', usage='%(prog)s [options]')
    subparsers = parser.add_subparsers(help='Choose a command', dest="command")

    expression_parser = subparsers.add_parser('expression',
                                              parents=[parent_parser],
                                              help='"Expression data"')
    mutation_parser = subparsers.add_parser('mutation',
                                            parents=[parent_parser],
                                            help='"Mutation data"')

    expression_parser.set_defaults(action=lambda: 'expression')
    mutation_parser.set_defaults(action=lambda: 'mutation')

    return parser.parse_args(sys.argv[1:])


def process_expression_data():
    logger.info(f'Process Expression Data')
    metric_column = 'expected_count'
    try:
        df = pandas.read_csv(arguments.input_file[0], sep='\t', lineterminator='\n')
        df = df[['GeneName', 'gene_id', metric_column]]
        df.rename(columns={metric_column: arguments.sample_id[0]}, inplace=True)
        print(df.head())
        columns_list = ["GeneName", "gene_id", arguments.sample_id[0]]
        df.to_csv(arguments.output_file[0], sep='\t', index=False, columns=columns_list)
    except Exception as exception:
        logger.error(f'Error while processing data in create_tsv_output tool -> {exception}')
        sys.exit(1)


def process_mutation_data():
    logger.info(f'Process Mutation Data')
    try:
        command = "sed 's/\\t/,/g' " + str(arguments.input_file[0]) + " > vcf_to_tsv.tsv"
        os.system(command)
        df = pandas.read_csv('vcf_to_tsv.tsv', lineterminator='\n')
        df.insert(0, 'sample_id', str(arguments.sample_id[0]))
        if list(df.columns)[0] == 'EFF[*].GENE':
            #df['sample_id'] = str(arguments.sample_id[0])
            df.rename(columns={'EFF[*].GENE': 'symbol', 'EFF[*].BIOTYPE': 'biotype', 'EFF[*].CODON': 'coding_sequence_change', 'SVTYPE': 'variant_class', 'EFF[*].AA': 'amino_acid_change','EFF[*].EFFECT': 'consequence', 'DP_HQ': 'read_depth', 'ALT_AF': 'allele_frequency', 'CHROM': 'chromosome', 'POS': 'seq_start_position', 'REF': 'ref_allele', 'ALT': 'alt_allele', 'ANN[*].GENEID': 'ensembl_gene_id', 'EFF[*].TRID': 'ensembl_transcript_id', 'ID': 'variation_id'}, inplace=True)
        else:
            df.rename(columns={'ANN[*].GENE': 'symbol', 'ANN[*].BIOTYPE': 'biotype', 'ANN[*].HGVS_C': 'coding_sequence_change', 'TYPE': 'variant_class', 'ANN[*].HGVS_P': 'amino_acid_change','ANN[*].EFFECT': 'consequence', 'DP': 'read_depth', 'AF': 'allele_frequency', 'CHROM': 'chromosome', 'POS': 'seq_start_position', 'REF': 'ref_allele', 'ALT': 'alt_allele', 'ANN[*].GENEID': 'ensembl_gene_id', 'ANN[*].FEATUREID': 'ensembl_transcript_id', 'ID': 'variation_id'}, inplace=True)
        print(df.head())
        df.to_csv(arguments.output_file[0], sep='\t', index=False)
    except Exception as exception:
        logger.error(f'Error while processing data in create_tsv_output tool -> {exception}')
        sys.exit(1)


arguments = parse_arguments()

# Process expression data
if arguments.command == "expression":
    process_expression_data()

# Process mutation data
if arguments.command == "mutation":
    process_mutation_data()

logger.info(f'Process completed.')
