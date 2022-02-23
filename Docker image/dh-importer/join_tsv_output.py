#! /usr/bin/env python3
"""
For generating single output table to export to DataHub.
Version: 1.2.0
"""

import argparse as ap
import logging
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
    parent_parser.add_argument('--input_file', nargs='+', required=True)
    parent_parser.add_argument('--output_file', nargs=1, required=True)

    parser = ap.ArgumentParser(prog='join_tsv_output', usage='%(prog)s [options]')
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


def merge_expression_data(df_expression, count):
    logger.info(f'Process Expression Data')
    column_name = df_expression.columns[2]
    logger.info(f'Sample_ID: {column_name}')
    df_expression = df_expression.drop('gene_id', 1)
    df_expression['sample_id'] = column_name
    df_expression.rename(columns={column_name: "rnaseq_count", "GeneName": "symbol"}, inplace=True)

    # Remove the header names from data frame except first data frame
    if count == 0:
        df_expression.to_csv(arguments.output_file[0], sep='\t', mode='a', index=False)
    else:
        df_expression.to_csv(arguments.output_file[0], sep='\t', mode='a', header=False, index=False)


def merge_mutation_data(df_mutation, count):
    logger.info(f'Process Mutation Data')
    # Remove the header names from data frame except first data frame
    if count == 0:
        df_mutation.to_csv(arguments.output_file[0], sep='\t', mode='a', index=False)
    else:
        df_mutation.to_csv(arguments.output_file[0], sep='\t', mode='a', header=False, index=False)


arguments = parse_arguments()
count = 0

for dataset_name in arguments.input_file:
    logger.info(f'Processing Dataset: {dataset_name}')
    try:
        df = pandas.read_csv(dataset_name, sep='\t', lineterminator='\n')
        df_initial = df.copy()

        # Process expression data
        if arguments.command == "expression":
            merge_expression_data(df, count)

        # Process mutation data
        if arguments.command == "mutation":
            merge_mutation_data(df, count)

        count = 1
    except Exception as exception:
        logger.error(f'Processing data in join_tsv_output tool -> {exception}')
        sys.exit(1)

logger.info(f'Process completed.')
