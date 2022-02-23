#!/usr/bin/env python3
"""
For preparing and uploading data to DataHub.
Version: 1.1.0
"""

import argparse as ap
import json
import logging
import os
import subprocess
import sys
import tarfile
import time

import pandas as pd
import requests as rq

# Disable warnings
rq.packages.urllib3.disable_warnings()

# Show all columns
pd.set_option('display.max_columns', None)

# create logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

VERIFY = True
BASE_DIR = "europdx"
DATA_DIR = BASE_DIR + "/data/UPDOG/PIPELINE"
MAPPING_DIR = BASE_DIR + "/mapping"
COMPRESSED_DATA_FILE = BASE_DIR + ".tgz"
METADATA_FILES_PREFIX = "PIPELINE_metadata-"
SAMPLEPLATFORM_FILES_PREFIX = "PIPELINE_sampleplatform-data"
FILES_SUFFIX = ".tsv"
METADATA_XSLX_WORKSHEETS = {"checklist",
                            "patient",
                            "sample",
                            "model",
                            "model_validation",
                            "sharing",
                            "loader"}
hdr = {}


def parse_arguments():
    """Read and parse commandline arguments"""
    parent_parser = ap.ArgumentParser(add_help=False)
    parent_parser.add_argument('--endpoint', nargs=1, required=True)
    parent_parser.add_argument('--insecure', action='store_const', const=True, default=True)
    parent_parser.add_argument('--password', nargs=1, required=False)
    parent_parser.add_argument('--reset', action='store_const', const=True, default=False)

    parser = ap.ArgumentParser(prog='dh-importer', usage='%(prog)s [options]')
    subparsers = parser.add_subparsers(help='Choose a command', dest="command")

    expression_parser = subparsers.add_parser('expression',
                                              parents=[parent_parser],
                                              help='"Expression data"')
    mutation_parser = subparsers.add_parser('mutation',
                                            parents=[parent_parser],
                                            help='"Mutation data"')
    featurecounts_parser = subparsers.add_parser('featurecounts',
                                              parents=[parent_parser],
                                              help='"Expression data"')

    expression_parser.add_argument('--metadata', nargs=1, required=True)
    expression_parser.add_argument('--sampleplatform', nargs=1, required=True)
    expression_parser.add_argument('--mappings', nargs=1, required=True)
    expression_parser.add_argument('--genome_assembly', nargs=1, required=True)
    expression_parser.add_argument('--datafile', nargs=1, required=True)
    expression_parser.add_argument('--outputfile', nargs=1, required=True)
    expression_parser.set_defaults(action=lambda: 'expression')

    featurecounts_parser.add_argument('--metadata', nargs=1, required=True)
    featurecounts_parser.add_argument('--sampleplatform', nargs=1, required=True)
    featurecounts_parser.add_argument('--mappings', nargs=1, required=True)
    featurecounts_parser.add_argument('--genome_assembly', nargs=1, required=True)
    featurecounts_parser.add_argument('--datafile', nargs=1, required=True)
    featurecounts_parser.add_argument('--outputfile', nargs=1, required=True)
    featurecounts_parser.set_defaults(action=lambda: 'expression')
    
    mutation_parser.add_argument('--metadata', nargs=1, required=True)
    mutation_parser.add_argument('--sampleplatform', nargs=1, required=True)
    mutation_parser.add_argument('--mappings', nargs=1, required=True)
    mutation_parser.add_argument('--genome_assembly', nargs=1, required=True)
    mutation_parser.add_argument('--datafile', nargs=1, required=True)
    mutation_parser.add_argument('--outputfile', nargs=1, required=True)
    mutation_parser.set_defaults(action=lambda: 'mutation')

    return parser.parse_args(sys.argv[1:])


def parse_status(resp):
    """Return SandBox status"""
    stat = json.loads(resp)
    return stat['status']


def wait_for_up(api_caller):
    """Check SandBox status"""
    global VERIFY
    retries = 100
    retry = retries
    api_response = sess.get(srv + '/api/status', verify=VERIFY, headers=hdr)
    status = parse_status(api_response.text)

    while retry > 0 and api_response.status_code == 200 and status != 'running':
        retry -= 1
        logger.info(f'Sand Box GET /status: {status}')
        time.sleep(15)
        api_response = sess.get(srv + '/api/status', verify=VERIFY, headers=hdr)
        status = parse_status(api_response.text)

    if api_response.status_code != 200:
        raise ConnectionAbortedError(api_response.text)

    if retry == 0:
        if api_caller == 'reset':
            api_response = sess.get(srv + '/api/log/reset', verify=VERIFY, headers=hdr)
            logger.error(api_response.text)
        if api_caller == 'import':
            api_response = sess.get(srv + '/api/log/import', verify=VERIFY, headers=hdr)
            logger.error(api_response.text)
        raise TimeoutError

    logger.info(f'Sand Box GET /status: {status}')


def process_expression_data():
    """Process expression data and augment with provided metadata"""
    expression_path = DATA_DIR + "/expression"

    create_directories(expression_path)

    # Load sampleplatform.xslx
    logger.info("Load sampleplatform.xslx data in Pandas")
    df_sampleplatform_origin = pd.read_excel(arguments.sampleplatform[0],
                                             skiprows=[1, 2, 3, 4],
                                             dtype='unicode')

    # Filter expression data
    df_sampleplatform_filtered = df_sampleplatform_origin[
        df_sampleplatform_origin["molecular_characterisation_type"] == "expression"]

    # Check if expression data exists
    if len(df_sampleplatform_filtered.index) < 1:
        logger.error("No expression data found in Sampleplatform.xslx.")
        sys.exit(1)

    # Select required columns
    df_sampleplatform_selected = df_sampleplatform_filtered[["model_id",
                                                             "sample_id",
                                                             "sample_origin",
                                                             "host_strain_nomenclature",
                                                             "passage",
                                                             "platform"]]
    # print(df_sampleplatform_selected.head(), end='\n' * 3)

    # Load pipeline output file
    logger.info("Load pipeline data in Pandas")
    df_pipeline_result = pd.read_csv(arguments.datafile[0], sep='\t', lineterminator='\n')

    crosscheck_sampleplatform_and_pipeline_data(df_pipeline_result, df_sampleplatform_selected)
    df_sampleplatform_selected["sample_id"] = df_sampleplatform_selected["sample_id"].astype(str)
    df_pipeline_result["sample_id"] = df_pipeline_result["sample_id"].astype(str)
    # Join results of two datasets
    df_result = pd.merge(df_pipeline_result, df_sampleplatform_selected,
                         left_on=['sample_id'],
                         right_on=['sample_id'],
                         how='left',
                         suffixes=('_left', '_right'))

    # Add missing columns
    ddf_result = df_result.assign(chromosome='',
                                 strand='',
                                 seq_start_position='',
                                 seq_end_position='',
                                 ucsc_gene_id='',
                                 ensembl_gene_id='',
                                 ensembl_transcript_id='',
                                 rnaseq_coverage='',
                                 rnaseq_fpkm='',
                                 rnaseq_tpm='',
                                 affy_hgea_probe_id='',
                                 affy_hgea_expression_value='',
                                 illumina_hgea_probe_id='',
                                 illumina_hgea_expression_value='',
                                 z_score='',
                                 genome_assembly=arguments.genome_assembly[0])

    # Reorder columns
    reordered_columns_list = ['model_id',
                              'sample_id',
                              'sample_origin',
                              'host_strain_nomenclature',
                              'passage',
                              'chromosome',
                              'strand',
                              'seq_start_position',
                              'seq_end_position', 'Length',
                              'symbol',
                              'ucsc_gene_id',
                              'ensembl_gene_id',
                              'ensembl_transcript_id',
                              'rnaseq_coverage',
                              'rnaseq_fpkm',
                              'rnaseq_tpm',
                              'rnaseq_count',
                              'affy_hgea_probe_id',
                              'affy_hgea_expression_value',
                              'illumina_hgea_probe_id',
                              'illumina_hgea_expression_value',
                              'z_score',
                              'genome_assembly',
                              'platform']
    
    df_result = df_result[reordered_columns_list]

    # Save result to csv
    df_result.to_csv(arguments.outputfile[0], sep='\t', index=False)
    os.link(arguments.outputfile[0], expression_path + "/" + arguments.outputfile[0])

    extract_excel_sheets()

    # Rename sampleplatform.xslx to sampleplatform.tsv
    # add extra tabs, and rename diagnosis_mappings.json
    df_sampleplatform_origin.to_csv(DATA_DIR + "/" + SAMPLEPLATFORM_FILES_PREFIX + FILES_SUFFIX,
                                    sep='\t',
                                    index=False)
    subprocess.call(["sed", "-i", "-e", 's/$/\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t/',
                     DATA_DIR + "/" + SAMPLEPLATFORM_FILES_PREFIX + FILES_SUFFIX])
    os.link(arguments.mappings[0], MAPPING_DIR + "/diagnosis_mappings.json")


def process_expression_data_featurecounts():
    """Process expression data and augment with provided metadata"""
    expression_path = DATA_DIR + "/expression"

    create_directories(expression_path)

    # Load sampleplatform.xslx
    logger.info("Load sampleplatform.xslx data in Pandas")
    df_sampleplatform_origin = pd.read_excel(arguments.sampleplatform[0],
                                             skiprows=[1, 2, 3, 4],
                                             dtype='unicode')

    # Filter expression data
    df_sampleplatform_filtered = df_sampleplatform_origin[
        df_sampleplatform_origin["molecular_characterisation_type"] == "expression"]

    # Check if expression data exists
    if len(df_sampleplatform_filtered.index) < 1:
        logger.error("No expression data found in Sampleplatform.xslx.")
        sys.exit(1)

    # Select required columns
    df_sampleplatform_selected = df_sampleplatform_filtered[["model_id",
                                                             "sample_id",
                                                             "sample_origin",
                                                             "host_strain_nomenclature",
                                                             "passage",
                                                             "platform"]]
    # print(df_sampleplatform_selected.head(), end='\n' * 3)

    # Load pipeline output file
    logger.info("Load pipeline data in Pandas")
    df_pipeline_result = pd.read_csv(arguments.datafile[0], sep='\t', lineterminator='\n')
    logger.info(df_pipeline_result.columns[-1])
    normalization_name = df_pipeline_result.columns[-1]

    crosscheck_sampleplatform_and_pipeline_data(df_pipeline_result, df_sampleplatform_selected)
    df_sampleplatform_selected["sample_id"] = df_sampleplatform_selected["sample_id"].astype(str)
    df_pipeline_result["sample_id"] = df_pipeline_result["sample_id"].astype(str)
    # Join results of two datasets
    df_result = pd.merge(df_pipeline_result, df_sampleplatform_selected,
                         left_on=['sample_id'],
                         right_on=['sample_id'],
                         how='left',
                         suffixes=('_left', '_right'))

    # Add missing columns
    df_result = df_result.assign(ucsc_gene_id='',
                                 ensembl_gene_id='',
                                 ensembl_transcript_id='',
                                 rnaseq_coverage='',
                                 affy_hgea_probe_id='',
                                 affy_hgea_expression_value='',
                                 illumina_hgea_probe_id='',
                                 illumina_hgea_expression_value='',
                                 z_score='',
                                 genome_assembly=arguments.genome_assembly[0])

    if normalization_name == 'rnaseq_tpm':
        df_result = df_result.assign(rnaseq_fpkm='')
    elif normalization_name == 'rnaseq_fpkm':
        df_result = df_result.assign(rnaseq_tpm='')
    else:
        df_result = df_result.assign(rnaseq_tpm='', 
                                    rnaseq_fpkm='')

    # Reorder columns
    if 'rpkm' not in df_result.columns:
        logger.info("one method")
        reordered_columns_list = ['model_id',
                              'sample_id',
                              'sample_origin',
                              'host_strain_nomenclature',
                              'passage',
                              'chromosome',
                              'strand',
                              'seq_start_position',
                              'seq_end_position',
                              'symbol',
                              'ucsc_gene_id',
                              'ensembl_gene_id',
                              'ensembl_transcript_id',
                              'rnaseq_coverage',
                              'rnaseq_fpkm',
                              'rnaseq_tpm',
                              'rnaseq_count',
                              'affy_hgea_probe_id',
                              'affy_hgea_expression_value',
                              'illumina_hgea_probe_id',
                              'illumina_hgea_expression_value',
                              'z_score',
                              'genome_assembly',
                              'platform']
    else:
        logger.info("all methods")
        reordered_columns_list = ['model_id',
                              'sample_id',
                              'sample_origin',
                              'host_strain_nomenclature',
                              'passage',
                              'chromosome',
                              'strand',
                              'seq_start_position',
                              'seq_end_position',
                              'symbol',
                              'ucsc_gene_id',
                              'ensembl_gene_id',
                              'ensembl_transcript_id',
                              'rnaseq_coverage',
                              'median_of_ratios',
                              'tmm', 'cpm', 'tpm','rpkm',
                              'affy_hgea_probe_id',
                              'affy_hgea_expression_value',
                              'illumina_hgea_probe_id',
                              'illumina_hgea_expression_value',
                              'z_score',
                              'genome_assembly',
                              'platform']
    df_result = df_result[reordered_columns_list]

    # Save result to csv
    df_result.to_csv(arguments.outputfile[0], sep='\t', index=False)
    os.link(arguments.outputfile[0], expression_path + "/" + arguments.outputfile[0])

    extract_excel_sheets()

    # Rename sampleplatform.xslx to sampleplatform.tsv
    # add extra tabs, and rename diagnosis_mappings.json
    df_sampleplatform_origin.to_csv(DATA_DIR + "/" + SAMPLEPLATFORM_FILES_PREFIX + FILES_SUFFIX,
                                    sep='\t',
                                    index=False)
    subprocess.call(["sed", "-i", "-e", 's/$/\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t/',
                     DATA_DIR + "/" + SAMPLEPLATFORM_FILES_PREFIX + FILES_SUFFIX])
    os.link(arguments.mappings[0], MAPPING_DIR + "/diagnosis_mappings.json")


def process_mutation_data():
    """Process mutation data and augment with provided metadata"""
    mutation_path = DATA_DIR + "/mut"

    create_directories(mutation_path)

    # Load sampleplatform.xslx
    logger.info("Load sampleplatform.xslx data in Pandas")
    df_sampleplatform_origin = pd.read_excel(arguments.sampleplatform[0],
                                             skiprows=[1, 2, 3, 4],
                                             dtype='unicode')
    print(df_sampleplatform_origin.head(), end='\n' * 3)
    # Filter mutation data
    df_sampleplatform_filtered = df_sampleplatform_origin[
        df_sampleplatform_origin["molecular_characterisation_type"] == "mutation"]

    # Check if mutation data exists
    logger.info(f'Number of Samples in Sampleplatform.xslx {len(df_sampleplatform_filtered.index)}')
    if len(df_sampleplatform_filtered.index) < 1:
        logger.error("No mutation data found in Sampleplatform.xslx.")
        sys.exit(1)

    # Select required columns
    df_sampleplatform_selected = df_sampleplatform_filtered[["model_id",
                                                             "sample_id",
                                                             "sample_origin",
                                                             "host_strain_nomenclature",
                                                             "passage",
                                                             "platform"]]

    # Load pipeline output file
    logger.info("Load pipeline data in Pandas")
    df_pipeline_result = pd.read_csv(arguments.datafile[0], sep='\t', lineterminator='\n')

    #to avoid errors with samples named with numbers
    df_sampleplatform_selected["sample_id"] = df_sampleplatform_selected["sample_id"].astype(str)
    df_pipeline_result["sample_id"] = df_pipeline_result["sample_id"].astype(str)

    crosscheck_sampleplatform_and_pipeline_data(df_pipeline_result, df_sampleplatform_selected)

    # Join results of two datasets
    df_result = pd.merge(df_pipeline_result, df_sampleplatform_selected,
                         left_on=['sample_id'],
                         right_on=['sample_id'],
                         how='left',
                         suffixes=('_left', '_right'))

    # Add missing columns
    if 'passage' not in df_result.columns:
        df_result = df_result.assign(passage='0.0')
    if 'biotype' not in df_result.columns:
        df_result = df_result.assign(biotype='')
    if 'coding_sequence_change' not in df_result.columns:
        df_result = df_result.assign(coding_sequence_change='')
    if 'variant_class' not in df_result.columns: 
        df_result = df_result.assign(variant_class='')
    if 'codon_change' not in df_result.columns:
        df_result = df_result.assign(codon_change='')
    if 'amino_acid_change' not in df_result.columns:
        df_result = df_result.assign(amino_acid_change='')
    if 'functional_prediction' not in df_result.columns:
        df_result = df_result.assign(functional_prediction='')
    if 'read_depth' not in df_result.columns:
        df_result = df_result.assign(read_depth='')
    if 'ucsc_gene_id' not in df_result.columns:
        df_result = df_result.assign(ucsc_gene_id='')
    if 'ncbi_gene_id' not in df_result.columns:
        df_result = df_result.assign(ncbi_gene_id='')
    if 'ncbi_transcript_id' not in df_result.columns:
        df_result = df_result.assign(ncbi_transcript_id='')
    if 'ensembl_gene_id' not in df_result.columns:
        df_result = df_result.assign(ensembl_gene_id='')
    if 'genome_assembly' not in df_result.columns:
        df_result = df_result.assign(genome_assembly=arguments.genome_assembly[0])
    #we do this just for now cause the parser cant handel too many symbols
    df_result = df_result.assign(variation_id='')
    
    # Reorder columns
    reordered_columns_list = ['model_id',
                              'sample_id',
                              'sample_origin',
                              'host_strain_nomenclature',
                              'passage',
                              'symbol',
                              'biotype',
                              'coding_sequence_change',
                              'variant_class',
                              'codon_change',
                              'amino_acid_change',
                              'consequence',
                              'functional_prediction',
                              'read_depth',
                              'allele_frequency',
                              'chromosome',
                              'seq_start_position',
                              'ref_allele',
                              'alt_allele',
                              'ucsc_gene_id',
                              'ncbi_gene_id',
                              'ncbi_transcript_id',
                              'ensembl_gene_id',
                              'ensembl_transcript_id',
                              'variation_id',
                              'genome_assembly',
                              'platform']
    df_result = df_result[reordered_columns_list]

    # Save result to csv
    df_result.to_csv(arguments.outputfile[0], sep='\t', index=False)
    os.link(arguments.outputfile[0], mutation_path + "/" + arguments.outputfile[0])

    extract_excel_sheets()

    # Rename sampleplatform.xslx to sampleplatform.tsv
    # add extra tabs, and rename diagnosis_mappings.json
    df_sampleplatform_origin.to_csv(DATA_DIR + "/" + SAMPLEPLATFORM_FILES_PREFIX + FILES_SUFFIX,
                                    sep='\t',
                                    index=False)
    subprocess.call(["sed", "-i", "-e", 's/$/\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t/',
                     DATA_DIR + "/" + SAMPLEPLATFORM_FILES_PREFIX + FILES_SUFFIX])
    os.link(arguments.mappings[0], MAPPING_DIR + "/diagnosis_mappings.json")


def create_directories(directory_path):
    """Create output and mapping directories for upload to SandBox"""
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    if not os.path.exists(MAPPING_DIR):
        os.makedirs(MAPPING_DIR)


def crosscheck_sampleplatform_and_pipeline_data(df_pipeline, df_sampleplatform_selected):
    """Check whether all the samples [sampleplatform.xslx and Pipeline Output] exist or not."""
    df_sampleplatform_unique = df_sampleplatform_selected["sample_id"].drop_duplicates().to_frame()
    df_pipeline_unique = df_pipeline["sample_id"].drop_duplicates().to_frame()
    df_sampleplatform_unique["sample_id"] = df_sampleplatform_unique["sample_id"].astype(str)
    df_pipeline_unique["sample_id"] = df_pipeline_unique["sample_id"].astype(str)
    df_all = df_sampleplatform_unique.merge(df_pipeline_unique,
                                            on=['sample_id'],
                                            how='left',
                                            indicator=True)
    df_left_only = df_all[df_all['_merge'] == 'left_only']
    unmatched_samples = len(df_left_only.index)
    if unmatched_samples > 0:
        logger.warning(f'sampleplatform.xslx has unmapped samples in pipeline: {unmatched_samples}')
        print(df_left_only, end='\n' * 3)
    df_all = df_pipeline_unique.merge(df_sampleplatform_unique,
                                      on=['sample_id'],
                                      how='left', indicator=True)
    df_left_only = df_all[df_all['_merge'] == 'left_only']
    unmatched_samples = len(df_left_only.index) - 1

    # 1 is subtracted because of column name in data
    if unmatched_samples > 0:
        logger.error(f'Pipeline output has unmapped samples in sampleplatform.xslx: {unmatched_samples}')
        print(df_left_only, end='\n' * 3)
        sys.exit(1)


def extract_excel_sheets():
    """Extract individual worksheets from metadata.xslx and create a separate file for each worksheet."""
    logger.info("Load metadata.xslx in Pandas")
    df_metadata_origin = pd.ExcelFile(arguments.metadata[0])
    for worksheet in METADATA_XSLX_WORKSHEETS:
        logger.info(f'Converted {worksheet} worksheet to tsv')
        temp_df = pd.read_excel(df_metadata_origin, worksheet)
        temp_df.to_csv(DATA_DIR + "/" + METADATA_FILES_PREFIX + worksheet + FILES_SUFFIX,
                       sep='\t',
                       index=False)


def create_tarfile():
    """Create compressed file for upload"""
    with tarfile.open(COMPRESSED_DATA_FILE, "w:gz") as tar:
        tar.add(BASE_DIR)


def upload_and_import():
    """Upload and import data (SandBox)"""
    api_response = sess.post(srv + '/api/upload',
                             verify=VERIFY,
                             headers=hdr,
                             files={'file': open(COMPRESSED_DATA_FILE, 'rb')})
    if api_response.status_code != 200 and api_response.status_code != 201:
        logger.error(f'POST /upload: {api_response.status_code}')
        raise ConnectionAbortedError(api_response.text)

    api_response = sess.post(srv + '/api/import', verify=VERIFY, headers=hdr)
    if api_response.status_code != 200:
        raise ConnectionAbortedError(api_response.text)

    wait_for_up(api_caller='import')


# Parse command line arguments
arguments = parse_arguments()

# Create Session
srv = arguments.endpoint[0].rstrip('/')
if arguments.password[0]:
    hdr = {"Authorization": "Bearer " + arguments.password[0]}
sess = rq.Session()
sess.verify = not arguments.insecure
logger.info(f'Disable TLS {arguments.insecure}')
logger.info(f'Reset Sand Box {arguments.reset}')

# Reset Sand Box
if arguments.reset:
    logger.info(f'POST /reset {hdr}')
    response = sess.post(srv + '/api/reset', verify=VERIFY, headers=hdr)
    if response.status_code != 200:
        raise ConnectionAbortedError(response.text)

    # Wait for software stack to start
    wait_for_up(api_caller='reset')

# Process expression data
if arguments.command == "expression":
    process_expression_data()

# Process expression data with featurecounts format
if arguments.command == "featurecounts":
    process_expression_data_featurecounts()

# Process mutation data
if arguments.command == "mutation":
    process_mutation_data()

create_tarfile()
upload_and_import()

logger.info('All the steps completed.')
