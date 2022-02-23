#! /usr/bin/env python
"""
Picard Alignment Metrics.
Version: 0.1.1
"""

from __future__ import print_function

import argparse
import os
import shutil
import sys


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('files', nargs='+',
                        help="The file[s] to use.")

    return parser.parse_args()


def main():
    args = parse_args()

    sample_name = sys.argv[4]
    RGID = "4"
    RGPL = "illumina"
    RGPU = "unit1"
    RGSM = "20"
    picard_stat_file = sample_name + "_picard_stat"
    interim_results_dir = '/mnt/volume/shared/reference-data/'

    # Check whether interim results directory exists or not. If not create it.
    if not os.path.exists(interim_results_dir):
        os.mkdir(interim_results_dir)

    # Remove existing files to avoid files exist exception.
    if os.path.exists(interim_results_dir + picard_stat_file):
        os.remove(interim_results_dir + picard_stat_file)

    # Adding readgroup information
    input_bam = args.files[0]

    try:
        os.system("picard AddOrReplaceReadGroups \
            INPUT=" + input_bam + " \
            OUTPUT=genome_bam_read_group.bam \
            SORT_ORDER=coordinate \
            RGID=" + RGID + " \
            RGLB=" + sample_name + " \
            RGPL=" + RGPL + " \
            RGPU=" + RGPU + " \
            RGSM=" + RGSM + " \
            CREATE_INDEX=true"
                  )
    except Exception as e:
        print('Error executing  picard AddOrReplaceReadGroups -> %s' % e)
        sys.exit(1)

    # Create dictionary for the reference file
    # SEQUENCE_DICTIONARY = os.path.splitext(args.files[3])[0] + ".dict"
    try:
        SEQUENCE_DICTIONARY = args.files[2] + ".dict"
        if not os.path.exists(SEQUENCE_DICTIONARY):
            symlink_to_fa = args.files[2] + ".fa"
            command = "ln -fs " + args.files[2] + " " + symlink_to_fa
            os.system(command)
            print("Command: " + command)
            os.system("picard CreateSequenceDictionary REFERENCE=" + symlink_to_fa)
    except Exception as e:
        print('Error executing  picard CreateSequenceDictionary -> %s' % e)
        sys.exit(1)

    # Picard Reorder Bam file
    try:
        os.system("picard ReorderSam \
            INPUT=genome_bam_read_group.bam \
            OUTPUT=genome_bam_read_group_reorder.bam \
            SEQUENCE_DICTIONARY=" + SEQUENCE_DICTIONARY + " \
            ALLOW_INCOMPLETE_DICT_CONCORDANCE=true \
            CREATE_INDEX=true"
                  )
    except Exception as e:
        print('Error executing  picard ReorderSam -> %s' % e)
        sys.exit(1)

    # Picard SortSam (generating sorted alignment bam file)
    try:
        os.system("picard SortSam \
            SORT_ORDER=coordinate \
            INPUT=genome_bam_read_group_reorder.bam \
            OUTPUT=genome_bam_read_group_reorder_sorted.bam \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=true"
                  )
    except Exception as e:
        print('Error executing  picard SortSam -> %s' % e)
        sys.exit(1)

    # CollectRnaSeqMetrics
    try:
        os.system("picard CollectRnaSeqMetrics \
            INPUT=genome_bam_read_group_reorder_sorted.bam \
            OUTPUT=" + picard_stat_file + " \
            REF_FLAT=" + args.files[1] + " \
            STRAND=NONE \
            VALIDATION_STRINGENCY=SILENT \
            CHART_OUTPUT=File"
                  )
    except Exception as e:
        print('Error executing  picard CollectRnaSeqMetrics -> %s' % e)
        sys.exit(1)

    #shutil.move(picard_stat_file, interim_results_dir)

    # Summary metrics compilation
    # try:
    #     os.system(
    #         "perl ../../../../../tools/pdx-analysis-workflows/RNA_PDX/summary_QC_metrics.pl " + interim_results_dir + quality_stat_file + " " + interim_results_dir + xenome_stat_file + " " + interim_results_dir + rsem_stat_file + " " + interim_results_dir + picard_stat_file + " > " + interim_results_dir + summary_stat_file)
    #     shutil.copy(interim_results_dir + summary_stat_file, os.getcwd())
    # except Exception as e:
    #     print('Error executing summary_QC_metrics.pl -> %s' % e)
    #     sys.exit(1)


if __name__ == '__main__':
    main()
