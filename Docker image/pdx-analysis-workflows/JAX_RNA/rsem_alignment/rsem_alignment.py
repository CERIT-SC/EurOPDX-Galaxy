#! /usr/bin/env python
from __future__ import print_function

"""
RSEM Alignment to transcriptome.
Version: 1.3.0
"""
import sys
import os
import shutil
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed-length', '--seed-length', default='25',
                        help="Seed length used by the read aligner. Providing the correct value is important for RSEM. [default: 25]")
    parser.add_argument('--forward-prob', '--forward-prob', default='0.5',
                        help="Probability of generating a read from the forward strand of a transcript. (Default: 0.5)")
    parser.add_argument('files', nargs='+',
                        help="The file[s] to use.")
    # parser.add_argument('-reference', '--directory', dest='odir', default='.',
    #                        help=
    #                        'Reference base file path with reference name e.g. /home/user/hg_38 '
    #                        '[current directory]')

    return parser.parse_args()


def main():
    # args = parse_args()
    sample_type = sys.argv[1]
    interim_results_dir = '/galaxy/reference-data/rsem/'
    print("[INFO] Sample Type is " + sample_type)

    try:
        if sample_type == 'single_end':
            #sample_name = sys.argv[6]
            sample_name = "sample"
            rsem_threads = sys.argv[7]
            rsem_stat = sys.argv[8]
            command = "rsem-calculate-expression -p " + rsem_threads + \
                      " --phred33-quals --seed-length " + sys.argv[2] + \
                      " --forward-prob " + sys.argv[3] + \
                      " --sort-bam-memory-per-thread 2G " \
                      "--time " \
                      "--output-genome-bam " \
                      "--sort-bam-by-coordinate " \
                      "--bowtie2 " + \
                      sys.argv[4] + " " + \
                      interim_results_dir + sys.argv[5] + " " + \
                      sample_name
        else:
            #sample_name = sys.argv[7]
            sample_name = "sample"
            rsem_threads = sys.argv[8]
            rsem_stat = sys.argv[9]
            command = "rsem-calculate-expression -p " + rsem_threads + \
                      " --phred33-quals --seed-length " + sys.argv[2] + \
                      " --forward-prob " + sys.argv[3] + \
                      " --sort-bam-memory-per-thread 2G " \
                      "--time " \
                      "--output-genome-bam " \
                      "--sort-bam-by-coordinate " \
                      "--bowtie2 " \
                      "--paired-end " + \
                      sys.argv[4] + " " + sys.argv[5] + " " + \
                      interim_results_dir + sys.argv[6] + " " + \
                      sample_name

        print('[INFO] Command: ' + command)
        os.system(command)
    except Exception as e:
        print('Error while executing rsem-calculate-expression -> %s %s' % command % e)
        sys.exit(1)

    try:
        shutil.copy(sample_name + ".stat/" + sample_name + ".cnt", rsem_stat)
    except Exception as e:
        print('Error saving the data in galaxy -> %s' % e)


if __name__ == '__main__':
    main()
    sys.exit(0)
