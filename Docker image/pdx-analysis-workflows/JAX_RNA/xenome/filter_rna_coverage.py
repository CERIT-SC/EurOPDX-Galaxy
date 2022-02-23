#! /usr/bin/env python
"""
For determining whether we have enough data to make meaningful RNA expression
estimations, we'll look at the number of human reads available after the
Xenome step.
Version: 1.1.0
"""

from __future__ import print_function

import argparse
import datetime
import os
import subprocess
import sys
import time


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--minimum-reads', default='50000',
                        help="Minimum number of human reads [default: 50000]")
    parser.add_argument('files', nargs='+',
                        help="The file[s] to test.")

    return parser.parse_args()


def process_file(fn, minimum, multiple):
    """
    Determine whether there was adequate coverage of bases in this file.
    NOTE: Returns True if the run was OK.
    :param fn:
    :param minimum:
    :param multiple:
    :return: True if meets criterion, False otherwise
    """
    line_found = False
    summary_section = False
    line = ""
    for line in open(fn):
        line = line.strip()
        if line == 'Summary':
            summary_section = True
            continue
        if summary_section and line.endswith('human'):
            line_found = True
            break
    if not line_found:
        print("Could not find coverage line in:", fn, file=sys.stderr)
        return False

    count = int(line.split()[0])
    if count < minimum:
        if multiple:
            print("{0}\t{1}".format(count, fn))
        else:
            print("Too low human read count: {0} {1}".format(
                count, fn))
        return False
    return True


def main():
    # args = parse_args()
    sample_type = sys.argv[1]
    interim_results_dir = '/galaxy/reference-data/'
    prefix_dir = interim_results_dir + "xenome/"

    delay_init = 180
    delay_normal = 90
    initial = True
    count = 0

    # if os.path.exists(interim_results_dir + xenome_stat_file):
    #     os.remove(interim_results_dir + xenome_stat_file)

    # if not os.path.exists(interim_results_dir):
    #     os.mkdir(interim_results_dir)

    # print("Initial delay is set to: " + str(delay_init))

    if not os.path.exists(prefix_dir):
        print(
            "Seems like Xenome Index Create Tool is not Run yet. Please run Xenome Index Create Tool first (It should "
            "be run only once when you started using Galaxy after download).")
        sys.exit(1)

    print("[INFO] Sample Type is " + sample_type)

    try:
        if sample_type == 'single_end':
            minimum_reads = int(sys.argv[4])
            sample_name = sys.argv[5]
            xenome_stat_file = sample_name + "_xenome_stat"
            xenome_threads = sys.argv[6]

            print("[INFO] Sample name is " + sample_name)
            sp = subprocess.Popen(
                ['xenome', 'classify',
                 '-T', xenome_threads,
                 '-P', prefix_dir + sys.argv[3],
                 # '--output-filename-prefix',  sample_name,
                 '-i', sys.argv[2]], close_fds=True)

        else:
            minimum_reads = int(sys.argv[5])
            xenome_stat_file = sys.argv[6] + "_xenome_stat"
            xenome_threads = sys.argv[7]
            sp = subprocess.Popen(
                ['xenome', 'classify',
                 '-T', xenome_threads,
                 '-P', prefix_dir + sys.argv[4],
                 '--pairs',
                 '--host-name', 'mouse',
                 '--graft-name', 'human',
                 '-i', sys.argv[2], '-i', sys.argv[3]], close_fds=True)
    except Exception as e:
        print('Error executing  xenome classify -> %s' % e)
        sys.exit(1)

    if sample_type == 'single_end':
        human_output = 'graft.fastq'
        mouse_output = 'host.fastq'
        both_output = 'both.fastq'
        neither_output = 'neither.fastq'
        ambiguous_output = 'ambiguous.fastq'
    else:
        human_output = 'human_1.fastq'
        mouse_output = 'mouse_1.fastq'
        both_output = 'both_1.fastq'
        neither_output = 'neither_1.fastq'
        ambiguous_output = 'ambiguous_1.fastq'

    while count == 0:
        if initial:
            time.sleep(delay_init)
            initial = False
        else:
            time.sleep(delay_normal)

        try:
            modified_time = datetime.datetime.fromtimestamp(os.path.getmtime(human_output))
            current_time = datetime.datetime.now()
            back_time = current_time - datetime.timedelta(minutes=5)
            print("[INFO] Current Time: %s Modified Time: %s Back Time: %s" % (current_time, modified_time, back_time))
        except OSError:
            modified_time = 0

        check_is_size_zero_bytes = os.stat(human_output).st_size

        if check_is_size_zero_bytes == 0:
            time.sleep(delay_normal)
        else:
            if modified_time < current_time - datetime.timedelta(minutes=5):
                print("[INFO] Check is matched. Going to kill the Xenome process.")
                sp.terminate()
                count = 1

    h1 = subprocess.Popen(['wc', '-l', human_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    h2 = subprocess.Popen(('cut', '-d', ' ', '-f', '1'), stdin=h1.stdout, stdout=subprocess.PIPE)
    human = int(h2.communicate()[0]) / 4

    m1 = subprocess.Popen(['wc', '-l', mouse_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    m2 = subprocess.Popen(('cut', '-d', ' ', '-f', '1'), stdin=m1.stdout, stdout=subprocess.PIPE)
    mouse = int(m2.communicate()[0]) // 4

    b1 = subprocess.Popen(['wc', '-l', both_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    b2 = subprocess.Popen(('cut', '-d', ' ', '-f', '1'), stdin=b1.stdout, stdout=subprocess.PIPE)
    both = int(b2.communicate()[0]) / 4

    n1 = subprocess.Popen(['wc', '-l', neither_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    n2 = subprocess.Popen(('cut', '-d', ' ', '-f', '1'), stdin=n1.stdout, stdout=subprocess.PIPE)
    neither = int(n2.communicate()[0]) / 4

    a1 = subprocess.Popen(['wc', '-l', ambiguous_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    a2 = subprocess.Popen(('cut', '-d', ' ', '-f', '1'), stdin=a1.stdout, stdout=subprocess.PIPE)
    ambiguous = int(a2.communicate()[0]) / 4

    total = human + mouse + both + neither + ambiguous

    print("[INFO] Total reads = " + str(total) + " Human: " + str(human) + " Mouse: " + str(mouse) + " Both: " + str(
        both) + " Neither: " + str(neither) + " Ambiguous: " + str(ambiguous))

    try:
        human_percentage = ((human * 100) // total)
        mouse_percentage = ((mouse * 100) // total)
        both_percentage = ((both * 100) // total)
        neither_percentage = ((neither * 100) // total)
        ambiguous_percentage = ((ambiguous * 100) // total)
    except Exception as e:
        print('Error finding percentage -> %s' % e)

    try:
        f = open(xenome_stat_file, "a")
        f.write("Summary\n")
        f.write("count\tpercent\tclass\n")
        f.write(str(human) + "\t" + str(human_percentage) + "\thuman\n")
        f.write(str(mouse) + "\t" + str(mouse_percentage) + "\tmouse\n")
        f.write(str(both) + "\t" + str(both_percentage) + "\tboth\n")
        f.write(str(neither) + "\t" + str(neither_percentage) + "\tneither\n")
        f.write(str(ambiguous) + "\t" + str(ambiguous_percentage) + "\tambiguous\n")
        f.close()
    except Exception as e:
        print('Error opening/saving results to the file -> %s' % e)
        f.close()

    multiple = True
    success = True

    if multiple:
        print("[INFO] Human reads\tRun")

    if human < minimum_reads:
        print("[INFO] Too low human read count: ", human)
        sys.exit(1)


if __name__ == '__main__':
    main()
