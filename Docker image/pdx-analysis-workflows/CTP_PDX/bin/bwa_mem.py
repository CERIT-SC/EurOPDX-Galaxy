#! /usr/bin/env python
from __future__ import print_function

"""
BWA-Mem Alignment to reference genome.
Version: 1.0.0
"""
import sys
import os


def averageLineLength(forward, reverse):
    fn1 = open(forward)
    fn2 = open(reverse)
    lines1 = fn1.readlines()
    lines2 = fn2.readlines()
    avg_forward = sum([len(line.strip('\n')) for line in lines1[1::4]]) / (len(lines1)/4)
    avg_reverse = sum([len(line.strip('\n')) for line in lines2[1::4]]) / (len(lines2)/4)
    return (avg_forward + avg_reverse)/2

def main():
    bwa_index = sys.argv[1]
    forward_reads = sys.argv[2]
    reverse_reads = sys.argv[3]
    types = sys.argv[4]
    try:
        file = open("read_group", "r")
        read_group = "'" + str(file.read()).strip() + "'"
        print(read_group)
        file.close()
    except Exception as e:
        print('Error while executing opening file %s' % e)
        sys.exit(1)
    try:
        #command = "seqtk mergepe " + forward_reads + " " + reverse_reads + " > out.fq"
        #print(command)
        #os.system(command)
        if types == "paired" :
            command = "bwa mem -t8 -R " + read_group + " " + bwa_index + " " + forward_reads + " " + reverse_reads + " > out.sam"
            print(command)
            os.system(command)
        else :
            command = "bwa mem -t8 -R " + read_group + " " + bwa_index + " " + forward_reads + " > out.sam"
            print(command)
            os.system(command)
        command = "samtools view -1 -S -b out.sam > aln.bam"
        print(command)
        os.system(command)
    except Exception as e:
        print('Error while executing bwamem -> %s %s' % command % e)
        sys.exit(1)

if __name__ == '__main__':
    main()
    sys.exit(0)
