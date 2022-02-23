#! /usr/bin/env python

import os
from sys import argv
from sys import stderr
import pandas as pd


def main():
    table=argv[1]
    filtername = argv[2]
    filtertext = open(filtername).read()
    df = pd.read_csv(table, sep="\t")
    df = df[df['GeneName'].isin(filtertext.split("\n"))]
    #df['sample_id'] = df['sample_id'].apply(lambda x : x.split('_')[0])
    f = open("filtered.csv", "w")
    df.to_csv(f , sep="\t", columns=['GeneName', 'sample_id', 'Quartile'])
    f.close()

if __name__ == '__main__':
    main()

