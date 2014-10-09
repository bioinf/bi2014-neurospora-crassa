#!/usr/bin/env python

import os
import sys
import tabulate
import argparse
from collections import namedtuple, defaultdict
from Bio import SearchIO

AlignmentData = namedtuple('AlignmentData', ['runid', 'score', 'range', 'exons'])

def get_run_id(run_path):
    return os.path.basename(run_path)[:-len('.parseable.out')]


def get_tool_id(run_path):
    return get_run_id(run_path).split('.')[-1]


def parse_n_fill_run_data(run_path, run_data):
    run_id = get_run_id(run_path)
    searchio_format = {'blast': 'blast-xml', 'exonerate': 'exonerate-vulgar'}[get_tool_id(run_path)]

    for query in SearchIO.parse(run_path, searchio_format):
        for hit in query.hits:
            for hsp in hit.hsps:
                exons = [x.hit_range for x in hsp.fragments]

                if get_tool_id(run_path) == 'blast':
                    score = hsp.bitscore
                elif get_tool_id(run_path) == 'exonerate':
                    score = hsp.score
                else:
                    score = 'N/A'

                alignment = AlignmentData(run_id, score, hsp.hit_range, exons)
                run_data[query.id][hit.id].append(alignment)


def print_run_data(run_data, output):
    for query_id, hits in run_data.items():
        for hit_id, alignments in hits.items():
            print('matching {} --> {}'.format(query_id, hit_id), file=output)
            data = map(list, sorted(alignments, key=lambda x: x.range[0]))
            print(tabulate.tabulate(data, headers=AlignmentData._fields), file=output)
            print('\n', file=output)

parser = argparse.ArgumentParser(description='Sweet comparator')
parser.add_argument('results', nargs='+', help='Results to compare.')
args = parser.parse_args()

all_run_data = defaultdict(lambda: defaultdict(list))
for path in args.results:
    parse_n_fill_run_data(path, all_run_data)
print_run_data(all_run_data, sys.stdout)


