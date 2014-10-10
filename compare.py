#!/usr/bin/env python

import os
import sys
import tabulate
import argparse
import genblastA_to_gff3
from collections import namedtuple, defaultdict
from Bio import SearchIO
from Bio import SeqIO

AlignmentData = namedtuple('AlignmentData', ['runid', 'score', 'matched', 'coverage', 'range', 'exons'])


def get_run_id(run_path):
    return os.path.basename(run_path)[:-len('.parseable.out')]


def get_tool_id(run_path):
    return get_run_id(run_path).split('.')[-1]

searchio_formats = {'blast': 'blast-xml', 'exonerate': 'exonerate-vulgar'}

def parse_n_fill_run_data_searchio(run_path, run_data, querydb):
    run_id = get_run_id(run_path)
    run_format = {'blast': 'blast-xml', 'exonerate': 'exonerate-vulgar'}[get_tool_id(run_path)]
    for query in SearchIO.parse(run_path, run_format):
        for hit in query.hits:
            for hsp in hit.hsps:
                exons = [x.hit_range for x in hsp.fragments]
                coverage = 'N/A'

                if querydb is not None:
                    total_matched = sum(x.query_span for x in hsp.fragments)
                    coverage = '{:.2f}%'.format(100 * total_matched / len(querydb[query.id]))

                if hasattr(hsp, 'score'):
                    score = hsp.score
                elif hasattr(hsp, 'bitscore'):
                    score = hsp.bitscore
                else:
                    score = 'N/A'

                if hasattr(hsp, 'ident_num') and hasattr(query, 'seq_len'):
                    matched = '{:.2f}%'.format(100 * hsp.ident_num / query.seq_len)
                else:
                    matched = 'N/A'

                alignment = AlignmentData(run_id, score, matched, coverage, hsp.hit_range, exons)
                run_data[query.id][hit.id].append(alignment)


def parse_n_fill_run_data_genblasta(run_path, run_data):
    run_id = get_run_id(run_path)
    with open(run_path, 'r') as runfile:
        for hit in genblastA_to_gff3.parse_genblastA(runfile):
                for hsp in hit['hsps'].values():
                    match = hit['match']
                    exons = [(hsp['match_start'] - 1, hsp['match_end'] - 1)]
                    matched = '{:.2f}%'.format(hsp['perc_id'])
                    alignment = AlignmentData(run_id, match['score'], matched, match['coverage_perc'] + '%',
                                              (int(match['match_start']), int(match['match_end'])), exons)
                    run_data[match['query_name']][match['match_name']].append(alignment)


def parse_n_fill_run_data(run_path, run_data, querydb):
    tool_id = get_tool_id(run_path)
    if tool_id in searchio_formats.keys():
        parse_n_fill_run_data_searchio(run_path, run_data, querydb)
    elif tool_id == 'genblasta':
        parse_n_fill_run_data_genblasta(run_path, run_data)


def print_run_data(run_data, output):
    for query_id, hits in run_data.items():
        for hit_id, alignments in hits.items():
            print('matching {} --> {}'.format(query_id, hit_id), file=output)
            data = map(list, sorted(alignments, key=lambda x: x.range[0]))
            print(tabulate.tabulate(data, headers=AlignmentData._fields), file=output)
            print('\n', file=output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sweet comparator')
    parser.add_argument('results', nargs='+', help='Results to compare.')
    parser.add_argument('-q', '--query', dest='query')
    args = parser.parse_args()

    all_run_data = defaultdict(lambda: defaultdict(list))

    querydb = None
    if args.query is not None:
        # should we use index or _index_db_?
        querydb = SeqIO.to_dict(SeqIO.parse(args.query, 'fasta'))
    for path in args.results:
        parse_n_fill_run_data(path, all_run_data, querydb)
    print_run_data(all_run_data, sys.stdout)