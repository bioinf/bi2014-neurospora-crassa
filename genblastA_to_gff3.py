#!/usr/bin/env python

import argparse
import sys
import re
import logging
import logging.config
import json
import os

START_STR = '//*****************START'
END_STR = '//******************END'
QUERY_NAME_STR = '//for query: '
HSP_STR = 'HSP_ID'
QUERY_NAME_RE = re.compile('^//for query: ([^/]*)//$')
GENOMIC_MATCH_RE = re.compile(
    '^(?P<query_name>[^|]*)\|(?P<match_name>[^:]*):(?P<match_start>\d+)\.\.(?P<match_end>\d+)\|(?P<strand>[+-])\|gene cover:(?P<coverage_num>\d+)\((?P<coverage_perc>[\d.]+)%\)\|score:(?P<score>[-\de.]+)\|rank:(?P<rank>\d+)$')
HSP_RE = re.compile(
    '^HSP_ID\[(?P<hsp_id>\d+)\]:\((?P<match_start>\d+)-(?P<match_end>\d+)\);query:\((?P<query_start>\d+)-(?P<query_end>\d+)\); pid: (?P<perc_id>[\d.]+)$')


def parse_genblastA(input_filename):
    in_record = False
    matches_by_query = dict()

    def dict_from_match_re(genomic_match):
        match_dict = genomic_match.groupdict()
        match_index = matches_by_query.get(match_dict['query_name'], 0) + 1
        match_dict['index'] = match_index
        matches_by_query[match_dict['query_name']] = match_index
        return match_dict

    # variables set during parsing
    hsp_dict = dict()
    genomic_match = None
    query_name = ''
    for line in input_filename:
        if not in_record:
            if line.startswith(START_STR):
                in_record = True
        else:
            # we're in a record
            if line.startswith(END_STR):
                # check that we've got a match to output (we might have NONE)
                if genomic_match:
                    match_dict = dict_from_match_re(genomic_match)
                    yield (dict(match=match_dict, hsps=hsp_dict))
                hsp_dict = dict()
                genomic_match = None
                query_name = ''
                in_record = False
            elif line.startswith(QUERY_NAME_STR):
                match = QUERY_NAME_RE.match(line.rstrip())
                if match == None:
                    logging.error('Query name regexp failed to match on line: {}'.format(line))
                    in_record = False
                else:
                    query_name = match.group(1)
            elif 'gene cover' in line:
                if genomic_match:
                    # we've already seen one match, need to output that
                    match_dict = dict_from_match_re(genomic_match)
                    yield (dict(match=match_dict, hsps=hsp_dict))
                    # and reset the hsp_dict, we're about to reset the genomic_match
                    hsp_dict = dict()
                genomic_match = GENOMIC_MATCH_RE.match(line.rstrip())
                if not genomic_match:
                    logging.error('Genomic match regexp failed to match on line: {}'.format(line))
                    in_record = False
                else:
                    # not much to do: we've got a genomic match saved in genomic_match, will use it once we've read the HSPs
                    logging.debug(
                        'Got match between {} and {} start: {} end: {}'.format(genomic_match.group('query_name'),
                                                                               genomic_match.group('match_name'),
                                                                               genomic_match.group('match_start'),
                                                                               genomic_match.group('match_end')))
            elif line.startswith(HSP_STR):
                match = HSP_RE.match(line.rstrip())
                if not match:
                    logging.error('HSP regexp failed to match on line: {}'.format(line))
                    in_record = False
                else:
                    # save the HSPs for this genomic match
                    hsp = dict(match_start=int(match.group('match_start')), match_end=int(match.group('match_end')),
                               query_start=int(match.group('query_start')), query_end=int(match.group('query_end')),
                               perc_id=float(match.group('perc_id')))
                    hsp_dict[int(match.group('hsp_id'))] = hsp


def write_gff_line(genomic_match, hsp_dict, query_name, output_file):
    # gff3 format
    # seq source type start end score strand phase attributes"
    first_hsp_id = min(hsp_dict.keys())
    last_hsp_id = max(hsp_dict.keys())
    target = 'Target={} {} {}'.format(query_name, hsp_dict[first_hsp_id]['query_start'],
                                      hsp_dict[last_hsp_id]['query_end'])
    attributes = 'ID={}_{};{}'.format(query_name, genomic_match['index'], target)
    gff_line = '\t'.join([genomic_match['match_name'], 'genBlastA', 'match',
                          genomic_match['match_start'], genomic_match['match_end'],
                          genomic_match['coverage_perc'], genomic_match['strand'],
                          '.', attributes]) + '\n'
    output_file.write(gff_line)


def write_bed_line(genomic_match, hsp_dict, query_name, output_file):
    # bed format
    # chrom chromStart chromEnd name score strand [other optional fields, see http://genome.ucsc.edu/FAQ/FAQformat.html#format1]
    # bed score is 0 to 1000 - here the coverage percentage is scaled to that range
    name = '{}_{}'.format(query_name, genomic_match['index'])
    bed_line = '\t'.join(
        [genomic_match['match_name'], int(genomic_match['match_start'] - 1), int(genomic_match['match_end'] - 1),
         name, str(float(genomic_match['coverage_perc']) * 10), genomic_match['strand']]) + '\n'
    output_file.write(bed_line)


def genblastA_process(input_file, output_file, output_format='gff3', min_perc_coverage=0.0, min_match_length=0,
                      min_perc_identity=0):
    if output_format == 'gff3':
        output_file.write('##gff-version 3\n')
    for match in parse_genblastA(input_file):
        genomic_match = match['match']
        hsp_dict = match['hsps']
        num_hsps = len(hsp_dict)
        match_length = abs(int(genomic_match['match_end']) - int(genomic_match['match_start']))
        avg_perc_identity = sum([hsp_dict[i]['perc_id'] for i in hsp_dict]) / num_hsps
        query_coverage_perc = float(genomic_match['coverage_perc'])
        if (avg_perc_identity >= min_perc_identity and
                    query_coverage_perc >= min_perc_coverage and
                    match_length >= min_match_length):
            if output_format == 'gff3':
                write_gff_line(match['match'], match['hsps'], match['match']['query_name'], output_file)
            elif output_format == 'bed':
                write_bed_line(match['match'], match['hsps'], match['match']['query_name'], output_file)
            else:
                sys.stderr.write('Unknown output format: {}\n'.format(args.output_format))
                sys.exit(1)
    input_file.close()
    output_file.close()


if __name__ == '__main__':
    log_config = os.getenv('LOG_CONFIG', None)
    if log_config:
        log_config_file = None
        try:
            log_config_file = open(log_config)
        except IOError as e:
            sys.stderr.write('Failed to load logging config from {}: {}\n'.format(log_config, str(e)))
        if log_config_file:
            config_dict = json.load(log_config_file)
            try:
                logging.config.dictConfig(config_dict)
            except (ValueError, TypeError, AttributeError, ImportError) as e:
                sys.stderr.write('Failed to parse log config dictionary: {}\n'.format(str(e)))
                logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description='parse genblastA output and produce GFF3')
    parser.add_argument('--min_perc_coverage', '-C', type=float, default=80.0,
                        help='Minimum coverage of the query sequence')
    parser.add_argument('--min_match_length', '-L', type=int, default=100, help='Shortest match length to accept')
    parser.add_argument('--min_perc_identity', '-I', type=float, default=80.0,
                        help='Minimum average % identity to accept')
    parser.add_argument('--output_format', '-F', type=str, choices=['gff3', 'bed'], default='gff3',
                        help='Output format: GFF3 or BED')
    parser.add_argument('genblastA_file', type=argparse.FileType(), help='genblastA format file')
    parser.add_argument('output_file', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='GFF3 output file')
    args = parser.parse_args()

    genblastA_process(args.genblastA_file, args.output_file, args.output_format, args.min_perc_coverage,
                      args.min_match_length, args.min_perc_identity)
