#!/usr/bin/env python

import argparse
import sweet_utils
import sys
from Bio import SeqIO
import re
from Bio import SearchIO
from Bio import pairwise2
from collections import namedtuple
from Bio.SubsMat import MatrixInfo as matlist
from sweet_utils import dbpath2seq

Exon = namedtuple('Exon', ['query_start','query_end', 'hit_start', 'hit_end'])
ExonerateInfo = namedtuple('ExonerateInfo', ['percentage'])

def printhsp(dna, exons):
    print(exons)
    prev = exons[0].hit_start - 50
    result = ''
    for _, _, s, e in exons:
        result += rna[prev:s].lower() + rna[s:e].upper()
    result += rna[exons[-1].hit_end:exons[-1].hit_end + 50].lower()
    print(result)


def getexons(aln, start_hit, start_query):
    i, j = start_hit, start_query

    isexon = False
    exons = []

    for a, b in zip(*aln):
        if not isexon and a != '-':
            start_hit, start_query = i, j
            isexon = True
        elif isexon and a == '-':
            exons.append(Exon(start_query, j, start_hit, i))
            isexon = False
        if b != '-':
            i += 3
        if a != '-':
            j += 1
        if b == '*':
            break
    if isexon:
        exons.append(Exon(start_query, j, start_hit, i))
    return exons

def retrieve_peptide(dna, exons):
    result = []
    for _, _, s, e in exons:
        result += list(dna[s:e].transcribe().translate())
    if result[-1] == '*':
        result.pop()
    return "".join(result)


def check_hsp(dna, exons):
    first_exon, last_exon, = exons[0], exons[-1]
    has_start_codon = (dna[first_exon.hit_start:first_exon.hit_end+3]
            .transcribe().translate()[0] == 'M')
    has_stop_codon = (dna[last_exon.hit_end - 3:last_exon.hit_end]
            .transcribe().translate()[0] == '*')
    if has_start_codon and has_stop_codon:
        return True
    else:
        return False

parser = argparse.ArgumentParser(
        description='Filtering genes with good structure properties.')
parser.add_argument(dest='found')
parser.add_argument(dest='proteins')
parser.add_argument(dest='dbname')
parser.add_argument('--perc', '-p', type=int, default=0, dest='percentage')

args = parser.parse_args()
genome = SeqIO.to_dict(SeqIO.parse(dbpath2seq(args.dbname), 'fasta'))
queries = SeqIO.to_dict(SeqIO.parse(args.proteins, 'fasta'))

good_entries = []
tool_id = sweet_utils.get_tool_id(args.found)

def exonerate_infos(filepath):
    ADDINFO = re.compile('^<< (\S+) >>$')
    with open(args.found, 'r') as foundfile:
        for line in foundfile:
            m = ADDINFO.match(line)
            if m is not None:
                yield ExonerateInfo(float(m.group(1)))

if tool_id == 'exonerate':
    exonerate_it = iter(exonerate_infos(args.found))

with open(args.found, 'r') as foundfile:
    for found in SearchIO.parse(foundfile, sweet_utils.get_run_format(args.found)):
        query = queries[found.id]
        good_found = SearchIO.QueryResult(id=found.id)
        for hit in found:
            if tool_id == 'exonerate':
                for hsp in hit.hsps:
                    hsp.ident_pct = next(exonerate_it).percentage

            if tool_id == 'blast':
                hit.sort(key=lambda hsp: hsp.ident_num,  reverse=True)
            elif tool_id == 'exonerate':
                hit.sort(key=lambda hsp: hsp.ident_pct, reverse=True)

            best_hsp = hit.hsps[0]
    
            hit_dna = genome[hit.id].seq[best_hsp.hit_start:best_hsp.hit_end]
            if best_hsp[0].hit_strand == -1:
                hit_dna = hit_dna.reverse_complement()

            if tool_id == 'blast':
                exons = getexons(best_hsp[0].aln, 0, best_hsp.query_start)
                percentage = best_hsp.ident_num * 100 / len(query)
            elif tool_id == 'exonerate':
                exons = []
                for fragment in best_hsp:
                    start_idx = best_hsp.hit_start
                    hit_start = fragment.hit_start - best_hsp.hit_start
                    hit_end = fragment.hit_end - best_hsp.hit_start
                    if fragment.hit_strand == -1:
                        hit_start, hit_end = hit_end - 1, hit_start
                        hit_start = len(hit_dna) - 1 - hit_start
                        hit_end = len(hit_dna) - hit_end
                    exons.append(Exon(
                        fragment.query_start, fragment.query_end,
                        hit_start, hit_end))
                percentage = 100
                percentage = best_hsp.ident_pct
    
            if check_hsp(hit_dna, exons) and percentage > args.percentage:
                good_hit = SearchIO.Hit(id=hit.id)
                good_hit.append(best_hsp)
                good_found.append(good_hit)
            
        if len(good_found) != 0:
            good_entries.append(good_found)
print(good_entries)
