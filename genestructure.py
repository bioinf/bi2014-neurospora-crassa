#!/usr/bin/env python

import argparse
import sweet_utils
import sys
from Bio import SeqIO
from Bio import SearchIO
from Bio import pairwise2
from collections import namedtuple
from Bio.SubsMat import MatrixInfo as matlist
from sweet_utils import dbpath2seq

Exon = namedtuple('Exon', ['query_start','query_end', 'hit_start', 'hit_end'])

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
    has_start_codon =(dna[first_exon.hit_start:first_exon.hit_end+3]
            .transcribe().translate()[0] == 'M')
    has_stop_codon = (dna[last_exon.hit_end - 3:last_exon.hit_end]
            .transcribe().translate()[0] == '*')
    if has_start_codon and has_stop_codon:
        return True
    else:
        return False

parser = argparse.ArgumentParser(
        description='Filtering genes with good structure properties.')
parser.add_argument(dest='found', type=argparse.FileType('r'))
parser.add_argument(dest='proteins')
parser.add_argument(dest='dbname')
parser.add_argument('--perc', '-p', type=int, default=0, dest='percentage')

args = parser.parse_args()
genome = SeqIO.to_dict(SeqIO.parse(dbpath2seq(args.dbname), 'fasta'))
queries = SeqIO.to_dict(SeqIO.parse(args.proteins, 'fasta'))

good_entries = []
for found in SearchIO.parse(args.found, 'blast-xml'):
    query = queries[found.id]
    good_found = SearchIO.QueryResult(id=found.id)
    for hit in found:
        hit.sort(key=lambda hsp: hsp.ident_num, reverse=True)
        best_hsp = hit.hsps[0]

        hit_dna = genome[hit.id].seq[best_hsp.hit_start:best_hsp.hit_end]
        if best_hsp.hit_strand == -1:
            hit_dna = hit_dna.reverse_complement()

        exons = getexons(best_hsp[0].aln, 0, best_hsp.query_start)
        percentage = best_hsp.ident_num * 100 / len(query)

        if check_hsp(hit_dna, exons) and percentage > args.percentage:
            good_hit = SearchIO.Hit(id=hit.id)
            good_hit.append(best_hsp)
            good_found.append(good_hit)
        
    if len(good_found) != 0:
        good_entries.append(good_found)
SearchIO.write(good_entries, sys.stdout, "blast-tab")
