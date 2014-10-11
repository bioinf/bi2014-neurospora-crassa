
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import re

start_codon = "ATG"

cinput = open('codons.txt', 'r')
all_codons = {}
for line in cinput.readlines():
	a = line.strip().split()
	all_codons[a[0]] = a[1]

parser = argparse.ArgumentParser(description='Checking alignment.')
parser.add_argument('-i', '--input', help="file with alignment")
parser.add_argument('-o', '--output', help="path to the output file")
args = parser.parse_args()

alignments = SeqIO.parse(args.input, "fasta")

output = open(args.output, 'w')

for alignment in alignments:
	ranges = [int(x) for x in re.findall("\d+", re.findall(r'\[.*\]', alignment.description)[0])]
	ranges.sort()
	seq = alignment.seq
	exon_union = ""
	for i in range(len(ranges) // 2):
		exon_union += str(seq[ranges[2*i]:ranges[2*i+1] + 1])
	print()

	if len(exon_union) < 7 or len(exon_union) % 3 != 0:
		continue
	
	if start_codon != str(exon_union[0:3]):
		continue

	good = True
	gen = ""
	for i in range(1, len(exon_union) // 3 - 1):
		codon = all_codons[exon_union[3*i:3*i+3]]
		if codon == 'X':
			good = False
			break
		gen += codon

	if not good and codons[exon_union[-3:]] != 'X':
		continue

	record = SeqRecord(
		Seq(gen),
		id=alignment.id,
		description = alignment.description
	)
	output.write(record.format("fasta"))