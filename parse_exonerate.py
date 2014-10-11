from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser(description='Parsing exonerate.')
parser.add_argument('-g', '--genome', help="path to fasta file with genome")
parser.add_argument('-e', '--exonerate', help="path to the result of exonerate")
parser.add_argument('-o', '--output', help="path to the output file")
args = parser.parse_args()
qresult = SearchIO.parse(args.exonerate, 'exonerate-vulgar')

genome = SeqIO.read(args.genome, "fasta").seq

output = open(args.output, 'w')

for query in qresult:
	minimum_score = float("inf")
	hit_range = (0, 0)
	exons = []
	hit_id = ""
	for hit in query.hits:
		for hsp in hit.hsps:
			if minimum_score > hsp.score:
				minimum_score = hsp.score
				exons = [x.hit_range for x in hsp.fragments]
				hit_id = hit.id
				hit_range = hsp.hit_range
	record = SeqRecord(
		Seq(genome[hit_range[0]:hit_range[1]]),
		id=query.id + "" + hit_id,
		description = str(minimum_score) + " " + str(hit_range) + " " + str(exons)
	)
	output.write(record.format("fasta"))
	
