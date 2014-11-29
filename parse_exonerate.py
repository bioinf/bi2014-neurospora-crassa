from Bio import SearchIO                       	
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse	

def get_substring(sequence, s, e, orientation):
	if orientation == 1:
		return sequence[s:e]
	else:
		return sequence[s:e].reverse_complement()	
	

parser = argparse.ArgumentParser(description='Parsing exonerate.')
parser.add_argument('-g', '--genome', help="path to fasta file with genome")
parser.add_argument('-e', '--exonerate', help="path to the result of exonerate")
parser.add_argument('-o', '--output', help="path to the output file")
args = parser.parse_args()
qresult = SearchIO.parse(args.exonerate, 'exonerate-vulgar')

sequences = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

output = open(args.output, 'w')

for query in qresult:
	best_score = 0
	hit_range = (0, 0)
	exons = []
	hit_id = ""
	orientation = 1				
	for hit in query.hits:
		for hsp in hit.hsps:
			if best_score < hsp.score:
				best_score = hsp.score
				exons = [x.hit_range for x in hsp.fragments]
				hit_id = hit.id
				hit_range = hsp.hit_range
				orientation =  hsp.hit_strand_all[0]
	substring = get_substring(sequences[hit_id].seq, hit_range[0], hit_range[1], orientation)
	print(query.id)
	if orientation == 1:	
		for i in range(len(exons)):
		        exons[i] = (exons[i][0] - hit_range[0], exons[i][1] - hit_range[0])
	if orientation == -1:
		for i in range(len(exons)):
			exons[i] = (hit_range[1] - exons[i][1], hit_range[1] - exons[i][0])
										 			
	record = SeqRecord(             		
		get_substring(sequences[hit_id].seq, hit_range[0], hit_range[1], orientation),
		id=query.id + "" + hit_id,
		description = str(best_score) + " " + str(hit_range) + " " + str(exons)
	)
	output.write(record.format("fasta"))
	
