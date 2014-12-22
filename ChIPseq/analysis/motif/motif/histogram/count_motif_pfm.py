# usage: count_motif.py motif.pwm sequence.fa
import sys
import numpy as np
from Bio import SeqIO
from Bio import Motif
from operator import itemgetter

def search_motif(mf,seq):
	"""search pwm for each motif in the motiflist form sequence"""
	result = [(score,pos) for pos,score in mf.search_pwm(seq, threshold=5.0)]
	if not result:
		return None
	sort_pos = sorted(result, key = itemgetter(0), reverse=True)
	return sort_pos[0][1]

mf = Motif.read(open(sys.argv[1]),'jaspar-pfm')
for record in SeqIO.parse(sys.argv[2],'fasta'):
	hit = search_motif(mf,record.seq) 
	if hit == None:
		continue
	else:
		record.id = record.id+'_'+str(hit)
		chr = record.id.split(':')[0]
		start = record.id.split(':')[1].split('-')[0]
		if hit >= 0:
			pos = int(start) + hit 
		else:
			pos = int(start) + len(record.seq) + hit
		print chr+'\t'+str(pos)+'\t'+str(pos+1)
