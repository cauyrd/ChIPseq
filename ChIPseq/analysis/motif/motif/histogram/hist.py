# usage: hist.py motiflist.txt sequence.fa control.fa
import sys,os
import numpy as np
from Bio import SeqIO
from Bio import Motif
from operator import itemgetter
import matplotlib
matplotlib.use('Agg')
import pylab

def search_motif(mf,seq):
	"""search pwm for each motif in the motiflist form sequence"""
	result = [(score,pos) for pos,score in mf.search_pwm(seq, threshold=5.0)]
	if not result:
		return None
	sort_pos = sorted(result, key = itemgetter(0), reverse=True)
	mid_pos = len(seq)/2
	pos = sort_pos[0][1]
	if pos >= 0:
		dist = pos - mid_pos
	else:
		dist = pos + mid_pos
	return dist

def draw_plot(motiffile):
	"""generating histogram"""
	count = []
	control = []
	mf = Motif.read(open(motiffile),'jaspar-pfm')
	for record in SeqIO.parse(sys.argv[2],'fasta'):
		hit = search_motif(mf,record.seq) 
		if hit == None:
			continue
		else:
			count.append(hit)
	for record in SeqIO.parse(sys.argv[3],'fasta'):
		hit = search_motif(mf,record.seq) 
		if hit == None:
			continue
		else:
			control.append(hit)
	# assume the sequence length is 201, center base +/- 100bp
	pylab.figure()
	pylab.hist(count, np.linspace(-100,100,101),color='g')
	num, bin = np.histogram(control, np.linspace(-100,100,101))
	pylab.plot(np.linspace(-100,100,100), num, color='r')
	pylab.xlabel('Distance relative to Stat5 motif')
	pylab.ylabel('No. Stat5 peaks')
	motifname = os.path.basename(motiffile)
	pylab.title(motifname.split('.')[0])
	pylab.savefig(motifname.split('.')[0]+'.png')

def main():
	loc = '/home/msistaff/yang4414/project/NGStoolbox/chipseq/motif/TRANSFACT/'
	ifp = open('motiflist.txt')
	for line in ifp:
		motiffile = loc+line.rstrip()
		draw_plot(motiffile)

if __name__ == '__main__':
	main()
