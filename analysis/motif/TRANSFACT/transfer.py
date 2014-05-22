# usage: count_motif.py motif.pwm sequence.fa
import sys
import numpy as np
from Bio import SeqIO
from Bio import Motif

count = 0
path = '/home/msistaff/yang4414/project/NGStoolbox/chipseq/motif/TRANSFACT/'
ifp = open(path+'matrixlist.txt')
for line in ifp:
	name = line.rstrip()
	mat = np.loadtxt(path+name)
	prefix = name.split('.')[0]
	pfm = np.transpose(mat)
	np.savetxt(path+prefix+'.pfm', pfm, fmt='%d')

