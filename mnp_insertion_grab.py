import sys
import os
import re
import numpy as np
from matplotlib import pyplot as plt
import collections
from collections import defaultdict
import os.path

sample_files = open("PER17-mnp.txt",'r').read().splitlines()

sample_nums = []

for line in sample_files:
	fds = line.split('\t')
	sample_nums.append(fds[0])

for this_sample in sample_nums:
	print ("stats from " + this_sample)

	if os.path.exists(this_sample+".fastq.txt"):
		infile = open(this_sample+".fastq.txt",'r').read().splitlines()
		## This is to keep track of insertions
		newfile = open(this_sample+".insertions.txt", 'w')
		ins_seqs = []
		for line in infile [0::4]:
			seq = line.split('\n')[0]
			ins_seq = str(seq[53:73])
			ins_seqs.append(ins_seq)
		newfile.writelines("%s\n" % ins_seq for ins_seq in ins_seqs)
		newfile.close()