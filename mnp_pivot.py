import sys
import os
import re
import numpy as np
from matplotlib import pyplot as plt
import collections
from collections import defaultdict
import os.path
import pandas as pd

sample_files = open("PER17-mnp.txt",'r').read().splitlines()

sample_nums = []

for line in sample_files:
	fds = line.split('\t')
	sample_nums.append(fds[0])

for this_sample in sample_nums:
	print ("stats from " + this_sample)

	if os.path.exists(this_sample+".insertions.txt"):
		infile = open(this_sample+".insertions.txt",'r').read().splitlines()
		## This is to keep track of insertions
		newfile = open(this_sample+".pivot.csv", 'w')
		ins_seqs = []
		val_seqs = []
		for line in infile:
			seq = line.split('\n')[0]
			ins_seqs.append(seq)
			val_seqs.append(1)
		df = pd.DataFrame({'foo':ins_seqs,
							'bar':val_seqs})
		table = pd.pivot_table(df,values='bar',index='foo',aggfunc=np.sum)
		table.to_csv(newfile,index=True)
		newfile.close()