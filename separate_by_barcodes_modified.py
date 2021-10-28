import sys
import os

###Just reads the barcode files
fwdbarcode = []
revbarcode = []
fnames = []
barcodefile = open(sys.argv[4],'r').read().splitlines()
for line in barcodefile:
	fds = line.split()
	fnames.append(fds[0])
	fwdbarcode.append(fds[1])
	revbarcode.append(fds[2])
	#print fds 
'''
fwdbarcode=["tatctca","tctcatc","agaggag","aactacg","tgcgac","gaaact","gactga","actgac","gtcagt","tgatccg","aggcatg","agatagg","tatctca","tctcatc","agaggag","aactacg","tgcgac","gaaact","gactga","actgac","gtcagt","tgatccg","actgac","gtcagt","tgatccg","aggcatg","agatagg","tatctca","tctcatc","agaggag","aactacg","tgcgac","gaaact", "actgac", "gtcagt", "tgatccg", "aggcatg"]
revbarcode=["agtttc","agtttc","agtttc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tcagtc","tttcaa","tttcaa","tttcaa","tttcaa","tttcaa","tttcaa","tttcaa","gtcgca","gtcgca","gtcgca","gtcgca","gtcgca","gtcgca","gtcgca","gtcgca","agtttc","agtttc","agtttc", "agtttc", "agtttc", "agtttc", "agtttc"]
fnames = ["34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","17","18","19","20","21","22","23","24","25","26","27", "29", "30", "31", "32"]
'''

barcodeset = list(zip(fwdbarcode,revbarcode))
foutfiles = []
routfiles = []
namesind = 0

if sys.argv[3]=="fastq":
	fastqtype = True
elif sys.argv[3]=="fasta":
	fastqtype = False

###makes a directory for the output, if that directory does not exist
directory = sys.argv[5]
# if not os.path.exists(directory):
#     os.makedirs(directory)

###generates the output files  fastq or fasta
for m in barcodeset:
	ftempfile = open(directory+"/"+fnames[namesind]+"_forward_"+m[0]+"_"+m[1]+"." + sys.argv[3], 'w+')
	rtempfile = open(directory+"/"+fnames[namesind]+"_reverse_"+m[0]+"_"+m[1]+"." + sys.argv[3], 'w+')
	foutfiles.append(ftempfile)
	routfiles.append(rtempfile)
	namesind+=1

###Reads the read files
reads1 = open(sys.argv[1],'r')
reads2 = open(sys.argv[2],'r')


###Starts going through the forward reads and making a table
print("Files are generated")
cntr1 = 0
fwddict = {}
templine=""
for l1 in reads1:
	if cntr1%4==0 and cntr1!=0:
		fds = templine.split("\t")
		name = fds[0] + "\t" +fds[2] + "\t" + fds[3]
		seq = fds[1].upper()
		flag = False
		fbarcode = []
		#print name.split()[0]
		for brcd in fwdbarcode:
			if seq.startswith(brcd.upper()): # +"ATC"):
				fbarcode.append(brcd)
				flag = True
				#break
		if flag:
			fwddict[name.split("\t")[0].split()[0]] = (name, seq, fbarcode)
		templine=l1.strip()
	elif cntr1==0:
		templine = l1.strip()
	else:
		templine+="\t"+l1.strip()
	cntr1+=1
	if cntr1%4000000 == 0:
		print(str(cntr1/4)+ " forward reads are checked")


###Reverse reads are going to be read, matched to forward ones, and write in output file
print("Forward reads are parsed")
cntr2 = 0
templine=""
for l2 in reads2:
	if cntr2%4==0 and cntr2!=0:
		fds=templine.split("\t")
		name = fds[0]+ "\t" +fds[2] + "\t" + fds[3]
		rseq = fds[1].upper()
		index = cntr2
		flag = False
		rbarcode = []
		if name.split("\t")[0].split()[0] in fwddict:
			#print "I'mhere"
			for brcd in revbarcode:
				if rseq.startswith(brcd.upper()): #+"CCC"):
					rbarcode.append(brcd)
					flag = True
					#break
			if flag:
				othername, fseq, fbarcode = fwddict[name.split("\t")[0].split()[0]]
				frb = set(zip(fbarcode, rbarcode))
				for pair in frb:
					fb, rb = pair
					if (fb,rb) in barcodeset:
						if fastqtype:
							foutfiles[barcodeset.index((fb,rb))].write(othername.split("\t")[0]  + "\n"+ fseq[len(fb):] + "\n" + othername.split("\t")[1] + "\n" + othername.split("\t")[2][len(fb):] + "\n" )
							routfiles[barcodeset.index((fb,rb))].write(name.split("\t")[0] + "\n"+ rseq[len(rb):] + "\n" + name.split("\t")[1]+"\n" + name.split("\t")[2][len(rb):] + "\n" )
						else:
							foutfiles[barcodeset.index((fb,rb))].write(">"+othername.split("\t")[0]  + "\n"+ fseq[len(fb):] + "\n" )
							routfiles[barcodeset.index((fb,rb))].write(">"+name.split("\t")[0] + "\n"+ rseq[len(rb):] + "\n" )
		templine=l2.strip()
	elif cntr2 == 0:
		templine=l2.strip()
	else:
		templine+="\t"+l2.strip()
	cntr2+=1
	if cntr2%4000000 == 0:
		print(str(cntr2/4)+ " reverse reads are checked")

print("Yah Yah Yah, we are done")
for f in foutfiles:
	f.close()
for f in routfiles:
	f.close()
