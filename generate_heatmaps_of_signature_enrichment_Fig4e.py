#This script was used to analyze the data in Figure 4e, corresponding to HTS data set PER37.

import sys
import os
import re
import numpy as np
import collections
from collections import defaultdict
import os.path
import pandas as pd
import math

### USER INPUTS: ###

ins_location_from_5prime = 34       #how far from the beginning of the F read is the 1st nt of the insertion?
upstream_seq = 'ATCAT'              #what's the sequence immediately upstream of the prime editing target site?
downstream_seq = 'TGATGG'           #what's the sequence immediately downstream of the prime editing target site?
HD_downstream_seq = 3               #what's an acceptable Hamming distance between the expected sequence downstream of the target site and the actual sequence (this allows the end of the insertion to be recognized despite sequencing errors)?
max_edits = 11                      #what's the maximum number of edits that can fit on your sequencing read?
max_length = 6                      #what's the maximum number of edits (array length) you want to analyze with heat maps?

output_file_name = '45fullrepeat' #this will be used to name the output files
F_sample_file = open("45repeatfull.fastq.assembled.fastq",'r').read().splitlines() #input file name here (fastq format; we analyzed fastq files that had F and R reads already assembled via PEAR)
signatures_file = open("PER37_signature_key.txt") #input file name here (key to convert signature to 1-letter symbols)


### DEFINING FUNCTIONS ###

def calc_H_dist(seq1,seq2): #defining a function to calculate Hamming distances between two sequences
    distance = 0
    for a,b in zip(seq1,seq2):
        if a != b:
            distance += 1
    return distance


def grab_dna_seqs(input_file): #reads in the DNA sequences from input file
    list_dna_seqs = []
    L = -1
    for line in input_file:
        L += 1
        if L%4 == 1:
            list_dna_seqs.append(line)
    return list_dna_seqs


def grab_insertions(input_seq, input_upstream_seq, input_ins_location, input_downstream_seq, input_HD_downstream_seq, input_max_edits): #determines how long insertion sequences are, then grabs the insertion sequences
    nts_upstream = len(input_upstream_seq) #checks how many nts upstream of the target site need to be considered to check for an insertion
    align_upstream = input_ins_location - nts_upstream 

    actual_upstream = input_seq[align_upstream:input_ins_location] #checks alignment upstream of the target site

    nts_downstream = len(input_downstream_seq) #will be used to check if the end of the insertion sequence has been found

    if actual_upstream == input_upstream_seq: #checking if an alignment was found upstream of the target site
        ins_location = input_ins_location
        found_alignment = 'true'  
    else:
        grabbed_ins = 'poor alignment'
        found_alignment = 'false'

    if found_alignment == 'true':   
        iterate = 'yes'
        end_edit = ins_location #this will hold the location corresponding to the end of the edited target
        n_edits = 0
        while iterate == 'yes': #iterate through a loop to find where the insertion sequence ends
            if n_edits == input_max_edits: #a maximum number of insertions can fit on the Illumina read, so insertions longer than this can't be grabbed
                iterate = 'no'
            elif calc_H_dist(input_seq[end_edit:end_edit+nts_downstream], input_downstream_seq) <= HD_downstream_seq: #checks if the expected sequence downstream of the target site is present; if so, this means the end of the insertion has been located
                iterate = 'no'
            else:
                n_edits += 1
                end_edit += 20

        if n_edits >= 1:
            length_ins = n_edits*20
            if n_edits <= (input_max_edits - 1):
                grabbed_ins = input_seq[ins_location: ins_location + length_ins]
            else:
                grabbed_ins = input_seq[ins_location: ins_location + length_ins - 17] #when you have the maximum number of edits, the read is too short to extract the full 20 nt ins, but we can still grab the signature   

        if n_edits == 0:
            grabbed_ins = 'no insertion'

        if grabbed_ins != 'no insertion' and grabbed_ins != 'poor alignment': #making sure the insertion sequences are as expected; if not, throw out the reads
            for i in range(1,n_edits):
                ins = grabbed_ins[3+20*(i-1):20*i]
                if i % 2 == 0:
                    if calc_H_dist(ins, 'GCCATCATCACCATCAT') >= 2:
                        grabbed_ins = 'wrong sequence'
                else:
                    if calc_H_dist(ins, 'ACTGGGCCATCTATCAT') >= 2:
                        grabbed_ins = 'wrong sequence'

    return grabbed_ins


def grab_signatures(input_ins, input_sig_key): 
        tot_length = len(input_ins)
        n_ins = math.ceil(tot_length/20)
        sig_3nt = 'no ins'
        converted_sigs = ''
        for i in range(0,n_ins): #grabbing the first 3 nt of each insertion (i.e., the signatures)
            sig_3nt = input_ins[0:3]
            input_ins = input_ins[20:] 
            if sig_3nt in input_sig_key.keys(): #converting the 3 nt signatures to 1 letter symbols
                sig_1letter = input_sig_key[sig_3nt]
                converted_sigs = converted_sigs + sig_1letter
            else:
                converted_sigs = converted_sigs + '.' #signatures that can't be found in the key are represented by .
                print(sig_3nt)
                print(tot_length)
   
        return converted_sigs



### GRABBING SIGNATURES AND CONVERTING THEM TO 1-LETTER SYMBOLS ###

F_dna_seqs = grab_dna_seqs(F_sample_file)

sig_dict = {} #making a dictionary where keys = 3-letter signatures and values = 1-letter symbols
for line in signatures_file:
    stripped_sigs = line.strip('\n').split('\t')
    if stripped_sigs[1].upper() not in sig_dict:	
        sig_dict[stripped_sigs[1].upper()] = stripped_sigs[0]
    else:
        print('Error: duplicate 3-letter signatures are present in the Signature Key')

F_ins_list = []
poor_alignments = 0
no_edits = 0
weird_edits = 0
for seq in F_dna_seqs: #grabbing insertion sequences from the fastq file
    grab_ins = grab_insertions(seq, upstream_seq, ins_location_from_5prime, downstream_seq, HD_downstream_seq, max_edits)
    if grab_ins == 'poor alignment':
        poor_alignments += 1
    elif grab_ins == 'no insertion':
        no_edits += 1
    elif grab_ins == 'wrong sequence':
        weird_edits += 1
    else:
        F_ins_list.append(grab_ins)


recorded_sigs = []
for ins in F_ins_list: #from the grabbed insertion sequences, extracting the 3 nt signatures, then converting them to 1 character symbols
    sig_symbols = grab_signatures(ins, sig_dict)
    recorded_sigs.append(sig_symbols)



### CREATING A PANDAS DATAFRAME TO REPORT HEATMAPS WITH % EACH SIGNATURE AT EACH POSITION, CALCULATED SEPARATELY FOR INSERTIONS OF DIFFERENT LENGTHS ###

heatmap_data = [] 
for length in range(1,max_length+1):
    list_of_length_L = []
    for record in recorded_sigs: #counting how many reads of each edit length we have
        if len(record) == length:
            list_of_length_L.append(record) 
    reads_length_L = len(list_of_length_L)

    for symbol in sig_dict.values(): #iterating through symbols, i.e., possible signature mutations
        symbol_fractions = []
        symbol_fractions.append(str(symbol)) #1st column in the dataframe will show the symbol (signature)
        symbol_fractions.append(length) #2nd column in the dataframe will show the insertion length
        symbol_fractions.append(reads_length_L) #3rd column in the dataframe reports the number of reads of each length
        
        if reads_length_L >= 1: #put this here to avoid a division by 0 error below
            for i in range(0,length): #looking at the insertion at each position -> what fraction of the reads (of ins length L) have the symbol of interest at each position?
                symbol_count = 0
                for record in list_of_length_L:
                    if record[i] == symbol:
                        symbol_count += 1
                symbol_frequency = symbol_count / reads_length_L
                symbol_fractions.append(symbol_frequency)
            for L in range(length, max_length):
                symbol_fractions.append('-') #fill in the heat map table with "-" for array lengths that aren't relevant for the given row
            heatmap_data.append(symbol_fractions)  
            

heatmap_dataframe = pd.DataFrame(heatmap_data)
heatmap_dataframe.sort_values(by=[0,1]).to_csv(output_file_name + '_heatmap_dataframe.csv')



### CREATING ANOTHER DATAFRAME, THIS TIME REPORTING 1 ROW OF DATA FOR EACH SIGNATURE (all reads used to calculate % signature at each position) ###

heatmap_data = [] 

for symbol in sig_dict.values(): #iterating through symbols, i.e., possible signature mutations
    symbol_fractions = []
    symbol_fractions.append(str(symbol)) #1st column in the dataframe will show the symbol (signature)

    for i in range(1,max_length+1): #looking at the insertion at each position -> what fraction of the reads (of ins length L) have the symbol of interest at each position?
        symbol_count = 0
        record_count = 0
        for record in recorded_sigs:
            if len(record) >= i:
                record_count += 1
                if record[i-1] == symbol: #record index of 0 corresponds to insertion number 1
                    symbol_count += 1
        if record_count >= 1:
            symbol_frequency = symbol_count / record_count #what fraction of records that are long enough to be considered have the symbol of interest at the position of interest?
        else:
            symbol_frequency = '-'
        symbol_fractions.append(symbol_frequency)
    heatmap_data.append(symbol_fractions)  

symbol_fractions = [] #adding a counter that will tally up n insertions used to calculate the % each symbol at each position
symbol_fractions.append('n insertions')
for i in range(1, max_length+1):
    record_count = 0
    for record in recorded_sigs:
        if len(record) >= i:
            record_count += 1
    symbol_fractions.append(record_count)
heatmap_data.append(symbol_fractions)

heatmap_2_dataframe = pd.DataFrame(heatmap_data)
heatmap_2_dataframe.sort_values(by=[0]).to_csv(output_file_name + '_heatmap_all_reads_for_%_each_position.csv')



### CREATING A DATAFRAME TO REPORT ALL OF THE RECORDS WITH 1-LETTER SYMBOLS REPRESENTING EACH SIGNATURE ###

signatures_data = []
for record in recorded_sigs:
    signatures_list = []
    signatures_list.append(len(record)) #1st column in dataframe will report the insertion length
    signatures_list.append(record) #2nd column in dataframe will report the inserted signatures (in 1-letter symbol format)
    signatures_list.append(int(1)) #3rd column will be used to count the number of reads of each insertion
    signatures_data.append(signatures_list)

columns = ['length', 'signatures (1 letter)', 'counts']
signatures_dataframe = pd.DataFrame(signatures_data, columns = columns)
grouped_signatures = signatures_dataframe.groupby(by=['signatures (1 letter)','length']).sum()
grouped_signatures.to_csv(output_file_name + '_signatures_dataframe.csv')


  
#REPORTING ALIGNMENT AND EDITING EFFICIENCY RESULTS                 
print('total n reads: ' + str(len(F_dna_seqs)))
print('total n insertions: ' + str(len(F_ins_list)))
print('poor alignments: ' + str(poor_alignments))
print('reads with no edits:' + str(no_edits))
print('reads with genomic sequences, big indels, etc.: ' + str(weird_edits))

signatures_file.close()