#This script was used to analyze the data in Figure 3. It quantifies the number of edits of each insertion length (including insertions that aren't increments of 20 bp), gives the frequency of each signature mutation across all reads (for Shannon entropy calculations), and reports every record extracted from all reads.

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

ins_location_from_5prime = 34       #how far from the beginning of the F read is 1st nt of the insertion? 
upstream_seq = 'TCACCATCAT'         #what's the sequence immediately upstream of the prime editing target site?
downstream_seq = 'TGATGGCAGA'       #what's the sequence immediately downstream of the prime editing target site?
HD_upstream_seq = 2                 #what's an acceptable Hamming Distance between the expected sequence upstream of the cutsite and the actual sequence? This will be used to check for alignment before extracting the insertion sequences.
HD_downstream_seq = 2               #what's an acceptable Hamming Distance between the expected sequence downstream of the cutsite and the actual sequence? This will be used to identify the end of each insertion sequence.
max_edits = 11                      #what's the maximum number of edits (20 bp insertions) that can fit on your sequencing read?

output_file_name = 'PER37_120' #this will be used to name the output files
F_sample_file = open("120.fastq.assembled.fastq",'r').read().splitlines() #input file name here (fastq F reads)
signatures_file = open("PER37_signature_key.txt") #input file name here (key to convert 3-bp signatures to 1-letter symbols)


### DEFINING FUNCTIONS ###

def calc_H_dist(seq1,seq2): #defining a function to calculate hamming distances between two sequences
    distance = 0
    for a,b in zip(seq1,seq2):
        if a != b:
            distance += 1
    return distance


def grab_dna_seqs(input_file): #reads in the DNA sequences from input files
    list_dna_seqs = []
    L = -1
    for line in input_file:
        L += 1
        if L%4 == 1:
            list_dna_seqs.append(line)
    return list_dna_seqs


def grab_insertions(input_seq, input_upstream_seq, input_ins_location, input_downstream_seq, input_HD_downstream_seq, input_HD_upstream_seq, input_max_edits): #determines how long insertion sequences are, then grabs the insertion sequences
    grabbed_ins = [] #will report 1) whether a good alignment was found and 2) the actual insertion sequence if applicable
    nts_upstream = len(input_upstream_seq) #checks how many nts upstream of the target site need to be considered to check for an insertion
    align_upstream = input_ins_location - nts_upstream 

    actual_upstream = input_seq[align_upstream:input_ins_location] #checks alignment upstream of the target site

    nts_downstream = len(input_downstream_seq) #will be used to check if the algorithm has found the end of the insertion sequence

    if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq: #checking if the expected sequence upstream of the edit site is present
        ins_location = input_ins_location
        found_alignment = 'true'  
    else:
        grabbed_ins = ['failed initial alignment', input_seq]
        found_alignment = 'false'

    if found_alignment == 'true':   
        iterate = 'yes'
        end_edit = ins_location #this will hold the location corresponding to the end of the edited target
        n_edits = 0
        while iterate == 'yes': #iterate through a loop to find where the insertion sequence ends
            if n_edits == input_max_edits: #a maximum number of insertions can fit on the Illumina read, so insertions longer than this can't be grabbed
                iterate = 'no'
            elif calc_H_dist(input_seq[end_edit:end_edit+nts_downstream], input_downstream_seq) <= input_HD_downstream_seq: #checks if the expected sequence downstream of the target site is present at the expected location for a 20-bp insertion
                iterate = 'no'
            else:
                downstream_match = 'scanning'
                for i in range(-19,19): #keep the range less than 20 so it doesn't recognize an insertion from the downstream editing event, which should be checked in the next iteration instead
                    if calc_H_dist(input_seq[end_edit+i:end_edit+nts_downstream+i], input_downstream_seq) <= HD_downstream_seq: #checks for the expected sequence downstream of the target site, using a scanning window to locate it
                        downstream_match = 'found'
                        iterate = 'no'
                        end_edit += i #this will capture insertions that aren't multiples of 20 bp
                        n_edits += 1
                if downstream_match != 'found': #if the expected sequence downstream of the target site hasn't been found yet, move down 20 nt and look again
                    n_edits += 1 
                    end_edit += 20

        if n_edits >= 1:
            length_ins = end_edit - ins_location #the end edit position was calculated above, so we can subtract the location of the first nt of the insertion to get the total ins length 
            if n_edits <= (input_max_edits - 1):
                grabbed_ins_seq = input_seq[ins_location: ins_location + length_ins]
            else:
                grabbed_ins_seq = input_seq[ins_location: ins_location + 20*(input_max_edits-1) + 3] #when you have the maximum number of edits, the read is too short to extract the full 20 nt for the last ins, but we can still grab the signature   
            grabbed_ins = ['good alignment', grabbed_ins_seq]

        if n_edits == 0:
            grabbed_ins = ['good alignment', ''] #indicates 0nt insertion

        grabbed_ins_description = 'assume correct sequence'
        if grabbed_ins[1] != '' and grabbed_ins[0] != 'poor alignment': #making sure the insertion sequences are as expected; if not, throw out the reads
            alignable_edits = math.floor(length_ins/20) #if there's an insertion that's less than 20 nt, it will be the terminal insertion.  We can check alignments for the early insertions and ignore the last insertion.
            for i in range(1,alignable_edits):
                ins = grabbed_ins[1][3+20*(i-1):20*i] #aligning the propagator sequence to the expected sequences
                if i % 2 == 0:
                    if calc_H_dist(ins, 'GCCATCATCACCATCAT') >= 2:
                        grabbed_ins_description = 'wrong sequence'
                else:
                    if calc_H_dist(ins, 'ACTGGGCCATCTATCAT') >= 2:
                        grabbed_ins_description = 'wrong sequence'
                    
        if grabbed_ins_description == 'wrong sequence':
            grabbed_ins = ['good alignment', 'wrong sequence']

    return grabbed_ins


def grab_signatures(input_ins, input_sig_key): 
        tot_length = len(input_ins)
        n_ins = math.ceil(tot_length/20)
        sig_3nt = 'no ins'
        converted_sigs = ''

        if len(input_ins)%20 == 0 or len(input_ins) == (max_edits-1)*20+3: #this will process signatures for normal length insertions or insertions that were too long to fully fit on the illumina read, but are presumably normal length
            for i in range(0,n_ins): #grabbing the first 3 nt of each insertion (i.e., the signatures)
                sig_3nt = input_ins[0:3]
                input_ins = input_ins[20:] 
                if sig_3nt in input_sig_key.keys(): #converting the 3 nt signatures to 1 letter symbols
                    sig_1letter = input_sig_key[sig_3nt]
                    converted_sigs = converted_sigs + sig_1letter
                else:
                    converted_sigs = converted_sigs + '.' #signatures that can't be found in the key are represented by .

        else: #this will process signatures for recording loci that have a non-20-bp insertion.  It will simply grab signatures for all insertions except the terminal, non-20-bp insertion
            for i in range(0,n_ins-1):
                sig_3nt = input_ins[0:3]
                input_ins = input_ins[20:] 
                if sig_3nt in input_sig_key.keys(): #converting the 3 nt signatures to 1 letter symbols
                    sig_1letter = input_sig_key[sig_3nt]
                    converted_sigs = converted_sigs + sig_1letter
                else:
                    converted_sigs = converted_sigs + '.' #signatures that can't be found in the key are represented by .

                
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
failed_alignments = []
irregular_size_ins = []
no_edits = 0
weird_edits = 0
for seq in F_dna_seqs: #grabbing insertion sequences from the fastq file
    grab_ins = grab_insertions(seq, upstream_seq, ins_location_from_5prime, downstream_seq, HD_downstream_seq, HD_upstream_seq, max_edits)
    if grab_ins[0] == 'failed initial alignment':
        failed_alignments.append(grab_ins[0])
    elif grab_ins[1] == '':
        no_edits += 1
    elif grab_ins[1] == 'wrong sequence':
        weird_edits += 1
    else:
        F_ins_list.append(grab_ins[1])
        if len(grab_ins[1])%20 != 0:
            irregular_size_ins.append(grab_ins[1])

recorded_sigs = []
for ins in F_ins_list: #from the grabbed insertion sequences, extracting the 3-bp signatures, then converting them to 1 character symbols
    sig_symbols = grab_signatures(ins, sig_dict)
    recorded_sigs.append(sig_symbols)



### TALLYING UP INSTANCES OF EACH INSERTION LENGTH ###

indel_histogram = [] #this will be a list of lists, used to make a pandas dataframe
for n in range (-19,0):
    ins_length_L = [] #this will hold the number of instances of the insertion length in question
    count = 0
    for ins in F_ins_list:
        if len(ins) == n:
            count += 1
    ins_length_L.append(n)
    ins_length_L.append(count)
    indel_histogram.append(ins_length_L)

ins_length_L = []
ins_length_L.append('0') #previously counted how many reads had 0 nt insertions.  Adding it to the dataframe.
ins_length_L.append(no_edits)
indel_histogram.append(ins_length_L)

for n in range(1,(max_edits-1)*20+3): #goes up to the maximum insertion length this pipeline can extract
    ins_length_L = [] #this will hold the number of instances of the insertion length in question
    count = 0
    for ins in F_ins_list:
        if len(ins) == n:
            count += 1
    ins_length_L.append(n)
    ins_length_L.append(count)
    indel_histogram.append(ins_length_L)


indel_columns = ['ins length (nt)','counts']
ins_size_dataframe = pd.DataFrame(indel_histogram, columns = indel_columns)
ins_size_dataframe.to_csv(output_file_name + '_indel_histogram.csv')



### CREATING A DATAFRAME TO REPORT ALL RECORDS, USING THE 1-LETTER SYMBOLS REPRESENTING EACH SIGNATURE ###

signatures_data = []
for record in recorded_sigs:
    if record != "": #this filters out records that consist of a single insertion that is less than 20 bp, causing no signature to be extracted above; hence, an empty string
        signatures_list = []
        signatures_list.append(len(record)) #1st column in dataframe will report the insertion length
        signatures_list.append(record) #2nd column in dataframe will report the inserted signatures (in 1-letter symbol format)
        signatures_list.append(int(1)) #3rd column will be used to count the number of reads of each insertion
        signatures_data.append(signatures_list)

columns = ['length', 'signatures (1 letter)', 'counts']
signatures_dataframe = pd.DataFrame(signatures_data, columns = columns)
grouped_signatures = signatures_dataframe.groupby(by=['signatures (1 letter)','length']).sum()
grouped_signatures.to_csv(output_file_name + '_signatures_dataframe.csv')



### COUNTING THE NUMBER OF OCCURENCES OF EACH SIGNATURE, TO BE USED FOR ENTROPY CALCULATIONS ###

entropy_data = []
for symbol in sig_dict.values():
    symbol_data = []
    odd_count = 0 
    even_count = 0
    for record in recorded_sigs:
        for i in range(0,len(record)): 
            if i%2 == 0: #considering signatures at odd and even steps separately, since some signatures are reused for alternating steps and are technically different for the two steps
                if record[i] == symbol:
                    odd_count += 1
        for i in range(1,len(record)): 
            if i%2 == 1:
                if record[i] == symbol:
                    even_count += 1
    
    symbol_data.append('odd-' + symbol) #adding the counts to the dataframe
    symbol_data.append(odd_count)
    entropy_data.append(symbol_data)

    symbol_data = []
    symbol_data.append('even-' + symbol)
    symbol_data.append(even_count)
    entropy_data.append(symbol_data)


un_ID_sig_count = 0 #also counting up the occurences of unidentified signatures
for record in recorded_sigs:
    for i in range(0, len(record)):
        if record[i] == ".":
            un_ID_sig_count += 1
symbol_data = []
symbol_data.append('.')
symbol_data.append(un_ID_sig_count)
entropy_data.append(symbol_data)


columns = ['signature','counts']
entropy_dataframe = pd.DataFrame(entropy_data, columns = columns)
entropy_dataframe.to_csv(output_file_name + '_entropy_dataframe.csv')
     
     
#REPORTING ALIGNMENT AND EDITING EFFICIENCY RESULTS 

print('total n reads: ' + str(len(F_dna_seqs)))
print('total n insertions: ' + str(len(F_ins_list)))
print('poor alignments: ' + str(len(failed_alignments)))
print('reads with no edits:' + str(no_edits))
print('reads with genomic sequences, big indels, etc.: ' + str(weird_edits))

signatures_file.close()