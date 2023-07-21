# This script was used to extract insertion sequences for Figure 3a-f, Figure 4, Figure 5, Extended Data Figure 1b, Extended Data Figure 3a, Extended Data Figure 4b and d, Extended Data Figure 5a-b, Extended Data Figure 6, Extended Data Figure 8, and Extended Data Figure 9.

# This script reads in fastq files (either from Illumina or Nanopore sequencing runs), then extracts the peCHYRON insertion sequences (records) from the raw reads.
# Insertion sequences are compared to correct peCHYRON edit sequences (3 nt randomized signature mutation + 17 nt constant propagator sequence). 
# If an insertion contains correct peCHYRON edits, the signature mutations are extracted, and each signature is converted to a 1-letter symbol to simplify interpretation of results. 
# The signature_key.txt input file defines the conversion between 3 nt signature mutations and 1-letter symbols. 

# The records with signatures converted to 1-letter symbols are outputted as a [sample name]_signatures_dataframe.csv file. 
# A count of how many times each signature appeared in the records is outputted as a [sample name]_entropy_dataframe.csv file, to enable Shannon entropy calculations.
# All insertion sizes (including those that aren't perfect inecrements of 20 nt) are tallied and outputted as an [sample name]_indel_histogram.csv file, to enable assessment of the abundance of correct vs. corrupted insertion sizes.
# Other output files ([sample name]_insertions_dataframe.csv, and [sample name]_overall_stats.csv) are also provided for quality control. 

import sys
import os
import re
import numpy as np
import collections
from collections import defaultdict
import os.path
import pandas as pd
import math
import Bio
from Bio.Seq import Seq


### USER INPUTS: ###

ins_location_from_5prime = 350       #how far from the beginning of the F read is 1st nt of the insertion? 
ins_location_sliding_window = 150    #how many nt upstream and downstream of the expected target site do you want to check for an alignment?
upstream_seq = 'CATCATCACCATCAT'     #what's the sequence immediately upstream of the prime editing target site?
downstream_seq = 'TGATGGCAGA'       #what's the sequence immediately downstream of the prime editing target site?
HD_upstream_seq = 3                 #what's an acceptable Hamming Distance between the expected sequence upstream of the target site and the actual sequence? This will be used to check for alignment before extracting the insertion sequences.
HD_downstream_seq = 2               #what's an acceptable Hamming Distance between the expected sequence downstream of the cutsite and the actual sequence? This will be used to identify the end of each insertion sequence.
max_edits = 11                      #what's the maximum number of edits (20 bp insertions) you want to consider?  This will determine the bounds of the indel histogram. 
base_input_file_name = 'PER69-1-'   #what's the constant part of the fastq file names?  For our fastq files, this string immediately precedes the sample number.  Here you are providing the names of the input files.
input_identifiers = [1,2,7,8,9]              #what sample numbers are you analyzing?  You can manually type in a list of sample numbers and comment out the For loop below, or leave this list empty and populate it with the For loop below.
#for num in range(1,9):
#    input_identifiers.append(num)

signatures_file = open("PER44_signature_key.txt") #input file name here (key to convert 3-bp signatures to 1-letter symbols)


### DEFINING FUNCTIONS ###

def calc_H_dist(seq1,seq2): 
    '''calculates the hamming distances between two sequences'''
    distance = 0
    for a,b in zip(seq1,seq2):
        if a != b:
            distance += 1
    return distance


def grab_dna_seqs(input_file): 
    '''reads in the DNA sequences from input files'''
    list_dna_seqs = []
    L = -1
    for line in input_file:
        L += 1
        if L%4 == 1:
            list_dna_seqs.append(line)
    return list_dna_seqs


def grab_insertions(input_seq, input_upstream_seq, input_ins_location, input_downstream_seq, input_HD_downstream_seq, input_HD_upstream_seq, input_max_edits): 
    '''determines the lengths of the insertion sequences, then grabs the insertion sequences'''
    
    grabbed_ins = [] #will report 1) whether a good alignment with the target site was found and 2) the actual insertion sequence if applicable
    
    nts_upstream = len(input_upstream_seq) #will be used to check for alignment alignment upstream of the target site
    align_upstream = input_ins_location - nts_upstream 
    actual_upstream = input_seq[align_upstream:input_ins_location] 

    nts_downstream = len(input_downstream_seq) #will be used to check for alignment downstream of the target site (i.e., to find the end of the insertion sequence)

    found_alignment = 'false'
    if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq: #checking if the expected sequence upstream of the edit site is present at the expected position, within an allowable Hamming distance
        ins_location = input_ins_location
        found_alignment = 'true'  

    if found_alignment == 'false': #checking for the expected sequence upstream of the target site, using a scanning window to detect if it is present at a different position than usual (this is especially important for Nanopore data, which is prone to indels due to sequencing error)
        for i in range(-ins_location_sliding_window, ins_location_sliding_window):
            align_upstream = input_ins_location - nts_upstream + i
            actual_upstream = input_seq[align_upstream:align_upstream + nts_upstream]
            if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq:
                found_alignment = 'true'
                ins_location = input_ins_location + i    
                break

    if found_alignment == 'false': #if an alignment wasn't found above, try aligning the reverse complement of the read, since Nanopore indiscriminately sequences both strands of DNA 
        seq = Seq(input_seq)
        rev_comp_seq = seq.reverse_complement()
        align_upstream = input_ins_location - nts_upstream

        actual_upstream = rev_comp_seq[align_upstream:input_ins_location] #checks alignment upstream of the target site
        if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq: #checking if the expected sequence upstream of the edit site is present at the expected position
            ins_location = input_ins_location
            found_alignment = 'true' 
            input_seq = rev_comp_seq

    if found_alignment == 'false': #checking for the expected sequence upstream of the target site in the reverse complement of the read, using a scanning window to detect if it is present at a different position than usual
        for i in range(-ins_location_sliding_window, ins_location_sliding_window):
            align_upstream = input_ins_location - nts_upstream + i
            actual_upstream = rev_comp_seq[align_upstream:align_upstream + nts_upstream]
            if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq:
                    found_alignment = 'true'
                    ins_location = input_ins_location + i
                    input_seq = rev_comp_seq
                    break
                   
    if found_alignment == 'false':
        grabbed_ins = ['failed initial alignment', input_seq] #report a failed alignment if all of the above alignment attempts were unsuccessful

    if found_alignment == 'true': #determine the length of the insertion for reads that successfully aligned  
        iterate = 'yes'
        end_edit = ins_location #this will hold the location corresponding to the end of the edited target
        n_edits = 0
        while iterate == 'yes': 
            if n_edits == input_max_edits: #a maximum number of insertions can fit on the Illumina read, so insertions longer than this can't be grabbed
                iterate = 'no'
            elif calc_H_dist(input_seq[end_edit:end_edit+nts_downstream], input_downstream_seq) <= input_HD_downstream_seq: #checks if the expected sequence downstream of the target site is present at the expected location for a 20-bp insertion, within an acceptable Hamming distance
                iterate = 'no'
            else:
                downstream_match = 'scanning' #if the expected sequence downstream of the target site wasn't found exactly where it's expected (an increment of 20 nt from the nick site), scan for it in every possible position
                for i in range(-19,19): #keep the scanning range less than 20 so it doesn't recognize an insertion from the downstream editing event, which should be checked in the next iteration instead
                    if calc_H_dist(input_seq[end_edit+i:end_edit+nts_downstream+i], input_downstream_seq) <= HD_downstream_seq: #checks for the expected sequence downstream of the target site
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
                grabbed_ins_seq = input_seq[ins_location: ins_location + 20*(input_max_edits-1) + 3] #when you have the maximum number of edits, an Illumina read is too short to extract the full 20 nt for the last insertion, but we can still grab the signature   
            grabbed_ins = ['good alignment', grabbed_ins_seq]

        if n_edits == 0:
            grabbed_ins = ['good alignment', ''] #indicates 0nt insertion

        grabbed_ins_description = 'assume correct sequence'
        if grabbed_ins[1] != '' and grabbed_ins[0] != 'failed initial alignment': #making sure the insertion sequences are as expected; if not, throw out the reads
            if len(grabbed_ins[1]) >= 20: #the total length of the insertion has to be at least 20 in order to be aligned
                alignable_edits = math.floor(length_ins/20) #if there's an insertion that's less than 20 nt, it will be the terminal insertion.  We can check alignments for the early insertions and ignore the last insertion.
                for i in range(0,alignable_edits):
                    ins = grabbed_ins[1][3+20*(i):20*(i+1)] #aligning the propagator sequence to the expected sequences
                    if i % 2 == 1:
                        if calc_H_dist(ins, 'GCCATCATCACCATCAT') >= 2: #B>A insertion
                            grabbed_ins_description = 'wrong sequence'
                    else:
                        if calc_H_dist(ins, 'ACTGGGCCATCTATCAT') >= 2: #A>B insertion
                            grabbed_ins_description = 'wrong sequence'
                        
            else:
                grabbed_ins_description = 'wrong sequence'
                    
        if grabbed_ins_description == 'wrong sequence':
            raw_seq = grabbed_ins[1]
            grabbed_ins = ['good alignment', 'wrong sequence', str(raw_seq)]

    return grabbed_ins


def grab_signatures(input_ins, input_sig_key): 
        '''reads in the full insertion sequences (including the propagator sequences), grabs only the 3-nt signatures, and converts them to 1-letter symbols'''
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

sig_dict = {} #making a dictionary where keys = 3-letter signatures and values = 1-letter symbols
for line in signatures_file:
    stripped_sigs = line.strip('\n').split('\t')
    if stripped_sigs[1].upper() not in sig_dict:	
        sig_dict[stripped_sigs[1].upper()] = stripped_sigs[0]
    else:
        print('Error: duplicate 3-letter signatures are present in the Signature Key')

for file in input_identifiers: #iterate through all of the input files in the directory, as long as their file names match the base_input_file_name and input_identifiers specified up top
    output_file_name = '{}{}'.format(base_input_file_name, file) #this will be used to name the output files
    F_sample_file = open("{}{}.fastq".format(base_input_file_name, file),'r').read().splitlines() #input file name here (fastq F reads)
    
    F_dna_seqs = grab_dna_seqs(F_sample_file)

    F_ins_list = []
    failed_alignments = []
    irregular_size_ins = []
    list_weird_edits = []
    no_edits = 0
    weird_edits = 0
    for seq in F_dna_seqs: #grabbing insertion sequences from the fastq file
        grab_ins = grab_insertions(seq, upstream_seq, ins_location_from_5prime, downstream_seq, HD_downstream_seq, HD_upstream_seq, max_edits)
        if grab_ins[0] == 'failed initial alignment':
            failed_alignments.append(grab_ins[1])
        elif grab_ins[1] == '':
            no_edits += 1
        elif grab_ins[1] == 'wrong sequence':
            weird_edits += 1
            list_weird_edits.append(grab_ins[2])
        else:
            F_ins_list.append(grab_ins[1])
            if len(grab_ins[1])%20 != 0:
                irregular_size_ins.append(grab_ins[1])

    recorded_sigs = []
    for ins in F_ins_list: #from the grabbed insertion sequences, extracting the 3-bp signatures, then converting them to 1 character symbols
        sig_symbols = grab_signatures(ins, sig_dict)
        recorded_sigs.append(sig_symbols)


    ### EXPORTING THE RAW LIST OF RAW INSERTIONS TO AN EXCEL FILE ###
    list_of_ins_list = []

    for entry in F_ins_list:
        list_of_ins_list.append([entry])

    for entry in list_weird_edits:
        list_of_ins_list.append([entry]) #this will add insertions that were <20 bp 

    insertions_dataframe = pd.DataFrame(list_of_ins_list)
    insertions_dataframe.to_csv(output_file_name + '_insertions_dataframe.csv')



    ### TALLYING UP INSTANCES OF EACH INSERTION LENGTH ###

    indel_histogram = [] #this will be a list of lists, used to make a pandas dataframe
    for n in range (-21,0):
        ins_length_L = [] #this will hold the number of instances of the insertion length in question
        count = 0
        for ins in F_ins_list:
            if len(ins) == n:
                count += 1
        for ins in list_weird_edits:
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
        for ins in list_weird_edits:
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

    signatures_list = []
    signatures_list.append(0) #adding number of reads with no edit to the dataframe
    signatures_list.append('')
    signatures_list.append(no_edits)
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

    report_headers = ['total n reads', 'total n insertions with aligned propagators', 'failed alignments', 'reads with no edits', 'reads without expected propagator sequence']
    report_results = [str(len(F_dna_seqs)), str(len(F_ins_list)), str(len(failed_alignments)), str(no_edits), str(weird_edits)]
    formatting_for_df = []
    formatting_for_df.append(report_results)
    reports_dataframe = pd.DataFrame(formatting_for_df, columns = report_headers)
    reports_dataframe.to_csv(output_file_name + '_overall_stats.csv')
    

signatures_file.close()