#This script was used to analyze the data in Figure 2c. 
#It reads in the data from the high throughput B->A pegRNA screen, then obtains counts of each insertion sequence at the genomic target site.
#Toward obtaining enrichment scores, we used two different scripts, mnp_insertion_grab.py and mnp_pivot.py, on HTS data for the raw hp-miniprep pool. This generated the data needed to make PER17_hpm_reads.txt (an input file below; lists each pegRNA in the hp-miniprep pool and the corresponding number of counts in the miniprep)
#Counts of each insertion sequence in the genome are tabulated alongside counts of the corresponding pegRNA in the hp-miniprep pool.  The resulting table enabled calculation of enrichment scores using equation (1) in the main text.

import sys
import os
import re
import numpy as np
import collections
from collections import defaultdict
import os.path
import pandas as pd
import Bio
from Bio.Seq import Seq

### USER INPUTS: ###

ins_location = 127                 #how far from the 5' end of the read is the 1st nt of the insertion?
upstream_seq = 'CATCA'             #what's the sequence immediately upstream of the prime editing target site?
min_distance1 = 7                  #what's the minimum allowable Hamming distance between a real insertion sequence and the site3 sequence?
min_distance2 = 7                  #what's the minimum allowable Hamming distance between a real insertion seequence and the 6xHis sequence (site3 and 6xHis sequences are both present in unedited cells because our starting cell line isn't homozygous for the 6xHis insertion)?
n_top_hits = 500000                #how many top hits (insertion seqs) do you want to report at the end?
lines_iterator = 10000000          #after how many lines should the memory be purged (to avoid memory errors while running the script)?


### DEFINING FUNCTIONS ###

def calc_H_dist(seq1,seq2): #calculates Hamming distance between two sequences
    distance = 0
    for a,b in zip(seq1,seq2):
        if a != b:
            distance += 1
    return distance


def check_distances (dist1,min_dist1,dist2,min_dist2,ins,hpm_dict): #compares the Hamming distances to their allowed values and prepares to make a dataframe of insertions
    list_for_dataframe = ['Hamming D too low', 1, 0, 0, 'N/A'] #reports to the dataframe if the insertion doesn't pass this filter
    if (dist1>=min_dist1) and (dist2>=min_dist2):
        list_for_df = []
        number = 1 #counter for n insertions of each sequence
        if ins in hpm_dict.keys(): #checks if the insertion sequence is in the hpm dictionary so we can calculate the enrichment scores
            hpm_counts = hpm_dict[ins]
        else:
            hpm_counts = 1
        list_for_dataframe = [insertion, number, distance1, distance2, hpm_counts]
    return(list_for_dataframe)



### GRABBING INSERTION SEQUENCES ###

site3 = 'CGTGCTCAGTCTGGGCCCCA' #sequence expected for unedited loci with wt site3
his6 = 'ATGATGGTGATGATGGCTTG' #sequence expected for unedited loci with 6xHis insertion only (no B->A edit)

ins_columns = ['ins seq', 'number', 'site3_HD', 'his6_HD', 'hpm_counts']  #creating a dataframe that will hold all the insertions
ins_df = pd.DataFrame(columns = ins_columns)

dna_seqs = [] #creating a list to hold all the DNA sequences (full reads)
ins_seqs = [] #creating a list to hold the 20-bp "insertion" sequences, including wt site3 and 1533 sequences
ins_for_df = [] #creating a list that will be used to make a dataframe to count hits


hpm_file = open("PER17_hpm_reads.txt") #input file name here (info about each pegRNA in the hp-miniprep pool and n reads from NGS on the miniprep pool)

hpm_dict = {} #making a dictionary where keys = hp-miniprep insertion sequences and values = number of reads in the hp-miniprep pool
for line in hpm_file:
    stripped_hpm_info = line.strip('\n').split('\t')
    hpm_upper = stripped_hpm_info[0].upper() #grabbing the 20-bp variable sequence from the pegRNA plasmid
    if hpm_upper not in hpm_dict:	
        hpm_dict[hpm_upper] = stripped_hpm_info[1]
    else:
        print('Error: duplicate hpm insertion sequences are present in the hpm file')

hpm_file.close()


with open("PER17_1_CKDL200162679-1a-D702-AK1681_HCJLGCCX2_L1_2.fq",'r') as sample_file: #input HTS data fastq file here
    counter = 0 #will keep track of total n reads
    L = -1  #will keep track of which lines in the fastq are DNA seqs and will tell the program when to purge the memory
    for line in sample_file:
        counter += 1
        L += 1
        if L%lines_iterator != (lines_iterator - 1): #the loop will start over after a user-defined n reads, to avoid memory errors
            if L%4 == 1:
                dna_seqs.append(line) #needs to cut off the line return at the end of the line

        else: #once the counter gets up to lines_iterator, tabulate the insertions, purge the memory, and restart the counter

            print('TOTAL n READS PARSED: ' + str(counter/4))
            print('dna seqs: ' + str(len(dna_seqs)))
            print('ins seqs: ' + str(len(ins_seqs)))
            print('ins_for_df: ' + str(len(ins_for_df)))

            nts_upstream = len(upstream_seq)
            for seq in dna_seqs:
                align_upstream = ins_location - nts_upstream - 1
                align_downstream = ins_location - 1

                if seq[align_upstream:align_downstream] == upstream_seq: #checks for the expected sequence upstream of the target site
                    insertion = seq[ins_location-1:ins_location+19] #grabs the 20-bp insertion
                    
                    distance1 = calc_H_dist(insertion,site3) #compares the 20-bp insertion to site3 seq and 6xHis seq
                    distance2 = calc_H_dist(insertion,his6)
                    list_for_df = check_distances (distance1,min_distance1,distance2,min_distance2,insertion,hpm_dict)
                    ins_for_df.append(list_for_df)

                elif seq[align_upstream + 1:align_downstream + 1] == upstream_seq: #adjusts target site position for reads that had a longer barcode on the illumina primer
                    insertion = seq[ins_location:ins_location+20]
                    
                    distance1 = calc_H_dist(insertion,site3) #compares the 20-bp insertion to site3 seq and his6 seq
                    distance2 = calc_H_dist(insertion,his6)
                    list_for_df = check_distances (distance1,min_distance1,distance2,min_distance2,insertion,hpm_dict)
                    ins_for_df.append(list_for_df)

                else:
                    continue

            sub_ins_df = pd.DataFrame(ins_for_df, columns = ins_columns) #makes a dataframe that will be purged repeatedly
            combined_df = pd.concat([ins_df, sub_ins_df]) #adds the dataframe that is about to be purged to the master dataframe
            counts_ins_df = combined_df.groupby(['ins seq','site3_HD','his6_HD','hpm_counts'], as_index=False)['number'].sum() #adds together n each insertion in the master dataframe, so it takes up less memory

            dna_seqs = [] #purging lists to prevent memory issues
            ins_seqs = [] 
            ins_for_df = [] 
            L = L%4 #restarts the iterator (to keep track of when memory should be purged)
            ins_df = counts_ins_df #updates the overall dataframe, so it holds the collapsed data set (identical insertions merged)


#this adds the "leftover" reads to the master dataframe once the memory purge iterator has gotten to the end of the input file
nts_upstream = len(upstream_seq)
for seq in dna_seqs: 
    align_upstream = ins_location - nts_upstream - 1
    align_downstream = ins_location - 1

    if seq[align_upstream:align_downstream] == upstream_seq: #checks for the expected sequence upstream of the target
        insertion = seq[ins_location-1:ins_location+19] #grabs the 20-bp insertion
        
        distance1 = calc_H_dist(insertion,site3) #compares the 20-bp insertion to site3 seq and his6 seq
        distance2 = calc_H_dist(insertion,his6)
        list_for_df = check_distances (distance1,min_distance1,distance2,min_distance2,insertion,hpm_dict)
        ins_for_df.append(list_for_df)

    elif seq[align_upstream + 1:align_downstream + 1] == upstream_seq: #adjusts target site position for reads that had a longer barcode on the illumina primer
        insertion = seq[ins_location:ins_location+20]
        
        distance1 = calc_H_dist(insertion,site3) #compares the 20-bp insertion to site3 seq and his6 seq
        distance2 = calc_H_dist(insertion,his6)
        list_for_df = check_distances (distance1,min_distance1,distance2,min_distance2,insertion,hpm_dict)
        ins_for_df.append(list_for_df)

    else:
        continue

sub_ins_df = pd.DataFrame(ins_for_df, columns = ins_columns) #makes a dataframe with the "leftover" data
combined_df = pd.concat([ins_df, sub_ins_df]) #adds the new dataframe to the master dataframe
counts_ins_df = combined_df.groupby(['ins seq','site3_HD','his6_HD','hpm_counts'], as_index=False)['number'].sum() #adds together n each insertion in the master dataframe

sorted_df = counts_ins_df.sort_values(by=['number','his6_HD'],ascending=[False,False])      
total_ins = sorted_df['number'].sum() #adds up total number of insertions in the dataframe

final_df = sorted_df.head(n_top_hits) #only take the top n hits and export them to an excel sheet
final_df.to_csv('PER17_dataframe_counted_insertions_GitHubtest.csv') 

print('total n reads inputted: ' + str(counter/4))
print('total n insertions grabbed: ' + str(total_ins))

