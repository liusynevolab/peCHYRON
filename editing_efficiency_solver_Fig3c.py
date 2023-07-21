# This script was used to analyze the data shown in Figure 3c.

# It doesn't read in any files. It is parameterized with experimental values under the USER INPUTS section.
# The extract_insertion_sequences_and_report_indel_sizes_and_signatures.py script was used to calculate fraction of reads with each number of insertions.

# The eta1 (A->B editing efficiency %/day) and eta2 (B->A editing efficiency %/day) are outputted as print statements.

import numpy as np
from scipy.optimize import root

##### USER INPUTS ####

eta_guess = [0.09, 0.01] #provide initial guesses for eta1 and eta2 (% editing / day)
Xviable0 = 0.824 #provide the experimentally-determined value for the fraction of loci that are viable for editing, but wildtype at the beginning of the experiment
t = 23 #how many days of editing have occurred at the timepoint of interest?
X1 = 0.605 #what fraction of reads have 1 edit at the timepoint of interest?
X2 = 0.054 #what fraction of reads have 2 edits at the timepoint of interest?

##### END USER INPUTS #####

##### DEFINE FUNCTIONS #####

def f(eta):
    # This stores the equations for eta1 and eta2, rearranged to their root form 
    eta1 = eta[0]
    eta2 = eta[1]
    
    X1_numerator = -(np.exp(-eta1*t-eta2*t)*(-np.exp(eta1*t)+np.exp(eta2*t))*eta1*Xviable0) #this equation is explained in the main text
    X1_denominator = eta1 - eta2
    
    X2_numerator = -(np.exp(-eta1*t-eta2*t)*eta1*eta2*(-np.exp(eta1*t)+np.exp(eta2*t)+np.exp(eta2*t)*eta1*t-np.exp(eta2*t)*eta2*t)*Xviable0) #this equation is explained in the main text
    X2_denominator = (eta1 - eta2)**2
    
    return [X1_numerator / X1_denominator - X1, X2_numerator / X2_denominator - X2]


# Use scipy's root function to solve the system of equations
result = root(f, eta_guess)

# Extract the values of eta1 and eta2 from the result
eta1_solution = result.x[0]
eta2_solution = result.x[1]

print("Solution:")
print("eta1 =", eta1_solution)
print("eta2 =", eta2_solution)