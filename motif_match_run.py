"""
match PWM to hg19
"""

################################################################################
# Modules

# import regular expression module
import re

# import sys module
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd

# import modules from motif_match.py
from motif_match import pwm_similarity_score
from motif_match import load_genome
from motif_match import max_similarity_score
from motif_match import max_score
from motif_match import max_score_multi

#########################################################################
# calculate the maximum of summed weights
pwm=pwm_similarity_score('MYC_matrix.txt')[0]
pwm_len=pwm_similarity_score('MYC_matrix.txt')[1]
pwm_order=pwm_similarity_score('MYC_matrix.txt')[2]
pwm_max_sum=0
for i in range(pwm_len):
	pwm_max_sum+=(max(pwm.loc[i,:]))

# create a dictionary for pwm
pwm_dic={}
for i in range(pwm_len):
	pwm_dic.update({(str(i)+'A'): pwm.loc[i,'A']})
	pwm_dic.update({str(i)+'C': pwm.loc[i,'C']})
	pwm_dic.update({str(i)+'G': pwm.loc[i,'G']})
	pwm_dic.update({str(i)+'T': pwm.loc[i,'T']})
	pwm_dic.update({str(i)+'N': pwm.loc[i,'N']})

chrs_clean=load_genome('hg19.genome.fa')

# Function main
def main():
    """
    This function sets the whole script will execute if run at the commandline.
    """
    # assume the input filename is the first argument
    train_filename = sys.argv[1]
    if len(sys.argv) == 2:
    	max_score(chrs_clean, train_filename, pwm_len, pwm, pwm_order, pwm_max_sum, pwm_dic)
    # assume the output filename is the second argument
    if len(sys.argv) == 3:
        number_lines = sys.argv[2]
        max_score_multi(chrs_clean, train_filename, number_lines, pwm_len, pwm, pwm_order, pwm_max_sum, pwm_dic)

# this will be executed when the script is run    
if __name__=='__main__':
    main()


"""
train_filename='MYC.train.labels.tsv'
import timeit
start = timeit.default_timer()
max_score(chrs_clean, train_filename, pwm_len, pwm, pwm_order, pwm_max_sum, pwm_dic)
stop = timeit.default_timer()
print stop - start
"""

