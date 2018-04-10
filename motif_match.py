"""
match PWM to a sequence
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

######################################################################
# functions

# Extracts the sequence in a FASTA file.
def load_FASTA (filename):
    """
    This function takes a filename (str) as input and returns
    the sequence (str) found in that file. 
    """
    # open file according to the filename given.
    with open (filename) as f:
       # reads what's in the file as a string, and then splits
       # the string by '>' into a list.
       lines = f.read().split('>')[1:]
       # close the file.
       f.close()
       # print what's after the first '\n'
       return ''.join(lines[0].split('\n')[1:])

def load_pwm (filename):
	with open (filename) as f:
		lines=f.read().split('\n')
		f.close()
		return lines


def pwm_similarity_score (motif_pwm_filename):
	# read in the pwm file, each position is the probability
	pwm = pd.read_csv(motif_pwm_filename,sep='\s+',header=None,engine='python')
	pwm_len=len(pwm)
	# insert the 5th column for N
	pwm.insert(4,4,0)
	# rename the columns
	pwm.columns = ['A', 'C', 'G', 'T', 'N']
	# create a dictionary for 5 letters
	pwm_order={'A':0,'C':1,'G':2,'T':3,'N':4}
	# return the matrix, length and dictionary
	return pwm,pwm_len,pwm_order



# load in hg19 or hg38 and extract only 24 chromosomes
def load_genome (genome_filename):
	# read in the genome.fa file using SeqIO
	genome=list(SeqIO.parse(genome_filename, 'fasta'))
	# only take chr1~22 + chrX + chrY
	chrs=[]
	for i in range(22):
		chrs.append('chr'+str(i+1))
	chrs.append('chrX')
	chrs.append('chrY')
	ids=[]
	for i in range(len(genome)):
		ids.append(genome[i].id)
	# put the 24 chrs into chrs_clean
	chrs_clean=[]
	for i in range(len(chrs)):
		chr=chrs[i]
		id_index=ids.index(chr)
		chrs_clean.append(genome[id_index])
	return chrs_clean

"""
# calcualte the maximum score for a given 200bp sequence
def max_similarity_score (sequence, pwm_len, pwm, pwm_order, pwm_max_sum):
	inp=sequence
	inp=inp.upper()
	max_ratio=0
	out=[0]*(200-pwm_len)
	for i in range(len(inp)-pwm_len):
		seq=inp[i:(i+pwm_len)]
		addup=0
		for k in range(pwm_len):
			addup+=pwm.loc[k,pwm_order[str(seq[k])]]
		ratio=addup/pwm_max_sum
		if ratio>max_ratio:
			max_ratio=ratio
	return max_ratio

"""
# calcualte the maximum score for a given sequence
# each sliding window, minus the first nucleotide, plus the last nucleotide
def max_similarity_score (sequence, pwm_len, pwm, pwm_order, pwm_max_sum, pwm_dic):
	# add this line if input is not in the Seq format
	# inp=Seq(sequence, generic_dna)
	# change the letters to the upper
	inp=sequence.upper()
	# the add_max starts with 0
	addup_max=0
	# calculate each score in a pwm_len window, step is 1bp
	for k in range(len(inp)-pwm_len):
		seq=inp[k:(k+pwm_len)]
		# find each bp probablity using dictionary
		keys=[str(m)+n for m,n in zip(range(pwm_len),[i for i in seq])]
		fe=[pwm_dic.get(key) for key in keys]
		# sum up the probability of each position
		addup=sum(fe)
		# replace addup_max if a bigger score is found
		if addup > addup_max:
			addup_max=addup
	## calculate scores of the reverse complement
	inp_rc=inp.reverse_complement()
	# calculate each score in a pwm_len window, step is 1bp
	for k in range(len(inp_rc)-pwm_len):
		seq=inp_rc[k:(k+pwm_len)]
		# find each bp probablity using dictionary
		keys=[str(m)+n for m,n in zip(range(pwm_len),[i for i in seq])]
		fe=[pwm_dic.get(key) for key in keys]
		# sum up the probability of each position
		addup=sum(fe)
		# replace addup_max if a bigger score is found
		if addup > addup_max:
			addup_max=addup
	# return the max ratio
	return addup_max/pwm_max_sum

# calculate the max score for each interval in a train_file
def max_score (chrs_clean, train_filename, pwm_len, pwm, pwm_order,pwm_max_sum, pwm_dic):
	# write every 10000 lines
	nlines = 51676736
	l=nlines/10000+1
	#nlines=20012
	#l=3
	nl=1
	while nl <= l:
		print nl
		ind_start=(nl-1)*10000+1
		ind_end=nl*10000
		out=[0]*10000
		if nl == l:
			ind_start=(nl-1)*10000+1
			ind_end=nlines
			out=[0]*(nlines-10000*(l-1))
		with open(train_filename) as fp:
			for i, line in enumerate(fp):
				if i >= (ind_start-1) and i < ind_end:
					chr=line.split('\t')[0]
					start=int(line.split('\t')[1])
					end=int(line.split('\t')[2])
					index=[m for m in range(24) if chrs_clean[m].id==chr][0]
					sequence=chrs_clean[index].seq[start:end]
					newi=i-10000*(nl-1)
					#out[newi-1]=max_similarity_score(sequence, pwm_len, pwm, pwm_order,pwm_max_sum,pwm_dic)
					out[newi]=str(chr)+'\t'+str(start)+'\t'+str(end)+'\t'+ str(max_similarity_score(sequence, pwm_len, pwm, pwm_order,pwm_max_sum,pwm_dic))
				elif i > ind_end:
					break
		fp.close()
		f=open('motif_scores.txt','a+')
		f.write('\n'.join(str(a) for a in out))
		f.write('\n')
		f.close()
		nl += 1

# use this for multi-tasking
# this calculate max_scores for each interval given by a file
# calculate the max score for each interval in a train_file
# NOTE!!! no header line in the train_file
def max_score_multi (chrs_clean, train_filename, number_lines, pwm_len, pwm, pwm_order,pwm_max_sum, pwm_dic):
	# read in filename and number of lines
	nlines=int(number_lines)
	# define an interval, run the data by each interval
	interval=10000
	# l is the number of iterations
	l=nlines/interval+1
	# nl starts at 1
	nl=1
	# loop when nl <= l
	while nl <= l:
		print nl
		ind_start=(nl-1)*interval+1
		ind_end=nl*interval
		out=[0]*interval
		# empty out list to store the scores
		if nl == l:
			ind_start=(nl-1)*interval+1
			ind_end=nlines
			if (nlines-interval*(l-1)) >= 1:
				out=[0]*(nlines-interval*(l-1))
			else:
				break
		# if start > end, break the loop
		if ind_start < ind_end:
		# read in interval lines from train_filename, calculate the maximum score of each interval
			with open(train_filename) as fp:
				for i, line in enumerate(fp):
					if i >= (ind_start-1) and i < ind_end:
						chr=line.split('\t')[0]
						start=int(line.split('\t')[1])
						end=int(line.split('\t')[2])
						index=[m for m in range(24) if chrs_clean[m].id==chr][0]
						sequence=chrs_clean[index].seq[start:end]
						newi=i-interval*(nl-1)
						#out[newi-1]=max_similarity_score(sequence, pwm_len, pwm, pwm_order,pwm_max_sum,pwm_dic)
						out[newi]=str(chr)+'\t'+str(start)+'\t'+str(end)+'\t'+ str(max_similarity_score(sequence, pwm_len, pwm, pwm_order,pwm_max_sum,pwm_dic))
						#print i
					elif i > ind_end:
						break
			fp.close()
		# write out the 100000 scores
			f=open('score_'+train_filename,'a+')
			f.write('\n'.join(str(a) for a in out))
			f.write('\n')
			f.close()
		nl += 1

"""
ind_start=10001
ind_end=20000
out=[0]*10000
with open("MYC.train.labels.tsv") as fp:
    for i, line in enumerate(fp):
        if i >= ind_start and i <= ind_end:
            chr=line.split('\t')[0]
            start=int(line.split('\t')[1])
            end=int(line.split('\t')[2])
            index=[m for m in range(24) if chrs_clean[m].id==chr][0]
            sequence=chrs_clean[index].seq[start:end]
            out[i-1]=max_similarity_score(sequence, pwm_len, pwm, pwm_order,pwm_max_sum,pwm_dic)
        elif i > ind_end:
            break
fp.close()


# calculate the maximum score for a given region file
def max_score (chrs_clean, train_filename, pwm_len, pwm, pwm_order,pwm_max_sum, pwm_dic):
	# train_data = pd.read_csv(train_filename,sep='\t',header=0,engine='python')
	# write every 1000 lines
	# nlines = 51676736
	# n=nlines/10000
	ind_start=1
	ind_end=1000
	out=[0]*1000
	with open(train_filename) as fp:
    	for i, line in enumerate(fp):
        	if i >= ind_start and i < ind_end:
            	chr=line.split('\t')[0]
            	start=line.split('\t')[1]
            	end=line.split('\t')[2]
            	index=[m for m in range(24) if chrs_clean[m].id==chr][0]
            	sequence=chrs_clean[index].seq[start:end]
            	out[i]=max_similarity_score(sequence, pwm_len, pwm, pwm_order,pwm_max_sum,pwm_dic)
        	elif i > ind_end:
            	break
	fp.close()
	return out
"""

"""
# create input sequences
for n in range(len(chrs_clean)):
	inp=chrs_clean[n].seq
	inp=inp.upper()
	inp_rc=inp.reverse_complement()
	# calculate similarity for each position
	f=open(chrs_clean[n].id+'.txt','a+')
	for i in range(len(inp)-pwm_len):
		seq=inp[i:(i+pwm_len)]
		addup=0
		for k in range(pwm_len):
			addup+=pwm.loc[k,pwm_order[str(seq[k])]]
			ratio=addup/pwm_max_sum
			pos=i+pwm_len/2
			out=str(pos)+'\t'+str(ratio)+'\n'
			f.write(out)
	f.close()



input_sequence='gggaggaggagccaagatggccgaataNNNNNNNNNNNNNNNNNctacagctcccagcgtgagcgacgcagaagacgggtgatttctgcatttccatctgaggtaccgggttcatctcactagggagtgccagacagtgggcgcaggccagtgtgtgtgcgcaccgtgcgcgagccgaagcagggcgaggcattgcctcacctgggaagcgcaaggggtcagggagttccctttccgagtcaaagaaaggggtgacggacgcacctggaaaatcgggtcactcccacccgaatattgcgcttttcagaccggcttaagaaacggcgcaccacgagactatatcccacacctggctcagagggtcctacgcccacggaatctcgctgattgctagcacagcagtctgagatcaaactgcaaggcggcaacgaggctgggggaggggcgcccgccattgcccaggcttgcttaggtaaacaaagcagccgggaagctcgaactgggtggagcccaccacagctcaaggaggcctgcctgcctctgtaggctccacctctgggggcagggcacagacaaacaaaaagacagcagtaacctctgcagacttaagtgtccctgtctgacagctttgaagagagcagtggttctcccagcacgcagctggagatctgagaacgggcagactgcctcctcaagtgggtccctgacccctgacccccgagcagcctaactgggaggcaccccccagcaggggcacactgacacctcacacggcagggtattccaacagacctgcagctgagggtcctgtctgttagaaggaaaactaacaaccagaaaggacatctacaccgaaaacccatctgtacatcaccatcatcaaagaccaaaagtagataaaaccacaaag'
inp=Seq(input_sequence, generic_dna)
inp=inp.upper()
inp_rc=inp.reverse_complement()


# sense strand scanning of the pwm, outputs a matrix with positions and similarity scores
out=''
for i in range(len(inp)-pwm_len):
	seq=inp[i:(i+pwm_len)]
	addup=0
	for k in range(pwm_len):
		addup+=pwm.loc[k,pwm_order[str(seq[k])]]
	ratio=addup/pwm_max_sum
	pos=i+pwm_len/2
	out+=str(pos)+'\t'+str(ratio)+'\n'




out=''
for i in range(len(inp)-pwm_len):
	seq=inp[i:(i+pwm_len)]
	addup=0
	for k in range(pwm_len):
		addup+=pwm.loc[k,pwm_order[str(seq[k])]]
	ratio=addup/pwm_max_sum
	pos=i+pwm_len/2
	out+=str(pos)+'\t'+str(ratio)+'\n'

	
f=open('test.txt','w')
f.write(out)
f.close()
"""












