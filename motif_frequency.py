"""
Motif Frequency Plots of mdscan results (.bed files)

"""

###############################################################################
#Modules

import numpy as np
import pylab as P
import optparse

###############################################################################
# Functions

def motif_freq_plot(motif_filename, shift_size):
    f = open(motif_filename, 'r')
    lines = f.readlines()
    f.close()
    line =[i.strip().split('\t') for i in lines]
    plot_matrix = []
    for i in line:
        plot_matrix.append([i[4], i[5]])
    for i in plot_matrix:
        i[0]=int(i[0])-shift_size
    vector_positions = [int(i[0]) for i in plot_matrix]
    P.figure()
    n,bins,patches=P.hist(vector_positions, 25, normed=True, histtype='bar')
    P.show()
    
    
