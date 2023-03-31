# imports
import argparse
import os
import scrublet as scr
import scipy as scipy
import numpy as np
import pandas as pd 
from scipy import sparse
from scipy import io
from matplotlib import pyplot as plt

def main():

   p = argparse.ArgumentParser(description='Calls scrublet on a sparse count matrix')

   p.add_argument('input', action='store', type=str,  
      help='sparse matrix values (from a count matrix) or dense matrix (from a cistopic matrix)')
   p.add_argument('cells', action='store', type=str,  
      help='cell names, column file of sparse matrix')
   p.add_argument('output', action='store', type=str,
      help='output prefix')
   p.add_argument('--expected_doublet_rate', action='store', type=float, default=0.06,
      help='estimated doublet rate, default 0.06')
   p.add_argument('--min_counts', action='store', type=int, default=2,
      help='minimum number of counts for scrublet before PCA, default 2')
   p.add_argument('--min_cells', action='store', type=int, default=3,
      help='minimum number of cells for scrublet before PCA, default 3')
   p.add_argument('--min_gene_variability_pctl', action='store', type=int, default=85,
      help='keeps only highest variabile genes before PCA, default 85')
   p.add_argument('--n_prin_comps', action='store', type=int, default=30,
      help='number of PCA comps to use before graph construction, default 30')

   args = p.parse_args()

   # load sparse matrix and transform to csc
   A = np.loadtxt(args.input)
   C = np.loadtxt(args.cells,dtype='str')
   counts_matrix = sparse.csc_matrix((A[:,2], (A[:,1]-1, A[:,0]-1)), shape=(int(max(A[:,1])), int(max(A[:,0])))) # subtracting 1 for 0 indexing in python
   del A

   # scrub
   scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate)
   doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts, min_cells=args.min_cells, min_gene_variability_pctl=args.min_gene_variability_pctl, n_prin_comps=args.n_prin_comps)

   # save
   outname1=args.output + ".doublet_distribution.png"
   outname2=args.output + ".doublet_scores.values"
   outname3=args.output + ".doublet_prediction.values"

   scrub.plot_histogram()
   plt.savefig(outname1)

   out1 = pd.DataFrame(np.vstack((C, doublet_scores))).transpose()
   out1.to_csv(outname2,sep='\t',header=None,index=None)

   out2 = pd.DataFrame(np.vstack((C, predicted_doublets))).transpose()
   out2.to_csv(outname3,sep='\t',header=None,index=None)

if __name__=='__main__':
    main()
