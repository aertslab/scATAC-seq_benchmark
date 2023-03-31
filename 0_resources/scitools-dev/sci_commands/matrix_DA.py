import argparse, re, sys, gzip, io
import scipy as scipy
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import dok_matrix
from scipy.stats import ranksums

def main():
   
   p = argparse.ArgumentParser(description='differential accessibility test on sparse matrix')

   p.add_argument('input', action='store', type=str,  
      help='sparse matrix values')
   p.add_argument('peaks', action='store', type=str,
      help='peak names, row file of sparse matrix')
   p.add_argument('cells', action='store', type=str,  
      help='cell names, column file of sparse matrix')
   p.add_argument('meta', action='store', type=str,
      help='metadata file, two column value file')
   p.add_argument('comparisons', action='store', type=str,
      help='comma separated list of strings to run pairwise comparisons, i.e. pg_1,pg_2,pg_3')
   p.add_argument('range1', action='store', type=int,
      help='0-based starting range of rows to look at')
   p.add_argument('range2', action='store', type=int,
      help='0-based starting range of rows to end at')

   p.add_argument('output', action='store', type=str,  
      help='output prefix')
   p.add_argument('--test', action='store', type=str, choices=["wilcox"], default="wilcox", 
      help='type of comparison: wilcox')
   p.add_argument('--fdr', action='store', type=str, choices=["BH"], default="BH",
      help='FDR test to perform, not implemented currently')

   args = p.parse_args()

   # read metadata
   meta = pd.read_csv(args.meta,sep="\t",header=None)
   meta = meta.rename(index=meta[0])

   # read string of pairwise comparisons
   comps = args.comparisons.split(",")

   # check existence of comparisons in metadata
   for comp in comps:
      if (meta[meta[1]==comp].shape[0] > 0):
         print("Comparison feature " + comp + " has " + str(meta[meta[1]==comp].shape[0]) + " elements found in metadata.")
      else:
         print("Comparison feature " + comp + " is not found in metadata. Exiting.")
         exit()

   # read sparse matrix into coo format
   A = np.loadtxt(args.input)
   rownames = np.loadtxt(args.peaks,dtype=object)
   colnames = np.loadtxt(args.cells,dtype=str)
   C = sparse.coo_matrix((A[:,2], (A[:,0]-1, A[:,1]-1)), shape=(int(max(A[:,0])), colnames.shape[0]))
   del A

   # reorder and potentially subset metadata by matrix colnames
   meta2 = meta.loc[colnames]
   print(str(meta.shape[0]) + " cells in metadata. " + str(colnames.shape[0]) + " cells in matrix. " + str(meta2.shape[0]) + " cells found in metadata. Continuing with elements found.")

   for comp in comps:
      if (meta2[meta2[1]==comp].shape[0] > 0):
         print("Comparison feature " + comp + " has " + str(meta2[meta2[1]==comp].shape[0]) + " elements found in matrix.")
      else:
         print("Comparison feature " + comp + " is not found in matrix. Exiting.")
         exit()

   # enumerate through comparsions
   for num1,comp1 in enumerate(comps):
      for num2,comp2 in enumerate(comps):
         if (num2 > num1):
            # create file
            filename = args.output + "." + comp1 + "_" + comp2 + "." + str(args.range1) + "_" + str(args.range2) + ".txt"
            with open(filename, 'a') as OUTFILE:
               # for each row
               for i in range(args.range1,args.range2):
                  name=rownames[i]
                  sub = C.tocsr()[i,:].todense()
                  a = sub[0,meta2.index.get_indexer_for((meta2[meta2[1] == comp1].index))]
                  b = sub[0,meta2.index.get_indexer_for((meta2[meta2[1] == comp2].index))]
                  meana = np.mean(np.squeeze(np.asarray(a)))
                  meanb = np.mean(np.squeeze(np.asarray(b)))
                  if (args.test == "wilcox"):
                     test = ranksums(np.squeeze(np.asarray(a)),np.squeeze(np.asarray(b)))
                  outstr = name + "\t" + str(meana) + "\t" + str(meanb) + "\t" + str(test[0]) + "\t" + str(test[1]) + "\n"
                  _ = OUTFILE.write(outstr)

if __name__=='__main__':
    main()
