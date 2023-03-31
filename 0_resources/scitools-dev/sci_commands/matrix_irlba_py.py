# tfidf in sparse format to umap

import argparse, re, sys, gzip, io
import scipy as scipy
import numpy as np
import pandas as pd
import irlb
from scipy import sparse
from scipy.sparse import dok_matrix

def main():
    p = argparse.ArgumentParser(description='IRLB on TFIDF in sparse format')
    p.add_argument('input', action='store', type=str,  
        help='input name')
    p.add_argument('cells', action='store', type=str,  
        help='cell names, column file of sparse matrix')
    p.add_argument('output', action='store', type=str,  
        help='output name')
    p.add_argument('-d', '--dim', action='store', type=int, default=50,
        help='dimensions (default 50)')

    args = p.parse_args()

    A = np.loadtxt(args.input)
    colnames = np.loadtxt(args.cells,dtype=str)
    rownames = ["DIM"+str(x) for x in range(1,args.dim+1)]

    C = sparse.coo_matrix((A[:,2], (A[:,0]-1, A[:,1]-1)), shape=(int(max(A[:,0])), int(max(A[:,1]))))
    del A
    S = irlb.irlb(C, args.dim)
    M = np.transpose(S[2])

    df = pd.DataFrame(M, index=rownames, columns=colnames)
    df.to_csv(args.output, index=True, header=True, sep='\t', index_label=False)

# end main function

if __name__=='__main__':
    main()
