import transformations as sp
import helpers as h
import scanpy as sc
from scipy.sparse import csr_matrix

if __name__ == "__main__":
    trans = sc.read_h5ad("tests/pseudo_shifted_log.h5ad")

    theirs = sc.read_h5ad("tests/pseudo_shifted_log_tgp.h5ad")
    
    print(csr_matrix(trans.X))
    print("")
    print(csr_matrix(theirs.X.T))
    print("")
    print(csr_matrix(trans.X / theirs.X.T))
    print("")
    print((trans.X - theirs.X.T).sum())
