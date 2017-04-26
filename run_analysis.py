import sys
import numpy as np
from pandas import *
from parse_pyvcf import get_matrices
from pprint import pprint

def run_softimpute(M):
    tol = 1e-10
    n,p=M.shape
    X = DataFrame(np.zeros((n,p)))
    M = DataFrame(M)
    for i in np.arange(1.0,0.01,-0.01):
        sftdiff = np.where(M.isnull(), X, M)

        U,S,V = np.linalg.svd(sftdiff,full_matrices=False)
        S = np.diag(S)
        Slmb = np.maximum(S-i,np.zeros(S.shape))
        X_new = np.dot(U, np.dot(S, V))
        diff = np.linalg.norm(X_new-X,'fro')**2/np.linalg.norm(X_new,'fro')**2 
        print "i:", i, "diff:", diff
        if diff < tol:
            return X_new
        X = X_new

def run_lmm(X,Y,U,G):
    
    Ulmm = X
    Umeans = Ulmm.mean(axis=0)
    Ulmm = Ulmm - Umeans
    Glmm = np.cov(Ulmm.transpose())
    Xlmm = U
    Xmeans = Xlmm.mean(axis=0)
    Xlmm = Xlmm - Xmeans
    Y = Y.values[:,[0]] # eye color
    Y = Y - Y.mean(axis=0)
    rho = 0.1
    V = np.dot(Ulmm, np.dot(Glmm,Ulmm.transpose())) + rho*np.cov(Ulmm)
 

    beta_inv = np.dot(Xlmm.transpose(),np.dot(np.linalg.inv(V), Xlmm))
    beta = np.linalg.inv(beta_inv)
    beta = np.dot(beta, Xlmm.transpose())
    beta = np.dot(beta, np.dot(np.linalg.inv(V),Y))
    pprint(beta)

    gamma = Glmm.dot(Ulmm.transpose()).dot(np.linalg.inv(V)).dot(Y-Xlmm.dot(beta))
    pprint(gamma)

    print "Ulmm:", np.max(Ulmm), np.min(Ulmm)
    print "Glmm:", np.max(Glmm), np.min(Glmm)
    print "Xlmm:", np.max(Xlmm), np.min(Xlmm)

def main(args):
    X,Y,U,G = get_matrices('Height')
    X_new = DataFrame(run_softimpute(X))


    pprint(X[:10][:10])
    pprint(X_new[:10][:10])
    pprint(Y)
    
    run_lmm(X_new, Y, U, G)

if __name__ == '__main__':
    main(sys.argv[1:])
