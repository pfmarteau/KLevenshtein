import numpy as np
import KLevenshtein
import Levenshtein

def gram(X, sigma, epsilon):
    G=np.zeros((len(X),len(X)))
    for i in range(len(X)):
        for j in range(i, len(X)):
            G[i,j]=KLevenshtein.similarity_int(X[i],X[j],sigma, epsilon)
            G[j,i]=G[i,j]
    return G

def checkPositiveDefiniteness():
    # check definiteness
    from numpy import linalg as LA
    sigma=5; epsilon = 1e-3
    N=100; L=100
    X=np.random.randint(26,size=(N,L))
    G=gram(X, sigma, epsilon)
    w, v = LA.eig(G)
    print('Eigen spectrum')
    print(w)

if __name__ == '__main__':
    import Levenshtein
    
    A = "the little mairmaid"
    B = "the little big man"
    C = "the big little man"
    D = "an old sailboat at sea"
    sigma=5
    epsilon=1
    L=[A,B,C,D]
    GLev=np.zeros((len(L),len(L)))
    GKLev=np.zeros((len(L),len(L)))
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            GLev[i,j]=Levenshtein.distance(L[i],L[j])
            GKLev[i,j]=KLevenshtein.similarity_str(L[i],L[j],sigma,epsilon)
            print("[%s] v.s. [%s] | lev= %d | klev= %1.4f" % (L[i], L[j], GLev[i,j], GKLev[i,j]))

    checkPositiveDefiniteness()


