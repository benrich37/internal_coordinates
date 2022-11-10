import numpy as np
import copy

def stretch_vec(vec):
    # assumes an n by 3 vec, as usually used to describe molecule geometry
    n = len(vec)
    output = list(np.zeros(3*n))
    for i in range(n):
        for j in range(3):
            output[3*i + j] = vec[i][j]
    return np.array(output)

def squish_vec(vec):
    # assumed a 1 by 3N vec
    n = len(vec)
    output = []
    for i in range(int(n/3)):
        output.append([0, 0, 0])
    for i in range(int(n/3)):
        for j in range(3):
            output[i][j] = vec[3*i+j]
    return np.array(output)

def stretch_vecs(vecs):
    output = []
    for v in vecs:
        output.append(stretch_vec(v))
    return output

def squish_vecs(vecs):
    output = []
    for v in vecs:
        output.append(squish_vec(v))
    return output

def normalize_basis(vecs):
    output = []
    for v in vecs:
        output.append(v/np.linalg.norm(v))
    return output

def get_vec_in_basis(vec, basis):
    cs = []
    for b in basis:
        cs.append(np.dot(vec, b))
    return cs

def get_vec_from_cs(cs, basis):
    vec = np.zeros(np.shape(basis[0]))
    for i in range(len(cs)):
        vec += basis[i]*cs[i]
    return vec

def get_cart_mom_n(posns, n):
    mom = np.zeros(np.shape(posns[0]))
    l = len(posns)
    for p in posns:
        mom += p*(np.linalg.norm(p)**n)
    mom = mom/l
    return mom

# def align_mom_n

def align_to_moments(posns, mom_n0, mom_n1, mom_n2, n0=0, n1=1, n2=2):
    cur_mom_n0 = get_cart_mom_n(posns, n0)