
from __future__ import division, print_function

from time import time

from Bio import Phylo as _Phylo

from BioExt.phylo import Phylo
from six.moves import cStringIO

import numpy as np

buf = cStringIO("((((((variant_2_57_copies:1e-10,variant_7_339_copies:1e-10):0.00135078,variant_8_118_copies:1e-10):0.0015468,variant_4_56_copies:1e-10):2.07765e-05,variant_5_233_copies:1e-10):4.05015e-05,variant_1_136_copies:1e-10):0.00173927,variant_3_67_copies:0.00916762,variant_6_418_copies:1e-10);")

# buf = cStringIO("(A,B,(C,D)E)F;")

N = 5000

def distance_matrix(tree):
    leaves = tree.get_terminals()
    dmat = np.zeros((len(leaves), len(leaves)), dtype=float)
    for i in range(len(leaves)):
        for j in range(len(leaves)):
            if j > i:
                dmat[i, j] = tree.distance(leaves[i], leaves[j])
    return dmat

def benchmark_read(cls):
    a = time()
    for _ in range(N):
        buf.seek(0)
        tree = cls.read(buf, 'newick')
    b = time()
    print('%d iterations took %g seconds' % (N, b - a))

def benchmark_dmat(cls):
    buf.seek(0)
    tree = cls.read(buf, 'newick')
    a = time()
    for _ in range(N):
        dmat = distance_matrix(tree)
    b = time()
    print('%d iterations took %g seconds' % (N, b - a))

def print_stuff(cls):
    buf.seek(0)
    tree = cls.read(buf, 'newick')
    print(tree)
    print(tree.get_terminals())
    print(distance_matrix(tree))

benchmark_read(_Phylo)
benchmark_dmat(_Phylo)
benchmark_read(Phylo)
benchmark_dmat(Phylo)

print_stuff(_Phylo)
print_stuff(Phylo)
