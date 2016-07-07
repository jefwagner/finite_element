from scipy.sparse import coo_matrix
import numpy as np
# import matplotlib.pyplot as plt

row = []
col = []
data = []

# def coo_matrix_fill(filename):
#     with open(filename) as f:
#         i = 0
#         for line in f:
#             row[i] = f.read(1)
#             f.read(1)
#             col[i] = f.reat(1)
#             f.read(1)
#             data[i] = f.read()
#             i += 1
#
def sparse_matrix_fill(filename):
    with open(filename, 'r') as f:
        i = 0
        for line in f:
            row[i] = f.read(1)
            f.read(1)
            col[i] = f.read(1)
            f.read(1)
            data[i] = f.read()




            print i

            i+=1

sparse_matrix_fill("unordered_mass_mat.txt")
# coo_matrix_fill("unordered_mass_mat.txt")
