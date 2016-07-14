from scipy.sparse import coo_matrix
import numpy as np
import matplotlib.pyplot as plt

def fileplot(filename):
    with open(filename, 'r') as f:
        row = []
        col = []
        data = []
        for line in f:
            count = 0
            for word in line.split(','):
                if (count == 0):
                    row.append(int(word))
                elif (count == 1):
                    col.append(int(word))
                else:
                    data.append(float(word))
                count += 1

        return (row, col, data)

data = fileplot("unordered_mass_mat.txt")

print data
