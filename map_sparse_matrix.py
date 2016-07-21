from scipy.sparse import coo_matrix
import numpy as np
import matplotlib.pyplot as plt

def file_plot_data(filename):
    with open(filename, 'r') as f:
        row = []
        col = []
        data = []
        num_points = f.readline()
        for line in f:
            count = 0
            for word in line.split():
                if (count == 0):
                    row.append(int(word))
                elif (count == 1):
                    col.append(int(word))
                else:
                    data.append(float(word))
                count += 1

        return (num_points, row, col, data)

data = file_plot_data("unordered_bound_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0],data[0]))

plt.spy(mat, precision=0.0001)
plt.title("Unordered Bound Matrix")
plt.show()
plt.close()

data = file_plot_data("ordered_bound_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0],data[0]))

plt.spy(mat, precision=0.0001)
plt.title("Ordered Bound Matrix")
plt.show()
plt.close()

data = file_plot_data("unordered_not_bound_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0], data[0]))

plt.spy(mat, precision=0.0001)
plt.title("Unordered Not Bound Matrix")
plt.show()
plt.close()

data = file_plot_data("ordered_not_bound_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0], data[0]))

plt.spy(mat, precision=0.0001)
plt.title("Ordered Not Bound Matrix")
plt.show()
plt.close()

print("Figures Printed!!")
