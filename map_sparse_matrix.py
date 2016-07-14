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
            for word in line.split(','):
                if (count == 0):
                    row.append(int(word))
                elif (count == 1):
                    col.append(int(word))
                else:
                    data.append(float(word))
                count += 1

        return (num_points, row, col, data)

data = file_plot_data("unordered_mass_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0],data[0]))

plt.spy(mat)
plt.title("Unordered Mass Matrix")
plt.show()
plt.close()

data = file_plot_data("ordered_mass_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0],data[0]))

plt.spy(mat)
plt.title("Ordered Mass Matrix")
plt.show()
plt.close()

data = file_plot_data("unordered_stiff_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0], data[0]))

plt.spy(mat)
plt.title("Unordered Stiffness Matrix")
plt.show()
plt.close()

data = file_plot_data("ordered_stiff_mat.txt")

mat = coo_matrix((data[3], (data[1], data[2])), shape=(data[0], data[0]))

plt.spy(mat)
plt.title("Ordered Stiffness Matrix")
plt.show()
plt.close()

print("Figures Printed!!")
