import numpy as np
import matplotlib.pyplot as plt

def get_data(filename):
    with open(filename, 'r') as f:
        x = []
        y = []
        lines = f.readlines()
        num_points = len(lines)
        for line in lines:
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))
        return (x,y)

data = get_data("Energy_plot.txt")

d = np.asarray(data[0])
energy = np.asarray(data[1])

plt.plot(d,energy)
plt.show()
